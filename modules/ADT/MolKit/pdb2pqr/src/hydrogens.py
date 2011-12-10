import os
import string
import math

from definitions import *
from utilities import *
from quatfit import *
from routines import *

__date__ = "16 May 2006"
__author__ = "Todd Dolinsky, Jens Erik Nielsen"

HDEBUG = 0
HYDPATH = "dat/HYDROGENS.xml"
ANGLE_CUTOFF = 20.0       # A - D - H(D) angle
DIST_CUTOFF = 3.3         # H(D) to A distance

class HydrogenHandler(sax.ContentHandler):
    """
        Extends the SAX XML Parser to parse the Hydrogens.xml
        class
    """
    def __init__(self):
        """
            Initalize the class.
        """
        
        self.curelement = ""
        self.curatom = None
        self.curobj = None
        self.curholder = None
        self.map = {}

    def startElement(self, name, attributes):
        """
            Create optimization holder objects or atoms
        """
        if name == "class":
            obj = OptimizationHolder()
            self.curholder = obj
            self.curobj = obj
        elif name == "atom":
            obj = DefinitionAtom()
            self.curatom = obj
            self.curobj = obj
        else:
            self.curelement = name
        return

    def endElement(self, name):
        """
            Complete whatever object is currently passed in
            by the name parameter
        """
        if name == "class": # Complete Residue object
            obj = self.curholder
            if not isinstance(obj, OptimizationHolder):
                raise ValueError, "Internal error parsing XML!"
          
            self.map[obj.name] = obj
            self.curholder = None
            self.curobj = None

        
        elif name == "atom": # Complete atom object
            atom = self.curatom
            if not isinstance(atom, DefinitionAtom):
                raise ValueError, "Internal error parsing XML!"
            atomname = atom.name
            if atomname == "":
                raise ValueError, "Atom name not set in XML!"
            else:
                self.curholder.map[atomname] = atom
                self.curatom = None
                self.curobj = self.curholder

        else: # Just free the current element namespace
            self.curelement = ""

        return self.map

    def characters(self, text):
        """
            Set a given attribute of the object to the text
        """
        
        if text.isspace(): return

        # If this is a float, make it so
        try:
            value = float(str(text))
        except ValueError:
            value = str(text)

        setattr(self.curobj, self.curelement, value)
      

class PotentialBond:
    """
        A small class containing the hbond structure
    """
    def __init__(self, atom1, atom2, dist):
        """
            Initialize the class

            Parameters
                atom1:  The first atom in the potential bond (Atom)
                atom2:  The second atom in the potential bond (Atom)
                dist:  The distance between the two atoms (float)
        """
        self.atom1 = atom1
        self.atom2 = atom2
        self.dist = dist

    def __str__(self):
        """
            String for debugging
        """
        txt = "%s %s" % (self.atom1.name, self.atom1.residue)
        txt += " to "
        txt += "%s %s" % (self.atom2.name, self.atom2.residue)
        txt += " (%.2f A)" % self.dist
        return txt
     

class Optimize:
    """
        The holder class for the hydrogen optimization
        routines. Individual optimization types inherit off of this
        class.  Any functions used by multiple types appear here.
    """
    def __init__(self):
        """
            Initialize the class
        """
        return

    def __str__(self):
        """
            String output for debugging
        """
        txt = "%s (%s)" % (self.residue, self.optinstance.opttype)
        return txt

    def debug(self, txt):
        """
            Easy way to turn on/off debugging
        """
        if HDEBUG: print txt

    def getHbondangle(self, atom1, atom2, atom3):
        """
            Get the angle between three atoms

            Parameters
                atom1:  The first atom (atom)
                atom2:  The second (vertex) atom (atom)
                atom3:  The third atom (atom)
            Returns
                angle:  The angle between the atoms (float)
        """
        angle = 0.0
        atom2Coords = atom2.getCoords()
        coords1 = subtract(atom3.getCoords(), atom2Coords)
        coords2 = subtract(atom1.getCoords(), atom2Coords)
        norm1 = normalize(coords1)
        norm2 = normalize(coords2)
        dotted = dot(norm1, norm2)
        if dotted > 1.0: # If normalized, this is due to rounding error
            dotted = 1.0
        rad = abs(math.acos(dotted))
        angle = rad*180.0/math.pi
        if angle > 180.0:
            angle = 360.0 - angle
        return angle

    def isHbond(self, donor, acc):
        """
            Determine whether this donor acceptor pair is a
            hydrogen bond
        """
        for donorhatom in donor.bonds:
            if not donorhatom.isHydrogen(): continue
        
            # Check the H(D)-A distance

            dist = distance(donorhatom.getCoords(), acc.getCoords())
          
            if dist > DIST_CUTOFF: continue

            # Ensure no conflicts if H(A)s if present

            flag = 1
            for acchatom in acc.bonds:
                if not acchatom.isHydrogen(): continue

                flag = 0

                # Check the H(D)-H(A) distance
                    
                hdist = distance(donorhatom.getCoords(), acchatom.getCoords())
                if hdist < 1.5: continue
                
                # Check the H(D)-H(A)-A angle
                
                angle = self.getHbondangle(donorhatom, acchatom, acc)
                if angle < 110.0: flag = 1

            if flag == 0: continue

            # Check the A-D-H(D) angle

            angle = self.getHbondangle(acc, donor, donorhatom)
                
            if angle <= ANGLE_CUTOFF: 
                self.debug("Found HBOND! %.4f %.4f" % (dist, angle))
                return 1

        # If we get here, no bond is formed

        return 0
        
    def getPairEnergy(self, donor, acceptor):
        """
            Get the energy between two atoms

            Parameters
                donor:    The first atom in the pair (Atom)
                acceptor: The second atom in the pair (Atom)
            Returns
                energy:   The energy of the pair (float)
        """

        # Initialize some variables
        
        max_hbond_energy = -10.0
        max_ele_energy = -1.0
        maxangle = ANGLE_CUTOFF
        max_dha_dist = DIST_CUTOFF
        max_ele_dist = 5.0
        energy = 0.0
        donorhs = []
        acceptorhs = []
        
        if not (donor.hdonor and acceptor.hacceptor): return energy

        # See if hydrogens are presently bonded to the acceptor and donor

        for bond in donor.bonds:
            if bond.isHydrogen(): donorhs.append(bond)
        for bond in acceptor.bonds:
            if bond.isHydrogen(): acceptorhs.append(bond)
            
        if donorhs == []: return energy

        for donorhatom in donorhs:

            dist = distance(donorhatom.getCoords(), acceptor.getCoords())
            if dist > max_dha_dist and dist < max_ele_dist: 
                energy += max_ele_energy/(dist*dist)
                continue

            if acceptorhs != []:

                # Case 1: Both donor and acceptor hydrogens are present

                for acceptorhatom in acceptorhs:

                    # Penalize if H(D) is too close to H(A)
                    
                    hdist = distance(donorhatom.getCoords(), acceptorhatom.getCoords())
                    if hdist < 1.5:
                        energy += -1 * max_hbond_energy
                        continue

                    # Assign energies based on angles
                    
                    angle1 = self.getHbondangle(acceptor, donor, donorhatom)
                    if angle1 <= maxangle:
                        angleterm = (maxangle - angle1)/maxangle
                        angle2 = self.getHbondangle(donorhatom, acceptorhatom, acceptor)
                        if angle2 < 110.0: angle2 = 1.0
                        else: angle2 = (110.0 - angle2)/110.0
                        energy += max_hbond_energy/pow(dist,3)*angleterm*angle2

            else:

                # Assign energies based on A-D-H(D) angle alone
              
                angle1 = self.getHbondangle(acceptor, donor, donorhatom)
                if angle1 <= maxangle:
                    angleterm = (maxangle - angle1)/maxangle
                    energy += max_hbond_energy/pow(dist,2)*angleterm
                  
            return energy
                
    def makeAtomWithNoBonds(self, atom, closeatom, addname):
        """
            Called for water oxygen atoms with no current bonds.
            Uses the closeatom to place the new atom directly
            colinear with the atom and the closeatom.

            Parameters
                atom:      The oxygen atom of the water
                closeatom: The nearby atom (donor/acceptor)
                addname:   The name of the atom to add
        """
        newcoords = []
        residue = atom.residue
        
        # Place along line, 1 A away
     
        vec = subtract(closeatom.getCoords(), atom.getCoords())
        dist = distance(atom.getCoords(), closeatom.getCoords())

        for i in range(3):
            newcoords.append(vec[i]/dist + atom.getCoords()[i])
            
        residue.createAtom(addname, newcoords)
        newatom = residue.getAtom(addname)
        self.routines.cells.addCell(newatom)

        # Set the bonds (since not in reference structure)

        if newatom not in atom.bonds: atom.bonds.append(newatom)
        if atom not in newatom.bonds: newatom.bonds.append(atom)

    def makeWaterWithOneBond(self, atom, addname):
        """
            Add an atom to a water residue that already has
            one bond.  Uses the water reference structure to
            align the new atom.
        """
        residue = atom.residue
        nextatom = atom.bonds[0]
        coords = [atom.getCoords(), nextatom.getCoords()]
        refcoords = [residue.reference.map[atom.name].getCoords(), \
                     residue.reference.map["H1"].getCoords()]
        refatomcoords = residue.reference.map["H2"].getCoords()

        # Make the atom

        newcoords = findCoordinates(2, coords, refcoords, refatomcoords)
        residue.createAtom(addname, newcoords)

        # Set the bonds (since not in reference structure)

        newatom = residue.getAtom(addname)
        if newatom not in atom.bonds: atom.bonds.append(newatom)
        if atom not in newatom.bonds: newatom.bonds.append(atom)


    def makeAtomWithOneBondH(self, atom, addname):
        """
            Add a hydrogen to an alcoholic donor with one
            existing bond.
        """
   
        residue = atom.residue
        nextatom = atom.bonds[0]
        coords = [atom.getCoords(), nextatom.getCoords()]
        refcoords = [residue.reference.map[atom.name].getCoords(), \
                     residue.reference.map[nextatom.name].getCoords()]
        refatomcoords = residue.reference.map[addname].getCoords()

        # Make the atom

        newcoords = findCoordinates(2, coords, refcoords, refatomcoords)
        residue.createAtom(addname, newcoords)

    def makeAtomWithOneBondLP(self, atom, addname):
        """
            Add a lone pair to an alcoholic donor with one
            existing bond.
        """

        # Initialize some variables
        
        residue = atom.residue
        for refname in atom.reference.bonds:
            if refname.startswith("H"): break

        nextatom = atom.bonds[0]
        coords = [atom.getCoords(), nextatom.getCoords()]
        refcoords = [residue.reference.map[atom.name].getCoords(), \
                     residue.reference.map[nextatom.name].getCoords()]
        refatomcoords = residue.reference.map[refname].getCoords()

        # Make the atom

        newcoords = findCoordinates(2, coords, refcoords, refatomcoords)
        residue.createAtom(addname, newcoords)     

        # Set the bonds (since not in reference structure)

        newatom = residue.getAtom(addname)
        if newatom not in atom.bonds: atom.bonds.append(newatom)
        if atom not in newatom.bonds: newatom.bonds.append(atom)

    def trySingleAlcoholicH(self, donor, acc, newatom):
        """
            After a new bond has been added using
            makeAtomWithOneBond*, try to find the best orientation
            by rotating to form a hydrogen bond.  If a bond
            cannot be formed, remove the newatom (thereby
            returning to a single bond).
        """
        # Initialize some variables

        besten = 999.99
        bestcoords = []
        residue = donor.residue
        pivot = donor.bonds[0]

        for i in range(72):
            residue.rotateTetrahedral(pivot, donor, 5.0)
            if self.isHbond(donor, acc):
                energy = self.getPairEnergy(donor, acc)
                if energy < besten:
                    bestcoords = newatom.getCoords()
                    besten = energy
        
        # If a hydrogen bond was made, set at best coordinates

        if bestcoords != []:
            newatom.x = bestcoords[0]
            newatom.y = bestcoords[1]
            newatom.z = bestcoords[2]
            self.routines.cells.addCell(newatom)
            return 1
        else:
            residue.removeAtom(newatom.name)
            return 0

    def trySingleAlcoholicLP(self, acc, donor, newatom):
        """
            After a new bond has been added using
            makeAtomWithOneBond*, ensure that a
            hydrogen bond has been made.  If so, try to
            minimze the H(D)-A-LP angle.  If that cannot
            be minimized, ignore the bond and remove the
            atom.
        """
     
        # Initialize some variables

        residue = acc.residue
        pivot = acc.bonds[0]
        bestangle = 180.00
        bestcoords = []
        
        # If a hydrogen bond was made, set at best distance
        
        if not self.isHbond(donor, acc):
            residue.removeAtom(newatom.name)
            return 0

        # Grab the H(D) that caused the bond

        for donorhatom in donor.bonds:
            if donorhatom.isHydrogen() and \
               self.getHbondangle(acc, donor, donorhatom) < ANGLE_CUTOFF: break
        
        for i in range(72):
            residue.rotateTetrahedral(pivot, acc, 5.0)
            angle = abs(self.getHbondangle(donorhatom, acc, newatom))
          
            if angle < bestangle:
                bestangle = angle
                bestcoords = newatom.getCoords()

        # Remove if geometry does not work
                            
        if bestangle > (ANGLE_CUTOFF * 2.0):
            self.debug("Removing due to geometry %.2f > %.2f" % (bestangle, ANGLE_CUTOFF*2.0))
            residue.removeAtom(newatom.name)
            return 0

        # Otherwise set to best coordinates
            
        newatom.x = bestcoords[0]
        newatom.y = bestcoords[1]
        newatom.z = bestcoords[2]
        self.routines.cells.addCell(newatom)
            
        return 1
       

    def getPositionsWithTwoBonds(self, atom):
        """
            Given a tetrahedral geometry with two
            existing bonds, return the two potential
            sets of coordinates that are possible for
            a new bond.
        """

        # Initialize some variables

        residue = atom.residue
        fixed = atom.bonds[0]
        rotate = atom.bonds[1]

        # Rotate by 120 degrees twice
      
        residue.rotateTetrahedral(fixed, atom, 120)
        loc1 = rotate.getCoords()  
        residue.rotateTetrahedral(fixed, atom, 120)
        loc2 = rotate.getCoords()
    
        # Set rotate back to original by one more rotation

        residue.rotateTetrahedral(fixed, atom, 120)

        return loc1, loc2

    def tryPositionsWithTwoBondsH(self, donor, acc, newname, loc1, loc2):
        """
            Try adding a new hydrogen two the two potential
            locations.  If both form hydrogen bonds, place at
            whatever returns the best bond as determined by
            getPairEnergy.
        """

        # Initialize some variables

        besten = 999.99
        bestcoords = []
        residue = donor.residue
        
        # Try the first position
        
        residue.createAtom(newname, loc1)
        if self.isHbond(donor, acc):
            besten = self.getPairEnergy(donor, acc)
            bestcoords = loc1

        # Try the second

        newatom = residue.getAtom(newname)
        newatom.x = loc2[0]
        newatom.y = loc2[1]
        newatom.z = loc2[2]
        if self.isHbond(donor, acc):
            energy = self.getPairEnergy(donor, acc)
            if energy < besten:
                bestcoords = loc2

        # Set at best coords

        if bestcoords != []:
            newatom.x = bestcoords[0]
            newatom.y = bestcoords[1]
            newatom.z = bestcoords[2]
            self.routines.cells.addCell(newatom)
            return 1
        else: 
            residue.removeAtom(newname)
            return 0
        
    def tryPositionsWithTwoBondsLP(self, acc, donor, newname, loc1, loc2):
        """
            Try placing an LP on a tetrahedral geometry with
            two existing bonds.  If this isn't a hydrogen bond
            it can return - otherwise ensure that the H(D)-A-LP
            angle is minimized.
        """

        # Initialize some variables

        bestangle = 180.00
        bestcoords = []
        residue = acc.residue

        # If the donor/acceptor pair is not an hbond return

        if not self.isHbond(donor, acc): return 0

        # Grab the H(D) that caused the bond

        for donorhatom in donor.bonds:
            if donorhatom.isHydrogen() and \
               self.getHbondangle(acc, donor, donorhatom) < ANGLE_CUTOFF: break

        # Try the first position
        
        residue.createAtom(newname, loc1)
        newatom = residue.getAtom(newname)
        angle = abs(self.getHbondangle(donorhatom, acc, newatom))
        if angle < bestangle:
            bestangle = angle
            bestcoords = loc1

        # Try the second

        newatom.x = loc2[0]
        newatom.y = loc2[1]
        newatom.z = loc2[2]
        angle = self.getHbondangle(donorhatom, acc, newatom)
        if angle < bestangle:
            bestcoords = loc2

        # Remove if geometry does not work

        if bestangle > (ANGLE_CUTOFF * 2.0):
            residue.removeAtom(newname)
            return 0
        
        # Otherwise set at best coords
        
        newatom.x = bestcoords[0]
        newatom.y = bestcoords[1]
        newatom.z = bestcoords[2]
        self.routines.cells.addCell(newatom)
     
        # Set the bonds (since not in reference structure)

        if newatom not in acc.bonds: acc.bonds.append(newatom)
        if acc not in newatom.bonds: newatom.bonds.append(acc)
        
        return 1

    def getPositionWithThreeBonds(self, atom):
        """
           If there's three bonds in a tetrahedral geometry,
           there's only one available position.  Find that
           position.
        """

        # Initialize some variables

        residue = atom.residue
        pivot = atom.bonds[0]
        rot1 = atom.bonds[1]
        rot2 = atom.bonds[2]

        # Find the two new positions

        residue.rotateTetrahedral(pivot, atom, 120)
        newcoords1 = rot1.getCoords()
        residue.rotateTetrahedral(pivot, atom, 120)
        newcoords2 = rot1.getCoords()
        residue.rotateTetrahedral(pivot, atom, 120)

        # Determine which is unoccupied

        if distance(rot2.getCoords(), newcoords1) > 0.1: return newcoords1
        else: return newcoords2    

    def tryPositionsWithThreeBondsH(self, donor, acc, newname, loc):
        """
            Try making a hydrogen bond with the lone available
            position.
        """
        residue = donor.residue
        residue.createAtom(newname, loc)
        if self.isHbond(donor, acc):
            newatom = residue.getAtom(newname)
            self.routines.cells.addCell(newatom)
            return 1
        else:
            residue.removeAtom(newname)
            return 0

    
    def tryPositionsWithThreeBondsLP(self, acc, donor, newname, loc):
        """
            Try making a hydrogen bond using the lone
            available hydrogen position.
        """
        residue = acc.residue
        if not self.isHbond(donor, acc): return 0

        # Grab the H(D) that caused the bond

        for donorhatom in donor.bonds:
            if donorhatom.isHydrogen() and \
               self.getHbondangle(acc, donor, donorhatom) < ANGLE_CUTOFF: break

        residue.createAtom(newname, loc)
        newatom = residue.getAtom(newname)

        # Remove if geometry does not work

        angle = abs(self.getHbondangle(donorhatom, acc, newatom))   
        if angle > (ANGLE_CUTOFF * 2.0):
            residue.removeAtom(newname)
            return 0

        # Otherwise keep it

        newatom = residue.getAtom(newname)
        self.routines.cells.addCell(newatom)
     
        # Set the bonds (since not in reference structure)

        if newatom not in acc.bonds: acc.bonds.append(newatom)
        if acc not in newatom.bonds: newatom.bonds.append(acc)

        return 1

class Flip(Optimize):
    """
        The holder for optimization of flippable residues.
    """
    def __init__(self, residue, optinstance, routines):
        """
            Initialize a potential flip.  Rather than flipping
            the given residue back and forth, take each atom
            that would be flipped and pre-flip it, making a
            new *FLIP atom in its place.

            Parameters
                residue:      The residue to flip (residue)
                optinstance:  The optimization instance containing
                              information about what to optimize
                
        """
        # Initialize some variables

        self.optinstance = optinstance
        self.residue = residue
        self.routines = routines
        self.atomlist = []
        self.hbonds = []

        map = {}

        # Get all moveable names for this angle/residue pair

        dihedral = optinstance.optangle
        pivot = dihedral.split()[2]
        moveablenames = self.routines.getMoveableNames(residue, pivot)

        # Cache current coordinates

        for name in moveablenames:
            atom = residue.getAtom(name)
            map[name] = atom.getCoords()
          
        # Flip the residue about the angle

        anglenum = residue.reference.dihedrals.index(dihedral)
        if anglenum == -1:
            raise ValueError, "Unable to find dihedral angle!"

        newangle = 180.0 + residue.dihedrals[anglenum]
        self.routines.setDihedralAngle(residue, anglenum, newangle) 

        # Create new atoms at cached positions

        for name in map:
            newname = "%sFLIP" % name
            residue.createAtom(newname, map[name])
            newatom = residue.getAtom(newname)
            self.routines.cells.addCell(newatom)

            # Set the bonds
            
            newatom.reference = residue.reference.map[name]
            for bond in newatom.reference.bonds:
                newbond = "%sFLIP" % bond
                if residue.hasAtom(newbond):
                    bondatom = residue.map[newbond]
                    if bondatom not in newatom.bonds: newatom.bonds.append(bondatom)
                    if newatom not in bondatom.bonds: bondatom.bonds.append(newatom)

                # And connect back to the existing structure

                newbond = bond
                if residue.hasAtom(newbond):
                    bondatom = residue.map[newbond]
                    if bondatom not in newatom.bonds: newatom.bonds.append(bondatom)
                    if newatom not in bondatom.bonds: bondatom.bonds.append(newatom)   

        residue.setDonorsAndAcceptors()

        # Add to the optimization list

        for name in moveablenames:

            # Get the atom
            
            atom = residue.getAtom(name)
            if not atom.isHydrogen() and \
               (atom.hdonor or atom.hacceptor): self.atomlist.append(atom)

            # And the FLIP

            atom = residue.getAtom("%sFLIP" % name)
            if not atom.isHydrogen() and \
               (atom.hdonor or atom.hacceptor): self.atomlist.append(atom)

        # Special case: Neutral unassigned HIS can be acceptors

        if isinstance(residue, HIS):
            if "HIS" == residue.name and \
               len(residue.patches) == 1:
                for atom in self.atomlist:
                    if atom.name.startswith("N"):
                        atom.hacceptor = 1

    def tryBoth(self, donor, acc, accobj):
        """
           Called when both the donor and acceptor are optimizeable;
           If one is fixed, we only need to try one side.  Otherwise
           first try to satisfy the donor - if that's succesful,
           try to satisfy the acceptor.  An undo may be necessary
           if the donor is satisfied and the acceptor isn't.
        """
        # If one residue if fixed, use other functions

        if donor.residue.fixed:
            if accobj.tryAcceptor(acc, donor): return 1
            else: return 0
        if acc.residue.fixed:
            if self.tryDonor(donor, acc): return 1
            else: return 0

    
        self.debug("Working on %s %s (donor) to %s %s (acceptor)" % \
                       (donor.residue, donor.name, acc.residue, acc.name))
        if self.isHbond(donor, acc):
            if accobj.tryAcceptor(acc, donor):
                self.fixFlip(donor)
                donor.hacceptor = 0
                self.debug("NET BOND SUCCESSFUL!")
                return 1
            else:
                return 0

    def tryDonor(self, donor, acc):
        """
           The main driver for adding a hydrogen to an
           optimizeable residue.
        """
        residue = self.residue

        # Do some error checking
        
        if not acc.hacceptor: return 0
   
        self.debug("Working on %s %s (donor) to %s %s (acceptor)" % \
                       (donor.residue, donor.name, acc.residue, acc.name))
            
        if self.isHbond(donor, acc):
            residue.fixed = donor.name
            self.fixFlip(donor)
            donor.hacceptor = 0
            return 1
        else:
            return 0

    def tryAcceptor(self, acc, donor):
        """
           The main driver for adding an LP to an optimizeable
           residue.
        """
        residue = acc.residue
        
        # Do some error checking

        if not donor.hdonor: return 0            
            
        self.debug("Working on %s %s (acceptor) to %s %s (donor)" % \
                   (acc.residue, acc.name, donor.residue, donor.name))
        if self.isHbond(donor, acc):
            residue.fixed = acc.name
            self.fixFlip(acc)
            acc.hdonor = 0
            return 1
        else: return 0

    def fixFlip(self, bondatom):
        """
           Called if a hydrogen bond has been found using
           the bondatom.  If bondatom is *FLIP, remove all *
           atoms, otherwise remove all *FLIP atoms. 
        """
   
        # Initialize some variables

        atomlist = []
        residue = bondatom.residue
        for atom in residue.getAtoms(): atomlist.append(atom)
         
        # Set a flag to see whether to delete the FLIPs or not
        
        flag = 0
        if bondatom.name.endswith("FLIP"): flag = 1
     
        # Delete the appropriate atoms
        
        for atom in atomlist:
            atomname = atom.name
            if atomname.endswith("FLIP") and flag: # Delete the other list
                if residue.hasAtom(atomname[:-4]):
                    self.routines.cells.removeCell(residue.getAtom(atomname[:-4]))
                    residue.removeAtom(atomname[:-4])  
            elif atomname.endswith("FLIP"):  # Delete the flip
                self.routines.cells.removeCell(atom)
                residue.removeAtom(atomname)
            else: continue

        residue.fixed = 1
            
    def finalize(self):
        """
            Finalizes a flippable back to its original state -
            since the original atoms are now *FLIP, it deletes
            the * atoms and renames the *FLIP atoms back to *.
        """
        residue = self.residue

        if residue.fixed: return
        atomlist = []
        for atom in residue.getAtoms(): atomlist.append(atom)
        for atom in atomlist:
            if atom.name.endswith("FLIP"):
                self.routines.cells.removeCell(atom)
                residue.removeAtom(atom.name[:-4])
                residue.renameAtom(atom.name, atom.name[:-4])
        residue.fixed = 1

    def complete(self):
        """
            Complete the flippable residue optimization.  Call the finalize
            function, and then rename all FLIP atoms back to their standard
            names.
        """
        residue = self.residue
        
        self.finalize()

        # Rename all *FLIP atoms
        for atom in residue.getAtoms():
            atomname = atom.name
            if atomname.endswith("FLIP"):
                residue.renameAtom(atomname, atomname[:-4]) 
        
class Alcoholic(Optimize):
    """
        The class for alcoholic residues
    """

    def __init__(self, residue, optinstance, routines):
        """
           Initialize the alcoholic class by removing
           the alcoholic hydrogen if it exists.
        """
        
        # Initialize some variables

        self.optinstance = optinstance
        self.residue = residue
        self.routines = routines
        self.atomlist = []
        self.hbonds = []

        name = optinstance.map.keys()[0]
        self.hname = name
        
        bondname = residue.reference.getAtom(name).bonds[0]
        self.atomlist.append(residue.getAtom(bondname))
        if residue.hasAtom(name):
            atom = residue.getAtom(name)
            self.routines.cells.removeCell(atom)
            residue.removeAtom(name)

    def tryBoth(self, donor, acc, accobj):
        """
           Called when both the donor and acceptor are optimizeable;
           If one is fixed, we only need to try one side.  Otherwise
           first try to satisfy the donor - if that's succesful,
           try to satisfy the acceptor.  An undo may be necessary
           if the donor is satisfied and the acceptor isn't.
        """
        # If one residue if fixed, use other functions

        residue = donor.residue
        if donor.residue.fixed:
            if accobj.tryAcceptor(acc, donor): return 1
            else: return 0
        if acc.residue.fixed:
            if self.tryDonor(donor, acc): return 1
            else: return 0

        if self.tryDonor(donor, acc):
            if accobj.tryAcceptor(acc, donor):
                self.debug("NET BOND SUCCESSFUL!")
                return 1
            else: # We need to remove the added H
                residue.removeAtom(self.hname)
                self.debug("REMOVED NET HBOND") 
                return 0
        else:
            return 0

    def tryDonor(self, donor, acc):
        """
           The main driver for adding a hydrogen to an
           optimizeable residue.
        """
        residue = self.residue

        # Do some error checking
        
        if not acc.hacceptor: return 0
            
        # Get the name of the atom to add

        newname = self.hname
        if residue.hasAtom(newname): return 0
     
        self.debug("Working on %s %s (donor) to %s %s (acceptor)" % \
                   (donor.residue, donor.name, acc.residue, acc.name))
        
        # Act depending on the number of bonds

        if len(donor.bonds) == 1: # No H or LP attached
            self.makeAtomWithOneBondH(donor, newname)
            newatom = donor.residue.getAtom(newname)
            return self.trySingleAlcoholicH(donor, acc, newatom)
        elif len(donor.bonds) == 2:
            loc1, loc2 = self.getPositionsWithTwoBonds(donor)
            return self.tryPositionsWithTwoBondsH(donor, acc, newname, loc1, loc2)
        elif len(donor.bonds) == 3:
            loc = self.getPositionWithThreeBonds(donor)
            return self.tryPositionsWithThreeBondsH(donor, acc, newname, loc)


    def tryAcceptor(self, acc, donor):
        """
           The main driver for adding an LP to an optimizeable
           residue.
        """
        residue = acc.residue
        
        # Do some error checking

        if not donor.hdonor: return 0            
     
        # Get the name of the LP to add

        if residue.hasAtom("LP2"): return 0
        elif residue.hasAtom("LP1"): newname = "LP2"
        else: newname = "LP1"
        
        self.debug("Working on %s %s (acceptor) to %s %s (donor)" % \
                   (acc.residue, acc.name, donor.residue, donor.name))

        # Act depending on the number of bonds

        if len(acc.bonds) == 1: # No H or LP attached
            self.makeAtomWithOneBondLP(acc, newname)
            newatom = acc.residue.getAtom(newname)
            return self.trySingleAlcoholicLP(acc, donor, newatom)      
        elif len(acc.bonds) == 2:
            loc1, loc2 = self.getPositionsWithTwoBonds(acc)
            return self.tryPositionsWithTwoBondsLP(acc, donor, newname, loc1, loc2)
        elif len(acc.bonds) == 3:
            loc = self.getPositionWithThreeBonds(acc)
            return self.tryPositionsWithThreeBondsLP(acc, donor, newname, loc)  
    
    def finalize(self):
        """
            Finalize an alcoholic residue.  Try to minimize
            conflict with nearby atoms by building away
            from them.  Called when LPs are still present
            so as to account for their bonds.
        """

        # Initialize some variables

        residue = self.residue
        atom = self.atomlist[0]

        # Conditions for return

        addname = self.hname
        if residue.fixed: return
        if residue.hasAtom(addname): return

        if len(atom.bonds) == 1:

            # Initialize variables
            
            pivot = atom.bonds[0]
            bestdist = 0.0
            bestcoords = []

            # Add atom and debump

            self.makeAtomWithOneBondH(atom, addname)
            newatom = residue.getAtom(addname)
            self.routines.cells.addCell(newatom)

            for i in range(18):
                residue.rotateTetrahedral(pivot, atom, 20.0)
                nearatom = self.routines.getClosestAtom(newatom)

                # If there is no closest conflict, return

                if nearatom == None: return
                
                dist = distance(nearatom.getCoords(), newatom.getCoords())
                if dist > bestdist:
                    bestdist = dist
                    bestcoords = newatom.getCoords()
                
            if bestcoords != []:
                newatom.x = bestcoords[0]
                newatom.y = bestcoords[1]
                newatom.z = bestcoords[2]
            
        elif len(atom.bonds) == 2:    
            loc1, loc2 = self.getPositionsWithTwoBonds(atom)
            residue.createAtom(addname, loc1)
            newatom = residue.getAtom(addname)
            self.routines.cells.addCell(newatom)
            
            # Debump residue if necessary by trying the other location
            
            nearatom = self.routines.getClosestAtom(newatom)    
            if nearatom == None: return     
            dist1 = distance(newatom.getCoords(), nearatom.getCoords())

            # Place at other location

            self.routines.cells.removeCell(atom)
            newatom.x = loc2[0]
            newatom.y = loc2[1]
            newatom.z = loc2[2]
            self.routines.cells.addCell(atom)

            nearatom = self.routines.getClosestAtom(newatom)
            if nearatom == None: return

            # If this is worse, switch back
            
            if distance(newatom.getCoords(), nearatom.getCoords()) < dist1:
                self.routines.cells.removeCell(atom)
                newatom.x = loc1[0]
                newatom.y = loc1[1]
                newatom.z = loc1[2]
                self.routines.cells.addCell(atom)
         
        elif len(atom.bonds) == 3:
          
            loc = self.getPositionWithThreeBonds(atom)
            residue.createAtom(addname, loc)
            self.routines.cells.addCell(residue.getAtom(addname))

    def complete(self):
        """
            Complete an alcoholic optimization.  Call finalize(), and then
            remove all extra LP atoms.
        """
        # Initialize some variables

        residue = self.residue
   
        self.finalize()
        residue.fixed = 1
        
        # Remove all LP atoms
        
        atomlist = []
        for atom in residue.getAtoms(): atomlist.append(atom)
        for atom in atomlist:
            if atom.name.startswith("LP"): residue.removeAtom(atom.name)

class Water(Optimize):
    """
        The class for water residues
    """
    def __init__(self, residue, optinstance, routines):
        """
             Initialize the water optimization class
        """
        self.optinstance = optinstance
        self.residue = residue
        self.routines = routines
        self.hbonds = []

        oxatom = residue.getAtom("O")
        if oxatom == None:
            raise ValueError, "Unable to find oxygen atom in %s!" % residue

        oxatom.hdonor = 1
        oxatom.hacceptor = 1
      
        self.atomlist = [oxatom]

    def tryBoth(self, donor, acc, accobj):
        """
           Called when both the donor and acceptor are optimizeable;
           If one is fixed, we only need to try one side.  Otherwise
           first try to satisfy the donor - if that's succesful,
           try to satisfy the acceptor.  An undo may be necessary
           if the donor is satisfied and the acceptor isn't.
        """
        # If one residue if fixed, use other functions

        residue = donor.residue

        if donor.residue.fixed:
            if accobj.tryAcceptor(acc, donor): return 1
            else: return 0
        if acc.residue.fixed:
            if self.tryDonor(donor, acc): return 1
            else: return 0

        if self.tryDonor(donor, acc):
            if accobj.tryAcceptor(acc, donor):
                self.debug("NET BOND SUCCESSFUL!")
                return 1
            else:  
                # We need to undo what we did to the donor

                self.debug("REMOVED NET HBOND")        
                if residue.hasAtom("H2"): residue.removeAtom("H2")      
                elif residue.hasAtom("H1"): residue.removeAtom("H1")
                return 0
        else:
            return 0

    def tryAcceptor(self, acc, donor):
        """
           The main driver for adding an LP to an optimizeable
           residue.
        """
        residue = acc.residue
        
        # Do some error checking

        if not donor.hdonor: return 0            
     

        # Get the name of the LP to add
        
        if residue.hasAtom("LP2"): return 0
        elif residue.hasAtom("LP1"): newname = "LP2"
        else: newname = "LP1"

        self.debug("Working on %s %s (acceptor) to %s %s (donor)" % \
        (acc.residue, acc.name, donor.residue, donor.name))

        # Act depending on the number of bonds

        if len(acc.bonds) == 0:

            if self.isHbond(donor, acc):

                # Find the best donor hydrogen and use that 

                besth = donor
                bestdist = distance(acc.getCoords(), donor.getCoords())
                for donorh in donor.bonds:
                    dist = distance(acc.getCoords(), donorh.getCoords())
                    if dist < bestdist:
                        besth = donorh
                        bestdist = dist
                        
                # Point the LP to the best H
                
                self.makeAtomWithNoBonds(acc, donorh, newname)
                self.debug("Added %s to %s" % (newname, acc.residue))
                return 1
                
            else: return 0
                
        elif len(acc.bonds) == 1: # No H or LP attached
            self.debug("Trying to add %s to %s with one bond" % (newname, acc.residue))
            self.makeWaterWithOneBond(acc, newname)
            newatom = acc.residue.getAtom(newname)
            return self.trySingleAlcoholicLP(acc, donor, newatom)
        elif len(acc.bonds) == 2:
            self.debug("Trying to add %s to %s with two bonds" % (newname, acc.residue))
            loc1, loc2 = self.getPositionsWithTwoBonds(acc)
            return self.tryPositionsWithTwoBondsLP(acc, donor, newname, loc1, loc2)
        elif len(acc.bonds) == 3:
            self.debug("Trying to add %s to %s with three bonds" % (newname, acc.residue))
            loc = self.getPositionWithThreeBonds(acc)
            return self.tryPositionsWithThreeBondsLP(acc, donor, newname, loc)  

    def tryDonor(self, donor, acc):
        """
           The main driver for adding a hydrogen to an
           optimizeable residue.
        """
        residue = self.residue

        # Do some error checking
        
        if not acc.hacceptor: return 0
   
        # Get the name of the atom to add
        
        if residue.hasAtom("H2"): return 0
        elif residue.hasAtom("H1"): newname = "H2"
        else: newname = "H1"
       
        self.debug("Working on %s %s (donor) to %s %s (acceptor)" % \
                   (donor.residue, donor.name, acc.residue, acc.name))

        # Act depending on the number of bonds
        if len(donor.bonds) == 0:
            self.makeAtomWithNoBonds(donor, acc, newname)
            if self.isHbond(donor, acc): return 1
            else:
                self.routines.cells.removeCell(residue.getAtom(newname))
                residue.removeAtom(newname)
                return 0
        if len(donor.bonds) == 1:
            self.makeWaterWithOneBond(donor, newname)
            newatom = donor.residue.getAtom(newname)
            return self.trySingleAlcoholicH(donor, acc, newatom)
        elif len(donor.bonds) == 2:
            loc1, loc2 = self.getPositionsWithTwoBonds(donor)
            return self.tryPositionsWithTwoBondsH(donor, acc, newname, loc1, loc2)
        elif len(donor.bonds) == 3:
            loc = self.getPositionWithThreeBonds(donor)
            return self.tryPositionsWithThreeBondsH(donor, acc, newname, loc)

    def finalize(self):
        """
            Finalize a water residue.  Try to minimize
            conflict with nearby atoms by building away
            from them.  Called when LPs are still present
            so as to account for their bonds.
        """

        residue = self.residue

        # Conditions for return

        if residue.fixed:  
            self.debug("Residue %s already fixed" % residue)
            return
        if residue.hasAtom("H2"): 
            self.debug("Residue %s already has H2" % residue)
            return

        atom = residue.getAtom("O")
        if not residue.hasAtom("H1"): addname = "H1"
        else:  addname = "H2"
       
        self.debug("Finalizing %s by adding %s (%i current O bonds)" % (residue, addname, len(atom.bonds)))
 
        if len(atom.bonds) == 0:

            newcoords = []

            # Build hydrogen away from closest atom
            
            closeatom = self.routines.getClosestAtom(atom)
            if closeatom != None:
                vec = subtract(atom.getCoords(), closeatom.getCoords())
                dist = distance(atom.getCoords(), closeatom.getCoords())

                for i in range(3):
                    newcoords.append(vec[i]/dist + atom.getCoords()[i])
                
            else:
                newcoords = add(atom.getCoords(), [1.0, 0.0, 0.0])      

            residue.createAtom(addname, newcoords)
            self.routines.cells.addCell(residue.getAtom(addname))

            self.finalize()
          
        elif len(atom.bonds) == 1:

            # Initialize variables
            
            pivot = atom.bonds[0]
            bestdist = 0.0
            bestcoords = []

            # Add atom and debump

            self.makeWaterWithOneBond(atom, addname)
            newatom = residue.getAtom(addname)
            self.routines.cells.addCell(newatom)

            for i in range(18):
                residue.rotateTetrahedral(pivot, atom, 20.0)
                nearatom = self.routines.getClosestAtom(newatom)

                # If there is no closest conflict, continue

                if nearatom == None: continue
    
                dist = distance(nearatom.getCoords(), newatom.getCoords())
                
                if dist > bestdist:
                    bestdist = dist
                    bestcoords = newatom.getCoords()
                
            if bestcoords != []:
                newatom.x = bestcoords[0]
                newatom.y = bestcoords[1]
                newatom.z = bestcoords[2]

            if addname == "H1":
                self.finalize()

            residue.fixed = 1

        elif len(atom.bonds) == 2:

            loc1, loc2 = self.getPositionsWithTwoBonds(atom)
            residue.createAtom(addname, loc1)
            newatom = residue.getAtom(addname)
            self.routines.cells.addCell(newatom)
            
            # Debump residue if necessary by trying the other location
            
            nearatom = self.routines.getClosestAtom(newatom)    
            if nearatom != None:    
                dist1 = distance(newatom.getCoords(), nearatom.getCoords())

                # Place at other location

                self.routines.cells.removeCell(atom)
                newatom.x = loc2[0]
                newatom.y = loc2[1]
                newatom.z = loc2[2]
                self.routines.cells.addCell(atom)

                nearatom = self.routines.getClosestAtom(newatom)
                if nearatom != None:

                    # If this is worse, switch back
            
                    if distance(newatom.getCoords(), nearatom.getCoords()) < dist1:
                        self.routines.cells.removeCell(atom)
                        newatom.x = loc1[0]
                        newatom.y = loc1[1]
                        newatom.z = loc1[2]
                        self.routines.cells.addCell(atom)
       
            if addname == "H1":
                self.finalize()    

        elif len(atom.bonds) == 3:

            loc = self.getPositionWithThreeBonds(atom)
            residue.createAtom(addname, loc)
            self.routines.cells.addCell(residue.getAtom(addname))            

          
    def complete(self):
        """
            Complete the water optimization class
        """
        self.finalize()

        residue = self.residue
        
        # Remove all LP atoms
        
        atomlist = []
        for atom in residue.getAtoms(): atomlist.append(atom)
        for atom in atomlist:
            if atom.name.startswith("LP"): residue.removeAtom(atom.name)

class Carboxylic(Optimize):
    """
        The class for carboxylic residues
    """
    def __init__(self, residue, optinstance, routines):
        """
            Initialize a case where the lone hydrogen atom
            can have four different orientations.  Works similar
            to initializeFlip by preadding the necessary atoms.

            This also takes into account that the carboxyl group
            has different bond lengths for the two C-O bonds -
            this is probably due to one bond being assigned
            as a C=O.  As a result hydrogens are only added to
            the C-O (longer) bond.
            
            Parameters
                residue:  The residue to flip (residue)
                dihedral: The angle to flip about
                hname:    The name of one of the hydrogens to add
            Returns
                optlist:  A list of optimizeable donors and
                          acceptors in the residue (list)
        """
        # Initialize some variables

        self.optinstance = optinstance
        self.residue = residue
        self.routines = routines
        self.atomlist = []
        self.hbonds = []
        self.hlist = []

        hname2 = ""
        hname = ""
        for name in optinstance.map.keys():
            if residue.hasAtom(name): hname2 = name
            else: hname = name
        bondatom1 = residue.getAtom(optinstance.map[hname].bond)
        bondatom2 = residue.getAtom(optinstance.map[hname2].bond)
        longflag = 0

        # If one bond in the group is significantly (0.05 A)
        # longer than the other, use that group only

        for pivotatom in bondatom1.bonds:
            if not pivotatom.isHydrogen(): break

        d1 = distance(pivotatom.getCoords(), bondatom1.getCoords())
        d2 = distance(pivotatom.getCoords(), bondatom2.getCoords())

        if d2 > d1 and abs(d1 - d2) > 0.05:
            longflag = 1
            
        elif d1 > d2 and abs(d1 - d2) > 0.05:          
            longflag = 1
            self.routines.cells.removeCell(residue.getAtom(hname2))
            residue.removeAtom(hname2)
        
            # Rebuild on bondatom1

            coords = [bondatom2.getCoords(), bondatom1.getCoords(), pivotatom.getCoords()]
            refcoords = [bondatom2.reference.getCoords(), bondatom1.reference.getCoords(), \
                         pivotatom.reference.getCoords()]
            refatomcoords = residue.reference.map[hname2].getCoords()
            newcoords = findCoordinates(3, coords, refcoords, refatomcoords)

            residue.createAtom(hname2, newcoords)
            self.routines.cells.addCell(residue.getAtom(hname2))

            # Switch the names so the new atom is on the longer bondatom2
            
            bondname1 = bondatom1.name
            bondname2 = bondatom2.name
            residue.renameAtom(bondname1, "TEMP")
            residue.renameAtom(bondname2, bondname1)
            residue.renameAtom("TEMP", bondname2)
       
        # First mirror the hydrogen about the same donor

        for di in residue.reference.dihedrals:
            if di.endswith(hname2): break

        anglenum = residue.reference.dihedrals.index(di)
        if anglenum == -1:
            raise ValueError, "Unable to find dihedral angle!"

        newangle = 180.0 + residue.dihedrals[anglenum]
        self.routines.setDihedralAngle(residue, anglenum, newangle) 

        hatom2 = residue.getAtom(hname2)
        newcoords = hatom2.getCoords()

        # Flip back to return original atom

        newangle = 180.0 + residue.dihedrals[anglenum]
        self.routines.setDihedralAngle(residue, anglenum, newangle) 

        # Rename the original atom and rebuild the new atom

        residue.renameAtom(hname2, "%s1" % hname2)
        newname = "%s2" % hname2
        residue.createAtom(newname, newcoords)
        newatom = residue.getAtom(newname)
        self.routines.cells.addCell(newatom)
        newatom.refdistance = hatom2.refdistance
        
        # Set the bonds for the new atom
            
        if bondatom2 not in newatom.bonds: newatom.bonds.append(bondatom2)
        if newatom not in bondatom2.bonds: bondatom2.bonds.append(newatom)

        # Return if this is the only atom to add

        if longflag:
            residue.setDonorsAndAcceptors()
            self.atomlist = [bondatom2]
            self.hlist = [residue.getAtom("%s1" % hname2)]
            self.hlist.append(residue.getAtom("%s2" % hname2))
            return 

        # Now rotate about the given dihedral angle, mirroring both

        anglenum = residue.reference.dihedrals.index(optinstance.optangle)
        if anglenum == -1:
            raise ValueError, "Unable to find dihedral angle!"

        newangle = 180.0 + residue.dihedrals[anglenum]
        self.routines.setDihedralAngle(residue, anglenum, newangle)

        newcoords = [hatom2.getCoords(), newatom.getCoords()]
   
        # Flip back to return original atom

        newangle = 180.0 + residue.dihedrals[anglenum]
        self.routines.setDihedralAngle(residue, anglenum, newangle) 

        # Make the new atoms

        for i in range(len(newcoords)):
            id = i + 1
            newname = "%s%i" % (hname, id)
            residue.createAtom(newname, newcoords[i])
            newatom = residue.getAtom(newname)
            self.routines.cells.addCell(newatom)

            # Set the bonds
            
            if bondatom1 not in newatom.bonds: newatom.bonds.append(bondatom1)
            if newatom not in bondatom1.bonds: bondatom1.bonds.append(newatom)

        residue.setDonorsAndAcceptors()

        self.hlist = [residue.getAtom("%s1" % hname2)]
        self.hlist.append(residue.getAtom("%s2" % hname2))
        self.hlist.append(residue.getAtom("%s1" % hname))
        self.hlist.append(residue.getAtom("%s2" % hname))

        self.atomlist = [bondatom2, bondatom1]

    def tryBoth(self, donor, acc, accobj):
        """
           Called when both the donor and acceptor are optimizeable;
           If one is fixed, we only need to try one side.  Otherwise
           first try to satisfy the donor - if that's succesful,
           try to satisfy the acceptor.  An undo may be necessary
           if the donor is satisfied and the acceptor isn't.
        """
        # If one residue if fixed, use other functions

        if donor.residue.fixed:
            if accobj.tryAcceptor(acc, donor): return 1
            else: return 0
        if acc.residue.fixed:
            if self.tryDonor(donor, acc): return 1
            else: return 0

        self.debug("Working on %s %s (donor) to %s %s (acceptor)" % \
                   (donor.residue, donor.name, acc.residue, acc.name))

        if self.isHbond(donor, acc):
            if accobj.tryAcceptor(acc, donor):
                self.fix(donor, acc)
                self.debug("NET BOND SUCCESSFUL!")
                return 1
            else:
                return 0

    def isCarboxylicHbond(self, donor, acc):
        """
            Determine whether this donor acceptor pair is a
            hydrogen bond
        """
        for donorhatom in donor.bonds:
            if not donorhatom.isHydrogen(): continue
        
            # Check the H(D)-A distance

            dist = distance(donorhatom.getCoords(), acc.getCoords())
            if dist > DIST_CUTOFF: continue

            # Check the A-D-H(D) angle

            angle = self.getHbondangle(acc, donor, donorhatom)
            if angle <= ANGLE_CUTOFF: 
                self.debug("Found HBOND! %.4f %.4f" % (dist, angle))
                return 1

        # If we get here, no bond is formed

        return 0
    
    def tryAcceptor(self, acc, donor):
        """
           The main driver for adding an LP to an optimizeable
           residue.
        """
        residue = acc.residue
        
        # Do some error checking

        if not donor.hdonor: return 0            

        self.debug("Working on %s %s (acceptor) to %s %s (donor)" % \
                   (acc.residue, acc.name, donor.residue, donor.name))

        # We want to ignore the Hs on the acceptor

        if self.isCarboxylicHbond(donor, acc):

            # Eliminate the closer hydrogen

            hyds = []
            dist = None
            for bond in acc.bonds:
                if bond.isHydrogen(): hyds.append(bond)
                    
            if len(hyds) < 2: return 1

            dist = distance(hyds[0].getCoords(), donor.getCoords())
            dist2 = distance(hyds[1].getCoords(), donor.getCoords())
            if dist < dist2: # Eliminate hyds[0]
                self.hlist.remove(hyds[0])
                self.routines.cells.removeCell(hyds[0])
                residue.removeAtom(hyds[0].name)
                donorhatom = residue.getAtom(hyds[1].name)
            else:
                self.hlist.remove(hyds[1])
                self.routines.cells.removeCell(hyds[1])
                residue.removeAtom(hyds[1].name)
                donorhatom = residue.getAtom(hyds[0].name)

            # If only one H is left, we're done

            if len(self.hlist) == 1:
                self.rename(donorhatom)
                residue.fixed = 1
            return 1
        
        else: return 0

     
    def tryDonor(self, donor, acc):
        """
           The main driver for adding a hydrogen to an
           optimizeable residue.
        """
        residue = self.residue

        # Do some error checking
        
        if not acc.hacceptor: return 0
 
        if self.isHbond(donor, acc):
            self.fix(donor, acc)
            return 1
        else: return 0

    def fix(self, donor, acc):
        """
            Fix the carboxylic residue.
        """

        self.debug("Fixing residue %s due to %s" % (donor.residue, donor.name))

        residue = donor.residue
        
        # Grab the H(D) that caused the bond

        for donorhatom in donor.bonds:
            if donorhatom.isHydrogen() and \
               self.getHbondangle(acc, donor, donorhatom) <= ANGLE_CUTOFF: break

        # Remove all the other available bonded hydrogens
        
        hydrogens = self.hlist[:]
        for atom in hydrogens:
            if atom != donorhatom:
                self.routines.cells.removeCell(atom)
                self.hlist.remove(atom)
                residue.removeAtom(atom.name)

        #for atom in self.atomlist:
        #    print "testing atom", atom
        #    for bond in atom.bonds:
        #        print "testing bond", bond
        #        if bond.isHydrogen() and bond != donorhatom:
        #            self.routines.cells.removeCell(bond)
        #            self.hlist.remove(bond)
        #            residue.removeAtom(bond.name)

        # Rename the atoms

        self.rename(donorhatom)

        residue.fixed = 1
        
    def finalize(self):
        """
            Finalize a protontated residue.  Try to minimize
            conflict with nearby atoms.
        """
        
        # Initialize some variables

        hydrogens = []
        bestdist = 0.0
        bestatom = None
        residue = self.residue

        if residue.fixed: return

        # For each atom, get the closest atom

        for hydatom in self.hlist:
            closeatom = self.routines.getClosestAtom(hydatom)
            dist = distance(hydatom.getCoords(), closeatom.getCoords())
            if dist > bestdist:
                bestdist = dist
                bestatom = hydatom

        # Keep the bestatom
        
        for hydatom in self.hlist:
            hydrogens.append(hydatom)

        for hydatom in hydrogens:
            if bestatom != hydatom:
                self.hlist.remove(hydatom)
                residue.removeAtom(hydatom.name)

        # Rename the atoms

        self.rename(bestatom)
        residue.fixed = 1

    def rename(self, hydatom):
        """
             Rename the optimized atoms appropriately.  This is done
             since the forcefields tend to require that the hydrogen is
             linked to a specific oxygen, and this atom may have different
             parameter values.

             Parameters
                 hydatom:  The hydrogen atom that was added. (atom)
        """
        residue = self.residue
        
        # Take off the extension

        hname = hydatom.name[:-1]
        residue.renameAtom(hydatom.name, hname)
        
        # PATCHES.xml expects *2 - if it's *1 that left, flip names
        
        if hydatom.name.endswith("1"):
            residue.renameAtom(hydatom.name,"%s2" %  hydatom.name[:-1])
            bondname0 = self.atomlist[0].name
            bondname1 = self.atomlist[1].name
            tempname = "FLIP"
            residue.renameAtom(self.atomlist[0].name, tempname)
            residue.renameAtom(self.atomlist[1].name, bondname0)
            residue.renameAtom(tempname, bondname1)
        elif hydatom.name.endswith("OXT"):
            residue.renameAtom(hydatom.name,"HO")
            bondname0 = self.atomlist[0].name
            bondname1 = self.atomlist[1].name
            tempname = "FLIP"
            residue.renameAtom(self.atomlist[0].name, tempname)
            residue.renameAtom(self.atomlist[1].name, bondname0)
            residue.renameAtom(tempname, bondname1)
   
   
    def complete(self):
        """
            If not already fixed, finalize
        """
        if not self.residue.fixed: self.finalize()
        
class hydrogenRoutines:
    """
        The main routines for hydrogen optimization.  This could
        potentially be extended from the routines object...
    """
    def __init__(self, routines):
        """
            Parse the XML file and store the data in a map
        """
        self.routines = routines
        self.protein = routines.protein
        self.optlist = []
        self.atomlist = []
        self.resmap = {}
    
        handler = HydrogenHandler()
        sax.make_parser()

        defpath = getDatFile(HYDPATH)
        if defpath == "":
            raise ValueError, "Could not find %s!" % HYDPATH 
     
        file = open(defpath)
        sax.parseString(file.read(), handler)
        file.close()

        self.map = handler.map

    def debug(self, text):
        """
            Print text to stdout for debugging purposes.

            Parameters
                text:  The text to output (string)
        """
        if HDEBUG: print text  

    def isOptimizeable(self, residue):
        """
            Check to see if the given residue is optimizeable
            There are three ways to identify a residue:

            1.  By name (i.e. HIS)
            2.  By reference name - a PDB file HSP has
                a HIS reference name
            3.  By patch - applied by PropKa, terminal selection

            Parameters
                residue:  The residue in question (Residue)
            Returns
                optinstance: None if not optimizeable, otherwise
                             the OptimizationHolder instance that
                             corresponds to the residue.
        """
        optinstance = None
        if not (isinstance(residue, Amino) or isinstance(residue, WAT)):
            return optinstance
        
        if residue.name in self.map: 
            optinstance = self.map[residue.name]
        elif residue.reference.name in self.map:
            optinstance = self.map[residue.reference.name]
        else:
            for patch in residue.patches:
                if patch in self.map:
                    optinstance = self.map[patch]
                    break

        # If alcoholic, make sure the hydrogen is present

        if optinstance != None:
            if optinstance.opttype == "Alcoholic":
                atomname = optinstance.map.keys()[0]
                if not residue.reference.hasAtom(atomname):
                    optinstance = None
                   
        return optinstance
                
    def setOptimizeableHydrogens(self):
        """
            Set any hydrogen listed in HYDROGENS.xml that
            is optimizeable.  Used BEFORE hydrogen optimization
            to label atoms so that they won't be debumped - i.e.
            if SER HG is too close to another atom, don't debump
            but wait for optimization.  This function should not
            be used if full optimization is not taking place.
        """
        for residue in self.protein.getResidues():
            optinstance = self.isOptimizeable(residue)
            if optinstance == None: continue
            for atom in residue.getAtoms():
                if atom.name in optinstance.map:
                    atom.optimizeable = 1

    def initializeFullOptimization(self):
        """
            Initialize the full optimization.  Detects all
            optimizeable donors and acceptors and sets the internal
            optlist.
        """
        self.routines.write("Initializing full optimization...\n")
        
        # Do some setup

        self.routines.cells = Cells(5)
        self.routines.cells.assignCells(self.protein)
        self.routines.calculateDihedralAngles()
        self.routines.setDonorsAndAcceptors()
        self.routines.updateInternalBonds()
        self.routines.setReferenceDistance()
        
        # First initialize the various types
   
        for residue in self.protein.getResidues():
            optinstance = self.isOptimizeable(residue)
            if optinstance == None: continue

            type = optinstance.opttype
            command = "%s(residue, optinstance, self.routines)" % type
            myobj = eval(command)
            self.atomlist += myobj.atomlist
            self.optlist.append(myobj)
            self.resmap[residue] = myobj
      
        self.routines.write("Done.\n")

    def initializeWaterOptimization(self):
        """
            Initialize optimization for waters only.  Detects all
            optimizeable donors and acceptors and sets the internal
            optlist.
        """

        self.routines.write("Initializing water bonding optimization...\n")
        
        # Do some setup

        self.routines.cells = Cells(5)
        self.routines.cells.assignCells(self.protein)
        self.routines.calculateDihedralAngles()
        self.routines.setDonorsAndAcceptors()
        self.routines.updateInternalBonds()
        self.routines.setReferenceDistance()
        
        # First initialize the various types
   
        for residue in self.protein.getResidues():
            optinstance = self.isOptimizeable(residue)
            if optinstance == None: continue

            type = optinstance.opttype
            if type == "Water":
                command = "%s(residue, optinstance, self.routines)" % type
                myobj = eval(command)
                self.atomlist += myobj.atomlist
                self.optlist.append(myobj)
                self.resmap[residue] = myobj
              
        self.routines.write("Done.\n")
    
    def optimizeHydrogens(self):
        """
            The main driver for the optimization.  Should be
            called only after the optlist has been initialized.
        """

        self.routines.write("Optimization progress:\n")
        
        optlist = self.optlist
        resmap = {}
        connectivity = {}

        # Initialize the detection progress

        if len(optlist) == 0: return

   
        self.routines.write("  Detecting potential hydrogen bonds:\n")
        self.routines.write("0% |                    | 100%\n", 1)
        self.routines.write("    ", 1)
        progress = 0.0
        increment = 1.0/len(optlist) 

        for obj in optlist:
            connectivity[obj] = []
            for atom in obj.atomlist:
                closeatoms = self.routines.cells.getNearCells(atom)
                for closeatom in closeatoms:
   
                    # Conditions for continuing
                
                    if atom.residue == closeatom.residue: continue
                    if not (closeatom.hacceptor or closeatom.hdonor): continue
                    if atom.hdonor and not atom.hacceptor \
                       and not closeatom.hacceptor: continue
                    if atom.hacceptor and not atom.hdonor \
                       and not closeatom.hdonor: continue
                
                    dist = distance(atom.getCoords(), closeatom.getCoords())
                    if dist < 4.3:
                        residue = atom.residue
                        hbond = PotentialBond(atom, closeatom, dist)

                        # Store the potential bond

                        obj.hbonds.append(hbond)

                        # Keep track of connectivity
                    
                        if closeatom in self.atomlist:
                            closeobj = self.resmap[closeatom.residue]
                            if closeobj not in connectivity[obj]:
                                connectivity[obj].append(closeobj)

            progress += increment
            while progress >= 0.0499:
                self.routines.write("*")
                progress -= 0.05

        if len(optlist) > 0: self.routines.write("\n")                   

        # Some residues might have no nearby hbonds - if so, place at
        #   default state

        for obj in optlist:
            if len(obj.hbonds) == 0:
                if obj.residue.fixed: continue
                self.debug("%s has no nearby partners - fixing." % obj.residue)
                obj.finalize()
            
        # Determine the distinct networks

        networks = []
        seen = []
        for obj in optlist:
            if obj.residue.fixed: continue
            if obj in seen: continue
            network = analyzeConnectivity(connectivity, obj)
            for obj in network:
                if obj not in seen: seen.append(obj)
             
            networks.append(network)

        # Initialize the output progress

        if len(networks) > 0:
            self.routines.write("  Optimizing hydrogen bonds:\n")
            self.routines.write("0% |                    | 100%\n", 1)
            self.routines.write("    ", 1)
            progress = 0.0
            increment = 1.0/len(networks)
            
        # Work on the networks

        for network in networks:
            txt = ""
            for obj in network:
                txt += "%s, " % obj
            self.debug("\nStarting network %s" % txt[:-2])

            ###  FIRST:  Only optimizeable to backbone atoms

            self.debug("* Optimizeable to backbone *")

            hbondmap = {}
            for obj in network:
                for hbond in obj.hbonds:
                    if hbond.atom2 not in self.atomlist:
                        hbondmap[hbond] = hbond.dist

            hbondlist = sortDictByValue(hbondmap)
            hbondlist.reverse()

            for hbond in hbondlist:
                atom = hbond.atom1
                atom2 = hbond.atom2
                obj = self.resmap[atom.residue]

                if atom.residue.fixed: continue
                if atom.hdonor: obj.tryDonor(atom, atom2)
                if atom.hacceptor: obj.tryAcceptor(atom, atom2)
              
            ### SECOND:  Non-dual water Optimizeable to Optimizeable

            self.debug("\n* Optimizeable to optimizeable *")

            hbondmap = {}
            seenlist = []
            for obj in network:
                for hbond in obj.hbonds:
                    if hbond.atom2 in self.atomlist \
                           and not (isinstance(hbond.atom1.residue, WAT) \
                           and isinstance(hbond.atom2.residue, WAT)):

                        # Only get one hbond pair
                     
                        if not (hbond.atom2, hbond.atom1) in seenlist:
                            hbondmap[hbond] = hbond.dist
                            seenlist.append((hbond.atom1, hbond.atom2))

            hbondlist = sortDictByValue(hbondmap)
            hbondlist.reverse()

            for hbond in hbondlist:
                atom = hbond.atom1
                atom2 = hbond.atom2
                obj1 = self.resmap[atom.residue]
                obj2 = self.resmap[atom2.residue]
                
                # Atoms may no longer exist if already optimized

                if not atom.residue.hasAtom(atom.name): continue
                if not atom2.residue.hasAtom(atom2.name): continue
            
                res = 0
                if atom.hdonor and atom2.hacceptor:
                    res = obj1.tryBoth(atom, atom2, obj2)
                  
                if atom.hacceptor and atom2.hdonor and res == 0:
                    obj2.tryBoth(atom2, atom, obj1)

            ### THIRD:  All water-water residues

            self.debug("\n* Water to Water *")

            hbondmap = {}
            seenlist = []
            for obj in network:
                for hbond in obj.hbonds:
                    residue = hbond.atom1.residue
                    if isinstance(residue, WAT) and \
                       isinstance(hbond.atom2.residue, WAT):
                        if not (hbond.atom2, hbond.atom1) in seenlist:
                            hbondmap[hbond] = hbond.dist
                            seenlist.append((hbond.atom1, hbond.atom2))
              
            
            hbondlist = sortDictByValue(hbondmap)
            hbondlist.reverse()

            for hbond in hbondlist:
                atom = hbond.atom1
                atom2 = hbond.atom2
                obj1 = self.resmap[atom.residue]
                obj2 = self.resmap[atom2.residue]
                
                res = 0
                if atom.hdonor and atom2.hacceptor:
                    res = obj1.tryBoth(atom, atom2, obj2)
                  
                if atom.hacceptor and atom2.hdonor and res == 0:
                    obj2.tryBoth(atom2, atom, obj1)
         
                   
            ### FOURTH: Complete all residues
    
            for obj in network:  obj.complete()
    
            # STEP 5:  Update progress meter

            progress += 100.0 * increment
            while progress >= 5.0:
                self.routines.write("*")
                progress -= 5.0

        if len(networks) > 0: self.routines.write("\n")

class OptimizationHolder:
    """
        A holder class for the XML parser.
    """
    def __init__(self):
        """
            Initialize the class.
        """
        self.name = ""
        self.map = {}
        self.opttype = ""
        self.optangle = ""

    def __str__(self):
        """
            A basic string representation for debugging
        """
        text = "%s\n" % self.name
        text += "Type: %s\n" % self.opttype
        if self.optangle != "":
            text += "Optimization Angle: %s\n" % self.optangle
        text += "Atoms: \n"
        for atomname in self.map:
            text += "\t%s\n" % str(self.map[atomname])
        return text
    
