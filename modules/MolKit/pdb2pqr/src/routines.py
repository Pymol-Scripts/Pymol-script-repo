"""
    Routines for PDB2PQR

    This module contains the protein object used in PDB2PQR and methods
    used to correct, analyze, and optimize that protein.
   
    ----------------------------
   
    PDB2PQR -- An automated pipeline for the setup, execution, and analysis of
    Poisson-Boltzmann electrostatics calculations

    Nathan A. Baker (baker@biochem.wustl.edu)
    Todd Dolinsky (todd@ccb.wustl.edu)
    Dept. of Biochemistry and Molecular Biophysics
    Center for Computational Biology
    Washington University in St. Louis

    Jens Nielsen (Jens.Nielsen@ucd.ie)
    University College Dublin

    Additional contributing authors listed in documentation and supporting
    package licenses.

    Copyright (c) 2003-2007.  Washington University in St. Louis.  
    All Rights Reserved.

    This file is part of PDB2PQR.

    PDB2PQR is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.
 
    PDB2PQR is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
   
    You should have received a copy of the GNU General Public License
    along with PDB2PQR; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA

    ----------------------------

"""

__date__ = "28 February 2006"
__author__ = "Jens Erik Nielsen, Todd Dolinsky"

CELL_SIZE = 2
BUMP_DIST = 2.0
BUMP_HDIST = 1.5
BONDED_SS_LIMIT = 2.5
PEPTIDE_DIST = 1.7
REPAIR_LIMIT = 10

import math
import copy
from pdb import *
from utilities import *
from quatfit import *
from forcefield import *
from structures import *
from protein import *
from definitions import *

class Routines:
    def __init__(self, protein, verbose):
        """
            Initialize the Routines class.  The class contains most
            of the main routines that run PDB2PQR

            Parameters
                protein:  The protein to run PDB2PQR on (Protein)
                verbose:  A flag to determine whether to write to
                          stdout
        """
        self.protein = protein
        self.verbose = verbose
        self.warnings = []
        self.cells = {}

            
    def write(self, message, indent=0):
        """
            Write a message to stderr for debugging if verbose

            Parameters
                message: The message to write (string)
                indent : The indent level (int, default=0)
        """
        out = ""
        if self.verbose:
           for i in range(indent):
               out += "\t"
           out += message
           sys.stderr.write(out)

    def getWarnings(self):
        """
            Get all warnings generated from routines
        """
        return self.warnings

    def applyNameScheme(self, forcefield):
        """
            Apply the naming scheme of the give forcefield to the atoms
            within the protein

            Parameters
                forcefield: The forcefield object (forcefield)
           
        """
        self.write("Applying the naming scheme to the protein...")
        for residue in self.protein.getResidues():
            if isinstance(residue, Amino) or isinstance(residue, WAT) \
               or isinstance(residue, Nucleic):
               resname = residue.ffname
            else: resname = residue.name

            for atom in residue.getAtoms():
                rname, aname = forcefield.getNames(resname, atom.name)
                if aname != None and rname != None:
                    atom.resName = rname
                    atom.name = aname
                    
        self.write("Done.\n")            
           
    def applyForcefield(self, forcefield):
        """
            Apply the forcefield to the atoms within the protein

            Parameters
                forcefield: The forcefield object (forcefield)
            Returns
                hitlist:    A list of atoms that were found in
                            the forcefield (list)
                misslist:   A list of atoms that were not found in
                            the forcefield (list)
        """
        self.write("Applying the forcefield to the protein...")
        misslist = []
        hitlist = []
        for residue in self.protein.getResidues():
            if isinstance(residue, Amino) or isinstance(residue, WAT) \
               or isinstance(residue, Nucleic):
               resname = residue.ffname
            else: resname = residue.name

            # Apply the parameters
           
            for atom in residue.getAtoms():
                atomname = atom.get("name")
                charge, radius = forcefield.getParams(resname, atomname)
                if charge != None and radius != None:
                    atom.set("ffcharge", charge)
                    atom.set("radius", radius)
                    hitlist.append(atom)
                else:
                    misslist.append(atom)
                    
        self.write("Done.\n")            
        return hitlist, misslist
   
            
    def updateSSbridges(self):
        """
            Check for SS-bridge partners, and if present, set appropriate
            partners
        """
        self.write("Updating SS bridges...\n")
        SGpartners = {}
        for residue in self.protein.getResidues():
            if isinstance(residue, CYS):
                atom = residue.getAtom("SG")
                if atom != None:
                    SGpartners[atom] = []

        for atom in SGpartners:
            for partner in SGpartners:
                if atom == partner or SGpartners[atom] != []: continue
                dist = distance(atom.getCoords(), partner.getCoords())
                if dist < BONDED_SS_LIMIT:
                    SGpartners[atom].append(partner)
                    SGpartners[partner].append(atom)
        
        for atom in SGpartners:
            res1 = atom.get("residue")
            numpartners = len(SGpartners[atom])
            if numpartners == 1:
                partner = SGpartners[atom][0]
                res2 = partner.get("residue")                  
                res1.set("SSbonded", 1)
                res1.set("SSbondedpartner", partner)
                self.applyPatch("CYX", res1)
                self.write("%s - %s\n" % (res1, res2), 1)
            elif numpartners > 1:
                error = "WARNING: %s has multiple potential " % res1
                error += "SS-bridge partners\n"
                self.write(error, 1)
                self.warnings.append(error)
            elif numpartners == 0:
                self.write("%s is a free cysteine\n" % res1, 1)
        self.write("Done.\n")

    def updateInternalBonds(self):
        """
            Update the internal bonding network using the reference
            objects in each atom.
        """
        for residue in self.protein.getResidues():
            if isinstance(residue, Amino) or isinstance(residue, WAT) or \
               isinstance(residue, Nucleic):
                for atom in residue.getAtoms():
                    if not atom.hasReference(): continue
                    for bond in atom.reference.bonds:
                        if not residue.hasAtom(bond): continue
                        bondatom = residue.getAtom(bond)
                        if bondatom not in atom.bonds:
                            atom.addBond(bondatom)
            
    def updateBonds(self):
        """
            Update the bonding network of the protein.  This happens
            in 3 steps:
              1.  Applying the PEPTIDE patch to all Amino residues
                  so as to add reference for the N(i+1) and C(i-1)
                  atoms
              2.  UpdateInternalBonds for inter-residue linking
              3.  Set the links to the N(i+1) and C(i-1) atoms
        """

        # Apply the peptide patch
       
        for residue in self.protein.getResidues():
            if isinstance(residue, Amino):
                if residue.isNterm or residue.isCterm: continue
                else: self.applyPatch("PEPTIDE", residue)

        # Update all internal bonds
                
        self.updateInternalBonds()

        # Set the peptide bond pointers

        for chain in self.protein.getChains():
            for i in range(chain.numResidues() - 1):
                res1 = chain.residues[i]
                res2 = chain.residues[i+1]
                if not isinstance(res1, Amino) or not isinstance(res2, Amino):
                    continue
                atom1 = res1.getAtom("C")
                atom2 = res2.getAtom("N")
                if atom1 != None: res2.peptideC = atom1
                if atom2 != None: res1.peptideN = atom2
                if atom1 == None or atom2 == None: continue
                if distance(atom1.getCoords(), atom2.getCoords()) > PEPTIDE_DIST:
                    text = "Gap in backbone detected between %s and %s!\n" % \
                           (res1, res2)
                    self.write(text, 1)
                    self.warnings.append(text)

    def applyPatch(self, patchname, residue):
        """
            Apply a patch to the given residue.  This is one of the key
            functions in PDB2PQR.  A similar function appears in
            definitions.py - that version is needed for residue level
            subtitutions so certain protonation states (i.e. CYM, HSE)
            are detectatble on input.

            This version looks up the particular patch name in the
            patchmap stored in the protein, and then applies the
            various commands to the reference and actual residue
            structures.

            See the inline comments for a more detailed explanation.

            Parameters
                patchname:  The name of the patch (string)
                residue:    The residue to apply the patch to (residue)
        """
        if patchname not in self.protein.patchmap:
            raise ValueError,"Unable to find patch %s!" % patchname
        
        # Make a copy of the reference, i.e. a new reference for
        # this patch.  Two examples:
        #     PEPTIDE is a special case, as it applies to
        #             every residue.
        #     CTERM only applies to one specific residue, so a
        #             deep copy is used.

        if patchname == "PEPTIDE":
            newreference = residue.reference
        else:
            newreference = copy.deepcopy(residue.reference)
            
        patch = self.protein.patchmap[patchname]

        # Add atoms from patch

        for atomname in patch.map: 
            newreference.map[atomname] = patch.map[atomname]
            for bond in patch.map[atomname].bonds:
                if bond not in newreference.map: continue
                if atomname not in newreference.map[bond].bonds:
                    newreference.map[bond].bonds.append(atomname)
                    
        # Remove atoms as directed by patch
            
        for remove in patch.remove:
            if remove in residue.map: residue.removeAtom(remove)
            if remove not in newreference.map: continue
            removebonds = newreference.map[remove].bonds
            del newreference.map[remove]
            for bond in removebonds:
                index = newreference.map[bond].bonds.index(remove)
                del newreference.map[bond].bonds[index]

        # Add the new dihedrals

        for dihedral in patch.dihedrals:
            newreference.dihedrals.append(dihedral)        

        # Point at the new reference

        residue.reference = newreference
        residue.patches.append(patchname)

        # Rename atoms as directed by patch

        for atom in residue.getAtoms():
            if atom.name in patch.altnames:
                residue.renameAtom(atom.name, patch.altnames[atom.name])

        # Replace each atom's reference with the new one

        for atomname in residue.map:
            if newreference.hasAtom(atomname):
                atom = residue.getAtom(atomname)
                atom.reference = newreference.map[atomname]

    def setStates(self):
        """
            Set the state of each residue.  This is the last step
            before assigning the forcefield, but is necessary so
            as to distinguish between various protonation states.

            See aa.py for residue-specific functions.
        """
        for residue in self.protein.getResidues():
            if isinstance(residue, Amino) or \
               isinstance(residue, Nucleic):
                residue.setState()

    def assignTermini(self, chain):
        """
            Assign the termini for the given chain by looking at
            the start and end residues.
        """

	if len(chain.residues) == 0: 
	    text = "Error: chain \"%s\" has 0 residues!" % chain.chainID
	    raise ValueError, text

        # Set the N-Terminus/ 5' Terminus

        res0 = chain.residues[0]
        if isinstance(res0, Amino):
            res0.set("isNterm",1)
            if isinstance(res0, PRO):
                self.applyPatch("NEUTRAL-NTERM", res0)
            else:
                self.applyPatch("NTERM",res0)
        elif isinstance(res0, Nucleic):
            res0.set("is5term",1)
            self.applyPatch("5TERM", res0)
                    
        # Set the C-Terminus/ 3' Terminus
            
        reslast = chain.residues[-1]
        if isinstance(reslast, Amino):
            reslast.set("isCterm",1)
            self.applyPatch("CTERM", reslast)
        elif isinstance(reslast, Nucleic):
            reslast.set("is3term",1)
            self.applyPatch("3TERM", reslast)
        else:
            for i in range(len(chain.residues)):
                resthis = chain.residues[-1 - i]
                if isinstance(resthis, Amino):
                    resthis.set("isCterm",1)
                    self.applyPatch("CTERM", resthis)
                    break
                elif resthis.name in ["NH2","NME"]: break
                elif isinstance(resthis, Nucleic):
                    resthis.set("is3term",1)
                    self.applyPatch("3TERM", resthis)
                    break
                    
    def setTermini(self):
        """
            Set the termini for the protein. First set all known
            termini by looking at the ends of the chain. Then
            examine each residue, looking for internal chain breaks.
        """
        self.write("Setting the termini... \n")

        # First assign the known termini

        for chain in self.protein.getChains():
            self.assignTermini(chain)

        # Now determine if there are any hidden chains

        letters = string.ascii_uppercase + string.ascii_lowercase
        numchains = len(self.protein.getChains())
        c = 0
        while c < len(self.protein.getChains()):

            chain = self.protein.chains[c]

            reslist = []
	    origlist = []	

	    # origlist holds the original residue list for the chain
	
	    for residue in chain.getResidues():
	        origlist.append(residue)
            
            for residue in origlist:
	        reslist.append(residue)
                oldid = residue.chainID
                
                # Look for ending termini
   		
		fixflag = 0
                if isinstance(residue, Amino):
		    if (residue.hasAtom("OXT") and not residue.isCterm):
			fixflag = 1

		elif isinstance(residue, Nucleic):
		    if ((residue.hasAtom("H3T") or residue.name.endswith("3"))\
		      and not residue.is3term):
			fixflag = 1

		if fixflag:

                    # Get an available chain ID
                    
                    numchains = len(self.protein.getChains())
                    chainid = letters[0]
                    id = 0
                    while chainid in self.protein.chainmap:
                        id += 1
                        chainid = letters[id]

                    # Make a new chain with these residues

                    newchain = Chain(chainid)
		  
                    self.protein.chainmap[chainid] = newchain
                    self.protein.chains.insert(c, newchain)		  

                    for res in reslist:
                        newchain.addResidue(res)
		        chain.residues.remove(res)
                        res.setChainID(chainid)
		
		
                    self.assignTermini(chain)
                    self.assignTermini(newchain)
                    
		    reslist = []
		    c += 1
		   
            c += 1

        # Update the final chain's chainID if it is "" unless it's all water

        if "" in self.protein.chainmap:

	    notwat = 0
	    for res in chain.residues:
		if not isinstance(res, WAT): 
		    notwat = 1
	            break

	    if notwat == 0: 
		self.write("Done.\n")
		return		

            chain = self.protein.chainmap[""]
            chainid = letters[0]
            id = 0
            while chainid in self.protein.chainmap:
                id += 1
                chainid = letters[id]

            # Use the new chainID

            self.protein.chainmap[chainid] = chain
            del self.protein.chainmap[""]

            for res in chain.residues:
                res.setChainID(chainid)
            
        self.write("Done.\n")

    def findMissingHeavy(self):
        """
            Repair residues that contain missing heavy (non-Hydrogen) atoms
        """
        self.write("Checking for missing heavy atoms... \n")
        misscount = 0
        heavycount = 0
        for residue in self.protein.getResidues():
            if not (isinstance(residue, Amino) or \
                    isinstance(residue, Nucleic)): continue
                      
            # Check for Missing Heavy Atoms

            for refatomname in residue.reference.map:
                if refatomname.startswith("H"): continue
                if refatomname in ["N+1","C-1"]: continue
                heavycount += 1
                if not residue.hasAtom(refatomname):
                    self.write("Missing %s in %s\n" % \
                               (refatomname, residue), 1)
                    misscount += 1
                    residue.addMissing(refatomname)                    
                            
            # Check for Extra Atoms

            atomlist = []
            for atom in residue.get("atoms"):
                atomlist.append(atom)

            for atom in atomlist:
                atomname = atom.get("name")                 
                if not residue.reference.hasAtom(atomname):
                    self.write("Extra atom %s in %s! - " % \
                               (atomname, residue), 1)
                    residue.removeAtom(atomname)
                    self.write("Deleted this atom.\n")

        if heavycount == 0:
            raise ValueError, "No heavy atoms found!"
                            
        misspct = 100.0 * float(misscount) / heavycount 
        if misspct > REPAIR_LIMIT:
            error = "This PDB file is missing too many (%i out of " % misscount
            error += "%i, %.2f%%) heavy atoms to accurately repair the file.  " % \
                     (heavycount, misspct)
            error += "The current repair limit is set at %i%%." % REPAIR_LIMIT
            raise ValueError, error
        elif misscount > 0:
            self.write("Missing %i out of %i heavy atoms (%.2f percent) - " %\
                       (misscount, heavycount, misspct))
            self.write("Will attempt to repair.\n")
            self.repairHeavy()
        else:
            self.write("No heavy atoms found missing - Done.\n")

    def rebuildTetrahedral(self, residue, atomname):
        """
            Rebuild a tetrahedral hydrogen group.  This is necessary
            due to the shortcomings of the quatfit routine - given a
            tetrahedral geometry and two existing hydrogens, the
            quatfit routines have two potential solutions.  This function
            uses basic tetrahedral geometry to fix this issue.

            Parameters
                residue:  The residue in question (residue)
                atomname: The atomname to add (string)
            Returns
                1 if successful, 0 otherwise
        """

        hcount = 0
        nextatomname = None
    
        atomref = residue.reference.map[atomname]
        bondname = atomref.bonds[0]

        # Return if the bonded atom does not exist
        
        if not residue.hasAtom(bondname): return 0

        # This group is tetrahedral if bondatom has 4 bonds,
        #  3 of which are hydrogens
        
        for bond in residue.reference.map[bondname].bonds:
            if bond.startswith("H"): hcount += 1
            else: nextatomname = bond

        # Check if this is a tetrahedral group

        if hcount != 3 or nextatomname == None: return 0

        # Now rebuild according to the tetrahedral geometry

        bondatom = residue.getAtom(bondname)
        nextatom = residue.getAtom(nextatomname)
        numbonds = len(bondatom.bonds)

        if numbonds == 1:

            # Place according to two atoms
            
            coords = [bondatom.getCoords(), nextatom.getCoords()]
            refcoords = [residue.reference.map[bondname].getCoords(), \
                         residue.reference.map[nextatomname].getCoords()]
            refatomcoords = atomref.getCoords()
            newcoords = findCoordinates(2, coords, refcoords, refatomcoords)
            residue.createAtom(atomname, newcoords)
        
            return 1
            
        elif numbonds == 2:

            # Get the single hydrogen coordinates

            hatom = None
            for bond in bondatom.reference.bonds:
                if residue.hasAtom(bond) and bond.startswith("H"):
                    hatom = residue.getAtom(bond)
                    break

            # Use the existing hydrogen and rotate about the bond

            residue.rotateTetrahedral(nextatom, bondatom, 120)
            newcoords = hatom.getCoords()
            residue.rotateTetrahedral(nextatom, bondatom, -120)
            residue.createAtom(atomname, newcoords)
            
            return 1
        
        elif numbonds == 3:

            # Find the one spot the atom can be

            hatoms = []
            for bond in bondatom.reference.bonds:
                if residue.hasAtom(bond) and bond.startswith("H"):
                    hatoms.append(residue.getAtom(bond))

            # If this is more than two something is wrong

            if len(hatoms) != 2: return 0

            # Use the existing hydrogen and rotate about the bond

            residue.rotateTetrahedral(nextatom, bondatom, 120)
            newcoords1 = hatoms[0].getCoords()
            residue.rotateTetrahedral(nextatom, bondatom, 120)
            newcoords2 = hatoms[0].getCoords()
            residue.rotateTetrahedral(nextatom, bondatom, 120)

            # Determine which one hatoms[1] is not in

            if distance(hatoms[1].getCoords(), newcoords1) > 0.1:
                residue.createAtom(atomname, newcoords1)
            else:
                residue.createAtom(atomname, newcoords2)
          
            return 1
            
    def addHydrogens(self):
        """
            Add the hydrogens to the protein.  This requires either
            the rebuildTetrahedral function for tetrahedral geometries
            or the standard quatfit methods.  These methods use three
            nearby bonds to rebuild the atom; the closer the bonds, the
            more accurate the results.  As such the peptide bonds are
            used when available.
        """
        count = 0
        self.write("Adding hydrogens to the protein...\n")
        for residue in self.protein.getResidues():
            if not (isinstance(residue, Amino) or \
                    isinstance(residue, Nucleic)): continue
            for atomname in residue.reference.map:
                if not atomname.startswith("H"): continue
                if residue.hasAtom(atomname): continue
                if isinstance(residue,CYS):
                    if residue.SSbonded and atomname == "HG": continue

                # If this hydrogen is part of a tetrahedral group,
                #  follow a different codepath

                if self.rebuildTetrahedral(residue, atomname):
                    count += 1
                    continue
         
                # Otherwise use the standard quatfit methods
                 
                coords = []
                refcoords = []
                
                refatomcoords = residue.reference.map[atomname].getCoords()
                bondlist = residue.reference.getNearestBonds(atomname)
                
                for bond in bondlist:
                    if bond == "N+1": atom = residue.peptideN
                    elif bond == "C-1": atom = residue.peptideC
                    else: atom = residue.getAtom(bond)

                    if atom == None: continue

                    # Get coordinates, reference coordinates

                    coords.append(atom.getCoords())
                    refcoords.append(residue.reference.map[bond].getCoords())

                    # Exit if we have enough atoms
                    
                    if len(coords) == 3: break

                if len(coords) == 3:
                    newcoords = findCoordinates(3, coords, refcoords, refatomcoords)
                    residue.createAtom(atomname, newcoords)
                    count += 1
                else:
                    self.write("Couldn't rebuild %s in %s!\n" % (atomname, residue),1)
                    
        self.write(" Added %i hydrogen atoms.\n" % count)

    def repairHeavy(self):
        """
            Repair all heavy atoms.  Unfortunately the first time we
            get to an atom we might not be able to rebuild it - it
            might depend on other atoms to be rebuild first (think side
            chains).  As such a 'seenmap' is used to keep track of what
            we've already seen and subsequent attempts to rebuild the
            atom.
        """
        self.write("Rebuilding missing heavy atoms... \n")
        for residue in self.protein.getResidues():
            if not (isinstance(residue, Amino) or \
                    isinstance(residue, Nucleic)): continue
            missing = residue.get("missing")
            if missing == []: continue

            # Initialize some variables

            seenmap = {}
            nummissing = len(missing)
              
            while len(missing) > 0:
                coords = []
                refcoords = []
                
                atomname = missing.pop(0)
                refatomcoords = residue.reference.map[atomname].getCoords()
                bondlist = residue.reference.getNearestBonds(atomname)

                for bond in bondlist:
                    if bond == "N+1": atom = residue.peptideN
                    elif bond == "C-1": atom = residue.peptideC
                    else: atom = residue.getAtom(bond)

                    if atom == None: continue

                    # Get coordinates, reference coordinates

                    coords.append(atom.getCoords())
                    refcoords.append(residue.reference.map[bond].getCoords())

                    # Exit if we have enough atoms
                    
                    if len(coords) == 3: break

                # We might need other atoms to be rebuilt first

                if len(coords) < 3:
                    try: seenmap[atomname] += 1
                    except KeyError: seenmap[atomname] = 1
                    if seenmap[atomname] > nummissing:
                        text = "Unable to rebuild atom %s in %s!" % \
                               (atomname, residue)
                        raise ValueError, text
                    else: missing.append(atomname)

                else: # Rebuild the atom
                    newcoords = findCoordinates(3, coords, refcoords, refatomcoords)
                    residue.createAtom(atomname, newcoords)
                    self.write("Added %s to %s at coordinates" % (atomname, residue),1)
                    self.write(" %.3f %.3f %.3f\n" % \
                           (newcoords[0], newcoords[1], newcoords[2]))
                    
        self.write("Done.\n")

    def setReferenceDistance(self):
        """
            Set the distance to the CA atom in the residue.
            This is necessary for determining which atoms are
            allowed to move during rotations.  Uses the
            shortestPath algorithm found in utilities.py.
        """
        for residue in self.protein.getResidues():
            if not isinstance(residue, Amino): continue

            # Initialize some variables
            
            map = {}
            caatom = residue.getAtom("CA")

            if caatom == None:
                text = "Cannot set references to %s without CA atom!\n"
                raise ValueError, text
            
            # Set up the linked map
            
            for atom in residue.getAtoms(): 
                map[atom] = atom.bonds
             
            # Run the algorithm
                
            for atom in residue.getAtoms():
                if atom.isBackbone():
                    atom.refdistance = -1 
                else:
                    atom.refdistance = len(shortestPath(map, atom, caatom)) - 1
        
    def debumpProtein(self):
        """
            Make sure that none of the added atoms were rebuilt
            on top of existing atoms.  See each called function
            for more information.
        """
        
        self.write("Checking if we must debump any residues... \n")

        # Do some setup

        self.cells = Cells(CELL_SIZE)
        self.cells.assignCells(self.protein)

        self.calculateDihedralAngles()
        self.setDonorsAndAcceptors()
        self.updateInternalBonds()
        self.setReferenceDistance()
    
        # Determine which residues to debump

        for residue in self.protein.getResidues():
            if not isinstance(residue, Amino): continue

            # Initialize variables

            conflictnames = []

            for atom in residue.getAtoms():
                atomname = atom.name
                if not atom.added: continue
                if atomname == "H": continue
                if atom.optimizeable: continue
               
                nearatoms = self.findNearbyAtoms(atom)
                
                # If something is too close, we must debump the residue
        
                if nearatoms != {}:
                    conflictnames.append(atomname)
                    for repatom in nearatoms:
                        self.write("%s %s is too close to %s %s\n" % \
                                  (residue, atomname, repatom.residue, repatom.name),1) 

            # If there are no conflicting atoms, move on

            if conflictnames == []: continue

            # Otherwise debump the residue

            self.write("Starting to debump %s...\n" % residue, 1)
            if self.debumpResidue(residue, conflictnames):
                 self.write("Debumping Successful!\n\n",1)
            else:
                 text = "WARNING: Unable to debump %s\n" % residue
                 self.write("********\n%s********\n\n" % text)
                 self.warnings.append(text)  

        self.write("Done.\n")

    def debumpResidue(self, residue, conflictnames):
        """
            Debump a specific residue.  Only should be called
            if the residue has been detected to have a conflict.
            If called, try to rotate about dihedral angles to
            resolve the conflict.

            Parameters
                residue:  The residue in question
                conflictnames:  A list of atomnames that were
                                rebuilt too close to other atoms
            Returns
                1 if successful, 0 otherwise
        """

        # Initialize some variables
        
        step = 0
        bestscore = 100
        anglenum = -1
        newcauses = []

        # Try (up to 10 times) to find a workable solution
        
        while step < 10:

            anglenum = self.pickDihedralAngle(residue, conflictnames, anglenum)
            
            if anglenum == -1: return 0
            
            self.write("Using dihedral angle number %i to debump the residue.\n" % anglenum, 1)

            for i in range(72):
                newangle = residue.dihedrals[anglenum] + 5.0
                self.setDihedralAngle(residue, anglenum, newangle)

                # Check for conflicts

                score = 0
                
                atomnames = residue.reference.dihedrals[anglenum].split()
                pivot = atomnames[2]
                moveablenames = self.getMoveableNames(residue, pivot)
                for name in moveablenames:
                    nearatoms = self.findNearbyAtoms(residue.getAtom(name))
                    for atom in nearatoms:
                        score += nearatoms[atom]
                
                if score == 0:
                    self.write("No conflicts found at angle %.2f.\n" % newangle, 1)
                    return 1

                # Set the best angle
                
                if score < bestscore:
                    bestangle = newangle

            self.setDihedralAngle(residue, anglenum, bestangle)
            step += 1
            

        # If we're here, debumping was unsuccessful

        return 0

    def calculateDihedralAngles(self):
        """
            Calculate the dihedral angle for every residue within the protein
        """
        for residue in self.protein.getResidues():
            if not isinstance(residue, Amino): continue
            residue.dihedrals = []

            refangles = residue.reference.dihedrals
            for di in refangles:
                coords = []
                atoms = di.split()
                for i in range(4):
                    atomname = atoms[i]
                    if residue.hasAtom(atomname):
                        coords.append(residue.getAtom(atomname).getCoords())
       
                if len(coords) == 4: angle = getDihedral(coords[0], coords[1], coords[2], coords[3])
                else: angle = None
                 
                residue.addDihedralAngle(angle)             

    def getClosestAtom(self, atom):
        """
            Get the closest atom that does not form a donor/acceptor pair.
            Used to detect potential conflicts.

            NOTE:  Cells must be set before using this function.

            Parameters
                atom:  The atom in question (Atom)
            Returns
                bestatom:  The closest atom to the input atom that does not
                           satisfy a donor/acceptor pair.
        """
        # Initialize some variables

        bestdist = 999.99
        bestatom = None
        residue = atom.residue

        # Get atoms from nearby cells
        
        closeatoms = self.cells.getNearCells(atom)

        # Loop through and see which is the closest
        
        for closeatom in closeatoms:
            closeresidue = closeatom.residue
            if closeresidue == residue: continue
            if not isinstance(closeresidue, Amino): continue
            if isinstance(residue, CYS):
                if residue.SSbondedpartner == closeatom: continue
                    
            # Also ignore if this is a donor/acceptor pair
                
            if atom.isHydrogen() and atom.bonds[0].hdonor \
               and closeatom.hacceptor: continue
            if closeatom.isHydrogen() and closeatom.bonds[0].hdonor \
                   and atom.hacceptor:
                continue

            dist = distance(atom.getCoords(), closeatom.getCoords())
            if dist < bestdist:
                bestdist = dist
                bestatom = closeatom
            
        return bestatom


    def findNearbyAtoms(self, atom):
        """
            Find nearby atoms for conflict-checking.  Uses
            neighboring cells to compare atoms rather than an all
            versus all O(n^2) algorithm, which saves a great deal
            of time.  There are several instances where we ignore
            potential conflicts; these include donor/acceptor pairs,
            atoms in the same residue, and bonded CYS bridges.

            Parameters
                atom:  Find nearby atoms to this atom (Atom)
            Returns
                nearatoms:  A list of atoms close to the atom.
        """
        # Initialize some variables

        nearatoms = {}
        residue = atom.residue
        cutoff = BUMP_DIST
        if atom.isHydrogen(): cutoff = BUMP_HDIST

        # Get atoms from nearby cells
        
        closeatoms = self.cells.getNearCells(atom)

        # Loop through and see if any are within the cutoff
        
        for closeatom in closeatoms:
            closeresidue = closeatom.residue
            if closeresidue == residue: continue
            if not isinstance(closeresidue, Amino): continue
            if isinstance(residue, CYS):
                if residue.SSbondedpartner == closeatom: continue
                    
            # Also ignore if this is a donor/acceptor pair
                
            if atom.isHydrogen() and atom.bonds[0].hdonor \
               and closeatom.hacceptor: continue
            if closeatom.isHydrogen() and closeatom.bonds[0].hdonor \
                   and atom.hacceptor:
                continue

            dist = distance(atom.getCoords(), closeatom.getCoords())
            if dist < cutoff:
                nearatoms[closeatom] = dist

        return nearatoms
        
    def pickDihedralAngle(self, residue, conflictnames, oldnum=None):
        """ 
            Choose an angle number to use in debumping

            Algorithm
                Instead of simply picking a random chiangle, this function
                uses a more intelligent method to improve efficiency.
                The algorithm uses the names of the conflicting atoms
                within the residue to determine which angle number
                has the best chance of fixing the problem(s). The method
                also insures that the same chiangle will not be run twice
                in a row.
            Parameters
                residue:    The residue that is being debumped (Residue)
                conflictnames: A list of atom names that are currently
                               conflicts (list)
                oldnum    : The old dihedral angle number (int)
            Returns
                bestnum    : The new dihedral angle number (int)
        """
        bestnum = -1
        best = 0
        for i in range(len(residue.dihedrals)):
            if i == oldnum: continue
            if residue.dihedrals[i] == None: continue

            score = 0
            atomnames = residue.reference.dihedrals[i].split()
            pivot = atomnames[2]
            
            moveablenames = self.getMoveableNames(residue, pivot)
         
            # If this pivot only moves the conflict atoms, pick it
            
            if conflictnames == moveablenames: return i 
            
            # Otherwise find the pivot with the most matches
            
            for name in conflictnames:
                if name in moveablenames:
                    score += 1
                    if score > best:
                        best = score
                        bestnum = i

        # Return the best angle.  If none were found, return -1.

        return bestnum

    def setDihedralAngle(self, residue, anglenum, angle):
        """
            Rotate a residue about a given angle. Uses the quatfit
            methods to perform the matrix mathematics.

            Parameters
                residue:   The residue to rotate
                anglenum:  The number of the angle to rotate as
                           listed in residue.dihedrals
                angle:     The desired angle.
                          
        """
        coordlist = []
        initcoords = []
        movecoords = []
        pivot = ""
        
        oldangle = residue.dihedrals[anglenum]
        diff = angle - oldangle

        atomnames = residue.reference.dihedrals[anglenum].split()

        pivot = atomnames[2]
        for atomname in atomnames:
            if residue.hasAtom(atomname):
                coordlist.append(residue.getAtom(atomname).getCoords())
            else:
                raise ValueError, "Error occurred while trying to debump!"

        initcoords = subtract(coordlist[2], coordlist[1])

        moveablenames = self.getMoveableNames(residue, pivot)

        for name in moveablenames:
            atom = residue.getAtom(name)
            movecoords.append(subtract(atom.getCoords(), coordlist[1]))
            
        newcoords = qchichange(initcoords, movecoords, diff)

        for i in range(len(moveablenames)):
            atom = residue.getAtom(moveablenames[i])
            self.cells.removeCell(atom)
            x = (newcoords[i][0] + coordlist[1][0])
            y = (newcoords[i][1] + coordlist[1][1])
            z = (newcoords[i][2] + coordlist[1][2])
            atom.set("x", x)
            atom.set("y", y)
            atom.set("z", z)
            self.cells.addCell(atom)
            

        # Set the new angle
        
        coordlist = []
        for atomname in atomnames:
            if residue.hasAtom(atomname):
                coordlist.append(residue.getAtom(atomname).getCoords())
            else:
                raise ValueError, "Error occurred while trying to debump!" 
        
        di = getDihedral(coordlist[0], coordlist[1], coordlist[2], coordlist[3])
        residue.dihedrals[anglenum] = di

    def getMoveableNames(self, residue, pivot):
        """
            Return all atomnames that are further away than the
            pivot atom.

            Parameters
                residue:  The residue to use
                pivot:    The pivot atomname
        """
        movenames = []
        refdist = residue.getAtom(pivot).refdistance
        for atom in residue.getAtoms():
            if atom.refdistance > refdist:
                movenames.append(atom.name)
        return movenames

    def setDonorsAndAcceptors(self):
        """
            Set the donors and acceptors within the protein
        """
        for residue in self.protein.getResidues():
            residue.setDonorsAndAcceptors()
        
    def runPROPKA(self, ph, ff, outname):
        """
            Run PROPKA on the current protein, setting protonation states to
            the correct values
            
            Parameters
               ph:  The desired pH of the system
               ff:  The forcefield name to be used
               outname: The name of the PQR outfile
        """
        self.write("Running propka and applying at pH %.2f... " % ph)

        # Initialize some variables

        linelen = 70
        txt = ""
        pkadic = {}
        warnings = []

        # Make sure PropKa has been installed.
        
        try:
            from propka.propkalib import runPKA
        except ImportError:
            text = "Couldn't find propka - make sure it has been installed!"
            raise ValueError, text

        # Reorder the atoms in each residue to start with N
   
        for residue in self.protein.getResidues():
            residue.reorder() 
        
        # Make a string with all non-hydrogen atoms

        for atom in self.protein.getAtoms():
            if not atom.isHydrogen():
                atomtxt = str(atom)
                if len(atomtxt) + 1 != linelen:
                    print "Atom line length (%i) does not match constant (%i)!" % \
                          ((len(atomtxt) +1), linelen)
                    sys.exit()
                txt += "%s\n" % atomtxt

        # The length of the overall text/line length ratio should be
        # the number of atoms without remainder

        txtlen = len(txt)
        if txtlen % linelen != 0:
            raise ValueError, "Extra characters in pka string!"


        # Run PropKa

        numatoms = int(txtlen) / linelen

        runPKA(numatoms, txt, outname)
        
        # Parse the results
        
        pkafile = open(outname)
        summary = 0
        while 1:
            line = pkafile.readline()
            if line == "": break
            if line.startswith("SUMMARY"): summary = 1
            elif line.startswith("-"): summary = 0
            elif summary:
                words = string.split(string.strip(line))
                key = ""
                for i in range(len(words) - 1):
                    key = "%s %s" % (key,words[i])
                key = string.strip(key)
                pkadic[key] = float(words[-1])
            
        if len(pkadic) == 0: return

        # Now apply each pka to the appropriate residue

        for residue in self.protein.getResidues():
            if not isinstance(residue, Amino): continue
            resname = residue.name
            resnum = residue.resSeq
            chainID = residue.chainID
    
            if residue.isNterm:
                key = "N+ %i %s" % (resnum, chainID)
                key = string.strip(key)
                if key in pkadic: 
                    value = pkadic[key]
                    del pkadic[key]
                    if ph >= value:
                        if ff in ["amber","charmm","tyl06"]:
                            warn = ("N-terminal %s" % key, "neutral")
                            warnings.append(warn)
                        else:
                            self.applyPatch("NEUTRAL-NTERM", residue)
                                
            if residue.isCterm:
                key = "C- %i %s" % (resnum, chainID)
                key = string.strip(key)
                if key in pkadic:
                    value = pkadic[key]
                    del pkadic[key]
                    if ph < value:
                        if ff in ["amber","charmm", "tyl06"]:
                            warn = ("C-terminal %s" % key, "neutral")
                            warnings.append(warn)
                        else:
                            self.applyPatch("NEUTRAL-CTERM", residue)
                                
            key = "%s %i %s" % (resname, resnum, chainID)
            key = string.strip(key)
            if key in pkadic:
                value = pkadic[key]
                del pkadic[key]
                if resname == "ARG" and ph >= value:
                    warn = (key, "neutral")
                    warnings.append(warn)
                elif resname == "ASP" and ph < value:
                    if residue.isCterm and ff in ["amber","tyl06"]:
                        warn = (key, "Protonated at C-Terminal")
                        warnings.append(warn)
                    elif residue.isNterm and ff in ["amber","tyl06"]:
                        warn = (key, "Protonated at N-Terminal")
                        warnings.append(warn)
                    else:
                        self.applyPatch("ASH", residue)
                elif resname == "CYS" and ph >= value:
                    if ff == "charmm":
                        warn = (key, "negative")
                        warnings.append(warn)
                    else:
                        self.applyPatch("CYM", residue)
                elif resname == "GLU" and ph < value:
                    if residue.isCterm and ff in ["amber","tyl06"]:
                        warn = (key, "Protonated at C-Terminal")
                        warnings.append(warn)
                    elif residue.isNterm and ff in ["amber","tyl06"]:
                        warn = (key, "Protonated at N-Terminal")
                        warnings.append(warn)
                    else:
                        self.applyPatch("GLH", residue)
                elif resname == "HIS" and ph < value:
                    self.applyPatch("HIP", residue)              
                elif resname == "LYS" and ph >= value:
                    if ff == "charmm":
                        warn = (key, "neutral")
                        warnings.append(warn)
                    elif ff in ["amber","tyl06"] and residue.get("isCterm"):
                        warn = (key, "neutral at C-Terminal")
                        warnings.append(warn)
                    elif ff == "tyl06" and residue.get("isNterm"):
                        warn = (key, "neutral at N-Terminal")
                        warnings.append(warn)  
                    else:
                        self.applyPatch("LYN", residue)
                elif resname == "TYR" and ph >= value:
                    if ff in ["charmm", "amber", "tyl06"]:
                        warn = (key, "negative")
                        warnings.append(warn)
                    else:
                        self.applyPatch("TYM", residue)
                        
        if len(warnings) > 0:
            init = "WARNING: Propka determined the following residues to be\n"
            self.warnings.append(init)
            init = "         in a protonation state not supported by the\n"
            self.warnings.append(init)
            init = "         %s forcefield!\n" % ff
            self.warnings.append(init)
            init = "         All were reset to their standard pH 7.0 state.\n"
            self.warnings.append(init)
            self.warnings.append("\n")
            for warn in warnings:
                text = "             %s (%s)\n" % (warn[0], warn[1])
                self.warnings.append(text)
            self.warnings.append("\n")
   
        if len(pkadic) > 0:
            warn = "         PDB2PQR could not identify the following residues\n"
            self.warnings.append(warn)
            warn = "         and residue numbers as returned by propka:\n"
            self.warnings.append(warn)
            self.warnings.append("\n")
            for item in pkadic:
                text = "             %s\n" % item
                self.warnings.append(text)
            self.warnings.append("\n")
        self.write("Done.\n")

class Cells:
    """
        The cells object provides a better way to search for nearby atoms. A
        pure all versus all search is O(n^2) - for every atom, every other atom
        must be searched.  This is rather inefficient, especially for large
        proteins where cells may be tens of angstroms apart.  The cell class
        breaks down the xyz protein space into several 3-D cells of desired
        size - then by simply examining atoms that fall into the adjacent
        cells one can quickly find nearby cells.

        NOTE:  Ideally this should be somehow separated from the routines
               object...
        """
    def __init__(self, cellsize):
        """
            Initialize the cells.

            Parameters
                cellsize:  The size of each cell (int)
        """
        self.cellmap = {}
        self.cellsize = cellsize

    def assignCells(self, protein):
        """
            Place each atom in a virtual cell for easy neighbor comparison
        """  
        for atom in protein.getAtoms():
            atom.cell = None
            self.addCell(atom)

    def addCell(self, atom):
        """
            Add an atom to the cell

            Parameters
                atom:  The atom to add (atom)
        """
        size = self.cellsize
        x = atom.get("x")
        if x < 0: x = (int(x)-1)/size*size
        else: x = int(x)/size*size
        y = atom.get("y")
        if y < 0: y = (int(y)-1)/size*size
        else: y = int(y)/size*size
        z = atom.get("z")
        if z < 0: z = (int(z)-1)/size*size
        else: z = int(z)/size*size
        key = (x,y,z)
        try:
            self.cellmap[key].append(atom)
        except KeyError:
            self.cellmap[key] = [atom]
        atom.set("cell", key)

    def removeCell(self, atom):
        """
             Remove the atom from a cell

             Parameters
                 atom:   The atom to add (atom)
        """
        oldcell = atom.get("cell")
        if oldcell == None: return
        atom.set("cell", None)
        self.cellmap[oldcell].remove(atom)

    def getNearCells(self, atom):
        """
            Find all atoms in bordering cells to an atom

            Parameters
                atom:  The atom to use (atom)
            Returns
                closeatoms:  A list of nearby atoms (list)
        """
        size = self.cellsize
        closeatoms = []
        cell = atom.get("cell")
        x = cell[0]
        y = cell[1]
        z = cell[2]
        for i in range(-1*size,2*size,size):
            for j in range(-1*size,2*size,size):
                for k in range(-1*size,2*size,size):
                    newkey = (x+i, y+j, z+k)
                    try:
                        newatoms = self.cellmap[newkey]
                        for atom2 in newatoms:
                            if atom == atom2: continue
                            closeatoms.append(atom2)
                    except KeyError: pass
                        
        return closeatoms
