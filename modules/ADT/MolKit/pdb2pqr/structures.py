"""
    Structures for PDB2PQR

    This module contains the structure objects used in PDB2PQR and their
    associated methods.

    Todd Dolinsky (todd@ccb.wustl.edu)
    Washington University in St. Louis

"""

__date__ = "22 October 2003"
__author__ = "Todd Dolinsky"

BACKBONE = ["N","CA","C","O","O2","HA","HN","H","tN"]

import string
from pdb import *

class Chain:
    """
        Chain class

        The chain class contains information about each chain within a given
        Protein object.  
    """

    def __init__(self, chainID):
        """
            Initialize the class

            Parameters
                chainID: The chainID for this chain as denoted in
                         the PDB file (string)
        """

        self.chainID = chainID
        self.residues = []

    def get(self, name):
        """
            Get a member of the Chain class

            Parameters
                name:     The name of the member
            Possible Values
                ID:       The ID of the chain
                Residues: The list of residues within the Chain
            Returns
                item:     The value of the member
        """
        if name == "atoms": self.getAtoms()
        else:
            try:
                item = getattr(self, name)
                return item
            except AttributeError:
                message = "Unable to get object \"%s\" in class Chain" % name
                raise ValueError, message

    def addResidue(self, residue):
        """
            Add a residue to the chain

            Parameters
                residue: The residue to be added (Residue)
        """
        self.residues.append(residue)

    def numResidues(self):
        """
            Get the number of residues for the chain

            Returns
                count:  Number of residues in the chain (int)
        """
        count = 0
        for residue in self.residues:
            count += 1
        return count

    def renumberResidues(self):
        """
            Renumber Atoms based on actual Residue number and not PDB resSeq
        """
        count = 1
        for residue in self.residues:
            residue.setResSeq(count)
            count += 1

    def numAtoms(self):
        """
            Get the number of atoms for the chain

            Returns
                count:  Number of atoms in the chain (int)
        """
        count = len(self.getAtoms())       
        return count
    
    def getAtoms(self):
        """
            Return a list of Atom objects contained in this chain

            Returns
                atomlist: List of Atom objects (list)
        """
        atomlist = []
        for residue in self.residues:
            myList = residue.get("atoms")
            for atom in myList:
                atomlist.append(atom)
        return atomlist


class Residue:
    """
        Residue class

        The residue class contains a list of Atom objects associated with that
        residue and other helper functions.
    """

    def __init__(self, atoms, sampleAtom):
        """
            Initialize the class

            Parameters
                atoms:      A list of Atom objects to be stored in this class
                            (list)
                sampleAtom: The final listed atom of the residue (Atom)
        """
        self.atoms = []
        self.name = sampleAtom.resName
        self.chainID = sampleAtom.chainID
        self.resSeq = sampleAtom.resSeq
        self.iCode = sampleAtom.iCode
        self.SSbonded = 0
        self.SSbondpartner = None
        self.type = 0
        self.map = {}
        self.chiangles = []
        self.isNterm = 0
        self.isCterm = 0
        self.is3term = 0
        self.is5term = 0
        self.missing = []
        self.debumpAtoms = []
        self.naname = None
        
        atomclass = ""
        for a in atoms:
            if isinstance(a,ATOM):
                atomclass = "ATOM"
            elif isinstance(a, HETATM):
                atomclass = "HETATM"
            atom = Atom(a, atomclass, self)
            atomname = atom.get("name")
            if atomname not in self.map:
                self.addAtom(atom)
            else: # Don't add duplicate atom              
                oldatom = self.getAtom(atomname)
                oldatom.set("altLoc","")

        if self.name == "HOH":
            self.name = "WAT"
            for atom in self.atoms:
                atom.set("resName","WAT")

    def get(self, name):
        """
            Get a member of the Residue class

            Parameters
                name:          The name of the member (string)
            Possible Values
                atoms:         The atoms in the residue
                name:          The name of the residue
                chainID:       The chainID associated with the residue
                resSeq:        The sequence number of the residue
                icode:         The iCode of the residue
                SSbonded:      1 if the residue has a SS bond, 0 otherwise
                SSbondpartner: The residue of the bond partner
                type:          The type associated with this residue
                isNterm:       1 if the residue is the N-Terminus, 0 otherwise
                isCterm:       1 if the residue is the C-Terminus, 0 otherwise
                missing:     List of missing atoms of the residue
            Returns
                item:          The value of the member
        """
        try:
            item = getattr(self, name)
            return item
        except AttributeError:
            message = "Unable to access object \"%s\" in class Residue" % name
            raise ValueError, message

    def set(self, name, value):
        """
            Set a member of the Residue class to a specific value 

            Parameters
                name:          The name of the object to set (string)
                value:         The object to append
            Possible Values
                atoms:         The atoms in the residue
                name:          The name of the residue
                chain:         The chainID associated with the residue
                resSeq:        The sequence number of the residue
                icode:         The iCode of the residue
                SSbonded:      1 if the residue has a SS bond, 0 otherwise
                SSbondpartner: The residue of the bond partner
                type:          The type associated with this residue
                isNterm:       1 if the residue is the N-Terminus, 0 otherwise
                isCterm:       1 if the residue is the C-Terminus, 0 otherwise
                isDirty:       1 if the residue is not missing atoms,
                               0 otherwise
            Notes
                resSeq points to the residue.setResSeq function
            Returns
                item:          The value of the member   
        """
        if name == "resSeq": self.setResSeq(value)
        else:
            try:
                setattr(self, name, value)
            except AttributeError:
                message = "Unable to set object \"%s\" in class Residue" % name
                raise ValueError, message

    def checkAtomNames(self):
        """
            Check to see if there are any misnamed hydrogens within the
            residue.  Converts hydrogens of type 1HH1 and 2H* to HH11 and
            H*2 to easily compare with the Amino Acid definition file. Also
            converts H*1 and H*2 to H*2 and H*3 when necessary. Rename the
            atom and update the residue.
        """
        resname = self.name
        names = {}
        newnames = {}
        atomlist = []

        # Make a copy of the atomlist

        for atom in self.atoms:
            atomlist.append(atom)

        # First make a list of names with int as last digit
        
        for atom in atomlist:
            atomname = atom.get("name")
            newname = atomname
            if atomname[0].isdigit() and len(atomname) > 1:
                if atomname[1] == "H":
                    newname = atomname[1:] + atomname[0]
            names[newname] = atomname

        # Now renumber the hydrogens as needed
            
        for n in names:
            newname = n
            if n == "HG11" and "HG13" not in names: newname = "HG12"
            elif n == "HG12" and "HG13" not in names: newname = "HG13"
            elif n == "HA1" and "HA2" in names: newname = "HA2"
            elif n == "HA2" and "HA1" in names: newname = "HA3"
            elif n == "HB1" and "HB3" not in names: newname = "HB2"
            elif n == "HB2" and "HB3" not in names: newname = "HB3"
            elif n == "HG1" and resname not in ["THR","SER","CYS"] and "HG3" not in names: newname = "HG2"
            elif n == "HG2" and "HG3" not in names: newname = "HG3"
            elif n == "HD1" and resname in ["ARG","LYS","PRO"]: newname = "HD2"
            elif n == "HD2" and resname in ["ARG","LYS","PRO"] and "HD1" in names: newname = "HD3"
            elif n == "HE1" and resname == "LYS": newname = "HE2"
            elif n == "HE2" and resname == "LYS" and "HE1" in names: newname = "HE3"
            elif n == "H1" and self.get("isNterm") and "H" not in names: newname = "H"
            elif n == "H1" and resname == "ACE": newname = "HH31"
            elif n == "H2" and resname == "ACE": newname = "HH32"
            elif n == "H3" and resname == "ACE": newname = "HH33"
            elif n == "HN": newname = "H"
            elif n == "HD1" and resname == "ILE": newname = "HD11"
            elif n == "HD2" and resname == "ILE": newname = "HD12"
            elif n == "HD3" and resname == "ILE": newname = "HD13"
            elif n == "HG1" and resname in ["SER","CYS"]: newname = "HG"
            
            newnames[names[n]] = newname

        # Now update the residue, being careful not to overwrite existing atoms

        old = {}
        for atom in atomlist:
            atomname = atom.get("name")
            newname = newnames[atomname]
            if atomname != newname:
                if atomname in old:
                    atom = old[atomname]
                    atom.set("name",newname)
                    self.map[newname] = atom
                    del old[atomname]
                elif newname in self.map:
                    oldatom = self.map[newname]
                    old[newname] = oldatom
                    self.renameAtom(atomname, newname)
                else:
                    self.renameAtom(atomname, newname)

        # There should be nothing left in old

        if len(old) != 0:
            raise ValueError, "Error Occurred when renaming hydrogens: %s" % old

    def updateIntraBonds(self, defresidue):
        """
            Update the IntraBonds for each atom in the residue

            Parameters
                defresidue:  The definition residue (DefinitionResidue)
        """
        for atom in self.atoms:
            atomname = atom.get("name")
            defatom = defresidue.getAtom(atomname)
            atom.set("intrabonds",[])
            if defatom == None:
                if self.isCterm and atomname == "OXT":
                    continue
                elif self.isNterm and atomname in ["H2","H3"]:
                    continue
                else:
                    raise ValueError, "Atom %s not found in updateIntraBonds!" % atomname
            for bondatomname in defatom.get("intrabonds"):
                if self.getAtom(bondatomname):
                    atom.addIntraBond(bondatomname)
    
    def numAtoms(self):
        """
            Get the number of atoms for the residue

            Returns
                count:  Number of atoms in the residue (int)
        """
        count = len(self.atoms)
        return count
                    
    def setResSeq(self, value):
        """
            Set the atom field resSeq to a certain value and
            change the residue's information.  The icode field is no longer
            useful.

            Parameters
                value:  The new value of resSeq (int)
        """
        self.iCode = ""
        self.resSeq = value
        for atom in self.atoms:
            atom.set("resSeq",value)
            atom.set("iCode","")
            
    def addAtom(self, atom):
        """
            Add the atom object to the residue.

            Parameters
                atom: The object to be added (ATOM)
        """
        self.atoms.append(atom)
        self.map[atom.get("name")] = atom

    def removeAtom(self, atomname):
        """
            Remove an atom from the residue object.

            Parameters
                atomname: The name of the atom to be removed (string)
        """
        atom = self.map[atomname]
        del self.map[atomname]
        for i in range(self.numAtoms()):
            a = self.atoms[i]
            if atom == a:
                del self.atoms[i]
                break
        for i in range(len(self.debumpAtoms)):
            a = self.debumpAtoms[i]
            if atom == a:
                del self.debumpAtoms[i]
                break
                

    def renameAtom(self, oldname, newname):
        """
            Rename an atom to a new name

            Parameters
                oldname: The old atom name (string)
                newname: The new atom name (string)
        """
        atom = self.map[oldname]
        atom.set("name",newname)
        self.map[newname] = atom
        del self.map[oldname]
        
    def createAtom(self, name, newcoords, type):
        """
            Add a new atom object to the residue. Uses an atom
            currently in the residue to seed the new atom
            object, then replaces the coordinates and name accordingly.

            Parameters
                name:      The name of the new atom (string)
                newcoords: The x,y,z coordinates of the new atom (list)
                type:      The type of atom, ATOM or HETATM
        """
        oldatom = self.atoms[0]
        newatom = Atom(oldatom, type, self)
        newatom.set("x",newcoords[0])
        newatom.set("y",newcoords[1])
        newatom.set("z",newcoords[2])
        newatom.set("name", name)
        newatom.set("occupancy",1.00)
        newatom.set("tempFactor",0.00)
        self.addAtom(newatom) 

    def addChiangle(self, value):
        """
            Add the value to the list of chiangles

            Parameters
                value: The value to be added (float)
        """
        self.chiangles.append(value)

    def addMissing(self, value):
        """
            Add the value to the list of missing atoms

            Parameters
                value: The name of the missing atom (string)
        """
        self.missing.append(value)

    def getAtom(self, name):
        """
            Retrieve an atom from the mapping

            Parameters
                resname: The name of the residue to retrieve (string)
        """
        try:
            return self.map[name]
        except KeyError:
            return None

    def getCharge(self):
        """
            Get the total charge of the residue.  In order to get rid
            of floating point rounding error, do the string
            transformation.

            Returns:
                charge: The charge of the residue (float)
        """
        charge = 0.0
        for atom in self.atoms:
            charge = charge + atom.get("ffcharge")

        charge = float("%.4f" % charge)
        return charge
    
    def addDebumpAtom(self, atom):
        """
            Add an atom to the check for debumping

            Parameters
                atom:  The atom to add to the list
        """
        self.debumpAtoms.append(atom)
        
class Atom(ATOM):
    """
        Class Atom

        The Atom class inherits off the ATOM object in pdb.py.  It is used
        for adding fields not found in the pdb that may be useful for analysis.
        Also simplifies code by combining ATOM and HETATM objects into a
        single class.
    """
    
    def __init__(self, atom, type, residue):
        """
            Initialize the new Atom object by using the old object.

            Parameters
                atom:    The original ATOM object (ATOM)
                type:    Either ATOM or HETATM (string)
                residue: A pointer back to the parent residue object (Residue)
        """
        if type == "ATOM" or type == "HETATM":
            self.type = type
        else:
            raise ValueError, "Invalid atom type %s!"
        self.serial = atom.serial
        self.name = atom.name
        self.altLoc = atom.altLoc
        self.resName = atom.resName
        self.chainID = atom.chainID
        self.resSeq = atom.resSeq
        self.iCode = atom.iCode
        self.x = atom.x
        self.y = atom.y
        self.z = atom.z
        self.occupancy = atom.occupancy
        self.tempFactor = atom.tempFactor
        self.segID = atom.segID
        self.element = atom.element
        self.charge = atom.charge
        
        self.intrabonds = []
        self.extrabonds = []
        self.residue = residue
        self.radius = 0.0
        self.ffcharge = 0.0
        self.hdonor = 0
        self.hacceptor = 0
        self.cell = None
        
    def __str__(self):
        """
            Returns a string of the new atom type.  Uses the ATOM string
            output but changes the first field to either by ATOM or
            HETATM as necessary.

            Returns
                out: String with ATOM/HETATM field set appropriately
        """
        orig = ATOM.__str__(self)
        type = string.ljust(self.type, 6)[:6]
        ffcharge = "%.4f" % self.ffcharge
        ffradius = "%.4f" % self.radius
        charge = string.rjust(ffcharge, 7)[:7]
        radius = string.ljust(ffradius, 6)[:6]
        out = "%s%s %s %s" % (type, orig[6:-20], charge, radius)
        out = "%s %s" % (out[:21], out[22:]) # Eliminate the chain ID
        return out

    def get(self, name):
        """
            Get a member of the Atom class

            Parameters
                name:       The name of the member (string)
            Possible Values
                type:       The type of Atom (either ATOM or HETATM)
                serial:     Atom serial number
                name:       Atom name
                altLoc:     Alternate location
                resName:    Residue name
                chainID:    Chain identifier
                resSeq:     Residue sequence number
                iCode:      Code for insertion of residues
                x:          Orthogonal coordinates for X in Angstroms.
                y:          Orthogonal coordinates for Y in Angstroms.
                z:          Orthogonal coordinates for Z in Angstroms.
                occupancy:  Occupancy
                tempFactor: Temperature Factor
                segID:      Segment identifier
                element:    Element symbol
                charge:     Charge on the atom
                bonds:      The bonds associated with the atom
                interbonds: The intrabonds associated with the atom
                extrabonds: The extrabonds assocaited with the atom
                residue:    The parent residue of the atom
                radius:     The radius of the atom
                ffcharge:   The forcefield charge on the atom
                hdonor:     Whether the atom is a hydrogen donor
                hacceptor:  Whether the atom is a hydrogen acceptor
            Returns
                item:       The value of the member
        """
        try:
            item = getattr(self, name)
            return item
        except AttributeError:
            message = "Unable to access object \"%s\" in class Atom" % name
            raise ValueError, message

    def set(self, name, value):
        """
            Set a member of the Atom class

            Parameters
                name:       The name of the member (string)
                value:      The value to set the member to
            Possible Values
                type:       The type of Atom (either ATOM or HETATM)
                serial:     Atom serial number
                name:       Atom name
                altLoc:     Alternate location
                resName:    Residue name
                chainID:    Chain identifier
                resSeq:     Residue sequence number
                iCode:      Code for insertion of residues
                x:          Orthogonal coordinates for X in Angstroms.
                y:          Orthogonal coordinates for Y in Angstroms.
                z:          Orthogonal coordinates for Z in Angstroms.
                occupancy:  Occupancy
                tempFactor: Temperature Factor
                segID:      Segment identifier
                element:    Element symbol
                charge:     Charge on the atom
                residue:    The parent residue of the atom
                radius:     The radius of the atom
                ffcharge:   The forcefield charge on the atom
                hdonor:     Whether the atom is a hydrogen donor
                hacceptor:  Whether the atom is a hydrogen acceptor
            Returns
                item:       The value of the member
        """
        try:
            setattr(self, name, value)
        except AttributeError:
            message = "Unable to set object \"%s\" in class Atom" % name
            raise ValueError, message   

    def getCoords(self):
        """
            Return the x,y,z coordinates of the atom in list form

            Returns
                List of the coordinates (list)
        """
        return [self.x, self.y, self.z]

    def addIntraBond(self, bondedatom):
        """
            Add a bond to the list of intrabonds

            Parameters:
                bondedatom: The atom to bond to (Atom)
        """
        self.intrabonds.append(bondedatom)

    def addExtraBond(self, bondedatom):
        """
            Add a bond to the list of extrabonds

            Parameters:
                bondedatom: The atom to bond to (Atom)
        """
        self.extrabonds.append(bondedatom)

    def isHydrogen(self):
        """
            Is this atom a Hydrogen atom?

            returns
                value: 1 if Atom is a Hydrogen, 0 otherwise
        """
        value = 0
        if self.name[0] == "H":
            value = 1
        return value

    def isBackbone(self):
        """
            Return true if atom name is in backbone, otherwise false

            Returns
                state: 1 if true, 0 if false
        """
        state = 0
        if self.name in BACKBONE:
            state = 1
        return state
