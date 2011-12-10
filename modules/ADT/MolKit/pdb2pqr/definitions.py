"""
    Definitions for PDB2PQR

    This file contains classes associated with Amino Acid and Rotamer
    definitions as used by PDB2PQR.

    Definition Files Created by
    Jens Erik Nielsen

    Ported to Python by Todd Dolinsky (todd@ccb.wustl.edu)
    Washington University in St. Louis
"""
    
__date__ = "8 September 2004"
__author__ = "Jens Erik Nielsen, Todd Dolinsky"

AAFILE = "AA.DAT"
NAFILE = "NA.DAT"
ROTAMERFILE = "ROTAMER.DAT"

import os
from pdb import *
from utilities import *
from structures import *
from routines import *

class Definition:
    """
        Definition class

        The Definition class contains the structured definitions found
        in the files and several mappings for easy access to the information.
    """
    def __init__(self):
        """
            Create a new Definition Object
        """
        self.AAdef = self.readDefinition(AAFILE)
        self.NAdef = self.readDefinition(NAFILE)

    def getAA(self):
        """
            Get the Amino Acid definition

            Returns
                The Amino Acid definition in self.chains[0]
        """
        return self.AAdef

    def getNA(self):
        """
            Get the Nucleic Acid definition

            Returns
                 The Nucleic Acid definition - a list of Nucleic Acid
                 residues
        """
        return self.NAdef
    
    def parseRotamer(self, reslines):
        """
            Parse ROTAMER.DAT, obtaining information about each atom and its
            position.

            Parameters
                reslines:  A list of lines containing each of the atoms of the
                           residue (list)
            Returns
                myResidue: The parsed residue (DefinitionResidue)
        """
        name = reslines[0][17:20]
        restype = 1
        myResidue = DefinitionResidue(name, restype)
        for i in range(len(reslines)):
            entries = string.split(reslines[i])
            atomname = entries[2]
            x = float(entries[5])
            y = float(entries[6])
            z = float(entries[7])

            atom = DefinitionAtom(i, atomname, name, x, y, z)
            myResidue.addAtom(atom)

        return myResidue

    def parseDefinition(self, reslines, name, restype):
        """
            Parse the definition file, obtaining information about each atom,
            its position, and its bonding information.

            Parameters
                reslines:  A list of lines containing each of the atoms of the
                           residue (list)
                name:      The name of the Residue (string)
                restype:   The type of the residue (int)
            Returns
                myResidue: The parsed residue (DefinitionResidue)
        """
        myResidue = DefinitionResidue(name, restype)
        refatom = -1 
    
        for i in range(0,len(reslines)-2):
            entries = string.split(reslines[i])
            atomname = entries[0]
            x = float(entries[1])
            y = float(entries[2])
            z = float(entries[3])
            
            atom = DefinitionAtom(i, atomname, name, x, y, z)
            myResidue.addAtom(atom)
        
            if atomname == "CA" and restype == 1: refatom = i
          
        line = reslines[-2]
        bonds = string.split(line)
        bondmap = {}
        for i in range(0,len(bonds),2):
            bondA = int(bonds[i])-1
            bondB = int(bonds[i+1])-1
            atomA = myResidue.get("atoms")[bondA]
            atomB = myResidue.get("atoms")[bondB]
            
            atomA.addIntraBond(atomB.get("name"))
            atomB.addIntraBond(atomA.get("name"))
            try:
                bondmap[bondA].append(bondB)
            except KeyError:
                bondmap[bondA] = [bondB]
            try:
                bondmap[bondB].append(bondA)
            except KeyError:
                bondmap[bondB] = [bondA]
                    
        line = reslines[-1]
        dihedrals = string.split(line)
        if int(dihedrals[0]) > 0:
            for i in range(1,len(dihedrals)):
                dihedralA = int(dihedrals[i]) - 1
                myResidue.addDihedral(myResidue.get("atoms")[dihedralA].get("name"))
  
        if len(myResidue.get("dihedralatoms")) != int(dihedrals[0]) * 4:
            raise ValueError, "Corrupt entry for torsion angles when parsing %s" % name

        if restype == 1:
            for i in range(myResidue.numAtoms()):
                atom = myResidue.get("atoms")[i]
                if atom.isBackbone():
                    atom.set("refdistance",-1)
                else:
                    atom.set("refdistance", len(shortestPath(bondmap, i, refatom)) - 1)
        return myResidue
    
    def readDefinition(self, defpath):
        """
            Read a definition file

            Parameters
                deffile: The path to the definition file (string)

            Returns
                def:  The definition chain (AADefinition)
        """
        if not os.path.isfile(defpath):
            for path in sys.path:
                testpath = "%s/%s" % (path, defpath)
                if os.path.isfile(testpath):
                    defpath = testpath
                    break
                    
        if not os.path.isfile(defpath):
            raise ValueError, "%s not found!" % defpath
        
        file = open(defpath)
        lines = file.readlines()
        file.close()
        reslines = []
        thisdef = DefinitionChain(defpath)

        for line in lines:
            if line.startswith("//"): pass
            elif line.startswith("*"):
                if len(reslines) > 0:
                    ids = string.split(reslines[0])
                    name = ids[0]
                    restype = int(ids[1])
                    reslines = reslines[1:]  
                    residue = self.parseDefinition(reslines, name, restype)
                    thisdef.addResidue(residue)
                    reslines = []
            else:
                reslines.append(string.strip(line))

        thisdef.renumberResidues()
        return thisdef        

    def readRotamerDefinition(self):
        """
            Read the Rotamer definitions
        """
        if os.path.isfile(ROTAMERFILE):
            file = open(ROTAMERFILE)
            lines =  file.readlines()
            file.close()
            reslines = []
            rotamerdef = Chain("ROTAMER")
            
            for line in lines:
                if line.startswith("*"): pass
                elif line.startswith("TER"):
                    if len(reslines) > 0:
                        residue = self.parseRotamer(reslines)
                        rotamerdef.addResidue(residue)
                        reslines = []
                else:
                    reslines.append(string.strip(line))

            rotamerdef.renumberResidues()
            self.chains.append(rotamerdef)
        else:
            raise ValueError, "%s not found!" % ROTAMERFILE

class DefinitionChain(Chain):
    """
        DefinitionChain class

        The DefinitionChain class extends the chain class to provide
        lookups for atom information.
    """
    def __init__(self, ID):
        """
           Initialize like the Chain constructor, but add necessary
           features

           Parameters
               ID: The ID of the chain
        """
        Chain.__init__(self, ID)
        self.residuemap = {}
        
    def addResidue(self, residue):
        """
            Add a residue to the chain

            Parameters
                residue: The residue to be added (Residue)
        """
        self.residues.append(residue)
        self.residuemap[residue.get("name")] = residue

    def getResidue(self, name):
        """
            Retrieve a residue from the mapping

            Parameters
                name: The name of the residue to retrieve (string)
        """
        try:
            return self.residuemap[name]
        except KeyError:
            return None

class DefinitionResidue(Residue):
    """
        DefinitionResidue class

        The DefinitionResidue class extends the Residue class to allow for a
        trimmed down initializing function.
    """
    def __init__(self, name, type):
        """
            Initialize the class using a few parameters

            Parameters:
                name: The abbreviated amino acid name of the DefinitionResidue
                type: The typecode associated with the residue
                      Available typecodes are:
                          1: Protein residue
                          2: Drug/small-molecule
                          3: Water
                number: ID number for residue
        """
        self.atoms = []
        self.dihedralatoms = []
        self.resSeq = 0
        self.type = type
        self.name = name
        self.map = {}

    def addDihedral(self, atom):
        """
            Add the atom to the list of dihedral bonds

            Parameters:
                atom: The atom to be added
        """
        self.dihedralatoms.append(atom)
        
    def makeBondList(self, residue, atomname):
        """
            For the given atomname, make a list of bonded atoms.
            First get all atoms present in the residue that are
            directly bonded to the atom - if this number is
            less than REFATOM_SIZE, take those atoms that are present and
            bonded to initial bond list and use them.

            Parameters
                residue:  The residue to check for present atoms (Residue)
                atomname: The atom name to sedd the list of bonds (string)
            Returns
                bonds:    A list of atomnames that are within two bonds of
                          the atom and present in residue (list)
        """
        bonds = []
        if atomname == "OXT":
            if "C" in residue.get("map"):
                bonds.append("C")
            else:
                return bonds
        else:
            defatom = self.getAtom(atomname)
            defbonds = defatom.get("intrabonds")
            for bondname in defbonds:
                if bondname in residue.get("map"):
                    bonds.append(bondname)
        if len(bonds) < REFATOM_SIZE:
            for bond in bonds:
                newatom = self.getAtom(bond)
                newbonds = newatom.get("intrabonds")
                for newname in newbonds:
                    if newname in residue.get("map") and newname not in bonds:
                        bonds.append(newname)
                        if len(bonds) == REFATOM_SIZE:
                            return bonds
            return bonds
        else:
            return bonds
        
class DefinitionAtom(Atom):

    """
        Class DefinitionAtom

        The DefinitionAtom class inherits off the Atom class.  It provides
        a trimmed down version of the initializating function from the Atom
        class for the definition files.
    """
    
    def __init__(self, serial, name, resName, x, y, z):
        """
            Initialize using a few basic parameters - set all other fields
            to null, which is necessary for debugging output by using the
            string function in the parent class.

            Parameters
                serial:  Atom serial number (int)
                name:    Atom name. (string)
                resName: Residue name. (string)
                resSeq:  Residue sequence number. (int)
                x:       Orthogonal coordinates for X in Angstroms. (float)
                y:       Orthogonal coordinates for Y in Angstroms. (float)
                z:       Orthogonal coordinates for Z in Angstroms. (float)
        """
        self.type = "ATOM"
        self.serial = serial
        self.name = name
        self.altLoc = ""
        self.resName = resName
        self.chainID = ""
        self.resSeq = 1
        self.iCode = ""
        self.x = x
        self.y = y
        self.z = z
        self.occupancy = 0.0
        self.tempFactor = 0.0
        self.segID = ""
        self.element = ""
        self.charge = ""
        self.ffcharge = 0.0
        self.radius = 0.0

        self.intrabonds = []
        self.residue = None
        self.refdistance = None

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




