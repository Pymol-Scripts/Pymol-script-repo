"""
    Structures for PDB2PQR

    This module contains the structure objects used in PDB2PQR and their
    associated methods.

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
__author__ = "Todd Dolinsky"

BACKBONE = ["N","CA","C","O","O2","HA","HN","H","tN"]

import string
from pdb import *
from utilities import *
from quatfit import *

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

    def getResidues(self):
        """
            Return a list of Residue objects in this chain
        """
        return self.residues
    
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

    def __init__(self, atoms):
        """
            Initialize the class

            Parameters
                atoms:      A list of Atom objects to be stored in this class
                            (list)
        """

        sampleAtom = atoms[-1]
        
        self.atoms = []
        self.name = sampleAtom.resName
        self.chainID = sampleAtom.chainID
        self.resSeq = sampleAtom.resSeq
        self.iCode = sampleAtom.iCode
       
        self.map = {}

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

    def __str__(self):
        """
            Basic string representation for debugging
        """
        text = "%s %s %i" % (self.name, self.chainID, self.resSeq)
        return text

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
                isNterm:       # of hydrogens if the residue is the N-Terminus, 0 otherwise
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
                isNterm:       # of hydrogens if the residue is the N-Terminus, 0 otherwise
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

    def setChainID(self, value):
        """
           Set the chainID field to a certain value
        """
        self.chainID = value
        for atom in self.atoms:
            atom.set("chainID", value)
        
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

        # Delete the atom from the map
        
        atom = self.map[atomname]
        bonds = atom.bonds
        del self.map[atomname]

        # Delete the atom from the list

        self.atoms.remove(atom)

        # Delete all instances of the atom as a bond
        
        for bondatom in bonds:
            if atom in bondatom.bonds:
                bondatom.bonds.remove(atom)

        del atom

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

    def getAtoms(self):
        return self.atoms

    def hasAtom(self, name):
        if name in self.map: return 1
        else: return 0

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
            atomcharge = atom.get("ffcharge")
            if atomcharge != None:
                charge = charge + atomcharge

        charge = float("%.4f" % charge)
        return charge

    def renameResidue(self, name):
        """
            Rename a given residue

            Parameters
                name:       The new name of the residue
        """
        self.name = name
        for atom in self.atoms:
            atom.resName = name

    def rotateTetrahedral(self, atom1, atom2, angle):
        """
            Rotate about the atom1-atom2 bond by a given angle
            All atoms connected to atom2 will rotate.

            Parameters:
                atom1:  The first atom of the bond to rotate about (atom)
                atom2:  The second atom of the bond to rotate about (atom)
                angle:  The number of degrees to rotate (float)
        """
        moveatoms = []
        movecoords = []
        
        initcoords = subtract(atom2.getCoords(), atom1.getCoords())

        # Determine which atoms to rotate
        
        for atom in atom2.bonds:
            if atom == atom1: continue
            moveatoms.append(atom)
            movecoords.append(subtract(atom.getCoords(), atom1.getCoords()))

        newcoords = qchichange(initcoords, movecoords, angle)
        for i in range(len(moveatoms)):
            atom = moveatoms[i]
            x = (newcoords[i][0] + atom1.get("x"))
            y = (newcoords[i][1] + atom1.get("y"))
            z = (newcoords[i][2] + atom1.get("z"))
            atom.set("x", x)
            atom.set("y", y)
            atom.set("z", z)
           

    def setDonorsAndAcceptors(self):
        """
            Set the donors and acceptors within the residue
        """
        if not hasattr(self, "reference"): return
        for atom in self.getAtoms():
            atomname = atom.get("name")
            resname = self.name

            atom.set("hdonor", 0)
            atom.set("hacceptor", 0)
         
            if atomname.startswith("N"):
                bonded = 0
                for bondedatom in atom.bonds:
                    if bondedatom.isHydrogen():
                        atom.set("hdonor",1)
                        bonded = 1
                        break
                if not bonded and self.reference.name == "HIS":
                    atom.set("hacceptor",1)
                    
            elif atomname.startswith("O") or \
                 (atomname.startswith("S") and self.reference.name == "CYS"):
                atom.set("hacceptor",1)
                for bondedatom in atom.bonds:
                    if bondedatom.isHydrogen():
                        atom.set("hdonor",1)
                        break     

    def reorder(self):
        """
            Reorder the atoms to start with N, CA, C, O if they exist
        """
        templist = []
        if self.hasAtom("N"): templist.append(self.getAtom("N"))
        if self.hasAtom("CA"): templist.append(self.getAtom("CA"))
        if self.hasAtom("C"): templist.append(self.getAtom("C"))
        if self.hasAtom("O"): templist.append(self.getAtom("O"))

        # Add remaining atoms
        for atom in self.atoms:
            if atom.name not in ["N", "CA", "C", "O"]:
                templist.append(atom)

        # Change the list pointer

        self.atoms = templist[:]

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
            raise ValueError, "Invalid atom type %s (Atom Class IN structures.py)!"
        self.serial = atom.serial
        self.name = atom.name
        self.altLoc = ""
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

        self.bonds = []
        self.reference = None
        self.residue = residue
        self.radius = None
        self.ffcharge = None
        self.hdonor = 0
        self.hacceptor = 0
        self.cell = None
        self.added = 0
        self.optimizeable = 0
        self.refdistance = 0
          
    def __str__(self):
        """
            Returns a string of the new atom type.  Uses the ATOM string
            output but changes the first field to either by ATOM or
            HETATM as necessary.

            Returns
                str: String with ATOM/HETATM field set appropriately
        """
        str = ""
        tstr = self.type
        str = str + string.ljust(tstr, 6)[:6]
        tstr = "%d" % self.serial
        str = str + string.rjust(tstr, 5)[:5]
        str = str + " "
        tstr = self.name
        if len(tstr) == 4:
            str = str + string.ljust(tstr, 4)[:4]
        else:
            str = str + " " + string.ljust(tstr, 3)[:3]

        tstr = self.resName
        if len(tstr) == 4:
            str = str + string.ljust(tstr, 4)[:4]
        else:
            str = str + " " + string.ljust(tstr, 3)[:3]
            
        str = str + " "
        tstr = self.chainID
        str = str + string.ljust(tstr, 1)[:1]
        tstr = "%d" % self.resSeq
        str = str + string.rjust(tstr, 4)[:4]
        str = str + "    "
        tstr = "%8.3f" % self.x
        str = str + string.ljust(tstr, 8)[:8]
        tstr = "%8.3f" % self.y
        str = str + string.ljust(tstr, 8)[:8]
        tstr = "%8.3f" % self.z
        str = str + string.ljust(tstr, 8)[:8]
        if self.ffcharge != None: ffcharge = "%.4f" % self.ffcharge
        else: ffcharge = "0.0000"
        str = str + string.rjust(ffcharge, 8)[:8]
        if self.radius != None: ffradius = "%.4f" % self.radius
        else: ffradius = "0.0000"
        str = str + string.rjust(ffradius, 7)[:7]
        return str
    
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

    def addBond(self, bondedatom):
        """
            Add a bond to the list of bonds

            Parameters:
                bondedatom: The atom to bond to (Atom)
        """
        self.bonds.append(bondedatom)

    def isHydrogen(self):
        """
            Is this atom a Hydrogen atom?

            Returns
                value: 1 if Atom is a Hydrogen, 0 otherwise
        """
        value = 0
        if self.name[0] == "H": value = 1
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

    def hasReference(self):
        """
            Determine if the atom object has a reference object or not.
            All known atoms should have reference objects.

            Returns
                1 if atom has a reference object, 0 otherwise.
        """

        if self.reference != None: return 1
        else: return 0

    
        
