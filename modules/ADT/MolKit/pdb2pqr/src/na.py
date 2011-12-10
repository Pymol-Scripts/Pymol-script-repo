"""
    Nucleic Acid Structures for PDB2PQR

    This module contains the base nucleic acid structures for
    pdb2pqr.

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

import string
from structures import *

class Nucleic(Residue):
    """
        Nucleic class

        This class provides standard features of the nucleic acids listed
        below

        Parameters
            atoms:  A list of Atom objects to be stored in this class
                     (list)
            ref:    The reference object for the amino acid.  Used to
                    convert from the alternate naming scheme to the
                    main naming scheme.
    """
    def __init__(self, atoms, ref):
        sampleAtom = atoms[-1]
        
        self.atoms = []
        self.name = sampleAtom.resName
        self.chainID = sampleAtom.chainID
        self.resSeq = sampleAtom.resSeq
        self.iCode = sampleAtom.iCode

        self.ffname = self.name
        self.map = {}
        self.dihedrals = []
        self.patches = []
        self.is3term = 0
        self.is5term = 0
        self.isCterm = 0
        self.isNterm = 0
        self.missing = []
        self.reference = ref
     
        # Create each atom

        for a in atoms:
            if a.name in ref.altnames: # Rename atoms
                a.name = ref.altnames[a.name]

            if a.name not in self.map:
                atom = Atom(a, "ATOM", self)
                self.addAtom(atom)

    def createAtom(self, atomname, newcoords):
        """
            Create an atom.  Overrides the generic residue's createAtom().

            Parameters
                atomname:  The name of the atom to add (string)
                newcoords: The coordinates of the atom (list)
        """
        oldatom = self.atoms[0]
        newatom = Atom(oldatom, "ATOM", self)
        newatom.set("x",newcoords[0])
        newatom.set("y",newcoords[1])
        newatom.set("z",newcoords[2])
        newatom.set("name", atomname)
        newatom.set("occupancy",1.00)
        newatom.set("tempFactor",0.00)
        newatom.added = 1
        self.addAtom(newatom) 

    def addAtom(self, atom):
        """
            Override the existing addAtom - include the link to the
            reference object
        """
        self.atoms.append(atom)
        atomname = atom.get("name")
        self.map[atomname] = atom
        try:
            atom.reference = self.reference.map[atomname]
            for bond in atom.reference.bonds:
                if self.hasAtom(bond):
                    bondatom = self.map[bond]
                    if bondatom not in atom.bonds: atom.bonds.append(bondatom)
                    if atom not in bondatom.bonds: bondatom.bonds.append(atom)
        except KeyError:
            atom.reference = None

    def addDihedralAngle(self, value):
        """
            Add the value to the list of chiangles

            Parameters
                value: The value to be added (float)
        """
        self.dihedrals.append(value)

    def setState(self):
        """ 
           Adds the termini for all inherited objects
        """
        if self.is5term: self.ffname = self.ffname + "5"
        if self.is3term: self.ffname = self.ffname + "3"
 
class A(Nucleic):
    """
        Adenosine class

        This class gives data about the Adenosine object, and inherits
        off the base residue class.
    """

    def __init__(self, atoms, ref):
        """
            Initialize the class

            Parameters
                atoms:      A list of Atom objects to be stored in this class
                            (list)
        """
        Nucleic.__init__(self, atoms, ref)
        self.reference = ref

    def setState(self):
        """
            Set the state to distinguish RNA from DNA.
        """
        if self.hasAtom("O2'"): self.ffname = "RA"
        else: self.ffname = "DA"
        Nucleic.setState(self)
     
class C(Nucleic):
    """
        Cytidine class

        This class gives data about the Cytidine object, and inherits
        off the base residue class.
    """

    def __init__(self, atoms, ref):
        """
            Initialize the class

            Parameters
                atoms:      A list of Atom objects to be stored in this class
                            (list)
        """
        Nucleic.__init__(self, atoms, ref)
        self.reference = ref
        
    def setState(self):
        """
            Set the state to distinguish RNA from DNA.
        """
        if self.hasAtom("O2'"): self.ffname = "RC"
        else: self.ffname = "DC"
        Nucleic.setState(self)
        
class G(Nucleic):
    """
        Guanosine class

        This class gives data about the Guanosine object, and inherits
        off the base residue class.
    """

    def __init__(self, atoms, ref):
        """
            Initialize the class

            Parameters
                atoms:      A list of Atom objects to be stored in this class
                            (list)
        """
        Nucleic.__init__(self, atoms, ref)
        self.reference = ref
        
    def setState(self):
        """
            Set the state to distinguish RNA from DNA.
        """
        if self.hasAtom("O2'"): self.ffname = "RG"
        else: self.ffname = "DG"
        Nucleic.setState(self)
        
class T(Nucleic):
    """
        Thymine class

        This class gives data about the Thymine object, and inherits
        off the base residue class.
    """

    def __init__(self, atoms, ref):
        """
            Initialize the class

            Parameters
                atoms:      A list of Atom objects to be stored in this class
                            (list)
        """
        Nucleic.__init__(self, atoms, ref)
        self.reference = ref
        
    def setState(self):
        """
            Set the state to distinguish RNA from DNA.  In this case it is
            always DNA.
        """
        self.ffname = "DT"
        Nucleic.setState(self)
        
class U(Nucleic):
    """
        Uridine class

        This class gives data about the Uridine object, and inherits
        off the base residue class.
    """

    def __init__(self, atoms, ref):
        """
            Initialize the class

            Parameters
                atoms:      A list of Atom objects to be stored in this class
                            (list)
        """
        Nucleic.__init__(self, atoms, ref)
        self.reference = ref
        
    def setState(self):
        """
            Set the state to distinguish RNA from DNA.  In this case it is
            always RNA.
        """
        self.ffname = "RU"
        Nucleic.setState(self)
