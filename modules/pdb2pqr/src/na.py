"""
    Nucleic Acid Structures for PDB2PQR

    This module contains the base nucleic acid structures for
    pdb2pqr.

    ----------------------------
   
    PDB2PQR -- An automated pipeline for the setup, execution, and analysis of
    Poisson-Boltzmann electrostatics calculations

    Copyright (c) 2002-2011, Jens Erik Nielsen, University College Dublin; 
    Nathan A. Baker, Battelle Memorial Institute, Developed at the Pacific 
    Northwest National Laboratory, operated by Battelle Memorial Institute, 
    Pacific Northwest Division for the U.S. Department Energy.; 
    Paul Czodrowski & Gerhard Klebe, University of Marburg.

	All rights reserved.

	Redistribution and use in source and binary forms, with or without modification, 
	are permitted provided that the following conditions are met:

		* Redistributions of source code must retain the above copyright notice, 
		  this list of conditions and the following disclaimer.
		* Redistributions in binary form must reproduce the above copyright notice, 
		  this list of conditions and the following disclaimer in the documentation 
		  and/or other materials provided with the distribution.
        * Neither the names of University College Dublin, Battelle Memorial Institute,
          Pacific Northwest National Laboratory, US Department of Energy, or University
          of Marburg nor the names of its contributors may be used to endorse or promote
          products derived from this software without specific prior written permission.

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
	ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
	WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
	IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
	INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
	BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
	DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
	LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
	OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED 
	OF THE POSSIBILITY OF SUCH DAMAGE.

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

class ADE(Nucleic):
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
     
class CYT(Nucleic):
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
        
class GUA(Nucleic):
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

class THY(Nucleic):
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
                
class URA(Nucleic):
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

