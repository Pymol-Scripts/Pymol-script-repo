"""
    Definitions for PDB2PQR

    This file contains classes associated with Amino Acid and Rotamer
    definitions as used by PDB2PQR.

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

"""
    
__date__ = "15 May 2008"
__author__ = "Jens Erik Nielsen, Todd Dolinsky, Yong Huang"

AAPATH = "dat/AA.xml"
NAPATH = "dat/NA.xml"
PATCHPATH = "dat/PATCHES.xml"

import os
import copy
import re
from xml import sax
from pdb import *
from utilities import *
from structures import *
from routines import *

class DefinitionHandler(sax.ContentHandler):
   
    def __init__(self):
        self.curelement = ""
        self.curatom = None
        self.curholder = None
        self.curobj = None
        self.map = {}
        self.patches = []
        return

    def startElement(self, name, attributes):
        if name == "residue":
            obj = DefinitionResidue()
            self.curholder = obj
            self.curobj = obj
        elif name == "patch":
            obj = Patch()
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
        if name == "residue": # Complete Residue object
            residue = self.curholder
            if not isinstance(residue, DefinitionResidue):
                raise ValueError, "Internal error parsing XML!"
            resname = residue.name
            if resname == "":
                raise ValueError, "Residue name not set in XML!"
            else:
                self.map[resname] = residue
                self.curholder = None
                self.curobj = None

        elif name == "patch": # Complete patch object
            patch = self.curholder
            if not isinstance(patch, Patch):
                raise ValueError, "Internal error parsing XML!"
            patchname = patch.name
            if patchname == "":
                raise ValueError, "Residue name not set in XML!"
            else:
                self.patches.append(patch)
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
        if text.isspace(): return

        # If this is a float, make it so
        try:
            value = float(str(text))
        except ValueError:
            value = str(text)

        # Special cases - lists and dictionaries
        if self.curelement == "bond":
            self.curobj.bonds.append(value)
        elif self.curelement == "dihedral":
            self.curobj.dihedrals.append(value)
        elif self.curelement == "altname":
            self.curholder.altnames[value] = self.curatom.name
        elif self.curelement == "remove":
            self.curobj.remove.append(value)
        else:
            setattr(self.curobj, self.curelement, value)
        return

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
        self.map = {}
        self.patches = {}

        handler = DefinitionHandler()
        sax.make_parser()

        for path in [AAPATH, NAPATH]:
            defpath = getDatFile(path)
            if defpath == "":
                raise ValueError, "%s not found!" % path

            file = open(defpath)
            sax.parseString(file.read(), handler)
            file.close()

            self.map.update(handler.map)
    
        # Now handle patches

        defpath = getDatFile(PATCHPATH)
        if defpath == "":
            raise ValueError, "%s not found!" % PATCHPATH
     
        handler.map = {}
        file = open(defpath)
        sax.parseString(file.read(), handler)
        file.close()

        # Apply specific patches to the reference object, allowing users
        #  to specify protonation states in the PDB file
        
        for patch in handler.patches:
            if patch.newname != "":

                # Find all residues matching applyto

                resnames = self.map.keys()
                for name in resnames:
                    regexp = re.compile(patch.applyto).match(name)
                    if not regexp: continue
                    newname = patch.newname.replace("*", name)
                    self.addPatch(patch, name, newname)
                  

            # Either way, make sure the main patch name is available
            
            self.addPatch(patch, patch.applyto, patch.name)

    def addPatch(self, patch, refname, newname):
        """
            Add a patch to a definition residue.

            Parameters
                patch:  The patch object to add (Patch)
                refname:  The name of the object to add the patch to (string)
                newname:  The name of the new (patched) object (string)
        """
        try:
            aadef = self.map[refname] # The reference
            patchResidue = copy.deepcopy(aadef)

            # Add atoms from patch

            for atomname in patch.map:
                patchResidue.map[atomname] = patch.map[atomname]
                for bond in patch.map[atomname].bonds:
                    if bond not in patchResidue.map: continue
                    if atomname not in patchResidue.map[bond].bonds:
                        patchResidue.map[bond].bonds.append(atomname)

            # Rename atoms as directed
            
            for key in patch.altnames:
                patchResidue.altnames[key] = patch.altnames[key]

            # Remove atoms as directed
                    
            for remove in patch.remove:
                if not patchResidue.hasAtom(remove): continue
                removebonds = patchResidue.map[remove].bonds
                del patchResidue.map[remove]
                for bond in removebonds:
                    if remove in patchResidue.map[bond].bonds:
                        patchResidue.map[bond].bonds.remove(remove)

            # Add the new dihedrals

            for dihedral in patch.dihedrals:
                patchResidue.dihedrals.append(dihedral)

            # Point at the new reference

            self.map[newname] = patchResidue

            # Store the patch

            self.patches[newname] = patch
                   
        except KeyError: # Just store the patch
            self.patches[newname] = patch

class Patch:
    """
        Patch the definitionResidue class
    """
    def __init__(self):
        """
            Initialize the Patch object.
        """
        self.name = ""
        self.applyto = ""
        self.map = {}
        self.remove = []
        self.altnames = {}
        self.dihedrals = []
        self.newname = ""
        
    def __str__(self):
        """
            A basic string representation for debugging
        """
        text = "%s\n" % self.name
        text += "Apply to: %s\n" % self.applyto
        text += "Atoms to add: \n"
        for atom in self.map:
            text += "\t%s\n" % str(self.map[atom])
        text += "Atoms to remove: \n"
        for remove in self.remove:
            text += "\t%s\n" % remove
        text += "Alternate naming map: \n"
        text += "\t%s\n" % self.altnames
        return text
        
class DefinitionResidue(Residue):
    """
        DefinitionResidue class

        The DefinitionResidue class extends the Residue class to allow for a
        trimmed down initializing function.
    """
    def __init__(self):
        """
            Initialize the class using a few parameters

            Parameters:
                name: The abbreviated amino acid name of the DefinitionResidue
        """
        self.name = ""
        self.dihedrals = []
        self.map = {}
        self.altnames = {}
        
    def __str__(self):
        """
            A basic string representation for debugging
        """
        text = "%s\n" % self.name
        text += "Atoms: \n"
        for atom in self.map:
            text += "\t%s\n" % str(self.map[atom])
        text += "Dihedrals: \n"
        for dihedral in self.dihedrals:
            text += "\t%s\n" % dihedral
        text += "Alternate naming map: \n"
        text += "\t%s\n" % self.altnames
        return text

    def addDihedral(self, atom):
        """
            Add the atom to the list of dihedral bonds

            Parameters:
                atom: The atom to be added
        """
        self.dihedralatoms.append(atom)

    def getNearestBonds(self, atomname):
        """
            Parameters
                number:   The number of bonds to get
            Returns
                bonds:    A list of atomnames that are within three bonds of
                          the atom and present in residue (list)
        """
        bonds = []
        lev2bonds = []
        atom = self.map[atomname]
        
        # Get directly bonded (length = 1) atoms
        
        for bondedatom in atom.bonds:
            if bondedatom not in bonds:
                bonds.append(bondedatom)

        # Get bonded atoms 2 bond lengths away
    
        for bondedatom in atom.bonds:
            for bond2 in self.map[bondedatom].bonds:
                if bond2 not in bonds and bond2 != atomname:
                    bonds.append(bond2)
                    lev2bonds.append(bond2)

        # Get bonded atoms 3 bond lengths away

        for lev2atom in lev2bonds:
            for bond3 in self.map[lev2atom].bonds:
                if bond3 not in bonds:
                    bonds.append(bond3)
         
        return bonds

class DefinitionAtom(Atom):
    """
        A trimmed down version of the Atom class
    """
    def __init__(self, name=None, x=None, y=None, z=None):
        """
            Initialize the class
        """
        self.name = name
        self.x = x
        self.y = y
        self.z = z
        if name == None:
            self.name = ""
        if x == None:
            self.x = 0.0
        if y == None:
            self.y = 0.0
        if z == None:
            self.z = 0.0
        self.bonds = []
     
    def __str__(self):
        """
            A basic string representation for debugging
        """
        text = "%s: %.3f %.3f %.3f" % (self.name, self.x, self.y, self.z)
        for bond in self.bonds:
            text += " %s" % bond
        return text
    
    def isBackbone(self):
        """
            Return true if atom name is in backbone, otherwise false

            Returns
                state: 1 if true, 0 if false
        """
        if self.name in BACKBONE: return 1
        else: return 0
