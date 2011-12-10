"""
    Definitions for PDB2PQR

    This file contains classes associated with Amino Acid and Rotamer
    definitions as used by PDB2PQR.

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

"""
    
__date__ = "28 February 2006"
__author__ = "Jens Erik Nielsen, Todd Dolinsky"

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
    def __init__(self):
        """
            Initialize the class
        """
        self.name = ""
        self.x = 0.0
        self.y = 0.0
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
        if self.name in BACKBONE: return 0
        else: return 1
