#!/usr/bin/python2 -O

"""
    Forcefield script

    This module takes a pdblist as input and replaces the occupancy and
    tempfactor fields with charge and radius fields, with values as defined
    by a particular forcefield.  The forcefield structure is modeled off of
    the structures.py file, where each forcefield is considered a chain of
    residues of atoms.

    Todd Dolinsky (todd@ccb.wustl.edu)
    Washington University in St. Louis
"""

__date__ = "6 October 2003"
__author__ = "Todd Dolinsky"

AMBER_FILE = "AMBER.DAT"
CHARMM_FILE = "CHARMM.DAT"
PARSE_FILE = "PARSE.DAT"

import string
import sys
import getopt
import os

class Forcefield:
    """
        Forcefield class

        The forcefield class contains definitions for a given forcefield.
        Each forcefield object contains a dictionary of residues, with each
        residue containing a dictionary of atoms.  Dictionaries are used
        instead of lists as the ordering is not important. The forcefield
        definition files are unedited, directly from the forcefield - all
        transformations are done within.

    """

    def __init__(self, ff):
        """
            Initialize the class by parsing the definition file

            Parameters
                ff: The name of the forcefield (string)
        """
        self.residues = {}
        self.name = ff

        defpath = ""
        if ff == "amber":
            defpath = AMBER_FILE
        elif ff == "charmm":
            defpath = CHARMM_FILE
        elif ff == "parse":
            defpath = PARSE_FILE
        else:
            raise ValueError, "Invalid forcefield %s!" % ff

        if not os.path.isfile(defpath):
            for path in sys.path:
                testpath = "%s/%s" % (path, defpath)
                if os.path.isfile(testpath):
                    defpath = testpath
                    break
        if not os.path.isfile(defpath):
            raise ValueError, "Unable to find forcefield %s!" % defpath

        file = open(defpath)

        lines = file.readlines()

        for line in lines:
            if not line.startswith("#"):
                fields = string.split(line)
                resname = fields[0]
                atomname = fields[1]
                charge = float(fields[2])
                radius = float(fields[3])

                atom = ForcefieldAtom(atomname, charge, radius)

                myResidue = self.getResidue(resname)
                if myResidue == None:
                    myResidue = ForcefieldResidue(resname)
                    self.residues[resname] = myResidue
                myResidue.addAtom(atom)

    def getResidue(self, resname):
        """
            Return the residue object with the given resname

            Parameters
                resname: The name of the residue (string)
            Returns
                residue: The residue object (ForcefieldResidue)
        """
        residue = None
        try:
            residue = self.residues[resname]
        except KeyError:
            pass
        return residue

    def getParams(self, residue, name):
        """
            Get the parameters associated with the input fields.
            The residue itself is needed instead of simply its name
            because  the forcefield may use a different residue name
            than the standard amino acid name.

            Parameters
                residue:  The residue (residue)
                name:     The atom name (string)
            Returns
                charge:   The charge on the atom (float)
                radius:   The radius of the atom (float)
        """
        charge = None
        radius = None
        resname = ""
        atomname = ""


        if self.name == "amber":
            resname, atomname = self.getAmberParams(residue, name)
        elif self.name == "charmm":
            resname, atomname = self.getCharmmParams(residue, name)
        elif self.name == "parse":
            resname, atomname = self.getParseParams(residue, name)

        defresidue = self.getResidue(resname)
        if defresidue == None:
            return charge, radius

        atom = defresidue.getAtom(atomname)
        if atom != None:
            charge = atom.get("charge")
            radius = atom.get("radius")

        return charge, radius

    def getAmberParams(self, residue, name):
        """
            Get the forcefield definitions from the Amber database

            Parameters
                residue:  The residue (residue)
                name:     The atom name (string)
            Returns
                resname:  The name of the amber residue
                atomname: The name of the amber atom
        """
        atomname = name
        type = residue.get("type")
        if type == 4:
            resname = residue.get("naname")
        else:
            resname = residue.get("name")
        
        # Residue Substitutions
            
        if residue.get("name") == "CYS" and "HG" not in residue.get("map"):
            resname = "CYX"
        elif residue.get("name") == "HIS":
            if "HD1" in residue.get("map") and "HE2" in residue.get("map"):
                resname = "HIP"
            elif "HD1" in residue.get("map"):
                resname = "HID"
            elif "HE2" in residue.get("map"):
                resname = "HIE"
            else:
                resname = "HID" # Default for no hydrogens
        elif residue.get("name") == "HSP":
            resname = "HIP"
        elif residue.get("name") == "HSE":
            resname = "HIE"
        elif residue.get("name") == "HSD":
            resname = "HID"
        elif residue.get("name") == "GLU" or residue.get("name") == "GLH":
            if "HE1" in residue.get("map"):
                resname = "GLH"
                if atomname == "HE1": atomname = "HE2"
                elif atomname == "OE1": atomname = "OE2"
                elif atomname == "OE2": atomname = "OE1"
            elif "HE2" in residue.get("map"): resname = "GLH"
        elif residue.get("name") == "ASP" or residue.get("name") == "ASH":
            if "HD1" in residue.get("map"):
                resname = "ASH"
                if atomname == "HD1": atomname = "HD2"
                elif atomname == "OD1": atomname = "OD2"
                elif atomname == "OD2": atomname = "OD1"
            elif "HD2" in residue.get("map"): resname = "ASH"
               
        if residue.get("isCterm") == 1:
            resname = "C" + resname
        elif residue.get("isNterm") == 1:
            resname = "N" + resname

        # Atom Substitutions
  
        if resname == "WAT":
            if atomname == "O": atomname = "OW"
            elif atomname == "H1": atomname = "HW"
            elif atomname == "H2": atomname = "HW"
        elif resname == "ILE":
            if atomname == "CD": atomname = "CD1"
        if resname[0] == "N": # N-terminal
            if atomname == "H": atomname = "H1"
        if (resname == "CCYS" or resname == "NCYS") and atomname == "HG": atomname = "HSG"
        return resname, atomname

    def getParseParams(self, residue, name):
        """
            Get the forcefield definitions from the Parse database

            Parameters
                residue:  The residue (residue)
                name:     The atom name (string)
            Returns
                resname:  The name of the amber residue
                atomname: The name of the amber atom
        """
        atomname = name
        resname = residue.get("name")

        # Terminal/Water Substitutions

        if residue.get("isNterm") and resname != "ACE":
            if resname == "PRO":
                resname = "PR+"
                if atomname == "H2": atomname = "HN1"
                elif atomname == "H3": atomname = "HN2"
            elif atomname in ["N","H","H2","H3","CA","HA","C","O"]:
                resname = "BK+"
                if atomname == "H": atomname = "H1"
        elif residue.get("isCterm"):
            if atomname in ["N","H","HA","CA","C","O","OXT"]:
                resname = "BK-"
                if atomname == "O": atomname = "O1"
                elif atomname == "OXT": atomname = "O2"

        elif residue.get("type") == 3:
            resname = "H2O"
            if atomname == "O": atomname = "OH"
            elif atomname == "H1": atomname = "HH1"
            elif atomname == "H2": atomname = "HH2"

        # Residue Substitutions
        if resname == "HSD": resname = "HID"
        elif resname in ["HIE","HSE"]: resname = "HIS"
        elif resname in ["HIP","HSP"]: resname = "HI+"
        elif resname == "ILE":
            if atomname == "HG12": atomname = "HG11"
            elif atomname == "HG13": atomname = "HG12"
            elif atomname == "CD": atomname = "CD1"
        elif resname == "CYS" and "HG" not in residue.get("map"):
            resname = "CSS"
        elif resname == "HIS":
            if "HD1" in residue.get("map") and "HE2" in residue.get("map"):
                resname = "HI+"
            elif "HD1" in residue.get("map"):
                resname = "HID"
            elif "HE2" in residue.get("map"):
                resname = "HIS"
        elif resname == "GLU" or resname == "GLH":
            if "HE1" in residue.get("map"):
                resname = "GL0"
                if atomname == "HE1": atomname = "HE2"
                elif atomname == "OE1": atomname = "OE2"
                elif atomname == "OE2": atomname = "OE1"
            elif "HE2" in residue.get("map"): resname = "GL0"
        elif resname == "ASP" or resname == "ASH":
            if "HD1" in residue.get("map"):
                resname = "AS0"
                if atomname == "HD1": atomname = "HD2"
                elif atomname == "OD1": atomname = "OD2"
                elif atomname == "OD2": atomname = "OD1"
            elif "HD2" in residue.get("map"): resname = "AS0"

        # Hydrogen Substitutions

        if atomname == "H": atomname = "HN"
        elif atomname == "HA2": atomname = "HA1"
        elif atomname == "HA3": atomname = "HA2"
        elif atomname == "HB2" and resname not in ["ALA"]: atomname = "HB1"
        elif atomname == "HB3" and resname not in ["ALA"]: atomname = "HB2"
        elif atomname == "HD2" and resname not in ["HIS","HI+","HID"]: atomname = "HD1"
        elif atomname == "HD3" and resname not in ["HIS","HI+","HID"]: atomname = "HD2"
        elif atomname == "HE2" and resname not in ["TRP","HIS","HI+","HID","GL0"]: atomname = "HE1"
        elif atomname == "HE3" and resname not in ["TRP","HIS","HI+","HID"]: atomname = "HE2"
        elif atomname == "HG2": atomname = "HG1"
        elif atomname == "HG3": atomname = "HG2"

        return resname, atomname

    def getCharmmParams(self, residue, name):
        """
            Get the forcefield definitions from the Charmm database

            Parameters
                residue:  The residue (residue)
                name:     The atom name (string)
            Returns
                resname:  The name of the Charmm residue
                atomname: The name of the Charmm atom
        """
        resname = residue.get("name")
        atomname = name

        #  Nucleic Acid Substitutions
        
        if residue.get("type") == 4:
            resname = resname[0]
            if resname == "A": resname = "ADE"
            elif resname == "C": resname = "CYT"
            elif resname == "G": resname = "GUA"
            elif resname == "T":
                resname = "THY"
                if atomname == "C7": atomname = "C5M"
                elif atomname == "H71": atomname = "H51"
                elif atomname == "H72": atomname = "H52"
                elif atomname == "H73": atomname = "H53" 
            elif resname == "U": resname = "URA"

            if atomname == "H5'1": atomname = "H5'"
            elif atomname == "H5'2": atomname = "H5''"
            elif atomname == "H2'1": atomname = "H2'"
            elif atomname in ["H2'2","HO'2"]: atomname = "H2''"
            
            if residue.getAtom("O2'") == None:
                if atomname in ["C2'","H2'","H2''"]: resname = "DEO1"

            if residue.getAtom("H5T") != None:
                if atomname in ["H5T","O5'","C5'"]: resname = "5TER"
            if residue.getAtom("H3T") != None:
                if atomname in ["H3T","O3'","C3'"]: resname = "3TER"
                
        # Terminal/Water Substitutions

        if residue.get("isNterm"):
            if resname == "GLY" and atomname in ["N","H","H2","H3","CA","HA2","HA3"]:
                resname = "GLYP"
                if atomname == "H": atomname = "HT1"
                elif atomname == "H2": atomname = "HT2"
                elif atomname == "H3": atomname = "HT3"
            elif resname == "PRO" and atomname in ["N","HN1","HN2","CD","CA","HD1","HD2","HA","H2","H3"]:
                resname = "PROP"
                if atomname == "H2": atomname = "HN1"
                elif atomname == "H3": atomname = "HN2"
            elif resname == "ACE":
                if atomname == "CH3": atomname = "CAY"
                elif atomname == "HH31": atomname = "HY1"
                elif atomname == "HH32": atomname = "HY2"
                elif atomname == "HH33": atomname = "HY3"
                elif atomname == "C": atomname = "CY"
                elif atomname == "O": atomname = "OY"
            else:
                if atomname in ["N","H","H2","H3","CA","HA"]:
                    resname = "NTER"
                    if atomname == "H": atomname = "HT1"
                    elif atomname == "H2": atomname = "HT2"
                    elif atomname == "H3": atomname = "HT3"               
        elif residue.get("isCterm"):
            if atomname in ["O","OXT","C"]:
                resname = "CTER"
                if atomname == "O":
                    atomname = "OT1"
                elif atomname == "OXT":
                    atomname = "OT2"
        elif residue.get("type") == 3:
            resname = "TP3M"
            if atomname == "O": atomname = "OH2"

        # Residue substitutions
            
        if resname == "ILE":
            if atomname == "CD1": atomname = "CD"
            elif atomname == "HD11": atomname = "HD1"
            elif atomname == "HD12": atomname = "HD2"
            elif atomname == "HD13": atomname = "HD3"
            elif atomname == "HG12": atomname = "HG11"
            elif atomname == "HG13": atomname = "HG12"
        elif resname == "CYS" and "HG" not in residue.get("map"):
            if atomname == "CB":
                resname = "DISU"
                atomname = "1CB"
            elif atomname == "SG":
                resname = "DISU"
                atomname = "1SG"
        elif resname == "HIS":
            if "HD1" in residue.get("map") and "HE2" in residue.get("map"):
                resname = "HSP"
            elif "HD1" in residue.get("map"):
                resname = "HSD"
            elif "HE2" in residue.get("map"):
                resname = "HSE"
        elif resname == "GLU" or resname == "GLH":
            if "HE1" in residue.get("map"):
                if atomname == "HE1": atomname = "HE2"
                elif atomname == "OE1": atomname = "OE2"
                elif atomname == "OE2": atomname = "OE1"
                if atomname in ["CG","HG3","HG1","HG2","CD","OE1","OE2","HE2"]: resname = "GLUP"
                else: resname == "GLU"
            elif "HE2" in residue.get("map"):
                if atomname in ["CG","HG3","HG1","HG2","CD","OE1","OE2","HE2"]: resname = "GLUP"
                else: resname == "GLU"
        elif resname == "ASP" or resname == "ASH":
            if "HD1" in residue.get("map"):
                if atomname == "HD1": atomname = "HD2"
                elif atomname == "OD1": atomname = "OD2"
                elif atomname == "OD2": atomname = "OD1"
                if atomname in ["CB","HB3","HB1","HB2","CG","OD1","OD2","HD2"]: resname = "ASPP"
                else: resname == "ASP"
            elif "HD2" in residue.get("map"):
                if atomname in ["CB","HB3","HB1","HB2","CG","OD1","OD2","HD2"]: resname = "ASPP"
                else: resname == "ASP"
                
        # HETATM Substitutions

        if resname == "ACE":
            if atomname == "CH3": atomname = "CAY"
            elif atomname == "HH31": atomname = "HY1"
            elif atomname == "HH32": atomname = "HY2"
            elif atomname == "HH33": atomname = "HY3"
            elif atomname == "C": atomname = "CY"
            elif atomname == "O": atomname = "OY"
        elif resname == "ADP":
            atomname = string.replace(atomname,"*","\'")
            
        # Hydrogen Substitutions

        if atomname == "H": atomname = "HN"
        elif atomname == "HA2": atomname = "HA1"
        elif atomname == "HA3": atomname = "HA2"
        elif atomname == "HB2" and resname not in ["ALA"]: atomname = "HB1"
        elif atomname == "HB3" and resname not in ["ALA"]: atomname = "HB2"
        elif atomname == "HD2" and resname not in ["HSP","HSE","HSD","ASPP"]: atomname = "HD1"
        elif atomname == "HD3" and resname not in ["HIS","HSE","HSD"]: atomname = "HD2"
        elif atomname == "HE2" and resname not in ["TRP","HSP","HSE","HSD","GLUP"]: atomname = "HE1"
        elif atomname == "HE3" and resname not in ["TRP","HSP","HSE","HSD"]: atomname = "HE2"
        elif atomname == "HG2": atomname = "HG1"
        elif atomname == "HG3": atomname = "HG2"
        elif atomname == "HG" and resname in ["SER","CYS"]: atomname = "HG1"
        
        return resname, atomname

class ForcefieldResidue:
    """
        ForcefieldResidue class

        The ForceFieldResidue class contains a mapping of all atoms within
        the residue for easy searching.
    """
    def __init__(self, name):
        """
            Initialize the ForceFieldResidue object

            Parameters
                name: The name of the residue (string)
        """
        self.name = name
        self.atoms = {}

    def addAtom(self, atom):
        """
            Add an atom to the ForcefieldResidue

            Parameters
                atom:  The atom to be added (atom)
        """
        atomname = atom.get("name")
        self.atoms[atomname] = atom

    def getAtom(self, atomname):
        """
            Return the atom object with the given atomname

            Parameters
                resname: The name of the atom (string)
            Returns
                residue: The atom object (ForcefieldAtom)
        """
        atom = None
        try:
            atom = self.atoms[atomname]
        except KeyError:
            pass
        return atom

class ForcefieldAtom:
    """
        ForcefieldAtom class

        The ForcefieldAtom object contains fields that are related to the
        forcefield at the atom level
    """
    
    def __init__(self, name, charge, radius):
        """
            Initialize the object

            Parameters
                name:    The atom name (string)
                charge:  The charge on the atom (float)
                radius:  The radius of the atom (float)
        """
        self.name = name
        self.charge = charge
        self.radius = radius

    def get(self, name):
        """
            Get a member of the ForcefieldAtom class

            Parameters
                name:       The name of the member (string)
            Possible Values
                name:    The atom name (string)
                charge:  The charge on the atom (float)
                radius:  The radius of the atom (float)
                epsilon: The epsilon assocaited with the atom (float)
            Returns
                item:       The value of the member
        """
        try:
            item = getattr(self, name)
            return item
        except AttributeError:
            message = "Unable to access object \"%s\" in class ForcefieldAtom" % name
            raise ValueError, message
