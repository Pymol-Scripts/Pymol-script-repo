"""
    Routines for PDB2PQR

    This module contains the protein object used in PDB2PQR and methods
    used to correct, analyze, and optimize that protein.

    Based on C code from Jens Erik Nielsen
    UCSD/HHMI
    
    Ported to Python by Todd Dolinsky (todd@ccb.wustl.edu)
    Washington University in St. Louis

"""

__date__ = "22 October 2003"
__author__ = "Jens Erik Nielsen, Todd Dolinsky"

CELL_SIZE = 2
BUMP_DIST = 2.0
BUMP_HDIST = 1.5
BONDED_SS_LIMIT = 2.5
LARGE_TORSION_ANGLE = 1000.0
PEPTIDE_DIST = 1.7
REPAIR_LIMIT = 10
REFATOM_SIZE = 3
HYDRO_BONDCOORDS = [[7.581,2.090,12.506],[6.458,2.162,13.159],[5.145,2.209,12.453]]
HYDRO_COORDS = [6.476, 2.186, 14.159]
NTERM_COORDS = [[-24.196, 48.790, -20.8], [-25.552, 49.881, -21.848], [-24.645, 49.491, -22.007]]
NTERM2_COORDS = [-24.001, 50.224, -22.226]
NTERM3_COORDS = [-24.869, 48.846, -22.770]
PEP_TRANS_N = [-1.252,1.877,0.883]
PEP_TRANS_CA = [-2.313,2.784,1.023]
OXT_COORDS = [-1.529,1.858,0.695]
AAS = ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLH","GLY","HIS",\
       "HID","HIE","HIP","HSD","HSE","HSP","ILE","LEU","LYS","MET",\
       "PHE","PRO","SER","THR","TRP","TYR","VAL"]
NAS = ["A","A5","A3","C","C5","C3","G","G5","G3","T","T5","T3","U",\
       "U5","U3","RA","RG","RC","RU","DA","DG","DC","DT"]
       

import random
from pdb import *
from utilities import *
from quatfit import *
from forcefield import *
from structures import *
from protein import *
from definitions import *

class Routines:
    def __init__(self, protein, verbose, definition=None):
        """
        """
        self.protein = protein
        self.definition = definition
        self.aadef = None
        self.verbose = verbose
        self.warnings = []
        self.cells = {}
        if definition != None:
            self.aadef = definition.getAA()
            self.nadef = definition.getNA()
            
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
        for chain in self.protein.getChains():
            for residue in chain.get("residues"):
                for atom in residue.get("atoms"):
                    atomname = atom.get("name")
                    charge, radius = forcefield.getParams(residue, atomname)
                    if charge != None and radius != None:
                        atom.set("ffcharge", charge)
                        atom.set("radius", radius)
                        hitlist.append(atom)
                    else:
                        misslist.append(atom)  
        self.write("Done.\n")            
        return hitlist, misslist

    def updateResidueTypes(self):
        """
            Find the type of residue as notated in the Amino Acid definition
        """
        self.write("Updating Residue Types... ")
        for chain in self.protein.getChains():
            for residue in chain.get("residues"):
                name = residue.get("name")
                if name in AAS:
                    residue.set("type",1)
                elif name == "WAT":
                    residue.set("type",3)
                elif name in NAS:
                    residue.set("type",4)
                else: # Residue is a ligand or unknown
                    residue.set("type",2)
                 
        self.write("Done\n")
            
    def updateSSbridges(self):
        """
            Check for SS-bridge partners, and if present, set appropriate
            partners
        """
        self.write("Updating SS bridges...\n")
        SGatoms = []
        SGpartners = []
        for atom in self.protein.getAtoms():
            if atom.name == "SG":
                SGatoms.append(atom)
                SGpartners.append([])

        for i in range(len(SGatoms)):
            for j in range(len(SGatoms)):
                dist = distance(SGatoms[i].getCoords(), SGatoms[j].getCoords())
                if i != j and dist < BONDED_SS_LIMIT:
                    SGpartners[i].append(j)
        
        for i in range(len(SGatoms)):
            res1 = SGatoms[i].get("residue")
            if len(SGpartners[i]) == 1:
                partner = SGpartners[i][0]
                if SGpartners[partner][0] == i:
                    res2 = SGatoms[partner].get("residue")
                    if i<partner:
                        self.write("CYS %4d - CYS %4d\n" % \
                                   (res1.get("resSeq"), res2.get("resSeq")), 1)
                    if res1.get("name") == "CYS":
                        res1.set("SSbonded", 1)
                        res1.set("SSbondpartner", SGatoms[partner])
                    else:
                        name = res1.get("name")
                        num = res1.get("resSeq")
                        error = "Tried to set SS bonding "
                        error += "for CYS %i, but the residue is a %s.  " % \
                                 (num, name)
                        error += "This should not occur - please contact the "
                        error += "author to report this bug."
                        raise ValueError, error
                
                else:
                    raise ValueError, "CYS %i unresolved!" % res1.get("resSeq")
            elif len(SGpartners[i]) > 1:
                error = "CYS %i has multiple potential " % res1.get("resSeq")
                error += "SS-bridge partners - PDB2PQR is unable to continue."
                raise ValueError, error
            elif len(SGpartners[i]) == 0:
                self.write("CYS %4d is a free cysteine\n" % res1.get("resSeq"), 1)
        self.write("Done.\n")

    def calculateChiangles(self):
        """
            Calculate the dihedral angle for every residue within the protein,
            using the Amino Acid definition.
        """
        self.write("Calculating all chiangles... ")
        for chain in self.protein.getChains():
            for residue in chain.get("residues"):
                residue.set("chiangles",[])
                name = residue.get("name")
                type = residue.get("type")
                definitionres = None
                if type == 1 or type == 3:
                    definitionres = self.aadef.getResidue(name)
                elif type == 4:
                    name = residue.get("naname")
                    definitionres = self.nadef.getResidue(name)
                if definitionres != None:
                    defdihedrals = definitionres.get("dihedralatoms")
                    for i in range(0, len(defdihedrals), 4):       
                        atom1 = residue.getAtom(defdihedrals[i])
                        atom2 = residue.getAtom(defdihedrals[i+1])
                        atom3 = residue.getAtom(defdihedrals[i+2])
                        atom4 = residue.getAtom(defdihedrals[i+3])
                        
                        if atom1 == None or atom2 == None \
                               or atom3 == None or atom4 == None:
                            residue.addChiangle(LARGE_TORSION_ANGLE)
                        else:
                            residue.addChiangle(getDihedral(atom1.getCoords(),\
                                                            atom2.getCoords(),\
                                                            atom3.getCoords(),\
                                                            atom4.getCoords()))
             
                else:
                    if residue.get("type") != 2:
                        error = "Unable to find Amino Acid definition for "
                        error += "%s!" % name
                        raise ValueError, error
        self.write("Done.\n")

    def updateExtraBonds(self):
        """
            Update peptide bonds between amino acids. Set the Termini.
        """
        self.write("Determining peptide bonds and termini... \n")
        for chain in self.protein.getChains():
            for i in range(chain.numResidues() - 1):
                residue1 = chain.get("residues")[i]
                residue2 = chain.get("residues")[i+1]
                res1type = residue1.get("type")
                res2type = residue2.get("type")
                
                if res1type == 1 and res2type == 1: 
                    atom1 = residue1.getAtom("C")
                    atom2 = residue2.getAtom("N")
                    if atom1 != None and atom2 != None:
                        if distance(atom1.getCoords(),
                                    atom2.getCoords()) < PEPTIDE_DIST:
                            atom1.addExtraBond(atom2)
                            atom2.addExtraBond(atom1)
                        else:
                            self.write("Gap in backbone detected in chain ",1)
                            self.write("%s between %s " % (chain.get("chainID"), \
                                                           residue1.get("name")))
                            self.write("%s and %s %s\n"%(residue1.get("resSeq"),\
                                                         residue2.get("name"),\
                                                         residue2.get("resSeq")))

                """ Set the appropriate termini """
                
          
                if res1type == 1 and i == 0:
                    residue1.set("isNterm",1)
                elif res1type == 1 and res2type != 1 and residue2.get("name") not in ["ACE","HMS"]:
                    # Check to make sure this is the last AA in the chain
                    if (i+2) > (chain.numResidues() - 1):
                         residue1.set("isCterm",1)
                    cterm = 1
                    for j in range(i+2, chain.numResidues()):
                        if chain.get("residues")[j].type == 1:
                            cterm = 0
                            break
                    if cterm == 1:
                        residue1.set("isCterm",1)
                elif res2type == 1 and i+2 == chain.numResidues():
                    residue2.set("isCterm",1)
                elif res1type == 4 and i == 0:
                    residue1.set("is5term",1) 
                elif res2type == 4 and i+2 == chain.numResidues():
                    residue2.set("is3term",1)
                elif res1type == 4 and residue1.getAtom("H3T") != None:      
                    residue1.set("is3term",1)
                    if res2type == 4: residue2.set("is5term",1)
                elif res2type == 4 and residue2.getAtom("H5T") != None:
                    residue2.set("is5term",1)
                    if res1type == 4: residue1.set("is3term",1)
                    
        self.write("Done.\n")

    def updateIntraBonds(self):
        """
            Update the bonds within a residue of the protein
        """
        for chain in self.protein.getChains():
            for residue in chain.get("residues"):
                type = residue.get("type")
                name = residue.get("name")
                if type == 1 or type == 3:
                    defresidue = self.aadef.getResidue(name)
                elif type == 4:
                    name = residue.get("naname")
                    defresidue = self.nadef.getResidue(name)
                else: continue
                if defresidue == None:
                    error = "Could not find definition for %s " % name
                    error += "even though it is type %i!" % type
                    raise ValueError, error
                residue.updateIntraBonds(defresidue)

    def correctNames(self):
        """
            Correct atom names so that they match those listed in
            the amino acid definition.  Handles C-Terminal Oxygens
            and various Hydrogen naming schemes.
        """
        self.write("Correcting all atom names... ")
        for chain in self.protein.getChains():
            for residue in chain.get("residues"):
                resname = residue.get("name")
                resnum = residue.get("resSeq")
                if residue.get("type") == 1:
                    residue.checkAtomNames()

                    if resname == "ILE":
                        atom = residue.getAtom("CD")
                        if atom !=  None:
                            residue.renameAtom("CD", "CD1")
           
                    if residue.get("isCterm") == 1:
                        for atomname in ["OT1","O\'"]:
                            atom = residue.getAtom(atomname)
                            if atom != None:
                                residue.renameAtom(atomname, "O")
                                self.write("\n")
                                self.write("Renaming %s to O" % atomname,1)
                                
                        for atomname in ["OT2", "O\'\'"]:
                            atom = residue.getAtom(atomname)
                            if atom != None:
                                residue.renameAtom(atomname, "OXT")
                                self.write("\n")
                                self.write("Renaming %s to OXT\n" % atomname,1)
                                
                    elif residue.get("isNterm") == 1:
                        if residue.getAtom("HT1") != None:
                            residue.renameAtom("HT1", "H")
                            self.write("\n")
                            self.write("Renaming HT1 to H",1)                  
                        if residue.getAtom("HT2") != None:
                            residue.renameAtom("HT2", "H2")
                            self.write("\n")
                            self.write("Renaming HT2 to H2",1)
                        if residue.getAtom("HT3") != None:
                            residue.renameAtom("HT3", "H3")
                            self.write("\n")
                            self.write("Renaming HT3 to H3",1)
                            
                elif residue.get("type") == 2:
                    if resname in ["ACE"]: # Acetyl N-Terminus
                        residue.checkAtomNames()
                        
                elif residue.get("type") == 3:
                    residue.checkAtomNames()
                    name = residue.get("name")                 
                    for atomname in ["OH2"]:
                        atom = residue.getAtom(atomname)
                        if atom != None:
                            residue.renameAtom(atomname, "O")
                    if residue.getAtom("O") == None:
                        error = "\tCannot Repair Water when " \
                                "Oxygen is missing!: See %s %i\n" % \
                                (resname, resnum)
                        raise ValueError, error

                elif residue.get("type") == 4:
                    id = ""

                    # Perform 3 Atom Naming Scheme Checks:
                    #   1. Replace all * with '
                    #   2. Convert 1H2' to H2'1
                    #   3. Replace all 5M with 7
                    
                    for atom in residue.get("atoms"):
                        atomname = atom.get("name")
                        newname = string.replace(atomname,"*","'")
                        if atomname != newname:
                            residue.renameAtom(atomname,newname)

                        try:
                            atomname = atom.get("name")
                            firstint = int(atomname[0])
                            newname = atomname[1:] + atomname[0]
                            residue.renameAtom(atomname,newname)
                        except ValueError: pass

                        atomname = atom.get("name")
                        if string.find(atomname,"5M") != -1:
                            newname = string.replace(atomname,"5M","7")
                            residue.renameAtom(atomname,newname)
                            self.write("Renaming %s to %s" % (atomname, newname),1)
                            self.write("\n")                     
           
                    # Determine if this is DNA/RNA and Definition name
                    
                    rna = 0
                    dna = 0 
                    name = residue.get("name")
                    if name[0] == "R" or name[0] == "D": name = name[1:]
                    if residue.getAtom("O2\'") != None: rna = 1
                    else: dna = 1

                    if rna and not dna: id = "R"
                    elif not rna and dna: id = "D"
                    else:
                        text = "Nucleic Acid %s %i" % (name, residue.resSeq)
                        text += "was found to be both DNA and RNA!"
                        raise ValueError, text
                    id += name[0]
                    if residue.get("is3term"): id += "3"
                    elif residue.get("is5term"): id += "5"
                    residue.set("naname",id)
                            
                else:   #residue is an unknown type
                    raise ValueError, "Unknown residue type!"

        self.write("Done.\n")

    def findMissingHeavy(self):
        """
            Repair residues that contain missing heavy (non-Hydrogen) atoms
        """
        self.write("Checking for missing heavy atoms... \n")
        misscount = 0
        for chain in self.protein.getChains():
            for residue in chain.get("residues"):
                resSeq = residue.get("resSeq")
                type = residue.get("type")
                if type == 1:
                    name = residue.get("name")
                    defresidue = self.aadef.getResidue(name)
                    if defresidue == None:
                        error = "Could not find definition for %s " % name
                        error += "even though it is an amino acid!"
                        raise ValueError, error
                    
                    # Check for Missing Heavy Atoms
                    
                    for defatom in defresidue.get("atoms"):
                        if not defatom.isHydrogen():
                            defname = defatom.get("name")
                            atom = residue.getAtom(defname)
                            if atom == None:
                                self.write("Missing %s in %s %i\n" % \
                                           (defname, name, resSeq), 1)
                                misscount += 1
                                residue.addMissing(defname)

                    if residue.get("isCterm") == 1:
                        atom = residue.getAtom("OXT")
                        if atom == None:
                            residue.addMissing("OXT")
                            misscount += 1
                            self.write("Missing OXT in %s %i\n" % (name, resSeq),1)
                            
                    # Check for Extra Atoms

                    atomlist = []
                    for atom in residue.get("atoms"):
                        atomlist.append(atom)

                    for atom in atomlist:
                        atomname = atom.get("name")
                        defatom = defresidue.getAtom(atomname)
                        if atomname == "OXT" and residue.get("isCterm"):
                            pass
                        elif atom.isHydrogen() and residue.get("isNterm"):
                            pass
                        elif defatom == None:
                            self.write("Extra atom %s in %s %i! - " % \
                                       (atomname, name, resSeq), 1)
                            residue.removeAtom(atomname)
                            self.write("Deleted this atom.\n")

                elif type == 4:
                    name = residue.get("naname")
                    defresidue = self.nadef.getResidue(name)
                    if defresidue == None:
                        error = "Could not find definition for %s " % name
                        error += "even though it is a nucleic acid!"
                        raise ValueError, error
                    
                    # Check for Missing Heavy Atoms
                    
                    for defatom in defresidue.get("atoms"):
                        if not defatom.isHydrogen():
                            defname = defatom.get("name")
                            atom = residue.getAtom(defname)
                            if atom == None:
                                resname = residue.get("name")
                                self.write("Missing %s in %s %i\n" % \
                                           (defname, resname, resSeq), 1)
                                misscount += 1
                                residue.addMissing(defname)
                            
                    # Check for Extra Atoms

                    atomlist = []
                    for atom in residue.get("atoms"):
                        atomlist.append(atom)

                    for atom in atomlist:
                        atomname = atom.get("name")
                        defatom = defresidue.getAtom(atomname)
                        if defatom == None:
                            self.write("Extra atom %s in %s %i! - " % \
                                       (atomname, name, resSeq), 1)
                            residue.removeAtom(atomname)
                            self.write("Deleted this atom.\n")
                            
        numatoms = self.protein.numAtoms()
        misspct = float(misscount) / (numatoms + misscount)
        if misspct > REPAIR_LIMIT / 100.0:
            error = "This PDB file is missing too many (%i out of " % misscount
            error += "%i, %i%%) heavy atoms to accurately repair the file.  " % \
                     ((numatoms + misscount), int(misspct*100))
            error += "The current repair limit is set at %i%%." % REPAIR_LIMIT
            raise ValueError, error
        elif misscount > 0:
            self.write("Missing %i out of %i heavy atoms (%.2f percent) - " %\
                       (misscount, (numatoms + misscount), (misspct*100)),1)
            self.write("Will attempt to repair.\n")
            sgadded = self.repairHeavy()
            if sgadded:
                self.updateSSbridges()
        else:
            self.write("No heavy atoms found missing - Done.\n")

    def rebuildMethyl(self, atomname, residue, defresidue):
        """
            Rebuild the final methyl hydrogen atom using equations from its
            tetrahedral geometry.  The normal quaternion/reference frame
            method does NOT work, since
                A.  It only works with 3 reference atoms, and
                B.  For methyl hydrogens, 3 reference atoms provide TWO
                    possible locations, where only one is potentially correct.

            Parameters
                atomname:   The name of the atom to rebuild (string)
                residue:    The residue (residue)
                defresidue: The definition residue (definitionResidue)
            Returns
                coords :  The new coords of the atom (list)
        """
        hyds = []
        rads = 109.5*math.pi/180.0
        bondname = None
        restname = None
        
        if not atomname.startswith("H"): return None
    
        # Get bonded atom that is one bond length away
    
        defatom = defresidue.getAtom(atomname)
        bondname = defatom.get("intrabonds")[0]
        if bondname not in residue.get("map"):
            return None
    
        # In methyl groups there are four atoms bonded to the bondatom -
        # the three hydrogens and the atom the bondatom is bonded to.
        
        bonds = defresidue.getAtom(bondname).get("intrabonds")
        if len(bonds) != 4: return None

        for bond in bonds:
            if bond.startswith("H") and bond in residue.get("map"):
                hyds.append(bond)
            elif restname == None:
                restname = bond
            elif bond.startswith("H"): pass
            else: return None

        if len(hyds) != 2: return None

        # We now have a methyl group - do the matrix math
        rows = []
        b = []
    
        bondatom = residue.getAtom(bondname)
        restatom  = residue.getAtom(restname)
        bonddist = distance(bondatom.getCoords(),residue.getAtom(hyds[0]).getCoords())
        restdist = distance(bondatom.getCoords(),restatom.getCoords())

        for hyd in hyds:
            hatom = residue.getAtom(hyd)
            lhs = cos(rads) * bonddist * bonddist
            rhs = bondatom.x*bondatom.x + bondatom.y*bondatom.y + bondatom.z*bondatom.z
            rhs = rhs - bondatom.x*hatom.x - bondatom.y*hatom.y - bondatom.z*hatom.z
            rhs = rhs - lhs
            rows.append([bondatom.x-hatom.x, bondatom.y-hatom.y, bondatom.z-hatom.z])
            b.append(rhs)
            
        lhs = cos(rads) * bonddist * restdist
        rhs = bondatom.x*bondatom.x + bondatom.y*bondatom.y + bondatom.z*bondatom.z
        rhs = rhs - bondatom.x*restatom.x - bondatom.y*restatom.y - bondatom.z*restatom.z
        rhs = rhs - lhs
        rows.append([bondatom.x-restatom.x, bondatom.y-restatom.y, bondatom.z-restatom.z])
        b.append(rhs)
        
        mat = Matrix(rows)
        return mat.LU(b)

    def addHydrogens(self):
        """
            Add hydrogens to the residue by using the definition.
        """
        count = 0
        self.write("Adding hydrogens to the protein...")
        for chain in self.protein.getChains():
            prevres = None
            for residue in chain.get("residues"):
                name = residue.get("name")
                type = residue.get("type")
                if type == 1:
                    if residue.get("isNterm") or len(chain.get("residues")) == 1:
                        prevres = residue
                    defresidue = self.aadef.getResidue(name)
                    for defatom in defresidue.get("atoms"):
                        refcoords = []
                        defcoords = []
                        if not defatom.isHydrogen(): continue
                        defname = defatom.get("name")
                        atom = residue.getAtom(defname)
                        if atom != None: continue
                        if atom == None and name == "HSN":
                            if defname == "HD1" and residue.getAtom("HE2"): continue
                            if defname == "HE2" and residue.getAtom("HD1"): continue
                        prevC = prevres.getAtom("C")

                        # For most backbone Hs, use the previous C atom and this residue's
                        #  N and CA atoms
                        
                        if defname == "H" and not residue.get("isNterm") and prevC != None:
                            refcoords.append(prevC.getCoords())
                            refcoords.append(residue.getAtom("N").getCoords())
                            refcoords.append(residue.getAtom("CA").getCoords())
                            defcoords = HYDRO_BONDCOORDS
                            defatomcoords = HYDRO_COORDS
                            newcoords = findCoordinates(REFATOM_SIZE, refcoords, defcoords, defatomcoords)
                            residue.createAtom(defname, newcoords, "ATOM")
                            residue.addDebumpAtom(residue.getAtom(defname))
                            count += 1
                        elif defname == "H" and residue.get("isNterm"): continue
                        elif residue.get("SSbonded") and defname == "HG": continue
                        else:
                            newcoords = self.rebuildMethyl(defname, residue, defresidue)
                            if newcoords != None:
                                residue.createAtom(defname, newcoords,"ATOM")
                                residue.addDebumpAtom(residue.getAtom(defname))
                                count += 1
                                continue
                            bonds = defresidue.makeBondList(residue,defname)
                            if len(bonds) < REFATOM_SIZE:
                                error = "Not enough bonds to remake hydrogen in %s %i" % \
                                        (name, residue.get("resSeq"))
                                raise ValueError, error
                            for i in range(REFATOM_SIZE):
                                refcoords.append(residue.getAtom(bonds[i]).getCoords())
                                defcoords.append(defresidue.getAtom(bonds[i]).getCoords())
                                defatomcoords = defresidue.getAtom(defname).getCoords()
                            newcoords = findCoordinates(REFATOM_SIZE, refcoords, defcoords, defatomcoords)
                            residue.createAtom(defname, newcoords, "ATOM")
                            residue.addDebumpAtom(residue.getAtom(defname))
                            count += 1
                        
                    # Add N-Terminal Hydrogens if Necessary
                    nterm = residue.get("isNterm")
                    cterm = residue.get("isCterm")
                    if nterm > 0:

                        if name != "PRO": hname = "H"
                        else: hname = "HA"

                        # First add the H at tetrahedral geometry
                        # See hydrogens.py for locations
                        if hname not in residue.map:
                            refcoords = []
                            defcoords = []
                            refcoords.append(residue.getAtom("N").getCoords())
                            refcoords.append(residue.getAtom("CA").getCoords())
                            defcoords.append([0,0,.3333])
                            defcoords.append([0,0,1.7963])
                            defatomcoords = [0.9428,0,0]
                            newcoords = findCoordinates(2, refcoords, defcoords, defatomcoords)
                            residue.createAtom(hname, newcoords, "ATOM")
                            residue.addDebumpAtom(residue.getAtom(hname))
                            count += 1

                        # Now add H2

                        refcoords = []
                        refcoords.append(residue.getAtom("CA").getCoords())
                        refcoords.append(residue.getAtom(hname).getCoords())
                        refcoords.append(residue.getAtom("N").getCoords())
                        defcoords = NTERM_COORDS

                        if nterm >= 2 and "H2" not in residue.get("map"):
                            defatomcoords = NTERM2_COORDS
                            newcoords = findCoordinates(REFATOM_SIZE, refcoords, defcoords, defatomcoords)
                            residue.createAtom("H2", newcoords, "ATOM")
                            residue.addDebumpAtom(residue.getAtom("H2"))
                            count += 1
                        
                        if nterm == 3 and "H3" not in residue.get("map"):
                            defatomcoords = NTERM3_COORDS
                            newcoords = findCoordinates(REFATOM_SIZE, refcoords, defcoords, defatomcoords)
                            residue.createAtom("H3", newcoords, "ATOM")
                            residue.addDebumpAtom(residue.getAtom("H3"))
                            count += 1
                            
                    elif cterm == 2: # Neutral C-terminus
                        refcoords = []
                        defcoords = []
                        refcoords.append(residue.getAtom("C").getCoords())
                        refcoords.append(residue.getAtom("O").getCoords())
                        refcoords.append(residue.getAtom("OXT").getCoords())
                        defcoords.append(defresidue.getAtom("C").getCoords())
                        defcoords.append(defresidue.getAtom("O").getCoords())
                        defcoords.append(OXT_COORDS)
                        defatomcoords = CTERM_COORDS
                        newcoords = findCoordinates(REFATOM_SIZE, refcoords, defcoords, defatomcoords)
                        residue.createAtom("HO", newcoords, "ATOM")
                        residue.addDebumpAtom(residue.getAtom("HO"))
                        count += 1
         
                elif type == 4:
                    name = residue.get("naname")
                    defresidue = self.nadef.getResidue(name)
                    for defatom in defresidue.get("atoms"):
                        refcoords = []
                        defcoords = []
                        if not defatom.isHydrogen(): continue
                        defname = defatom.get("name")
                        atom = residue.getAtom(defname)
                        if atom != None: continue
                        bonds = defresidue.makeBondList(residue,defname)
                        if len(bonds) < REFATOM_SIZE:
                            error = "Not enough bonds to remake hydrogen in %s %i" % \
                                    (name, residue.get("resSeq"))
                            raise ValueError, error
                        for i in range(REFATOM_SIZE):
                            refcoords.append(residue.getAtom(bonds[i]).getCoords())
                            defcoords.append(defresidue.getAtom(bonds[i]).getCoords())
                        defatomcoords = defresidue.getAtom(defname).getCoords()
                        newcoords = findCoordinates(REFATOM_SIZE, refcoords, defcoords, defatomcoords)
                        residue.createAtom(defname, newcoords, "ATOM")
                        residue.addDebumpAtom(residue.getAtom(defname))
                        count += 1
                        
                prevres = residue
        self.write(" Added %i hydrogen atoms.\n" % count)


    def repairAA(self, residues, resnum):
        """
            Repair heavy atoms in Amino Acid (type 1) residues

            Parameters
                residues: The list of residues in the chain (list)
                resnum:   The index of the residue to fix (int)
            Returns
                sgadded:  1 if an CYS SG is added, 0 otherwise
                          Used to detect if updateSSbridges must be
                          called
        """
        sgadded = 0
        residue = residues[resnum]
        seen = {}
        resname = residue.get("name")
        resSeq = residue.get("resSeq")
        missing = residue.get("missing")
        origlen = len(missing)
        defresidue = self.aadef.getResidue(resname)
        while len(missing) > 0:
            bonds = []
            refcoords = []
            defcoords = []
            missing.reverse()
            atomname = missing.pop()
            missing.reverse()

            if atomname == "O" or atomname == "C":
                N = None
                if residue.getAtom("CA") != None:
                    bonds.append("CA")
                if atomname == "O" and residue.getAtom("C") != None:
                    bonds.append("C")
                elif atomname == "C" and residue.getAtom("O") != None:
                    bonds.append("O")

                try:
                    if len(bonds) != 2: raise IndexError
                    if residue.get("isCterm") and residue.getAtom("N") != None:
                        bonds.append("N")
                    else:
                        nextres = residues[resnum + 1]
                        N = nextres.getAtom("N")
                except IndexError:
                    text = "\tUnable to repair %s %i %s\n" % (resname, resSeq, atomname)
                    raise ValueError, text
                for i in range(len(bonds)):
                    refcoords.append(residue.getAtom(bonds[i]).getCoords())
                    defcoords.append(defresidue.getAtom(bonds[i]).getCoords())
                defatomcoords = defresidue.getAtom(atomname).getCoords()
                if N != None:
                    refcoords.append(N.getCoords())
                    defcoords.append(PEP_TRANS_N)
                    bonds.append("N")

            elif atomname == "N" and not residue.get("isNterm"):# and resname == "GLY":
                try:
                    prevres = residues[resnum - 1]
                except IndexError:
                    text = "\tUnable to repair %s %i\n" % (resname, resSeq)
                    raise ValueError, text
                if prevres.getAtom("C") != None:
                    bonds.append("C")
                    refcoords.append(prevres.getAtom("C").getCoords())
                    defcoords.append(defresidue.getAtom("C").getCoords())
                if prevres.getAtom("CA") != None:
                    bonds.append("CA")
                    refcoords.append(prevres.getAtom("CA").getCoords())
                    defcoords.append(defresidue.getAtom("CA").getCoords())  
                elif prevres.getAtom("O") != None:
                    bonds.append("O")
                    refcoords.append(prevres.getAtom("O").getCoords())
                    defcoords.append(defresidue.getAtom("O").getCoords()) 

                if residue.getAtom("CA") != None:
                    bonds.append("CA")
                    refcoords.append(residue.getAtom("CA").getCoords())
                    defcoords.append(PEP_TRANS_CA)
                        
                defatomcoords = PEP_TRANS_N
                        
            else:
                bonds = defresidue.makeBondList(residue, atomname)
                if atomname == "OXT":
                    defatomcoords = OXT_COORDS
                else:
                    defatomcoords = defresidue.getAtom(atomname).getCoords()
                for i in range(len(bonds)):
                    refcoords.append(residue.getAtom(bonds[i]).getCoords())
                    defcoords.append(defresidue.getAtom(bonds[i]).getCoords())
                            
            # Now refcoords has the reference atom coordinates from the PDB
            #     defcoords has the reference atom frame coordinates from AA.DAT
            #     defatomcoords has the reference frame new atom coordinates
                    
            if len(bonds) < REFATOM_SIZE:
                if atomname not in seen and origlen > 1:
                    seen[atomname] = 1
                    missing.append(atomname)
                elif seen[atomname] < origlen - 1:
                    seen[atomname] += 1
                    missing.append(atomname)
                else:
                    text = "\tUnable to repair %s %i\n" % (resname, resSeq)
                    raise ValueError, text
            else:
                newcoords = findCoordinates(REFATOM_SIZE, refcoords, defcoords, defatomcoords)
                residue.createAtom(atomname, newcoords,"ATOM")
                self.write("Added %s to %s %i at coordinates " % (atomname, resname, resSeq),1)
                self.write("%.3f %.3f %.3f\n" % (newcoords[0], newcoords[1], newcoords[2]))
                residue.addDebumpAtom(residue.getAtom(atomname))
                if atomname == "SG" and resname == "CYS":
                    sgadded = 1

        return sgadded

    def repairNA(self, residue):
        """
            Repair heavy atoms in Nucleic Acid (type 4) residues

            Parameters
                residue:  The residue to repair (residue)
        """
        seen = {}
        resname = residue.get("name")
        resSeq = residue.get("resSeq")
        missing = residue.get("missing")
        origlen = len(missing)
        defresidue = self.nadef.getResidue(residue.get("naname"))
        while len(missing) > 0:
            bonds = []
            refcoords = []
            defcoords = []
            missing.reverse()
            atomname = missing.pop()
            missing.reverse()

            bonds = defresidue.makeBondList(residue, atomname)
            defatomcoords = defresidue.getAtom(atomname).getCoords()
            for i in range(len(bonds)):
                refcoords.append(residue.getAtom(bonds[i]).getCoords())
                defcoords.append(defresidue.getAtom(bonds[i]).getCoords())
                            
            # Now refcoords has the reference atom coordinates from the PDB file
            #     defcoords has the reference atom frame coordinates from NA.DAT
            #     defatomcoords has the reference frame new atom coordinates
                    
            if len(bonds) < REFATOM_SIZE:
                if atomname not in seen and origlen > 1:
                    seen[atomname] = 1
                    missing.append(atomname)
                elif seen[atomname] < origlen - 1:
                    seen[atomname] += 1
                    missing.append(atomname)
                else:
                    text = "\tUnable to repair %s %i\n" % (resname, resSeq)
                    raise ValueError, text
            else:
                newcoords = findCoordinates(REFATOM_SIZE, refcoords, defcoords, defatomcoords)
                residue.createAtom(atomname, newcoords,"ATOM")
                self.write("Added %s to %s %i at coordinates " % (atomname, resname, resSeq),1)
                self.write("%.3f %.3f %.3f\n" % (newcoords[0], newcoords[1], newcoords[2]))
                residue.addDebumpAtom(residue.getAtom(atomname))

    def repairHeavy(self):
        """
            Repair all heavy atoms in residues with missing atoms

            Returns
                sgadded:  1 if an CYS SG is added, 0 otherwise
                          Used to detect if updateSSbridges must be
                          called.
        """
        self.write("Attempting to repair heavy atoms...\n")
        sgadded = 0
        for chain in self.protein.getChains():
            residues = chain.get("residues")
            for resnum in range(len(residues)):
                residue = residues[resnum]
                type = residue.get("type")
                if type == 1:
                    sgadded = self.repairAA(residues, resnum)
                elif type == 4:
                    self.repairNA(residue)
            
        self.write("Done.\n")

        return sgadded

    def createCells(self):
        """
            Place each atom in a virtual cell for easy neighbor comparison
        """       
        self.cells = {}
        for chain in self.protein.getChains():
            for residue in chain.get("residues"):
                for atom in residue.get("atoms"):
                    self.addCell(atom)

    def addCell(self, atom):
        """
            Add an atom to the cell

            Parameters
                atom:  The atom to add (atom)
        """
        size = CELL_SIZE
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
            self.cells[key].append(atom)
        except KeyError:
            self.cells[key] = [atom]
        atom.set("cell", key)

    def removeCell(self, atom):
        """
             Remove the atom from a cell

             Parameters
                 atom:   The atom to add (atom)
        """
        oldcell = atom.get("cell")
        atom.set("cell", None)
        try:
            for i in range(len(self.cells[oldcell]) - 1):
                if self.cells[oldcell][i] == atom:
                    self.cells[oldcell].pop(i)
                    return
        except KeyError: # Shouldn't occur, but it is okay
            pass

    def getNearCells(self, atom):
        """
            Find all atoms in bordering cells to an atom

            Parameters
                atom:  The atom to use (atom)
            Returns
                closeatoms:  A list of nearby atoms (list)
        """
        size = CELL_SIZE
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
                        newatoms = self.cells[newkey]
                        for atom2 in newatoms:
                            if atom == atom2: continue
                            closeatoms.append(atom2)
                    except KeyError: pass
                        
        return closeatoms
    
    def debumpProtein(self):
        """
             Ensure that an added atom is not on top of another atom.  If it
             is, debump the residue.  A threshold is used to determine
             all nearby atoms so the entire protein need not be searched for
             every new dihedral angle.
        """
        self.write("Checking if we must debump any residues...\n")
        self.createCells()
        cells = self.cells  
        bumpresidues = []
        debumpmap = {}
        debumpAtoms = []
        for chain in self.protein.getChains():
            for residue in chain.get("residues"):
                for atom in residue.get("debumpAtoms"):
                    if atom != None: debumpAtoms.append(atom)
       
        for atom1 in debumpAtoms:
            atomname = atom1.get("name")
            if atomname == "H": continue
            residue1 = atom1.get("residue")
            if residue1.get("isNterm") or residue1.get("isCterm"): continue

            # NOTE: For now, disable NA debumping
            if residue1.get("type") == 4: continue
            
            closeatoms = []
            coords1 = atom1.getCoords()

            closeatoms = self.getNearCells(atom1)
                
            for atom2 in closeatoms:
                if atom2.get("residue") == residue1: continue
                elif residue1.get("SSbondpartner") == atom2: continue
                coords2 = atom2.getCoords()
                residue2 = atom2.get("residue")
                dist = distance(coords1, coords2)
                compdist = BUMP_DIST
                if atom1.isHydrogen(): compdist = BUMP_HDIST
                if dist < compdist:
                    if residue1 not in bumpresidues:
                        self.write("Must debump %s %i due to %s in %s %i\n" % \
                                   (residue1.name, residue1.resSeq, atom2.name, \
                                    residue2.name, residue2.resSeq), 1)
                        debumpmap[residue1] = [atomname]
                        bumpresidues.append(residue1)
                    else:
                        debumpmap[residue1].append(atomname)
                    break
            
        self.write("Done.\n")

        for i in range(len(bumpresidues)):
            residue = bumpresidues[i]
            causenames = debumpmap[residue]
            type = residue.get("type")
            if type == 1:
                defresidue = self.aadef.getResidue(residue.get("name"))
            elif type == 4:
                defresidue = self.nadef.getResidue(residue.get("naname"))
            self.write("Starting to debump %s %i...\n" % \
                       (residue.get("name"), residue.get("resSeq")))
            value = self.debumpResidue(residue, causenames, defresidue)
            if value == 1:
                self.write("Debumping Successful.\n")
            else:
                text = "WARNING: Unable to debump %s %i\n" % \
                       (residue.get("name"), residue.get("resSeq"))
                self.write("********\n")
                self.write(text)
                self.write("********\n")
                self.warnings.append(text)
  
 
    def debumpResidue(self, residue, causenames, defresidue):
        """
             Ensure that an added atom was not added on top of another
             atom by rotating certain atoms about a dihedral angle.  Finds
             a working angle with no conflicts from other residues and sets
             the new atom coordinates.

             Parameters
                 residue:    The residue to be debumped (Residue)
                 causenames: A list of atom names that must be moved (list)
                 defresidue: The definition for this residue (DefinitionResidue)
        """

        # Get the first chi number

        oldchinum = -1
        chinum = self.pickChiangle(residue, defresidue, causenames, oldchinum)
        if chinum == -1:
            return 0
        step = 0
        bestangle = LARGE_TORSION_ANGLE
        bestvalue = LARGE_TORSION_ANGLE
        value = LARGE_TORSION_ANGLE
        newcauses = []
        
        while bestvalue > 0.0 and step < 10:
            if oldchinum != -1:
                causenames = newcauses
                chinum = self.pickChiangle(residue, defresidue, newcauses, oldchinum)
            if chinum == -1:
                return 0   
            self.write("Using dihedral angle number %i to debump the residue.\n" % chinum , 1)
            for angle in range(-180,180,5):
                value, newcauses = self.setChiangle(residue, chinum, angle, defresidue, causenames)
                if value < bestvalue:
                    bestangle = angle
                    bestvalue = value
                    if value == 0.0:
                        self.write("No conflicts found at angle %.2f.\n" % angle, 1)
                        return 1
                        
            if value > 0.0:    
                value, newcauses = self.setChiangle(residue, chinum, bestangle, defresidue, causenames)
                step += 1

            oldchinum = chinum

        return 0

    def pickChiangle(self, residue, defresidue, causenames, oldnum):
        """
            Choose a chiangle number to use in debumping

            Algorithm
                Instead of simply picking a random chiangle, this function
                uses a more intelligent method to improve efficiency.
                The algorithm uses the names of the conflicting atoms
                within the residue to determine which chiangle number
                has the best chance of fixing the problem(s).  If more than
                one chiangle number resolves the same number of atoms,
                the new number is picked randomly.  The method also
                insures that the same chiangle will not be run twice
                in a row.
            Parameters
                residue:    The residue that is being debumped (Residue)
                defresidue: The definition of the residue (DefinitionResidue)
                causenames: A list of atom names that are currently
                            conflicts (list)
                oldnum    : The old chiangle number (int)
            Returns
                chinum    : The new chiangle number (int)
        """
        scores = []
        best = 0
        for i in range(len(residue.get("chiangles"))):
            if i == oldnum: continue
            if residue.get("chiangles")[i] == LARGE_TORSION_ANGLE: continue
            score = 0
            rootname = defresidue.get("dihedralatoms")[i*4 + 2]
            movenames = self.getMovenames(residue, defresidue, rootname)
            for causename in causenames:
                if causename in movenames:
                    score += 1
                    if score > best:
                        best = score
                        scores = [i]
                    elif score == best:
                        scores.append(i)

        # If no chinums move the problem atom, we can't debump

        if len(scores) == 0:
            return -1
        
        chinum = scores[random.randint(0, len(scores) - 1)]

        return chinum
    
    def setChiangle(self, residue, chinum, chiangle, defresidue, causenames=[]):
        """
             Set the chiangle and move the appropriate atoms

             Parameters
                 residue :   The residue that is being debumped (Residue)
                 chinum  :   The chi number (int)
                 chiangle:   The new angle to set (float)
                 defresidue: The definition of the residue (DefinitionResidue)
                 causenames: The atom names that must be moved (list)
             Returns
                 value     : A value indicating whether any atoms are too
                             close or not - any non-zero value means that
                             there is still a conflict.
        """
        BASIS = 1
        POINT = 2

        torsatoms = []
        rootname = ""
        
        oldchi = residue.get("chiangles")[chinum]
        if oldchi > 180.0:
            raise ValueError,"Invalid dihedral angle size %.3f!" % oldchi

        difchi = chiangle - oldchi
        
        for i in range(4):
            defatomname = defresidue.get("dihedralatoms")[chinum*4 + i]
            atom = residue.getAtom(defatomname)
            if atom != None:
                torsatoms.append(atom.getCoords())
                if i == POINT:
                    rootname = defatomname
            else:
                raise ValueError, "Error occurred while trying to debump!"

        initcoords = subtract(torsatoms[POINT], torsatoms[BASIS])

        movenames = self.getMovenames(residue, defresidue, rootname)
        movecoords = []
        for name in movenames:
            atom = residue.getAtom(name)
            movecoords.append(subtract(atom.getCoords(), torsatoms[BASIS]))
            
        newcoords = qchichange(initcoords, movecoords, difchi)

        for i in range(len(movenames)):
            name = movenames[i]
            atom = residue.getAtom(name)
            self.removeCell(atom)
            x = (newcoords[i][0] + torsatoms[BASIS][0])
            y = (newcoords[i][1] + torsatoms[BASIS][1])
            z = (newcoords[i][2] + torsatoms[BASIS][2])
            atom.set("x", x)
            atom.set("y", y)
            atom.set("z", z)
            self.addCell(atom)
            
            #print "%s is moved to coords %5.3f %5.3f %5.3f" % (name, atom.x, atom.y, atom.z)

        torsatoms = []
        for i in range(4):
            defatomname = defresidue.get("dihedralatoms")[chinum*4 + i]
            atom = residue.getAtom(defatomname)
            if atom != None:
                torsatoms.append(atom.getCoords())
            else:
                raise ValueError, "Error occurred while trying to set chiangle!"
        
        di = getDihedral(torsatoms[0], torsatoms[1], torsatoms[2], torsatoms[3])
        #print "New Chiangle: %5.3f" % di
        residue.get("chiangles")[chinum] = di

        if causenames == []:
            return

        newcauses = []
        value = 0.0
        for atom1 in residue.get("atoms"):
            atomname = atom1.get("name")
            if not residue.get("isNterm") and (atomname in movenames or atomname in causenames):
                coords1 = atom1.getCoords()
                closeatoms = self.getNearCells(atom1)
                for atom2 in closeatoms:
                    if atom2.get("residue") == residue: continue
                    elif residue.get("SSbondpartner") == atom2: continue
                    coords2 = atom2.getCoords()
                    dist = distance(coords1, coords2)
                    compdist = BUMP_DIST
                    if atom1.isHydrogen(): compdist = BUMP_HDIST
                    if dist < compdist:
                        if atomname not in newcauses:
                            newcauses.append(atomname)
                        value += (BUMP_DIST - dist)                    
    
        return value, newcauses

    def getMovenames(self, residue, defresidue, rootname):
        """
            Determine the names of the atoms that are to be moved based
            on the chiangle's root atom.

            Parameters
                residue:    The residue that is being debumped (Residue)
                defresidue: The definition of the residue (DefinitionResidue)
                rootname:   The name of the root atom used to calculate
                            distance.
            Returns
                movenames:  A list of atom names to move (list)
        """
        movenames = []
        rootatom = defresidue.getAtom(rootname)
        initdist = rootatom.get("refdistance")
        for atom in residue.get("atoms"):
            atomname = atom.get("name")
            defatom = defresidue.getAtom(atomname)
            if defatom == None and (residue.get("isCterm") or residue.get("isNterm")):
                continue
           
            if defatom.get("refdistance") > initdist:
                if residue.get("name") == "ILE" and rootname == "CG1":
                    if atomname in ["CD1","HD11","HD12","HD13","HG12","HG13"]:
                        movenames.append(atomname)
                elif residue.get("name") == "THR" and rootname == "OG1":
                    if atomname == "HG1":
                        movenames.append(atomname)
                else:
                    movenames.append(atomname)
        return movenames

    def setDonorsAndAcceptors(self):
        """
            Set the donors and acceptors within the protein
        """
        self.updateIntraBonds()
        for atom in self.protein.getAtoms():
            atomname = atom.get("name")
            resname = atom.residue.get("name")
            if atomname.startswith("N"):
                for bond in atom.get("intrabonds"):
                    if bond[0] == "H":
                        atom.set("hdonor",1)
                        break
                if not ((atomname == "NZ" and resname == "LYS") or \
                       atom.residue.get("isNterm") or atomname == "N"):
                    atom.set("hacceptor",1)
            elif atomname.startswith("O") or atomname.startswith("S"):
                atom.set("hacceptor",1)
                for bond in atom.get("intrabonds"):
                    if bond[0] == "H":
                        atom.set("hdonor",1)
                        break
                    
    def printHbond(self):
        """
            Print a list of all hydrogen bonds to stdout.  A hydrogen bond
            is defined when a donor has a hydrogen within 3.3 Angstroms, and
            the Hyd-Donor-Accepor angle is less than 20 degrees.
        """
        self.write("Printing hydrogen bond list...\n")
        from hydrogens import hydrogenRoutines
        hydRoutines = hydrogenRoutines(self)
        dlist = []
        alist = []
        self.setDonorsAndAcceptors()
        for atom in self.protein.getAtoms():
            if atom.get("hdonor"): dlist.append(atom)
            if atom.get("hacceptor"):
                alist.append(atom)
        
        for donor in dlist:
            donorhs = []
            for bond in donor.get("intrabonds"):
                if bond[0] == "H":
                    donorhs.append(bond)
            for acc in alist:
                if acc == donor: continue
                for donorh in donorhs:
                    donorhatom = donor.get("residue").getAtom(donorh)
                    dist = distance(donorhatom.getCoords(), acc.getCoords())
                    if dist > 3.3: continue
                    angle = hydRoutines.getHbondangle(acc, donor, donorhatom)
                    if angle <= 20.0:
                        print "Donor: %s %s %i  \tAcceptor: %s %s %i\tHdist: %.2f\tAngle: %.2f" % \
                              (donor.resName, donor.name, donor.residue.resSeq, acc.resName, \
                               acc.name, acc.residue.resSeq, dist, angle)
                        
    def optimizeHydrogens(self):
        """
            Wrapper function for hydrogen optimizing routines.  The routines
            were too extensive to properly fit within this file.
        """
        self.write("Beginning to optimize hydrogens...\n")
        from hydrogens import hydrogenRoutines
        self.updateIntraBonds()
        self.calculateChiangles()
        myhydRoutines = hydrogenRoutines(self)
        myhydRoutines.readHydrogenDefinition()
        myhydRoutines.optimizeHydrogens()

    def optimizeWaters(self):
        """
            Wrapper function for water optimizing routines.
        """
        run = 0
        for atom in self.protein.getAtoms():
            res = atom.get("residue")
            if res.get("type") == 3:
                if res.getAtom("H1") == None or \
                   res.getAtom("H2") == None:
                    run = 1
                    break
        if run == 0: return
        from hydrogens import hydrogenRoutines
        self.write("Optimizing water hydrogens.\n")
        mywatRoutines = hydrogenRoutines(self)
        mywatRoutines.readHydrogenDefinition()
        mywatRoutines.optimizeWaters()
        self.write("Done optimizing hydrogens.\n")

    def randomizeWaters(self):
        """
            Instead of optimizing, find each WAT O and place H1 and H2
            while giving it a random orientation
        """
        run = 0
        for atom in self.protein.getAtoms():
            res = atom.get("residue")
            if res.get("type") == 3:
                if res.getAtom("H1") == None or \
                   res.getAtom("H2") == None:
                    run = 1
                    break
        if run == 0: return
        from hydrogens import hydrogenRoutines
        self.write("Randomizing water hydrogens.\n")
        myrandRoutines = hydrogenRoutines(self)
        myrandRoutines.readHydrogenDefinition()
        myrandRoutines.randomizeWaters()
        self.write("Done randomizing hydrogens.\n")
