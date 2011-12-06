"""
    Routines for PDB2PQR

    This module contains the protein object used in PDB2PQR and associated
    methods
    
    Todd Dolinsky (todd@ccb.wustl.edu)
    Washington University in St. Louis

"""

__date__ = "14 November 2003"
__author__ = "Todd Dolinsky"

from pdb import *
from structures import *

class Protein:
    """
        Protein class

        The protein class represents the parsed PDB, and provides a
        hierarchy of information - each Protein contains a list of Chain
        objects as provided in the PDB file.  Each Chain then contains its
        associated list of Residue objects, and each Residue contains a list
        of Atom objects, completing the hierarchy.
    """

    def __init__(self, pdblist):
        """
            Initialize using parsed PDB file

            Parameters
                pdblist: List of Classes of PDB lines as created
                         by pdb.py->readPDB
        """

        self.chainmap, self.chains = self.createProtein(pdblist)
        #for chain in self.chains:
            #chain.renumberResidues()

    def createProtein(self, pdblist):
        """
            Fill the Protein with chains, residues, and atoms

            Parameters
                pdblist: List of Classes of PDB lines as created
                         by pdb.py->readPDB (list)
            Returns
                dict:    Mapping of chain ID to chain object
                list:    List of chain objects sorted by chain ID (dict)
        """

        dict = {}
        list = []
        
        previousAtom = None
        residue = []
        numModels = 0

        for record in pdblist:
            if isinstance(record, ATOM) or isinstance(record, HETATM):
                chainID = record.chainID
                resSeq = record.resSeq
                resName = record.resName
                iCode = record.iCode

                if previousAtom == None:
                    previousAtom = record
                
                if chainID not in dict:
                    myChain = Chain(chainID)
                    dict[chainID] = myChain
                        
                if resSeq != previousAtom.resSeq or \
                       iCode != previousAtom.iCode:
                    myResidue = Residue(residue, previousAtom)
                    dict[previousAtom.chainID].addResidue(myResidue)
                    residue = []

                residue.append(record)
                previousAtom = record

            elif isinstance(record, END):
                myResidue = Residue(residue, previousAtom)
                dict[previousAtom.chainID].addResidue(myResidue)
                residue = []

            elif isinstance(record, MODEL):
                numModels += 1
                if residue == []: continue
                if numModels > 1:
                    myResidue = Residue(residue, previousAtom)
                    dict[previousAtom.chainID].addResidue(myResidue)
                    break

        if residue != [] and numModels <= 1:
            myResidue = Residue(residue, previousAtom)
            dict[previousAtom.chainID].addResidue(myResidue)

        chainmap = dict.copy()
        if dict.has_key(""):
            dict["ZZ"] = dict[""]
            del dict[""]

        keys = dict.keys()
        keys.sort()

        for key in keys:
            list.append(dict[key])
        
        return chainmap, list

    def printAtoms(self, atomlist):
        """
            Get the text for the entire protein
            Parameters
                atomlist: The list of atoms to include (list)
            Returns
                text:     The list of (stringed) atoms (list)
        """
        self.reSerialize()
        text = []
        for atom in atomlist:
            text.append("%s\n" % str(atom))
        return text

    def reSerialize(self):
        """
            Generate new serial numbers for atoms in the protein
        """
        count = 1
        for atom in self.getAtoms():
            atom.set("serial", count)
            count += 1
    
    def numResidues(self):
        """
            Get the number of residues for the entire protein (including
            multiple chains)

            Returns
                count:  Number of residues in the protein (int)
        """
        count = 0
        for chain in self.chains:
            count += chain.numResidues()
        return count

    def numAtoms(self):
        """
            Get the number of atoms for the entire protein(including
            multiple chains)

            Returns
                count:  Number of atoms in the protein (int)
        """
        count = len(self.getAtoms())
        return count

    def getAtoms(self):
        """
            Return all Atom objects in list format

            Returns
                atomlist:  List of Atom objects in the protein (list)
        """

        atomlist = []
        for chain in self.chains:
            for atom in chain.getAtoms():
                atomlist.append(atom)
        return atomlist

    def getCharge(self):
        """
            Get the total charge on the protein
            NOTE:  Since the misslist is used to identify incorrect
                   charge assignments, this routine does not list the
                   3 and 5 termini of nucleic acid chains as having
                   non-integer charge even though they are (correctly)
                   non-integer.
            Returns:
                misslist: List of residues with non-integer
                          charges (list)
                charge:   The total charge on the protein (float)
        """
        charge = 0.0
        misslist = []
        for chain in self.chains:
            for residue in chain.get("residues"):
                rescharge = residue.getCharge()
                charge = charge + rescharge
                if residue.get("is3term") or residue.get("is5term"):
                    continue
                if float("%i" % rescharge) != rescharge:
                    misslist.append(residue)
        return misslist, charge

    def getChains(self):
        """
            Get the chains object

            Returns
                chains: The list of chains in the protein (chain)
        """
        return self.chains
