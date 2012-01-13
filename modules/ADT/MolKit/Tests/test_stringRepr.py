#
#
#
#
# $Id: test_stringRepr.py,v 1.2 2005/06/08 18:06:51 rhuey Exp $
#

import unittest
from string import split
from MolKit.protein import ProteinSet, ChainSet, ResidueSet
from MolKit.protein import Protein, Chain, Residue
from MolKit.molecule import MoleculeSet, AtomSet
from MolKit.molecule import Atom
from MolKit.tree import TreeNodeSet
from MolKit.stringSelector import StringSelector
from MolKit import Read



class StringReprBaseTests(unittest.TestCase):

    def setUp(self):
        self.mols = Read('Data/stringSel.pdb')
        self.mol = self.mols[0]
        self.stringSel = StringSelector()
    

    def tearDown(self):
        """
        clean-up
        """
        for mol in self.mols:
            del mol


    #tests for correction stringRepr initialization below sets
    def test_correct_MoleculeSet_stringRepr(self):
        """
         check moleculeSet stringRepr is correct
        """
        self.assertEquals(self.mols.stringRepr, "stringSel")


    def test_correct_ChainSet_stringRepr(self):
        """
         check ChainSet stringRepr is correct
        """
        self.assertEquals(self.mols.chains.stringRepr, "stringSel:")


    def test_correct_ResidueSet_stringRepr(self):
        """
         check ResidueSet stringRepr is correct
        """
        self.assertEquals(self.mols.chains.residues.stringRepr, "stringSel::")


    def test_correct_AtomSet_stringRepr(self):
        """
         check AtomSet stringRepr is correct
        """
        self.assertEquals(self.mols.chains.residues.atoms.stringRepr, "stringSel:::")


    #tests for correction stringRepr initialization below individual
    def test_correct_Molecule_ChainSet_stringRepr(self):
        """
         check single Molecule's chainSet stringRepr is correct
        """
        self.assertEquals(self.mols[0].chains.stringRepr, "stringSel:")


    def test_correct_Chain1_ResidueSet_stringRepr(self):
        """
         check Chain1's residues stringRepr is correct
        """
        self.assertEquals(self.mols.chains[0].residues.stringRepr, "stringSel:A:")


    def test_correct_Chain2_ResidueSet_stringRepr(self):
        """
         check two chain's residues stringRepr is correct
        """
        self.assertEquals(self.mols.chains[1].residues.stringRepr, "stringSel:W:")


    def Xtest_correct_2Chain_ResidueSet_stringRepr(self):
        """
         check two Chains residues stringRepr is correct
        """
        self.assertEquals(self.mols.chains[0:].residues.stringRepr, "stringSel::")


    def test_correct_last_Chain_ResidueSet_stringRepr(self):
        """
         check last Chain's residues stringRepr is correct
        """
        self.assertEquals(self.mols.chains[-1].residues.stringRepr, "stringSel:W:")


    def test_correct_Mol0_Chain0_ResidueSet_stringRepr(self):
        """
         check single Mol's single Chain's residues stringRepr 
        """
        self.assertEquals(self.mols[0].chains[0].residues.stringRepr, "stringSel:A:")


    def Xtest_correct_Mol0_Chains1_2_ResidueSet_stringRepr(self):
        """
         check single Mol's two Chain's residues stringRepr 
        """
        self.assertEquals(self.mols[0].chains[0:].residues.stringRepr, "stringSel::")


    def test_correct_Mols_Chains1_residues_atoms_stringRepr(self):
        """
         check single Chain's residues's atoms stringRepr 
        """
        self.assertEquals(self.mols.chains[0].residues.atoms.stringRepr, "stringSel:A::")


    def test_select_with_empty_string(self):
        """
         test result with empty string returns all mols
        """
        result, msg = self.stringSel.select(self.mols, "")
        self.assertEquals(result.stringRepr, "stringSel")



    def test_setLevel_protein_protein(self):
        """
         test setting protein level to protein
        """
        mols = self.mols.setLevel(Protein)
        self.assertEqual(mols.stringRepr, "stringSel")


    def test_setLevel_protein_chain(self):
        """
         test setting protein level to chain
        """
        chains = self.mols.setLevel(Chain)
        self.assertEqual(chains.stringRepr, "stringSel:")


    def test_setLevel_protein_residue(self):
        """
         test setting protein level to residue
        """
        res = self.mols.setLevel(Residue)
        self.assertEqual(res.stringRepr, "stringSel::")


    def test_setLevel_protein_residue(self):
        """
         test setting protein level to Atom
        """
        ats = self.mols.setLevel(Atom)
        self.assertEqual(ats.stringRepr, "stringSel:::")




    def test_setLevel_chain_protein(self):
        """
         test setting chain level to protein
        """
        mols = self.mols.chains.setLevel(Protein)
        self.assertEqual(mols.stringRepr, "stringSel")


    def test_setLevel_chain_chain(self):
        """
         test setting chain level to chain
        """
        chains = self.mols.chains.setLevel(Chain)
        self.assertEqual(chains.stringRepr, "stringSel:")


    def test_setLevel_chain_residue(self):
        """
         test setting chain level to residue
        """
        res = self.mols.chains.setLevel(Residue)
        self.assertEqual(res.stringRepr, "stringSel::")


    def test_setLevel_chain_atom(self):
        """
         test setting chain level to Atom
        """
        ats = self.mols.chains.setLevel(Atom)
        self.assertEqual(ats.stringRepr, "stringSel:::")



    def test_setLevel_residues_protein(self):
        """
         test setting residue level to protein
        """
        mols = self.mols.chains.residues.setLevel(Protein)
        self.assertEqual(mols.stringRepr, "stringSel")


    def test_setLevel_residues_chain(self):
        """
         test setting residues level to chain
        """
        chains = self.mols.chains.residues.setLevel(Chain)
        self.assertEqual(chains.stringRepr, "stringSel:")


    def test_setLevel_residues_residue(self):
        """
         test setting residues level to residue
        """
        res = self.mols.chains.residues.setLevel(Residue)
        self.assertEqual(res.stringRepr, "stringSel::")


    def test_setLevel_residues_atom(self):
        """
         test setting residue level to Atom
        """
        ats = self.mols.chains.residues.setLevel(Atom)
        self.assertEqual(ats.stringRepr, "stringSel:::")


class TwoMoleculeStringReprTests(StringReprBaseTests):

    def setUp(self):
        self.mols1 = Read('Data/stringSel.pdb')
        self.mol1 = self.mols1[0]
        self.mols2 = Read('Data/protease.pdb')
        self.mol2 = self.mols2[0]
        self.mols = self.mols1 + self.mols2
        self.stringSel = StringSelector()
    

    #tests for correction stringRepr initialization below sets
    def test_correct_MoleculeSet_stringRepr(self):
        """
         check moleculeSet stringRepr is correct
        """
        self.assertEquals(self.mols.stringRepr, "stringSel/+/protease")


    def test_correct_ChainSet_stringRepr(self):
        """
         check ChainSet stringRepr is correct
        """
        self.assertEquals(self.mols.chains.stringRepr, "stringSel:/+/protease:")


    def test_correct_ResidueSet_stringRepr(self):
        """
         check ResidueSet stringRepr is correct
        """
        self.assertEquals(self.mols.chains.residues.stringRepr, "stringSel::/+/protease::")


    def test_correct_AtomSet_stringRepr(self):
        """
         check AtomSet stringRepr is correct
        """
        self.assertEquals(self.mols.chains.residues.atoms.stringRepr, "stringSel:::/+/protease:::")


    #tests for correction stringRepr initialization below individual
    def test_correct_Molecule_ChainSet_stringRepr(self):
        """
         check single Molecule's chainSet stringRepr is correct
        """
        self.assertEquals(self.mol1.chains.stringRepr, "stringSel:")


    def test_correct_Chain1_ResidueSet_stringRepr(self):
        """
         check Chain1's residues stringRepr is correct
        """
        self.assertEquals(self.mols.chains[0].residues.stringRepr, "stringSel:A:")


    def test_correct_Chain2_ResidueSet_stringRepr(self):
        """
         check two chain's residues stringRepr is correct
        """
        self.assertEquals(self.mols.chains[1].residues.stringRepr, "stringSel:W:")


    def Qtest_correct_2Chain_ResidueSet_stringRepr(self):
        """
         check two Chains residues stringRepr is correct
        """
        self.assertEquals(self.mols.chains[0-2].residues.stringRepr, "stringSel::")


    def test_correct_last_Chain_ResidueSet_stringRepr(self):
        """
         check last Chain's residues stringRepr is correct
        """
        self.assertEquals(self.mols.chains[-1].residues.stringRepr, "protease:B:")


    def test_correct_Mol0_Chain0_ResidueSet_stringRepr(self):
        """
         check single Mol's single Chain's residues stringRepr 
        """
        self.assertEquals(self.mols[0].chains[0].residues.stringRepr, "stringSel:A:")


    def Xtest_correct_Mol0_Chains1_2_ResidueSet_stringRepr(self):
        """
         check single Mol's two Chain's residues stringRepr 
        """
        self.assertEquals(self.mols[0].chains[0:].residues.stringRepr, "stringSel::")


    def test_correct_Mols_Chains1_residues_atoms_stringRepr(self):
        """
         check single Chain's residues's atoms stringRepr 
        """
        self.assertEquals(self.mols.chains[0].residues.atoms.stringRepr, "stringSel:A::")


    def test_select_with_empty_string(self):
        """
         test result with empty string returns all mols
        """
        result, msg = self.stringSel.select(self.mols, "")
        self.assertEqual(self.mols.stringRepr, "stringSel/+/protease")



    def test_setLevel_protein_protein(self):
        """
         test setting protein level to protein
        """
        mols = self.mols.setLevel(Protein)
        self.assertEqual(mols.stringRepr, "stringSel/+/protease")


    def test_setLevel_protein_chain(self):
        """
         test setting protein level to chain
        """
        chains = self.mols.setLevel(Chain)
        self.assertEqual(chains.stringRepr, "stringSel:/+/protease:")


    def test_setLevel_protein_residue(self):
        """
         test setting protein level to residue
        """
        res = self.mols.setLevel(Residue)
        self.assertEqual(res.stringRepr, "stringSel::/+/protease::")


    def test_setLevel_protein_residue(self):
        """
         test setting protein level to Atom
        """
        ats = self.mols.setLevel(Atom)
        self.assertEqual(ats.stringRepr, "stringSel:::/+/protease:::")


    def test_setLevel_chain_protein(self):
        """
         test setting chain level to protein
        """
        mols = self.mols.chains.setLevel(Protein)
        self.assertEqual(mols.stringRepr, "stringSel/+/protease")


    def test_setLevel_chain_chain(self):
        """
         test setting chain level to chain
        """
        chains = self.mols.chains.setLevel(Chain)
        self.assertEqual(chains.stringRepr, "stringSel:/+/protease:")


    def test_setLevel_chain_residue(self):
        """
         test setting chain level to residue
        """
        res = self.mols.chains.setLevel(Residue)
        self.assertEqual(res.stringRepr, "stringSel::/+/protease::")


    def test_setLevel_chain_atom(self):
        """
         test setting chain level to Atom
        """
        ats = self.mols.chains.setLevel(Atom)
        self.assertEqual(ats.stringRepr, "stringSel:::/+/protease:::")



    def test_setLevel_residues_protein(self):
        """
         test setting residue level to protein
        """
        mols = self.mols.chains.residues.setLevel(Protein)
        self.assertEqual(mols.stringRepr, "stringSel/+/protease")


    def test_setLevel_residues_chain(self):
        """
         test setting residues level to chain
        """
        chains = self.mols.chains.residues.setLevel(Chain)
        self.assertEqual(chains.stringRepr, "stringSel:/+/protease:")


    def test_setLevel_residues_residue(self):
        """
         test setting residues level to residue
        """
        res = self.mols.chains.residues.setLevel(Residue)
        self.assertEqual(res.stringRepr, "stringSel::/+/protease::")


    def test_setLevel_residues_atom(self):
        """
         test setting residue level to Atom
        """
        ats = self.mols.chains.residues.setLevel(Atom)
        self.assertEqual(ats.stringRepr, "stringSel:::/+/protease:::")





if __name__ == '__main__':
    unittest.main()
