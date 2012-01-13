#
#
#
#
# $Id: test_stringSelector.py,v 1.19 2008/11/24 21:03:09 rhuey Exp $
#

import unittest
from string import split
from MolKit.stringSelector import StringSelector, CompoundStringSelector
from MolKit.molecule import Atom, AtomSet
from MolKit.protein import Protein, Chain, Residue, ProteinSet, ChainSet, ResidueSet
from MolKit.sets import Sets


molecules = None
molecules2 = None

class StringSelectorBaseTests(unittest.TestCase):

    def initMolecules(self):
        global molecules
        if molecules is None:
            from MolKit import Read
            molecules = Read('Data/stringSel.pdb')
        self.molecules = molecules

    def setUp(self):
        if not hasattr(self, 'molecules'):
            self.initMolecules()


    def tearDown(self):
        """
        clean-up
        """
        for m in self.molecules:
            del(m)
        del(self.molecules)


    def test_constructor(self):
        """
        instantiate a StringSelector
        """
        stringSel = StringSelector()
        self.assertEquals(stringSel.__class__, StringSelector)


    def test_select_with_empty_string(self):
        """
         test result with empty string returns all molecules
        """
        stringSel = StringSelector()
        result, msg = stringSel.select(self.molecules, "")
        self.assertEquals(result, self.molecules)
        self.assertEquals(result.stringRepr, self.molecules[0].name)

    def test_select_with_lambda_atoms_expr(self):
        """
         test result with lambda x:x.name==CA returns 3 atoms
        """
        selString = ":::lambda x:x.name=='CA'"
        #sanity check
        stringSel = StringSelector()
        result, msg = stringSel.select(self.molecules, selString)
        self.assertEquals(len(result), 3)
        self.assertEquals(result.stringRepr, selString)


    def test_select_with_lambda_residue_expr(self):
        """
         test result with lambda x:x.name==CA returns 9 atoms
        """
        selString = "::lambda x:x.type=='GLN':"
        #sanity check
        stringSel = StringSelector()
        result, msg = stringSel.select(self.molecules, selString)
        self.assertEquals(len(result), 9)
        self.assertEquals(result.stringRepr, selString)


    def test_select_with_lambda_chain_expr(self):
        """
         test result with lambda x:x.id=='W' returns 5 atoms
        """
        selString = ":lambda x:x.id=='W'::"
        #sanity check
        stringSel = StringSelector()
        result, msg = stringSel.select(self.molecules, selString)
        self.assertEquals(len(result), 5)
        self.assertEquals(result.stringRepr, selString)


    def test_select_with_lambda_molecule_expr(self):
        """
         test result with lambda x:x.name=='stringSel' returns 29 atoms
        """
        selString = "lambda x:x.name=='stringSel':::"
        #sanity check
        stringSel = StringSelector()
        result, msg = stringSel.select(self.molecules, selString)
        self.assertEquals(len(result), 29)
        self.assertEquals(result.stringRepr, selString)

    def test_select_with_CA_returns_3atoms(self):
        """
         test result with CA returns 3 atoms
        """
        selString = ":::CA"
        #sanity check
        stringSel = StringSelector()
        result, msg = stringSel.select(self.molecules, selString)
        self.assertEquals(len(result), 3)
        self.assertEquals(result.stringRepr, selString)


    def test_select_with_backbone_returns_12atoms(self):
        """
         test result with backbone returns 12 atoms
        """
        selString = ":::backbone"
        #sanity check
        stringSel = StringSelector()
        result, msg = stringSel.select(self.molecules, selString)
        self.assertEquals(len(result), 12)
        self.assertEquals(result.stringRepr, selString)



    def test_select_with_hetatm_returns_5atoms(self):
        """
         test result with hetatm returns 5 atoms
        """
        selString = ":::hetatm"
        #sanity check
        stringSel = StringSelector()
        result, msg = stringSel.select(self.molecules, selString)
        self.assertEquals(len(result), 5)
        self.assertEquals(result.stringRepr, selString)


    def test_select_with_backbone_returns_no_waters(self):
        """
         test result with backbone returns no waters
        """
        selString = ":::backbone"
        stringSel = StringSelector()
        result, msg = stringSel.select(self.molecules, selString)
        self.assertEquals(len(result), 12)
        self.assertEqual('HOH' not in result.parent.type, True)
        self.assertEquals(result.stringRepr, selString)


    def test_select_with_bad_key_returns_EmptySet(self):
        """
         test result with bad_key returns EmptySet
        """
        selString = ":::'backbone'"
        stringSel = StringSelector()
        result, msg = stringSel.select(self.molecules, selString)
        self.assertEquals(len(result), 0)
        #previously
        #self.assertEquals(result, None)


    def test_select_with_simple_Atom_full_name(self):
        """
         test result with simple Atom full_name
        """
        subset = self.molecules.allAtoms[:2]
        selString = subset.full_name()
        stringSel = StringSelector()
        result, msg = stringSel.select(self.molecules, selString)
        self.assertEqual(result, subset)
        self.assertEquals(result.stringRepr, selString)


    def test_select_with_simple_Residue_full_name(self):
        """
         test result with simple Residue full_name
        """
        subset = self.molecules.chains.residues[:2]
        selString = subset.full_name()
        stringSel = StringSelector()
        result, msg = stringSel.select(self.molecules, selString)
        self.assertEqual(result, subset)
        self.assertEquals(result.stringRepr, selString)


    def test_select_with_compound_Atom_full_name(self):
        """
         test result with compound Atom full_name
        """
        subset = self.molecules.allAtoms[:2]
        selString = subset[0].full_name() + ';' + subset[1].full_name()
        stringSel = StringSelector()
        result, msg = stringSel.select(self.molecules, selString)
        self.assertEqual(result, subset)
        self.assertEquals(result.stringRepr, selString)


    def test_select_with_5_Atom_full_names(self):
        """
         test result with >2 Atom full_names joined with ';'
        """
        subset = self.molecules.allAtoms[:5]
        selString = ''
        for a in subset:
            selString += a.full_name() + ';'
        selString = selString[:-1]
        stringSel = StringSelector()
        result, msg = stringSel.select(self.molecules, selString)
        self.assertEqual(result, subset)
        self.assertEquals(result.stringRepr, selString)


    def test_select_with_compound_Residue_full_name(self):
        """
         test result with compound Residue full_name
        """
        subset = self.molecules.chains.residues[:2]
        selString = subset[0].full_name() + ';' + subset[1].full_name()
        stringSel = StringSelector()
        result, msg = stringSel.select(self.molecules, selString)
        self.assertEqual(result, subset)
        self.assertEquals(result.stringRepr, selString)


    def test_select_with_5_Residue_full_names(self):
        """
         test result with 5 Residue full_names joined with ';'
        """
        subset = self.molecules.chains.residues[:5]
        selString = ''
        for a in subset:
            selString += a.full_name() + ';'
        selString = selString[:-1]
        stringSel = StringSelector()
        result, msg = stringSel.select(self.molecules, selString)
        self.assertEqual(result, subset)
        self.assertEquals(result.stringRepr, selString)


    def test_select_all_atoms_from_compound_string(self):
        """
         test result using all the atoms' full_name()
        """
        selString = ''
        subset = self.molecules.allAtoms
        for a in subset:
            selString += a.full_name() + ';'
        selString = selString[:-1]
        stringSel = StringSelector()
        result, msg = stringSel.select(self.molecules, selString)
        self.assertEqual(result, subset)
        self.assertEquals(result.stringRepr, selString)
        

    def test_select_set_from_string(self):
        """
         test result using string which is key in Sets dict
        """
        sets = Sets()
        subset = self.molecules[0].allAtoms
        set_name = 'first'
        sets.add(set_name, subset)
        selString = ':::' + set_name
        stringSel = StringSelector()
        result, msg = stringSel.select(self.molecules, selString, sets=sets)
        self.assertEqual(result, subset)
        self.assertEquals(result.stringRepr, selString)
        



class CompoundStringSelectorBaseTests(unittest.TestCase):

    def initMolecules(self):
        global molecules, molecules2
        from MolKit import Read
        if molecules is None:
            molecules = Read('Data/stringSel.pdb')
        if molecules2 is None:
            molecules2 = Read('Data/protease.pdb')
        self.mols1 = molecules
        self.mols2 = molecules2
        self.mols = molecules + molecules2

    def setUp(self):
        if not hasattr(self, 'molecules'):
            self.initMolecules()


#    def setUp(self):
#        from MolKit import Read
#        if not hasattr(self, 'mols'):
#            self.mols1 = Read('Data/stringSel.pdb')
#            self.mols2 = Read('Data/protease.pdb')
#            self.mols = self.mols1 + self.mols2


    def test_constructor(self):
        """
        instantiate a StringSelector
        """
        stringSel = CompoundStringSelector()
        self.assertEquals(stringSel.__class__, CompoundStringSelector)


class AddStringSelectorTests(CompoundStringSelectorBaseTests):
    #add
    def test_select_add_mols(self):
        """
        test selecting with stringRepr of mols of added mols
        """
        stringSel = CompoundStringSelector()
        selString = self.mols.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, self.mols)
        self.assertEquals(selString, "stringSel/+/protease")


    def test_select_add_mols_chains(self):
        """
        test selecting with stringRepr of chains of added mols
        """
        stringSel = CompoundStringSelector()
        selString = self.mols.chains.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, self.mols.chains)
        self.assertEquals(selString, "stringSel:/+/protease:")


    def test_select_add_mols_residues(self):
        """
        test selecting with stringRepr of residues of added mols
        """
        stringSel = CompoundStringSelector()
        selString = self.mols.chains.residues.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, self.mols.chains.residues)
        self.assertEquals(selString, "stringSel::/+/protease::")



    def test_select_add_mols_atoms(self):
        """
        test selecting with stringRepr of atoms of added mols
        """
        stringSel = CompoundStringSelector()
        selString = self.mols.chains.residues.atoms.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, self.mols.chains.residues.atoms)
        #print "selString=", selString
        self.assertEquals(selString, "stringSel:::/+/protease:::")


    def test_select_add_mols_allAtoms(self):
        """
        test selecting with stringRepr of atoms of added mols
        """
        stringSel = CompoundStringSelector()
        selString = self.mols.allAtoms.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, self.mols.allAtoms)
        self.assertEquals(selString, "stringSel:::/+/protease:::")



class UnionStringSelectorTests(CompoundStringSelectorBaseTests):
    #union: '|'
    def test_select_union_mols(self):
        """
        test selecting with stringRepr of union mols
        """
        stringSel = CompoundStringSelector()
        diff_mols = self.mols | self.mols1
        selString = diff_mols.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_mols)
        self.assertEquals(selected, self.mols)
        self.assertEquals(selString, "stringSel/+/protease/|/stringSel")


    def test_select_union_chains(self):
        """
        test selecting with stringRepr of union chains
        """
        stringSel = CompoundStringSelector()
        diff_chains = self.mols.chains | self.mols1.chains
        selString = diff_chains.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_chains)
        self.assertEquals(selected, self.mols.chains)
        self.assertEquals(selString, "stringSel:/+/protease:/|/stringSel:")


    def test_select_union_residues(self):
        """
        test selecting with stringRepr of union residues
        """
        stringSel = CompoundStringSelector()
        diff_residues = self.mols.chains.residues | self.mols1.chains.residues
        selString = diff_residues.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_residues)
        self.assertEquals(selected, self.mols.chains.residues)
        self.assertEquals(selString, "stringSel::/+/protease::/|/stringSel::")


    def test_select_union_atoms(self):
        """
        test selecting with stringRepr of union atoms
        """
        stringSel = CompoundStringSelector()
        diff_atoms = self.mols.chains.residues.atoms | self.mols1.chains.residues.atoms
        selString = diff_atoms.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_atoms)
        self.assertEquals(selected, self.mols.chains.residues.atoms)
        self.assertEquals(selString, "stringSel:::/+/protease:::/|/stringSel:::")


    def test_select_union_allAtoms(self):
        """
        test selecting with stringRepr of union allAtoms
        """
        stringSel = CompoundStringSelector()
        diff_atoms = self.mols1.allAtoms | self.mols2.allAtoms
        selString = diff_atoms.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_atoms)
        self.assertEquals(selected, self.mols.allAtoms)
        self.assertEquals(selString, "stringSel:::/|/protease:::")


    #union: '|' returning single copy
    def test_select_union_mols_single(self):
        """
        test selecting with stringRepr of union 2X same mols returning single
        """
        stringSel = CompoundStringSelector()
        diff_mols = self.mols1 | self.mols1
        selString = diff_mols.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_mols)
        self.assertEquals(selected, self.mols1)
        self.assertEquals(selString, "stringSel")


    def test_select_union_chains_single(self):
        """
        test selecting with stringRepr of union 2X same chains returning
        single
        """
        stringSel = CompoundStringSelector()
        diff_chains = self.mols1.chains | self.mols1.chains
        selString = diff_chains.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_chains)
        self.assertEquals(selected, self.mols1.chains)
        self.assertEquals(selString, "stringSel:")


    def test_select_union_residues_all(self):
        """
        test selecting with stringRepr of union 2X same residues returning
        single 
        """
        stringSel = CompoundStringSelector()
        diff_residues = self.mols1.chains.residues | self.mols1.chains.residues
        selString = diff_residues.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_residues)
        self.assertEquals(selected, self.mols1.chains.residues)
        self.assertEquals(selString, "stringSel::")


    def test_select_union_atoms_all(self):
        """
        test selecting with stringRepr of union 2X atoms returning single copy
        """
        stringSel = CompoundStringSelector()
        diff_atoms = self.mols1.chains.residues.atoms | self.mols1.chains.residues.atoms
        selString = diff_atoms.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_atoms)
        self.assertEquals(selected, self.mols1.chains.residues.atoms)
        self.assertEquals(selString, "stringSel:::")


    def test_select_union_allAtoms_all(self):
        """
        test selecting with stringRepr of union 2X allAtoms returning single
        copy
        """
        stringSel = CompoundStringSelector()
        diff_atoms = self.mols1.allAtoms | self.mols1.allAtoms
        selString = diff_atoms.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_atoms)
        self.assertEquals(selected, self.mols1.allAtoms)
        self.assertEquals(selString, "stringSel:::")



class SubtractStringSelectorTests(CompoundStringSelectorBaseTests):
    #subtract
    def test_select_subtract_mols(self):
        """
        test selecting with stringRepr of subtracted mols
        """
        stringSel = CompoundStringSelector()
        diff_mols = self.mols - self.mols1
        selString = diff_mols.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_mols)
        self.assertEquals(selString, "stringSel/+/protease/-/stringSel")


    def test_select_subtract_chains(self):
        """
        test selecting with stringRepr of subtracted chains
        """
        stringSel = CompoundStringSelector()
        diff_chains = self.mols.chains - self.mols1.chains
        selString = diff_chains.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_chains)
        self.assertEquals(selString, "stringSel:/+/protease:/-/stringSel:")


    def test_select_subtract_residues(self):
        """
        test selecting with stringRepr of subtracted residues
        """
        stringSel = CompoundStringSelector()
        diff_residues = self.mols.chains.residues - self.mols1.chains.residues
        selString = diff_residues.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_residues)
        self.assertEquals(selString, "stringSel::/+/protease::/-/stringSel::")


    def test_select_subtract_atoms(self):
        """
        test selecting with stringRepr of subtracted atoms
        """
        stringSel = CompoundStringSelector()
        diff_atoms = self.mols.chains.residues.atoms - self.mols1.chains.residues.atoms
        selString = diff_atoms.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_atoms)
        self.assertEquals(selString, "stringSel:::/+/protease:::/-/stringSel:::")


    def test_select_subtract_allAtoms(self):
        """
        test selecting with stringRepr of subtracted allAtoms
        """
        stringSel = CompoundStringSelector()
        diff_atoms = self.mols.allAtoms - self.mols1.allAtoms
        selString = diff_atoms.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_atoms)
        self.assertEquals(selString, "stringSel:::/+/protease:::/-/stringSel:::")


class IntersectionStringSelectorTests(CompoundStringSelectorBaseTests):
    #intersection: '&'
    def test_select_intersect_mols(self):
        """
        test selecting with stringRepr of intersect mols
        """
        stringSel = CompoundStringSelector()
        diff_mols = self.mols & self.mols1
        selString = diff_mols.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_mols)
        self.assertEquals(selString, "stringSel/+/protease/&/stringSel")
        self.assertEquals(selected, self.mols1)


    def test_select_intersect_chains(self):
        """
        test selecting with stringRepr of intersected chains
        """
        stringSel = CompoundStringSelector()
        diff_chains = self.mols.chains & self.mols1.chains
        selString = diff_chains.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_chains)
        self.assertEquals(selString, "stringSel:/+/protease:/&/stringSel:")
        self.assertEquals(selected, self.mols1.chains)


    def test_select_intersect_residues(self):
        """
        test selecting with stringRepr of intersected residues
        """
        stringSel = CompoundStringSelector()
        diff_residues = self.mols.chains.residues & self.mols1.chains.residues
        selString = diff_residues.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_residues)
        self.assertEquals(selected, self.mols1.chains.residues)
        self.assertEquals(selString, "stringSel::/+/protease::/&/stringSel::")


    def test_select_intersect_atoms(self):
        """
        test selecting with stringRepr of intersected atoms
        """
        stringSel = CompoundStringSelector()
        diff_atoms = self.mols.chains.residues.atoms & self.mols1.chains.residues.atoms
        selString = diff_atoms.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_atoms)
        self.assertEquals(selected, self.mols1.chains.residues.atoms)
        self.assertEquals(selString, "stringSel:::/+/protease:::/&/stringSel:::")


    def test_select_intersect_allAtoms(self):
        """
        test selecting with stringRepr of intersected allAtoms
        """
        stringSel = CompoundStringSelector()
        diff_atoms = self.mols.allAtoms & self.mols1.allAtoms
        selString = diff_atoms.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_atoms)
        self.assertEquals(selected, self.mols1.allAtoms)
        self.assertEquals(selString, "stringSel:::/+/protease:::/&/stringSel:::")


    #intersection: '&' returning []
    def test_select_intersect_mols_empty(self):
        """
        test selecting with stringRepr of empty intersect mols
        """
        stringSel = CompoundStringSelector()
        diff_mols = self.mols1 & self.mols2
        selString = diff_mols.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_mols)
        self.assertEquals(selString, "stringSel/&/protease")
        self.assertEquals(selected, ProteinSet())


    def test_select_intersect_chains_empty(self):
        """
        test selecting with stringRepr of empty intersected chains
        """
        stringSel = CompoundStringSelector()
        diff_chains = self.mols1.chains & self.mols2.chains
        selString = diff_chains.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_chains)
        self.assertEquals(selected, ChainSet())
        self.assertEquals(selString, "stringSel:/&/protease:")


    def test_select_intersect_residues_empty(self):
        """
        test selecting with stringRepr of empty intersected residues
        """
        stringSel = CompoundStringSelector()
        diff_residues = self.mols1.chains.residues & self.mols2.chains.residues
        selString = diff_residues.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_residues)
        self.assertEquals(selected, ResidueSet())
        self.assertEquals(selString, "stringSel::/&/protease::")


    def test_select_intersect_atoms_empty(self):
        """
        test selecting with stringRepr of empty intersected atoms
        """
        stringSel = CompoundStringSelector()
        diff_atoms = self.mols1.chains.residues.atoms & self.mols2.chains.residues.atoms
        selString = diff_atoms.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_atoms)
        self.assertEquals(selected, AtomSet())
        self.assertEquals(selString, "stringSel:::/&/protease:::")


    def test_select_intersect_allAtoms(self):
        """
        test selecting with stringRepr of intersected allAtoms
        """
        stringSel = CompoundStringSelector()
        diff_atoms = self.mols1.allAtoms & self.mols2.allAtoms
        selString = diff_atoms.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_atoms)
        self.assertEquals(selected, AtomSet())
        self.assertEquals(selString, "stringSel:::/&/protease:::")



class XorStringSelectorTests(CompoundStringSelectorBaseTests):
    #xor: '^'
    def test_select_xor_mols(self):
        """
        test selecting with stringRepr of xor mols
        """
        stringSel = CompoundStringSelector()
        diff_mols = self.mols ^ self.mols1
        selString = diff_mols.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_mols)
        self.assertEquals(selected, self.mols2)
        self.assertEquals(selString, "stringSel/+/protease/^/stringSel")


    def test_select_xor_chains(self):
        """
        test selecting with stringRepr of xor chains
        """
        stringSel = CompoundStringSelector()
        diff_chains = self.mols.chains ^ self.mols1.chains
        selString = diff_chains.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_chains)
        self.assertEquals(selected, self.mols2.chains)
        self.assertEquals(selString, "stringSel:/+/protease:/^/stringSel:")


    def test_select_xor_residues(self):
        """
        test selecting with stringRepr of xor residues
        """
        stringSel = CompoundStringSelector()
        diff_residues = self.mols.chains.residues ^ self.mols1.chains.residues
        selString = diff_residues.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_residues)
        self.assertEquals(selected, self.mols2.chains.residues)
        self.assertEquals(selString, "stringSel::/+/protease::/^/stringSel::")


    def test_select_xor_atoms(self):
        """
        test selecting with stringRepr of xor atoms
        """
        stringSel = CompoundStringSelector()
        diff_atoms = self.mols.chains.residues.atoms ^ self.mols1.chains.residues.atoms
        selString = diff_atoms.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_atoms)
        self.assertEquals(selected, self.mols2.chains.residues.atoms)
        self.assertEquals(selString, "stringSel:::/+/protease:::/^/stringSel:::")


    def test_select_xor_allAtoms(self):
        """
        test selecting with stringRepr of xor allAtoms
        """
        stringSel = CompoundStringSelector()
        diff_atoms = self.mols.allAtoms ^ self.mols1.allAtoms
        selString = diff_atoms.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_atoms)
        self.assertEquals(selected, self.mols2.allAtoms)
        self.assertEquals(selString, "stringSel:::/+/protease:::/^/stringSel:::")


    #xor: '^' returning all
    def test_select_xor_mols_all(self):
        """
        test selecting with stringRepr of xor mols returning all
        """
        stringSel = CompoundStringSelector()
        diff_mols = self.mols1 ^ self.mols2
        selString = diff_mols.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_mols)
        self.assertEquals(selected, self.mols)
        self.assertEquals(selString, "stringSel/^/protease")


    def test_select_xor_chains_all(self):
        """
        test selecting with stringRepr of xor chains returning all
        """
        stringSel = CompoundStringSelector()
        diff_chains = self.mols1.chains ^ self.mols2.chains
        selString = diff_chains.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_chains)
        self.assertEquals(selected, self.mols.chains)
        self.assertEquals(selString, "stringSel:/^/protease:")


    def test_select_xor_residues_all(self):
        """
        test selecting with stringRepr of xor residues returning all
        """
        stringSel = CompoundStringSelector()
        diff_residues = self.mols1.chains.residues ^ self.mols2.chains.residues
        selString = diff_residues.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_residues)
        self.assertEquals(selected, self.mols.chains.residues)
        self.assertEquals(selString, "stringSel::/^/protease::")


    def test_select_xor_atoms_all(self):
        """
        test selecting with stringRepr of xor atoms returning all
        """
        stringSel = CompoundStringSelector()
        diff_atoms = self.mols1.chains.residues.atoms ^ self.mols2.chains.residues.atoms
        selString = diff_atoms.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_atoms)
        self.assertEquals(selected, self.mols.chains.residues.atoms)
        self.assertEquals(selString, "stringSel:::/^/protease:::")


    def test_select_xor_allAtoms_all(self):
        """
        test selecting with stringRepr of xor allAtoms returning all
        """
        stringSel = CompoundStringSelector()
        diff_atoms = self.mols1.allAtoms ^ self.mols2.allAtoms
        selString = diff_atoms.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, diff_atoms)
        self.assertEquals(selected, self.mols.allAtoms)
        self.assertEquals(selString, "stringSel:::/^/protease:::")


class SubsetStringSelectorTests(CompoundStringSelectorBaseTests):

    #subset: 's'
    def test_select_subset_mols(self):
        """
        test selecting with stringRepr for subset
        """
        stringSel = CompoundStringSelector()
        diff_mols = self.mols[-1:]
        selString = "(stringSel/+/protease\\s\\protease)"
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        #selected, msg = stringSel.select(self.mols, selString)
        self.assertEquals(selected, diff_mols)
        self.assertEquals(selected, self.mols[-1:])
        self.assertEquals(selected.stringRepr, selString)


    def test_select_subset_chains(self):
        """
        test selecting with stringRepr of subset of chains
        """
        stringSel = CompoundStringSelector()
        diff_chains = self.mols.chains.get('B')
        selString = "(stringSel:/+/protease:\\s\\B)"
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        #selected, msg = stringSel.select(self.mols, selString)
        self.assertEquals(selected, diff_chains)
        self.assertEquals(selected, self.mols.chains[-1:])
        self.assertEquals(selected.stringRepr, selString)


    def test_select_subset_residues(self):
        """
        test selecting with stringRepr of subset residues
        """
        stringSel = CompoundStringSelector()
        #this is the name of the 7th residue in the first molecule
        diff_residues = self.mols.chains.residues.get("HOH617")
        selString = "(stringSel::/+/protease::\\s\\HOH617)"
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        #selected, msg = stringSel.select(self.mols, selString)
        self.assertEquals(selected, diff_residues)
        self.assertEquals(selected, self.mols.chains.residues[7:8])
        self.assertEquals(selected.stringRepr, selString)


    def test_select_subset_atoms(self):
        """
        test selecting with stringRepr of subset atoms
        select the first atom in protease from mols.chains.residues.atoms
        """
        stringSel = CompoundStringSelector()
        diff_atoms = self.mols.chains.residues.atoms.get('0')
        selString = "(stringSel:::/+/protease:::\\s\\0)"
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        #selected, msg = stringSel.select(self.mols, selString)
        self.assertEquals(selected, diff_atoms)
        self.assertEquals(selected, self.mols.chains.residues.atoms[0:1])
        self.assertEquals(selected.stringRepr,  selString)


    def test_select_subset_allAtoms(self):
        """
        test selecting with stringRepr of subset allAtoms
        select the first atom in protease from mols.allAtoms
        """
        stringSel = CompoundStringSelector()
        diff_atoms = self.mols.allAtoms.get('0')
        selString = "(stringSel:::/+/protease:::\\s\\0)"
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        #selected, msg = stringSel.select(self.mols, selString)
        self.assertEquals(selected, diff_atoms)
        self.assertEquals(selected, self.mols.allAtoms[0:1])
        self.assertEquals(selected.stringRepr,  selString)


class MixedSelectorTests(CompoundStringSelectorBaseTests):

    def test_select_complicated_sidechain_atoms(self):
        """
        test selecting w/stringRepr of complicated sidechain atoms(bug 672)
        ARG8 atoms excluding backbone atoms and CB atoms

        """
        stringSel = CompoundStringSelector()
        #arg is ARG8 in chain A
        arg = self.mols2.chains.residues[7]
        sidechain = arg.atoms-arg.atoms.get('backbone')-arg.atoms.get('CB')
        #diff_atoms = self.mols1.allAtoms | self.mols1.allAtoms
        selString = sidechain.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, sidechain)
        self.assertEquals(selected, self.mols2.allAtoms[67:73])
        self.assertEquals(selString, 'protease:A:ARG8:/-/(protease:A:ARG8:\\s\\backbone)/-/(protease:A:ARG8:\\s\\CB)')
        self.assertEquals(selString, selected.stringRepr)


    def test_select_complicated_sidechain_atoms2(self):
        """
        test selecting w/stringRepr of complicated sidechain atoms(bug 673)

        """
        stringSel = CompoundStringSelector()
        #arg is ARG8 in chain A
        arg = self.mols2.chains.residues[7]
        sidechain = arg.atoms.get('sidechain')-arg.atoms.get('CB')
        #diff_atoms = self.mols1.allAtoms | self.mols1.allAtoms
        selString = sidechain.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        self.assertEquals(selected, sidechain)
        self.assertEquals(selected, self.mols2.allAtoms[67:73])
        self.assertEquals(selString, selected.stringRepr)


    def test_select_get_from_get_result(self):
        """
        test selecting w/stringRepr of get from the result of a get

        """
        stringSel = CompoundStringSelector()
        m = self.mols2
        cats = m.allAtoms.get('C.?')
        cbats = cats.get('CB')
        expectedString = "((protease:::\\s\\C.?)\\s\\CB)"
        selString = cbats.stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        #selected, msg = stringSel.select(self.mols, selString)
        self.assertEquals(selected, cbats)
        self.assertEquals(selString, expectedString)
        self.assertEquals(selected.stringRepr, expectedString)


    def test_select_set_from_string(self):
        """
         test result using key into Sets
        """
        sets = Sets()
        subset = self.mols.allAtoms[:3]
        set_name = 'first_three_atoms'
        sets.add(set_name, subset)
        selString = ':::' + set_name
        css = CompoundStringSelector()
        result, msg = css.select(self.mols, selString, sets=sets)
        self.assertEqual(result, subset)
        self.assertEquals(result.stringRepr, selString)
        
    
    def test_select_complicated_sidechain_atoms3(self):
        """
        test selecting w/stringRepr of complicated sidechain atoms
        -> carbon atoms in sidechains of ARG residues in protease except CB atoms
        """
        stringSel = CompoundStringSelector()
        #arg is ARG8 in chain A
        args = self.mols2.chains.residues.get('ARG*')
        sidechains = args.atoms.get('sidechain')
        sc_carbons = sidechains.get('C*')
        #sc_carbons = sidechains.get('C.?')
        target = sc_carbons - args.atoms.get('CB')
        selString = target.stringRepr
        #check that built set has expected stringRepr
        targetStr = '(((protease::\\s\\ARG*):\\s\\sidechain)\\s\\C*)/-/((protease::\\s\\ARG*):\\s\\CB)'
        self.assertEquals(selString, targetStr)
        # try to get same atoms by selecting using built set's stringRepr
        selected, msg = stringSel.select(self.mols, selString, returnMsg=True)
        #sanity check that 24 atoms are found
        self.assertEquals(len(selected), 24)
        #check that the same atoms are found
        self.assertEquals(selected, target)
        #check that selected set has expected stringRepr
        self.assertEquals(targetStr, selected.stringRepr)


    def test_select_with_lambda_atoms_expr_with_parentheses(self):
        """
         test result with lambda x:len(x.bonds)==1 returns 468 atoms for hsg1.pdb
        """
        selString = ":::lambda x:len(x.bonds)==1"
        stringSel = StringSelector()
        from MolKit import Read
        mol = Read("Data/hsg1.pdb")
        mol[0].buildBondsByDistance()
        result, msg = stringSel.select(mol, selString)
        self.assertEquals(len(result), 468)
        self.assertEquals(result.stringRepr, selString)



if __name__ == '__main__':
    test_cases = [
        'StringSelectorBaseTests',
        'CompoundStringSelectorBaseTests',
        'AddStringSelectorTests',
        'UnionStringSelectorTests',
        'SubtractStringSelectorTests',
        'IntersectionStringSelectorTests',
        'XorStringSelectorTests',
        'SubsetStringSelectorTests',
        'MixedSelectorTests',
    ]
    unittest.main( argv=([__name__ ,] + test_cases) )
    #unittest.main()
