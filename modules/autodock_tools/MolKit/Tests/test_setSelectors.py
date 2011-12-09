#
#
#
#
# $Id: test_setSelectors.py,v 1.6 2006/03/29 17:45:38 rhuey Exp $
#

import unittest
from string import split
from MolKit.protein import ProteinSetSelector, ChainSetSelector, ChainSet
from MolKit.protein import ResidueSetSelector, ResidueSet
from MolKit.molecule import AtomSetSelector, AtomSet, BondSetSelector, BondSet
from MolKit.tree import TreeNodeSet
from MolKit.stringSelector import StringSelector
from MolKit import Read
from MolKit.sets import Sets



class ProteinSetSelectorBaseTests(unittest.TestCase):

    def test_constructor(self):
        """
        instantiate an ProteinSetSelector
        """
        stringSel = ProteinSetSelector()
        self.assertEquals(stringSel.__class__, ProteinSetSelector)


#    def test_constructorOptions(self):
#        """
#         test possible constructor options
#            options = [userPref='cS']
#        """
#        stringSel = ProteinSetSelector(self.mols, self.selString, userPref='cI')
#        self.assertEquals(stringSel.__class__, ProteinSetSelector)



class ProteinSetSelectorTests(ProteinSetSelectorBaseTests):


    def setUp(self):
        self.mols = Read('Data/stringSel.pdb')
        self.mol = self.mols[0]
        self.stringSel = ProteinSetSelector()
    

    def tearDown(self):
        """
        clean-up
        """
        del(self.mol)


    def test_select_with_empty_string(self):
        """
         test result with empty string returns all mols
        """
        result, msg = self.stringSel.select(self.mols, "")
        self.assertEquals(result, self.mols)


    def test_select_end(self):
        """
         test select with '$'  returns last item
        """
        new_mols = Read("Data/1crn.pdb")
        self.mols +=new_mols
        result, msg = self.stringSel.select(self.mols, "$")
        self.assertEquals(result[-1], self.mols[-1])


    def test_select_with_valid_index(self):
        """
         test string with valid_index returns set with 1 item 
        """
        selString = "0"
        #selString = "1"
        result, msg = self.stringSel.select(self.mols, selString)
        self.assertEquals(len(result), 1)


    def test_select_with_invalid_index_returns_empty_set(self):
        """
         test string with invalid_index returns empty set
        """
        selString = "2"
        result, msg = self.stringSel.select(self.mols, selString)
        self.assertEqual(result.__class__, self.mols.__class__)
        self.assertEqual(len(result), 0)
        self.assertEqual(msg[0], selString)


    def test_select_with_valid_range(self):
        """
         test string with valid_range returns set with 2 items
        """
        new_mols = Read("Data/1crn.pdb")
        self.mols +=new_mols
        selString = "0-1"
        #selString = "1-2"
        result, msg = self.stringSel.select(self.mols, selString)
        self.assertEquals(len(result), 2)


    def test_select_with_invalid_range(self):
        """
         test string with invalid_range returns set with 0 items 
        """
        new_mols = Read("Data/1crn.pdb")
        self.mols +=new_mols
        selString = "4-6"
        result, msg = self.stringSel.select(self.mols, selString)
        self.assertEquals(len(result), 0)


    def test_select_with_valid_regex(self):
        """
         test string with valid_regex returns set with 1 items
         <this regex is intended to match 1, 2,or 3 followed by anything>
        """
        new_mols = Read("Data/1crn.pdb")
        self.mols +=new_mols
        selString = "[1-3]*"
        result, msg = self.stringSel.select(self.mols, selString)
        #print "result=", result
        self.assertEquals(len(result), 1)


    def test_select_with_valid_regex_2(self):
        """
         test string with valid_regex returns set with 1 items
         <this regex is intended to match anything in range s-z
         followed by anything>
        """
        new_mols = Read("Data/1crn.pdb")
        self.mols +=new_mols
        selString = "[s-z]*"
        result, msg = self.stringSel.select(self.mols, selString)
        #print "result=", result
        self.assertEquals(len(result), 1)


    def test_select_valid_set(self):
        """
         test string with valid set name returns 1 set
        """
        new_mols = Read("Data/1crn.pdb")
        self.mols +=new_mols
        these_sets = Sets()
        key = 'first'
        these_sets.add(key, self.mols[1:])
        selString = key
        result, msg = self.stringSel.select(self.mols, selString, 
                                        sets=these_sets)
        #print "result=", result
        self.assertEquals(len(result), 1)
        self.assertEquals(result, self.mols[1:])



class ChainSetSelectorBaseTests(unittest.TestCase):

    def test_constructor(self):
        """
        instantiate an ChainSetSelector
        """
        stringSel = ChainSetSelector()
        self.assertEquals(stringSel.__class__, ChainSetSelector)



class ChainSetSelectorTests(ChainSetSelectorBaseTests):


    def setUp(self):
        self.mols = Read('Data/stringSel.pdb')
        new_mols = Read("Data/1crn.pdb")
        self.mols +=new_mols
        self.chains = self.mols.chains
        self.stringSel = ChainSetSelector()
    

    def tearDown(self):
        """
        clean-up
        """
        del(self.chains)
        del(self.mols)


    def test_select_with_empty_string(self):
        """
         test result with empty string returns all chains
        """
        result, msg = self.stringSel.select(self.chains, "")
        self.assertEquals(result, self.chains)


    def test_select_end(self):
        """
        test select with '$'  returns last item
        """
        result, msg = self.stringSel.select(self.chains, "$")
        self.assertEquals(result, self.chains[-1:])


    def test_select_with_valid_index(self):
        """
         test string with valid_index returns set with 1 item 
        """
        selString = "1"
        result, msg = self.stringSel.select(self.chains, selString)
        self.assertEquals(len(result), 1)


    def test_select_with_invalid_index_returns_empty_set(self):
        """
         test string with invalid_index returns empty set
        """
        selString = "7"
        result, msg = self.stringSel.select(self.chains, selString)
        self.assertEqual(result.__class__, self.chains.__class__)
        self.assertEqual(len(result), 0)
        self.assertEqual(msg[0], selString)


    def test_select_with_valid_range(self):
        """
         test string with valid_range returns set with 2 items
        """
        selString = "1-2"
        result, msg = self.stringSel.select(self.chains, selString)
        self.assertEquals(len(result), 2)


    def test_select_with_invalid_range(self):
        """
         test string with invalid_range returns set with 0 items 
        """
        selString = "4-6"
        result, msg = self.stringSel.select(self.chains, selString)
        self.assertEquals(len(result), 0)


    def test_select_with_valid_regex(self):
        """
         test string with valid_regex returns set with 1 items
         <this regex is intended to match anything in range A-Z
         followed by anything>
        """
        selString = "A"
        result, msg = self.stringSel.select(self.chains, selString)
        #print "result=", result
        self.assertEquals(len(result), 1)


    def test_select_valid_set(self):
        """
         test string with valid set name returns 1 set
        """
        these_sets = Sets()
        key = 'first'
        these_sets.add(key, self.mols.chains[1:])
        selString = key
        result, msg = self.stringSel.select(self.chains, selString, 
                                        sets=these_sets)
        #print "result=", result
        self.assertEquals(len(result), 2)
        self.assertEquals(result, self.chains[1:])




class ResidueSetSelectorBaseTests(unittest.TestCase):

    def test_constructor(self):
        """
        instantiate an ResidueSetSelector
        """
        stringSel = ResidueSetSelector()
        self.assertEquals(stringSel.__class__, ResidueSetSelector)



class ResidueSetSelectorTests(ResidueSetSelectorBaseTests):


    def setUp(self):
        self.mols = Read('Data/stringSel.pdb')
        new_mols = Read("Data/1crn.pdb")
        self.mols +=new_mols
        self.residues = self.mols.chains.residues
        self.stringSel = ResidueSetSelector()
    

    def tearDown(self):
        """
        clean-up
        """
        del(self.residues)
        del(self.mols)


    def test_select_with_empty_string(self):
        """
         test result with empty string returns all residues
        """
        result, msg = self.stringSel.select(self.residues, "")
        self.assertEquals(result, self.residues)


    def test_select_end(self):
        """
        test select with '$'  returns last item
        """
        result, msg = self.stringSel.select(self.residues, "$")
        self.assertEquals(result, self.residues[-1:])


    def test_select_with_valid_index(self):
        """
         test string with valid_index returns set with 1 item 
        """
        selString = "1"
        result, msg = self.stringSel.select(self.residues, selString)
        self.assertEquals(len(result), 1)


    def test_select_with_invalid_index_returns_empty_set(self):
        """
         test string with invalid_index returns empty set
        """
        selString = "777"
        result, msg = self.stringSel.select(self.residues, selString)
        self.assertEqual(result.__class__, self.residues.__class__)
        self.assertEqual(len(result), 0)
        self.assertEqual(msg[0], selString)


    def test_select_with_valid_range(self):
        """
         test string with valid_range returns set with 2 items
        """
        selString = "1-2"
        result, msg = self.stringSel.select(self.residues, selString)
        self.assertEquals(len(result), 2)


    def test_select_with_invalid_range(self):
        """
         test string with invalid_range returns set with 0 items 
        """
        selString = "64-71"
        result, msg = self.stringSel.select(self.residues, selString)
        self.assertEquals(len(result), 0)


    def test_select_with_valid_regex(self):
        """
         test string with valid_regex returns set with 1 items
         <this regex is intended to match anything in range A-Z
         followed by anything>
        """
        selString = "A*"
        result, msg = self.stringSel.select(self.residues, selString)
        #print "result=", result
        self.assertEquals(len(result), 11)


    def test_select_valid_set(self):
        """
         test string with valid set name returns 1 set
        """
        these_sets = Sets()
        new_set = self.residues[:10]
        key = 'first_residues'
        these_sets.add(key, new_set)
        selString = key
        result, msg = self.stringSel.select(self.residues, selString, 
                                        sets=these_sets)
        self.assertEquals(len(result), len(new_set))
        self.assertEquals(result, new_set)


class AtomSetSelectorBaseTests(unittest.TestCase):

    def test_constructor(self):
        """
        instantiate an AtomSetSelector
        """
        stringSel = AtomSetSelector()
        self.assertEquals(stringSel.__class__, AtomSetSelector)



class AtomSetSelectorTests(AtomSetSelectorBaseTests):


    def setUp(self):
        self.mols = Read('Data/stringSel.pdb')
        new_mols = Read("Data/1crn.pdb")
        self.mols +=new_mols
        self.atoms = self.mols.chains.residues.atoms
        self.stringSel = AtomSetSelector()
    

    def tearDown(self):
        """
        clean-up
        """
        del(self.atoms)
        del(self.mols)


    def test_select_with_empty_string(self):
        """
         test result with empty string returns all atoms
        """
        result, msg = self.stringSel.select(self.atoms, "")
        self.assertEquals(result, self.atoms)


    def test_select_end(self):
        """
        test select with '$'  returns last item
        """
        result, msg = self.stringSel.select(self.atoms, "$")
        self.assertEquals(result, self.atoms[-1:])


    def test_select_with_valid_index(self):
        """
         test string with valid_index returns set with 1 item 
        """
        selString = "1"
        result, msg = self.stringSel.select(self.atoms, selString)
        self.assertEquals(len(result), 1)


    def test_select_with_invalid_index_returns_empty_set(self):
        """
         test string with invalid_index returns empty set
        """
        selString = "777"
        result, msg = self.stringSel.select(self.atoms, selString)
        self.assertEqual(result.__class__, self.atoms.__class__)
        self.assertEqual(len(result), 0)
        self.assertEqual(msg[0], selString)


    def test_select_with_valid_range(self):
        """
         test string with valid_range returns set with 2 items
        """
        selString = "1-2"
        result, msg = self.stringSel.select(self.atoms, selString)
        #FIX THIS
        self.assertEquals(len(result), 2)


    def test_select_with_invalid_range(self):
        """
         test string with invalid_range returns set with 0 items 
        """
        selString = "649-713"
        result, msg = self.stringSel.select(self.atoms, selString)
        self.assertEquals(len(result), 0)


    def test_select_with_valid_regex(self):
        """
         test string with valid_regex returns set with 1 items
         <this regex is intended to match anything in range A-Z
         followed by anything>
        """
        selString = "CA"
        result, msg = self.stringSel.select(self.atoms, selString)
        #print "result=", result
        self.assertEquals(len(result), 49)


    def test_select_valid_set(self):
        """
         test string with valid set name returns 1 set
        """
        these_sets = Sets()
        new_set = self.atoms[:10]
        key = 'first_atoms'
        these_sets.add(key, new_set)
        selString = key
        result, msg = self.stringSel.select(self.atoms, selString, 
                                        sets=these_sets)
        self.assertEquals(len(result), len(new_set))
        self.assertEquals(result, new_set)



class StringSelectorBaseTests(unittest.TestCase):

    def test_constructor(self):
        """
        instantiate an StringSelector()
        """
        stringSel = StringSelector()
        self.assertEquals(stringSel.__class__, StringSelector)


class StringSelectorTests(StringSelectorBaseTests):

    def setUp(self):
        self.mols = Read('Data/stringSel.pdb')
        self.mol = self.mols[0]
        self.mol.buildBondsByDistance()
        self.stringSel = StringSelector()
    

    def tearDown(self):
        """
        clean-up
        """
        del(self.mol)


    #tests with 1 change of level
    def test_1level_select_with_empty_string(self):
        """
        test result with empty string ":" returns all chains
        """
        #stringSel = MVProteinSetSelector(self.mols, self.selString)
        #result, msg =  stringSel.go()
        result, msg = self.stringSel.select(self.mols, ":")
        self.assertEquals(result, self.mols.chains)
        self.assertEquals(result.__class__, self.mols.chains.__class__)



    def test_1level_select_end(self):
        """
         test select with '$'  returns last molecule.chains
        """
        new_mols = Read("Data/1crn.pdb")
        self.mols +=new_mols
        result, msg = self.stringSel.select(self.mols, "$:")
        self.assertEquals(result, self.mols[-1:].chains)


    def test_1level_select_end_2(self):
        """
         test select with ':$'  returns molecule's last chain
        """
        new_mols = Read("Data/1crn.pdb")
        self.mols +=new_mols
        result, msg = self.stringSel.select(self.mols, ":$")
        self.assertEquals(result[0], self.mols[-1:].chains[-1])


    def test_1level_select_with_valid_index(self):
        """
         test string with valid_index returns set with 1 item 
        """
        selString = "0:"
        #selString = "1:"
        result, msg = self.stringSel.select(self.mols, selString)
        self.assertEquals(len(result), len(self.mols[0].chains))
        self.assertEquals(result.__class__, self.mols.chains.__class__)


    def test_1level_select_with_invalid_index_returns_empty_set(self):
        """
         test string with invalid_index returns empty ChainSet
       FIX THIS: should it be an empty ChainSet?

        """
        selString = "2:"
        result, msg = self.stringSel.select(self.mols, selString)
        self.assertEqual(result.__class__, ChainSet)
        #self.assertEqual(result.__class__, self.mols.chains.__class__)
        self.assertEqual(len(result), 0)
        self.assertEqual(msg[0], selString[0])


    def test_1level_select_with_valid_range(self):
        """
         test string with valid_range returns ChainSet with 2 items
        """
        new_mols = Read("Data/1crn.pdb")
        self.mols +=new_mols
        selString = "0-1:"
        #selString = "1-2:"
        result, msg = self.stringSel.select(self.mols, selString)
        self.assertEquals(len(result), len(self.mols.chains))


    def test_1level_select_with_invalid_range(self):
        """
         test string with invalid_range returns set with 0 items 
       FIX THIS: should it be an empty ChainSet?
        """
        new_mols = Read("Data/1crn.pdb")
        self.mols +=new_mols
        selString = "4-6:"
        result, msg = self.stringSel.select(self.mols, selString)
        self.assertEquals(len(result), 0)
        self.assertEquals(result.__class__, ChainSet)
        #self.assertEquals(result.__class__, self.mols.chains.__class__)


    def test_1level_select_with_valid_regex(self):
        """
         test string with valid_regex returns set with 1 items
         <this regex is intended to match 1, 2,or 3 followed by anything>
        """
        new_mols = Read("Data/1crn.pdb")
        self.mols +=new_mols
        selString = "[1-3]*:"
        result, msg = self.stringSel.select(self.mols, selString)
        #print "result=", result
        #1crn has 1 chain
        self.assertEquals(len(result), 1)
        self.assertEquals(result.__class__, self.mols.chains.__class__)


    def test_1level_select_with_valid_regex_2(self):
        """
         test string with valid_regex returns set with 1 items
         <this regex is intended to match anything in range s-z
         followed by anything>
        """
        new_mols = Read("Data/1crn.pdb")
        self.mols +=new_mols
        selString = "[s-z]*:"
        result, msg = self.stringSel.select(self.mols, selString)
        #print "result=", result
        self.assertEquals(len(result), len(self.mols[0].chains))
        self.assertEquals(result.__class__, self.mols.chains.__class__)


    #tests with 2 changes of level
    def test_2level_select_with_empty_string(self):
        """
        test result with empty string "::" returns all residues
        """
        #stringSel = MVProteinSetSelector(self.mols, self.selString)
        #result, msg =  stringSel.go()
        result, msg = self.stringSel.select(self.mols, "::")
        self.assertEquals(result, self.mols.chains.residues)
        self.assertEquals(result.__class__, self.mols.chains.residues.__class__)


    def test_2level_select_end(self):
        """
         test select with '$'  returns last molecule.residues
        """
        new_mols = Read("Data/1crn.pdb")
        self.mols +=new_mols
        result, msg = self.stringSel.select(self.mols, "$::")
        self.assertEquals(result, self.mols[-1:].chains.residues)


    def test_2level_select_end_2(self):
        """
         test select with ':$'  returns molecule's last residue
        """
        new_mols = Read("Data/1crn.pdb")
        self.mols +=new_mols
        result, msg = self.stringSel.select(self.mols, "::$")
        self.assertEquals(result[0], self.mols[-1:].chains.residues[-1])


    def test_2level_select_with_valid_index(self):
        """
        test string with valid_index:: returns residue set of 54 
        """
        selString = "0::"
        #selString = "1::"
        result, msg = self.stringSel.select(self.mols, selString)
        self.assertEquals(len(result), len(self.mols.chains.residues))
        self.assertEquals(result.__class__, self.mols.chains.residues.__class__)


    def test_2level_select_with_invalid_index_returns_empty_set(self):
        """
        test string with invalid_index:: returns empty ResidueSet
       FIX THIS: should it be an empty residueSet?

        """
        selString = "2::"
        result, msg = self.stringSel.select(self.mols, selString)
        self.assertEqual(result.__class__, ResidueSet)
        #self.assertEqual(result.__class__, self.mols.chains.residues.__class__)
        self.assertEqual(len(result), 0)
        self.assertEqual(msg[0], selString[0])


    def test_2level_select_with_valid_range(self):
        """
         test string with valid_range returns all residues 
        """
        new_mols = Read("Data/1crn.pdb")
        self.mols +=new_mols
        selString = "0-1::"
        #selString = "1-2::"
        result, msg = self.stringSel.select(self.mols, selString)
        self.assertEquals(len(result), len(self.mols.chains.residues))


    def test_2level_select_with_invalid_range(self):
        """
        test string with invalid_range returns set with 0 items 
       FIX THIS: should it be an empty ResidueSet?
        """
        new_mols = Read("Data/1crn.pdb")
        self.mols +=new_mols
        selString = "4-6::"
        result, msg = self.stringSel.select(self.mols, selString)
        self.assertEquals(len(result), 0)
        self.assertEquals(result.__class__, ResidueSet)
        #self.assertEquals(result.__class__, self.mols.chains.residues.__class__)


    def test_2level_select_with_valid_regex(self):
        """
         test string with valid_regex:: returns set with 46 residues
         <this regex is intended to match 1, 2,or 3 followed by anything>
        """
        new_mols = Read("Data/1crn.pdb")
        self.mols +=new_mols
        selString = "[1-3]*::"
        result, msg = self.stringSel.select(self.mols, selString)
        #print "result=", result
        #result matches 1crn only
        self.assertEquals(len(result), len(self.mols[1].chains.residues))
        self.assertEquals(result.__class__, self.mols.chains.residues.__class__)


    def test_2level_select_with_valid_regex_2(self):
        """
         test string with valid_regex:: returns residueSet of 8
         <this regex is intended to match anything in range s-z
         followed by anything>
        """
        new_mols = Read("Data/1crn.pdb")
        self.mols +=new_mols
        selString = "[s-z]*::"
        result, msg = self.stringSel.select(self.mols, selString)
        #print "result=", result
        #result matches mols[0] only
        self.assertEquals(len(result), len(self.mols[0].chains.residues))
        self.assertEquals(result.__class__, self.mols.chains.residues.__class__)


    #tests with 3 changes of level
    def test_3level_select_with_empty_string(self):
        """
        test result with empty string ":::" returns all atoms
        """
        #stringSel = MVProteinSetSelector(self.mols, self.selString)
        #result, msg =  stringSel.go()
        result, msg = self.stringSel.select(self.mols, ":::")
        self.assertEquals(result, self.mols.allAtoms)
        self.assertEquals(result.__class__, self.mols.allAtoms.__class__)


    def test_3level_select_end(self):
        """
         test select with '$'  returns last molecule.atoms
        """
        new_mols = Read("Data/1crn.pdb")
        self.mols +=new_mols
        result, msg = self.stringSel.select(self.mols, "$:::")
        self.assertEquals(result, self.mols[-1:].allAtoms)


    def test_3level_select_end_2(self):
        """
         test select with ':$'  returns molecule's last atom
        """
        new_mols = Read("Data/1crn.pdb")
        self.mols +=new_mols
        result, msg = self.stringSel.select(self.mols, ":::$")
        self.assertEquals(result[0], self.mols[-1:].allAtoms[-1])


    def test_3level_select_with_valid_index(self):
        """
         test string with valid_index returns set of atoms
        """
        selString = "0:::"
        #selString = "1:::"
        result, msg = self.stringSel.select(self.mols, selString)
        #print "result=", result
        self.assertEquals(len(result), len(self.mols[0].allAtoms))
        self.assertEquals(result.__class__, self.mols.allAtoms.__class__)


    def test_3level_select_with_invalid_index_returns_empty_set(self):
        """
         test string with invalid_index returns empty AtomSet
       FIX THIS: should it be an empty AtomSet?

        """
        selString = "2:::"
        result, msg = self.stringSel.select(self.mols, selString)
        self.assertEqual(result.__class__, AtomSet)
        self.assertEqual(len(result), 0)
        self.assertEqual(msg[0], selString[0])


    def test_3level_select_with_valid_range(self):
        """
         test string with valid_range returns 1crn.allAtoms 
        """
        new_mols = Read("Data/1crn.pdb")
        self.mols +=new_mols
        selString = "0-1:::"
        #selString = "1-2:::"
        result, msg = self.stringSel.select(self.mols, selString)
        self.assertEquals(len(result), len(self.mols.allAtoms))


    def test_3level_select_with_invalid_range(self):
        """
         test string with invalid_range returns set with 0 items 
       FIX THIS: should it be an empty AtomSet?
        """
        new_mols = Read("Data/1crn.pdb")
        self.mols +=new_mols
        selString = "4-6:::"
        result, msg = self.stringSel.select(self.mols, selString)
        self.assertEquals(len(result), 0)
        self.assertEquals(result.__class__, AtomSet)
        #self.assertEquals(result.__class__, self.mols.allAtoms.__class__)


    def test_3level_select_with_valid_regex(self):
        """
         test string with valid_regex returns AtomSet 1crn.allAtoms
         <this regex is intended to match 1, 2,or 3 followed by anything>
        """
        new_mols = Read("Data/1crn.pdb")
        self.mols +=new_mols
        selString = "[1-3]*:::"
        result, msg = self.stringSel.select(self.mols, selString)
        #print "result=", result
        #result matches 1crn only
        self.assertEquals(len(result), len(self.mols[1].allAtoms))
        self.assertEquals(result.__class__, self.mols.allAtoms.__class__)


    def test_3level_select_with_valid_regex_2(self):
        """
         test string with valid_regex returns set with 1 items
         <this regex is intended to match anything in range s-z
         followed by anything>
        """
        new_mols = Read("Data/1crn.pdb")
        self.mols +=new_mols
        selString = "[s-z]*:::"
        result, msg = self.stringSel.select(self.mols, selString)
        #print "result=", result
        #result matches mols[0] only
        self.assertEquals(len(result), len(self.mols[0].allAtoms))
        self.assertEquals(result.__class__, self.mols.allAtoms.__class__)


class BondSetSelectorBaseTests(unittest.TestCase):

    def test_constructor(self):
        """
        instantiate an BondSetSelector
        """
        stringSel = BondSetSelector()
        self.assertEquals(stringSel.__class__, BondSetSelector)



class BondSetSelectorTests(BondSetSelectorBaseTests):


    def setUp(self):
        self.mols = Read('Data/stringSel.pdb')
        for m in self.mols:
            m.buildBondsByDistance()
        self.bonds = self.mols.chains.residues.atoms.bonds[0]
        self.stringSel = BondSetSelector()
    

    def tearDown(self):
        """
        clean-up
        """
        del(self.bonds)
        del(self.mols)


    def test_select_with_empty_string(self):
        """
         test result with empty string returns all bonds
        """
        result, msg = self.stringSel.select(self.bonds, "")
        self.assertEquals(result, self.bonds)


    def test_select_end(self):
        """
         test select with '$'  returns last item
        """
        result, msg = self.stringSel.select(self.bonds, "$")
        self.assertEquals(result[-1], self.bonds[-1])


    def test_select_with_valid_index(self):
        """
         test string with valid_index returns set with 1 item 
        """
        selString = "0"
        result, msg = self.stringSel.select(self.bonds, selString)
        self.assertEquals(len(result), 1)
        self.assertEquals(result[0], self.bonds[0])


    def test_select_with_invalid_index_returns_empty_set(self):
        """
         test string with invalid_index returns empty set
        """
        selString = "400"
        result, msg = self.stringSel.select(self.bonds, selString)
        self.assertEqual(result.__class__, self.bonds.__class__)
        self.assertEqual(len(result), 0)
        self.assertEqual(msg[0], selString)


    def test_select_with_valid_range(self):
        """
         test string with valid_range returns set with 2 items
        """
        selString = "0-1"
        result, msg = self.stringSel.select(self.bonds, selString)
        self.assertEquals(len(result), 2)


    def test_select_with_invalid_range(self):
        """
         test string with invalid_range returns set with 0 items 
        """
        selString = "400-623"
        result, msg = self.stringSel.select(self.bonds, selString)
        self.assertEquals(len(result), 0)


    def test_select_with_valid_lambda_exp(self):
        """
         test string with valid_lambda_exp returns set with items
         <this regex is intended to select all bonds involving NE2>
        """
        selString = "lambda x: x.atom1.name=='NE2' or x.atom2.name=='NE2'"
        result, msg = self.stringSel.select(self.bonds, selString)
        #print "result=", result
        self.assertEquals(len(result), 1)


    def test_no_select_with_valid_regex_2(self):
        """
         bondSelector does not support regex
         test string with valid_regex returns set with 0 items
        """
        selString = "[1-2]"
        result, msg = self.stringSel.select(self.bonds, selString)
        self.assertEquals(len(result), 0)


    def test_select_valid_set(self):
        """
         test string with valid set name returns 1 set
        """
        these_sets = Sets()
        key = 'first'
        these_sets.add(key, self.bonds[1:])
        selString = key
        result, msg = self.stringSel.select(self.bonds, selString, 
                                        sets=these_sets)
        #print "result=", result
        self.assertEquals(len(result), 23)
        self.assertEquals(result, self.bonds[1:])




if __name__ == '__main__':
    #unittest.main()
    test_cases = [
        'ProteinSetSelectorTests',
        'ChainSetSelectorTests',
        'ResidueSetSelectorTests',
        'AtomSetSelectorTests',
        'StringSelectorTests',
        'BondSetSelectorTests',
    ]
    unittest.main( argv=([__name__,] + test_cases) )
