
#
#
#
#
# $Id: test_selection.py,v 1.10 2010/09/10 19:31:49 sanner Exp $
#

import unittest
from string import split
from MolKit.stringSelector import StringSelector, CompoundStringSelector
from MolKit.molecule import Atom, AtomSet, MoleculeSet
from MolKit.protein import Protein, Chain, Residue, ProteinSet, ChainSet, ResidueSet


molecules1 = None
molecules2 = None
molecules3 = None

class SelectionBaseTest(unittest.TestCase):

    def initMolecules(self):
        global molecules1, molecules2, molecules3
        if molecules1 is None:
            #print "reading hsg1"
            from MolKit import Read
            molecules1 = Read('Data/hsg1.pdb')
        self.molecules1 = molecules1
        if molecules2 is None:
            #print "reading 2plv"
            from MolKit import Read
            molecules2 = Read('Data/2plv.pdb')
        self.molecules2 = molecules2
        if molecules3 is None:
            #print "reading 1gyc"
            from MolKit import Read
            molecules3 = Read('Data/1gyc.pdb')
        self.molecules3 = molecules3
        self.molecules = molecules1+molecules2+molecules3

    def setUp(self):
        if not hasattr(self, 'molecules'):
            self.initMolecules()


class MoleculeSetSelectionTests(SelectionBaseTest):
    #singleton tests: regexp, relative, index, NamedResSet, NamedAtomSet

    def setUp(self):
        if not hasattr(self, 'molecules'):
            self.initMolecules()

    #simple/singleton
    #   regexp: 
    def test_MoleculeSet_regexp(self):
        """
         test regexp used in Mols.get(regexp)
        """
        result = self.molecules.get("*g*")
        self.assertEquals(result, self.molecules1+self.molecules3)
        self.assertEquals(result.stringRepr,"(hsg1/+/2plv/+/1gyc\\s\\*g*)")


    #   relative: 
    def test_MoleculeSet_relative_index(self):
        """
         test relative index used in Mols.get(relative_index)
        """
        result = self.molecules.get("#1")
        self.assertEquals(result, self.molecules2)
        self.assertEquals(result.stringRepr,"(hsg1/+/2plv/+/1gyc\\s\\#1)")


    #   index: 
    def test_MoleculeSet_index(self):
        """
         test index used in Mols.get(index)
        """
        result = self.molecules.get("1")
        self.assertEquals(result, self.molecules2)
        self.assertEquals(result.stringRepr,"(hsg1/+/2plv/+/1gyc\\s\\1)")


    #   NamedChainSet:  n/a
    def test_MoleculeSet_NamedChainSet(self):
        """
         test NamedResSet used in Mols.get(NamedChainSet)
        """
        result = self.molecules.get("proteic")
        self.assertEquals(result, MoleculeSet())

    #   NamedResSet:  n/a
    def test_MoleculeSet_NamedResSet(self):
        """
         test NamedResSet used in Mols.get(NamedResSet)
        """
        result = self.molecules.get("acidic")
        self.assertEquals(result, MoleculeSet())


    #   NamedAtomSet: n/a
    def test_MoleculeSet_NamedAtomSet(self):
        """
         test NamedAtomSet used in Mols.get(NamedAtomSet)
        """
        result = self.molecules.get("backbone")
        self.assertEquals(result, MoleculeSet())


    #compound
    # !!!!NOTE!!!! RANGE includes both ends!!!
    # RANGE: REGEXP-REGEXP
    #       regexp -  regexp
    #       regexp -  regexp(fail)
    def test_MoleculeSet_range_regexp_regexp(self):
        """
         test range regexp-regexp used in Mols.get(regexp-regexp)
        """
        result = self.molecules.get("*1-1*")
        self.assertEquals(result, self.molecules)
        self.assertEquals(result.stringRepr,"(hsg1/+/2plv/+/1gyc\\s\\*1-1*)")

    def test_MoleculeSet_range_regexp_regexp_w_spaces(self):
        """
         test range regexp-regexp_w_spaces used in Mols.get(regexp-regexp)
        """
        result = self.molecules.get("*1 - 1*")
        self.assertEquals(result, self.molecules)
        self.assertEquals(result.stringRepr,"(hsg1/+/2plv/+/1gyc\\s\\*1 - 1*)")

    def test_MoleculeSet_range_regexp_regexp_fail_second(self):
        """
         test range regexp-regexp_fail_second used in Mols.get(regexp-regexp)
        """
        result = self.molecules.get("*1-6*")
        self.assertEquals(result, MoleculeSet())

    #       regexp(fail) -  regexp
    def test_MoleculeSet_range_regexp_fail_first_regexp(self):
        """
         test range regexp_fail_first-regexp used in Mols.get(regexp-regexp)
        """
        result = self.molecules.get("*6-1*")
        self.assertEquals(result, MoleculeSet())

    #       regexp(fail) -  regexp(fail)
    def test_MoleculeSet_range_regexp_regexp_fail_both(self):
        """
         test range regexp-regexp_fail_both used in Mols.get(regexp-regexp)
        """
        result = self.molecules.get("*6-8*")
        self.assertEquals(result, MoleculeSet())

    # RANGE: REGEXP-RELATIVE
    #       regexp -  relative
    #       regexp -  relative(fail)
    def test_MoleculeSet_range_regexp_relative(self):
        """
         test range regexp-relative used in Mols.get(regexp-relative)
        """
        result = self.molecules.get("hs*-#2")
        self.assertEquals(result, self.molecules)
        self.assertEquals(result.stringRepr,"(hsg1/+/2plv/+/1gyc\\s\\hs*-#2)")

    def test_MoleculeSet_range_regexp_relative_w_spaces(self):
        """
         test range regexp-relative_w_spaces used in Mols.get(regexp-relative)
        """
        result = self.molecules.get("hs* - #2")
        self.assertEquals(result, self.molecules)
        self.assertEquals(result.stringRepr,"(hsg1/+/2plv/+/1gyc\\s\\hs* - #2)")

    def test_MoleculeSet_range_regexp_relative_fail_second(self):
        """
         test range regexp-relative_fail_second used in Mols.get(regexp-relative)
        """
        result = self.molecules.get("*1-#1*")
        self.assertEquals(result, MoleculeSet())

    #       regexp(fail) -  relative
    def test_MoleculeSet_range_regexp_fail_first_relative(self):
        """
         test range regexp_fail_first-relative used in Mols.get(regexp-relative)
        """
        result = self.molecules.get("*6-#1")
        self.assertEquals(result, MoleculeSet())

    #       regexp(fail) -  relative(fail)
    def test_MoleculeSet_range_regexp_relative_fail_both(self):
        """
         test range regexp-relative_fail_both used in Mols.get(regexp-relative)
        """
        result = self.molecules.get("*6-#8")
        self.assertEquals(result, MoleculeSet())

    # RANGE: REGEXP-INDEX
    #       regexp -  index
    #       regexp -  index(fail)
    def test_MoleculeSet_range_regexp_index(self):
        """
         test range regexp-index used in Mols.get(regexp-index)
        """
        result = self.molecules.get("hs*-2")
        self.assertEquals(result, self.molecules)
        self.assertEquals(result.stringRepr,"(hsg1/+/2plv/+/1gyc\\s\\hs*-2)")

    def test_MoleculeSet_range_regexp_index_w_spaces(self):
        """
         test range regexp-index_w_spaces used in Mols.get(regexp-index)
        """
        result = self.molecules.get("hs* - 2")
        self.assertEquals(result, self.molecules)
        self.assertEquals(result.stringRepr,"(hsg1/+/2plv/+/1gyc\\s\\hs* - 2)")

    def test_MoleculeSet_range_regexp_index_fail_second(self):
        """
         test range regexp-index_fail_second used in Mols.get(regexp-index)
        """
        result = self.molecules.get("*1-3")
        self.assertEquals(result, MoleculeSet())

    #       regexp(fail) -  index
    def test_MoleculeSet_range_regexp_fail_first_index(self):
        """
         test range regexp_fail_first-index used in Mols.get(regexp-index)
        """
        result = self.molecules.get("*6-1")
        self.assertEquals(result, MoleculeSet())

    #       regexp(fail) -  index(fail)
    def test_MoleculeSet_range_regexp_index_fail_both(self):
        """
         test range regexp-index_fail_both used in Mols.get(regexp-index)
        """
        result = self.molecules.get("*6-8")
        self.assertEquals(result, MoleculeSet())

    # RANGE: REGEXP-NamedChainSet
    #       regexp -  NamedChainSet
    #       regexp -  NamedChainSet(fail)
    def test_MoleculeSet_range_regexp_NamedChainSet(self):
        """
         test range regexp-NamedChainSet used in Mols.get(regexp-NamedChainSet)
        """
        result = self.molecules.get("hs*-dna")
        self.assertEquals(result, MoleculeSet())

    def test_MoleculeSet_range_regexp_NamedChainSet_w_spaces(self):
        """
         test range regexp-NamedChainSet_w_spaces used in Mols.get(regexp-NamedChainSet)
        """
        result = self.molecules.get("hs* - dna")
        self.assertEquals(result, MoleculeSet())

    def test_MoleculeSet_range_regexp_NamedChainSet_fail_second(self):
        """
         test range regexp-NamedChainSet_fail_second used in Mols.get(regexp-NamedChainSet)
        """
        result = self.molecules.get("*1-mistake")
        self.assertEquals(result, MoleculeSet())

    #       regexp(fail) -  NamedChainSet
    def test_MoleculeSet_range_regexp_fail_first_NamedChainSet(self):
        """
         test range regexp_fail_first-NamedChainSet used in Mols.get(regexp-NamedChainSet)
        """
        result = self.molecules.get("*6-dna")
        self.assertEquals(result, MoleculeSet())

    #       regexp(fail) -  NamedChainSet(fail)
    def test_MoleculeSet_range_regexp_NamedChainSet_fail_both(self):
        """
         test range regexp-NamedChainSet_fail_both used in Mols.get(regexp-NamedChainSet)
        """
        result = self.molecules.get("*6-mistake")
        self.assertEquals(result, MoleculeSet())


    # RANGE: REGEXP-NamedResSet
    #       regexp -  NamedResSet
    #       regexp -  NamedResSet(fail)
    def test_MoleculeSet_range_regexp_NamedResSet(self):
        """
         test range regexp-NamedResSet used in Mols.get(regexp-NamedResSet)
        """
        result = self.molecules.get("hs*-buried")
        self.assertEquals(result, MoleculeSet())

    def test_MoleculeSet_range_regexp_NamedResSet_w_spaces(self):
        """
         test range regexp-NamedResSet_w_spaces used in Mols.get(regexp-NamedResSet)
        """
        result = self.molecules.get("hs* - buried")
        self.assertEquals(result, MoleculeSet())

    def test_MoleculeSet_range_regexp_NamedResSet_fail_second(self):
        """
         test range regexp-NamedResSet_fail_second used in Mols.get(regexp-NamedResSet)
        """
        result = self.molecules.get("*1-mistake")
        self.assertEquals(result, MoleculeSet())

    #       regexp(fail) -  NamedResSet
    def test_MoleculeSet_range_regexp_fail_first_NamedResSet(self):
        """
         test range regexp_fail_first-NamedResSet used in Mols.get(regexp-NamedResSet)
        """
        result = self.molecules.get("*6-buried")
        self.assertEquals(result, MoleculeSet())

    #       regexp(fail) -  NamedResSet(fail)
    def test_MoleculeSet_range_regexp_NamedResSet_fail_both(self):
        """
         test range regexp-NamedResSet_fail_both used in Mols.get(regexp-NamedResSet)
        """
        result = self.molecules.get("*6-mistake")
        self.assertEquals(result, MoleculeSet())


    # RANGE: REGEXP-NamedAtomSet
    #       regexp -  NamedAtomSet
    #       regexp -  NamedAtomSet(fail)
    def test_MoleculeSet_range_regexp_NamedAtomSet(self):
        """
         test range regexp-NamedAtomSet used in Mols.get(regexp-NamedAtomSet)
        """
        result = self.molecules.get("hs*-backbone")
        self.assertEquals(result, MoleculeSet())

    def test_MoleculeSet_range_regexp_NamedAtomSet_w_spaces(self):
        """
         test range regexp-NamedAtomSet_w_spaces used in Mols.get(regexp-NamedAtomSet)
        """
        result = self.molecules.get("hs* - backbone")
        self.assertEquals(result, MoleculeSet())

    def test_MoleculeSet_range_regexp_NamedAtomSet_fail_second(self):
        """
         test range regexp-NamedAtomSet_fail_second used in Mols.get(regexp-NamedAtomSet)
        """
        result = self.molecules.get("*1-mistake")
        self.assertEquals(result, MoleculeSet())

    #       regexp(fail) -  NamedAtomSet
    def test_MoleculeSet_range_regexp_fail_first_NamedAtomSet(self):
        """
         test range regexp_fail_first-NamedAtomSet used in Mols.get(regexp-NamedAtomSet)
        """
        result = self.molecules.get("*6-backbone")
        self.assertEquals(result, MoleculeSet())

    #       regexp(fail) -  NamedAtomSet(fail)
    def test_MoleculeSet_range_regexp_NamedAtomSet_fail_both(self):
        """
         test range regexp-NamedAtomSet_fail_both used in Mols.get(regexp-NamedAtomSet)
        """
        result = self.molecules.get("*6-mistake")
        self.assertEquals(result, MoleculeSet())

    #
    # RANGE: RELATIVE-REGEXP
    #       relative -  regexp
    #       relative -  regexp(fail)
    def test_MoleculeSet_range_relative_regexp(self):
        """
         test range relative-regexp used in Mols.get(relative-regexp)
        """
        result = self.molecules.get("#1-1*")
        self.assertEquals(result, self.molecules[1:])
        self.assertEquals(result.stringRepr,"(hsg1/+/2plv/+/1gyc\\s\\#1-1*)")

    def test_MoleculeSet_range_relative_regexp_w_spaces(self):
        """
         test range relative-regexp_w_spaces used in Mols.get(relative-regexp)
        """
        result = self.molecules.get("#1 - 1*")
        self.assertEquals(result, self.molecules[1:])
        self.assertEquals(result.stringRepr,"(hsg1/+/2plv/+/1gyc\\s\\#1 - 1*)")

    def test_MoleculeSet_range_relative_regexp_fail_second(self):
        """
         test range relative-regexp_fail_second used in Mols.get(relative-regexp)
        """
        result = self.molecules.get("#1-6*")
        self.assertEquals(result, MoleculeSet())

    #       relative(fail) -  regexp
    def test_MoleculeSet_range_relative_fail_first_regexp(self):
        """
         test range relative_fail_first-regexp used in Mols.get(relative-regexp)
        """
        result = self.molecules.get("#6-1*")
        self.assertEquals(result, MoleculeSet())

    #       relative(fail) -  regexp(fail)
    def test_MoleculeSet_range_relative_regexp_fail_both(self):
        """
         test range relative-regexp_fail_both used in Mols.get(relative-regexp)
        """
        result = self.molecules.get("#6-8*")
        self.assertEquals(result, MoleculeSet())

    # RANGE: RELATIVE-RELATIVE
    #       relative -  relative
    #       relative -  relative(fail)
    def test_MoleculeSet_range_relative_relative(self):
        """
         test range relative-relative used in Mols.get(relative-relative)
        """
        result = self.molecules.get("#1-#2")
        self.assertEquals(result, self.molecules[1:])
        self.assertEquals(result.stringRepr,"(hsg1/+/2plv/+/1gyc\\s\\#1-#2)")

    def test_MoleculeSet_range_relative_relative_w_spaces(self):
        """
         test range relative-relative_w_spaces used in Mols.get(relative-relative)
        """
        result = self.molecules.get("#1 - #2")
        self.assertEquals(result, self.molecules[1:])
        self.assertEquals(result.stringRepr,"(hsg1/+/2plv/+/1gyc\\s\\#1 - #2)")

    def test_MoleculeSet_range_relative_relative_fail_second(self):
        """
         test range relative-relative_fail_second used in Mols.get(relative-relative)
        """
        result = self.molecules.get("#1-#5")
        self.assertEquals(result, MoleculeSet())

    #       relative(fail) -  relative
    def test_MoleculeSet_range_relative_fail_first_relative(self):
        """
         test range relative_fail_first-relative used in Mols.get(relative-relative)
        """
        result = self.molecules.get("#6-#1")
        self.assertEquals(result, MoleculeSet())

    #       relative(fail) -  relative(fail)
    def test_MoleculeSet_range_relative_relative_fail_both(self):
        """
         test range relative-relative_fail_both used in Mols.get(relative-relative)
        """
        result = self.molecules.get("#6-#8")
        self.assertEquals(result, MoleculeSet())

    # RANGE: RELATIVE-INDEX
    #       relative -  index
    #       relative -  index(fail)
    def test_MoleculeSet_range_relative_index(self):
        """
         test range relative-index used in Mols.get(relative-index)
        """
        result = self.molecules.get("#1-2")
        self.assertEquals(result, self.molecules[1:])
        self.assertEquals(result.stringRepr,"(hsg1/+/2plv/+/1gyc\\s\\#1-2)")

    def test_MoleculeSet_range_relative_index_w_spaces(self):
        """
         test range relative-index_w_spaces used in Mols.get(relative-index)
        """
        result = self.molecules.get("#1 - 2")
        self.assertEquals(result, self.molecules[1:])
        self.assertEquals(result.stringRepr,"(hsg1/+/2plv/+/1gyc\\s\\#1 - 2)")

    def test_MoleculeSet_range_relative_index_fail_second(self):
        """
         test range relative-index_fail_second used in Mols.get(relative-index)
        """
        result = self.molecules.get("#1-3")
        self.assertEquals(result, MoleculeSet())

    #       relative(fail) -  index
    def test_MoleculeSet_range_relative_fail_first_index(self):
        """
         test range relative_fail_first-index used in Mols.get(relative-index)
        """
        result = self.molecules.get("#6-1")
        self.assertEquals(result, MoleculeSet())

    #       relative(fail) -  index(fail)
    def test_MoleculeSet_range_relative_index_fail_both(self):
        """
         test range relative_index-fail_both used in Mols.get(relative-index)
        """
        result = self.molecules.get("#6-8")
        self.assertEquals(result, MoleculeSet())

    # RANGE: RELATIVE-NamedResSet
    #       relative -  NamedResSet
    #       relative -  NamedResSet(fail)
    def test_MoleculeSet_range_relative_NamedResSet(self):
        """
         test range relative-NamedResSet used in Mols.get(relative-NamedResSet)
        """
        result = self.molecules.get("#1-buried")
        self.assertEquals(result, MoleculeSet())

    def test_MoleculeSet_range_relative_NamedResSet_w_spaces(self):
        """
         test range relative-NamedResSet_w_spaces used in Mols.get(relative-NamedResSet)
        """
        result = self.molecules.get("#1 - buried")
        self.assertEquals(result, MoleculeSet())

    def test_MoleculeSet_range_relative_NamedResSet_fail_second(self):
        """
         test range relative-NamedResSet_fail_second used in Mols.get(relative-NamedResSet)
        """
        result = self.molecules.get("#1-mistake")
        self.assertEquals(result, MoleculeSet())

    #       relative(fail) -  NamedResSet
    def test_MoleculeSet_range_relative_fail_first_NamedResSet(self):
        """
         test range relative_fail-first_NamedResSet used in Mols.get(relative-NamedResSet)
        """
        result = self.molecules.get("#6-buried")
        self.assertEquals(result, MoleculeSet())

    #       relative(fail) -  NamedResSet(fail)
    def test_MoleculeSet_range_relative_NamedResSet_fail_both(self):
        """
         test range relative-NamedResSet_fail_both used in Mols.get(relative-NamedResSet)
        """
        result = self.molecules.get("#6-mistake")
        self.assertEquals(result, MoleculeSet())


    # RANGE: RELATIVE-NamedAtomSet
    #       relative -  NamedAtomSet
    #       relative -  NamedAtomSet(fail)
    def test_MoleculeSet_range_relative_NamedAtomSet(self):
        """
         test range relative-NamedAtomSet used in Mols.get(relative-NamedAtomSet)
        """
        result = self.molecules.get("#1-backbone")
        self.assertEquals(result, MoleculeSet())

    def test_MoleculeSet_range_relative_NamedAtomSet_w_spaces(self):
        """
         test range relative-NamedAtomSet_w_spaces used in Mols.get(relative-NamedAtomSet)
        """
        result = self.molecules.get("#1 - backbone")
        self.assertEquals(result, MoleculeSet())

    def test_MoleculeSet_range_relative_NamedAtomSet_fail_second(self):
        """
         test range relative-NamedAtomSet_fail_second used in Mols.get(relative-NamedAtomSet)
        """
        result = self.molecules.get("#1-mistake")
        self.assertEquals(result, MoleculeSet())

    #       relative(fail) -  NamedAtomSet
    def test_MoleculeSet_range_relative_fail_first_NamedAtomSet(self):
        """
         test range relative_fail_first-NamedAtomSet used in Mols.get(relative-NamedAtomSet)
        """
        result = self.molecules.get("#6-backbone")
        self.assertEquals(result, MoleculeSet())

    #       relative(fail) -  NamedAtomSet(fail)
    def test_MoleculeSet_range_relative_NamedAtomSet_fail_both(self):
        """
         test range relative-NamedAtomSet_fail_both used in Mols.get(relative-NamedAtomSet)
        """
        result = self.molecules.get("#6-mistake")
        self.assertEquals(result, MoleculeSet())

    #
    # RANGE: INDEX-REGEXP
    #       index -  regexp
    #       index -  regexp(fail)
    def test_MoleculeSet_range_index_regexp(self):
        """
         test range index-regexp used in Mols.get(index-regexp)
        """
        result = self.molecules.get("1-1*")
        self.assertEquals(result, self.molecules[1:])
        self.assertEquals(result.stringRepr,"(hsg1/+/2plv/+/1gyc\\s\\1-1*)")

    def test_MoleculeSet_range_index_regexp_w_spaces(self):
        """
         test range index-regexp_w_spaces used in Mols.get(index-regexp)
        """
        result = self.molecules.get("1 - 1*")
        self.assertEquals(result, self.molecules[1:])
        self.assertEquals(result.stringRepr,"(hsg1/+/2plv/+/1gyc\\s\\1 - 1*)")

    def test_MoleculeSet_range_index_regexp_fail_second(self):
        """
         test range index-regexp_fail_second used in Mols.get(index-regexp)
        """
        result = self.molecules.get("1-6*")
        self.assertEquals(result, MoleculeSet())

    #       index(fail) -  regexp
    def test_MoleculeSet_range_index_fail_first_regexp(self):
        """
         test range index_fail_first-regexp used in Mols.get(index-regexp)
        """
        result = self.molecules.get("6-1*")
        self.assertEquals(result, MoleculeSet())

    #       index(fail) -  regexp(fail)
    def test_MoleculeSet_range_index_regexp_fail_both(self):
        """
         test range index-regexp_fail_both used in Mols.get(index-regexp)
        """
        result = self.molecules.get("6-8*")
        self.assertEquals(result, MoleculeSet())

    # RANGE: INDEX-RELATIVE
    #       index -  relative
    #       index -  relative(fail)
    def test_MoleculeSet_range_index_relative(self):
        """
         test range index-relative used in Mols.get(index-relative)
        """
        result = self.molecules.get("1-#2")
        self.assertEquals(result, self.molecules[1:])
        self.assertEquals(result.stringRepr,"(hsg1/+/2plv/+/1gyc\\s\\1-#2)")

    def test_MoleculeSet_range_index_relative_w_spaces(self):
        """
         test range index-relative_w_spaces used in Mols.get(index-relative)
        """
        result = self.molecules.get("1 - #2")
        self.assertEquals(result, self.molecules[1:])
        self.assertEquals(result.stringRepr,"(hsg1/+/2plv/+/1gyc\\s\\1 - #2)")

    def test_MoleculeSet_range_index_relative_fail_second(self):
        """
         test range index-relative_fail_second used in Mols.get(index-relative)
        """
        result = self.molecules.get("1-#5")
        self.assertEquals(result, MoleculeSet())

    #       index(fail) -  relative
    def test_MoleculeSet_range_index_fail_first_relative(self):
        """
         test range index_fail_first-relative used in Mols.get(index-relative)
        """
        result = self.molecules.get("6-#1")
        self.assertEquals(result, MoleculeSet())

    #       index(fail) -  relative(fail)
    def test_MoleculeSet_range_index_relative_fail_both(self):
        """
         test range index-relative_fail_both used in Mols.get(index-relative)
        """
        result = self.molecules.get("6-#8")
        self.assertEquals(result, MoleculeSet())

    # RANGE: INDEX-INDEX
    #       index -  index
    #       index -  index(fail)
    def test_MoleculeSet_range_index_index(self):
        """
         test range index-index used in Mols.get(index-index)
        """
        result = self.molecules.get("1-2")
        self.assertEquals(result, self.molecules[1:])
        self.assertEquals(result.stringRepr,"(hsg1/+/2plv/+/1gyc\\s\\1-2)")

    def test_MoleculeSet_range_index_index_w_spaces(self):
        """
         test range index-index_w_spaces used in Mols.get(index-index)
        """
        result = self.molecules.get("1 - 2")
        self.assertEquals(result, self.molecules[1:])
        self.assertEquals(result.stringRepr,"(hsg1/+/2plv/+/1gyc\\s\\1 - 2)")

    def test_MoleculeSet_range_index_index_fail_second(self):
        """
         test range index-index_fail_second used in Mols.get(index-index)
        """
        result = self.molecules.get("1-3")
        self.assertEquals(result, MoleculeSet())

    #       index(fail) -  index
    def test_MoleculeSet_range_index_fail_first_index(self):
        """
         test range index_fail_first-index used in Mols.get(index-index)
        """
        result = self.molecules.get("6-1")
        self.assertEquals(result, MoleculeSet())

    #       index(fail) -  index(fail)
    def test_MoleculeSet_range_index_index_fail_both(self):
        """
         test range index_index-fail_both used in Mols.get(index-index)
        """
        result = self.molecules.get("6-8")
        self.assertEquals(result, MoleculeSet())

    # RANGE: INDEX-NamedResSet
    #       index -  NamedResSet
    #       index -  NamedResSet(fail)
    def test_MoleculeSet_range_index_NamedResSet(self):
        """
         test range index-NamedResSet used in Mols.get(index-NamedResSet)
        """
        result = self.molecules.get("1-buried")
        self.assertEquals(result, MoleculeSet())

    def test_MoleculeSet_range_index_NamedResSet_w_spaces(self):
        """
         test range index-NamedResSet_w_spaces used in Mols.get(index-NamedResSet)
        """
        result = self.molecules.get("1 - buried")
        self.assertEquals(result, MoleculeSet())

    def test_MoleculeSet_range_index_NamedResSet_fail_second(self):
        """
         test range index-NamedResSet_fail_second used in Mols.get(index-NamedResSet)
        """
        result = self.molecules.get("1-mistake")
        self.assertEquals(result, MoleculeSet())

    #       index(fail) -  NamedResSet
    def test_MoleculeSet_range_index_fail_first_NamedResSet(self):
        """
         test range index_fail-first_NamedResSet used in Mols.get(index-NamedResSet)
        """
        result = self.molecules.get("6-buried")
        self.assertEquals(result, MoleculeSet())

    #       index(fail) -  NamedResSet(fail)
    def test_MoleculeSet_range_index_NamedResSet_fail_both(self):
        """
         test range index-NamedResSet_fail_both used in Mols.get(index-NamedResSet)
        """
        result = self.molecules.get("6-mistake")
        self.assertEquals(result, MoleculeSet())


    # RANGE: INDEX-NamedAtomSet
    #       index -  NamedAtomSet
    #       index -  NamedAtomSet(fail)
    def test_MoleculeSet_range_index_NamedAtomSet(self):
        """
         test range index-NamedAtomSet used in Mols.get(index-NamedAtomSet)
        """
        result = self.molecules.get("1-backbone")
        self.assertEquals(result, MoleculeSet())

    def test_MoleculeSet_range_index_NamedAtomSet_w_spaces(self):
        """
         test range index-NamedAtomSet_w_spaces used in Mols.get(index-NamedAtomSet)
        """
        result = self.molecules.get("1 - backbone")
        self.assertEquals(result, MoleculeSet())

    def test_MoleculeSet_range_index_NamedAtomSet_fail_second(self):
        """
         test range index-NamedAtomSet_fail_second used in Mols.get(index-NamedAtomSet)
        """
        result = self.molecules.get("1-mistake")
        self.assertEquals(result, MoleculeSet())

    #       index(fail) -  NamedAtomSet
    def test_MoleculeSet_range_index_fail_first_NamedAtomSet(self):
        """
         test range index_fail_first-NamedAtomSet used in Mols.get(index-NamedAtomSet)
        """
        result = self.molecules.get("6-backbone")
        self.assertEquals(result, MoleculeSet())

    #       index(fail) -  NamedAtomSet(fail)
    def test_MoleculeSet_range_index_NamedAtomSet_fail_both(self):
        """
         test range index-NamedAtomSet_fail_both used in Mols.get(index-NamedAtomSet)
        """
        result = self.molecules.get("6-mistake")
        self.assertEquals(result, MoleculeSet())



class ChainSetSelectionTests(SelectionBaseTest):
    #singleton tests: regexp, relative, index, NamedChainSet, NamedResSet, NamedAtomSet

    def setUp(self):
        if not hasattr(self, 'chains'):
            self.initMolecules()
            self.chains = self.molecules.chains

    #simple/singleton
    #   regexp: 
    def test_ChainSet_regexp(self):
        """
         test regexp used in chains.get(regexp)
        """
        result = self.chains.get("A")
        self.assertEquals(result, self.chains[0:1] + self.chains[8:9])
        self.assertEquals(result.stringRepr,"(hsg1:/+/2plv:/+/1gyc:\\s\\A)")

    #       regexp decimal digits
    def test_ChainSet_regexp_decimal_digits(self):
        """
         test range regexp-decimal_digits used in chains.get("\d")
        """
        result = self.chains.get("\d")
        self.assertEquals(len(result), 4)
        self.assertEquals(result.stringRepr,"(hsg1:/+/2plv:/+/1gyc:\\s\\\\d)")

    #       regexp non decimal digits
    def test_ChainSet_regexp_non_decimal_digits(self):
        """
         test range regexp-non_decimal_digits used in chains.get("\D")
        """
        result = self.chains.get("\D")
        self.assertEquals(len(result), 6)
        self.assertEquals(result.stringRepr,"(hsg1:/+/2plv:/+/1gyc:\\s\\\\D)")

    #       regexp alphabetic range
    def test_ChainSet_regexp_alphabetic_range(self):
        """
         test range regexp-alphabetic_range used in chains.get("[A-Z]")
        """
        result = self.chains.get("[A-Z]")
        self.assertEquals(len(result), 5)
        self.assertEquals(result.stringRepr,"(hsg1:/+/2plv:/+/1gyc:\\s\\[A-Z])")

    #       regexp numeric range
    def test_ChainSet_regexp_numeric_range(self):
        """
         test range regexp-numeric_range used in chains.get("[0-9]")
        """
        result = self.chains.get("[0-9]")
        self.assertEquals(len(result), 4)
        self.assertEquals(result.stringRepr,"(hsg1:/+/2plv:/+/1gyc:\\s\\[0-9])")

    #       regexp compound range
    def test_ChainSet_regexp_compound_range(self):
        """
         test range regexp-compound_range used in chains.get("[0-9A-Z]")
        """
        result = self.chains.get("[0-9A-Z]")
        self.assertEquals(len(result), 9)
        self.assertEquals(result.stringRepr,"(hsg1:/+/2plv:/+/1gyc:\\s\\[0-9A-Z])")

    #       regexp or
    def test_ChainSet_regexp_or_regexp(self):
        """
         test range regexp-or-regexp used in chains.get("W|B")
        """
        result = self.chains.get("W|A")
        self.assertEquals(len(result), 3)
        self.assertEquals(result.stringRepr,"(hsg1:/+/2plv:/+/1gyc:\\s\\W|A)")


    #       regexp set  NB "(?B2)" didnot work here
    def test_ChainSet_regexp_set(self):
        """
         test range regexp-set used in chains.get("[B2]")
        """
        result = self.chains.get("[B2]")
        self.assertEquals(len(result), 2)
        self.assertEquals(result.stringRepr,"(hsg1:/+/2plv:/+/1gyc:\\s\\[B2])")


    #       regexp with space
    def test_ChainSet_regexp_with_space(self):
        """
         test range regexp-with_space used in chains.get("\s")
        """
        result = self.chains.get("\s")
        self.assertEquals(len(result), 1)
        self.assertEquals(result.stringRepr,"(hsg1:/+/2plv:/+/1gyc:\\s\\\\s)")


    #   relative: 
    def test_ChainSet_relative_index(self):
        """
         test relative index used in chains.get(relative_index)
        """
        result = self.chains.get("#1")
        self.assertEquals(result, self.molecules[0].chains[0:1]+ self.molecules[1].chains[0:1]+ self.molecules[2].chains[0:1])
        self.assertEquals(result.stringRepr,"(hsg1:/+/2plv:/+/1gyc:\\s\\#1)")


    #   index: 
    def test_ChainSet_index_which_matches_id(self):
        """
        NOTE: for chains, IDS have precedence over indices
         test index used in chains.get(index)
        """
        result = self.chains.get("1")
        #@@
        #index matches a chain id which is matched first!!!
        self.assertEquals(result, self.molecules[1].chains[0:1])
        self.assertEquals(result.stringRepr,"(hsg1:/+/2plv:/+/1gyc:\\s\\1)")

    def test_ChainSet_index(self):
        """
         test index used in chains.get(index)
        """
        result = self.chains.get("0")
        #@@
        #index matches a chain id which is matched first!!!
        self.assertEquals(result, self.chains[0:1])
        self.assertEquals(result.stringRepr,"(hsg1:/+/2plv:/+/1gyc:\\s\\0)")


    #   NamedChainSet:  
    def test_ChainSet_NamedChainSet_proteic(self):
        """
         test NamedResSet used in chains.get(NamedResSet)
        """
        result = self.chains.get("proteic")
        #result.id = ['A','B','1','2','3']
        self.assertEquals(len(result), 5)
        self.assertEquals(result.stringRepr,"(hsg1:/+/2plv:/+/1gyc:\\s\\proteic)")


    def test_ChainSet_NamedChainSet_dna(self):
        """
         test NamedResSet used in chains.get(dna)
        """
        from MolKit import Read
        mol = Read("Data/bdna_HS.pdb")
        result = mol.chains.get("dna")
        #print result.id  #['A', 'B']
        self.assertEquals(len(result), 2)
        self.assertEquals(result.stringRepr,"(bdna_HS:\\s\\dna)")


    #   NamedAtomSet: n/a
    def test_ChainSet_NamedAtomSet(self):
        """
         test NamedAtomSet used in chains.get(NamedAtomSet)
        """
        result = self.chains.get("backbone")
        self.assertEquals(result, ChainSet())


    #compound
    # !!!!NOTE!!!! RANGE includes both ends!!!
    # RANGE: REGEXP-REGEXP
    #       regexp -  regexp
    #       regexp -  regexp(fail)
    def test_ChainSet_range_regexp_regexp(self):
        """
         test range regexp-regexp used in chains.get(regexp-regexp)
        """
        result = self.chains.get("W*-*Z")
        #self.assertEquals(result, self.chains)
        self.assertEquals(result.stringRepr,"(hsg1:/+/2plv:/+/1gyc:\\s\\W*-*Z)")

    def test_ChainSet_range_regexp_regexp_w_spaces(self):
        """
         test range regexp-regexp_w_spaces used in chains.get(regexp-regexp)
        """
        result = self.chains.get("W* - *Z")
        #self.assertEquals(result, self.chains)
        self.assertEquals(result.stringRepr,"(hsg1:/+/2plv:/+/1gyc:\\s\\W* - *Z)")

    def test_ChainSet_range_regexp_regexp_fail_second(self):
        """
         test range regexp-regexp_fail_second used in chains.get(regexp-regexp)
        """
        result = self.chains.get("W*-6*")
        self.assertEquals(result, ChainSet())

    #       regexp(fail) -  regexp
    def test_ChainSet_range_regexp_fail_first_regexp(self):
        """
         test range regexp_fail_first-regexp used in chains.get(regexp-regexp)
        """
        result = self.chains.get("*6-*Z")
        self.assertEquals(result, ChainSet())

    #       regexp(fail) -  regexp(fail)
    def test_ChainSet_range_regexp_regexp_fail_both(self):
        """
         test range regexp-regexp_fail_both used in chains.get(regexp-regexp)
        """
        result = self.chains.get("*6-*6")
        self.assertEquals(result, ChainSet())

    # RANGE: REGEXP-RELATIVE
    #       regexp -  relative
    #       regexp -  relative(fail)
    def test_ChainSet_range_regexp_relative(self):
        """
         test range regexp-relative used in chains.get(regexp-relative)
        """
        result = self.chains.get("W*-#2")
        #self.assertEquals(result, self.chains)
        self.assertEquals(result.stringRepr,"(hsg1:/+/2plv:/+/1gyc:\\s\\W*-#2)")

    def test_ChainSet_range_regexp_relative_w_spaces(self):
        """
         test range regexp-relative_w_spaces used in chains.get(regexp-relative)
        """
        result = self.chains.get("W* - #2")
        #self.assertEquals(result, self.chains)
        self.assertEquals(result.stringRepr,"(hsg1:/+/2plv:/+/1gyc:\\s\\W* - #2)")

    def test_ChainSet_range_regexp_relative_fail_second(self):
        """
         test range regexp-relative_fail_second used in chains.get(regexp-relative)
        """
        result = self.chains.get("*1-#11")
        self.assertEquals(result, ChainSet())

    #       regexp(fail) -  relative
    def test_ChainSet_range_regexp_fail_first_relative(self):
        """
         test range regexp_fail_first-relative used in chains.get(regexp-relative)
        """
        result = self.chains.get("*6-#1")
        self.assertEquals(result, ChainSet())

    #       regexp(fail) -  relative(fail)
    def test_ChainSet_range_regexp_relative_fail_both(self):
        """
         test range regexp-relative_fail_both used in chains.get(regexp-relative)
        """
        result = self.chains.get("*6-#8")
        self.assertEquals(result, ChainSet())

    # RANGE: REGEXP-INDEX
    #       regexp -  index
    #       regexp -  index(fail)
    def test_ChainSet_range_regexp_index(self):
        """
         test range regexp-index used in chains.get(regexp-index)
        """
        result = self.chains.get("W*-2")
        #self.chains.id= ['W', 'A', 'B', '1', '2', '4', '3', ' ', 'A', 'Z']
        #self.assertEquals(result, self.chains[:3])
        self.assertEquals(result.stringRepr,"(hsg1:/+/2plv:/+/1gyc:\\s\\W*-2)")

    def test_ChainSet_range_regexp_index_w_spaces(self):
        """
         test range regexp-index_w_spaces used in chains.get(regexp-index)
        """
        result = self.chains.get("W* - 2")
        #self.chains.id= ['W', 'A', 'B', '1', '2', '4', '3', ' ', 'A', 'Z']
        #self.assertEquals(result, self.chains[:3])
        self.assertEquals(result.stringRepr,"(hsg1:/+/2plv:/+/1gyc:\\s\\W* - 2)")

    def test_ChainSet_range_regexp_index_fail_second(self):
        """
         test range regexp-index_fail_second used in chains.get(regexp-index)
        """
        result = self.chains.get("W*-13")
        self.assertEquals(result, ChainSet())

    #       regexp(fail) -  index
    def test_ChainSet_range_regexp_fail_first_index(self):
        """
         test range regexp_fail_first-index used in chains.get(regexp-index)
        """
        result = self.chains.get("T*-1")
        self.assertEquals(result, ChainSet())

    #       regexp(fail) -  index(fail)
    def test_ChainSet_range_regexp_index_fail_both(self):
        """
         test range regexp-index_fail_both used in chains.get(regexp-index)
        """
        result = self.chains.get("T*-18")
        self.assertEquals(result, ChainSet())

    # RANGE: REGEXP-NamedChainSet
    #       regexp -  NamedChainSet
    def test_ChainSet_range_regexp_NamedChainSet(self):
        """
         test range regexp-NamedResSet used in chains.get(regexp-NamedChainSet)
        """
        result = self.chains.get("W*-proteic")
        self.assertEquals(len(result), 4)
        self.assertEquals(result.stringRepr,"(hsg1:/+/2plv:/+/1gyc:\\s\\W*-proteic)")

    # RANGE: REGEXP-NamedResSet
    #       regexp -  NamedResSet
    def test_ChainSet_range_regexp_NamedResSet(self):
        """
         test range regexp-NamedResSet used in chains.get(regexp-NamedResSet)
        """
        result = self.chains.get("hs*-buried")
        self.assertEquals(result, ChainSet())

    #       regexp -  NamedResSet(fail)
    def test_ChainSet_range_regexp_NamedResSet_fail_second(self):
        """
         test range regexp-NamedResSet_fail_second used in chains.get(regexp-NamedResSet)
        """
        result = self.chains.get("*1-mistake")
        self.assertEquals(result, ChainSet())

    #       regexp(fail) -  NamedResSet
    def test_ChainSet_range_regexp_fail_first_NamedResSet(self):
        """
         test range regexp_fail_first-NamedResSet used in chains.get(regexp-NamedResSet)
        """
        result = self.chains.get("*6-buried")
        self.assertEquals(result, ChainSet())

    #       regexp(fail) -  NamedResSet(fail)
    def test_ChainSet_range_regexp_NamedResSet_fail_both(self):
        """
         test range regexp-NamedResSet_fail_both used in chains.get(regexp-NamedResSet)
        """
        result = self.chains.get("*6-mistake")
        self.assertEquals(result, ChainSet())


    # RANGE: REGEXP-NamedAtomSet
    #       regexp -  NamedAtomSet
    def test_ChainSet_range_regexp_NamedAtomSet(self):
        """
         test range regexp-NamedAtomSet used in chains.get(regexp-NamedAtomSet)
        """
        result = self.chains.get("hs*-backbone")
        self.assertEquals(result, ChainSet())

    #       regexp -  NamedAtomSet(fail)
    def test_ChainSet_range_regexp_NamedAtomSet_fail_second(self):
        """
         test range regexp-NamedAtomSet_fail_second used in chains.get(regexp-NamedAtomSet)
        """
        result = self.chains.get("*1-mistake")
        self.assertEquals(result, ChainSet())

    #       regexp(fail) -  NamedAtomSet
    def test_ChainSet_range_regexp_fail_first_NamedAtomSet(self):
        """
         test range regexp_fail_first-NamedAtomSet used in chains.get(regexp-NamedAtomSet)
        """
        result = self.chains.get("*6-backbone")
        self.assertEquals(result, ChainSet())

    #       regexp(fail) -  NamedAtomSet(fail)
    def test_ChainSet_range_regexp_NamedAtomSet_fail_both(self):
        """
         test range regexp-NamedAtomSet_fail_both used in chains.get(regexp-NamedAtomSet)
        """
        result = self.chains.get("*6-mistake")
        self.assertEquals(result, ChainSet())

    #
    # RANGE: RELATIVE-REGEXP
    #       relative -  regexp
    def test_ChainSet_range_relative_regexp(self):
        """
         test range relative-regexp used in chains.get(relative-regexp)
        """
        result = self.chains.get("#1-Z*")
        self.assertEquals(result, self.chains)
        self.assertEquals(result.stringRepr,"(hsg1:/+/2plv:/+/1gyc:\\s\\#1-Z*)")

    #       relative -  regexp(fail)
    def test_ChainSet_range_relative_regexp_fail_second(self):
        """
         test range relative-regexp_fail_second used in chains.get(relative-regexp)
        """
        result = self.chains.get("#1-Q*")
        self.assertEquals(result, ChainSet())

    #       relative(fail) -  regexp
    def test_ChainSet_range_relative_fail_first_regexp(self):
        """
         test range relative_fail_first-regexp used in chains.get(relative-regexp)
        """
        result = self.chains.get("#16-Z*")
        self.assertEquals(result, ChainSet())

    #       relative(fail) -  regexp(fail)
    def test_ChainSet_range_relative_regexp_fail_both(self):
        """
         test range relative-regexp_fail_both used in chains.get(relative-regexp)
        """
        result = self.chains.get("#16-8*")
        self.assertEquals(result, ChainSet())

    # RANGE: RELATIVE-RELATIVE
    #       relative -  relative
    def test_ChainSet_range_relative_relative(self):
        """
         test range relative-relative used in chains.get(relative-relative)
        """
        result = self.chains.get("#1-#2")
        #.flat #1 matches the first chain, #2 the last, hence all of them
        self.assertEquals(result, self.chains)
        self.assertEquals(result.stringRepr,"(hsg1:/+/2plv:/+/1gyc:\\s\\#1-#2)")

    def test_ChainSet_range_relative_relative_2(self):
        """
         test range relative-relative used in chains.get(relative-relative)
        """
        result = self.chains.get("#1-#3")
        #.flat #1 matches the first chain, #3 the last of second molecule
        self.assertEquals(result, self.chains[:6])
        self.assertEquals(result, self.molecules[0].chains + self.molecules[1].chains[:3])
        self.assertEquals(result.stringRepr,"(hsg1:/+/2plv:/+/1gyc:\\s\\#1-#3)")

    #       relative -  relative(fail)
    def test_ChainSet_range_relative_relative_fail_second(self):
        """
         test range relative-relative_fail_second used in chains.get(relative-relative)
        """
        result = self.chains.get("#1-#15")
        self.assertEquals(result, ChainSet())

    #       relative(fail) -  relative
    def test_ChainSet_range_relative_fail_first_relative(self):
        """
         test range relative_fail_first-relative used in chains.get(relative-relative)
        """
        result = self.chains.get("#16-#1")
        self.assertEquals(result, ChainSet())

    #       relative(fail) -  relative(fail)
    def test_ChainSet_range_relative_relative_fail_both(self):
        """
         test range relative-relative_fail_both used in chains.get(relative-relative)
        """
        result = self.chains.get("#16-#18")
        self.assertEquals(result, ChainSet())

    # RANGE: RELATIVE-INDEX
    #       relative -  index
    def test_ChainSet_range_relative_index(self):
        """
         test range relative-index used in chains.get(relative-index)
        """
        result = self.chains.get("#1-2")
        self.assertEquals(result, self.molecules[0].chains)
        self.assertEquals(result.stringRepr,"(hsg1:/+/2plv:/+/1gyc:\\s\\#1-2)")

    #       relative -  index(fail)
    def test_ChainSet_range_relative_index_fail_second(self):
        """
         test range relative-index_fail_second used in chains.get(relative-index)
        """
        result = self.chains.get("#1-13")
        self.assertEquals(result, ChainSet())


    #       relative(fail) -  index
    def test_ChainSet_range_relative_fail_first_index(self):
        """
         test range relative_fail_first-index used in chains.get(relative-index)
        """
        result = self.chains.get("#16-1")
        self.assertEquals(result, ChainSet())

    #       relative(fail) -  index(fail)
    def test_ChainSet_range_relative_index_fail_both(self):
        """
         test range relative_index-fail_both used in chains.get(relative-index)
        """
        result = self.chains.get("#16-18")
        self.assertEquals(result, ChainSet())

    # RANGE: RELATIVE-NamedResSet
    #       relative -  NamedResSet
    def test_ChainSet_range_relative_NamedResSet(self):
        """
         test range relative-NamedResSet used in chains.get(relative-NamedResSet)
        """
        result = self.chains.get("#1-buried")
        self.assertEquals(result, ChainSet())

    #       relative -  NamedResSet(fail)
    def test_ChainSet_range_relative_NamedResSet_fail_second(self):
        """
         test range relative-NamedResSet_fail_second used in chains.get(relative-NamedResSet)
        """
        result = self.chains.get("#1-mistake")
        self.assertEquals(result, ChainSet())

    #       relative(fail) -  NamedResSet
    def test_ChainSet_range_relative_fail_first_NamedResSet(self):
        """
         test range relative_fail-first_NamedResSet used in chains.get(relative-NamedResSet)
        """
        result = self.chains.get("#16-buried")
        self.assertEquals(result, ChainSet())

    #       relative(fail) -  NamedResSet(fail)
    def test_ChainSet_range_relative_NamedResSet_fail_both(self):
        """
         test range relative-NamedResSet_fail_both used in chains.get(relative-NamedResSet)
        """
        result = self.chains.get("#16-mistake")
        self.assertEquals(result, ChainSet())


    # RANGE: RELATIVE-NamedAtomSet
    #       relative -  NamedAtomSet
    def test_ChainSet_range_relative_NamedAtomSet(self):
        """
         test range relative-NamedAtomSet used in chains.get(relative-NamedAtomSet)
        """
        result = self.chains.get("#1-backbone")
        self.assertEquals(result, ChainSet())

    #       relative -  NamedAtomSet(fail)
    def test_ChainSet_range_relative_NamedAtomSet_fail_second(self):
        """
         test range relative-NamedAtomSet_fail_second used in chains.get(relative-NamedAtomSet)
        """
        result = self.chains.get("#1-mistake")
        self.assertEquals(result, ChainSet())

    #       relative(fail) -  NamedAtomSet
    def test_ChainSet_range_relative_fail_first_NamedAtomSet(self):
        """
         test range relative_fail_first-NamedAtomSet used in chains.get(relative-NamedAtomSet)
        """
        result = self.chains.get("#16-backbone")
        self.assertEquals(result, ChainSet())

    #       relative(fail) -  NamedAtomSet(fail)
    def test_ChainSet_range_relative_NamedAtomSet_fail_both(self):
        """
         test range relative-NamedAtomSet_fail_both used in chains.get(relative-NamedAtomSet)
        """
        result = self.chains.get("#16-mistake")
        self.assertEquals(result, ChainSet())

    #
    # RANGE: INDEX-REGEXP
    #       index -  regexp
    def test_ChainSet_range_index_regexp(self):
        """
         test range index-regexp used in chains.get(index-regexp)
        """
        result = self.chains.get("1-Z*")
        self.assertEquals(result, self.chains[1:])
        self.assertEquals(result.stringRepr,"(hsg1:/+/2plv:/+/1gyc:\\s\\1-Z*)")

    #       index -  regexp(fail)
    def test_ChainSet_range_index_regexp_fail_second(self):
        """
         test range index-regexp_fail_second used in chains.get(index-regexp)
        """
        result = self.chains.get("1-6*")
        self.assertEquals(result, ChainSet())

    #       index(fail) -  regexp
    def test_ChainSet_range_index_fail_first_regexp(self):
        """
         test range index_fail_first-regexp used in chains.get(index-regexp)
        """
        result = self.chains.get("6-1*")
        self.assertEquals(result, ChainSet())

    #       index(fail) -  regexp(fail)
    def test_ChainSet_range_index_regexp_fail_both(self):
        """
         test range index-regexp_fail_both used in chains.get(index-regexp)
        """
        result = self.chains.get("6-8*")
        self.assertEquals(result, ChainSet())

    # RANGE: INDEX-RELATIVE
    #       index -  relative
    def test_ChainSet_range_index_relative(self):
        """
         test range index-relative used in chains.get(index-relative)
        """
        result = self.chains.get("1-#2")
        self.assertEquals(result, self.chains[1:])
        self.assertEquals(result.stringRepr,"(hsg1:/+/2plv:/+/1gyc:\\s\\1-#2)")

    #       index -  relative(fail)
    def test_ChainSet_range_index_relative_fail_second(self):
        """
         test range index-relative_fail_second used in chains.get(index-relative)
        """
        result = self.chains.get("1-#15")
        self.assertEquals(result, ChainSet())

    #       index(fail) -  relative
    def test_ChainSet_range_index_fail_first_relative(self):
        """
         test range index_fail_first-relative used in chains.get(index-relative)
        """
        result = self.chains.get("16-#1")
        self.assertEquals(result, ChainSet())

    #       index(fail) -  relative(fail)
    def test_ChainSet_range_index_relative_fail_both(self):
        """
         test range index-relative_fail_both used in chains.get(index-relative)
        """
        result = self.chains.get("16-#18")
        self.assertEquals(result, ChainSet())

    # RANGE: INDEX-INDEX
    #       index -  index
    def test_ChainSet_range_index_index(self):
        """
         test range index-index used in chains.get(index-index)
        """
        result = self.chains.get("1-2")
        self.assertEquals(result, self.chains[1:3])
        self.assertEquals(result.stringRepr,"(hsg1:/+/2plv:/+/1gyc:\\s\\1-2)")

    #       index -  index(fail)
    def test_ChainSet_range_index_index_fail_second(self):
        """
         test range index-index_fail_second used in chains.get(index-index)
        """
        result = self.chains.get("1-13")
        self.assertEquals(result, ChainSet())

    #       index(fail) -  index
    def test_ChainSet_range_index_fail_first_index(self):
        """
         test range index_fail_first-index used in chains.get(index-index)
        """
        result = self.chains.get("16-1")
        self.assertEquals(result, ChainSet())

    #       index(fail) -  index(fail)
    def test_ChainSet_range_index_index_fail_both(self):
        """
         test range index_index-fail_both used in chains.get(index-index)
        """
        result = self.chains.get("16-18")
        self.assertEquals(result, ChainSet())

    # RANGE: INDEX-NamedResSet
    #       index -  NamedResSet
    def test_ChainSet_range_index_NamedResSet(self):
        """
         test range index-NamedResSet used in chains.get(index-NamedResSet)
        """
        result = self.chains.get("1-buried")
        self.assertEquals(result, ChainSet())

    #       index -  NamedResSet(fail)
    def test_ChainSet_range_index_NamedResSet_fail_second(self):
        """
         test range index-NamedResSet_fail_second used in chains.get(index-NamedResSet)
        """
        result = self.chains.get("1-mistake")
        self.assertEquals(result, ChainSet())

    #       index(fail) -  NamedResSet
    def test_ChainSet_range_index_fail_first_NamedResSet(self):
        """
         test range index_fail-first_NamedResSet used in chains.get(index-NamedResSet)
        """
        result = self.chains.get("6-buried")
        self.assertEquals(result, ChainSet())

    #       index(fail) -  NamedResSet(fail)
    def test_ChainSet_range_index_NamedResSet_fail_both(self):
        """
         test range index-NamedResSet_fail_both used in chains.get(index-NamedResSet)
        """
        result = self.chains.get("6-mistake")
        self.assertEquals(result, ChainSet())


    # RANGE: INDEX-NamedAtomSet
    #       index -  NamedAtomSet
    def test_ChainSet_range_index_NamedAtomSet(self):
        """
         test range index-NamedAtomSet used in chains.get(index-NamedAtomSet)
        """
        result = self.chains.get("1-backbone")
        self.assertEquals(result, ChainSet())

    #       index -  NamedAtomSet(fail)
    def test_ChainSet_range_index_NamedAtomSet_fail_second(self):
        """
         test range index-NamedAtomSet_fail_second used in chains.get(index-NamedAtomSet)
        """
        result = self.chains.get("1-mistake")
        self.assertEquals(result, ChainSet())

    #       index(fail) -  NamedAtomSet
    def test_ChainSet_range_index_fail_first_NamedAtomSet(self):
        """
         test range index_fail_first-NamedAtomSet used in chains.get(index-NamedAtomSet)
        """
        result = self.chains.get("6-backbone")
        self.assertEquals(result, ChainSet())

    #       index(fail) -  NamedAtomSet(fail)
    def test_ChainSet_range_index_NamedAtomSet_fail_both(self):
        """
         test range index-NamedAtomSet_fail_both used in chains.get(index-NamedAtomSet)
        """
        result = self.chains.get("6-mistake")
        self.assertEquals(result, ChainSet())


class ResidueSetSelectionTests(SelectionBaseTest):
    #singleton tests: regexp, relative, index, NamedResSet, NamedAtomSet

    def setUp(self):
        if not hasattr(self, 'residues'):
            self.initMolecules()
            self.residues = self.molecules.chains.residues

    #simple/singleton
    #   regexp: 
    def test_ResidueSet_regexp(self):
        """
         test regexp used in residues.get(regexp)
        """
        result = self.residues.get("ALA*")
        self.assertEquals(len(result), 121)
        self.assertEquals(result.stringRepr,"(hsg1::/+/2plv::/+/1gyc::\\s\\ALA*)")


    #   relative: 
    def test_ResidueSet_relative_index(self):
        """
         test relative index used in residues.get(relative_index)
        """
        result = self.residues.get("#1")
        self.assertEquals(len(result),len(self.molecules.chains))
        self.assertEquals(result.stringRepr,"(hsg1::/+/2plv::/+/1gyc::\\s\\#1)")


    #   index: 
    def test_ResidueSet_index_1(self):
        """
         test index used in residues.get(index)
        """
        result = self.residues.get("1")
        self.assertEquals(result, self.molecules.chains.residues[1:2])
        self.assertEquals(result.stringRepr,"(hsg1::/+/2plv::/+/1gyc::\\s\\1)")

    def test_ResidueSet_index_0(self):
        """
         test index used in residues.get(index)
        """
        result = self.residues.get("0")
        self.assertEquals(result, self.residues[0:1])
        self.assertEquals(result.stringRepr,"(hsg1::/+/2plv::/+/1gyc::\\s\\0)")


    #   NamedResSet:  
    def test_ResidueSet_NamedResSet(self):
        """
         test NamedResSet used in residues.get(NamedResSet)
        """
        result = self.residues.get("acidic")
        self.assertEquals(result.stringRepr,"(hsg1::/+/2plv::/+/1gyc::\\s\\acidic)")
        self.assertEquals(len(result), 132)


    #   NamedAtomSet: n/a
    def test_ResidueSet_NamedAtomSet(self):
        """
         test NamedAtomSet used in residues.get(NamedAtomSet)
        """
        result = self.residues.get("backbone")
        self.assertEquals(result, ResidueSet())


    #compound
    # !!!!NOTE!!!! RANGE includes both ends!!!
    # RANGE: REGEXP-REGEXP
    #       regexp -  regexp
    def test_ResidueSet_range_regexp_regexp(self):
        """
         test range regexp-regexp used in residues.get(regexp-regexp)
        """
        result = self.residues.get("ALA*-TRP*")
        self.assertEquals(result.stringRepr,"(hsg1::/+/2plv::/+/1gyc::\\s\\ALA*-TRP*)")
        #with correction which limited relative ranges to #num-#num only
        #self.assertEquals(len(result), 1153)
        #self.assertEquals(len(result), 1996)

    #       regexp -  regexp(fail)
    def test_ResidueSet_range_regexp_regexp_fail_second(self):
        """
         test range regexp-regexp_fail_second used in residues.get(regexp-regexp)
        """
        result = self.residues.get("ALA*-FOO*")
        self.assertEquals(result, ResidueSet())

    #       regexp(fail) -  regexp
    def test_ResidueSet_range_regexp_fail_first_regexp(self):
        """
         test range regexp_fail_first-regexp used in residues.get(regexp-regexp)
        """
        result = self.residues.get("FOO*-TRP*")
        self.assertEquals(result, ResidueSet())

    #       regexp(fail) -  regexp(fail)
    def test_ResidueSet_range_regexp_regexp_fail_both(self):
        """
         test range regexp-regexp_fail_both used in residues.get(regexp-regexp)
        """
        result = self.residues.get("FOO*-*BAR")
        self.assertEquals(result, ResidueSet())

    # RANGE: REGEXP-RELATIVE
    #       regexp -  relative
    def test_ResidueSet_range_regexp_relative(self):
        """
         test range regexp-relative used in residues.get(regexp-relative)
        """
        result = self.residues.get("PRO*-#2")
        self.assertEquals(result.stringRepr,"(hsg1::/+/2plv::/+/1gyc::\\s\\PRO*-#2)")
        #with correction which limited relative ranges to #num-#num only
        #self.assertEquals(len(result),  4)
        #self.assertEquals(len(result),  2049)

    #       regexp -  relative(fail)
    def test_ResidueSet_range_regexp_relative_fail_second(self):
        """
         test range regexp-relative_fail_second used in residues.get(regexp-relative)
        """
        result = self.residues.get("ALA*-#11111")
        self.assertEquals(result, ResidueSet())

    #       regexp(fail) -  relative
    def test_ResidueSet_range_regexp_fail_first_relative(self):
        """
         test range regexp_fail_first-relative used in residues.get(regexp-relative)
        """
        result = self.residues.get("FOO*-#1")
        self.assertEquals(result, ResidueSet())

    #       regexp(fail) -  relative(fail)
    def test_ResidueSet_range_regexp_relative_fail_both(self):
        """
         test range regexp-relative_fail_both used in residues.get(regexp-relative)
        """
        result = self.residues.get("FOO*-#88888")
        self.assertEquals(result, ResidueSet())

    # RANGE: REGEXP-INDEX
    #       regexp -  index
    def test_ResidueSet_range_regexp_index(self):
        """
         test range regexp-index used in residues.get(regexp-index)
        """
        result = self.residues.get("PRO*-2")
        self.assertEquals(result.stringRepr,"(hsg1::/+/2plv::/+/1gyc::\\s\\PRO*-2)")
        #with correction which limited relative ranges to #num-#num only
        #self.assertEquals(len(result), 7)
        #self.assertEquals(len(result), 0)

    #       regexp -  index(fail)
    def test_ResidueSet_range_regexp_index_fail_second(self):
        """
         test range regexp-index_fail_second used in residues.get(regexp-index)
        """
        result = self.residues.get("PRO*-13333")
        self.assertEquals(result, ResidueSet())

    #       regexp(fail) -  index
    def test_ResidueSet_range_regexp_fail_first_index(self):
        """
         test range regexp_fail_first-index used in residues.get(regexp-index)
        """
        result = self.residues.get("FOO*-1")
        self.assertEquals(result, ResidueSet())

    #       regexp(fail) -  index(fail)
    def test_ResidueSet_range_regexp_index_fail_both(self):
        """
         test range regexp-index_fail_both used in residues.get(regexp-index)
        """
        result = self.residues.get("FOO*-18")
        self.assertEquals(result, ResidueSet())

    # RANGE: REGEXP-NamedResSet
    #       regexp -  NamedResSet
    def test_ResidueSet_range_regexp_NamedResSet(self):
        """
         test range regexp-NamedResSet used in residues.get(regexp-NamedResSet)
        """
        result = self.residues.get("PRO*-buried")
        #with correction which limited relative ranges to #num-#num only
        #self.assertEquals(len(result), 1449)
        #self.assertEquals(len(result), 2030)
        self.assertEquals(result.stringRepr,"(hsg1::/+/2plv::/+/1gyc::\\s\\PRO*-buried)")

    #       regexp -  NamedResSet(fail)
    def test_ResidueSet_range_regexp_NamedResSet_fail_second(self):
        """
         test range regexp-NamedResSet_fail_second used in residues.get(regexp-NamedResSet)
        """
        result = self.residues.get("PRO*-mistake")
        self.assertEquals(result, ResidueSet())

    #       regexp(fail) -  NamedResSet
    def test_ResidueSet_range_regexp_fail_first_NamedResSet(self):
        """
         test range regexp_fail_first-NamedResSet used in residues.get(regexp-NamedResSet)
        """
        result = self.residues.get("FOO*-buried")
        self.assertEquals(result, ResidueSet())

    #       regexp(fail) -  NamedResSet(fail)
    def test_ResidueSet_range_regexp_NamedResSet_fail_both(self):
        """
         test range regexp-NamedResSet_fail_both used in residues.get(regexp-NamedResSet)
        """
        result = self.residues.get("FOO*-mistake")
        self.assertEquals(result, ResidueSet())


    # RANGE: REGEXP-NamedAtomSet
    #       regexp -  NamedAtomSet
    def test_ResidueSet_range_regexp_NamedAtomSet(self):
        """
         test range regexp-NamedAtomSet used in residues.get(regexp-NamedAtomSet)
        """
        result = self.residues.get("PRO*-backbone")
        self.assertEquals(result, ResidueSet())

    #       regexp -  NamedAtomSet(fail)
    def test_ResidueSet_range_regexp_NamedAtomSet_fail_second(self):
        """
         test range regexp-NamedAtomSet_fail_second used in residues.get(regexp-NamedAtomSet)
        """
        result = self.residues.get("PRO*-mistake")
        self.assertEquals(result, ResidueSet())

    #       regexp(fail) -  NamedAtomSet
    def test_ResidueSet_range_regexp_fail_first_NamedAtomSet(self):
        """
         test range regexp_fail_first-NamedAtomSet used in residues.get(regexp-NamedAtomSet)
        """
        result = self.residues.get("FOO*-backbone")
        self.assertEquals(result, ResidueSet())

    #       regexp(fail) -  NamedAtomSet(fail)
    def test_ResidueSet_range_regexp_NamedAtomSet_fail_both(self):
        """
         test range regexp-NamedAtomSet_fail_both used in residues.get(regexp-NamedAtomSet)
        """
        result = self.residues.get("FOO*-atomSetmistake")
        self.assertEquals(result, ResidueSet())

    #
    # RANGE: RELATIVE-REGEXP
    #       relative -  regexp
    def test_ResidueSet_range_relative_regexp(self):
        """
         test range relative-regexp used in residues.get(relative-regexp)
        """
        result = self.residues.get("#1-PRO*")
        self.assertEquals(result.stringRepr,"(hsg1::/+/2plv::/+/1gyc::\\s\\#1-PRO*)")
        #with correction which limited relative ranges to #num-#num only
        #self.assertEquals(len(result), 1461)
        self.assertEquals(len(result), 2149)

    #       relative -  regexp(fail)
    def test_ResidueSet_range_relative_regexp_fail_second(self):
        """
         test range relative-regexp_fail_second used in residues.get(relative-regexp)
        """
        result = self.residues.get("#1-Q*")
        self.assertEquals(result, ResidueSet())

    #       relative(fail) -  regexp
    def test_ResidueSet_range_relative_fail_first_regexp(self):
        """
         test range relative_fail_first-regexp used in residues.get(relative-regexp)
        """
        result = self.residues.get("#16666-PRO*")
        self.assertEquals(result, ResidueSet())

    #       relative(fail) -  regexp(fail)
    def test_ResidueSet_range_relative_regexp_fail_both(self):
        """
         test range relative-regexp_fail_both used in residues.get(relative-regexp)
        """
        result = self.residues.get("#166666-FOO*")
        self.assertEquals(result, ResidueSet())

    # RANGE: RELATIVE-RELATIVE
    #       relative -  relative
    def test_ResidueSet_range_relative_relative(self):
        """
         test range relative-relative used in residues.get(relative-relative)
        """
        result = self.residues.get("#1-#2")
        #.flat #1 matches the first reisude, #2 the second
        self.assertEquals(len(result), 20)
        self.assertEquals(result.stringRepr,"(hsg1::/+/2plv::/+/1gyc::\\s\\#1-#2)")

    def test_ResidueSet_range_relative_relative_2(self):
        """
         test range relative-relative used in residues.get(relative-relative)
        """
        result = self.residues.get("#1-#3")
        self.assertEquals(len(result), 30)
        self.assertEquals(result.stringRepr,"(hsg1::/+/2plv::/+/1gyc::\\s\\#1-#3)")

    #       relative -  relative(fail)
    def test_ResidueSet_range_relative_relative_fail_second(self):
        """
         test range relative-relative_fail_second used in residues.get(relative-relative)
        """
        result = self.residues.get("#1-#15555")
        self.assertEquals(result, ResidueSet())

    #       relative(fail) -  relative
    def test_ResidueSet_range_relative_fail_first_relative(self):
        """
         test range relative_fail_first-relative used in residues.get(relative-relative)
        """
        result = self.residues.get("#166666-#1")
        self.assertEquals(result, ResidueSet())

    #       relative(fail) -  relative(fail)
    def test_ResidueSet_range_relative_relative_fail_both(self):
        """
         test range relative-relative_fail_both used in residues.get(relative-relative)
        """
        result = self.residues.get("#166666-#1888888")
        self.assertEquals(result, ResidueSet())

    # RANGE: RELATIVE-INDEX
    #       relative -  index
    def test_ResidueSet_range_relative_index(self):
        """
         test range relative-index used in residues.get(relative-index)
        """
        result = self.residues.get("#1-2") #matches 0,1,2
        self.assertEquals(result.stringRepr,"(hsg1::/+/2plv::/+/1gyc::\\s\\#1-2)")
        #with correction which limited relative ranges to #num-#num only
        #self.assertEquals(len(result), 30)
        self.assertEquals(len(result), 3)

    #       relative -  index(fail)
    def test_ResidueSet_range_relative_index_fail_second(self):
        """
         test range relative-index_fail_second used in residues.get(relative-index)
        """
        result = self.residues.get("#1-133333")
        self.assertEquals(result, ResidueSet())


    #       relative(fail) -  index
    def test_ResidueSet_range_relative_fail_first_index(self):
        """
         test range relative_fail_first-index used in residues.get(relative-index)
        """
        result = self.residues.get("#166666-1")
        self.assertEquals(result, ResidueSet())

    #       relative(fail) -  index(fail)
    def test_ResidueSet_range_relative_index_fail_both(self):
        """
         test range relative_index-fail_both used in residues.get(relative-index)
        """
        result = self.residues.get("#166666-1888888")
        self.assertEquals(result, ResidueSet())

    # RANGE: RELATIVE-NamedResSet
    #       relative -  NamedResSet
    def test_ResidueSet_range_relative_NamedResSet(self):
        """
         test range relative-NamedResSet used in residues.get(relative-NamedResSet)
        """
        result = self.residues.get("#1-buried")
        self.assertEquals(result.stringRepr,"(hsg1::/+/2plv::/+/1gyc::\\s\\#1-buried)")
        #with correction which limited relative ranges to #num-#num only
        #self.assertEquals(len(result), 1545)
        self.assertEquals(len(result), 2157)

    #       relative -  NamedResSet(fail)
    def test_ResidueSet_range_relative_NamedResSet_fail_second(self):
        """
         test range relative-NamedResSet_fail_second used in residues.get(relative-NamedResSet)
        """
        result = self.residues.get("#1-mistake")
        self.assertEquals(result, ResidueSet())

    #       relative(fail) -  NamedResSet
    def test_ResidueSet_range_relative_fail_first_NamedResSet(self):
        """
         test range relative_fail-first_NamedResSet used in residues.get(relative-NamedResSet)
        """
        result = self.residues.get("#166666-buried")
        self.assertEquals(result, ResidueSet())

    #       relative(fail) -  NamedResSet(fail)
    def test_ResidueSet_range_relative_NamedResSet_fail_both(self):
        """
         test range relative-NamedResSet_fail_both used in residues.get(relative-NamedResSet)
        """
        result = self.residues.get("#166666-mistake")
        self.assertEquals(result, ResidueSet())


    # RANGE: RELATIVE-NamedAtomSet
    #       relative -  NamedAtomSet
    def test_ResidueSet_range_relative_NamedAtomSet(self):
        """
         test range relative-NamedAtomSet used in residues.get(relative-NamedAtomSet)
        """
        result = self.residues.get("#1-backbone")
        self.assertEquals(result, ResidueSet())

    #       relative -  NamedAtomSet(fail)
    def test_ResidueSet_range_relative_NamedAtomSet_fail_second(self):
        """
         test range relative-NamedAtomSet_fail_second used in residues.get(relative-NamedAtomSet)
        """
        result = self.residues.get("#1-mistake")
        self.assertEquals(result, ResidueSet())

    #       relative(fail) -  NamedAtomSet
    def test_ResidueSet_range_relative_fail_first_NamedAtomSet(self):
        """
         test range relative_fail_first-NamedAtomSet used in residues.get(relative-NamedAtomSet)
        """
        result = self.residues.get("#1666666-backbone")
        self.assertEquals(result, ResidueSet())

    #       relative(fail) -  NamedAtomSet(fail)
    def test_ResidueSet_range_relative_NamedAtomSet_fail_both(self):
        """
         test range relative-NamedAtomSet_fail_both used in residues.get(relative-NamedAtomSet)
        """
        result = self.residues.get("#1666666-mistake")
        self.assertEquals(result, ResidueSet())

    #
    # RANGE: INDEX-REGEXP
    #       index -  regexp
    def test_ResidueSet_range_index_regexp(self):
        """
         test range index-regexp used in residues.get(index-regexp)
        """
        result = self.residues.get("1-PRO*")
        self.assertEquals(result.stringRepr,"(hsg1::/+/2plv::/+/1gyc::\\s\\1-PRO*)")
        #with correction which limited relative ranges to #num-#num only
        #self.assertEquals(len(result), 1454)
        self.assertEquals(len(result), 2148)

    #       index -  regexp(fail)
    def test_ResidueSet_range_index_regexp_fail_second(self):
        """
         test range index-regexp_fail_second used in residues.get(index-regexp)
        """
        result = self.residues.get("1-666666*")
        self.assertEquals(result, ResidueSet())

    #       index(fail) -  regexp
    def test_ResidueSet_range_index_fail_first_regexp(self):
        """
         test range index_fail_first-regexp used in residues.get(index-regexp)
        """
        result = self.residues.get("66666-1*")
        self.assertEquals(result, ResidueSet())

    #       index(fail) -  regexp(fail)
    def test_ResidueSet_range_index_regexp_fail_both(self):
        """
         test range index-regexp_fail_both used in residues.get(index-regexp)
        """
        result = self.residues.get("6666-FOO*")
        self.assertEquals(result, ResidueSet())

    # RANGE: INDEX-RELATIVE
    #       index -  relative
    def test_ResidueSet_range_index_relative(self):
        """
         test range index-relative used in residues.get(index-relative)
        """
        result = self.residues.get("1-#2")
        self.assertEquals(result.stringRepr,"(hsg1::/+/2plv::/+/1gyc::\\s\\1-#2)")
        #with correction which limited relative ranges to #num-#num only
        #self.assertEquals(len(result), 10)
        self.assertEquals(len(result), 2175)

    #       index -  relative(fail)
    def test_ResidueSet_range_index_relative_fail_second(self):
        """
         test range index-relative_fail_second used in residues.get(index-relative)
        """
        result = self.residues.get("1-#15555")
        self.assertEquals(result, ResidueSet())

    #       index(fail) -  relative
    def test_ResidueSet_range_index_fail_first_relative(self):
        """
         test range index_fail_first-relative used in residues.get(index-relative)
        """
        result = self.residues.get("16-#1")
        #with correction which limited relative ranges to #num-#num only
        #self.assertEquals(result, ResidueSet())
        self.assertEquals(len(result), 2159)

    #       index(fail) -  relative(fail)
    def test_ResidueSet_range_index_relative_fail_both(self):
        """
         test range index-relative_fail_both used in residues.get(index-relative)
        """
        result = self.residues.get("16666-#188888")
        self.assertEquals(result, ResidueSet())

    # RANGE: INDEX-INDEX
    #       index -  index
    def test_ResidueSet_range_index_index(self):
        """
         test range index-index used in residues.get(index-index)
        """
        result = self.residues.get("1-2")
        self.assertEquals(result.stringRepr,"(hsg1::/+/2plv::/+/1gyc::\\s\\1-2)")
        #with correction which limited relative ranges to #num-#num only
        #self.assertEquals(len(result), len(self.molecules.chains)*2)
        self.assertEquals(len(result), 2)

    #       index -  index(fail)
    def test_ResidueSet_range_index_index_fail_second(self):
        """
         test range index-index_fail_second used in residues.get(index-index)
        """
        result = self.residues.get("1-133333")
        self.assertEquals(result, ResidueSet())

    #       index(fail) -  index
    def test_ResidueSet_range_index_fail_first_index(self):
        """
         test range index_fail_first-index used in residues.get(index-index)
        """
        result = self.residues.get("166666-1")
        self.assertEquals(result, ResidueSet())

    #       index(fail) -  index(fail)
    def test_ResidueSet_range_index_index_fail_both(self):
        """
         test range index_index-fail_both used in residues.get(index-index)
        """
        result = self.residues.get("166666-188888")
        self.assertEquals(result, ResidueSet())

    # RANGE: INDEX-NamedResSet
    #       index -  NamedResSet
    def test_ResidueSet_range_index_NamedResSet(self):
        """
         test range index-NamedResSet used in residues.get(index-NamedResSet)
        """
        result = self.residues.get("1-buried")
        self.assertEquals(result.stringRepr,"(hsg1::/+/2plv::/+/1gyc::\\s\\1-buried)")
        #with correction which limited relative ranges to #num-#num only
        #self.assertEquals(len(result),  1538)
        self.assertEquals(len(result),  2156)

    #       index -  NamedResSet(fail)
    def test_ResidueSet_range_index_NamedResSet_fail_second(self):
        """
         test range index-NamedResSet_fail_second used in residues.get(index-NamedResSet)
        """
        result = self.residues.get("1-mistake")
        self.assertEquals(result, ResidueSet())

    #       index(fail) -  NamedResSet
    def test_ResidueSet_range_index_fail_first_NamedResSet(self):
        """
         test range index_fail-first_NamedResSet used in residues.get(index-NamedResSet)
        """
        result = self.residues.get("666666-buried")
        self.assertEquals(result, ResidueSet())

    #       index(fail) -  NamedResSet(fail)
    def test_ResidueSet_range_index_NamedResSet_fail_both(self):
        """
         test range index-NamedResSet_fail_both used in residues.get(index-NamedResSet)
        """
        result = self.residues.get("6-mistake")
        self.assertEquals(result, ResidueSet())


    # RANGE: INDEX-NamedAtomSet
    #       index -  NamedAtomSet
    def test_ResidueSet_range_index_NamedAtomSet(self):
        """
         test range index-NamedAtomSet used in residues.get(index-NamedAtomSet)
        """
        result = self.residues.get("1-backbone")
        self.assertEquals(result, ResidueSet())

    #       index -  NamedAtomSet(fail)
    def test_ResidueSet_range_index_NamedAtomSet_fail_second(self):
        """
         test range index-NamedAtomSet_fail_second used in residues.get(index-NamedAtomSet)
        """
        result = self.residues.get("1-mistake")
        self.assertEquals(result, ResidueSet())

    #       index(fail) -  NamedAtomSet
    def test_ResidueSet_range_index_fail_first_NamedAtomSet(self):
        """
         test range index_fail_first-NamedAtomSet used in residues.get(index-NamedAtomSet)
        """
        result = self.residues.get("666666-backbone")
        self.assertEquals(result, ResidueSet())

    #       index(fail) -  NamedAtomSet(fail)
    def test_ResidueSet_range_index_NamedAtomSet_fail_both(self):
        """
         test range index-NamedAtomSet_fail_both used in residues.get(index-NamedAtomSet)
        """
        result = self.residues.get("666666-mistake")
        self.assertEquals(result, ResidueSet())


class AtomSetSelectionTests(SelectionBaseTest):
    #singleton tests: regexp, relative, index, NamedResSet, NamedAtomSet

    def setUp(self):
        if not hasattr(self, 'atoms'):
            self.initMolecules()
            self.atoms = self.molecules.chains.residues.atoms

    #simple/singleton
    #   regexp: 

    def test_AtomSet_regexp(self):
        """
         test regexp used in atoms.get(regexp)
        """
        result = self.atoms.get("C*")
        self.assertEquals(result.stringRepr,"(hsg1:::/+/2plv:::/+/1gyc:::\\s\\C*)")
        self.assertEquals(len(result), 7781)


    def test_AtomSet_regexp_sanity_no_wild_cards(self):
        """
         sanity check test regexp no wild cards
        """
        result = self.atoms.get("C")
        one_letter = self.atoms.get(lambda x:x.element[0]=='C' and len(x.name)==1)
        self.assertEquals(result, one_letter)


    def test_AtomSet_regexp_sanity_one_dot_wild_card(self):
        """
         sanity check test regexp one dot wild card
        """
        result = self.atoms.get("C.")
        two_letters = self.atoms.get(lambda x:x.element[0]=='C' and len(x.name)==2)
        self.assertEquals(result, two_letters)


    def test_AtomSet_regexp_sanity_combo_wild_cards(self):
        """
         sanity check test regexp combo wild cards
        """
        result = self.atoms.get("C*")
        two_or_more_result = self.atoms.get("C.*")
        one_result = self.atoms.get("C")
        self.assertEquals(len(result), len(two_or_more_result + one_result))


    #   relative: 
    def test_AtomSet_relative_index(self):
        """
         test relative index used in atoms.get(relative_index)
        """
        result = self.atoms.get("#1")
        self.assertEquals(len(result),len(self.molecules.chains.residues))
        self.assertEquals(result.stringRepr,"(hsg1:::/+/2plv:::/+/1gyc:::\\s\\#1)")


    #   index: 
    def test_AtomSet_index_1(self):
        """
         test index used in atoms.get(index)
        """
        result = self.atoms.get("1")
        self.assertEquals(result, self.molecules.chains.residues.atoms[1:2])
        self.assertEquals(result.stringRepr,"(hsg1:::/+/2plv:::/+/1gyc:::\\s\\1)")

    def test_AtomSet_index_0(self):
        """
         test index used in atoms.get(index)
        """
        result = self.atoms.get("0")
        self.assertEquals(result, self.atoms[0:1])
        self.assertEquals(result.stringRepr,"(hsg1:::/+/2plv:::/+/1gyc:::\\s\\0)")


    #   NamedResSet:  n/a
    def test_AtomSet_NamedResSet(self): 
        """
         test NamedResSet used in atoms.get(NamedResSet)
        """
        result = self.atoms.get("acidic")
        self.assertEquals(result, AtomSet())


    #   NamedAtomSet: 
    def test_AtomSet_NamedAtomSet(self):
        """
        test NamedAtomSet used in atoms.get(NamedAtomSet)
        """
        result = self.atoms.get("backbone")
        self.assertEquals(result.stringRepr,"(hsg1:::/+/2plv:::/+/1gyc:::\\s\\backbone)")
        self.assertEquals(len(result), 6165)


    #compound
    # !!!!NOTE!!!! RANGE includes both ends!!!
    # RANGE: REGEXP-REGEXP
    #       regexp -  regexp
    def test_AtomSet_range_regexp_regexp(self):
        """
         test range regexp-regexp used in atoms.get(regexp-regexp)
        """
        result = self.atoms.get("C*-N*")
        self.assertEquals(result.stringRepr,"(hsg1:::/+/2plv:::/+/1gyc:::\\s\\C*-N*)")
        self.assertEquals(len(result), 12799)

    #       regexp -  regexp(fail)
    def test_AtomSet_range_regexp_regexp_fail_second(self):
        """
         test range regexp-regexp_fail_second used in atoms.get(regexp-regexp)
        """
        result = self.atoms.get("C*-FOO*")
        self.assertEquals(result, AtomSet())

    #       regexp(fail) -  regexp
    def test_AtomSet_range_regexp_fail_first_regexp(self):
        """
         test range regexp_fail_first-regexp used in atoms.get(regexp-regexp)
        """
        result = self.atoms.get("FOO*-C*")
        self.assertEquals(result, AtomSet())

    #       regexp(fail) -  regexp(fail)
    def test_AtomSet_range_regexp_regexp_fail_both(self):
        """
         test range regexp-regexp_fail_both used in atoms.get(regexp-regexp)
        """
        result = self.atoms.get("FOO*-*BAR")
        self.assertEquals(result, AtomSet())

    # RANGE: REGEXP-RELATIVE
    #       regexp -  relative
    def test_AtomSet_range_regexp_relative(self):
        """
         test range regexp-relative used in atoms.get(regexp-relative)
        """
        result = self.molecules[0].chains[1].residues[0:2].atoms.get("C*-#5")
        #with correction which limited relative ranges to #num-#num only
        self.assertEquals(len(result),  11)
        self.assertEquals(result.stringRepr,"(hsg1:B:0-1:\s\\C*-#5)")

    #       regexp -  relative(fail)
    def test_AtomSet_range_regexp_relative_fail_second(self):
        """
         test range regexp-relative_fail_second used in atoms.get(regexp-relative)
        """
        result = self.molecules[0].chains[0].residues[0].atoms.get("C*-#11111")
        self.assertEquals(result, AtomSet())

    #       regexp(fail) -  relative
    def test_AtomSet_range_regexp_fail_first_relative(self):
        """
         test range regexp_fail_first-relative used in atoms.get(regexp-relative)
        """
        result = self.molecules[0].chains[0].residues[0].atoms.get("FOO*-#1")
        self.assertEquals(result, AtomSet())

    #       regexp(fail) -  relative(fail)
    def test_AtomSet_range_regexp_relative_fail_both(self):
        """
         test range regexp-relative_fail_both used in atoms.get(regexp-relative)
        """
        result = self.molecules[0].chains[0].residues[0].atoms.get("FOO*-#88888")
        self.assertEquals(result, AtomSet())

    # RANGE: REGEXP-INDEX
    #       regexp -  index
    def test_AtomSet_range_regexp_index(self):
        """
         test range regexp-index used in atoms.get(regexp-index)
        """
        result = self.molecules[0].chains[1].residues[0:2].atoms.get("C*-7")
        self.assertEquals(result.stringRepr,"(hsg1:B:0-1:\\s\\C*-7)")
        self.assertEquals(len(result), 7)

    #       regexp -  index(fail)
    def test_AtomSet_range_regexp_index_fail_second(self):
        """
         test range regexp-index_fail_second used in atoms.get(regexp-index)
        """
        result = self.molecules[0].chains[2].residues[0:2].atoms.get("C*-7777")
        self.assertEquals(result, AtomSet())

    #       regexp(fail) -  index
    def test_AtomSet_range_regexp_fail_first_index(self):
        """
         test range regexp_fail_first-index used in atoms.get(regexp-index)
        """
        result = self.molecules[0].chains[2].residues[0:2].atoms.get("FOO*-7")
        self.assertEquals(result, AtomSet())

    #       regexp(fail) -  index(fail)
    def test_AtomSet_range_regexp_index_fail_both(self):
        """
         test range regexp-index_fail_both used in atoms.get(regexp-index)
        """
        result = self.molecules[0].chains[2].residues[0:2].atoms.get("FOO*-778")
        self.assertEquals(result, AtomSet())

    # RANGE: REGEXP-NamedResSet
    #       regexp -  NamedResSet
    def test_AtomSet_range_regexp_NamedResSet(self):
        """
         test range regexp-NamedResSet used in atoms.get(regexp-NamedResSet)
        """
        result = self.atoms.get("PRO*-buried")
        self.assertEquals(result, AtomSet())

    #       regexp -  NamedResSet(fail)
    def test_AtomSet_range_regexp_NamedResSet_fail_second(self):
        """
         test range regexp-NamedResSet_fail_second used in atoms.get(regexp-NamedResSet)
        """
        result = self.atoms.get("PRO*-mistake")
        self.assertEquals(result, AtomSet())

    #       regexp(fail) -  NamedResSet
    def test_AtomSet_range_regexp_fail_first_NamedResSet(self):
        """
         test range regexp_fail_first-NamedResSet used in atoms.get(regexp-NamedResSet)
        """
        result = self.atoms.get("FOO*-buried")
        self.assertEquals(result, AtomSet())

    #       regexp(fail) -  NamedResSet(fail)
    def test_AtomSet_range_regexp_NamedResSet_fail_both(self):
        """
         test range regexp-NamedResSet_fail_both used in atoms.get(regexp-NamedResSet)
        """
        result = self.atoms.get("FOO*-mistake")
        self.assertEquals(result, AtomSet())


    # RANGE: REGEXP-NamedAtomSet
    #       regexp -  NamedAtomSet
    def test_AtomSet_range_regexp_NamedAtomSet(self):
        """
         test range regexp-NamedAtomSet used in atoms.get(regexp-NamedAtomSet)
        """
        result = self.molecules[0].chains[1].residues[0:2].atoms.get("C*-backbone")
        self.assertEquals(result.stringRepr,"(hsg1:B:0-1:\\s\\C*-backbone)")
        self.assertEquals(len(result), 10)

    #       regexp -  NamedAtomSet(fail)
    def test_AtomSet_range_regexp_NamedAtomSet_fail_second(self):
        """
         test range regexp-NamedAtomSet_fail_second used in atoms.get(regexp-NamedAtomSet)
        """
        result = self.atoms.get("C*-mistake")
        self.assertEquals(result, AtomSet())

    #       regexp(fail) -  NamedAtomSet
    def test_AtomSet_range_regexp_fail_first_NamedAtomSet(self):
        """
         test range regexp_fail_first-NamedAtomSet used in atoms.get(regexp-NamedAtomSet)
        """
        result = self.atoms.get("FOO*-backbone")
        self.assertEquals(result, AtomSet())

    #       regexp(fail) -  NamedAtomSet(fail)
    def test_AtomSet_range_regexp_NamedAtomSet_fail_both(self):
        """
         test range regexp-NamedAtomSet_fail_both used in atoms.get(regexp-NamedAtomSet)
        """
        result = self.atoms.get("FOO*-mistake")
        self.assertEquals(result, AtomSet())

    #
    # RANGE: RELATIVE-REGEXP
    #       relative -  regexp
    def test_AtomSet_range_relative_regexp(self):
        """
         test range relative-regexp used in atoms.get(relative-regexp)
        """
        result = self.molecules[0].chains[1].residues[0:2].atoms.get("#1-N*")
        self.assertEquals(result.stringRepr,"(hsg1:B:0-1:\\s\\#1-N*)")
        #with correction which limited relative ranges to #num-#num only
        self.assertEquals(len(result), 16)

    #       relative -  regexp(fail)
    def test_AtomSet_range_relative_regexp_fail_second(self):
        """
         test range relative-regexp_fail_second used in atoms.get(relative-regexp)
        THIS was TOO SLOW to use all the atoms
        """
        #result = self.atoms.get("#1-Q*")
        result = self.molecules[0].allAtoms.get("#1-Q*")
        self.assertEquals(result, AtomSet())

    #       relative(fail) -  regexp
    def test_AtomSet_range_relative_fail_first_regexp(self):
        """
         test range relative_fail_first-regexp used in atoms.get(relative-regexp)
        THIS was TOO SLOW to use all the atoms
        """
        #result = self.atoms.get("#16666666-C*")
        result = self.molecules[0].allAtoms.get("#1666666-C*")
        self.assertEquals(result, AtomSet())

    #       relative(fail) -  regexp(fail)
    def test_AtomSet_range_relative_regexp_fail_both(self):
        """
         test range relative-regexp_fail_both used in atoms.get(relative-regexp)
        THIS was TOO SLOW to use all the atoms
        """
        #result = self.atoms.get("#166666-FOO*")
        result = self.molecules[0].allAtoms.get("#1666666-FOO*")
        self.assertEquals(result, AtomSet())

    # RANGE: RELATIVE-RELATIVE
    #       relative -  relative
    def test_AtomSet_range_relative_relative(self):
        """
         test range relative-relative used in atoms.get(relative-relative)
        """
        result = self.molecules[0].chains[1].residues[0:2].atoms.get("#1-#3")
        self.assertEquals(result.stringRepr,"(hsg1:B:0-1:\\s\\#1-#3)")
        #.flat #1 matches the first atom, #2 the second
        self.assertEquals(len(result), 6)

    #       relative -  relative(fail)
    def test_AtomSet_range_relative_relative_fail_second(self):
        """
         test range relative-relative_fail_second used in atoms.get(relative-relative)
        """
        result = self.molecules[0].chains[2].residues[0:2].atoms.get("#1-#333333")
        self.assertEquals(result, AtomSet())

    #       relative(fail) -  relative
    def test_AtomSet_range_relative_fail_first_relative(self):
        """
         test range relative_fail_first-relative used in atoms.get(relative-relative)
        """
        result = self.molecules[0].chains[2].residues[0:2].atoms.get("#1111-#3")
        self.assertEquals(result, AtomSet())

    #       relative(fail) -  relative(fail)
    def test_AtomSet_range_relative_relative_fail_both(self):
        """
         test range relative-relative_fail_both used in atoms.get(relative-relative)
        THIS was TOO SLOW 
        """
        result = self.molecules[0].chains[2].residues[0:2].atoms.get("#1111-#33333")
        self.assertEquals(result, AtomSet())

    # RANGE: RELATIVE-INDEX
    #       relative -  index
    def test_AtomSet_range_relative_index(self):
        """
         test range relative-index used in atoms.get(relative-index)
        """
        result = self.molecules[0].chains[1].residues[0:2].atoms.get("#1-3")
        self.assertEquals(result.stringRepr,"(hsg1:B:0-1:\\s\\#1-3)")
        #with correction which limited relative ranges to #num-#num only
        self.assertEquals(len(result), 4)


    #       relative -  index(fail)
    def test_AtomSet_range_relative_index_fail_second(self):
        """
         test range relative-index_fail_second used in atoms.get(relative-index)
        THIS was TOO SLOW 
        """
        result = self.molecules[0].chains[2].residues[0:2].atoms.get("#1-333333")
        self.assertEquals(result, AtomSet())


    #       relative(fail) -  index
    def test_AtomSet_range_relative_fail_first_index(self):
        """
         test range relative_fail_first-index used in atoms.get(relative-index)
        THIS was TOO SLOW 
        """
        result = self.molecules[0].chains[2].residues[0:2].atoms.get("#1666666-#3")
        self.assertEquals(result, AtomSet())

    #       relative(fail) -  index(fail)
    def test_AtomSet_range_relative_index_fail_both(self):
        """
         test range relative_index-fail_both used in atoms.get(relative-index)
        THIS was TOO SLOW 
        """
        result = self.molecules[0].chains[2].residues[0:2].atoms.get("#1666666-#333333")
        self.assertEquals(result, AtomSet())

    # RANGE: RELATIVE-NamedResSet
    #       relative -  NamedResSet
    def test_AtomSet_range_relative_NamedResSet(self):
        """
         test range relative-NamedResSet used in atoms.get(relative-NamedResSet)
        """
        result = self.molecules[0].chains[2].residues[0:2].atoms.get("#1-buried")
        self.assertEquals(result, AtomSet())

    #       relative -  NamedResSet(fail)
    def test_AtomSet_range_relative_NamedResSet_fail_second(self):
        """
         test range relative-NamedResSet_fail_second used in atoms.get(relative-NamedResSet)
        THIS was TOO SLOW 
        """
        result = self.molecules[0].chains[2].residues[0:2].atoms.get("#1-mistake")
        self.assertEquals(result, AtomSet())

    #       relative(fail) -  NamedResSet
    def test_AtomSet_range_relative_fail_first_NamedResSet(self):
        """
         test range relative_fail-first_NamedResSet used in atoms.get(relative-NamedResSet)
        THIS was TOO SLOW 
        """
        result = self.molecules[0].chains[2].residues[0:2].atoms.get("#1666-buried")
        self.assertEquals(result, AtomSet())

    #       relative(fail) -  NamedResSet(fail)
    def test_AtomSet_range_relative_NamedResSet_fail_both(self):
        """
        test range relative-NamedResSet_fail_both used in atoms.get(relative-NamedResSet)
        THIS was TOO SLOW 
        """
        result = self.molecules[0].chains[2].residues[0:2].atoms.get("#1666-mistake")
        self.assertEquals(result, AtomSet())


    # RANGE: RELATIVE-NamedAtomSet
    #       relative -  NamedAtomSet
    def test_AtomSet_range_relative_NamedAtomSet(self):
        """
         test range relative-NamedAtomSet used in atoms.get(relative-NamedAtomSet)
        THIS was TOO SLOW 
        """
        result = self.molecules[0].chains[1].residues[0:2].atoms.get("#1-backbone")
        self.assertEquals(result.stringRepr,"(hsg1:B:0-1:\\s\\#1-backbone)")
        #with correction which limited relative ranges to #num-#num only
        self.assertEquals(len(result), 11)


    #       relative -  NamedAtomSet(fail)
    def test_AtomSet_range_relative_NamedAtomSet_fail_second(self):
        """
         test range relative-NamedAtomSet_fail_second used in atoms.get(relative-NamedAtomSet)
        """
        result = self.molecules[0].chains[2].residues[0:2].atoms.get("#1-mistake")
        self.assertEquals(result, AtomSet())

    #       relative(fail) -  NamedAtomSet
    def test_AtomSet_range_relative_fail_first_NamedAtomSet(self):
        """
         test range relative_fail_first-NamedAtomSet used in atoms.get(relative-NamedAtomSet)
        """
        result = self.molecules[0].chains[2].residues[0:2].atoms.get("#166666-backbone")
        self.assertEquals(result, AtomSet())

    #       relative(fail) -  NamedAtomSet(fail)
    def test_AtomSet_range_relative_NamedAtomSet_fail_both(self):
        """
         test range relative-NamedAtomSet_fail_both used in atoms.get(relative-NamedAtomSet)
         #THIS was too slow
        """
        result = self.molecules[0].chains[2].residues[0:2].atoms.get("#166666-mistake")
        self.assertEquals(result, AtomSet())

    #
    # RANGE: INDEX-REGEXP
    #       index -  regexp
    def test_AtomSet_range_index_regexp(self):
        """
         test range index-regexp used in atoms.get(index-regexp)
        """
        result = self.atoms.get("1-C*")
        self.assertEquals(result.stringRepr,"(hsg1:::/+/2plv:::/+/1gyc:::\\s\\1-C*)")
        self.assertEquals(len(result), 12811)

    #       index -  regexp(fail)
    def test_AtomSet_range_index_regexp_fail_second(self):
        """
         test range index-regexp_fail_second used in atoms.get(index-regexp)
        """
        result = self.atoms.get("1-666666*")
        self.assertEquals(result, AtomSet())

    #       index(fail) -  regexp
    def test_AtomSet_range_index_fail_first_regexp(self):
        """
         test range index_fail_first-regexp used in atoms.get(index-regexp)
        """
        result = self.atoms.get("66666-1*")
        self.assertEquals(result, AtomSet())

    #       index(fail) -  regexp(fail)
    def test_AtomSet_range_index_regexp_fail_both(self):
        """
         test range index-regexp_fail_both used in atoms.get(index-regexp)
        """
        result = self.atoms.get("6666-FOO*")
        self.assertEquals(result, AtomSet())

    # RANGE: INDEX-RELATIVE
    #       index -  relative
    def test_AtomSet_range_index_relative(self):
        """
         test range index-relative used in atoms.get(index-relative)
         THIS WAS TOO SLOW
        """
        result = self.molecules[0].chains[0].residues[0:2].atoms.get("1-#2")
        #self.assertEquals(len(result), 2)
        #with correction which limited relative ranges to #num-#num only
        self.assertEquals(len(result), 8)

    #       index -  relative(fail)
    def test_AtomSet_range_index_relative_fail_second(self):
        """
         test range index-relative_fail_second used in atoms.get(index-relative)
         THIS WAS TOO SLOW
        """
        result = self.molecules[0].chains[2].residues[0:2].atoms.get("1-#222222")
        self.assertEquals(result, AtomSet())

    #       index(fail) -  relative
    def test_AtomSet_range_index_fail_first_relative(self):
        """
         test range index_fail_first-relative used in atoms.get(index-relative)
         THIS WAS TOO SLOW
        """
        result = self.molecules[0].chains[2].residues[0:2].atoms.get("1666666-#2")
        self.assertEquals(result, AtomSet())

    #       index(fail) -  relative(fail)
    def test_AtomSet_range_index_relative_fail_both(self):
        """
         test range index-relative_fail_both used in atoms.get(index-relative)
         THIS WAS TOO SLOW
        """
        result = self.molecules[0].chains[2].residues[0:2].atoms.get("1666666-#22222")
        self.assertEquals(result, AtomSet())

    # RANGE: INDEX-INDEX
    #       index -  index
    def test_AtomSet_range_index_index(self):
        """
         test range index-index used in atoms.get(index-index)
        """
        result = self.atoms.get("1-2")
        self.assertEquals(result.stringRepr,"(hsg1:::/+/2plv:::/+/1gyc:::\\s\\1-2)")
        self.assertEquals(len(result), 2)

    #       index -  index(fail)
    def test_AtomSet_range_index_index_fail_second(self):
        """
         test range index-index_fail_second used in atoms.get(index-index)
        """
        result = self.atoms.get("1-133333")
        self.assertEquals(result, AtomSet())

    #       index(fail) -  index
    def test_AtomSet_range_index_fail_first_index(self):
        """
         test range index_fail_first-index used in atoms.get(index-index)
        """
        result = self.atoms.get("166666-1")
        self.assertEquals(result, AtomSet())

    #       index(fail) -  index(fail)
    def test_AtomSet_range_index_index_fail_both(self):
        """
         test range index_index-fail_both used in atoms.get(index-index)
        """
        result = self.atoms.get("166666-188888")
        self.assertEquals(result, AtomSet())

    # RANGE: INDEX-NamedResSet
    #       index -  NamedResSet
    def test_AtomSet_range_index_NamedResSet(self):
        """
         test range index-NamedResSet used in atoms.get(index-NamedResSet)
        """
        result = self.atoms.get("1-buried")
        self.assertEquals(result, AtomSet())


    #       index -  NamedResSet(fail)
    def test_AtomSet_range_index_NamedResSet_fail_second(self):
        """
         test range index-NamedResSet_fail_second used in atoms.get(index-NamedResSet)
        """
        result = self.atoms.get("1-mistake")
        self.assertEquals(result, AtomSet())

    #       index(fail) -  NamedResSet
    def test_AtomSet_range_index_fail_first_NamedResSet(self):
        """
         test range index_fail-first_NamedResSet used in atoms.get(index-NamedResSet)
        """
        result = self.atoms.get("666666-buried")
        self.assertEquals(result, AtomSet())


    #       index(fail) -  NamedResSet(fail)
    def test_AtomSet_range_index_NamedResSet_fail_both(self):
        """
         test range index-NamedResSet_fail_both used in atoms.get(index-NamedResSet)
        """
        result = self.atoms.get("666666-mistake")
        self.assertEquals(result, AtomSet())


    # RANGE: INDEX-NamedAtomSet
    #       index -  NamedAtomSet
    def test_AtomSet_range_index_NamedAtomSet(self):
        """
         test range index-NamedAtomSet used in atoms.get(index-NamedAtomSet)
        """
        result = self.atoms.get("1-backbone")
        self.assertEquals(result.stringRepr,"(hsg1:::/+/2plv:::/+/1gyc:::\\s\\1-backbone)")
        self.assertEquals(len(result),  12668)


    #       index -  NamedAtomSet(fail)
    def test_AtomSet_range_index_NamedAtomSet_fail_second(self):
        """
         test range index-NamedAtomSet_fail_second used in atoms.get(index-NamedAtomSet)
        """
        result = self.atoms.get("1-mistake")
        self.assertEquals(result, AtomSet())

    #       index(fail) -  NamedAtomSet
    def test_AtomSet_range_index_fail_first_NamedAtomSet(self):
        """
         test range index_fail_first-NamedAtomSet used in atoms.get(index-NamedAtomSet)
        """
        result = self.atoms.get("666666-backbone")
        self.assertEquals(result, AtomSet())

    #       index(fail) -  NamedAtomSet(fail)
    def test_AtomSet_range_index_NamedAtomSet_fail_both(self):
        """
         test range index-NamedAtomSet_fail_both used in atoms.get(index-NamedAtomSet)
        """
        result = self.atoms.get("666666-mistake")
        self.assertEquals(result, AtomSet())



if __name__ == '__main__':
    test_cases = [
        'MoleculeSetSelectionTests',
        'ChainSetSelectionTests',
        'ResidueSetSelectionTests',
        'AtomSetSelectionTests',
    ]
    unittest.main( argv=([__name__ , ] + test_cases) )
    #unittest.main( argv=([__name__ , '-v' ] + test_cases) )
    #unittest.main()
