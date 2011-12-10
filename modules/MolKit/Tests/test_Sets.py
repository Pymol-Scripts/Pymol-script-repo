## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#
#
#
#
# $Header: /opt/cvs/python/packages/share1.5/MolKit/Tests/test_Sets.py,v 1.2 2007/07/24 17:30:40 vareille Exp $
#

import unittest
import numpy.oldnumeric as Numeric
from string import split
from MolKit.sets import Sets
from MolKit import Read
from MolKit.molecule import AtomSet
from MolKit.protein import ResidueSet, ChainSet


class SetsConstructorTests(unittest.TestCase):
    """
    tests of Sets class constructor options
"""    
    def test_constructor(self):
        """
        test creating an empty sets instance
        """
        sets = Sets()
        self.assertEquals(sets.__class__, Sets)



class SetsTests(SetsConstructorTests):
    """
    tests of Sets class methods
"""    

    def setUp(self):
        self.mol = Read('Data/1crn.pdb')[0]
        self.sets = Sets()


    def tearDown(self):
        """
        clean-up
        """
        del(self.mol)


    def test_add(self):
        """
        test adding an AtomSet
        """
        k = 'myfirst set'
        atSet = self.mol.allAtoms
        sets = self.sets
        sets.add(k, atSet)
        self.assertEquals(len(sets.keys()), 1)
        self.assertEquals(sets.keys()[0], k)
        self.assertEquals(sets.values()[0], atSet)


    def test_add_residues(self):
        """
        test adding a ResidueSet
        """
        sets = self.sets
        k = 'myfirst set'
        resSet = self.mol.chains.residues
        sets.add(k, resSet)
        self.assertEquals(len(sets.keys()), 1)
        self.assertEquals(sets.keys()[0], k)
        self.assertEquals(sets.values()[0], resSet)


    def test_remove(self):
        """
        test removing an AtomSet by key
        """
        sets = self.sets
        k = 'myfirst set'
        atSet = self.mol.allAtoms
        sets.add(k, atSet)
        sets.remove(k)
        self.assertEquals(sets, {})


    def test_remove_by_instance(self):
        """
        test removing an AtomSet by instance
        """
        sets = self.sets
        k = 'myfirst set'
        atSet = self.mol.allAtoms
        sets.add(k, atSet)
        sets.removeByInstance(atSet)
        self.assertEquals(sets, {})


    def test_add_twice_overwrites(self):
        """
        test adding with the same key overwrites;
        """
        sets = self.sets
        k = 'myfirst set'
        atSet = self.mol.allAtoms
        sets.add(k, atSet)
        atSet10 = self.mol.allAtoms[:10]
        sets.add(k, atSet10)
        self.assertEquals(len(sets), 1)
        self.assertEquals(sets.values()[0],  atSet10)


    def test_add_twice_overwrite_False(self):
        """
        test adding with the same key with overwrite=False 
        raises AssertionError;
        """
        sets = self.sets
        k = 'myfirst set'
        atSet = self.mol.allAtoms
        sets.add(k, atSet)
        self.assertRaises(AssertionError, sets.add, (k, atSet), {'overwrite':False})


    def test_get(self):
        """
        test get sets of specified TreeNodeSet type
        """
        sets = self.sets
        k = 'myfirst set'
        atSet = self.mol.allAtoms
        sets.add(k, atSet)
        kr = 'allResidues'
        resSet = self.mol.chains.residues
        sets.add(kr, resSet)
        kr2 = '10Residues'
        resSet10 = self.mol.chains.residues[:10]
        sets.add(kr2, resSet10)
        #check get returns 1 item
        self.assertEquals(sets.get(AtomSet), {k:atSet})
        #check get returns >1 item
        self.assertEquals(len(sets.get(ResidueSet)), 2)
        #check what get returns 
        #when there are no items of specified TreeNodeSet type
        self.assertEquals(sets.get(ChainSet), {})



if __name__ == '__main__':
    unittest.main()
