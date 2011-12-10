## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#
#
#
#
# $Id: test_distanceSelector.py,v 1.2 2007/07/24 17:30:40 vareille Exp $
#

import unittest
import numpy.oldnumeric as Numeric
from MolKit.distanceSelector import CloserThanVDWSelector, CloserThanVDWPlusConstantSelector
from MolKit import Read


class BaseTests(unittest.TestCase):

    def setUp(self):
        self.lig = Read('Data/ind_crystal.pdb')[0]
        self.lig.buildBondsByDistance()
        self.rec = Read('Data/hsg1.pdbqt')[0]
        self.rec.buildBondsByDistance()
    

    def tearDown(self):
        """
        clean-up
        """
        del(self.lig)
        del(self.rec)


    def test_CloserThanVDWSelector(self):
        """
        instantiate a CloserThanVDWSelector
        """
        d_sel = CloserThanVDWSelector()
        self.assertEquals(d_sel.__class__, CloserThanVDWSelector)


    def test_CloserThanVDWSelector_setupCutoff(self):
        """
         CloserThanVDWSelector cutoff is expected shape
        """
        cutoff = CloserThanVDWSelector().setupCutoff(self.rec.allAtoms, self.lig.allAtoms, 3.0)
        self.assertEquals(cutoff.shape[0], len(self.lig.allAtoms))
        self.assertEquals(cutoff.shape[1], len(self.rec.allAtoms))


    def test_CloserThanVDWSelector_select(self):
        """
        CloserThanVDWSelector.select returns expected number of atoms 
        """
        d_sel = CloserThanVDWSelector()
        resultD, distD = d_sel.select(self.lig.allAtoms, self.rec.allAtoms)
        #check that 14 atoms in the ligand are near 1 or more atoms in the receptor
        self.assertEquals(len(resultD.keys()), 14)
        d = {}
        for k in resultD.values():
            for item in k:
                d[item] = 0
        #check that 16 atoms in the receptor are near an atom in the ligand
        self.assertEquals(len(d.keys()), 16)
                

    def test_CloserThanVDWPlusConstantSelector(self):
        """
        instantiate a CloserThanVDWPlusConstantSelector
        """
        d_sel = CloserThanVDWPlusConstantSelector()
        self.assertEquals(d_sel.__class__, CloserThanVDWPlusConstantSelector)


    def test_CloserThanVDWPlusConstantSelector_setupCutoff(self):
        """
         CloserThanVDWPlusConstantSelector cutoff is expected shape
        """
        cutoff = CloserThanVDWPlusConstantSelector().setupCutoff(self.rec.allAtoms, self.lig.allAtoms, 3.0)
        self.assertEquals(cutoff.shape[0], len(self.lig.allAtoms))
        self.assertEquals(cutoff.shape[1], len(self.rec.allAtoms))


    def test_CloserThanVDWPlusConstantSelector_select(self):
        """
        CloserThanVDWPlusConstantSelector.select returns expected number of atoms 
        """
        d_sel = CloserThanVDWPlusConstantSelector(constant=1)
        resultD, distD = d_sel.select(self.lig.allAtoms, self.rec.allAtoms)
        #check that 45 atoms in the ligand are near 1 or more atoms in the receptor
        self.assertEquals(len(resultD.keys()), 45)
        d = {}
        for k in resultD.values():
            for item in k:
                d[item] = 0
        #check that 16 atoms in the receptor are near an atom in the ligand
        # when using sum of vdw PLUS 1.0 angstrom
        self.assertEquals(len(d.keys()), 93)
                



if __name__ == '__main__':
    unittest.main()
