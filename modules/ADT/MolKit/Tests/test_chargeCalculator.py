## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#
#
#
#
# $Id: test_chargeCalculator.py,v 1.2 2007/07/24 17:30:40 vareille Exp $
#

import unittest
import numpy.oldnumeric as Numeric
from string import split
from MolKit.chargeCalculator import GasteigerChargeCalculator
from MolKit.chargeCalculator import KollmanChargeCalculator
from MolKit import Read


class BaseTests(unittest.TestCase):

    def setUp(self):
        self.mol = Read('Data/1crn.pdb')[0]
        self.mol.buildBondsByDistance()
    

    def tearDown(self):
        """
        clean-up
        """
        del(self.mol)


    def test_gasteiger_calc_constructor(self):
        """
        instantiate an GasteigerChargeCalculator
        """
        g_calc = GasteigerChargeCalculator()
        self.assertEquals(g_calc.__class__, GasteigerChargeCalculator)


    def test_gasteiger_calc_addCharges(self):
        """
         GasteigerChargeCalculator addCharges to each atom, giving expected total
        """
        atlen = len(self.mol.allAtoms)
        g_calc = GasteigerChargeCalculator()
        g_calc.addCharges(self.mol.allAtoms)
        g_ats = self.mol.allAtoms.get(lambda x: x.chargeSet=='gasteiger') 
        self.assertEquals(atlen, len(g_ats))
        total_charge = Numeric.add.reduce(self.mol.allAtoms.charge)
        self.assertAlmostEquals(total_charge, 1.221245e-15, 6)


    def test_kollman_calc_constructor(self):
        """
        instantiate an KollmanChargeCalculator
        """
        k_calc = KollmanChargeCalculator()
        self.assertEquals(k_calc.__class__, KollmanChargeCalculator)


    def test_Kollman_calc_addCharges(self):
        """
         KollmanChargeCalculator addCharges to each atom, giving expected total
        """
        atlen = len(self.mol.allAtoms)
        k_calc = KollmanChargeCalculator()
        k_calc.addCharges(self.mol.allAtoms)
        k_ats = self.mol.allAtoms.get(lambda x: x.chargeSet=='Kollman') 
        self.assertEquals(atlen, len(k_ats))
        total_charge = Numeric.add.reduce(self.mol.allAtoms.charge)
        self.assertAlmostEquals(total_charge, -19.554, 4)


    def test_Kollman_calc_from_bug(self):
        """
         test KollmanChargeCalculator can add charges to 1a30
        """
        mol = Read("Data/1a30.pdb")[0]
        atlen = len(mol.allAtoms)
        k_calc = KollmanChargeCalculator()
        k_calc.addCharges(mol.allAtoms)
        k_ats = mol.allAtoms.get(lambda x: x.chargeSet=='Kollman') 
        self.assertEquals(atlen, len(k_ats))
        total_charge = Numeric.add.reduce(mol.allAtoms.charge)
        self.assertAlmostEquals(total_charge, 5.7060, 4)

if __name__ == '__main__':
    unittest.main()
