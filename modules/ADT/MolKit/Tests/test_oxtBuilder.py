#
#
#
#
# $Id: test_oxtBuilder.py,v 1.1 2005/09/15 16:12:38 rhuey Exp $
#

import unittest
from string import split
from MolKit.oxtBuilder import OxtBuilder


class BaseTests(unittest.TestCase):
    def setUp(self):
        from MolKit import Read
        self.mol = Read('Data/2plv_no_oxt.pdb')[0]
        self.mol.buildBondsByDistance()
    

    def tearDown(self):
        """
        clean-up
        """
        del(self.mol)


    def test_constructor(self):
        """
        instantiate an HydrogenBuilder
        """
        oxt_builder = OxtBuilder()
        self.assertEquals(oxt_builder.__class__, OxtBuilder)


#    def test_constructorOptions(self):
#        """
#         test possible constructor options
#            options = NONE!!
#        """
#    
#        oxt_builder = OxtBuilder(???=???)
#        self.assertEquals(oxt_builder.__class__, OxtBuilder)


    def test_add_oxt(self):
        """
         test add_oxt 
        """
        beforeLen = len(self.mol.allAtoms)
        oxt_builder = OxtBuilder()
        #the last chain is all waters so skip it
        for ch in self.mol.chains[:-1]:
            catom = ch.residues[-1].atoms[2]
            new_at = oxt_builder.add_oxt(catom)
            #print "added ", new_at.full_name()
        afterLen = len(self.mol.allAtoms)
        #print "beforeLen=", beforeLen, ' afterLen=', afterLen
        self.assertEquals(beforeLen+len(self.mol.chains[:-1]), afterLen)






if __name__ == '__main__':
    unittest.main()
