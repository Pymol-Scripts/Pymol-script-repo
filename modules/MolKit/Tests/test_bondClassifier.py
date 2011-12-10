#
#
#
#
# $Id: test_bondClassifier.py,v 1.3 2003/11/03 19:20:09 rhuey Exp $
#

import unittest
from string import split
from MolKit.bondClassifier import BondClassifier
from MolKit.bondSelector import AmideBondSelector
from MolKit.bondSelector import PeptideBackBoneBondSelector
from MolKit.bondSelector import CycleBondSelector
from MolKit.bondSelector import LeafBondSelector



class BaseTests(unittest.TestCase):
    def setUp(self):
        from MolKit import Read
        self.mol = Read('Data/1crn_2.pdb')[0]
        self.mol.buildBondsByDistance()
        
    
    def tearDown(self):
        """
        clean-up
        """
        del self.mol


    def compareBondSets(self, set1, set2):
        """
        check that two bond sets have same bonds
        """
        self.assertEquals(len(set1), len(set2))
        idlist1 = [id(x) for x in set1]
        idlist2 = [id(x) for x in set2]
        idlist1.sort()
        idlist2.sort()
        for b1,b2 in zip(idlist1, idlist2):
            self.assertEquals(b1,b2)


    def test_constructor(self):
        """
        instantiate an BondClassifier
        """
        bndClassifier = BondClassifier()
        self.assertEquals(bndClassifier.__class__, BondClassifier)


    def test_constructorOptions(self):
        """
         test possible constructor options
            options = {key: bondSelector, key2:bondSelector2,....}
        """
        bndClassifier = BondClassifier({'amide': AmideBondSelector()})
        self.assertEquals(bndClassifier.__class__, BondClassifier)


    def test_inputParameters(self):
        """
         test nul input to classify
            
        """
        bndClassifier = BondClassifier({'amide': AmideBondSelector()})
        self.assertRaises(AssertionError, bndClassifier.classify)


    def test_classifierVSselectors(self):
        """
         make sure that the classifier returns the same results as its
         selectors
        """
        bndClassifier = BondClassifier({'amide': AmideBondSelector(),
                                        'cycle': CycleBondSelector(),
                                        'leaf': LeafBondSelector(),
                                        'peptide': PeptideBackBoneBondSelector(),
                                        })
        #all the bonds in the molecule as a BondSet
        bnds = self.mol.allAtoms.bonds[0]

        #get the specified bonds as a BondSet
        localDict = {}
        localDict['amide'] = AmideBondSelector().select(bnds)
        localDict['cycle'] = CycleBondSelector().select(bnds)
        localDict['leaf'] =  LeafBondSelector().select(bnds)
        localDict['peptide'] = PeptideBackBoneBondSelector().select(bnds)

        #make the classifier do the same thing
        resultDict = bndClassifier.classify(bnds)

        for k in resultDict.keys():
            self.compareBondSets(localDict[k], resultDict[k])


if __name__ == '__main__':
    unittest.main()

