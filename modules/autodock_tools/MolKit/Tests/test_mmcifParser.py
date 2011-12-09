#$Header: /opt/cvs/python/packages/share1.5/MolKit/Tests/test_mmcifParser.py,v 1.14 2008/09/12 17:33:56 sargis Exp $
#
#$Id: test_mmcifParser.py,v 1.14 2008/09/12 17:33:56 sargis Exp $ 
import unittest
from string import split
from MolKit.mmcifParser import MMCIFParser


mols = None
class MMCIFParserBaseTests(unittest.TestCase):

    def initParser(self):
        global mols
        if mols is None:
            parser = MMCIFParser(filename='./Data/2PLV.cif')
            mols = parser.parse()
        self.mols = mols

    def setUp(self):
        if not hasattr(self, 'mols'):
            self.initParser() 

    def tearDown(self):
        pass
           
    def test_number_of_chains(self):
        """Check that there are 6 chains"""
        self.assertEqual(len(self.mols.chains),6)

    def test_number_of_residues(self):
        """Check that there are 1335 residues"""
        self.assertEqual(len(self.mols.chains.residues),857)
        
    def test_number_of_atoms(self):
        """Check that there are 7162 atoms"""
        self.assertEqual(len(self.mols.chains.residues.atoms),7162)

    def test_chain_id(self):
        """Check that chain id are in following order 'A','B','C','D','E','F'"""
        self.assertEqual(self.mols.chains.id,['A','B','C','D','E','F'])

    def test_number_of_atoms(self):
        """Check that there are 7162 atoms"""
        self.assertEqual(len(self.mols.chains.residues.atoms),7162)
        
    def test_parseSSData(self):
        """Check that parseSSData() parses secondary structure"""
        parser1 = MMCIFParser(filename='Data/1CRN.cif')
        mols1 = parser1.parse()
        ssDataForMol = parser1.parseSSData(mols1)   

    
    def test_parseFileWrongFormat(self):
        parser1 = MMCIFParser("Data/fakePDB.pdb")
        mols1 = parser1.parse()
        assert mols1 is None or len(mols1)==0
    
    def test_parseNAD(self):
        parser1 = MMCIFParser("Data/1GGA_.cif")
        mol1 = parser1.parse()
        P = mol1.chains[-1].residues.atoms[0].chemElem
        assert P is 'P' 
        
if __name__ == '__main__':
    unittest.main()
