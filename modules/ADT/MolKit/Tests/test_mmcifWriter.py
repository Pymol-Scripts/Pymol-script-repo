#$Header: /opt/cvs/python/packages/share1.5/MolKit/Tests/test_mmcifWriter.py,v 1.3 2006/04/19 00:52:07 sargis Exp $
#
#$Id: test_mmcifWriter.py,v 1.3 2006/04/19 00:52:07 sargis Exp $ 
import unittest, os
from string import split
from MolKit.mmcifParser import MMCIFParser
from MolKit.mmcifWriter import MMCIFWriter
class MMCIFParserBaseTests(unittest.TestCase):
               
    def test_writer_1(self):
        """
        Checks that MMCIFParser can read the file written by MMCIFWriter
        """
        parser = MMCIFParser(filename='Data/1CRN.cif')
        mol = parser.parse()
        writer = MMCIFWriter()
        writer.write('Data/1CRN_.cif',mol)
        new_parser = MMCIFParser(filename='Data/1CRN_.cif')
        mol = parser.parse()
        os.remove('Data/1CRN_.cif')
    

    def test_write_2(self):
        """
        Tests that PdbWriter can write molecule parsed by MMCIFParser
        """
        # read a molecule
        parser = MMCIFParser(filename='Data/1CRN.cif')
        mol = parser.parse()
        # instanciate a PdbWriter and call the write method with the
        # default arguments
        from MolKit.pdbWriter import PdbWriter
        writer = PdbWriter()
        writer.write('Data/1crn_mmcifwriter.pdb', mol)
        os.remove('Data/1crn_mmcifwriter.pdb')

if __name__ == '__main__':
    unittest.main()
