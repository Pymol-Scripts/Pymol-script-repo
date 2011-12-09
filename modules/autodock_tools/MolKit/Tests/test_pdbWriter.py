#############################################################################
#
# Author: Sophie I. COON, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################
#
# $Header: /opt/cvs/python/packages/share1.5/MolKit/Tests/test_pdbWriter.py,v 1.8 2010/08/25 18:43:14 annao Exp $
#
# $Id: test_pdbWriter.py,v 1.8 2010/08/25 18:43:14 annao Exp $
#

import sys, string, os
import unittest
from MolKit.pdbWriter import PdbWriter, PdbqWriter, PdbqsWriter, PdbqtWriter
try:
    import stride
    haveStride = True
except:
    haveStride = False

class PDBWriterTests(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        if hasattr(self, 'mol'): del self.mol


    def test_Write_1(self):
        """
        Test the default option of the write method of a PdbWriter.
        Just write the ATOM, HETATM and CONECT records and bondOrigin
        is File.
        """
        # read a molecule
        from MolKit import Read
        self.mol = Read("Data/1crn.pdb")[0]

        # instanciate a PdbWriter and call the write method with the
        # default arguments
        from MolKit.pdbWriter import PdbWriter
        writer = PdbWriter()
        writer.write('Data/1crn_writer.pdb', self.mol)
        # This should write only the ATOM and CONECT records
        # The TER records are automatically created with the ATOM records

        # 1- Make sure that the file has been created
        import os
        self.failUnless(os.path.exists('Data/1crn_writer.pdb'))
        # 2- Make sure that the file created the proper records
        # Get the default records from the new pdb files
        from MolKit.pdbParser import PdbParser
        nmol = Read("Data/1crn_writer.pdb")[0]

        # COMPARE 
        ncoords = nmol.allAtoms.coords
        ocoords = self.mol.allAtoms.coords
        self.assertEqual(ncoords, ocoords)

        nname = nmol.allAtoms.name
        oname = nmol.allAtoms.name
        self.assertEqual(nname, oname)

        nbonds = nmol.allAtoms.bonds[0]
        obonds = self.mol.allAtoms.bonds[0]
        self.assertEqual(len(nbonds), len(obonds))
        nbatms = nbonds.atom1+nbonds.atom2
        nbatms = nbatms.uniq()
        nbatms.sort()
        
        nbatms = nbonds.atom1+nbonds.atom2
        nbatms = nbatms.uniq()
        nbatms.sort()
        os.system("rm Data/1crn_writer.pdb")

    def test_Write_2(self):
        """
        Test writing out CONECT records for the bonds described in
        the pdb file and built by distance.
        """
        from MolKit import Read
        self.mol = Read("Data/1crn.pdb")[0]
        self.mol.buildBondsByDistance()
        # instanciate a PdbWriter and call the write method with the
        # default arguments
        from MolKit.pdbWriter import PdbWriter
        writer = PdbWriter()
        writer.write('Data/1crn_writer.pdb', self.mol)
        from MolKit.pdbParser import PdbParser
        nmol = Read("Data/1crn_writer.pdb")[0]
        os.system("rm Data/1crn_writer.pdb")

    def test_Write_3(self):
        """
        Test writing out CONECT records for the bonds described in
        the pdb file eventhough the bonds built by distance exist.
        """
        from MolKit import Read
        self.mol = Read("Data/1crn.pdb")[0]
        self.mol.buildBondsByDistance()
        # instanciate a PdbWriter and call the write method with the
        # default arguments
        from MolKit.pdbWriter import PdbWriter
        writer = PdbWriter()
        writer.write('Data/1crn_writer.pdb', self.mol)
        from MolKit.pdbParser import PdbParser
        nmol = Read("Data/1crn_writer.pdb")[0]
    
        os.system("rm Data/1crn_writer.pdb")

    def test_Write_4(self):
        """
        Test writing out the ATOM records and the CONECT records
        without the HETATM
        """
        from MolKit import Read
        self.mol = Read("Data/1bsr.pdb")[0]
        # instanciate a PdbWriter and call the write method with the
        # default arguments
        from MolKit.pdbWriter import PdbWriter
        writer = PdbWriter()
        writer.write('Data/1bsr_writer.pdb', self.mol,
                     records=['ATOM', 'CONECT'], bondOrigin=('File',))
        nmol = Read("Data/1bsr_writer.pdb")[0]

        nhoh = nmol.chains.residues.get(lambda x: x.type=="HOH")
        ohoh = self.mol.chains.residues.get(lambda x: x.type=="HOH")
        self.failUnless(len(nhoh)==0)
        os.system("rm Data/1bsr_writer.pdb")
        

    def test_Write_5(self):
        from MolKit import Read
        self.mol = Read("Data/1crn.pdb")[0]
        self.mol.secondaryStructureFromFile()
        from MolKit.pdbWriter import PdbWriter
        writer = PdbWriter()
        writer.write('Data/1crn_writer.pdb', self.mol,
                     records=['ATOM', 'HETATM', 'CONECT',
                              'TURN','HELIX', 'SHEET'],
                     bondOrigin='all', ssOrigin='File')
        # Make sure that the ss information has been written out
        # properly.
        nmol = Read("Data/1crn_writer.pdb")[0]
        nmol.secondaryStructureFromFile()
        nsset = nmol.chains[0].secondarystructureset
        osset = self.mol.chains[0].secondarystructureset
        self.assertEqual(len(nsset), len(osset))
        for nss, oss in map(None, nsset, osset):
            self.assertEqual(nss.name, oss.name)
            self.assertEqual(nss.start.name,oss.start.name)
            self.assertEqual(nss.end.name, oss.end.name)
            self.assertEqual(len(nss.residues), len(oss.residues))
        os.system("rm Data/1crn_writer.pdb")

    def test_Write_6(self):
        
        if not haveStride: return
        from MolKit import Read
        self.mol = Read("Data/1crn.pdb")[0]
        self.mol.secondaryStructureFromStride()
        from MolKit.pdbWriter import PdbWriter
        writer = PdbWriter()
        writer.write('Data/1crn_writer.pdb', self.mol,
                     records=['ATOM', 'HETATM', 'CONECT',
                              'TURN','HELIX', 'SHEET'],
                     bondOrigin='all', ssOrigin='Stride')
        # Make sure that the ss information has been written out
        # properly.
        nmol = Read("Data/1crn_writer.pdb")[0]
        nmol.secondaryStructureFromFile()
        nsset = nmol.chains[0].secondarystructureset
        osset = self.mol.chains[0].secondarystructureset
        self.assertEqual(len(nsset), len(osset))
        for nss, oss in map(None, nsset, osset):
            self.assertEqual(nss.name, oss.name)
            self.assertEqual(nss.start.name,oss.start.name)
            self.assertEqual(nss.end.name, oss.end.name)
            self.assertEqual(len(nss.residues), len(oss.residues))
        os.system("rm Data/1crn_writer.pdb")
        
    def test_Write_7(self):
        """
        Test the default option of the write method of a PdbWriter.
        Just write the ATOM, HETATM and CONECT records and bondOrigin
        is File.
        """
        # read a molecule
        from MolKit import Read
        self.mol = Read("Data/1crn.pdb")[0]

        # instanciate a PdbWriter and call the write method with the
        # default arguments
        from MolKit.pdbWriter import PdbWriter
        writer = PdbWriter()
        self.mol.write('Data/1crn_writer.pdb', self.mol, writer)
        # This should write only the ATOM and CONECT records
        # The TER records are automatically created with the ATOM records

        # 1- Make sure that the file has been created
        import os
        self.failUnless(os.path.exists('Data/1crn_writer.pdb'))
        # 2- Make sure that the file created the proper records
        # Get the default records from the new pdb files
        from MolKit.pdbParser import PdbParser
        nmol = Read("Data/1crn_writer.pdb")[0]

        # COMPARE 
        ncoords = nmol.allAtoms.coords
        ocoords = self.mol.allAtoms.coords
        self.assertEqual(ncoords, ocoords)

        nname = nmol.allAtoms.name
        oname = nmol.allAtoms.name
        self.assertEqual(nname, oname)

        nbonds = nmol.allAtoms.bonds[0]
        obonds = self.mol.allAtoms.bonds[0]
        self.assertEqual(len(nbonds), len(obonds))
        nbatms = nbonds.atom1+nbonds.atom2
        nbatms = nbatms.uniq()
        nbatms.sort()
        
        nbatms = nbonds.atom1+nbonds.atom2
        nbatms = nbatms.uniq()
        nbatms.sort()
        os.system("rm Data/1crn_writer.pdb")


    def test_Write_4Char_AtomNames(self):
        """ check that atoms with names of 4+char are written properly"""
        from MolKit import Read
        self.mol = Read("Data/HTet.mol2")[0]
        from MolKit.pdbWriter import PdbWriter
        writer = PdbWriter()
        writer.write('Data/test_HTet.pdb', self.mol)
        testmol = Read("Data/test_HTet.pdb")[0]
        self.assertEqual(len(self.mol.allAtoms), len(testmol.allAtoms))
        #check the 4char atoms have correct element types
        self.assertEqual(testmol.allAtoms[0].element, 'C')
        self.assertEqual(testmol.allAtoms[1].element, 'Cl')
        self.assertEqual(testmol.allAtoms[22].element, 'N')
        self.assertEqual(testmol.allAtoms[23].element, 'H')
        #check the 4char atoms have reasonable names
        self.assertEqual(testmol.allAtoms[0].name, 'C11')
        self.assertEqual(testmol.allAtoms[1].name, 'Cl22')
        self.assertEqual(testmol.allAtoms[22].name, 'N61')
        self.assertEqual(testmol.allAtoms[23].name, 'HN61')
        os.system("rm Data/test_HTet.pdb")
        




class PdbWriterTest(unittest.TestCase):

    def setUp(self):
        pass
    

    def tearDown(self):
        """
        clean-up
        """
        pass


    def compare(self, file1, file2):
        """
        check that two written files have the same content
        NB: RECORDS ARE NOT THE SAME SO SKIP CHECKING THEM
        ALSO: CONECT records are not the same so skip them
        """
        lines1 = open(file1).readlines()
        atomlines1 = []
        for l in lines1:
            if string.find(l, 'ATOM')==0:
                atomlines1.append(l)
        lines2 = open(file2).readlines()
        atomlines2 = []
        for l in lines2:
            if string.find(l, 'ATOM')==0:
                atomlines2.append(l)
        if len(atomlines1)!=len(atomlines2):
            return False, 'lengths=', len(atomlines1), ' ', len(atomlines2)
        for i in range(len(atomlines1)):
            l1 = atomlines1[i]
            l2 = atomlines2[i]
            #ll_1 = string.split(l1)
            #ll_2 = string.split(l2)
            #if len(ll_1)!=len(ll_2):
                #return False, 'number of items'
            #if atomlines1[i][:75]!=atomlines2[i][:75]:
            if l1[:30]!=l2[:30]:
                print l1[:30],' v. ',  l2[:30]
                return False, 'start ' + str(i)
            if float(l1[30:38])!=float(l2[30:38]):
                return False, 'xcoord ' + str(i)
            if float(l1[38:46])!=float(l2[38:46]):
                return False, 'ycoord ' + str(i)
            if float(l1[46:54])!=float(l2[46:54]):
                return False, 'zcoord ' + str(i)
            if l1[54:70]!=l2[54:70]:
                return False, 'end ' + str(i)
        return True, None

class PdbWriterInitTest(PdbWriterTest):


    def test_constructor(self):
        """
        instantiate an PdbWriter
        """
        writer = PdbWriter()
        self.assertEquals(writer.__class__, PdbWriter)


class PdbWriterAtomLinesTest(PdbWriterTest):

    def setUp(self):
        from MolKit import Read
        self.mol = Read('Data/1crn.pdb')[0]
        self.mol.buildBondsByDistance()

    
    def tearDown(self):
        del(self.mol)


    def test_write(self):
        """
        test writing a pdb file
        """
        writer = PdbWriter()
        writer.write('test_pdbWriter.pdb', self.mol, bondOrigin=('File',))
        ans, errors = self.compare('Data/1crn.pdb', 'test_pdbWriter.pdb') 
        self.assertEquals(errors, None)
        self.assertEquals(ans, True)



class PdbqWriterAtomLinesTest(PdbWriterTest):


    def setUp(self):
        from MolKit import Read
        self.mol = Read('Data/1crn_Hs.pdbq')[0]
        self.mol.buildBondsByDistance()
    

    def tearDown(self):
        del(self.mol)


    def test_write(self):
        """
        test writing a pdbq file
        """
        writer = PdbqWriter()
        writer.write('test_pdbqWriter.pdbq', self.mol, bondOrigin=('File',))
        ans, errors = self.compare('Data/1crn_Hs.pdbq',
                                   'test_pdbqWriter.pdbq') 
        self.assertEquals(errors, None)
        self.assertEquals(ans, True)



class PdbqsWriterAtomLinesTest(PdbWriterTest):

    def setUp(self):
        from MolKit import Read
        self.mol = Read('Data/1crn.pdbqs')[0]
        self.mol.buildBondsByDistance()
    

    def tearDown(self):
        del(self.mol)


    def test_write(self):
        """
        test writing a pdbqs file
        """
        writer = PdbqsWriter()
        writer.write('test_pdbqsWriter.pdbqs', self.mol, bondOrigin=('File',))
        ans, errors = self.compare('Data/1crn.pdbqs', 'test_pdbqsWriter.pdbqs') 
        self.assertEquals(errors, None)
        self.assertEquals(ans, True)


class PdbqtWriterAtomLinesTest(PdbWriterTest):

    def setUp(self):
        from MolKit import Read
        self.mol = Read('Data/hsg1.pdbqt')[0]
        self.mol.buildBondsByDistance()
    

    def tearDown(self):
        del(self.mol)


    def test_write(self):
        """
        test writing a pdbqs file
        """
        writer = PdbqtWriter()
        writer.write('test_pdbqtWriter.pdbqt', self.mol, bondOrigin=('File',))
        ans, errors = self.compare('Data/hsg1.pdbqt', 'test_pdbqtWriter.pdbqt') 
        self.assertEquals(errors, None)
        self.assertEquals(ans, True)




if __name__ == '__main__':
    unittest.main()
