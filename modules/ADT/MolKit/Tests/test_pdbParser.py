#############################################################################
#
# Author: Sophie I. COON, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################
#
# $Header: /opt/cvs/python/packages/share1.5/MolKit/Tests/test_pdbParser.py,v 1.12 2009/07/10 22:22:46 rhuey Exp $
#
# $Id: test_pdbParser.py,v 1.12 2009/07/10 22:22:46 rhuey Exp $
#


import sys
from mglutil.regression import testplus

def test_parseMol_BasicRec():
    from MolKit.pdbParser import PdbParser
    parser = PdbParser('Data/1crn.pdb')
    mols = parser.parse()


def test_parseMol_HelixRec():
    from MolKit.pdbParser import PdbParser
    parser = PdbParser('Data/1crn.pdb')
    parser.SetReadOptions('HELIX')
    assert 'HELIX' in parser.pdbRecordParser.keys()
    mols = parser.parse()
    assert parser.recordInfo.has_key('HELIX')

def test_parseSSData():
    from MolKit.pdbParser import PdbParser
    parser = PdbParser('Data/1crn.pdb')
    mols = parser.parse()
    ssDataForMol = parser.parseSSData(mols[0])


def test_parseSSData_2():
    """ Calling parseSSData after setting the readOption to parse HELIX,
    STRAND, TURN, etc..."""
    from MolKit.pdbParser import PdbParser
    parser = PdbParser('Data/1crn.pdb')
    parser.SetReadOptions('HELIX')
    parser.SetReadOptions('SHEET')
    parser.SetReadOptions('TURN')
    mols = parser.parse()
    assert parser.recordInfo.has_key('HELIX')
    assert parser.recordInfo.has_key('SHEET')
    assert parser.recordInfo.has_key('TURN')
    ssDataForMol = parser.parseSSData(mols[0])

def test_parseHYDBND_1():
    """
    """
    from MolKit.pdbParser import PdbParser
    parser = PdbParser("Data/smfx_badHYDBND.pdb")
    parser.PDBtags.remove("HYDBND")
    del parser.pdbRecordParser['HYDBND']
    parser.parse()
    failure = False
    hdbndLines = parser.getRecords(parser.allLines, "HYDBND")
    try:
        parser.parse_PDB_HYDBND(hbndLines)
    except:
        failure = True
    assert failure 

def test_parseHYDBND_2():
    """
    """
    from MolKit.pdbParser import PdbParser
    parser = PdbParser("Data/smfx_goodHYDBND.pdb")
    parser.parse()
    mol = parser.mol
    atmsInHBond = mol.allAtoms.get(lambda x: hasattr(x, 'hbonds'))
    hbonds = atmsInHBond.hbonds
    assert len(hbonds) == 21
    

def test_parseHYDBND_3():
    """
    check that HYDBND records are read correctly
    when names start with a space (eg residue name "  A1")
    """
    from MolKit.pdbParser import PdbParser
    parser = PdbParser("Data/bdna_HS.pdb")
    m = parser.parse()[0]
    assert len(m.allAtoms.get(lambda x: hasattr(x, 'hbonds')))==110



def test_parsePDBQ():
    """
    check that pdbq format files can be read
    """
    from MolKit.pdbParser import PdbqParser
    parser = PdbqParser("Data/1ent.pdbqs")
    m = parser.parse()[0]
    assert len(m.allAtoms)==2901


def test_parsePDBQS():
    """
    check that pdbqs format files can be read
    """
    from MolKit.pdbParser import PdbqsParser
    parser = PdbqsParser("Data/1ent.pdbqs")
    m = parser.parse()[0]
    assert len(m.allAtoms)==2901


def test_parsePDBQT():
    """
    check that pdbqt format files can be read
    """
    from MolKit.pdbParser import PdbqtParser
    parser = PdbqtParser("Data/1ent_rec.pdbqt")
    m = parser.parse()[0]
    assert len(m.allAtoms)==2901

## commented out the test because it does not work properly on older version of MacOS (ppcDarwin7)
## def test_parseFileWithNoRead():
##     from MolKit.pdbParser import PdbParser
##     import os
##     # Change the read permission on the file before reading it
##     os.system("chmod 220 Data/1crn.pdb")
##     parser = PdbParser("Data/1crn.pdb")
##     mols = parser.parse()
##     assert mols is None
##     # put back the read permission on the file.
##     os.system("chmod 644 Data/1crn.pdb")


def test_parseFileWrongFormat():
    from MolKit.pdbParser import PdbParser
    parser = PdbParser("Data/fakePDB.pdb")
    mols = parser.parse()
    assert mols is None or len(mols)==0


    parser = PdbParser("Data/hpi1s.mol2")
    mols = parser.parse()
    assert mols is None or len(mols)==0

    from MolKit.mol2Parser import Mol2Parser
    parser = Mol2Parser("Data/1crn.pdb")
    mols=parser.parse()
    assert mols is None or len(mols)==0

def test_parsegetRecords():
    # Test the change made in the lambda function of getRecords method
    from MolKit.pdbParser import PdbParser
    parser = PdbParser("Data/1crn.pdb")
    parser.readFile()
    parser.getKeys()
    
    atmRec = parser.getRecords(parser.allLines,"ATOM")
    assert len(atmRec) == 327

def test_parseNMRModel():
    from MolKit import Read
    mols = Read("./Data/2EZNsmall.pdb")
    assert len(mols)==3, "Expecting 3 NMR models got %s"%len(mols)
    assert mols.chains.id == [' ', ' ', ' '], "Expecting chain id ' ' got %s"%mols.chains.id

def test_parseNMRModel_Conformations():
    from MolKit import Read
    mols = Read("./Data/2EZNsmall.pdb", modelsAs='conformations')
    assert len(mols)==1, "Expecting 1 NMR models got %s"%len(mols)
    assert len(mols[0].allAtoms[0]._coords)==3, "Expecting 3 NMR atom _coords got %s"%len(mols[0].allAtoms[0]._coords)
    assert mols.chains.id == [' '], "Expecting chain id ' ' got %s"%mols.chains.id

def test_parseTERRecord():
    from MolKit import Read

    # Protein with TER records specifying chain A and ChainB
    print "Testing 1bsr: protein with normal TER records"
    mol = Read("./Data/1bsr.pdb")[0]
    assert len(mol.chains)==2
    assert mol.chains.name == ['A', 'B']
    print mol.chains.name

    # Protein with TER records specifying, 1- several chains and 2-
    # Gaps in chain E.
    # If chain ID is the same on each side of a TER records
    # this doesn't create a new chain but a gap information.
#    print "Testing 1ACB.pdb: protein with TER specifying GAPS in chain E"
#    mol = Read("./Data/1ACB.pdb")[0]
#    assert len(mol.chains)==3, "Expecting 3 chains got %s"%len(mol.chains)
#    assert mol.chains.gaps == [[('LEU13', 'ILE16'),('TYR146', 'ALA149')], [], []],"Expecting two gaps ('LEU13', 'ILE16'), ('TYR146', 'ALA149') in first chain and got %s"%mol.chains.gaps
#    print mol.chains.name
#    print mol.chains.gaps

    # Protein with incorrect TER record...Only TER record no information
    # on the chain ID.
    # TER record doesn't have the right information
    print "Testing ind.pdb: protein with INCORRECT TER record"
    mol = Read("./Data/ind.pdb")[0]
    assert len(mol.chains)==1, "Expecting 1 chain I got %s %s"%(len(mol.chains),mol.chains.name)
    print mol.chains.name

    # Protein with No TER records at all and 1 chain
    print "Testing 1crn_Hs.pdbq: protein with no TER Record and 1 chain"
    mol = Read("Data/1crn_Hs.pdbq")[0]
    assert len(mol.chains)==1, "Expecting 1 chain got %s"%len(mol.chains)
    assert mol.chains[0].name == " ", "Expecting chain name ' ', got %s"%mol.chains[0].name
    print mol.chains.name

    # Protein with 1 TER records and 2 chains A and B
    print "Testing 1crn2Chains.pdb: protein with only 1 TER record and 2 chains"
    mol = Read("Data/1crn2Chains.pdb")[0]
    assert len(mol.chains)==2, "Expecting 1 chain got %s"%len(mol.chains)
    assert mol.chains.name == ['A','B'], "Expecting chain name ' ', got %s"%mol.chains.name
    print mol.chains.name
    
    # Protein with No TER records at all and 2 chains
    print "Testing 1crn2Chains.pdb: protein with only 1 TER record at the end of the ATOM"
    mol = Read("Data/1crn2ChainsNOTER.pdb")[0]
    assert len(mol.chains)==2, "Expecting 1 chain got %s"%len(mol.chains)
    assert mol.chains.name == ['A','B'], "Expecting chain name ' ', got %s"%mol.chains.name
    print mol.chains.name

def test_parsePQR():
    from MolKit.pdbParser import PQRParser
    parser = PQRParser("./Data/mead.pqr")
    mols = parser.parse()
    assert len(mols)==1, "Expecting 1 molecule mead got %s"%mols.name
    mol = mols[0]
    assert len(mol.chains)==1, "Expecting 1 chain got %s"%len(mol.chains)
    

    
harness = testplus.TestHarness( __name__,
                                funs=testplus.testcollect( globals()),
                                #funs=[test_parsegetRecords],
                                )

if __name__ == '__main__':
    print harness
    sys.exit( len( harness))

