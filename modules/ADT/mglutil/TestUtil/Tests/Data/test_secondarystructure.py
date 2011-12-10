#############################################################################
#
# Author: Sophie I. COON, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################
#
# $Header: /opt/cvs/python/packages/share1.5/mglutil/TestUtil/Tests/Data/test_secondarystructure.py,v 1.3 2009/03/12 20:52:33 vareille Exp $
#
# $Id: test_secondarystructure.py,v 1.3 2009/03/12 20:52:33 vareille Exp $
#
import sys
"""
This module implements a set of functions to test the commands of the
secondarystructurecommands module
"""
mv= None
setUp='startMoleculeViewer'
tearDown = 'quitMoleculeViewer'

# We should create a test protein directory ....
def startMoleculeViewer():
    global mv
    from Pmv.moleculeViewer import MoleculeViewer
    mv = MoleculeViewer(customizer = './.empty', logMode = 'overwrite', withShell=0)
    mv.setUserPreference(('trapExceptions', '0'), log = 0)
    mv.setUserPreference(('warningMsgFormat', 'printed'), log = 0)
    mv.loadCommand('fileCommands', 'readMolecule', 'Pmv')
    mv.loadCommand('deleteCommands','deleteMol', 'Pmv')
    mv.loadCommand("bondsCommands", "buildBondsByDistance", "Pmv")
    mv.setOnAddObjectCommands(['buildBondsByDistance','displayLines'], log=0)
    mv.loadModule("interactiveCommands", 'Pmv')
    mv.loadModule('secondaryStructureCommands', 'Pmv')

def quitMoleculeViewer():
    mv.Exit(0)


############################################################
# ComputeSecondaryStructureCommand Tests                   #
############################################################
def test_computeSecondaryStructure_emptyViewer():
    """
    Test if the computeSecondaryStructure behaves properly when no molecule has
    been loaded in the viewer and with the default values
    """
    mv.computeSecondaryStructure(mv.getSelection())

def test_computeSecondaryStructure_Normal():
    """
    Test the normal behavior of the computeSecondaryStructure on 1crn
    """
    # When computing the secondary structure of a 1crn we expect:
    #  - a secondary structure set containing 10 secondary structure elements
    #  - the selectionCommands.sets__ dictionnary contains 10 extra keys
    #  - a new attribute of the mol has been created called hasSS its value
    # should be ['From File']
    #  - resWithSS ?
    mv.readMolecule("1crn.pdb")
    mv.computeSecondaryStructure("1crn", molModes = {'1crn':'From File'})
    mol = mv.Mols[0]
    assert hasattr(mol, 'hasSS') and mol.hasSS == ['From File']
    assert hasattr(mol.chains[0], 'secondarystructureset') and \
           len(mol.chains[0].secondarystructureset)==10
    assert not hasattr(mol,'_ExtrudeSecondaryStructureCommand__hasSSGeom') or \
           mol._ExtrudeSecondaryStructureCommand__hasSSGeom == 0
    mv.deleteMol("1crn")
    assert len(mv.Mols)==0
    
def test_colorSecondaryStructure_Normal():
    """
    Test the normal behavior of the computeSecondaryStructure on 1crn
    """
    # When computing the secondary structure of a 1crn we expect:
    #  - a secondary structure set containing 10 secondary structure elements
    #  - the selectionCommands.sets__ dictionnary contains 10 extra keys
    #  - a new attribute of the mol has been created called hasSS its value
    # should be ['From File']
    #  - resWithSS ?
    mv.readMolecule("1crn.pdb")
    mv.computeSecondaryStructure("1crn", molModes = {'1crn':'From File'})
    mv.colorBySecondaryStructure("1crn")
    mol = mv.Mols[0]
    #assert hasattr(mol, 'hasSS') and mol.hasSS == ['From File']
    #assert hasattr(mol.chains[0], 'secondarystructureset') and \
    #       len(mol.chains[0].secondarystructureset)==10
    #assert not hasattr(mol,'_ExtrudeSecondaryStructureCommand__hasSSGeom') or \
    #       mol._ExtrudeSecondaryStructureCommand__hasSSGeom == 0
    mv.deleteMol("1crn")
    assert len(mv.Mols)==0
    

def test_computeSecondaryStructure_nosheet2D():
    """
    Test the computeSecondaryStructure on a molecule (1crnnosheet2D) created
    using 1crn.pdb (res1-res29) removing all the O atoms, but leaving the
    Secondary Structure information in the file.
    We expect to:
      - create the proper secondarystructureset containing 5 SS elements
      - but no sheet2D
    """
    mv.readMolecule('1crnnosheet2D.pdb')
    mv.select('1crnnosheet2D')
    mv.computeSecondaryStructure(mv.getSelection())
    chain = mv.Mols[0].chains[0]
    mol = mv.Mols[0]
    assert mol.hasSS == ['From File']
    assert hasattr(chain, 'secondarystructureset')
    assert len(chain.secondarystructureset) == 5
    mv.deleteMol('1crnnosheet2D')
    assert len(mv.Mols) == 0

def test_computeSecondaryStructure_fileThenStrideOnMol():
    """
    ComputeSecondaryStructure on 1crn using the info in the file then
    using stride."""
    mv.readMolecule("1crn.pdb")
    mv.select("1crn")
    mv.computeSecondaryStructure(mv.getSelection(), {'1crn':'From File'})
    mol = mv.Mols[0]
    assert mol.hasSS == ['From File']
    assert hasattr(mol.chains[0], 'secondarystructureset') and \
           len(mol.chains[0].secondarystructureset)==10
    mv.computeSecondaryStructure(mv.getSelection(), {'1crn':'From Stride'})
    assert mol.hasSS == ['From Stride']
    assert hasattr(mol.chains[0], 'secondarystructureset') and \
           len(mol.chains[0].secondarystructureset)==11

    mv.deleteMol('1crn')
    assert len(mv.Mols)==0
    
def test_computeSecondaryStructure_twochains():
    """
    Test computeSecondaryStructure on a molecule protease.pdb with two chains
    """
    mv.readMolecule('protease.pdb')
    mv.select('protease')
    mv.computeSecondaryStructure(mv.getSelection())
    mol = mv.Mols[0]
    assert mol.hasSS == ['From File']
    assert hasattr(mol.chains[0], 'secondarystructureset')
    assert hasattr(mol.chains[1], 'secondarystructureset')
    mv.deleteMol('protease')
    assert len(mv.Mols) == 0
    
def test_computeSecondaryStructure_nofileinfo():
    """
    Test computeSecondaryStructure on 7ins.pdb with no information in the file
    using stride by default
    """
    from MolKit.protein import SecondaryStructureSet
    mv.readMolecule('7ins.pdb')
    mv.select('7ins')
    mv.computeSecondaryStructure(mv.getSelection())
    mol = mv.Mols[0]
    assert mol.hasSS == ['From Stride']
    c0 = mol.chains[0]
    assert (isinstance(c0.secondarystructureset, SecondaryStructureSet) and \
            len(c0.secondarystructureset)==6)
    c1 = mol.chains[1]
    assert (isinstance(c1.secondarystructureset, SecondaryStructureSet) and \
            len(c1.secondarystructureset)==3)
    c2 = mol.chains[2]
    assert (isinstance(c2.secondarystructureset, SecondaryStructureSet) and \
            len(c2.secondarystructureset)==4)
    c3 = mol.chains[3]
    assert (isinstance(c3.secondarystructureset, SecondaryStructureSet) and \
            len(c3.secondarystructureset)==3)
    c4 = mol.chains[4]
    assert (isinstance(c4.secondarystructureset, SecondaryStructureSet) and \
            len(c4.secondarystructureset)==4)
    c5 = mol.chains[5]
    assert (isinstance(c5.secondarystructureset, SecondaryStructureSet) and \
            len(c5.secondarystructureset)==2)
    c6 = mol.chains[6]
    assert not hasattr(c6, 'secondarystructureset')
    c7 = mol.chains[7]
    assert not hasattr(c7, 'secondarystructureset')

    mv.deleteMol('7ins')
    assert len(mv.Mols) == 0

def test_computeSecondaryStructure_fileOnMolAndStrideOnMol():
    """
    Test computeSecondaryStructure command on two proteins protease.pdb with
    File information and 7ins with no File information.
    The molMode are specified protease : From File and 7ins From Stride
    """
    from MolKit.protein import SecondaryStructureSet
    mv.readMolecule('protease.pdb')
    mv.readMolecule('7ins.pdb')
    mv.computeSecondaryStructure(mv.getSelection(),
                                 molModes={'protease':'From File',
                                          '7ins':'From Stride'})
    prot = mv.Mols[0]
    ins = mv.Mols[1]
    assert prot.hasSS == ['From File']
    assert isinstance(prot.chains[0].secondarystructureset,
                      SecondaryStructureSet)
    
    assert isinstance(prot.chains[1].secondarystructureset,
                      SecondaryStructureSet)

    assert ins.hasSS == ['From Stride']
    c0 = ins.chains[0]
    assert (isinstance(c0.secondarystructureset, SecondaryStructureSet) and \
            len(c0.secondarystructureset)==6)
    c1 = ins.chains[1]
    assert (isinstance(c1.secondarystructureset, SecondaryStructureSet) and \
            len(c1.secondarystructureset)==3)
    c2 = ins.chains[2]
    assert (isinstance(c2.secondarystructureset, SecondaryStructureSet) and \
            len(c2.secondarystructureset)==4)
    c3 = ins.chains[3]
    assert (isinstance(c3.secondarystructureset, SecondaryStructureSet) and \
            len(c3.secondarystructureset)==3)
    c4 = ins.chains[4]
    assert (isinstance(c4.secondarystructureset, SecondaryStructureSet) and \
            len(c4.secondarystructureset)==4)
    c5 = ins.chains[5]
    assert (isinstance(c5.secondarystructureset, SecondaryStructureSet) and \
            len(c5.secondarystructureset)==2)
    c6 = ins.chains[6]
    assert not hasattr(c6, 'secondarystructureset')
    c7 = ins.chains[7]
    assert not hasattr(c7, 'secondarystructureset')
    mv.deleteMol('7ins')
    mv.deleteMol('protease')
    assert len(mv.Mols) == 0

def test_computeSecondaryStructure_forceFileWhenNoInfo():
    """
    Test computeSecondaryStructure command on 7ins.pdb which doesn't have
    information on the secondarystructure in its file. Forces to compute
    from the file information. What do we expect ?
    """
    from MolKit.protein import SecondaryStructureSet
    mv.readMolecule('7ins.pdb')
    mv.select('7ins')
    mv.computeSecondaryStructure(mv.getSelection(),
                                 molModes={'7ins':'From File'})
    ins = mv.Mols[0]
    assert ins.hasSS == []
    
    c0 = ins.chains[0]
    assert not hasattr(c0, 'secondarystructureset')
    c1 = ins.chains[1]
    assert not hasattr(c1, 'secondarystructureset')
    c2 = ins.chains[2]
    assert not hasattr(c2, 'secondarystructureset')
    c3 = ins.chains[3]
    assert not hasattr(c3, 'secondarystructureset')
    c4 = ins.chains[4]
    assert not hasattr(c4, 'secondarystructureset')
    c5 = ins.chains[5]
    assert not hasattr(c5, 'secondarystructureset')
    c6 = ins.chains[6]
    assert not hasattr(c6, 'secondarystructureset')
    c7 = ins.chains[7]
    assert not hasattr(c7, 'secondarystructureset')

    mv.deleteMol('7ins')
    assert len(mv.Mols) == 0
                                 

############################################################
# ExtrudeSecondaryStructureCommand Tests                   #
############################################################
def test_extrudeSecondaryStructure_emptyViewer():
    """ Test if extrudeSecondaryStructure behaves properly when no molecule has
    been loaded in the viewer and with the default values"""
    mv.extrudeSecondaryStructure(mv.getSelection())

def test_extrudeSecondaryStructure_noDheet2D():
    """
    Testing extrudeSecondaryStructure for a protein 1crnnosheet2D with 1 chain
    having a secondarystructureset holding 5 secondary structure elements
    but no sheet2D because no residues have an O.
    """
    mv.readMolecule('1crnnosheet2D.pdb')
    mv.computeSecondaryStructure(mv.getSelection())
    chain = mv.Mols[0].chains[0]
    mv.extrudeSecondaryStructure(mv.getSelection(), display=1)
    assert hasattr(chain, 'sheet2D')
    assert chain.sheet2D['ssSheet2D'] is None
    mv.deleteMol('1crnnosheet2D')
    assert len(mv.Mols)==0

def test_extrudeSecondaryStructure_defaultParam():
    """
    Function to test the extrudeSecondaryStructure command for 1crn with 1
    chain with the default parameters
    """
    mv.readMolecule('1crn.pdb')
    mv.computeSecondaryStructure(mv.getSelection())
    chain = mv.Mols[0].chains[0]
    mv.extrudeSecondaryStructure(mv.getSelection(), display=1)
    
    mv.deleteMol('1crn')
    assert len(mv.Mols)==0

   
################################################################
# DisplayExtrudedSSCommand Tests                               #
################################################################
def test_displaySecondaryStructure_emptyViewer():
    """ Test if displaySecondaryStructure behaves properly when no
    molecule has been loaded in the viewer and with the default values"""
    mv.displayExtrudedSS(mv.getSelection())

def test_displaySecondaryStructure_beforeCompute():
    """
    Test the display secondarystructurecommands before the
    computing and extruding the secondary structure information
    """
    mv.readMolecule('1crn.pdb')
    mv.displayExtrudedSS(mv.getSelection())
    mv.deleteMol('1crn')
    assert len(mv.Mols) == 0
################################################################
# OTHER secondaryStructure  Tests                              #
################################################################

def test_secondaryStructure_fileThenStrideOnMol():
    """
    ComputeSecondaryStructure on 1crn using the info in the file then
    using stride."""
    mv.readMolecule("1crn.pdb")
    mv.select("1crn")
    mv.computeSecondaryStructure(mv.getSelection(), {'1crn':'From File'})
    mv.extrudeSecondaryStructure(mv.getSelection())
    mol = mv.Mols[0]
    assert mol.hasSS == ['From File']
    assert hasattr(mol.chains[0], 'secondarystructureset') and \
           len(mol.chains[0].secondarystructureset)==10

    mv.computeSecondaryStructure(mv.getSelection(), {'1crn':'From Stride'})
    mv.extrudeSecondaryStructure(mv.getSelection())
    assert mol.hasSS == ['From Stride']
    assert hasattr(mol.chains[0], 'secondarystructureset') and \
           len(mol.chains[0].secondarystructureset)==11

    mv.deleteMol('1crn')
    assert len(mv.Mols)==0

def test_secondaryStructure_cleanBeforeExtrude():
    """ This test makes sure that everything is cleaned by the clean method
    The clean method is called when computing the SS using stride after using
    the info in the file.
    """
    mv.readMolecule('1crn.pdb')
    mv.select('1crn')
    mv.computeSecondaryStructure(mv.getSelection(), {'1crn':'From File'})
    mol = mv.Mols[0]
    mv.computeSecondaryStructure.clean(mol)
    mv.deleteMol('1crn')
    assert len(mv.Mols) == 0
    
############################################################
# RibbonCommand Tests                                      #
############################################################
def test_ribbon_emptyViewer():
    """ Test if the ribbon behaves properly when no molecule has
    been loaded in the viewer and with the default values"""
    mv.ribbon(mv.getSelection())

def test_ribbon_2tbv():
    """ Test if ribbon on the 2tbv with the default param"""
    mv.readMolecule('./2tbv.pdb')
    assert len(mv.Mols) == 1 and mv.Mols[0].name == '2tbv'
    mv.select('2tbv')
    mv.ribbon(mv.getSelection())
    mv.deleteMol('2tbv')


############################################################
# OTHER BUGS                                               #
############################################################
def test_DanielBug():
    mv.readMolecule("./fx.pdb")
    mv.loadCommand("selectionCommands", ['selectFromString',], 'Pmv')
    mv.loadModule("editCommands",'Pmv')
    mv.loadModule("deleteCommands", 'Pmv')

    mv.selectFromString('','','','H*',1)
    mv.deleteAtomSet(mv.getSelection())
    mv.add_hGC("fx:::", polarOnly = 1, renumber = 1,
               method = 'noBondOrder', log = 0)
    mv.computeS


## if __name__ == '__main__':
##     testplus.chdir()
##     print harness
##     sys.exit( len( harness))
