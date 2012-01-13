#
# $Header: /opt/cvs/python/packages/share1.5/mglutil/TestUtil/Tests/Data/test_displayCommands.py,v 1.3 2009/03/12 20:52:33 vareille Exp $
#
# $Id: test_displayCommands.py,v 1.3 2009/03/12 20:52:33 vareille Exp $
#

import sys
from exceptions import ValueError, AssertionError
#from mglutil.regression import testplus
import unittest
"""
Test module for the Pmv module: displayCommands
"""
## ############################################################################
## ### BASE TEST CLASS FOR THE DISPLAY MODULE
## ############################################################################

class DisplayBaseTest(unittest.TestCase):

    def setUp(self):
        from Pmv.moleculeViewer import MoleculeViewer
        self.mv = MoleculeViewer(customizer = './.empty',
                                 logMode = 'overwrite', withShell=0)
        self.mv.loadCommand('fileCommands', 'readMolecule', 'Pmv')
        self.mv.loadCommand('deleteCommands','deleteMol', 'Pmv')
        self.mv.loadCommand("bondsCommands", "buildBondsByDistance", "Pmv")
        self.mv.setOnAddObjectCommands(['buildBondsByDistance',
                                        'displayLines'], log=0)
        self.mv.loadModule("interactiveCommands", 'Pmv')
        # Don't want to trap exceptions and errors...
        # the user pref is set to 1 by
        # default
        self.mv.setUserPreference(('trapExceptions', '0'), log = 0)
        self.mv.setUserPreference(('warningMsgFormat', 'printed'), log = 0)
        
        self.mv.loadModule('displayCommands', 'Pmv')
        
    def tearDown(self):
        self.mv.Exit(0)

## ############################################################################
## ### DISPLAYLINES TESTS
## ############################################################################

class DisplayLinesTest(DisplayBaseTest):
    def test_displaylines_emptyviewer(self):
        """
        Testing __call__ on empty viewer
        """
        # here we test if the command can be called with nodes arg only
        if self.mv.displayLines.flag & 1:
            self.mv.displayLines(self.mv.getSelection())

    def test_displaylines_defaultargs(self):
        """
        Testing __call__ on 7ins with default arguments not specified
        """
        # Reading the test molecule
        self.mv.readMolecule('7ins.pdb')
        assert len(self.mv.Mols) == 1 and self.mv.Mols[0].name == '7ins'
        # Test body
        self.mv.displayLines(self.mv.getSelection())
        self.assertEqual(self.mv.Mols[0].geomContainer.atoms['bonded'],
                         self.mv.Mols[0].allAtoms)

        
    def test_displaylines_negate(self):
        """
        Testing __call__ on 7ins with only = 0 and negate = 1
        """
        # Reading the test molecule
        self.mv.readMolecule('7ins.pdb')
        # Test body
        self.mv.displayLines(self.mv.getSelection(), negate=1)
        self.assertEqual(len(self.mv.Mols[0].geomContainer.atoms['bonded']),
                         0)

    def test_displaylines_only(self):
        """
        Testing __call__ on 7ins with only = 1 and negate = 0
        """
        # Reading the test molecule
        self.mv.readMolecule('7ins.pdb')
        # Test body
        self.mv.select(self.mv.Mols[0].chains[0])
        self.mv.displayLines(self.mv.getSelection(), only = 1)
        from MolKit.molecule import Atom
        # need to assert something about the lines geom as well.
        self.assertEqual(len(self.mv.Mols[0].geomContainer.atoms['bonded']),
                         len(self.mv.Mols[0].chains[0].findType(Atom)))

    def test_displaylines_bug33(self):
        """
        Test written for the bug #33 reported in MGLBUGZILLA by Daniel.
        only does not work when multiple molecule and selection only on
        one mol.
        """
        # Reading the test molecules
        self.mv.readMolecule('7ins.pdb')
        self.mv.readMolecule('1crn.pdb')

        # Test body
        self.mv.select(self.mv.Mols[0].chains[0].residues[1:15])
        nodes = self.mv.getSelection()
        self.mv.displayLines(nodes, only=1)
        self.assertEqual(self.mv.Mols[0].geomContainer.atoms['bonded'],
                         nodes.atoms.uniq())
        self.assertEqual(len(self.mv.Mols[1].geomContainer.atoms['bonded']),
                         0)
        



############################################################################
### TEST SHOWMOLECULE
############################################################################
class ShowMoleculeTest(DisplayBaseTest):

    def test_showmolecules_emptyviewer(self):
        """
        __call__ on empty viewer
        """
        from exceptions import ValueError
        if self.mv.showMolecules.flag & 1:
            self.mv.showMolecules(self.mv.getSelection())
        else:
            print "Cannot be called with empty selection"
    def test_showmolecules_defaultargs(self):
        """
        Testing __call__ with default argument on the 7ins.
        negate = 0
        """
        # Reading the test molecule
        self.mv.readMolecule('7ins.pdb')
        # Test body
        self.mv.showMolecules(['7ins',])
        self.assertEqual(self.mv.Mols[0].geomContainer.geoms['master'].visible,
                         1)
        # Deleting the test molecule.
        # self.mv.deleteMol('7ins')

    def test_showmolecules_withargs(self):
        """
        Testing __call__ with the following arguments:
        molName = ['7ins',]
        negate = 1
        """
        # Reading the test molecule
        self.mv.readMolecule('7ins.pdb')
        # Test body
        self.mv.showMolecules(['7ins',], negate=1)
        self.assertEqual(self.mv.Mols[0].geomContainer.geoms['master'].visible,
                         0)
        # Deleting the test molecule.
        # self.mv.deleteMol('7ins')

    def test_showmolecules_withargs2(self):
        """
        Testing __call__ with the following arguments:
        molName = '7ins' bad input the molName should a list of molecule name.
        negate = 1

        """
        # Reading the test molecule
        self.mv.readMolecule('7ins.pdb')
        # Test body
        oldvisible = self.mv.Mols[0].geomContainer.geoms['master'].visible
        self.mv.showMolecules('7ins', negate=1)
        self.assertEqual(self.mv.Mols[0].geomContainer.geoms['master'].visible,
                         oldvisible)
        # Deleting the test molecule.
        # self.mv.deleteMol('7ins')

############################################################################
### TEST DISPLAYCPK
############################################################################
class DisplayCPKTest(DisplayBaseTest):

    def test_displaycpk_emptyviewer(self):
        """
        __call__ on empty viewer
        """

        # here we test if the command can be called with nodes arg only
        if self.mv.displayCPK.flag & 1   :
            self.mv.displayCPK(self.mv.getSelection())
        else:
            raise ValueError("WARNING: self.mv.displayCPK cannot be called with only self.mv.getSelection()")

    def test_displaycpk_bug33(self):
        """
        Test written for the bug #33 reported in MGLBUGZILLA by Daniel. only does
        not work when multiple molecule and selection only on one mol.
        """
        # Reading the test molecules
        self.mv.readMolecule('7ins.pdb')
        self.mv.readMolecule('1crn.pdb')
        self.mv.displayCPK(self.mv.getSelection())

        # Test body
        self.mv.select(self.mv.Mols[0].chains[0].residues[1:15])
        nodes = self.mv.getSelection()
        self.mv.displayCPK(nodes, only=1)
        self.assertEqual(self.mv.Mols[0].geomContainer.atoms['cpk'].sort(),
                         nodes.atoms.uniq().sort())
        self.assertEqual(len(self.mv.Mols[1].geomContainer.atoms['cpk']),0)

        # delete molecules.
##         self.mv.deleteMol('7ins')
##         self.mv.deleteMol('1crn')

    def test_dispalyCPKWithVariousScaleFactors(self):
        self.mv.loadCommand("selectionCommands","selectFromString", "Pmv")
        self.mv.loadCommand("selectionCommands","clearSelection",  "Pmv")
        self.mv.readMolecule('1crn.pdb')
        self.mv.selectFromString('', '', '1-10', 'CA', 1)
        self.mv.displayCPK(self.mv.getSelection(), scaleFactor=2.0, quality=10)
        self.mv.clearSelection()
        self.mv.selectFromString('', '', '11-20', 'CA', 1)
        self.mv.displayCPK(self.mv.getSelection(), scaleFactor=.5, quality=10)
        self.mv.displayCPK('1crn', negate=1)
        self.mv.displayCPK('1crn')
        self.mv.deleteMol('1crn')

    def test_dispalyCPKAssignRadii_1(self):
        self.mv.readMolecule('1crn.pdb')
        # no radius assigned
        nodes = self.mv.getSelection()
        wrad = filter(lambda x: hasattr(x, 'radius'), self.mv.allAtoms)
        self.assertEqual(len(wrad), 0)
        self.mv.displayCPK(nodes)
        wrad = filter(lambda x: hasattr(x, 'radius'), self.mv.allAtoms)
        self.assertEqual(len(wrad),len(self.mv.allAtoms))
        self.assertEqual(self.mv.Mols[0].unitedRadii,1)
        self.mv.displayCPK(nodes, negate=1)

        self.mv.assignAtomsRadii("1crn", united=0, overwrite=1)
        self.mv.displayCPK(nodes)
        self.assertEqual(self.mv.Mols[0].unitedRadii,0)

##         self.mv.deleteMol("1crn")
##         assert len(self.mv.Mols)==0

############################################################################
### TEST DISPLAYSTICKSANDBALLS
############################################################################
class DisplaySticksAndBallsTest(DisplayBaseTest):

    def test_displaysticksballs_emptyviewer(self):
        """
        __call__ on empty viewer
        """
        # here we test if the command can be called with nodes arg only
        if self.mv.displaySticksAndBalls.flag & 1:
            self.mv.displaySticksAndBalls(self.mv.getSelection())
        else:
            raise ValueError("WARNING: self.mv.displaySticksAndBalls cannot be called with only self.mv.getSelection()")

    def test_displaysticksballs_bug33(self):
        """
        Test written for the bug #33 reported in MGLBUGZILLA by Daniel.
        only does not work when multiple molecule and selection only on
        one mol.
        """
        # Reading the test molecules
        self.mv.readMolecule('7ins.pdb')
        self.mv.readMolecule('1crn.pdb')
        self.mv.displaySticksAndBalls(self.mv.getSelection())
        # not bondsBOTest body
        self.mv.select(self.mv.Mols[0].chains[0].residues[1:15])
        nodes = self.mv.getSelection()
        self.mv.displaySticksAndBalls(nodes, only=1)
        self.assertEqual(self.mv.Mols[0].geomContainer.atoms['sticks'],
                         nodes.atoms.uniq())
        self.assertEqual(self.mv.Mols[0].geomContainer.atoms['balls'],
                         nodes.atoms.uniq())

        self.assertEqual(len(self.mv.Mols[1].geomContainer.atoms['sticks']), 0)
        self.assertEqual(len(self.mv.Mols[1].geomContainer.atoms['balls']),0)

##         # delete molecules.
##         self.mv.deleteMol('7ins')
##         self.mv.deleteMol('1crn')
##         assert len(self.mv.Mols) == 0



if __name__ == '__main__':
    unittest.main()
