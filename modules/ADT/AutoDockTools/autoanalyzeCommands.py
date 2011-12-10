## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#############################################################################
#
# Author: Ruth HUEY, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################


# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/autoanalyzeCommands.py,v 1.260.2.7 2011/05/09 20:45:22 rhuey Exp $
#
# $Id: autoanalyzeCommands.py,v 1.260.2.7 2011/05/09 20:45:22 rhuey Exp $
#
#
#
#
#
#
#

"""
This Module facilitates analyzing results of autodock jobs. 

    * The first step is 'Read Docking Log'  The selected file is parsed 
    which sets self.docked to a new Docking instance.  The Docking class 
    has attributes:

        o dlgParser

            x 'dlg': full pathname of dlg

        o dpo 

        o ch:a conformation handler.

            x 'clusterNum': 

            x 'clusterList': 

            x 'modelList': a list of docked conformations 

        o macroFile: the Macromolecule file used

        o 'macro': filename of macromolecule (eg '1hvrCorr.pdbqt')

        o 'macroStem': name of macromolecule up to last '.' (eg '1hvrCorr')

        o ligand:  the original ligand 

        o output: lines containing summary of docking

    The new Docking is also entered in the dictionary 'dockings' as a separate item 
    whose key is the file and whose value is the Docking.


After the selected docking log file is parsed, the user can:

    * select a displayed docked conformation using the 'Choose A Docked Conformation' menubutton.  This opens a DockingChooser widget which is a ListChooser allowing selection either in the widget or in the viewer of any of the displayed docking. Information about each docked conformation is displayed in the information window of the DockingChooser as different entries are high-lighted.  

    * display the macromolecule via the "Show Macromolecule" menubutton.  This menubutton is linked to a file browsers in case the macromolecule whose name is parsed from the docking log file is not in the current directory. (FIX THIS: what if the macromolecule is in a different directory but there is a molecule with the same name here???). The user can change the visibility, sampling, isovalue, renderMode and visibility of bounding box  for each of  the displayed grids

    * display the autogrids used in the docking via the "Show Grids Used For Calc" menubutton.  This menubutton is linked to a ListChooser which lets the user select whether to load all or some of the grids. The user can interactively change the visibility of each grid's isosurface, its sampling value, its isovalue, its rendermode (LINE or FILL) and the visibility of its bounding box. 

    * The user is able to visualize extra grid maps using the "Show Grid" button. 

    * If the current docking has clusters, the user is able to visualize a results histogram for it with 'Show Histogram'. The histogram can be printed.

    * Result Summaries for docking(s) can be viewed, edited and saved with 'Get Output'

    * Dockings can be deleted via 'Delete Docking Log'

"""
import types, Tkinter, commands, os, copy, glob, popen2, time
from string import find, join, replace, split, rfind, strip
import numpy.oldnumeric as Numeric, math
import re
import Pmw

from DejaVu.Box import Box
from DejaVu.Geom import Geom
from DejaVu.Spheres import Spheres
from DejaVu.Cylinders import Cylinders
from DejaVu.colorTool import Map, RGBRamp, array2DToImage
from DejaVu.bitPatterns import patternList
from DejaVu import colorTool

from ViewerFramework.VFCommand import CommandGUI
from Pmv.moleculeViewer import EditAtomsEvent
#from ViewerFramework.VF import ModificationEvent

from mglutil.gui.InputForm.Tk.gui import InputFormDescr
from mglutil.gui.BasicWidgets.Tk.colorWidgets import ColorChooser
from mglutil.gui.BasicWidgets.Tk.thumbwheel import ThumbWheel

#from mglutil.gui.InputForm.Tk.gui import InputFormDescr, ListChooser
from mglutil.gui.BasicWidgets.Tk.customizedWidgets import ExtendedSliderWidget, ListChooser
from mglutil.util.callback import CallBackFunction
from mglutil.math.statetocoords import StateToCoords
from mglutil.math.rmsd import RMSDCalculator

from MolKit import molecule, protein, Read
from MolKit.tree import TreeNode, TreeNodeSet
from MolKit.molecule import  MoleculeSet, Molecule, AtomSet, Atom, Bond, HydrogenBond, HydrogenBondSet
from MolKit.protein import ProteinSet, Protein, ResidueSet, Residue
from MolKit.pdbParser import PdbParser, PdbqParser, PdbqsParser, PdbqtParser
from MolKit.stringSelector import CompoundStringSelector
from AutoDockTools.DockingParameters import DockingParameters
from AutoDockTools.Conformation import AutodockState, ConformationHandler
from AutoDockTools.ConfPlayer import ConformationPlayer, PopulationPlayer
from AutoDockTools.Docking import Docking, FoxResultProcessor
from AutoDockTools.DlgParser import DlgParser
from AutoDockTools.EntropiaParser import EntropiaParser
from AutoDockTools.TrajParser import TrajParser
from AutoDockTools.EpdbParser import EpdbParser
from AutoDockTools.histogram import HistogramRI
from AutoDockTools.cluster import Clusterer
from AutoDockTools import interactiveHistogramGraph
from AutoDockTools.InteractionDetector import InteractionDetector

from MolKit.torTree import TorTree


from Pmv.mvCommand import MVCommand
from Pmv.guiTools import MoleculeChooser, BarButton, Kill
#from Pmv.Grid import AutoGrid, AutoGridSurfaceGui

from SimpleDialog import SimpleDialog

#these are the texts on menubuttons, menu entries etc:
menuText = {}
#menuText['AnalyzeMB'] = ' Analyze Docking Results '
menuText['AnalyzeMB'] = 'Analyze'
#A Dockings
#1-6, 5 is a separator
menuText['DockingLogMB'] = 'Dockings'
menuText['readDLG'] = 'Open...'
menuText['readVinaResult'] = 'Open AutoDock vina result...'
menuText['readVSResult'] = 'Open AutoDock Virtual Screening Result ligand...'
menuText['deleteDLG'] = 'Clear...'
menuText['selectDLG'] = 'Select...'
menuText['readDirDLG'] = 'Open All...'
menuText['seeSpots'] = 'Show as Spheres...'
menuText['showBindingSite'] = 'Show Interactions'

#B Macromolecule
#1-2
menuText['MoleculesMB'] = 'Macromolecule'
menuText['readMacro'] = 'Open...'
menuText['chooseMacro'] = 'Choose...'

#C 
#null

#D Grids
#1-3
menuText['GridsMB'] = 'Grids'
menuText['showGridsMB'] = 'Open...'
menuText['addGridMB'] = 'Open Other...'
menuText['epdbMol'] = 'Calculate Energy...'

#E Conformations
#1-3
menuText['StatesMB'] = 'Conformations'
menuText['showStatesMB'] = 'Play...'
menuText['showStatesByEnergyMB'] = 'Play, ranked by energy...'
menuText['readStatesMB'] = 'Open...'
menuText['writeResultMB']  = 'Save...'
menuText['chooseConfMB'] = 'Load...'
menuText['showPopulationMB'] = 'View initial population...'

#F Clusterings
menuText['ClusteringMB'] = 'Clusterings'
menuText['showStatesCLUSTERINGMB'] = 'Show...'
menuText['readStatesCLUSTERINGMB'] = 'Read...'
menuText['makeCLUSTERINGMB'] = 'Recluster...'
menuText['makeSubsetCLUSTERINGMB'] = 'Recluster on Subset...'
menuText['writeClusteringMB']  = 'Save...'

#ZULU Results
menuText['ResultsMB'] = 'Other Results'
menuText['showHistogramMB'] = 'Plot Histogram...'
menuText['showStatesHISTOGRAMMB'] = 'Show Histogram'
menuText['showChartMB'] = 'Tabulate Energies...'
menuText['getOutputMB'] = 'Extract Histogram...'
menuText['writeVSResult'] = 'Write AutoDock Virtual Screening Result ligand...'

#used to determine which player widget to show

def checkHasInitializedDockings(vf):
    if not hasattr(vf, 'dockings'):
        vf.dockings = {}
        vf.docked = None


def hideShowHide(vf, mol):
    #nb mol is a molecule
    if hasattr(vf,'showMolecules'):
        cbvar=vf.showMolecules.molList[mol.name]
        cbvar.set(0)
        vf.showMolecules.checkButton_cb(mol.name)
    else:
        vf.displayLines(mol,only=0,negate=1)


def toggleShowHide(vf,molname):
    #if vf.commands.has_key('showMolecules'):
    if hasattr(vf,'showMolecules'):
        cbvar=vf.showMolecules.molList[molname]
        if cbvar.get(): cbvar.set(0)
        else:cbvar.set(1)
        vf.showMolecules.checkButton_cb(molname)
    else:
        vf.displayLines(mol,only=0,negate=1)


def checkNameStr(molname):
    spChar=['?','*','.','$','#',':','-',',']
    for item in spChar:
        if find(molname,item)>-1:
            molname = replace(molname,item,'_')
    return molname


class ADChooseMacro(MVCommand):
    """This class is used to choose the macromolecule from molecules already present in the viewer
    \nPackage : AutoDockTools
    \nModule  : autoanalyzeCommands
    \nClass   : ADChooseMacro
    \nCommand : ADanalyze_chooseMacromolecule
    \nSynopsis:\n    
        None<---ADanalyze_chooseMacromolecule(macroMol)
    \nRequired Arguments:\n    
        macroMol --- molecule present in viewer
    """

            
    def chooseMolecule_cb(self, event = None):
        """
        invoked from guiCallback"""
        mol = self.chooser.getMolSet()
        self.chooser.ipf.form.withdraw()
        if mol:
            self.doitWrapper(mol, log=1, redraw=0)
        else:
            return "ERROR"


    def __call__(self, mol, **kw):
        mols = self.vf.expandNodes(mol)
        if len(mols)==0:
            return 'ERROR'
        #need to construct a list of the current molecules connected to dockings
        self.buildMolList()
        apply(self.doit, (mols[0],), kw)


    def __call__(self, macroMol=None, **kw):
        """None<-ADanalyze_chooseMacromolecule(macroMol)
        \nmacroMol: molecule present in viewer
        """
        mols = self.vf.expandNodes(macroMol)
        if len(mols)==0:
            return 'ERROR'
        macroMol = mols[0]
        if not self.vf.docked:
            print "no current docking"
            return 'ERROR'
        #if hasattr(self.vf.docked, 'macroMol'):
        #    self.vf.warningMsg('current docking already has a macromolecule')
        #    return 'ERROR'
        apply(self.doitWrapper,(macroMol,), kw)


    def doit(self, macroMol):
        if not self.vf.docked:
            print "no current docking"
            return 'ERROR'
        self.vf.docked.macroMol = macroMol


    def guiCallback(self):
        """called each time the 'choose macromolecule' button is pressed"""
        if not hasattr(self.vf, 'docked'):
            msg = "no current docking: you must read a docking log first"
            self.vf.warningMsg(msg)
            return 'ERROR'
        if not len(self.vf.Mols):
            msg = 'no molecules in viewer'
            self.vf.warningMsg(msg)
            return 'ERROR'
        if hasattr(self.vf.docked, 'macroMol'):
            t = 'Current docking already has a macromolecule:'+ self.vf.docked.macroMol.name+"\nDo you want to replace it?"
            d = SimpleDialog(self.vf.GUI.ROOT, text = t,
            buttons = ["Yes", "No", "Cancel"], default = 0, 
                    title = "Replace current macromolecule? ")
            replace = d.go()
            if replace>0:
                return 'ERROR'
        self.chooser = MoleculeChooser(self.vf, mode='single', 
                        title='Choose Macromolecule for Current Docking')
        self.chooser.ipf.append({'name':'Select Button',
                                 'widgetType':Tkinter.Button,
                                 'text':'Select Macromolecule',
                                 'wcfg':{'bd':6},
                                 'gridcfg':{'sticky':Tkinter.E+Tkinter.W},
                                 'command': self.chooseMolecule_cb})
        self.form = self.chooser.go(modal=0, blocking=0)
        lb = self.chooser.ipf.entryByName['Molecule']['widget'].lb
        lb.bind("<Double-Button-1>",self.chooseMolecule_cb)


ADChooseMacroGUI=CommandGUI()
ADChooseMacroGUI.addMenuCommand('AutoToolsBar', menuText['AnalyzeMB'], 
        menuText['chooseMacro'], cascadeName = menuText['MoleculesMB'])


        
class ADReadVSResult(MVCommand):
    """This class is used to load a pdbqt file containing a docked VS result
    \nPackage : AutoDockTools
    \nModule  : autoanalyzeCommands
    \nClass   : ADReadVSResult
    \nCommand : ADanalyze_readVSResult
    \nSynopsis:\n 
        None<-ADanalyze_readVSResult(resultFile=None)
    \nArguments:\n   
        resultFile: pdbqt file from Virtual Screening AutoDock result including information on
                        clustering, hydrogen-bonds to receptor, close-contacts ...
    """


    def onAddCmdToViewer(self):
        checkHasInitializedDockings(self.vf)
        self.first = 1
        if self.vf.hasGui:
            for item in ['colorByMolecules','colorByAtomType', 'color']:
                if not hasattr(self.vf, item):
                    self.vf.loadCommand('colorCommands', item, 'Pmv')
            for item in ['displayLines','showMolecules']:
                if not hasattr(self.vf, item):
                    self.vf.loadCommand('displayCommands', item, 'Pmv')
        if 'VSCloseContacts_displayMode' not in self.vf.userpref.keys():
            doc = """Choose display mode for atoms in close contact in a VirtualScreening result:\nasSelection, asCPK,  asSticks or None\n. Default value is 'asSelection'"""
            self.vf.userpref.add('VSCloseContacts_displayMode', 'asSelection', 
                    ['asSelection', 'asCPK', 'asSticks','None'],
                    callbackFunc = [self.set_VSCloseContacts_displayMode], doc=doc ,
                    category='AutoDockTools')
            self.vf.userpref['VSCloseContacts_displayMode']['value'] = 'None'
        if 'VSHydrogenBond_displayMode' not in self.vf.userpref.keys():
            doc = """Choose display mode for hydrogen bonds in a VirtualScreening result:\nasSpheres, asCylinders  or None\n. Default value is 'asSpheres'"""
            self.vf.userpref.add('VSHydrogenBond_displayMode', 'asSpheres', 
                    ['asSpheres', 'asCylinders', 'asLines', 'None'],
                    callbackFunc = [self.set_VSHydrogenBond_displayMode], doc=doc ,
                    category='AutoDockTools')
            self.vf.userpref['VSHydrogenBond_displayMode']['value'] = 'None'
        if 'VSPiInteractions_displayMode' not in self.vf.userpref.keys():
            doc = """Choose display mode for atoms with PiInteractions in a VirtualScreening result:\nasCylindersAndCones, asSelection or None\n. Default value is 'asCylindersAndCones'"""
            self.vf.userpref.add('VSPiInteractions_displayMode', 'asCylindersAndCones', 
                    ['asCylindersAndCones', 'asSelection', 'None'],
                    callbackFunc = [self.set_VSPiInteractions_displayMode], doc=doc ,
                    category='AutoDockTools')
            self.vf.userpref['VSPiInteractions_displayMode']['value'] = 'None'
        if 'VS_ShowClusteringHistogram' not in self.vf.userpref.keys():
            doc = """Choose to show or not show Clustering Histogram for a VirtualScreening result:\n'Yes' or 'No'\n. Default value is 'No'"""
            self.vf.userpref.add('VS_ShowClusteringHistogram', 'Yes', 
                    ['Yes', 'No', 'None'],
                    callbackFunc = [self.set_VS_ShowClusteringHistogramMode], doc=doc ,
                    category='AutoDockTools')
            self.vf.userpref['VS_ShowClusteringHistogram']['value'] = 'No'
                    
                    

    def set_VS_ShowClusteringHistogramMode(self, name, oldval, newval): 
        self.vf.userpref['VS_ShowClusteringHistogram']['value'] = newval


    def set_VSCloseContacts_displayMode(self, name, oldval, newval): 
        self.vf.userpref['VSCloseContacts_displayMode']['value'] = newval


    def set_VSHydrogenBond_displayMode(self, name, oldval, newval): 
        self.vf.userpref['VSHydrogenBond_displayMode']['value'] = newval


    def set_VSPiInteractions_displayMode(self, name, oldval, newval): 
        self.vf.userpref['VSPiInteractions_displayMode']['value'] = newval



    # CLEAN UP THE receptor here
    def onRemoveObjectFromViewer(self, obj):
        #remove the hbonds to this atoms in this molecule from the receptor..
        other = None
        has_hbonds = 0
        displayMode = self.vf.userpref['VSCloseContacts_displayMode']['value']
        if obj.__dict__.has_key('LE_close_ats'): #isV2
            if len(self.vf.Mols): other = self.vf.Mols[-1]
            #print "V2: other.name=", other.name
            kkeys = ['LE_close_ats', 'LC_close_ats']
            for k in kkeys:
                other_ats = obj.__dict__[k]
                if len(other_ats):
                    other = other_ats[0].top
                    if displayMode=='asSelection':
                        self.vf.deselect(other_ats)
                    elif displayMode=='asCPK':
                        self.vf.undisplayCPK(other_ats)    
                    elif displayMode=='asSticks':
                        self.vf.undisplaySticksAndBalls(other_ats)
            hats = obj.allAtoms.get(lambda x: hasattr(x, 'hbonds'))
            if len(hats):
                has_hbonds = 1
                for at in hats:
                    for b in at.hbonds:
                        b.donAt.hbonds.remove(b)
                        b.accAt.hbonds.remove(b)
                        if b.hAt is not None:
                            b.hAt.hbonds.remove(b)
                        del(b)
            other_hats = other.allAtoms.get(lambda x: hasattr(x, 'hbonds'))
            if len(other_hats):
                for at in other_hats:
                    for b in at.hbonds:
                        try:
                            b.donAt.hbonds.remove(b)
                        except: pass
                        try:
                            b.accAt.hbonds.remove(b)
                        except: pass
                        if b.hAt is not None:
                            try:
                                b.hAt.hbonds.remove(b)
                            except: pass
                        del(b)
                    try:
                        delattr(at, 'hbonds')
                        delattr(at, '_hbonds')
                    except:
                        pass
        else: #isV1
            #?who is other
            if len(self.vf.Mols):
                other = self.vf.Mols[-1]
                #print "V1: other.name=", other.name
                hats = other.allAtoms.get(lambda x: hasattr(x, 'hbonds'))
                if len(hats):
                    has_hbonds=1
                    for at in hats: delattr(at, 'hbonds')
                    #msg = other.name + " had hbonds"
                    #self.vf.warningMsg(msg)
                if hasattr(other, 'close_ats') and len(other.close_ats):
                    if displayMode=='asSelection':
                        self.vf.deselect(other.close_ats)
                    elif displayMode=='asCPK':
                        self.vf.undisplayCPK(other.close_ats)    
                    elif displayMode=='asSticks':
                        self.vf.undisplaySticksAndBalls(other.close_ats)
                    delattr(other, 'close_ats')
        #@@ DO APPROPRIATE THING
        if has_hbonds:
            if self.vf.userpref['VSHydrogenBond_displayMode']['value']=="asSpheres":
                try:
                    self.vf.hbondsAsSpheres.dismiss_cb()
                except:
                    self.warningMsg("Error calling hbondsAsSpheres.dismiss_cb")
            if self.vf.userpref['VSHydrogenBond_displayMode']['value']=="asCylinders":
                try:
                    self.vf.hbondsAsCylinders.dismiss_cb()
                except:
                    self.warningMsg("Error calling hbondsAsCylinders.dismiss_cb")
            if self.vf.userpref['VSHydrogenBond_displayMode']['value']=="asLines":
                try:
                    self.vf.showHBonds.dismiss_cb()
                except:
                    self.warningMsg("Error calling showHBonds.dismiss_cb")
        #if self.vf.userpref['VS_ShowClusteringHistogram']['value']=="Yes":
        #    try:
        #         print "do something here"
        #    except:
        #        self.warningMsg("Error calling showHBonds.dismiss_cb")


        #is it necessary to do anything to these? reset??
        #receptor.hbond_ct = 0 receptor.hb = AtomSet() #???  receptor.cc = AtomSet() receptor.pi = AtomSet()
        #@@ DO APPROPRIATE THING
        if self.vf.userpref['VSPiInteractions_displayMode']['value']=="asCylindersAndCones":
            pi_pi_geom = self.vf.GUI.VIEWER.findGeomsByName("pi_pi")
            if len(pi_pi_geom):
                pi_pi_geom[0].Set(visible=0)
            pi_cation_geom = self.vf.GUI.VIEWER.findGeomsByName("pi_cation")
            if len(pi_cation_geom):
                pi_cation_geom[0].Set(visible=0)
            cation_pi_geom = self.vf.GUI.VIEWER.findGeomsByName("cation_pi")
            if len(cation_pi_geom):
                cation_pi_geom[0].Set(visible=0)
        if other is None:
            for mol in self.vf.Mols:   
                if mol!=obj:
                    other = mol
        if not other:
            return
        if not hasattr(other, 'close_ats'):
            return
        close_ats = other.close_ats
        if len(self.vf.Mols):
            name = self.vf.Mols[0].name
        if len(self.vf.Mols)>1:
            for z in range(len(self.vf.Mols)):
                name += ',' +self.vf.Mols[z].name
            #print "resettting name to", name
            self.vf.Mols.setStringRepr(name)
            #for LC-case ligand display cleanup
        self.vf.GUI.VIEWER.Redraw()


    def guiCallback(self):
        """called each time the 'open virtual screening result' menubutton is pressed"""
        resultFile = self.vf.askFileOpen(types=[('select virtual screening filename:', '*.pdbqt'),
                            ('all files','*')],
                            title = 'VirtualScreening Result File:')
        if resultFile:
            kw = {}
            #ask = 1
            apply( self.doitWrapper, (resultFile,),  kw)


    def processArrowEvent(self, event):
        #detect molecules with conformations as well as vina_results
        # update VSCloseContacts + VSHydrogenBonds
        displayMode = self.vf.userpref['VSCloseContacts_displayMode']['value']
        hbond_displayMode = self.vf.userpref['VSHydrogenBond_displayMode']['value'] 
        from Pmv.moleculeViewer import EditAtomsEvent
        mols = self.vf.Mols.get(lambda x: hasattr(x, 'isVSResult') and len(x.allAtoms[0]._coords[0])>1)
        if not len(mols):
            return
        if hasattr(event, 'keysym') and 'Right' in event.keysym:
            #print "Right"
            for m in mols:
                #print m.name
                #geom = m.geomContainer.geoms['ProteinLabels']
                ind = m.allAtoms[0].conformation
                nconf = len(m.allAtoms[0]._coords)
                if ind+1 < nconf:
                    ind = ind+1
                    m.allAtoms.setConformation(ind)
                    event2 = EditAtomsEvent('coords', m.allAtoms)
                    self.vf.dispatchEvent(event2)
                    #@@ 4-12-11: TODO @@
                    if hasattr(m, 'bindingSite'):
                        delattr(m, 'bindingSite')
                        c = self.vf.ADanalyze_showBindingSite
                        c.build()
                #if hasattr(m, 'hbonds') and len(m.hbonds):
                HATS = AtomSet()
                HATS += m.allAtoms.get(lambda x: hasattr(x, 'hbonds') and len(x.hbonds))
                HATS += self.receptor.allAtoms.get(lambda x: hasattr(x, 'hbonds') and len(x.hbonds))
                if len(HATS):
                    if hbond_displayMode=='asSpheres':
                        self.vf.hbondsAsSpheres(HATS)
                        #self.vf.hbondsAsSpheres(m.allAtoms,self.receptor.allAtoms)
                    elif hbond_displayMode=='asLines':
                        self.vf.showHBonds(HATS)
                    elif hbond_displayMode=='asCylinders':
                        self.vf.hbondsAsCylinders(HATS)
                    elif hbond_displayMode=='asCylinders':
                        self.vf.extrudeHBonds(HATS)
                if ind==0 and hasattr(m, 'LE_close_ats') and len(m.LE_close_ats): 
                    if displayMode=='asSelection': 
                        self.vf.select(m.LE_close_ats, log=0, only=1)    
                    elif displayMode=='asCPK':
                        self.vf.displayCPK(m.LE_close_ats, log=0, only=1, scaleFactor=0.6)    
                    elif displayMode=='asSticks':
                        self.vf.displaySticksAndBalls(m.LE_close_ats, log=0, sticksBallsLicorice="Licorice", only=1)    
                    if m.clustNB.winfo_ismapped():
                        m.clustNB.draw.itemconfig((m.le_index,), fill='red')
                        m.clustNB.draw.itemconfig((m.lc_index,), fill='blue')
                if ind==1 and hasattr(m, 'LC_close_ats') and len(m.LC_close_ats):
                    if displayMode=='asSelection':
                        self.vf.select(m.LC_close_ats, log=0, only=1)    
                    elif displayMode=='asCPK':
                        self.vf.displayCPK(m.LC_close_ats, log=0, only=1, scaleFactor=0.6)    
                    elif displayMode=='asSticks':
                        self.vf.displaySticksAndBalls(m.LC_close_ats, log=0, sticksBallsLicorice="Licorice", only=1)    
                    if hasattr(m, 'clustNB'):
                        if m.clustNB.winfo_ismapped():
                            m.clustNB.draw.itemconfig((m.lc_index,), fill='red')
                            m.clustNB.draw.itemconfig((m.le_index,), fill='blue')
                else: #hbonds
                    #print "OFF THE END!!"
                    print 'end: last virtual screening result for ', m.name
                    #event = EditAtomsEvent('coords', m.allAtoms)
                    #self.vf.dispatchEvent(event)
                    #print " end "#, ind, "  vina_energy =", m.vina_results[-1][0]
        if hasattr(event, 'keysym') and 'Left' in event.keysym:
            #print "Left"
            for m in mols:
                #print m.name,
                ind = m.allAtoms[0].conformation
                if ind > 0:
                    ind = ind - 1
                    m.allAtoms.setConformation(ind)
                    event = EditAtomsEvent('coords', m.allAtoms)
                    self.vf.dispatchEvent(event)
                    #g = m.geomContainer.geoms['ProteinLabels']
                    #g.Set(labels = (str(m.vina_results[ind][0]),))
                    #print "LEFT: vina_results[",ind,"][0]=", m.vina_results[ind][0]
                    #if ind==0:
                    #    print "over-all best vina_energy=", m.vina_results[ind][0]
                    #else:
                    #    print ind+1," vina_energy=", m.vina_results[ind][0]
                    if hasattr(m, 'bindingSite'):
                        delattr(m, 'bindingSite')
                        c = self.vf.ADanalyze_showBindingSite
                        c.build()
                else:
                    #print "LESS THAN ZERO!!"
                    print 'back to first: best VSresult for ', m.name
                    #print " at first"
                    #print " over-all best vina_energy=", m.vina_results[0][0]
                event2 = EditAtomsEvent('coords', m.allAtoms)
                self.vf.dispatchEvent(event2)
                #if hasattr(m, 'LE_close_ats') and len(m.LE_close_ats):
                HATS = AtomSet()
                HATS += m.allAtoms.get(lambda x: hasattr(x, 'hbonds') and len(x.hbonds))
                HATS += self.receptor.allAtoms.get(lambda x: hasattr(x, 'hbonds') and len(x.hbonds))
                if len(HATS):
                    if hbond_displayMode=='asSpheres':
                        self.vf.hbondsAsSpheres(HATS)
                        #self.vf.hbondsAsSpheres(m.allAtoms,self.receptor.allAtoms)
                    elif hbond_displayMode=='asLines':
                        self.vf.showHBonds(HATS)
                    elif hbond_displayMode=='asCylinders':
                        self.vf.hbondsAsCylinders(HATS)
                    elif hbond_displayMode=='asCylinders':
                        self.vf.extrudeHBonds(HATS)
                if ind==0 and hasattr(m, 'LE_close_ats') and len(m.LE_close_ats): 
                    if displayMode=='asSelection':
                        self.vf.select(m.LE_close_ats, log=0, only=1)    
                    elif displayMode=='asCPK':
                        self.vf.displayCPK(m.LE_close_ats, log=0, only=1, scaleFactor=0.6)    
                    elif displayMode=='asSticks':
                        self.vf.displaySticksAndBalls(m.LE_close_ats, log=0, sticksBallsLicorice="Licorice", only=1)    
                    if hasattr(m, 'clustNB') and m.clustNB.winfo_ismapped():
                        m.clustNB.draw.itemconfig((m.le_index,), fill='red')
                        m.clustNB.draw.itemconfig((m.lc_index,), fill='blue')
                #if hasattr(m, 'LC_close_ats') and len(m.LC_close_ats):
                if ind==1 and hasattr(m, 'LC_close_ats') and len(m.LC_close_ats):
                    if displayMode=='asSelection':
                        self.vf.select(m.LC_close_ats, log=0, only=1)    
                    elif displayMode=='asCPK':
                        self.vf.displayCPK(m.LC_close_ats, log=0, only=1, scaleFactor=0.6)    
                    elif displayMode=='asSticks':
                        self.vf.displaySticksAndBalls(m.LC_close_ats, log=0, sticksBallsLicorice="Licorice", only=1)    
                    if hasattr(m, 'clustNB'):
                        if m.lc_index!=m.le_index:
                            if m.clustNB.winfo_ismapped():
                                m.clustNB.draw.itemconfig((m.lc_index,), fill='blue')
                                m.clustNB.draw.itemconfig((m.le_index,), fill='red')
        #hbonds
        #update closeContacts and ???
        


    def __call__(self, resultFile=None, **kw):
        """None<-ADanalyze_readVSResult(resultFile=None)
        \nresultFile: augmented pdbqt file containing docked coordinates plus
                      various information about clustering and interactions 
                      between docked ligand and receptor...
        """
        if resultFile is None:
            return 'ERROR'
        apply(self.doitWrapper,(resultFile,), kw)


    def doit(self, resultFile, **kw):
        #retrieve the receptor name first  and check that it's in the viewer...
        #from MolKit import Read
        #prelim = Read(resultFile)[0]
        newparser = PdbqtParser(resultFile, modelsAs='conformations')
        prelim = newparser.parse()
        if len(prelim): 
            ligand=prelim[0]
            ligand.isVSResult = 1
            parser = ligand.parser
        else: 
            msg = "ERROR: unable to open %s" %resultFile 
            return msg
        if parser.isV2:
            rec_name = parser.__dict__['AD_rec']
            rms = str(parser.__dict__['rmstol'])
            num_runs = str(parser.__dict__['AD_runs'])
        elif hasattr(parser, 'macro_close_ats') and len(parser.macro_close_ats):
            #['xJ1_xtal:B:GLY16:N,HN~ZINC02025973_vs:d:<0>:O3\n','xJ1_xtal:B:LEU63:N,HN~ZINC02025973_vs:d:<0>:O5\n']    
            first = parser.macro_close_ats[0]
            #'xJ1_xtal:B:GLY16:N,HN~ZINC02025973_vs:d:<0>:O3\n'
            rec_name = first.split(':')[0]
        else:
            msg = "Improperly formated pdbqt+ file! receptor not specified in " + resultFile
            self.warningMsg(msg)
            return 'ERROR'
        if rec_name in self.vf.Mols.name:
            index = self.vf.Mols.name.index(rec_name)
            #print rec_name, " already in viewer"
            #receptor = self.vf.Mols[index]
        else:
            #msg = "Please read in receptor " + rec_name + '.pdbqt' + " first"
            msg = rec_name + " not in viewer! Please read it in ..."
            self.warningMsg(msg)
            receptorFile = self.vf.askFileOpen(types=[('select virtual screening receptor filename:', '*.pdbqt'),
                        ('all files','*')],
                        title = 'VirtualScreening Receptor File:')
            if receptorFile:
                self.vf.readMolecule(receptorFile)
                index = -1
            else:
                return 'ERROR'
        receptor = self.receptor = self.vf.Mols[index]
        # The receptor is in the viewer! rec_name NOW read in the ligand and process...
        # reset receptor for new ligand
        receptor.hbond_ct = 0
        receptor.hb = AtomSet() #???
        #@@ ??? ok to reset self.vf.docked ???
        if parser.isV2:
            if self.first:
                #check for preexisting callback from MoleculeReader
                #eM = self.vf.GUI.VIEWER.currentCamera.eventManager
                #for key in ["<Right>", "<Left>"]:
                #    if eM.eventHandlers.has_key(key):
                #        eM.RemoveCallback(key,self.vf.readMolecule.processArrowEvent)
                #        print "removed readMolecule.processArrowEvent for ", key
                self.vf.GUI.addCameraCallback("<Right>", self.processArrowEvent)
                self.vf.GUI.addCameraCallback("<Left>", self.processArrowEvent)
                self.first = 0
            if len(ligand.allAtoms[0]._coords)>1:
                msg = "Added %d docked conformations to %s\nUse keyboard arrow keys to view" %( len(ligand.allAtoms[0]._coords), ligand.name)
                self.vf.warningMsg(msg)
        d = self.vf.docked = Docking()
        top = None
        if  self.vf.userpref['VS_ShowClusteringHistogram']['value']=='Yes':
            #provide toplevel for showing the histogram
            top = Tkinter.Toplevel(master=self.vf.GUI.ROOT)
            try:
                tstr = ligand.name + ': rms=' + rms + ' ' + num_runs + ' runs'
            except:
                pass  #FIX THIS @@
            if parser.isV2:
                hist_list = ligand.parser.__dict__['AD_histogram']
                #['-7.720:4:**', '-7.430:3', '-7.140:2', '-7.090:11:*', 
                # '-7.070:5', '-7.000:1', '-6.990:8', '-6.950:2', 
                #'-6.810:3', '-6.770:9', '-6.770:7', '-6.720:2', 
                #'-6.600:3', '-6.570:4', '-6.570:3', '-6.560:1', 
                #'-6.460:1', '-6.420:4', '-6.420:3', '-6.370:4', 
                #'-6.370:6', '-6.330:2', '-6.330:5', '-6.310:1', 
                #'-6.290:3', '-6.250:5', '-6.240:1', '-6.180:1', 
                #'-6.170:1', '-6.150:1', '-6.130:1', '-6.080:1', 
                #'-6.050:1', '-5.990:1', '-5.980:1', '-5.980:1', 
                #'-5.980:1', '-5.900:2', '-5.900:1', '-5.820:1', 
                #'-5.730:1', '-5.600:1', '-5.540:1']
                #@@ line 382 Docking.py
                dataList = []
                reverseList = []
                ctr = 2
                cl_ctr = 0
                for item in hist_list:
                    item_list = item.split(':')
                    dataList.append(map(float,item_list[:2]))
                    if len(item_list)>2:
                        if item_list[-1]=='**':
                            LE = item_list
                        if item_list[-1]=='*':
                            LC = item_list
                    num_confs = int(item_list[1])
                    if num_confs>1: 
                        reverseList.append(range(ctr, ctr+num_confs+1))
                    ctr += 1
                    cl_ctr += 1
                top.title(tstr)
                xlabel = 'BINDING ENERGY'
                ligand.clustNB = interactiveHistogramGraph.InteractiveHistogramGraph(ligand.name,
                        master=top, nodeList = dataList, reverseIndex=reverseList,
                        label_text=ligand.name + ':' + rms + ' rms', xlabel_text=xlabel, 
                        ylabel_text='#\nC\nO\nN\nF\nO\nR\nM\nA\nT\nI\nO\nN\nS')
                le = 10000
                max_height = -1
                max_index = None
                for i,j in ligand.clustNB.geoms.items():
                    if j.point[0]<le:
                        le = j.point[0]
                        le_index = i
                    if j.height>max_height:
                        max_height = j.height
                        lc_index = i
                ligand.clustNB.draw.itemconfig((le_index,), fill='red')
                #always color le_index red to start...
                #remember lc_index to support showing LC
                ligand.le_index = le_index
                ligand.lc_index = lc_index

        if parser.isV1:
            resultlist = d.readVSResult(resultFile, receptor, top=top)
            if resultlist is None:
                return "ERROR in reading VSResult (ligand) ", resultFile
            ligand = resultlist[0]
        #mols = self.vf.readMolecule(pdbqtFile, modelsAs=modelsAs, setupUpdates=setupUpdates)
        self.vf.addMolecule(ligand, ask=True)
        ligand.buildBondsByDistance()
        parser = ligand.parser
        #PROCESS INFORMATION FROM PARSER-hb,cc and pi-pi
        #close_contacts and hydrogen bonds
        isV2 = parser.isV2
        css = CompoundStringSelector()
        LE_close_ats = AtomSet()
        LC_close_ats = AtomSet()
        if not isV2:
            LE_close_ats = ligand.close_ats+receptor.close_ats
        else:
            LE_close_ats_list = parser.__dict__['AD_LE_vdw']
            for cal_str in LE_close_ats_list:
                more_ats, more_ats_str = css.select(self.vf.Mols, strip(cal_str))
                if len(more_ats):
                    LE_close_ats+=more_ats
            ligand.LE_close_ats = LE_close_ats
            LC_close_ats_list = parser.__dict__['AD_LC_vdw']
            for cal_str in LC_close_ats_list:
                more_ats, more_ats_str = css.select(self.vf.Mols, strip(cal_str))
                if len(more_ats):
                    LC_close_ats+=more_ats
            ligand.LC_close_ats = LC_close_ats
            ligand.hbond_ct = 0
            #build hbonds here
            le_hb_list = []
            for kk in ['AD_LE_hbd', 'AD_LE_hba']: 
                le_hb_list.extend(parser.__dict__[kk])
            lc_hb_list = []
            for kk in ['AD_LC_hbd', 'AD_LC_hba']: 
                lc_hb_list.extend(parser.__dict__[kk])
            ssStr = "~~"
            for hlist, hbs in zip([le_hb_list, lc_hb_list],[[],[]]):
                for hb_str in hlist:
                    #VALIDITY CHECK
                    ok = False
                    #['xJ1_xtal:B:GLY16:N,HN~ZINC02025973_vs:d:<0>:O3\n','xJ1_xtal:B:LEU63:N,HN~ZINC02025973_vs:d:<0>:O5\n']    
                    # v2
                    #d:<0>:N1~~B:ASN83:N,d:<0>:O3~~B:LYS20:NZ
                    str_1, str_2 = hb_str.split(ssStr)
                    if not len(str_1) or not len(str_2):
                        print "skipping improper hb_str ", hb_str
                        continue
                    ligStr = ligand.chains[0].id #@@FIX THIS
                    ok = True
                    if ligStr in str_1 and rec_name in str_2:
                        lig_str = str_1
                        rec_str = str_2
                    elif ligStr in str_2 and rec_name in str_1:
                        lig_str = str_2
                        rec_str = str_1
                    else:
                        print "ligand name: ", ligand.name, " not found in hb_str ", hb_str
                        ok = False                    
                        print "not ok!: lig_str=", lig_str, ' rec_str=', rec_str
                        raise 'ligand name!'
                    if ok:
                        #print "ok!: lig_str=", lig_str, ' rec_str=', rec_str
                        #rec_str xJ1_xtal:B:GLY16:N,HN
                        # try to retrieve receptor atoms
                        rec_ats = css.select(self.vf.Mols, strip(rec_str))
                        #rec_ats = self.vf.expandNodes(strip(rec_str))
                        if not len(rec_ats):
                            print "skipping hb_str because ",  rec_str, " did not match atoms in the receptor"
                        #(<AtomSet instance> holding 2 Atom, "xJ1_xtal:B:ARG57:NE,HE", '')
                        rec_ats = rec_ats[0]
                        #add ligand name to lig_str
                        lig_ats = css.select(self.vf.Mols, strip(lig_str))[0]
                        #(<AtomSet instance> holding 1 Atom, "ZINC02025973_vs_le:d:<0>:O3", '')
                        #lig_ats = self.vf.expandNodes(strip(lig_str))
                        if not len(lig_ats):
                            print "skipping hb_str because ",  lig_str, " did not match atoms in the ligand"
                        #(<AtomSet instance> holding 1 Atom, "ZINC02026663_vs:d:<0>:O3", '')
                        #lig_ats = lig_ats[0]
                        donor_ats = rec_ats
                        accAt = lig_ats[0]
                        if "," in lig_str: 
                            donor_ats = lig_ats
                            accAt = rec_ats[0]
                        if not hasattr(accAt, 'hbonds'):
                            #accAt.hbonds = HydrogenBondSet()
                            accAt._hbonds = [HydrogenBondSet()]
                            accAt.hbonds = accAt._hbonds[0]
                        #?? or ??
                        donAt = donor_ats[0]
                        if not hasattr(donAt, 'hbonds'):
                            donAt._hbonds = [HydrogenBondSet()]
                            donAt.hbonds = donAt._hbonds[0]
                        hAt = None
                        if len(donor_ats)==2:
                            hAt = donor_ats[1]
                            if not hasattr(hAt, 'hbonds'):
                                #hAt.hbonds = HydrogenBondSet()
                                hAt._hbonds = [HydrogenBondSet()]
                                hAt.hbonds = hAt._hbonds[0]
                        try:
                            hb = HydrogenBond(donAt, accAt, hAt)
                            #add them to these atoms
                            donAt.hbonds.append(hb)
                            accAt.hbonds.append(hb)
                        except:
                            msg = "%s exception on creating hydrogen bond between %s and %s " %(hb_str,donAt.full_name(), accAt.full_name())
                            self.vf.warningMsg(msg)
                            continue
                        if hAt is not None:
                            if not hasattr(hAt, 'hbonds'):
                                #hAt.hbonds = HydrogenBondSet()
                                hAt._hbonds = [HydrogenBondSet()]
                                hAt.hbonds = hAt._hbonds[0]
                            hAt.hbonds.append(hb)
                        receptor.hbond_ct += 1
                        ligand.hbond_ct += 1
                        #whew!
        if len(LE_close_ats):
            displayMode = self.vf.userpref['VSCloseContacts_displayMode']['value']
            if displayMode=='asSelection':
                self.vf.select(LE_close_ats, log=0, only=1)    
            elif displayMode=='asCPK':
                self.vf.displayCPK(LE_close_ats, log=0, only=1, scaleFactor=0.6)    
            elif displayMode=='asSticks':
                self.vf.displaySticksAndBalls(LE_close_ats, log=0, sticksBallsLicorice="Licorice", only=1)    
        #hbonds
        # 4-15-2011: @@ TODO @@
        #if isV2:
            #receptor.hbond_ct = 0
            #for k in ['AD_LE_hba', 'AD_LE_hbd', 'AD_LC_hba', 'AD_LC_hbd']:
            #    try:
            #        num_hbonds = len(parser.__dict__[k])
            #        ligand.hbond_ct += num_hbonds
            #        receptor.hbond_ct += num_hbonds
            #        #build the hbonds
            #        for sss in parser.__dict[k]:
            #            #'*:d:<0>:O3~~xJ1_xtal:B:LEU63:N' 
            #            if sss.find("``")>-1:
            #                ligStr,recStr = sss.split('~~')
            #    except:
            #        pass
        if  receptor.hbond_ct + ligand.hbond_ct>0: 
            if self.vf.userpref['VSHydrogenBond_displayMode']['value'] == 'asSpheres':
                self.vf.hbondsAsSpheres(receptor.allAtoms+ligand.allAtoms, log=0)
            if self.vf.userpref['VSHydrogenBond_displayMode']['value'] == 'asCylinders':
                self.vf.hbondsAsCylinders(receptor.allAtoms+ligand.allAtoms, log=0)
            if self.vf.userpref['VSHydrogenBond_displayMode']['value'] == 'asLines':
                self.vf.showHBonds(receptor.allAtoms+ligand.allAtoms, log=0)
        #@@ pi @@
        if self.vf.userpref['VSPiInteractions_displayMode']['value'] != 'None':
            #pi-cation
            pi_cation_list = []
            if hasattr(ligand.parser, 'pi_cation') and len(ligand.parser.pi_cation):
                pi_centers = []
                cation_centers = []
                for pi_c_str in ligand.parser.pi_cation:
                #['xJ1_xtal:B:LYS55:NZ~~ZINC02026663_vs:d:<0>:C6,C7,C3,C4,C8,C5']    
                    pi_str, cation_str = pi_c_str.split("~~")
                    pi_ats = self.vf.expandNodes(pi_str)
                    if len(pi_ats)==1:
                        #pi_ats_coords= pi_ats[0].coords 
                        cation_ats = self.vf.expandNodes(cation_str)
                        if len(cation_ats)>1:
                            pi_cation_list.append([pi_ats[0], cation_ats])
                    elif len(pi_ats)>1:
                        if len(cation_ats):
                            pi_cation_list.append([cation_ats[0], pi_ats])

            #pi-pi
            pi_pi_list = []
            if hasattr(ligand.parser, 'pi_pi') and len(ligand.parser.pi_pi):
                for pi_pi_str in ligand.parser.pi_pi:
                    pi_str1, pi_str2 = pi_pi_str.split("~~")
                    pi_ats1 = self.vf.expandNodes(pi_str1)
                    pi_ats2 = self.vf.expandNodes(pi_str2)
                    pi_pi_list.append([pi_ats1, pi_ats2])
                    print 'do something here'
            if not len(pi_cation_list) and not len(pi_pi_list):
                return
            if self.vf.userpref['VSPiInteractions_displayMode']['value'] == 'asCylindersAndCones':
                if len(pi_pi_list): #pi-pi
                    vertices = []
                    faces = []
                    radii = []
                    centers = []
                    ctr = 0
                    for first,second in pi_pi_list:
                        center1 = (Numeric.add.reduce(first.coords)/len(first)).tolist() 
                        radii.append(0.7)
                        vertices.append(center1)
                        center2 = (Numeric.add.reduce(second.coords)/len(second)).tolist() 
                        radii.append(0.7)
                        vertices.append(center2)
                        faces.append((ctr,ctr+1))
                        ctr += 2
                    pi_pi_geom = self.vf.GUI.VIEWER.findGeomsByName("pi_pi")
                    if len(pi_pi_geom): 
                        pi_pi_geom[0].Set(faces=faces, radii=radii, vertices=vertices, visible=1)
                    else:
                        return "no pi_pi geom in viewer"
                if len(pi_cation_list): #pi-cation
                    vertices = []
                    faces = []
                    radii = []
                    ctr = 0
                    for first,second in pi_cation_list:
                        if hasattr(first, '__len__'):
                            center1 = (Numeric.add.reduce(first.coords)/len(first)).tolist()
                            radii.append(0.7)
                        else:
                            center1 = (first.coords)
                            radii.append(0.1)
                        if hasattr(second, '__len__'): #<AtomSet instance> holding 6 Atom, "ZINC02026663_vs:d:<0>:C6,C7,C3..."
                            center2 = (Numeric.add.reduce(second.coords)/len(second)).tolist()
                            radii.append(0.7)
                        else:
                            center2 = (second.coords)
                            radii.append(0.1)
                        faces.append((ctr,ctr+1))
                        vertices.append(center1)
                        vertices.append(center2)
                        ctr += 2
                    pi_cation_geom = self.vf.GUI.VIEWER.findGeomsByName("pi_cation")
                    if len(pi_cation_geom): 
                        pi_cation_geom[0].Set(radii=radii, faces=faces, 
                                                vertices=vertices, visible=1)
                    else:
                        return "no pi_geom in viewer"
            if self.vf.userpref['VSPiInteractions_displayMode']['value']=="asSelection":
                for pi_c_str in ligand.parser.pi_cation:
                #['xJ1_xtal:B:LYS55:NZ~~ZINC02026663_vs:d:<0>:C6,C7,C3,C4,C8,C5']    
                    pi_str, cation_str = pi_c_str.split("~~")
                    pi_ats = self.vf.expandNodes(pi_str)
                    cation_at = self.vf.expandNodes(cation_str)
                    self.vf.select(cation_at)
                    self.vf.select(pi_ats)
                for pi_pi_str in ligand.parser.pi_pi:
                    pi_str1, pi_str2 = pi_pi_str.split("~~")
                    pi_ats1 = self.vf.expandNodes(pi_str1)
                    pi_ats2 = self.vf.expandNodes(pi_str2)
                    if len(pi_ats1): self.vf.select(pi_ats1)
                    if len(pi_ats2): self.vf.select(pi_ats2)
                    self.vf.select(pi_ats2)
            self.vf.GUI.VIEWER.Redraw()


    def describeClust(self, mol, event=None):
        print "in VS describeClust"
        ##raise 'describeClust'
        return

        

ADReadVSResultGUI=CommandGUI()
ADReadVSResultGUI.addMenuCommand('AutoToolsBar', menuText['AnalyzeMB'], 
        menuText['readVSResult'], cascadeName = menuText['MoleculesMB'])


        
class ADWriteVSResult(MVCommand):
    """This class is used to write a pdbqt+ file from the current conformation of the current docking,  
    \nPackage : AutoDockTools
    \nModule  : autoanalyzeCommands
    \nClass   : ADWriteVSResult
    \nCommand : ADanalyze_writeVSResult
    \nSynopsis:\n 
        None<-ADanalyze_writeVSResult(resultFile=None)
    \nArguments:\n   
        resultFile: pdbqt+  filename for current AutoDock docking to include information on
                        clustering, hydrogen-bonds to receptor, close-contacts ...
    """
    def onAddCmdToViewer(self):
        checkHasInitializedDockings(self.vf)
        if self.vf.hasGui:
            for item in ['colorByMolecules','colorByAtomType', 'color']:
                if not hasattr(self.vf, item):
                    self.vf.loadCommand('colorCommands', item, 'Pmv')
            for item in ['displayLines','showMolecules']:
                if not hasattr(self.vf, item):
                    self.vf.loadCommand('displayCommands', item, 'Pmv')
                    

    def guiCallback(self):
        """called each time the 'save virtual screening result' menubutton is pressed"""
        docked = self.vf.docked
        if docked is None:
            msg= "There is no current docking in the viewer. Please load one first using Analyze->Dockings->Open..."
            self.warningMsg(msg)
            return 'ERROR'
        if isinstance(docked, Molecule):
            msg = "'Write as VSResult' is not available because current docking is an AutoDock Vina result." 
            self.warningMsg(msg)
            return 'ERROR'
        conf = None
        #if hasattr(docked,'ch'):
        #    conf = docked.ch.current_conf
        #if conf is None:
        #    msg= "There is no current docked conformation.\nPlease use Analyze->Clusterings->Show..\nand click on a bar in the histogram to select a conformation"
        #    self.warningMsg(msg)
        #    return 'ERROR'    
        dlo = docked.dlo_list[0]
        macroStem = dlo.macroStem
        if macroStem not in self.vf.Mols.name:
            msg= "Please use 'Analyze->Macromolecule->Open...' to add %s to the viewer first"%(macroStem+".pdbqt")
            self.warningMsg(msg)
            return 'ERROR'
        #build the interactions
        #@@ 5/6/2011 pi-pi not supported yet
        #t = 'Detect pi-pi interactions?'
        #d = SimpleDialog(self.vf.GUI.ROOT, text = t,
        #buttons = ["No", "Yes", "Cancel"], default = 1, 
        #        title = "Pi-pi interactions?")
        #detect_pi = d.go()
        #if detect_pi==2:
        #    return 'ERROR'
        detect_pi = 0
        resultFile = self.vf.askFileSave(types=[('docked result file','*_VS.pdbqt')], title = 'Docked Ligand File:')
        if resultFile:
            kw = {'detect_pi': detect_pi}
            apply( self.doitWrapper, (resultFile,),  kw)


    def __call__(self, resultFile=None, **kw):
        """None<-ADanalyze_writeVSResult(resultFile=None)
        \nresultFile: augmented pdbqt file containing docked coordinates plus
                      various information about clustering and interactions 
                      between docked ligand and receptor...
        """
        apply(self.doitWrapper,(resultFile,), kw)


    def doit(self, resultFile, **kw):
        if not self.vf.docked:
            msg= "There is no current docking in the viewer. Please load a docked result first"
            self.warningMsg(msg)
            return 'ERROR'
        docked = self.vf.docked
        #conf = docked.ch.current_conf
        #if conf is None:
        #    msg= "There is no current docked conformation.\nPlease use Analyze->Clustering->Show..\nand click on a bar in the histogram to select a conformation"
        #    self.warningMsg(msg)
        #    return 'ERROR'
        #ind = str(conf.run)
        if resultFile is None:
            resultFile = docked.ligMol.name + "_VS" + ".pdbqt"    
        new_ligand_name = os.path.splitext(os.path.basename(resultFile))[0]
        #retrieve the receptor name first  and check that it's in the viewer...
        macroMol = None
        if hasattr(docked,'macroMol'):
            macroMol = docked.macroMol
        #verify that macromolecule is in the viewer
        if macroMol is None:
            msg= " No macromolecule for current docking.\nPlease use Analyze->Macromolecule->Open.. or Analyze->Macromolecule->Choose... first"
            self.vf.warningMsg(msg)
            return 'ERROR'
        receptor_filename = macroMol.parser.filename     
        dlg_list = []
        for dlo in docked.dlo_list:
            dlg_list.append(dlo.filename)
        frp = FoxResultProcessor(receptor_filename) 
        frp.process(dlg_list=dlg_list, rms_tolerance=2.0, outputfilename=resultFile)


###        #build the interactions
###        detect_pi = kw.get('detect_pi', False)
###        self.intF = InteractionDetector(receptor_filename, detect_pi=detect_pi) #someday, detect_pi=True
###        #setup to write output 
###        max_cl_to_write = kw.get("max_cl_to_write", -1)
###        nconfs = len(docked.ch.conformations)
###        #docking should already have clusterer + clustering
###        rmsd = kw.get('rmsd', None)
###        ind = str(conf.run)
###        if resultFile is None:
###            resultFile = docked.ligMol.name + "_vs_" + ind + ".pdbqt"    
###        new_ligand_name = os.path.splitext(os.path.basename(resultFile))[0]
###        #setup to write output 
###        max_cl_to_write = kw.get("max_cl_to_write", -1)
###        nconfs = len(docked.ch.conformations)
###        #docking should already have clusterer + clustering
###        rmsd = kw.get('rmsd', None)
###        if rmsd is None and len(docked.clusterer.clustering_dict.keys()):
###            rmsd = docked.clusterer.clustering_dict.keys()[0]
###        if rmsd is None:
###            msg = "Current docking has no clustering! Please cluster results first!"
###            self.warningMsg(msg)
###            return "ERROR"
###        cl_dict = docked.clusterer.clustering_dict[rmsd]
###        nclusters = len(cl_dict)
###        cl_lengths = map(len, cl_dict)
###        # adjust the number of clusters to write down if necessary
###        if max_cl_to_write<0 or max_cl_to_write > nclusters:
###            max_cl_to_write = nclusters
###            #print 'set max_cl_to_write to nclusters ', nclusters
###        #WRITE THE CURRENT CONFORMATION!!
###        index_LC = cl_lengths.index(max(cl_lengths))
###        largest_equals_best = True
###        if index_LC!=0:
###            largest_equals_best = False
###        #WRITE THIS CONFORMATION!
###        origname = docked.ligMol.name
###        docked.ligMol.name = new_ligand_name
###        outfilename = new_ligand_name + ".pdbqt"
###        # number of atoms used for ligand efficiency is NOW correctly heavy atoms ONLY!!
###        total_num_lig_ats = len(docked.ligMol.allAtoms)
###        num_h_ats = len(docked.ligMol.allAtoms.get(lambda x: x.element=='H'))
###        num_lig_ats = total_num_lig_ats-num_h_ats
###        new_coords = docked.ligMol.allAtoms.coords
###        # setup file for writing
###        optr = open(outfilename, 'w')
###        sss, dr_str = docked.clusterer.getInfoStr(comment='USER  AD> ', ind=ind, rms=rmsd,ncl_to_write=max_cl_to_write, include_dlgfilename_run=True)
###        sss_list = sss.split('\n')
###        dr_str_list = dr_str.split('\n') #dlgfilename and run info
###        #look for omitted clusters
###        omitted = ""
###        if sss_list[0].find("omitted")>-1:
###            omitted = sss_list[0] 
###            sss_list = sss_list[1:]
###        num_sss = len(sss_list)
###        if num_sss == 1:
###            num_to_write = 1
###        else:
###            num_to_write = min(num_sss-1, max_cl_to_write)
###        b_efficiency = conf.binding_energy/num_lig_ats
###        optr.write( "%s %4.2f %3d runs %3d clusters\n"%(sss_list[0],rmsd, nconfs, nclusters))
###        optr.write( "USER  AD> ligand efficiency  %6.4f\n"%(b_efficiency))
###        optr.write( "USER  AD> rmsd, LE, clu_size, clu_e_range, dlgfilename, run#, b_curCRDs, b_LE, b_LC\n") #?CHANGE?
###        num_clusters = len(cl_dict)
###        if num_clusters < num_to_write:
###            num_to_write = num_clusters
###        for i in range(0, num_to_write):  #starts with '#binding' so numbers are 1-based
###            cl = cl_dict[i]
###            len_cl = len(cl)
###            #add energy_range
###            b_cl = cl[0]   #best in cluster
###            w_cl = cl[-1]  #worst in cluster
###            e_range = w_cl.binding_energy - b_cl.binding_energy
###            b_LE = int(docked.ch.current_conf==cl_dict[0][0])
###            b_LC = int(docked.ch.current_conf==cl_dict[index_LC][0])
###            b_CURRENT = int(docked.ch.current_conf==cl_dict[i][0])
###            if i==0: #BEST ENERGY
###                if index_LC!=0:
###                    optr.write( "%s,%d,%4.2f,%s,%d,1,0\n"%(sss_list[i+1], len_cl,  e_range, dr_str_list[i], b_LE))
###                else:
###                    optr.write( "%s,%d,%4.2f,%s,%d,1,1\n"%(sss_list[i+1], len_cl,  e_range, dr_str_list[i], b_LE))
###            elif i==index_LC: #here Largest Cluster is NOT also Lowest Energy
###                optr.write( "%s,%d,%4.2f,%s,0,0,%d\n"%(sss_list[i+1], len_cl, e_range, dr_str_list[i], b_LC))
###            elif b_CURRENT: # these coordinates are neither LE nor LE_LC
###                optr.write( "%s,%d,%4.2f,%s,1,0,0\n"%(sss_list[i+1], len_cl, e_range, dr_str_list[i]))
###            else: 
###                optr.write( "%s,%d,%4.2f,%s,0,0,0\n"%(sss_list[i+1], len_cl, e_range, dr_str_list[i]))
###        if len(omitted): #if some clusters have been skipped
###            optr.write(omitted + "\n")

###        #now write the rest of the ligand from its parser
###        atm_ct = 0
###        for line in docked.ligMol.parser.allLines:
###            #update the coordinates here...
###            if line.find("HETATM")==0 or line.find("ATOM")==0:
###                nl = line[:30] +  "%8.3f%8.3f%8.3f"%(docked.ligMol.allAtoms[atm_ct].coords[0], docked.ligMol.allAtoms[atm_ct].coords[1], docked.ligMol.allAtoms[atm_ct].coords[2])  + line[54:]
###                atm_ct += 1
###                optr.write(nl + "\n")        
###            else:    
###                optr.write(line + "\n")
###        optr.close()
###        docked.ligMol.name = origname 
###        #now write the stuff from intF
###        sss2 = self.intF.processLigand(outfilename, outputfilename=outfilename) 
###        #add the VSResult header line
###        optr = open(outfilename, 'r')
###        lines = optr.readlines()
###        optr.close()
###        optr = open(outfilename, 'w')
###        optr.write("REMARK VirtualScreeningResult %s\n"%(time.asctime( time.localtime(time.time()))))
###        for l in lines: 
###            optr.write(l)
###        optr.close()


ADWriteVSResultGUI=CommandGUI()
ADWriteVSResultGUI.addMenuCommand('AutoToolsBar', menuText['AnalyzeMB'], 
        menuText['writeVSResult'], cascadeName = menuText['MoleculesMB'])



class ADReadMacro(MVCommand):
    """This class is used to load the macromolecule specified in the docking log
    \nPackage : AutoDockTools
    \nModule  : autoanalyzeCommands
    \nClass   : ADReadMacro
    \nCommand : ADanalyze_readMacromolecule
    \nSynopsis:\n 
        None<-ADanalyze_readMacromolecule(macroFile=None)
    \nArguments:\n   
        macroFile: pdbqt file for macromolecule
    """


    def onAddCmdToViewer(self):
        checkHasInitializedDockings(self.vf)
        if self.vf.hasGui:
            for item in ['colorByMolecules','colorByAtomType', 'color']:
                if not hasattr(self.vf, item):
                    self.vf.loadCommand('colorCommands', item, 'Pmv')
            for item in ['displayLines','showMolecules']:
                if not hasattr(self.vf, item):
                    self.vf.loadCommand('displayCommands', item, 'Pmv')
                    

    def onRemoveObjectFromViewer(self, obj):
        if self.vf.docked and hasattr(self.vf.docked, 'macroMol'):
            d = self.vf.docked
            if d.macroMol==obj:
                delattr(d, 'macroMol')
            for v in self.vf.dockings.values():
                if hasattr(v,'macroMol') and v.macroMol==obj:
                    delattr(v,'macroMol')


    def guiCallback(self):
        """called each time the 'show macro' button is pressed"""
        if not self.vf.docked:
            self.vf.warningMsg('Please Read a Docking Log First')
            return
        if hasattr(self.vf.docked, 'macroMol'):
            #self.vf.warningMsg('current docking already has a macromolecule')
            #return
            t = 'Current docking already has a macromolecule:'+ self.vf.docked.macroMol.name+"\nDo you want to replace it?"
            d = SimpleDialog(self.vf.GUI.ROOT, text = t,
            buttons = ["Yes", "No", "Cancel"], default = 0, 
                    title = "Replace current macromolecule? ")
            replace = d.go()
            if replace>0:
                return 'ERROR'
        if self.vf.docked.__class__==Protein:
            macroFile = self.vf.askFileOpen(types=[('select receptor filename:', '*.pdbqt'),
                                ('all files','*')],
                                title = "AutoDock Vina Receptor File" )
            if not macroFile:
                return 'ERROR'
        elif len(self.vf.docked.dlo_list):
            macroFile = self.vf.docked.dlo_list[0].macroFile
            if macroFile==".pdbqs" or macroFile=='.pdbqt':
                msgStr = "unable to read macromolecule:\nmacroFile not specified in this docking log"
                self.vf.warningMsg(msgStr)
                return 'ERROR'
        else:
            msgStr = "current docking has no macromolecule!!!"
            self.vf.warningMsg(msgStr)
            return 'ERROR'
        fileList = glob.glob(macroFile)
        if macroFile not in fileList:
            #show warning
            msgStr = macroFile + " not found in this directory"
            self.vf.warningMsg(msgStr)
            msgName = macroFile + ':'
            macroFile = self.vf.askFileOpen(types=[('select macromolecule filename:', '*.pdbq*'),
                                ('all files','*')],
                                title = msgName)
            #if macroFile is None:
            if not macroFile:
                return 'ERROR'
        #self.doitWrapper(self.vf.docked.dlo_list[0].macroFile, log=1,redraw=1)
        self.doitWrapper(macroFile, log=1,redraw=1)

        

    def __call__(self, macroFile=None, **kw):
        """None<-ADanalyze_readMacromolecule(macroFile=None)
        \nmacroFile: pdbq* file for macromolecule
        """
        if macroFile is None:
            return 'ERROR'
        if not self.vf.docked:
            print "no current docking"
            return 'ERROR'
        if hasattr(self.vf.docked, 'macroMol'):
            self.vf.warningMsg('current docking already has a macromolecule')
            return 'ERROR'
        apply(self.doitWrapper,(macroFile,), kw)


    def doit(self, macroFile):
        if not self.vf.docked:
            print "no current docking"
            return 'ERROR'
        macrolist = self.vf.readMolecule(macroFile,topCommand=0)
        if macrolist is None:
            return "ERROR"
        macro = macrolist[0]
        macro.buildBondsByDistance()
        self.vf.docked.macroMol = macro
        self.vf.displayLines(macro, topCommand=0, redraw=1)


ADReadMacroGUI=CommandGUI()
ADReadMacroGUI.addMenuCommand('AutoToolsBar', menuText['AnalyzeMB'], 
        menuText['readMacro'], cascadeName = menuText['MoleculesMB'])


        
class ADEPDBMol(MVCommand):
    """This class is used to run AutoDock in the cmd mode, epdb a molecule in
    the context of the grids and attach energy values to its atoms
    \nPackage : AutoDockTools
    \nModule  : autoanalyzeCommands
    \nClass   : ADEPDBMol
    \nCommand : ADanalyze_epdbMolecule
    \nSynopsis:\n
        epdbList<-ADanalyze_epdbMolecule
    \nRequired Arguments:\n       
        epdbList  --- list of tuples of energies for the atoms
                (atomType, vdW_energy, estat_energy)
        \ndpfFile --- docking parameter file for evaluation
        \nmolFile --- pdbq file to be evaluated
    """


    def onAddCmdToViewer(self):
        self.parser = EpdbParser()
        for item in ['colorByMolecules','colorByAtomType', 'colorByProperty']:
            if not hasattr(self.vf, item):
                self.vf.loadCommand('colorCommands', item, 'Pmv')
        if not self.vf.colorMaps.has_key('rgb256'):
            mod = __import__("ViewerFramework")
            VFPath = mod.__path__[0]
            self.vf.loadColorMap(os.path.join(VFPath, "ColorMaps","rgb256_map.py"), 
                                    topCommand=0)         


    def guiCallback(self):
        """called each time the 'epdb mol' button is pressed"""
        #need to get dpf and molecule files...
        dpfFile= self.vf.askFileOpen(types=[('select dpf:', '*.dpf'), 
                                ('all files','*')],
                                title = 'Select Docking Parameter File')
        if dpfFile:
            molFile= self.vf.askFileOpen(types=[('select molecule filename:', '*.pdbq'), 
                                ('all files','*')],
                                title = 'Select molecule')
            if molFile:
                kw = {}
                kw['log'] = 1
                kw['redraw'] = 1
                apply(self.doitWrapper, (dpfFile, [molFile],),kw)


    def __call__(self, dpfFile=None, molList=None, **kw):
        """epdbList<-ADanalyze_epdbMolecule
        \nepdbList is list of tuples of energies for the atoms
                (atomType, vdW_energy, estat_energy)
        \ndpfFile: docking parameter file for evaluation
        \nmolFile: pdbq file to be evaluated
        """
        if dpfFile and molList:
            return apply(self.doitWrapper, (dpfFile, molList,), kw)
        else:
            return 'ERROR'


    def doit(self, dpfFile, molList, **kw):
        #FIX THIS: it assumes autodock3 is name of executable
        cmd = "autodock3 -p "+ dpfFile + " -c"
        #cmd = "echo '" + molList + "\n'|autodock3 -p "+ dpfFile + "-c"
        outptr, inptr, errptr = popen2.popen3(cmd,bufsize=0)
        inptr.flush()
        #should these be kept for any reason?
        lineList = []
        line = outptr.readline()
        ct = 0
        while find(line, 'NOW IN COMMAND MODE')<0  and find(line, 'sorry')<0:
            lineList.append(line)
            line = outptr.readline()
            ct = ct + 1
        if find(line,'sorry')>-1:
            print 'command mode error: ', line
            return 'ERROR'
        #now build + enter epdb cmd
        for molFile in molList:
            cmd1 = 'epdb '+ molFile + ' \n'
            inptr.write(cmd1)
            inptr.flush()

            ll2 = [line]
            line = outptr.readline()
            #FIX THIS: process could fail to get this far...
            while find(line, 'Torsional Free Energy')<0:
                ll2.append(line)
                line = outptr.readline()
            ll2.append(line)
            outptr.flush()
            inptr.flush()

            #try to do it a second time
            #inptr.write(cmd1)
            #inptr.flush()
            #ll3 = [line]
            #line = outptr.readline()
            #FIX THIS: process could fail to get this far...
            #while find(line, 'Torsional Free Energy')<0:
            #    ll3.append(line)
            #    line = outptr.readline()
            #ll3.append(line)

            #to close popen object
            #cmd2 = 'stop\n'
            ###TEST THIS!!!!
            #inptr.write(cmd2)
            #inptr.close()
            parser = self.parser
            parser._parse(ll2)
            #parse lines and build dictionary IEA 
            # keys= AtomType, NB+EE, NBE, EEE, PC
            #also totals
            #Estimated Free Energy of Binding, FDE, FIE, FIEL, TFE
            #OUTPUT parser info to python shell:
            print 'Estimated Free Energy:', parser.estFreeEnergy
            print 'Final Docked Energy:', parser.finalDockedEnergy
            print 'Final Intermolecular Energy:', parser.finalIntermolEnergy
            print 'Final Internal Energy:', parser.finalInternalEnergy
            print 'Torsional Free Energy:', parser.torsionalFreeEnergy
            mols = self.vf.readMolecule(molFile)
            mol = mols[0]
            mol.estFreeEnergy = parser.estFreeEnergy
            mol.finalDockedEnergy = parser.finalDockedEnergy
            mol.finalIntermolEnergy = parser.finalIntermolEnergy
            mol.finalInternalEnergy = parser.finalInternalEnergy
            mol.torsionalFreeEnergy = parser.torsionalFreeEnergy

            for i in range(len(mol.allAtoms)):
                at = mol.allAtoms[i]
                at.vdw_energy = parser.vdw_energies[i]
                at.estat_energy = parser.estat_energies[i]
                at.total_energy = parser.total_energies[i]
            maxi = max(parser.vdw_energies)
            mini = min(parser.vdw_energies)
            self.vf.colorByProperty(self.vf.Mols[-1].allAtoms, ('lines',),
                                    'vdw_energy',
                                    maxi=max(parser.vdw_energies),
                                    mini=min(parser.vdw_energies),
                                    propertyLevel='Atom',
                                    colormap='rgb256', log=0)
        cmd2 = 'stop\n'
        inptr.write(cmd2)
        inptr.flush()
        inptr.close()
        outptr.flush()
        outptr.close()
        return parser.estFreeEnergy
        #return map(None,parser.atTypes,parser.vdWs, parser.estats)


ADEPDBMolGUI=CommandGUI()
ADEPDBMolGUI.addMenuCommand('AutoToolsBar', menuText['AnalyzeMB'], 
        menuText['epdbMol'], cascadeName = menuText['GridsMB'])

        
class ADSeeSpots(MVCommand):
    """This class is used to display all docked conformations as spheres
    \nPackage : AutoDockTools
    \nModule  : autoanalyzeCommands
    \nClass   : ADSeeSpots
    \nCommand : ADanalyze_showDockingsAsSpheres 
    \nSynopsis:\n
        None<-ADanalyze_showDockingsAsSpheres
    \nRequired Arguments:\n
        docklist: list of keys into self.vf.dockings
    """

    def onAddCmdToViewer(self):
        # Create a ColorChooser.
        self.palette = ColorChooser(commands=self.color_cb,
                    exitFunction = self.hidePalette_cb)
        # pack
        self.palette.pack(fill='both', expand=1)
        # hide it 
        self.palette.hide()
        self.quality = 10
        self.managedGeometries = []


    def reset(self, event=None):
        #this function is called when docking shown as spheres
        #is in process of being cleared from viewer
        if hasattr(self, 'cmg'):
            self.cmg.master.withdraw()
            delattr(self, 'cmg')
            self.dismiss_cb()


    def hidePalette_cb(self, event=None):
        self.palette.hide()


    def guiCallback(self):
        """called each time the 'Visualize dockings as Spheres' button is pressed"""
        docklist = self.vf.dockings.keys()
        if len(docklist)==0:
            self.vf.warningMsg('Please Read a Docking Log First')
            return
        elif len(docklist)==1 and  hasattr(docklist[0],'ligMol') and hasattr(docklist[0].ligMol, 'vina_energy'):
            self.warningMsg("'Show As Spheres' is not implemented for vina results")
            return 
        else:
            self.doitWrapper(docklist,log=1,redraw=0)



    def __call__(self, docklist, **kw):
        """None<-ADanalyze_showDockingsAsSpheres
        \ndocklist: list of keys into self.vf.dockings
        """
        if docklist and len(docklist):
            return apply(self.doitWrapper, (docklist,), kw)
        else:
            return 'ERROR'


    def doit(self, docklist, **kw):
        self.keys = docklist
        res = self.addSpheres(docklist)
        if res=='ERROR':
            return res
        if not hasattr(self, 'form'):
            ifd = self.getIfd(docklist)
            self.form = self.vf.getUserInput(self.ifd, modal=0, blocking=0,
                 width=380, height=390)
            self.form.root.protocol('WM_DELETE_WINDOW',self.dismiss_cb)
            self.setUpWidgets()
        else:
            self.lb.delete(0,'end')
            for item in docklist:
                self.lb.insert('end', item)
                gC = self.vf.dockings[item].ligMol.geomContainer
                g = gC.geoms['dockedSpheres']
                if g.visible:
                    self.lb.select_set('end')
            self.form.deiconify()
        self.updateVis()
                

    def setUpWidgets(self):
        self.radii_esw = self.ifd.entryByName['radii']['widget']
        self.lb = self.ifd.entryByName['dockedLC']['widget'].lb
        self.quality_sl = self.ifd.entryByName['quality']['widget']
        self.quality_sl.set(self.quality)


    def getIfd(self, docklist):
        #print 'docklist=', docklist
        dl2 = []
        for n in docklist:
            dl2.append((n, None))
        ifd = self.ifd = InputFormDescr(title='Show Docked Conformations as Spheres')
        ifd.append({'name': 'dockingLabel',
            'widgetType': Tkinter.Label,
            'wcfg':{'text': str(len(docklist)) + ' Current Docking(s)\n(centers=average of coords)'},
            'gridcfg':{'sticky': 'wens', 'columnspan':2}})
        ifd.append({'name': 'dockedLC',
            'widgetType':ListChooser,
            'wcfg':{
                'entries': dl2,
                'mode': 'multiple',
                'title': '',
                'command': self.updateVis,
                #'command': CallBackFunction(self.addSpheres, docklist),
                'lbwcfg':{'height':5, 
                    'selectforeground': 'red',
                    'exportselection': 0,
                    'width': 30},
            },
            'gridcfg':{'sticky':'wens', 'row':2,'column':0, 'columnspan':2}})
        ifd.append( {'name':'radii',
            'widgetType': ExtendedSliderWidget,
            'wcfg':{'label':'radii',
                 'minval':.01, 'maxval':1.0,
                'immediate':1,
                'init':.2,
                'width':250,
                'command':self.updateRadii,
                'sliderType':'float',
                'entrypackcfg':{'side':'right'}},
            'gridcfg':{'sticky':'we', 'columnspan':2}})
            #'gridcfg':{'sticky':'wens', 'columnspan':2}})
        ifd.append( {'name':'quality',
            'widgetType': Tkinter.Scale,
            'wcfg':{'label':'quality',
                'troughcolor':'green',
                 'from_':2, 'to_':20,
                'orient':'horizontal',
                'length':'2i',
                'command':self.updateQuality },
            'gridcfg':{'sticky':'we', 'columnspan':2}})
            #'gridcfg':{'sticky':'wens'}})
        ifd.append({'name':'updateColorBut',
            'widgetType':Tkinter.Button,
            'wcfg': { 'text':'Choose color',
                'command': self.updateColor},
            'gridcfg':{'sticky':'we' }})
            #'gridcfg':{'sticky':'wes', 'row':-1, 'column':1}})
        ifd.append({'name':'colorByPropBut',
            'widgetType':Tkinter.Button,
            'wcfg': { 'text':'Color by Prop',
                'command': self.colorByProp},
            'gridcfg':{'sticky':'we', 'row':-1, 'column':1}})
        ifd.append({'name':'closeBut',
            'widgetType':Tkinter.Button,
            'wcfg': { 'text':'Dismiss',
                'command': self.dismiss_cb},
            'gridcfg':{'sticky':'wens','columnspan':2 }})
        return ifd


    def colorByProp(self, event=None):
        curSel = []
        props = []
        for item in self.lb.curselection():
            curSel.append(self.lb.get(item))
        for k in self.keys:
            if k in curSel:
                docking = self.vf.dockings[k]
                gC = docking.ligMol.geomContainer
                props = docking.sphere_energies
                cm = self.vf.colorMaps['rgb256']
                minval = min(props)
                maxval = max(props)
                cm.configure(mini=minval, maxi=maxval)
                cols = cm.Map(props, mini=minval, maxi=maxval)
                g = gC.geoms['dockedSpheres']
                if g.visible:
                    #cols has shape (#dockings,4)
                    g.Set(radii=self.radii_esw.get(), materials=cols, 
                            inheritMaterial=0)
        if len(props) and not hasattr(self, 'cmg'):
            self.vf.showCMGUI('rgb256')
            self.cmg = self.vf.showCMGUI.cmg
            self.cmg.addCallback(self.colorByProp)
            self.cmg.configure(mini=minval, maxi=maxval)
            self.cmg.docking = docking
        self.vf.GUI.VIEWER.Redraw()
 

    def updateVis(self, event=None):
        curSel = []
        for item in self.lb.curselection():
            curSel.append(self.lb.get(item))
        for k in self.keys:
            gC = self.vf.dockings[k].ligMol.geomContainer
            if gC.geoms.has_key('dockedSpheres'):
                g = gC.geoms['dockedSpheres']
            else:
                res = self.addSpheres(self.keys)
                if res=='ERROR':
                    return res
            g.Set(visible = k in curSel)
            #g.RedoDisplayList()
        self.vf.GUI.VIEWER.Redraw()


    def updateRadii(self, event=None):
        curSel = []
        for item in self.lb.curselection():
            curSel.append(self.lb.get(item))
        for k in self.keys:
            if k in curSel:
                gC = self.vf.dockings[k].ligMol.geomContainer
                g = gC.geoms['dockedSpheres']
                g.Set(radii= self.radii_esw.get())
        self.vf.GUI.VIEWER.Redraw()
 

    def updateQuality(self, event=None):
        self.quality = self.quality_sl.get()
        curSel = []
        for item in self.lb.curselection():
            curSel.append(self.lb.get(item))
        for k in self.keys:
            if k in curSel:
                gC = self.vf.dockings[k].ligMol.geomContainer
                g = gC.geoms['dockedSpheres']
                g.Set(quality=self.quality)
        self.vf.GUI.VIEWER.Redraw()


    def color_cb(self, colors):
        for k in self.keys:
            gC = self.vf.dockings[k].ligMol.geomContainer
            g = gC.geoms['dockedSpheres']
            if g.visible:
                g.Set(radii = self.radii_esw.get(), materials = (colors[:3],),
                            inheritMaterial=0)
        self.vf.GUI.VIEWER.Redraw()
        col = colorTool.TkColor(colors[:3])
        self.quality_sl.config(troughcolor = col)


    def updateColor(self, event=None):
        self.palette.master.deiconify()


    def dismiss_cb(self, event=None):
        self.form.withdraw()
        if hasattr(self, 'palette'):
            self.palette.hide()


    def pickedVerticesToDockings(self, geom, vertIndList):
        """Function called to convert picked vertices into dockings"""
        #print "in pickedVerticesToDockings: geom=", geom, ' + vertIndList=', vertIndList
        ind_list = self.ifd.entryByName['dockedLC']['widget'].lb.curselection()
        key_ind = int(ind_list[0])
        print 'highlighted dockinglog filename is ', self.keys[key_ind]
        docking = self.vf.dockings[self.keys[key_ind]]
        geomCont = docking.ligMol.geomContainer
        msg = "Picked sphere(s) corresponds to docking run(s) with energy :\n"
        for ind in vertIndList:
            lig = docking.ligMol
            if not hasattr(lig,'spw'):
                lig.spw = ConformationPlayer(lig, docking, self.vf, 
                        lig.name, docking.ch.conformations, form2=1)
                        #buttonMask={'recordB':False})
            conf = docking.ch.conformations[ind]
            docking.ligMol.spw.nextFrame(ind+1)
            msg += "%4d  % 6.4f\n" %(ind+1, docking.sphere_energies[ind])
        print msg
        if self.cmg.legend is None: self.cmg.createCML()
        self.cmg.showLegendVar.set(1)
        self.cmg.configure({},updateGui=True)
        self.cmg.legend.Set(visible=1, unitsString=str(conf.energy)+ "kcal",
                            name=str(ind+1)+'th run') 
        return AtomSet()

    
    def addSpheres(self, docklist, event=None):
        for k in self.keys:
            vis = k in docklist
            docking = self.vf.dockings[k]
            if hasattr(docking.ligMol, 'vina_energy'):
                self.warningMsg("'Show as Spheres' is not implemented for vina results")
                return 'ERROR'
            geomCont = docking.ligMol.geomContainer
            if geomCont.geoms.has_key('dockedSpheres'):
                geomCont.geoms['dockedSpheres'].Set(visible=vis)
                self.vf.GUI.VIEWER.Redraw()
                continue
            verts = self.getVertices(docking)
            docking.sphere_energies = self.getEnergies(docking)
            sphs = Spheres('dockedSpheres', 
                    materials=((0,1,0),), shape=(0,3), visible=vis,
                    radii=0.1, quality=15, 
                    inheritMaterial=0)
                    #pickable=0, inheritMaterial=0)
            self.managedGeometries.append(sphs)
            geomCont.addGeom(sphs)
            geomCont.geomPickToAtoms['dockedSpheres'] = self.pickedVerticesToDockings
            geomCont.geomPickToBonds['dockedSpheres'] = None
            #geomCont.addGeom(sphs, self)
            sphs.replace = True
            self.managedGeometries.append(sphs)
            geomCont.VIEWER.AddObject(sphs, parent=geomCont.masterGeom, redo=0)
            sphs.Set(vertices = verts, visible=0)


    def getVertices(self, docking):
        verts = []
        num = len(docking.ligMol.allAtoms)*1.0
        if hasattr(docking.ligMol, 'vina_energy'):
            self.warningMsg("Showing results as spheres is not implemented for vina results")
            return 'ERROR'
        if docking.flex_res:
            ind = len(docking.flex_res.atoms)
            num = num - ind
            ind = int(ind)
            for c in docking.ch.conformations:
                center = Numeric.add.reduce(Numeric.array(c.getCoords()[:-ind])).astype('f')
                this_vert = (center/(num)).tolist()
                verts.append(this_vert)
        else:
            for c in docking.ch.conformations:
                center = Numeric.add.reduce(Numeric.array(c.getCoords())).astype('f')
                this_vert = (center/(num)).tolist()
                verts.append(this_vert)
        return verts


    def getEnergies(self, docking):
        if hasattr(docking.ligMol, 'vina_energy'):
            self.warningMsg("Showing energies is not implemented for vina results")
            return 'ERROR'
        energies = []
        for c in docking.ch.conformations:
            if hasattr(c, 'energy'):
                energies.append(c.energy)
            elif hasattr(c, 'total_energies'):
                this_e = Numeric.add.reduce(c.total_energies)
                energies.append(this_e)
                c.energy = this_e
            else:
                self.warningMsg("NO ENERGIES AVAILABLE")
                energies.append(0)
        return energies


ADSeeSpotsGUI=CommandGUI()
ADSeeSpotsGUI.addMenuCommand('AutoToolsBar', menuText['AnalyzeMB'], 
        menuText['seeSpots'], cascadeName = menuText['DockingLogMB'], separatorAbove=1)


        
class ADShowBindingSite(MVCommand):
    """This class is used to display docked ligand and surrounding atoms in macromolecule
    \nPackage : AutoDockTools
    \nModule  : autoanalyzeCommands
    \nClass   : ADShowBindingSite
    \nCommand : ADanalyze_showBindingSite 
    \nSynopsis:\n
        None<-ADanalyze_showBindingSite
    \nRequired Arguments:\n
        docking: autodock result to display
    """

    def onAddCmdToViewer(self):
        self.vf.loadModule('labelCommands', 'Pmv')
        self.vf.loadModule('displayCommands', 'Pmv')
        self.vf.loadModule('hbondCommands', 'Pmv')
        self.percentCutoff = 1.0
        self.updateDisplay = Tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.setBKGColor = Tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.showHbatSpheres = Tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.showCloseContactSpheres = Tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.showSS = Tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.showMM = Tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.showMsms = Tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.showResLabels = Tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.showPiPi = Tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.showPiCation = Tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.infoType = Tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.infoType.set("")
        self.infoTypeList = ["ligand contacts", "receptor contacts", "hydrogen bonds"]
        self.sphereType = Tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.sphereType.set("wireframe")
        self.sphereTypeList = ["wireframe", "solid"]
        self.backgroundcolor=[1.0,1.0,1.0]
        self.managed_sets = {}  #current sets to display
        for v in [self.updateDisplay, self.showCloseContactSpheres, self.showSS, \
                    self.showMsms, self.showResLabels]:
            v.set(1)
        self.showSS.set(0)  #turned off 7/2009 due to problems with vina results
        self.masterGeom = Geom('AD_displayIntGeom', shape=(0,0), 
                                pickable=0, protected=True)
        self.masterGeom.isScalable = 0
        display_geoms_list  = self.vf.GUI.VIEWER.findGeomsByName('AD_display_geoms')
        if display_geoms_list==[]:
            display_geoms = Geom("AD_display_geoms", shape=(0,0), protected=True)
            self.vf.GUI.VIEWER.AddObject(display_geoms, parent=self.vf.GUI.miscGeom)
            display_geoms_list = [display_geoms]
        display_geoms = display_geoms_list[0]
        display_int_geom_list  = self.vf.GUI.VIEWER.findGeomsByName('AD_displayInteractions_geom')
        if display_int_geom_list==[]:
        	display_int_geom = Geom('AD_displayInteractions_geoms', shape=(0,0), 
                                pickable=0, protected=True)
        	display_int_geom.isScalable = 0
        	self.vf.GUI.VIEWER.AddObject(display_int_geom, parent = display_geoms)
            	display_int_geom_list = [display_int_geom]
        display_int_geom = display_int_geom_list[0]
        pi_pi_geom_list = self.vf.GUI.VIEWER.findGeomsByName("pi_pi")
        if pi_pi_geom_list==[]:
            pi_pi_geom = Cylinders('pi_pi',quality=50,
                            inheritLineWidth=0, lineWidth=10,
                            radii=(0.7), pickable=0,
                            materials = ((1,1,0),),
                            inheritMaterial=False)
            pi_pi_geom_list = [pi_pi_geom]
        self.pi_pi_geom = pi_pi_geom_list[0]
        t_shaped_geom_list = self.vf.GUI.VIEWER.findGeomsByName("t_shaped")
        if t_shaped_geom_list==[]:
            t_shaped_geom = Cylinders('t_shaped',quality=50,
                            inheritLineWidth=0, lineWidth=10,
                            radii=(0.7), pickable=0,
                            materials = ((1,0,0),),
                            inheritMaterial=False)
            t_shaped_geom_list = [t_shaped_geom]
        self.t_shaped_geom = t_shaped_geom_list[0]
        cation_pi_geom_list = self.vf.GUI.VIEWER.findGeomsByName("cation_pi")
        if cation_pi_geom_list ==[]: 
            cation_pi_geom = Cylinders('cation_pi',quality=50,
                            inheritLineWidth=0, lineWidth=10,
                            radii=(0.2, 0.7), pickable=0,
                            materials = ((1,1,0),),
                            inheritMaterial=False)
            cation_pi_geom_list ==[cation_pi_geom] 
        self.cation_pi_geom = cation_pi_geom_list[0]
        pi_cation_geom_list = self.vf.GUI.VIEWER.findGeomsByName("pi_cation")
        if pi_cation_geom_list ==[]: 
            pi_cation_geom = Cylinders('pi_cation',quality=50,
                            inheritLineWidth=0, lineWidth=10,
                            radii=((0.7, 0.01),), pickable=0,
                            materials = ((1,1,0),),
                            inheritMaterial=False)
            pi_cation_geom_list ==[pi_cation_geom] 
        self.pi_cation_geom = pi_cation_geom_list[0]
        from opengltk.OpenGL import GL
        for obj in [self.pi_pi_geom, self.t_shaped_geom, self.cation_pi_geom,
                self.pi_cation_geom]:
            try:
            	self.vf.GUI.VIEWER.AddObject(obj, parent=display_int_geom)
            except:
                pass
            obj.Set(inheritFrontPolyMode = 0, 
                    frontPolyMode=GL.GL_LINE,
                    opacity=(.3,.3), inheritLineWidth=0,
                    inheritCulling=0, lineWidth=1)


    def update_managed_sets(self, event=None):
        rDict = self.intDescr.results
        key_list = [('macro_close_res','mResLab'),   ('ss_res','mSSRes'),\
                    ('macro_close_non_hb','mClose'), ('macro_hb_atoms','mhbnds'), \
                    ('lig_close_non_hb','lClose'), ('lig_hb_atoms','lhbnds'), \
                    ('lig_close_carbons', 'lcloseCs')]
        for k1,k2 in key_list:
            nodes = rDict[k1]
            if k2 in self.vf.sets.keys():
                self.vf.sets.pop(k2)
                self.vf.saveSet(nodes, k2, addToDashboard=False)
            else:
                self.vf.saveSet(nodes, k2, addToDashboard=True)


    def onRemoveObjectFromViewer(self, obj):
        if self.vf.docked and hasattr(self.vf.docked, 'ligMol') and obj==self.vf.docked.ligMol:
            self.reset()
            for v in [self.updateDisplay, self.showCloseContactSpheres, self.showSS, \
                    self.showMsms, self.showResLabels]:
                v.set(1)
            if hasattr(self, 'displayform'):
                self.close()
            try:
                delattr(self, 'docking')
            except:
                pass



    def getSSResidues(self, mol):
        from MolKit.protein import ResidueSet
        res = ResidueSet()
        for k, v in mol.geomContainer.atoms.items():
	        if k[:3] in ['Coi', 'Hel', 'Str', 'Tur']:
		        res += v
        return res


    def setupUndoBefore(self, docking):
        #???color???
        if type(docking)==types.StringType:
            docking = self.vf.dockings[docking]
        if not hasattr(docking, 'macroMol'):
            return 
        if hasattr(docking, 'has_ribbon'):
            delattr(docking, 'has_ribbon')
        lig = docking.ligMol
        macro = docking.macroMol
        kw = {'log':0,'topCommand':0, 'redraw':True}
        kw['lineWidth'] = self.vf.displayLines.lastUsedValues['default']['lineWidth']
        macro_gca = macro.geomContainer.atoms
        geomSet = macro_gca['bonded']
        boSet =  macro_gca['bondorder'].uniq()
        if len(boSet) == 0:
            kw['displayBO']=False
        else:
            kw['displayBO']=True
        # The undo of a display command is to display ONLY what was
        # displayed before, if something was already displayed
        self.addUndoCall( (geomSet,), kw, self.vf.displayLines.name)
        kw = {'log':0,'topCommand':0, 'redraw':True}
        old_colors = macro.allAtoms.colors['lines']
        self.addUndoCall( (geomSet,old_colors,), kw, self.vf.color.name)
        geomSet = macro_gca['ResidueLabels']
        kw = self.vf.labelByProperty.lastUsedValues['default']
        if len(geomSet):
            kw['only'] = 1
        else:
            geomSet = macro.findType(Residue)
            kw['negate'] = 1
        self.addUndoCall( (geomSet,), kw, self.vf.labelByProperty.name)
        kw = {'log':0,'topCommand':0, 'redraw':True}
        geomSet = macro_gca['cpk']
        if len(geomSet):
            kw['only'] = 1
        else:
            geomSet = macro.allAtoms
            kw['negate'] = 1
        #self.addUndoCall( (geomSet,), kw, self.vf.displayCPK.name)
        #MACRO SS
        res = self.getSSResidues(macro)
        if len(res):
            kw = {'log':0,'topCommand':0, 'redraw':True, 'negate':1}
            self.addUndoCall( (res,), kw, self.vf.ribbon.name)
        #kw = {'log':0,'topCommand':0, 'redraw':True, 'negate':1}
        #self.addUndoCall( (macro.allAtoms,), kw, self.vf.ribbon.name)
        #ligand
        lgca = lig.geomContainer.atoms
        defaultValues = self.vf.displaySticksAndBalls.lastUsedValues['default']
        kw={}
        kw['bRad'] = defaultValues['bRad']
        kw['bScale'] = defaultValues['bScale']
        kw['bquality'] = defaultValues['bquality']
        kw['cradius'] = defaultValues['cradius']
        kw['cquality'] = defaultValues['cquality']
        kw['sticksBallsLicorice'] = defaultValues['sticksBallsLicorice']
        ballset = lgca['balls']
        stickset = lgca['sticks']
        if len(ballset)==0:
            # no balls displayed
            if len(stickset) == 0: # negate
                # no sticks displayed 
                kw['negate'] = True
                kw['redraw'] = True
                self.addUndoCall( (lig.allAtoms,), kw,
                                  self.vf.displaySticksAndBalls.name )
            else:
                # noballs was on
                kw['negate'] = False
                kw['redraw'] = True
                kw['only'] = True
                self.addUndoCall( (stickset,), kw,
                                  self.vf.displaySticksAndBalls.name )
        else:
            kw['redraw'] = True
            kw['only'] = True
            self.addUndoCall( (stickset,), kw,
                              self.vf.displaySticksAndBalls.name )
        lgcg = lig.geomContainer.geoms
        lgca = lig.geomContainer.atoms
        surfName = lig.name + '_msms'
        self.addUndoCall( (lig,), {'surfName':surfName, 'negate':1}, self.vf.displayMSMS.name)
        kw = {'log':0,'topCommand':0, 'redraw':True}
        geomSet = lgca['cpk']
        if len(geomSet):
            #undisplay all
            kw['negate'] = 1
            self.addUndoCall( (lig.allAtoms,), kw, self.vf.displayCPK.name)
            #redisplay current only...
            kw['negate'] = 0
            self.addUndoCall( (geomSet,), kw, self.vf.displayCPK.name)
        else:
            geomSet = lig.allAtoms
            kw['negate'] = 1
            self.addUndoCall( (geomSet,), kw, self.vf.displayCPK.name)
        #extrudedSS
        keys_to_exclude = ['sticks', 'AtomLabels', 'CAsticks', 'bonded', 'cpk', 'ChainLabels',\
                'ProteinLabels', 'ResidueLabels', 'nobnds', 'balls', 'bondorder', 'lines', 'CAballs']
        displayed_res = ResidueSet()
        gc = macro.geomContainer
        for k in gc.atoms.keys():
            if not k in keys_to_exclude:
                displayed_res += gc.atoms[k].findType(Residue)
        if len(displayed_res):
            kw = {'log':0,'topCommand':0, 'only':True, 'redraw':True}
            self.addUndoCall( (displayed_res,), kw, self.vf.displayExtrudedSS.name)
        else:
            kw = {'log':0,'topCommand':0, 'redraw':True, 'negate':True}
            self.addUndoCall( (macro.chains.residues,), kw, self.vf.displayExtrudedSS.name)
        #        self.vf.displayExtrudedSS(macro_atoms, topCommand=0 )
        #flex_res
        #try to undo flex_res display also
        if hasattr(docking, 'flex_res') and len(docking.flex_res):
            ind = len(lig.chains.residues) - len(docking.flex_res)
            flex_res = lig.chains.residues[ind:]
            ##try to set the label fontScales here
            #g = lig.geomContainer.geoms['ResidueLabels']
            #g.Set(fontScales=(.3,.3,.3))
            kw = {'log':0,'topCommand':0, 'redraw':True}
            kw['lineWidth'] = self.vf.displayLines.lastUsedValues['default']['lineWidth']
            geomSet = lgca['bonded']
            boSet =  lgca['bondorder'].uniq()
            if len(boSet) == 0:
                kw['displayBO']=False
            else:
                kw['displayBO']=True
            # The undo of a display command is to display ONLY what was
            # displayed before, if something was already displayed
            self.addUndoCall( (geomSet,), kw, self.vf.displayLines.name)
            kw = {'log':0,'topCommand':0, 'redraw':True}
            old_colors = lig.allAtoms.colors['lines']
            self.addUndoCall( (geomSet,old_colors,), kw, self.vf.color.name)
            geomSet = lgca['ResidueLabels']
            kw = self.vf.labelByProperty.lastUsedValues['default']
            if len(geomSet):
                kw['only'] = 1
            else:
                geomSet = lig.findType(Residue)
                kw['negate'] = 1
            self.addUndoCall( (geomSet,), kw, self.vf.labelByProperty.name)
            kw = {'log':0,'topCommand':0, 'redraw':True}
            geomSet = lgca['cpk']
            if len(geomSet):
                #undisplay all
                kw['negate'] = 1
                self.addUndoCall( (lig.allAtoms,), kw, self.vf.displayCPK.name)
                #redisplay current only...
                kw['negate'] = 0
                self.addUndoCall( (geomSet,), kw, self.vf.displayCPK.name)
                #kw['only'] = 1
            else:
                geomSet = lig.allAtoms
                kw['negate'] = 1
            self.addUndoCall( (geomSet,), kw, self.vf.displayCPK.name)
            #LIG SS
            res = self.getSSResidues(lig)
            if len(res):
                kw = {'log':0,'topCommand':0, 'redraw':True, 'negate':1}
                self.addUndoCall( (res,), kw, self.vf.ribbon.name)
        #self.vf.setbackgroundcolor
        old_color = self.vf.GUI.VIEWER.currentCamera.backgroundColor
        if len(old_color)>3: old_color = list(old_color[:3])
        self.addUndoCall( [old_color,], {}, self.vf.setbackgroundcolor.name)
        self.addUndoCall( (), {}, 'ADanalyze_showBindingSite.restoreCpk')
        self.addUndoCall( (), {}, 'hbondsAsSpheres.dismiss_cb')
        self.addUndoCall( (), {}, 'ADanalyze_showBindingSite.close')
        self.addUndoCall( (), {}, 'setbackgroundcolor.dismiss_cb')


    def setBackGroundColor(self, event=None):
        self.vf.setbackgroundcolor.guiCallback()
        self.backgroundcolor = self.vf.GUI.VIEWER.currentCamera.backgroundColor
        self.setBKGColor.set(0)


    def updateCpk(self):
        docking = self.docking
        macro = docking.macroMol
        lig = docking.ligMol
        m_cpk = macro.geomContainer.geoms['cpk']
        l_cpk = lig.geomContainer.geoms['cpk']
        #self.vf.GUI.VIEWER.CurrentCameraBackgroundColor((1.00,1.00,1.00))
        #self.vf.setbackgroundcolor(self.backgroundcolor, topCommand=0)
        #self.vf.setbackgroundcolor([1.00,1.00,1.00], topCommand=0)
        from opengltk.OpenGL import GL
        for g in [m_cpk, l_cpk]:
            g.inheritFrontPolyMode = 0
            if self.sphereType.get()=='wireframe':
                g.Set(frontPolyMode=GL.GL_LINE)
                g.Set(opacity=[.3,.3])
                g.Set(lineWidth=1, inheritLineWidth=0)
                g.Set(stacks=25, slices=25)
            else:
                g.Set(frontPolyMode=GL.GL_FILL)
                #g.Set(opacity=[1., 1.])


    def restoreCpk(self, **kw ):
        docking = self.vf.docked
        macro = docking.macroMol
        lig = docking.ligMol
        from opengltk.OpenGL import GL
        m_cpk = macro.geomContainer.geoms['cpk']
        l_cpk = lig.geomContainer.geoms['cpk']
        for g in [m_cpk, l_cpk]:
            g.inheritFrontPolyMode = 0
            if self.sphereType.get()=='solid':
                g.Set(frontPolyMode=GL.GL_FILL)
                #g.Set(opacity=[1., 1.])
            else:
                g.Set(frontPolyMode=GL.GL_LINE)
                g.Set(opacity=[.3,.3])
                g.Set(lineWidth=1, inheritLineWidth=0)
                g.Set(stacks=25, slices=25)


    def guiCallback(self):
        """called each time the 'Show Interactions' button is pressed"""
        docklist = self.vf.dockings.keys()
        if len(docklist)==0:
            self.vf.warningMsg("Currently no dockings are available. Please use Analyze->Open to read in a '.dlg' file first...")
            return
        elif len(docklist)==1:
            self.doitWrapper(docklist[0])
        else:
            ifd2= InputFormDescr(title='Select Docking: ')
            ifd2.append({'widgetType':'ListChooser',
                'name':'dlgFiles',
                'entries':docklist,
                'mode':'single',
                'wcfg':{'title':'select docking'},
                'lbwcfg':{'height':4},
                'gridcfg':{'sticky':'w', 'column':100,
                'rowspan':10}})
            val = self.vf.getUserInput(ifd2)
            if len(val)>0 and len(val['dlgFiles'])>0:
                docking = val['dlgFiles'][0]
                kw = {}
                apply(self.doitWrapper, val['dlgFiles'][0], kw)
            else:
                return "ERROR"


    def __call__(self, docking_filename, **kw):
        """None<-ADanalyze_showBindingSite
        \ndocking: 
        """
        if docking_filename:
            return apply(self.doitWrapper, (docking_filename,), kw)
        else:
            return 'ERROR'


    def revert(self, event=None, **kw):
        #try:
        #    self.vf.undo()
        #    delattr(self, 'docking')
        #except:
        #    pass
        if hasattr(self.docking.ligMol, 'spw'):
            self.docking.ligMol.spw.updateBindingSite = 0
        count = 0
        try:
            s = self.vf.undoCmdStack[-1]
            while len(s)>1 and s[1]!=self.name and count<10:
                self.vf.undo()
                count = count + 1
                s = self.vf.undoCmdStack[-1]
            self.vf.undo()
            self.form.withdraw()
            self.close()
        except:
            pass


    def close(self, event=None, **kw):
        for g in [self.pi_pi_geom, self.cation_pi_geom, 
                        self.pi_cation_geom]:
                g.Set(visible=0)
        self.displayform.destroy()
        delattr(self, 'displayform')


    def save(self, filename=None, event=None):
        self.vf.saveImage.guiCallback()


    def updatePercentCutoff(self, event=None):
        percentCutoff = self.ifd.entryByName['vdwPercentCutoff']['widget'].get()
        msg = 'resetting percentCutoff to ' + str(percentCutoff)
        self.warningMsg(msg)
        self.build(percentCutoff)
        
     
    def toggleUpdate(self, event=None):
        if hasattr(self.docking.ligMol, 'spw'):
            self.docking.ligMol.spw.updateBindingSite = self.updateDisplay.get()


    def buildDisplayForm(self, event=None):
        ifd = self.ifd= InputFormDescr(title='Binding Site Interactions display options')
        ifd.append({'name': 'updateCB',
                    'widgetType':Tkinter.Checkbutton,
                    'tooltip': 'Each time you change the conformation using the player,\nthe display will be rebuilt. To undo use\nRevert or Edit->Undo ADanalyze_showBindingSite.\nIf you have displayed multiple conformations,\npossibly you may need to Undo displayExtrudedSS multiple times..',
                    'wcfg':{'text':"update display for each new conformation",
                            'variable':self.updateDisplay},
                    'command': self.toggleUpdate,
                    'gridcfg':{'sticky':Tkinter.W, 'columnspan':2}})
        ifd.append({'name':'setBGcolorCB',
                    'widgetType':Tkinter.Checkbutton,
                    'tooltip': 'customize the Viewer background color...',
                    'wcfg':{'text':"set background color",
                            'variable':self.setBKGColor},
                    'command': self.setBackGroundColor,
                    'gridcfg':{'sticky':Tkinter.W, 'row':-1, 'column':2}})
        ifd.append({'name': 'msmsCB',
                    'widgetType':Tkinter.Checkbutton,
                    'tooltip': 'display solvent-excluded-surface geometry for ligand',
                    'wcfg':{'text':"display msms",
                            'variable':self.showMsms},
                    'command': self.updateMsms,
                    'gridcfg':{'sticky':Tkinter.W, 'columnspan':3}})
        ifd.append({'widgetType':Tkinter.Label,
                    'tooltip': 'set which spheres to display and whether to use wireframe or solid',
                    'wcfg': {'text':"display spheres as"},
                    'gridcfg':{'sticky':'w'} })
        ifd.append({'widgetType':Pmw.ComboBox,
            'name':'sphere_type',
            'wcfg':{'entryfield_value':self.sphereType.get(),
                    'labelpos':'w',
                    'listheight':'50',
                    'scrolledlist_items': self.sphereTypeList,
                    'selectioncommand': self.updateSphereDisplay,
                    },
            'gridcfg':{'sticky':'w', 'row':-1, 'column':1}})
        ifd.append({'widgetType':Tkinter.Label,
                    'tooltip': 'set which spheres to display ',
                    'wcfg': {'text':"on atoms in:"},
                    'gridcfg':{'sticky':'e'} })
        ifd.append({'name': 'closeAtsCB',
                    'widgetType':Tkinter.Checkbutton,
                    'tooltip': "display spheres on pairs of atoms closer than sum of vdw radii*Scaling Factor.\nWhen msms surface for 'ligand' is displayed, spheres are displayed only on the 'receptor'",
                    'wcfg':{'text':"close contact",
                            'variable':self.showCloseContactSpheres},
                    'command': self.updateCCCpk,
                    'gridcfg':{'sticky':Tkinter.W, 'row':-1, 'column':1}})
        ifd.append({'name': 'hbondsCB',
                    'widgetType':Tkinter.Checkbutton,
                    'tooltip': "display spheres on atoms involved in hydrogen bonds\nbetween ligand atoms and receptor atoms.\nWhen msms surface for ligand is displayed, spheres are displayed only on the receptor\n(Hbonds are displayed as small green spheres)",
                    'wcfg':{'text':"hydrogen bonds",
                            'variable':self.showHbatSpheres},
                    'command': self.updateHBonds,
                    'gridcfg':{'sticky':Tkinter.W, 'row':-1, 'column':2}})
        ifd.append({'name':'vdwPercentCutoff',
                    'widgetType':ThumbWheel,
                    'tooltip': "used in finding 'close atoms' which are closer than:\ndistance = (atom1.vdw+atom2.vdw)*scaling_factor\n(smaller values find fewer contacts; larger values more)",
                    'wcfg':{'labCfg':{'text': 'VDW Scaling Factor'},
                             'showLabel':1, 'width':150,
                             'width':75, 'height':20,
                             'min':0.1, 'type':float, 'precision':2,
                             'callback': self.updatePercentCutoff,
                             'value':1.0, 'continuous':0, 'oneTurn':2,
                             'wheelPad':2, 'height':20},
                    'gridcfg':{'sticky':'e', 'columnspan':3}})
        ifd.append({'widgetType':Tkinter.Label,
                    'tooltip': 'set residues for which to display extruded secondary structure',
                    'wcfg': {'text':"display ribbon for:"},
                    'gridcfg':{'sticky':'w'} })
        ifd.append({'name': 'ssCB',
                    'widgetType':Tkinter.Checkbutton,
                    'tooltip': 'display ribbon secondary structure for sequences of >3\nresidues in receptor with atoms close to ligand atoms.\ngaps of 1 residue are allowed',
                    'wcfg':{'text':"near residues",
                            'variable':self.showSS},
                    'command': self.showHideSecondaryStructure, 
                    'gridcfg':{'sticky':'e'}})
        ifd.append({'name': 'ssALLCB',
                    'widgetType':Tkinter.Checkbutton,
                    'tooltip': 'display ribbon secondary structure/nfor all residues in receptor',
                    'wcfg':{'text':"for all residues",
                            'variable':self.showMM},
                    'command': self.updateMM, 
                    'gridcfg':{'sticky':'e', 'row':-1, 'column':1}})
        ifd.append({'name': 'ssOtherCB',
                    'widgetType':Tkinter.Button,
                    'tooltip': "read in a different molecule for receptor secondary structure\n(useful for AD4 dockings 'missing' flexible residues in receptor)",
                    'wcfg':{'text':"open other file", 'relief':'flat'},
                    'command': self.openOtherMM, 
                    'gridcfg':{'sticky':Tkinter.EW, 'row':-1, 'column':2}})
        ifd.append({'name': 'resLabCB',
                    'widgetType':Tkinter.Checkbutton,
                    'tooltip': "display labels on all residues in 'receptor' which have\none or more atoms close to an atom in the 'ligand'",
                    'wcfg':{'text':"display labels on residues",
                            'variable':self.showResLabels},
                    'command': self.updateResLabels,
                    'gridcfg':{'sticky':Tkinter.W, 'columnspan':3}})
        ifd.append({'widgetType':Tkinter.Label,
                    'tooltip': 'print list of residue contacts in python shell',
                    'wcfg': {'text':"output list of:"},
                    'gridcfg':{'sticky':'w'} })
        ifd.append({'widgetType':Pmw.ComboBox,
            'name':'info_type',
            'wcfg':{'entryfield_value':self.infoType.get(),
                    'labelpos':'w',
                    'listheight':'80',
                    'scrolledlist_items': self.infoTypeList,
                    'selectioncommand': self.printInfo,
                    },
            'gridcfg':{'sticky':'w', 'row':-1, 'column':1, 'columnspan':2}})
        ifd.append({'name': 'piLabCB',
                    'widgetType':Tkinter.Checkbutton,
                    'tooltip': "display pi-pi interactions between 'ligand' and 'receptor'. Note: initial calculation is slow...",
                    'wcfg':{'text':"display pi-pi interactions",
                            'variable':self.showPiPi},
                    'command': self.updatePiPi,
                    'gridcfg':{'sticky':Tkinter.W}})
        ifd.append({'name': 'piCationLabCB',
                    'widgetType':Tkinter.Checkbutton,
                    'tooltip': "display cation-pi and pi-cation interactions between 'ligand' and 'receptor'. Note: initial calculation is slow...",
                    'wcfg':{'text':"display pi-cation interactions",
                            'variable':self.showPiCation},
                    'command': self.updatePiCation,
                    'gridcfg':{'sticky':Tkinter.W, 'row':-1, 'column':2}})
        ifd.append({'name': 'revertB',
                    'widgetType': Tkinter.Button,
                    'text':'Revert',
                    'tooltip': "Revert to previous display and destroy this form....",
                    'wcfg':{'bd':4,'state':'disabled'},
                    'gridcfg':{'sticky':'nesw'},
                    'command':self.revert})
        ifd.append({'name': 'closeB',
                    'widgetType': Tkinter.Button,
                    'text':'Close',
                    'tooltip': "Destroy this form.\nIt will be rebuilt if you change the conformation\nUse 'Edit->Undo ADanalyze_showBindingSiteInteractions' to restore...",
                    'wcfg':{'bd':4},
                    'gridcfg':{'sticky':'nesw', 'row':-1, 'column':1},
                    'command':self.close})
        ifd.append({'name': 'saveButton',
                    'widgetType':Tkinter.Button,
                    'wcfg':{'text':'Save Image'},
                    'command': self.save,
                    'gridcfg':{'sticky':Tkinter.W+Tkinter.E,'row':-1,'column':2}})
        self.displayform = self.vf.getUserInput(self.ifd, modal=0, blocking=0)
        self.displayform.root.protocol('WM_DELETE_WINDOW',self.close)
        e = self.ifd.entryByName['sphere_type']['widget']._entryfield
        e._entryFieldEntry.configure(width=8)
        e2 = self.ifd.entryByName['info_type']['widget']._entryfield
        e2._entryFieldEntry.configure(width=12)
        self.displayform.autoSize()
          

    def updateSphereDisplay(self, event=None):
        docking = self.docking
        lig = docking.ligMol
        macro = docking.macroMol
        oldVal = self.sphereType.get()
        if hasattr(self, 'ifd'):
            t = self.ifd.entryByName['sphere_type']['widget'].get()
            self.sphereType.set(t)
            if oldVal!=t:
                if t =='wireframe':
                    self.updateCpk()
                elif t=='solid':
                    self.restoreCpk()
        else:
            if self.sphereType.get()=='solid':
                self.restoreCpk()
            else:
                self.updateCpk()


    def openOtherMM(self, event=None):    
        mmFile = self.vf.askFileOpen(types=[('select filename for secondary structure calculation:', '*.pdb*'),
                            ('all files','*')],
                            title = 'Alternative Receptor File:')
        if mmFile:
            alt_receptor = self.vf.readMolecule(mmFile)
            if not alt_receptor:
                msg = "Error in reading " + mmFile
                self.warningMsg(msg)
                return 
            self.vf.displayLines(alt_receptor, negate=True, topCommand=0)
            self.vf.displaySticksAndBalls(alt_receptor, sticksBallsLicorice='Licorice',
                                negate=True, cradius=0.1, topCommand=0)
            self.docking.ss_macro = alt_receptor[0]
            

    def printInfo(self, event=None):
        t = self.ifd.entryByName['info_type']['widget'].get()
        if t =='ligand contacts':
            self.print_ligand_residue_contacts()
        elif t=='receptor contacts':
            self.print_macro_residue_contacts()
        if t=='hydrogen bonds':
            self.print_hydrogen_bonds()
        self.infoType.set(t)


    def print_macro_residue_contacts(self, event=None):
        print "\n\nresidues in 'receptor'-> 'ligand' residues in close contact"
        self.intDescr.print_macro_residue_contacts()
        print "\n"


    def print_ligand_residue_contacts(self, event=None):
        print "\n\nresidues in 'ligand'-> 'receptor' residues in close contact"
        self.intDescr.print_ligand_residue_contacts()
        print "\n"


    def print_hydrogen_bonds(self, event=None):
        print "\n\nhydrogen bonds (donor residue->acceptor residue(s))"
        self.intDescr.print_hb_residue()
        print "\n"


    def updateMsms(self):
        docking = self.docking
        lig = docking.ligMol
        macro = docking.macroMol
        rDict = self.intDescr.results
        lig_all_ats = rDict['lig_close_atoms']
        lig_hb_ats = rDict['lig_hb_atoms']
        lig_close_non_hb = rDict['lig_close_non_hb']
        surfName = lig.name + '_msms'
        docking.surfName = surfName
        if self.showMsms.get():
            #when display msms, hide ligand cpk
            self.vf.displayCPK(lig_all_ats, negate=1, topCommand=0)
            surfName = lig.name + '_msms'
            self.vf.displayMSMS(lig.allAtoms, surfName = surfName, topCommand=0)
            g = lig.geomContainer.geoms[surfName]
            from opengltk.OpenGL import GL
            g.inheritFrontPolyMode = 0
            g.Set(sharpColorBoundaries=False)
            g.Set(inheritSharpColorBoundaries=False)
            g.Set(opacity=[.75,.75])
            try:
                g.sortPoly()
            except:
                pass
        else:
            #when hide msms, show ligand cpk
            cpk_ats = AtomSet()
            if self.showCloseContactSpheres.get():
                cpk_ats += lig_close_non_hb[:]
            if self.showHbatSpheres.get():
                cpk_ats += lig_hb_ats[:]
            self.vf.displayCPK(cpk_ats, scaleFactor=1.0, quality=25, topCommand=0)
            self.vf.colorByAtomType(cpk_ats, geomsToColor=['cpk'], topCommand=0)
            self.vf.displayMSMS(lig.allAtoms, negate=1, topCommand=0)
        self.updateCpk()
        self.adjustMsms()


    def updateHBonds(self):
        docking = self.docking
        lig = docking.ligMol
        macro = docking.macroMol
        c = self.vf.hbondsAsSpheres
        rDict = self.intDescr.results
        macro_hb_atoms = rDict['macro_hb_atoms']
        lig_hb_atoms =  rDict['lig_hb_atoms']
        if not len(macro_hb_atoms + lig_hb_atoms):
            return
        h_ats = rDict['macro_hb_atoms'][:]
        if not self.showMsms.get():
            h_ats += rDict['lig_hb_atoms'][:]
        if self.showHbatSpheres.get(): 
            self.vf.displayCPK(h_ats, scaleFactor=1.0, quality=25, topCommand=0)
            self.vf.colorByAtomType(h_ats, geomsToColor=['cpk'], topCommand=0)
        else:
            self.vf.displayCPK(h_ats, negate=1, topCommand=0)
        self.updateCpk()


    def updateCCCpk(self, event=None):
        docking = self.docking
        lig = docking.ligMol
        macro = docking.macroMol
        rDict = self.intDescr.results
        cpk_atoms = rDict['macro_close_non_hb'][:]
        if not self.showMsms.get():
            cpk_atoms+=rDict['lig_close_non_hb'][:]
        if self.showCloseContactSpheres.get(): 
            self.vf.displayCPK(cpk_atoms, scaleFactor=1.0, quality=25, topCommand=0)
            self.vf.colorByAtomType(cpk_atoms, geomsToColor=['cpk'], topCommand=0)
        else:
            self.vf.displayCPK(cpk_atoms, negate=1, topCommand=0)
        self.updateCpk()


    def showHideSecondaryStructure(self):
        docking = self.docking
        macro = docking.ss_macro
        #reference macro molecule for AD4 dockings with flexible sidechains 
        lig = docking.ligMol
        rDict = self.intDescr.results
        ss_res = rDict['ss_res']
        if self.showSS.get():
            if len(ss_res)>3:
                if not self.showMM.get():
                    self.vf.displayExtrudedSS(ss_res.atoms, only=True, topCommand=0)
                else:
                    self.vf.displayExtrudedSS(macro.allAtoms, only=True, topCommand=0)
        else:
            #just hide the ribbon for the adjacent-close-residue sequences
            self.vf.displayExtrudedSS(ss_res.atoms, negate=True, topCommand=0)


    def updateSecondaryStructure(self):
        docking = self.docking
        macro = docking.macroMol
        #reference macro molecule for AD4 dockings with flexible sidechains 
        ss_macro = docking.ss_macro
        lig = docking.ligMol
        from opengltk.OpenGL import GL
        #attempt to show ribbon for contiguous residues in macromolecule
        rDict = self.intDescr.results
        ss_res = rDict['ss_res']
        if self.showSS.get():
            if len(ss_res)>3:
                if not self.showMM.get():
                    self.vf.ribbon(ss_res.atoms, only=True, topCommand=0)
                else:
                    if not hasattr(docking, 'has_ribbon'):
                        docking.has_ribbon = 1
                        self.vf.ribbon(ss_macro.allAtoms, topCommand=0)
                        #self.vf.ribbon(macro.allAtoms, topCommand=0)
                macro.geomContainer.geoms['secondarystructure'].Set(culling=GL.GL_NONE)
                if ss_macro!=macro:
                    if len(ss_macro.geomContainer.atoms['secondarystructure']):
                        ss_macro.geomContainer.geoms['secondarystructure'].Set(culling=GL.GL_NONE)
            else:
                #self.showSS.set(0)
                if not self.showMM.get():
                    res = self.getSSResidues(macro)
                    if len(res):
                        self.vf.displayExtrudedSS(res, negate=True, topCommand=0)
        elif self.showMM.get(): #and  not hasattr(docking, 'has_ribbon'):
            if not hasattr(docking, 'has_ribbon'):
                docking.has_ribbon = 1
                #self.vf.ribbon(macro.allAtoms, topCommand=0)
                self.vf.ribbon(ss_macro.allAtoms, topCommand=0)
                ss_macro.geomContainer.geoms['secondarystructure'].Set(culling=GL.GL_NONE)
        else:
            res = self.getSSResidues(macro)
            if len(res):
                self.vf.displayExtrudedSS(res, negate=True, topCommand=0)
                


    def updateMM(self, event=None):
        docking = self.docking
        macro = docking.macroMol
        ss_macro = docking.ss_macro
        if self.showMM.get():
            keys_to_exclude = ['sticks', 'AtomLabels', 'CAsticks', 'bonded', 'cpk', 'ChainLabels',\
                    'ProteinLabels', 'ResidueLabels', 'nobnds', 'balls', 'bondorder', 'lines', 'CAballs']
            displayed_atoms = AtomSet()
            gc = ss_macro.geomContainer
            for k in gc.atoms.keys():
                if not k in keys_to_exclude:
                    displayed_atoms += gc.atoms[k].findType(Atom)
            if len(displayed_atoms)!= len(ss_macro.allAtoms):
                if not hasattr(docking, 'has_ribbon'):
                    self.vf.ribbon(ss_macro.allAtoms, topCommand=0)
                    docking.has_ribbon=1
                else:
                    self.vf.displayExtrudedSS(ss_macro.allAtoms, topCommand=0)
            else:
                self.vf.displayExtrudedSS(ss_macro.allAtoms, topCommand=0)
                #self.vf.displayExtrudedSS(macro.allAtoms, topCommand=0)
        else:
            ss_res = self.intDescr.results['ss_res']
            if self.showSS.get() and len(ss_res)>3:
                self.vf.displayExtrudedSS(ss_res.atoms, only=True, topCommand=0)
                if ss_macro!=macro:
                    self.vf.displayExtrudedSS(ss_macro.allAtoms, negate=True, topCommand=0)
            else:
                self.vf.displayExtrudedSS(ss_macro.allAtoms, negate=True, topCommand=0)


    def reset(self):
        try:
            docking = self.docking
            macro = docking.macroMol
            ss_macro = docking.ss_macro
            lig = docking.ligMol
            self.vf.displayCPK(macro.allAtoms, negate=True, topCommand=0)
            residues_to_label = macro.chains.residues + self.flex_res
            self.vf.labelByProperty(residues_to_label, textcolor=(0.15,0.15,0.15), negate=True, topCommand=0)
            if not self.showMM.get() and hasattr(docking, 'has_ribbon'):
                self.vf.displayExtrudedSS(ss_macro.chains.residues, negate=True, topCommand=0)
            self.vf.displayCPK(lig.allAtoms, negate=True, topCommand=0)
            for item in ['msmsatoms','non_hbas','non_hbas_lig','lig_cs',\
            'hbas', 'hb_macro_ats', 'hb_lig_ats','reslabres', 'ssatoms','close_ats', \
            'close_lig_ats', 'close_macro_ats', 'close_only_macro', 'closeMacroRes']:
                try:
                    delattr(lig, item)
                except:
                    pass
                try:
                    delattr(macro, item)
                except:
                    pass
        except:
            pass


    def buildLigandDisplay(self):
        docking = self.docking
        lig = docking.ligMol
        #1. displayligand as sticksNballs, labelled by name
        #check whether the ligand is already displayed
        lgc = lig.geomContainer
        if len(lgc.atoms['sticks'])!=len(lig.allAtoms):
            self.vf.displaySticksAndBalls(lig.allAtoms, cquality=25, bquality=25, 
                        cradius=0.1, bRad=0.10, topCommand=0)
            if len(self.flex_res):
                self.vf.displaySticksAndBalls(self.flex_res.atoms, cquality=25, bquality=25, 
                                            cradius=0.1, bRad=0.1, topCommand=0)
            self.vf.colorByAtomType(lig, geomsToColor=['sticks', 'balls'],topCommand=0) #byAtom
        #1b. displayligand as msms, with low opacity
        surfName = lig.name+ '_msms'
        docking.surfName = surfName
        #self.surfName = surfName
        if self.flex_res_index:
            ind = self.flex_res_index
            self.vf.computeMSMS(lig.chains.residues[:ind], surfName, perMol=0, topCommand=0, display=self.showMsms.get())
        else:
            self.vf.computeMSMS(lig, surfName, perMol=1, topCommand=0, display=self.showMsms.get())
        self.adjustMsms()


    def adjustMsms(self):
        docking = self.docking
        lig = docking.ligMol
        lgc = lig.geomContainer
        rDict = self.intDescr.results
        surfName = docking.surfName
        #g = lig.geomContainer.geoms[surfName]
        #g.Set(sharpColorBoundaries=False)
        #g.Set(inheritSharpColorBoundaries=False)
        l_msms = lgc.geoms[surfName]
        l_msms.Set(inheritMaterial=0, sharpColorBoundaries=False,
                   inheritSharpColorBoundaries=False, opacity=[.75,.75])
        l_msms._setTransparent(1)
        lig.msmsatoms = lgc.atoms[surfName][:]
        self.vf.color(lig.allAtoms, [(0.5,0.5,0.5),], [surfName], topCommand=0)
        lig_close_atoms = rDict['lig_close_atoms']
        if len(lig_close_atoms):
            self.vf.colorByAtomType(lig_close_atoms, geomsToColor=[surfName, 'balls', 'sticks'],topCommand=0 ) # per atom type
        l_msms.sortPoly()
        #rDict = self.intDescr.results
        #macro_close_res = rDict['macro_close_res'] + rDict['macro_hb_res']
        #residues_to_label = macro_close_res + self.flex_res
        #self.vf.displaySticksAndBalls(residues_to_label.atoms, 
        #                        only=True,
        #                        sticksBallsLicorice='Licorice', cradius=0.1,topCommand=0)

    def buildCloseAtoms(self):
        docking = self.docking
        macro = docking.macroMol
        lig = docking.ligMol
        mgc = macro.geomContainer
        self.vf.displayLines(macro.allAtoms, negate=1,topCommand=0)
        rDict = self.intDescr.results
        macro_close_res = rDict['macro_close_res']
        macro_hbond_res = rDict['macro_hb_res']
        macro_res = macro_close_res + macro_hbond_res
        self.vf.displaySticksAndBalls(macro,
                        negate=True,
                        sticksBallsLicorice='Licorice', 
                        cradius=0.1)
        if len(macro_res):
            self.vf.displaySticksAndBalls(macro_res.atoms,
                        sticksBallsLicorice='Licorice', 
                        cradius=0.1)
            #self.vf.displayLines(macro_res.atoms, topCommand=0)
            residues_to_label = macro_res + self.flex_res
            self.vf.labelByProperty(residues_to_label, textcolor=(0.15, 0.15, 0.15),  properties=['name'],topCommand=0)
            reslabels = macro.geomContainer.geoms['ResidueLabels']
            reslabels.Set(fontScales =(0.3, 0.3, 0.3))
            self.vf.colorByAtomType(macro_close_res.atoms, geomsToColor=['balls','sticks','lines', 'cpk'], topCommand=0)
        surfName = lig.name+ '_msms'
        self.vf.color(lig.allAtoms, [(0.5,0.5,0.5),], [surfName],topCommand=0)
        self.vf.color(lig.allAtoms.get(lambda x:x.element=='C'), [(0.5,0.5,0.5),], ['balls', 'sticks'],topCommand=0)  #gray
        lig_close_atoms = rDict['lig_close_atoms']
        if len(lig_close_atoms):
            self.vf.colorByAtomType(lig_close_atoms, geomsToColor=[surfName, 'balls', 'sticks'],topCommand=0)
            
 
    def buildHbondDisplay(self):
        docking = self.docking
        lig = docking.ligMol
        macro = docking.macroMol
        rDict = self.intDescr.results
        c = self.vf.hbondsAsSpheres
        c.dismiss_cb()
        hbas = rDict['macro_hbas']
        if len(hbas):
            #self.vf.displayLines(hbas, topCommand=0)
            self.vf.displaySticksAndBalls(hbas, sticksBallsLicorice='Licorice',cradius=0.1, topCommand=0)
            #display hbonds as spheres 
            all = lig.allAtoms + macro.allAtoms
            c(all, topCommand=0)
            c.spacing_esw.set(.35)
            c.radii_esw.set(.06)
            self.vf.colorByAtomType(hbas, geomsToColor=['sticks', 'balls', 'lines'], topCommand=0)#atom
        surfName = lig.name+ '_msms'
        lig_cs = rDict['lig_close_carbons']
        if len(lig_cs): 
            self.vf.color(lig_cs, [(1,1,1),], [surfName,'balls', 'sticks'], topCommand=0)
        cpk_ats = AtomSet()
        if self.showCloseContactSpheres.get():
            cpk_ats += rDict['macro_close_non_hb'][:]
            if not self.showMsms.get():
                cpk_ats += rDict['lig_close_non_hb'][:]
        if self.showHbatSpheres.get(): 
            cpk_ats += rDict['macro_hb_atoms'][:]
            if not self.showMsms.get():
                cpk_ats += rDict['lig_hb_atoms'][:]
        if len(cpk_ats):
            self.vf.displayCPK(cpk_ats, scaleFactor=1.0, quality=25, topCommand=0)
            self.vf.colorByAtomType(cpk_ats, geomsToColor=['cpk'], topCommand=0)
            self.updateCpk()


    def build(self, percentCutoff=1.):
        #print "in build"
        self.showPiPi.set(0)
        self.showPiCation.set(0)
        self.pi_pi_geom.Set(visible=0)
        self.cation_pi_geom.Set(visible=0)
        self.pi_cation_geom.Set(visible=0)
        docking = self.docking
        lig = docking.ligMol
        macro = docking.macroMol
        self.flex_res = ResidueSet()
        self.flex_res_index = 0
        if hasattr(docking, 'flex_res') and len(docking.flex_res):
            self.flex_res_index = ind = len(lig.chains.residues) - len(docking.flex_res)
            self.flex_res = lig.chains.residues[ind:]
            #try to set the label fontScales here
            g = lig.geomContainer.geoms['ResidueLabels']
            g.Set(fontScales=(.3,.3,.3))
        from MolKit.interactionDescriptor import InteractionDescriptor
        #(1)
        if len(self.flex_res):
            flex_res_atoms = self.flex_res.atoms
            lig_atoms = lig.allAtoms - flex_res_atoms
            macro_atoms = macro.allAtoms + flex_res_atoms
            self.intDescr = InteractionDescriptor(lig_atoms, macro_atoms, percentCutoff)
        else:
            self.intDescr = InteractionDescriptor(lig, macro, percentCutoff)
        self.reset()
        self.vf.GUI.VIEWER.stopAutoRedraw()
        # build geoms for  the various displays...
        #rDict = self.intDescr.results
        #macro_close_res = rDict['macro_close_res'] + rDict['macro_hb_res']
        #residues_to_label = macro_close_res + self.flex_res
        #self.vf.displaySticksAndBalls(residues_to_label.atoms, 
        #                        only=True,
        #                        sticksBallsLicorice='Licorice', cradius=0.1,topCommand=0)
        self.buildLigandDisplay()
        self.buildCloseAtoms()
        self.buildHbondDisplay()
        self.updateCpk()
        self.updateSecondaryStructure()
        # add an input form here to toggle 
        # various displayed geometries on and off
        if not hasattr(self, 'displayform'):
            self.buildDisplayForm()
        docking.bindingSite = True
        self.vf.GUI.VIEWER.startAutoRedraw()
        self.vf.setbackgroundcolor(self.backgroundcolor)
        #c = self.vf.setbackgroundcolor
        #if len(c.cmdForms):
        #    f = c.cmdForms.values()[0]
        #    try:
        #        if not f.f.master.master.winfo_ismapped():
        #            self.vf.setbackgroundcolor.guiCallback()
        #    except:
        #        pass
        #else:
        #    self.vf.setbackgroundcolor.guiCallback()
        self.update_managed_sets()

    def updatePiCation(self, event=None):
        self.intDescr.detectPiInteractions()
        rDict = self.intDescr.results
        #add_geom = False
        #if not 'cation_pi' in rDict.keys():
        #    self.intDescr.detectPiInteractions()
        cation_pi = rDict['cation_pi']
        pi_cation = rDict['pi_cation']
        if len(cation_pi)==0 and len(pi_cation)==0:
            self.warningMsg('no cation_pi or pi_cation interactions detected')
            self.showPiCation.set(0)
            return
        if self.showPiCation.get():
            verts = []
            faces = []
            radii = []
            ctr = 0
            #these are ordered ligand_receptor 
            #and can be ligand_cation, receptor_pi OR ligand_pi, receptor_cation
            if len(cation_pi):
                for first,second in cation_pi: 
                    #either (atom, list of atoms) OR (list, atom)
                    if type(first)==types.ListType:
                        center1 = (Numeric.add.reduce(AtomSet(first).coords)/len(first)).tolist()
                        radii.append(0.7)
                    else:
                        center1 = (first.coords)
                        radii.append(0.1)
                    if type(second)==types.ListType:
                        center2 = (Numeric.add.reduce(AtomSet(second).coords)/len(second)).tolist()
                        radii.append(0.7)
                    else:
                        center2 = (second.coords)
                        radii.append(0.1)
                    faces.append((ctr,ctr+1))
                    verts.append(center1)
                    verts.append(center2)
                    ctr += 2
                self.cation_pi_geom.Set(radii=radii, faces=faces, 
                                        vertices=verts, visible=1)
            verts = []
            faces = []
            radii = []
            ctr = 0
            if len(pi_cation):
                for first,second in pi_cation: #either (atom, list of atoms)
                    if type(first)==types.ListType:
                        center1 = (Numeric.add.reduce(AtomSet(first).coords)/len(first)).tolist()
                        radii.append(0.7)
                    else:
                        center1 = (first.coords)
                        radii.append(0.1)
                    if type(second)==types.ListType:
                        center2 = (Numeric.add.reduce(AtomSet(second).coords)/len(second)).tolist()
                        radii.append(0.7)
                    else:
                        center2 = (second.coords)
                        radii.append(0.1)
                    faces.append((ctr,ctr+1))
                    verts.append(center1)
                    verts.append(center2)
                    ctr += 2
                self.pi_cation_geom.Set(radii=radii, faces=faces, 
                                            vertices=verts, visible=1)
        else:
            self.cation_pi_geom.Set(visible=0)
            self.pi_cation_geom.Set(visible=0)
    

    def updatePiPi(self, event=None):
        self.intDescr.detectPiInteractions()
        rDict = self.intDescr.results
        #add_geom = False
        #if not 'pi_pi' in rDict.keys():
        #    self.intDescr.detectPiInteractions()
        #    add_geom = True
        pi_pi = rDict['pi_pi']
        if not len(pi_pi):
            self.warningMsg('no pi_pi interactions detected')
        t_shaped = rDict['t_shaped']
        if not len(t_shaped):
            self.warningMsg('no t_shaped interactions detected')
        if len(pi_pi)==0 and len(t_shaped)==0:
            self.showPiPi.set(0)
            return
        if self.showPiPi.get():
            verts = []
            faces = []
            ctr = 0
            for first, second in pi_pi: #either (atom, list of atoms)
                #print "processing pi_pi"
                center1 = (Numeric.add.reduce(first.coords)/len(first)).tolist()
                center2 = (Numeric.add.reduce(second.coords)/len(second)).tolist()
                faces.append((ctr,ctr+1))
                verts.append(center1)
                verts.append(center2)
                ctr += 2
            for first, second in t_shaped: #either (atom, list of atoms)
                #print "processing t_shaped"
                center1 = (Numeric.add.reduce(first.coords)/len(first)).tolist()
                center2 = (Numeric.add.reduce(second.coords)/len(second)).tolist()
                faces.append((ctr,ctr+1))
                verts.append(center1)
                verts.append(center2)
                ctr += 2
            self.pi_pi_geom.Set(faces=faces, vertices=verts, visible=1)
        else:
            self.pi_pi_geom.Set(visible=0)

    
    def updateResLabels(self, event=None):
        rDict = self.intDescr.results
        macro_close_res = rDict['macro_close_res'] + rDict['macro_hb_res']
        residues_to_label = macro_close_res + self.flex_res
        show = self.showResLabels.get()
        negate = not show
        #self.vf.labelByProperty(macro_close_res, negate=negate, 
        self.vf.displaySticksAndBalls(residues_to_label[0].top, negate=True, 
                                sticksBallsLicorice='Licorice'), 
        self.vf.labelByProperty(residues_to_label, negate=negate, 
                                textcolor=(0.15, 0.15, 0.15),  
                                properties=['name'],topCommand=0)
        self.vf.displaySticksAndBalls(residues_to_label.atoms, 
                                #only=True,
                                sticksBallsLicorice='Licorice', cradius=0.1,topCommand=0)


    def doit(self, docking_filename, **kw):
        try:
            docking = self.vf.dockings[docking_filename]
            if not hasattr(docking, 'macroMol'):
                msg = 'Please choose the macromolecule for this docking first!'
                self.warningMsg(msg)
                return "ERROR"
        except:
            msg = "no docking for dlg filename", docking_filename, " currently in viewer!"
            self.warningMsg(msg)
            return "ERROR"
        self.docking = docking
        if not hasattr(self.docking, 'ss_macro'):
            self.docking.ss_macro = self.docking.macroMol
        self.build()
        docking.bindingSite = True
        if hasattr(docking.ligMol, 'spw'):
            docking.ligMol.spw.updateBindingSite = self.updateDisplay.get()


ADShowBindingSiteGUI=CommandGUI()
ADShowBindingSiteGUI.addMenuCommand('AutoToolsBar', menuText['AnalyzeMB'], 
        menuText['showBindingSite'], cascadeName = menuText['DockingLogMB'], separatorAbove=1)



class ADMakeAllGrids(MVCommand):
    """ toplevel GUI with one row per grid_type allowing:loading a grid (first, then->),toggling the visiblity of its isosurface,changing sampling rate,changing isocontour isovalue,changing the render mode between LINE and FILL visibility of the grid's bounding box
    \nPackage : AutoDockTools
    \nModule  : autoanalyzeCommands
    \nClass   : ADMakeAllGrids
    \nCommand : ADanalyze_showGridIsocontours 
    \nSynopsis:\n     
        None<---ADanalyze_showGridIsocontours(gridFile)
    \nRequired Arguments:\n
        gridFile --- AutoGrid map filename
    \nOptional Argument:\n
        ask --- whether to open a filebrowser if specified filename does not exist
    """


    #def checkDependencies(self):
    #    import isocontour
    #    import Pmv.Grid


    def onAddCmdToViewer(self):
        checkHasInitializedDockings(self.vf)
        if not hasattr(self.vf, 'grids'):
            self.vf.grids = {}
        #self.vf.loadModule('displayCommands','Pmv')
        if self.vf.hasGui:
            for item in ['colorByMolecules','colorByAtomType', 'color']:
                if not hasattr(self.vf, item):
                    self.vf.loadCommand('colorCommands', item, 'Pmv')
            for item in ['displayLines','showMolecules']:
                if not hasattr(self.vf, item):
                    self.vf.loadCommand('displayCommands', item, 'Pmv')


    def showIFD(self):
        d = self.vf.docked
        if not d:
            self.vf.warningMsg('Please Read a Docking Log First')
            return
        if hasattr(d,'ligMol') and hasattr(d.ligMol, 'vina_energy'):
            self.vf.warningMsg('Current docking is a vina result. Visualizing grids is not possible.')
            return "ERROR"
        if not hasattr(self, 'ifd'):
            self.hasGrid = 0
            if len(d.dlo_list):
                dlo = d.dlo_list[-1]
                if hasattr(dlo, 'macroStem'):
                    macroLab = '%s AUTOGRIDS:'%dlo.macroStem
            else:
                macroLab = '%s AUTOGRIDS:'
            dlgstr = ': ' + d.dlo_list[-1].filename
            ifd = self.ifd= InputFormDescr(title='Visualize AutoGrids'+dlgstr)
            ifd.append({'name': 'macroLab',
                        'widgetType':Tkinter.Label,
                        'text':macroLab,
                        'gridcfg':{'sticky':Tkinter.W, 'columnspan':3}})
            ifd.append({'name': 'centerLab',
                        'widgetType':Tkinter.Label,
                        'text':'center: ',
                        'gridcfg':{'sticky':Tkinter.W, 'row':-1,'column':3,
                            'columnspan':2}})
            ifd.append({'name': 'nptsLab',
                        'widgetType':Tkinter.Label,
                        'text':'number of points: ',
                        'gridcfg':{'sticky':Tkinter.W, 'columnspan':3}})
            ifd.append({'name': 'spaceLab',
                        'widgetType':Tkinter.Label,
                        'text':'spacing: ',
                        #'gridcfg':{'sticky':Tkinter.W, 'columnspan':5}})
                        'gridcfg':{'sticky':Tkinter.W, 'row':-1,'column':3,
                            'columnspan':2}})
            ifd.append({'widgetType': Tkinter.Label,
                        'text':'___________________________________________________',
                        'wcfg':{'bd':6},
                        'gridcfg':{'sticky':Tkinter.W, 'columnspan':4}})
            ifd.append({'widgetType': Tkinter.Label,
                        'text':'Display\nMap',
                        'wcfg':{'bd':6},
                        'gridcfg':{'sticky':Tkinter.W}})
            ifd.append({'widgetType': Tkinter.Label,
                        'text':'Sampling',
                        'gridcfg':{'sticky':Tkinter.W, 'row':-1,'column':1}})
            ifd.append({'widgetType': Tkinter.Label,
                        'text':'IsoValue',
                        'wcfg':{'bd':6},
                        'gridcfg':{'sticky':Tkinter.W+Tkinter.E,
                             'row':-1,'column':2,'columnspan':2}})
            ifd.append({'widgetType': Tkinter.Label,
                        'text':'RenderMode',
                        'wcfg':{'bd':6},
                        'gridcfg':{'sticky':Tkinter.W, 'row':-1,'column':4}})
            ifd.append({'name': 'closeB',
                'widgetType': Tkinter.Button,
                'text':'Close',
                'wcfg':{'bd':4},
                'gridcfg':{'sticky':Tkinter.E+Tkinter.W, 'columnspan':6},
                'command':self.Close_cb})
            self.form = self.vf.getUserInput(self.ifd, modal=0,blocking=0)
            self.form.root.protocol('WM_DELETE_WINDOW',self.Close_cb)
        else:
            if hasattr(self,'form'):
                self.form.deiconify()


    def guiCallback(self):
        d = self.vf.docked
        if not d:
            self.vf.warningMsg('Please Read a Docking Log First')
            return
        if hasattr(d.ligMol, 'vina_energy'):
            self.vf.warningMsg('Current docking is a vina result. Visualizing grids is not possible.')
            return "ERROR"
        #build or open ifd
        self.showIFD()
        gridObjs = []
        if not hasattr(d, 'gridObjects'):
            dlo = d.dlo_list[0]
            p = dlo.parser
            if p.version>=4.0:
                at_types = dlo.dpo['ligand_types']['value'].split()
            else:
                at_types =  dlo.dpo['types']['value']
            for item in at_types:
                gridFile = dlo.macroStem + '.' + item + '.map'
                gridObjs.append(gridFile)
            d.gridObjs = gridObjs
        ifd2= InputFormDescr(title='Load Grid')
        ifd2.append({'widgetType':'ListChooser',
            'name':'gridObjs',
            'entries':d.gridObjs,
            'mode':'multiple',
            'wcfg':{'title':'Pick grids '},
            'lbwcfg':{'height':4},
            'gridcfg':{'sticky':'w', 'column':100,
            'rowspan':10}})
        val = self.vf.getUserInput(ifd2)
        if len(val)>0 and len(val['gridObjs'])>0:
            filelist = val['gridObjs']
            for item in filelist:
                if os.path.exists(item):
                    self.doitWrapper(item, log=1, redraw=0)
                else:
                    fileList = glob.glob(item)
                    gridFile = self.vf.askFileOpen(types=[('select mapfile','*.map'),
                            ('all files','*')],
                            title=item)
                    #if gridFile is not None:
                    if gridFile:
                        self.doitWrapper(gridFile, log=1, redraw=0)
                    else: return 'ERROR'


    def Close_cb(self, event=None):
        self.form.root.withdraw()


    def __call__(self, gridFile, ask=1, **kw):
        """None<---ADanalyze_showGridIsocontours(gridFile)
        \ngridFile --- AutoGrid map filename
        \nask --- whether to open a filebrowser if specified filename does not exist
        """
        if gridFile is None:
            return 'ERROR'
        elif not len(strip(gridFile)):
            #this is (probably) only from a script
            raise IOError
        else:
            kw['ask'] = ask
            return apply(self.doitWrapper, (gridFile,), kw)


    def doit(self, gridFile, **kw):
        """gridFile macroMol.name.atType.map"""
        d = self.vf.docked
        if hasattr(d,'ligMol') and hasattr(d.ligMol, 'vina_energy'):
            self.vf.warningMsg('Current docking is a vina result. Visualizing grids is not possible.')
            return "ERROR"
        if self.vf.grids.has_key(gridFile): 
            newGrid = self.vf.grids[gridFile]
        else:
            if os.path.exists(gridFile):
                newGrid = AutoGrid(gridFile)
            else:
                ask = kw.get('ask',1)
                if ask:
                    #if the file doesn't exist ask: 
                    gridFile = self.vf.askFileOpen(types=[('select mapfile','*.map'),
                        ('all files','*')],
                        title=gridFile)
                    if gridFile:
                        newGrid = AutoGrid(gridFile)
                    else: 
                        return 'ERROR'
                else:
                    #script tries to read a non-existent file
                    #so we want to raise a particular error here
                    raise IOError

                #self.lastCmdLog = [self.logString(gridFile)]
            atType = split(gridFile,'.')[-2]
            nname = os.path.basename(gridFile)
            if hasattr(d,'gridObjs') and nname in d.gridObjs: 
                d.gridObjs.remove(nname)
            if os.path.dirname(gridFile)==os.getcwd():
                newGrid.name = nname
            else:
                newGrid.name = gridFile
            newGrid.atType = atType
            self.vf.grids[gridFile] = newGrid
            newGrid.srf = None
        #now make the ifd_entry
        if newGrid.srf:
            t = 'this grid already has a srf'
            self.vf.warningMsg(t)
            return 'ERROR'
        #open the form for isocontours
        self.showIFD()
        newGrid.surfGUI = AutoGridSurfaceGui(newGrid,self.vf)
        newGrid.surfGUI.makeIfdEntry(self.ifd,['closeB'])
        if self.hasGrid==0:
            #update the CENTER, npts and spacing labels
            cenStr = ""
            for c in newGrid.CENTER:
                cenStr = cenStr + str(round(c,4)) + ' '
            #self.ifd.entryByName['centerLab']['widget'].config(text='center: '+str(newGrid.CENTER))
            entD = self.ifd.entryByName
            entD['centerLab']['widget'].config(text='center: '+ cenStr)
            entD['nptsLab']['widget'].config(text='number of points: '+str(newGrid.NELEMENTS))
            entD['spaceLab']['widget'].config(text='spacing: '+str(newGrid.SPACING))
            self.hasGrid = 1
        if not hasattr(d,'grids'):
            d.grids = {}
        d.grids[newGrid.name] = newGrid
        return newGrid


ADMakeAllGridsGUI=CommandGUI()
ADMakeAllGridsGUI.addMenuCommand('AutoToolsBar', menuText['AnalyzeMB'], 
            menuText['showGridsMB'], cascadeName=menuText['GridsMB'])
                


class ADGetOutput(MVCommand):
    """When the docking log is parsed, summary lines are detected and saved in  
    \nPackage : AutoDockTools
    \nModule  : autoanalyzeCommands
    \nClass   : ADGetOutput
    \nCommand : ADanalyze_showResultsOutput
    \nSynopsis:\n    
        None<---ADanalyze_showResultsOutput(dockStr)
    \nRequired Arguments:\n    
        dockStr --- list of dockingfiles, one for each docking present in the
        viewer
    """

    def guiCallback(self):
        """called each time the 'get output' button is pressed"""
        #THIS SHOULD BE A LIST BOX OF CURRENT DOCKING DICTIONARIES
        entries = self.vf.dockings.keys()
        if len(entries)==0:
            self.vf.warningMsg('Please Read a Docking Log First')
            return
        elif len(entries)>1:
            ifd = self.ifd=InputFormDescr(title = 'Select Dockings')
            ifd.append({'name': 'dockings',
                'title':' ',
                'widgetType':'ListChooser',
                'entries': entries,
                'mode': 'multiple',
                'gridcfg':{'row':0,'column':0},
                'lbwcfg':{'height':10, 'width':60}})
            vals = self.vf.getUserInput(ifd)
            if len(vals)>0 and len(vals['dockings'])>0:
                dockstr=vals['dockings']
                self.doitWrapper(dockstr,log=1,redraw=0)
        else:
            self.doitWrapper(entries,log=1,redraw=0)


    def doit(self, dockStr):
        """list of selected keys into self.vf.dockings"""
        ##self.lastCmdLog = self.logString(dockStr)

        #ss="""
	        #CLUSTERING HISTOGRAM
	        #____________________
#
#________________________________________________________________________________
         #|           |     |           |     |                                    
 #Clus | Lowest   | Run | Mean      | Num | Histogram                          
  #-ter | Docked  |     | Docked    | in  |                                    
 #Rank | Energy    |     | Energy    | Clus|    5    10   15   20   25   30   35
#_____|________|_____|______|_____|____:____|____:____|____:____|____:___
   #"""
        ss =        "            CLUSTERING HISTOGRAM     \n"
        ss = ss +   "            ____________________     \n\n"
        ss = ss +   "_____|___________|_____|___________|_____|____:____|____:____|____:____|____:___\n"
        ss = ss +   "Clus | Lowest    | Run | Mean      | Num | Histogram                          \n"
        ss = ss +   "-ter | Docked    |     | Docked    | in  |                                    \n" 
        ss = ss +   "Rank | Energy    |     | Energy    | Clus|    5    10   15   20   25   30   35\n"
        ss = ss +   "_____|___________|_____|___________|_____|____:____|____:____|____:____|____:___\n"



        for item in dockStr:
            d = self.vf.dockings[item]
            #FIX THIS: does Docking know its output????
            for dlo in d.dlo_list:
                for l in dlo.output:
                    ss = ss + str(l)    
                    ss = ss + '\n\n'
        ifd = InputFormDescr(title='Docking Summary')
        ifd.append({'name':'summaryText',
                'size':[100,40],
                'label':'',
                'defaultValue':ss,
                'widgetType':'ScrolledText',
                'readButton':1,
                'writeButton':1,
                'readFileType':[('DockingSummary Files','*.summary')],
                'writeFileType':[('DockingSummary Files','*.summary')]})
        vals= self.vf.getUserInput(ifd)
        if len(vals)==0:
            return 'ERROR'


ADGetOutputGUI=CommandGUI()
ADGetOutputGUI.addMenuCommand('AutoToolsBar', menuText['AnalyzeMB'], 
        menuText['getOutputMB'] , cascadeName=menuText['StatesMB'])



class ADGetAGrid(MVCommand):
    """Adds a row to the toplevel GUI of ADMakeAllGrids enabling:loading a USER grid (first, then->),toggling its visiblity,changing sampling rate,changing isocontour value,changing the render mode between lines and fill.
    \nPackage : AutoDockTools
    \nModule  : autoanalyzeCommands
    \nClass   : ADGetAGrid
    \nCommand : ADanalyze_addExtraGridIsocontour
    \nSynopsis:\n    
        None<---ADanalyze_addExtraGridIsocontour(gridFile)
    \nRequired Arguments:\n    
        gridFile --- extra autogrid output file for which an isocontour will be displayed
        viewer
    """
    def guiCallback(self):
        """called each time the 'select grid' button is pressed"""
        d = self.vf.docked
        if not d:
            self.vf.warningMsg('Please Read a Docking Log First')
            return
        if hasattr(d.ligMol, 'vina_energy'):
            self.warningMsg("Current docking is a vina result. Visualizing grids is not possible")
            return "ERROR"
        gridFile = self.vf.askFileOpen(types=[('select grid filename:', '*.map'),
                            ('all files','*')],
                            title = 'Grid File:')
        if gridFile:
            self.doitWrapper(gridFile,log=1,redraw=0)


    def doit(self, gridFile):
        ##self.lastCmdLog = self.logString(gridFile)
        if not hasattr(self.vf.ADanalyze_showGridIsocontours, 'ifd'):
            self.vf.ADanalyze_showGridIsocontours.guiCallback()
        self.vf.ADanalyze_showGridIsocontours(gridFile,topCommand=0)

ADGetAGridGUI=CommandGUI()
ADGetAGridGUI.addMenuCommand('AutoToolsBar', menuText['AnalyzeMB'], 
        menuText['addGridMB'], cascadeName=menuText['GridsMB'])



class ADSelectDLG(MVCommand):
    """Allows the user to change current docking 
    \nPackage : AutoDockTools
    \nModule  : autoanalyzeCommands
    \nClass   : ADSelectDLG
    \nCommand : ADanalyze_selectDLG
    \nSynopsis:\n    
        None<---ADanalyze_selectDLG
    \nRequired Arguments:\n    
        dockstr --- selected key into self.vf.dockings
    """
    def onAddCmdToViewer(self):
        checkHasInitializedDockings(self.vf)
        #self.vf.loadModule('displayCommands','Pmv')
        if self.vf.hasGui:
            for item in ['colorByMolecules','colorByAtomType', 'color']:
                if not hasattr(self.vf, item):
                    self.vf.loadCommand('colorCommands', item, 'Pmv')
            for item in ['displayLines','showMolecules']:
                if not hasattr(self.vf, item):
                    self.vf.loadCommand('displayCommands', item, 'Pmv')
        if 'Player GUI' not in self.vf.userpref.keys():
            doc = """Use StatesPlayerWidget or ConformationPlayer to
            visualized docked structures. Default value is 'StatesPlayerWidget'"""
            self.vf.userpref.add('Player GUI', 'ConformationPlayer', 
                    ['StatesPlayerWidget', 'ConformationPlayer'],
                    callbackFunc = [self.set_playerGUI], doc=doc ,
                    category='AutoDockTools')
        try:
            from AutoDockTools import WebServices
            self.vf.browseCommands('WebServices', package='AutoDockTools')
        except:
            pass
            

    def set_playerGUI(self, name, oldval, newval):
        self.vf.userpref['Player GUI']['value'] = newval

    def guiCallback(self):
        """called each time the 'select docking logfile ' button is pressed"""
        entries = self.vf.dockings.keys()
        if not len(entries):
            self.vf.warningMsg('No docking logs have been read')
            return
        elif len(entries)==1:
            self.vf.warningMsg('Current docking is only docking in viewer')
            return

        ifd = self.ifd=InputFormDescr(title = 'Select a Docking')
        ifd.append({'name': 'dockings',
            'title':'Pick Docking ',
            'widgetType':'ListChooser',
            'entries': entries,
            'mode': 'single',
            'gridcfg':{'row':0,'column':0},
            'lbwcfg':{'height':5, 'width':40}})
        vals = self.vf.getUserInput(ifd)
        if len(vals)>0 and len(vals['dockings'])>0:
            dockstr=vals['dockings'][0]
            self.doitWrapper(dockstr,log=1,redraw=0)


    def doit(self, dockstr):
        """selected key into self.vf.dockings"""
        #msg = 'Current Docking: '+dockstr
        #self.vf.GUI.pickLabel.configure(text=msg)
        ##self.lastCmdLog = self.logString(dockstr)
        self.vf.docked = self.vf.dockings[dockstr]
        #if hasattr(self.vf.ADanalyze_showHistogram, 'canvas') and self.vf.ADanalyze_showHistogram.canvas.winfo_ismapped():
        #    l= self.vf.ADanalyze_showHistogram.buildIt()

ADSelectDLGGUI=CommandGUI()
ADSelectDLGGUI.addMenuCommand('AutoToolsBar', menuText['AnalyzeMB'], 
        menuText['selectDLG'], cascadeName = menuText['DockingLogMB'])



class ADDeleteDLG(MVCommand):
    """Allows the user to delete docking 
    \nPackage : AutoDockTools
    \nModule  : autoanalyzeCommands
    \nClass   : ADDeleteDLG
    \nCommand : ADanalyze_deleteDLG
    \nSynopsis:\n    
        None<---ADanalyze_deleteDLG
    \nRequired Arguments:\n    
        dockingKey --- selected key into self.vf.dockings
    """
    def onAddCmdToViewer(self):
        checkHasInitializedDockings(self.vf)
        if self.vf.hasGui:
            for item in ['colorByMolecules','colorByAtomType', 'color']:
                if not hasattr(self.vf, item):
                    self.vf.loadCommand('colorCommands', item, 'Pmv')
            for item in ['displayLines','showMolecules']:
                if not hasattr(self.vf, item):
                    self.vf.loadCommand('displayCommands', item, 'Pmv')


    def guiCallback(self):
        """called each time the 'delete docking logfile ' button is pressed"""
        entries = self.vf.dockings.keys()
        if not len(entries):
            self.vf.warningMsg('No docking logs have been read')
            return
        ifd = self.ifd=InputFormDescr(title = 'Delete a Docking')
        ifd.append({'name': 'dockings',
            'title':'Pick Docking To Delete',
            'widgetType':'ListChooser',
            'entries': entries,
            #THINK ABOUT THIS!
            'mode': 'single',
            'gridcfg':{'row':0,'column':0},
            'lbwcfg':{'height':5, 'width':40}})
        vals = self.vf.getUserInput(ifd)
        if len(vals)>0 and len(vals['dockings'])>0:
            dockingKey = vals['dockings'][0]
            self.doitWrapper(dockingKey,log=1,redraw=0)


    def doit(self, dockingKey):
        """selected key into self.vf.dockings"""
        if not dockingKey in self.vf.dockings.keys():
            raise KeyError, dockingKey + " is not a docking in viewer"
        docking = self.vf.dockings[dockingKey]
        if hasattr(self.vf.ADanalyze_showDockingsAsSpheres, 'keys'):
            for k in self.vf.ADanalyze_showDockingsAsSpheres.keys:
                if dockingKey in k:
                    self.vf.ADanalyze_showDockingsAsSpheres.reset()
        ##DOES IT MATTER WHETHER selected dlg is current:
        #if self.vf.docked==docking:
            #msg='Current Docking: None'
            #self.vf.GUI.pickLabel.configure(text=msg)
        dlgMols = []
        self.vf.deleteMol(docking.ligMol, log=0)
        ####THIS DELETES THE DOCKING RIGHT HERE!!!
        if hasattr(docking,'chooser'):
            docking.chooser.form.root.destroy()    
        if hasattr(self.vf.ADanalyze_showGridIsocontours, 'ifd'):
            docking.grids = {}
            self.vf.GUI.VIEWER.Redraw()
        #if docking==self.vf.docked and hasattr(self.vf.ADanalyze_showHistogram,'canvas'):
        #    l = self.vf.ADanalyze_showHistogram.canvas.find_all()
        #    for item in l:
        #        self.vf.ADanalyze_showHistogram.canvas.delete(item)
        #        self.vf.ADanalyze_showHistogram.close()
        #        self.vf.ADanalyze_showHistogram.canvas.config({'height':450,'width':660})
        if self.vf.dockings.has_key(dockingKey):
            del self.vf.dockings[dockingKey]
        if docking==self.vf.docked: 
            if len(self.vf.dockings.keys()):
                key = self.vf.dockings.keys()[0]
                self.vf.docked = self.vf.dockings[key]
            else:
                self.vf.docked = None
        elif self.vf.docked is None and len(self.vf.dockings.keys()):
            self.vf.docked = self.vf.dockings.values()[0]
            


ADDeleteDLGGUI=CommandGUI()
ADDeleteDLGGUI.addMenuCommand('AutoToolsBar', menuText['AnalyzeMB'], 
        menuText['deleteDLG'], cascadeName = menuText['DockingLogMB'])



class ADGetDirDLGs(MVCommand):
    """Allows the user to read all the docking logs in a chosen directory
    \nPackage : AutoDockTools
    \nModule  : autoanalyzeCommands
    \nClass   : ADGetDirDLGs
    \nCommand : ADanalyze_readAllDLGinDirectory
    \nSynopsis:\n
        None<---ADanalyze_readAllDLGinDirectory(dlgDir, makeOneDocking)
    \nRequired Arguments:\n
        dlgDir --- directory containing any number of docking log files
    \nOptional Arguments:\n    
        \nmakeOneDocking --- whether to compile all dlg results into one Docking object with 1 list of conformations and ability to do clustering if not, a separate Docking is created for each docking logfile.
    """


    def guiCallback(self):
        """called each time the 'read all docking logfiles' button is pressed"""
        #need to have a listchooser here with entries of dlg currently read-in
        dlgFile = self.vf.askFileOpen(types=[('select any docking log filename:', '*.dlg'),
                            ('all files','*')],
                            title = 'To Set Directory:Select any dlg file')
        if dlgFile:
            #ask whether to make it one Docking or many
            t = 'Should all the docking logfiles in this directory be read into one Docking?'
            d = SimpleDialog(self.vf.GUI.ROOT, text = t,
            buttons = ["No", "Yes", "Cancel"], default = 1, 
                    title = "Create a Single Docking? ")
            makeOneDocking = d.go()
            if makeOneDocking==2:
                return 'ERROR'
            dlgDir=os.path.dirname(dlgFile)
            self.doitWrapper(dlgDir, makeOneDocking, log=1, redraw=0)


    def __call__(self, dlgDir, makeOneDocking=0, **kw):
        """None<---ADanalyze_readAllDLGinDirectory(dlgDir, makeOneDocking)
        \ndlgDir --- directory containing any number of docking log files
        makeOneDocking --- whether to compile all dlg results into one Docking object with 1 list of conformations and ability to do clustering if not, a separate Docking is created for each docking logfile.
        """
        apply(self.doitWrapper, (dlgDir, makeOneDocking,), kw)


    def doit(self, dlgDir, makeOneDocking, **kw):
        oldPref = self.vf.userpref['Warning Message Format']['value']
        self.vf.setUserPreference(('Warning Message Format','printed'))
        dlgList = glob.glob(dlgDir + '/*.dlg')
        if len(dlgList):
            dlgList.sort()
            kw = {'log':0}
            for dlgName in dlgList:
                apply(self.vf.ADanalyze_readDLG,(dlgName, makeOneDocking,), kw)
        self.vf.setUserPreference(('Warning Message Format', oldPref))
                    

ADGetDirDLGsGUI=CommandGUI()
ADGetDirDLGsGUI.addMenuCommand('AutoToolsBar', menuText['AnalyzeMB'],
        menuText['readDirDLG'], cascadeName = menuText['DockingLogMB'])



class ADGetDLG(MVCommand):
    """Allows the user to select a filename for the docking log
    \nPackage : AutoDockTools
    \nModule  : autoanalyzeCommands
    \nClass   : ADGetDLG
    \nCommand : ADanalyze_readDLG
    \nSynopsis:\n
        None<-ADanalyze_readDLG(dlgFile, addToPrevious=0)
    \nRequired Arguments:\n    
        dlgFile --- docking log file
    \nOptional Arguments:\n    
        addToPrevious --- if 0, start a new Docking otherwise add docked conformations to a previous docking
    """


    def onRemoveObjectFromViewer(self, obj):
        if hasattr(obj, 'docking') and hasattr(obj, 'torTree'):
            #clean-up self.vf.dockings:
            remove_keys = []
            for k, v in self.vf.dockings.items():
                if v==obj.docking:
                    remove_keys.append(k)
            for k in remove_keys:
                d = self.vf.dockings[k] 
                if hasattr(d, 'ch') and hasattr(d.ch, 'mol') and d.ch.mol==obj:
                    for c in d.ch.conformations:
                        delattr(c, 'mol')
                    delattr(d.ch, 'mol')
                del self.vf.dockings[k] 
            #clean-up self.vf.docked:
            if obj.docking==self.vf.docked:
                d = self.vf.docked
                #if hasattr(d, 'ch') and hasattr(d.ch, 'mol') and d.ch.mol==obj:
                if hasattr(d, 'ch'): 
                    if hasattr(d.ch, 'conformations'):
                        for conf in d.ch.conformations:
                            if hasattr(conf, 'subset'):
                                delattr(conf, 'subset')
                    if hasattr(d.ch, 'mol'):
                        delattr(d.ch, 'mol')
                if hasattr(d, 'clusterer') and hasattr(d.clusterer, 'subset'):
                    delattr(d.clusterer,'subset')
                self.vf.docked = None

            #remove reference to obj.ligMol
            if hasattr(obj.docking, 'ligMol'):
                delattr(obj.docking, 'ligMol')

            #clean-up obj.docking
            del obj.docking

            #try to clean-out torTree reference
            for key in obj.torTree.rootNode.__dict__.keys():
                obj.torTree.rootNode.__dict__[key] = None
            for item in obj.torTree.torsionMap:
                for key in item.__dict__.keys():
                    item.__dict__[key] = None
            del obj.torTree


    def checkMolLines(self, mol):
        #if len(mol.geomContainer.geoms['lines'].faceSet.faces.array)==0:
        if len(mol.geomContainer.geoms['bonded'].faceSet.faces.array)==0:
            cbvar = self.vf.showMolecules.molList[mol.name]
            n = 1-cbvar.get()
            #CHECK FOR BONDS HERE:
            if len(mol.allAtoms.bonds[0])==0:
                mol.buildBondsByDistance()
            self.vf.displayLines(mol, negate=n, topCommand=0)
            if not hasattr(mol, 'colored') or mol.colored: return
            #FIX THIS
            #THIS IS STUPID
            ligName = self.vf.docked.ligMol.name
            if mol.name==ligName and not mol.colored:
                self.vf.colorByMolecules(mol,['lines'], topCommand=0, redraw=1)
            else:
                self.vf.colorByAtomType(mol,['lines'], topCommand=0)
                aromaticCs = AtomSet(filter(lambda x: x.autodock_element=='A',mol.allAtoms))
                if len(aromaticCs):
                    self.vf.color(aromaticCs,((0.,1.,0.),),['lines'],topCommand=0, redraw=1)
                mol.colored=1


    def onAddCmdToViewer(self):
        checkHasInitializedDockings(self.vf)
        if self.vf.hasGui:
            for item in ['colorByMolecules','colorByAtomType', 'color']:
                if not hasattr(self.vf, item):
                    self.vf.loadCommand('colorCommands', item, 'Pmv')
            for item in ['displayLines','showMolecules']:
                if not hasattr(self.vf, item):
                    self.vf.loadCommand('displayCommands', item, 'Pmv')
            #add callback to showMolecules
            self.vf.showMolecules.addCallback(self.checkMolLines)
            self.hasShowHide=1
        else:
            self.hasShowHide=0


    def guiCallback(self):
        """called each time the 'select docking logfile' button is pressed"""
        dlgFile = self.vf.askFileOpen(types=[('select docking log filename:', '*.dlg'),
                            ('all files','*')],
                            title = 'Docking Log File:')

        if dlgFile:
            addToPrevious = 0
            if self.vf.docked and len(self.vf.docked.dlo_list):
                dlo = self.vf.docked.dlo_list[-1]
                if dlo.filename:
                    t = ' Add to current docking containing ' + dlo.filename + '?'
                else:
                    t = ' Add to current docking containing ' + \
                            self.vf.docked.ligMol.name + '?'
                d = SimpleDialog(self.vf.GUI.ROOT, text = t,
                buttons = ["No", "Yes", "Cancel"], default = 0, 
                        title = "Add To Previous Docking: ")
                addToPrevious = d.go()
                if addToPrevious==2:
                    return 'ERROR'
            kw = {}
            ask = 1
            apply( self.doitWrapper, (dlgFile, addToPrevious,ask),  kw)


    def __call__(self, dlgFile, addToPrevious=0, ask=1, **kw):
        """None<-ADanalyze_readDLG(dlgFile, addToPrevious=0)
        \ndlgFile --- docking log file
        \naddToPrevious --- if 0, start a new Docking otherwise add docked conformations to a previous docking
        """
        apply(self.doitWrapper, (dlgFile, addToPrevious, ask,), kw) 


    def doit(self, dlgFile, addToPrevious, ask):
        """dlgFile: log of a docking result"""
        #CHECK FOR XML extension
        ext = os.path.splitext(os.path.basename(dlgFile))[-1]
        if ext=='.xml':
            from AutoDockTools.XMLParser import XMLParser
            p = XMLParser()
            if not addToPrevious or not self.vf.docked:
                d = self.vf.docked = Docking(parser=p)
                print isinstance(d.parser, XMLParser)
            if not isinstance(self.vf.docked.parser, XMLParser):
                print "unable to add xml result to previous result"
                return "ERROR"
            d.readXMLResults(dlgFile)
        else:
            #d will be a new docking instance
            if not addToPrevious or not self.vf.docked:
                d = self.vf.docked = Docking()
            #set up rmsTool for the clusterer here
            #this is used by the spw in showStates
            d = self.vf.docked
            #now readDLG creates a new dlo and adds it to d.dlo_list
            d.readDlg(dlgFile)
        if hasattr(d, 'clusterer'):
            coords = d.ligMol.allAtoms.coords[:]
            d.clusterer.rmsTool = RMSDCalculator(coords)
            d.clusterer.rmsToolRef = '0'
            d.clusterer.usesSubset = 0
        #NB: thisDLO is a docking log object which should know its confs
        if not d.ligMol in self.vf.Mols:
            mol = self.vf.addMolecule(d.ligMol, ask=ask)
            if mol==None:
                oldPref = self.vf.userpref['Warning Message Format']['value']
                self.vf.setUserPreference(('Warning Message Format','pop-up'))
                msg = '\n' + dlgFile + ' skipped:\nUser chose not to build ligand molecule!!!'
                self.warningMsg(msg)
                self.vf.setUserPreference(('Warning Message Format', oldPref))
                return "ERROR"
            mol.colored = 0
            mol.docking = d

        #FIX THIS:
        d2 = self.vf.dockings
        if d not in d2.values():
            filename = os.path.basename(dlgFile)
            if filename in d2.keys():
                #fix this: could be >1 duplicate
                filename = filename + '_' + str(len(d2.keys()))
            d2[filename] = d
        else:
            ind = d2.values().index(d)
            filename = d2.keys()[ind]

        #if hasattr(self.vf.ADanalyze_showHistogram, 'canvas'):
        #    l = self.vf.ADanalyze_showHistogram.canvas.find_all()
        #    for item in l:
        #        self.vf.ADanalyze_showHistogram.canvas.delete(item)
        #        self.vf.ADanalyze_showHistogram.buildIt()

        modelLen = str(len(d.ch.conformations))

        if not addToPrevious:
            msg = 'Read '+ modelLen + ' docked conformations from '  + filename +'\n' 
            if modelLen>1: msg = msg + '\nUse Analyze->Conformations->Play... to view'
        else:
            msg = 'read ' + os.path.basename(dlgFile) +'\n\t->'+ modelLen
        if hasattr(d, 'warnings') and len(d.warnings):
            msg = msg +  "\n\nType 'self.docked.warnings' in python shell to list " + str(len(d.warnings))+ " warning message(s) from dlg"
        self.vf.warningMsg(msg)



ADGetDLGGUI=CommandGUI()
ADGetDLGGUI.addMenuCommand('AutoToolsBar', menuText['AnalyzeMB'],
        menuText['readDLG'], cascadeName = menuText['DockingLogMB'])
###ADGetDLGGUI.menuBarCfg.update({'background':'tan','relief':'sunken'})



class ADLoadVinaResult(MVCommand):
    """Allows the user to read a multi-model vina result 
    \nPackage : AutoDockTools
    \nModule  : autoanalyzeCommands
    \nClass   : ADLoadVinaResult
    \nCommand : ADanalyze_readVinaResult
    \nSynopsis:\n
        None<-ADanalyze_readVinaResult(pdbqtFile, asMolecules=0)
    \nRequired Arguments:\n    
        pdbqtFile --- vina result file
    \nOptional Arguments:\n    
        modelsAs --- if 'conformations', add conformations to a single molecule, else build separate molecules
    """

    def onAddCmdToViewer(self):
        if self.vf.hasGui:
            #print "in ADLVR: onAddCmdToViewer"
            self.first = 1
            

    def onRemoveObjectFromViewer(self, obj):
        if self.vf.docked and self.vf.docked==obj:
            self.vf.docked = None


    def processArrowEvent(self, event):
        #detect molecules with conformations and vina_results
        from Pmv.moleculeViewer import EditAtomsEvent
        mols = self.vf.Mols.get(lambda x: hasattr(x, 'vina_results') and len(x.allAtoms[0]._coords)>1)
        #print "pAE: mols=", mols
        if not len(mols):
            return
        #c = self.vf.GUI.VIEWER.currentCamera
        #c.focus_set()
        #c.Enter_cb()
        if hasattr(event, 'keysym') and event.keysym=='Right':
            #print "Right"
            for m in mols:
                #print m.name
                geom = m.geomContainer.geoms['ProteinLabels']
                ind = m.allAtoms[0].conformation
                nconf = len(m.allAtoms[0]._coords)
                ###if ind+1 < nconf:
                #print "R:ind=", ind
                if ind < nconf:
                    ###ind = ind+1
                    ###m.allAtoms.setConformation(ind)
                    ###event = EditAtomsEvent('coords', m.allAtoms)
                    ###self.vf.dispatchEvent(event)
                    #print "RIGHT: comparing ", str(m.vina_results[ind][0]), " with ", geom.labels
                    #if len(geom.vertexSet)==1 and float(geom.labels[0])>=m.vina_results[ind][0]:
                    #    print "in right if"
                        #geom.Set(labels = (str(m.vina_results[ind+1][0]),)) 
                        #self.vf.labelByProperty(mols, properties = ['vina_energy'], font={'fontScales':(.5,.5,.5)})
                    #self.vf.labelByProperty(mols, properties = ['vina_energy'], font={'fontScales':(.5,.5,.5)})
                    g = m.geomContainer.geoms['ProteinLabels']
                    energy = m.vina_results[ind][0]
                    #g.Set(labels = (str(energy),))
                    #print "RIGHT: vina_results[",ind,"][0]=", m.vina_results[ind][0]
                    g.Set(labels = (str(m.vina_results[ind][0]),))
                    if ind==0:
                        print "over-all best vina_energy=", m.vina_results[ind][0]
                    else:
                        print ind+1," vina_energy=", m.vina_results[ind][0]
                    if hasattr(m, 'bindingSite'):
                        delattr(m, 'bindingSite')
                        c = self.vf.ADanalyze_showBindingSite
                        c.build()
                else:
                    #print "OFF THE END!!"
                    print 'end: last vina result for ', m.name
                event = EditAtomsEvent('coords', m.allAtoms)
                self.vf.dispatchEvent(event)
                    #print " end "#, ind, "  vina_energy =", m.vina_results[-1][0]
        if hasattr(event, 'keysym') and event.keysym=='Left':
            #print "Left"
            for m in mols:
                #print m.name,
                geom = m.geomContainer.geoms['ProteinLabels']
                ind = m.allAtoms[0].conformation
                #print "L: ind=", ind
                ###if ind > 0:
                if ind >= 0:
                    ###ind = ind - 1
                    ###m.allAtoms.setConformation(ind)
                    ###event = EditAtomsEvent('coords', m.allAtoms)
                    ###self.vf.dispatchEvent(event)
                    g = m.geomContainer.geoms['ProteinLabels']
                    g.Set(labels = (str(m.vina_results[ind][0]),))
                    #print "LEFT: vina_results[",ind,"][0]=", m.vina_results[ind][0]
                    if ind==0:
                        print "over-all best vina_energy=", m.vina_results[ind][0]
                    else:
                        print ind+1," vina_energy=", m.vina_results[ind][0]
                    if hasattr(m, 'bindingSite'):
                        delattr(m, 'bindingSite')
                        c = self.vf.ADanalyze_showBindingSite
                        c.build()
                else:
                    #print "LESS THAN ZERO!!"
                    print 'first: best vina result for ', m.name
                    #print " at first"
                    #print " over-all best vina_energy=", m.vina_results[0][0]
                event = EditAtomsEvent('coords', m.allAtoms)
                self.vf.dispatchEvent(event)
        

    def guiCallback(self):
        """called each time the 'select docking logfile' button is pressed"""
        if not hasattr(self, 'modelsAs'):
            self.modelsAs = Tkinter.IntVar(master=self.vf.GUI.ROOT)

        pdbqtFile = self.vf.askFileOpen(types=[('select vina result filename:', '*.pdbqt'),
                            ('all files','*')],
                            title = 'AutoDock Vina Result File:')
        if pdbqtFile:
            ext = os.path.splitext(os.path.basename(pdbqtFile))[-1]
            if ext!='.pdbqt': 
                msg = "invalid file:" + pdbqtFile + ".  extension must be .pdbqt"
                self.vf.warningMsg(msg)
                return
            ifd = InputFormDescr(title="Load MODELS as:")
            ifd.append({'name': 'ModelsAsMols',
                'text': 'separate molecules',
                'widgetType':Tkinter.Radiobutton,
                'tooltip':'Check this button to add a separate molecule for each model',
                'variable': self.modelsAs,
                'value': '1',
                'text': 'Multiple molecules            ',
                'gridcfg': {'sticky':'w'}})
            ifd.append({'name': 'ModelsAsConfs',
                'widgetType':Tkinter.Radiobutton,
                'tooltip':'Check this button to add a single molecule\nwith a separate conformation for each model',
                'variable': self.modelsAs,
                'value': '0',
                'text': 'Single molecule with multiple conformations       ',
                'gridcfg': {'sticky':'w'}}) 
            d = self.vf.getUserInput(ifd)
            modelsAs = 'conformations' 
            ans = d['ModelsAsMols']
            if ans>0: modelsAs = 'molecules'
            #self.vf.readMolecule(pdbqtFile)
            apply( self.doitWrapper, (pdbqtFile, modelsAs), {})


    def __call__(self, pdbqtFile, modelsAs='conformations', **kw):
        """None<-ADanalyze_readVinaResult(pdbqtFile, modelsAs='conformations')
        \npdbqtFile --- vina result file
        \nmodelsAs --- if 'molecules', a new molecule is added for each result otherwise a single molecule with multiple conformations is built
        """
        apply(self.doitWrapper, (pdbqtFile, modelsAs), kw) 


    def doit(self, pdbqtFile, modelsAs='conformations', setupUpdates=1):
        """pdbqtFile: output of an AutoDock Vina docking
           modelsAs: build a single molecule with 'conformations' or multiple 'molecules'
        """
        #print "calling readMolecule with setupUpdates=", setupUpdates
        mols = self.vf.readMolecule(pdbqtFile, modelsAs=modelsAs, setupUpdates=setupUpdates)
        m = mols[0]
        if modelsAs=='conformations':
            mols.current_index = 0
            self.vf.dockings = {}
            self.vf.dockings[m.name] = m
            self.vf.docked = m
            m.dlo_list = []
            m.ligMol = m
            if self.first:
                #check for preexisting callback from MoleculeReader
                #eM = self.vf.GUI.VIEWER.currentCamera.eventManager
                #for key in ["<Right>", "<Left>"]:
                #    if eM.eventHandlers.has_key(key):
                #        eM.RemoveCallback(key,self.vf.readMolecule.processArrowEvent)
                #        print "removed readMolecule.processArrowEvent for ", key
                self.vf.GUI.addCameraCallback("<Right>", self.processArrowEvent)
                self.vf.GUI.addCameraCallback("<Left>", self.processArrowEvent)
                self.first = 0
            msg = "Added %d docked conformations to %s\nUse keyboard arrow keys to view" %( len(m.allAtoms[0]._coords), m.name)
            #msg = " Added %d conformations to %s" %( len(m.allAtoms[0]._coords) - 1, m.name)
            self.warningMsg(msg)
        #set up label by vina_energy etc
        self.vf.labelByProperty(mols, properties = ['vina_energy'], font={'fontScales':(.5,.5,.5)})
        print "over-all best vina_energy=", m.vina_results[0][0]

        



ADLoadVinaResultGUI=CommandGUI()
ADLoadVinaResultGUI.addMenuCommand('AutoToolsBar', menuText['AnalyzeMB'],
        menuText['readVinaResult'], cascadeName = menuText['DockingLogMB'])
###ADGetDLGGUI.menuBarCfg.update({'background':'tan','relief':'sunken'})



class ClusterDockingChooser(MoleculeChooser):


    def __init__(self, viewer, mode='single', title= 'Choose Clustered Docking'):
        MoleculeChooser.__init__(self,viewer, mode, title)
        self.ifd =self.ipf= InputFormDescr(title=title)


    def showConf(self, event=None):
        if not hasattr(self.vf, 'showMolecules'): 
            return
        ind = self.lb.curselection()
        if len(ind) and hasattr(self.vf, 'showMolecules'):
            curStr = self.lb.get(ind)
            if find(curStr, 'no conformation')>-1:
                return
            name = split(curStr)[0]
            #THIS won't WORK UNLESS 
            #mol previously had displayLines called on it
            #self.vf.showMolecules([name])
            mol = self.vf.Mols.NodesFromName(name)[0]
            self.vf.ADanalyze_readDLG.checkMolLines(mol)

        
    def go(self, modal=1, blocking=0, event='<ButtonRelease-1>'):
        """start the form"""
        d = self.vf.docked
        if not hasattr(d, 'clusterer'):
            self.vf.warningMsg('current docking has no clusterer')
            return
        if not len(d.clusterer.clustering_dict.keys()):
            self.vf.warningMsg('current clusterer has no clusterings')
            return
        #FIX THIS: add choice of possible keys
        #thisDLO = d.dlo_list[-1]
        #rmstol = self.vf.dpo['rmstol']['value']
        #rmstol = thisDLO.dpo['rmstol']['value']
        #FIX THIS:???
        rmstol = d.clusterer.clustering_dict.keys()[0]
        cluList = d.clusterer.clustering_dict[rmstol]
        #cluList = d.ch.clusterer.clustering_dict[rmstol]
        #cluList = d.ch.conformations
        runs = len(d.ch.conformations)
        #runs = thisDLO.dpo['runs']['value']
        #runs = d.dpo['runs']['value']
        entries=[]
        titleStr = 'select from ' + str(runs) + ' dockings:\n(double click to update coords)\n(Rank_SubRank              docked energy)'
        for i in range(len(cluList)):
            #self.headLine = "RANK "+str(i)+": mean energy ="+str(cluList[i].lowEnergy)
            for j in range(len(cluList[i])):
                rankStr = str(i+1)+'_'+str(j+1)
                docked=cluList[i][j]
                if self.vf.docked.version<4.0:
                    if hasattr(docked, 'mol') and hasattr(docked,'intermol_energy') and \
                                hasattr(docked, 'torsional_energy'):
                        freeEnergy = docked.intermol_energy + docked.torsional_energy
                        if i==0 and j==0: entries.append((docked.mol.name + ' ' + 'input', 'n/a'))
                        dockedMol = docked.mol.name + " " +rankStr + "                               "+str(docked.energy)
                        dockedStr= "Rank:                          "+rankStr+\
                                 "\nDocked Energy:         " +str(docked.energy) + \
                                 "\nCluster RMS:                " + str(docked.clRMS) + \
                                 "\nRef RMS:                     " + str(docked.refRMS)+ \
                                 "\nfreeEnergy :               "+str(freeEnergy)+\
                                 "\nkI :                           " +str(docked.inhib_constant) + \
                                 "\nInterMolecularEnergy: " + str(docked.intermol_energy) + \
                                 "\nInternal Energy :          " + str(docked.internal_energy)
                    else:
                        dockedMol = docked.mol.name + ' (only has docked energy)      '+ str(docked.energy)
                        dockedStr = "Rank:           "+str(i+1)+"\nDocked Energy: " +str(docked.energy)
                        if docked.clRMS is not None:
                            dockedStr = dockedStr + "\nCluster RMS:   " + str(docked.clRMS)
                            if docked.refRMS is not None:
                                dockedStr = dockedStr + "\nRef RMS:       " + str(docked.refRMS)
                else:   #version 4.0
                    if i==0 and j==0: entries.append((docked.mol.name + ' ' + 'input', 'n/a'))
                    freeEnergy = docked.intermol_energy + docked.torsional_energy
                    dockedMol = docked.mol.name + " " +rankStr + "                               "+str(docked.binding_energy)
                    dockedStr = "Rank: " + rankStr
                    try:
                        dockedStr = dockedStr + "\nBinding Energy:            " +str(docked.binding_energy)
                    except:
                        dockedStr = dockedStr + "\nBinding Energy:           unavailable"
                    try:
                        dockedStr = dockedStr + "\nkI :                               " +str(docked.inhib_constant) + docked.inhib_constant_units 
                    except:
                        dockedStr = dockedStr + "\nkI :                   unavailable" 
                    try:
                        dockedStr = dockedStr + "\nIntermolecular Energy :    "+str(docked.intermol_energy)
                    except:
                        dockedStr = dockedStr + "\nIntermolecular Energy :    unavailable"
                    try:
                        dockedStr = dockedStr + "\nInternal Energy :           " + str(docked.total_internal)
                    except:
                        dockedStr = dockedStr + "\nInternal Energy :    unavailable"
                    try:
                        dockedStr = dockedStr + "\nTorsional Energy :          " + str(docked.torsional_energy)
                    except:
                        dockedStr = dockedStr + "\nTorsional Energy :    unavailable" 
                    try:
                        dockedStr = dockedStr + "\nUnbound Extended Energy:  " + str(docked.unbound_energy)
                    except:
                        dockedStr = dockedStr + "\nUnbound Extended Energy:    unavailable"
                    try:
                        dockedStr = dockedStr + "\nCluster RMS:                 " + str(docked.clRMS)
                    except:
                        dockedStr = dockedStr + "\nCluster RMS:        unavailable"
                    try:
                        dockedStr = dockedStr + "\nRef RMS:                      " + str(docked.refRMS) 
                    except:
                        dockedStr = dockedStr + "\nRef RMS:     unavailable"
                entries.append((dockedMol,dockedStr))

        self.ifd.insert(0,{'name':'Molecule',
                            'widgetType': 'ListChooser',
                            'title':titleStr,
                            'command': self.showConf,
                            'commandEvent': "<ButtonPress-1>",
                            'mode':self.mode,
                            'entries':entries})
        if not (modal or blocking):
            #self.ifd.append({'widgetType':Tkinter.Button,
                            #'text':'Show Conformation',
                            #'gridcfg':{'sticky': Tkinter.E+Tkinter.W,'columnspan':2},
                            #'command':self.show_cb})
            self.ifd.append({'widgetType':Tkinter.Button,
                            'text':'Write Current Coords',
                            'gridcfg':{'sticky': Tkinter.E+Tkinter.W,'columnspan':2},
                            'command':self.save_cb})
            self.ifd.append({'widgetType':Tkinter.Button,
                            'text':'Dismiss',
                            'gridcfg':{'sticky': Tkinter.E+Tkinter.W,'columnspan':2},
                            'command':self.done_cb})
        #self.bindPickingCallback(event)
        from Pmv.picker import AtomPicker
        self.ap=AtomPicker(self.vf, None,0,callbacks=[self.onPick],immediate=1)
        self.ap.go(modal=0)
        val=self.vf.getUserInput(self.ifd, modal=modal, blocking=blocking)
        self.ifd.entryByName['Molecule']['widget'].comments.config(height=9)
        if val:
            if modal or blocking:
                ##self.unbindPickingCallback()
                molNames=''
                for m in val['Molecule']: 
                    mlist=split(m)
                    if len(mlist)==2:
                        molNames=molNames+mlist[0]
                mols = self.vf.Mols.NodesFromName(molNames)
                self.ap.stop()
                if self.mode=='single': return mols[0]
                return mols
            else:
                self.form = val
                self.lb = self.ifd.entryByName['Molecule']['widget'].lb
                #this won't work because it preempts curselection...
                #self.lb.bind("<Button-1>",self.show_cb, add='+')
                self.lb.bind("<Double-Button-1>",self.show_cb)
                return val
        else:
            return val


    def done_cb(self,event=None):
        self.ap.stop()
        self.ifd.form.root.withdraw()


    def show_cb(self, event=None):
        c=self.vf.ADanalyze_chooseDockedConformations
        lc=c.chooser.ifd.entryByName['Molecule']['widget']
        curselection = lc.lb.curselection()
        if len(curselection)==0:
            t="Select Docked Conformation First"
            self.vf.warningMsg(t)
            return
        ind = int(curselection[0])
        #ind = int(curselection[0]) - 1
        dockedName = split(lc.entries[ind][0])[0]
        dockedMol = self.vf.Mols.NodesFromName(dockedName)[0]
        if ind == 0:
            dockedMol.allAtoms.setConformation(0)
        else:
            confName = split(lc.entries[ind][0])[1]
            #names are '1_1' etc
            clu, rank = map(int, split(confName,'_'))
            clusterer = self.vf.docked.clusterer
            #FIX THIS: it lets you only look at first rms clustering
            clu_keys = clusterer.clustering_dict.keys()
            if not len(clu_keys):
                t = 'No Clusterings available:Make a clustering first!'
                self.vf.warningMsg(t)
                return
            rmstol = clusterer.clustering_dict.keys()[0]
            clustering = clusterer.clustering_dict[rmstol]
            conf = clustering[clu-1][rank-1]
            dockedMol.docking.ch.set_conformation(conf)

        if not self.vf.hasGui: return
        event = EditAtomsEvent('coords', dockedMol.allAtoms)
        self.vf.dispatchEvent(event)
        #modEvent = ModificationEvent('edit','coords', dockedMol.allAtoms)
        #dockedMol.geomContainer.updateGeoms(modEvent)
        #self.vf.GUI.VIEWER.Redraw()
        
        #if find(dockedName, 'noconformation')>-1:
        #    msg='no conformation for this curselection'
        #    self.vf.warningMsg(msg)
        #    return
        #dockedMol = self.vf.Mols.NodesFromName(dockedName)[0]
        #for item in self.vf.Mols:
        #    if item != dockedMol:
        #        #undisplay the rest EXCEPT for macromolecule, maybe
        #        if item.name != self.vf.docked.macroStem:
        #            ##cbvar=self.vf.showMolecules.molList[item.name]
        #            ##cbvar.set(0)
        #            self.vf.showMolecules([item.name], negate = 1)
        #    else:
        #        #if len(item.geomContainer.atoms['lines'])==0:
        #        if len(item.geomContainer.atoms['bonded'])==0:
        #            item.buildBondsByDistance()
        #            self.vf.displayLines(item)
        #        #cbvar=self.vf.showMolecules.molList[item.name]
        #        #cbvar.set(0)
        #        self.vf.showMolecules([item.name], negate = 0)
        

    def save_cb(self, event=None):
        ligandlines = self.vf.docked.ligMol.parser.allLines
        sel = self.vf.getSelection()
        if len(sel)==0:
            t="Select Docked Conformation First"
            self.vf.warningMsg(t)
            return

        #docked = sel[0]
        ##FIX THIS
        #dockedlines = docked.allLines
        ###get output file, open it and write
        savefile=self.vf.askFileSave(types=[('docked ligand file','*.docked.pdbq'),('autodock file','*.pdbq'),('any','*.*')],
                title = 'Docked Ligand File:')
        if savefile:
            outf=open(savefile, 'w')
            ctr=0
            for i in range(len(ligandlines)):
                item = ligandlines[i]
                print item
                if find(item,'ATOM')!=0 and find(item, 'HETATM') !=0:
                    outstring = item + '\n'
                else:
                    #have to build coord str
                    #SHOULD THESE BE KEPT FOR SOME REASON???
                    #dockeditem = dockedlines[ctr]
                    #thisDLO = self.vf.docked.dlo_list[-1]
                    #coords = thisDLO.ligMol.allAtoms[ctr].coords
                    coords = self.vf.docked.ligMol.allAtoms[ctr].coords
                    cStr = '%8.3f%8.3f%8.3f'%(coords[0],coords[1],coords[2])
                    outstring = item[:30] + cStr + item[54:] + '\n'
                    #outstring = item[:30] + dockeditem[30:54]+item[54:]
                    ctr=ctr+1
                outf.write(outstring)
            outf.close()


    def getMolSet(self):
        """method to get currently selected molecules when the chooser is used
        in modal=0 and blocking=0 mode"""
        
        val = self.form.checkValues()
        molNames = ""
        for m in val['Molecule']: 
            mlist=split(m)
            if len(mlist)==2:
                molNames=molNames+mlist[0]
        if len(molNames)>0:
            mols = self.vf.Mols.NodesFromName( molNames )
            if self.mode=='single': return mols[0]
            return mols


    def onPick(self,atoms):
        listChooser = self.ifd.entryByName['Molecule']['widget']
        tkListBox = listChooser.lb
        if len(atoms):
            pickedMol = atoms[0].top
            #then need to make pickedMol the selection in self.lc
            for i in range(len(listChooser.entries)):
                listChooserlist=split(listChooser.entries[i][0])
                if pickedMol.name == listChooserlist[0]:
                    self.pickedMolIndex= i
                    tkListBox.select_clear(0,'end')
                    listChooser.select(i)
                    return
            t="error: %s not a docked conformation" %pickedMol.name
            self.vf.warningMsg(t)



class ModelDockingChooser(MoleculeChooser):


    def __init__(self, viewer, mode='single', title= 'Choose Docking Model'):
        MoleculeChooser.__init__(self,viewer,mode,title)
        self.ifd =self.ipf= InputFormDescr(title=title)


    def go(self, modal=1, blocking=0, event='<ButtonRelease-1>'):
        """start the form"""
        d = self.vf.docked.dlo_list[-1]
        #d = self.vf.docked
        modelList=d.ch.modelList
        runs = d.dlo_list[0].dpo['runs']['value']
        entries=[]
        titleStr = 'select from ' + str(runs) + ' dockings:\n(double click to update coords)'
        for i in range(len(modelList)):
            docked=modelList[i]
            dockedMol = docked.mol.name + "                                  "+str(docked.freeEnergy)
            dockedStr= "freeEnergy :           "+str(docked.freeEnergy)+"\nkI :                     " +str(docked.kI) + "\nInterMolecularEnergy:   " + str(docked.interMolEnergy) + "\nInternal Energy :       " + str(docked.internalEnergy)
            entries.append((dockedMol,dockedStr))
            ##if i!=len(cluList)-1:
                ##entries.append('__________________________________________')

        self.ifd.insert(0,{'name':'Molecule',
                            'widgetType': 'ListChooser',
                            'title':titleStr,
                            'mode':self.mode,
                            'entries':entries})
        if not (modal or blocking):
            self.ifd.append({'widgetType':Tkinter.Button,
                            'text':'Dismiss',
                            'command':self.done_cb})
        self.bindPickingCallback(event)
        val=self.vf.getUserInput(self.ifd, modal=modal, blocking=blocking)
        self.ifd.entryByName['Molecule']['widget'].comments.config(height=9)
        if val:
            if modal or blocking:
                self.unbindPickingCallback()
                molNames=''
                for m in val['Molecule']: 
                    mlist=split(m)
                    if len(mlist)==2:
                        molNames=molNames+mlist[0]
                mols = self.vf.Mols.NodesFromName(molNames)
                if self.mode=='single': return mols[0]
                return mols
            else:
                self.form = val
                return val
        else:
            return val


    def done_cb(self,event=None):
        self.ifd.form.root.withdraw()


    def getMolSet(self):
        """method to get currently selected molecules when the chooser is used
        in modal=0 and blocking=0 mode"""
        
        val = self.form.checkValues()
        molNames = ""
        for m in val['Molecule']: 
            mlist=split(m)
            if len(mlist)==2:
                molNames=molNames+mlist[0]
        if len(molNames)>0:
            mols = self.vf.Mols.NodesFromName( molNames )
            if self.mode=='single': return mols[0]
            return mols


    def onPick(self,pick):
        listChooser = self.ifd.entryByName['Molecule']['widget']
        tkListBox = listChooser.lb
        atoms = self.vf.findPickedAtoms(pick)
        if len(atoms):
            pickedMol = atoms[0].top
            #then need to make pickedMol the selection in self.lc
            for i in range(len(listChooser.entries)):
                listChooserlist=split(listChooser.entries[i][0])
                if pickedMol.name == listChooserlist[0]:
                    self.pickedMolIndex= i
                    tkListBox.select_clear(0,'end')
                    listChooser.select(i)
                    return
            t= "error: %s not a docked conformation" %pickedMol.name
            self.vf.warningMsg(t)


colorsList =['black','brown','blue','red']



class ADDrawHistogram(MVCommand):

    def onAddCmdToViewer(self):
        checkHasInitializedDockings(self.vf)


    def guiCallback(self):
        #could have info attached to a molecule
        if not self.vf.docked:
            self.vf.warningMsg('Please Read a Docking Log First')
            return
        if not hasattr(self, 'ifd'):
            self.buttons = {}
            ifd = self.ifd= InputFormDescr(title='Visualize Results Histogram')
            ifd.append({'name': 'histCanvas',
                        'widgetType':Tkinter.Canvas,
                        'wcfg':{'width':450,
                            'height':450,
                            'bg':'white'},
                        'gridcfg':{'sticky':Tkinter.W, 'columnspan':4}})
            ifd.append({'name': 'quitButton',
                        'widgetType':Tkinter.Button,
                        'text':'Close',
                        'command': self.close,
                        'gridcfg':{'sticky':Tkinter.W+Tkinter.E}})
            ifd.append({'name': 'resetButton',
                        'widgetType':Tkinter.Button,
                        'text':'Reset Ligand',
                        'command': self.reset,
                        'gridcfg':{'sticky':Tkinter.W+Tkinter.E,'row':-1,'column':1}})
            ifd.append({'name': 'infoButton',
                        'widgetType':Tkinter.Button,
                        'text':'Info...',
                        'command': self.info,
                        'gridcfg':{'sticky':Tkinter.W+Tkinter.E,'row':-1,'column':2}})
            ifd.append({'name': 'saveButton',
                        'widgetType':Tkinter.Button,
                        'text':'Write as postscript ',
                        'command': self.savePS,
                        'gridcfg':{'sticky':Tkinter.W+Tkinter.E,'row':-1,'column':3}})
            self.form = self.vf.getUserInput(self.ifd, modal=0,blocking=0)
            self.form.root.protocol('WM_DELETE_WINDOW',self.close)
            self.canvas=self.ifd.entryByName['histCanvas']['widget']
            Tkinter.Widget.bind(self.canvas,"<Button-1>", self.MouseDown)
            Tkinter.Widget.bind(self.canvas,"<Button1-ButtonRelease>", self.MouseUp)
            Tkinter.Widget.bind(self.canvas,"<ButtonPress-2>", self.showInfo)
            Tkinter.Widget.bind(self.canvas,"<ButtonRelease-2>", self.hideInfo)
            ifd2=self.ifd2=InputFormDescr(title='Docked Conformation Info')
            ifd2.append({'name': 'histText',
                        'widgetType':Tkinter.Text,
                        'wcfg':{'width':60,
                            'height':10, 
                            'bg':'white'},
                        'gridcfg':{'sticky':Tkinter.W, 'columnspan':3}})
            self.form2 = self.vf.getUserInput(self.ifd2, modal=0,blocking=0)
            self.histText=self.ifd2.entryByName['histText']['widget']
            self.form2.root.withdraw()
        else:
            self.form.deiconify()
        self.buildIt()


    def savePS(self,event=None):
        thisDLO = self.vf.docked.dlo_list[-1]
        dockingName = os.path.basename( thisDLO.filename)
        #dockingName = os.path.basename( self.vf.docked.dlgFile)
        tstr = dockingName +' Histogram Postscript File'
        savefile= self.vf.askFileSave(types=[('postscript files','*.ps')],
                title = tstr)
        if savefile:
            self.canvas.postscript({'file':savefile,'colormode':'mono'})


    def info(self, event=None):
        msg = "Click rectangle with B-2 for docked stats...\nClick with B-1 to set ligand coords \n(Color changes for delta Lowest Energy>2Kcal)"
        self.vf.warningMsg(msg)
        self.form.lift()


    def getDocking(self,dname):
        ans = None
        l = split(dname, '_')
        i = int(l[0])
        j = int(l[1])
        #FIX THIS
        #thisDLO = self.vf.docked.dlo_list[-1]
        cluster_list = self.vf.docked.clusterer.clustering_dict.values()[0]
        return cluster_list[i][j]
        #return self.vf.docked.cluster_list[i][j]
        #for item in self.vf.docked.ch.conformations:
            #for m in item.members:
                #if m.name==dname:
                    #return m
        #if not ans:
            #msg = dname + ' not a docked conformation'
            #self.vf.warningMsg(msg)
            #return None


    def hideInfo(self, event):
        self.form2.root.withdraw()


    def showInfo(self, event):
        self.lastx = self.startx = self.canvas.canvasx(event.x)
        self.lasty = self.starty = self.canvas.canvasy(event.y)
        self.currentObject= self.canvas.find_closest(self.startx,self.starty)[0]
        if self.currentObject in self.buttons.keys():
            nstr=self.buttons[self.currentObject]
            docked = self.getDocking(nstr)
            if docked:
                self.histText.delete(0.0,'end')
                if not hasattr(docked, 'run'):
                    docked.run = None
                msg = 'name: '+docked.mol.name+'\nenergy: '
                msg = msg + str(docked.energy)+'\nrefRMS: '
                msg = msg + str(docked.refRMS)+'\nclusterRMS: '
                msg = msg + str(docked.clRMS) + '\nrun:' + str(docked.run)
                self.histText.insert('end',msg)
                self.form2.deiconify()
                top1=self.histText.winfo_toplevel()
                oldgeom=top1.geometry()
                geomlist=split(oldgeom,'+')
                newgeom=geomlist[0]+'+100+50'
                top1.geometry(newgeom)
        

    def MouseDown(self, event):
        self.currentObject = None
        self.lastx = self.startx = self.canvas.canvasx(event.x)
        self.lasty = self.starty = self.canvas.canvasy(event.y)
        self.selObj= self.canvas.find_closest(self.startx,self.starty)[0]
        #FIX THIS:
        docked = self.vf.docked
        ch = docked.ch
        ligMol = docked.ligMol
        if self.selObj in self.buttons.keys():
            self.canvas.itemconfig(self.selObj, width=2)
            nstr=self.buttons[self.selObj]
            nl=map(lambda x:x.name, self.vf.Mols)
            cluNum, rank = map(int,split(nstr, '_'))
            cluster = docked.clusterer.clustering_dict.values()[0][cluNum]
            #cluster = docked.cluster_list.data[cluNum]
            conf = cluster[rank]
            allAtoms = ligMol.allAtoms
            allAtoms.setConformation(docked.ch.confIndex)
            allAtoms.updateCoords(conf.getCoords())
            #docked.ch.set(docked.ch.conformations.index(conf))
            event = EditAtomsEvent('coords', ligMol.allAtoms)
            self.vf.dispatchEvent(event)
            #modEvent = ModificationEvent('edit','coords', ligMol.allAtoms)
            #ligMol.geomContainer.updateGeoms(modEvent)
            #self.vf.GUI.VIEWER.Redraw()
            
            # eg 1_1
            #if nstr in nl: 
                #self.vf.directSelect(self.buttons[self.selObj])
            #docked = self.getDocking(nstr)
            #if docked:
            #msg = 'energy: '+str(conf.energy)+'\nrefRMS: '+str(conf.refRMS)+'\nclRMS: ' + str(conf.clRMS)
            #print msg
            #self.vf.warningMsg(msg)
            self.canvas.itemconfig(self.selObj, width=1)
            self.form.lift()
            

    def MouseUp(self, event):
        self.lastx = self.canvas.canvasx(event.x)
        self.lasty = self.canvas.canvasy(event.y)
        if self.selObj and self.selObj in self.buttons.keys():
            self.canvas.itemconfig(self.selObj, width=1)


    def getEnergyText(self,i):
        #FIX THIS
        cluster = self.vf.docked.clusterer.clustering_dict.values()[0][i]
        #cluster = self.vf.docked.cluster_list[i]
        #cluster = self.vf.docked.ch.conformations[i]
        if i==0:
            self.clCtr=0
            self.curEnergy=cluster[0].docking_energy
            #self.curEnergy=cluster[0].energy
        elif cluster[0].docking_energy-self.curEnergy>2:
            #elif cluster[0].energy-self.curEnergy>2:
            self.curEnergy=cluster[0].docking_energy
            #self.curEnergy=cluster[0].energy
            self.clCtr=(self.clCtr+1)%4
        colorText=colorsList[self.clCtr]
        cluster.colorText=colorText
        xcoord=self.xoffset+(i+0.5)*self.clusterWidth
        if self.clusterWidth<16:
            if i%10==0:
                self.canvas.create_text(xcoord,self.yoffset+10, text=i, fill=colorText)
        elif self.clusterWidth <50:
            self.canvas.create_text(xcoord,self.yoffset+10, text=str(i+1), fill=colorText)
        else:
            self.canvas.create_text(xcoord,self.yoffset+10, text=str(cluster[0].docking_energy), fill=colorText)
            #self.canvas.create_text(xcoord,self.yoffset+10, text=str(cluster[0].energy), fill=colorText)
        

    def getButton(self,docked,i,j):
        #docked is cluster[i].member[j]
        #FIX THIS
        cluster = self.vf.docked.clusterer.clustering_dict.values()[0][i]
        colorText = cluster.colorText
        #colorText = self.vf.docked.ch.conformations[docked.rank-1].colorText
        y0=self.yoffset-j*self.indHeight
        y1=y0-self.indHeight
        x0=self.xoffset+i*self.clusterWidth
        x1=x0+self.clusterWidth
        key =self.canvas.create_rectangle(x0,y0,x1,y1,fill = 'white',outline=colorText)
        self.buttons[key] = str(i)+'_'+str(j)
        tstr=str(i+1) + "_"+str(j+1)
        if self.clusterWidth>50:
            self.canvas.create_text(x0+0.5*self.clusterWidth,y0-0.5*self.indHeight, text=tstr,fill=colorText)


    def buildIt(self):
        if not self.vf.docked:
            self.vf.warningMsg('Please Read a Docking Log First')
            return
        #FIX THIS
        if hasattr(self.vf.docked, 'hasClusters') and not \
                    self.vf.docked.hasClusters:
            self.vf.warningMsg('Please cluster Docking First')
            return
        if hasattr(self.vf.docked, 'clusterer') and \
                len(self.vf.docked.clusterer.clustering_dict.keys())==0:
            return
        #delete oldstuff first??
        l= self.canvas.find_all()
        for item in l:
            self.canvas.delete(item)
        colorCtr=0
        #FIX THIS
        cList=self.vf.docked.clusterer.clustering_dict.values()[0]
        #cList=self.vf.docked.cluster_list
        if len(cList):
            #determine how long canvas has to be...
            clen=len(cList)
            self.xoffset=10
            self.yoffset=360
            self.clusterWidth = int(400./float(clen))
            #determine how wide canvas has to be...
            #num_members = self.getArray(cList)
            #maxNumMembers = Numeric.maximum.reduce(num_members)
            num_members = map(lambda x: len(x), cList)
            maxNumMembers = max(num_members)
            num_membersSum = Numeric.add.reduce(Numeric.array(num_members))
            nWidth = 100 + maxNumMembers*40
            self.indHeight=int(330./float(maxNumMembers))
            if self.clusterWidth<30:
                self.clusterWidth=30
            if self.indHeight>30:
                self.indHeight=30
            if self.indHeight<10:
                self.indHeight=15
                newheight=15*maxNumMembers+70
                self.yoffset=newheight-40    
                self.canvas.config({'height':newheight})
            totalwidth=self.clusterWidth*clen+20
            if totalwidth>400:
                self.canvas.config({'width':totalwidth})
            thisDLO = self.vf.docked.dlo_list[-1]
            msg = "Current Docking:\n" + thisDLO.filename
            #msg = "Current Docking:\n" + self.vf.docked.dlgFile
            self.canvas.create_text(200,15,text=msg)
            rmstol = thisDLO.dpo['rmstol']['value']
            #rmstol = self.vf.docked.dpo['rmstol']['value']
            if rmstol:
                msg = "Clusters (rmstol-tolerance="+str(rmstol)+")"
                self.canvas.create_text(200,385,text=msg)
            for i in range(len(cList)):
                self.getEnergyText(i)
                for j in range(len(cList[i])):
                    self.getButton(cList[i][j],i,j)
        else:
            valList=[(0,0,100,200),(100,0,150,175)]
            for item in valList:
                self.canvas.create_rectangle(item[0],item[1],item[2],item[3], fill ='green')


    def reset(self):
        if self.vf.docked and hasattr(self.vf.docked,'ligMol'):
            mol = self.vf.docked.ligMol
            allAtoms = mol.allAtoms
            allAtoms.setConformation(0)
            event = EditAtomsEvent('coords', allAtoms)
            self.vf.dispatchEvent(event)
            #modEvent = ModificationEvent('edit','coords', allAtoms)
            #mol.geomContainer.updateGeoms(modEvent)
            #self.vf.GUI.VIEWER.Redraw()


    def close(self):
        self.form.root.withdraw()


ADDrawHistogramGUI=CommandGUI()
ADDrawHistogramGUI.addMenuCommand('AutoToolsBar', menuText['AnalyzeMB'],
        menuText['showHistogramMB'], cascadeName=menuText['StatesMB'])
            


class ADMacroLigandChart(MVCommand):

    def onAddCmdToViewer(self):
        self.macros = []
        self.ligands = []
        self.labels={}
        self.dockingsShown=[]


    def close(self,event=None):
        self.form.root.withdraw()


    def showOptions(self,ifd, event=None):
        row_w=ifd.entryByName['rowChoices']['widget']
        row_w.gridcfg=ifd.entryByName['rowChoices']['gridcfg']
        column_w=ifd.entryByName['columnChoices']['widget']
        column_w.gridcfg=ifd.entryByName['columnChoices']['gridcfg']
        #if sorting by row, use a column as the key and vice versa
        if self.sortChoice=='row':
            column_w.grid(column_w.gridcfg)
            row.grid_forget()
        else:
            row_w.grid(row_w.gridcfg)
            column_w.grid_forget()


    def sortInfo_cb(self,event=None):
        #allow user to specify row or column
        #for sort: use rows' 1st entry
        #print "calling sortByRow w/column1"
        if not hasattr(self, 'sortChoice'):
            self.sortChoice=Tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.sortChoice.set('row')
        if not hasattr(self, 'rowChoice'):
            self.rowChoice=Tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.rowChoice.set('1')
            self.columnChoice=Tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.columnChoice.set('1')
        ifd = InputFormDescr(title='Sort Options')
        ifd.append({ 'widgetType':Tkinter.Label,
                'wcfg':{'text':'   '},
                #'wcfg':{'text':'On Row or Column'},
                'gridcfg':{'sticky':Tkinter.W+Tkinter.E}})
        ifd.append({'name': 'rowBut',
                'widgetType':Tkinter.Radiobutton,
                'wcfg':{'text':'row', 
                    'value':'row',
                    'variable':self.sortChoice},
                'gridcfg':{'sticky':Tkinter.E}})
        #ifd.append({'name': 'colBut',
                #'widgetType':Tkinter.Radiobutton,
                #'wcfg':{'text':'column', 
                    ##'command':CallBackFunction(self.showOptions,ifd), 
                    #'value':'column',
                    #'variable':self.sortChoice},
                #'gridcfg':{'sticky':Tkinter.W, 'row': -1, 'column':1}})
        ifd.append({ 'widgetType':Tkinter.Label,
                'wcfg':{'text':'Using as Keys:'},
                'gridcfg':{'sticky':Tkinter.W+Tkinter.E}})
        columnList=map(lambda x: 'column' + str(x),
                    range(1,len(self.ifd.entryByColumn)))
        ifd.append({'name': 'columnChoices',
                'widgetType':Tkinter.Radiobutton,
                'listtext':columnList,
                'defaultValue':self.columnChoice.get(),
                'groupedBy':10,
                'wcfg':{'variable':self.columnChoice},
                'gridcfg':{'sticky':Tkinter.E}})
        #rowList=map(lambda x: 'row'+str(x),range(1,len(self.ifd.entryByRow)))
        #ifd.append({'name': 'rowChoices',
                #'widgetType':Tkinter.Radiobutton,
                #'listtext':rowList,
                #'defaultValue':self.rowChoice.get(),
                #'groupedBy':10,
                #'wcfg':{'variable':self.rowChoice},
                #'gridcfg':{'sticky':Tkinter.W,'row':3, 'column':1}})
        vals = self.vf.getUserInput(ifd)
        if self.sortChoice.get()=='row':
            direction='row'
            p=int(self.columnChoice.get()[-1])
        else:
            #t= "column-based sort under developement"
            #self.vf.warningMsg(t)
            #return
            direction='column'
            p=int(self.rowChoice.get()[-1])
        self.sortByRowColumn(direction,p)


    def sortByRowColumn(self,direction,  keyNum):
        if direction=='row':
            dict1=self.ifd.entryByRow
            dict2=self.ifd.entryByColumn
        else:
            dict1=self.ifd.entryByColumn
            dict2=self.ifd.entryByRow
        d1len=len(dict1.keys())
        d2len=len(dict2.keys())
        d=len(dict1.keys())-len(dict2[keyNum])
        for i in range(d): 
            dict2[keyNum].append([])
        d2len=len(dict2.keys())
        wList=dict2[keyNum]
        if len(wList)<d2len:wList.append([])
        keyList=[]
        sList=[]
        minVal='999999'
        for item in wList[1:]:
            if len(item)>0: 
                key=item[0].cget('text')
                try:
                    key1=float(key)
                except:
                    key1=key
                keyList.append(key1)
                sList.append(key1)
            else:
                keyList.append(minVal)
                sList.append(minVal)
                minVal=minVal+'9'
        keyList.insert(0,wList[0])
        sList.sort()
        sList.insert(0,wList[0])
        #mapList : each entry at position 'i' is position in origianl of current entry
        #eg [x,a,e,i]->[x,i,a,e] mapList=[0,3,1,2]
        mapList=[]
        for item in sList:mapList.append(keyList.index(item))
        #now index of keyList in sList shows how to repack widgets...
        self.newDict1=self.updateDict1(dict1,mapList)
        self.newDict2=self.updateDict2(dict2,mapList)
        #for i in range(len(dict1.keys())):
            #self.newDict1[i]=[]
        #self.newDict1[0]=dict1[0]
        for i in range(len(dict1.keys())):
            newRow=sList.index(keyList[i])
            for item in dict1[i]:
                if item!=[]:
                    item[0].oldgridcfg[direction]=newRow
                    item[0].grid(item[0].oldgridcfg)
        #if direction=='row':
            #self.ifd.entryByRow=self.newDict1
            #self.ifd.entryByColumn=self.newDict2
        #else:
            #self.ifd.entryByColumn=self.newDict1
            #self.ifd.entryByRow=self.newDict2


    def updateDict1(self, dict1, mapList):
        newDict={}
        for i in range(len(dict1.keys())):
            newindex=mapList[i]
            newDict[i]=dict1[newindex]
        return newDict


    def updateDict2(self, dict2, mapList):
        newDict={}
        for i in range(len(dict2.keys())):
            itemList=dict2[i]
            cutoff=len(itemList)-1
            newitemList=[]
            for j in range(len(itemList)):
                newindex= mapList[j]
                if newindex>cutoff:
                    newitemList.append([])
                else:
                    newitemList.append(itemList[newindex])
            newDict[i]=newitemList
        return newDict
            

    def changeInfo(self,event=None):
        t="Select new label for buttons:"
        d = SimpleDialog(self.vf.GUI.ROOT, text = t,
              buttons = ["lowest energy","RMSD", "number in lowest energy cluster"], default = 0, title = "Possible Labels: ")
        ok = d.go()
        llist=['lowEnergy','RMSD','cluNum']
        newLabel = llist[ok]
        for docking in self.dockingsShown:
            if hasattr(docking, 'clusterer') and hasattr(docking.clusterer, 'clustering_dict'):
                key = docking.clusterer.clustering_dict.keys()[0]
                conf = len(docking.clusterer.clustering_dict[key][0][0])
                clLen = len(docking.clusterer.clustering_dict[key][0])
            else:
                continue
            if ok == 0: 
                label=str(conf.docking_energy)
                self.ifd.entryByName['infoButton']['widget'].config(text="'Lowest Energy' labels...")
            elif ok==1:
                label=str(conf.refRMS)
                self.ifd.entryByName['infoButton']['widget'].config(text="'reference RMSD' labels...")
            else:
                label=str(clLen)
                self.ifd.entryByName['infoButton']['widget'].config(text="Number in Cluster 1")
            docking.button.config(text=label)
        

    def showInfo(self, docked, event=None):
        self.infoText.delete(0.0,'end')
        conf0 = docked.ch.conformations[0]
        msg = 'name: '+docked.dlo_list[0].filename+'\nenergy: '+str(conf0.lowEnergy)+'\nrefRMSD: '+str(conf0.refRMSD)+'\nclusterRMSD: ' + str(conf0.clRMSD)
        self.infoText.insert('end',msg)
        self.form2.deiconify()
        

    def hideInfo(self, docking, event=None):
        self.form2.root.withdraw()
        

    def savePS(self,filename=None,event=None):
        if not filename:
            filename= self.vf.askFileSave(types=[('docked ligand file','*.docked.pdbq')],
                title = 'Chart File:')
        if filename:
            outf=open(filename, 'w')
        else: return
        maxrow=len(self.labels.keys())-1
        #rowKeys=self.labels.keys()
        #rowKeys.sort()
        #maxrow=rowKeys[-1]
        maxcolumn=0
        for item in self.labels.values():
            #each item is a dictionary w/ column keys
            itemKeys=item.keys()
            itemKeys.sort()
            maxitem=itemKeys[-1]
            if maxitem>maxcolumn: 
                maxcolumn=maxitem
        #now build print out maxrow of maxcolumn+1 
        spaceStr="                      "
        for i in range(maxrow+1):
            columnLabels=self.labels[i]
            outstr=''
            for j in range(maxcolumn+1):
                if columnLabels.has_key(j):
                    outstr="%s%18.18s"%(outstr,columnLabels[j])
                else:
                    outstr="%s%18.18s"%(outstr,spaceStr)
            outstr=outstr+"\n"
            outf.write(outstr)
        outf.close()
                

    def buildForm(self):
        self.buttons = {}
        #self.labels={}
        self.labels[0]={}
        self.labels[0][0]='Ligands\Macros'
        #self.dockingsShown=[]
        self.rowctr=1
        self.columnctr=1
        #self.macros=[]
        #self.ligands=[]
        ifd = self.ifd= InputFormDescr(title='Results Chart')
        self.bButtons=[['quitButton','infoButton','sortButton','saveButton']]
        ifd.append({'name': 'cornerLabel',
                    'widgetType':Tkinter.Label,
                    'wcfg':{'text':'Ligands\Macros',
                            'width':16},
                    'gridcfg':{'sticky':Tkinter.W+Tkinter.E}})
        ifd.append({'name': 'quitButton',
                    'widgetType':Tkinter.Button,
                    'wcfg':{'text':'Close'},
                    'command': self.close,
                    'gridcfg':{'sticky':Tkinter.W+Tkinter.E}})
        ifd.append({'name': 'infoButton',
                    'widgetType':Tkinter.Button,
                    'wcfg':{'text':"'Lowest Energy' labels"},
                    'command': self.changeInfo,
                    'gridcfg':{'sticky':Tkinter.W+Tkinter.E,'row':-1,'column':1}})
        ifd.append({'name': 'sortButton',
                    'widgetType':Tkinter.Button,
                    'wcfg':{'text':"Sort"},
                    'command': self.sortInfo_cb,
                    'gridcfg':{'sticky':Tkinter.W+Tkinter.E,'row':-1,'column':2}})
        ifd.append({'name': 'saveButton',
                    'widgetType':Tkinter.Button,
                    'wcfg':{'text':'Write as postscript '},
                    'command': self.savePS,
                    'gridcfg':{'sticky':Tkinter.W+Tkinter.E,'row':-1,'column':3}})
        #These will be dictionaries w/ keys 0,1,2,3...
        self.ifd.entryByRow={}
        self.ifd.entryByColumn={}
        self.form = self.vf.getUserInput(self.ifd, modal=0,\
                            blocking=0,scrolledFrame = 1)
        self.form.root.protocol('WM_DELETE_WINDOW',self.close)

        w1=self.ifd.entryByName['cornerLabel']['widget']
        self.ifd.entryByRow[0]=[[w1]]
        self.ifd.entryByColumn[0]=[[w1]]
        w1.oldgridcfg=w1.grid_info()
        self.bgColor=self.ifd.entryByName['cornerLabel']['widget'].cget('bg')


    def guiCallback(self):
        if not self.vf.docked:
            self.vf.warningMsg('Please Read a Docking Log First')
            return
        #FIX THIS
        thisDLO = self.vf.docked.dlo_list[-1]
        if not len(thisDLO.conformations):
            self.vf.warningMsg('Current docking has no docked models ')
            return

        if not hasattr(thisDLO.dpo, 'receptor_stem') or \
                    thisDLO.dpo.receptor_stem=='':
            self.vf.warningMsg('Current docking has no specified macromolecule stem ')
            return

        if not hasattr(self, 'ifd'):
            self.buildForm()
        elif not self.form.root.winfo_ismapped():
            self.form.root.deiconify()
        ifd2=self.ifd2=InputFormDescr(title='Select Dockings')
        entries = self.vf.dockings.keys()
        if len(entries)>1:
            ifd2.append({'name': 'dockings',
                'title':'Docking Logs:',
                'widgetType':'ListChooser',
                'entries': entries,
                'mode': 'multiple',
                'gridcfg':{'row':0,'column':0},
                'lbwcfg':{'height':5, 'width':40}})
            vals = self.vf.getUserInput(ifd2)
            if len(vals)>0 and len(vals['dockings'])>0:
                dockstr=vals['dockings']
                self.doitWrapper(dockstr,log=1, redraw=0)
        else:
            self.doitWrapper(entries,log=1, redraw=0)
        ifd3=self.ifd3=InputFormDescr(title='Docked Conformation Info')
        ifd3.append({'name': 'infoText',
                    'widgetType':Tkinter.Text,
                    'wcfg':{'width':60,
                        'height':10, 
                        'bg':'white'},
                    'gridcfg':{'sticky':Tkinter.W, 'columnspan':3}})
        self.form2 = self.vf.getUserInput(self.ifd3, modal=0,blocking=0)
        self.infoText=self.ifd3.entryByName['infoText']['widget']
        self.form2.root.withdraw()


    def doit(self, dockingList):
        ##self.lastCmdLog = self.logString(dockingList)
        for item in dockingList:
            docking=self.vf.dockings[item]
            if hasattr(docking.dlo_list[0].dpo,'receptor_stem'):
                if not hasattr(docking,'button'):
                    if not hasattr(self, 'ifd'):
                        self.buildForm()
                    self.addButton(docking,'dlg')
                    self.dockingsShown.append(docking)
               

    def addButton(self, docking,labelKey):
        if hasattr(docking, 'button'):
            t = docking.dlo_list[0].filename + ' already read in'
            self.vf.warningMsg(t)
            return
        #FIX THIS
        macro = docking.dlo_list[0].dpo.receptor_stem
        ligand = docking.ligMol.name    
        #DON'T change . to _ (?)
        #9/18: convent: mol in viewer has corrected name, button will have filename
        label=str(docking.ch.conformations[0].binding_energy)
        #label=str(docking.ch.conformations[0].lowEnergy)
        print "checking ", macro , " vs ", self.macros
        if macro in self.macros:
            if ligand in self.ligands:
                rowNum=self.ligands.index(ligand)+1
                columnNum=self.macros.index(macro)+1
                newButton= self.makeIfdEntry(self.ifd,self.bButtons,rowNum,columnNum,label,docking)
                self.updateEntryByRow(newButton)
                self.updateEntryByColumn(newButton)
            else:
                #here add a button under macro for ligand, in a new row labeled ligand
                newButton= self.addRowButton(macro,ligand,label,docking)
        else:
            if ligand in self.ligands:
                #here add a button in ligand row, in a new column labeled macro
                newButton= self.addColumnButton(macro,ligand,label,docking)
            else:
                #here add a newmacrocolumn and ligand row button
                newButton = self.addColumnRowButton(macro,ligand,label,docking)
        docking.button = newButton


    def updateEntryByRow(self, newButton):
        dict=self.ifd.entryByRow
        rnum=int(newButton.row)
        cnum=int(newButton.column)
        if not dict.has_key(rnum):
            dict[rnum]=[]
        while len(dict[rnum])<cnum+1:
            dict[rnum].append([])
        dict[rnum][cnum].append(newButton)


    def updateEntryByColumn(self, newButton):
        dict=self.ifd.entryByColumn
        #FIX THIS
        rnum=int(newButton.row)
        cnum=int(newButton.column)
        if not dict.has_key(cnum):
            dict[cnum]=[]
        while len(dict[cnum])<rnum+1:
            dict[cnum].append([])
        dict[cnum][rnum].append(newButton)


    def addRowButton(self,macro,ligand,label,docking):
        #in this case, column already exists so add new row + button in column
        #set up ifd.entryByRow with maxColumn spaces (?)
        columnNum=self.macros.index(macro)+1
        self.ligands.append(ligand)
        rowNum=self.ligands.index(ligand)+1
        newWidget=self.makeIfdLabelButton(self.ifd, self.bButtons,rowNum,0,ligand)
        self.updateEntryByRow(newWidget)
        self.updateEntryByColumn(newWidget)
        newButton=self.makeIfdEntry(self.ifd,self.bButtons,rowNum,columnNum,label,docking)
        self.updateEntryByRow(newButton)
        self.updateEntryByColumn(newButton)
        self.rowctr=self.rowctr+1
        return newButton
            

    def addColumnButton(self,macro,ligand,label,docking):
        #in this case, add column + button in column
        #set up ifd.entryByColumn with maxRow spaces (?)
        rowNum=self.ligands.index(ligand)+1
        self.macros.append(macro)
        columnNum=self.macros.index(macro)+1
        newWidget=self.makeIfdLabelButton(self.ifd, self.bButtons,0,columnNum,macro)
        self.updateEntryByColumn(newWidget)
        self.updateEntryByRow(newWidget)
        newButton=self.makeIfdEntry(self.ifd,self.bButtons,rowNum,columnNum,label,docking)
        self.updateEntryByColumn(newButton)
        self.updateEntryByRow(newButton)
        self.columnctr=self.columnctr+1
        return newButton
            

    def addColumnRowButton(self,macro,ligand,label,docking):
        #in this case, add new row + new column + button in row-column
        #set up ifd.entryByRow with maxColumn spaces (?)
        #set up ifd.entryByColumn with maxRow spaces (?)
        self.ligands.append(ligand)
        self.macros.append(macro)
        rowNum=self.ligands.index(ligand)+1
        columnNum=self.macros.index(macro)+1
        newWidget=self.makeIfdLabelButton(self.ifd, self.bButtons,0,columnNum,macro)
        self.updateEntryByColumn(newWidget)
        self.updateEntryByRow(newWidget)
        newWidget=self.makeIfdLabelButton(self.ifd, self.bButtons,rowNum,0,ligand)
        self.updateEntryByRow(newWidget)
        self.updateEntryByColumn(newWidget)
        newButton=self.makeIfdEntry(self.ifd,self.bButtons,rowNum,columnNum,label,docking)
        self.rowctr=self.rowctr+1
        self.columnctr=self.columnctr+1
        self.updateEntryByRow(newButton)
        self.updateEntryByColumn(newButton)
        return newButton


    def chooseConformation(self,bestConf,label,event=None):
        #NB bestConf is a name
        currentButton=self.ifd.entryByName[label]['widget']
        mol = self.vf.Mols.NodesFromName(bestConf)[0]
        self.hasShowHide = self.vf.commands.has_key('showMolecules')
##         if len(mol.geomContainer.geoms['lines'].faceSet.faces.array)==0:
        if len(mol.geomContainer.geoms['bonded'].faceSet.faces.array)==0:
            mol.buildBondsByDistance()
            self.vf.displayLines(mol,topCommand=0)
            self.vf.colorByAtomType(mol,['lines'],topCommand=0,redraw=1)
            if self.hasShowHide: 
                self.vf.showMolecules.molList[bestConf].set(1)
                self.vf.showMolecules([bestConf])
        elif self.hasShowHide:
            toggleShowHide(self.vf,bestConf)
        else:
##             mol.geomContainer.geoms['lines'].Set(visible=1)
            mol.geomContainer.geoms['bonded'].Set(visible=1)
            self.vf.GUI.VIEWER.Redraw()

        self.colorButton(currentButton, mol.geomContainer.geoms['master'].visible)


    def showMolecule(self, molFile, event=None):
        #should the user be asked about whether to load molecule?
        self.hasShowHide=self.vf.commands.has_key('showMolecules')
        rnum=rfind(molFile, '.')
        if rnum<0:
            print "illegal name", molFile
            return
        molName=molFile[:rnum]
        currentButton=self.ifd.entryByName[molFile]['widget']
        #molType=molFile[rnum+1:]
        if molName in self.vf.Mols.name:
            mol = self.vf.Mols.NodesFromName(molName)[0]
##             if len(mol.geomContainer.geoms['lines'].faceSet.faces.array)==0:
            if len(mol.geomContainer.geoms['bonded'].faceSet.faces.array)==0:
                mol.buildBondsByDistance()
                self.vf.displayLines(mol,topCommand=0)
                self.vf.colorByAtomType(mol,['lines'],topCommand=0,redraw=1)
                self.vf.showMolecules.molList[mol.name].set(1)
                self.vf.showMolecules([mol.name])
            #check whether it's currently visible
            elif self.hasShowHide:
                self.vf.colorByAtomType(mol,['lines'],topCommand=0,redraw=1)
                toggleShowHide(self.vf,molName)
            else:
##                 mol.geomContainer.geoms['lines'].Set(visible=1)
                mol.geomContainer.geoms['bonded'].Set(visible=1)
                self.vf.GUI.VIEWER.Redraw()
        else:
            #if molType=='pdbqs':
            #    molParser=self.vf.readPDBQS
            #    fileStr='*.pdbqs'
            #elif molType=='pdbq':
            #    molParser=self.vf.readPDBQ
            #    fileStr='*.pdbq'
            #else:
            #    print "error!"
            fileList = glob.glob(molFile)
            if molFile in fileList:
                molList = self.vf.readMolecule(molFile,topCommand=0)
                #molList = apply(molParser,(molFile,),{'log':0,'redraw':1})
                mol = molList[0]
            else:
                #show warning
                msgStr = molFile + " not found in this directory"
                self.vf.warningMsg(msgStr)
                msgName = molFile + ':'
                newFile = self.vf.askFileOpen(types=[('select filename:', fileStr), 
                        ('all files','*')],
                        title = msgName)
                if newFile:
                    molList = self.vf.readMolecule(molFile,topCommand=0)
                    #molList = apply(molParser,(newFile,),{'log':0,'redraw':1})
                    mol = molList[0]
                else: return
            mol.buildBondsByDistance()
            self.vf.displayLines(mol,log=0,redraw=0)
        #if the molecule hooked to button is visible, make bg='white'
        self.colorButton(currentButton, mol.geomContainer.geoms['master'].visible)


    def colorButton(self, button, displayed):
        if displayed:
            button.config(bg='white')
        else:
            button.config(bg=self.bgColor)
            

    def saveLabel(self,row,column,label):
        if self.labels.has_key(row):
            self.labels[row][column]=label
        else:
            self.labels[row]={}
            self.labels[row][column]=label
            

    def makeIfdLabelButton(self,ifd,bButtons,rowNum,columnNum,label):
        shortLabel = os.path.split(label)[-1]
        if not hasattr(ifd, 'form'):
            t= 'Error: form not initialized'
            self.vf.warningMsg(t)
            return
        form = ifd.form
        #this is probably unnecessary(?):how dynamic is this???
        oldRow=None
        for itemList in bButtons:
            for item in itemList: 
                if oldRow==None:
                    oldRow= int(ifd.entryByName[item]['widget'].grid_info()['row'])
                ifd.entryByName[item]['widget'].grid_forget()
        entry= {'name': '%s'%label,
            'widgetType':Tkinter.Button,
            'wcfg':{'width':16, 'relief':'ridge','text':shortLabel,
                'command':CallBackFunction(self.showMolecule,label)}, 
            'gridcfg':{'sticky':Tkinter.E,
                'row':rowNum,'column':columnNum}}
        ifd.entryByName[entry['name']]=entry
        form.addEntry(entry)
        rowctr=2
        newWidget=ifd.entryByName[entry['name']]['widget']
        newWidget.row=int(newWidget.grid_info()['row'])
        newWidget.column=int(newWidget.grid_info()['column'])
        self.saveLabel(newWidget.row, newWidget.column,shortLabel)
        for itemList in bButtons:
            columnctr=0
            for item in itemList:
                #FIX THIS: it isn't general and assumes rows of 1 item
                gridcfg=ifd.entryByName[item]['gridcfg']
                gridcfg['row']= oldRow+1
                gridcfg['column']=columnctr
                columnctr=columnctr+1
                ifd.entryByName[item]['widget'].grid(gridcfg)
            #rowctr=rowctr+1
        form.lift()
        newWidget.oldgridcfg=newWidget.grid_info()
        return newWidget


    def makeIfdEntry(self,ifd,bButtons, rowNum,columnNum,label,docking):
        #this is where the button is made and put in self.ifd2
        if not hasattr(ifd, 'form'):
            t= 'Error: form not initialized'
            self.vf.warningMsg(t)
            return
        form = ifd.form
        #this is probably unnecessary(?):how dynamic is this???
        for itemList in bButtons:
            for item in itemList: 
                ifd.entryByName[item]['widget'].grid_forget()
        #should it be ligMol or macromol
        bestConf=docking.ligMol.name
        #bestConf=docking.ch.conformations[0].members[0].mol.name
        entry= {'name': '%s'%label,
            'widgetType':Tkinter.Button,
            'wcfg':{'width':16,'relief':'ridge',
            'command': CallBackFunction(self.chooseConformation,bestConf,label),    
                'text':label},
            'gridcfg':{'sticky':Tkinter.W+Tkinter.E,
                'row':rowNum,'column':columnNum}} 
        ifd.entryByName[entry['name']]=entry
        form.addEntry(entry)
        newWidget=ifd.entryByName[entry['name']]['widget']
        newWidget.bind('<ButtonPress-2>', CallBackFunction(self.showInfo,docking))
        newWidget.bind('<ButtonRelease-2>', CallBackFunction(self.hideInfo,docking))
        newWidget.row=int(newWidget.grid_info()['row'])
        newWidget.column=int(newWidget.grid_info()['column'])
        self.saveLabel(newWidget.row, newWidget.column,label)
        rowctr=2
        for itemList in bButtons:
            columnctr=0
            for item in itemList:
                #FIX THIS: it isn't general and assumes rows of 1 item
                gridcfg=ifd.entryByName[item]['gridcfg']
                gridcfg['row']=gridcfg['row']+rowctr+2
                gridcfg['column']=columnctr
                columnctr=columnctr+1
                ifd.entryByName[item]['widget'].grid(gridcfg)
            rowctr=rowctr+1
        form.lift()
        self.buttons[docking.dlo_list[0].filename] = ifd.entryByName[entry['name']]['widget']
        #self.buttons[docking.dlgFile] = ifd.entryByName[entry['name']]['widget']
        newWidget.oldgridcfg=newWidget.grid_info()
        return newWidget
        
ADMacroLigandChartGUI=CommandGUI()
ADMacroLigandChartGUI.addMenuCommand('AutoToolsBar', menuText['AnalyzeMB'],
            menuText['showChartMB'], cascadeName=menuText['StatesMB'])
            


class ADDockingChooser(MVCommand):


    def onAddCmdToViewer(self):
        checkHasInitializedDockings(self.vf)
        if self.vf.hasGui:
            if not hasattr(self.vf,'select'):
                self.vf.loadCommand('selectionCommands','select', 'Pmv')
            if not hasattr(self.vf,'displayLines'):
                self.vf.loadCommand('displayCommands', 'displayLines', 'Pmv')
            if not hasattr(self.vf,'showMolecules'):
                self.vf.loadCommand('displayCommands', 'showMolecules', 'Pmv')
        #self.vf.loadModule('displayCommands','Pmv')


    def __init__(self, mode='single', title='Choose a Docked Conformation'):
        MVCommand.__init__(self)
        self.mode=mode
        self.title=title
        

    def guiCallback(self, event=None):
        if not self.vf.docked:
            self.vf.warningMsg('Please Read a Docking Log First')
            return
        d = self.vf.docked
        if hasattr(d.ligMol, 'vina_energy'):
            self.warningMsg("This is a vina result. Please use keyboard arrow keys to view docked results.")
            return "ERROR"
        #d = self.vf.docked.dlo_list[-1]
        if hasattr(d,'chooser'):
            if d.chooser.form.root.winfo_ismapped(): return
            else:
                from Pmv.picker import AtomPicker
                self.chooser.ap=AtomPicker(self.vf,None,0,callbacks=[self.chooser.onPick],immediate=1)
                self.chooser.ap.go(modal=0)
                d.chooser.form.deiconify()
                return

        #FIX THIS
        #if d.hasClusters:
        if hasattr(d, 'clusterer'):
            self.chooser= ClusterDockingChooser(self.vf, self.mode, self.title)
            #self.chooser.ifd.append({'name':'Select Cluster Button',
                                    #'widgetType':Tkinter.Button,
                                    #'text': 'Select',
                                    #'wcfg':{'bd':6},
                                    #'gridcfg':{'sticky':Tkinter.E+Tkinter.W},
                                    #'command': CallBackFunction(self.selectDockingCluster_cb,self.chooser)})                
            self.form = self.chooser.go(modal=0, blocking=0)
            if not self.form:
                return
        #elif d.ch.modelList:
        #    self.chooser= ModelDockingChooser(self.vf, self.mode, self.title)
        #    #self.chooser.ifd.append({'name':'Select Model Button',
        #                            #'widgetType':Tkinter.Button,
        #                            #'text': 'Select',
        #                            #'wcfg':{'bd':6},
        #                            #'gridcfg':{'sticky':Tkinter.E+Tkinter.W},
        #                            #'command': CallBackFunction(self.selectDockingModel_cb,self.chooser)})                
        #    self.form = self.chooser.go(modal=0, blocking=0)
        else:
            self.vf.warningMsg('No clusters found')
            return
        d.chooser = self.chooser
        d.chooser_callback = self.selectDockingCluster_cb
        self.chooser.dlg = d.dlo_list[-1].filename
        ns = os.path.splitext(os.path.basename(self.chooser.dlg))[0]
        tstr = ns + ' Conformation Chooser'
        self.form.root.title(tstr)


    def selectDockingCluster_cb(self,chooser,event=None):
        """called each time the Select Docking button is pressed"""
        ##THIS ISN'T ENOUGH:self.vf.ADanalyze_selectDLG(self.chooser.dlg)
        confs=chooser.getMolSet()
        if confs: self.doitWrapper(confs,log=1,redraw=0)


    def selectDockingModel_cb(self,event=None):
        """called each time the Select Docking button is pressed"""
        ##THIS ISN'T ENOUGH:self.vf.ADanalyze_selectDLG(self.chooser.dlg)
        confs=self.chooser.getMolSet()
        if confs: self.doitWrapper(confs)


    def doit(self, nodes):
        ##self.lastCmdLog = self.logString(nodes)

        if type(nodes)==types.StringType:
            nodes = self.vf.Mols.NodesFromName(nodes)

        if issubclass(nodes.__class__, TreeNode):
            nodes = nodes.setClass([nodes])

        assert isinstance(nodes, MoleculeSet)
        self.vf.setIcomLevel(Molecule)

        for mol in nodes: 
            #check whether mol is visible
            if hasattr(self.vf, 'showMolecules'):
                self.vf.showMolecules([mol.name])
            self.vf.select(mol)
        #need to redo the info string
        #msg = 'Current Docking: '+self.vf.docked.dlgFile
        #self.vf.GUI.pickLabel.configure(text=msg)


ADDockingChooserGUI=CommandGUI()
ADDockingChooserGUI.addMenuCommand('AutoToolsBar', menuText['AnalyzeMB'],
        menuText['chooseConfMB'], cascadeName = menuText['StatesMB'])



class ReadAutoDockStates(MVCommand):

    def onRemoveObjectFromViewer(self, obj):
        #try to clean up a circular reference
        if hasattr(obj, 'docking'):
            d = obj.docking
            delattr(d, 'ligMol')
            delattr(obj, 'docking')


    def guiCallback(self):
        """called each time the 'read states' button is pressed"""
        dpfFile = self.vf.askFileOpen(types=[('select docking parameter filename:', '*.dpf'),
                        ('all files','*')],
                        title = 'Docking Parameter File:')
        if dpfFile:
            trjFile = self.vf.askFileOpen(types=[('select results filename:', '*.res'),
                        ('select trajectory filename:', '*.trj'), 
                        ('all files','*')],
                        title = 'Trajectory or Results File:')
            if trjFile:
                if self.vf.docked:
                    dlo = self.vf.docked.dlo_list[-1]
                    if dlo.filename:
                        t = ' Add to current docking containing ' + dlo.filename + '?'
                    else:
                        t = ' Add to current docking containing ' + \
                                    self.vf.docked.ligMol.name + '?'
                    d = SimpleDialog(self.vf.GUI.ROOT, text = t,
                    buttons = ["No", "Yes", "Cancel"], default = 0, 
                                title = "Add To Previous Docking: ")
                    addToPrevious = d.go()
                    if addToPrevious==2:
                        return 'ERROR'
                    self.doitWrapper(dpfFile, trjFile, addToPrevious, ask=1)
                else:
                    self.doitWrapper(dpfFile, trjFile, 0, ask=1)

        
    def __call__(self, dpfFile, trjFile, addToPrevious=0, ask=1, **kw):
        """None<-ADanalyze_readStates(dpfFile, trjFile)"""
        if not (dpfFile and trjFile): return 'ERROR'
        #print "call: addToPrevious=", addToPrevious, "ask=", ask
        apply(self.doitWrapper, (dpfFile, trjFile, addToPrevious, ask), kw)   


    def getMolecule(self, molFileName, ask=1):
        msgName = molFileName + ':'
        ligFile= self.vf.askFileOpen(types=[('select ligand filename:', '*.pdbq'),
                        ('all files','*')],
                        title = msgName)
        if ligFile:
            lig = self.vf.readPDBQ(ligFile, ask=ask)[0]
        else: return 'ERROR'
        lig.buildBondsByDistance()
        self.vf.displayLines(lig, topCommand=0,redraw=1)
        return lig


    def doit(self, dpfFile, resultsFile, addToPrevious, ask):
        """
        FIX THIS BY USING CORRECT FORM of parser for trjFile
        """ 
        #print "doit: addToPrevious=", addToPrevious, "ask=", ask
        ext = os.path.splitext(resultsFile)[-1]
        if ext=='.dlg':
            if addToPrevious:
                d = self.vf.docked
            else:
                d = Docking()
            d.readDlg(resultsFile)
        elif ext=='.res' or ext=='.results':
            if addToPrevious and isinstance(self.vf.docked.parser,EntropiaParser):
                d = self.vf.docked
            else:
                parser = EntropiaParser()
                d = Docking(parser)
            d.readEntropiaResults(dpfFile, resultsFile)
        elif ext=='.trj':
            if addToPrevious and isinstance(self.vf.docked.parser,TrajParser):
                d = self.vf.docked
            else:
                parser = TrajParser()
                d = Docking(parser)
            d.readEntropiaResults(dpfFile, resultsFile)
        else:
            self.vf.warningMsg('unknown filetype')
            return 'ERROR'
            
        if not d.ligMol:
            #which happens if d.dpo['move']['value'] isn't found
            ligFile= self.vf.askFileOpen(types=[('select ligand filename:', '*.pdbq'), 
                            ('all files','*')],
                            title = 'Select Ligand FILE')
            if ligFile:
                mol = d.ligMol = self.vf.readPDBQ(ligFile,ask=ask,topCommand=0)[0]
            else:
                return 'ERROR'
            #NEED TO make a conformation handler etc etc
            dlo = d.dlo_list[-1]
            d.ch = ConformationHandler(mol, dlo.dpo['about']['value'])
            d.addConformations(parser)
            #add the clusterer if there isn't one at this point
        else:
            # read dlg or MolKit.Read found file
            if d.ligMol not in self.vf.Mols:
                #print "calling addMolecule with ask=", ask
                d.ligMol = self.vf.addMolecule(d.ligMol, ask=ask)

        if d.ligMol=='ERROR': 
            return 'ERROR'

        #FIX THIS: should molecule know its docking???
        d.ligMol.docking = d
        self.vf.dockings[resultsFile] = d
        if not self.vf.docked:
            self.vf.docked = d

        #is this WISE???
        assert len(d.ch.conformations), 'no conformations found '
        if not hasattr(d, 'clusterer'):
            if hasattr(d, 'version') and d.version>=4.0:
                d.clusterer = Clusterer(d.ch.conformations, sort='energy')
            else:
                d.clusterer = Clusterer(d.ch.conformations, sort='docking')
            # all written results are sorted by energy
            d.clusterer.sorted_data = 1
        d.clusterer.usesSubset = 0
        coords = d.ligMol.allAtoms.coords[:]
        d.clusterer.rmsTool = RMSDCalculator(coords)
        d.clusterer.rmsToolRef = '0'
        d.clusterer.usesSubset = 0


ReadAutoDockStatesGUI = CommandGUI()
ReadAutoDockStatesGUI.addMenuCommand('AutoToolsBar', menuText['AnalyzeMB'],
        menuText['readStatesMB'],cascadeName=menuText['StatesMB'])



################ Start of StatesPlayerWidget ###############
 


class StatesPlayerWidget:

    def __init__(self, mol, docking, vf, titleStr=None, sequenceList=None,\
                        idList=None, wname='spw', delta=0, 
                        ask=1, **kw):
        self.mol = mol
        #think about this:
        #assert hasattr(mol, 'docking'), 'mol has no docking'
        #assert hasattr(mol.docking,'ch'), 'mol.docking has no conformation handler'
        #assert hasattr(mol.docking.ch,'conformations'), 'conformation handler has no conformations'
        #check that mol has docking which has ch with conformations???
        self.docking = docking
        self.vf = vf
        for item in ['colorByMolecules','colorByAtomType', 'colorByProperty']:
            if not hasattr(self.vf, item):
                self.vf.loadCommand('colorCommands', item, 'Pmv')
        if not self.vf.colorMaps.has_key('rgb256'):
            mod = __import__("ViewerFramework")
            VFPath = mod.__path__[0]
            self.vf.loadColorMap(os.path.join(VFPath, "ColorMaps","rgb256_map.py"), 
                                    topCommand=0)         
        self.update(sequenceList, idList)

        if self.vf.hasGui:
            self.doTorsionsOnly = Tkinter.IntVar(master=self.vf.GUI.ROOT)
            self.form = self.buildFORM(titleStr)
        else :
            self.doTorsionsOnly = False
        # create a coordinate slot in the molecule for spw use
        mol.allAtoms.addConformation(mol.allAtoms.coords)
        self.coordSlot = len(mol.allAtoms[0]._coords) - 1
        self.ask = ask


    def update(self, sequenceList=None, idList=None):
        #sequenceList is a (possibly) ordered list of conformations
        if not sequenceList:
            #FIX THIS:!!!
            self.sequenceList = self.docking.ch.conformations
            #self.sequenceList = self.mol.docking.ch.conformations
        else:
            self.sequenceList = sequenceList

#        for item in self.sequenceList:
#            item.getCoords()

        #idList is a list of strings to associate with sequenceList
        #eg could be '1-1','1-2','1-3',...,'1-7','2-1', etc
        # or '1,'15', '24', '2' etc etc

        ##indicies in idList are 1 more than indices in sequenceList
        if not idList:
            #insert 0 for original state
            idL = range(0, len(self.sequenceList) + 1)
            self.idList = map(str, idL)
        else:
            #or it could be explicit
            #SHOULD always be incremented by 1 so 0 can be original coords
            self.idList = idList
        #ALSO: close form3 and clear counter
        #these maynot exist yet:
        if hasattr(self, 'form'):
            if hasattr(self.form, 'form3'):
                self.closeform3()
            if hasattr(self.form, 'ent2'):
                self.form.ent2.delete(0,'end')
                self.form.ent2.insert(0, str(self.idList[0]))
                #this calls applyState with reset flag
                self.applyState(-1)


    def custom_validate(self, text):
        if text in self.idList:
            return 1
        else:
            return -1


    def custom_counter(self, text, factor, increment, **kw):
        # text is current content of entry
        # factor is 1 for increment and -1 for decrement
        # increment is value of increment megawidget option

        if not text in self.idList:
            raise ValueError, text + ' not in current idList'
        #check whether ind+factor is in range
        newval = self.idList.index(text) + factor
        if newval<0 or newval>(len(self.idList)-1):
            return text
        else:
            return self.idList[newval]


    def buildFORM(self, titleStr):
        #??FIX THIS:
        mol = self.mol
        at0 = mol.allAtoms[0]
        if not titleStr:
            titleStr = 'Show ' + mol.name + ' Sequence'
        self.stop = 0
        if hasattr(self, 'form'):
            self.form.deiconify()
            return
        maxval = len(self.sequenceList)
        self.rmsVar = Tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.rmsVar.set('rms(ref=0) 0.0000')
        self.energyVar = Tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.colorType = Tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.colorType.set('atom')
        ifd = mol.ifd2 = InputFormDescr(title=titleStr)
        ifd.append({'name':'energyLab',
            'widgetType': Tkinter.Label,
            'textvariable': self.energyVar,
            'wcfg':{'bd':4},
            'gridcfg':{'sticky':'nesw', 'columnspan':2}})
            #'gridcfg':{'sticky':'nesw','row':-1, 'column':1}}),
        ifd.append({'name':'rmsLab',
            'widgetType': Tkinter.Label,
            'textvariable': self.rmsVar,
            'wcfg':{'bd':4},
            'gridcfg':{'sticky':'nesw'}}),
            #'gridcfg':{'sticky':'nesw','row':-1, 'column':1}}),
        colorTypeList = ['atom','molecule', 'no change']
        #colorTypeList = ['atom','vdw','elec_stat','total','molecule']
        ifd.append({'widgetType':Pmw.ComboBox,
                            'name':'colorType',
                            #'defaultValue':self.colorType.get(),
                            'wcfg':{'label_text':'color by',
		 		                    'entryfield_value':self.colorType.get(),
                                    'labelpos':'w',
                                    'listheight':'80',
                                    'scrolledlist_items': colorTypeList,
				                    'selectioncommand': self.updateColor,
                                    },
                            #'gridcfg':{'sticky':'w'}})
                            'gridcfg':{'sticky':'nesw','row':-1, 'column':1}}),
        ifd.append({'widgetType':Pmw.Counter,
			    'name':'statesCounter',
			    'required':1,
			    'wcfg':{#'labelpos': 'n,
		 		    #'label_text':'conformation:  ',
                    'autorepeat':0,
		 		    'entryfield_value':self.idList[0],
                    'datatype': self.custom_counter,
		 		    'entry_width':9,
		 		    'entryfield_validate': self.custom_validate },
		 	    'gridcfg':{'sticky':'nesw', 'columnspan':2}})
        ifd.append({'name':'selectCB',
            'widgetType': Tkinter.Checkbutton,
            'text':'Show IdList',
            'wcfg':{'bd':4},
            'gridcfg':{'sticky':'nesw', 'columnspan':1},
            #'gridcfg':{'sticky':'nesw', 'columnspan':2},
            'command': self.showStatesList})
        ifd.append({'name': 'doTransCB',
            'widgetType': Tkinter.Checkbutton,
            'text':'torsions only',
            #'text':'transform torsions only',
            'variable': self.doTorsionsOnly,
            'wcfg':{'bd':4},
            'gridcfg':{'sticky':'nesw', 'row':-1,'column':1}})
            #'gridcfg':{'sticky':'nesw', 'columnspan':2}})
        ifd.append({'name': 'playB',
            'widgetType': Tkinter.Button,
            'text':'Play Sequence',
            'wcfg':{'bd':4},
            'gridcfg':{'sticky':'nesw','columnspan':1},
            'command':self.Play_cb})
        ifd.append({'name': 'playReverseB',
            'widgetType': Tkinter.Button,
            'text':'Play it Reverse',
            'wcfg':{'bd':4},
            'gridcfg':{'sticky':'nesw','row':-1, 'column':1},
            'command':self.PlayRev_cb})
        ifd.append({'name': 'stopB',
            'widgetType': Tkinter.Button,
            'text':'Stop',
            'wcfg':{'bd':4},
            'gridcfg':{'sticky':'nesw'},
            'command':self.Stop_cb})
        ifd.append({'name': 'pauseB',
            'widgetType': Tkinter.Button,
            'text':'Pause',
            'wcfg':{'bd':4},
            'gridcfg':{'sticky':'nesw','row':-1, 'column':1},
            'command':self.Pause_cb})
        ifd.append({'name': 'closeB',
            'widgetType': Tkinter.Button,
            'text':'Make rms refcoords',
            'wcfg':{'bd':4},
            'gridcfg':{'sticky':'nesw'},
            'command':self.MakeRef_cb})
        ifd.append({'name': 'buildB',
            'widgetType': Tkinter.Button,
            'text':'Build',
            'wcfg':{'bd':4},
            'gridcfg':{'sticky':'nesw','row':-1, 'column':1},
            'command':self.Build_cb})
        ifd.append({'name': 'closeB',
            'widgetType': Tkinter.Button,
            'text':'Close',
            'wcfg':{'bd':4},
            'gridcfg':{'sticky':'nesw', 'columnspan':2},
            'command':self.Close_cb})
        form = self.vf.getUserInput(ifd, modal=0,blocking=0)
        form.root.protocol('WM_DELETE_WINDOW',CallBackFunction(self.Close_cb,mol))

        form.ifd = ifd
        ctr = ifd.entryByName['statesCounter']['widget']
        entF = ctr.component('entryfield')
        form.ent2 = entF._entryFieldEntry
        da = ctr.component('downarrow')
        ua = ctr.component('uparrow')
        for item in [da,ua]:
            item.bind('<ButtonPress-1>', self.SetState_cb, '+')
        form.ent2.bind('<Return>', self.SetState_cb, '+')
        form.counter = form.ifd.entryByName['statesCounter']['widget']
        form.showList = form.ifd.entryByName['selectCB']['widget']
        form.showVar = form.ifd.entryByName['selectCB']['variable']
        ctr2 = ifd.entryByName['colorType']['widget']
        entF2 = ctr2.component('entryfield')
        #all this to get a handle to the Tkinter.Entry buried in the pww
        form.ent3 = entF2._entryFieldEntry
        form.ent3.bind('<Return>', self.updateColor, '+')
        form.ent3.config(width=8)
        return form


    def MakeRef_cb(self, event=None):
        """None<-MakeRef_cb(mol, event=None)

        makes current molecule coordinates rmstool refCoords
        """
        clusterer = self.vf.docked.clusterer
        rmsTool = clusterer.rmsTool
        if clusterer.usesSubset:
            rmsTool.setRefCoords(clusterer.subset.coords[:])
        else:
            rmsTool.setRefCoords(self.mol.allAtoms.coords[:])
        clusterer.rmsToolRef = self.form.ent2.get()


    def updateColor(self, event=None):
        """None<-updateColor(mol, event=None)

        makes current molecule coordinates rmstool refCoords
        """
        confNum = self.form.ent2.get()
        colorType = self.form.ent3.get()
        #this is already taken care of before
        at0 = self.mol.allAtoms[0]
        has_elec = hasattr(at0, 'estat_energy')
        has_vdw = hasattr(at0, 'vdw_energy')
        has_tot = hasattr(at0, 'total_energy')
        if colorType==None or colorType=='no change':
            return
        colormap = self.vf.colorMaps['rgb256']
        if confNum=='0': 
            if colorType=='molecule':
                self.vf.colorByMolecules(self.mol.allAtoms, ('lines',), topCommand=0)
            else:
                self.vf.colorByAtomType(self.mol.allAtoms, ('lines',), topCommand=0)
        elif colorType=='atom':
            self.vf.colorByAtomType(self.mol.allAtoms, ('lines',), topCommand=0)
        elif colorType=='molecule':
            self.vf.colorByMolecules(self.mol.allAtoms, ('lines',), topCommand=0)
        elif colorType=='elec_stat' and has_elec:
            val_list = Numeric.array(self.mol.allAtoms.estat_energy)
            mini = min(val_list)
            maxi = max(val_list)
            self.vf.colorByProperty(self.mol.allAtoms, ('lines',),
                'estat_energy',colormap='rgb256',mini=mini, maxi=maxi, 
                propertyLevel='Atom', topCommand=0)
        elif colorType=='vdw' and has_vdw:
            val_list = Numeric.array(self.mol.allAtoms.vdw_energy)
            mini = min(val_list)
            maxi = max(val_list)
            self.vf.colorByProperty(self.mol.allAtoms, ('lines',),
                'vdw_energy', colormap='rgb256',mini=mini, maxi=maxi, propertyLevel='Atom',topCommand=0)
        elif colorType=='total' and has_tot:
            val_list = Numeric.array(self.mol.allAtoms.total_energy)
            mini = min(val_list)
            maxi = max(val_list)
            self.vf.colorByProperty(self.mol.allAtoms, ('lines',),
                'total_energy',colormap='rgb256',mini=mini, maxi=maxi,propertyLevel='Atom',topCommand=0)


    def Build_cb(self, event=None):
        """None<-Build_cb(mol, event=None)

        builds new molecule with current coordinates and adds it to the viewer
        """
        #FIRST CHECK THAT THIS HASN'T already been built
        #get the current counter content for name of new molecule
        #w = self.form.ifd.entryByName['statesCounter']['widget']
        numStr = self.form.counter.get()
        #numStr = w.get()
        newname = self.mol.name + '_conf_' + numStr
        if newname in self.vf.Mols.name:
            msg = newname + ' already in viewer. Not building a second copy'
            self.vf.warningMsg(msg)
            return 'ERROR'

        allLines = self.mol.parser.allLines
        newLines = []
        #coords = self.mol.chains.residues.atoms.coords
        coords = self.mol.allAtoms.coords
        c = self.docking.ch.current_conf
        estat_energies = []
        hasESTAT = 0
        vdw_energies = []
        hasVDW = 0
        if hasattr(c, 'estat_energies'):
            estat_energies = c.estat_energies
        if len(estat_energies):
            hasESTAT = 1
        if hasattr(c, 'vdw_energies'):
            vdw_energies = c.vdw_energies
        if len(vdw_energies):
            hasVDW = 1
        #does this matter here?
        #coords = mol.allAtoms.coords
        ctr = 0
        #foundCONECT = 0
        for l in allLines:
            if find(l, 'ATOM')==0 or find(l, 'HETA')==0:
                cc = coords[ctr]
                restLine = l[54:]
                if hasESTAT and hasVDW:
                    #occupancy is rec[54:60]
                    #temperatureFactor is rec[60:66]
                    #charge is rec[70:76]
                    restLine = '%6.3f%6.3f'%(vdw_energies[ctr],estat_energies[ctr]) 
                else:
                    restLine = l[54:66]
                newLines.append(l[:30] +'%8.3f%8.3f%8.3f'%(cc[0],cc[1],cc[2])+ \
                                    restLine + l[66:])
                ctr = ctr+1
            #this code builds new CONECT records but doesn't appear
            # to be necessary
            #elif find(l,'CONECT')==0:
                #if not foundCONECT:
                ##have to redo the CONECT records here
                    #foundCONECT = 1
                    #for at in self.mol.allAtoms:
                        #l = 'CONECT%5d'%at.number 
                        #for b in at.bonds:
                            #at2 = b.atom1
                            #if at2==at:
                                #at2 = b.atom2
                            #l = l + '%5d'%at2.number
                        #newLines.append(l)
                        ##print 'for ', at.full_name(),':',l
            else:
                newLines.append(l)

        pdbqParser = PdbqParser() 
        pdbqParser.allLines = newLines
        filename = pdbqParser.filename = self.mol.parser.filename + '_conf_' + numStr
        newMol = pdbqParser.parse()[0]          
        newMol.name = newname
        #newMol.name = mol.name + '_conf_' + numStr
        newMol = self.vf.addMolecule(newMol, ask=self.ask)
        newMol.colored = 0
        if hasESTAT:
            newMol.allAtoms.estat_energy = estat_energies[:]
            newMol.allAtoms.temperatureFactor = estat_energies[:]
        if hasVDW:
            newMol.allAtoms.vdw_energy = vdw_energies[:]
            newMol.allAtoms.occupancy = vdw_energies[:]
        if hasVDW and hasESTAT:
            total_energy = Numeric.add(newMol.allAtoms.estat_energy, newMol.allAtoms.vdw_energy)
            newMol.allAtoms.total_energy = total_energy.tolist()
        #FIX THIS: this function doesn't log itself...


    def Close_cb(self,  event=None):
        """None<-Close_cb( event=None)
           withdraws self.form
        """
        self.form.withdraw()


    def SetState_cb(self, event = None):
        """None<-SetState_cb(mol, event=None)
           calls applyState to mol with contents of statesCounter
        """
        idStr = self.form.counter.get()
        #val = int(self.form.counter.get())
        #NB: -1 because of offset problem
        confInd = self.idList.index(idStr) - 1
        #always call applyState with a zero-based integer index
        self.applyState(confInd)
        #self.applyState(val)


    def Play_cb(self,  event = None):
        cam = self.vf.GUI.VIEWER.cameras[0]
        self.stop = 0
        #this command plays the states in self.sequenceList 
        # via its matching idList
        curIdStr = self.form.counter.get()
        curConfInd = self.idList.index(curIdStr)
        for id in self.idList[curConfInd:]:
            if self.stop: break
            #this is for display in counter, so 1-based
            self.form.ent2.delete(0,'end')
            self.form.ent2.insert(0, str(id))
            #NB: -1 because of offset problem
            newInd = self.idList.index(id) - 1
            self.applyState(newInd)
            cam.update()
            

    def PlayRev_cb(self, event = None):
        cam = self.vf.GUI.VIEWER.cameras[0]
        self.stop = 0
        #FIX THIS
        #HAVE TO START AT CURRENT VALUE AND GO DOWN
        curIdStr = self.form.counter.get()
        #NB: reserve 0 as reset to original cooords
        if curIdStr==0:
            return
        #this assumes that 0 may be in list and played as '1'
        curConfInd = self.idList.index(curIdStr)
        revList =[]
        for v in self.idList[:curConfInd]:
            revList.append(v)
        revList.reverse()
        for id in revList:
            if self.stop: break
            #this is for display in counter, so 1-based
            self.form.ent2.delete(0,'end')
            self.form.ent2.insert(0, id)
            #NB: -1 because of offset problem
            newInd = self.idList.index(id) - 1
            #FIX THIS ????:
            #self.form.counter.decrement()
            self.applyState(newInd)
            cam.update()


    def Stop_cb(self, event=None):
        self.stop = 1
        self.applyState(-1)
        self.form.ent2.delete(0,'end')
        self.form.ent2.insert(0,'0')
        self.form.ent3.delete(0, 'end')
        self.form.ent3.insert(0,'atom')
        self.updateColor()
            

    def Pause_cb(self, event=None):
        self.stop = 1


    def applyState(self, confInd):
        """None<-applyState(mol, confInd)"""
        mol = self.mol
        ch = self.docking.ch
        clust = self.docking.clusterer
        if clust.usesSubset:
            allAts = clust.subset
        else:
            allAts = mol.allAtoms
        
        # -1 is key for go back to original
        if int(confInd)==-1:
            mol.allAtoms.setConformation(0)
            conf = None
        else:
            #in this case want to get conformation
            conf = self.sequenceList[confInd]
            if self.vf.hasGui:
                doTorsion = self.doTorsionsOnly.get()
            else :
                doTorsion = self.doTorsionsOnly
            if doTorsion :
                coords = conf.getTorsionOnlyCoords()
                mol.allAtoms.updateCoords(coords[:], self.coordSlot)
            else:
                #use conformation handler so if there is extra
                # stuff like estat_energies or vdw_energies it
                # will get updated also
                ch.set_conformation(conf, self.coordSlot)
        
                
        #NB: ch.current_conf could be None if reset (-1)
        #ch.current_conf = conf
        # 9/10: changed to set allAts above
        #t = 'rms(ref='+ self.rmsToolRef +') %8.4f'%(round(self.rmsTool.computeRMSD(mol.allAtoms.coords[:]),3))
        if hasattr(clust, 'rmsToolRef'):
            t = 'rms(ref='+ clust.rmsToolRef +') %8.4f'%(round(clust.rmsTool.computeRMSD(allAts.coords[:]),3))
        else:
            t = 'no rms available'
        if conf:
            s = 'Energies: binding = %-5.5s : docking = %-5.5s' %(conf.binding_energy, conf.docking_energy)
        else:
            s = 37*' '
        #if not self.vf.hasGui: return
        event = EditAtomsEvent('coords', mol.allAtoms)
        self.vf.dispatchEvent(event)
        #modEvent = ModificationEvent('edit','coords', mol.allAtoms)
        #mol.geomContainer.updateGeoms(modEvent)
        #self.vf.measureDistanceGC.update()
        #self.vf.measureAngleGC.update()
        #self.vf.measureTorsionGC.update()
        if self.vf.hasGui: 
            self.rmsVar.set(t)
            self.energyVar.set(s)
            self.updateColor()
            self.vf.GUI.VIEWER.Redraw()
        #NB: in ONE test: modEvent was slightly slower (.40 vs .34)
        #self.vf.displayLines(mol, negate=1, topCommand=0)
        #self.vf.displayLines(mol, topCommand=0, redraw=1)

    def showStatesList(self, event=None):
        ##THIS WON'T WORK WITHOUT MORE INFO:
        ##IF just pass list of conformations, don't know their indices
        if self.form.showVar.get()=='0':
            self.form.showList.config(text = 'Show Conformation List')
            self.closeform3()
            return
        else:
            self.form.showList.config(text = 'Hide Conformation List')
        #FIXT THIS@@
        #statesToList = self.idList
        #if not len(statesToList): 
            #self.vf.warningMsg('no bin currently selected')
            #return
        #sTL_1based = map(lambda x: x+1, statesToList)
        mol = self.mol
        entries = self.idList
        #entries = []
        #for i in self.idList:
            ##FIX THIS @@
            ##make the conformation labels 1 based to go with counter
            #entries.append(int(i) + 1)
        if hasattr(self.form, 'form3'): return
        ifd2 = InputFormDescr(title='Choose Conformation')
        ifd2.append({'widgetType':'ListChooser',
            'name':'stateObjs',
            'entries':entries,
            'wcfg':{'title':'Pick state',
                    'mode':'single'},
            'lbwcfg':{'height':10},
            'gridcfg':{'sticky':'nsew', 'column':100,
                     'rowspan':10}})
        self.form.form3 = self.vf.getUserInput(ifd2, modal=0, blocking=0)
        self.form.form3.root.protocol('WM_DELETE_WINDOW',self.form.form3.root.withdraw)

        self.form.lb = ifd2.entryByName['stateObjs']['widget'].lb
        ifd2.entryByName['stateObjs']['widget'].title.config(text=
            'Double click to set state')
        self.form.lb.bind('<Double-Button-1>', self.showSpecConf)


    def showSpecConf(self, event=None):
        lb = self.form.lb
        if len(lb.curselection()):
            #state = lb.get(lb.curselection()[0])
            idStr = lb.get(lb.curselection()[0])
            #NB: -1 because of offset problem
            confInd = self.idList.index(idStr) -1
            #state = int(lb.get(lb.curselection()[0]))
            #FIX THIS??? applyState takes a 1-based number
            #ind = self.idList.index(int(state)) + 1
            self.applyState(confInd)
            self.form.ent2.delete(0,'end')
            self.form.ent2.insert(0, idStr)
            self.vf.GUI.VIEWER.cameras[0].update()


    def closeform3(self, event=None):
        if hasattr(self.form, 'form3'):
            self.form.form3.destroy()
            delattr(self.form, 'form3')

######################################END OF WIDGET###################


class ShowAutoDockStatesBaseCmd(MVCommand):
    """base cmd for cmds allowing user to show different docked states of a ligand"""


    def onAddCmdToViewer(self):
        checkHasInitializedDockings(self.vf)
        from AutoDockTools.autoanalyzeCommands import StatesPlayerWidget
        self.mode = 'single'
        self.title = 'Choose Molecule'
        self.molList = []
        self.dockings = []
        for item in ['colorByProperty','colorByAtomType']:
            if not hasattr(self.vf, item):
                self.vf.loadCommand('colorCommands', item, 'Pmv')
        if self.vf.hasGui:
            self.doTorsionsOnly = Tkinter.IntVar(master=self.vf.GUI.ROOT)
            #self.vf.loadModule('measureCommands')


    def onRemoveObjectFromViewer(self, obj):
        if hasattr(obj, 'ifd'):
            delattr(obj, 'ifd')
        if hasattr(obj, 'form'):
            obj.form.destroy()
        if hasattr(obj, 'form2'):
            obj.form2.destroy()
        if hasattr(obj, 'spw'):
            if hasattr(obj.spw, 'sequenceList'):
                for c in obj.spw.sequenceList:
                    if hasattr(c, 'mol'):
                        delattr(c, 'mol')
            if hasattr(obj.spw, 'mol'):
                delattr(obj.spw, 'mol')
        if obj in self.molList:
            ind = self.molList.index(obj)
            self.molList.remove(obj)
            d = self.dockings[ind]
            self.dockings.remove(d)
        if hasattr(obj, 'torTree'):
            for key in obj.torTree.rootNode.__dict__.keys():
                obj.torTree.rootNode.__dict__[key] = None
            for item in obj.torTree.torsionMap:
                for key in item.__dict__.keys():
                    item.__dict__[key] = None


            
    def chooseMolecule_cb(self, event = None):
        """
        invoked from guiCallback"""
        mol = self.chooser.getMolSet()
        if mol in self.molList:
            self.chooser.form.withdraw()
            self.doitWrapper(mol, log=1, redraw=0)
            #else:
                #if not hasattr(mol, 'docking'):
                    #msg = mol.name + ' has no docking'
                    #self.vf.warningMsg(msg)
                    #return 'ERROR'
        else:
            msg = mol.name + ' is not connected to a docking'
            self.vf.warningMsg(msg)
            return "ERROR"



    def __call__(self, mol, log=1, **kw):
        mols = self.vf.expandNodes(mol)
        if len(mols)==0:
            return 'ERROR'
        #need to construct a list of the current molecules connected to dockings
        self.buildMolList()
        kw['log'] = log
        #kw.setdefault('log', 1)
        apply(self.doitWrapper, (mols[0],), kw)


    def doit(self, mol, **kw):
        pass


    def buildMolList(self):
        for k, d in self.vf.dockings.items():
            self.molList.append(d.ligMol)
            self.dockings.append(d)


    def guiCallback(self):
        """called each time the 'show states' button is pressed"""
        if not len(self.vf.Mols):
            msg = 'no molecules in viewer'
            self.vf.warningMsg(msg)
            return 'ERROR'
        self.dockings = []
        molList = self.molList = []
        #FIX THIS!!!
        self.buildMolList()
        #for k, d in self.vf.dockings.items():
            #self.molList.append(d.ligMol)
            #self.dockings.append(d)
        if not len(molList):
            msg = 'no molecules with dockings in viewer'
            self.vf.warningMsg(msg)
            return 'ERROR'
        elif len(molList)==1:
            self.doitWrapper(molList[0], log=1, redraw=0)
        else:
            self.chooser = MoleculeChooser(self.vf, self.mode, self.title)
            self.chooser.ipf.append({'name':'Select Button',
                                     'widgetType':Tkinter.Button,
                                     'text':'Select Molecule',
                                     'wcfg':{'bd':6},
                                     'gridcfg':{'sticky':Tkinter.E+Tkinter.W},
                                     'command': self.chooseMolecule_cb})
            self.form = self.chooser.go(modal=0, blocking=0)
            lb = self.chooser.ipf.entryByName['Molecule']['widget'].lb
            lb.bind("<Double-Button-1>",self.chooseMolecule_cb)


#################################END OF SHOWSTATESBASECMD ###################



class ShowAutoDockStates(ShowAutoDockStatesBaseCmd):
    """allows user to show different docked states of a ligand"""


    def doit(self, mol, ask=1):
        #by not specifying sequenceList and idList
        #spw is built with mol.docking.ch.conformations and a range (see class).
        mols = self.vf.expandNodes(mol)
        if len(mols)==0:
            return 'ERROR'
        mol = mols[0]
        if hasattr(mol, 'vina_energy'):
            self.warningMsg("This is vina result. Please use keyboard arrow keys to view docked results")
            return 'ERROR'
        #need to build a molList
        index = self.molList.index(mol)
        docking = self.dockings[index]
        #mol = docking.ligMol
        mol.stop = 0
        if hasattr(mol, 'spw'):
            #10/26/04:
            confs = docking.ch.conformations
            mol.spw.update(confs)
            mol.spw.nextFrame(1)
            #mol.spw.update()
            if not mol.spw.form.f.winfo_ismapped():
                mol.spw.form.deiconify()
        else:
            #TRACE whether ask is necessary!!!
            if self.vf.userpref['Player GUI']['value']=='ConformationPlayer':
                mol.spw = ConformationPlayer(mol, docking, self.vf, 
                            mol.name, docking.ch.conformations, form2=1,
                            ask=ask)
                mol.spw.nextFrame(1)
            else:
                mol.spw = StatesPlayerWidget(mol, docking, self.vf, mol.name, ask=ask)


ShowAutoDockStatesGUI = CommandGUI()
ShowAutoDockStatesGUI.addMenuCommand('AutoToolsBar', 
    menuText['AnalyzeMB'], menuText['showStatesMB'],
        cascadeName=menuText['StatesMB'])


#################################END OF SHOWSTATES###################

class ShowAutoDockStatesByEnergy(ShowAutoDockStatesBaseCmd):
    """allows user to show docked states of a ligand ordered by energy"""


    def doit(self, mol, ask=1):
        #by not specifying sequenceList and idList
        #spw is built with mol.docking.ch.conformations and a range (see class).
        mols = self.vf.expandNodes(mol)
        if len(mols)==0:
            return 'ERROR'
        mol = mols[0]
        if hasattr(mol, 'vina_energy'):
            self.warningMsg("This is a vina result. Please use keyboard arrow keys to view docked results ranked by energy")
            return 'ERROR'
        #need to build a molList
        index = self.molList.index(mol)
        docking = self.dockings[index]
        #mol = docking.ligMol
        mol.stop = 0
        if not hasattr(docking, 'clusterer'):
            docking.clusterer = Clusterer(docking.ch.conformations, sort='energy')
        clist = []
        for i in docking.clusterer.argsort:
            clist.append(docking.clusterer.data[int(i)])
            
        if hasattr(mol, 'spw'):
            mol.spw.update(clist)
            mol.spw.nextFrame(1)
            if not mol.spw.form.f.winfo_ismapped():
                mol.spw.form.deiconify()
        else:
            #TRACE whether ask is necessary!!!
            mol.spw = ConformationPlayer(mol, docking, self.vf, 
                            mol.name, clist, form2=1,
                            ask=ask)
							# buttonMask={'recordB':False})
            mol.spw.nextFrame(1)


ShowAutoDockStatesByEnergyGUI = CommandGUI()
ShowAutoDockStatesByEnergyGUI.addMenuCommand('AutoToolsBar', 
    menuText['AnalyzeMB'], menuText['showStatesByEnergyMB'],
        cascadeName=menuText['StatesMB'])


#################################END OF SHOWSTATES BY ENERGY###################


class ShowAutoDockPopulation(ShowAutoDockStatesBaseCmd):
    """allows user to show population created in autodock """


    def doit(self, mol, ask=1):
        #by not specifying sequenceList and idList
        #spw is built with mol.docking.ch.conformations and a range (see class).
        mols = self.vf.expandNodes(mol)
        if len(mols)==0:
            return 'ERROR'
        mol = mols[0]
        if hasattr(mol, 'vina_energy'):
            self.warningMsg("Viewing initial population not available for vina results")
            return 'ERROR'
        #need to build a molList
        index = self.molList.index(mol)
        docking = self.dockings[index]
        #mol = docking.ligMol
        mol.stop = 0
        if docking.ph:
            if hasattr(mol, 'ppw'):
                #10/26/04:
                confs = docking.ph.individuals
                mol.ppw.update(confs)
                #mol.ppw.update()
                if not mol.ppw.form.f.winfo_ismapped():
                    mol.ppw.form.deiconify()
            else:
                #TRACE whether ask is necessary!!!
                idList = range(len(docking.ph.all_populations))
                titleStr = mol.name + " population player"
                mol.ppw = PopulationPlayer(mol, docking, self.vf, 
                                titleStr, docking.ph.individuals, form2=1,
                                ask=ask, idList=idList)
        else:
            msg = "It is not possible to visualize the initial population for the current docking. Only AutoDock4 results computed with 'outlev' set to 4 in the dpf include the necessary information about the initial population"
            self.warningMsg(msg)
            return "ERROR"


ShowAutoDockPopulationGUI = CommandGUI()
ShowAutoDockPopulationGUI.addMenuCommand('AutoToolsBar', 
    menuText['AnalyzeMB'], menuText['showPopulationMB'],
        cascadeName=menuText['StatesMB'])



#################################END OF SHOWSTATES###################



class ShowAutoDockStatesHISTOGRAM(ShowAutoDockStatesBaseCmd):
    """allows user to show histogram of docked states of a ligand"""


    def onRemoveObjectFromViewer(self, obj):
        if hasattr(obj, 'ifd'):
            delattr(obj, 'ifd')
        if hasattr(obj, 'form'):
            obj.form.destroy()
        if hasattr(obj, 'form2'):
            obj.form2.destroy()
        if hasattr(obj, 'histNB'):
            obj.histNB.master.destroy()
        for item in ['ehist','elist','min','max','lb','ent','ent2',\
                     'histform', 'histNB','ifd','ifd2','nbins']:
            if hasattr(obj, item):
                delattr(obj, item)
        if obj in self.molList:
            ind = self.molList.index(obj)
            self.molList.remove(obj)
            d = self.dockings[ind]
            self.dockings.remove(d)
        if hasattr(obj, 'spw'):
            if hasattr(obj.spw, 'sequenceList'):
                for c in obj.spw.sequenceList:
                    if hasattr(c, 'mol'):
                        delattr(c, 'mol')
            if hasattr(obj.spw, 'mol'):
                delattr(obj.spw, 'mol')
        if hasattr(obj, 'torTree'):
            for key in obj.torTree.rootNode.__dict__.keys():
                obj.torTree.rootNode.__dict__[key] = None
            for item in obj.torTree.torsionMap:
                for key in item.__dict__.keys():
                    item.__dict__[key] = None


    def buildInteractiveGraph(self, mol, docking):
        if hasattr(mol, 'histNB'):
            mol.histNB.exit()
        r = (float(mol.min.get()), float(mol.max.get()))
        mol.ehist = HistogramRI(mol.elist,mol.nbins.get(),range=r)
        mol.ehist.createReverseIndex()
        nodeList = mol.ehist.array
        tstr = mol.name + ' histogram'
        top = Tkinter.Toplevel(master=self.vf.GUI.ROOT)
        top.title(tstr)
        mol.histNB = interactiveHistogramGraph.InteractiveHistogramGraph(mol.name, 
                master=top, nodeList = nodeList,
                reverseIndex=mol.ehist.reverseIndex)
        mol.histNB.draw.bind('<Button-1>', CallBackFunction(self.showHist,mol,docking),'+')
        mol.histform.withdraw()


    def showHist(self, mol, docking, event=None):
        #histogram
        titleStr = 'Show ' + mol.name + ' Histogram Sequence'
        #THIS IS WHERE spw is built or deiconified
        #NB: currentInfo is 0-based
        #this is a list [0,3,8,9]
        #currentInfo is a list of indices
        cI = mol.histNB.currentInfo
        cI0 = map(lambda x:x+1, cI)
        idList = map(str, cI0)
        idList.insert(0,'0')
        clist = []
        confs = docking.ch.conformations
        #clist is a constructed list of conformations 
        #corresponding to indices in cI
        for n in cI:
            #CHECK THIS are these 0-based or 1-based
            clist.append(confs[n])

        if hasattr(mol, 'spw'):
            mol.spw.update(clist, idList)
            mol.spw.nextFrame(1)
        else:
            #the first time
            if self.vf.userpref['Player GUI']['value']=='StatesPlayerWidget':
                mol.spw = StatesPlayerWidget(mol, docking, self.vf, mol.name, 
                                            clist, idList, ask=self.ask)
            else:
                titleStr = mol.name + " player"
                mol.spw = ConformationPlayer(mol, docking, self.vf, titleStr,
                                            clist, idList=idList, form2=1, ask=self.ask)
                mol.spw.nextFrame(1)
        # make sure widget is visible: de-iconify if iconified.
        if not mol.spw.form.f.winfo_ismapped():
            mol.spw.form.deiconify()


    def doit(self, mol, ask=1):
        mols = self.vf.expandNodes(mol)
        if len(mols)==0:
            return 'ERROR'
        mol = mols[0]
        #need this for testing purposes (sigh)
        self.ask = ask
        ind = self.molList.index(mol)
        #mol = docking.ligMol
        docking = self.dockings[ind]
        mol.stop2 = 0
        maxval = len(docking.ch.conformations)
        if not hasattr(mol, 'elist'):
            s = docking.ch.conformations
            elist = []
            for c in s: 
                if c.docking_energy:
                    #print "using 'docking_energy'"
                    elist.append(c.docking_energy)
                elif c.binding_energy: 
                    #print "using 'binding_energy'"
                    elist.append(c.binding_energy)
                elif c.energy: 
                    #print "using 'energy'"
                    elist.append(c.energy)
            if not len(elist): 
                self.vf.warningMsg('No energies available for docking.ch.conformations')
                return 'ERROR'
            mol.elist = Numeric.array(elist)
            mol.r = [Numeric.minimum.reduce(mol.elist), 
                        Numeric.maximum.reduce(mol.elist)]
            mol.nbins = Tkinter.IntVar(master=self.vf.GUI.ROOT)
            mol.nbins.set(10)
            mol.min = Tkinter.StringVar(master=self.vf.GUI.ROOT)
            mol.min.set(str(mol.r[0]))
            mol.max = Tkinter.StringVar(master=self.vf.GUI.ROOT)
            mol.max.set(str(mol.r[1]))

        titleStr = 'Show ' + mol.name + ' Conformations'
        ifd =  InputFormDescr(title=titleStr)
        ifd.append({'name': 'nbinEnt',
            'widgetType': Tkinter.Entry,
            'wcfg':{'bd':4,
                'label': 'Number of Bins:',
                'textvariable': mol.nbins
                },
            'gridcfg':{'sticky':'nesw','columnspan':2} })
        ifd.append({'name': 'minEnt',
            'widgetType': Tkinter.Entry,
            'wcfg':{'bd':4,
                'label': 'Energy Minimum:',
                'textvariable': mol.min
                },
            'gridcfg':{'sticky':'nesw','columnspan':2} })
        ifd.append({'name': 'maxEnt',
            'widgetType': Tkinter.Entry,
            'wcfg':{'bd':4,
                'label': 'Energy Maximum:',
                'textvariable': mol.max
                },
            'gridcfg':{'sticky':'nesw','columnspan':2} })
        ifd.append({'name': 'buildBut',
            'widgetType': Tkinter.Button,
            'text':'Build Histogram:',
            'wcfg':{'bd':4},
            'command': CallBackFunction(self.buildInteractiveGraph, mol, docking),
            'gridcfg':{'sticky':'nesw','columnspan':2} })
        ifd.append({'name': 'resetRangeBut',
            'widgetType': Tkinter.Button,
            'text':'Reset Range',
            'wcfg':{'bd':4},
            'command': CallBackFunction(self.resetRange, mol),
            'gridcfg':{'sticky':'nesw'} })
        ifd.append({'name': 'closeB',
            'widgetType': Tkinter.Button,
            'text':'Close',
            'wcfg':{'bd':4},
            'gridcfg':{'sticky':'nesw','row':-1,'column':1},
            'command':CallBackFunction(self.CloseHIST_cb, mol)})
        mol.histform = self.vf.getUserInput(ifd, modal=0,blocking=0)
        mol.histform.root.protocol('WM_DELETE_WINDOW',CallBackFunction(self.CloseHIST_cb,mol))
        mol.histform.ifd = ifd


    def resetRange(self, mol, event=None):
        mol.min.set(mol.r[0])
        mol.max.set(mol.r[1])


    def CloseHIST_cb(self, mol, event=None):
        """None<-CloseHIST_cb(mol, event=None)
           destroys mol.histform
        """
        if hasattr(mol, 'histform'):
            mol.histform.destroy()
            delattr(mol, 'histform')


ShowAutoDockStatesHISTOGRAMGUI = CommandGUI()
ShowAutoDockStatesHISTOGRAMGUI.addMenuCommand('AutoToolsBar', 
    menuText['AnalyzeMB'], menuText['showStatesHISTOGRAMMB'], 
    cascadeName=menuText['StatesMB'])

#################################END OF SHOWSTATESHISTOGRAM###################


class ShowAutoDockClusteringStates(ShowAutoDockStatesBaseCmd):
    """allows user to show states in a CLUSTERING of docked states of a ligand"""


    def onRemoveObjectFromViewer(self, obj):
        if hasattr(obj, 'ifd'):
            delattr(obj, 'ifd')
        if hasattr(obj, 'form'):
            obj.form.destroy()
        if hasattr(obj, 'clustNB'):
            obj.clustNB.master.destroy()
        for item in [ 'clustform', 'clustNB','ifd','ifd2','nbins']:
            if hasattr(obj, item):
                delattr(obj, item)
        if hasattr(obj, 'spw'):
            if hasattr(obj.spw, 'sequenceList'):
                for c in obj.spw.sequenceList:
                    if hasattr(c, 'mol'):
                        delattr(c, 'mol')
            if hasattr(obj.spw, 'mol'):
                delattr(obj.spw, 'mol')
        if obj in self.molList:
            ind = self.molList.index(obj)
            self.molList.remove(obj)
            d = self.dockings[ind]
            if hasattr(d, 'clusterer'):
                if hasattr(d.clusterer, 'subset'):
                    delattr(d.clusterer, 'subset')
            self.dockings.remove(d)
        if hasattr(obj, 'torTree'):
            for key in obj.torTree.rootNode.__dict__.keys():
                obj.torTree.rootNode.__dict__[key] = None
            for item in obj.torTree.torsionMap:
                for key in item.__dict__.keys():
                    item.__dict__[key] = None


    def showClust(self, mol, docking, event=None):
        if not hasattr(mol.clustNB, 'currentNode'):
            return
        titleStr = mol.name + ' rms '+str(mol.cur_rms)+' clustering:'+ \
                mol.clustNB.currentNode.name + ' cluster'
        #THIS IS WHERE spw is built, updated, and/or deiconified

        #currentInfo is a list of indices
        cI = mol.clustNB.currentInfo
        curNode = mol.clustNB.currentNode
        #name of node is '0' etc... represents which cluster
        #in clustering is active
        curName = curNode.name
        ind = int(curName)

        # create a list of conformations for this cluster
        clist = []
        l = docking.clusterer.clustering_dict[mol.cur_rms]

        # create a list of ids for the conformations
        #start with a zero to enable reset to original state
        idList = ['0']
        cluName = str(int(curName)+1)
        for i in range(len(l[ind])):
            clist.append(l[ind][i])
            idList.append(cluName+'_'+str(i+1))

        if hasattr(mol, 'spw'):
            mol.spw.update(clist, idList)
            mol.spw.nextFrame(1)
        else:
            #the first time
            if self.vf.userpref['Player GUI']['value']=='StatesPlayerWidget':
                mol.spw = StatesPlayerWidget(mol, docking, self.vf, mol.name, clist, 
                            idList, ask=self.ask)
            else:
                mol.spw = ConformationPlayer(mol, docking, self.vf, mol.name, clist, 
                            idList=idList, form2=1, 
                            ask=self.ask)
							#buttonMask={'recordB':False})
                mol.spw.nextFrame(1)
                            #mol.name, docking.ch.conformations, idList=idList, 
        # make sure widget is visible: de-iconify if iconified.
        if not mol.spw.form.f.winfo_ismapped():
            mol.spw.form.deiconify()


    def sortData(self, clusterer):
        #this method should not be called
        if hasattr(clusterer, 'sorted_data'):
            return
        newL = []
        for i in clusterer.argsort:
            newL.append(clusterer.data[int(i)])
        clusterer.data = newL
        clusterer.sorted_data = 1


    def setClusterer(self, docking):
        ifd= InputFormDescr(title='Select Clusterer')
        ifd.append({'widgetType':'ListChooser',
            'name':'clusterSel',
            'entries':docking.clusterer_dict.keys(),
            'mode':'single',
            'title':'energy used for clustering\n(if subset, name_energy)',
            'lbwcfg':{'height':4},
            'gridcfg':{'sticky':'w', 'column':100,
            'rowspan':10}})
        val = self.vf.getUserInput(ifd)
        if len(val)>0 and len(val['clusterSel'])>0:
            #print 'set clusterer to ', val['clusterSel'][0]
            clustKey = val['clusterSel'][0]
            if clustKey not in ['binding', 'docking']:
                #need to setup the subset stuff
                #changed so that set names with underscores would be ok
                ind = rfind(clustKey, '_')
                key = clustKey[:ind]
                #key = split(clustKey, '_')[0]
                docking.clusterer = docking.clusterer_dict[clustKey]
                #update conformations-their subset and subset_coords
                nodes = self.vf.ADanalyze_makeSubsetClustering.vf.sets[key]
                atoms = nodes.findType(Atom)
                for conf in self.vf.docked.ch.conformations:
                    if hasattr(conf, 'subset_coords'):
                        delattr(conf, 'subset_coords')
                    conf.subset = atoms
                    junk = conf.getCoords_subset()
                docking.clusterer.usesSubset = key
                docking.clusterer.subset = atoms
            else:
                docking.clusterer = docking.clusterer_dict[clustKey]
                docking.clusterer.usesSubset = 0
        

    def doit(self, mol, ask=1):
        """None<-doit(mol, ask)
            this is where clustering sequences are set up 
        self.ADanalyze_showClusteringStates
        """
        mols = self.vf.expandNodes(mol)
        if len(mols)==0:
            return 'ERROR'
        mol = mols[0]
        if hasattr(mol, 'vina_energy'):
            self.warningMsg("This is a vina result. Clustering is not available ")
            return "ERROR"
        self.ask = ask
        #in inherited guiCallback method: 
        # self.molList and self.dockings are built like this
        #    for k, d in self.vf.dockings.items():
        #       self.molList.append(d.ligMol)
        #       self.dockings.append(d)
        ind = self.molList.index(mol)
        docking = self.dockings[ind]
        mol.stop2 = 0
        docking_has_clusterer = hasattr(docking, 'clusterer')
        if not docking_has_clusterer or len(docking.clusterer.clustering_dict.keys())==0:
            #check whether a file of clusterings is to be read
            msg = 'Clusterer of current docking has no clusterings to show. Recluster or read in a previous clustering first...'
            self.warningMsg(msg)
            return "ERROR"
        #check how many clusterers the docking has
        num_clusterers = len(docking.clusterer_dict.keys())
        #if more than one, choose which one to use
        if num_clusterers >1:
            self.setClusterer(docking)
        #check how many clusterings the clusterer has
        entries = docking.clusterer.clustering_dict.keys()
        num_entries = len(entries)
        if num_entries==0:
            msg = "no available clusterings to show"
            self.warningMsg(msg)
            return 'ERROR'
        elif num_entries==1:
            #if there is only one available rms, just use it
            self.showClustering(entries[0], mol, docking)
            return 
        entries = map(lambda x:'%6.3f'%x, entries)
        entries.sort()
        #sort and format keys here
        titleStr = 'Show ' + mol.name + ' Clusterings'
        ifd =  InputFormDescr(title=titleStr)
        ifd.append({'name': 'lc',
            'widgetType': 'ListChooser',
            'entries': entries,
            'command': CallBackFunction(self.setClustering, mol, docking),
            'title': 'Select RMS Tolerance',
            'wcfg': { 'mode':self.mode},
            'gridcfg':{'sticky':'nesw','columnspan':2} })
        ifd.append({'name': 'closeB',
            'widgetType': Tkinter.Button,
            'text':'Close',
            'wcfg':{'bd':4},
            'gridcfg':{'sticky':'nesw'},
            'command':CallBackFunction(self.CloseCF_cb, mol)})
        mol.clustform = self.vf.getUserInput(ifd, modal=0,blocking=0)
        mol.clustform.root.protocol('WM_DELETE_WINDOW', CallBackFunction(self.CloseCF_cb,mol))
        mol.clustform.ifd = ifd
        w = mol.clustform.ifd.entryByName['lc']['widget']
        self.lb = w.lb


    def setClustering(self, mol, docking, event=None):
        #FIX THIS!!!!
        #get the list of confs corresponding to lb.curselection()
        ind = self.lb.curselection()
        if len(ind):
            rms = float(self.lb.get(ind))
            self.showClustering(rms, mol, docking)


    def showClustering(self, rms, mol, docking):
        mol.cluSEQ = docking.clusterer.clustering_dict[rms]
        #cluSEQ is a list of lists of conformations
        mol.cur_rms = rms
        dataList = []
        reverseList = []
        #rLctr is 0...n, used by interactiveHistogram
        rLctr = 0
        confL = docking.ch.conformations
        e = docking.clusterer.energy_used
        for l in mol.cluSEQ:
            if docking.version>=4.0:
                e = str(docking.version)
                try:
                    energy = l[0].energy
                except:
                    energy = l[0].binding_energy
            elif e=='binding':
                energy = l[0].binding_energy
            else:
                energy = l[0].docking_energy
            dataList.append([float(energy), len(l)])
            #this ends up being 0-6,7,8-9
            reverseList.append(range(rLctr,rLctr+len(l)))
            rLctr = rLctr + len(l)
        tstr = mol.name + ': rms =  ' + str(rms) + '  clustering'
        top = Tkinter.Toplevel(master=self.vf.GUI.ROOT)
        top.title(tstr)
        #NOW build the histogram
        if e[0]=='d':
            xlabel = 'DOCKING ENERGY'
        elif e=='4.0':
            xlabel = 'ENERGY'
        else:
            xlabel = 'BINDING ENERGY'
        mol.clustNB = interactiveHistogramGraph.InteractiveHistogramGraph(mol.name,
            master=top, nodeList = dataList, reverseIndex=reverseList,
            label_text=mol.name + ':' + str(rms) + ' rms', xlabel_text=xlabel, 
            ylabel_text='#\nC\nO\nN\nF\nO\nR\nM\nA\nT\nI\nO\nN\nS')
        mol.clustNB.draw.bind('<Button-1>', CallBackFunction(self.showClust, mol, docking),'+')
        try:
            mol.clustform.withdraw()
        except:
            pass


    def CloseCF_cb(self, mol, event=None):
        """None<-CloseCF_cb(mol, event=None)
           destroys mol.clustform
        """
        if hasattr(mol, 'clustform'):
            mol.clustform.destroy()
            delattr(mol, 'clustform')


ShowAutoDockStatesCLUSTERINGGUI = CommandGUI()
ShowAutoDockStatesCLUSTERINGGUI.addMenuCommand('AutoToolsBar', 
    menuText['AnalyzeMB'], menuText['showStatesCLUSTERINGMB'], 
    cascadeName=menuText['ClusteringMB'])



class ReadAutoDockClusteringStates(ShowAutoDockStatesBaseCmd):
    """allows user to read in a CLUSTERING from a 'clust' file"""

    def doit(self, mol, clustFile=None, ask=1):
        """None<-doit(mol, clustFile, ask)
            this is where clustering sequences are read in from a file 
        self.ADanalyze_readClusteringStates
        """
        mols = self.vf.expandNodes(mol)
        if len(mols)==0:
            return 'ERROR'
        mol = mols[0]
        self.ask = ask
        #NB: in guiCallback: self.molList and self.dockings are built like this
        #    for k, d in self.vf.dockings.items():
        #       self.molList.append(d.ligMol)
        #       self.dockings.append(d)
        ind = self.molList.index(mol)
        docking = self.dockings[ind]
        mol.stop2 = 0
        #FIX THIS!!!!
        #if not docking_has_clusterer or len(docking.clusterer.clustering_dict.keys())==0:
        #    #check whether a file of clusterings is to be read
        if hasattr(docking, 'version') and docking.version>=4.0:
            sort = 'energy'
        else:
            sort = 'docking'
        docking_has_clusterer = hasattr(docking, 'clusterer')
        if not clustFile:
            clustFile= self.vf.askFileOpen(types=[('select clustering filename:', '*.clust'), 
                            ('all files','*')],
                            title = 'Select clustering file')
            #if the user specifies a clust file, create clusterer to try to read it
            if clustFile:
                if not docking_has_clusterer:
                    docking.clusterer = Clusterer(docking.ch.conformations, 
                                            sort=sort)
                    docking.clusterer.usesSubset = 0
                    docking.clusterer_dict[sort] = clusterer
                try:
                    docking.clusterer.read(clustFile)
                except:
                    #if clustfile is on binding energy, need to have a
                    #clusterer, whose argsort is on binding_energy
                    if sort=='docking' and not 'binding' in docking.clusterer_dict.keys():
                        clusterer = Clusterer(docking.ch.conformations, 
                                            sort='binding')
                        docking.clusterer = clusterer
                        clusterer.usesSubset = 0
                        docking.clusterer_dict['binding'] = clusterer
                    docking.clusterer.read(clustFile)
                self.vf.docked = docking
            else:
                #if read is cancelled, return msg + quit
                return 'ERROR'


ReadAutoDockStatesCLUSTERINGGUI = CommandGUI()
ReadAutoDockStatesCLUSTERINGGUI.addMenuCommand('AutoToolsBar', 
    menuText['AnalyzeMB'], menuText['readStatesCLUSTERINGMB'], 
    cascadeName=menuText['ClusteringMB'])



class WriteAutoDockStates(MVCommand):
    """allows user to write Entropia type clustered docking results
    """

    def onAddCmdToViewer(self):
        checkHasInitializedDockings(self.vf)
        from AutoDockTools.autoanalyzeCommands import StatesPlayerWidget
        self.mode = 'single'
        self.title = 'Choose Molecule'
        self.molList = []
        self.dockings = []

    def setClusterer(self, docking):
        ifd= InputFormDescr(title='Select Clusterer')
        ifd.append({'widgetType':'ListChooser',
            'name':'clusterSel',
            'entries':docking.clusterer_dict.keys(),
            'mode':'single',
            'title':'energy used for clustering',
            'lbwcfg':{'height':4},
            'gridcfg':{'sticky':'w', 'column':100,
            'rowspan':10}})
        val = self.vf.getUserInput(ifd)
        if len(val)>0 and len(val['clusterSel'])>0:
            #print 'set clusterer to ', val['clusterSel'][0]
            docking.clusterer = docking.clusterer_dict[val['clusterSel'][0]]


    def buildMolList(self):
        molList = self.molList = []
        for k, d in self.vf.dockings.items():
            self.molList.append(d.ligMol)
            self.dockings.append(d)


    def guiCallback(self):
        """called each time the 'Write Results File' button is pressed"""
        if not len(self.vf.Mols):
            msg = 'no molecules in viewer'
            self.vf.warningMsg(msg)
            return 'ERROR'
        self.dockings = []
        self.buildMolList()
        #molList = self.molList = []
        ##FIX THIS!!!
        #for k, d in self.vf.dockings.items():
            #self.molList.append(d.ligMol)
            #self.dockings.append(d)
        molList = self.molList
        if not len(molList):
            msg = 'no molecules with dockings in viewer'
            self.vf.warningMsg(msg)
            return 'ERROR'
        elif len(molList)==1:
            mol = molList[0]
            #if there are more than one clusterers available, make user pick
            if len(mol.docking.clusterer_dict.keys())>1:
                self.setClusterer(mol.docking)
            filename = self.vf.askFileSave(types=[('docked results file','*.res')], title = 'Docked Result File:')
            if filename:
                self.doitWrapper(molList[0], filename, log=1, redraw=0)
            #self.doitWrapper(molList[0], log=1, redraw=0)
        else:
            self.chooser = MoleculeChooser(self.vf, self.mode, self.title)
            self.chooser.ipf.append({'name':'Select Button',
                                     'widgetType':Tkinter.Button,
                                     'text':'Select Molecule',
                                     'wcfg':{'bd':6},
                                     'gridcfg':{'sticky':Tkinter.E+Tkinter.W},
                                     'command': self.chooseMolecule_cb})
            self.form = self.chooser.go(modal=0, blocking=0)
            lb = self.chooser.ipf.entryByName['Molecule']['widget'].lb
            lb.bind("<Double-Button-1>",self.chooseMolecule_cb)

            
    def chooseMolecule_cb(self, event = None):
        """
        invoked from guiCallback"""
        mol = self.chooser.getMolSet()
        if mol in self.molList:
            if len(mol.docking.clusterer_dict.keys())>1:
                self.setClusterer(mol.docking)
            self.chooser.form.withdraw()
            filename = self.vf.askFileSave(types=[('docked results file','*.res')], title = 'Docked Result File:')
            if filename:
                self.doitWrapper(mol, filename, log=1, redraw=0)
        else:
            msg = mol.name + ' is not connected to a docking'
            self.vf.warningMsg(msg)
            return "ERROR"


    def __call__(self, mol, filename, **kw):
        mols = self.vf.expandNodes(mol)
        if len(mols)==0:
            return 'ERROR'
        self.buildMolList()
        apply(self.doit, (mols[0], filename), kw)


    def doit(self, mol, filename):
        """None<-doit(mol)
NB:conformations are written in order of the argsort matrix built in the
clusterer so that they will correspond with the written clust file
        """
        ind = self.molList.index(mol)
        #mol = docking.ligMol
        docking = self.dockings[ind]
        if not len(docking.ch.conformations):
            return 'ERROR'
        fptr = open(filename, 'w')
        #for c in docking.ch.conformations:
            #c.writeRes(fptr)
        #at this point, use the order in the clusterer argsort
        cl = docking.clusterer
        for i in range(len(cl.argsort)):
            conf = cl.data[int(cl.argsort[i])]
            if not hasattr(conf, 'mol'):
                conf.mol = mol
                #FIX THIS!!!
                conf.mol.allAtoms = mol.chains.residues.atoms
            conf.writeRes(fptr)
        fptr.close()



WriteAutoDockStatesGUI = CommandGUI()
WriteAutoDockStatesGUI.addMenuCommand('AutoToolsBar', 
    menuText['AnalyzeMB'], menuText['writeResultMB'], 
    cascadeName=menuText['StatesMB'])



class WriteAutoDockClustering(MVCommand):
    """allows user to write Clustering results
    """

    def onAddCmdToViewer(self):
        checkHasInitializedDockings(self.vf)
        from AutoDockTools.autoanalyzeCommands import StatesPlayerWidget
        self.mode = 'single'
        self.title = 'Choose Molecule'
        self.molList = []
        self.dockings = []


    def onRemoveObjectFromViewer(self, obj):
        if obj in self.molList:
            ind = self.molList.index(obj)
            self.molList.remove(obj)
            d = self.dockings[ind]
            self.dockings.remove(d)


    def guiCallback(self):
        """called each time the 'Write Clust' button is pressed"""
        if not len(self.vf.Mols):
            msg = 'no molecules in viewer'
            self.vf.warningMsg(msg)
            return 'ERROR'
        self.dockings = []
        molList = self.molList = []
        #FIX THIS!!!
        for k, d in self.vf.dockings.items():
            self.molList.append(d.ligMol)
            self.dockings.append(d)
        if not len(molList):
            msg = 'no molecules with dockings in viewer'
            self.vf.warningMsg(msg)
            return 'ERROR'
        elif len(molList)==1:
            mol = molList[0]
            if len(mol.docking.clusterer_dict.keys())>1:
                self.setClusterer(mol.docking)
            filename = self.vf.askFileSave(types=[('clustering results file','*.clust')], title = 'Clustering Result File:')
            if filename:
                self.doitWrapper(molList[0], filename, log=1, redraw=0)
            #self.doitWrapper(molList[0], log=1, redraw=0)
        else:
            self.chooser = MoleculeChooser(self.vf, self.mode, self.title)
            self.chooser.ipf.append({'name':'Select Button',
                                     'widgetType':Tkinter.Button,
                                     'text':'Select Molecule',
                                     'wcfg':{'bd':6},
                                     'gridcfg':{'sticky':Tkinter.E+Tkinter.W},
                                     'command': self.chooseMolecule_cb})
            self.form = self.chooser.go(modal=0, blocking=0)
            lb = self.chooser.ipf.entryByName['Molecule']['widget'].lb
            lb.bind("<Double-Button-1>",self.chooseMolecule_cb)

            
    def chooseMolecule_cb(self, event = None):
        """
        invoked from guiCallback"""
        mol = self.chooser.getMolSet()
        if mol in self.molList:
            self.chooser.form.withdraw()
            if len(mol.docking.clusterer_dict.keys())>1:
                self.setClusterer(mol.docking)
            filename = self.vf.askFileSave(types=[('clustering results file','*.clust')], title = 'Clustering Result File:')
            if filename:
                self.doitWrapper(mol, filename, log=1, redraw=0)
        else:
            msg = mol.name + ' is not connected to a docking'
            self.vf.warningMsg(msg)
            return "ERROR"



    def __call__(self, mol, filename, **kw):
        mols = self.vf.expandNodes(mol)
        if len(mols)==0:
            return 'ERROR'
        apply(self.doit, (mols[0], filename), kw)


    def setClusterer(self, docking):
        ifd= InputFormDescr(title='Select Clusterer')
        ifd.append({'widgetType':'ListChooser',
            'name':'clusterSel',
            'entries':docking.clusterer_dict.keys(),
            'mode':'single',
            'title':'energy used for clustering',
            'lbwcfg':{'height':4},
            'gridcfg':{'sticky':'w', 'column':100,
            'rowspan':10}})
        val = self.vf.getUserInput(ifd)
        if len(val)>0 and len(val['clusterSel'])>0:
            #print 'set clusterer to ', val['clusterSel'][0]
            docking.clusterer = docking.clusterer_dict[val['clusterSel'][0]]


    def doit(self, mol, filename):
        """None<-doit(mol)
        """
        mols = self.vf.expandNodes(mol)
        if len(mols)==0:
            return 'ERROR'
        mol = mols[0]
        ind = self.molList.index(mol)
        #mol = docking.ligMol
        docking = self.dockings[ind]
        if not len(docking.clusterer.clustering_dict.keys()):
            return 'ERROR'
        docking.clusterer.write(filename)


WriteAutoDockClusteringGUI = CommandGUI()
WriteAutoDockClusteringGUI.addMenuCommand('AutoToolsBar', 
    menuText['AnalyzeMB'], menuText['writeClusteringMB'], 
    cascadeName=menuText['ClusteringMB'])



class MakeAutoDockCLUSTERING(MVCommand):
    """allows user to make clustering of docking results at specified rms tolerances"""


    def guiCallback(self):
        """called each time the 'make clustering ' button is pressed"""
        #if there are more than one dockings, allow user to choose
        if len(self.vf.dockings.keys())>1:
            self.vf.ADanalyze_selectDLG.guiCallback()
        #need to get at the clusterer, an outputfile name and a tolerance list:
        if not self.vf.docked:
            self.vf.warningMsg('Please Read a Docking Log First')
            return
        if hasattr(self.vf.docked.ligMol, 'vina_energy'):
            self.warningMsg("Current docking is a vina result. Clustering is not available ")
            return "ERROR"
        if not hasattr(self.vf.docked, 'clusterer'):
            self.vf.warningMsg('current docking has no clusterer')
            return 'ERROR'
        clusterer = self.vf.docked.clusterer
        if not len(clusterer.data):
            self.vf.warningMsg('current clusterer has no conformations')
            return 'ERROR'
        #check that docked.molecule is in its input conformation...
        if not hasattr(self.vf.docked, 'ligMol'):
            self.vf.warningMsg('current docking has no ligand!')
            return 'ERROR'
        a0 = self.vf.docked.ligMol.allAtoms[0]
        conf = a0.conformation
        crds = a0._coords
        current_crds = a0.coords
        #check that the ligand is in its input conformation
        #if len(crds)>1 and (crds[0][0]!=current_crds[0] or \
        #    crds[0][1]!=current_crds[1] or \
        #    crds[0][2]!=current_crds[2]) :
            #t = 'Warning!\nLigand not in input conformation.\nDo you want to cluster anyway?' 
            #d = SimpleDialog(self.vf.GUI.ROOT, text = t,
            #        buttons = ["Yes", "Cancel"], default = 0, 
            #        title = "Cluster on non-input coordinates?")
            #dontCluster = d.go()
            #if dontCluster>0:
            #    return 'ERROR'
        if not hasattr(self, 'form'):
            self.buildFORM()
            self.form = self.vf.getUserInput(self.ifd, modal=0,blocking=0)
            self.form.root.protocol('WM_DELETE_WINDOW',self.Cancel_cb)
        else:
            self.form.deiconify()


    def buildFORM(self):
        # entry for tol
        # button to designate file name
        # button to go or cancel
        d = self.vf.docked
        version = d.version
        mol = d.ligMol
        self.tolStr = Tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.tolStr.set('0.5 2.0')
        self.sort = Tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.sort.set('docking')
        self.outputFile = Tkinter.StringVar(master=self.vf.GUI.ROOT)
        titleStr = 'Cluster ' + mol.name + ' Conformations'
        ifd =  self.ifd = InputFormDescr(title = titleStr)
        ifd.append({'name': 'tolEnt',
            'widgetType': Tkinter.Entry,
            'wcfg':{'bd':4,
                'label': 'Tolerances:',
                'textvariable': self.tolStr,
                },
            'gridcfg':{'sticky':'nesw','columnspan':2} })
        if version>=4.0:
            ifd.append({'name':'dsortbut',
                    'widgetType':Tkinter.Label,
                    'wcfg':{'text':'Clustering on AD4 Energy'}, 
                    'gridcfg':{'sticky':Tkinter.E}})
        else:
            ifd.append({'name':'dsortbut',
                    'widgetType':Tkinter.Radiobutton,
                    'wcfg':{'text':'docking energy', 
                        'value':'docking',
                        'variable':self.sort},
                    'gridcfg':{'sticky':Tkinter.E}})
            ifd.append({'name':'bsortbut',
                    'widgetType':Tkinter.Radiobutton,
                    'wcfg':{'text':'binding energy', 
                        'value':'binding',
                        'variable':self.sort},
                    'gridcfg':{'sticky':Tkinter.E}})
        ifd.append({'name': 'fileEnt',
            'widgetType': Tkinter.Entry,
            'wcfg':{'bd':4,
                'label': 'Outputfile Name:',
                'textvariable': self.outputFile,
                },
            'gridcfg':{'sticky':'nesw','columnspan':2} })
        ifd.append({'name': 'closeB',
            'widgetType': Tkinter.Button,
            'wcfg':{
                'text':'OK',
                'command':self.OK_cb,
                'bd':4},
            'gridcfg':{'sticky':Tkinter.E+Tkinter.W}}),
        ifd.append({'name': 'cancelB',
            'widgetType': Tkinter.Button,
            'wcfg':{
                'bd':4,
                'text':'Cancel',
                'command':self.Cancel_cb},
            'gridcfg':{'sticky':Tkinter.E+Tkinter.W, 'row':-1, 'column':1}})


    def OK_cb(self, event=None):
        tolList = map(float, split(self.tolStr.get()))
        filename = strip(self.outputFile.get())
        sort = self.sort.get()
        self.form.withdraw()
        if not len(filename):
            self.doitWrapper(tolList, sort, None)
        else:
            self.doitWrapper(tolList, sort, filename)


    def Cancel_cb(self, event=None):
        self.form.withdraw()


    def __call__(self, toleranceList, sort='docking', outputFile=None, **kw):
        """None<-ADanalyze_makeClustering(toleranceList, sort, outputFile)
operates on current self.vf.docked object's clusterer object with data 
which are conformations
toleranceList: list of rms tolerances for clusterings
sort: type of energy for clustering, either docking or binding.
outputFile: file to get written clustering results, generally *.clust
        """
        if not self.vf.docked:
            return 'ERROR'
        if not hasattr(self.vf.docked, 'clusterer'):
            return 'ERROR'
        if not len(self.vf.docked.clusterer.data):
            return 'ERROR'
        apply(self.doitWrapper, (toleranceList, sort, outputFile,), kw)



    def doit(self, toleranceList, sort, outputFile):
        """None<-ADanalyze_makeClustering.doit(toleranceList, sort, outputFile)
toleranceList: space separated list of floats 
sort: type of energy for clustering, either docking or binding.
outputFile: optional name of file for written output
        """
        #if no conformations, return
        if not self.vf.docked:
            return 'ERROR'
        docked = self.vf.docked
        #is this necessary?
        #if not hasattr(docked, 'clusterer'):
        #    return 'ERROR'

        if hasattr(docked, 'version') and docked.version>=4.0:
            docked.clusterer = Clusterer(docked.ch.conformations, sort='energy')
            clusterer = docked.clusterer
            clusterer.usesSubset = 0
        elif sort[0]=='d':
            if docked.clusterer_dict.has_key('docking'):
                clusterer = docked.clusterer_dict['docking']
            else:
                clusterer = Clusterer(docked.ch.conformations,
                                        sort='docking')
                docked.clusterer_dict['docking'] = clusterer
                clusterer.usesSubset = 0
        else:
            if docked.clusterer_dict.has_key('binding'):
                clusterer = docked.clusterer_dict['binding']
            else:
                clusterer = Clusterer(docked.ch.conformations,
                                        sort='binding')
                docked.clusterer_dict['binding'] = clusterer
                clusterer.usesSubset = 0
        if not len(clusterer.data):
            print 'returning lack of data error'
            return 'ERROR'
        for tol in toleranceList:
            print 'clustering at ', tol
            clusterer.make_clustering(tol)
        if outputFile:
            clusterer.write(outputFile)
        #make docked's clusterer the latest clusterer used
        docked.clusterer = clusterer

        #set up rmsTool
        coords = docked.ligMol.allAtoms.coords[:]
        docked.clusterer.rmsTool = RMSDCalculator(coords)
        docked.clusterer.rmsToolRef = '0'
        ###WHAT IS THIS????
        if hasattr(docked.ligMol, 'spw') and \
                hasattr(docked.ligMol.spw, 'updateSortList'):
            docked.ligMol.spw.updateSortList()

MakeAutoDockCLUSTERINGGUI = CommandGUI()
MakeAutoDockCLUSTERINGGUI.addMenuCommand('AutoToolsBar', 
    menuText['AnalyzeMB'], menuText['makeCLUSTERINGMB'], 
    cascadeName=menuText['ClusteringMB'])



def getCoords_subset(conf):
    """Return coordinates of conf's current subset, if there is one.

    If the coordinates haven't been computed yet,
        then compute, save, and return them.
    Otherwise, return the previously-computed coordinates.
    """
    if not conf.subset:
        return conf.getCoords()
    if not hasattr(conf, 'subset_coords'):
        oldCoords = conf.mol.allAtoms.coords[:]
        oldConf = conf.mol.allAtoms[0].conformation
        conf.mol.allAtoms.updateCoords(conf.getCoords(), 
                    conf.mol.stoc.confIndex)
        conf.subset_coords = conf.subset.coords[:]
        conf.mol.allAtoms.updateCoords(oldCoords, oldConf)
    return conf.subset_coords


def getRMSD_subset(conf, refCoords=None):
    """Return RMSD of conf's subset relative to refCoords.

    If refCoords is not given, the original coordinates for the
    subset will be used as the reference.
    """
    if not refCoords:
        refCoords = getCoords_subset(conf)
    rmsd_calc = RMSDCalculator(refCoords)            
    return rmsd_calc.computeRMSD(getCoords_subset(conf))


def get_distance_subset(clusterer, a, b):
    """return RMSD between subsets of atoms in two conformations, a and b.
    """
    ax = clusterer.data.index(a)
    bx = clusterer.data.index(b)

    if clusterer.dist_matrix[ax][bx] >= 0.0:
        # return previously saved distance
        return clusterer.dist_matrix[ax][bx]
    else:
        # compute, save, and return distance
        dist = getRMSD_subset(a, getCoords_subset(b))
        clusterer.dist_matrix[ax][bx] = clusterer.dist_matrix[bx][ax] = dist
        return dist



class MakeAutoDockSubsetCLUSTERING(MakeAutoDockCLUSTERING):
    """allows user to make clustering of docking results at specified rms tolerances using a subset of the atoms in the docked ligand"""

    def onRemoveObjectFromViewer(self, obj):
        if self.vf.docked and hasattr(self.vf.docked, 'ligMol'):
            if obj!=self.vf.docked.ligMol:
                return
            d = self.vf.docked
            if hasattr(d, 'ch') and hasattr(d.ch, 'conformations'):
                for c in d.ch.conformations:
                    if hasattr(c, 'subset'):
                        delattr(c, 'subset')


    def guiCallback(self):
        """called each time the 'make clustering ' button is pressed"""
        #need to get at the clusterer, an outputfile name and a tolerance list:
        if not self.vf.docked:
            self.vf.warningMsg('Please Read a Docking Log First')
            return
        if hasattr(self.vf.docked.ligMol, 'vina_energy'):
            self.warningMsg("Current docking is a vina result. Clustering is not available ")
            return "ERROR"
        if not hasattr(self.vf.docked, 'clusterer'):
            self.vf.warningMsg('current docking has no clusterer')
            return 'ERROR'
        clusterer = self.vf.docked.clusterer
        if not len(clusterer.data):
            self.vf.warningMsg('current clusterer has no conformations')
            return 'ERROR'
        if not len(self.vf.sets.values()):
            self.vf.warningMsg('currently no saved sets')
            return 'ERROR'
        lig_ids = map(id, self.vf.docked.ligMol.allAtoms)
        set_entries = []
        for k, val in self.vf.sets.items():
            #check that the sets are of atoms in docked.ligMol
            ok = 1
            vset = val
            #vset = val[0]
            atset = vset.findType(Atom)
            ats_ids = map(id, atset)
            for vv in ats_ids:
                if vv not in lig_ids:
                    ok = 0
                    break
            if ok:
                set_entries.append(k)
        if not len(set_entries):
            self.vf.warningMsg('currently no saved sets of docked ligand')
            return 'ERROR'
        self.set_entries = set_entries
        if not hasattr(self, 'form'):
            self.buildFORM()
            self.form = self.vf.getUserInput(self.ifd, modal=0,blocking=0)
            self.form.root.protocol('WM_DELETE_WINDOW',self.Cancel_cb)
        else:
            self.form.deiconify()
            #update the set entries in ComBox here
            lb = self.ifd.entryByName['set']['widget']._list
            lb.delete(0,'end')
            for k in self.set_entries:
                lb.insert('end', k)


    def buildFORM(self):
        # entry for tol
        # button to designate file name
        # button to go or cancel
        d = self.vf.docked
        mol = d.ligMol
        self.tolStr = Tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.tolStr.set('0.5 2.0')
        self.sort = Tkinter.StringVar(master=self.vf.GUI.ROOT)
        version = d.version
        self.sort.set('docking')
        if version>=4.0:
            self.sort.set('energy')
        self.outputFile = Tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.keyName = Tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.keyName.set(self.set_entries[0])
        titleStr = 'Cluster ' + mol.name + ' Conformations'
        ifd =  self.ifd = InputFormDescr(title = titleStr)
        ifd.append({'name':'set',
                    'widgetType':Pmw.ComboBox,
                    'wcfg':{'label_text': 'available subsets',
                            'labelpos':'w',
                            'listheight':'80',
                            'entryfield_value': self.keyName.get(),
                            'scrolledlist_items':self.set_entries,},
                            #'selectioncommand':self.update, }
                    'gridcfg':{'sticky':'nesw','columnspan':2} })
        ifd.append({'name': 'tolEnt',
            'widgetType': Tkinter.Entry,
            'wcfg':{'bd':4,
                'label': 'Tolerances:',
                'textvariable': self.tolStr,
                },
            'gridcfg':{'sticky':'nesw','columnspan':2} })
        if version>=4.0:
            ifd.append({'name':'dsortbut',
                    'widgetType':Tkinter.Radiobutton,
                    'wcfg':{'text':'energy', 
                        'value':'energy',
                        'variable':self.sort},
                    'gridcfg':{'sticky':Tkinter.E}})
        else:
            ifd.append({'name':'dsortbut',
                    'widgetType':Tkinter.Radiobutton,
                    'wcfg':{'text':'docking energy', 
                        'value':'docking',
                        'variable':self.sort},
                    'gridcfg':{'sticky':Tkinter.E}})
            ifd.append({'name':'bsortbut',
                    'widgetType':Tkinter.Radiobutton,
                    'wcfg':{'text':'binding energy', 
                        'value':'binding',
                        'variable':self.sort},
                    'gridcfg':{'sticky':Tkinter.E}})
        ifd.append({'name': 'fileEnt',
            'widgetType': Tkinter.Entry,
            'wcfg':{'bd':4,
                'label': 'Outputfile Name:',
                'textvariable': self.outputFile,
                },
            'gridcfg':{'sticky':'nesw','columnspan':2} })
        ifd.append({'name': 'closeB',
            'widgetType': Tkinter.Button,
            'wcfg':{
                'text':'OK',
                'command':self.OK_cb,
                'bd':4},
            'gridcfg':{'sticky':Tkinter.E+Tkinter.W}}),
        ifd.append({'name': 'cancelB',
            'widgetType': Tkinter.Button,
            'wcfg':{
                'bd':4,
                'text':'Cancel',
                'command':self.Cancel_cb},
            'gridcfg':{'sticky':Tkinter.E+Tkinter.W, 'row':-1, 'column':1}})


    def OK_cb(self, event=None):
        tolList = map(float, split(self.tolStr.get()))
        filename = strip(self.outputFile.get())
        key = self.ifd.entryByName['set']['widget'].get()
        self.form.withdraw()
        sort = self.sort.get()
        if not len(filename):
            self.doitWrapper(key, tolList, sort, None)
        else:
            self.doitWrapper(key, tolList, sort, filename)


    def __call__(self, key, toleranceList, sort, outputFile=None, **kw):
        """None<-ADanalyze_makeSubsetClustering(key, toleranceList, outputFile)
key:  name of subset of atoms to be basis of clustering
toleranceList: list of rms tolerances for clusterings
sort: type of energy for clustering, either docking or binding.
outputFile: file to get written clustering results, generally *.clust
        """
        if key not in self.vf.sets.keys():
            msg = 'invalid set name:'+ key 
            self.vf.warningMsg(msg)
            return 'ERROR'
        if not self.vf.docked:
            msg = 'no current docking in viewer'
            self.vf.warningMsg(msg)
            return 'ERROR'
        if not len(toleranceList):
            msg = 'no rms tolerances specified'
            self.vf.warningMsg(msg)
            return 'ERROR'
        apply(self.doitWrapper, (key, toleranceList, sort, outputFile,), kw)



    def doit(self, key, toleranceList, sort, outputFile):
        """None<-ADanalyze_makeSubsetClustering.doit(key, 
                        toleranceList, sort, outputFile)
key:  name of subset to be basis of clustering
toleranceList: space separated list of floats 
sort: type of energy for clustering, either docking or binding.
outputFile: optional name of file for written output
        """
        nodes = self.vf.sets[key]
        #NB: atoms = subset used for this clustering
        atoms = nodes.findType(Atom)
        if not len(atoms):
            print 'no atoms found in key specified as basis for custom clustering'
            return 'ERROR'
        else:
            print 'making clustering with subset of ', len(atoms), ' atoms'
        if not self.vf.docked:
            print 'no current docking'
            return 'ERROR'
        #if no conformations, return
        if not len(self.vf.docked.ch.conformations):
            print 'no current conformations'
            return 'ERROR'
        docked = self.vf.docked
        if hasattr(docked, 'version') and docked.version>=4.0:
            clusterer = Clusterer(docked.ch.conformations, 
                            sort='energy')
        elif sort[0]=='d':
            clusterer = Clusterer(docked.ch.conformations,
                                    sort='docking')
        else:
            clusterer = Clusterer(docked.ch.conformations,
                                    sort='binding')
        keystr = key + '_' + sort
        docked.clusterer_dict[keystr] = clusterer
        docked.clusterer = clusterer
        clusterer.usesSubset = key
        clusterer.subset = atoms
        #set subset member of all conformations:
        for conf in docked.ch.conformations:
            conf.subset = atoms
            if hasattr(conf, 'subset_coords'):
                delattr(conf, 'subset_coords')
        #replace default _get_distance with _get_distance_subset
        #clusterer.set_get_distance(get_distance_subset)
        clusterer.set_get_distance(clusterer._get_distance_subset)
        for tol in toleranceList:
            print 'clustering at ', tol
            clusterer.make_clustering(tol)
        if outputFile:
            clusterer.write(outputFile)

        #set up rmsTool for subset
        docked.clusterer.rmsTool = RMSDCalculator(atoms.coords[:])
        docked.clusterer.rmsToolRef = '0'


MakeAutoDockSubsetCLUSTERINGGUI = CommandGUI()
MakeAutoDockSubsetCLUSTERINGGUI.addMenuCommand('AutoToolsBar', 
    menuText['AnalyzeMB'], menuText['makeSubsetCLUSTERINGMB'], 
    cascadeName=menuText['ClusteringMB'])



commandList = [
    {'name':'ADanalyze_readVSResult','cmd':ADReadVSResult(),'gui':ADReadVSResultGUI},
    {'name':'ADanalyze_readVinaResult','cmd':ADLoadVinaResult(),'gui':ADLoadVinaResultGUI},
    {'name':'ADanalyze_readDLG','cmd':ADGetDLG(),'gui':ADGetDLGGUI},
    {'name':'ADanalyze_readAllDLGInDirectory','cmd':ADGetDirDLGs(),'gui':ADGetDirDLGsGUI},
    {'name':'ADanalyze_selectDLG','cmd':ADSelectDLG(),'gui':ADSelectDLGGUI},
    {'name':'ADanalyze_deleteDLG','cmd':ADDeleteDLG(),'gui':ADDeleteDLGGUI},
    {'name':'ADanalyze_readMacromolecule','cmd':ADReadMacro(),'gui':ADReadMacroGUI},
    {'name':'ADanalyze_chooseMacromolecule','cmd':ADChooseMacro(),'gui':ADChooseMacroGUI},
    {'name':'ADanalyze_showDockingsAsSpheres','cmd':ADSeeSpots(),'gui':ADSeeSpotsGUI},
    {'name':'ADanalyze_showBindingSite','cmd':ADShowBindingSite(),'gui':ADShowBindingSiteGUI},
    {'name':'ADanalyze_readStates','cmd':ReadAutoDockStates(),'gui':ReadAutoDockStatesGUI},
    {'name':'ADanalyze_showStates','cmd':ShowAutoDockStates(),'gui':ShowAutoDockStatesGUI},
    {'name':'ADanalyze_showStatesByEnergy','cmd':ShowAutoDockStatesByEnergy(),'gui':ShowAutoDockStatesByEnergyGUI},
    #{'name':'ADanalyze_showPopulation','cmd':ShowAutoDockPopulation(),'gui':ShowAutoDockPopulationGUI},
    {'name':'ADanalyze_chooseDockedConformations','cmd':ADDockingChooser(),'gui':ADDockingChooserGUI},
    {'name':'ADanalyze_showStatesHISTOGRAM','cmd':ShowAutoDockStatesHISTOGRAM(),'gui':ShowAutoDockStatesHISTOGRAMGUI},
    {'name':'ADanalyze_showResultsOutput','cmd':ADGetOutput(),'gui':ADGetOutputGUI},
    #{'name':'ADanalyze_showHistogram','cmd':ADDrawHistogram(),'gui':ADDrawHistogramGUI},
    {'name':'ADanalyze_getChart','cmd':ADMacroLigandChart(),'gui':ADMacroLigandChartGUI},
    #{'name':'ADanalyze_writeStates','cmd':WriteAutoDockStates(),'gui':WriteAutoDockStatesGUI},
    {'name':'ADanalyze_showClusteringStates','cmd':ShowAutoDockClusteringStates(),'gui':ShowAutoDockStatesCLUSTERINGGUI},
    {'name':'ADanalyze_readClusteringStates','cmd':ReadAutoDockClusteringStates(),'gui':ReadAutoDockStatesCLUSTERINGGUI},
    {'name':'ADanalyze_makeClustering','cmd':MakeAutoDockCLUSTERING(),'gui':MakeAutoDockCLUSTERINGGUI},
    {'name':'ADanalyze_makeSubsetClustering','cmd':MakeAutoDockSubsetCLUSTERING(),'gui':MakeAutoDockSubsetCLUSTERINGGUI},
    {'name':'ADanalyze_writeClustering','cmd':WriteAutoDockClustering(),'gui':WriteAutoDockClusteringGUI},
    {'name':'ADanalyze_writeVSResult','cmd':ADWriteVSResult(),'gui':ADWriteVSResultGUI},
    ]

try:
    from Pmv.Grid import AutoGrid, AutoGridSurfaceGui
    for i in [ #{'name':'ADanalyze_epdbMolecule',   'cmd':ADEPDBMol(), 'gui':ADEPDBMolGUI},
    {'name':'ADanalyze_addExtraGridIsocontour','cmd':ADGetAGrid(),'gui':ADGetAGridGUI}, {'name':'ADanalyze_showGridIsocontours','cmd':ADMakeAllGrids(),'gui':ADMakeAllGridsGUI}]:
        commandList.insert(7,i)
except:
    print 'skipping the isocontour-dependent commands'


def initModule(vf):
    for dict in commandList:
        vf.addCommand(dict['cmd'],dict['name'],dict['gui'])

    if hasattr(vf, 'GUI'):
        for item in vf.GUI.menuBars['AutoToolsBar'].menubuttons.values():
            item.configure(background = 'tan')
        if not hasattr(vf.GUI, 'adtBar'):
            vf.GUI.adtBar = vf.GUI.menuBars['AutoToolsBar']
            vf.GUI.adtFrame = vf.GUI.adtBar.menubuttons.values()[0].master
            





