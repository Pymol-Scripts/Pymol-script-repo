## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#############################################################################
#
# Author: Ruth HUEY, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################
#
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/autoflexCommands.py,v 1.66.2.1 2009/03/19 19:07:30 rhuey Exp $
#
# $Id: autoflexCommands.py,v 1.66.2.1 2009/03/19 19:07:30 rhuey Exp $
#
#
#
#
#
#
#

"""
This Module facilitates producing a  files for AutoDock. The steps in this process are:

    * Set the macromolecule: 

        o Read a PDBQT Macromolecule 

        o Choose Macromol...

    * Select which residues are to be flexible in macromolecule using Pmv selection tools:

        o ICOM Select 

        o SelectFromString

        o Select Spherical Region

    * The results of the previous steps are written to a file. The user selects a filename via a filebrowser.  
    
"""
import numpy.oldnumeric as Numeric

from DejaVu import viewerConst
from ViewerFramework.VFCommand import CommandGUI
##  from ViewerFramework.gui import InputFormDescr
from mglutil.gui.InputForm.Tk.gui import InputFormDescr
from mglutil.gui.InputForm.Tk.gui import InputFormDescr
from mglutil.gui.BasicWidgets.Tk.thumbwheel import ThumbWheel
from mglutil.gui.BasicWidgets.Tk.customizedWidgets import ListChooser,\
                                                ExtendedSliderWidget

from Pmv.mvCommand import MVCommand, MVBondICOM, MVAtomICOM
from MolKit.tree import TreeNode, TreeNodeSet
from MolKit.molecule import Atom, AtomSet, BondSet
from MolKit.protein import Residue, ResidueSet, Chain
from MolKit.pdbWriter import PdbqtWriter
from MolKit.bondSelector import RotatableBondSelector, AmideBondSelector
from MolKit.bondSelector import GuanidiniumBondSelector, LeafBondSelector

from Pmv.guiTools import MoleculeChooser
import types, string, Tkinter, os, Pmw
from AutoDockTools.autotorsCommands import MAXTORS, SetRotatableBonds
from AutoDockTools.atomTypeTools import AutoDock4_AtomTyper

from SimpleDialog import SimpleDialog

menuText = {}
menuText['AutoFlexMB'] = 'Flexible Residues'
menuText['InputMB'] = 'Input'
menuText['Read Macro'] = 'Open Macromolecule...'
menuText['Choose Macro'] = 'Choose Macromolecule...'
menuText['Set Residues'] = 'Choose Torsions in Currently Selected Residues...'
menuText['Set Hinge'] = 'Set up Hinge...'
menuText['Edit Hinge'] = 'Edit Hinge...'
menuText['Step Back'] = 'Redisplay Macromolecule'
menuText['WriteMB'] = 'Output'
menuText['writeFlexible'] = 'Save Flexible PDBQT...'
menuText['writeRigid'] = 'Save Rigid PDBQT...'
menuText['writeDir'] = 'Save Multiple Flexible PDBQTS...'



class AF_MacroReader(MVCommand):
    """ allows user to select a filename for the macromolecule"""


    def onRemoveObjectFromViewer(self, obj):
        if hasattr(self.vf, 'flexDict'):
            dict = self.vf.flexDict
            if dict.has_key('macrofilename'): 
                ok = False
                macrofilename = dict['macrofilename']
                for m in self.vf.Mols:
                    if m.parser.filename==macrofilename:
                        ok = True
                        break
                if not ok:
                    del dict['macrofilename']
            if dict.has_key('macroname') and dict['macroname'] not in self.vf.Mols.name:
                del dict['macroname']
                if dict.has_key('macromol') and dict['macromol']:
                    del dict['macromol']


    def onAddCmdToViewer(self):
        if not hasattr(self.vf, 'flexDict'):
            self.vf.flexDict={}
        if not self.vf.commands.has_key('readPDBQT'):
            self.vf.loadCommand('fileCommands', 'readPDBQT', 'Pmv')


    def guiCallback(self):
        """called each time the 'select pdbqt macromolecule' button is pressed"""
        macroFile = self.vf.askFileOpen(types=[('PDBQT files', '*.pdbqt')],
                title = 'PDBQT Macromolecule File:')
        if macroFile:
            filename=os.path.split(macroFile)[-1]
            ext = os.path.splitext(filename)[1]
            if ext!='.pdbqt':
                msg = 'File can only be PDBQT format'
                self.warningMsg(msg)
                return 'ERROR'
            self.doitWrapper(macroFile)


    def __call__(self, macroFile, **kw):
        """None<-ADflex_readMacro(macroFile)"""
        if not macroFile: return 'ERROR'
        ext = os.path.splitext(macroFile)[1]
        if ext!='.pdbqt':
            msg = 'File must be PDBQT format'
            self.warningMsg(msg)
            return 'ERROR'
        apply(self.doitWrapper, (macroFile,), kw)


    def doit(self, macroFile):
        mollist = self.vf.readPDBQT(macroFile, topCommand=0)
        if not len(mollist): return 'ERROR'
        mol = mollist[0]
        mol.allAtoms.used = 0
        dict = self.vf.flexDict
        dict['macrofilename'] = macroFile
        dict['macroname'] = mol.name
        dict['macromol'] = mol


AF_MacroReaderGUI=CommandGUI()
AF_MacroReaderGUI.addMenuCommand('AutoToolsBar', menuText['AutoFlexMB'], \
        menuText['Read Macro'], cascadeName = menuText['InputMB'])



class AF_MacroChooser(MVCommand):
    """ allows user to choose a molecule already present for the macromolecule"""


    def __init__(self, mode='single', title = 'Choose Macromolecule'):
        MVCommand.__init__(self)
        self.mode = mode
        self.title = title


    def onRemoveObjectFromViewer(self, obj):
        if hasattr(self.vf, 'flexDict'):
            dict = self.vf.flexDict
            if dict.has_key('macrofilename'): 
                ok = False
                macrofilename = dict['macrofilename']
                for m in self.vf.Mols:
                    if m.parser.filename==macrofilename:
                        ok = True
                        break
                if not ok:
                    del dict['macrofilename']
            if dict.has_key('macroname') and dict['macroname'] not in self.vf.Mols.name:
                del dict['macroname']
                if dict.has_key('macromol') and dict['macromol']:
                    del dict['macromol']



    def chooseMolecule_cb(self, event = None):
        """called each time the 'choose Molecule' button is pressed"""
        mols = self.chooser.getMolSet()
        kw = {'redraw':0}
        if mols: apply(self.doitWrapper, (mols,), kw)
        self.chooser.form.withdraw()


    def guiCallback(self):
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


    def __call__(self, nodes, **kw):
        """None<-ADflex_chooseMacro(nodes)"""
        nodes = self.vf.expandNodes(nodes)
        if not len(nodes): return 'ERROR'
        apply(self.doitWrapper, (nodes,), kw)


    def doit(self, nodes, **kw):
        nodes = self.vf.expandNodes(nodes)
        if not len(nodes): return 'ERROR'
        mol = nodes[0]
        mol.allAtoms.used=0
        #if mol is from a pdbqt file, do not need to add 'q' or 't'
        filetype = os.path.splitext(os.path.basename(mol.parser.filename))[1]
        msg = ""
        chargeMsg = ""
        typeMsg = ""
        if filetype!='.pdbqt':
            #make sure all atoms have charge: 'q'
            ats = mol.allAtoms.get(lambda x: x.chargeSet==None)
            if len(ats):
                mol.buildBondsByDistance()
                self.vf.computeGasteiger(mol, topCommand=0)
                chargeMsg = "added gasteiger charges "
            #make sure that all atoms have autodock_element: 't'
            ats = mol.allAtoms.get(lambda x: hasattr(x, 'autodock_element'))
            if len(ats)!= len(mol.allAtoms):
                ad4_typer = AutoDock4_AtomTyper()
                mol.buildBondsByDistance()
                ad4_typer.setAutoDockElements(mol)
                typeMsg = " added autodock4 atom types "
        if len(chargeMsg):
            msg += chargeMsg
            if len(typeMsg):
                msg = msg + ' and ' + typeMsg  + ' to ' + mol.name
        elif len(typeMsg):
            msg = typeMsg + ' to ' + mol.name
        hs = mol.allAtoms.get(lambda x: x.element=='H' and len(x.bonds))
        nphs = hs.get(lambda x: x.bonds[0].atom1.element=='C' or x.bonds[0].atom2.element=='C')
        if len(nphs):
            lenNPHS = 0
            beforeLen = len(mol.allAtoms)
            if 'automerge_nphs' in kw.keys():
                self.vf.mergeNPHSGC(mol.allAtoms)
            else:
                nphs_msg="There appear to be some nonpolar hydrogen in "+ mol.name+ "  Do you wish to merge them to conform to AutoDock4 atom types? "
            d = SimpleDialog(self.vf.GUI.ROOT, text=nphs_msg,
                buttons=['No', 'Yes'], default=1,
                title="Merge Non-polar Hydrogens?")
            mergeNPHS = d.go()
            if mergeNPHS:
                self.vf.mergeNPHSGC(mol.allAtoms)
                lenNPHS = beforeLen - len(mol.allAtoms)
            if lenNPHS:
                msg = msg + " and merged " + str(lenNPHS) + " non-polar hydrogens"
        if len(msg):
            self.warningMsg(msg)
        dict = self.vf.flexDict
        dict['macroname'] = mol.name
        dict['macrofilename'] = mol.parser.filename
        dict['macromol'] = mol
        if hasattr(self.vf, 'gpo') and hasattr(self.vf.gpo, 'receptor') and self.vf.gpo.receptor==mol:
            msg = mol.name + " is currently the 'macromolecule'\nin the Grid menu. Make sure the macromolecule specified in the gpf does not include the flexible residues!!"
            self.warningMsg(msg)


    def onPick(self,event):
        listChooser = self.ipf.entryByName['Molecule']['widget']
        tkListBox = listChooser.lb
        atom,geom = self.vf.findPickedAtom(event)
        if atom:
            pickedMol = atom.top
            #then need to make pickedMol the selection in self.lc
            for i in range(len(listChooser.entries)):
                listChooserlist=string.split(listChooser.entries[i][0])
                if pickedMol.name == listChooserlist[0]:
                    self.pickedMolIndex= i
                    tkListBox.select_clear(0,'end')
                    listChooser.select(i)
                    return
            t= "error: %s not in mv.Mols" %pickedMol.name
            self.vf.warningMsg(t)


AF_MacroChooserGUI=CommandGUI()
AF_MacroChooserGUI.addMenuCommand('AutoToolsBar', menuText['AutoFlexMB'],
            menuText['Choose Macro'], cascadeName = menuText['InputMB'])



class AF_SelectResidues(MVCommand):
    """ allows user to set up a set of residues in macromolecule whose sidechains are to be flexed in an autodock run"""


    def onAddCmdToViewer(self):
        if not hasattr(self.vf, 'flexDict'):
            self.vf.flexDict={}
        if not hasattr(self.vf, 'colorByAtomType'):
            self.vf.loadCommand('colorCommands', 'colorByAtomType', 'Pmv')
        self.torscount = 0


    def guiCallback(self):
        """called each time the 'Set Selected Residues' button is pressed"""
        if not self.vf.flexDict.has_key('macroname'): 
            t='select protein first'
            self.vf.warningMsg(t)
            return 'ERROR'
        else: macroname=self.vf.flexDict['macroname']
        if len(self.vf.selection)==0:
            msg = "Please select residues to be modelled as flexible first!"
            self.warningMsg(msg)
            return 'ERROR'
        nodes = self.vf.getSelection()
        nodes = nodes.findType(Residue).uniq()
        if not len(nodes):
            t='no current residues selected'
            self.vf.warningMsg(t)
            return 'ERROR'
        #warn if there is there is not a specific subselection
        mol = nodes[0].top
        title = "Process ALL residues in %s?" %mol.name 
        if len(nodes)==len(mol.chains.residues):
            msg="CAUTION: currently processing all the residues in "+mol.name+ " will be very time consuming.  Do you wish to continue? "
            d = SimpleDialog(self.vf.GUI.ROOT, text=msg,
                buttons=['No', 'Yes'], default=1,
                title=title)
            useAll = d.go()
            if not useAll:
                return "ERROR"
        if not nodes.__class__==ResidueSet:
            self.vf.setIcomLevel(Residue, topCommand=0)
            nodes = nodes.findType(Residue).uniq()
        mol = nodes[0].top
        mol.allAtoms.used=0
        kw = {'redraw':0}
        self.torscount = 0
        apply(self.doitWrapper, (nodes,), kw)


    def __call__(self, nodes=None, **kw):
        nodes = self.vf.expandNodes(nodes)
        if not len(nodes): return "ERROR"
        if nodes.__class__ != ResidueSet:
            nodes = nodes.findType(Residue).uniq()
        apply(self.doitWrapper, (nodes,), kw)


    def doit(self, nodes):
        flex_residues = self.vf.expandNodes(nodes)
        if not len(flex_residues): return 'ERROR'

        #remove all prolines
        proList = ResidueSet(filter(lambda x: x.type=='PRO', flex_residues))
        if proList:
            flex_residues = flex_residues - proList

        #remove all waters
        h20List = ResidueSet(filter(lambda x: x.type=='HOH', flex_residues))
        if h20List:
            flex_residues = flex_residues - h20List

        if not len(flex_residues):
            t='No non-water and non-proline Residues selected!'
            self.vf.warningMsg(t)
            return 'ERROR'

        map(self.setAutoFlexFields, flex_residues)
        #map(self.setSideChain, flex_residues)
        #map(self.setTorsionFields,flex_residues)
        #for item in flex_residues:
        #    self.setTorsionFields(item)
        #map(self.getAmideBonds,flex_residues)
        #remove any residues with no possible torsions
        flexList = filter(lambda x: x.torscount!=0, flex_residues)
        flexList = ResidueSet(flexList)
        if not len(flexList):
            t='Current Selected Residues have no active torsions!'
            self.vf.warningMsg(t)
            return 'ERROR'
        self.torscount = Numeric.add.reduce(flexList.torscount)
        mol = flex_residues[0].top
        allResidues = mol.findType(Residue)
        rigidResidues = allResidues-flex_residues
        dict = self.vf.flexDict
        dict['flex_residues'] = flex_residues
        dict['flex_residues_number'] = len(flex_residues)
        dict['rigidResidues'] = rigidResidues
        self.vf.ADflex_processResidues.guiCallback()


    def setAutoFlexFields(self, res):
        #process residues
        if hasattr(res, 'setup'): 
            return
        res.setup = 1
        res.atoms.used = 0
        res.atoms.bonds[0].possibleTors = 0
        res.atoms.bonds[0].activeTors = 0
        backbone_names = ['C','N','O','HN','HN1','HN2', 'HA', 
                    'H1','H2','H3','HO', 'H']
        #includes CA
        sidechain = res.atoms.get(lambda x: x.name not in backbone_names)
        res.sideChain = sidechain
        bondlist = res.bondlist = sidechain.bonds[0]
        bondlist.leaf = 0
        bondlist.possibleTors = 0
        bondlist.activeTors = 0
        rbs = RotatableBondSelector()
        rotatables = rbs.select(bondlist)
        for b in rotatables:
            b.possibleTors = 1
            b.activeTors = 1
        amides = AmideBondSelector().select(bondlist)
        for b in amides:
            b.activeTors = 0
            b.possibleTors = 1
        guanidiniums = GuanidiniumBondSelector().select(bondlist)
        for b in guanidiniums:
            b.activeTors = 0
            b.possibleTors = 1
        leaves = LeafBondSelector().select(bondlist)
        for b in leaves:
            b.activeTors = 0
            b.possibleTors = 0
        res.torscount = len(bondlist.get(lambda x: x.activeTors==1))
        #this field is not used in AutoDock4
        res.torsdof = res.torscount
        res.torscount = len(bondlist.get(lambda x: x.activeTors==1))
        res.torsdof = res.torscount
        
        caAtoms = res.atoms.get(lambda x: x.name=='CA')
        #get returns an AtomSet
        if caAtoms:  #this checks for len(caAtoms)
            res.rootlist = caAtoms
        elif self.vf.hasGui:
            #use inputForm to get root atom and atoms for flexible part
            rootname = Tkinter.StringVar(master=self.vf.GUI.ROOT)
            if hasattr(res, 'rootlist'):
                rootname.set(res.rootlist[0].name)
            else:
                at0 = res.atoms.get(lambda x: x._uniqIndex==0)[0]
                rootname.set(at0.name)
            s = 'Set Root Atom for ' + res.name
            ifd = InputFormDescr(title=s)
            ifd.append({'name':'sideChainLC',
                'widgetType': 'ListChooser',
                'mode': 'single',
                'entries':res.atoms.name,
                'title': 'Select Root Atom',
                'lbwcfg':{'height':20,'selectforeground':'red','exportselection':0},
                'gridcfg':{'sticky':Tkinter.W +Tkinter.E}}),
            vals= self.vf.getUserInput(ifd, modal=1, blocking=1)
            if vals:
                try:
                    atList = vals['sideChainLC']
                    res.rootlist = res.atoms.get(lambda x, atList=atList: x.name==atList[0])
                    res.sideChain = res.atoms
                except:
                    msg = 'rootatom not in res, using defaults'
                    self.vf.warningMsg(msg)
                    res.rootlist = AtomSet([res.atoms.get(lambda x: x._uniqIndex == 0)[0]])
                    res.sideChain = res.atoms
            else:
                msg = 'rootatom not in res, using defaults'
                self.vf.warningMsg(msg)
                res.rootlist = AtomSet([res.atoms.get(lambda x: x._uniqIndex == 0)[0]])
                res.sideChain = res.atoms

        else:
            res.rootlist = AtomSet([res.atoms.get(lambda x: x._uniqIndex == 0)[0]])
            res.sideChain = res.atoms


AF_SelectResiduesGUI = CommandGUI()
AF_SelectResiduesGUI.addMenuCommand('AutoToolsBar', menuText['AutoFlexMB'],menuText['Set Residues'])



class AF_ProcessResidues(SetRotatableBonds, MVBondICOM):
    """ allows user to process interactively a set of residues in macromolecule whose sidechains are to be flexed in an autodock run"""


    def __init__(self, func=None):
        SetRotatableBonds.__init__(self)
        MVBondICOM.__init__(self)
        self.save = None
        self.guiUp = 0
        self.pickLevel = 'parts'


    def onAddCmdToViewer(self):
        if not hasattr(self.vf, 'flexDict'):
            self.vf.flexDict = {}
        if not hasattr(self.vf, 'colorByAtomType'):
            self.vf.loadCommand('colorCommands', 'colorByAtomType', 'Pmv')
        if not hasattr(self.vf, 'setICOM'):
            self.vf.loadCommand('interactiveCommands', 'setICOM', 'Pmv')
        if self.vf.hasGui and not self.torsStr:
            self.torsStr = Tkinter.StringVar(master=self.vf.GUI.ROOT)


    def guiCallback(self):
        """called INDIRECTLY each time the 'Choose Torsions in Currently Selected Residues...' button is pressed"""
        if not self.vf.flexDict.has_key('flex_residues'): 
            t='select residues first'
            self.vf.warningMsg(t)
            return 'ERROR'
        nodes = self.vf.flexDict['flex_residues']
        self.mol = nodes[0].top
        if not len(nodes):
            t='no residues selected'
            self.vf.warningMsg(t)
            return 'ERROR'
        if not nodes.__class__==ResidueSet:
            t='selection must be of type ResidueSet'
            self.vf.warningMsg(t)
            return 'ERROR'
        self.changeFlexResCt(0)
        if not hasattr(self, 'ifd'):
            self.torsStr= Tkinter.StringVar(master=self.vf.GUI.ROOT)
            s = 'Number of rotatable bonds ='+str(0)+ ' / '+str(MAXTORS)+'\n'
            self.torsStr.set(s)
            self.renameAromatic= Tkinter.IntVar(master=self.vf.GUI.ROOT)
            infoStr = 'Pick or drag-&-pick bonds.  \nGreen = rotatable, \nMagenta = non-rotatable, \nRed = unrotatable.\n\n'
            ifd = self.ifd=InputFormDescr(title='Torsion Count')
            ifd.append({'name': 'maxTorsLab',
                'widgetType':Tkinter.Label,
                'wcfg':{'text':infoStr},
                'gridcfg':{'sticky':Tkinter.W + Tkinter.E}}),
            ifd.append({'name':'torsEntryLab',
                'widgetType':Tkinter.Label,
                'wcfg':{'textvariable':self.torsStr},
                'gridcfg':{'sticky':Tkinter.W +Tkinter.E}}),
            ifd.append({'name': 'noAmideBut',
                'widgetType':Tkinter.Checkbutton,
                'wcfg':{'text':'amide torsions are allowed',
                    #'activebackground':'white',
                    'selectcolor':'white',
                    'indicatoron':0,
                    'command':self.setNoAmideTors_cb},
                'gridcfg':{'sticky':Tkinter.W+Tkinter.E, 'columnspan':2}}),
            ifd.append({'name': 'closeBut',
                'widgetType':Tkinter.Button,
                'wcfg':{'text':'Close','command':self.close_cb},
                'gridcfg':{'sticky':Tkinter.W+Tkinter.E, 'columnspan':2}})
            self.form= self.vf.getUserInput(self.ifd, modal=0, blocking=0)
        else:
            if hasattr(self,'form'):
                self.form.root.deiconify()
        ## !@@! Is THIS APPROPRIATE?!?!?
        if hasattr(self.vf.ADflex_setResidues,'torscount'):
            self.torscount = self.vf.ADflex_setResidues.torscount
        else:
            self.torscount = 0
        self.currentNodes=nodes.sideChain
        self.vf.displayLines(self.currentNodes, topCommand=0)
        self.vf.displayLines(self.currentNodes, only=1, topCommand=0, redraw=1)
        pTbndset = self.currentNodes.bonds[0].get(lambda x:x.activeTors==1)
        pTatomset = (pTbndset.atom1 + pTbndset.atom2).uniq()
        if len(pTatomset):
            geom = self.mol.geomContainer.geoms['AtomLabels']
            geom.Set(billboard=True, fontStyle='solid', fontScales=(.3,.3,.3,))
            self.vf.labelByProperty(pTatomset, ('name', ), topCommand=0,location='Last', redraw=1)
        self.buildCol(self.mol, self.torscount)
        self.pTatomset = pTatomset
        self.save = self.vf.ICmdCaller.commands.value[None]
        self.vf.setICOM(self, topCommand = 0)
        self.vf.setIcomLevel( Atom )
        self.vf.flexDict['torscount'] = self.torscount
        #need to set up display by torsion activity here
        

    def close_cb(self, event=None):
        #at this point, remove any residues which have no active torsions left
        dict = self.vf.flexDict
        flex_residues = dict['flex_residues']
        badList = []
        for item in flex_residues:
            if not hasattr(item,'torscount'):
                badList.append(item)
            if not item.torscount:
                print "eliminating", item.name
                badList.append(item)
        badSet = ResidueSet(badList)
        flex_residues = flex_residues - badSet
        dict['flex_residues'] = flex_residues
        dict['rigidResidues'] = dict['rigidResidues']+badSet
        dict['flex_residues_number'] = len(flex_residues)
        if len(badSet):
            badAtoms = badSet.atoms
            self.currentNodes = self.currentNodes-badAtoms
        self.changeFlexResCt(0)
        self.vf.displayLines(self.currentNodes, topCommand=0, only=1)
        self.vf.colorByAtomType(self.currentNodes)
        self.vf.labelByProperty(self.pTatomset, ('name', ), negate=True, topCommand=0, redraw=1)
        self.pTatomset = AtomSet()
        #set the PCOM to something else
        self.vf.setICOM(self.save, topCommand = 0)
        self.vf.GUI.VIEWER.Redraw()
        self.form.root.withdraw()
        self.save = None
        

    def setNoAmideTors_cb(self, event=None, log=1, redraw=0):
        self.setNoAmideTors()
        self.buildCol(self.mol, self.torscount)


    def setNoAmideTors(self, log=0):
        self.vf.flexDict['noAmides']=self.hasAmide
        if self.hasAmide:
            self.hasAmide = 0
            self.turnOffAmides()
            self.ifd.entryByName['noAmideBut']['widget'].config(text='amide torsions are not allowed')
        else:
            self.hasAmide = 1
            self.turnOnAmides()
            self.ifd.entryByName['noAmideBut']['widget'].config(text='amide torsions are allowed')


    def turnOffAmides(self):
        map(self.turnOffAmide, self.vf.flexDict['flex_residues'])
        

    def turnOnAmides(self):
        map(self.turnOnAmide, self.vf.flexDict['flex_residues'])


    def turnOnAmide(self, res):
        if not hasattr(res, 'amidebonds'):
            return
        for item in res.amidebonds:
            #only turn on the ones which were turned off
            if item.possibleTors and not item.activeTors:
                item.activeTors = 1
                self.torscount = self.torscount + 1


    def turnOffAmide(self, res):
        if not hasattr(res, 'amidebonds'):
            return
        for item in res.amidebonds:
            #only turn off the ones which were turned on
            if item.possibleTors and item.activeTors:
                item.activeTors = 0
                self.torscount = self.torscount - 1


    def buildCol(self, mol, torscount):
        #process residues
        s = 'Number of rotatable bonds ='+str(torscount) + ' / '+str(MAXTORS)+'\n'
        self.torsStr.set(s)
        currentbonds=mol.geomContainer.atoms['bonded'].bonds[0]
        col = []
        for b in currentbonds:
            if b.possibleTors:
                if b.activeTors: col.append((0,1,0))
                else: col.append((1,0,1))
            else:
                col.append((1,0,0))
        mol.geomContainer.geoms['bonded'].Set(materials=col,
                                              inheritMaterial=False,
                                              matBind=viewerConst.PER_PART)
        self.vf.GUI.VIEWER.Redraw()
        

    def changeFlexResCt(self,delta):
        n=self.vf.flexDict['flex_residues_number']
        n=n+delta
        self.vf.flexDict['flex_residues_number']=n
        msg = 'current %d flexible residues'%n
        self.vf.GUI.pickLabel.configure(text=msg)


    def checkAromatic_cb(self, event=None):
        if self.renameAromatic.get():
            self.ifd.entryByName['checkAromBut']['widget'].config(text='aromatic carbons named A..')
            map(self.nameA,self.vf.flexDict['flex_residues'])
        else:
            self.ifd.entryByName['checkAromBut']['widget'].config(text='aromatic carbons named C..')
            map(self.renameC,self.vf.flexDict['flex_residues'])


    def nameA(self, res):
        print 'changing the names of aromatic carbons is no longer supported'
        return


    def renameC(self,res):
        print 'changing the names of aromatic carbons is no longer supported'
        return


    def __call__(self, bonds,  **kw):
        kw['topCommand'] = 0
        kw['busyIdle'] = 1
        kw['log'] = 0
        #self.setUpDisplay()
        apply(self.doitWrapper, (bonds,), kw)


    def doit(self, bonds):
        if not self.currentNodes:
            msg='No residues currently processed'
            self.vf.warningMsg(msg)
            return 'ERROR'
        for bond in bonds:
            if not bond.possibleTors: continue
            atoms = AtomSet([bond.atom1, bond.atom2])
            self.vf.ADflex_setBondRotatableFlag(atoms, not bond.activeTors,
                topCommand = 0, redraw = 1, log =1, setupUndo = 1)


    def stop(self):
        self.done_cb()


    def getObjects(self,pick):
        flexMol = self.vf.flexDict['flex_residues'].top.uniq()[0]
        flexMolGeom = flexMol.geomContainer.geoms['bonded']
        for o, val in pick.hits.items(): #loop over geometries
            primInd = map(lambda x: x[0], val)
            if o != flexMolGeom: continue
            else: g = o.mol.geomContainer
            if g.geomPickToBonds.has_key(o.name):
                func = g.geomPickToBonds[o.name]
                if func: return func(o, primInd)
            else:
                l = []
                bonds = g.atoms[o.name].bonds[0]
                for i in range(len(primInd)):
                    l.append(bonds[int(primInd[i])])
                return BondSet(l)


    def done_cb(self):
        self.guiUp = 0
        if self.form: self.form.withdraw()
        self.vf.colorByAtomType(self.vf.flexDict['flex_residues'],
                    topCommand=0, redraw=1)


    def dismiss(self):
        self.vf.setICOM(self.save, topCommand=0)
        self.save = None
        self.done_cb()


AF_ProcessResiduesGUI = CommandGUI()



class AF_ProcessHingeResidues(SetRotatableBonds, MVBondICOM):
    """ allows user to process interactively a set of hinge residues in macromolecule whose sidechains are to be flexed in an autodock run"""


    def __init__(self, func=None):
        SetRotatableBonds.__init__(self)
        MVBondICOM.__init__(self)
        self.save = None
        self.guiUp = 0
        self.pickLevel = 'parts'


    def onAddCmdToViewer(self):
        if not hasattr(self.vf, 'flexDict'):
            self.vf.flexDict = {}
        if not hasattr(self.vf, 'colorByAtomType'):
            self.vf.loadCommand('colorCommands', 'colorByAtomType', 'Pmv')
        if not hasattr(self.vf, 'setICOM'):
            self.vf.loadCommand('interactiveCommands', 'setICOM', 'Pmv')
        if self.vf.hasGui and not self.torsStr:
            self.torsStr = Tkinter.StringVar(master=self.vf.GUI.ROOT)


    def setAutoFlexFieldsHinge(self, res, atomOne, atomTwo, atoms):
        if hasattr(res, 'processed'): 
            print res.full_name(), ' already has autoflexFields'
            return
        #PROCESS HINGE RESIDUES
        #backbone_names = ['CA', 'C','N','O','HN','HN1','HN2', 'HA', 
        backbone_names = ['C','N','O','HN','HN1','HN2', 'HA', 
                    'H1','H2','H3','HO', 'H']
        #only process the sidechains of atoms which are designated to be moved by hinge
        #resatoms = res.atoms - AtomSet([atomOne, atomTwo])
        resatoms = res.atoms.get(lambda x: x in atoms)
        sidechain = resatoms.get(lambda x: x.name not in backbone_names)
        res.sideChain = sidechain
        bondlist = res.bondlist = sidechain.bonds[0]
        bondlist.leaf = 0
        rbs = RotatableBondSelector()
        rotatables = rbs.select(bondlist)
        for b in rotatables:
            b.possibleTors = 1
            b.activeTors = 1
        amides = AmideBondSelector().select(bondlist)
        for b in amides:
            b.activeTors = 0
            b.possibleTors = 1
        guanidiniums = GuanidiniumBondSelector().select(bondlist)
        for b in guanidiniums:
            b.activeTors = 0
            b.possibleTors = 1
        leaves = LeafBondSelector().select(bondlist)
        for b in leaves:
            b.activeTors = 0
            b.possibleTors = 0
        res.torscount = len(bondlist.get(lambda x: x.activeTors==1))
        #this field is not used in AutoDock4
        res.torsdof = res.torscount
        #res.rootlist = AtomSet([resatoms[0]]) ??is this necessary?
        #res.rootlist = AtomSet([resatoms.get(lambda x: x._uniqIndex == 0)[0]])
        res.processed = True


    def changeHingeFlexResCt(self, delta):
        n = self.vf.flexDict.get('flex_residues_number',0)
        n = n + delta
        self.vf.flexDict['flex_residues_number'] = n
        msg = 'current %d flexible residues'%n
        self.vf.GUI.pickLabel.configure(text=msg)


    def guiCallback(self):
        """called INDIRECTLY each time the 'Set Selected Residues' button is pressed"""
        if not self.vf.flexDict.has_key('hinge_list'): 
            t='set hinge first'
            self.vf.warningMsg(t)
            return 'ERROR'
        h_list = self.vf.flexDict['hinge_list']
        if not len(h_list):
            t='currently hinge list is empty!'
            self.vf.warningMsg(t)
            return 'ERROR'
        if not hasattr(self, 'hinge'):
            self.hinge_to_process = self.vf.flexDict['hinge_list'][-1]
        hinge = self.hinge_to_process
        (atomOne,atomTwo), atoms = hinge
        self.mol = atoms[0].top
        #self.mol.allAtoms.bonds[0].possibleTors = 0
        #self.mol.allAtoms.bonds[0].activeTors = 0
        if not len(atoms):
            t='no hinge atoms specified'
            self.vf.warningMsg(t)
            return 'ERROR'
        self.currentNodes = atoms
        hinge_resSet = atoms.parent.uniq()
        for res in hinge_resSet:
            self.setAutoFlexFieldsHinge(res, atomOne, atomTwo, atoms)
        self.torscount = Numeric.add.reduce(hinge_resSet.torscount)
        self.changeHingeFlexResCt(len(hinge_resSet))
        if not hasattr(self, 'ifdX'):
            self.torsStr= Tkinter.StringVar(master=self.vf.GUI.ROOT)
            s = 'Number of rotatable bonds ='+str(0)+ ' / '+str(MAXTORS)+'\n'
            self.torsStr.set(s)
            self.renameAromatic= Tkinter.IntVar(master=self.vf.GUI.ROOT)
            infoStr = 'Pick or drag-&-pick bonds.  \nGreen = rotatable, \nMagenta = non-rotatable, \nRed = unrotatable.\n\n'
            ifdX = self.ifdX=InputFormDescr(title='Torsion Count')
            ifdX.append({'name': 'maxTorsLab',
                'widgetType':Tkinter.Label,
                'wcfg':{'text':infoStr},
                'gridcfg':{'sticky':Tkinter.W + Tkinter.E}}),
            ifdX.append({'name':'torsEntryLab',
                'widgetType':Tkinter.Label,
                'wcfg':{'textvariable':self.torsStr},
                'gridcfg':{'sticky':Tkinter.W +Tkinter.E}}),
            ifdX.append({'name': 'closeBut',
                'widgetType':Tkinter.Button,
                'wcfg':{'text':'Close','command':self.closeX_cb},
                'gridcfg':{'sticky':Tkinter.W+Tkinter.E, 'columnspan':2}})
            self.formX= self.vf.getUserInput(self.ifdX, modal=0, blocking=0)
        else:
            if hasattr(self,'formX'):
                self.formX.root.deiconify()
        self.save = self.vf.ICmdCaller.commands.value[None]
        self.vf.setICOM(self, topCommand = 0)
        self.vf.setIcomLevel( Atom )
####    #need to set up display by torsion activity here
        #could have preexisting flexible residues... or not
        self.vf.flexDict.setdefault('torscount', 0)
        if 'flex_residues' in self.vf.flexDict.keys() and \
                len(self.vf.flexDict['flex_residues']):
            self.vf.flexDict['torscount'] += self.torscount
        else:
            self.vf.flexDict['torscount'] = self.torscount
        self.vf.displayLines(self.currentNodes, only=1, topCommand=0, redraw=1)
        self.buildCol(self.mol, self.torscount)


    def closeX_cb(self, event=None):
        dict = self.vf.flexDict
        self.vf.displayLines(self.currentNodes, topCommand=0, only=1)
        self.vf.colorByAtomType(self.currentNodes)
        self.vf.GUI.VIEWER.Redraw()
        self.formX.root.withdraw()
        self.vf.flexDict.setdefault('torscount', 0)
        #this is dangerous
        #!@@! fix this: what if two hinges... or if remove one... or all
        flex_residues_torscount = 0
        if hasattr(self.vf.ADflex_setResidues, 'torscount') and \
                    self.vf.ADflex_setResidues.torscount>0:
            flex_residues_torscount = self.vf.ADflex_setResidues.torscount
        self.vf.flexDict['torscount'] = self.torscount + flex_residues_torscount
        self.vf.setICOM(self.save, topCommand = 0)
        self.save = None
        self.dismiss()
        

    def buildCol(self, mol, torscount):
        #process hinge
        s = 'Number of rotatable bonds ='+str(torscount) + ' / '+str(MAXTORS)+'\n'
        self.torsStr.set(s)
        currentbonds = mol.geomContainer.atoms['bonded'].bonds[0]
        col = []
        for b in currentbonds:
            if hasattr(b, 'possibleTors') and b.possibleTors:
                if hasattr(b, 'activeTors') and b.activeTors: col.append((0,1,0))
                else: col.append((1,0,1))
            else:
                col.append((1,0,0))
        mol.geomContainer.geoms['bonded'].Set(materials=col,
                                              inheritMaterial=False,
                                              matBind=viewerConst.PER_PART)
        self.vf.GUI.VIEWER.Redraw()
        

    def changeFlexResCt(self, delta):
        self.vf.flexDict.setdefault('flex_residues_number',0)
        n = self.vf.flexDict['flex_residues_number']
        n = n + delta
        self.vf.flexDict['flex_residues_number'] = n
        msg = 'current %d flexible residues'%n
        self.vf.GUI.pickLabel.configure(text=msg)


    def __call__(self, bonds,  **kw):
        kw['topCommand'] = 0
        kw['busyIdle'] = 1
        kw['log'] = 0
        #kw['flexRes'] = False
        #self.setUpDisplay()
        apply(self.doitWrapper, (bonds,), kw)


    def doit(self, bonds, **kw):
        if not self.currentNodes:
            msg='No residues currently processed'
            self.vf.warningMsg(msg)
            return 'ERROR'
        for bond in bonds: 
            if not bond.possibleTors: 
                continue
            atoms = AtomSet([bond.atom1, bond.atom2])
            self.vf.ADflex_setBondRotatableFlag(atoms, not bond.activeTors,
                topCommand = 0, redraw = 1, log =1, setupUndo = 1,
                flexRes=True)


    def stop(self):
        self.done_cb()


    def getObjects(self,pick):
        flexMol = self.mol
        flexMolGeom = flexMol.geomContainer.geoms['bonded']
        for o, val in pick.hits.items(): #loop over geometries
            primInd = map(lambda x: x[0], val)
            if o != flexMolGeom: continue
            else: g = o.mol.geomContainer
            if g.geomPickToBonds.has_key(o.name):
                func = g.geomPickToBonds[o.name]
                if func: return func(o, primInd)
            else:
                l = []
                bonds = g.atoms[o.name].bonds[0]
                for i in range(len(primInd)):
                    l.append(bonds[int(primInd[i])])
                return BondSet(l)


    def done_cb(self):
        self.guiUp = 0
        if self.form: self.form.withdraw()
        #self.vf.ADflex_setHinge.spheres.visible=False
        self.vf.GUI.VIEWER.Redraw()


    def dismiss(self):
        self.vf.setICOM(self.save, topCommand=0)
        self.save = None
        self.done_cb()

AF_ProcessHingeResiduesGUI = CommandGUI()



class AF_EditHinge(MVCommand, MVAtomICOM):
    """ allows user to remove atoms to be moved from a hinge interactively"""


    def __init__(self, func=None):
        MVCommand.__init__(self, func)
        MVAtomICOM.__init__(self)
        self.save = None
        self.mode = None
        self.atomList = AtomSet([])
        #self.undoAtList = AtomSet([])


    def onAddCmdToViewer(self):
        if not hasattr(self.vf, 'flexDict'):
            self.vf.flexDict = {}
        if not hasattr(self.vf, 'colorByAtomType'):
            self.vf.loadCommand('colorCommands', 'colorByAtomType', 'Pmv')
        if not hasattr(self.vf, 'setICOM'):
            self.vf.loadCommand('interactiveCommands', 'setICOM', 'Pmv')
        if not hasattr(self.vf, 'ADflex_setHinge'):
            self.vf.loadCommand('autoflexCommands', 'ADflex_setHinge', 'AutoDockTools')
        if self.vf.hasGui:
            self.mode = Tkinter.StringVar(master=self.vf.GUI.ROOT)
    

    def continuousUpdate_cb(self, name, oldval, newval):
        if newval == 'yes':
            self.continuousUpdate = 1
            for event in ['<B2-Motion>', '<B3-Motion>', '<Shift-B3-Motion>']:
                self.vf.GUI.addCameraCallback(event, self.update_cb)
        else:
            self.continuousUpdate = 0
            for event in ['<B2-Motion>', '<B3-Motion>', '<Shift-B3-Motion>']:
                self.vf.GUI.removeCameraCallback(event, self.update_cb)
                

    def update_cb(self, event=None):
        #print "in update_cb"
        self.update()


    def update(self, forward=1, event=None):
        #print 'in update: self.atomList=', self.atomList
        if not len(self.atomList):
            print 'Currently no atoms: resetting geoms'
            return
        for at in self.atomList:
            c1 = self.getTransformedCoords(at)
            self.vf.ADflex_setHinge.lineVertices.append(tuple(c1))
        self.vf.ADflex_setHinge.update()


    def guiCallback(self):
        """  """
        if not hasattr(self,'form'):
            ifd = self.ifd = InputFormDescr(title="Edit Atoms to be moved by current Hinge by:")
            #ifd.append({'name': 'testLab',
            #    'widgetType':Tkinter.Label,
            #    'wcfg':{'text':'Mode:'},
            #    'gridcfg':{'columnspan':3,'sticky':'w'}})
            ifd.append({'widgetType':Tkinter.Radiobutton,
                'wcfg':{'text':'adding Atoms',
                    'variable':self.mode,
                    'value': 'add'},
                'gridcfg':{'sticky':'we'}})
            ifd.append({'widgetType':Tkinter.Radiobutton,
                'wcfg':{'text':'removing Atoms', 
                    'variable':self.mode,
                    'value': 'remove'},
                'gridcfg':{'row':-1,'sticky':'we'}})
            #ifd.append({'widgetType':Tkinter.Button,
            #    'wcfg':{'text':'Close', 'command':self.close_cb},
            #    'gridcfg':{'sticky':'we', 'columnspan':2}})
            self.form = self.vf.getUserInput(ifd, modal=0, blocking=0)
            self.form.root.protocol('WM_DELETE_WINDOW',self.close_cb)
        else:
            self.form.root.deiconify()
        self.save = self.vf.ICmdCaller.commands.value[None]
        self.vf.setICOM(self)
        

    def close_cb(self, event=None):
        #at this point, just cleanup
        self.stop()
        self.form.withdraw()
        self.vf.GUI.VIEWER.Redraw()
        self.mode.set("")
        #update hinge_atoms and non_hinge_atoms here?
        flexDict = self.vf.flexDict
        all = 'all_hinge_atoms'
        opp = 'non_hinge_atoms'
        #print 'before call to getAllHingeAtoms, ha=', len(flexDict[all]),
        #print 'nha=', len(flexDict[opp]),
        flexDict = self.vf.flexDict
        flexDict[all], flexDict[opp] = self.vf.ADflex_setHinge.getAllHingeAtoms()
        #print '2:after call to getAllHingeAtoms, ha=', len(flexDict.get(all,[])),
        #print '2:nha=', len(flexDict.get(opp, []))



    def doit(self, atoms):
        if not hasattr(self.vf, 'flexDict'):
            self.warningMsg('no flexDIct!')
            return "ERROR"
        dict = self.vf.flexDict
        if not 'hinge_list' in dict.keys():
            self.warningMsg('no hinge_list in flexDict!')
            return "ERROR"
        hinges = self.vf.flexDict['hinge_list']
        if not len(hinges):
            self.warningMsg('current hinges in hinge_list !')
            return "ERROR"
        if len(hinges):
            current_hinge = hinges[-1]
            if not len(current_hinge)==2:
                self.warningMsg('last hinge in hinge_list is ill-formed!')
                return 'ERROR'
        atoms_to_move = current_hinge[1]
        if not len(atoms_to_move):
            self.warningMsg("no atoms to move in current hinge!")
            return "ERROR"
        if self.mode.get()=='remove':
            new_atoms_to_move = atoms_to_move - atoms
            #remember to put these atoms into rigid 
        else:
            new_atoms_to_move = atoms_to_move + atoms
            #remember to remove these atoms from rigid 
        #print 'setting atoms to move to ', new_atoms_to_move
        new_ats = new_atoms_to_move.uniq()
        self.vf.flexDict['hinge_list'][-1][1] = new_ats
        self.vf.ADflex_setHinge.atoms = new_ats
        if self.vf.hasGui:
            self.vf.ADflex_setHinge.hingeAtoms.Set(vertices = new_ats.coords)
            self.vf.ADflex_setHinge.hingeAtoms.visible = 1
            self.vf.GUI.VIEWER.Redraw()


    def startICOM(self):
        self.vf.setIcomLevel( Atom )


    def dismiss(self):
        self.vf.setICOM(self.save, topCommand=0)
        self.save = None
        self.vf.GUI.VIEWER.Redraw()
        self.done_cb()

    def stop(self):
        self.dismiss()

    def done_cb(self):
        self.guiUp = 0
        if self.form: self.form.withdraw()



AF_EditHingeGUI = CommandGUI()
AF_EditHingeGUI.addMenuCommand('AutoToolsBar', menuText['AutoFlexMB'],\
        menuText['Edit Hinge'])




class AF_SetHinge(MVCommand, MVAtomICOM):
    """ allows user to setup a hinge interactively to be flexed in an autodock run"""


    def __init__(self, func=None):
        MVCommand.__init__(self, func)
        MVAtomICOM.__init__(self)
        self.save = None
        self.guiUp = 0
        self.atoms = []
        self.atomList = []
        self.atomCenters = []
        self.atomOne = None   #hinge pt 1
        self.atomTwo = None   #hinge pt 2
        self.coordSlot = None
        self.labelStrs = []
        self.labelCenters = []
        self.lineVertices = []
        self.snakeLength = 2
        self.old_hinge_method = ""
        self.old_atoms_method = ""
        #self.pickLevel = 'parts'


    def getTransformedCoords(self, atom):
        # when there is no viewer, the geomContainer is None
        if not atom.top.geomContainer:
            return atom.coords
        g = atom.top.geomContainer.geoms['master']
        c = self.vf.transformedCoordinatesWithInstances(AtomSet([atom]))
        return  Numeric.array(c[0], 'f')


    #def setupUndoBefore(self, ats):
    #    self.addUndoCall((),{}, self.name+'.undo')


    #def undo(self,*args, **kw):
    #    if len(self.lineVertices):
    #        self.atomList = self.atomList[:-1]
    #        self.update(forward=0)


    def onAddCmdToViewer(self):
        if not hasattr(self.vf, 'flexDict'):
            self.vf.flexDict = {}
        if not hasattr(self.vf, 'colorByAtomType'):
            self.vf.loadCommand('colorCommands', 'colorByAtomType', 'Pmv')
        if not hasattr(self.vf, 'setICOM'):
            self.vf.loadCommand('interactiveCommands', 'setICOM', 'Pmv')
        if self.vf.hasGui:
            from DejaVu.Spheres import Spheres
            from DejaVu.Geom import Geom
            miscGeom = self.vf.GUI.miscGeom
            self.masterGeom = Geom("setAutoFlexHingeGeoms", 
                                    shape=(0,0), protected=True)
            self.vf.GUI.VIEWER.AddObject(self.masterGeom, parent=miscGeom)
            self.spheres = Spheres(name='AutoFlexHinge_spheres', 
                materials=((1.,1.,0),), shape=(0,3), radii=0.2, 
                quality=15, inheritMaterial=0, protected=True)
            from opengltk.OpenGL import GL
            self.spheres.frontPolyMode = GL.GL_LINE
            #marker for end points of hinge
            self.spheres.Set(visible=1, tagModified=False)
            self.spheres.pickable = 0
            from DejaVu.IndexedPolylines import IndexedPolylines
            #marker for axis of hinge
            self.line = IndexedPolylines('distLine', materials = ((1,1,0),),
                                      inheritMaterial=0, lineWidth=3, 
                                      stippleLines=1, protected=True,
                                      visible=True)
            self.vf.GUI.VIEWER.AddObject(self.line, parent=self.masterGeom)
            self.line.pickable = 0
            from DejaVu.Points import CrossSet
            #markers for atoms moved by hinge motion
            self.hingeAtoms = CrossSet('AutoFlexHinge_hingeAtoms', 
                inheritMaterial=0, materials=((1.,0.3,0),),
                offset=0.1,lineWidth=2, protected=True, visible=True)
            self.hingeAtoms.Set(visible=0, tagModified=False)
            self.hingeAtoms.pickable = 0
            self.selLevel = Tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.selectionBase = Tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.selectHingeMethod = Tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.selectionCenter = Tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.adjustPt = Tkinter.IntVar(master=self.vf.GUI.ROOT)
            self.adjustPt.set(0)
            self.editAtomsToMove = Tkinter.IntVar(master=self.vf.GUI.ROOT)
            self.editAtomsToMove.set(0)
            self.vf.GUI.VIEWER.AddObject(self.spheres, redo=0,
                                parent=self.masterGeom)
            self.vf.GUI.VIEWER.AddObject(self.hingeAtoms, redo=0,
                                parent=self.masterGeom)
    

    def continuousUpdate_cb(self, name, oldval, newval):
        if newval == 'yes':
            self.continuousUpdate = 1
            for event in ['<B2-Motion>', '<B3-Motion>', '<Shift-B3-Motion>']:
                self.vf.GUI.addCameraCallback(event, self.update_cb)
        else:
            self.continuousUpdate = 0
            for event in ['<B2-Motion>', '<B3-Motion>', '<Shift-B3-Motion>']:
                self.vf.GUI.removeCameraCallback(event, self.update_cb)
                

    def update_cb(self, event=None):
        if not len(self.atomList):
            return
        vi = self.vf.GUI.VIEWER
        if vi.redirectTransformToRoot:
            return
        if vi.currentObject==vi.rootObject:
            return 
        self.update()


    def update(self, forward=1, event=None):
        if not len(self.atomList):
            print 'Currently no atoms: resetting geoms'
            self.spheres.Set(vertices=[])
            self.labels.Set(vertices=[])
            self.line.Set(vertices=[])
            #self.hingeAtoms.Set(vertices=[])  #???
            return
        self.lineVertices=[]
        for at in self.atomList:
            c1 = self.getTransformedCoords(at)
            self.lineVertices.append(tuple(c1))
        self.spheres.Set(vertices=self.lineVertices)
        if len(self.lineVertices)==1:
            #this fixes case of stepping back over 1st label
            self.labels.Set(vertices=[])
            self.line.Set(vertices=[])
        elif len(self.atomList)>1:
            self.line.Set(vertices=self.lineVertices, faces=[(0,1),] )
            self.atomOne = self.atomList[0]
            if self.coordSlot==None:
                mol = self.atomOne.top
                self.coordSlot = len(mol.allAtoms[0]._coords)
                mol.allAtoms.addConformation(mol.allAtoms.coords[:])
                mol.allAtoms.setConformation(self.coordSlot)
            self.atomTwo = self.atomList[1]
            if self.atomOne==self.atomTwo:
                self.vf.warningMsg('ERROR:two hinge atoms are identical')
                self.clear_cb()
                return
        #setting spheres doesn't trigger redraw so do it explicitly
        self.vf.GUI.VIEWER.Redraw()


    def buildForm(self):
        #build ifd and form
        ifd = self.ifd = InputFormDescr(title='Set Up Hinge')
        ifd.append({'widgetType': Tkinter.Label,
                    'name':'set hinge atoms label',
                    'wcfg':{'text':'Set Hinge Atoms by:',
                            'bd':2,'relief':'ridge'},
                    'gridcfg':{'sticky':Tkinter.W+Tkinter.E,
                               'columnspan':8}})
        ifd.append({'widgetType':Tkinter.Radiobutton,
                    'name':'pickingRB',
                    'wcfg':{'variable':self.selectHingeMethod,
                            'text':'picking',
                            'value':'picking',
                            'command':self.updateHinge},  #1
                    'gridcfg':{'sticky':Tkinter.W}}),
        ifd.append({'widgetType':Tkinter.Radiobutton,
                    'name':'selectionRB',
                    'wcfg':{'variable':self.selectHingeMethod,
                            'text':'current selection',
                            'value':'cursel',
                            'command':self.updateHinge},  #1
                    'gridcfg':{'sticky':Tkinter.W, 'row':-1, 'column':3}}),
#        ifd.append({'widgetType':Tkinter.Checkbutton,
#                    'name':'pt1CB',
#                    'wcfg':{ 'text':'adjust  x y z coords of hingePts ',
#                            'command':self.adjustPt_cb},
#                    'gridcfg':{'sticky':Tkinter.W}}),
#        ifd.append({'widgetType':Tkinter.Checkbutton,
#                    'name':'editLastHingeCB',
#                    'wcfg':{'variable':self.editAtomsToMove,
#                            'text':'edit atoms to move',
#                            'command':self.editAtomsToMove_cb},
#                    'gridcfg':{'sticky':Tkinter.W}}),
##                    'gridcfg':{'sticky':Tkinter.W, 'row':-1, 'column':3}}),
#        ifd.append({'name':'xval_tw',
#            'widgetType':ThumbWheel,
#            'wType':ThumbWheel,
#            'wcfg':{ 'labCfg':{ 
#                        'text': 'x',
#                    },
#                    'type':'float',
#                    'precision':2,
#                    'width':30,
#                    'continuous':1,
#                    'wheelPad':2,
#                    'height':20,
#                    'value':0.,
#                    'callback':self.updatePt1,
#                    'oneTurn':50.,},
#             'gridcfg':{'columnspan':1,'sticky':'e'}})
#        ifd.append({'name':'yval_tw',
#            'widgetType':ThumbWheel,
#            'wType':ThumbWheel,
#            'wcfg':{ 'labCfg':{ 
#                        'text': 'y',
#                    },
#                    'type':'float',
#                    'precision':2,
#                    'width':30,
#                    'continuous':1,
#                    'wheelPad':2,
#                    'height':20,
#                    'value':0.,
#                    'callback':self.updatePt1,
#                    'oneTurn':50.,},
#             'gridcfg':{'columnspan':1,'sticky':'e', 'row':-1, 'column':1}})
#        ifd.append({'name':'zval_tw',
#            'widgetType':ThumbWheel,
#            'wType':ThumbWheel,
#            'wcfg':{ 'labCfg':{ 
#                        'text': 'z',
#                    },
#                    'type':'float',
#                    'precision':2,
#                    'width':30,
#                    'continuous':1,
#                    'wheelPad':2,
#                    'height':20,
#                    'value':0.,
#                    'callback':self.updatePt1,
#                    'oneTurn':50.,},
#             'gridcfg':{'columnspan':1,'sticky':'e', 'row':-1, 'column':2}})
#        ifd.append({'name':'xval2_tw',
#            'widgetType':ThumbWheel,
#            'wType':ThumbWheel,
#            'wcfg':{ 'labCfg':{ 
#                        'text': 'x',
#                    },
#                    'type':'float',
#                    'precision':2,
#                    'width':30,
#                    'continuous':1,
#                    'wheelPad':2,
#                    'height':20,
#                    'value':0.,
#                    'callback':self.updatePt2,
#                    'oneTurn':50.,},
#             'gridcfg':{'columnspan':1,'sticky':'e', 'row':-1, 'column':3}})
#        ifd.append({'name':'yval2_tw',
#            'widgetType':ThumbWheel,
#            'wType':ThumbWheel,
#            'wcfg':{ 'labCfg':{ 
#                        'text': 'y',
#                    },
#                    'type':'float',
#                    'precision':2,
#                    'width':30,
#                    'continuous':1,
#                    'wheelPad':2,
#                    'height':20,
#                    'value':0.,
#                    'callback':self.updatePt2,
#                    'oneTurn':50.,},
#             'gridcfg':{'columnspan':1,'sticky':'e', 'row':-1, 'column':4}})
#        ifd.append({'name':'zval2_tw',
#            'widgetType':ThumbWheel,
#            'wType':ThumbWheel,
#            'wcfg':{ 'labCfg':{ 
#                        'text': 'z',
#                    },
#                    'type':'float',
#                    'precision':2,
#                    'width':30,
#                    'continuous':1,
#                    'wheelPad':2,
#                    'height':20,
#                    'value':0.,
#                    'callback':self.updatePt2,
#                    'oneTurn':50.,},
#             'gridcfg':{'columnspan':1,'sticky':'e', 'row':-1, 'column':5}})
        ifd.append({'name': 'atoms to move label',
                    'widgetType': Tkinter.Label,
                    'wcfg':{'bd':2,'relief':'ridge',
                            'text':'Atoms to move:'},
                    'gridcfg':{'sticky':Tkinter.W+Tkinter.E,
                               'columnspan':8}})
        ifd.append({'widgetType':Tkinter.Radiobutton,
                    'name':'betweenRB',
                    'wcfg':{'variable':self.selectionBase,
                            'text':'between hinge points',
                            'value':'between',
                            'command':self.updateBase},
                    'gridcfg':{'sticky':Tkinter.W, 'columnspan':2}})
        ifd.append({'widgetType':Tkinter.Radiobutton,
                    'name':'selectionRB',
                    'wcfg':{'text':'current selection',
                            'variable':self.selectionBase,
                            'value':'selection',
                            'command':self.updateBase},
                    'gridcfg':{'sticky':Tkinter.W, 'row':-1, 'column':3}})
        ifd.append({'widgetType':Tkinter.Radiobutton,
                    'name':'savedSetRB',
                    'wcfg':{'text':'a saved set',
                            'variable':self.selectionBase,
                            'value':'set',
                            'command':self.updateBase},
                    'gridcfg':{'sticky':Tkinter.W, 'row':-1, 'column':5}})
                    #'gridcfg':{'sticky':Tkinter.W}})
        ifd.append({'widgetType':Tkinter.Checkbutton,
                    'name':'editLastHingeCB',
                    'wcfg':{'variable':self.editAtomsToMove,
                            'text':'edit atoms to move',
                            'command':self.editAtomsToMove_cb},
                    'gridcfg':{'sticky':Tkinter.W}})
                    #'gridcfg':{'sticky':Tkinter.W, 'row':-1, 'column':2}}),
        ifd.append({'name':'selectB',
                    'widgetType':Tkinter.Button,
                    'wcfg':{'text':'Set Hinge',
                            'command':self.setHinge_cb},
                    'gridcfg':{'sticky':Tkinter.W+Tkinter.E, 'columnspan':2}})
        ifd.append({'name':'selectTorsionsB',
                    'widgetType':Tkinter.Button,
                    'wcfg':{'text':'Select Torsions in Last Hinge',
                            'command':self.selectTorsions_cb},
                    'gridcfg':{'sticky':Tkinter.W+Tkinter.E,'row':-1, 'column':2, 'columnspan':2}})
        ifd.append({'name':'clearLastB',
                    'widgetType':Tkinter.Button,
                    'wcfg':{'text':'Remove Last Hinge',
                            'command':self.removeLastHinge_cb},
                    'gridcfg':{'sticky':Tkinter.W+Tkinter.E, 'row':-1,'columnspan':2, 'column':4 }})
        ifd.append({'name':'clearB',
                    'widgetType':Tkinter.Button,
                    'wcfg':{'text':'Clear',
                            'command':self.clear_cb},
                    'gridcfg':{'sticky':Tkinter.W+Tkinter.E, 'columnspan':2}})
        ifd.append({'name':'clearAllB',
                    'widgetType':Tkinter.Button,
                    'wcfg':{'text':'Clear Hinge List',
                            'command':self.clearAll_cb},
                    'gridcfg':{'sticky':Tkinter.W+Tkinter.E,
                               'row':-1, 'column':2, 'columnspan':2}})
        ifd.append({'name':'closeB',
                    'widgetType':Tkinter.Button,
                    'wcfg':{'text':'Close',
                            'command':self.close_cb},
                    'gridcfg':{'sticky':Tkinter.W+Tkinter.E,
                               'row':-1,'column':4, 'columnspan':2}})

        self.form = self.vf.getUserInput(ifd, modal=0,blocking=0)
        self.form.root.protocol('WM_DELETE_WINDOW',self.close_cb)
        #self.selection_base = self.ifd.entryByName['atoms to move']['widget']
        #self.lb = self.ifd.entryByName['baseMols']['widget'].lb
        #self.lb.config({'selectmode':'multiple'})
        #self.lb.bind('<Enter>',self.updateBase)
        #self.lb.bind('<Leave>',self.updateBase)
#        self.xval = self.ifd.entryByName['xval_tw']['widget']
#        self.yval = self.ifd.entryByName['yval_tw']['widget']
#        self.zval = self.ifd.entryByName['zval_tw']['widget']
#        self.xval2 = self.ifd.entryByName['xval2_tw']['widget']
#        self.yval2 = self.ifd.entryByName['yval2_tw']['widget']
#        self.zval2 = self.ifd.entryByName['zval2_tw']['widget']
#        self.hide_xyz()
        self.old_val = ""


    def adjustPt_cb(self, event=None):
        if self.adjustPt.get():
            self.repack_xyz()
        else:
            self.hide_xyz()


    def editAtomsToMove_cb(self, event=None):
        flexDict = self.vf.flexDict
        if 'hinge_list' not in flexDict.keys() or len(self.vf.flexDict['hinge_list'])==0:
            msg = 'no hinges to edit'
            self.warningMsg(msg)
            self.editAtomsToMove.set(0)
            return
        if self.editAtomsToMove.get():
            self.vf.ADflex_editHinge.guiCallback()
        else:
            self.vf.ADflex_editHinge.close_cb()
        #reset the flexDict entries here
        all = 'all_hinge_atoms'
        opp = 'non_hinge_atoms'
        flexDict[all], flexDict[opp] = self.getAllHingeAtoms()


#    def updatePt1(self, event=None):
#        print 'setting atomOne.x to', self.xval.get()
#        self.atomOne.coords[0] = self.xval.get()
#        print 'setting atomOne.y to', self.yval.get()
#        self.atomOne.coords[1] = self.yval.get()
#        print 'setting atomOne.z to', self.zval.get()
#        self.atomOne.coords[2] = self.zval.get()
#        #self.lineVertices.append(tuple(c1))
#        self.spheres.Set(vertices=[self.atomOne.coords,])


#    def updatePt2(self, event=None):
#        print 'in updatePt2'
#        print 'setting atomTwo.x to', self.xval2.get()
#        self.atomTwo.coords[0] = self.xval2.get()
#        print 'setting atomTwo.y to', self.yval2.get()
#        self.atomTwo.coords[1] = self.yval2.get()
#        print 'setting atomTwo.z to', self.zval2.get()
#        self.atomTwo.coords[2] = self.zval2.get()
#        self.spheres.Set(vertices=[self.atomOne.coords,])



#    def update_xyz(self, val):
#        if self.old_val=="":
#            if val=='xyz':
#                self.repack_xyz()
#        elif val!=self.old_val:
#            if self.old_val=='xyz':
#                self.hide_xyz()
#            elif val=='xyz': 
#                self.repack_xyz()


    def getReverseSubtree(self, at1, at2, debug=False):
        m1 = at1.top
        m2 = at2.top
        assert m1==m2
        m1.allAtoms._subtree_selector = 0
        for b in at1.bonds:
            a2 = b.neighborAtom(at1)
            if a2.number<at1.number:
                a2._subtree_selector = 1
        at1._subtree_selector = at2._subtree_selector = 1
        non_moving_ats = m1.allAtoms.get(lambda x: x._subtree_selector==1)
        m1._MarkTree(at1)
        non_moving_ats._subtree_selector = 0
        results = m1.allAtoms.get(lambda x: x._subtree_selector==1)
        if not debug:
            delattr(m1.allAtoms, '_subtree_selector')
        return results
        

    def updateBase(self, val=None):
        base =  self.selectionBase.get()
        if base == 'between':
            if self.atomOne is not None and self.atomTwo is not None:
                #print 'calling getReverseSubtree with ', self.atomOne.full_name(), '+', self.atomTwo.full_name()
                self.atoms = self.getReverseSubtree(self.atomOne, self.atomTwo)
                self.hingeAtoms.Set(vertices = self.atoms.coords)
                self.hingeAtoms.visible = 1
        elif base=='selection':
            self.atoms = self.vf.getSelection()[:].findType(Atom)
            self.hingeAtoms.Set(vertices = self.atoms.coords)
            self.hingeAtoms.visible = 1
        elif base=='set':
            self.atoms = self.getSet()
            if self.atoms:
                self.hingeAtoms.Set(vertices = self.atoms.coords)
                self.hingeAtoms.visible = 1
        self.vf.GUI.VIEWER.Redraw()


    def getSet(self, event=None):
        idf = InputFormDescr(title ='')
        entries = []
        names = self.vf.sets.keys()
        names.sort()
        for key, value in self.vf.sets.items():
            entries.append( (key, value.comments) )
        idf.append({'name':'set',
                    'widgetType':ListChooser,
                    'wcfg':{'mode':'single',
                            'entries': entries,
                            'title':'Choose an item: '}})
        vals = self.vf.getUserInput(idf)
        if len(vals['set'])> 0:
            setNames = vals['set']
            print 'returning ', setNames
            return self.vf.sets[setNames[0]]
 

    def selectTorsions_cb(self, event=None):
        hinge_list = self.vf.flexDict['hinge_list']
        if len(hinge_list):
            hinges = []
            for item in hinge_list:
                entry = str(item[0])[1:-1] # + AtomSet(item[1])?
                hinges.append(entry)
            s = 'Select Hinge to process '
            ifd = InputFormDescr(title=s)
            ifd.append({'name':'hingeLC',
                'widgetType': 'ListChooser',
                'mode': 'single',
                'entries':hinges,
                'title': s,
                'lbwcfg':{'height':20,'selectforeground':'red','exportselection':0},
                'gridcfg':{'sticky':Tkinter.W +Tkinter.E}}),
            vals= self.vf.getUserInput(ifd, modal=1, blocking=1)
            if vals:
                try:
                    hinge = vals['hingeLC'][0]
                    index = hinges.index(hinge) 
                    self.hinge_to_process = hinge_list[index]
                except:
                    msg = "Unable to process hinge designated: ", hinge
                    self.warningMsg(msg)
                    return 'ERROR'
            last_hinge = self.vf.flexDict['hinge_list'][-1]
            if len(last_hinge)==2 and len(last_hinge[0])==2:
                self.vf.ADflex_processHingeResidues.guiCallback()
            else:
                self.warningMsg("last hinge is ill-formed")
                return
        else:
            self.warningMsg("currently there are no specified hinges")
            return 
        #if not len(self.atoms):
        #    print "currently no hinge atoms have been designated"
        #    return
        #self.vf.ADflex_processHingeResidues.guiCallback()


    def clear_cb(self, event=None):
        #callback for clear button
        self.selectionBase.set("")
        self.atomOne = None
        self.atomTwo = None
        if len(self.atoms):
            resSet = self.atoms.parent.uniq()
            tors_ct = 0
            for r in resSet:
                if hasattr(r, 'torscount'):
                    tors_ct += r.torscount
            #remove these torsions from self.torscount
            if hasattr(self, 'torscount'):
                self.torscount = self.torscount - tors_ct
        self.atoms = []
        self.lineVertices = []
        self.hingeAtoms.Set(vertices=[])
        self.spheres.Set(vertices=[])
        self.line.Set(vertices=[])
        #self.selection_base.getcurselection()
        self.selectionBase.set("")
        self.selectHingeMethod.set("")
        self.vf.GUI.VIEWER.Redraw()


    def getAllHingeAtoms(self):
        #rebuild hinge_residues
        d = {}
        for h in self.vf.flexDict['hinge_list']:
            d[h[0][0]] = 1
            d[h[0][1]] = 1
            for at in h[1]: d[at] = 1
            all_hinge_atoms = AtomSet(d.keys())
            resSet = all_hinge_atoms.parent.uniq()
            non_hinge_atoms = resSet.atoms - all_hinge_atoms
        #caution:  should this reorder the residues?
        #resSet.sort()
        return all_hinge_atoms, non_hinge_atoms


    def removeLastHinge_cb(self, event=None):
        flexDict = self.vf.flexDict
        if len(flexDict['hinge_list'])>1:
            flexDict['hinge_list'] = flexDict['hinge_list'][:-1]
            flexDict['all_hinge_atoms'], flexDict['non_hinge_atoms'] = self.getAllHingeAtoms()
        else:
            flexDict['hinge_list']=[]
            flexDict['all_hinge_atoms'] = AtomSet()
            flexDict['non_hinge_atoms'] = AtomSet()
        self.clear_cb()


    def clearAll_cb(self, event=None):
        self.vf.flexDict['hinge_list']=[]
        #self.vf.flexDict['hinge_residues'] = ResidueSet()
        self.vf.flexDict['all_hinge_atoms'] = AtomSet()
        self.vf.flexDict['non_hinge_atoms'] = AtomSet()
        self.torscount = 0
        #reset gui components
        self.clear_cb()


    def hide_xyz(self, event=None):
        for b in [self.xval, self.yval, self.zval,
                    self.xval2, self.yval2, self.zval2]:
            b.grid_forget()


    def repack_xyz(self, event=None):
        cfgstrings = ['xval_tw', 'yval_tw', 'zval_tw', 'xval2_tw', 'yval2_tw', 'zval2_tw']
        buttons = [self.xval, self.yval, self.zval, self.xval2, self.yval2, self.zval2]
        if hasattr(self, 'atomOne'):
            self.xval.set(self.atomOne.coords[0])
            self.yval.set(self.atomOne.coords[1])
            self.zval.set(self.atomOne.coords[2])
        if hasattr(self, 'atomTwo'):
            self.xval2.set(self.atomTwo.coords[0])
            self.yval2.set(self.atomTwo.coords[1])
            self.zval2.set(self.atomTwo.coords[2])
        for cfg_name, b in zip(cfgstrings, buttons):
            b.grid(self.ifd.entryByName[cfg_name]['gridcfg'])


    def updateHinge(self, event=None):
        #callback for changes to radiobuttons for hinge
        val = self.selectHingeMethod.get()
        #if different method for setting hinge, clear geoms
        if val!=self.old_hinge_method:
            self.clear_cb()
            self.old_hinge_method = val
        #here manage picking an atom
        #self.update_xyz(val)
        if hasattr(self, 'ap'):  # if there is an atom picker delete it
            self.ap.stop()
            del self.ap
        if val == 'picking':
            from Pmv.picker import AtomPicker
            self.ap=AtomPicker(self.vf, None, gui=0,
                               callbacks=[self.setHingePt_cb], immediate=1)
            self.ap.go(modal=0)
        elif val == 'cursel':
            #use current selection here
            selnodes = self.vf.getSelection()[:]
            len_nodes = len(selnodes)
            if len_nodes>=2:
                self.atomOne=selnodes[0]
                c1 = self.getTransformedCoords(self.atomOne)
                self.lineVertices.append(tuple(c1))
                if selnodes[1]!=self.atomOne:
                    self.atomTwo=selnodes[1]
                else:
                    return 'ERROR'
                c2 = self.getTransformedCoords(self.atomTwo)
                self.lineVertices.append(tuple(c2))
                self.line.Set(vertices=self.lineVertices, faces=[(0,1),])
                self.line.visible=1
                self.spheres.Set(vertices=self.lineVertices)
                self.vf.GUI.VIEWER.Redraw()
                if len_nodes>2:
                    self.atoms = selnodes[2:]
                    c2 = self.getTransformedCoords(self.atomTwo)
                    self.setHingePt_cb(selnodes)
        elif val == 'xyz':
            pt = [float(self.xval.get()), float(self.yval.get()), float(self.zval.get())]
            self.spheres.Set(vertices = (pt,))
            if hasattr(self, 'ap'):
                self.ap.stop()
                del self.ap
            #self.hingeAtoms.Set(visible=0, tagModified=False)
            self.centerList = [map(float,[self.xval.get(),self.yval.get(),self.zval.get()]),]
            self.spheres.Set(vertices = self.centerList)
        else:
            msg = val + "is not a recognized center for SelectInSphere"
            self.vf.warningMsg(msg)
            return 'ERROR'
        #self.drawHinge()
        if hasattr(self, 'form') and self.form:
            self.form.lift()
        self.old_val=val
        self.vf.GUI.VIEWER.Redraw()


    def get_state(self, event=None):
        val = 'no first atom'
        if self.atomOne is not None:
            if self.atomTwo is not None: 
                if len(self.atoms):
                    val = 'ready'
                else:
                    val = 'no atoms to move'
            else:
                val = 'no second atom'
        return val


    def check_torscount(self, atoms):
        res = atoms.parent.uniq()
        for r in res:
            if not hasattr(r, 'torscount'):
                r.atoms.bonds[0].possibleTors = 0
                r.atoms.bonds[0].activeTors = 0
                r.torscount = 0


    def saveHinge(self, newval, event=None):
        if len(newval)!=2:
            print newval, ' not required format which is [(atom1, atom2), atoms_to_move]'
            return
        if len(newval[0])!=2:
            print newval, ' not required format which is [(atom1, atom2), atoms_to_move]'
            return
        (atomOne, atomTwo), atoms = newval
        dict = self.vf.flexDict
        dict.setdefault('hinge_list', [])
        self.check_torscount(atoms)
        #print "flexDict['hinge_list']=", dict['hinge_list']
        kw = {'log': 1}
        apply(self.doitWrapper,((atomOne, atomTwo), atoms), kw)
        

    def setHinge_cb(self, event=None):
        #callback for Set Hinge button
        #print "in setHinge_cb with selectionBase= ", self.selectionBase.get()
        state = self.get_state()
        if state=='ready':
            newval = [(self.atomOne, self.atomTwo), self.atoms]
            #print "calling saveHinge with ", newval
            self.saveHinge(newval)
        elif state=='no atoms to move': 
            self.warningMsg('atoms to be moved by hinge must be specified first')
        elif state=='no second atom':
            self.warningMsg('second hinge atom for hinge must be specified first')
        else:
            self.warningMsg('two hinge points must be specified first!')
                


    def setHingePt_cb(self, atoms, event=None):
        if not len(atoms):
            return
        state = self.get_state()
        if state=='no first atom':
            self.atomOne=atoms[0]
            c1 = self.getTransformedCoords(self.atomOne)
            self.lineVertices.append(tuple(c1))
        elif state=='no second atom':
            self.atomTwo=atoms[0]
            c2 = self.getTransformedCoords(self.atomTwo)
            self.lineVertices.append(tuple(c2))
            self.line.Set(vertices = self.lineVertices, faces=[(0,1),])
        elif state=='no atoms to move':
            self.warningMsg('atoms to be moved by hinge must be specified first')
        else:
            newval = [(self.atomOne, self.atomTwo), self.atoms]
            self.saveHinge(newval)
            self.lineVertices = []
        self.spheres.Set(vertices=self.lineVertices)
        self.hingeAtoms.Set(vertices=atoms.coords,) 
        self.vf.GUI.VIEWER.Redraw()


    def guiCallback(self):
        """  """
        if not hasattr(self,'form'):
            self.buildForm()
        else:
            self.form.root.deiconify()
        

    def close_cb(self, event=None):
        #at this point, just cleanup
        self.line.Set(visible=False)
        self.spheres.Set(visible=False)
        self.hingeAtoms.Set(visible=False)
        self.vf.GUI.VIEWER.Redraw()
        self.form.root.withdraw()
        self.stop()


    def changeFlexResCt(self,delta):
        n=self.vf.flexDict.get('flex_residues_number',0)
        n=n+delta
        self.vf.flexDict['flex_residues_number']=n
        msg = 'current %d flexible residues'%n
        self.vf.GUI.pickLabel.configure(text=msg)


    def __call__(self, (atomOne, atomTwo), atoms,  **kw):
        """
        atomOne: first atom in hinge
        atomTwo: second atom in hinge
        atoms: all atoms which will be moved by this hinge

        """
        kw['topCommand'] = 0
        kw['busyIdle'] = 1
        kw['log'] = 0
        #self.setUpDisplay()
        atomOne = self.vf.expandNodes(atomOne)
        if not len(atomOne):
            return "ERROR"
        atomOne = atomOne[0]
        atomTwo = self.vf.expandNodes(atomTwo)
        if not len(atomTwo):
            return "ERROR"
        atomTwo = atomTwo[0]
        atoms = self.vf.expandNodes(atoms)
        if not len(atoms):
            return "ERROR"
        apply(self.doitWrapper, ((atomOne, atomTwo), atoms), kw)


    def doit(self, (atomOne, atomTwo), atoms, **kw):
        """
        atomOne
        """
        dict = self.vf.flexDict
        #if 'hinge_list' is not in dict.keys, add it with [] as its value
        hinge_list = dict.setdefault('hinge_list', [])
        newval = ((atomOne, atomTwo), atoms)
        if not newval in hinge_list:
            dict['hinge_list'].append(newval)
            dict['all_hinge_atoms'], dict['non_hinge_atoms'] = self.getAllHingeAtoms()
            #created a log string here
        else:
            print newval, " already saved in flexDict['hinge_list']"
            return "ERROR"


    def stop(self):
        self.dismiss()


    def getObjects(self,pick):
        print "in get Objects with ", pick
        if 'flex_residues' in self.vf.flexDict.keys():
            flexMol = self.vf.flexDict['flex_residues'].top.uniq()[0]
        elif 'hinge_list' in self.vf.flexDict.keys():
            #hinge_list= [[(at1,at2),atoms]]
            #get the first atom's top
            flexMol = self.vf.flexDict['hinge_list'][0][0][0].top
        else:
            return 'ERROR'
        flexMolGeom = flexMol.geomContainer.geoms['bonded']
##         flexMolGeom = flexMol.geomContainer.geoms['lines']
        for o, val in pick.hits.items(): #loop over geometries
            primInd = map(lambda x: x[0], val)
            if o != flexMolGeom: continue
            else: g = o.mol.geomContainer
            if g.geomPickToBonds.has_key(o.name):
                func = g.geomPickToBonds[o.name]
                if func: return func(o, primInd)
            else:
                l = []
                bonds = g.atoms[o.name].bonds[0]
                for i in range(len(primInd)):
                    l.append(bonds[int(primInd[i])])
                return BondSet(l)


    def done_cb(self):
        self.guiUp = 0
        if self.form: self.form.withdraw()


    def dismiss(self):
        self.vf.setICOM(self.save, topCommand=0)
        self.save = None
        self.done_cb()


AF_SetHingeGUI = CommandGUI()
AF_SetHingeGUI.addMenuCommand('AutoToolsBar', menuText['AutoFlexMB'],\
        menuText['Set Hinge'])



class AF_SetBondRotatableFlag(SetRotatableBonds):
    """set the flag that tells whether a bond is rotatable in aan AutoDock
    ligand"""


    def onAddCmdToViewer(self):
        if not hasattr(self.vf, 'setICOM'):
            self.vf.loadCommand('interactiveCommands', 'setICOM', 'Pmv') 
        self.torsStr = None


    def buildCol(self, mol, torscount):
        self.vf.ADflex_processResidues.torsStr('Current Active Torsions: ' + str(torscount))
        currentbonds=mol.geomContainer.atoms['bonded'].bonds[0]
##         currentbonds=mol.geomContainer.atoms['lines'].bonds[0]
        col = []
        for b in currentbonds:
            if b.possibleTors:
                if b.activeTors: col.append((0,1,0))
                else: col.append((1,0,1))
            else:
                col.append((1,0,0))
        mol.geomContainer.geoms['bonded'].Set(materials=col,
                                              inheritMaterial=False,
                                              matBind=viewerConst.PER_PART)
##         mol.geomContainer.geoms['lines'].Set(materials=col,
##             matBind=viewerConst.PER_PART)
        self.vf.GUI.VIEWER.Redraw()
        

    def setupUndoBefore(self, atoms, rotatable, flexRes=False):
        self.addUndoCall( (atoms, not rotatable), {'redraw':1,'flexRes':flexRes},
                  self.name )
        
        
    def doit(self, atoms, rotatable, flexRes=False):
        """arguments are the numbers of the two atoms"""

        assert rotatable in [ 0, 1 ]
        if len(atoms) < 2: 
            return 'ERROR'
        
        bonds = atoms[:2].bonds
        if len(bonds[0])==0:
            print 'ERROR: no bond between ...'
            return 'ERROR'

        bond = bonds[0][0]

        if not hasattr(bond, 'possibleTors') or bond.possibleTors==0: 
            return 'ERROR'

        # update torsion count
        res = atoms[0].parent
        mol = res.top
        torscount = res.torscount

        if bond.activeTors!=rotatable:
            if rotatable==0: inc = -1
            else: inc = 1
            torscount = torscount + inc
            bond.activeTors = rotatable
            res.torscount = torscount
        else: 
            return 'ERROR'

        ttc = self.vf.flexDict['torscount']
        ttc = ttc + inc
        if not flexRes:
            self.vf.ADflex_setResidues.torscount = ttc
            self.vf.ADflex_processResidues.torscount = ttc
            self.vf.flexDict['torscount'] = ttc
            if self.vf.hasGui:
                self.vf.ADflex_processResidues.buildCol(mol, ttc)
        else:
            #print 'ttc=', ttc, ' inc=', inc, 'newtorscount =', ttc
            self.vf.ADflex_processHingeResidues.torscount = ttc
            self.vf.flexDict['torscount'] = ttc
            if self.vf.hasGui:
                self.vf.ADflex_processHingeResidues.buildCol(mol, ttc)


    def __call__(self, atoms, rotatable, **kw):
        """None <- setBondRotatableFlag(atoms, rotatable, **kw)
        rotatable can be either 1 or 0
        """
        atoms = self.vf.expandNodes(atoms)
        if len(atoms)<2: return
        assert isinstance( atoms[0], Atom )
        kw.setdefault('flexRes', False)
        apply( self.doitWrapper, (atoms, rotatable), kw ) 



class AF_StepBack(MVCommand):

    def onAddCmdToViewer(self):
        if not hasattr(self.vf, 'flexDict'):
            self.vf.flexDict = {}
        if not hasattr(self.vf, 'colorByAtomType'):
            self.vf.loadCommand('colorCommands', 'colorByAtomType', 'Pmv')


    def guiCallback(self):
        """called each time the 'Redisplay Macromolecule' button is pressed"""
        dict=self.vf.flexDict
        if not dict.has_key('macroname') or not dict['macroname']:
            msg='No Protein designated'
            self.vf.warningMsg(msg)
            return
        if not dict.has_key('flex_residues') or not dict['flex_residues']:
            msg='No residues designated'
            self.vf.warningMsg(msg)
            return
        self.doitWrapper(dict['macroname'])


    def doit(self, molName):
        mols = self.vf.Mols.NodesFromName(molName)
        if not len(mols)==1:
            msg='error in finding ' + molName+ ' in self.vf.Mols'
            self.vf.warningMsg(msg)
            return 'ERROR'
        mol=mols[0]
        #need to redraw the molecule and close the ProcessResidue panel
        self.vf.displayLines(mol, topCommand=0)
        self.vf.colorByAtomType(mol.allAtoms, topCommand=0, redraw=1)
        #need to reset total torscounts:
        self.vf.flexDict['torscount'] = 0
        self.vf.ADflex_setResidues.torscount = 0
        c = self.vf.ADflex_processResidues
        c.torscount = 0
        c.dismiss()
            

AF_StepBackGUI = CommandGUI()
AF_StepBackGUI.addMenuCommand('AutoToolsBar', menuText['AutoFlexMB'], menuText['Step Back'])

        

class AF_FlexFileWriter(MVCommand):
    """ allows user to choose an output filename and write autoflex stuff"""


    def onAddCmdToViewer(self):
        self.writer = PdbqtWriter()
        self.writer.recType = 'none'
        self.debug = False


    def guiCallback(self):
        """called each time the 'writeFlexFile' button is pressed"""
        outfile = self.vf.askFileSave(types=[('autoflex file:', '*.pdbqt')],
                title = 'AutoFlex File:')
        if outfile:
            self.doitWrapper(outfile)


    def __call__(self, outfile=None, **kw):
        if not outfile: return 'ERROR'
        apply(self.doitWrapper, (outfile,), kw)


    def doit(self, outfile):
        self.writtenAtoms = []
        outfileptr=open(outfile,'w')
        dict=self.vf.flexDict
        if not dict.has_key('macrofilename'):
            self.vf.warningMsg("No Protein File Selected!")
            return 'ERROR'
        macrofilename = dict['macrofilename']
        if not dict.has_key('flex_residues') and not dict.has_key('hinge_list'):
            self.vf.warningMsg("Must Specify Residues for Flexible Side Chains or set Hinge First")
            return 'ERROR'
        flex_residues=dict.get('flex_residues',[])
        flexResTorsCt = 0
        for r in flex_residues:
            flexResTorsCt = flexResTorsCt + r.torscount
        flex_types = {}
        dict = self.vf.flexDict
        self.outatom_counter = 1
        self.torsionCtr = 0
        for item in flex_residues:
            self.writeResidue(item, outfileptr)
            ###don't add residue's torsdof
            ###self.torsdof += item.torsdof
            for a in item.atoms:
                flex_types[a.autodock_element] = 1
        if 'hinge_list' in dict.keys():
            for item in dict['hinge_list']:
                self.writeHinge(item, outfileptr)
                #include hinge atoms' types, also
                for a in item[0]:
                    flex_types[a.autodock_element] = 1
                for a in item[1]:
                    flex_types[a.autodock_element] = 1
        outfileptr.close()
        dict['flex_filename'] = outfile
        dict['flex_types'] = flex_types.keys()


    def writeHinge(self, item, outfileptr):
        #item format is [(atomOne, atomTwo), atoms_to_move]
        ###print "in writeHinge with ",
        ###for a in item[0]:
        ###    print a.full_name(), 
        ###print
        ###print len(item[1])
        ###for a in item[1]:
        ###    print a.parent.parent.id,
        ###print
        (atomOne, atomTwo), atoms_to_move = item
        atomOne.rnum = 0
        atomTwo.rnum = 0 
        atoms_to_move.rnum = -1
        res = atomOne.parent
        outfileptr.write("BEGIN_RES %s %s%4s\n" %(res.type, res.parent.id, res.number))
        #change this when start to support nested torsions
        # get total number of active torsions plus one for the hinge itsel
        parents_atoms_to_move = atoms_to_move.parent.uniq()
        for res in parents_atoms_to_move:
            #@@WHEN DOES THIS HAPPEN???
            if not hasattr(res, 'torscount'): 
                res.torscount = 0
        ntors = Numeric.add.reduce(parents_atoms_to_move.torscount) + 1
        oStr = "REMARK    %d active torsions:\nREMARK  status: ('A' for Active; 'I' for Inactive)\n" %(ntors)
        outfileptr.write(oStr)
        #!@@! these are original numbers...
        atomOne.number = self.outatom_counter
        self.outatom_counter = self.outatom_counter + 1
        atomTwo.number = self.outatom_counter 
        remark_counter = 1
        oStr = "REMARK    %d  A    between atoms: %s_%d  and %s_%d\n" %(remark_counter, atomOne.name, atomOne.number, atomTwo.name, atomTwo.number)
        outfileptr.write(oStr)
        remark_counter += 1
        #need to write remarks about all the nested torsions here...
        all_moving = AtomSet([atomOne , atomTwo]) + atoms_to_move
        bondset = AtomSet(all_moving).bonds[0].get(lambda x: hasattr(x, 'activeTors') and x.activeTors)
        bondset.written = 0
        bTors = 'I'
        for b in bondset:
            if not (b.atom1 in all_moving and b.atom2 in all_moving):
                continue
            if hasattr(b, 'activeTors') and b.activeTors:
                outstring = "REMARK " +" %3d   A    between atoms: %-3s  and  %-3s \n" %(remark_counter, b.atom1.name + '_' + b.atom1.parent.name+'_'+b.atom1.name, b.atom2.name + '_' + b.atom2.parent.name + '_'+b.atom2.name)
                outfileptr.write(outstring)
                remark_counter += 1
        outfileptr.write("ROOT\n")
        self.writer.write_atom(outfileptr, atomOne)
        atomOne.used = 1
        outfileptr.write("ENDROOT\n")
        # this is already set here...atomTwo.number = self.outatom_counter
        oStr = "BRANCH %d %d\n" %(atomOne.number, atomTwo.number)
        outfileptr.write(oStr)
        self.writer.write_atom(outfileptr, atomTwo)
        self.outatom_counter = self.outatom_counter + 1
        atomTwo.used = 1
        self.writtenAtoms = []
        #start with atomOne in order to proceed by writing atoms bonded to it..
        atoms_to_write = AtomSet([atomOne])
        atoms_to_write.extend(atoms_to_move)
        for at1 in atoms_to_write:
            if not at1.used:
                at1.number = self.outatom_counter
                self.writer.write_atom(outfileptr, at1)
                self.outatom_counter = self.outatom_counter + 1
                at1.used = 1
            for bond in at1.bonds:
                at2 = bond.atom1
                if at2==at1:
                    at2 = bond.atom2
                if at2.used or at2 not in atoms_to_move:
                    continue
                if not hasattr(bond, 'activeTors') or not bond.activeTors:
                    # don't do anything for atoms bonded by non-rotatable bonds
                    continue
                else:
                    #here have to write subtree
                    marker = self.outatom_counter
                    outstring = "BRANCH %3d %3d\n"%(at1.number, marker)
                    outfileptr.write(outstring)
                    #print "calling writesubtree with ", at1.parent.name, ' ', at1.name, '-', at2.name
                    self.writeSubtree(outfileptr, at1, at2)
                    outstring = "ENDBRANCH %3d %3d\n"%(at1.number, marker)
                    outfileptr.write(outstring)
        outfileptr.write( "ENDBRANCH %3d %3d\n"%(atomOne.number, atomTwo.number))
        outfileptr.write("END_RES %s %s%4s\n" %(res.type, res.parent.id, res.number))


    def writeResidue(self, res, outfileptr):
        # should have ALREADY eliminated residues with no active torsions
        # need to have already done autotors stuff to this residue's sidechain:
        # atoms should know which charge to use
        # also must have these fields: torscount, torsdof, and bonds know if active
        # also res.bondlist
        #start output
        outfileptr.write("BEGIN_RES %s %s%4s\n" %(res.type, res.parent.id, res.number))
        #first write out remarks about torsions
        outfileptr.write("REMARK  " + "%d" %(res.torscount) + " active torsions:\n")
        outfileptr.write("REMARK  status: ('A' for Active; 'I' for Inactive)\n")
        #only want to process bonds pertaining to sideChain atoms
        for b in res.bondlist:
            if b.activeTors == 1: bTors = 'A'
            else:   bTors = 'I'
            if b.possibleTors:
                if bTors=='A':
                    self.torsionCtr = self.torsionCtr + 1
                    outstring = "REMARK " +" %3d  %s    between atoms: %-3s  and  %-3s \n" %(self.torsionCtr,bTors,b.atom1.name,b.atom2.name)
                else:
                    outstring = "REMARK " +"      %s    between atoms: %-3s  and  %-3s \n" %(bTors,b.atom1.name,b.atom2.name)
                outfileptr.write(outstring)
        #next write out  root which is always the CA atom
        outfileptr.write("ROOT\n")
        assert hasattr(res, 'rootlist')
        #reset used field to serve as visited flag
        res.atoms.used = 0
        #rootlist grows to include atoms up to first active tors in each subtree
        for at in res.rootlist:
            at.used = 1
            for bond in at.bonds:
                if bond.activeTors and bond.possibleTors: continue
                at2 = bond.atom1
                if at2==at: at2 = bond.atom2
                #only track and write out bonds to the sideChain atoms
                if at2 not in res.sideChain: continue
                if at2.used: continue
                if at2 not in res.rootlist:
                    res.rootlist.append(at2)
                    at2.rnum = len(res.rootlist)
        #remove all atoms which have no rotatable bonds
        #from BOTH the rootlist and the sideChain atoms
        #this means they will be written in the rigid portion
        badList = AtomSet([])
        for at in res.rootlist:
            hasTorsion=0
            for b in at.bonds:
                if b.activeTors:
                    hasTorsion=1
                    break
            if not hasTorsion:
                badList.append(at)
        if len(badList) and len(badList)!=len(res.rootlist):
            res.rootlist = res.rootlist - badList
            res.sideChain = res.sideChain - badList
        #now visit atoms connected to expanded rootlist
        for at in res.rootlist:
            at.number = self.outatom_counter
            self.writer.write_atom(outfileptr, at)
            at.newindex = self.outatom_counter
            at.used = 1
            self.outatom_counter = self.outatom_counter + 1
        outfileptr.write("ENDROOT\n")
        #last write out the rest of the stuff, using writeSubtree.....
        for at in res.rootlist:
            for b in at.bonds:
                at2 = b.atom1
                if at2 == at: at2 = b.atom2
                if at2 in res.rootlist:
                    continue
                if at2 not in res.sideChain:
                    continue
                if at2.used:
                    continue
                self.process(at, at2, outfileptr)
        outfileptr.write("END_RES %s %s%4s\n" %(res.type, res.parent.id, res.number))


    def process(self, fromAtom, nextAtom, outfileptr):
        startIndex = fromAtom.number
        endIndex = self.outatom_counter
        outstring = "BRANCH %3d %3d\n"%(startIndex, endIndex)
        outfileptr.write(outstring)
        queue = self.writeBreadthFirst(outfileptr, fromAtom, nextAtom)
        if self.debug: print fromAtom.name, ':', nextAtom.name, ': queue=', queue
        if len(queue):
            for fromAtom, nextAtom in queue:
                if self.debug: print " processing queue entry: ", fromAtom.name, '-', nextAtom.name
                self.process(fromAtom, nextAtom, outfileptr)
        outstring = "ENDBRANCH %3d %3d\n"%(startIndex, endIndex)
        outfileptr.write(outstring)


    def writeLevel(self, atom, outfptr):
        """
        write all atoms bonded to atoms bonded to this atom by non-rotatable
        bonds
        """
        if self.debug: 
            print "\n\nin writeLevel with ", atom.name, " outatom_counter=", self.outatom_counter
            print "len(", atom.name, ").bonds=", len(atom.bonds)
        queue = []
        nextAts = []
        for b in atom.bonds:
            if self.debug:
                print "processing b=", b.atom1.name, '-', b.atom2.name, ' activeTors=', b.activeTors
                print 'atom1 in writtenAtoms=', b.atom1 in self.writtenAtoms
                print 'atom2 in writtenAtoms=', b.atom2 in self.writtenAtoms
            if b.activeTors: 
                at2 = b.atom1
                if at2==atom: at2=b.atom2
                queue.append((atom, at2))
                if self.debug: print atom.name, 'wL: queue=', queue
                continue
            a2 = b.atom1
            if a2==atom:
                a2 = b.atom2
            if a2.used: 
                if self.debug: print "!!a2 is already used!!", a2.name
                continue    
            if a2 not in self.writtenAtoms:
                a2.number = a2.newindex = self.outatom_counter
                if self.debug: print "writeLevel: wrote bonded atom named=", a2.name, 'a2.used=', a2.used
                self.writer.write_atom(outfptr, a2)
                self.writtenAtoms.append(a2)
                a2.used = 1
                self.outatom_counter+=1
                nextAts.append(a2)
        for a2 in nextAts:
            if self.debug: 
                print 'in for nextAts loop with a2=', a2.name
                print 'calling wL'
            nq = self.writeLevel(a2, outfptr)
            if len(nq):
                if self.debug: print "extending queue with", nq
                queue.extend(nq)
        if self.debug:
            print 'returning queue=', queue
        return queue
            

    def writeBreadthFirst(self, outfptr, fromAtom, startAtom):
        """
            None <-writeBreadthFirst(outfptr, fromAtom, startAtom)
            writeBreadthFirst visits all the atoms in the current level
            then the first level down etc in a Breadth First Order traversal. 
                            1                <-1
                        5 6   7 8            <-3
                     9 10   11 12            <-4
            It is used to write out the molecule with the correct format 
            for AutoDock. Specifically, appropriate BRANCH/ENDBRANCH 
            statements are added. 
        """
        if self.debug: 
            print "in wBF with fromAtom=", fromAtom.name, '+ startAtom=', startAtom.name, 'startAtom.used=', startAtom.used
        queue = []
        if startAtom.used==0:
            startAtom.used = 1
            startAtom.number = startAtom.newindex = self.outatom_counter
            self.writer.write_atom(outfptr,startAtom)
            if self.debug: print 'wBF: wrote ', startAtom.name
            self.writtenAtoms.append(startAtom)
            self.outatom_counter += 1
            if self.debug: print "self.outatom_counter=", self.outatom_counter
            activeTors = []
            #outfptr.write(outstring)
            for bond in startAtom.bonds:
                if not hasattr(bond, 'activeTors'):
                    continue
                at2 = bond.atom1
                if at2==startAtom: at2 = bond.atom2
                if at2==fromAtom: continue  #this is current bond
                elif not at2.used:
                    if bond.activeTors:
                        queue.append((startAtom,at2))
                    else:
                        at2.number = at2.newindex = self.outatom_counter
                        if self.debug: 
                            print "\n\nwriting and calling wL with nA=", at2.name, '-', at2.number
                        self.writer.write_atom(outfptr, at2)
                        if self.debug: print 'wBF2: wrote ', at2.name
                        at2.written = 1
                        self.writtenAtoms.append(at2)
                        at2.newindex = self.outatom_counter
                        self.outatom_counter = self.outatom_counter + 1
                        if self.debug: print '!!!2:calling wL'
                        newQ = self.writeLevel(at2, outfptr)
                        if self.debug: print "newQ=", newQ
                        at2.used = 1
                        if len(newQ): 
                            if self.debug: print "@@@@len(newq)=", len(newQ)
                            queue.extend(newQ)
                            if self.debug: print "queue=", queue
            if self.debug: 
                print " currently queue=",
                for atom1, atom2 in queue: 
                    print atom1.name, '-', atom2.name, ',',
                print
        return  queue


    def writeSubtree(self, outfptr, fromAtom, startAtom):
        """
            None <-writeSubtree(outfptr, fromAtom, startAtom)
            writeSubtree recursively visits the atoms in the current 
            'subtree' of the molecule in a Depth First Order traversal. 
            It is used to write out the molecule with the correct format 
            for AutoDock. Specifically, appropriate BRANCH/ENDBRANCH 
            statements are added. 
        """
        if startAtom.used==0:
            startAtom.used = 1
            at = startAtom
            for bond in startAtom.bonds:
                if bond.activeTors: 
                    continue
                marker = self.outatom_counter
                nextAtom = bond.atom1
                if nextAtom==startAtom: 
                    nextAtom = bond.atom2
                if nextAtom==fromAtom: 
                    continue
                if not nextAtom.used:
                    if hasattr(bond,'incycle'):
                        if not hasattr(nextAtom, 'cycleout'):
                            nextAtom.cycleout = 1
                            nextAtom.newindex = self.outatom_counter
                            nextAtom.number = self.outatom_counter
                            self.writer.write_atom(outfptr,nextAtom)
                            self.writtenAtoms.append(nextAtom)
                            self.outatom_counter = self.outatom_counter+1
                    else:
                        nextAtom.newindex = self.outatom_counter
                        nextAtom.number = self.outatom_counter
                        self.writer.write_atom(outfptr,nextAtom)
                        self.writtenAtoms.append(nextAtom)
                        self.outatom_counter = self.outatom_counter+1
            for bond in startAtom.bonds:
                marker = self.outatom_counter
                nextAtom = bond.atom1
                if nextAtom==startAtom: 
                    nextAtom = bond.atom2
                if nextAtom==fromAtom: 
                    continue
                if not nextAtom.used:
                    testcond = len(nextAtom.bonds)
                    if bond.activeTors and bond.possibleTors:
                        if testcond >0:
                            outstring = "BRANCH %3d %3d\n"%(at.newindex,marker)
                            outfptr.write(outstring)
                        nextAtom.newindex = self.outatom_counter
                        nextAtom.number = self.outatom_counter
                        self.writer.write_atom(outfptr,nextAtom)
                        self.writtenAtoms.append(nextAtom)
                        self.outatom_counter = self.outatom_counter+1
                    self.WriteSubtree(startAtom, nextAtom)
                    if bond.activeTors and bond.possibleTors:
                        if testcond >0:
                            outstring = "ENDBRANCH %3d %3d\n"%(at.newindex,marker)
                            outfptr.write(outstring)
        return


    def writeFile_cb(self):
        """This callback function allows the user to select an output file, if none has
        been designated before. It opens that file and calls the write method of 
        AutoTors to write pdbqt-formatted for AutoDock4, including appropriate keywords
        such as ROOT/ENDROOT, BRANCH/ENDBRANCH etc."""
        if self.rootnum > 0:
            if self.outfile == None:
                self.message( "no output file specified, please select:\n")
                from Dialog import Dialog
                from FileDialog import SaveFileDialog
                fd2 = SaveFileDialog(self.ROOT,title="File for AutoTors Formatted Output")
                self.outfile = fd2.go (key = "test")
            if self.outfile != None:
                f=open(self.outfile, 'w')
                self.write()
        else: 
            self.message( "WRITE ERROR:no root atoms specified yet!\n")


AF_FlexFileWriterGUI = CommandGUI()
AF_FlexFileWriterGUI.addMenuCommand('AutoToolsBar', menuText['AutoFlexMB'], \
            menuText['writeFlexible'], cascadeName = menuText['WriteMB'])



class AF_RigidFileWriter(MVCommand):
    """ allows user to choose an output filename and write reduced macromolecule"""


    def onAddCmdToViewer(self):
        if not hasattr(self.vf, 'writePDBQT'):
            self.vf.loadCommand('fileCommands','writePDBQT', 'Pmv')
        self.PdbqtWriter = PdbqtWriter()
        self.PdbqtWriter.recType = 'none'


    def __call__(self, outfile=None, **kw):
        """None<-ADflex_writeRigidFile(outfile)"""
        if not outfile: return 'ERROR'
        apply(self.doitWrapper, (outfile,), kw)


    def doit(self, outfile):
        #write the reduced macromolecule here
        d = self.vf.flexDict
        if not (d.has_key('flex_residues') or  d.has_key('all_hinge_atoms')) :
            msg = 'flexible residues or hinge must be specified first!'
            self.vf.warningMsg(msg)
            return 'ERROR'
        flex_residues = d.get('flex_residues', ResidueSet())
        all_hinge_atoms = d.get('all_hinge_atoms', AtomSet())
        mol = d['macromol']
        #mol = self.vf.Mols.NodesFromName(d['macroname'])
        rigidResidues = d.get('rigidResidues', mol.chains.residues)
        #build AtomSet of rigidResidue atoms PLUS N,C and O from flex
        rigidAtoms = rigidResidues.findType(Atom) + d.get('non_hinge_atoms', AtomSet())
        #flexAtoms = flex_residues.findType(Atom)
        outatoms = mol.allAtoms
        if len(flex_residues):
            flexSideChain = flex_residues.sideChain
            outatoms = outatoms - flexSideChain
        if len(all_hinge_atoms):
            #remove the hinge atoms, also...
            outatoms = outatoms - all_hinge_atoms
        #outatoms = outatoms - hinge_atoms
        #renumber this atoms?
        outatoms.number = range(1, len(outatoms)+1)
        outatoms_type = {}
        for a in outatoms:
            outatoms_type[a.autodock_element] = 1
        outptr = open(outfile, 'w')
        for at in outatoms:
            self.PdbqtWriter.write_atom(outptr, at)
        outptr.close()
        d = self.vf.flexDict
        d['rigid_filename'] = outfile
        d['rigid_types'] = outatoms_type.keys()


    def guiCallback(self):
        outfile = self.vf.askFileSave(types=[('pdbqt file', '*.pdbqt')],
                title = 'Autoflex Non-Flexible Residue Output File:')
        if outfile: 
            self.doitWrapper(outfile)


AF_RigidFileWriterGUI = CommandGUI()
AF_RigidFileWriterGUI.addMenuCommand('AutoToolsBar',  menuText['AutoFlexMB'], \
            menuText['writeRigid'], cascadeName = menuText['WriteMB'])



class AF_LigandDirectoryWriter(MVCommand):
    """ allows user to choose a directory of formatted ligands and write flexible output for each ligand"""


    def onAddCmdToViewer(self):
        if not hasattr(self.vf, 'writePDBQ'):
            self.vf.loadCommand('fileCommands','writePDBQ', 'Pmv')


    def __call__(self, ligDir=None, **kw):
        """None<-ADflex_writeFlexDir(ligDir)"""
        if not ligDir: return 'ERROR'
        apply(self.doitWrapper, (ligDir,), kw)


    def doit(self, ligDir):
        fileList= os.listdir(ligDir)
        for item in fileList:
            xx=string.split(item,'.')
            if len(xx)>1 and xx[-1]=='pdbqt':
                ligFileFullName=string.join((ligDir, item),'/')
                self.vf.ADflex_readLigand(item,1)
                ##9/11 FIX THIS STUPID NAME
                outputFile='autoflex_'+item
                self.vf.ADflex_writeFlexFile(outputFile)
        

    def guiCallback(self):
        """called each time the 'writeFlexDir' button is pressed"""
        outfile = self.vf.askFileOpen(types=[('select any formatted autotors file:', '*.out.pdbqt'),('formatted pdbqt file:', '*.pdbqt')],
                title = 'To Set Directory: Select any autotors formatted file:')
        if outfile:
            ligDir=os.path.split(outfile)[0]
            self.doitWrapper(ligDir)


AF_LigandDirectoryWriterGUI = CommandGUI()
AF_LigandDirectoryWriterGUI.addMenuCommand('AutoToolsBar', menuText['AutoFlexMB'], \
            menuText['writeDir'], cascadeName = menuText['WriteMB'])



commandList = [
    {'name':'ADflex_readMacro','cmd':AF_MacroReader(),'gui':AF_MacroReaderGUI},
    {'name':'ADflex_chooseMacro','cmd':AF_MacroChooser(),'gui':AF_MacroChooserGUI},
    {'name':'ADflex_setResidues','cmd':AF_SelectResidues(),'gui':AF_SelectResiduesGUI},
    {'name':'ADflex_processResidues','cmd':AF_ProcessResidues(),'gui':None},
    #{'name':'ADflex_processHingeResidues','cmd':AF_ProcessHingeResidues(),'gui':None},
    {'name':'ADflex_setBondRotatableFlag','cmd':AF_SetBondRotatableFlag(),'gui':None},
    #{'name':'ADflex_setHinge','cmd':AF_SetHinge(),'gui':AF_SetHingeGUI},
    #{'name':'ADflex_editHinge','cmd':AF_EditHinge(),'gui':None},
    {'name':'ADflex_stepBack','cmd':AF_StepBack(),'gui':AF_StepBackGUI},
    {'name':'ADflex_writeFlexFile','cmd':AF_FlexFileWriter(),'gui':AF_FlexFileWriterGUI},
    {'name':'ADflex_writeRigidFile','cmd':AF_RigidFileWriter(),'gui':AF_RigidFileWriterGUI},
    #{'name':'ADflex_writeFlexDir','cmd':AF_LigandDirectoryWriter(),'gui':AF_LigandDirectoryWriterGUI}
]



def initModule(vf):
    for dict in commandList:
        vf.addCommand(dict['cmd'], dict['name'], dict['gui'])
    
    if vf.hasGui:
        vf.GUI.menuBars['AutoToolsBar'].menubuttons[menuText['AutoFlexMB']].config(bg='tan',underline='-1')    
        if not hasattr(vf.GUI, 'adtBar'):
            vf.GUI.adtBar = vf.GUI.menuBars['AutoToolsBar']
            vf.GUI.adtFrame = vf.GUI.adtBar.menubuttons.values()[0].master
#        if hasattr(vf.GUI, 'ligandLabel'):
#            vf.GUI.ligandLabelLabel.pack_forget()
#            vf.GUI.ligandLabelLabel.pack(side='left')
#            vf.GUI.ligandLabel.pack_forget()
#            vf.GUI.ligandLabel.pack(side='left')
#        else:
#            vf.GUI.ligandLabelLabel = Tkinter.Label(vf.GUI.adtFrame, \
#                            text="Ligand:", bg='tan')
#            vf.GUI.ligandLabelLabel.pack(side='left')
#            vf.GUI.ligandLabel=Tkinter.Label(vf.GUI.adtFrame, text="None", width=4,
#                                     relief='sunken', borderwidth=1,
#                                     anchor='w' )
#            vf.GUI.ligandLabel.pack(side='left')
#        if hasattr(vf.GUI, 'receptorLabel'):
#            vf.GUI.receptorLabelLabel.pack_forget()
#            vf.GUI.receptorLabelLabel.pack(side='left')
#            vf.GUI.receptorLabel.pack_forget()
#            vf.GUI.receptorLabel.pack(side='left')
#        else:
#            vf.GUI.receptorLabelLabel = Tkinter.Label(vf.GUI.adtFrame, \
#                            text="Receptor:", bg='tan')
#            vf.GUI.receptorLabelLabel.pack(side='left')
#            vf.GUI.receptorLabel=Tkinter.Label(vf.GUI.adtFrame, text="None", width=4,
#                                     relief='sunken', borderwidth=1,
#                                     anchor='w' )
#            vf.GUI.receptorLabel.pack(side='left')

        
