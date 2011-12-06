## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#############################################################################
#
# Author: Ruth HUEY, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2005
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/pyAutoDockCommands.py,v 1.21 2008/07/15 22:24:11 rhuey Exp $
#
# $Id: pyAutoDockCommands.py,v 1.21 2008/07/15 22:24:11 rhuey Exp $
#
#
#
#
#
#
#
"""
This Module contains commands which allow the user to compute various
component pair-wise energies of the autodock305 forcefield
and assign the results to individual atoms. 

"""

import Tkinter, numpy.oldnumeric as Numeric, Pmw

from PyAutoDock.MolecularSystem import MolecularSystem
from PyAutoDock.scorer import WeightedMultiTerm, Distance
from PyAutoDock.vanDerWaals import VanDerWaals 
from PyAutoDock.vanDerWaals import HydrogenBonding, HydrogenBonding12_10
from PyAutoDock.vanDerWaals import NewHydrogenBonding, NewHydrogenBonding12_10
from PyAutoDock.InternalEnergy import InternalEnergy
from PyAutoDock.electrostatics import Electrostatics
from PyAutoDock.desolvation import Desolvation, NewDesolvation
from PyAutoDock.AutoDockScorer import AutoDock305Scorer, AutoDock4Scorer
from PyAutoDock.AutoDockScorer import AutoDockTermWeights305, AutoDockTermWeights4

from Pmv.mvCommand import MVCommand, MVAtomICOM, MVBondICOM
from Pmv.stringSelectorGUI import StringSelectorGUI
from MolKit import Read
from MolKit.molecule import AtomSet, Atom, BondSet
from ViewerFramework.VFCommand import Command, CommandGUI
from DejaVu.Geom import Geom
from mglutil.gui.InputForm.Tk.gui import InputFormDescr
from mglutil.gui.BasicWidgets.Tk.thumbwheel import ThumbWheel
from mglutil.gui.BasicWidgets.Tk.customizedWidgets import ExtendedSliderWidget
from mglutil.util.misc import ensureFontCase


parentGeom = Geom("pyADGeoms", shape=(0,0), pickable=0)
parentGeom.isScalable = 0
parentGeom.in_viewer = False

pep_aromList  = ['PHE_CD1', 'PHE_CG', 'PHE_CD2', 'PHE_CE1',\
             'PHE_CE2', 'PHE_CZ', 'TYR_CD1', 'TYR_CG', 'TYR_CD2', 'TYR_CE1',\
             'TYR_CE2', 'TYR_CZ', 'HIS_CD2', 'HIS_CE1', 'HIS_CG', 'TRP_CD1',\
             'TRP_CG', 'TRP_CD2', 'TRP_CE2', 'TRP_CZ2', 'TRP_CH2', 'TRP_CZ3',\
             'TRP_AE3', 'PHE_AD1', 'PHE_AG', 'PHE_AD2', 'PHE_AE1',\
             'PHE_AE2', 'PHE_AZ', 'TYR_AD1', 'TYR_AG', 'TYR_AD2', 'TYR_AE1',\
             'TYR_AE2', 'TYR_AZ', 'HIS_AD2', 'HIS_AE1', 'HIS_AG', 'TRP_AD1',\
             'TRP_AG', 'TRP_AD2', 'TRP_AE2', 'TRP_AZ2', 'TRP_AH2', 'TRP_AZ3',\
             'TRP_AE3']


class PyADCalcClosestPairs(MVCommand):
    """For each atom in one AtomSet, determine the index of closest atom in a second
AtomSet and the distance between them.  A display option exists which draws
lines between each atom in one set and the closest atom to it in the second
set which are closer than 'cutoff' angstrom distance.  By default, cutoff is
set to 8 angstrom but can be set by the user.  The color of the lines drawn
from the atoms in first set to their closest atom in the second set can be
customized. Similarly, the color of the lines drawn from the atoms in the
second set to their closest atom in the first set can be customized."""

    def onAddCmdToViewer(self):
        from DejaVu.IndexedPolylines import IndexedPolylines
        miscGeom = self.vf.GUI.miscGeom
        # we don't need that that check anymore,
        # as miscGeom is always added when we instanciate ViewerFrameworkGUI
        #if miscGeom not in self.vf.GUI.VIEWER.rootObject.children:
        #    self.vf.GUI.VIEWER.AddObject(miscGeom)
        
        
        if not parentGeom.in_viewer:
            self.vf.GUI.VIEWER.AddObject(parentGeom, parent=miscGeom)
            parentGeom.in_viewer = True
        self.first_lines = IndexedPolylines('first_CloseAts', 
                    materials = ((0,1,1),), lineWidth=4, 
                    stippleLines=1, inheritMaterial=0)
        self.second_lines = IndexedPolylines('second_CloseAts', 
                    materials = ((0,1,0),), lineWidth=4, 
                    stippleLines=1, inheritMaterial=0)
        self.vf.GUI.VIEWER.AddObject(self.first_lines, parent=parentGeom)
        self.vf.GUI.VIEWER.AddObject(self.second_lines, parent=parentGeom)


    def onRemoveObjectFromViewer(self, obj):
        ats = obj.findType(Atom)
        ats_with_closest_at = ats.get(lambda x: hasattr(x, 'closest_at'))
        #remove the atoms linked to these atoms
        if len(ats_with_closest_at):
            delattr(ats_with_closest_at, 'closest_at')
            self.first_lines.Set(visible=False)
            self.second_lines.Set(visible=False)
        #remove these atoms from any linked to them
        ats_with_closest_at = self.vf.allAtoms.get(lambda x: hasattr(x, 'closest_at'))
        #just remove all of them for simplicity's sake... instead of looping
        #over them checking whether to remove each one...
        if len(ats_with_closest_at):
            delattr(ats_with_closest_at, 'closest_at')

            
    def buildForm(self):
        ifd = self.ifd = InputFormDescr(title = "For Closest Pairs: Specify 'first' and 'second' group of atoms:")
        self.hideDSel = Tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.hideASel = Tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.displayClosePairsVar = Tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.useSelection = Tkinter.IntVar(master=self.vf.GUI.ROOT)
        ifd.append({'name':'DkeyLab',
                    'widgetType':Tkinter.Label,
                    'text':'For First atoms:',
                    'gridcfg':{'sticky':'w','columnspan':3}})
        ifd.append({'name':'selDRB0',
            'widgetType':Tkinter.Radiobutton,
            'wcfg':{'text':'Use all atoms',
                'variable': self.hideDSel,
                'value':1,
                'command':self.hideDSelector
                },
            'gridcfg':{'sticky':'w'}})
        ifd.append({'name':'selDRB1',
            'widgetType':Tkinter.Radiobutton,
            'wcfg':{'text':'Set atoms to use',
                'variable': self.hideDSel,
                'value':0,
                'command':self.hideDSelector
                },
            'gridcfg':{'sticky':'w', 'row':-1, 'column':1}})
        ifd.append({'name': 'keyDSelector',
            'wtype':StringSelectorGUI,
            'widgetType':StringSelectorGUI,
            'wcfg':{ 'molSet': self.vf.Mols,
                    'vf': self.vf,
                    'all':1,
                    'crColor':(0.,1.,.2),
            },
            'gridcfg':{'sticky':'we', 'columnspan':3 }})
        #now the acceptors
        ifd.append({'name':'AkeyLab',
                    'widgetType':Tkinter.Label,
                    'text':'For Second atoms:',
                    'gridcfg':{'sticky':'w','columnspan':3}})
        ifd.append({'name':'selARB0',
            'widgetType':Tkinter.Radiobutton,
            'wcfg':{'text':'Use all atoms',
                'variable': self.hideASel,
                'value':1,
                'command':self.hideASelector
                },
            'gridcfg':{'sticky':'w'}})
        ifd.append({'name':'selARB1',
            'widgetType':Tkinter.Radiobutton,
            'wcfg':{'text':'Set atoms to use',
                'variable': self.hideASel,
                'value':0,
                'command':self.hideASelector
                },
            'gridcfg':{'sticky':'w', 'row':-1, 'column':1}})
        ifd.append({'name': 'keyASelector',
            'wtype':StringSelectorGUI,
            'widgetType':StringSelectorGUI,
            'wcfg':{ 'molSet': self.vf.Mols,
                    'vf': self.vf,
                    'all':1,
                    'crColor':(0.,1.,.2),
            },
            'gridcfg':{'sticky':'we', 'columnspan':3 }})
        ifd.append({'name': 'showClosePairs',
            'widgetType':Tkinter.Radiobutton,
            'wcfg':{'text': 'Display Close Pairs',
                    'variable': self.displayClosePairsVar,
                    'value':1,
                    'command':self.hideDisplay,
                    },
            'gridcfg':{'sticky':Tkinter.W}})
        ifd.append({'name': 'hideClosePairs',
            'widgetType':Tkinter.Radiobutton,
            'wcfg':{'text': 'Hide Close Pairs',
                    'variable': self.displayClosePairsVar,
                    'value':0,
                    'command':self.hideDisplay,
                    },
            'gridcfg':{'sticky':'w', 'row':-1, 'column':1}})
        ifd.append( {'name':'cutoff',
            'widgetType': ExtendedSliderWidget,
            'wcfg':{'label':'cutoff',
                'minval':1.0, 'maxval':20.0,
                'immediate':0,
                'init':8.,
                'labelsCursorFormat':'%3.1f',
                'entrywcfg':{'width':4},
                'width':250,
                'sliderType':'float',
                'entrypackcfg':{'side':'right'}},
            'gridcfg':{'sticky':'wens', 'columnspan':3}})
        ifd.append({'name':'colorLab',
                    'widgetType':Tkinter.Label,
                    'text':'Set color of lines from:',
                    'gridcfg':{'sticky':'w','columnspan':3}})
        ifd.append({'name':'firstColorLab',
                    'widgetType':Tkinter.Label,
                    'text':'First atoms to closest atom in Second set:',
                    'gridcfg':{'sticky':'w','columnspan':3}})
        ifd.append( {'name':'rvalCtr1',
            'widgetType': Pmw.Counter,
            'wcfg':{
            'labelpos' : 'w',
            'label_text' : 'R:', 
            'entry_width' : 4,
            'entryfield_value' : '1.0',
            'datatype' :  {'counter':'real',
                        'separator': '.'},
            'entryfield_validate' : {'validator' : 'real',
                                           'min' : 0., 'max' : 1.0,
                                           'separator': '.'},
            'increment':0.1},
            'gridcfg':{'sticky':'w',}})
        ifd.append( {'name':'gvalCtr1',
            'widgetType': Pmw.Counter,
            'wcfg':{
            'labelpos' : 'w',
            'label_text' : 'G:', 
            'entry_width' : 4,
            'entryfield_value' : '1.0',
            'datatype' :  {'counter':'real',
                        'separator': '.'},
            'entryfield_validate' : {'validator' : 'real',
                                           'min' : 0., 'max' : 1.0,
                                           'separator': '.'},
            'increment':0.1},
            'gridcfg':{'sticky':'w','row':-1, 'column':1}})
        ifd.append( {'name':'bvalCtr1',
            'widgetType': Pmw.Counter,
            'wcfg':{
            'labelpos' : 'w',
            'label_text' : 'B:', 
            'entry_width' : 4,
            'entryfield_value' : '1.0',
            'datatype' :  {'counter':'real',
                        'separator': '.'},
            'entryfield_validate' : {'validator' : 'real',
                                           'min' : 0., 'max' : 1.0,
                                           'separator': '.'},
            'increment':0.1},
            'gridcfg':{'sticky':'w','row':-1, 'column':2}})
        ifd.append({'name':'secondColorLab',
                    'widgetType':Tkinter.Label,
                    'text':'Second atoms to closest atom in First set:',
                    'gridcfg':{'sticky':'w','columnspan':3}})
        ifd.append( {'name':'rvalCtr2',
            'widgetType': Pmw.Counter,
            'wcfg':{
            'labelpos' : 'w',
            'label_text' : 'R:', 
            'entry_width' : 4,
            'entryfield_value' : '1.0',
            'datatype' :  {'counter':'real',
                        'separator': '.'},
            'entryfield_validate' : {'validator' : 'real',
                                           'min' : 0., 'max' : 1.0,
                                           'separator': '.'},
            'increment':0.1},
            'gridcfg':{'sticky':'w',}})
        ifd.append( {'name':'gvalCtr2',
            'widgetType': Pmw.Counter,
            'wcfg':{
            'labelpos' : 'w',
            'label_text' : 'G:', 
            'entry_width' : 4,
            'entryfield_value' : '1.0',
            'datatype' :  {'counter':'real',
                        'separator': '.'},
            'entryfield_validate' : {'validator' : 'real',
                                           'min' : 0., 'max' : 1.0,
                                           'separator': '.'},
            'increment':0.1},
            'gridcfg':{'sticky':'w','row':-1, 'column':1}})
        ifd.append( {'name':'bvalCtr2',
            'widgetType': Pmw.Counter,
            'wcfg':{
            'labelpos' : 'w',
            'label_text' : 'B:', 
            'entry_width' : 4,
            'entryfield_value' : '1.0',
            'datatype' :  {'counter':'real',
                        'separator': '.'},
            'entryfield_validate' : {'validator' : 'real',
                                           'min' : 0., 'max' : 1.0,
                                           'separator': '.'},
            'increment':0.1},
            'gridcfg':{'sticky':'w','row':-1, 'column':2}})
        ifd.append({'widgetType': Tkinter.Button,
            'text':'Ok',
            'wcfg':{'bd':6},
            'gridcfg':{'sticky':'ew', 'columnspan':2 },
            'command':self.Accept_cb})
        ifd.append({'widgetType': Tkinter.Button,
            'text':'Cancel',
            'wcfg':{'bd':6},
            'gridcfg':{'sticky':'ew', 'column':2,'row':-1},
            'command':self.Close_cb})
        self.form = self.vf.getUserInput(self.ifd, modal=0, blocking=0)
        self.hideDSelector()
        self.hideASelector()
        self.hideDisplay()


    def guiCallback(self):
        if not len(self.vf.Mols):
            self.warningMsg('no molecules in viewer')
            return 
        if not hasattr(self, 'ifd'):
            self.buildForm()
        else:
            self.form.deiconify()


    def Close_cb(self, event=None):
        self.form.withdraw()


    def Accept_cb(self, event=None):
        self.form.withdraw()
        firstNodesToCheck = self.ifd.entryByName['keyDSelector']['widget'].get()
        if not len(firstNodesToCheck):
            self.warningMsg('ERROR! no atoms in first group specified')
            return "ERROR"
        else:
            firstAts = firstNodesToCheck.findType(Atom)
            if self.useSelection.get():
                curSel = self.vf.getSelection()
                if len(curSel):
                    curAts = curSel.findType(Atom)
                    if curAts != self.vf.allAtoms:
                        firstAts = curAts.inter(firstAts)
                        if firstAts is None:
                            self.warningMsg('ERROR! no atoms in first group specified')
                            return "ERROR"
                else:
                    self.warningMsg('ERROR! no atoms in first group specified')
                    return "ERROR"
        secondNodesToCheck = self.ifd.entryByName['keyASelector']['widget'].get()
        if not len(secondNodesToCheck): 
            self.warningMsg('ERROR! no atoms in second group specified')
            return "ERROR"
        else:
            secondAts = secondNodesToCheck.findType(Atom)
            if self.useSelection.get():
                curSel = self.vf.getSelection()
                if len(curSel):
                    curAts = curSel.findType(Atom)
                    if curAts != self.vf.allAtoms:
                        secondAts = curAts.inter(secondAts)
                        if not secondAts:  #NB inter returns empty set!
                            self.warningMsg('ERROR! no atoms in second group specified')
                            return "ERROR"
                else:
                    self.warningMsg('ERROR! no atoms in second group specified')
                    return "ERROR"
        self.Close_cb()
        display = self.displayClosePairsVar.get()
        cutoff = self.ifd.entryByName['cutoff']['widget'].get()
        mat1 = (float(self.ifd.entryByName['rvalCtr1']['widget'].get()),
                float(self.ifd.entryByName['gvalCtr1']['widget'].get()),
                float(self.ifd.entryByName['bvalCtr1']['widget'].get()),)
        mat2 = (float(self.ifd.entryByName['rvalCtr2']['widget'].get()),
                float(self.ifd.entryByName['gvalCtr2']['widget'].get()),
                float(self.ifd.entryByName['bvalCtr2']['widget'].get()),)
        return self.doitWrapper(firstAts, secondAts, display, cutoff, mat1, mat2)


    def hideDSelector(self, event=None):
        e = self.ifd.entryByName['keyDSelector']
        if self.hideDSel.get():
            e['widget'].grid_forget()
        else:
            e['widget'].grid(e['gridcfg'])
        self.form.autoSize() 


    def hideASelector(self, event=None):
        e = self.ifd.entryByName['keyASelector']
        if self.hideASel.get():
            e['widget'].grid_forget()
        else:
            e['widget'].grid(e['gridcfg'])
        self.form.autoSize() 


    def hideDisplay(self, event=None):
        names = ['cutoff', 'colorLab', 'firstColorLab', 'rvalCtr1', 'gvalCtr1', 'bvalCtr1', \
                    'secondColorLab','rvalCtr2', 'gvalCtr2', 'bvalCtr2'] 
        for name in names:
            e = self.ifd.entryByName[name]
            if not self.displayClosePairsVar.get():
                e['widget'].grid_forget()
                self.first_lines.Set(visible=0)
                self.second_lines.Set(visible=0)
            else:
                apply(e['widget'].grid, (), e['gridcfg'])
                self.first_lines.Set(visible=1)
                self.second_lines.Set(visible=1)
        self.form.autoSize() 
        self.vf.GUI.VIEWER.Redraw()


    def __call__(self, first, second, display=False, cutoff=8.0, \
                    mat1=(0,1,1), mat2 = (0,1,0), **kw):
        """None<-PyAD_calcClosestPairs(first, second)
first:one set of atoms
second:second set of atoms
        """
        apply(self.doitWrapper, (first, second, display, cutoff, mat1, mat2), kw)


    def doit(self, first, second, display, cutoff, mat1, mat2):
        firstAts = self.vf.expandNodes(first)
        secondAts = self.vf.expandNodes(second)
        ms = self.ms = MolecularSystem()
        ms.add_entities(firstAts) 
        ms.add_entities(secondAts) 
        dist = ms.get_dist_mat(0,1)
        #label each first atom by distance to closest second atom
        distances = map(min, dist) #this returns list of lowest vals
        #warning: min distance is always from an atom to itself...
        if display:
            vertices = []
            faces = []
            ct = 0
        for i in range(len(dist)):
            a = firstAts[i]
            min_dist = a.min_dist = distances[i]
            min_index = a.min_index = dist[i].index(min_dist)
            b = secondAts[min_index]
            a.closest_at = b
            if display:
                if min_dist < cutoff:
                    vertices.extend([a.coords, b.coords])
                    faces.append((ct, ct+1))
                    ct=ct+2
        if display:
            self.firstVerts = vertices
            self.firstFaces = faces
            self.first_lines.Set(vertices=vertices, faces=faces, materials=(mat1,))
            self.first_lines.Set(visible=True)
            vertices = []
            faces = []
            ct = 0
        #label each second atom by distance to closest first atom
        swap_dist = Numeric.swapaxes(dist, 0,1)
        column_distances = map(min, swap_dist)
        #now do the second ones
        for i in range(len(swap_dist)):
            b = secondAts[i]
            min_dist = b.min_dist = column_distances[i]
            min_index = b.min_index = swap_dist[i].tolist().index(min_dist)
            a = firstAts[min_index]
            b.closest_at = a
            if display:
                if min_dist < cutoff:
                    vertices.extend([b.coords, a.coords])
                    faces.append((ct, ct+1))
                    ct=ct+2
        if display:
            self.secondVerts = vertices
            self.secondFaces = faces
            self.second_lines.Set(vertices=vertices, faces=faces, materials=(mat2,))
            self.second_lines.Set(visible=True)
            self.vf.GUI.VIEWER.Redraw()


PyADCalcClosestPairsGUI=CommandGUI()
PyADCalcClosestPairsGUI.addMenuCommand('AutoToolsBar', 'PyAutoDock', 'calculate closest pairs')


class PyADCalcVDWEnergies(MVCommand):
    """For each atom in one AtomSet, determine the vdw energy vs all the atoms in a second
    AtomSet"""

    def __init__(self, func=None):
        MVCommand.__init__(self, func)
        self.flag = self.flag | self.objArgOnly
        self.weight = .1485
        self.scorer = VanDerWaals()
        self.prop = 'vdw_energy'
        self.title = "For VanDerWaals Energies: Specify 'first' and 'second' group of atoms:"
        self.weightLabelTxt = "VDW weight"


    def onAddCmdToViewer(self):
        self.hideDSel = Tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.hideASel = Tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.weightLabel = Tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.weightLabel.set(self.weightLabelTxt)


    def buildForm(self):
        ifd = self.ifd = InputFormDescr(title = "Specify 'first' and 'second' group of atoms:")
        self.useSelection = Tkinter.IntVar(master=self.vf.GUI.ROOT)
        ifd.append({'name':'DkeyLab',
                    'widgetType':Tkinter.Label,
                    'text':'For First atoms:',
                    'gridcfg':{'sticky':'w','columnspan':3}})
        ifd.append({'name':'selDRB0',
            'widgetType':Tkinter.Radiobutton,
            'wcfg':{'text':'Use all atoms',
                'variable': self.hideDSel,
                'value':1,
                'command':self.hideDSelector
                },
            'gridcfg':{'sticky':'w'}})
        ifd.append({'name':'selDRB1',
            'widgetType':Tkinter.Radiobutton,
            'wcfg':{'text':'Set atoms to use',
                'variable': self.hideDSel,
                'value':0,
                'command':self.hideDSelector
                },
            'gridcfg':{'sticky':'w', 'row':-1, 'column':1}})
        ifd.append({'name': 'keyDSelector',
            'wtype':StringSelectorGUI,
            'widgetType':StringSelectorGUI,
            'wcfg':{ 'molSet': self.vf.Mols,
                    'vf': self.vf,
                    'all':1,
                    'crColor':(0.,1.,.2),
            },
            'gridcfg':{'sticky':'we', 'columnspan':2 }})
        #now the second atoms
        ifd.append({'name':'AkeyLab',
                    'widgetType':Tkinter.Label,
                    'text':'For Second atoms:',
                    'gridcfg':{'sticky':'w','columnspan':3}})
        ifd.append({'name':'selARB0',
            'widgetType':Tkinter.Radiobutton,
            'wcfg':{'text':'Use all atoms',
                'variable': self.hideASel,
                'value':1,
                'command':self.hideASelector
                },
            'gridcfg':{'sticky':'w'}})
        ifd.append({'name':'selARB1',
            'widgetType':Tkinter.Radiobutton,
            'wcfg':{'text':'Set atoms to use',
                'variable': self.hideASel,
                'value':0,
                'command':self.hideASelector
                },
            'gridcfg':{'sticky':'w', 'row':-1, 'column':1}})
        ifd.append({'name': 'keyASelector',
            'wtype':StringSelectorGUI,
            'widgetType':StringSelectorGUI,
            'wcfg':{ 'molSet': self.vf.Mols,
                    'vf': self.vf,
                    'all':1,
                    'crColor':(0.,1.,.2),
            },
            'gridcfg':{'sticky':'we', 'columnspan':2 }})
        ifd.append({'name':'weight',
            'wtype':ThumbWheel,
            'widgetType':ThumbWheel,
            'wcfg':{ 'labCfg':{'text':self.weightLabel.get(), 
                        'font':(ensureFontCase('helvetica'),12,'bold')},
                'showLabel':2, 'width':100,
                'min':0.0, 'max':2, 'type':float, 'precision':4,
                'value':self.weight,
                'showLabel':1,
                'continuous':1, 'oneTurn':2, 'wheelPad':2, 'height':20},
            'gridcfg':{'sticky':'we', 'columnspan':2}})
        ifd.append({'widgetType': Tkinter.Button,
            'text':'Ok',
            'wcfg':{'bd':6},
            'gridcfg':{'sticky':'ew', 'columnspan':2 },
            'command':self.Accept_cb})
        ifd.append({'widgetType': Tkinter.Button,
            'text':'Cancel',
            'wcfg':{'bd':6},
            'gridcfg':{'sticky':'ew', 'column':2,'row':-1},
            'command':self.Close_cb})
        self.form = self.vf.getUserInput(self.ifd, modal=0, blocking=0)
        self.hideDSelector()
        self.hideASelector()


    def guiCallback(self):
        if not len(self.vf.Mols):
            self.warningMsg('no molecules in viewer')
            return 
        if not hasattr(self, 'ifd'):
            self.buildForm()
        else:
            self.form.deiconify()


    def Close_cb(self, event=None):
        self.form.withdraw()


    def Accept_cb(self, event=None):
        self.form.withdraw()
        firstNodesToCheck = self.ifd.entryByName['keyDSelector']['widget'].get()
        if not len(firstNodesToCheck):
            self.warningMsg('no first atoms specified')
            firstAts = AtomSet([])
        else:
            firstAts = firstNodesToCheck.findType(Atom)
            if self.useSelection.get():
                curSel = self.vf.getSelection()
                if len(curSel):
                    curAts = curSel.findType(Atom)
                    if curAts != self.vf.allAtoms:
                        firstAts = curAts.inter(firstAts)
                        if firstAts is None:
                            msg = 'no specified first atoms in current selection'
                            self.warningMsg(msg)
                            firstAts = AtomSet([])
                else:
                    msg = 'no current selection'
                    self.warningMsg(msg)
                    firstAts = AtomSet([])
        ###FIX THIS!!!
        secondNodesToCheck = self.ifd.entryByName['keyASelector']['widget'].get()
        if not len(secondNodesToCheck):
            self.warningMsg('no second group of atoms specified')
            secondAts = AtomSet([])
            #return
        else:
            secondAts = secondNodesToCheck.findType(Atom)
            if self.useSelection.get():
                curSel = self.vf.getSelection()
                if len(curSel):
                    curAts = curSel.findType(Atom)
                    if curAts != self.vf.allAtoms:
                        secondAts = curAts.inter(secondAts)
                        if not secondAts:  #NB inter returns empty set!
                            msg = 'no specified second group of atoms in current selection'
                            self.warningMsg(msg)
                            secondAts = AtomSet([])
                else:
                    msg = 'no current selection'
                    self.warningMsg(msg)
                    secondAts = AtomSet([])
        weight = self.ifd.entryByName['weight']['widget'].get()
        self.Close_cb()
        return self.doitWrapper(firstAts, secondAts, weight)


    def hideDSelector(self, event=None):
        e = self.ifd.entryByName['keyDSelector']
        if self.hideDSel.get():
            e['widget'].grid_forget()
        else:
            e['widget'].grid(e['gridcfg'])
        self.form.autoSize() 


    def hideASelector(self, event=None):
        e = self.ifd.entryByName['keyASelector']
        if self.hideASel.get():
            e['widget'].grid_forget()
        else:
            e['widget'].grid(e['gridcfg'])
        self.form.autoSize() 


    def __call__(self, first, second, weight=.1485, **kw):
        """None<-PyAD_calcVDWEnergies(first, second, weight=.1485)
first:one set of atoms
second:second set of atoms
weight: vdw weight to use
        """
        apply(self.doitWrapper, (first, second, weight), kw)


    def doit(self, first, second, weight):
        firstAts = self.vf.expandNodes(first)
        secondAts = self.vf.expandNodes(second)
        ms = self.ms = MolecularSystem()
        ms.add_entities(firstAts) 
        ms.add_entities(secondAts) 
        wmt = WeightedMultiTerm()
        wmt.set_molecular_system(ms)
        wmt.add_term( self.scorer, weight)
        result = wmt.get_score_array()
        #label each first atom by sum of its interaction energies 
        for i in range(len(first)):
            setattr(firstAts[i], self.prop, Numeric.add.reduce(result[i]))
        #label each second atom by sum of its vdw interaction energies
        swap_result = Numeric.swapaxes(result, 0,1)
        for i in range(len(swap_result)):
            setattr(secondAts[i], self.prop, Numeric.add.reduce(swap_result[i]))


PyADCalcVDWEnergiesGUI=CommandGUI()
PyADCalcVDWEnergiesGUI.addMenuCommand('AutoToolsBar', 'PyAutoDock', 'vdW', cascadeName = 'calculate intermolecular energies')



class PyAD4CalcVDWEnergies(PyADCalcVDWEnergies):
    """For each atom in one AtomSet, determine the vdw energy vs all the atoms in a second
    AtomSet"""

    def __init__(self, func=None):
        PyADCalcVDWEnergies.__init__(self, func)
        self.flag = self.flag | self.objArgOnly
        self.weight = .1560
        self.scorer = VanDerWaals()
        self.prop = 'ad4_vdw_energy'
        self.title = "For VanDerWaals Energies: Specify 'first' and 'second' group of atoms:"
        self.weightLabelTxt = "AD4 VDW weight"


    def __call__(self, first, second, weight=.1560, **kw):
        """None<-PyAD4_calcVDWEnergies(first, second, weight=.1560)
first:one set of atoms
second:second set of atoms
weight: vdw weight to use
        """
        apply(self.doitWrapper, (first, second, weight), kw)


PyAD4CalcVDWEnergiesGUI=CommandGUI()
PyAD4CalcVDWEnergiesGUI.addMenuCommand('AutoToolsBar', 'PyAutoDock', 'AD4_vdW', cascadeName = 'calculate intermolecular energies')



class PyADCalcHBONDEnergies(PyADCalcVDWEnergies):
    """For each atom in one AtomSet, determine the hydrogen bonding energy vs all the atoms in a second
    AtomSet"""


    def __init__(self, func=None):
        PyADCalcVDWEnergies.__init__(self, func)
        self.flag = self.flag | self.objArgOnly
        self.weight = .0656
        self.scorer = HydrogenBonding()
        self.prop = 'hb_energy'
        self.title = "For HydrogenBonding Energies: Specify 'first(receptor)' and 'second(ligand)' group of atoms:"
        self.weightLabelTxt = "HydrogenBonding weight"


    def __call__(self, first, second, weight=.0656, **kw):
        """None<-PyAD_calcHBONDEnergies(first, second, weight=.0656)
first:one set of atoms
second:second set of atoms
weight: hbonding weight to use
        """
        apply(self.doitWrapper, (first, second, weight), kw)


PyADCalcHBONDEnergiesGUI=CommandGUI()
PyADCalcHBONDEnergiesGUI.addMenuCommand('AutoToolsBar', 'PyAutoDock', 'hydrogen bonding', cascadeName = 'calculate intermolecular energies')



class PyAD4CalcHBONDEnergies(PyADCalcVDWEnergies):
    """For each atom in one AtomSet, determine the hydrogen bonding energy vs all the atoms in a second
    AtomSet"""


    def __init__(self, func=None):
        PyADCalcVDWEnergies.__init__(self, func)
        self.flag = self.flag | self.objArgOnly
        self.weight = .0974
        self.scorer = NewHydrogenBonding()
        self.prop = 'ad4_hb_energy'
        self.title = "For AD4 HydrogenBonding Energies: Specify 'first(receptor)' and 'second(ligand)' group of atoms:"
        self.weightLabelTxt = "AD4 HydrogenBonding weight"


    def __call__(self, first, second, weight=.0974, **kw):
        """None<-PyAD4_calcHBONDEnergies(first, second, weight=.0974)
first:one set of atoms
second:second set of atoms
weight: hbonding weight to use
        """
        apply(self.doitWrapper, (first, second, weight), kw)


    def doit(self, first, second, weight):
        firstAts = self.vf.expandNodes(first)
        secondAts = self.vf.expandNodes(second)
        ms = self.ms = MolecularSystem()
        ms.add_entities(firstAts) 
        ms.add_entities(secondAts) 
        wmt = WeightedMultiTerm()
        wmt.set_molecular_system(ms)
        wmt.add_term( self.scorer, weight)
        result = wmt.get_score_array()
        #label each first atom by sum of its vdw interaction energies 
        for i in range(len(first)):
            hbondmin = min(result[i])
            hbondmax = max(result[i])
            setattr(firstAts[i], self.prop, hbondmin+hbondmax)
            setattr(firstAts[i], 'hbondmin', hbondmin)
            setattr(firstAts[i], 'hbondmax', hbondmax)
            setattr(firstAts[i], 'sum_hb_energy', Numeric.add.reduce(result[i]))
        #label each second atom by sum of its vdw interaction energies
        swap_result = Numeric.swapaxes(result, 0,1)
        for i in range(len(swap_result)):
            hbondmin = min(swap_result[i])
            hbondmax = max(swap_result[i])
            setattr(secondAts[i], self.prop, hbondmin+hbondmax)
            setattr(secondAts[i], 'hbondmin', hbondmin)
            setattr(secondAts[i], 'hbondmax', hbondmax)
            setattr(secondAts[i], 'sum_hb_energy', Numeric.add.reduce(swap_result[i]))



PyAD4CalcHBONDEnergiesGUI=CommandGUI()
PyAD4CalcHBONDEnergiesGUI.addMenuCommand('AutoToolsBar', 'PyAutoDock', 'AD4_hydrogen bonding', cascadeName = 'calculate intermolecular energies')



class PyADCalcESTATEnergies(PyADCalcVDWEnergies):
    """For each atom in one AtomSet, determine the electrostatics energy vs all the atoms in a second
    AtomSet
    """

    def __init__(self, func=None):
        PyADCalcVDWEnergies.__init__(self, func)
        self.flag = self.flag | self.objArgOnly
        self.weight = .1146
        self.scorer = Electrostatics()
        self.prop = 'estat_energy'
        self.title = "For Electrostatics Energies: Specify 'first' and 'second' group of atoms:"
        self.weightLabelTxt = "Electrostatics weight"


    def __call__(self, first, second, weight=.1146, **kw):
        """None<-PyAD_calcESTATEnergies(first, second, weight=.1146)
first:one set of atoms
second:second set of atoms
weight: electrostatics weight to use
        """
        apply(self.doitWrapper, (first, second, weight), kw)


PyADCalcESTATEnergiesGUI=CommandGUI()
PyADCalcESTATEnergiesGUI.addMenuCommand('AutoToolsBar', 'PyAutoDock', 'electrostatics', cascadeName = 'calculate intermolecular energies')



class PyAD4CalcESTATEnergies(PyADCalcVDWEnergies):
    """For each atom in one AtomSet, determine the electrostatics energy vs all the atoms in a second
    AtomSet
    """


    def __init__(self, func=None):
        PyADCalcVDWEnergies.__init__(self, func)
        self.flag = self.flag | self.objArgOnly
        self.weight = .1465
        self.scorer = Electrostatics()
        self.prop = 'ad4_estat_energy'
        self.title = "For Electrostatics Energies: Specify 'first' and 'second' group of atoms:"
        self.weightLabelTxt = "AD4 Electrostatics weight"


    def __call__(self, first, second, weight=.1465, **kw):
        """None<-PyAD4_calcESTATEnergies(first, second, weight=.1465)
first:one set of atoms
second:second set of atoms
weight: electrostatics weight to use
        """
        apply(self.doitWrapper, (first, second, weight), kw)


PyAD4CalcESTATEnergiesGUI=CommandGUI()
PyAD4CalcESTATEnergiesGUI.addMenuCommand('AutoToolsBar', 'PyAutoDock', 'AD4_electrostatics', cascadeName = 'calculate intermolecular energies')



class PyADCalcDSOLVEnergies(PyADCalcVDWEnergies):
    """For each atom in one AtomSet, determine the hydrogen bonding energy vs all the atoms in a second
    AtomSet"""

    def __init__(self, func=None):
        PyADCalcVDWEnergies.__init__(self, func)
        self.flag = self.flag | self.objArgOnly
        self.weight = .1711
        self.scorer = Desolvation()
        self.prop = 'dsolv_energy'
        self.title = "For Desolvation Energies: Specify 'first(receptor)' and 'second(ligand)' group of atoms:"
        self.weightLabelTxt = "Desolvation weight"


    def __call__(self, first, second, weight=.1711, **kw):
        """None<-PyAD_calcDSOLVEnergies(first, second, DSOLVwt=.1711)
first:one set of atoms
second:second set of atoms
weight: dsolv weight to use
        """
        apply(self.doitWrapper, (first, second, weight), kw)


    def doit(self, first, second, weight):
        firstAts = self.vf.expandNodes(first)
        secondAts = self.vf.expandNodes(second)
        #KLUG here to fix AtSolPar and AtVol for protein
        #
        bothAts = firstAts + secondAts
        for a in bothAts:
            if a.parent.type + '_' + a.name in pep_aromList:
                a.autodock_element=='A'
                a.AtSolPar = .1027
            elif a.autodock_element=='A':
                a.AtSolPar = .1027
            elif a.autodock_element=='C':
                a.AtSolPar = .6844
            else:
                a.AtSolPar = 0.0
        ms = self.ms = MolecularSystem()
        ms.add_entities(firstAts) 
        ms.add_entities(secondAts) 
        wmt = WeightedMultiTerm()
        wmt.set_molecular_system(ms)
        wmt.add_term( self.scorer, weight)
        result = wmt.get_score_array()
        #fix O/H in the first row: arbitrarily set to solvation constant
        for i, n in zip(range(len(result)), result[0]):
            if n==.236 or n==.118:
                #print "1st row: zeroing ", i
                result[0][i]=0

        #label each first atom by its desolvation interaction energies 
        for i in range(len(first)):
            a = firstAts[i]
            if a.element=='O':
                setattr(a, self.prop, .236)   #check this
            elif a.element=='H':
                setattr(a, self.prop, .118)   #check this
            else:
                setattr(a, self.prop, Numeric.add.reduce(result[i]))
        #label each second atom by sum of its vdw interaction energies
        swap_result = Numeric.swapaxes(result, 0,1)
        #label each first atom by its desolvation interaction energies 
        for i in range(len(swap_result)):
            #secondAts[i].vdw_energy = Numeric.add.reduce(swap_result[i])
            a = secondAts[i]
            if a.element=='O':
                setattr(a, self.prop, .236)   #check this
            elif a.element=='H':
                setattr(a, self.prop, .118)   #check this
            else:
                setattr(a, self.prop, Numeric.add.reduce(swap_result[i]))
            #setattr(secondAts[i], self.prop, Numeric.add.reduce(swap_result[i]))


PyADCalcDSOLVEnergiesGUI=CommandGUI()
PyADCalcDSOLVEnergiesGUI.addMenuCommand('AutoToolsBar', 'PyAutoDock', 'desolvation', cascadeName = 'calculate intermolecular energies')



class PyAD4CalcDSOLVEnergies(PyADCalcVDWEnergies):
    """For each atom in one AtomSet, determine the hydrogen bonding energy vs all the atoms in a second
    AtomSet"""

    def __init__(self, func=None):
        PyADCalcVDWEnergies.__init__(self, func)
        self.flag = self.flag | self.objArgOnly
        self.weight = .1159
        self.scorer = NewDesolvation()
        self.prop = 'ad4_dsolv_energy'
        self.title = "For Desolvation Energies: Specify 'first' and 'second' group of atoms:"
        self.weightLabelTxt = "AD4 Desolvation weight"


    def __call__(self, first, second, weight=.1159, **kw):
        """None<-PyAD4_calcDSOLVEnergies(first, second, DSOLVwt=.1711)
first:one set of atoms
second:second set of atoms
weight: dsolv weight to use
        """
        apply(self.doitWrapper, (first, second, weight), kw)


PyAD4CalcDSOLVEnergiesGUI=CommandGUI()
PyAD4CalcDSOLVEnergiesGUI.addMenuCommand('AutoToolsBar', 'PyAutoDock', 
                            'AD4_desolvation', cascadeName = 'calculate intermolecular energies')



class PyADCalcAD3Energies(PyADCalcVDWEnergies):
    """For each atom in one AtomSet, determine the autodock3 energy vs all the atoms in a second
    AtomSet"""

    def __init__(self, func=None):
        PyADCalcVDWEnergies.__init__(self, func)
        self.flag = self.flag | self.objArgOnly
        self.weight = None
        self.weightLabel = None
        self.scorer = AutoDock305Scorer()
        self.prop = 'ad305_energy'
        self.title = "For AutoDock305 Energies: Specify 'first(receptor)' and 'second(ligand)' group of atoms:"
        self.weightLabelTxt = ""


    def onAddCmdToViewer(self):
        self.hideDSel = Tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.hideASel = Tkinter.IntVar(master=self.vf.GUI.ROOT)


    def buildForm(self):
        ifd = self.ifd = InputFormDescr(title = self.title)
        self.hideDSel = Tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.hideASel = Tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.useSelection = Tkinter.IntVar(master=self.vf.GUI.ROOT)
        ifd.append({'name':'DkeyLab',
                    'widgetType':Tkinter.Label,
                    'text':'For First atoms:',
                    'gridcfg':{'sticky':'w','columnspan':3}})
        ifd.append({'name':'selDRB0',
            'widgetType':Tkinter.Radiobutton,
            'wcfg':{'text':'Use all atoms',
                'variable': self.hideDSel,
                'value':1,
                'command':self.hideDSelector
                },
            'gridcfg':{'sticky':'w'}})
        ifd.append({'name':'selDRB1',
            'widgetType':Tkinter.Radiobutton,
            'wcfg':{'text':'Set atoms to use',
                'variable': self.hideDSel,
                'value':0,
                'command':self.hideDSelector
                },
            'gridcfg':{'sticky':'w', 'row':-1, 'column':1}})
        ifd.append({'name': 'keyDSelector',
            'wtype':StringSelectorGUI,
            'widgetType':StringSelectorGUI,
            'wcfg':{ 'molSet': self.vf.Mols,
                    'vf': self.vf,
                    'all':1,
                    'crColor':(0.,1.,.2),
            },
            'gridcfg':{'sticky':'we', 'columnspan':2 }})
        #now the acceptors
        ifd.append({'name':'AkeyLab',
                    'widgetType':Tkinter.Label,
                    'text':'For Second atoms:',
                    'gridcfg':{'sticky':'w','columnspan':3}})
        ifd.append({'name':'selARB0',
            'widgetType':Tkinter.Radiobutton,
            'wcfg':{'text':'Use all atoms',
                'variable': self.hideASel,
                'value':1,
                'command':self.hideASelector
                },
            'gridcfg':{'sticky':'w'}})
        ifd.append({'name':'selARB1',
            'widgetType':Tkinter.Radiobutton,
            'wcfg':{'text':'Set atoms to use',
                'variable': self.hideASel,
                'value':0,
                'command':self.hideASelector
                },
            'gridcfg':{'sticky':'w', 'row':-1, 'column':1}})
        ifd.append({'name': 'keyASelector',
            'wtype':StringSelectorGUI,
            'widgetType':StringSelectorGUI,
            'wcfg':{ 'molSet': self.vf.Mols,
                    'vf': self.vf,
                    'all':1,
                    'crColor':(0.,1.,.2),
            },
            'gridcfg':{'sticky':'we', 'columnspan':2 }})
        ifd.append({'widgetType': Tkinter.Button,
            'text':'Ok',
            'wcfg':{'bd':6},
            'gridcfg':{'sticky':'ew', 'columnspan':2 },
            'command':self.Accept_cb})
        ifd.append({'widgetType': Tkinter.Button,
            'text':'Cancel',
            'wcfg':{'bd':6},
            'gridcfg':{'sticky':'ew', 'column':2,'row':-1},
            'command':self.Close_cb})
        self.form = self.vf.getUserInput(self.ifd, modal=0, blocking=0)
        self.hideDSelector()
        self.hideASelector()


    def Accept_cb(self, event=None):
        self.form.withdraw()
        firstNodesToCheck = self.ifd.entryByName['keyDSelector']['widget'].get()
        if not len(firstNodesToCheck):
            self.warningMsg('no first atoms specified')
            firstAts = AtomSet([])
        else:
            firstAts = firstNodesToCheck.findType(Atom)
            if self.useSelection.get():
                curSel = self.vf.getSelection()
                if len(curSel):
                    curAts = curSel.findType(Atom)
                    if curAts != self.vf.allAtoms:
                        firstAts = curAts.inter(firstAts)
                        if firstAts is None:
                            msg = 'no specified first atoms in current selection'
                            self.warningMsg(msg)
                            firstAts = AtomSet([])
                else:
                    msg = 'no current selection'
                    self.warningMsg(msg)
                    firstAts = AtomSet([])
        secondNodesToCheck = self.ifd.entryByName['keyASelector']['widget'].get()
        if not len(secondNodesToCheck):
            self.warningMsg('no second group of atoms specified')
            secondAts = AtomSet([])
            #return
        else:
            secondAts = secondNodesToCheck.findType(Atom)
            if self.useSelection.get():
                curSel = self.vf.getSelection()
                if len(curSel):
                    curAts = curSel.findType(Atom)
                    if curAts != self.vf.allAtoms:
                        secondAts = curAts.inter(secondAts)
                        if not secondAts:  #NB inter returns empty set!
                            msg = 'no specified second group of atoms in current selection'
                            self.warningMsg(msg)
                            secondAts = AtomSet([])
                else:
                    msg = 'no current selection'
                    self.warningMsg(msg)
                    secondAts = AtomSet([])
        self.Close_cb()
        return self.doitWrapper(firstAts, secondAts)


    def __call__(self, first, second, **kw):
        """None<-PyAD_calcAD3Energies(first, second)
first:one set of atoms
second:second set of atoms
        """
        apply(self.doitWrapper, (first, second,), kw)


    def doit(self, first, second ):
        firstAts = self.vf.expandNodes(first)
        secondAts = self.vf.expandNodes(second)
        bothAts = firstAts + secondAts
        for a in bothAts:
            if a.parent.type + '_' + a.name in pep_aromList:
                a.autodock_element=='A'
                a.AtSolPar = .1027
            elif a.autodock_element=='A':
                a.AtSolPar = .1027
            elif a.autodock_element=='C':
                a.AtSolPar = .6844
            else:
                a.AtSolPar = 0.0
        ms = self.ms = MolecularSystem()
        ms.add_entities(firstAts) 
        ms.add_entities(secondAts) 
        self.scorer.set_molecular_system(ms)
        score_array = self.scorer.get_score_array()
        self.scorer.labels_atoms_w_nrg(score_array)

PyADCalcAD3EnergiesGUI=CommandGUI()
PyADCalcAD3EnergiesGUI.addMenuCommand('AutoToolsBar', 'PyAutoDock', 'AutoDock305', cascadeName = 'calculate intermolecular energies')



class PyAD4CalcAD4Energies(PyADCalcAD3Energies):
    """For each atom in one AtomSet, determine the autodock4 energy vs all the atoms in a second
    AtomSet"""


    def __init__(self, func=None):
        PyADCalcAD3Energies.__init__(self, func)
        self.flag = self.flag | self.objArgOnly
        self.weight = None
        self.weightLabel = None
        self.scorer = AutoDock4Scorer()
        self.prop = 'ad4_energy'
        self.title = "For AutoDock4 Energies: Specify 'first(receptor)' and 'second(ligand)' group of atoms:"


    def __call__(self, first, second, **kw):
        """None<-PyAD4_calcAD4Energies(first, second)
first:one set of atoms
second:second set of atoms
        """
        apply(self.doitWrapper, (first, second,), kw)


    def doit(self, first, second):
        #FIX THIS: some scorers are symmetric and others are not
        firstAts = self.vf.expandNodes(first)
        secondAts = self.vf.expandNodes(second)
        ms = self.ms = MolecularSystem()
        ms.add_entities(firstAts) 
        ms.add_entities(secondAts) 
        self.scorer.set_molecular_system(ms)
        score_array = self.scorer.get_score_array()
        self.scorer.labels_atoms_w_nrg(score_array)
##         #NEED TO CORRECT hbond
##         #result = self.scorer.get_score_array()
##         sal = score_array_list = []
##         for t, w in self.scorer.terms:
##             sal.append(t.get_score_array()*w)
##         hbond_array = sal[1]
##         result = Numeric.add(sal[0], sal[2])
##         result = Numeric.add(result, sal[3])

##         #label each first atom by sum of its vdw interaction energies 
##         for i in range(len(first)):
##             #firstAts[i].vdw_energy = Numeric.add.reduce(result[i])
##             hbond_val =  min(hbond_array[i])+max(hbond_array[i])
##             setattr(firstAts[i], self.prop, Numeric.add.reduce(result[i])+hbond_val)
##         #label each second atom by sum of its vdw interaction energies
##         swap_result = Numeric.swapaxes(result, 0,1)
##         swap_hbond_array = Numeric.swapaxes(hbond_array, 0, 1)
##         for i in range(len(swap_result)):
##             #secondAts[i].vdw_energy = Numeric.add.reduce(swap_result[i])
##             hbond_val =  min(swap_hbond_array[i])+max(swap_hbond_array[i])
##             setattr(secondAts[i], self.prop, Numeric.add.reduce(swap_result[i])+hbond_val)


PyAD4CalcAD4EnergiesGUI=CommandGUI()
PyAD4CalcAD4EnergiesGUI.addMenuCommand('AutoToolsBar', 'PyAutoDock', 'AutoDock4', cascadeName = 'calculate intermolecular energies')



class PyAD4CalcInternalEnergies(PyAD4CalcAD4Energies):
    """For one AtomSet, determine the autodock4 energy of all the non-bond interactions:exclude 1-1, 1-2, 1-3, 1-4 and weedbonds in the same rigid piece.  This scorer uses AutoDock4Scorer for estat, vdw and dsolv energies PLUS a hbond1210 scorer for hbond"""

    def __init__(self, func=None):
        PyAD4CalcAD4Energies.__init__(self, func)
        self.flag = self.flag | self.objArgOnly
        iec = self.scorer = InternalEnergy(exclude_one_four=True, weed_bonds=True)
        ad4s = AutoDock4Scorer()
        self.weight = AutoDockTermWeights4() 
        #add the electrostatics term
        iec.add_term(ad4s.terms[0][0], self.weight.vdw_weight)
        #add a HydrogenBonding1210
        iec.add_term(NewHydrogenBonding12_10(), self.weight.hbond_weight)
        iec.add_term(ad4s.terms[2][0], self.weight.vdw_weight)
        iec.add_term(ad4s.terms[3][0], self.weight.dsolv_weight)
        self.prop = 'ad4_internal_energy'
        self.title = "For AD4 Internal Energies: Specify a group of atoms:"
        self.weightLabelTxt = ""


    def buildForm(self):
        ifd = self.ifd = InputFormDescr(title = self.title)
        self.hideDSel = Tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.hideASel = Tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.useSelection = Tkinter.IntVar(master=self.vf.GUI.ROOT)
        ifd.append({'name':'DkeyLab',
                    'widgetType':Tkinter.Label,
                    'text':'Atoms for internal energy calculation:',
                    'gridcfg':{'sticky':'w','columnspan':3}})
        ifd.append({'name':'selDRB0',
            'widgetType':Tkinter.Radiobutton,
            'wcfg':{'text':'Use all atoms',
                'variable': self.hideDSel,
                'value':1,
                'command':self.hideDSelector
                },
            'gridcfg':{'sticky':'w'}})
        ifd.append({'name':'selDRB1',
            'widgetType':Tkinter.Radiobutton,
            'wcfg':{'text':'Set atoms to use',
                'variable': self.hideDSel,
                'value':0,
                'command':self.hideDSelector
                },
            'gridcfg':{'sticky':'w', 'row':-1, 'column':1}})
        ifd.append({'name': 'keyDSelector',
            'wtype':StringSelectorGUI,
            'widgetType':StringSelectorGUI,
            'wcfg':{ 'molSet': self.vf.Mols,
                    'vf': self.vf,
                    'all':1,
                    'crColor':(0.,1.,.2),
            },
            'gridcfg':{'sticky':'we', 'columnspan':2 }})
        ifd.append({'widgetType': Tkinter.Button,
            'text':'Ok',
            'wcfg':{'bd':6},
            'gridcfg':{'sticky':'ew', 'columnspan':2 },
            'command':self.Accept_cb})
        ifd.append({'widgetType': Tkinter.Button,
            'text':'Cancel',
            'wcfg':{'bd':6},
            'gridcfg':{'sticky':'ew', 'column':2,'row':-1},
            'command':self.Close_cb})
        self.form = self.vf.getUserInput(self.ifd, modal=0, blocking=0)
        self.hideDSelector()



    def Accept_cb(self, event=None):
        self.form.withdraw()
        firstNodesToCheck = self.ifd.entryByName['keyDSelector']['widget'].get()
        if not len(firstNodesToCheck):
            self.warningMsg('no first atoms specified')
            firstAts = AtomSet([])
        else:
            firstAts = firstNodesToCheck.findType(Atom)
            if self.useSelection.get():
                curSel = self.vf.getSelection()
                if len(curSel):
                    curAts = curSel.findType(Atom)
                    if curAts != self.vf.allAtoms:
                        firstAts = curAts.inter(firstAts)
                        if firstAts is None:
                            msg = 'no specified first atoms in current selection'
                            self.warningMsg(msg)
                            firstAts = AtomSet([])
                else:
                    msg = 'no current selection'
                    self.warningMsg(msg)
                    firstAts = AtomSet([])
        self.Close_cb()
        return self.doitWrapper(firstAts)


    def doit(self, first):
        firstAts = self.vf.expandNodes(first)
        ms = self.ms = MolecularSystem()
        ms.add_entities(firstAts) 
        self.scorer.set_molecular_system(ms)
        result = self.scorer.get_score_array()
        for i in range(len(first)):
            setattr(firstAts[i], self.prop, Numeric.add.reduce(result[i]))



PyAD4CalcInternalEnergiesGUI=CommandGUI()
PyAD4CalcInternalEnergiesGUI.addMenuCommand('AutoToolsBar', 'PyAutoDock', 'AD4', cascadeName = 'calculate internal energies') 



class PyAD305CalcInternalEnergies(PyAD4CalcInternalEnergies):
    """For one AtomSet, determine the autodock3 energy of all the non-bond interactions:exclude 1-1, 1-2, 1-3, 1-4 and weedbonds in the same rigid piece.  This scorer uses AutoDock4Scorer for estat, vdw and dsolv energies PLUS a hbond1210 scorer for hbond"""

    def __init__(self, func=None):
        PyAD4CalcInternalEnergies.__init__(self, func)
        self.flag = self.flag | self.objArgOnly
        iec = self.scorer = InternalEnergy(exclude_one_four=True, weed_bonds=True)
        ad3s = AutoDock305Scorer()
        self.weight = AutoDockTermWeights305() 
        #ad305 internal energy comprised of estat, vdw and hb1210
        #add the electrostatics term
        iec.add_term(ad3s.terms[0][0], self.weight.vdw_weight)
        #add a HydrogenBonding1210
        iec.add_term(HydrogenBonding12_10(), self.weight.hbond_weight)
        #add vanDerWaals
        iec.add_term(ad3s.terms[2][0], self.weight.vdw_weight)
        #ad305 does not include desolvation
        #iec.add_term(ad3s.terms[3][0], self.weight.dsolv_weight)
        self.prop = 'ad3_internal_energy'
        self.title = "For AD3 Internal Energies: Specify a group of atoms:"
        self.weightLabelTxt = ""


PyAD305CalcInternalEnergiesGUI=CommandGUI()
PyAD305CalcInternalEnergiesGUI.addMenuCommand('AutoToolsBar', 'PyAutoDock', 'AD305', cascadeName = 'calculate internal energies')





commandList = [
    
    {'name':'PyAD_calcClosestPairs','cmd':PyADCalcClosestPairs(),
        'gui':PyADCalcClosestPairsGUI},
    {'name':'PyAD_calcVDWEnergies','cmd':PyADCalcVDWEnergies(),
        'gui':PyADCalcVDWEnergiesGUI},
    {'name':'PyAD_calcHBONDEnergies','cmd':PyADCalcHBONDEnergies(),
        'gui':PyADCalcHBONDEnergiesGUI},
    {'name':'PyAD_calcESTATEnergies','cmd':PyADCalcESTATEnergies(),
        'gui':PyADCalcESTATEnergiesGUI},
    {'name':'PyAD_calcDSOLVEnergies','cmd':PyADCalcDSOLVEnergies(),
        'gui':PyADCalcDSOLVEnergiesGUI},
    {'name':'PyAD_calcAD3Energies','cmd':PyADCalcAD3Energies(),
        'gui':PyADCalcAD3EnergiesGUI},
    {'name':'PyAD4_calcVDWEnergies','cmd':PyAD4CalcVDWEnergies(),
        'gui':PyAD4CalcVDWEnergiesGUI},
    {'name':'PyAD4_calcHBONDEnergies','cmd':PyAD4CalcHBONDEnergies(),
        'gui':PyAD4CalcHBONDEnergiesGUI},
    {'name':'PyAD4_calcESTATEnergies','cmd':PyAD4CalcESTATEnergies(),
        'gui':PyAD4CalcESTATEnergiesGUI},
    {'name':'PyAD4_calcDSOLVEnergies','cmd':PyAD4CalcDSOLVEnergies(),
        'gui':PyAD4CalcDSOLVEnergiesGUI},
    {'name':'PyAD4_calcAD4Energies','cmd':PyAD4CalcAD4Energies(),
        'gui':PyAD4CalcAD4EnergiesGUI},
    {'name':'PyAD305_calcInternalEnergies','cmd':PyAD305CalcInternalEnergies(),
        'gui':PyAD305CalcInternalEnergiesGUI},
    {'name':'PyAD4_calcInternalEnergies','cmd':PyAD4CalcInternalEnergies(),
        'gui':PyAD4CalcInternalEnergiesGUI},
    ]

def initModule(vf):

    for dict in commandList:
        vf.addCommand(dict['cmd'],dict['name'],dict['gui'])


    if hasattr(vf, 'GUI'):
        for item in vf.GUI.menuBars['AutoToolsBar'].menubuttons.values():
            item.configure(background = 'tan')
            item.configure(underline = '-1')


