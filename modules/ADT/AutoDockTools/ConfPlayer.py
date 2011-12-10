
#############################################################################
#
# Author: Ruth HUEY
#
# Copyright: M. Sanner TSRI 2003
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/ConfPlayer.py,v 1.48.4.1 2011/01/07 19:39:42 rhuey Exp $
#
# $Id: ConfPlayer.py,v 1.48.4.1 2011/01/07 19:39:42 rhuey Exp $
#
#

"""
Graphical User Interface to Play AutoDock Conformations which are results of
docking experiment

"""
from mglutil.gui.BasicWidgets.Tk.player import Player
import Tkinter, Pmw, tkMessageBox
import types, time, os
import numpy.oldnumeric as Numeric
from string import find
from mglutil.util.callback import CallBackFunction
from mglutil.gui.InputForm.Tk.gui import InputFormDescr, InputForm, evalString
from mglutil.util.callback import CallbackManager
from mglutil.gui.BasicWidgets.Tk.thumbwheel import ThumbWheel
from mglutil.util.packageFilePath import findFilePath
from mglutil.util.misc import ensureFontCase
import tkMessageBox
from Pmv.moleculeViewer import EditAtomsEvent
#from ViewerFramework.VF import ModificationEvent
from MolKit.pdbParser import PdbqParser, PdbqtParser
from MolKit.molecule import AtomSet, Molecule
from mglutil.math.rmsd import RMSDCalculator
from Pmv.guiTools import MoleculeChooser



class ConformationPlayer(Player):

    def __init__(self, mol, docking, vf, titleStr=None, sequenceList=None,
                        idList = None, wname='spw', delta=0, form2=1,
                        ask=1, **kw):

        self.mol = mol
        #ask is used when newly built molecule is added to viewer (for testing)
        self.ask = ask
        self.coordSlot = len(mol.allAtoms[0]._coords) - 1
        self.docking = docking
        self.updateBindingSite = True
        self.vf = vf
        self.update(sequenceList, idList)
        for item in ['colorByMolecules','colorByAtomType', 'colorByProperty']:
            if not hasattr(self.vf, item):
                self.vf.loadCommand('colorCommands', item, 'Pmv')
        self.vf.loadModule('hbondCommands')
        try:
            if not self.vf.colorMaps.has_key('rgb256'):
                mod = __import__("ViewerFramework")
                VFPath = mod.__path__[0]
                self.vf.loadColorMap(os.path.join(VFPath, "ColorMaps","rgb256_map.py"), 
                                        topCommand=0) 
        except:
            print 'cannot load rgb256'
        self.startFrame = 0
        self.endFrame = len(self.sequenceList)
        kw['titleStr'] = titleStr
        kw['master'] = vf.GUI.ROOT
        kw['endFrame'] = self.endFrame
        kw['maxFrame'] = self.endFrame
        kw['form2'] = form2
        apply(Player.__init__, (self,), kw)
        try:
            self.form.ifd.entryByName['setanimB']['widget'].grid_forget()
            #print 'withdrew SetAnim button'
            self.form.autoSize()
        except:
            pass
        try:
            self.form.ifd.entryByName['recordB']['widget'].grid_forget()
            self.form.autoSize()
        except:
            pass        
        #Player.__init__self, master=vf.GUI.ROOT, endFrame=self.endFrame, 
            #maxFrame=self.endFrame )

    
    def update(self, sequenceList=None, idList=None):
        #sequenceList is a (possibly) ordered list of conformations
        #print 'in cp update: len(sL)=', len(sequenceList)
        if not sequenceList:
            #FIX THIS:!!!
            self.sequenceList = self.docking.ch.conformations
            #self.sequenceList = self.mol.docking.ch.conformations
        else:
            self.sequenceList = sequenceList
        #NB: played sequences have '0' at beginning
        #check to be sure that endFrame and startFrame are ok for use
        # with this sequenceList (that is, not too large)
        lenSeq = len(self.sequenceList)
        self.maxFrame = lenSeq
        if hasattr(self, 'endFrame'):
            self.endFrame = lenSeq
        if hasattr(self, 'startFrame'):
            if self.startFrame>lenSeq:
                print 'startFrame> len(sequenceList)\nresetting startFrame to 0'
                self.startFrame = 0
        #update the thumbwheels here
        if hasattr(self, 'playModeForm'):
            e = self.playModeForm.descr.entryByName
            endTW = e['endFrameTW']['widget']
            #print 'setting endTW to ', self.endFrame
            endTW.max = self.endFrame
            endTW.set(self.endFrame)
            startTW = e['startFrameTW']['widget']
            startTW.set(self.startFrame)
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
            self.idList = map(str, idList)
        #ALSO: close form3 and clear counter
        #these maynot exist yet:
        if hasattr(self, 'playModeForm'):
            if hasattr(self.playModeForm, 'form3'):
                self.closeform3()
        if hasattr(self, 'form'):
            if hasattr(self.form, 'ent2'):
                newLen = max(map(len, self.idList))
                if newLen>3:
                    self.form.ent2.config(width=newLen)
                self.form.ent2.delete(0,'end')
                #could use startFrame if it is valid here:
                if self.startFrame<=len(self.sequenceList) and self.startFrame>0:
                    next_val = str(self.idList[self.startFrame])
                    self.form.ent2.insert(0, next_val)
                    self.currentFrameIndex = self.startFrame
                    self.applyState(self.startFrame-1)
                else:
                    #print self.startFrame, ": index out of range for ", sequenceList, "; resetting startFrame to  0"
                    self.form.ent2.insert(0, str(self.idList[0]))
                    #this calls applyState with reset flag
                    self.applyState(-1)


    def nextFrame(self, id):
        #FRAMES, currentFrameIndex etc ARE INTEGERS!!!
        strID = self.idList[id]
        if self.hasCounter and self.gui:
            self.form.ent2.delete(0,'end')
            self.form.ent2.insert(0, str(strID))
        self.currentFrameIndex = int(id)
        #self.currentFrameIndex = self.idList.index(str(id))
        self.applyState(self.currentFrameIndex-1)


    #methods to set Frame to a specific frame
    #SetState_cb, GoToStart_cb, GoToEnd_cb, setCurrentFrameIndex
    def SetState_cb(self,  event=None):
        #do nothing if no counter
        if self.hasCounter:
            id = self.form.counter.get()
            #check whether id is in current idList
            if id in self.idList:
                ind = self.idList.index(id)
            #elif hasattr(self, 'orderedIdList') and id in self.orderedIdList:
            #    #try to find it in orderedIdList
            #    ind = self.orderedIdList.index(id)
            else:
                #in this case can only reset to 0
                ind = 0
            self.nextFrame(ind)


    def applyState(self, confInd):
        """None<-applyState(mol, confInd)"""
        mol = self.mol
        ch = self.docking.ch
        clust = self.docking.clusterer
        allAts = mol.allAtoms
        if hasattr(clust, 'usesSubset') and clust.usesSubset:
            allAts = clust.subset
        
        # -1 is key for go back to original
        if int(confInd)==-1:
            mol.allAtoms.setConformation(0)
            conf = None
        else:
            #in this case want to get conformation
            conf = self.sequenceList[confInd]
            ch.set_conformation(conf, self.coordSlot)

        if hasattr(clust, 'rmsToolRef'):
            t = 'rms(ref='+ clust.rmsToolRef +') %8.4f'%(round(clust.rmsTool.computeRMSD(allAts.coords[:]),3))
        else:
            t = 'no rms available'
        if not self.vf.hasGui: return
        if hasattr(self.docking, 'bindingSite'):
            self.vf.displayMSMS(self.docking.ligMol, surfName=[self.docking.surfName], negate=True, topCommand=0)
            gc = self.docking.ligMol.geomContainer
            try:
                del gc.atoms[self.docking.surfName]
                geom = gc.geoms[self.docking.surfName]
                geom.Set(protected=False)
                del gc.geoms[self.docking.surfName]
                self.vf.GUI.VIEWER.RemoveObject(geom)
            except:
                pass
            gc.msms={}
            gc.msmsAtoms = {}
        event = EditAtomsEvent('coords', mol.allAtoms)
        self.vf.dispatchEvent(event)
        self.updateColor()
        if hasattr(self, 'playModeForm'):
            if self.buildHBondVar.get():
                self.buildHBonds()
        if self.updateBindingSite and hasattr(self.docking, 'bindingSite'):
            delattr(self.docking, 'bindingSite')
            cmd = self.vf.ADanalyze_showBindingSite
            try:
                percentCutoff = cmd.ifd.entryByName['vdwPercentCutoff']['widget'].get()
            except:
                percentCutoff = 1.0
            cmd.build(percentCutoff)
        self.vf.GUI.VIEWER.Redraw()
        self.showStats()


    def custom_validate(self, text):
        if not len(text):
            return -1
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
        newval = self.idList.index(text) + factor*self.stepSize
        if newval<0 or newval>(len(self.idList)-1):
            return text
        else:
            return self.idList[newval]


    def updateColor(self,event=None):
        """None<-updateColor(mol, event=None)

        makes current molecule coordinates rmstool refCoords
        """
        #print 'in updateColor'
        confNum = self.form.ent2.get()
        #print 'confNum = ', confNum
        if hasattr(self, 'playModeForm'):
            colorType = self.playModeForm.ent3.get()
        else:
            colorType = 'atom'
        #this is already taken care of before
        at0 = self.mol.allAtoms[0]
        has_elec = hasattr(at0, 'estat_energy')
        has_vdw = hasattr(at0, 'vdw_energy')
        has_tot = hasattr(at0, 'total_energy')
        if colorType==None or colorType=='no change':
            return
        try:
            colormap = self.vf.colorMaps['rgb256']
        except:
            print 'setting colormap failed'
            return
        if confNum=='0': 
            if colorType=='molecule':
                self.vf.colorByMolecules(self.mol.allAtoms, ('lines',), topCommand = 0)
            else:
                self.vf.colorByAtomType(self.mol.allAtoms, ('lines',), topCommand = 0)
        elif colorType=='atom':
            self.vf.colorByAtomType(self.mol.allAtoms, ('lines',), topCommand = 0)
        elif colorType=='molecule':
            self.vf.colorByMolecules(self.mol.allAtoms, ('lines',), topCommand = 0)
        elif colorType=='elec_stat' and has_elec:
            val_list = Numeric.array(self.mol.allAtoms.estat_energy)
            mini = min(val_list)
            maxi = max(val_list)
            self.vf.colorByProperty(self.mol.allAtoms, ('lines',),
                'estat_energy',colormap='rgb256',propertyLevel='Atom',
                mini=mini, maxi=maxi, topCommand=0)
        elif colorType=='vdw' and has_vdw:
            val_list = Numeric.array(self.mol.allAtoms.vdw_energy)
            mini = min(val_list)
            maxi = max(val_list)
            self.vf.colorByProperty(self.mol.allAtoms, ('lines',),
                'vdw_energy', colormap='rgb256',propertyLevel='Atom',
                mini=mini, maxi=maxi, topCommand=0)
        elif colorType=='total' and has_tot:
            val_list = Numeric.array(self.mol.allAtoms.total_energy)
            mini = min(val_list)
            maxi = max(val_list)
            self.vf.colorByProperty(self.mol.allAtoms, ('lines',),
                'total_energy',colormap='rgb256',propertyLevel='Atom',
                mini=mini, maxi=maxi, topCommand=0)


#    def buildForm(self, titleStr):
#        #??FIX THIS:
#        mol = self.mol
#        at0 = mol.allAtoms[0]
#        if not titleStr:
#            titleStr = 'Show ' + mol.name + ' Sequence'
#        self.stop = 0
#        if hasattr(self, 'form'):
#            self.form.deiconify()
#            return
#        maxval = len(self.sequenceList)
#        self.doTorsionsOnly = Tkinter.IntVar()
#        self.rmsVar = Tkinter.StringVar()
#        self.rmsVar.set('rms(ref=0) 0.0000')
#        self.energyVar = Tkinter.StringVar()
#        ifd = mol.ifd2 = InputFormDescr(title=titleStr)
#        ifd.append({'name':'energyLab',
#            'widgetType': Tkinter.Label,
#            'textvariable': self.energyVar,
#            'tooltip':'docking and binding energies of current conf',
#            'wcfg':{'bd':4},
#            'gridcfg':{'sticky':'ew', 'columnspan':2}})
#            #'gridcfg':{'sticky':'ew','row':-1, 'column':1}}),
#        ifd.append({'name':'rmsLab',
#            'widgetType': Tkinter.Label,
#            'tooltip':'rms of current conf vs current rms ref conf',
#            'textvariable': self.rmsVar,
#            'wcfg':{'bd':4},
#            'gridcfg':{'sticky':'ew'}}),
#            #'gridcfg':{'sticky':'ew','row':-1, 'column':1}}),
#        colorTypeList = ['atom','molecule', 'no change']
#        #colorTypeList = ['atom','vdw','elec_stat','total','molecule']
#        ifd.append({'widgetType':Pmw.ComboBox,
#                            'name':'colorType',
#                            'tooltip':'used to set coloring scheme for confs',
#                            'wcfg':{'label_text':'Color by',
#		 		                    'entryfield_value':'atom',
#                                    'labelpos':'ew',
#                                    'listheight':'80',
#                                    'scrolledlist_items': colorTypeList,
#				                    'selectioncommand': self.updateColor,
#                                    },
#                            #'gridcfg':{'sticky':'ew'}})
#                            'gridcfg':{'sticky':'ew','row':-1, 'column':1}}),
#        ifd.append({'widgetType':Pmw.Counter,
#			    'name':'statesCounter',
#			    'required':1,
#                'tooltip':'used to show frames via random access',
#			    'wcfg':{#'labelpos': 'n,
#		 		    #'label_text':'conformation:  ',
#                    'autorepeat':0,
#		 		    'entryfield_value':self.idList[0],
#                    'datatype': self.custom_counter,
#		 		    'entry_width':9,
#		 		    'entryfield_validate': self.custom_validate },
#		 	    'gridcfg':{'sticky':'ew', 'columnspan':2}})
#        #ifd.append({'name':'selectCB',
#        #    'widgetType': Tkinter.Checkbutton,
#        #    'tooltip':"open a list showing ids of current conf sequence",
#        #    'text':'show id list',
#        #    'wcfg':{'bd':4},
#        #    'gridcfg':{'sticky':'ew', 'columnspan':1},
#        #    #'gridcfg':{'sticky':'ew', 'columnspan':2},
#        #    'command': self.showStatesList})
#        ifd.append({'name': 'doTransCB',
#            'widgetType': Tkinter.Checkbutton,
#            'text':'torsions only',
#            'tooltip':"show each conf's torsional transformations only",
#            #'text':'transform torsions only',
#            'variable': self.doTorsionsOnly,
#            'wcfg':{'bd':4},
#            'gridcfg':{'sticky':'ew', 'row':-1,'column':1}})
#            ##'gridcfg':{'sticky':'ew', 'columnspan':2}})
#        ifd.append({'name': 'playB',
#            'widgetType': Tkinter.Button,
#            'tooltip':'play from current sequence to last conf',
#            'text':'Play Sequence',
#            'wcfg':{'bd':4},
#            'gridcfg':{'sticky':'ew','columnspan':1},
#            'command':self.Play_cb})
#        ifd.append({'name': 'playRevB',
#            'widgetType': Tkinter.Button,
#            'tooltip':'play from current sequence to input conf',
#            'text':'Play it Reverse',
#            'wcfg':{'bd':4},
#            'gridcfg':{'sticky':'ew','row':-1, 'column':1},
#            'command':self.PlayRev_cb})
#        ifd.append({'name': 'stopB',
#            'widgetType': Tkinter.Button,
#            'text':'Stop',
#            'tooltip':'stop and reset to input conf',
#            'wcfg':{'bd':4},
#            'gridcfg':{'sticky':'ew'},
#            'command':self.Stop_cb})
#        ifd.append({'name': 'pauseB',
#            'widgetType': Tkinter.Button,
#            'tooltip':'pause at current conf',
#            'text':'Pause',
#            'wcfg':{'bd':4},
#            'gridcfg':{'sticky':'ew','row':-1, 'column':1},
#            'command':self.Pause_cb})
#        ifd.append({'name': 'rmsB',
#            'widgetType': Tkinter.Button,
#            'tooltip':'makes current conf rms ref',
#            'text':'Make rms refcoords',
#            'wcfg':{'bd':4},
#            'gridcfg':{'sticky':'ew'},
#            'command':self.MakeRef_cb})
#        ifd.append({'name': 'buildB',
#            'widgetType': Tkinter.Button,
#            'tooltip':'builds new molecule with current conf coords',
#            'text':'Build',
#            'wcfg':{'bd':4},
#            'gridcfg':{'sticky':'ew','row':-1, 'column':1},
#            'command':self.Build_cb})
#        ifd.append({'name':'selectCB',
#            'widgetType': Tkinter.Checkbutton,
#            'text':'Show IdList',
#            'tooltip':"open a list showing ids of current conf sequence",
#            'wcfg':{'bd':4},
#            'gridcfg':{'sticky':'ew', 'columnspan':1},
#            #'gridcfg':{'sticky':'ew', 'columnspan':2},
#            'command': self.showStatesList})

#        ifd.append({'name': 'closeB',
#            'widgetType': Tkinter.Button,
#            'tooltip':'closes player',
#            'text':'Close',
#            'wcfg':{'bd':4},
#            'gridcfg':{'sticky':'ew', 'columnspan':2},
#            'command':self.Close_cb})
#        form = self.vf.getUserInput(ifd, modal=0,blocking=0)
#        form.root.protocol('WM_DELETE_WINDOW',CallBackFunction(self.Close_cb,mol))

#        form.ifd = ifd
#        ctr = ifd.entryByName['statesCounter']['widget']
#        entF = ctr.component('entryfield')
#        form.ent2 = entF._entryFieldEntry
#        da = ctr.component('downarrow')
#        ua = ctr.component('uparrow')
#        for item in [da,ua]:
#            item.bind('<ButtonPress-1>', self.SetState_cb, '+')
#        form.ent2.bind('<Return>', self.SetState_cb, '+')
#        form.counter = form.ifd.entryByName['statesCounter']['widget']
#        form.showList = form.ifd.entryByName['selectCB']['widget']
#        form.showVar = form.ifd.entryByName['selectCB']['variable']
#        ctr2 = ifd.entryByName['colorType']['widget']
#        entF2 = ctr2.component('entryfield')
#        #all this to get a handle to the Tkinter.Entry buried in the pww
#        form.ent3 = entF2._entryFieldEntry
#        form.ent3.bind('<Return>', self.updateColor, '+')
#        form.ent3.config(width=8)

#        #copied lines:
#        form.playB = form.ifd.entryByName['playB']['widget']
#        form.playRevB = form.ifd.entryByName['playRevB']['widget']
#        #set up link to balloon help which needs to change, also
#        form.playTT = form.ifd.entryByName['playB']['balloon']
#        form.playRevTT = form.ifd.entryByName['playRevB']['balloon']
#        if self.hasCounter:
#            ctr = ifd.entryByName['statesCounter']['widget']
#            entF = ctr.component('entryfield')
#            form.ent2 = entF._entryFieldEntry
#            da = ctr.component('downarrow')
#            ua = ctr.component('uparrow')
#            for item in [da,ua]:
#                item.bind('<ButtonPress-1>', self.SetState_cb, '+')
#            form.ent2.bind('<Return>', self.SetState_cb, '+')
#            form.counter = form.ifd.entryByName['statesCounter']['widget']
#        #print 'returning form'
#        return form


    def MakeRef_cb(self, event=None):
        clusterer = self.docking.clusterer
        rmsTool = clusterer.rmsTool
        if not hasattr(clusterer, 'rmsTool'):
            print 'Not available: Make a clustering first!'
            return
        if clusterer.usesSubset:
            rmsTool.setRefCoords(clusterer.subset.coords[:])
        else:
            rmsTool.setRefCoords(self.mol.allAtoms.coords[:])
        clusterer.rmsToolRef = self.form.ent2.get()
        for conf in self.sequenceList:
            if hasattr(conf, 'clRMS'):
                delattr(conf, 'clRMS')
        #force a redo of the info panel here
        if self.showStatsVar.get()==1:
            self.showStats()


    #methods for changing playMode
    def SetMode_cb(self, event=None):
        #print 'SetMode'
        #playMode options:
        #   0   play once and stop
        #   1   play continuously in 1 direction
        #   2   play once in 2 directions
        #   3   play continuously in 2 directions
        #play framerate is frame/per second
        if not hasattr(self, 'playModeForm'):
            self.showPlayMode = Tkinter.IntVar(master=self.vf.GUI.ROOT)
            self.showFrameParmWidgets = Tkinter.IntVar(master=self.vf.GUI.ROOT)
            self.playModeVar = Tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.playModeVar.set('once and stop')
            self.playModeList=[ 'once and stop', 
                                'continuously in 1 direction',
                                'once in 2 directions', 
                                'continuously in 2 directions']
            self.frameParmsList=[ 'framerateLabel','framerateTW','startFrameLabel', 
                                  'startFrameTW', 'endFrameLabel', 'endFrameTW', 
                                  'stepSizeLabel', 'stepSizeTW']

            self.showStatsVar = Tkinter.IntVar(master=self.vf.GUI.ROOT)
            self.showListVar = Tkinter.IntVar(master=self.vf.GUI.ROOT)
            self.doTorsionsOnly = Tkinter.IntVar(master=self.vf.GUI.ROOT)
            self.buildHBondVar = Tkinter.IntVar(master=self.vf.GUI.ROOT)
            ifd2 = InputFormDescr(title='Set Play Options')    
            ifd2.append({'name':'showStatsCB',
                'widgetType': Tkinter.Checkbutton,
                'tooltip':'show binding energy, docking energy etc\nfor current conf in a separate window',
                'wcfg':{ 'text':'Show Info',
                        'command': self.showStats,
                        'variable': self.showStatsVar,
                       },
                'gridcfg':{'sticky':'we'}})
            ifd2.append({'name':'hbondsCB',
                'widgetType': Tkinter.Checkbutton,
                'tooltip':'build and display hydrogen bonds for current conf',
                'wcfg':{ 'text':'Build  H-bonds',
                        'command': self.buildHBonds,
                        'variable': self.buildHBondVar,
                       },
                'gridcfg':{'sticky':'ew', 'row':-1, 'column':1}})
            #colorTypeList = ['atom', 'molecule','no change']
            colorTypeList = ['atom','vdw','elec_stat','total','molecule']
            ifd2.append({'widgetType':Pmw.ComboBox,
                            'name':'colorType',
                            'tooltip':'used to set coloring scheme for confs',
                            'wcfg':{'label_text':'Color by',
 		                            'entryfield_value':'atom',
                                    'labelpos':'w',
                                    'listheight':'80',
                                    'scrolledlist_items': colorTypeList,
                                    'selectioncommand': self.updateColor,
                                    },
                            'gridcfg':{'sticky':'nesw'}}),
            ifd2.append({'name':'selectCB',
                'widgetType': Tkinter.Checkbutton,
                'tooltip':'show ids of current ordered conformation list',
                'wcfg':{ 'text':'Show Conf List',
                        'command': self.showStatesList,
                        'variable': self.showListVar,
                       },
                'gridcfg':{'sticky':'ew', 'row':-1, 'column':1}})
                #'gridcfg':{'sticky':'ew'}})
            ifd2.append({'name': 'rmsB',
                'widgetType': Tkinter.Button,
                'tooltip':'makes current conf rms reference',
                'text':'Make clust RMS ref',
                'wcfg':{ 'command':self.MakeRef_cb,
                        },
                'gridcfg':{'sticky':'we'}})
            ifd2.append({'name':'setRMSB',
                'widgetType': Tkinter.Button,
                'tooltip':'choose a molecule as reference for rms calculation',
                'wcfg':{ 'text':'Choose mol for RMS ref',
                        'command': self.SetRMSRef_cb,
                       },
                'gridcfg':{'sticky':'nesw', 'row':-1, 'column':1}})
            ifd2.append({'name':'playModeMb',
                'widgetType': Tkinter.Menubutton,
                'tooltip':'set play mode choice',
                'wcfg':{ 'text':'Play Mode',
                       },
                'gridcfg':{'sticky':'we'}})
                #'gridcfg':{'sticky':'w', 'columnspan':2}})
            ifd2.append({'name':'adjustFrameParmsMb',
                'widgetType': Tkinter.Checkbutton,
                'tooltip':'opens panel to set play rate, start conf number, end conf number \nand step size for playing conf sequence',
                'wcfg':{ 'text':'Play Parameters',
                        'command': self.showFrameParms_cb,
                        'variable': self.showFrameParmWidgets,
                       },
                'gridcfg':{'sticky':'w', 'row':-1, 'column':1}})
            ifd2.append( {'name': 'framerateLabel',
                    'widgetType':Tkinter.Label,
                    'wcfg':{'text':'frame rate:',
                        'font':(ensureFontCase('helvetica'),12,'bold')},
                    'gridcfg':{'sticky':'w'}})
            ifd2.append({'name': 'framerateTW',
                    'wtype':ThumbWheel,
                    'widgetType':ThumbWheel,
                    'tooltip':'set max num of confs to be displayed per second',
                    'wcfg':{
                        'labCfg':{
                            'fg':'black',
                            'side':'left',
                            'text':''
                            },
                        'showLabel':1, 'width':100,
                        'min':0,
                        'max':100,
                        'lockBMin':1,
                        'lockBMax':0,
                        'lockBIncrement':1,
                        'value':self.framerate,
                        'oneTurn':100,
                        'type':'float',
                        'increment':.1,
                        'callback':self.setMode_cb,
                        'canvascfg':{'bg':'red'},
                        'wheelLabcfg1':{'font':(ensureFontCase('times'),14,'bold')},
                        'wheelLabcfg2':{'font':(ensureFontCase('times'),14,'bold')},
                        'continuous':1, 'wheelPad':1, 'height':20},
                    'gridcfg':{'sticky':'nesw', 'row':-1, 'column':1}})
            ifd2.append( {'name': 'startFrameLabel',
                    'widgetType':Tkinter.Label,
                    'wcfg':{'text':'start frame:',
                        'font':(ensureFontCase('helvetica'),12,'bold')},
                    'gridcfg':{'sticky':'w'}})
            ifd2.append({'name': 'startFrameTW',
                    'wtype':ThumbWheel,
                    'widgetType':ThumbWheel,
                    'tooltip':'set number of first conf to be displayed',
                    'wcfg':{
                        'labCfg':{
                            'fg':'black',
                            'side':'left',
                            'text':''
                            },
                        'showLabel':1, 'width':100,
                        'min':0,
                        'max':self.endFrame,
                        'lockBMin':0,
                        'lockBMax':1,
                        'lockBIncrement':1,
                        'value':self.startFrame,
                        'oneTurn':10,
                        'type':'int',
                        'increment':1,
                        'callback':self.setMode_cb,
                        'canvascfg':{'bg':'green'},
                        'wheelLabcfg1':{'font':(ensureFontCase('times'),14,'bold')},
                        'wheelLabcfg2':{'font':(ensureFontCase('times'),14,'bold')},
                        'continuous':1, 'wheelPad':1, 'height':20},
                    'gridcfg':{'sticky':'ew', 'row':-1,  'column':1}})
            ifd2.append( {'name': 'endFrameLabel',
                    'widgetType':Tkinter.Label,
                    'wcfg':{'text':'end frame:',
                        'font':(ensureFontCase('helvetica'),12,'bold')},
                    'gridcfg':{'sticky':'w'}})
            ifd2.append({'name': 'endFrameTW',
                    'wtype':ThumbWheel,
                    'widgetType':ThumbWheel,
                    'tooltip':'set number of last conf to be displayed',
                    'wcfg':{
                        'labCfg':{
                            'fg':'black',
                            'side':'left',
                            'text':''
                            },
                        'showLabel':1, 'width':100,
                        'min':self.startFrame,
                        'max':self.maxFrame,
                        'lockBMin':1,
                        'lockBMax':0,
                        'lockBIncrement':1,
                        'value':self.endFrame,
                        'oneTurn':10,
                        'type':'int',
                        'increment':1,
                        'callback':self.setMode_cb,
                        'canvascfg':{'bg':'green'},
                        'wheelLabcfg1':{'font':(ensureFontCase('times'),14,'bold')},
                        'wheelLabcfg2':{'font':(ensureFontCase('times'),14,'bold')},
                        'continuous':1, 'wheelPad':1, 'height':20},
                    'gridcfg':{'sticky':'nesw', 'row':-1, 'column':1}})
            ifd2.append( {'name': 'stepSizeLabel',
                    'widgetType':Tkinter.Label,
                    'wcfg':{'text':'step size:',
                        'font':(ensureFontCase('helvetica'),12,'bold')},
                    'gridcfg':{'sticky':'w'}})
            ifd2.append({'name': 'stepSizeTW',
                    'wtype':ThumbWheel,
                    'widgetType':ThumbWheel,
                    'tooltip':'set step before next conf number: default is 1',
                    'wcfg':{
                        'labCfg':{
                            'fg':'black',
                            'side':'left',
                            'text':''
                            },
                        'showLabel':1, 'width':100,
                        'min':1,
                        'max':1000,
                        'lockBMin':1,
                        'lockBMax':0,
                        'lockBIncrement':1,
                        'value':self.stepSize,
                        'oneTurn':10,
                        'type':'int',
                        'increment':1,
                        'callback':self.setMode_cb,
                        'canvascfg':{'bg':'blue'},
                        'wheelLabcfg1':{'font':(ensureFontCase('times'),14,'bold')},
                        'wheelLabcfg2':{'font':(ensureFontCase('times'),14,'bold')},
                        'continuous':1, 'wheelPad':1, 'height':20},
                    'gridcfg':{'sticky':'nesw', 'row':-1, 'column':1}})
            ifd2.append({'name':'buildB',
                'widgetType': Tkinter.Button,
                'tooltip':'build a new molecule with current conf coords\nand add it to viewer',
                'wcfg':{ 'text':'Build Current',
                        'command': self.Build_cb,
                       },
                'gridcfg':{'sticky':'we'}})
                #'gridcfg':{'sticky':'ew', 'row':-1, 'column':1}})
            ifd2.append({'name':'buildAllB',
                'widgetType': Tkinter.Button,
                'tooltip':'build a new molecule with each conf coords\nand add them to viewer',
                'wcfg':{ 'text':'Build All',
                        'command': self.BuildAll_cb,
                       },
                'gridcfg':{'sticky':'nesw', 'row':-1, 'column':1}})
            ifd2.append({'name':'writeB',
                'widgetType': Tkinter.Button,
                'tooltip':'write a new file with current conf coords',
                'wcfg':{ 'text':'Write Current',
                        'command': self.Write_cb,
                       },
                'gridcfg':{'sticky':'we'}})
            ifd2.append({'name':'writeAllB',
                'widgetType': Tkinter.Button,
                'tooltip':'write a new file with each conf coords',
                'wcfg':{ 'text':'Write All',
                        'command': self.WriteAll_cb,
                       },
                'gridcfg':{'sticky':'nesw', 'row':-1, 'column':1}})
            ifd2.append({'name':'cancelB',
                        'widgetType': Tkinter.Button,
                        'wcfg':{
                            'text': 'Close',
                            'command': self.cancelPlayMode_cb,
                        },
                        'gridcfg':{'sticky':'ew','columnspan':1}}),
            ifd2.append({'name':'writeComplexB',
                'widgetType': Tkinter.Button,
                'tooltip':'write a new file with receptor and current conf coords',
                'wcfg':{ 'text':'Write Complex',
                        'command': self.Write_Complex_cb,
                       },
                'gridcfg':{'sticky':'nesw', 'row':-1, 'column':1}})
            self.playModeForm = InputForm(self.master, self.root,
                        descr = ifd2,
                        modal = 0, blocking = 0)
            self.framerateWidget = ifd2.entryByName['framerateTW']['widget']
            self.startFrameWidget = ifd2.entryByName['startFrameTW']['widget']
            self.endFrameWidget = ifd2.entryByName['endFrameTW']['widget']
            self.stepSizeWidget = ifd2.entryByName['stepSizeTW']['widget']
            self.frameParmCfgs = []
            self.frameParmWidgets = []
            for i in self.frameParmsList:
                ent = ifd2.entryByName[i]
                self.frameParmCfgs.append(ent['gridcfg'])
                self.frameParmWidgets.append(ent['widget'])
            self.playModeMb = ifd2.entryByName['playModeMb']['widget']
            self.playModeMb.bind('<ButtonPress>', self.buildPlayModeMenu, add='+')
            self.showFrameParms_cb()
            #self.showStatsVar = ifd2.entryByName['showStatsCB']['variable']
            self.showStatsList = ifd2.entryByName['showStatsCB']['widget']
            self.showList = ifd2.entryByName['selectCB']['widget']
            #COMMENT OUT THESE LINES 11/4/2004
            #ctr1 = ifd2.entryByName['sortType']['widget']
            #entF1 = ctr1.component('entryfield')
            ##all this to get a handle to the Tkinter.Entry buried in the pww
            #self.playModeForm.ent1 = entF1._entryFieldEntry
            #self.playModeForm.ent1.bind('<Return>', self.updateSort, '+')
            #self.playModeForm.ent1.config(width=8)

            ctr2 = ifd2.entryByName['colorType']['widget']
            entF2 = ctr2.component('entryfield')
            #all this to get a handle to the Tkinter.Entry buried in the pww
            self.playModeForm.ent3 = entF2._entryFieldEntry
            self.playModeForm.ent3.bind('<Return>', self.updateColor, '+')
            self.playModeForm.ent3.config(width=8)
        else:
            #Should the widgets be reset here?
            self.playModeVar.set(self.playModeList[self.playMode])
            self.framerateWidget.set(self.framerate)
            self.startFrameWidget.set(self.startFrame)
            self.endFrameWidget.set(self.endFrame)
            self.stepSizeWidget.set(self.stepSize)
            self.playModeForm.deiconify()
        self.playModeForm.autoSize()
        

    def updateSortList(self, event=None):
        #print 'in updateSortList'
        if not hasattr(self.docking, 'clusterer'):
            print 'self.docking has no clusterer'
            return
        clustVals = self.docking.clusterer.clustering_dict.keys()
        if not len(clustVals):
            print 'self.docking.clusterer has no clusterings'
            return
        #if hasattr(self, 'playModeForm'):
            ##sortList = self.sortList = ['run']
            #lb = self.playModeForm.descr.entryByName['sortType']['widget']._list._listbox
            #entries = lb.get(0, 'end')
            #for cV in clustVals:
                #if not cV in self.sortList:
                    #self.sortList.append(cV)
                    #lb.insert('end', cV)


    def updateSort(self, event=None):
        #print 'in updateSort'
        if not hasattr(self.docking, 'clusterer'):
            print 'docking has no clusterer'
            return
        #update Choose Conformation widget if nec
        #FIX THIS: what if started with cluster and then went to all confs
        if not hasattr(self, 'origSequenceList'):
            self.origSequenceList = self.sequenceList[:][:]
        #for case of going from a cluster to all confs
        if len(self.origSequenceList)<len(self.sequenceList):
            self.origSequenceList = self.sequenceList[:][:]
        if not hasattr(self, 'origIdList'):
            self.origIdList = self.idList[:]
        #for case of going from a cluster to all confs
        if len(self.origIdList)<len(self.idList):
            self.origIdList = self.idList[:]

        updateform3 = 0
        k = self.playModeForm.ent1.get()
        if k=='run':
            self.sequenceList = self.origSequenceList[:][:]
            self.idList = self.origIdList[:][:]
            self.showStatesList()
            return
        cl = self.docking.clusterer
        vals = cl.clustering_dict[float(k)]
        seqList = []
        #nb: ids are all strings!!!
        idList = ['0']
        #what if there are more than 1 clusterings?
        for i in range(len(vals)):
            #add the conformations in order
            l = vals[i]
            seqList.extend(l)
            #add the corresponding ids
            for j in range(len(l)):
                newId = str(i+1)+'_'+str(j+1)
                idList.append(newId)
        newLen = max(map(len, self.idList))
        if newLen>3:
            print 'updating ent2 width'
            self.form.ent2.config(width=newLen)
        self.orderedSequenceList = seqList[:][:]
        self.orderedIdList = idList[:][:]
        self.sequenceList = self.orderedSequenceList[:][:]
        self.idList = self.orderedIdList[:][:]
        self.showStatesList()


    def buildPlayModeMenu(self, event=None):
        mB = self.playModeMb
        keyList = self.playModeList
        if not self.showPlayMode.get():
            #mB.config(bg='white')
            if not hasattr(mB, 'menu'):
                mB.menu = Tkinter.Menu(mB)
                mB['menu'] = mB.menu
            else:
                mB.menu.delete(1, 'end')
            for i in range(len(keyList)):
                mB.menu.add_radiobutton(label=keyList[i], var=self.playModeVar, 
                            value=keyList[i], command=self.setMode_cb)
            self.showPlayMode.set(1)
        else:
            mB.menu.unpost()
            self.showPlayMode.set(0)
            
           

    def buildHBonds(self, event=None):
        if not hasattr(self.docking, 'macroMol'):
            self.vf.warningMsg( 'must read docking macromolecule first')
            self.buildHBondVar.set(0)
            return
        c = self.vf.showHBonds
        if not self.buildHBondVar.get():
            print 'hiding hbonds'
            if hasattr(c, 'ifd') and hasattr(c, 'form') and \
                        c.form.root.winfo_exists():
                c.dismiss_cb()
            return
        #print 'building hbonds ' for ligMol and macroMol
        hbs = self.vf.buildHBonds(self.docking.ligMol, self.docking.macroMol, topCommand=0)
        if hbs=="ERROR":
            msg = "1:unable to build hbonds for docking constructed from %s"%self.docking.dlo_list[0].parser.filename
            self.vf.warningMsg(msg)
            #return
        hbs2 = self.vf.buildHBonds(self.docking.macroMol, self.docking.ligMol, topCommand=0)
        if hbs2=="ERROR":
            msg = "2:unable to build hbonds for docking constructed from %s"%self.docking.dlo_list[0].parser.filename
            self.vf.warningMsg(msg)
            #return
        else:
            if len(hbs2):
                hbs.update(hbs2)
        if self.docking.ch.current_conf:
            conf = self.docking.ch.current_conf
            if len(hbs.keys()): 
                #link repr of these hbonds to current conformation
                ss = '%d hydrogen bonds formed: \n'%(len(hbs.keys()))
                for h, v in hbs.items():
	                ss = ss + '% 24s :  % 24s\n'%(h.full_name(),AtomSet(v).full_name())
                conf.hbstr = ss
            else:
                conf.hbstr = "no hydrogen bonds formed\n"
        if len(hbs.keys()):    #in this case, display them
            if hasattr(c, 'ifd') and  hasattr(c, 'form') and \
                        c.form.root.winfo_exists():
                #in this case have to destroy previous ifd and rebuild form
                #this is twisted...
                c.dismiss_cb()
            self.vf.showHBonds(self.vf.allAtoms, topCommand=0)
            #show hbonds only for ligand and macromol???
            #allAtoms = self.docking.ligMol.allAtoms + \
                        #self.docking.macroMol.allAtoms
            #self.vf.showHBonds(allAtoms)
        

    def showFrameParms_cb(self, event=None):
        if not self.showFrameParmWidgets.get():
            for w in self.frameParmWidgets:
                w.grid_forget()
        else:
            for i in range(len(self.frameParmWidgets)):
                w = self.frameParmWidgets[i]
                cfg = self.frameParmCfgs[i]
                w.grid(cfg)
        self.playModeForm.autoSize()


    def setMode_cb(self, event=None):
        curVal = self.playModeVar.get()
        #print 'setting playMode to ', curVal
        self.playMode = self.playModeList.index(curVal)
        #print 'setting playMode to ', curVal
        self.framerate = round(self.framerateWidget.get(),4)
        #print 'setting self.framerate ->', self.framerate
        self.timestamp= 1./self.framerate
        self.startFrame = self.startFrameWidget.get()
        self.endFrame = self.endFrameWidget.get()
        #print 'set endFrame to', self.endFrame
        self.stepSize = self.stepSizeWidget.get()
        self.playMode = self.playModeList.index(curVal)
        #i think this restarts the player's trip memory
        #that is, hasn't gone in any direction yet
        self.oneDirection = 0


    def setPlayMode_cb(self, event=None):
        self.setMode_cb()
        self.cancelPlayMode_cb()
        self.oneDirection = 0
        

    def SetRMSRef_cb(self, event=None):
        if not len(self.vf.Mols):
            self.vf.warningMsg('no molecules currently in viewer')
            return 'ERROR'
        if hasattr(self, 'chooser') and self.chooser.form.root.winfo_exists():
            return
            
        self.chooser = MoleculeChooser(self.vf,mode = 'extended',
                                       title='Choose RMS Ref Molecule' )
        self.chooser.ipf.append({'name':'Rms Reference Button',
                                 'widgetType':Tkinter.Button,
                                 'text':'Select RMS Ref Molecule',
                                 'wcfg':{'bd':6},
                                 'gridcfg':{'sticky':Tkinter.E+Tkinter.W},
                                 'command': self.rmsRefMolecule_cb})
        self.rmsRefForm = self.chooser.go(modal=0, blocking=0)
        self.rmsRefForm.autoSize()
        

    def rmsRefMolecule_cb(self, event=None):
        mols = self.chooser.getMolSet()
        if mols is not None and len(mols):
            mol = mols[0]
            if len(mol.allAtoms)!=len(self.docking.ligMol.allAtoms):
                self.vf.warningMsg('reference molecule must have same number of atoms as docked ligand molecule')
                return
            self.rmsRef = mol
            self.rmsTool = RMSDCalculator(mol.allAtoms.coords[:])
            print 'set rmsRef to ', mol.name
            #need to attach  this to clusterer + to conformations (?)
            for conf in self.sequenceList:
                if hasattr(conf, 'refRMS'):
                    delattr(conf, 'refRMS')
            self.chooser.done_cb()
        else:
            print 'nothing selected'
           


    def Build_cb(self, event=None):
        #print building current
        """None<-Build_cb(mol, event=None)

        builds new molecule with current coordinates and adds it to the viewer
        """
        #FIRST CHECK THAT THIS HASN'T already been built
        #get the current counter content for name of new molecule
        #w = self.form.ifd.entryByName['statesCounter']['widget']
        numStr = self.form.counter.get()
        #numStr = w.get()
        #remember idList has '0' always added at the beginning for input conf
        if numStr!='0':
            confInd =  self.idList.index(numStr) - 1
        else:
            confInd = 0
        #CHECK THIS IS CORRECT!!!
        conf = self.sequenceList[confInd]
        self.buildConf(conf, numStr)


    def BuildAll_cb(self, event=None):
        """None<-BuildAll_cb(mol, event=None)

        builds new molecule with  coordinates of each conf and adds them to the viewer
        """
        #print building allConfs
        for i in range(len(self.sequenceList)):
            conf = self.sequenceList[i]
            nStr = self.idList[i+1]
            self.buildConf(conf, nStr)

    def buildConf(self, conf, nameStr):
        newname = self.mol.name + '_conf_' + nameStr
        if newname in self.vf.Mols.name:
            msg = newname + ' already in viewer. Not building a second copy'
            self.vf.warningMsg(msg)
            return 'ERROR'
        allLines = self.mol.parser.allLines
        newLines = []
        #see if molecule is in conformation conf
        cur_conf = self.docking.ch.current_conf
        if conf!=cur_conf:
            self.docking.ch.set_conformation(conf)
        #coords = self.mol.chains.residues.atoms.coords
        coords = self.mol.allAtoms.coords
        #c = self.docking.ch.current_conf
        estat_energies = []
        hasESTAT = 0
        vdw_energies = []
        hasVDW = 0
        if hasattr(conf, 'estat_energies'):
            estat_energies = self.docking.ch.current_conf.estat_energies
        if len(estat_energies):
            hasESTAT = 1
        if hasattr(conf, 'vdw_energies'):
            vdw_energies = self.docking.ch.current_conf.vdw_energies
        if len(vdw_energies):
            hasVDW = 1
        ctr = 0
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
            else:
                newLines.append(l)
        if isinstance(self.mol.parser, PdbqParser):
            pdbqParser = PdbqParser() 
        else:
            pdbqParser = PdbqtParser() 
        pdbqParser.allLines = newLines
        filename = pdbqParser.filename = self.mol.parser.filename + '_conf_' + nameStr
        newMol = pdbqParser.parse()[0]          
        newMol.name = newname
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


    def Write_Complex_cb(self, event=None):
        #print 'writing complex...'
        """None<-Write_complex_cb(mol, event=None)
        writes a new file with current coordinates AND receptor
        """
        if not hasattr(self.docking, 'macroMol'):
            self.vf.warningMsg('Please set the Macromolecule for this docking first via:\nAnalyze->Macromolecule->Open\nor\nAnalyze->Macromolecule->Choose')
            return "ERROR"
        file = self.vf.askFileSave(types=[('pdbqt files', '*.pdbqt'), ("all", "*")], 
                        title="write complex of receptor + current conf:")
        if file is not None:
            #step 1: write all the lines from the receptor's parser into the file
            macroFile = self.docking.macroMol.parser.filename
            #print 'macroFile=', macroFile
            optr = open(macroFile)
            #print "receptor ptr=", optr
            lines = optr.readlines()
            optr.close()
            #print "read lines for receptor=",  len(lines)
            #open the new file to write these lines
            fptr = open(file, 'w')
            #print 'opened ', file
            for l in lines:
                fptr.write(l)
            liglines = self.docking.ligMol.parser.allLines
            ctr = 0
            #print 'writing ', len(liglines), ' lines for the ligand atoms'
            for l in liglines:
                if l.find("ATOM")==0 or l.find("HETATM")==0:
                    crds = self.docking.ligMol.allAtoms[ctr].coords
                    rec = "%s%8.3f%8.3f%8.3f%s\n"%(l[:30],crds[0], crds[1], crds[2],l[54:] ) 
                    fptr.write(rec)
                    ctr += 1
            fptr.close()


    def Write_cb(self, event=None):
        #print writing current
        """None<-Write_cb(mol, event=None)
        writes a new file with current coordinates
        """
        file = self.vf.askFileSave(types=[('pdbqt files', '*.pdbqt'), ("all", "*")], 
                        title="write current conf:")
        if file is not None:
            self.docking.write_current_conformation(filename=file)


    def WriteAll_cb(self, event=None):
        """None<-WriteAll_cb(mol, event=None)
        writes a new file with  coordinates of each conf 
        """
        old_conf = self.docking.ch.current_conf
        for conf in self.docking.ch.conformations:
            self.docking.ch.set_conformation(conf)
            self.docking.write_current_conformation()
        self.docking.ch.set_conformation(old_conf)


    def showStats(self, event=None, force_update=False):
        if not hasattr(self, 'playModeForm'):
            return
        if not hasattr(self, 'showStatsVar'):
            return
        if self.showStatsVar.get()==0:
            if hasattr(self, 'statsForm'):
                self.statsForm.withdraw()
            return
        conf = self.docking.ch.current_conf
        if conf==None or self.form.counter.get()=='0':
            t2 = 'Input Conformation'
            ss = '            \n            \n            \n     No Stats     \n            \n\n\n\n'
        #elif len(self.sequenceList)==1:
        #    t2 = 'only one run: no clustering possible'
        #    ss = '            \n            \n            \n     No Stats     \n            \n\n\n\n'
        else:
            ind = self.sequenceList.index(conf)
            id = self.idList[ind+1]
            t =  'Conformation ' + id
            if len(self.sequenceList)==1:
                t2 = 'only one run: no clustering'
            else:
                t2 = t + ' Info\n'
            ss = self.getStatString(conf)
        if not hasattr(self, 'statsForm'):
            self.statsString = Tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.statsString.set(ss)
            ifd9 = InputFormDescr(title=t2)
            ifd9.append({'widgetType':Tkinter.Label,
                'name':'statsLabel',
                'wcfg':{'text':ss,
                        'textvariable':self.statsString},
                'gridcfg':{'sticky':'nsew',
                     'rowspan':10}})
            #ifd9.append({ 'name': 'closeB',
                #'widgetType': Tkinter.Button,
                #'tooltip':'closes stats',
                #'text':'Close',
                #'wcfg':{'bd':4},
                #'gridcfg':{'sticky':'ew'},
                #'command':self.CloseStats_cb})
            self.statsForm = self.vf.getUserInput(ifd9, modal=0, blocking=0)
            self.statsForm.autoSize()
            self.statsString = ifd9.entryByName['statsLabel']['wcfg']['textvariable'] 
        else:
            self.statsString.set(ss) 
            self.statsForm.f.master.master.wm_title(t2)
            self.statsForm.deiconify()
            self.statsForm.autoSize()


    def CloseStats_cb(self, event=None):
        self.statsForm.withdraw()


    def getStatString(self, conf):
        ss = ''
        found = 0
        #if ref has been set to something else, use self.rmsTool
        # get energy range for cluster containing this conf
        if hasattr(self.mol, 'cluSEQ'):
            for clu in self.mol.cluSEQ:
                if conf in clu:
                    found = 1
                    break
        if found:
            energy_range = round((clu[-1].binding_energy - clu[0].binding_energy), 3)
            ss = ss + 'binding energy range in this cluster of %d is %.2f\n'%(len(clu), energy_range)
        for a in ['binding_energy', 'docking_energy', 'ad4_energy', 'ligand_efficiency',
                  'inhib_constant', 'inhib_constant_units', 'intermol_energy',
                  'vdw_hb_desolv_energy', 'electrostatic_energy',
                  'moving_ligand_fixed_receptor',
                  'moving_ligand_moving_receptor', 'total_internal',
                  'ligand_internal', 'receptor_internal', 'torsional_energy',
                  'unbound_energy', 'internal_energy', 'vdw_energy',
                  'estat_energy', 'filename']:
            if hasattr(conf, a) and getattr(conf,a) is not None:
                newS = a + '=' + str(getattr(conf,a)) + '\n'
                ss = ss + newS
        #if ref has been set to something else, use self.rmsTool
        if hasattr(self, 'rmsTool') and not hasattr(conf, 'refRMS'):
            #have to calc cluster rms specially
            if hasattr(conf, 'clRMS'):
                ss = ss + 'clRMS=' + str(conf.clRMS) + '\n'
            else:
                ss = ss + 'clRMS=n/a\n'
            newV = round(self.rmsTool.computeRMSD(conf.getCoords()[:]), 3)
            ss = ss + 'refRMS='+str(newV) + '\n'
            conf.refRMS = newV
        elif hasattr(conf, 'clRMS'):
            ss = ss + 'clRMS=' + str(conf.clRMS) + '\n'
            #if only 1 cluster_dict.key then clRMS isn't changed
            ss = ss + 'refRMS='
            if hasattr(conf, 'refRMS'):
                ss = ss + str(conf.refRMS) + '\n'
            else:
                ss = ss + 'n/a\n'
        else:   #add a new clRMS here using the clusterer's rmsTool
            newV = round(self.docking.clusterer.rmsTool.computeRMSD(conf.getCoords()[:]), 3)
            ss  = ss + 'clRMS='+str(newV) + '\n'  
            conf.clRMS = newV
            #print "set conf.clRMS to ", newV
            #ss = ss + 'refRMS=' + str(conf.refRMS) + '\n'
            ss = ss + 'refRMS='
            if hasattr(conf, 'refRMS'):
                ss = ss + str(conf.refRMS) + '\n'
            else:
                ss = ss + 'n/a\n'
        for a in [ 'rseed1', 'rseed2']:
            if hasattr(conf, a):
                newS = a + '=' + str(getattr(conf,a)) + '\n'
                ss = ss + newS
        if hasattr(conf, 'hbstr'):
            ss = ss + conf.hbstr
        if hasattr(conf, 'ss'):
            print "conf already has ss=", ss
        else:
            conf.ss = ss
        return ss


    def showStatesList(self, event=None):
        if hasattr(self.playModeForm, 'form3'):
            self.closeform3()
        if self.showListVar.get()==0:
            self.showList.config(text = 'Show Conf List')
            return
        else:
            self.showList.config(text = 'Hide Conf List')
        mol = self.mol
        entries = self.idList
        if hasattr(self.form, 'form3'): return
        ifd4 = InputFormDescr(title='Choose Conformation')
        ifd4.append({'widgetType':'ListChooser',
            'name':'stateObjs',
            'entries':entries,
            'wcfg':{'title':'Pick state',
                    'mode':'single'},
            'lbwcfg':{'height':10},
            'gridcfg':{'sticky':'nsew', 'column':100,
                     'rowspan':10}})
        self.playModeForm.form3 = self.vf.getUserInput(ifd4, modal=0, blocking=0)
        self.playModeForm.form3.autoSize()
        self.playModeForm.form3.root.protocol('WM_DELETE_WINDOW',self.playModeForm.form3.root.withdraw)

        self.form.lb = ifd4.entryByName['stateObjs']['widget'].lb
        ifd4.entryByName['stateObjs']['widget'].title.config(text=
            'Double click to set state')
        self.form.lb.bind('<Double-Button-1>', self.showSpecConf)
        self.form.autoSize()


    def showSpecConf(self, event=None):
        lb = self.form.lb
        if len(lb.curselection()):
            #state = lb.get(lb.curselection()[0])
            idStr = lb.get(lb.curselection()[0])
            #NB: -1 because of offset problem
            confInd = self.idList.index(idStr) -1
            self.applyState(confInd)
            self.form.ent2.delete(0,'end')
            self.form.ent2.insert(0, idStr)
            self.vf.GUI.VIEWER.cameras[0].update()


    def closeform3(self, event=None):
        if hasattr(self.playModeForm, 'form3'):
            self.playModeForm.form3.destroy()
            delattr(self.playModeForm, 'form3')


class PopulationPlayer(ConformationPlayer):

    def __init__(self, mol, docking, vf, titleStr=None, sequenceList=None,
                        idList = None, wname='spw', delta=0, form2=1,
                        ask=1, **kw):
        kw['titleStr'] = titleStr
        kw['sequenceList'] = sequenceList
        kw['idList'] = idList
        kw['wname'] = wname
        kw['delta'] = delta
        kw['form2'] = form2
        kw['ask'] = ask
        apply(ConformationPlayer.__init__, (self, mol, docking, vf,),kw) 
    

    def update(self, sequenceList=None, idList=None):
        #sequenceList is a (possibly) ordered list of conformations
        #print 'in cp update: len(sL)=', len(sequenceList)
        if not sequenceList:
            #FIX THIS:!!!
            ph = self.docking.ph
            self.sequenceList = ph.all_populations[ph.current_pop_ind]
            self.sequence_id = 0
        else:
            self.sequenceList = sequenceList
            self.sequence_id = self.docking.ph.all_populations.index(sequenceList)
        #NB: played sequences have '0' at beginning
        #check to be sure that endFrame and startFrame are ok for use
        # with this sequenceList (that is, not too large)
        lenSeq = len(self.sequenceList)
        self.maxFrame = lenSeq
        if hasattr(self, 'endFrame'):
            self.endFrame = lenSeq
        if hasattr(self, 'startFrame'):
            if self.startFrame>lenSeq:
                print 'startFrame> len(sequenceList)\nresetting startFrame to 0'
                self.startFrame = 0
        #update the thumbwheels here
        if hasattr(self, 'playModeForm'):
            self.ww = self.playModeForm.descr.entryByName['selectCB']['widget']
            self.ww.config({'text': 'Show List of Populations'})
            e = self.playModeForm.descr.entryByName
            endTW = e['endFrameTW']['widget']
            #print 'setting endTW to ', self.endFrame
            endTW.max = self.endFrame
            endTW.set(self.endFrame)
            startTW = e['startFrameTW']['widget']
            startTW.set(self.startFrame)
        ##indicies in idList are 1 more than indices in sequenceList
        if not idList:
            #insert 0 for original state
            idL = range(0, len(self.sequenceList) + 1)
            self.idList = map(str, idL)
        else:
            #or it could be explicit
            #SHOULD always be incremented by 1 so 0 can be original coords
            self.idList = map(str, idList)
        #ALSO: close form3 and clear counter
        #these maynot exist yet:
        ##if hasattr(self, 'playModeForm'):
        ##    if hasattr(self.playModeForm, 'form3'):
        ##        self.closeform3()
        if hasattr(self, 'form'):
            if hasattr(self.form, 'ent2'):
                newLen = max(map(len, self.idList))
                if newLen>3:
                    self.form.ent2.config(width=newLen)
                self.form.ent2.delete(0,'end')
                #could use startFrame if it is valid here:
                if self.startFrame<=len(self.sequenceList) and self.startFrame>0:
                    next_val = str(self.idList[self.startFrame])
                    self.form.ent2.insert(0, next_val)
                    self.currentFrameIndex = self.startFrame
                    self.applyState(self.startFrame-1)
                else:
                    #print self.startFrame, ": index out of range for ", sequenceList, "; resetting startFrame to  0"
                    self.form.ent2.insert(0, str(self.idList[0]))
                    #this calls applyState with reset flag
                    self.applyState(-1)


    def applyState(self, confInd):
        """None<-applyState(mol, confInd)"""
        mol = self.mol
        ph = self.docking.ph
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
            ph.set_conformation(conf, self.coordSlot)
                
        if hasattr(clust, 'rmsToolRef'):
            t = 'rms(ref='+ clust.rmsToolRef +') %8.4f'%(round(clust.rmsTool.computeRMSD(allAts.coords[:]),3))
        else:
            t = 'no rms available'
        if not self.vf.hasGui: return
        event = EditAtomsEvent('coords', mol.allAtoms)
        self.vf.dispatchEvent(event)
        #modEvent = ModificationEvent('edit','coords', mol.allAtoms)
        #mol.geomContainer.updateGeoms(modEvent)
        self.updateColor()
        if hasattr(self, 'playModeForm'):
            if self.buildHBondVar.get():
                self.buildHBonds()
        self.vf.GUI.VIEWER.Redraw()
        self.showStats()


    def updateColor(self,event=None):
        """None<-updateColor(mol, event=None)

        makes current molecule coordinates rmstool refCoords
        """
        #print 'in updateColor'
        confNum = self.form.ent2.get()
        #print 'confNum = ', confNum
        if hasattr(self, 'playModeForm'):
            colorType = self.playModeForm.ent3.get()
        else:
            colorType = 'atom'
        #this is already taken care of before
        at0 = self.mol.allAtoms[0]
        has_elec = hasattr(at0, 'estat_energy')
        has_vdw = hasattr(at0, 'vdw_energy')
        has_tot = hasattr(at0, 'total_energy')
        if colorType==None or colorType=='no change':
            return
        try:
            colormap = self.vf.colorMaps['rgb256']
        except:
            print 'setting colormap failed'
            return
        if confNum=='0': 
            if colorType=='molecule':
                self.vf.colorByMolecules(self.mol.allAtoms, ('lines',), topCommand = 0)
            else:
                self.vf.colorByAtomType(self.mol.allAtoms, ('lines',), topCommand = 0)
        elif colorType=='atom':
            self.vf.colorByAtomType(self.mol.allAtoms, ('lines',), topCommand = 0)
        elif colorType=='molecule':
            self.vf.colorByMolecules(self.mol.allAtoms, ('lines',), topCommand = 0)
        elif colorType=='elec_stat' and has_elec:
            val_list = Numeric.array(self.mol.allAtoms.estat_energy)
            mini = min(val_list)
            maxi = max(val_list)
            self.vf.colorByProperty(self.mol.allAtoms, ('lines',),
                'estat_energy',colormap='rgb256',propertyLevel='Atom',
                mini=mini, maxi=maxi, topCommand=0)
        elif colorType=='vdw' and has_vdw:
            val_list = Numeric.array(self.mol.allAtoms.vdw_energy)
            mini = min(val_list)
            maxi = max(val_list)
            self.vf.colorByProperty(self.mol.allAtoms, ('lines',),
                'vdw_energy', colormap='rgb256',propertyLevel='Atom',
                mini=mini, maxi=maxi,topCommand=0)
        elif colorType=='total' and has_tot:
            val_list = Numeric.array(self.mol.allAtoms.total_energy)
            mini = min(val_list)
            maxi = max(val_list)
            self.vf.colorByProperty(self.mol.allAtoms, ('lines',),
                'total_energy',colormap='rgb256',propertyLevel='Atom',
                mini=mini, maxi=maxi,topCommand=0)



    def updateSortList(self, event=None):
        pass

    def updateSort(self, event=None):
        pass



    def buildHBonds(self, event=None):
        pass
        

    def SetRMSRef_cb(self, event=None):
        pass
        

    def rmsRefMolecule_cb(self, event=None):
        pass


    def buildConf(self, conf, nameStr):
        newname = self.mol.name + '_'+str(self.docking.ph.current_pop_ind)+ '_' + str(self.sequence_id) + '_ind_' + nameStr
        if newname in self.vf.Mols.name:
            msg = newname + ' already in viewer. Not building a second copy'
            self.vf.warningMsg(msg)
            return 'ERROR'
        allLines = self.mol.parser.allLines
        newLines = []
        #see if molecule is in conformation conf
        cur_conf = self.docking.ph.current_conf
        if conf!=cur_conf:
            self.docking.ph.set_conformation(conf)
        #coords = self.mol.chains.residues.atoms.coords
        coords = self.mol.allAtoms.coords
        ctr = 0
        for l in allLines:
            if find(l, 'ATOM')==0 or find(l, 'HETA')==0:
                cc = coords[ctr]
                cc = [round(cc[0],3), round(cc[1], 3), round(cc[2],3)]
                restLine = l[54:66]
                newLines.append(l[:30] +'%8.3f%8.3f%8.3f'%(cc[0],cc[1],cc[2])+ \
                                    restLine + l[66:])
                ctr = ctr+1
            else:
                newLines.append(l)
        if isinstance(self.mol.parser, PdbqParser):
            pdbqParser = PdbqParser() 
        else:
            pdbqParser = PdbqtParser() 
        pdbqParser.allLines = newLines
        filename = pdbqParser.filename = self.mol.parser.filename + '_conf_' + nameStr
        newMol = pdbqParser.parse()[0]          
        newMol.name = newname
        newMol = self.vf.addMolecule(newMol, ask=self.ask)
        newMol.colored = 0


    def WriteAll_cb(self, event=None):
        """None<-WriteAll_cb(mol, event=None)
        writes a new file with  coordinates of each conf 
        """
        #print writing allConfs
        old_conf = self.docking.ph.current_conf
        for conf in self.docking.ph.conformations:
            self.docking.ph.set_conformation(conf)
            self.docking.write_current_conformation()
        self.docking.ph.set_conformation(old_conf)


    def showStats(self, event=None, force_update=False):
        pass

    def CloseStats_cb(self, event=None):
        self.statsForm.withdraw()


    def getStatString(self, conf):
        pass


    def showStatesList(self, event=None):
        #this is used to set population
        #if hasattr(self.playModeForm, 'form3'):
        #    self.closeform3()
        if self.showListVar.get()==0:
            self.showList.config(text = 'Show Populations List')
            return
        else:
            self.showList.config(text = 'Hide Populations List')
        mol = self.mol
        entries = self.idList
        #entries = []
        #for i in self.idList:
            ##FIX THIS @@
            ##make the conformation labels 1 based to go with counter
            #entries.append(int(i) + 1)
        #CHECK THIS!!!
        if hasattr(self.form, 'form3'): return
        ifd4 = InputFormDescr(title='Choose Population')
        ifd4.append({'widgetType':'ListChooser',
            'name':'stateObjs',
            'entries':entries,
            'wcfg':{'title':'Pick Population',
                    'mode':'single'},
            'lbwcfg':{'height':10},
            'gridcfg':{'sticky':'nsew', 'column':100,
                     'rowspan':10}})
        self.playModeForm.form3 = self.vf.getUserInput(ifd4, modal=0, blocking=0)
        self.playModeForm.form3.autoSize()
        self.playModeForm.form3.root.protocol('WM_DELETE_WINDOW',self.playModeForm.form3.root.withdraw)

        self.form.lb = ifd4.entryByName['stateObjs']['widget'].lb
        ifd4.entryByName['stateObjs']['widget'].title.config(text=
            'Double click to set population')
        self.form.lb.bind('<Double-Button-1>', self.showSpecPop)
        self.form.autoSize()


    def showSpecPop(self, event=None):
        #set the population here->
        d = self.mol.docking
        lb = self.form.lb
        if len(lb.curselection()):
            idStr = lb.get(lb.curselection()[0])
            #NB: -1 because of offset problem
            confInd = self.idList.index(idStr) -1
            d.ph.set_current_pop(confInd)
            self.update(d.ph.conformations)


    def closeform3(self, event=None):
        if hasattr(self.playModeForm, 'form3'):
            self.playModeForm.form3.destroy()
            delattr(self.playModeForm, 'form3')




