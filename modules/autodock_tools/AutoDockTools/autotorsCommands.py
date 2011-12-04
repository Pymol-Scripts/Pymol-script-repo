#############################################################################
#
# Author: Ruth HUEY, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/autotorsCommands.py,v 1.101.2.5 2009/09/18 18:01:41 rhuey Exp $
#
# $Id: autotorsCommands.py,v 1.101.2.5 2009/09/18 18:01:41 rhuey Exp $
#
#
#
#
#
#
#
"""
This Module facilitates selecting and formatting a ligand for a subsequent 
AutoDock run.  The steps in this process are:

    * The user selects the small molecule from a list of molecules 
already in the moleculeViewer OR as a PDBQ file, a PDB file or
a MOL2 file from a fileBrowser.  

    * The user selects the ROOT atom of the ligand either: 

        o     by picking it or 

        o     by autoroot which sets the root to be the atom in the 
            molecule which has the smallest 'largest sub-tree.'

    * Next the user decides which possible and active torsions he wants 
to disallow, changing them from active to inactive. This is done by picking 
an active 'green' bond which turns it inactive or 'purple'. This is 
reversible. The user can also disallow all peptide backbone torsions and/or 
all torsions of amide bonds.

    * Carbons in cycles can be tested for aromaticity.  If the angle 
between the normals to adjacent atoms in the cycle is less than 7.5 Degrees, 
the cycle is considered aromatic: its carbons are renamed "A.." and their 
element type set to 'A'. (This is for the force-field calculations done 
in AutoDock.) This Module does this conversion reversibly. Also, the user 
is able to select a carbon to convert (reversibly) and he can change the
the value of the aromaticity cut-off.

    * Non-polar hydrogens and lone pairs are merged which means that the charge of 
each is added to its heavy atom and the hydrogen atoms themselves are not written 
in the output file, thus in some sense 'removing' them from the molecule. 
'Fewer' atoms simplifies the AutoDock run.

    * The last function of this Module is to write a file which contains 
the correctly formatted ligand atoms.  The ROOT section of the molecule 
expands from the selected ROOT atom out to include all atoms adjacent to it 
up to the first active torsion.  The active torsions set the position of 
BRANCH key words in the output pdbq file (and their corresponding 
ENDBRANCH  key words). These keywords are nested to set up  a 
Breadth-First Order Traversal.  Autotors also calculates the torsional degrees 
of freedom (TORSDOF) which is the number of possible torsions less the number of 
symmetry-equivalent torsions (such as a bond to a NH3). This key word is the 
last line of the pdbq file. 
"""



from MoleculePreparation import LigandPreparation, AD4LigandPreparation
from DejaVu import viewerConst
from DejaVu.Geom  import Geom
from PyBabel.cycle import RingFinder
from Pmv.mvCommand import MVCommand, MVAtomICOM, MVBondICOM
#from Pmv.selectionCommands import MVSelectAtomCommand, MVSelectCommand
from MolKit import Read
from MolKit.torTree import TorTree
from MolKit.tree import TreeNode, TreeNodeSet
from MolKit.molecule import AtomSet, Atom, BondSet
from MolKit.protein import Protein, Chain, ChainSet
from MolKit.pdbWriter import PdbqWriter, PdbqtWriter

from MolKit.mol2Parser import Mol2Parser

from ViewerFramework.VFCommand import Command, CommandGUI
##  from ViewerFramework.gui import InputFormDescr
from mglutil.gui.InputForm.Tk.gui import InputFormDescr
from Pmv.guiTools import MoleculeChooser
from Pmv.qkollua import q
from SimpleDialog import SimpleDialog

import numpy.oldnumeric as Numeric
import types, Tkinter, math, os, Pmw, tkMessageBox
from string import split, find

#create global geometries for this module
from DejaVu.Spheres import Spheres
rootSph = Spheres(name='autotors_rootSph', materials=((0.,1.,0),),
          shape = (0,3), radii = 0.3, quality = 15, inheritMaterial=0,
          vertices=((0.,0.,0.),), visible=0, pickable=0)
markSph = Spheres(name='autotors_markSph', materials=((0.,1.,0),),
          shape = (0,3), radii = 0.15, quality = 15,inheritMaterial=0,
          vertices=((0.,0.,0.),), visible=0, pickable=0)

def check_autotors_geoms(VFGUI):
    autotors_geoms_list = VFGUI.VIEWER.findGeomsByName('autotors_geoms')
    if autotors_geoms_list==[]:
        autotors_geoms = Geom("autotors_geoms", shape=(0,0), protected=True)
        VFGUI.VIEWER.AddObject(autotors_geoms, parent=VFGUI.miscGeom)
        autotors_geoms_list = [autotors_geoms]
    return autotors_geoms_list[0]


#FIX THIS: this should be the same value as in the source code of autodock
MAXTORS = 32

def enter_cb(event=None):
    try:
        print event.widget
    except:
        pass



#these are the texts on menubuttons, menu entries, etc.
menuText = {}
menuText['AutoTorsMB'] = 'Ligand'
menuText['Input Molecule'] = 'Input'
menuText['Read Molecule'] = 'Open...(AD3)'
menuText['Choose Molecule'] = 'Choose...(AD3)'
menuText['Rigid Molecule'] = 'Open as Rigid...(AD3)'
menuText['Read Molecule4'] = 'Open...'
menuText['Choose Molecule4'] = 'Choose...'
menuText['Rigid Molecule4'] = 'Open as Rigid...'
menuText['Ref Molecule'] = 'Write Reference File'

menuText['DefineRigidRootMB'] = 'Torsion Tree'
menuText['ByPicking'] = 'Choose Root...'
menuText['PickChain'] = 'Add a chain to root'
menuText['UnPickChain'] = 'Remove a chain from root'
menuText['Automatically'] = 'Detect Root...'
menuText['SRA1'] = 'Show Root Expansion'
menuText['SRA2'] = 'Hide Root Expansion'
menuText['ShowRootAtoms'] = menuText['SRA1']
menuText['ShowAutotorsRootSphMB'] = 'Show/Hide Root Marker'

menuText['DefineRotatableBondsMB'] = 'Rotatable Bonds'
menuText['DefineRotatableBonds'] = 'Choose Torsions...'
menuText['SetTorsionNumber'] = 'Set Number of Torsions...'

menuText['MActive1'] = 'Make all rotatable bonds non-rotatable'
menuText['MActive2'] = 'Make all rotatable bonds rotatable'
menuText['MAmide1'] = 'Make amide bonds non-rotatable'
menuText['MAmide2'] = 'Make amide bonds rotatable'
menuText['MGuan1'] = 'Make guanidinium bonds non-rotatable'
menuText['MGuan2'] = 'Make guanidinium bonds rotatable'
menuText['MPeptide1'] = 'Make peptide backbone bonds non-rotatable'
menuText['MPeptide2'] = 'Make peptide backbone bonds rotatable'
menuText['MAmide'] = menuText['MAmide1']
menuText['MPeptide'] = menuText['MPeptide1']
menuText['MSelected1'] = 'Make bonds between selected atoms non-rotatable'
menuText['MSelected2'] = 'Make rotatable bonds between selected atoms rotatable'
menuText['MSelected'] = menuText['MSelected1']

menuText['AromaticCarbonsMB'] = 'Aromatic Carbons'
menuText['RenameAromaticCarbons'] = 'Rename (C>A)'
menuText['RestoreAliphaticCarbons'] = 'Restore (A>C)'
menuText['SetCarbonNames'] = 'Set Names...'
menuText['ChangeAromaticityCriteria'] = 'Aromaticity Criterion...'

menuText['NonPolarHydrogensMB'] = 'Merge NonPolar Hydrogens'
menuText['Merge'] = 'Merge'
menuText['Restore'] = 'Restore'

menuText['WriteMB'] = 'Output'
menuText['WritePDBQMB'] = 'Save as PDBQ...(AD3)'
menuText['WritePDBQTMB'] = 'Save as PDBQT...'


menuText['AutomaticAutotorsSetupMB'] = 'Quick Setup...'

#other text strings
#label in define rotatable bonds gui
menuText['torStr1'] = 'Number of rotatable bonds = '
menuText['torStr2'] = ' / ' + str(MAXTORS) + '\n'


#these are the warning msgs:
warningText= {}

warningText['noMolecule'] = 'Sorry, you need to read in a molecule first.\n\nUse the menu\n\n"%s:\n %s\n  :%s".' %(menuText['AutoTorsMB'],menuText['Input Molecule'], menuText['Read Molecule4'])

warningText['noAtorsMol'] = 'Sorry, you need to read or choose a molecule first.\n\nUse either:\n\n%s:\n    %s:\n      %s\n  or\n      %s' % (menuText['AutoTorsMB'],menuText['Input Molecule'], menuText['Read Molecule4'], menuText['Choose Molecule4'])

charge_errorfile = None

def checkMolCharges(mol, vf):
    #charge_errorfile must be opened to append
    if not hasattr(vf, 'checkResCharges'):
        vf.loadCommand('editCommands',['checkResCharges'])
    totalCharge,resList = vf.checkResCharges(mol, topCommand=0)
    errCharge = round(totalCharge,0) - totalCharge
    #if errCharge < -.005 or errCharge > .005:
    if errCharge < -.07 or errCharge > .07:
        msg = 'Non-integral charge on '+ mol.name +': '+ \
                    str(totalCharge) + '\n\n'
        lenres = len(resList)
        if lenres:
            msg = msg + 'correct %d residues:\n' %len(resList)
            truncated = 0
            if lenres>15:
                truncated=1
                msg= msg + '(truncated at 15..)\n'
                resList = resList[:15]
            ss = ''
            for item in resList:
                ss = ss + item.name + '     '+ str(item.err)+'\n'
            msg = msg + ss
            if truncated:
                msg= msg + '....\n'
        msg = msg + '\nCharges should be corrected in written output file!'
        vf.warningMsg(msg)
    if charge_errorfile is not None: charge_errorfile.write(msg)
    return errCharge, resList


autoMergeNPHS = 1

def set_autoMergeNPHS(name, oldval, newval):
    ##FIX THIS: how to let autotors see it?
    global autoMergeNPHS
    autoMergeNPHS = newval
    ##print 'in set_autoMergeNPHS', autoMergeNPHS


def initLPO(mol, mode='interactive',repairs="", root=0, outputfilename=None,
                    cleanup='nphs_lps'):
    hs = AtomSet(filter(lambda x:x.element=='H', mol.allAtoms))
    #do NOT add hydrogens here
    repairs = 'bonds'
    #if len(hs):
        ##this should be a userpreference... and (?) exposed in gui??
        #repairs = 'bonds'
    chargeType = mol.allAtoms[0].chargeSet
    num_zero_charge = len(mol.allAtoms.get(lambda x: hasattr(x, 'charge') and x.charge==0))
    if num_zero_charge==len(mol.allAtoms):
        print 'forcing addition of gasteiger charges to molecule with all zero charges'
        charges_to_add = 'gasteiger'
    elif chargeType!=None:
        charges_to_add = None
    else:
        charges_to_add = 'gasteiger'
    #FIX THIS: set from user preference or something
    #cleanup = 'nphs_lps'
    #FIX THIS: set from user preference or something
    #allowd_bonds = None
    #may never need autoRoot so don't automatically determine it

    mol.LPO = LigandPreparation(mol, mode, charges_to_add=charges_to_add, 
                repairs=repairs, root=root, outputfilename=outputfilename,
                cleanup=cleanup)




def initLPO4(mol, mode='interactive',repairs="", root=0, outputfilename=None,
                    cleanup='nphs_lps'):
    hs = AtomSet(filter(lambda x:x.element=='H', mol.allAtoms))
    #do NOT add hydrogens here
    repairs = 'bonds'
    #if len(hs):
        ##this should be a userpreference... and (?) exposed in gui??
        #repairs = 'bonds'
    chargeType = mol.allAtoms[0].chargeSet
    num_zero_charge = len(mol.allAtoms.get(lambda x: hasattr(x, 'charge') and x.charge==0))
    if num_zero_charge==len(mol.allAtoms):
        charges_to_add = 'gasteiger'
        print 'forcing addition of gasteiger charges to molecule with all zero charges'
    elif chargeType!=None:
        charges_to_add = None
    else:
        charges_to_add = 'gasteiger'
    #FIX THIS: set from user preference or something
    #cleanup = 'nphs_lps'
    #FIX THIS: set from user preference or something
    #allowd_bonds = None
    #may never need autoRoot so don't automatically determine it
    mol.LPO = AD4LigandPreparation(mol, mode, repairs, charges_to_add,
                root=root, outputfilename=outputfilename, cleanup=cleanup)

class AdtSetMode(MVCommand):

    def __init__(self):
        MVCommand.__init__(self)
        self.levels=["AD4.2","AD4.0","AD3.05"]
        self.levelColors = {
            'AD4.2':'DarkGreen', #'#19B219'
            'AD4.1':'DarkGreen', #'#19B219'
            'AD4.0':'DarkCyan',
            'AD3.05':'DarkBlue',
            }
        self.levelBarNames = {
            'AD4.2':'AutoTools41Bar',
            'AD4.1':'AutoTools41Bar',
            'AD4.0':'AutoTools4Bar',
            'AD3.05':'AutoTools3Bar',
            }


    def onAddCmdToViewer(self):
        self.first = 1
        import Tkinter, types
        from AutoDockTools import autotors4Commands, autotors41Commands, autotors3Commands
        from AutoDockTools import autogpf4Commands, autogpf41Commands, autogpf3Commands
        from AutoDockTools import autodpf4Commands, autodpf41Commands, autodpf3Commands
        from AutoDockTools import autoflex4Commands, autoflex41Commands
        from AutoDockTools import autostart4Commands, autostart41Commands, autostart3Commands
        from AutoDockTools import autoanalyze4Commands, autoanalyze41Commands, autoanalyze3Commands

        self.cmdDict = {}
        self.cmdDict['AD3.05'] = ['autotors3Commands', 'autogpf3Commands', 'autodpf3Commands',\
                                        'autostart3Commands', 'autoanalyze3Commands']
        self.cmdDict['AD4.0'] = ['autotors4Commands', 'autoflex4Commands', 'autogpf4Commands',\
                                 'autodpf4Commands','autostart4Commands', 'autoanalyze4Commands']
        self.cmdDict['AD4.2'] = ['autotors41Commands', 'autoflex41Commands', 'autogpf41Commands', \
                                 'autodpf41Commands', 'autostart41Commands', 'autoanalyze41Commands']

        self.barNames = {
            'AD4.2': 'AutoTools41Bar',
            'AD4.0': 'AutoTools4Bar',
            'AD3.05': 'AutoTools3Bar',
            }

        self.frameNames = {
            'AD4.2': 'adt41Frame',
            'AD4.0': 'adt4Frame',
            'AD3.05': 'adt3Frame',
            }

        self.modeLabelNames = {
            'AD4.2': 'adt41ModeLabel',
            'AD4.0': 'adt4ModeLabel',
            'AD3.05': 'adt3ModeLabel',
            }
        self.gpfSetGrid = {
            }


        if self.vf.hasGui:
            import Tkinter
            self.modeVar = Tkinter.StringVar(master=self.vf.GUI.ROOT)
            #for the moment, the default should be 4.0
            self.modeVar.set(self.levels[1])
            self.oldModeVar = Tkinter.StringVar(master=self.vf.GUI.ROOT)
            if 'AutoTools4Bar' in self.vf.GUI.menuBars.keys():
                self.oldModeVar.set('AD4.0')
            elif 'AutoTools41Bar' in self.vf.GUI.menuBars.keys():
                self.oldModeVar.set('AD4.2')
            elif 'AutoTools3Bar' in self.vf.GUI.menuBars.keys():
                self.oldModeVar.set('AD3.05')
            if hasattr(self.vf.GUI, 'adt4ModeLabel'):
                self.vf.GUI.adt4ModeLabel.bind("<Double-Button-1>", self.guiCallback)
    
    
    def Close_cb(self, event=None):
        self.form.withdraw()


    def guiCallback(self, event=None):
        if not hasattr(self, 'ifd'):
            self.buildForm()
        else:
            self.form.deiconify()


    def setMode_cb(self, event=None):
        self.doitWrapper(self.modeVar.get())


    def __call__(self, ModeStr, **kw):
        """ADMode <- ADTSetMode(ModeStr, **kw)
        set the current AutoDock MODE and activate commands for this level and inactivate others.
        """
        return apply( self.doitWrapper, (ModeStr,), kw)


    def doit(self, ModeStr, **kw):
        if type(ModeStr)!=types.StringType:
            return "ERROR"
        if ModeStr not in self.levels:
            msg = ModeStr + "string does not map to a valid level"
            self.warningMsg(msg)
            return "ERROR"
        oldModeStr = self.oldModeVar.get()
        if ModeStr==oldModeStr:
            if self.first: 
                self.updateCmds(ModeStr)
                msg = "updateCmds, then return"
            else:
                msg = ModeStr + " already in use so return"
                self.warningMsg(msg)
            return "ERROR"
        #hide old toolbar
        if oldModeStr!='':
            oldBarName = self.levelBarNames[oldModeStr]
            if self.vf.GUI.menuBars.has_key(oldBarName):
                self.vf.GUI.menuBars[oldBarName].pack_forget()
        #pack new toolbar
        barName = self.levelBarNames[ModeStr]
        if self.vf.GUI.menuBars.has_key(barName):
            self.vf.GUI.menuBars[barName].pack(fill='x',expand=1)
        else:
            #load its cmds and pack it here
            for modName in self.cmdDict[ModeStr]:
                self.vf.browseCommands(modName, commands=None, package='AutoDockTools')
        self.vf.GUI.currentADTBar = barName
        self.vf.GUI.menuBars[barName]._frame.master.config({'bg':'tan','height':25,'relief':'flat'})
        import Tkinter
        col = self.levelColors[ModeStr]
        frame = self.frameNames[ModeStr]
        if not hasattr(self.vf.GUI, frame):
            setattr(self.vf.GUI, frame, self.vf.GUI.menuBars[barName].menubuttons.values()[0].master)
        frameInst = getattr(self.vf.GUI, frame)
        modeLabelName = self.modeLabelNames[ModeStr]
        if not hasattr(self.vf.GUI, modeLabelName):
            setattr(self.vf.GUI, modeLabelName,  Tkinter.Label(frameInst, text=ModeStr, width=len(ModeStr),
                                 relief='sunken', borderwidth=1, fg='DarkGreen',
                                 bg = 'ivory',anchor='w' ))
            getattr(self.vf.GUI, modeLabelName).pack(side='left')
            getattr(self.vf.GUI, modeLabelName).bind("<Double-Button-1>", self.vf.ADTSetMode.guiCallback)
        else:
            getattr(self.vf.GUI, modeLabelName).configure(bg='ivory', fg=col, text=ModeStr,
                                            font='Helvetica 10 bold',
                                            width=len(ModeStr))
            getattr(self.vf.GUI, modeLabelName).bind("<Double-Button-1>", self.vf.ADTSetMode.guiCallback)
        self.oldModeVar.set(ModeStr)
        self.updateCmds(ModeStr)
        self.first = 0
        self.modeVar.set(ModeStr)
        

    def updateCmds(self, ModeStr ):
        #remember the current modeStr as 'ADmode'
        self.vf.ADmode = ModeStr
        if ModeStr=="AD4.2":
            self.vf.ADgpf4_setAtomTypes = self.vf.AD41gpf_setAtomTypes
            self.vf.ADtors_defineRotBonds = self.vf.AD41tors_defineRotBonds
            self.vf.ADtors_markRoot = self.vf.AD41tors_markRoot
            self.vf.ADtors_setCarbonNames = self.vf.AD41tors_setCarbonNames
            self.vf.ADgpf_setGrid = self.vf.AD41gpf_setGrid
            self.vf.ADdpf4_initLigand = self.vf.AD41dpf_initLigand
            self.vf.ADflex_setResidues = self.vf.AD41flex_setResidues
            if hasattr(self.vf, 'AD41analyze_showGridIsocontours'):
                self.vf.ADanalyze_showGridIsocontours = self.vf.AD41analyze_showGridIsocontours
            #self.vf.ADanalyze_showHistogram = self.vf.AD41analyze_showHistogram
            self.vf.ADanalyze_showDockingsAsSpheres = self.vf.AD41analyze_showDockingsAsSpheres
            self.vf.ADanalyze_showBindingSite = self.vf.AD41analyze_showBindingSite
            self.vf.ADanalyze_readDLG = self.vf.AD41analyze_readDLG
            self.vf.ADanalyze_selectDLG = self.vf.AD41analyze_selectDLG
            self.vf.ADanalyze_chooseDockedConformations = self.vf.AD41analyze_chooseDockedConformations
            self.vf.ADanalyze_makeSubsetClustering = self.vf.AD41analyze_makeSubsetClustering
        elif ModeStr=="AD3.05":
            self.vf.ADtors_defineRotBonds = self.vf.AD3tors_defineRotBonds
            self.vf.ADtors_markRoot = self.vf.AD3tors_markRoot
            self.vf.ADtors_setCarbonNames = self.vf.AD3tors_setCarbonNames
            self.vf.ADgpf_setGrid = self.vf.AD3gpf_setGrid
            if hasattr(self.vf, 'AD3analyze_showGridIsocontours'):
                self.vf.ADanalyze_showGridIsocontours = self.vf.AD3analyze_showGridIsocontours
            #self.vf.ADanalyze_showHistogram = self.vf.AD3analyze_showHistogram
            self.vf.ADanalyze_showDockingsAsSpheres = self.vf.AD3analyze_showDockingsAsSpheres
            self.vf.ADanalyze_showBindingSite = self.vf.AD3analyze_showBindingSite
            self.vf.ADanalyze_readDLG = self.vf.AD3analyze_readDLG
            self.vf.ADanalyze_selectDLG = self.vf.AD3analyze_selectDLG
            self.vf.ADanalyze_chooseDockedConformations = self.vf.AD3analyze_chooseDockedConformations
            self.vf.ADanalyze_makeSubsetClustering = self.vf.AD3analyze_makeSubsetClustering
        elif ModeStr=="AD4.0":
            from AutoDockTools.WebServices import WebServices, WebServices4GUI
            if not hasattr(self.vf, 'ADweb_services'):
                self.vf.addCommand(WebServices(), 'ADweb_services', WebServices4GUI)
            self.vf.ADgpf4_setAtomTypes = self.vf.AD4gpf_setAtomTypes
            self.vf.ADtors_markRoot = self.vf.AD4tors_markRoot
            self.vf.ADtors_setCarbonNames = self.vf.AD4tors_setCarbonNames
            self.vf.ADtors_defineRotBonds = self.vf.AD4tors_defineRotBonds
            self.vf.ADgpf_setGrid = self.vf.AD4gpf_setGrid
            self.vf.ADdpf4_initLigand = self.vf.AD4dpf_initLigand
            self.vf.ADflex_setResidues = self.vf.AD4flex_setResidues
            if hasattr(self.vf, 'AD4analyze_showGridIsocontours'):
                self.vf.ADanalyze_showGridIsocontours = self.vf.AD4analyze_showGridIsocontours
            #self.vf.ADanalyze_showHistogram = self.vf.AD4analyze_showHistogram
            self.vf.ADanalyze_showDockingsAsSpheres = self.vf.AD4analyze_showDockingsAsSpheres
            self.vf.ADanalyze_showBindingSite = self.vf.AD4analyze_showBindingSite
            self.vf.ADanalyze_readDLG = self.vf.AD4analyze_readDLG
            self.vf.ADanalyze_selectDLG = self.vf.AD4analyze_selectDLG
            self.vf.ADanalyze_chooseDockedConformations = self.vf.AD4analyze_chooseDockedConformations
            self.vf.ADanalyze_makeSubsetClustering = self.vf.AD4analyze_makeSubsetClustering


    def buildForm(self):
        if not hasattr(self, 'ifd'):
            from mglutil.gui.InputForm.Tk.gui import InputFormDescr
            import Tkinter
            ifd = self.ifd = InputFormDescr(title = "AutoDock Mode:")
            levelLabels = ['AutoDock 4.2', 'AutoDock 4.0', 'AutoDock 3.05']
            levels = ['AD4.2', 'AD4.0', 'AD3.05']
            for level, levlabel in zip(levels, levelLabels):
                ifd.append({'name':level, 
                            'widgetType': Tkinter.Radiobutton,
                            'wcfg':{'text':levlabel,
                                    'variable':self.modeVar,
                                    'value':level,
                                    'justify':'left',
                                    'activebackground':self.levelColors[level],
                                    'selectcolor':self.levelColors[level],
                                    'command':self.setMode_cb},
                            'gridcfg':{'sticky':'we'}})
            ifd.append({'name':'dismiss',
                        'widgetType':Tkinter.Button,
                        'defaultValue':1,
                        'wcfg':{'text':'Dismiss',
                                'command':self.Close_cb},
                        'gridcfg':{'sticky':'we'}
                        })
            self.form = self.vf.getUserInput(self.ifd, modal=0, blocking=0)
            self.form.root.protocol('WM_DELETE_WINDOW',self.Close_cb)
        else:
            self.form.deiconify()


class AtorsMoleculeChooser(MVCommand):
    """allows user to choose as ligand  a molecule already in the viewer"""


    def onAddCmdToViewer(self):
        if not hasattr(self.vf, 'atorsDict'):
            self.vf.atorsDict={}
        if not hasattr(self.vf,'readMolecule'):
            self.vf.loadCommand('fileCommands', 'readMolecule', 'Pmv',
                                topCommand=0)


    def __init__(self, mode='single', title = 'Choose Molecule for AutoDock3'):
        MVCommand.__init__(self)
        self.mode = mode
        self.title = title


    def chooseLigand_cb(self, event = None):
        """called each time the 'choose Ligand' button is pressed"""
        mols = self.chooser.getMolSet()
        if mols is not None: 
            self.chooser.form.withdraw()
            self.doitWrapper(mols, redraw=0)
            #self.doitWrapper(mols, log=0, redraw=0)


    def guiCallback(self):
        self.chooser = MoleculeChooser(self.vf, self.mode, self.title)
        self.chooser.ipf.append({'name':'Select Button',
                                 'widgetType':Tkinter.Button,
                                 'text':'Select Autotors Molecule',
                                 'wcfg':{'bd':6},
                                 'gridcfg':{'sticky':'we'},
                                 'command': self.chooseLigand_cb})
        self.form = self.chooser.go(modal=0, blocking=0)
        lb = self.chooser.ipf.entryByName['Molecule']['widget'].lb
        lb.bind("<Double-Button-1>",self.chooseLigand_cb)


    def __call__(self, nodes, **kw):
        """None<-ADtors_chooseLigand(nodes)
nodes:ligand for autodock
        """
        apply(self.doitWrapper, (nodes,), kw)


    def doit(self, nodes, **kw):
        mol = self.vf.expandNodes(nodes)[0]
        filename = mol.parser.filename
        ftype=split(filename,'.')[-1]
        #FIX THIS: could it be something else?
        #do we care???
        if ftype not in ['pdbq', 'pdb', 'pdbqs','pdbqt', 'mol2', 'dlg']:
            msg =  "unknown filetype: "+ ftype
            self.vf.warningMsg(msg)
            #????
            return
        useTorTree=kw.get('useTorTree', 0)
        if not useTorTree and hasattr(mol, 'torTree'):
            msg = mol.name + ' already has a a torsion tree. Do you want to use it to set the activity of rotatable bonds?'
            d = SimpleDialog(self.vf.GUI.ROOT, text=msg, 
                buttons=['No','Yes'], default=1, 
                title='Use Previous Torsion Tree?')
            useTorTree=d.go()
            print 'set useTorTree to', useTorTree

        #include user_preferences here...
        #what about reprocessing?
        self.vf.atorsDict['molecule'] = mol
        if not hasattr(mol, 'LPO') or mol.LPO.version>=4:
            if useTorTree:
                root = mol.ROOT
            cleanup = "nphs_lps"
            if self.vf.userpref['Automerge NPHS']['value']==0:
                cleanup = "lps"
            initLPO(mol, cleanup=cleanup)
            title = "summary for " + mol.name
            self.vf.warningMsg(mol.LPO.summarize(), title=title)
            self.vf.allAtoms = self.vf.Mols.chains.residues.atoms
            #put in useTorTree stuff here:
            if useTorTree:
                self.rebuildTorTree(mol, root)

        if self.vf.hasGui:
            self.vf.colorByAtomType(mol, ['lines'], topCommand=0)
            #nb aromCs could be an empty AtomSet
            self.vf.color(mol.LPO.aromCs,((0,1,0,),), ['lines'], topCommand=0)
            self.vf.centerScene(topCommand=0)
            self.vf.displayLines(mol,topCommand=0)
            self.vf.GUI.VIEWER.Redraw()
            #self.vf.GUI.ligandLabelLabel.config(text='AD3 Ligand:')
            #self.vf.GUI.ligandLabel.config(text=mol.name, width=len(mol.name))


    def rebuildTorTree(self, mol, root):
        #print '2:rebuilding torTree(flexibility pattern) from file'
        #have to rebuild it to capture which bonds are referenced
        #all bonds off to start
        allAts = mol.allAtoms
        for b in allAts.bonds[0]:
            if b.activeTors:
                b.activeTors = 0
        torscount = 0
        tM = mol.torTree.torsionMap
        for i in range(len(tM)):
            bnum0, bnum1 = tM[i].bond
            a0 = allAts.get(lambda x: x.number==bnum0+1)[0]
            #print "a0=", a0.name
            #a0 = allAts[bnum0]
            a0.tt_ind = bnum0
            a1 = allAts.get(lambda x: x.number==bnum1+1)[0]
            #print "a1=", a1.name
            #a1 = allAts[bnum1]
            a1.tt_ind = bnum1
            b = AtomSet([a0,a1]).bonds[0]
            assert b is not None
            if hasattr(b, 'possibleTors'):
                assert b.possibleTors
            else:
                b.possibleTors = 1
            b.activeTors = 1
            torscount = torscount + 1
        #this is also done in AtorsInitMol FIX THIS!!!
        #mol.torscount = torscount
        mol.torscount = len(mol.allAtoms.bonds[0].get(lambda x: x.activeTors==1))
        mol.ROOT = root
        mol.ROOT.rnum0 = 0


    def onPick(self,event):
        listChooser = self.ipf.entryByName['Molecule']['widget']
        tkListBox = listChooser.lb
        atom,geom = self.vf.findPickedAtom(event)
        if atom is not None:
            pickedMol = atom.top
            #then need to make pickedMol the selection in self.lc
            for i in range(len(listChooser.entries)):
                listChooserlist=split(listChooser.entries[i][0])
                if pickedMol.name == listChooserlist[0]:
                    self.pickedMolIndex= i
                    tkListBox.select_clear(0,'end')
                    listChooser.select(i)
                    return
            t= "error: on  %s " %pickedMol.name
            self.vf.warningMsg(t)


AtorsMoleculeChooserGUI=CommandGUI()
AtorsMoleculeChooserGUI.addMenuCommand('AutoToolsBar', menuText['AutoTorsMB'], menuText['Choose Molecule'], cascadeName = menuText['Input Molecule'])



class AtorsReader(MVCommand):
    """allows user to select a file for the ligand via a file browser"""
 

    def onAddCmdToViewer(self):
        if not hasattr(self.vf, 'atorsDict'):
            self.vf.atorsDict={}
        if not hasattr(self.vf, 'readMolecule'):
            self.vf.loadCommand('fileCommands', 'readMolecule', 'Pmv')


    def guiCallback(self):
        """called each time the 'select molecule' button is pressed"""
        molFile = self.vf.askFileOpen(types=[('PDBQ files:', '*.pdbq'),\
                ('PDB files:', '*.pdb'), ('MOL2 files:','*.mol2'),\
                ('all files:', '*')],\
                title = 'Ligand file for AutoDock3:')
        if not molFile: return
        self.doitWrapper(molFile, ask=1, redraw=1)


    def __call__(self, filename, log=1, **kw):
        """None<-ADtors_readLigand(filename)
filename:file to read to get ligand for autodock
        """
        kw['log'] = log
        apply(self.doitWrapper, (filename,),kw)


    def doit(self, filename, **kw):
        #FIX THIS:
        #do we care???
        #could it be something else?
        ftype = split(filename,'.')[-1]
        if ftype not in ['pdbq','pdb','pdbqt', 'pdbqs','mol2']:
            msg = "unknown ligand file type " + ftype
            self.vf.warningMsg(msg)
            return
        mols = self.vf.readMolecule(filename, log=0)
        if not mols:
            return 'ERROR'
        if len(mols)>1:
            msg = str(len(mols)) + ' molecules in ', filename
            self.vf.warningMsg(msg)
            maxAts = 0
            for m in mols:
                numAts = len(m.allAtoms)
                if numAts>maxAts:
                    mol = m
                    maxAts = numAts
        else:
            mol = mols[0]
        if not mol.chains[0].hasBonds: 
            mol.buildBondsByDistance()

        cleanup = "nphs_lps"
        if self.vf.userpref['Automerge NPHS']['value']==0:
            cleanup = "lps"
        initLPO(mol, cleanup=cleanup)
        title = "summary for " + mol.name
        self.vf.warningMsg(mol.LPO.summarize(), title=title)
        #ALWAYS update vf.allAtoms
        self.vf.allAtoms = self.vf.Mols.chains.residues.atoms
        self.vf.atorsDict['molecule'] = mol

        if self.vf.hasGui:
            self.vf.colorByAtomType(mol, ['lines'], topCommand=0)
            #nb aromCs could be an empty AtomSet
            self.vf.color(mol.LPO.aromCs,((0,1,0,),), ['lines'], topCommand=0)
            self.vf.centerScene(topCommand=0)
            self.vf.displayLines(mol,topCommand=0)
            self.vf.GUI.VIEWER.Redraw()
            #self.vf.GUI.ligandLabelLabel.config(text='AD3 Ligand:')
            #self.vf.GUI.ligandLabel.config(text=mol.name, width=len(mol.name))
        #what to do about this stuff?
        #if hasattr(self.vf.ADtors_defineRotBonds, 'noATBut'):
        #    menuText['MAmide'] = menuText['MAmide1']
        #    menuText['MGuan'] = menuText['MGuan1']
        #    menuText['MPeptide'] = menuText['MPeptide1']
        #    menuText['MActive'] = menuText['MActive1']
        #    menuText['MSelected'] = menuText['MSelected1']
        #    self.vf.ADtors_defineRotBonds.noATBut.config(text=menuText['MAmide'])
        #    self.vf.ADtors_defineRotBonds.noGTBut.config(text=menuText['MGuan'])
        #    self.vf.ADtors_defineRotBonds.noPBTBut.config(text=menuText['MPeptide'])
        #    self.vf.ADtors_defineRotBonds.noACTBut.config(text=menuText['MActive'])
        #    self.vf.ADtors_defineRotBonds.noSELBut.config(text=menuText['MSelected'])


AtorsReaderGUI = CommandGUI()
AtorsReaderGUI.addMenuCommand('AutoToolsBar', menuText['AutoTorsMB'], menuText['Read Molecule'], cascadeName = menuText['Input Molecule'])



class Ators4MoleculeChooser(MVCommand):
    """allows user to choose as ligand  a molecule already in the viewer"""


    def onAddCmdToViewer(self):
        if not hasattr(self.vf, 'atorsDict'):
            self.vf.atorsDict={}
        if not hasattr(self.vf,'readMolecule'):
            self.vf.loadCommand('fileCommands', 'readMolecule', 'Pmv',
                                topCommand=0)


    def __init__(self, mode='single', title = 'Choose Molecule for AutoDock4'):
        MVCommand.__init__(self)
        self.mode = mode
        self.title = title


    def chooseLigand_cb(self, event = None):
        """called each time the 'choose Ligand4' button is pressed"""
        mols = self.chooser.getMolSet()
        if mols is not None: 
            self.chooser.form.withdraw()
            self.doitWrapper(mols, redraw=0)
            #self.doitWrapper(mols, log=0, redraw=0)


    def guiCallback(self):
        self.chooser = MoleculeChooser(self.vf, self.mode, self.title)
        self.chooser.ipf.append({'name':'Select Button',
                                 'widgetType':Tkinter.Button,
                                 'text':'Select Molecule for AutoDock4',
                                 'wcfg':{'bd':6},
                                 'gridcfg':{'sticky':'we'},
                                 'command': self.chooseLigand_cb})
        self.form = self.chooser.go(modal=0, blocking=0)
        lb = self.chooser.ipf.entryByName['Molecule']['widget'].lb
        lb.bind("<Double-Button-1>",self.chooseLigand_cb)


    def __call__(self, nodes, **kw):
        """None<-ADtors4_chooseLigand(nodes)
nodes:ligand for autodock4
        """
        apply(self.doitWrapper, (nodes,), kw)


    def doit(self, nodes, **kw):
        mol = self.vf.expandNodes(nodes)[0]
        filename = mol.parser.filename
        ftype=split(filename,'.')[-1]
        #FIX THIS: could it be something else?
        #do we care???
        if ftype not in ['pdbq', 'pdb', 'pdbqs','pdbqt', 'mol2', 'dlg']:
            msg =  "unknown filetype: "+ ftype
            self.vf.warningMsg(msg)
            #????
            return
        useTorTree=0
        if hasattr(mol, 'torTree'):
            if 'useTorTree' in kw.keys():
                useTorTree = kw['useTorTree']
            else:
                msg = mol.name + ' already has a a torsion tree. Do you want to use it to set the activity of rotatable bonds?'
                d = SimpleDialog(self.vf.GUI.ROOT, text=msg, 
                    buttons=['No','Yes'], default=1, 
                    title='Use Previous Torsion Tree?')
                useTorTree=d.go()
                print 'set useTorTree to', useTorTree

        #include user_preferences here...
        #what about reprocessing?
        self.vf.atorsDict['molecule'] = mol
        if not hasattr(mol, 'LPO') or mol.LPO.version==3:
            cleanup = "nphs_lps"
            if self.vf.userpref['Automerge NPHS']['value']==0:
                cleanup = "lps"
            initLPO4(mol, cleanup=cleanup)
            title = "summary for " + mol.name
            self.vf.warningMsg(mol.LPO.summarize(), title=title)
            #warn if there are any atoms with zero charge besides carbons
            zero_charge_atoms = mol.allAtoms.get(lambda x: x.charge==0 and x.element!='C')
            if zero_charge_atoms is not None and len(zero_charge_atoms):
                msg = "These atoms have zero charge:\n"
                for a in zero_charge_atoms:
                    msg = msg + a.name + " "
                self.vf.warningMsg(msg)
            self.vf.allAtoms = self.vf.Mols.chains.residues.atoms
            #put in useTorTree stuff here:
            if useTorTree:
                self.rebuildTorTree(mol, mol.ROOT)
                #self.rebuildTorTree(mol, root)

        if self.vf.hasGui:
            self.vf.colorByAtomType(mol, ['lines'], topCommand=0)
            #nb aromCs could be an empty AtomSet
            self.vf.color(mol.LPO.aromCs,((0,1,0,),), ['lines'], topCommand=0)
            self.vf.centerScene(topCommand=0)
            self.vf.displayLines(mol,topCommand=0)
            self.vf.GUI.VIEWER.Redraw()
            #self.vf.GUI.ligandLabelLabel.config(text='Ligand:')
            #self.vf.GUI.ligandLabel.config(text=mol.name, width=len(mol.name))


    def rebuildTorTree(self, mol, root):
        #print '1:rebuilding torTree(flexibility pattern) from file'
        #have to rebuild it to capture which bonds are referenced
        #all bonds off to start
        allAts = mol.allAtoms
        for b in allAts.bonds[0]:
            if b.activeTors:
                b.activeTors = 0
        torscount = 0
        tM = mol.torTree.torsionMap
        for i in range(len(tM)):
            bnum0, bnum1 = tM[i].bond
            a0 = allAts.get(lambda x: x.number==bnum0 + 1)[0]
            #a0 = allAts[bnum0]
            a0.tt_ind = bnum0
            #a1 = allAts[bnum1]
            a1 = allAts.get(lambda x: x.number==bnum1 + 1)[0]
            a1.tt_ind = bnum1
            b = AtomSet([a0,a1]).bonds[0]
            if hasattr(b, 'possibleTors'):
                assert b.possibleTors
            else:
                b.possibleTors = 1
            b.activeTors = 1
            torscount = torscount + 1

        #this is also done in AtorsInitMol FIX THIS!!!
        mol.torscount = torscount
        mol.ROOT = root
        mol.ROOT.rnum0 = 0


    def onPick(self,event):
        listChooser = self.ipf.entryByName['Molecule']['widget']
        tkListBox = listChooser.lb
        atom,geom = self.vf.findPickedAtom(event)
        if atom is not None:
            pickedMol = atom.top
            #then need to make pickedMol the selection in self.lc
            for i in range(len(listChooser.entries)):
                listChooserlist=split(listChooser.entries[i][0])
                if pickedMol.name == listChooserlist[0]:
                    self.pickedMolIndex= i
                    tkListBox.select_clear(0,'end')
                    listChooser.select(i)
                    return
            t= "error: on  %s " %pickedMol.name
            self.vf.warningMsg(t)


Ators4MoleculeChooserGUI=CommandGUI()
Ators4MoleculeChooserGUI.addMenuCommand('AutoToolsBar', menuText['AutoTorsMB'], menuText['Choose Molecule4'], cascadeName = menuText['Input Molecule'])



class Ators4Reader(MVCommand):
    """allows user to select a file for the ligand via a file browser"""
 

    def onAddCmdToViewer(self):
        if not hasattr(self.vf, 'atorsDict'):
            self.vf.atorsDict={}
        if not hasattr(self.vf, 'readMolecule'):
            self.vf.loadCommand('fileCommands', 'readMolecule', 'Pmv')


    def guiCallback(self):
        """called each time the 'select molecule' button is pressed"""
        molFile = self.vf.askFileOpen(types=[ ('PDBQT files:', '*.pdbqt'),\
                ('PDBQ files:', '*.pdbq'),\
                ('MOL2 files:','*.mol2'), ('PDB files:', '*.pdb'),\
                ('all files:', '*')],\
                title = 'Ligand File for AutoDock4:')
        if not molFile: return
        self.doitWrapper(molFile, ask=1, redraw=1)


    def __call__(self, filename, log=1, **kw):
        """None<-ADtors4_readLigand(filename)
filename:file to read to get ligand for autodock
        """
        kw['log'] = log
        apply(self.doitWrapper, (filename,),kw)


    def doit(self, filename, **kw):
        #FIX THIS:
        #do we care???
        #could it be something else?
        ftype = split(filename,'.')[-1]
        if ftype not in ['pdbq','pdb','pdbqt', 'pdbqs','mol2']:
            msg = "unknown ligand file type " + ftype
            self.vf.warningMsg(msg)
            return
        mols = self.vf.readMolecule(filename, log=0)
        if not mols:
            return 'ERROR'
        if len(mols)>1:
            msg = str(len(mols)) + ' molecules in ', filename
            self.vf.warningMsg(msg)
            maxAts = 0
            for m in mols:
                numAts = len(m.allAtoms)
                if numAts>maxAts:
                    mol = m
                    maxAts = numAts
        else:
            mol = mols[0]
        if not mol.chains[0].hasBonds: 
            mol.buildBondsByDistance()

        cleanup = "nphs_lps"
        if self.vf.userpref['Automerge NPHS']['value']==0:
            cleanup = "lps"
        initLPO4(mol, cleanup=cleanup)
        title = "summary for " + mol.name
        self.vf.warningMsg(mol.LPO.summarize(), title=title)
        #ALWAYS update vf.allAtoms
        self.vf.allAtoms = self.vf.Mols.chains.residues.atoms
        self.vf.atorsDict['molecule'] = mol

        if self.vf.hasGui:
            self.vf.colorByAtomType(mol, ['lines'], topCommand=0)
            #nb aromCs could be an empty AtomSet
            self.vf.color(mol.LPO.aromCs,((0,1,0,),), ['lines'], topCommand=0)
            self.vf.centerScene(topCommand=0)
            self.vf.displayLines(mol,topCommand=0)
            self.vf.GUI.VIEWER.Redraw()
            #self.vf.GUI.ligandLabelLabel.config(text='Ligand:')
            #self.vf.GUI.ligandLabel.config(text=mol.name, width=len(mol.name))
        #what to do about this stuff?
        #if hasattr(self.vf.ADtors_defineRotBonds, 'noATBut'):
        #    menuText['MAmide'] = menuText['MAmide1']
        #    menuText['MGuan'] = menuText['MGuan1']
        #    menuText['MPeptide'] = menuText['MPeptide1']
        #    menuText['MActive'] = menuText['MActive1']
        #    menuText['MSelected'] = menuText['MSelected1']
        #    self.vf.ADtors_defineRotBonds.noATBut.config(text=menuText['MAmide'])
        #    self.vf.ADtors_defineRotBonds.noGTBut.config(text=menuText['MGuan'])
        #    self.vf.ADtors_defineRotBonds.noPBTBut.config(text=menuText['MPeptide'])
        #    self.vf.ADtors_defineRotBonds.noACTBut.config(text=menuText['MActive'])
        #    self.vf.ADtors_defineRotBonds.noSELBut.config(text=menuText['MSelected'])


Ators4ReaderGUI = CommandGUI()
Ators4ReaderGUI.addMenuCommand('AutoToolsBar', menuText['AutoTorsMB'], menuText['Read Molecule4'], cascadeName = menuText['Input Molecule'])



class AtorsRefWriter(MVCommand):
    """allows user to prepare a reference file for the ligand"""
 

    def onAddCmdToViewer(self):
        if not hasattr(self.vf, 'atorsDict'):
            self.vf.atorsDict={}
        if not hasattr(self.vf, 'readMolecule'):
            self.vf.loadCommand('fileCommands', 'readMolecule', 'Pmv')
        self.PdbqWriter = PdbqWriter()
        ##self.PdbqtWriter = PdbqtWriter()


    def guiCallback(self):
        """called each time the 'select reference' button is pressed"""
        #check that a ligand molecule has been written
        #ligfile = self.vf.atorsDict['outfile']
        ligfile = self.vf.atorsDict['molecule'].LPO.outputfilename
        if not ligfile:
            ligfile = self.vf.askFileOpen(types=[ ('formatted ligand files:', '*.pdbq'),\
                    ('all files:', '*')],\
                    title = 'Formatted Ligand File:')
        #if ligfile is not None:
        if ligfile:
            reffile = self.vf.askFileOpen(types=[ ('PDB files:', '*.pdb'),\
                    ('all files:', '*')],\
                    title = 'Autotors Reference File:')
            #if reffile is not None:
            if reffile:
                self.doitWrapper(ligfile, reffile, ask=1, redraw=1)


    def __call__(self, ligfile, reffile, log=1, **kw):
        """None<-ADtors_prepRef(ligfile, reffile)
ligfile: written output file
reffile:input file to reorder to be rms reference for autodock
        """
        kw['log'] = log
        apply(self.doitWrapper, (ligfile, reffile,),kw)


    def doit(self, ligfile, reffile, **kw):
        #FIX THIS:
        #do we care???
        #could it be something else?
        ligs = self.vf.readMolecule(ligfile)
        if not ligs:
            return 'ERROR'
        lig = ligs[0]
        if not lig.chains[0].hasBonds: 
            lig.buildBondsByDistance()
        refs = self.vf.readMolecule(reffile, log=0)
        if not refs:
            return 'ERROR'
        ref = refs[0]
        if not ref.chains[0].hasBonds: 
            ref.buildBondsByDistance()
        #check for same number of atoms and same atom names
        assert len(lig.allAtoms)==len(ref.allAtoms)
        refAts = ref.allAtoms
        refAtNames = refAts.name
        ref.allAtoms.written = 0
        lig.allAtoms.written = 0
        changedAts = []
        for ligAt in lig.allAtoms:
            #ligAtName = ligAt.name
            if ligAt.name[0]=='A':
                if len(ligAt.name)==1:
                    ligAt.name = 'C' 
                else:
                    ligAt.name = 'C' + ligAt.name[1:]
                changedAts.append(ligAt)
            #deal with hydrogens
            if ligAt.element!='H':
                assert ligAt.name in refAtNames
        outfilename = os.path.splitext(os.path.basename(reffile))[0] + '.ref.pdb'
        fptr = open(outfilename, 'w')
        ctr = 1
        for a in lig.allAtoms:
            #don't trust hydrogen names
            if a.element!='H' and a.name in refAtNames:
                #a2 is in ref with same name as a
                a2 = refAts.get(lambda x, n=a.name:x.name==n)[0]
            else:
                #try to match a neighbor atom's name (esp for hydrogens)
                foundAt = 0
                for b in a.bonds:
                    #neighbor is in ligand
                    neighbor = a.bonds[0].neighborAtom(a) 
                    #n2 is in ref
                    n2s = refAts.get(lambda x, n=neighbor.name: x.name==n)
                    if n2s is not None:
                        #n2 is in ref
                        for n2 in n2s:
                            for b in n2.bonds:
                                #a2 is in ref
                                a2 = b.neighborAtom(n2)
                                if a2.written: 
                                    continue
                                #take the first one of same element 
                                elif a2.element==a.element:
                                    #is this enough???
                                    foundAt = 1
                                    break
                    if not foundAt:
                        print 'could not match ', a.full_name(), ' ', a.number
            #mark this pair of atoms as written
            a.written = 1
            a2.written = 1
            a2.number = ctr 
            ctr = ctr + 1
            self.PdbqWriter.write_atom(fptr,a2)
            ##self.PdbqtWriter.write_atom(fptr,a2)
        fptr.close()
        for a in changedAts:
            if len(a.name)==1:
                a.name = 'A'
            else:
                a.name = 'A' + a.name[1:]
            


AtorsRefWriterGUI = CommandGUI()
AtorsRefWriterGUI.addMenuCommand('AutoToolsBar', menuText['AutoTorsMB'], menuText['Ref Molecule'], cascadeName = menuText['Input Molecule'])



class RigidMolecule(MVCommand):
    """allows user to write molecule with ROOT, ENDROOT + TORSDOF 0 added"""


    def guiCallback(self):
        molFile = self.vf.askFileOpen(types=[('select file:', '*.pdbq')],
                title = 'Read PDBQ for Rigid Autotors Output:')
        if not molFile: return
        outFile = self.vf.askFileSave(types=[('outputfile:', '*.pdbq')],
               title = 'Rigid Autotors Outputfile:')
        #if outFile is not None:
        if outFile: 
            self.doitWrapper(molFile, outFile, log=1, redraw=0)



    def __call__(self, molFile, outFile, **kw):
        """None<-ADtors_rigidLigand(molFile, outFile)
molFile:file to read to get ligand
outFile: file to write formatted rigid ligand
        """
        if not molFile:
            return 'ERROR'
        if not outFile:
            return 'ERROR'
        apply(self.doitWrapper, (molFile, outFile,),kw)



    def doit(self, molFile, outFile):
        molFileptr = open(molFile, 'r')
        allLines = molFileptr.readlines()
        molFileptr.close()
        outfptr = open(outFile, 'w')
        outstring = 'REMARK  0 active torsions:\n'
        outfptr.write(outstring)
        outstring = 'ROOT\n'
        outfptr.write(outstring)
        for l in allLines:
            startLine = l[:4]
            if startLine == 'ATOM' or startLine == 'HETA':
                outfptr.write(l)
        outstring = 'ENDROOT\n'
        outfptr.write(outstring)
        outstring = 'TORSDOF 0\n'
        outfptr.write(outstring)
        outfptr.close()


RigidMoleculeGUI = CommandGUI()
RigidMoleculeGUI.addMenuCommand('AutoToolsBar', menuText['AutoTorsMB'], menuText['Rigid Molecule'], cascadeName = menuText['Input Molecule'], separatorBelow=1)



class RigidMolecule4(MVCommand):
    """allows user to write molecule with ROOT, ENDROOT + TORSDOF 0 added"""


    def guiCallback(self):
        molFile = self.vf.askFileOpen(types=[('select file:', '*.pdbqt')],
                title = 'Read PDBQT for Rigid Autotors Output:')
        if not molFile: return
        outFile = self.vf.askFileSave(types=[('outputfile:', '*.pdbqt')],
               title = 'Rigid Autotors4 Outputfile:')
        #if outFile is not None:
        if outFile: 
            self.doitWrapper(molFile, outFile, log=1, redraw=0)



    def __call__(self, molFile, outFile, **kw):
        """None<-ADtors4_rigidLigand(molFile, outFile)
molFile:file to read to get ligand
outFile: file to write formatted rigid ligand
        """
        if not molFile:
            return 'ERROR'
        if not outFile:
            return 'ERROR'
        apply(self.doitWrapper, (molFile, outFile,),kw)



    def doit(self, molFile, outFile):
        molFileptr = open(molFile, 'r')
        allLines = molFileptr.readlines()
        molFileptr.close()
        outfptr = open(outFile, 'w')
        outstring = 'REMARK  0 active torsions:\n'
        outfptr.write(outstring)
        outstring = 'ROOT\n'
        outfptr.write(outstring)
        for l in allLines:
            startLine = l[:4]
            if startLine == 'ATOM' or startLine == 'HETA':
                outfptr.write(l)
        outstring = 'ENDROOT\n'
        outfptr.write(outstring)
        outstring = 'TORSDOF 0\n'
        outfptr.write(outstring)
        outfptr.close()


RigidMolecule4GUI = CommandGUI()
RigidMolecule4GUI.addMenuCommand('AutoToolsBar', menuText['AutoTorsMB'], menuText['Rigid Molecule4'], cascadeName = menuText['Input Molecule'], separatorBelow=1)


class AUTOTORSWriter(MVCommand):
    """allows user to select and write an output file for the formatted ligand
    for AutoDock3
"""


    def onAddCmdToViewer(self):
        if not hasattr(self.vf, 'atorsDict'):
            self.vf.atorsDict={}


    def doit(self, filename):
        #need to be sure filename is a string
        dict = self.vf.atorsDict
        if not dict.has_key('molecule'):
            self.vf.warningMsg(warningText['noAtorsMol'])
            return 'ERROR'
        mol = dict['molecule']
        #print "atcommands: calling write with ", filename
        #FIX THIS: if extension is pdbqt, write pdbqt
        if mol.LPO.version>=4:
            msg = mol.name + " currently formatted for AutoDock4!\nUnable to write AutoDock3 file"
            self.vf.warningMsg(msg)
            return "ERROR"
        mol.LPO.write(filename)
        rootSph_list = self.vf.GUI.VIEWER.findGeomsByName('rootSph')
        rootSph = rootSph_list[0]
        markSph_list = self.vf.GUI.VIEWER.findGeomsByName('markSph')
        markSph = rootSph_list[0]
        if self.vf.hasGui:
            rootSph.Set(visible=0)
            #self.vf.GUI.ligandLabelLabel.config(text='AD3 Ligand:')
            #self.vf.GUI.ligandLabel.config(text=filename, width=len(filename))


    def __call__(self, filename, **kw):
        """None<-ADtors_writeFormattedPDBQ(filename)
            filename: file to write formatted ligand
        """
        apply(self.doitWrapper, (filename,), kw)


    def guiCallback(self):
        hasGui = self.vf.hasGui
        dict = self.vf.atorsDict
        if not dict.has_key('molecule'):
            self.vf.warningMsg(warningText['noAtorsMol'])
            return
        mol = dict['molecule']
        if not hasattr(mol, 'ROOT'):
            self.vf.warningMsg('Must select root before writing file')
            return
        #newfile = self.vf.askFileSave(types=[('PDBQ files:', '*.out.pdbq',)],
        #    title = 'Formatted Autotors Molecule File:')
        ##if newfile is not None:
        #if newfile:
        #    self.doitWrapper(newfile, log=1, redraw=0)
        currentPath = os.getcwd()
        defaultFilename = os.path.join(currentPath, mol.name) + '.pdbq'
        newfile = self.vf.askFileSave(idir=currentPath, ifile=defaultFilename, 
            types=[('PDBQ files:', '*.pdbq',)],
            title = 'Formatted Autotors Molecule File:')
        if newfile:
            self.doitWrapper(newfile, log=1, redraw=0)


AUTOTORSWriterGUI=CommandGUI()
#AUTOTORSWriterGUI.addMenuCommand('AutoToolsBar', menuText['AutoTorsMB'], 'Write PDBQ ...')
AUTOTORSWriterGUI.addMenuCommand('AutoToolsBar', menuText['AutoTorsMB'], menuText['WritePDBQMB'],
            cascadeName = menuText['WriteMB'])



class AUTOTORS4Writer(MVCommand):
    """allows user to select and write an output file for the formatted ligand
for AutoDock4
ADD PDBQTWriter switch here-->>>
"""


    def onAddCmdToViewer(self):
        if not hasattr(self.vf, 'atorsDict'):
            self.vf.atorsDict={}


    def doit(self, filename):
        #need to be sure filename is a string
        dict = self.vf.atorsDict
        if not dict.has_key('molecule'):
            self.vf.warningMsg(warningText['noAtorsMol'])
            return 'ERROR'
        mol = dict['molecule']
        #print "atcommands: calling write with ", filename
        #FIX THIS: if extension is pdbqt, write pdbqt
        if mol.LPO.version==3:
            msg = mol.name + " currently formatted for AutoDock3!\nUnable to write AutoDock4 file"
            self.vf.warningMsg(msg)
            return "ERROR"
        mol.LPO.write(filename)
        rootSph_list = self.vf.GUI.VIEWER.findGeomsByName('rootSph')
        rootSph = rootSph_list[0]
        markSph_list = self.vf.GUI.VIEWER.findGeomsByName('markSph')
        markSph = rootSph_list[0]
        if self.vf.hasGui:
            rootSph.Set(visible=0)
            #self.vf.GUI.ligandLabelLabel.config(text='Ligand:')
            #self.vf.GUI.ligandLabel.config(text=filename, width=len(filename))


    def __call__(self, filename, **kw):
        """None<-ADtors4_writeFormattedPDBQT(filename)
            filename: file to write formatted ligand
        """
        apply(self.doitWrapper, (filename,), kw)


    def guiCallback(self):
        hasGui = self.vf.hasGui
        dict = self.vf.atorsDict
        if not dict.has_key('molecule'):
            self.vf.warningMsg(warningText['noAtorsMol'])
            return
        mol = dict['molecule']
        if not hasattr(mol, 'ROOT'):
            self.vf.warningMsg('Must select root before writing file')
            return
        currentPath = os.getcwd()
        defaultFilename = os.path.join(currentPath, mol.name) + '.pdbqt'
        newfile = self.vf.askFileSave(idir=currentPath, ifile=defaultFilename, 
            types=[('PDBQT files:', '*.pdbqt',)],
            title = 'Formatted Autotors Molecule File:')
        #if newfile is not None:
        if newfile:
            self.doitWrapper(newfile, log=1, redraw=0)


AUTOTORS4WriterGUI=CommandGUI()
AUTOTORS4WriterGUI.addMenuCommand('AutoToolsBar', menuText['AutoTorsMB'], menuText['WritePDBQTMB'],
            cascadeName = menuText['WriteMB'])



class MarkRoot(MVCommand):
    """shows current extent of root portion of the molecule:it includes all contiguous atoms starting with those adjacent to the designated root atom out to first active torsion """
    

    def onAddCmdToViewer(self):
        if not hasattr(self.vf, 'atorsDict'):
            self.vf.atorsDict={}
        self.markonoff = 0
        if self.vf.hasGui:
            self.markOn_Off = Tkinter.IntVar(master=self.vf.GUI.ROOT)
            self.markOn_Off.set(0)
            #initialize the geometries here:
            
            miscGeom = self.vf.GUI.miscGeom
            # we don't need this check anymore because miscGeom 
            #       is always added when we instantiate ViewerFrameworkGUI
            #if miscGeom not in self.vf.GUI.VIEWER.rootObject.children:
            #    self.vf.GUI.VIEWER.AddObject(miscGeom, redo=0)
            parentGeom_list = self.vf.GUI.VIEWER.findGeomsByName('autotors_geoms')
            if parentGeom_list==[]:
                parentGeom = Geom("autotors_geoms", shape=(0,0))
                #if parentGeom not in self.vf.GUI.VIEWER.rootObject.children:
                #parentGeom.replace = True
                self.vf.GUI.VIEWER.AddObject(parentGeom, redo=0, parent=miscGeom)
            else:
                parentGeom = parentGeom_list[0]
            markSph_list = self.vf.GUI.VIEWER.findGeomsByName('markSph')
            if markSph_list==[]:
                markSph = Spheres(name='autotors_markSph', materials=((0.,1.,0),),\
                        shape = (0,3), radii = 0.15, quality = 15,inheritMaterial=0,\
                        vertices=((0.,0.,0.),), visible=0, pickable=0)
                self.vf.GUI.VIEWER.AddObject(markSph, redo=0, parent=parentGeom)
            else:
                markSph = markSph_list[0]
            if not hasattr(self.vf, 'setICOM'):
                self.vf.loadCommand('interactiveCommands', 'setICOM', 'Pmv') 


    def __call__(self,**kw):
        """None<-ADtors_markRoot()
starts or stops marking root expansion display:
all atoms in root portion of ligand are marked with small sphere"""
        apply(self.doitWrapper,(), kw)


    def doit(self):
        if self.markonoff==0:
            return
        if not self.vf.atorsDict.has_key('molecule'):
            self.vf.warningMsg("you must select a ligand molecule first!")
            return "ERROR"
        mol = self.vf.atorsDict['molecule']
        rootSph_list = self.vf.GUI.VIEWER.findGeomsByName('rootSph')
        rootSph = rootSph_list[0]
        markSph_list = self.vf.GUI.VIEWER.findGeomsByName('markSph')
        markSph = rootSph_list[0]
        if not hasattr(mol, 'ROOT'):
            rootSph.Set(visible = 0)
            return
        allBonds = mol.allAtoms.bonds[0]
        for b in allBonds:
            b.marked = 0
        #allBonds.marked = 0
        self.neighborList = AtomSet()
        self.getNeighbors(mol.ROOT)
        #make neighbor geometry mark all those vertices
        #remove the b.marked attribute
        delattr(allBonds, 'marked')
        activeNeighbors = self.neighborList
        markSph_list = self.vf.GUI.VIEWER.findGeomsByName('markSph')
        markSph = markSph_list[0]
        if len(activeNeighbors):
            markSph.Set(visible=1)
            markSph.Set(vertices=activeNeighbors.coords)
        else:
            markSph.Set(visible=0)
        del self.neighborList


    def getNeighbors(self, at):
        numBonds = len(at.bonds)
        notActBonds = filter(lambda x: x.activeTors!=1, at.bonds)
        numMarkedBonds = filter(lambda x: x.marked==1, at.bonds)
        numNotActive = len(notActBonds)
        for b in notActBonds:
            if b.marked: 
                continue
            if numBonds==numMarkedBonds: 
                continue
            b.marked = 1
            if b.atom1!=at: 
                self.neighborList.append(b.atom1)
                self.getNeighbors(b.atom1)
            else:
                self.neighborList.append(b.atom2)
                self.getNeighbors(b.atom2)
                    

    def guiCallback(self):
        #to start or stop marking root expansion:
        markSph_list = self.vf.GUI.VIEWER.findGeomsByName('markSph')
        markSph = markSph_list[0]
        menuBarKey = self.vf.GUI.currentADTBar
        if self.markOn_Off.get():
            self.markonoff = 0
            markSph.Set(visible=0)
            
            menu = self.vf.GUI.menuBars[menuBarKey].menubuttons[menuText['AutoTorsMB']].menu
            #menu = self.vf.GUI.menuBars['AutoToolsBar'].menubuttons[menuText['AutoTorsMB']].menu
            children = menu.children[menuText['DefineRigidRootMB']]
            ind = children.index(menuText['ShowRootAtoms'])
            children.entryconfig(ind,{'label':menuText['SRA1']})
            menuText['ShowRootAtoms'] = menuText['SRA1']
            if self.vf.hasGui: 
                self.markOn_Off.set(0)
        else:
            self.markOn_Off.set(1)
            self.markonoff = 1
            if not self.vf.atorsDict.has_key('molecule'):
                self.vf.warningMsg(warningText['noAtorsMol'])
                return
            molecule = self.vf.atorsDict['molecule']
            if not hasattr(molecule, 'ROOT'):
                self.vf.warningMsg('select root FIRST!')
                return
            menu = self.vf.GUI.menuBars[menuBarKey].menubuttons[menuText['AutoTorsMB']].menu
            #menu = self.vf.GUI.currentADTBar.menubuttons[menuText['AutoTorsMB']].menu
            #menu = self.vf.GUI.menuBars['AutoToolsBar'].menubuttons[menuText['AutoTorsMB']].menu
            children = menu.children[menuText['DefineRigidRootMB']]
            ind = children.index(menuText['ShowRootAtoms'])
            children.entryconfig(ind,{'label':menuText['SRA2']})
            menuText['ShowRootAtoms'] = menuText['SRA2']
            self.doitWrapper(log=1,redraw=1)


MarkRootGUI=CommandGUI()
MarkRootGUI.addMenuCommand('AutoToolsBar', menuText['AutoTorsMB'], menuText['SRA1'],
            cascadeName = menuText['DefineRigidRootMB'])



class SelectRoot(MVCommand, MVAtomICOM):
    """allows user to pick an atom to be ROOT, the rigid portion 
of ligand which has rotatable BRANCHES """


    def onAddCmdToViewer(self):
        if not hasattr(self.vf, 'atorsDict'):
            self.vf.atorsDict={}
        if self.vf.hasGui: 
            #initialize the geometries here:
            parent = check_autotors_geoms(self.vf.GUI)
            #if rootSph not in self.vf.GUI.VIEWER.rootObject.children:
            rootSph_list = self.vf.GUI.VIEWER.findGeomsByName('rootSph')
            if rootSph_list==[]:
                rootSph = Spheres(name='autotors_rootSph', materials=((0.,1.,0),),\
                        shape = (0,3), radii = 0.3, quality = 15, inheritMaterial=0,\
                        vertices=((0.,0.,0.),), visible=0, pickable=0)
                rootSph.replace = True
                self.vf.GUI.VIEWER.AddObject(rootSph, 
                            redo=0, parent=parent)
            if not hasattr(self.vf, 'setICOM'):
                self.vf.loadCommand('interactiveCommands', 'setICOM', 'Pmv') 




    def __init__(self, func=None):
        MVCommand.__init__(self, func)
        MVAtomICOM.__init__(self)
        self.save = None


    def doit(self, atoms=None):
        if len(atoms)==0: return 
        atom = atoms[0]
        mol = atom.top
        if not hasattr(mol, 'LPO'):
            self.vf.warningMsg("picked atom not in a formatted molecule")
            return 'ERROR'
        if mol!=self.vf.atorsDict['molecule']:
            #FIX THIS: is this what should happen??
            self.vf.warningMsg("picked atom not in current ligand molecule\nSetting this molecule as current ligand molecule")
            self.vf.atorsDict['molecule'] = mol
        index = mol.allAtoms.index(atom)
        mol.LPO.setroot(index)
        rootSph_list = self.vf.GUI.VIEWER.findGeomsByName('rootSph')
        rootSph = rootSph_list[0]
        if self.vf.hasGui:
            rootSph.Set(vertices=(atom.coords,), visible=1)
            self.vf.GUI.VIEWER.Redraw()
            self.vf.ADtors_markRoot(topCommand=0,redraw=1)
            if self.save:
                self.vf.setICOM(self.save)
                self.save = None

            

    def guiCallback(self):
        if not self.vf.atorsDict.has_key("molecule"):
            self.vf.warningMsg(warningText['noAtorsMol'])
            return
        self.save = self.vf.ICmdCaller.commands.value[None]
        self.vf.setICOM(self, topCommand=0)
        self.vf.setIcomLevel( Atom )


    def __call__(self, atom, **kw):
        """None <- selectRoot(atom, **kw) 
set the root atom by setting mv.atorsDict['rootlist'] to [atom]"""
        if not self.vf.atorsDict.has_key("molecule"):
            self.vf.warningMsg(warningText['noAtorsMol'])
            return 'ERROR'
        if not atom:
            return 'ERROR'
        atoms = self.vf.expandNodes(atom)
        if not atoms:
            return 'ERROR'
        atoms = atoms.findType(Atom)
        if not atoms:
            return 'ERROR'
        apply( self.doitWrapper, (atoms,), kw)
        

SelectRootGUI=CommandGUI()
SelectRootGUI.addMenuCommand('AutoToolsBar', menuText['AutoTorsMB'],\
 menuText['ByPicking'], cascadeName = menuText['DefineRigidRootMB'])



class SetTorsionNumberGUICommand(MVCommand):
    """provides gui to ADtors_setTorsionNumber to specified number to inactivate and whether 
to inactive those which move the fewest atoms or those which move the most
    """


    def onRemoveObjectFromViewer(self, obj):
        dict = self.vf.atorsDict
        if dict.has_key('molecule') and obj==dict['molecule'] and hasattr(self, 'ifd'):
            #print 'deleting settorsion ifd'
            delattr(self, 'ifd')
        

    def onAddCmdToViewer(self):
        if not hasattr(self.vf, 'atorsDict'):
            self.vf.atorsDict={}
        if self.vf.hasGui:
            self.typeVar = Tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.typeVar.set('fewest')
            self.numTorsions = Tkinter.IntVar(master=self.vf.GUI.ROOT)


    def __call__(self, numTors, type='fewest',**kw):
        """None<-ADtors_setTorsionNumberGC(numTors, type)
        numTors number of activeTorsions at end
        type whether to inactive most or fewest movers
        """
        apply(self.doitWrapper,(numTors, type), kw)


    def doit(self, numTors, type):
        #print 'do something here'
        pass


    def guiCallback(self):
        #check that there's an ators molecule, root etc
        dict = self.vf.atorsDict
        if not dict.has_key('molecule'):
            self.vf.warningMsg("No ligand molecule selected")
            return 'ERROR'

        mol = dict['molecule']
        if not hasattr(mol, 'ROOT'):
            self.vf.warningMsg("Must define root before setting torsions")
            return 'ERROR'

        self.vf.ADtors_defineRotBonds.buildCol()
        #ok to build form etc
        if not hasattr(self, 'ifd'):
            self.buildForm()
            self.maxtors = mol.torscount
        else:
            self.form.deiconify()
        self.ctr.setentry(mol.torscount)
        posTors = len(mol.possible_tors_bnds)
        if self.maxtors<32:
            self.ctr._counterEntry._validationInfo['max'] = self.maxtors


    def slider_cb(self, event=None):
        #Pmw.Counter callback
        dict = self.vf.atorsDict
        mol = dict['molecule']
        posTors = mol.possible_tors
        val = int(self.ctr.get())
        self.numTorsions.set(val)
        #print "calling mol.LPO.limit_torsions with ", val, ' and ', self.typeVar.get()
        mol.LPO.limit_torsions(val, self.typeVar.get())
        #print "mol.activeTors=", len(filter(lambda x: x.activeTors==1, mol.allAtoms.bonds[0]))
        self.vf.ADtors_defineRotBonds.buildCol()


    def buildForm(self):
        ifd = self.ifd = InputFormDescr(title = 'Set Number of Active Torsions')
        ifd.append({'name': 'typeLab',
            'widgetType':Tkinter.Label,
            'text':'set number of active torsions moving:',
            'gridcfg':{'sticky':'we', 'columnspan':2}})
        ifd.append({'name':    'fewestRB',
            'widgetType':Tkinter.Radiobutton,
            'wcfg':{'text':'fewest atoms',
                    'variable':self.typeVar,
                    'value':'fewest'},
            'gridcfg':{'sticky':'w'}})
        ifd.append({'name':    'mostRB',
            'widgetType':Tkinter.Radiobutton,
            'wcfg':{'text':'most atoms',
                    'variable':self.typeVar,
                    'value':'most'}, 
            'gridcfg':{'sticky':'w','row':-1,'column':1}})
        ifd.append({'name':    'dividerLab',
            'widgetType':Tkinter.Label,
            'wcfg':{'text':'________________________________' }, 
            'gridcfg':{'sticky':'we','columnspan':2}})
        ifd.append({'widgetType':Pmw.Counter,
                'name':'numTorsCounter',
                'required':1,
                'wcfg':{'labelpos': 'n',
                    'label_text':'number of active torsions:  ',
                    'autorepeat':0,
                    'entryfield_value':0,
                    'entry_width':9,
                    'entryfield_validate':{'validator' : 'integer',
                               'min' : '0',
                               'max' : 32 },
                    'increment':1},
                'gridcfg':{'sticky':'nesw', 'columnspan':2}})
        ifd.append({'name':    'closeBut',
            'widgetType':Tkinter.Button,
            'wcfg':{'text':'Dismiss',
                    'command':self.Dismiss_cb}, 
            'gridcfg':{'sticky':'we','columnspan':4}})
        self.form = self.vf.getUserInput(self.ifd, modal=0, blocking=0)
        self.ctr = self.ifd.entryByName['numTorsCounter']['widget']
        da = self.ctr.component('downarrow')
        ua = self.ctr.component('uparrow')
        for item in [da, ua]:
            item.bind('<ButtonPress-1>', self.slider_cb, '+')
        entF = self.ctr.component('entryfield')._entryFieldEntry
        entF.bind('<Return>', self.slider_cb, '+')
        self.form.root.protocol('WM_DELETE_WINDOW',self.Dismiss_cb)
            

    def Dismiss_cb(self, event=None):
        mol = self.vf.atorsDict['molecule']
        self.vf.colorByAtomType(mol,
                    topCommand=0, redraw=1)
        #aromaticCs = mol.allAtoms.get(lambda x: x.autodock_element=='A')
        #FIX THIS
        aromaticCs = mol.LPO.aromCs
        if len(aromaticCs):
            self.vf.color(aromaticCs,((0.,1.,0.),),['lines'],topCommand=0)    
        self.form.withdraw()


SetTorsionNumberGUICommandGUI=CommandGUI()
SetTorsionNumberGUICommandGUI.addMenuCommand('AutoToolsBar', menuText['AutoTorsMB'],\
 menuText['SetTorsionNumber'], cascadeName = menuText['DefineRigidRootMB'] )



class SetTorsionNumber(MVCommand):
    """sets number of torsions to specified number by inactivating either those which
move the fewest atoms or those which move the most. if number is > than
current but less than possible, torsions are reactivated
    """


    def __call__(self, numTors, type='fewest',simpleModel=1, **kw):
        """None<-ADtors_setTorsionNumber(numTors, type='fewest', simpleModel=1)
simpleModel: 
    numTorsions torsions of type type are set active, rest inactive
non simpleModel:
    torsions in current ators molecule are activated or inactivated until total active 
is equal to numTors. Method decides whether to inactive those which move the fewest 
atoms or those which move the most.
        """
        dict = self.vf.atorsDict
        if not dict.has_key('molecule'):
            msg = 'no current autotors ligand molecule selected'
            return 'ERROR'
        mol = dict['molecule']
        assert hasattr(mol, 'LPO')
        #if not mol.processed_bonds:
            ##print 'calling processBonds from SetTorsionNumber'
            #self.vf.ADtors_processBonds(mol, topCommand=0)
        if not hasattr(mol, 'ROOT'):
            msg = 'must set root before limiting torsions'
            return 'ERROR'
        apply(self.doitWrapper,(numTors, type, simpleModel,), kw)


    def doit(self, numTors, type, simpleModel):
        dict = self.vf.atorsDict
        mol = dict['molecule']
        mol.LPO.limit_torsions(numTors, type)
        #if simpleModel:
        #    self.setTorsions(mol, numTors, type)
        #    return

#        #???update + keep other model???
#        if torscount==numTors:
#            msg = 'specified number==number present: no adjustment'
#            self.vf.warningMsg(msg)
#            if self.vf.hasGui:
#                self.vf.ADtors_defineRotBonds.buildCol()
#            return 'ERROR'
#        elif torscount<numTors:
#            if torscount==torsionMapNum:
#                msg = 'specified number > number possible: no adjustment'
#                self.vf.warningMsg(msg)
#                return 'ERROR'
#            else:
#                #in this case turn on as many as possible
#                if numTors>=torsionMapNum:
#                    #turn on everything
#                    delta = torsionMapNum - torscount
#                else:
#                    delta = numTors - torscount
#                self.turnOnTorsions(delta, type)
#        else:
#            #torscount>numTors
#            #in this case turn them off 
#            delta = torscount - numTors
#            self.turnOffTorsions(delta, type)


    def setTorsions(self, mol, numTors, type):
        dict = self.vf.atorsDict
        tNum = len(mol.torTree.torsionMap)
        if numTors>tNum:
            msg='too many torsions specified! '+ str(numTors)+  ' reducing to'+str(tNum)
            self.vf.warningMsg(msg)
            numTors = tNum
        if type=='fewest':
            rangeList = range(numTors)
        else:
            rangeList = []
            for k in range(1, numTors+1):
                rangeList.append(-k)
        #turn them all off
        torsionMap = mol.torTree.torsionMap
        for i in range(len(torsionMap)):
            node = torsionMap[i]
            b = mol.allAtoms.get(lambda x, node=node: x.tt_ind in node.bond).bonds[0][0]
            b.activeTors = 0
        #turn on the right number at correct end
        for i in rangeList:
            node = torsionMap[i]
            b = mol.allAtoms.get(lambda x, node=node: x.tt_ind in node.bond).bonds[0][0]
            b.activeTors = 1
        mol.torscount = numTors


    def turnOnTorsions(self, delta, type = 'fewest'):
        dict = self.vf.atorsDict
        mol = dict['molecule']
        allAts = mol.allAtoms
        torsionMap = mol.torTree.torsionMap
        torscount = mol.torscount
        #turn on delta torsions + adjust torscount in dict 
        if type=='fewest':
            rangeList = range(delta)
        else:
            rangeList = []
            for k in range(1, delta+1):
                rangeList.append(-k)
        for i in rangeList:
            node = torsionMap[i]
            b = allAts.get(lambda x, node=node: x.tt_ind in node.bond).bonds[0][0]
            if not b.activeTors:
                b.activeTors = 1
            else:
                lastInd = rangeList[-1]
                if type=='fewest':
                    rangeList.append(lastInd+1)
                else:
                    rangeList.append(lastInd-1)
        #torscount should be torscount + delta here
        mol.torscount = mol.torscount + numTors
        
        
    def turnOffTorsions(self, delta, type = 'fewest'):
        dict = self.vf.atorsDict
        mol = dict['molecule']
        allAts = mol.allAtoms
        torsionMap = mol.torTree.torsionMap
        torscount = mol.torscount
        #turn on delta torsions + adjust torscount in dict 
        if type=='fewest':
            rangeList = range(delta)
        else:
            rangeList = []
            for k in range(1, delta+1):
                rangeList.append(-k)
        for i in rangeList:
            node = torsionMap[i]
            if node.bond==(None,None):
                print 'error in turnOff torsions with ', rangeList
                break
            b = allAts.get(lambda x, node=node: x.tt_ind in node.bond).bonds[0][0]
            if b.activeTors:
                b.activeTors = 0
            else:
                lastInd = rangeList[-1]
                if type=='fewest':
                    rangeList.append(lastInd+1)
                else:
                    rangeList.append(lastInd-1)
        mol.torscount = mol.torscount - delta



class AutoRoot(MVCommand):
    """causes program to pick an atom to be ROOT: the one which has
the smallest 'largest sub-tree'"""


    def __call__(self, **kw):
        """None<-ADtors_autoRoot()
sets mv.atorsDict['rootlist'] to a list containing the atom with the smallest
largest subtree.
        """
        if not self.vf.atorsDict.has_key("molecule"):
            self.vf.warningMsg(warningText['noAtorsMol'])
            return 'ERROR'
        apply(self.doitWrapper,(), kw)



    def doit(self):
        """root atom is selected as follows:
        all atoms are evaluated for the number of atoms in subtrees.(The counting
        process is cut-off if the number counted is greater than the current smallest
        number counted).The atom with the smallest largest subtree is selected 
        for root.  Ties are resolved as follows: if only one is in a cycle, 
        it is selected, else arbitrarily the first found is selected.
        The selected atom is set to root as in setRoot"""
        #for each atom in atomlist: for all its bonds, get len(mol.subTree)
        #keep longest as member: maxbranch
        #then get the atom with the smallest maxbranch
        #but if getAutoRoot has been already called, keep previous center
        #if non-polar hydrogens have been merged, a correction in maxbranch is made..
        ###self.log()
        dict = self.vf.atorsDict
        if not dict.has_key('molecule'):
            self.vf.warningMsg(warningText['noAtorsMol'])
            return
        mol= dict['molecule']
        if len(mol.chains)>1:
            msg = "AutoRoot not implemented for molecules with >1 chain"
            self.vf.warningMsg(msg)
            return 'ERROR'
        mol.LPO.autoroot()
        if self.vf.hasGui:
            rootSph_list = self.vf.GUI.VIEWER.findGeomsByName('rootSph')
            rootSph = rootSph_list[0]
            rootSph.Set(vertices=(mol.autoRoot.coords,), visible=1)
            self.vf.ADtors_markRoot(topCommand=0,redraw=1)
            self.vf.GUI.message('autoRoot set to:' + mol.autoRoot.full_name())
            self.vf.GUI.VIEWER.Redraw()
        

    def guiCallback(self):
        self.doitWrapper(log=1,redraw=0)


AutoRootGUI=CommandGUI()
AutoRootGUI.addMenuCommand('AutoToolsBar', menuText['AutoTorsMB'],\
menuText['Automatically'], cascadeName = menuText['DefineRigidRootMB'])



############################################################################
############################################################################
#
#  Rotatable bonds definition
#
############################################################################
############################################################################



class SetRotatableBonds(MVCommand):

    form = None
    running = 0
    hasAmide = 1
    hasGuan = 1
    hasPeptide = 1
    hasSelected = 1
    torsStr = None
    success = 1


    def __init__(self, func=None):
        MVCommand.__init__(self, func)
        self.hasActive = 1
        

    def onAddCmdToViewer(self):
        if not hasattr(self.vf, 'atorsDict'):
            self.vf.atorsDict={}
        if self.vf.hasGui:
            if not hasattr(self.vf, 'setICOM'):
                self.vf.loadCommand('interactiveCommands', 'setICOM', 'Pmv') 
            if not self.torsStr:
                SetRotatableBonds.torsStr = Tkinter.StringVar(master=self.vf.GUI.ROOT)
            if not hasattr(self.vf,'labelByProperty'):
                self.vf.loadCommand('labelCommands', 'labelByProperty', 'Pmv',
                                    topCommand=0)


    def buildCol(self):
        mol = self.vf.atorsDict['molecule']
        #torscount = mol.torscount
        torscount = len(mol.allAtoms.bonds[0].get(lambda x: x.activeTors))
        #currentbonds=mol.geomContainer.atoms['lines'].bonds[0]
        currentbonds=mol.geomContainer.atoms['bonded'].bonds[0]
        col = []
        for b in currentbonds:
            if b.possibleTors:
                if b.activeTors: col.append((0,1,0))
                else: col.append((1,0,1))
            else:
                col.append((1,0,0))
        #mol.geomContainer.geoms['lines'].Set(materials = col,
        mol.geomContainer.geoms['bonded'].Set(materials = col,
                                              inheritMaterial=False,
                                              matBind =viewerConst.PER_PART)
        self.torsStr.set(menuText['torStr1'] + str(torscount) + menuText['torStr2'])
        self.vf.ADtors_markRoot(topCommand=0,redraw=1)


    def setNoActiveTors(self,log=0):
        mol = self.vf.atorsDict['molecule']
        if self.hasActive==1:
            self.hasActive = 0
            mol.LPO.set_all_torsions(0)  #force off here
        else:
            self.hasActive = 1
            mol.LPO.set_all_torsions(1)  #force on here


    def setNoActiveTors_cb(self, event=None, log=1, redraw=0):
        mol = self.vf.atorsDict['molecule']
        possiblebonds = filter(lambda x: x.possibleTors==1, mol.allAtoms.bonds[0])
        #if not len(pTatomset): return
        #activebonds = pTatomset.bonds[0]
        if not len(possiblebonds): return 
        if log:
            if self.vf.hasGui:
                msg='self.ADtors_defineRotBonds.setNoActiveTors_cb(log=0,redraw='+str(redraw)+')'
            else:
                msg='self.ADtors_defineRotBonds.setNoActiveTors(log=0)'
            self.vf.log(msg)
        if hasattr(self, 'noACTBut'):
            if not self.hasActive:
                menuText['MActive'] = menuText['MActive1']
                menuText['MGuan'] = menuText['MGuan1']
                menuText['MSelected'] = menuText['MSelected1']
                menuText['MAmide'] = menuText['MAmide1']
                menuText['MPeptide'] = menuText['MPeptide1']
            else:
                menuText['MActive'] = menuText['MActive2']
                menuText['MGuan'] = menuText['MGuan2']
                menuText['MSelected'] = menuText['MSelected2']
                menuText['MAmide'] = menuText['MAmide2']
                menuText['MPeptide'] = menuText['MPeptide2']
            self.noSELBut.config(text=menuText['MSelected'])
            self.noGTBut.config(text=menuText['MGuan'])
            self.noACTBut.config(text=menuText['MActive'])
            self.noATBut.config(text=menuText['MAmide'])
            self.noPBTBut.config(text=menuText['MPeptide'])
        self.setNoActiveTors()
        self.buildCol()


    def setNoGuanidiniumTors_cb(self, event= None, log=1, redraw=0):
        mol = self.vf.atorsDict['molecule']
        guanbonds = mol.guanidiniumbnds
        if not len(guanbonds): return
        if log:
            if self.vf.hasGui:
                msg='self.ADtors_defineRotBonds.setNoGuanidiniumTors_cb(log=0,redraw='+str(redraw)+')'
            else:
                msg='self.ADtors_defineRotBonds.setNoGuanidiniumTors(log=0)'
            self.vf.log(msg)
        if hasattr(self, 'noGTBut'):
            if not self.hasGuan:
                menuText['MGuan'] = menuText['MGuan1']
            else:
                menuText['MGuan'] = menuText['MGuan2']
            self.noGTBut.config(text=menuText['MGuan'])
        self.setNoGuanidiniumTors()
        self.buildCol()


    def setNoGuanidiniumTors(self, log=0):
        mol = self.vf.atorsDict['molecule']
        guanbonds = mol.guanidiniumbnds
        if not len(guanbonds): return 
        torscount = mol.torscount
        if self.hasGuan==1:
            self.hasGuan = 0
            mol.LPO.set_guanidinium_torsions(0) #force off 
        else:
            self.hasGuan = 1
            mol.LPO.set_guanidinium_torsions(1) #force off 
        if log:
            msg = 'self.ADtors_defineRotBonds.setNoGuanidiniumTors(log=0)'
            self.vf.log(msg)
            self.vf.message(msg)
        #if turned off: self.hasGuan==0, else: self.hasGuan==1
        return



    def setNoAmideTors(self, log=0):
        mol = self.vf.atorsDict['molecule']
        amidebonds = mol.amidebnds
        if not len(amidebonds):
            return
        if self.hasAmide==1:
            self.hasAmide = 0
            mol.LPO.set_amide_torsions(0) #force off 
        else:
            self.hasAmide = 1
            mol.LPO.set_amide_torsions(1) #force on 
        if log:
            msg = 'self.ADtors_defineRotBonds.setNoAmideTors(log=0)'
            self.vf.log(msg)
            self.vf.message(msg)
        #if turned off: self.hasAmide==0, else: self.hasAmide==1
        return


    def setNoSelected_cb(self, event= None, log=1,redraw=0):
        sel = self.vf.getSelection()
        if not sel:
            return
        selAts = sel.findType(Atom)
        if not selAts:
            return
        selbonds = selAts.bonds
        if not selbonds:
            return
        selbonds = selbonds[0]
        if not len(selbonds): 
            return
        if log:
            if self.vf.hasGui:
                msg='self.ADtors_defineRotBonds.setNoSelected_cb(log=0,redraw='+str(redraw)+')'
            else:
                msg='self.ADtors_defineRotBonds.setNoSelected(log=0)'
            self.vf.log(msg)
        if hasattr(self, 'noSELBut'):
            print 'self.hasSelected=',self.hasSelected
            if not self.hasSelected:
                menuText['MSelected'] = menuText['MSelected1']
            else:
                menuText['MSelected'] = menuText['MSelected2']
            self.noSELBut.config(text=menuText['MSelected'])
        self.setNoSelected()
        self.buildCol()


    def setNoSelected(self, log=0):
        mol = self.vf.atorsDict['molecule']
        #dict = self.vf.atorsDict
        selbonds = self.vf.getSelection().findType(Atom).bonds[0]
        if not len(selbonds): return
        torscount = mol.torscount
        if self.hasSelected==1:
            self.hasSelected=0
            for item in selbonds:
                ind1 = mol.allAtoms.index(item.atom1)
                ind2 = mol.allAtoms.index(item.atom2)
                mol.LPO.toggle_torsion(ind1, ind2) #force off here
        else:
            self.hasSelected = 1
            for item in selbonds:
                ind1 = mol.allAtoms.index(item.atom1)
                ind2 = mol.allAtoms.index(item.atom2)
                mol.LPO.toggle_torsion(ind1, ind2) #force on here
        if log:
            msg='self.ADtors_defineRotBonds.setNoSelected(log=0)'
            self.vf.log(msg)


    def setNoAmideTors_cb(self, event= None, log=1, redraw=0):
        amidebonds = self.vf.atorsDict['molecule'].amidebnds
        if not len(amidebonds): return
        if log:
            if self.vf.hasGui:
                msg='self.ADtors_defineRotBonds.setNoAmideTors_cb(log=0,redraw='+str(redraw)+')'
            else:
                msg='self.ADtors_defineRotBonds.setNoAmideTors(log=0)'
            self.vf.log(msg)
        if hasattr(self, 'noATBut'):
            if not self.hasAmide:
                menuText['MAmide'] = menuText['MAmide1']
            else:
                menuText['MAmide'] = menuText['MAmide2']
            self.noATBut.config(text=menuText['MAmide'])
        self.setNoAmideTors()
        self.buildCol()


    def setNoPeptideTors(self,log=0):
        mol = self.vf.atorsDict['molecule']
        pepbackbonds = mol.ppbbbnds
        if not len(pepbackbonds): 
            return
        torscount = mol.torscount
        if self.hasPeptide==1: #this is about the text on the button
            self.hasPeptide=0
            mol.LPO.set_ppbb_torsions(0) #force off 
        else:
            self.hasPeptide = 1
            mol.LPO.set_ppbb_torsions(1) #force on 
        if log:
            msg = 'self.ADtors_defineRotBonds.setNoPeptideTors(log=0)'
            self.vf.log(msg)
            self.vf.message(msg)


    def setNoPeptideTors_cb(self, event=None,log=1,redraw=0):
        mol = self.vf.atorsDict['molecule']
        pepbackbonds = mol.ppbbbnds
        if not len(pepbackbonds): return
        torscount = mol.torscount
        if log:
            if self.vf.hasGui:
                msg='self.ADtors_defineRotBonds.setNoPeptideTors_cb(log=0,redraw='+str(redraw)+')'
            else:
                msg='self.ADtors_defineRotBonds.setNoPeptideTors_cb(log=0)'
            self.vf.log(msg)
        if hasattr(self, 'noPBTBut'):
            if not self.hasPeptide:
                menuText['MPeptide'] = menuText['MPeptide1']
            else:
                menuText['MPeptide'] = menuText['MPeptide2']
            self.noPBTBut.config(text=menuText['MPeptide'])
        self.setNoPeptideTors()
        self.buildCol()



class DefiningRotatableBonds(SetRotatableBonds, MVBondICOM):
    """ """


    def __init__(self, func=None):
        SetRotatableBonds.__init__(self)
        MVBondICOM.__init__(self)
        self.save = None
        self.guiUp = 0
        #Overwrite default value inherited from ICOM
        self.pickLevel = 'parts'



    def onRemoveObjectFromViewer(self, obj):
        dict = self.vf.atorsDict
        #possibly window was open but 'molecule' removed from d before this
        #is called
        if self.form and hasattr(self.form,'root') and \
                self.form.root.winfo_ismapped() and not dict.has_key('molecule'):
            self.form.withdraw()
        


    def __call__(self, bonds,  **kw):
        """None<-ADtors_defineRotBonds(bonds)
bonds: possible torsions whose rotability is 
switched:1->0 and 0->1
            """
        kw['topCommand'] = 0
        kw['busyIdle'] = 1
        kw['log'] = 0
        self.setUpDisplay()
        apply(self.doitWrapper, (bonds,), kw)


    def doit(self, bonds):
        for bond in bonds:
            if not bond.possibleTors: continue
            mol = bond.atom1.top
            if not hasattr(mol, 'LPO'):
                print "bond in a non-formatted molecule"
                return "ERROR"
            ind1 = mol.allAtoms.index(bond.atom1)
            ind2 = mol.allAtoms.index(bond.atom2)
            mol.LPO.toggle_torsion(ind1, ind2)
            self.buildCol()
                                        

    def setUpDisplay(self):
        dict = self.vf.atorsDict
        if not dict.has_key('molecule'):
            self.vf.warningMsg(warningText['noAtorsMol'])
            return 
        mol = dict['molecule']
        if not hasattr(mol, 'LPO'):
            cleanup = "nphs_lps"
            if self.vf.userpref['Automerge NPHS']['value']==0:
                cleanup = "lps"
            initLPO(mol, cleanup=cleanup)
            title = "summary for " + mol.name
            self.vf.warningMsg(mol.LPO.summarize(),title=title)
            self.vf.allAtoms = self.vf.Mols.chains.residues.atoms
        if self.vf.hasGui and not self.guiUp:
            self.guiUp = 1
            pTbndset = mol.possible_tors_bnds
            pTatomset = (pTbndset.atom1 + pTbndset.atom2).uniq()
            if len(pTatomset):
                geom = mol.geomContainer.geoms['AtomLabels']
                geom.Set(billboard=True, fontStyle='solid', fontScales=(.3,.3,.3,))
                self.vf.labelByProperty(pTatomset, ('name', 'number'),
                        topCommand=0, redraw=1)
            self.displayForm()
            self.buildCol()


    def guiCallback(self):
        self.setUpDisplay()
        self.save = self.vf.ICmdCaller.commands.value[None]
        self.vf.setICOM(self, topCommand = 0)
        self.vf.setIcomLevel( Atom )


    def stop(self):
        self.done_cb()


    def getObjects(self, pick):
        atorsMolSet = self.vf.expandNodes(self.vf.atorsDict['molecule'])
        atorsMolGeom = atorsMolSet[0].geomContainer.geoms['bonded']
##         atorsMolGeom = atorsMolSet[0].geomContainer.geoms['lines']

        for o, val in pick.hits.items(): #loop over geometries
            primInd = map(lambda x: x[0], val)
            if o != atorsMolGeom: continue
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


    def dismiss(self):
        self.vf.setICOM(self.save, topCommand = 0)
        self.save = None
        self.done_cb()


    def done_cb(self):
        self.guiUp = 0
        dict = self.vf.atorsDict
        mol = dict['molecule']
        if not hasattr(mol, 'LPO'):
            print ' formatting ', mol.name, ' for autotors'
            cleanup = "nphs_lps"
            if self.vf.userpref['Automerge NPHS']['value']==0:
                cleanup = "lps"
            initLPO(mol, cleanup=cleanup)
            print "initLPO #4"
            self.vf.warningMsg(mol.LPO.summarize(),title=title)
            self.vf.allAtoms = self.vf.Mols.chains.residues.atoms
        if self.form:
            self.form.withdraw()
        pTbndset = mol.possible_tors_bnds
        pTatomset = (pTbndset.atom1 + pTbndset.atom2).uniq()
        if pTatomset is not None and len(pTatomset):
            self.vf.labelByProperty(pTatomset,('name','number',),
                  negate=1, topCommand=0, redraw=1)
        #FIX THIS: couldn't get here w/o molecule
        self.vf.colorByAtomType(mol, topCommand=0, redraw=1)
        aromaticCs = mol.LPO.aromCs
        #FIX THIS
        if len(aromaticCs):
            self.vf.color(aromaticCs,((0.,1.,0.),),['lines'],topCommand=0)    



    def displayForm(self):
        if self.form is not None:
            self.form.deiconify()
            return
        mol = self.vf.atorsDict['molecule']  
        self.hasAmide = mol.has_amide
        self.hasPeptide = mol.has_backbone
        self.hasGuanidinium = mol.has_guanidinium
        self.hasActive = 1
        self.hasSelected = 1
        torscount = len(mol.allAtoms.bonds[0].get(lambda x: x.activeTors))
        #self.hasAmide = 1
        #self.hasActive = 1
        #self.hasPeptide = 1
        #self.hasSelected = 1
        self.torsStr.set(menuText['torStr1'] + str(torscount) + menuText['torStr2'])
        ifd = self.ifd= InputFormDescr(title = 'Torsion Count')
        ifd.append({'name':'instrText',
                'widgetType':Tkinter.Label,
                'wcfg': {'text':'Pick or drag-&-pick bonds.  \nGreen = rotatable, \nMagenta = non-rotatable, \nRed = unrotatable.\n\n'},
                'gridcfg':{'sticky':'we'}}),
        ifd.append({'name':'torsEntryLab',
                'widgetType':Tkinter.Label,
                'wcfg':{'textvariable':self.torsStr},
                'gridcfg':{'sticky':'we'}}),
        ifd.append({'name':'noPepBut',
                'widgetType':Tkinter.Button,
              'wcfg':{'text':'Make peptide backbone bonds non-rotatable',
                  'command':self.setNoPeptideTors_cb},
                'gridcfg':{'sticky':'we',
                       'columnspan':2}}),
        ifd.append({'name':'noAmideBut',
                'widgetType':Tkinter.Button,
               'wcfg':{'text':'Make amide bonds non-rotatable',
                   'command':self.setNoAmideTors_cb},
                'gridcfg':{'sticky':'we', 'columnspan':2}}),
        ifd.append({'name':'noGuanBut',
                'widgetType':Tkinter.Button,
               'wcfg':{'text':'Make guanidinium bonds non-rotatable',
                   'command':self.setNoGuanidiniumTors_cb},
                'gridcfg':{'sticky':'we', 'columnspan':2}}),
        ifd.append({'name':'noSelectedBut',
                'widgetType':Tkinter.Button,
               'wcfg':{'text':'Make bonds between selected atoms non-rotatable',
                   'command':self.setNoSelected_cb},
                'gridcfg':{'sticky':'we', 'columnspan':2}}),
        ifd.append({'name':'noActiveBut',
                'widgetType':Tkinter.Button,
               'wcfg':{'text':'Make all active bonds non-rotatable',
                   'command':self.setNoActiveTors_cb},
                'gridcfg':{'sticky':'we', 'columnspan':2}}),
        ##ifd.append({'name':'forceTorsLab',
            ##'widgetType':Tkinter.Label,
            ##'wcfg':{'text':"(to force possible activity:\npress 'Shift' while picking a red bond\nNB:must not be in a ring)"},
            ##'gridcfg':{'sticky':'w' + 'e'}}),
        ifd.append({'name':'done',
                'widgetType':Tkinter.Button,
                'wcfg':{'text':'Done','command':self.dismiss},
                'gridcfg':{'sticky':'we',
                       'columnspan':2}})
        self.form = self.vf.getUserInput(self.ifd, modal=0, blocking=0)
        self.form.root.protocol('WM_DELETE_WINDOW',self.dismiss)
        self.noATBut = self.ifd.entryByName['noAmideBut']['widget']
        self.noGTBut = self.ifd.entryByName['noGuanBut']['widget']
        self.noPBTBut = self.ifd.entryByName['noPepBut']['widget']
        self.noACTBut = self.ifd.entryByName['noActiveBut']['widget']
        self.noSELBut = self.ifd.entryByName['noSelectedBut']['widget']
        #self.hasPeptide = mol.has_backbone
        if not self.hasPeptide:
            menuText['MPeptide'] = menuText['MPeptide2']
        else:
            menuText['MPeptide'] = menuText['MPeptide1']
        self.noPBTBut.config(text=menuText['MPeptide'])
        #self.hasGuanidinium = mol.has_guanidinium
        if not self.hasGuanidinium:
            menuText['MGuan'] = menuText['MGuan2']
        else:
            menuText['MGuan'] = menuText['MGuan1']
        self.noGTBut.config(text=menuText['MGuan'])
        #self.hasAmide = mol.has_amide
        if not self.hasAmide:
            menuText['MAmide'] = menuText['MAmide2']
        else:
            menuText['MAmide'] = menuText['MAmide1']
        self.noATBut.config(text=menuText['MAmide'])

#menuText['MActive1'] = 'Make all rotatable bonds non-rotatable'
#menuText['MActive2'] = 'Make all rotatable bonds rotatable'
        

DefiningRotatableBondsGUI = CommandGUI()
DefiningRotatableBondsGUI.addMenuCommand('AutoToolsBar', menuText['AutoTorsMB'],\
    menuText['DefineRotatableBonds'], cascadeName = menuText['DefineRigidRootMB'])



class SetBondRotatableFlag(SetRotatableBonds):
    """set the flag that tells whether a bond is rotatable in aan AutoDock
    ligand"""


    def setupUndoBefore(self, atoms, rotatable):
        self.addUndoCall( (atoms, not rotatable), {'redraw':1},
                  self.name )
        

    def __call__(self, atoms, rotatable, **kw):
        """None <- setBondRotatableFlag(atoms, rotatable, **kw)
rotatable can be either 1 or 0
        """
        atoms = self.vf.expandNodes(atoms)
        if len(atoms)<2: return
        assert isinstance( atoms[0], Atom )
        apply( self.doitWrapper, (atoms, rotatable), kw )

        
    def doit(self, atoms, rotatable):
        """atoms are two atoms in bond whose rotability is to be toggled"""
        dict = self.vf.atorsDict

        assert rotatable in [ 0, 1 ]
        if len(atoms) < 2: 
            return 'ERROR'
        
        bonds = atoms[:2].bonds
        if len(bonds[0])==0:
            print 'ERROR: no bond between ...'
            return 'ERROR'

        bond = bonds[0][0]

        mol = dict['molecule']
        if not hasattr(mol, 'LPO'):
            msg = 'ERROR: ', mol.name, ' not formatted! Choose as ligand first!'
            self.vf.warningMsg(msg)
            #print 'calling processBonds from setrotatable doit'
            #self.vf.ADtors_processBonds(mol, topCommand=0)
            return 'ERROR'

        if bond.possibleTors==0: 
            return 'ERROR'

        ind1 = mol.allAtoms.index(bond.atom1)
        ind2 = mol.allAtoms.index(bond.atom2)

        if bond.activeTors!=rotatable:
            mol.LPO.toggle_torsion(ind1, ind2)
        else: 
            return 'ERROR'
        if self.vf.hasGui:
            self.buildCol()



class CheckAromatic(MVCommand):
    """checks rings for planarity and changes carbons in planar rings to
element type and name: 'A'"""


    def onAddCmdToViewer(self):
        if not hasattr(self.vf, 'atorsDict'):
            self.vf.atorsDict={}
        d = self.vf.atorsDict
        self.vf.atorsDict['init_aCs'] = 0
        bondD = self.aromDict = {}
        bondD['TRP'] = ['NE1','CD1','CG','CD2','CE2','CZ2',\
                            'CH2','CZ3','CE3']
        bondD['TYR'] = ['CD1','CG','CD2','CE1','CE2','CZ']
        bondD['PHE'] = ['CD1','CG','CD2','CE1','CE2','CZ']
        bondD['HIS'] = ['CD2','CE1','CG','ND1','NE2']
        bondD['PRO'] = ['N','CD','CA','CB','CG']
        self.pep_aromList = ['PHE_CD1', 'PHE_CG', 'PHE_CD2', 'PHE_CE1',\
            'PHE_CE2', 'PHE_CZ', 'TYR_CD1', 'TYR_CG', 'TYR_CD2', 'TYR_CE1',\
            'TYR_CE2', 'TYR_CZ', 'HIS_CD2', 'HIS_CE1', 'HIS_CG', 'TRP_CD1',\
            'TRP_CG', 'TRP_CD2', 'TRP_CE2', 'TRP_CZ2', 'TRP_CH2', 'TRP_CZ3',\
            'TRP_CE3']
        self.useProteinAromaticList = 1
        doc = """For proteins, use standard look-up dictionary for aromatic bonds. Valid values are 1 or 0"""
        self.vf.userpref.add('Autotors: Protein Aromatic List', 1, [0,1],
                             callbackFunc = [self.set_useProteinAromaticList], 
                             doc=doc, category='AutoDockTools')


    def set_useProteinAromaticList(self, name, oldval, newval):
        self.useProteinAromaticList = newval


    def peptideA_init(self, mol):
        dkeys = self.aromDict.keys()
        resSet = mol.chains.residues.get(lambda x, dkeys=dkeys: x.type in dkeys)
        bondDict = {}
        if not resSet or len(resSet)==0:
            mol.cyclecount = 0
            return bondDict
        mol.cyclecount = numres = len(resSet)
        #NB: cyclenum is 1 based because 0 means not numbered
        for i in range(numres):
            res = resSet[i]
            ats = res.atoms
            #bondDict keys are 1 based too
            keys = self.aromDict[res.type]
            bnds = bondDict[i+1] = ats.get(lambda x, \
                                keys=keys:x.name in keys).bonds[0]
            for b in bnds:
                bnds.cyclenum = i + 1
        return bondDict


    def getPeptideAromatics(self, mol):
        atSet = AtomSet([])
        allAts = mol.allAtoms

        arom_ats = allAts.get(lambda x, l=self.pep_aromList:\
            x.parent.type+'_'+x.name in l)

        if not arom_ats or len(arom_ats)==0:
            return AtomSet([])
        for at in arom_ats:
            at.name = 'A' + at.name[1:]
            at.autodock_element = 'A'
        #print 'returning ', len(arom_ats), ' arom_ats'
        return arom_ats


    def doit(self):
        if not self.vf.atorsDict.has_key('molecule'):
            self.vf.warningMsg(warningText['noAtorsMol'])
            return
        dict = self.vf.atorsDict
        mol = dict['molecule']
        if not hasattr(mol, 'LPO'):
            msg = mol.name + ' not selected as Ligand molecule! ' 
            return 'ERROR'
        #re-detect there here
        aromaticCs = mol.LPO.setAromaticCarbons()
        if not len(aromaticCs):
            msg = mol.name + ' currently has no aromaticCs'
        if len(aromaticCs)>0 and self.vf.hasGui:
            self.vf.color(aromaticCs,((0.,1.,0.),),['lines'],topCommand=0)    
    

    def __call__(self, **kw):
        """None<- ADtors_changePlanarCarbonsToA
changes names of carbon atoms in planar rings from C* to A*"""
        apply(self.doitWrapper, (), kw)


CheckAromaticGUI = CommandGUI()
CheckAromaticGUI.addMenuCommand('AutoToolsBar', menuText['AutoTorsMB'],\
menuText['RenameAromaticCarbons'], cascadeName = menuText['AromaticCarbonsMB'])



class StopCheckAromatic(MVCommand):
    """restores carbons in aromatic rings to element type 'C' and name: 'C...'"""


    def doit(self):
        dict = self.vf.atorsDict
        mol = dict['molecule']
        aromaticCs = mol.LPO.aromCs
        if not len(aromaticCs):
            msg = mol.name + ' currently has no aromaticCs'
            return 'ERROR'
        if len(aromaticCs)>0:
            mol.LPO.set_carbon_names(aromaticCs, 'C')
            if self.vf.hasGui and len(aromaticCs):
                self.vf.color(aromaticCs,((1.,1.,1.),),['lines'],topCommand=0)    
    

    def __call__(self, **kw):
        """None<- ADtors_changeAromaticCarbonsToC
changes names of carbon atoms in planar rings from A* to C*"""
        apply(self.doitWrapper, (), kw)


StopCheckAromaticGUI = CommandGUI()
StopCheckAromaticGUI.addMenuCommand('AutoToolsBar', menuText['AutoTorsMB'],\
menuText['RestoreAliphaticCarbons'], cascadeName = menuText['AromaticCarbonsMB'])



class SetCarbonNames(MVCommand, MVAtomICOM):
    """toggles Carbon names between 'A' + 'C'"""


    def onAddCmdToViewer(self):
        if not hasattr(self.vf, 'atorsDict'):
            self.vf.atorsDict={}
        if self.vf.hasGui and not hasattr(self.vf, 'setICOM'):
            self.vf.loadCommand('interactiveCommands', 'setICOM', 'Pmv') 


    def __init__(self, func=None):
        MVCommand.__init__(self, func)
        MVAtomICOM.__init__(self)
        self.save = None


    def setupUndoBefore(self, atoms):
        oldname = atoms.name
        oldae = atoms.autodock_element
        for at in atoms:
            if at.name[0]=='C':
                at.name = 'A' + at.name[1:]
                at.autodock_element = 'A'
            elif at.name[0] == 'A':
                at.name = 'C' + at.name[1:]
                at.autodock_element = 'C'
        self.addUndoCall( (atoms,), {'redraw':1},
                  self.vf.ADtors_setCarbonNames.name)
        for i in range(len(atoms)):
            atoms[i].name = oldname[i]
            atoms[i].autodock_element = oldae[i]


    def __call__(self, nodes, **kw):
        """None<- ADtors_setCarbonNames
        
allows user to interactively toggle names of carbons from A* to C* and vice
versa"""
        if not self.vf.atorsDict.has_key('molecule'):
            self.vf.warningMsg(warningText['noAtorsMol'])
            return "ERROR"
        if not len(nodes):
            return "ERROR"
        atoms = self.vf.expandNodes(nodes)
        if not len(atoms):
            return "ERROR"
        aSet = atoms.findType(Atom)
        if not len(aSet): 
            return 'ERROR'
        atorsAtoms = None
        molecules, ats = self.vf.getNodesByMolecule(aSet)
        for mol, atomSets in map(None, molecules, ats):
            if mol == self.vf.atorsDict['molecule']:
                atorsAtoms = atomSets
                break
        if not atorsAtoms:
            return "ERROR"
        apply(self.doitWrapper,(atorsAtoms,), kw)


    def doit(self, atoms):
        #As turn into Cs
        #Cs turn into As
        dict = self.vf.atorsDict
        mol = dict['molecule']
        aromCs = mol.LPO.aromCs
        # WARNING: this forces all carbons in atoms to 'A' type
        As = AtomSet(filter(lambda x: (x.element=='C' and x.autodock_element=='A'), atoms))
        Cs = AtomSet(filter(lambda x: (x.element=='C' and x.autodock_element=='C'), atoms))
        mol.LPO.set_carbon_names(As, 'C')
        mol.LPO.set_carbon_names(Cs, 'A')
        if self.vf.hasGui:
            #changed Cs turn green and changed As turn white
            self.vf.color(Cs,((0.,1.,0.,),),['lines'],topCommand=0, redraw=1)
            self.vf.color(As,((1.,1.,1.,),),['lines'],topCommand=0, redraw=1)
            self.vf.GUI.VIEWER.Redraw()


    def guiCallback(self):
        if not self.vf.atorsDict.has_key('molecule'):
            self.vf.warningMsg(warningText['noAtorsMol'])
        self.save = self.vf.ICmdCaller.commands.value[None]
        self.vf.setICOM(self, topCommand=0)
        self.vf.setIcomLevel( Atom )
        ifd = self.ifd = InputFormDescr(title = 'STOP')
        ifd.append({'name':'stopBut',
            'widgetType':Tkinter.Button,
            'wcfg':{'text':'Stop Setting Carbon Names', 
                'command':self.dismiss_cb},
            'gridcfg':{'sticky':'we'}})
        self.form = self.vf.getUserInput(ifd, modal=0, blocking=0)
        self.form.root.protocol('WM_DELETE_WINDOW',self.dismiss_cb)


    def dismiss_cb(self, event = None):
        if hasattr(self, 'form'):
            self.form.root.withdraw()
        if self.save is not None:
            self.vf.setICOM(self.save)
            self.save = None


    def stop(self):
        self.dismiss_cb()


SetCarbonNamesGUI = CommandGUI()
SetCarbonNamesGUI.addMenuCommand('AutoToolsBar', menuText['AutoTorsMB'],\
menuText['SetCarbonNames'], cascadeName = menuText['AromaticCarbonsMB'])



class ChangeAromaticCutOff(MVCommand):
    """allows user to change requirement for aromaticity of cycles:
default is < 5 degrees between normals to adjacent atoms.  User enters a new angle in degrees"""


    def onAddCmdToViewer(self):
        if not hasattr(self.vf, 'atorsDict'):
            self.vf.atorsDict={}
        if self.vf.hasGui: 
            self.cutVal = Tkinter.StringVar(master=self.vf.GUI.ROOT)


    def __call__(self, val,  **kw):
        """None<-ADtors_changePlanarityCriteria
        val: new angle to use as criteria for planarity detection
if angle between normals to adjacent carbons in cycle is less than val, 
cycle is declared to be planar"""
        apply(self.doitWrapper, (val,), kw)


    def doit(self,val):
        if self.vf.hasGui:self.cutVal.set(str(val))
        dict = self.vf.atorsDict
        mol = dict['molecule']
        #LigandPreparationObject uses AromaticCarbonManager,ACM
        mol.LPO.changePlanarityCriteria(val)
        if self.vf.hasGui:
            Cs = AtomSet(filter(lambda x: x.autodock_element=='C',mol.allAtoms))
            As = AtomSet(filter(lambda x: x.autodock_element=='A',mol.allAtoms))
            if self.vf.hasGui:
                self.vf.color(As,((0.,1.,0.,),),['lines'],topCommand=0, redraw=1)
                self.vf.color(Cs,((1.,1.,1.,),),['lines'],topCommand=0, redraw=1)
                self.vf.GUI.VIEWER.Redraw()


    def guiCallback(self):
        mol = self.vf.atorsDict['molecule']
        val = mol.LPO.ACM.cutoff
        self.cutVal.set(str(val))
        ifd = self.ifd = InputFormDescr(title='CutOff Angle')
        ifd.append({'name':'aromVal',
            'widgetType':Tkinter.Entry,
            'wcfg':{
                'label': 'Enter angle in Degrees:',
                'textvariable':self.cutVal},
            'gridcfg':{'sticky':'w'}})
        val = self.vf.getUserInput(self.ifd)
        if val is not None:
            newval = float(val['aromVal'])
            self.doitWrapper(newval,log=1,redraw=0)
        

ChangeAromaticCutOffGUI=CommandGUI()
ChangeAromaticCutOffGUI.addMenuCommand('AutoToolsBar', menuText['AutoTorsMB'],\
menuText['ChangeAromaticityCriteria'], cascadeName = menuText['AromaticCarbonsMB'])



class TogglerootSphere(MVCommand):
    """lets user toggle rootSph visibility"""

    def onRemoveObjectFromViewer(self, obj):
        dict = self.vf.atorsDict
        rootSph_list = self.vf.GUI.VIEWER.findGeomsByName('rootSph')
        rootSph = rootSph_list[0]
        markSph_list = self.vf.GUI.VIEWER.findGeomsByName('markSph')
        markSph = rootSph_list[0]
        if dict.has_key('molecule') and obj==dict['molecule'] and rootSph.visible:
            #print 'hiding root sphere'
            #self.vf.ADtors_addChainToRootGC.chainSph.Set(visible=0)
            markSph.Set(visible=0)
            rootSph.Set(visible=0)
        

    def __call__(self, event=None, **kw):
        """None<- ADtors_showRootSphere
Allows the user to toggle the visibility of the sphere marking the atom
designated as the root atom in the ligand molecule"""
        kw['event'] = event
        apply(self.doitWrapper, (), kw)


    def doit(self, event=None):
        rootSph_list = self.vf.GUI.VIEWER.findGeomsByName('rootSph')
        rootSph = rootSph_list[0]
        if rootSph.visible:
            rootSph.Set(visible=0)
            #self.vf.ADtors_addChainToRootGC.chainSph.Set(visible=0)
            markSph.Set(visible=0)
        else:
            rootSph.Set(visible=1)
            #if len(self.vf.atorsDict['chain_rootlist']):
                #self.vf.ADtors_addChainToRootGC.chainSph.Set(visible=1)
            self.vf.ADtors_markRoot(topCommand=0,redraw=1)


    def guiCallback(self):
        self.doitWrapper(log=1,redraw=1)


TogglerootSphereGUI=CommandGUI()
TogglerootSphereGUI.addMenuCommand('AutoToolsBar', menuText['AutoTorsMB'],\
menuText['ShowAutotorsRootSphMB'], cascadeName = menuText['DefineRigidRootMB'], separatorBelow=1)



class AutoAutoTors(MVCommand):
    """performs default actions on designated molecule"""


    def onAddCmdToViewer(self):
        if not hasattr(self.vf, 'atorsDict'):
            self.vf.atorsDict={}
        self.ligand = None
        if self.vf.hasGui:
            self.inviewer = Tkinter.IntVar(master=self.vf.GUI.ROOT)
            self.outfileVar = Tkinter.StringVar(master=self.vf.GUI.ROOT)


    def onRemoveObjectFromViewer(self, obj):
        if obj==self.ligand:
            self.ligand = None


    def guiCallback(self):
        if not hasattr(self, 'ifd'):
            self.inviewer.set(2)
            self.buildForm()
        else:
            self.inviewer.set(2)
            self.outfileVar.set('')
            self.form.deiconify()


    def buildForm(self):
        ifd = self.ifd = InputFormDescr(title = 'AutoAutotors Parameters')
        ifd.append({'name': 'typeMolLab',
            'widgetType':Tkinter.Label,
            'text':'molecule location:',
            'gridcfg':{'sticky':'w'}})
        ifd.append({'name':    'pickMol',
            'widgetType':Tkinter.Radiobutton,
            'wcfg':{'text':'in viewer',
                    'variable':self.inviewer,
                    'value':0, 
                    'command':self.chooseMol},
            'gridcfg':{'sticky':'w','row':-1,'column':1}})
        ifd.append({'name':    'readMol',
            'widgetType':Tkinter.Radiobutton,
            'wcfg':{'text':'from file',
                    'variable':self.inviewer,
                    'value':1, 
                    'command':self.readMol},
            'gridcfg':{'sticky':'w','row':-1,'column':2}})
        ifd.append({'name':'outputFile',
            'widgetType':Tkinter.Entry,
            'wcfg':{
                'label': 'Output filename:',
                'width': 50,
                'textvariable':self.outfileVar},
            'gridcfg':{'sticky':'w','columnspan':4}})
        ifd.append({'name':    'goBut',
            'widgetType':Tkinter.Button,
            'wcfg':{'text':'Accept',
                    'command':self.go_cb}, 
            'gridcfg':{'sticky':'we', 'columnspan':2}})
        ifd.append({'name':    'closeBut',
            'widgetType':Tkinter.Button,
            'wcfg':{'text':'Cancel',
                    'command':self.cancel_cb}, 
            'gridcfg':{'sticky':'we','row':-1,'column':2}})
        self.form = self.vf.getUserInput(self.ifd, modal=0, blocking=0)
        self.form.root.protocol('WM_DELETE_WINDOW',self.cancel_cb)
            

    def cancel_cb(self, event=None):
        self.form.withdraw()


    def go_cb(self, event=None):
        if not self.ligand:
            if not self.inviewer.get():
                self.readMol()
            else:
                self.chooseMol()
        #print "self.outfileVar.get()=", self.outfileVar.get()
        self.doitWrapper(self.ligand,
                outfile=self.outfileVar.get(), 
                ask_outfile=0)
        

    def chooseMol(self, event=None):
        self.chooser = MoleculeChooser(self.vf, 'single', 'Choose Ligand')
        self.chooser.ipf.append({'name':'Select Button',
             'widgetType':Tkinter.Button,
             'text':'Select Molecule',
             'wcfg':{'bd':6},
             'gridcfg':{'sticky':'e'+'w'},
             'command': self.chooseMolecule_cb})
        self.form2 = self.chooser.go(modal=0, blocking=0)
        lb = self.chooser.ipf.entryByName['Molecule']['widget'].lb
        lb.bind("<Double-Button-1>", self.chooseMolecule_cb)


    def readMol(self, event=None):
        ligand = self.vf.askFileOpen(types=[('MOL2 files', '*.mol2'),
            ('PDBQ files', '*.pdbq'),('PDB files','*.pdb'),\
            ('all files', '*')], title = 'Ligand File:')
        #if ligand is not None:
        if ligand:
            #mols = self.vf.readMolecule(ligand)
            mols = Read(ligand)
            if not mols:
                return "ERROR"
            else:
                mol = mols[0]
            if not mol.chains[0].hasBonds: 
                mol.buildBondsByDistance()
            self.ligand = mol
            filename=os.path.split(ligand)[-1]
            ftype = split(ligand,'.')[-1]
            outfile = None
            if self.vf.hasGui and hasattr(self, 'outfileVar'):
                nameStr = mol.name + '.out.pdbqt'
                self.outfileVar.set(nameStr)
                outfile = nameStr
            #if ftype=='pdb':
                #self.processPDBLigand(mol, outfile)
            #FIX THIS:
            #could it be some other type???
            if ftype not in ['pdbq','mol2','pdb','pqr','pdbqt']:
                t= "file must be of type: pdb, pdb, mol2, pqr, pdbq, pdbqt"
                self.vf.warningMsg(t)
                return 'ERROR'
        else: self.ligand=None


    def chooseMolecule_cb(self, event = None):
        """called each time the molecule to be processed is already in viewer"""
        mols = self.chooser.getMolSet()
        if not mols: return
        if issubclass(mols.__class__, TreeNode):
            mols = mols.setClass([mols])
        mol = mols[0]
        if self.vf.hasGui and hasattr(self, 'outfileVar'):
            nameStr = mol.name + '.out.pdbqt'
            self.outfileVar.set(nameStr)
        self.ligand = mol
        self.chooser.form.withdraw()
        stem = mol.name
        #filename = os.path.split(mol.parser.filename)[1]
        filename = stem + ".out.pdbqt"
        self.outfileVar.set(filename)
        #ftype = split(filename,'.')[-1]
        #if ftype=='pdb': 
            #self.processPDBLigand(mol, outfile=self.outfileVar.get())

        #self.doitWrapper(self.ligand,
        #        CtoA=self.CtoA.get(),
        #        outfile=self.outfileVar.get(), 
        #        ask_outfile=1)

        #self.doitWrapper(self.ligand, outfile=self.outfileVar.get(), 
        #        ask_outfile=1)


    def __call__(self, mol, outfile=None, ask_outfile=1, **kw):
        """None<- ADtors_automaticLigandFormatting
outfile: name of file for formatted pdbqt output
ask_outfile: whether to ask the user about output filename
By default ask is equal to 1 which opens a dialog box to get 
user input for outfilename
if ask_outfile=0, outputfile name is name of molecule + '.out.pdbqt'
"""
        kw['outfile'] = outfile
        kw['ask_outfile'] = ask_outfile
        apply(self.doitWrapper, (mol,), kw)


    def doit(self, mol, **kw):
        outfile = kw['outfile']
        ask_outfile = kw['ask_outfile']
        #print 'in doit, ask_outfile=', ask_outfile
        if hasattr(self, 'form'):
            self.form.withdraw()
        dict = self.vf.atorsDict
        if type(mol)==types.StringType:
            #if there are any Mols, try to get this one
            if len(self.vf.Mols) and mol in self.vf.Mols.name:
                mol = self.vf.Mols.NodesFromName(mol)
            else:
                mol = self.vf.readMolecule(mol)
        if issubclass(mol.__class__, TreeNode):
            mol = mol.setClass([mol])
        assert isinstance(mol, TreeNodeSet)
        mol = mol[0]
        cleanup = "nphs_lps"
        if self.vf.userpref['Automerge NPHS']['value']==0:
            cleanup = "lps"
        initLPO4(mol, mode='automatic', root='auto', outputfilename=outfile,
                    cleanup=cleanup)
        title = "Automatic AutoDock4 Ligand Setup Summary for " + mol.name
        self.vf.warningMsg(mol.LPO.summarize(),title=title)
         

AutoAutoTorsGUI=CommandGUI()
AutoAutoTorsGUI.addMenuCommand('AutoToolsBar', menuText['AutoTorsMB'], \
menuText['AutomaticAutotorsSetupMB'], cascadeName = menuText['Input Molecule'])



class StopAutoTors(MVCommand):
    """hides AutoToolsBar and makes rootSph invisible"""

    def __call__(self, **kw):
        """None<-ADtors_stop
general cleanup of atorsDict and geometries used by autotorsCommands"""
        apply(self.doitWrapper,(), kw)

        
    def doit(self):
        #ASK MICHEL ABOUT THIS!!!!
        dict = self.vf.atorsDict
        if self.vf.hasGui:
            if hasattr(self.vf.GUI.VIEWER,'lastPickedObject') and hasattr(self.vf.GUI.VIEWER.lastPickedObject, 'mol'):
                if self.vf.GUI.VIEWER.lastPickedObject.mol not in self.vf.Mols:
                    self.vf.GUI.VIEWER.lastPickedObject.mol=None
                elif self.vf.GUI.VIEWER.lastPickedObject.mol == dict['molecule']:
                    self.vf.GUI.VIEWER.lastPickedObject.mol=None
            c=self.vf.ADtors_defineRotBonds
            if hasattr(c,'form') and c.form and c.form.root.winfo_viewable():
                c.quitsubPanel1()
            rootSph.Set(vertices=((0.,0.,0.),))
            rootSph.Set(visible=0)
            markSph.Set(vertices=((0.,0.,0.),))
            markSph.Set(visible=0)
        for k in ['amidebonds','aromaticCs','autoRoot',\
                  'pTatomset','pepbackbonds','rootlist','rootnum',\
                  'molecule','chain_rootlist']:
            val = dict.get(k, None)
            if val is not None:
                del dict[k]
        dict['aromaticCs'] = AtomSet([])


    def guiCallback(self):
        self.doitWrapper(log=1,redraw=0)



class AtorsInit(MVCommand):
    """ADtors_init is DEPRECATED. Remove from script"""


    def __call__(self, **kw):
        """ADtors_init is DEPRECATED. Remove from script"""
        apply(self.doitWrapper,(),kw)


    def doit(self):
        msg = "Deprecated command: ADtors_init\nRemove from script!!!"
        self.vf.warningMsg(msg)
   


class AtorsInitMol(MVCommand):
    """ADtors_initLigand is DEPRECATED. Remove from script"""


    def __call__(self, **kw):
        """ADtors_initLigand is DEPRECATED. Remove from script"""
        apply(self.doitWrapper,(),kw)


    def doit(self):
        msg = "Deprecated command: ADtors_initLigand\nRemove from script!!!"
        self.vf.warningMsg(msg)



class ProcessCharges(MVCommand):
    """ADtors_processCharges is DEPRECATED. Remove from script"""


    def __call__(self, **kw):
        """ADtors_processCharges is DEPRECATED. Remove from script"""
        apply(self.doitWrapper,(),kw)


    def doit(self):
        msg = "Deprecated command: ADtors_processCharges\nRemove from script!!!"
        self.vf.warningMsg(msg)
   


class ProcessBonds(MVCommand):
    """ADtors_processBonds is DEPRECATED. Remove from script"""


    def __call__(self, **kw):
        """ADtors_processBonds is DEPRECATED. Remove from script"""
        apply(self.doitWrapper,(),kw)


    def doit(self):
        msg = "Deprecated command: ADtors_processBonds\nRemove from script!!!"
        self.vf.warningMsg(msg)
   



commandList = [
    
    {'name':'ADtors4_readLigand','cmd':Ators4Reader(),
        'gui':Ators4ReaderGUI},
    {'name':'ADtors4_chooseLigand','cmd':Ators4MoleculeChooser(),
        'gui':Ators4MoleculeChooserGUI},
    {'name':'ADtors4_rigidLigand','cmd':RigidMolecule4(),
        'gui':RigidMolecule4GUI},
    {'name':'ADtors_readLigand','cmd':AtorsReader(),
        'gui':AtorsReaderGUI},
    {'name':'ADtors_chooseLigand','cmd':AtorsMoleculeChooser(),
        'gui':AtorsMoleculeChooserGUI},
    {'name':'ADtors_rigidLigand','cmd':RigidMolecule(),
        'gui':RigidMoleculeGUI},
    {'name':'ADtors_automaticLigandFormatting','cmd':AutoAutoTors(),
        'gui':AutoAutoTorsGUI},
#    {'name':'ADtors_writeRef','cmd':AtorsRefWriter(),
#        'gui':AtorsRefWriterGUI},
    {'name':'ADtors_setRoot','cmd':SelectRoot(),'gui':SelectRootGUI},
    {'name':'ADtors_autoRoot','cmd':AutoRoot(),'gui':AutoRootGUI},
#    {'name':'ADtors_addChainToRootGC','cmd':AddChainToRootGUICommand(),
#        'gui':AddChainToRootGUICommandGUI},
#    {'name':'ADtors_addChainToRoot','cmd':AddChainToRoot(),'gui':None},
#    {'name':'ADtors_removeChainFromRootGC','cmd':RemoveChainFromRootGUICommand(),
#        'gui':RemoveChainFromRootGUICommandGUI},
#    {'name':'ADtors_removeChainFromRoot','cmd':RemoveChainFromRoot(),'gui':None},
    {'name':'ADtors_markRoot','cmd':MarkRoot(),'gui':MarkRootGUI},
    {'name':'ADtors_showRootSphere','cmd':TogglerootSphere(),
        'gui':TogglerootSphereGUI},

    {'name':'ADtors_defineRotBonds', 'cmd':DefiningRotatableBonds(),
         'gui':DefiningRotatableBondsGUI },
    {'name':'ADtors_setBondRotatableFlag', 'cmd':SetBondRotatableFlag(),
         'gui':None },
    {'name':'ADtors_limitTorsionsGC','cmd':SetTorsionNumberGUICommand(),
        'gui':SetTorsionNumberGUICommandGUI},
    {'name':'ADtors_limitTorsions','cmd':SetTorsionNumber(),'gui':None},

#    {'name':'ADtors_changePlanarCarbonsToA','cmd':CheckAromatic(),
#        'gui':CheckAromaticGUI},
#    {'name':'ADtors_changeAromaticCarbonsToC','cmd':StopCheckAromatic(),
#        'gui':StopCheckAromaticGUI},
    {'name':'ADtors_setCarbonNames','cmd':SetCarbonNames(),
        'gui':SetCarbonNamesGUI},
    {'name':'ADtors_changePlanarityCriteria','cmd':ChangeAromaticCutOff(),
        'gui':ChangeAromaticCutOffGUI},
    {'name':'ADtors4_writeFormattedPDBQT','cmd':AUTOTORS4Writer(),
        'gui':AUTOTORS4WriterGUI},
    {'name':'ADtors_writeFormattedPDBQ','cmd':AUTOTORSWriter(),
        'gui':AUTOTORSWriterGUI},
    ]

def initModule(vf):
    vf.addCommand(AtorsInit(), 'ADtors_init')
    vf.addCommand( AtorsInitMol(), 'ADtors_initLigand')
    vf.addCommand( ProcessCharges(), 'ADtors_processCharges')
    vf.addCommand( ProcessBonds(), 'ADtors_processBonds')
    if not hasattr(vf, 'ADTSetMode'):
        vf.addCommand(AdtSetMode(), 'ADTSetMode')

    for dict in commandList:
        vf.addCommand(dict['cmd'],dict['name'],dict['gui'])

    vf.addCommand(StopAutoTors(), 'ADtors_stop')

    if vf.hasGui:
        vf.GUI.menuBars['AutoToolsBar']._frame.config( {'background':'tan'})
        for item in vf.GUI.menuBars['AutoToolsBar'].menubuttons.values():
            item.configure(background = 'tan')
        if not hasattr(vf.GUI, 'adtBar'):
            vf.GUI.adtBar = vf.GUI.menuBars['AutoToolsBar']
            vf.GUI.adtFrame = vf.GUI.adtBar.menubuttons.values()[0].master
#        if not hasattr(vf.GUI, 'ligandLabel'):
#            vf.GUI.ligandLabelLabel = Tkinter.Label(vf.GUI.adtFrame, \
#                            text="Ligand:", bg='tan')
#            vf.GUI.ligandLabelLabel.pack(side='left')
#            vf.GUI.ligandLabel=Tkinter.Label(vf.GUI.adtFrame, text="None", width=4,
#                                     relief='sunken', borderwidth=1,
#                                     anchor='w' )
#            vf.GUI.ligandLabel.pack(side='left')


