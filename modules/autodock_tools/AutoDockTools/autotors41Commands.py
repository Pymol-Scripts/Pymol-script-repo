#############################################################################
#
# Author: Ruth HUEY, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/autotors41Commands.py,v 1.2.4.1 2009/02/17 23:48:12 rhuey Exp $
#
# $Id: autotors41Commands.py,v 1.2.4.1 2009/02/17 23:48:12 rhuey Exp $ 
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


import Tkinter
from ViewerFramework.VFCommand import CommandGUI
from AutoDockTools.autotorsCommands import rootSph, markSph,\
menuText, AtorsMoleculeChooser, MAXTORS, AdtSetMode,\
AtorsReader, Ators4MoleculeChooser, Ators4Reader, AtorsRefWriter, \
RigidMolecule, RigidMolecule4, AUTOTORSWriter, AUTOTORS4Writer, \
MarkRoot, SelectRoot, SetTorsionNumberGUICommand, SetTorsionNumber, \
AutoRoot, SetRotatableBonds, DefiningRotatableBonds, SetBondRotatableFlag,\
CheckAromatic, StopCheckAromatic, SetCarbonNames, ChangeAromaticCutOff, \
TogglerootSphere, AutoAutoTors, StopAutoTors, AtorsInit, AtorsInitMol, \
ProcessCharges, ProcessBonds, rootSph, markSph, check_autotors_geoms,\
MAXTORS, menuText, warningText, checkMolCharges,\
autoMergeNPHS, set_autoMergeNPHS


Ators41MoleculeChooserGUI=CommandGUI()
Ators41MoleculeChooserGUI.addMenuCommand('AutoTools41Bar', menuText['AutoTorsMB'], menuText['Choose Molecule4'], cascadeName = menuText['Input Molecule'])



Ators41ReaderGUI = CommandGUI()
Ators41ReaderGUI.addMenuCommand('AutoTools41Bar', menuText['AutoTorsMB'], menuText['Read Molecule4'], cascadeName = menuText['Input Molecule'])


AtorsRefWriterGUI = CommandGUI()
AtorsRefWriterGUI.addMenuCommand('AutoTools41Bar', menuText['AutoTorsMB'], menuText['Ref Molecule'], cascadeName = menuText['Input Molecule'])


RigidMolecule41GUI = CommandGUI()
RigidMolecule41GUI.addMenuCommand('AutoTools41Bar', menuText['AutoTorsMB'], menuText['Rigid Molecule4'], cascadeName = menuText['Input Molecule'])


AUTOTORS41WriterGUI=CommandGUI()
AUTOTORS41WriterGUI.addMenuCommand('AutoTools41Bar', menuText['AutoTorsMB'], menuText['WritePDBQTMB'],
            cascadeName = menuText['WriteMB'])


MarkRootGUI=CommandGUI()
MarkRootGUI.addMenuCommand('AutoTools41Bar', menuText['AutoTorsMB'], menuText['SRA1'],
            cascadeName = menuText['DefineRigidRootMB'])


SelectRootGUI=CommandGUI()
SelectRootGUI.addMenuCommand('AutoTools41Bar', menuText['AutoTorsMB'],\
 menuText['ByPicking'], cascadeName = menuText['DefineRigidRootMB'])


SetTorsionNumberGUICommandGUI=CommandGUI()
SetTorsionNumberGUICommandGUI.addMenuCommand('AutoTools41Bar', menuText['AutoTorsMB'],\
 menuText['SetTorsionNumber'], cascadeName = menuText['DefineRigidRootMB'] )



AutoRootGUI=CommandGUI()
AutoRootGUI.addMenuCommand('AutoTools41Bar', menuText['AutoTorsMB'],\
menuText['Automatically'], cascadeName = menuText['DefineRigidRootMB'])


DefiningRotatableBondsGUI = CommandGUI()
DefiningRotatableBondsGUI.addMenuCommand('AutoTools41Bar', menuText['AutoTorsMB'],\
    menuText['DefineRotatableBonds'], cascadeName = menuText['DefineRigidRootMB'])



CheckAromaticGUI = CommandGUI()
CheckAromaticGUI.addMenuCommand('AutoTools41Bar', menuText['AutoTorsMB'],\
menuText['RenameAromaticCarbons'], cascadeName = menuText['AromaticCarbonsMB'])



StopCheckAromaticGUI = CommandGUI()
StopCheckAromaticGUI.addMenuCommand('AutoTools41Bar', menuText['AutoTorsMB'],\
menuText['RestoreAliphaticCarbons'], cascadeName = menuText['AromaticCarbonsMB'])



SetCarbonNamesGUI = CommandGUI()
SetCarbonNamesGUI.addMenuCommand('AutoTools41Bar', menuText['AutoTorsMB'],\
menuText['SetCarbonNames'], cascadeName = menuText['AromaticCarbonsMB'])


ChangeAromaticCutOffGUI=CommandGUI()
ChangeAromaticCutOffGUI.addMenuCommand('AutoTools41Bar', menuText['AutoTorsMB'],\
menuText['ChangeAromaticityCriteria'], cascadeName = menuText['AromaticCarbonsMB'])



TogglerootSphereGUI=CommandGUI()
TogglerootSphereGUI.addMenuCommand('AutoTools41Bar', menuText['AutoTorsMB'],\
menuText['ShowAutotorsRootSphMB'], cascadeName = menuText['DefineRigidRootMB'])


AutoAutoTorsGUI=CommandGUI()
AutoAutoTorsGUI.addMenuCommand('AutoTools41Bar', menuText['AutoTorsMB'], \
menuText['AutomaticAutotorsSetupMB'], cascadeName = menuText['Input Molecule'])





commandList = [
    
    {'name':'AD41tors_readLigand','cmd':Ators4Reader(),
        'gui':Ators41ReaderGUI},
    {'name':'AD41tors_chooseLigand','cmd':Ators4MoleculeChooser(),
        'gui':Ators41MoleculeChooserGUI},
    {'name':'AD41tors_rigidLigand','cmd':RigidMolecule4(),
        'gui':RigidMolecule41GUI},
    {'name':'AD41tors_automaticLigandFormatting','cmd':AutoAutoTors(),
        'gui':AutoAutoTorsGUI},
#    {'name':'AD41tors_writeRef','cmd':AtorsRefWriter(),
#        'gui':AtorsRefWriterGUI},
    {'name':'AD41tors_setRoot','cmd':SelectRoot(),'gui':SelectRootGUI},
    {'name':'AD41tors_autoRoot','cmd':AutoRoot(),'gui':AutoRootGUI},
#    {'name':'AD41tors_addChainToRootGC','cmd':AddChainToRootGUICommand(),
#        'gui':AddChainToRootGUICommandGUI},
#    {'name':'AD41tors_addChainToRoot','cmd':AddChainToRoot(),'gui':None},
#    {'name':'AD41tors_removeChainFromRootGC','cmd':RemoveChainFromRootGUICommand(),
#        'gui':RemoveChainFromRootGUICommandGUI},
#    {'name':'AD41tors_removeChainFromRoot','cmd':RemoveChainFromRoot(),'gui':None},
    {'name':'AD41tors_markRoot','cmd':MarkRoot(),'gui':MarkRootGUI},
    {'name':'AD41tors_showRootSphere','cmd':TogglerootSphere(),
        'gui':TogglerootSphereGUI},
    {'name':'AD41tors_defineRotBonds', 'cmd':DefiningRotatableBonds(),
         'gui':DefiningRotatableBondsGUI },
    {'name':'AD41tors_limitTorsionsGC','cmd':SetTorsionNumberGUICommand(),
        'gui':SetTorsionNumberGUICommandGUI},
##    {'name':'AD41tors_changePlanarCarbonsToA','cmd':CheckAromatic(),
##        'gui':CheckAromaticGUI},
##    {'name':'AD41tors_changeAromaticCarbonsToC','cmd':StopCheckAromatic(),
##        'gui':StopCheckAromaticGUI},
    {'name':'AD41tors_setCarbonNames','cmd':SetCarbonNames(),
        'gui':SetCarbonNamesGUI},
    {'name':'AD41tors_changePlanarityCriteria','cmd':ChangeAromaticCutOff(),
        'gui':ChangeAromaticCutOffGUI},
    {'name':'AD41tors_writeFormattedPDBQT','cmd':AUTOTORS4Writer(),
        'gui':AUTOTORS41WriterGUI},
    ]

def initModule(vf):

    for dict in commandList:
        vf.addCommand(dict['cmd'],dict['name'],dict['gui'])
    if not hasattr(vf, 'ADtors_limitTorsions'):
        vf.addCommand(SetTorsionNumber(), 'ADtors_limitTorsions')
    if not hasattr(vf, 'ADtors_setBondRotatableFlag'):
        vf.addCommand(SetBondRotatableFlag(), 'ADtors_setBondRotatableFlag')
    if not hasattr(vf, 'ADtors_stop'):
        vf.addCommand(StopAutoTors(), 'ADtors_stop')



    if not hasattr(vf, 'ADTSetMode'):
        vf.addCommand(AdtSetMode(), 'ADTSetMode')
    if vf.hasGui:
        vf.GUI.menuBars['AutoTools41Bar']._frame.config( {'background':'tan'})
        for item in vf.GUI.menuBars['AutoTools41Bar'].menubuttons.values():
            item.configure(background = 'tan')
        if not hasattr(vf.GUI, 'adt41Bar'):
            vf.GUI.adt41Bar = vf.GUI.menuBars['AutoTools41Bar']
            vf.GUI.adt41Frame = vf.GUI.adt41Bar.menubuttons.values()[0].master
        if not hasattr(vf.GUI, 'adt41ModeLabel'):
            mbs = {}
            packing_list = []
            for c in vf.GUI.adt41Frame.children.values():
                if isinstance(c, Tkinter.Menubutton):
                    mbs[c.cget('text')] = c
                    packing_list.append(c.cget('text'))
                    c.pack_forget()
            vf.GUI.adt41ModeLabel=Tkinter.Label(vf.GUI.adt41Frame, text="ADT4.1", width=6,
                             relief='sunken', borderwidth=1, fg='DarkGreen',
                             bg = 'ivory',anchor='w' )
            vf.GUI.adt41ModeLabel.pack(side='left')
            vf.GUI.adt41ModeLabel.bind("<Double-Button-1>", vf.ADTSetMode.guiCallback())
            for t in packing_list:
                try:
                    c = mbs[t]
                    c.pack(side='left')
                except:
                    pass

#        if not hasattr(vf.GUI, 'ligandLabel'):
#            vf.GUI.ligandLabelLabel = Tkinter.Label(vf.GUI.adtFrame, \
#                            text="Ligand:", bg='tan')
#            vf.GUI.ligandLabelLabel.pack(side='left')
#            vf.GUI.ligandLabel=Tkinter.Label(vf.GUI.adtFrame, text="None", width=41,
#                                     relief='sunken', borderwidth=1,
#                                     anchor='w' )
#            vf.GUI.ligandLabel.pack(side='left')


