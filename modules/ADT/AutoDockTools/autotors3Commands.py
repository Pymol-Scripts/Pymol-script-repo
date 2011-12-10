#############################################################################
#
# Author: Ruth HUEY, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/autotors3Commands.py,v 1.3 2009/02/10 17:44:05 rhuey Exp $
#
# $Id: autotors3Commands.py,v 1.3 2009/02/10 17:44:05 rhuey Exp $
#
#
#
#
#
#
#
"""
This Module is used to setup files for autodock3 calculations.
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



#???Set RotatableBonds???
from AutoDockTools.autotorsCommands import AtorsMoleculeChooser


AtorsMoleculeChooserGUI=CommandGUI()
AtorsMoleculeChooserGUI.addMenuCommand('AutoTools3Bar', menuText['AutoTorsMB'],\
menuText['Choose Molecule'], cascadeName = menuText['Input Molecule'])


AtorsReaderGUI = CommandGUI()
AtorsReaderGUI.addMenuCommand('AutoTools3Bar', menuText['AutoTorsMB'],\
menuText['Read Molecule'], cascadeName = menuText['Input Molecule'])


AtorsRefWriterGUI = CommandGUI()
AtorsRefWriterGUI.addMenuCommand('AutoTools3Bar', menuText['AutoTorsMB'],\
menuText['Ref Molecule'], cascadeName = menuText['Input Molecule'])


RigidMoleculeGUI = CommandGUI()
RigidMoleculeGUI.addMenuCommand('AutoTools3Bar', menuText['AutoTorsMB'],\
menuText['Rigid Molecule'], cascadeName = menuText['Input Molecule'])


AUTOTORSWriterGUI=CommandGUI()
AUTOTORSWriterGUI.addMenuCommand('AutoTools3Bar', menuText['AutoTorsMB'],\
menuText['WritePDBQMB'], cascadeName = menuText['WriteMB'])


MarkRootGUI=CommandGUI()
MarkRootGUI.addMenuCommand('AutoTools3Bar', menuText['AutoTorsMB'],\
menuText['SRA1'], cascadeName = menuText['DefineRigidRootMB'])


SelectRootGUI=CommandGUI()
SelectRootGUI.addMenuCommand('AutoTools3Bar', menuText['AutoTorsMB'],\
menuText['ByPicking'], cascadeName = menuText['DefineRigidRootMB'])


SetTorsionNumberGUICommandGUI=CommandGUI()
SetTorsionNumberGUICommandGUI.addMenuCommand('AutoTools3Bar', menuText['AutoTorsMB'],\
 menuText['SetTorsionNumber'], cascadeName = menuText['DefineRigidRootMB'] )


AutoRootGUI=CommandGUI()
AutoRootGUI.addMenuCommand('AutoTools3Bar', menuText['AutoTorsMB'],\
menuText['Automatically'], cascadeName = menuText['DefineRigidRootMB'])


DefiningRotatableBondsGUI = CommandGUI()
DefiningRotatableBondsGUI.addMenuCommand('AutoTools3Bar', menuText['AutoTorsMB'],\
menuText['DefineRotatableBonds'], cascadeName = menuText['DefineRigidRootMB'])


CheckAromaticGUI = CommandGUI()
CheckAromaticGUI.addMenuCommand('AutoTools3Bar', menuText['AutoTorsMB'],\
menuText['RenameAromaticCarbons'], cascadeName = menuText['AromaticCarbonsMB'])



StopCheckAromaticGUI = CommandGUI()
StopCheckAromaticGUI.addMenuCommand('AutoTools3Bar', menuText['AutoTorsMB'],\
menuText['RestoreAliphaticCarbons'], cascadeName = menuText['AromaticCarbonsMB'])



SetCarbonNamesGUI = CommandGUI()
SetCarbonNamesGUI.addMenuCommand('AutoTools3Bar', menuText['AutoTorsMB'],\
menuText['SetCarbonNames'], cascadeName = menuText['AromaticCarbonsMB'])



ChangeAromaticCutOffGUI=CommandGUI()
ChangeAromaticCutOffGUI.addMenuCommand('AutoTools3Bar', menuText['AutoTorsMB'],\
menuText['ChangeAromaticityCriteria'], cascadeName = menuText['AromaticCarbonsMB'])



TogglerootSphereGUI=CommandGUI()
TogglerootSphereGUI.addMenuCommand('AutoTools3Bar', menuText['AutoTorsMB'],\
menuText['ShowAutotorsRootSphMB'], cascadeName = menuText['DefineRigidRootMB'])


AutoAutoTors3GUI=CommandGUI()
AutoAutoTors3GUI.addMenuCommand('AutoTools3Bar', menuText['AutoTorsMB'], \
menuText['AutomaticAutotorsSetupMB'], cascadeName = menuText['Input Molecule'])


commandList = [
    
    {'name':'AD3tors_readLigand','cmd':AtorsReader(),
        'gui':AtorsReaderGUI},
    {'name':'AD3tors_chooseLigand','cmd':AtorsMoleculeChooser(),
        'gui':AtorsMoleculeChooserGUI},
    {'name':'AD3tors_rigidLigand','cmd':RigidMolecule(),
        'gui':RigidMoleculeGUI},
    {'name':'AD3tors_automaticLigandFormatting','cmd':AutoAutoTors(),
        'gui':AutoAutoTors3GUI},
#    {'name':'AD3tors_writeRef','cmd':AtorsRefWriter(),
#        'gui':AtorsRefWriterGUI},
    {'name':'AD3tors_setRoot','cmd':SelectRoot(),'gui':SelectRootGUI},
    {'name':'AD3tors_autoRoot','cmd':AutoRoot(),'gui':AutoRootGUI},
#    {'name':'AD3tors_addChainToRootGC','cmd':AddChainToRootGUICommand(),
#        'gui':AddChainToRootGUICommandGUI},
#    {'name':'AD3tors_addChainToRoot','cmd':AddChainToRoot(),'gui':None},
#    {'name':'AD3tors_removeChainFromRootGC','cmd':RemoveChainFromRootGUICommand(),
#        'gui':RemoveChainFromRootGUICommandGUI},
#    {'name':'AD3tors_removeChainFromRoot','cmd':RemoveChainFromRoot(),'gui':None},
    {'name':'AD3tors_markRoot','cmd':MarkRoot(),'gui':MarkRootGUI},
    {'name':'AD3tors_showRootSphere','cmd':TogglerootSphere(),
        'gui':TogglerootSphereGUI},

    {'name':'AD3tors_defineRotBonds', 'cmd':DefiningRotatableBonds(),
         'gui':DefiningRotatableBondsGUI },
    {'name':'AD3tors_limitTorsionsGC','cmd':SetTorsionNumberGUICommand(),
        'gui':SetTorsionNumberGUICommandGUI},
##    {'name':'AD3tors_changePlanarCarbonsToA','cmd':CheckAromatic(),
##        'gui':CheckAromaticGUI},
##    {'name':'AD3tors_changeAromaticCarbonsToC','cmd':StopCheckAromatic(),
##        'gui':StopCheckAromaticGUI},
    {'name':'AD3tors_setCarbonNames','cmd':SetCarbonNames(),
        'gui':SetCarbonNamesGUI},
    {'name':'AD3tors_changePlanarityCriteria','cmd':ChangeAromaticCutOff(),
        'gui':ChangeAromaticCutOffGUI},
    {'name':'AD3tors_writeFormattedPDBQ','cmd':AUTOTORSWriter(),
        'gui':AUTOTORSWriterGUI},
    ]


def initModule(vf):

    for dict in commandList:
        vf.addCommand(dict['cmd'],dict['name'],dict['gui'])
    
    if not hasattr(vf, 'ADTSetMode'):
        vf.addCommand(AdtSetMode(), 'ADTSetMode')
    if not hasattr(vf, 'ADtors_limitTorsions'):
        vf.addCommand(SetTorsionNumber(), 'ADtors_limitTorsions')
    if not hasattr(vf, 'ADtors_setBondRotatableFlag'):
        vf.addCommand(SetBondRotatableFlag(), 'ADtors_setBondRotatableFlag')
    if not hasattr(vf, 'ADtors_stop'):
        vf.addCommand(StopAutoTors(), 'ADtors_stop')

    if vf.hasGui:
        vf.GUI.menuBars['AutoTools3Bar']._frame.config( {'background':'tan'})
        for item in vf.GUI.menuBars['AutoTools3Bar'].menubuttons.values():
            item.configure(background = 'tan')
        if not hasattr(vf.GUI, 'adt3Bar'):
            vf.GUI.adt3Bar = vf.GUI.menuBars['AutoTools3Bar']
            vf.GUI.adt3Frame = vf.GUI.adt3Bar.menubuttons.values()[0].master
        if not hasattr(vf.GUI, 'adt3ModeLabel'):
            mbs = {}
            packing_list = []
            for c in vf.GUI.adt3Frame.children.values():
                if isinstance(c, Tkinter.Menubutton):
                    mbs[c.cget('text')] = c
                    packing_list.append(c.cget('text'))
                    c.pack_forget()
            vf.GUI.adt3ModeLabel=Tkinter.Label(vf.GUI.adt3Frame, text="ADT3.0", width=6,
                             relief='sunken', borderwidth=1, fg='DarkBlue',
                             bg = 'ivory',anchor='w' )
            vf.GUI.adt3ModeLabel.pack(side='left')
            vf.GUI.adt3ModeLabel.bind("<Double-Button-1>", vf.ADTSetMode.guiCallback())
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
#            vf.GUI.ligandLabel=Tkinter.Label(vf.GUI.adtFrame, text="None", width=4,
#                                     relief='sunken', borderwidth=1,
#                                     anchor='w' )
#            vf.GUI.ligandLabel.pack(side='left')


