#############################################################################
#
# Author: Ruth HUEY, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/autodpf3Commands.py,v 1.1 2008/06/02 23:55:43 rhuey Exp $
#
# $Id: autodpf3Commands.py,v 1.1 2008/06/02 23:55:43 rhuey Exp $
#
#
#
#
#
#

"""
This Module facilitates producing a docking parameter file for AutoDock. The steps in this process are:

    * Selecting the macromolecule filename: The user can select the macromolecule for autodpf in three ways: it can be chosen from molecules previously added to the moleculeViewer, it can be picked as a PDB file,  or it can be picked as a MOL2 file:

        o Choose Macromol...

        o Select PDB Macromolecule 

        o Select MOL2 Macromolecule

    * Selecting the small molecule which has been previously formatted by AutoTors: 

        o Via Reading a PDBQ-File which adds the ligand to the viewer

    * The user sets parameters pertaining to the small molecule 

        o Checking that a grid map exists for each of the ligand atom types 

        o Indicating whether a floating grid map exists

        o Setting the initial translation of the small molecule

            - by choosing  the 'random' option which sets a random starting position for the ligand

            - by entering the desired coordinates in the entry

        o Setting the initial quaternion of the small molecule

            - by choosing  the 'random' option which sets a random starting quaternion.

            - by entering the desired initial quaternion -Qx,Qy,Qz,Qw in the entry.  Qx, Qy, Qz define the unit vector in the direction of rigid body rotation and Qw the angle of rotation about this unit vector.

        o Setting the coefficient of the torsional DOF

        o By choosing to set the initial dihedrals for the small molecule or not: If not, AutoDock assumes that the chi1, chi2, chi3 etc are all zero and does not change the initial ligand torsion angles. If the user chooses to set the initial dihedrals, he further chooses:

            - for them to be randomly assigned 

            - an initial relative dihedral angle for each active torsion in the ligand.

        o The user can specify two types of torsion constraints for the ligand:

            -  Gaussian constraints which use an inverted Gaussian bell curve to calculate the energy function input of the constraint.  This type of constraint is specified by two floating point numbers: the perferred angle in the range -180-+180decreeds and the half-width which is the difference between two angles at which the energy is half the barrier PLUS an integer which identifies the torsion according to the list at the top of the AutoTors-generated input ligand PDBQ file. More than one constraint of this type may be specified for a single torsion.

            - Hard torsion constraints may also be specified. These differ from the previous type in that the torsion is never allowed to take values bewond the range defined and in that the second parameter is the full width of the allowed range of torsion angles. Moreover, only one constraint of this type is allowed per torsion.

        o If the user specifies torsion constraints, he may also specify the height of the energy barrier to be applied to these constraints.

        o If the user specifies Gaussian torsion constraints, he may also specify whether to store and output the torsion energies

    * The user sets parameters pertaining to docking algorithm(s) he wishes to use
:
        o Setting Simulated Annealing parameters.

        o Setting Genetic Algorithm parameters (GA).

        o Setting Local Search parameters (LS).

    It is important to remember that any of these may be used alone but only GA and LS may be used together


    * The user adjusts these additional parameters: 
    
        o the step sizes of translation, quaternion rotation and dihedral torsion change.
        o  energy parameters including energy assigned to atoms outside the grid volume, the maximum allowable initial energy and the maximum number of retries.

        o output format parameters including the level of detail for the output, the rms cluster tolerance, the reference file for rms calculations and whether to do symmetry checking in the rms calculations.


    * The user selects which kind of docking parameter file to write : 
    
        o Simulated Annealing 

        o GA

        o LS

        o GALS


    * The results of the previous steps are written to a file. The user selects a filename via a filebrowser.  By convention, the file should have a .dpf extension. If no macromolecule has been selected, it is not possible to write a grid parameter file and the user gets a warning message to that effect. Likewise, the types of the maps to be calculated must be set before the grid parameter file is written and a warning message to this effect appears if the types have not been set.
(A checkbutton, "DONE", allows the user to withdraw the autoTools menuBar)
    
"""
from ViewerFramework.VFCommand import CommandGUI


from AutoDockTools.autodpfCommands import DpfSetDpo, DpfLoadDefaults,\
DpfMacroSelector, DpfMacroChooser, DpfInitLigand, DpfLigandChooser,\
DpfLigPDBQReader, DpfEditor, DpfSAWriter, DpfGAWriter, DpfLSWriter, DpfGALSWriter, \
DpfClusterWriter, SimAnneal, GA, LS, SetDockingRunParms, StopAutoDpf,\
menuText, checkHasDpo, sortKeyList, setDpoFields



DpfLoadDefaultsGUI = CommandGUI()
DpfLoadDefaultsGUI.addMenuCommand('AutoTools3Bar', menuText['AutoDpfMB'], menuText['ReadDpfMB'])


DpfMacroSelectorGUI=CommandGUI()
DpfMacroSelectorGUI.addMenuCommand('AutoTools3Bar', menuText['AutoDpfMB'],\
menuText['ReadMacro'], cascadeName = menuText['MacromoleculeMB'])


DpfMacroChooserGUI=CommandGUI()
DpfMacroChooserGUI.addMenuCommand('AutoTools3Bar', menuText['AutoDpfMB'],\
menuText['ChooseMacro'], cascadeName = menuText['MacromoleculeMB'])


DpfInitLigandGUI=CommandGUI()
DpfInitLigandGUI.addMenuCommand('AutoTools3Bar', menuText['AutoDpfMB'], menuText['AdjustLigand'], cascadeName = menuText['SetLigandParmsMB'])


DpfLigandChooserGUI=CommandGUI()
DpfLigandChooserGUI.addMenuCommand('AutoTools3Bar', menuText['AutoDpfMB'],\
menuText['ChooseLigand'], cascadeName = menuText['SetLigandParmsMB'])


        
DpfLigPDBQReaderGUI = CommandGUI()
DpfLigPDBQReaderGUI.addMenuCommand('AutoTools3Bar', menuText['AutoDpfMB'],\
menuText['ReadLigand'],  cascadeName = menuText['SetLigandParmsMB'])


DpfEditorGUI=CommandGUI()
DpfEditorGUI.addMenuCommand('AutoTools3Bar', menuText['AutoDpfMB'], menuText['EditDpfMB'])



DpfSAWriterGUI=CommandGUI()
DpfSAWriterGUI.addMenuCommand('AutoTools3Bar', menuText['AutoDpfMB'],\
menuText['WriteSA'], cascadeName = menuText['WriteDpfMB'])



DpfGAWriterGUI=CommandGUI()
DpfGAWriterGUI.addMenuCommand('AutoTools3Bar', menuText['AutoDpfMB'], menuText['WriteGA'],\
cascadeName = menuText['WriteDpfMB'])



DpfLSWriterGUI=CommandGUI()
DpfLSWriterGUI.addMenuCommand('AutoTools3Bar', menuText['AutoDpfMB'], menuText['WriteLS'],\
cascadeName = menuText['WriteDpfMB'])



DpfGALSWriterGUI=CommandGUI()
DpfGALSWriterGUI.addMenuCommand('AutoTools3Bar', menuText['AutoDpfMB'],\
menuText['WriteGALS'], cascadeName = menuText['WriteDpfMB'])



DpfClusterWriterGUI=CommandGUI()
DpfClusterWriterGUI.addMenuCommand('AutoTools3Bar', menuText['AutoDpfMB'],\
menuText['WriteCluster'], cascadeName = menuText['WriteDpfMB'])



SimAnnealGUI=CommandGUI()
SimAnnealGUI.addMenuCommand('AutoTools3Bar', menuText['AutoDpfMB'],menuText['SA'], cascadeName = menuText['SetSearchParmsMB'])



GAGUI=CommandGUI()
GAGUI.addMenuCommand('AutoTools3Bar', menuText['AutoDpfMB'],menuText['GA'], cascadeName = menuText['SetSearchParmsMB'])



LSGUI=CommandGUI()
LSGUI.addMenuCommand('AutoTools3Bar', menuText['AutoDpfMB'],'Local Search Parameters ', cascadeName = menuText['SetSearchParmsMB'])



SetDockingRunParmsGUI=CommandGUI()
SetDockingRunParmsGUI.addMenuCommand('AutoTools3Bar', menuText['AutoDpfMB'],
    menuText['SetDockingRunParmsMB'])



commandList = [
    {'name':'AD3dpf_read','cmd':DpfLoadDefaults(),'gui':DpfLoadDefaultsGUI},
    {'name':'AD3dpf_chooseMacromolecule','cmd':DpfMacroChooser(),'gui':DpfMacroChooserGUI},
    {'name':'AD3dpf_readMacromolecule','cmd':DpfMacroSelector(),'gui':DpfMacroSelectorGUI},
    #ligand
    {'name':'AD3dpf_chooseFormattedLigand','cmd':DpfLigandChooser(),'gui':DpfLigandChooserGUI},
    {'name':'AD3dpf_readFormattedLigand','cmd':DpfLigPDBQReader(),'gui':DpfLigPDBQReaderGUI},
    {'name':'ADdpf_initLigand','cmd':DpfInitLigand(),'gui':DpfInitLigandGUI},

    {'name':'AD3dpf_setGAparameters','cmd':GA(),'gui':GAGUI},
    {'name':'AD3dpf_setSAparameters','cmd':SimAnneal(),'gui':SimAnnealGUI},
    {'name':'AD3dpf_setLSparameters','cmd':LS(),'gui':LSGUI},
    {'name':'AD3dpf_setDockingParameters','cmd':SetDockingRunParms(),'gui':SetDockingRunParmsGUI},
    {'name':'AD3dpf_writeGALS','cmd':DpfGALSWriter(),'gui':DpfGALSWriterGUI},
    {'name':'AD3dpf_writeGA','cmd':DpfGAWriter(),'gui':DpfGAWriterGUI},
    {'name':'AD3dpf_writeSA','cmd':DpfSAWriter(),'gui':DpfSAWriterGUI},
    {'name':'AD3dpf_writeLS','cmd':DpfLSWriter(),'gui':DpfLSWriterGUI},
    #{'name':'AD3dpf_writeCluster','cmd':DpfClusterWriter(),'gui':DpfClusterWriterGUI},
    {'name':'AD3dpf_edit','cmd':DpfEditor(),'gui':DpfEditorGUI},
]



def initModule(vf):

    for dict in commandList:
        vf.addCommand(dict['cmd'], dict['name'], dict['gui'])
    if not hasattr(vf, 'ADdpf_setDpo'):
        vf.addCommand(DpfSetDpo(), 'ADdpf_setDpo', None)

    if vf.hasGui:
        for item in vf.GUI.menuBars['AutoTools3Bar'].menubuttons.values():
            item.configure(background = 'tan')
        if not hasattr(vf.GUI, 'adtBar'):
            vf.GUI.adtBar = vf.GUI.menuBars['AutoTools3Bar']
            vf.GUI.adtFrame = vf.GUI.adtBar.menubuttons.values()[0].master




        
