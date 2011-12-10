#############################################################################
#
# Author: Ruth HUEY, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################


# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/autoanalyze4Commands.py,v 1.7 2009/02/26 22:14:47 rhuey Exp $
#
# $Id: autoanalyze4Commands.py,v 1.7 2009/02/26 22:14:47 rhuey Exp $
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
from ViewerFramework.VFCommand import CommandGUI
from AutoDockTools.autoanalyzeCommands import menuText,\
checkHasInitializedDockings, hideShowHide, toggleShowHide,\
checkNameStr, ADChooseMacro, ADReadMacro, ADEPDBMol,\
ADSeeSpots, ADShowBindingSite, ADMakeAllGrids, ADGetOutput,\
ADGetAGrid, ADSelectDLG, ADDeleteDLG, ADGetDirDLGs, ADGetDLG,\
ClusterDockingChooser, ModelDockingChooser, ADDrawHistogram,\
ADMacroLigandChart, ADDockingChooser, ReadAutoDockStates,\
StatesPlayerWidget, ShowAutoDockStatesBaseCmd, ShowAutoDockStates,\
ShowAutoDockStatesByEnergy, ShowAutoDockPopulation,\
ShowAutoDockStatesHISTOGRAM, ShowAutoDockClusteringStates,\
ReadAutoDockClusteringStates, WriteAutoDockStates,\
WriteAutoDockClustering, MakeAutoDockCLUSTERING,\
MakeAutoDockSubsetCLUSTERING


ADChooseMacroGUI=CommandGUI()
ADChooseMacroGUI.addMenuCommand('AutoTools4Bar', menuText['AnalyzeMB'], 
        menuText['chooseMacro'], cascadeName = menuText['MoleculesMB'])


ADReadMacroGUI=CommandGUI()
ADReadMacroGUI.addMenuCommand('AutoTools4Bar', menuText['AnalyzeMB'], 
        menuText['readMacro'], cascadeName = menuText['MoleculesMB'])


ADEPDBMolGUI=CommandGUI()
ADEPDBMolGUI.addMenuCommand('AutoTools4Bar', menuText['AnalyzeMB'], 
        menuText['epdbMol'], cascadeName = menuText['GridsMB'])


ADSeeSpotsGUI=CommandGUI()
ADSeeSpotsGUI.addMenuCommand('AutoTools4Bar', menuText['AnalyzeMB'], 
        menuText['seeSpots'], cascadeName = menuText['DockingLogMB'])


ADShowBindingSiteGUI=CommandGUI()
ADShowBindingSiteGUI.addMenuCommand('AutoTools4Bar', menuText['AnalyzeMB'], 
        menuText['showBindingSite'], cascadeName = menuText['DockingLogMB'])


ADMakeAllGridsGUI=CommandGUI()
ADMakeAllGridsGUI.addMenuCommand('AutoTools4Bar', menuText['AnalyzeMB'], 
            menuText['showGridsMB'], cascadeName=menuText['GridsMB'])
                

ADGetOutputGUI=CommandGUI()
ADGetOutputGUI.addMenuCommand('AutoTools4Bar', menuText['AnalyzeMB'], 
        menuText['getOutputMB'] , cascadeName=menuText['StatesMB'])



ADGetAGridGUI=CommandGUI()
ADGetAGridGUI.addMenuCommand('AutoTools4Bar', menuText['AnalyzeMB'], 
        menuText['addGridMB'], cascadeName=menuText['GridsMB'])



ADSelectDLGGUI=CommandGUI()
ADSelectDLGGUI.addMenuCommand('AutoTools4Bar', menuText['AnalyzeMB'], 
        menuText['selectDLG'], cascadeName = menuText['DockingLogMB'])


ADDeleteDLGGUI=CommandGUI()
ADDeleteDLGGUI.addMenuCommand('AutoTools4Bar', menuText['AnalyzeMB'], 
        menuText['deleteDLG'], cascadeName = menuText['DockingLogMB'])


ADGetDirDLGsGUI=CommandGUI()
ADGetDirDLGsGUI.addMenuCommand('AutoTools4Bar', menuText['AnalyzeMB'],
        menuText['readDirDLG'], cascadeName = menuText['DockingLogMB'])


ADGetDLGGUI=CommandGUI()
ADGetDLGGUI.addMenuCommand('AutoTools4Bar', menuText['AnalyzeMB'],
        menuText['readDLG'], cascadeName = menuText['DockingLogMB'])
###ADGetDLGGUI.menuBarCfg.update({'background':'tan','relief':'sunken'})



ADDrawHistogramGUI=CommandGUI()
ADDrawHistogramGUI.addMenuCommand('AutoTools4Bar', menuText['AnalyzeMB'],
        menuText['showHistogramMB'], cascadeName=menuText['StatesMB'])



ADMacroLigandChartGUI=CommandGUI()
ADMacroLigandChartGUI.addMenuCommand('AutoTools4Bar', menuText['AnalyzeMB'],
            menuText['showChartMB'], cascadeName=menuText['StatesMB'])
            

ADDockingChooserGUI=CommandGUI()
ADDockingChooserGUI.addMenuCommand('AutoTools4Bar', menuText['AnalyzeMB'],
        menuText['chooseConfMB'], cascadeName = menuText['StatesMB'])


ReadAutoDockStatesGUI = CommandGUI()
ReadAutoDockStatesGUI.addMenuCommand('AutoTools4Bar', menuText['AnalyzeMB'],
        menuText['readStatesMB'],cascadeName=menuText['StatesMB'])


ShowAutoDockStatesGUI = CommandGUI()
ShowAutoDockStatesGUI.addMenuCommand('AutoTools4Bar', 
    menuText['AnalyzeMB'], menuText['showStatesMB'],
        cascadeName=menuText['StatesMB'])


ShowAutoDockStatesByEnergyGUI = CommandGUI()
ShowAutoDockStatesByEnergyGUI.addMenuCommand('AutoTools4Bar', 
    menuText['AnalyzeMB'], menuText['showStatesByEnergyMB'],
        cascadeName=menuText['StatesMB'])


ShowAutoDockPopulationGUI = CommandGUI()
ShowAutoDockPopulationGUI.addMenuCommand('AutoTools4Bar', 
    menuText['AnalyzeMB'], menuText['showPopulationMB'],
        cascadeName=menuText['StatesMB'])



ShowAutoDockStatesHISTOGRAMGUI = CommandGUI()
ShowAutoDockStatesHISTOGRAMGUI.addMenuCommand('AutoTools4Bar', 
    menuText['AnalyzeMB'], menuText['showStatesHISTOGRAMMB'], 
    cascadeName=menuText['StatesMB'])



ShowAutoDockStatesCLUSTERINGGUI = CommandGUI()
ShowAutoDockStatesCLUSTERINGGUI.addMenuCommand('AutoTools4Bar', 
    menuText['AnalyzeMB'], menuText['showStatesCLUSTERINGMB'], 
    cascadeName=menuText['ClusteringMB'])


ReadAutoDockStatesCLUSTERINGGUI = CommandGUI()
ReadAutoDockStatesCLUSTERINGGUI.addMenuCommand('AutoTools4Bar', 
    menuText['AnalyzeMB'], menuText['readStatesCLUSTERINGMB'], 
    cascadeName=menuText['ClusteringMB'])



WriteAutoDockStatesGUI = CommandGUI()
WriteAutoDockStatesGUI.addMenuCommand('AutoTools4Bar', 
    menuText['AnalyzeMB'], menuText['writeResultMB'], 
    cascadeName=menuText['StatesMB'])



WriteAutoDockClusteringGUI = CommandGUI()
WriteAutoDockClusteringGUI.addMenuCommand('AutoTools4Bar', 
    menuText['AnalyzeMB'], menuText['writeClusteringMB'], 
    cascadeName=menuText['ClusteringMB'])



MakeAutoDockCLUSTERINGGUI = CommandGUI()
MakeAutoDockCLUSTERINGGUI.addMenuCommand('AutoTools4Bar', 
    menuText['AnalyzeMB'], menuText['makeCLUSTERINGMB'], 
    cascadeName=menuText['ClusteringMB'])



MakeAutoDockSubsetCLUSTERINGGUI = CommandGUI()
MakeAutoDockSubsetCLUSTERINGGUI.addMenuCommand('AutoTools4Bar', 
    menuText['AnalyzeMB'], menuText['makeSubsetCLUSTERINGMB'], 
    cascadeName=menuText['ClusteringMB'])



commandList = [
    {'name':'AD4analyze_readDLG','cmd':ADGetDLG(),'gui':ADGetDLGGUI},
    {'name':'AD4analyze_readAllDLGInDirectory','cmd':ADGetDirDLGs(),'gui':ADGetDirDLGsGUI},
    {'name':'AD4analyze_selectDLG','cmd':ADSelectDLG(),'gui':ADSelectDLGGUI},
    {'name':'AD4analyze_deleteDLG','cmd':ADDeleteDLG(),'gui':ADDeleteDLGGUI},
    {'name':'AD4analyze_readMacromolecule','cmd':ADReadMacro(),'gui':ADReadMacroGUI},
    {'name':'AD4analyze_chooseMacromolecule','cmd':ADChooseMacro(),'gui':ADChooseMacroGUI},
    {'name':'AD4analyze_showDockingsAsSpheres','cmd':ADSeeSpots(),'gui':ADSeeSpotsGUI},
    {'name':'AD4analyze_showBindingSite','cmd':ADShowBindingSite(),'gui':ADShowBindingSiteGUI},
    #{'name':'AD4analyze_readStates','cmd':ReadAutoDockStates(),'gui':ReadAutoDockStatesGUI},
    {'name':'AD4analyze_showStates','cmd':ShowAutoDockStates(),'gui':ShowAutoDockStatesGUI},
    {'name':'AD4analyze_showStatesByEnergy','cmd':ShowAutoDockStatesByEnergy(),'gui':ShowAutoDockStatesByEnergyGUI},
    {'name':'AD4analyze_showPopulation','cmd':ShowAutoDockPopulation(),'gui':ShowAutoDockPopulationGUI},
    {'name':'AD4analyze_chooseDockedConformations','cmd':ADDockingChooser(),'gui':ADDockingChooserGUI},
    #{'name':'AD4analyze_showStatesHISTOGRAM','cmd':ShowAutoDockStatesHISTOGRAM(),'gui':ShowAutoDockStatesHISTOGRAMGUI},
    #{'name':'AD4analyze_showResultsOutput','cmd':ADGetOutput(),'gui':ADGetOutputGUI},
    #{'name':'AD4analyze_showHistogram','cmd':ADDrawHistogram(),'gui':ADDrawHistogramGUI},
    #{'name':'AD4analyze_getChart','cmd':ADMacroLigandChart(),'gui':ADMacroLigandChartGUI},
    #{'name':'AD4analyze_writeStates','cmd':WriteAutoDockStates(),'gui':WriteAutoDockStatesGUI},
    {'name':'AD4analyze_showClusteringStates','cmd':ShowAutoDockClusteringStates(),'gui':ShowAutoDockStatesCLUSTERINGGUI},
    #{'name':'AD4analyze_readClusteringStates','cmd':ReadAutoDockClusteringStates(),'gui':ReadAutoDockStatesCLUSTERINGGUI},
    {'name':'AD4analyze_makeClustering','cmd':MakeAutoDockCLUSTERING(),'gui':MakeAutoDockCLUSTERINGGUI},
    {'name':'AD4analyze_makeSubsetClustering','cmd':MakeAutoDockSubsetCLUSTERING(),'gui':MakeAutoDockSubsetCLUSTERINGGUI},
    #{'name':'AD4analyze_writeClustering','cmd':WriteAutoDockClustering(),'gui':WriteAutoDockClusteringGUI},
    ]

try:
    from Pmv.Grid import AutoGrid, AutoGridSurfaceGui
    for i in [ #{'name':'AD4analyze_epdbMolecule', 'cmd':ADEPDBMol(), 'gui':ADEPDBMolGUI},
    {'name':'AD4analyze_addExtraGridIsocontour','cmd':ADGetAGrid(),'gui':ADGetAGridGUI}, {'name':'AD4analyze_showGridIsocontours','cmd':ADMakeAllGrids(),'gui':ADMakeAllGridsGUI}]:
        commandList.insert(7,i)
except:
    print 'skipping the isocontour-dependent commands'


def initModule(vf):
    for dict in commandList:
        vf.addCommand(dict['cmd'],dict['name'],dict['gui'])
    #if not hasattr(vf, 'ADanalyze_showHistogram') and hasattr(vf, 'AD4analyze_showHistogram'):
    #	vf.ADanalyze_showHistogram = vf.AD4analyze_showHistogram
    if not hasattr(vf, 'ADanalyze_showDockingsAsSpheres') and hasattr(vf, 'AD4analyze_showDockingsAsSpheres'):
    	vf.ADanalyze_showDockingsAsSpheres = vf.AD4analyze_showDockingsAsSpheres
    if not hasattr(vf, 'ADanalyze_showGridIsocontours') and hasattr(vf, 'AD4analyze_showGridIsocontours'):
    	vf.ADanalyze_showGridIsocontours = vf.AD4analyze_showGridIsocontours
    if not hasattr(vf, 'ADanalyze_showBindingSite') and hasattr(vf, 'AD4analyze_showBindingSite'):
    	vf.ADanalyze_showBindingSite = vf.AD4analyze_showBindingSite
    if not hasattr(vf, 'ADanalyze_chooseDockedConformations') and hasattr(vf, 'AD4analyze_chooseDockedConformations'):
    	vf.ADanalyze_chooseDockedConformations = vf.AD4analyze_chooseDockedConformations
    if not hasattr(vf, 'ADanalyze_readDLG') and hasattr(vf, 'AD4analyze_readDLG'):
    	vf.ADanalyze_readDLG = vf.AD4analyze_readDLG
    if not hasattr(vf, 'ADanalyze_selectDLG') and hasattr(vf, 'AD4analyze_selectDLG'):
    	vf.ADanalyze_selectDLG = vf.AD4analyze_selectDLG
    if not hasattr(vf, 'ADanalyze_makeSubsetClustering') and hasattr(vf, 'AD4analyze_makeSubsetClustering'):
    	vf.ADanalyze_makeSubsetClustering = vf.AD4analyze_makeSubsetClustering

    if hasattr(vf, 'GUI'):
        for item in vf.GUI.menuBars['AutoTools4Bar'].menubuttons.values():
            item.configure(background = 'tan')
        if not hasattr(vf.GUI, 'adtBar'):
            vf.GUI.adtBar = vf.GUI.menuBars['AutoTools4Bar']
            vf.GUI.adtFrame = vf.GUI.adtBar.menubuttons.values()[0].master
            





