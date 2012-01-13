modlist = [
("AutoDockTools","autoanalyzeCommands","""
This Module facilitates analyzing results of autodock jobs. 

    * The first step is 'Read Docking Log'  The selected file is parsed which sets items in the dictionary 'docked', an attribute of the molecular viewer. The selected file is also entered in the dictionary 'dockings' as a separate item. The selected file is appended to mv.dockings['dlgFiles'] and its short filename to mv.dockings['dlgEntries']. These keys are set in mv.docked:
            
        o 'macroFile': the Macromolecule file used

        o 'ligand':  the original ligand 

        o 'types': the kinds of grid files used 

        o 'runs': the number of docking runs 

        o 'clusterNum': the number of clusters produced

        o 'clusterList': a list of ADClusters which have members of ADDocked instances. 

        o 'rmsd-tolerance': rmsd value used in this docking for clustering

        o 'output': lines containing summary of docking

        o 'ligResList': ordered list of residues found in input-ligand (obsolete)

        o 'printStart': line number in dlg for beginning of output

        o 'printEnd': line number in dlg for end of output

        o 'modelList': a list of docked conformtions when no clustering was done (??) 

        o 'macro': filename of macromolecule (eg '1hvrCorr.pdbqs')

        o 'macroStem': name of macromolecule up to last '.' (eg '1hvrCorr')

        o 'dlg': full pathname of dlg

        o 'ligAllLines': list of all lines in dlg for input ligand including autotors keyword lines

After the selected docking log file is parsed, the user can:

    * select a displayed docked conformation using the 'Choose A Docked Confomration' menubutton.  This opens a DockingChooser widget which is a ListChooser allowing selection either in the widget or in the viewer of any of the displayed docking. Information about each docked conformation is displayed in the information window of the DockingChooser as different entries are high-lighted.  

    * display the macromolecule via the "Show Macromolecule" menubutton.  This menubutton is linked to a file browsers in case the macromolecule whose name is parsed from the docking log file is not in the current directory. (FIX THIS: what if the macromolecule is in a different directory but there is a molecule with the same name here???). The user can change the visibility, sampling, isovalue, renderMode and visibility of bounding box  for each of  the displayed grids

    * display the autogrids used in the docking via the "Show Grids Used For Calc" menubutton.  This menubutton is linked to a ListChooser which lets the user select whether to load all or some of the grids. The user can interactively change the visibility of each grid's isosurface, its sampling value, its isovalue, its rendermode (LINE or FILL) and the visibility of its bounding box. 

    * The user is able to visualize extra grid maps using the "Show Grid" button. 

    * If the current docking has clusters, the user is able to visualize a results histogram for it with 'Show Histogram'. The histogram can be printed.

    * Result Summaries for docking(s) can be viewed, edited and saved with 'Get Output'

    * Dockings can be deleted via 'Delete Docking Log'

"""),
("AutoDockTools","autodpfCommands","""
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
    
"""),
("AutoDockTools","autogpfCommands","""
This Module facilitates producing a grid parameter file for AutoGrid. The steps in this process are:

    * 'Macromolecule': Selecting the macromolecule: 
        The user can select the macromolecule for autogpf in two ways: 
           - it can be chosen from molecules previously added to the moleculeViewer  
           - it can be read in from a file:

        o Choose Macromol...

        o Read Macromolecule 


    * 'Set Map Types': Setting the types of maps to generate: 

        o Set Map Types Directly

        o By Choosing Ligand

        o By Reading Formatted File

The user can change the types of maps to be calculated.
He decides which types of possible hydrogen bonding he wishes to model. 
For instance, IF hydrogens are present  AND nitrogens, oxygens and /or sulfurs, 
the user can decide to model N-H bonds, O-H bonds and/or S-H bonds.  
He sets which type of dielectric to use:
    -distance-dependent dielectric  
    -constant dielectric  
(Other ligand-related commands allow the user to set energy parameters for new 
atom types or to set up a specialized 'covalent' grid-map.)


    * 'Set Grid': The user positions the grid and sets its dimensions by:

        o Setting the center of the grid maps: 

            - by picking an atom or

            - by entering the full-name of an atom or 

            - by entering the desired coordinates in entries 'x center', 'y center', 
'z center' (NB: ALL entries must be 'activated' by a 'Return')

            - by choosing  the 'Center on Macromolecule' option which sets the 
center of the grid to the geometric center of the macromolecule (obtained by 
averaging all its coordinates)

            - by choosing  the 'Center on Ligand' option which sets the center of 
the grid to the geometric center of the ligand (obtained by averaging all its 
coordinates)

        o Setting the number of grid points in each direction (which has to be an 
even number) and the spacing between the points. This is done by using the 
corresponding scale widgets.

        o Adjusting the position of the grid using scales for x-offset, y-offset 
and z-offset.  These scales allow the user to move the grid box up to 10 angstroms 
in any direction along any of the three axes. 
(NOTE that the units of these scales are tenths of Angstroms and the new coordinates 
of the center are reflected in the x-center, y-center, z-center entries)

    * 'Set Other Options': The user adjusts these additional parameters: 
    
        o the smoothing factor can be changed from its default 0.5Angstrom value.  
This changes the radius of the area within which the minimum energy is stored.
        o  electrostatic potential map may or may not be generated by AutoGrid

        o floating point potential map may or may not be generated 

        o the user may decide whether or not to use the default distance dependent 
dielectric constant.  If not, he can enter his desired dielectric constant or use 
the default value, 40. It should be noted that this entered value is multiplied 
by 0.1146 by the program for input to AutoGrid.

    * 'Write GPF': The results of the previous steps are written to a file. 
The user selects a filename via a filebrowser.  By convention, the file should 
have a .gpf extension. If no macromolecule has been selected, it is not possible 
to write a grid parameter file and the user gets a warning message to that effect. 
Likewise, the types of the maps to be calculated must be set before the grid 
parameter file is written and a warning message to this effect appears if the 
types have not been set.

    * 'Edit GPF': Allows user to edit a grid parameter file.  If one has been
written, it is automatically loaded. Otherwise, the user can select any *.gpf
file to edit from a file browser.
    
"""),
("AutoDockTools","autostartCommands","""
This Module facilitates starting autogrid and autodock jobs and managing them

"""),
("AutoDockTools","autotorsCommands","""
This Module facilitates selecting and formatting a ligand for a subsequent 
AutoDock run.  The steps in this process are:

    * The user selects the small molecule from a list of molecules 
already in the moleculeViewer OR as a PDBQ file or as a MOL2 file from 
a fileBrowser.  

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
is able to select a carbon to convert (reversibly).  

    * Non-polar hydrogens can be merged which means that the charge of 
each is added to its carbon and the hydrogen atoms themselves are not written 
in the output file, thus in some sense 'removing' them from the molecule. 
'Fewer' atoms simplifies the AutoDock run.

    * The last function of this Module is to write a file which contains 
the correctly formatted ligand atoms.  The ROOT section of the molecule 
expands from the selected ROOT atom out to include all atoms adjacent to it 
up to the first active torsion.  The active torsions set the position of 
BRANCH and TORS key words in the output pdbq file (and their corresponding 
ENDBRANCH and ENDTORS key words). These keywords are nested to set up  a 
Depth-First Order Traversal.  Autotors also calculates the torsional degrees 
of freedom (TORSDOF) which is the number of possible torsions less the number of 
symmetry-equivalent torsions (such as a bond to a NH3). This key word is the 
last line of the pdbq file. 
"""),
("AutoDockTools", "autopilotCommands", ""),
]












