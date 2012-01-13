#!/usr/bin/env python
#
# This is a module in development for making conformationally averaged PBE maps
#
debug=False
import sys, os

print __file__
import os
try:
    file_name=__file__
    if file_name[:2]=='./':
        scriptpath=os.getcwd()
    else:
        scriptpath=os.path.join(os.getcwd(),os.path.split(file_name)[0])
        if scriptpath[-1] == "/":
            scriptpath=scriptpath[:-1]
except:
    scriptpath=os.path.split(sys.argv[0])[0]
    if scriptpath=='.':
        scriptpath=os.getcwd()
#
# Add to import path
#
pdb2pqr_path=os.path.split(scriptpath)[0]
sys.path.append(pdb2pqr_path)

import string
import math
import string
import getopt
import time
from src.pdb import *
from src.utilities import *
from src.structures import *
from src.definitions import *
from src.forcefield import *  
from src.routines import *
from src.protein import *
from src.server import *
from StringIO import *
from src.hydrogens import *

class conf_avg:

    def __init__(self,options):
        """Initialize class and decide which kind of job to do"""
        self.options=options # Store options so we can access them anywhere
        potentials=[]
		# If directoryPath is specified then use that, otherwise use pdbfilename
        if options.directoryPath!='':
            listOfFiles=os.listdir(options.directoryPath)
            for currentPDB in listOfFiles:
                currentPDB=os.path.join(options.directoryPath,currentPDB) # Jens added this fix
                pots=self.process_one_pdb(currentPDB)
                potentials.append(pots)
        else:
            # Single file
            potentials.append(self.process_one_pdb(os.path.join(os.getcwd(),options.pdbfilename)))
		#
		# Average potentials
		#
        avg_pots=self.average_potentials(potentials)
        return
        
    #
    # ------
    #
    
    def process_one_pdb(self,pdbfilename):
        """Do everything for one input file"""
        print "Working on: %s" %pdbfilename
        pdbfile = getPDBFile(pdbfilename)
        
        if self.options.MD:
            #
            # Run an MD simulation for this PDB file and calculate potentials for all the snapshots
            #
            snapshots=self.run_MD(pdbfilename)
        else:
            snapshots=[pdbfilename]
        #
        # Get the potentials for everything
        #
        potentials=[]
        for pdbname in snapshots:
            pots=self.get_potentials(pdbfilename)
            potentials.append(pdbname)
        return potentials
        
    #
    # ------
    #
    
    def run_MD(self,pdbfilename):
        """Run an MD simulation and return a number of snapshots"""
        files=os.listdir(os.getcwd())
        addfiles=[pdbfilename]
        #for file in files:
            #addfiles.append(os.path.join(os.getcwd(),file))
        import Gromacs_class as Gclass
        G=Gclass.GROMACS(addfiles)
        #
        # Create Gromacs input file
        # 
        pdbfile=os.path.split(pdbfilename)[1]
        G.pdb2gmx(pdbfile,forcefield=1,ignore_Hs=True,auto_select_his=True)
        #
        # Set up the simulation box. The argument gives the distance between the box edges
        # and the protein in nm
        #
        G.center(1.5)
        G.solvate()
        #
        # Energy minimise and do the pre-MD
        #
        G.EM(2000,1000)
        params={}
        params['ref_t']='%d' %options.temperature
        G.PR_MD(250)
        # ----------------------------------
        # 500 ps timestep is 2 fs
        G.MD(self.options.MDtime*500,params)
        #
        # Get the snapshots
        #
        filenames=G.get_snapshots(self.options.numsnapshots)
        return filenames
        
        
    def get_potentials(self,currentPDB):
		"""Get the potentials by first running pdb2pqr and then apbs"""
		myProtein,apbs_inputfile=self.run_pdb2pqr(currentPDB)
		potentials=self.run_apbs(myProtein,apbs_inputfile)
		return potentials
        
        
    def run_pdb2pqr(self,currentPDB):
        """Run pdb2pqr, prepare input for apbs"""
        pdbfile = getPDBFile(currentPDB)
        pdblist, errlist = readPDB(pdbfile)
        #
        # Instantiate pdb2pqr
        #
        myDefinition = Definition()
        myProtein = Protein(pdblist, myDefinition)

        #
        # Setup everything
        #
        myRoutines = Routines(myProtein, verbose)
        myRoutines.updateResidueTypes()
        myRoutines.updateSSbridges()
        myRoutines.updateBonds()
        myRoutines.setTermini()
        myRoutines.updateInternalBonds()

        myforcefield=Forcefield(ff, myDefinition, None)
        myRoutines.applyNameScheme(myforcefield)

        myRoutines.findMissingHeavy()
        myRoutines.addHydrogens()
        myRoutines.debumpProtein()
        myProtein.reSerialize()
        #
        # Add and optimze hydrogens:
        # 
        from src.hydrogens import hydrogenRoutines
        myRoutines.updateInternalBonds()
        myRoutines.calculateDihedralAngles()
        myhydRoutines = hydrogenRoutines(myRoutines)
        #
        # Now optimize hydrogens
        #
        myhydRoutines.setOptimizeableHydrogens()
        myhydRoutines.initializeFullOptimization()
        myhydRoutines.optimizeHydrogens()
        myhydRoutines.cleanup()
        myRoutines.setStates()

        print "Created protein object (after processing myRoutines) -"
        print "\tNumber of residues in protein: %s" % myProtein.numResidues()
        print "\tNumber of atoms in protein   : %s" % myProtein.numAtoms()

        #
        # Assign charges
        #
        for chain in myProtein.getChains():
            for residue in chain.get("residues"):
                for atom in residue.get("atoms"):
                    atomname = atom.get("name")
                    charge, radius = myforcefield.getParams1(residue, atomname)
                    atom.set("radius", radius)
                    atom.set("ffcharge", charge)
        #
        #
		method=""
		async=0
		split=0
		import pdb2pka.inputgen_pKa as IP
		igen = IP.inputGen(currentPDB)
		igen.maps=None
		igen.set_type('background')
		igen.pdie=8.0
		igen.sdie=80.0
		all_center,extent=igen.getCenter()
		igen.setfineCenter(all_center)
		print 'Center: %5.1fA %5.1fA %5.1fA' %(all_center[0],all_center[1],all_center[2])
		print 'Extent: %5.1fA %5.1fA %5.1fA'  %(extent[0],extent[1],extent[2])

		apbs_inputfile=igen.printInput()
		return myProtein, apbs_inputfile
        
    def run_apbs(self,myProtein,apbs_inputfile):
		"""runs apbs"""
		import pdb2pka.apbs 
		APBS=pdb2pka.apbs.runAPBS()
		potentials = APBS.runAPBS(myProtein, apbs_inputfile)
		APBS.cleanup()
		return potentials
        
    def average_potentials(self,potentials):
        """This function averages many potential maps"""
        avg_pots=[]
        for i in range(0,len(potentials[0])):
            currSum=0
            for j in range(0,len(potentials)):
                currSum+=potentials[j][i]
            currAvg=currSum/len(potentials)
            avg_pots.append(currAvg)

        print avg_pots
        return avg_pots

#
# ----
#

if __name__=='__main__':
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options]',version='%prog 1.0')
    parser.add_option('-p','--pdb',dest='pdbfilename',action='store',type='string',default='2lzt.pka.pdb',help='The PDB file. Default: %default')
    parser.add_option('-d','--dir',dest='directoryPath',action='store',type='string',default='',
                  help='Direcotry of the PDB files/snapshots. Default: %default')
    #
    # Flags
    #
    parser.add_option('--MD',dest='MD',action='store_true',default=False,help='Perform an MD simulation and use snapshots for calculating electrostatic potential')
    parser.add_option('-s','--MDsnapshots',dest='numsnapshots',action='store',type='int',default=100,help='Number of MD snapshots to use. Default: %default')
    parser.add_option('--MDtime',dest='MDtime',action='store',type='int',default=100,help='Time in picoseconds that MD should be run for. Default: %default')
    parser.add_option('-t','--temp',dest='temperature',action='store',type='float',default=310.15,help='Temperature for the MD run. Default: %default')
    #
    # We can think about adding flags for not solvating the structure etc here
    #



    (options, args) = parser.parse_args()

    verbose=True
    ff='parse'

    I=conf_avg(options)

