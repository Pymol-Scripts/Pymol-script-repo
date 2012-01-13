"""
    APBS interface for PDB2PQR

    Todd Dolinsky (todd@ccb.wustl.edu)
    Washington University in St. Louis

    Jens Erik Nielsen

"""

__date__  = "16 August 2005"
__author__ = "Todd Dolinsky, Jens Erik Nielsen"

import sys
import time
import string
from src import psize
from src import inputgen

from apbslib import *


Python_kb = 1.3806581e-23
Python_Na = 6.0221367e+23
NOSH_MAXMOL = 20
NOSH_MAXCALC = 20


class APBSError(Exception):
    """ APBSError class

        The APBSError class inherits off the Exception module and returns
        a string defining the nature of the error. 
    """
    
    def __init__(self, value):
        """
            Initialize with error message

            Parameters
                value:  Error Message (string)
        """
        self.value = value
        
    def __str__(self):
        """
            Return the error message
        """
        return `self.value`

def getUnitConversion():
    """
        Get the unit conversion from kT to kJ/mol

        Returns
            factor: The conversion factor (float)
    """
    temp = 298.15
    factor = Python_kb/1000.0 * temp * Python_Na
    return factor

def runAPBS(PQR, INPUT):
    """
        Run APBS using PQR and INPUT strings!
    """


    # Initialize the MALOC library
    startVio()

    # Initialize variables, arrays
    com = Vcom_ctor(1)
    rank = Vcom_rank(com)
    size = Vcom_size(com)
    mgparm = MGparm()
    pbeparm = PBEparm()
    mem = Vmem_ctor("Main")
    pbe = new_pbelist(NOSH_MAXMOL)
    pmg = new_pmglist(NOSH_MAXMOL)
    pmgp = new_pmgplist(NOSH_MAXMOL)
    realCenter = double_array(3)
    totEnergy = []
    x = []
    y = []
    z = []
    chg = []
    rad = []
   
    
    # Start the main timer
    main_timer_start = time.clock()



    # Parse the input file
    nosh = NOsh_ctor(rank, size)

    if not parseInputFromString(nosh, INPUT):
        stderr.write("main:  Error while parsing input file.\n")
        raise APBSError, "Error occurred!"

    # Load the molecules using Valist_load routine

    alist = new_valist(NOSH_MAXMOL)
    #atoms = protein.getAtoms()
    #protsize = len(atoms)
    atoms = string.split(PQR,"\n")
    for i in range(len(atoms)):
        atom = atoms[i]
        params = string.split(atom)
        x.append(float(params[5]))
        y.append(float(params[6]))
        z.append(float(params[7]))
        chg.append(float(params[8]))
        rad.append(float(params[9]))
  
    myAlist = make_Valist(alist,0)
    Valist_load(myAlist, len(atoms), x,y,z,chg,rad) 

    # Initialize the energy holders

    for i in range(nosh.ncalc): totEnergy.append(0.0)
    potList = []
    
    # Initialize the force holders
    forceList = []
  
    # Load the dieletric maps

    dielXMap = new_gridlist(NOSH_MAXMOL)
    dielYMap = new_gridlist(NOSH_MAXMOL)
    dielZMap = new_gridlist(NOSH_MAXMOL)
  
 
 
    # Load the kappa maps
    kappaMap = new_gridlist(NOSH_MAXMOL)
  

    # Load the charge maps
    chargeMap = new_gridlist(NOSH_MAXMOL)
   
    
    # Do the calculations
 
    sys.stdout.write("Preparing to run %d PBE calculations. \n" % nosh.ncalc)
   
    for icalc in xrange(nosh.ncalc):
        sys.stdout.write("---------------------------------------------\n")
        calc = NOsh_getCalc(nosh, icalc)
        mgparm = calc.mgparm
        pbeparm = calc.pbeparm
        if calc.calctype != 0:
            sys.stderr.write("main:  Only multigrid calculations supported!\n")
            raise APBSError, "Only multigrid calculations supported!"

        for k in range(0, nosh.nelec):
            if NOsh_elec2calc(nosh,k) >= icalc:
                break

        name = NOsh_elecname(nosh, k+1)
        if name == "":
            sys.stdout.write("CALCULATION #%d:  MULTIGRID\n" % (icalc+1))
        else:
            sys.stdout.write("CALCULATION #%d (%s): MULTIGRID\n" % ((icalc+1),name))
        sys.stdout.write("Setting up problem...\n")
	
        # Routine initMG
	
        if initMG(icalc, nosh, mgparm, pbeparm, realCenter, pbe, 
              alist, dielXMap, dielYMap, dielZMap, kappaMap, chargeMap, 
              pmgp, pmg) != 1:
            sys.stderr.write("Error setting up MG calculation!\n")
            raise APBSError, "Error setting up MG calculation!"
	
        # Print problem parameters 

        printMGPARM(mgparm, realCenter)
        printPBEPARM(pbeparm)
        # Solve the problem : Routine solveMG
	
        thispmg = get_Vpmg(pmg,icalc)

        if solveMG(nosh, thispmg, mgparm.type) != 1:
            stderr.write("Error solving PDE! \n")
            raise APBSError, "Error Solving PDE!"

        # Set partition information : Routine setPartMG

        if setPartMG(nosh, mgparm, thispmg) != 1:
            sys.stderr.write("Error setting partition info!\n")
            raise APBSError, "Error setting partition info!"
	
        ret, totEnergy[icalc] = energyMG(nosh, icalc, thispmg, 0,
                                         totEnergy[icalc], 0.0, 0.0, 0.0)
	
        # Set partition information
	
        # Write out data from MG calculations : Routine writedataMG	
        writedataMG(rank, nosh, pbeparm, thispmg)
	
        # Write out matrix from MG calculations	
        writematMG(rank, nosh, pbeparm, thispmg)

        # GET THE POTENTIALS
              
        potentials = getPotentials(nosh, pbeparm, thispmg, myAlist)
        potList.append(potentials)
        
    # Handle print statements

    if nosh.nprint > 0:
        sys.stdout.write("---------------------------------------------\n")
        sys.stdout.write("PRINT STATEMENTS\n")
    for iprint in xrange(nosh.nprint):
        if NOsh_printWhat(nosh, iprint) == NPT_ENERGY:
            printEnergy(com, nosh, totEnergy, iprint)
        elif NOsh_printWhat(nosh, iprint) == NPT_FORCE:
            printForce(com, nosh, nforce, atomforce, iprint)
        else:
            sys.stdout.write("Undefined PRINT keyword!\n")
            break

    sys.stdout.write("----------------------------------------\n")
    sys.stdout.write("CLEANING UP AND SHUTTING DOWN...\n")

    # Clean up APBS structures
    
    #killForce(mem, nosh, nforce, atomforce)
    killEnergy()
    killMG(nosh, pbe, pmgp, pmg)
    killChargeMaps(nosh, chargeMap)
    killKappaMaps(nosh, kappaMap)
    killDielMaps(nosh, dielXMap, dielYMap, dielZMap)
    killMolecules(nosh, alist)
    delete_Nosh(nosh)

    # Clean up Python structures

    #ptrfree(nfor)
    delete_double_array(realCenter)
    #delete_int_array(nforce)
    #delete_atomforcelist(atomforce)
    delete_valist(alist)
    delete_gridlist(dielXMap)
    delete_gridlist(dielYMap)
    delete_gridlist(dielZMap)
    delete_gridlist(kappaMap)
    delete_gridlist(chargeMap)
    delete_pmglist(pmg)
    delete_pmgplist(pmgp)
    delete_pbelist(pbe)
    
    
    # Clean up MALOC structures
    delete_Com(com)
    delete_Mem(mem)
    
    sys.stdout.write("\n")
    sys.stdout.write("Thanks for using APBS!\n\n")

    # Stop the main timer
    main_timer_stop = time.clock()
    sys.stdout.write("Total execution time:  %1.6e sec\n" % (main_timer_stop - main_timer_start))

    #Return potentials

    return potList

if __name__ == "__main__":

    # This would really be something like
    #   atoms = protein.getAtoms()
    #   PQR = protein.printAtoms(atoms)
    PQR = "ATOM      1  I   ION     1       0.000   0.000  0.000  1.00  3.00"

    size = psize.Psize()
    size.parseString(PQR)
    size.setAll()
   
    # The actual name doesn't matter since we're loading from the string!
    input = inputgen.Input("dummy.pqr", size, "mg-auto", 0)

    # Print out the number of elec statements

    print "Number of elecs: ", len(input.elecs)

    # Let's set the dielectric in the second elec statement

    input.elecs[1].sdie = 55.55

    # Now run APBS
    
    potentials = runAPBS(PQR,str(input))

    # And print the results!

    print "Now we have: ", potentials
