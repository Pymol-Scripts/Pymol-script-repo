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
try:
    from apbslib import *
except:
    #
    # We need apbslib.so, apbslib.so and apbslib.py
    #
    print
    print 'Missing libraries for interfacing with APBS'
    print
    print 'You need to find _apbslib.so and apbslib.py and symlink into the pdb2pqr/pdb2pka source code directory'
    print 'The files can be found in the share/tools/python dir of your apbs installation'
    print
    sys.exit(0)

Python_kb = 1.3806581e-23
Python_Na = 6.0221367e+23
NOSH_MAXMOL = 20
NOSH_MAXCALC = 20
ACD_ERROR = 2    # < Error setting up calculation>

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

class runAPBS:

    def __init__(self):
        return
    

    def getUnitConversion(self):
        """
            Get the unit conversion from kT to kJ/mol

            Returns
                factor: The conversion factor (float)
        """
        temp = 298.15
        factor = Python_kb/1000.0 * temp * Python_Na
        return factor

    def runAPBS(self,protein, inputpath, CM=None):
        """
            Run APBS, using the protein instead of a pqr file

            Parameters
                protein:    The protein object (protein)
                inputpath:  The path to the APBS input file (string)
            Returns
                potentials: A list of lists of potentials at atom
                            locations - one list for each APBS
                            calculation
        """
        #
        # Initialize the MALOC library
        startVio()

        # Initialize variables, arrays
        self.com = Vcom_ctor(1)
        self.rank = Vcom_rank(self.com)
        self.size = Vcom_size(self.com)
        self.mgparm = MGparm()
        self.pbeparm = PBEparm()
        self.mem = Vmem_ctor("Main")
        self.pbe = new_pbelist(NOSH_MAXMOL)
        self.pmg = new_pmglist(NOSH_MAXMOL)
        self.pmgp = new_pmgplist(NOSH_MAXMOL)
        self.realCenter = double_array(3)
        self.totEnergy = []
        self.x = []
        self.y = []
        self.z = []
        self.chg = []
        self.rad = []
        #nforce = int_array(NOSH_MAXCALC)
        #atomforce = new_atomforcelist(NOSH_MAXCALC)
        #nfor = ptrcreate("int",0)

        # Start the main timer
        self.main_timer_start = time.clock()

        # Check invocation
        #stdout.write(getHeader())

        # Parse the input file
        self.nosh = NOsh_ctor(self.rank, self.size)
        #nosh = NOsh()
        #NOsh_ctor2(nosh, rank, size)
        sys.stdout.write("Parsing input file %s...\n" % inputpath)
        if NOsh_parseInputFile(self.nosh, inputpath) != 1:
            sys.stderr.write("main:  Error while parsing input file.\n")
            raise APBSError, "Error while parsing input file!"

        # Load the molecules using Valist_load routine

        self.alist = new_valist(NOSH_MAXMOL)
        self.atoms = protein.getAtoms()
        self.protsize = len(self.atoms)

        # SETUP CALCULATIONS

        if NOsh_setupElecCalc(self.nosh, self.alist) != 1:
            sys.stderr.write("Error setting up ELEC calculations\n")
            raise APBSError, "Error while setting up calculations!"

        if NOsh_setupApolCalc(self.nosh, self.alist) == ACD_ERROR:
            sys.stderr.write("Error setting up APOL calculations\n")
            raise APBSError, "Error while setting up calculations!"

        #
        # DEBUGGING
        #
        self.cm_list=[]

        #print 'These are the charges in the PQR file'
        #print 'atom name\tresnum\tresname\tcharge\tradius'
        #res_charge=0.0
        #old_res=-1
        #old_res_name=''
        for i in range(len(self.atoms)):
            atom = self.atoms[i]
            self.x.append(atom.get("x"))
            self.y.append(atom.get("y"))
            self.z.append(atom.get("z"))
            self.chg.append(atom.get("ffcharge"))
            self.rad.append(atom.get("radius"))
            self.cm_list.append([atom.resSeq,atom.name,atom.get("ffcharge")])
            #if atom.resSeq!=old_res:
                #if old_res!=-1:
                #    print '%4i %4s %6.4f' %(old_res,old_res_name,res_charge)
            #    res_charge=atom.get('ffcharge')
            #    old_res=atom.resSeq
            #    old_res_name=atom.resName
            #else:
            #    res_charge=res_charge+atom.get('ffcharge')
                
            #print '%5s\t%4i\t%5s\t%6.4f\t%6.4f' %(atom.name,atom.resSeq,atom.resName,atom.get("ffcharge"),atom.get('radius'))
        #print '%4i %4s %6.4f' %(old_res,old_res_name,res_charge)
        #
        # DEBUG
        #
        if CM:
            CM.display_charges(self.cm_list)
        #

        self.myAlist = make_Valist(self.alist,0)
        Valist_load(self.myAlist, self.protsize, self.x,self.y,self.z,self.chg,self.rad) 

        # Initialize the energy holders
        self.totEnergy = [ 0.0 ] * int( self.nosh.ncalc )
        
        potList = []

        # Initialize the force holders
        forceList = []

        # Load the dieletric maps

        self.dielXMap = new_gridlist(NOSH_MAXMOL)
        self.dielYMap = new_gridlist(NOSH_MAXMOL)
        self.dielZMap = new_gridlist(NOSH_MAXMOL)

        if loadDielMaps(self.nosh, self.dielXMap, self.dielYMap, self.dielZMap) != 1:
            sys.stderr.write("Error reading dielectric maps!\n")
            raise APBSError, "Error reading dielectric maps!"

        # Load the kappa maps
        self.kappaMap = new_gridlist(NOSH_MAXMOL)
        if loadKappaMaps(self.nosh, self.kappaMap) != 1:
            sys.stderr.write("Error reading kappa maps!\n")
            raise APBSError, "Error reading kappa maps!"

        # Load the potential maps
        self.potMap = new_gridlist(NOSH_MAXMOL)
        if loadPotMaps(self.nosh, self.potMap) != 1:
            sys.stderr.write("Error reading potential maps!\n")
            raise APBSError, "Error reading potential maps!"

        # Load the charge maps
        self.chargeMap = new_gridlist(NOSH_MAXMOL)
        if loadChargeMaps(self.nosh, self.chargeMap) != 1:
            sys.stderr.write("Error reading charge maps!\n")
            raise APBSError, "Error reading charge maps!"

        # Do the calculations

        sys.stdout.write("Preparing to run %d PBE calculations. \n" % self.nosh.ncalc)

        for icalc in xrange(self.nosh.ncalc):
            sys.stdout.write("---------------------------------------------\n")
            self.calc = NOsh_getCalc(self.nosh, icalc)
            self.mgparm = self.calc.mgparm
            self.pbeparm = self.calc.pbeparm
            if self.calc.calctype != 0:
                sys.stderr.write("main:  Only multigrid calculations supported!\n")
                raise APBSError, "Only multigrid calculations supported!"

            for k in range(0, self.nosh.nelec):
                if NOsh_elec2calc(self.nosh,k) >= icalc:
                    break

            name = NOsh_elecname(self.nosh, k+1)
            #if name == "":
            #    sys.stdout.write("CALCULATION #%d:  MULTIGRID\n" % (icalc+1))
            #else:
            #    sys.stdout.write("CALCULATION #%d (%s): MULTIGRID\n" % ((icalc+1),name))
            #sys.stdout.write("Setting up problem...\n")

            # Routine initMG

            if initMG(icalc, self.nosh, self.mgparm, self.pbeparm, self.realCenter, self.pbe, 
                  self.alist, self.dielXMap, self.dielYMap, self.dielZMap, self.kappaMap, self.chargeMap, 
                  self.pmgp, self.pmg, self.potMap) != 1:
                sys.stderr.write("Error setting up MG calculation!\n")
                raise APBSError, "Error setting up MG calculation!"

            # Print problem parameters 

            #printMGPARM(self.mgparm, self.realCenter)
            #printPBEPARM(self.pbeparm)

            # Solve the problem : Routine solveMG

            self.thispmg = get_Vpmg(self.pmg,icalc)

            if solveMG(self.nosh, self.thispmg, self.mgparm.type) != 1:
                stderr.write("Error solving PDE! \n")
                raise APBSError, "Error Solving PDE!"

            # Set partition information : Routine setPartMG

            if setPartMG(self.nosh, self.mgparm, self.thispmg) != 1:
                sys.stderr.write("Error setting partition info!\n")
                raise APBSError, "Error setting partition info!"

            ret, self.totEnergy[icalc] = energyMG(self.nosh, icalc, self.thispmg, 0,
                                                  self.totEnergy[icalc], 0.0, 0.0, 0.0)

            # Set partition information

            #aforce = get_AtomForce(atomforce, icalc)
            #forceMG(mem, nosh, pbeparm, mgparm, thispmg, nfor, aforce, alist)
            #ptrset(nforce,ptrvalue(nfor), icalc)

            # Write out data from MG calculations : Routine writedataMG	
            writedataMG(self.rank, self.nosh, self.pbeparm, self.thispmg)

            # Write out matrix from MG calculations	
            writematMG(self.rank, self.nosh, self.pbeparm, self.thispmg)

            # GET THE POTENTIALS

            potentials = getPotentials(self.nosh, self.pbeparm, self.thispmg, self.myAlist)
            potList.append(potentials)

        #
        # Cleanup
        #
        return potList

    #
    # ------
    #

    def get_potentials(self,protein):

        import sys
        import copy

        sys.setrecursionlimit(10000)

        delete_valist(self.alist)
        self.alist = new_valist(NOSH_MAXMOL)
        self.atoms = protein.getAtoms()
        self.protsize = len(self.atoms)
        proteincopy = copy.copy(protein)

        for i in range(len(self.atoms)):
            atom = self.atoms[i]
            self.x.append(atom.get("x"))
            self.y.append(atom.get("y"))
            self.z.append(atom.get("z"))
            self.chg.append(atom.get("ffcharge"))
            self.rad.append(atom.get("radius"))
            #self.cm_list.append([atom.resSeq,atom.name,atom.get("ffcharge")])
        #
        # DEBUG
        #
        #if CM:
        #    CM.display_charges(self.cm_list)
        #
        self.myAlist = make_Valist(self.alist,0)

        xlist = self.x[-1*(self.protsize):]
        ylist = self.y[-1*(self.protsize):] 
        zlist = self.z[-1*(self.protsize):]
        chglist = self.chg[-1*(self.protsize):]
        radlist = self.rad[-1*(self.protsize):]

        Valist_load(self.myAlist, self.protsize, xlist, ylist, zlist, chglist, radlist)

        potentials=[]
        potentials = getPotentials(self.nosh,self.pbeparm,self.thispmg,self.myAlist)

        protein = copy.copy(proteincopy)

        # Free up the memory allocated for self.myAlist
        remove_Valist(self.myAlist)

        return potentials


    #
    # ------
    #
        
    def cleanup(self):

        # Handle print statements

        #if self.nosh.nprint > 0:
        #    sys.stdout.write("---------------------------------------------\n")
        #    sys.stdout.write("PRINT STATEMENTS\n")
        #for iprint in xrange(self.nosh.nprint):
        #    if NOsh_printWhat(self.nosh, iprint) == NPT_ENERGY:
        #        printEnergy(self.com,self.nosh, self.totEnergy, iprint)
        #    elif NOsh_printWhat(self.nosh, iprint) == NPT_FORCE:
        #        printForce(self.com, self.nosh, self.nforce, self.atomforce, self.iprint)
        #    else:
        #        sys.stdout.write("Undefined PRINT keyword!\n")
        #        break

        #sys.stdout.write("----------------------------------------\n")
        #sys.stdout.write("CLEANING UP AND SHUTTING DOWN...\n")

        # Clean up APBS structures

        #killForce(mem, nosh, nforce, atomforce)
        killEnergy()
        killMG(self.nosh, self.pbe, self.pmgp, self.pmg)
        killChargeMaps(self.nosh, self.chargeMap)
        killKappaMaps(self.nosh, self.kappaMap)
        killPotMaps(self.nosh, self.potMap)
        killDielMaps(self.nosh, self.dielXMap, self.dielYMap, self.dielZMap)

        if self.myAlist.number == 0:
            self.myAlist = make_Valist(self.alist,0)
            Valist_load(self.myAlist, self.protsize, self.x,self.y,self.z,self.chg,self.rad) 
        killMolecules(self.nosh, self.alist)


        del self.nosh

        # Clean up Python structures

        #ptrfree(nfor)
        delete_double_array(self.realCenter)
        #delete_int_array(nforce)
        #delete_atomforcelist(atomforce)
        delete_valist(self.alist)
        delete_gridlist(self.dielXMap)
        delete_gridlist(self.dielYMap)
        delete_gridlist(self.dielZMap)
        delete_gridlist(self.kappaMap)
        delete_gridlist(self.potMap)
        delete_gridlist(self.chargeMap)
        delete_pmglist(self.pmg)
        delete_pmgplist(self.pmgp)
        delete_pbelist(self.pbe)

        # Clean up MALOC structures
        del self.com
        del self.mem

        #sys.stdout.write("\n")
        #sys.stdout.write("Thanks for using APBS!\n\n")

        # Stop the main timer
        #main_timer_stop = time.clock()
        #sys.stdout.write("Total execution time:  %1.6e sec\n" % (main_timer_stop - self.main_timer_start))

        #Return

        return 
