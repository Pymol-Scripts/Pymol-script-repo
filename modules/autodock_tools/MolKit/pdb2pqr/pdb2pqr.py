#!/usr/bin/python
"""
    Driver for PDB2PQR

    This module takes a PDB file as input and performs optimizations
    before yielding a new PDB-style file as output.

    Ported to Python by Todd Dolinsky (todd@ccb.wustl.edu)
    Washington University in St. Louis

    Parsing utilities provided by Nathan A. Baker (baker@biochem.wustl.edu)
    Washington University in St. Louis
"""

__date__  = "17 March 2007"
__author__ = "Todd Dolinsky, Nathan Baker"
__version__ = "1.2.1"

import string
import sys
import getopt
import os
#Sargis: This is required to fix xml bug on linux with no ssl.so
#TODO: remove this when it fixed
try:
    import _hashlib
except ImportError:
    site = os.path.split(os.__file__)[0]                    
    sys.path.remove(os.path.join(site,'site-packages'))
import time
from src import pdb
from src import utilities
from src import structures
from src import routines
from src import protein
from src import server
from src.pdb import *
from src.utilities import *
from src.structures import *
from src.definitions import *
from src.forcefield import *
from src.routines import *
from src.protein import *
from src.server import *
from src.hydrogens import *
from StringIO import *

def usage(rc):
    """
        Print usage for this script to stdout.

        Parameters
            rc:  Exit status (int)
    """

    str = "\n"
    str = str + "pdb2pqr  (Version %s)\n" % __version__
    str = str + "\n"
    str = str + "This module takes a PDB file as input and performs\n"
    str = str + "optimizations before yielding a new PDB-style file as\n"
    str = str + "output\n"
    str = str + "\n"
    str = str + "Usage: pdb2pqr.py [options] --ff=<forcefield> <path> <output-path>\n"
    str = str + "    Required Arguments:\n"
    str = str + "        <forcefield>  :  The forcefield to use - currently amber\n"
    str = str + "                         charmm, parse, and tyl06  are supported.\n"
    str = str + "        <path>        :  The path to the PDB file or an ID\n"
    str = str + "                         to obtain from the PDB archive\n"
    str = str + "        <output-path> :  The desired output name of the PQR file\n"
    str = str + "                         to be generated\n"
    str = str + "    Optional Arguments:\n"
    str = str + "        --nodebump    :  Do not perform the debumping operation\n"
    str = str + "        --noopt       :  Do not perform hydrogen optimization\n"
    str = str + "        --chain       :  Keep the chain ID in the output PQR file\n" 
    str = str + "        --assign-only :  Only assign charges and radii - do not add\n"
    str = str + "                         atoms, debump, or optimize.\n"
    str = str + "        --clean       :  Do no optimization, atom addition, or\n"
    str = str + "                         parameter assignment, just return the\n"
    str = str + "                         original PDB file in aligned format.\n"
    str = str + "        --ffout=<name>:  Instead of using the standard canonical\n"
    str = str + "                         naming scheme for residue and atom names,\n"
    str = str + "                         use the names from the given forcefield.\n"
    str = str + "        --with-ph=<ph>:  Use propka to calculate pKas and apply them\n"
    str = str + "                         to the molecule given the pH value. Actual\n"
    str = str + "                         PropKa results will be output to \n"
    str = str + "                         <output-path>.propka.\n"
    str = str + "        --apbs-input  :  Create a template APBS input file based on\n"
    str = str + "                         the generated PQR file.\n"
    str = str + "        --ligand=<path>: Calculate the parameters for the ligand in\n"
    str = str + "                         mol2 format at the given path. Pdb2pka must\n"
    str = str + "                         be compiled\n"
    str = str + "        --verbose (-v):  Print information to stdout\n"
    str = str + "        --help    (-h):  Display the usage information\n"

    # Check to see if there are usage statements from the
    # extensions directory

    extensions = getAvailableExtensions()
    if len(extensions) > 0:
        str = str + "\n    Optional Arguments from Extensions Directory:\n"
        for ext in extensions:
            str += extensions[ext].usage()
    
    str = str + "\n"
    sys.stderr.write(str)
    sys.exit(rc)

def printHeader(atomlist, reslist, charge, ff, warnings, options):
    """
        Print the header for the PQR file

        Parameters:
            atomlist: A list of atoms that were unable to have
                      charges assigned (list)
            reslist:  A list of residues with non-integral charges
                      (list)
            charge:   The total charge on the protein (float)
            ff:       The forcefield name (string)
            warnings: A list of warnings generated from routines (list)
            options:  A dictionary of command lnie options (float)
        Returns
            header:   The header for the PQR file (string)
    """
    header = "REMARK   1 PQR file generated by PDB2PQR (Version %s)\n" % __version__
    header = header + "REMARK   1\n"
    header = header + "REMARK   1 Forcefield Used: %s\n" % ff
    if "ffout" in options:
        header = header + "REMARK   1 Naming Scheme Used: %s\n" % options["ffout"]
    header = header + "REMARK   1\n"
    
    if "ph" in options:
        header = header + "REMARK   1 pKas calculated by propka and assigned using pH %.2f\n" % options["ph"]
        header = header + "REMARK   1\n"

    for warning in warnings:
        header = header + "REMARK   5 " + warning 
    header = header + "REMARK   5\n"
    
    if len(atomlist) != 0:
        header += "REMARK   5 WARNING: PDB2PQR was unable to assign charges\n"
        header += "REMARK   5          to the following atoms (omitted below):\n"
        for atom in atomlist:
            header += "REMARK   5              %i %s in %s %i\n" % \
                      (atom.get("serial"), atom.get("name"), \
                       atom.get("residue").get("name"), \
                       atom.get("residue").get("resSeq"))
        header += "REMARK   5\n"
    if len(reslist) != 0:
        header += "REMARK   5 WARNING: Non-integral net charges were found in\n"
        header += "REMARK   5          the following residues:\n"
        for residue in reslist:
            header += "REMARK   5              %s - Residue Charge: %.4f\n" % \
                      (residue, residue.getCharge())
        header += "REMARK   5\n"
    header += "REMARK   6 Total charge on this protein: %.4f e\n" % charge
    header += "REMARK   6\n"

    return header
  
def runPDB2PQR(pdblist, ff, options):
    """
        Run the PDB2PQR Suite

        Parameters
            pdblist: The list of objects that was read from the PDB file
                     given as input (list)
            ff:      The name of the forcefield (string)
            options: A dictionary of PDB2PQR options, including:
                     verbose: When 1, script will print information to stdout
                              When 0, no detailed information will be printed (int)
                     debump:  When 1, debump heavy atoms (int)
                     opt:     When 1, run hydrogen optimization (int)
                     ph:      The desired ph of the system (float)
                     outname: The name of the desired output file
        Returns
            header:  The PQR file header (string)
            lines:   The PQR file atoms (list)
            missedligandresidues:  A list of ligand residue names whose charges could
                     not be assigned (ligand)
    """
    ph = None
    pkaname = ""
    outname = ""
    outroot = ""
    typemapname = ""
    lines = []

    # userff is CGI-based User Forcefield file object

    if "userff" in options: userff = options["userff"]
    else: userff = None
    
    if "verbose" in options: verbose = 1
    else: verbose = 0

    if "opt" in options: optflag = 1
    else: optflag = 0

    if "chain" in options: chainflag = 1
    else: chainflag = 0

    if "outname" not in options or options["outname"] == None:
        text = "Error: Output name not set!"
        raise ValueError, text
    else:
        outname = options["outname"]
        period = string.find(outname,".")
        if period > 0: outroot = outname[0:period]
        else: outroot = outname

    if "ph" in options:
        pka = 1
        ph = options["ph"]
        pkaname = outroot + ".propka"
        if os.path.isfile(pkaname): os.remove(pkaname)
    else: pka = 0

    typemapname = "%s-typemap.html" % outroot

    extmap = options["extensions"]
    
    start = time.time()

    if verbose:
        print "Beginning PDB2PQR...\n"

    myDefinition = Definition()
    if verbose:
        print "Parsed Amino Acid definition file."   

    # Check for the presence of a ligand!  This code is taken from pdb2pka/pka.py

    if "ligand" in options:
        from pdb2pka.ligandclean import ligff
        myProtein, myDefinition, Lig = ligff.initialize(myDefinition, options["ligand"], pdblist, verbose)        
    else:
        myProtein = Protein(pdblist, myDefinition)

    if verbose:
        print "Created protein object -"
        print "\tNumber of residues in protein: %s" % myProtein.numResidues()
        print "\tNumber of atoms in protein   : %s" % myProtein.numAtoms()
        
    myRoutines = Routines(myProtein, verbose)

    myRoutines.setTermini()
    myRoutines.updateBonds()

    if "clean" in options:
        header = ""
        lines = myProtein.printAtoms(myProtein.getAtoms(), chainflag)
      
        # Process the extensions
        for ext in extmap:
            module = extmap[ext]
            call = "module.%s(myRoutines, outroot)" % ext
            eval(call)  
    
        if verbose:
            print "Total time taken: %.2f seconds\n" % (time.time() - start)
        return header, lines

    if not "assign-only" in options:

        myRoutines.findMissingHeavy()
        myRoutines.updateSSbridges()

        if "debump" in options:
            myRoutines.debumpProtein()  

        if pka:
            myRoutines.runPROPKA(ph, ff, pkaname)

        myRoutines.addHydrogens()

        if optflag:
            myhydRoutines = hydrogenRoutines(myRoutines)
            myhydRoutines.setOptimizeableHydrogens()

        if "debump" in options:
            myRoutines.debumpProtein()  

        if optflag:
            myhydRoutines.initializeFullOptimization()
            myhydRoutines.optimizeHydrogens()
        else:
            myhydRoutines = hydrogenRoutines(myRoutines)
            myhydRoutines.initializeWaterOptimization()
            myhydRoutines.optimizeHydrogens()

    else:  # Special case for HIS if using assign-only
        for residue in myProtein.getResidues():
            if isinstance(residue, HIS):
                myRoutines.applyPatch("HIP", residue)

    myRoutines.setStates()

    myForcefield = Forcefield(ff, myDefinition, userff)
    hitlist, misslist = myRoutines.applyForcefield(myForcefield)
  
    ligsuccess = 0
    if "ligand" in options:

        # If this is independent, we can assign charges and radii here
 
        for residue in myProtein.getResidues():
            if isinstance(residue, LIG):
                templist = []
                Lig.make_up2date(residue)
                for atom in residue.getAtoms():
                    atom.ffcharge = Lig.ligand_props[atom.name]["charge"]
                    atom.radius = Lig.ligand_props[atom.name]["radius"]
                    if atom in misslist:
                        misslist.pop(misslist.index(atom))
                        templist.append(atom)

                charge = residue.getCharge()
                if abs(charge - int(charge)) > 0.001:
                    # Ligand parameterization failed
                    myRoutines.warnings.append("WARNING: PDB2PQR could not successfully parameterize\n")
                    myRoutines.warnings.append("         the desired ligand; it has been left out of\n")
                    myRoutines.warnings.append("         the PQR file.\n")
                    myRoutines.warnings.append("\n")
                    
                    # remove the ligand
                    myProtein.residues.remove(residue) 
                    for chain in myProtein.chains:
                        if residue in chain.residues: chain.residues.remove(residue)
                else:
                    ligsuccess = 1
                    # Mark these atoms as hits
                    hitlist = hitlist + templist
    
    # Temporary fix; if ligand was successful, pull all ligands from misslist
    if ligsuccess:
        templist = misslist[:]
        for atom in templist:
            if isinstance(atom.residue, Amino) or isinstance(atom.residue, Nucleic): continue
            misslist.remove(atom)

    # Creat the Typemap
    myProtein.createHTMLTypeMap(myDefinition, typemapname)

    # Grab the protein charge

    reslist, charge = myProtein.getCharge()

    # If we want a different naming scheme, use that

    if "ffout" in options:
        scheme = options["ffout"]
        userff = None # Currently not supported
        if scheme != ff: myNameScheme = Forcefield(scheme, myDefinition, userff)
        else: myNameScheme = myForcefield
        myRoutines.applyNameScheme(myNameScheme)

    header = printHeader(misslist, reslist, charge, ff, myRoutines.getWarnings(), options)
    lines = myProtein.printAtoms(hitlist, chainflag)

    # Determine if any of the atoms in misslist were ligands
    missedligandresidues = []
    for atom in misslist:
        if isinstance(atom.residue, Amino) or isinstance(atom.residue, Nucleic): continue
        if atom.resName not in missedligandresidues:
            missedligandresidues.append(atom.resName)

    # Process the extensions
 
    for ext in extmap:
        module = extmap[ext]
        call = "module.%s(myRoutines, outroot)" % ext
        eval(call)

    if verbose:
        print "Total time taken: %.2f seconds\n" % (time.time() - start)

    return header, lines, missedligandresidues

def getAvailableExtensions(displayflag=0):
    """
        Grab available extensions from the extensions directory

        Parameters
            displayflag: Display the error message if 1
        Returns
            extensions: A map containing the extensions name and
                        the module instance.
    """
    extensions = {}
    dir = "%s" % os.path.dirname(sys.argv[0])
    if dir == "": extdir = "extensions"
    else: extdir = "%s/extensions" % dir
    for filename in os.listdir(extdir):
        if filename.endswith(".py"):
            if filename == "__init__.py": continue
            
            # Test to see if we can find the function

            name = filename[:-3]
            try:
                e = __import__("extensions.%s" % name, globals(), locals(), name)
                if callable(eval("e.%s" % name)) and \
                   callable(eval("e.usage")):
                    extensions[name] = e
            except (AttributeError, ImportError):
                txt = "\nWarning: Missing either \"%s\" or \"usage\" functions in %s!" %\
                      (name, filename)
                txt += "\nThis extension will not be included.\n\n"
                if displayflag:
                    sys.stderr.write(txt)

    return extensions

def mainCommand():
    """
        Main driver for running program from the command line.
    """
    shortOptlist = "h,v"
    longOptlist = ["help","verbose","ff=","ffout=","nodebump","noopt","with-ph=","apbs-input","chain","clean","assign-only", "ligand="]

    extensions = getAvailableExtensions(1)
    longOptlist += extensions.keys()

    try: opts, args = getopt.getopt(sys.argv[1:], shortOptlist, longOptlist)
    except getopt.GetoptError, details:
        sys.stderr.write("GetoptError:  %s\n" % details)
        usage(2)

    if len(args) != 2:
        sys.stderr.write("Incorrect number (%d) of arguments!\n" % len(args))
        usage(2)

    options = {"debump":1,"opt":1,"extensions":{}}
 
    outpath = None
    ff = None
    for o,a in opts:
        undashed = o[2:]
        if o in ("-v","--verbose"):
            options["verbose"] = 1
        elif o in ("-h","--help"):
            usage(2)
            sys.exit()
        elif o == "--nodebump":  del options["debump"]
        elif o == "--noopt":    del options["opt"]
        elif o == "--apbs-input": options["input"] = 1
        elif o == "--with-ph":
            try:
                ph = float(a)
                options["ph"] = ph
                if ph < 0.0 or ph > 14.0: raise ValueError
            except ValueError:
                text = "%s is not a valid pH!  " % a
                text += "Please choose a pH between 0.0 and 14.0."
                raise ValueError, text
        elif o == "--assign-only":
            del options["debump"]
            del options["opt"]
            options["assign-only"] = 1
        elif o == "--clean":
            del options["debump"]
            del options["opt"]
            options["clean"] = 1
        elif o == "--ff":      
            ff = a
            
            # Check to make sure forcefield file is available

            defpath = getFFfile(ff)
            if defpath == "":
                raise ValueError, "Unable to find parameter files for forcefield %s!" % ff
            
        elif o == "--chain": options["chain"] = 1
        elif o == "--ffout":
            if a in ["amber","AMBER","charmm","CHARMM","parse","PARSE","tyl06","TYL06"]:
                options["ffout"] = a
            else:
                raise ValueError, "Invalid forcefield naming scheme %s!" % a
        elif o == "--ligand":
            if os.path.isfile(a):
                options["ligand"] = open(a)
            else:
                raise ValueError, "Unable to find ligand file %s!\n" % a
        elif undashed in extensions.keys():
            options["extensions"][undashed] = extensions[undashed]
            
    if ff == None and "clean" not in options:
        raise ValueError, "Forcefield not specified!"

    text =  "\n--------------------------\n"
    text += "PDB2PQR - a Python-based structural conversion utility\n"
    text += "--------------------------\n"
    text += "Please cite your use of PDB2PQR as:\n"
    text += "  Dolinsky TJ, Nielsen JE, McCammon JA, Baker NA.\n"
    text += "  PDB2PQR: an automated pipeline for the setup, execution,\n"
    text += "  and analysis of Poisson-Boltzmann electrostatics calculations.\n"
    text += "  Nucleic Acids Research 32 W665-W667 (2004).\n\n"
    sys.stderr.write(text)
            
    path = args[0]
    file = getPDBFile(path)
    pdblist, errlist = readPDB(file)
    
    if len(pdblist) == 0 and len(errlist) == 0:
        try: os.remove(path)
        except OSError: pass
        raise ValueError, "Unable to find file %s!\n" % path

    if len(errlist) != 0 and "verbose" in options:
        print "Warning: %s is a non-standard PDB file.\n" % path
        print errlist

    outpath = args[1]
    options["outname"] = outpath

    header, lines, missedligands = runPDB2PQR(pdblist, ff, options)

    # Print the PQR file
    outfile = open(outpath,"w")
    outfile.write(header)
    for line in lines:
        outfile.write(line)
    outfile.close()

    if "input" in options:
        from src import inputgen
        from src import psize
        method = "mg-auto"
        size = psize.Psize()
        size.parseInput(outpath)
        size.runPsize(outpath)
        async = 0 # No async files here!
        input = inputgen.Input(outpath, size, method, async)
        input.printInputFiles()
   
def mainCGI():
    """
        Main driver for running PDB2PQR from a web page
    """
    import cgi
    import cgitb

    cgitb.enable()
    form = cgi.FieldStorage()

    options = {"extensions":{}}
 
    ff = form["FF"].value 
    input = 0
  
    if form.has_key("DEBUMP"): options["debump"] = 1
    if form.has_key("OPT"): options["opt"] = 1
    if form.has_key("PROPKA"):
        try:
            ph = float(form["PH"].value)
            if ph < 0.0 or ph > 14.0: raise ValueError
            options["ph"] = ph
        except ValueError:
             text = "The entered pH of %.2f is invalid!  " % form["PH"].value
             text += "Please choose a pH between 0.0 and 14.0."
             print "Content-type: text/html\n"
             print text
             sys.exit(2)
    if form.has_key("PDBID"):
        file = getPDBFile(form["PDBID"].value)
    elif form.has_key("PDB"):
        file = StringIO(form["PDB"].value)
    if form.has_key("INPUT"):
        input = 1
        options["apbs"] = 1
    if form.has_key("USERFF"):
        userff = StringIO(form["USERFF"].value)
        ff = "user-defined"
        options["userff"] = userff
    if form.has_key("FFOUT"):
        if form["FFOUT"].value != "internal":
            options["ffout"] = form["FFOUT"].value
    if form.has_key("CHAIN"):
        options["chain"] = 1
    if form.has_key("LIGAND"):
        options["ligand"] = StringIO(form["LIGAND"].value)    
  
    pdblist, errlist = readPDB(file)    
    if len(pdblist) == 0 and len(errlist) == 0:
        text = "Unable to find PDB file - Please make sure this is "
        text += "a valid PDB file ID!"
        print "Content-type: text/html\n"
        print text
        sys.exit(2)
    elif len(pdblist) > 10000 and "opt" in options:
        text = "<HTML><HEAD>"
        text += "<TITLE>PDB2PQR Error</title>"
        text += "<link rel=\"stylesheet\" href=\"%s\" type=\"text/css\">" % STYLESHEET
        text += "</HEAD><BODY><H2>PDB2PQR Error</H2><P>"
        text += "Due to server limits, we are currently unable to optimize "
        text += "proteins of greater than 10000 atoms on the server.  If you "
        text += "want to forgo optimization please try the server again.<P>"
        text += "Otherwise you may use the standalone version of PDB2PQR that "
        text += "is available from the <a href=\"http://pdb2pqr.sourceforge.net\">"
        text += "PDB2PQR SourceForge project page</a>."
        text += "</BODY></HTML>"
        print "Content-type: text/html\n"
        print text
        sys.exit(2)
        
    try:
        starttime = time.time()
        name = setID(starttime)
 
        pqrpath = startServer(name)
        options["outname"] = pqrpath
        header, lines, missedligands = runPDB2PQR(pdblist, ff, options)
        file = open(pqrpath, "w")
        file.write(header)
        for line in lines:
            file.write("%s\n" % string.strip(line))
        file.close()
                
        if input:
            from src import inputgen
            from src import psize
            method = "mg-auto"
            size = psize.Psize()
            size.parseInput(pqrpath)
            size.runPsize(pqrpath)
            async = 0 # No async files here!
            myinput = inputgen.Input(pqrpath, size, method, async)
            myinput.printInputFiles()
                    
        endtime = time.time() - starttime
        createResults(header, input, name, endtime, missedligands)
        logRun(options, endtime, len(lines), ff, os.environ["REMOTE_ADDR"])

    except StandardError, details:
        print "Content-type: text/html\n"
        print details
        createError(name, details)
    
if __name__ == "__main__":
    """ Determine if called from command line or CGI """
    
    if not os.environ.has_key("REQUEST_METHOD"):
        mainCommand()    
    else:
        mainCGI()
