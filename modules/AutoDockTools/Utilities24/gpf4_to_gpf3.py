#!/usr/bin/env python

import os
from string import find, strip

from MolKit import Read
from MolKit.pdbWriter import PdbqtWriter
from AutoDockTools.atomTypeTools import AutoDock4_AtomTyper
from AutoDockTools.GridParameters import GridParameters, grid_parameter_list




if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: gpf4_to_gpf3.py -s gpf4_stem -r receptor_filename -l ligand_filename"
        print
        print "    Description of command..."
        print "        [-s]    gpf4_stem"
        print "        [-r]    receptor_filename"
        print "        [-l]    ligand_filename"
        print "    Optional parameters:"
        print "        [-o]    outputfilename"
        print "        [-v]    verbose output"


    # process command arguments
    try:
        #opt_list, args = getopt.getopt(sys.argv[1:], 'r:o:bs:ampvHG')
        opt_list, args = getopt.getopt(sys.argv[1:], 's:r:l:o:v')

    except getopt.GetoptError, msg:
        print 'gpf4_to_gpf3.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    #-s: pdbq_stem
    gpf4_filename =  None

    #-r: receptor_filename
    receptor_filename =  None

    #-l: ligand_filename
    ligand_filename =  None

    # optional parameters
    # optional parameters
    verbose = None
    outputfilename = None

    #'s:r:l:v'
    for o, a in opt_list:
        if o in ('-s', '--s'):
            gpf4_stem = a
            gpf4_filename = a + '.gpf'
            if verbose: print 'set gpf4_filename to ', gpf4_filename
        if o in ('-r', '--r'):
            receptor_filename = a
            if verbose: print 'set receptor_filename to ', receptor_filename
        if o in ('-l', '--l'):
            ligand_filename = a
            if verbose: print 'set ligand_filename to ', ligand_filename
        if o in ('-o', '--o'):
            outputfilename = a
            if verbose: print 'set outputfilename to ', outputfilename
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print 'set verbose to ', True
        if o in ('-h', '--'):
            usage()
            sys.exit()


    if not ligand_filename:
        print 'gpf4_to_gpf3: ligand_filename must be specified.'
        usage()
        sys.exit()
    if not receptor_filename:
        receptor_filename = gpf4_stem + '.pdbqs'
        print "gpf4_to_gpf3: using default receptor_filename:", receptor_filename
    if not gpf4_filename:
        print 'gpf4_to_gpf3: no gpf4_filename specified using GridParameter defaults'

    #what about nucleic acids???
    GPO3 = GridParameters()
    GPO3.read4(gpf4_filename)
    GPO3.read(gpf4_filename)
    if verbose: print 'read ', gpf4_filename
    #have to set the ligand_types
    GPO3.set_ligand(ligand_filename)
    #build an autogrid3 'types' string from autogrid4 ligand_types
    GPO3.set_ligand_types3(GPO3['ligand_types']['value'])
    if verbose: print 'set ligand to  ', ligand_filename
    GPO3.set_receptor(receptor_filename)
    if verbose: print 'set receptor to  ', receptor_filename
    if outputfilename is None:
        outputfilename = gpf4_stem + '_3.gpf'
    GPO3.write(outputfilename, grid_parameter_list)
    if verbose: 
        print "wrote ", outputfilename, ' using:'
        for p in grid_parameter_list:
            print p
    

# To execute this command type:
# gpf4_to_gpf3.py -s gpf4_stem -v




