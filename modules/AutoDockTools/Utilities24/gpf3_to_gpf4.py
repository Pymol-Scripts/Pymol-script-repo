#!/usr/bin/env python

import os
from string import find, strip

from MolKit import Read
from MolKit.pdbWriter import PdbqtWriter
from AutoDockTools.atomTypeTools import AutoDock4_AtomTyper
from AutoDockTools.GridParameters import GridParameters, grid_parameter_list4




if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: gpf3_to_gpf4.py -s gpf3_stem -r receptor_filename -l ligand_filename"
        print
        print "    Description of command..."
        print "        [-s]    gpf3_stem"
        print "        [-r]    receptor_filename (pdbqt)"
        print "        [-l]    ligand_filename (pdbqt)"
        print "    Optional parameters:"
        print "        [-o]    outputfilename"
        print "        [-v]    verbose output"


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 's:r:l:o:v')

    except getopt.GetoptError, msg:
        print 'gpf3_to_gpf4.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    #-s: pdbqt_stem
    gpf3_filename =  None

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
            gpf3_stem = a
            gpf3_filename = a + '.gpf'
            if verbose: print 'set gpf3_filename to ', gpf3_filename
        if o in ('-r', '--r'):
            receptor_filename = a
            if verbose: print 'set receptor_filename to ', receptor_filename
            extension = os.path.splitext(receptor_filename)[1]
            if extension!='.pdbqt':
                print receptor_filename, " not in pdbqt format"
                usage()
                sys.exit()
        if o in ('-l', '--l'):
            ligand_filename = a
            #check that the ligand_filename type is pdbqt
            extension = os.path.splitext(ligand_filename)[1]
            if extension!='.pdbqt':
                print ligand_filename, " not in pdbqt format"
                usage()
                sys.exit()
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
        print 'gpf3_to_gpf4: ligand_filename must be specified.'
        usage()
        sys.exit()
    if not receptor_filename:
        receptor_filename = gpf3_stem + '.pdbqt'
        print "gpf3_to_gpf4: using default receptor_filename:", receptor_filename
    if not gpf3_filename:
        print 'gpf3_to_gpf4: no gpf3_filename specified using GridParameter defaults'

    #what about nucleic acids???
    GPO4 = GridParameters()
    if gpf3_filename:
        GPO4.read(gpf3_filename)
        if verbose: print 'read ', gpf3_filename
    else:
        if verbose: print 'using gpo defaults'
    GPO4.set_ligand4(ligand_filename)
    if verbose: print 'set ligand to  ', ligand_filename
    GPO4.set_receptor4(receptor_filename)
    if verbose: print 'set receptor to  ', receptor_filename
    if outputfilename is None:
        outputfilename = gpf3_stem + '_4.gpf'
    GPO4.write4(outputfilename, grid_parameter_list4)
    if verbose: 
        print "wrote ", outputfilename, ' using:'
        for p in grid_parameter_list4:
            print p
    

# To execute this command type:
# gpf3_to_gpf4.py -s gpf3_stem -v




