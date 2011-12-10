#!/usr/bin/env python

import os
from string import find

from MolKit import Read
from MolKit.pdbWriter import PdbWriter




if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: pdbqt_to_pdb.py -s filename"
        print
        print "    Description of command..."
        print "        [-f]    pdbqt_filename"
        print "    Optional parameters:"
        print "        [-o]    alternative pdb_filename"
        print "        [-v]    verbose output"


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'f:o:v')

    except getopt.GetoptError, msg:
        print 'pdbqt_to_pdb.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    #-s: pdbqt_filename_stem
    pdbqt_filename =  None

    # optional parameters
    verbose = None
    pdb_filename =  None

    #'f:o:v'
    for o, a in opt_list:
        if o in ('-f', '--f'):
            pdbqt_filename = a
            stem = pdbqt_filename.split('.')[0]
            pdb_filename = stem + ".pdb"
            if verbose: print 'set pdbqt_filename to ', pdbqt_filename
        if o in ('-o', '--o'):
            pdb_filename = a 
            if verbose: 
                print 'set output pdb_filename to ', a
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print 'set verbose to ', True
        if o in ('-h', '--'):
            usage()
            sys.exit()


    if not pdbqt_filename:
        print 'pdbqt_to_pdb: pdbqt_filename must be specified.'
        usage()
        sys.exit()

    mols = Read(pdbqt_filename)
    if verbose: print 'read ', pdbqt_filename
    mol = mols[0]
    mol.buildBondsByDistance()
    #fix number for problem files with alternative positions
    mol.allAtoms.number = range(1, len(mol.allAtoms)+1)

    #pdb_filename = mol.name + '.pdb'
    writer = PdbWriter()
    writer.write(pdb_filename, mol.allAtoms, records=['ATOM', 'HETATM'])

    if verbose:
        print "wrote ", ctr, " atoms to", pdb_filename
    

# To execute this command type:
# pdbqt_to_pdb.py -f pdbqt_filename_stem [-o outputfilename] -v




