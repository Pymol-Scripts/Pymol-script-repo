#!/usr/bin/env python

import os
from string import find

from MolKit import Read
from MolKit.pdbWriter import PdbqsWriter
from AutoDockTools.atomTypeTools import SolvationParameterizer




if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: pdbqt_to_pdbqs.py -s filename"
        print
        print "    Description of command..."
        print "        [-s]    receptor_filename_stem"
        print "    Optional parameters:"
        print "        [-o]    alternative pdbqs_filename"
        print "        [-v]    verbose output"


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 's:o:v')

    except getopt.GetoptError, msg:
        print 'pdbqt_to_pdbqs.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    #-s: pdbqt_filename_stem
    pdbqt_filename =  None

    # optional parameters
    verbose = None
    pdbqs_filename =  None

    #'s:o:v'
    for o, a in opt_list:
        if o in ('-s', '--s'):
            pdbqt_filename = a + ".pdbqt"
            pdbqs_filename = a + ".pdbqs"
            if verbose: print 'set AD4 receptor_filename to ', pdbqt_filename
        if o in ('-o', '--o'):
            pdbqs_filename = a 
            if verbose: 
                print 'set AD3 pdbqs_filename to ', a
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print 'set verbose to ', True
        if o in ('-h', '--'):
            usage()
            sys.exit()


    if not pdbqt_filename:
        print 'pdbqt_to_pdbqs: pdbqt_filename_stem must be specified.'
        usage()
        sys.exit()

    #what about nucleic acids???

    mols = Read(pdbqt_filename)
    if verbose: print 'read ', pdbqt_filename
    mol = mols[0]
    #fix number for problem files with alternative positions
    mol.allAtoms.number = range(1, len(mol.allAtoms)+1)

    #possible clean-up???
    #need to type atoms + assign babel_types
    solparer = SolvationParameterizer()
    solparer.addParameters(mol.allAtoms)
    if verbose:
        print "set autodock3 solvation parameters for ", mol.name

    #pdbqs_filename = mol.name + '.pdbqs'
    writer = PdbqsWriter()
    fptr = open(pdbqs_filename, 'w')
    ctr = 0
    for line in mol.parser.allLines:
        if find(line, 'ATOM')<0 and find(line, "HETA")<0:
            fptr.write(line)
        else:
            this_atom = mol.allAtoms[ctr]
            writer.write_atom(fptr, this_atom)
            ctr = ctr + 1
    fptr.close()
    if verbose:
        print "wrote ", ctr, " atoms to", pdbqs_filename
    

# To execute this command type:
# pdbqt_to_pdbqs.py -s pdbqt_filename_stem  -v




