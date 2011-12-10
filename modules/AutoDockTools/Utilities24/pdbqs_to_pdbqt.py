#!/usr/bin/env python

import os
from string import find

from MolKit import Read
from MolKit.pdbWriter import PdbqtWriter
from MolKit.chargeCalculator import GasteigerChargeCalculator
from AutoDockTools.atomTypeTools import AutoDock4_AtomTyper




if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: pdbqs_to_pdbqt.py -s filename"
        print
        print "    Description of command..."
        print "        [-s]    receptor_filename_stem"
        print "    Optional parameters:"
        print "        [-C]    preserve input charges, (default is addition of gasteiger charges)"   
        print "        [-p]    preserve input charges on atom type, eg -p Zn"
        print "        [-o]    alternative pdbqt_filename"
        print "        [-v]    verbose output"


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 's:Cp:o:v')

    except getopt.GetoptError, msg:
        print 'pdbqs_to_pdbqt.py: %s' %msg
        usage()
        sys.exit(2)

    charges_to_add = 'gasteiger'
    # initialize required parameters
    #-s: pdbqs_filename_stem
    pdbqs_filename =  None

    # optional parameters
    verbose = None
    pdbqt_filename =  None
    #-C do not add charges
    #-p preserve charges on specific atom types
    preserve_charge_types=''
    #'s:Cp:o:v'
    for o, a in opt_list:
        if o in ('-s', '--s'):
            pdbqs_filename = a + ".pdbqs"
            pdbqt_filename = a + ".pdbqt"
            if verbose: print 'set receptor_filename to ', pdbqs_filename
        if o in ('-C', '--C'):
            charges_to_add = None
            if verbose: print 'do not add charges'
        if o in ('-p', '--p'):
            preserve_charge_types+=a
            preserve_charge_types+=','
            if verbose: print 'preserve initial charges on ', preserve_charge_types
        if o in ('-o', '--o'):
            pdbqt_filename = a 
            if verbose: 
                print 'set pdbqt_filename to ', a
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print 'set verbose to ', True
        if o in ('-h', '--'):
            usage()
            sys.exit()


    if not pdbqs_filename:
        print 'pdbqs_to_pdbqt: pdbqs_filename_stem must be specified.'
        usage()
        sys.exit()

    #what about nucleic acids???
    mols = Read(pdbqs_filename)
    if verbose: print 'read ', pdbqs_filename
    mol = mols[0]
    #fix number for problem files with alternative positions
    mol.allAtoms.number = range(1, len(mol.allAtoms)+1)
    #possible clean-up???
    #need to type atoms + assign babel_types
    AD4_typer = AutoDock4_AtomTyper(set_aromatic_carbons=False)
    AD4_typer.setAutoDockElements(mol)
    if verbose:
        print "set autodock4 autodock_element for ", mol.name
    if charges_to_add is not None:
        preserved = {}
        preserved_types = preserve_charge_types.split(',') 
        for t in preserved_types:
            if not len(t): continue
            ats = mol.allAtoms.get(lambda x: x.autodock_element==t)
            for a in ats:
                if a.chargeSet is not None:
                    preserved[a] = [a.chargeSet, a.charge]

    if charges_to_add=='gasteiger':
        GasteigerChargeCalculator().addCharges(mol.allAtoms)

    if charges_to_add is not None:
        #restore any previous charges
        for atom, chargeList in preserved.items():
            atom._charges[chargeList[0]] = chargeList[1]
            atom.chargeSet = chargeList[0]

    #pdbqt_filename = mol.name + '.pdbqt'
    writer = PdbqtWriter()
    fptr = open(pdbqt_filename, 'w')
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
        print "wrote ", ctr, " atoms to", pdbqt_filename
    

# To execute this command type:
# pdbqs_to_pdbqt.py -s pdbqs_filename_stem  -v




