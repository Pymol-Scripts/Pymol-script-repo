
#!/usr/bin/env python

import os
from string import find

from MolKit import Read
from MolKit.pdbWriter import PdbqWriter
from MolKit.molecule import AtomSet

from AutoDockTools.MoleculePreparation import AD4LigandPreparation


if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: repair_ligand4.py -s filename"
        print
        print "    Description of command..."
        print "        [-f]    pdbqt_filename"
        print "    Optional parameters:"
        print "        [-o]    alternative pdbqt_filename"
        print "        (default is 'repaired_' +pdbqt_filename)"
        print "        [-v]    verbose output"


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'f:o:v')

    except getopt.GetoptError, msg:
        print 'repair_ligand4.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    #-f: pdbqt_filename_stem
    pdbqt_filename =  None

    # optional parameters
    verbose = None
    outputfilename =  None

    #'f:o:v'
    for o, a in opt_list:
        if o in ('-f', '--f'):
            pdbqt_filename = a
            if verbose: print 'set pdbqt_filename to ', pdbqt_filename
            outputfilename =  'repaired_' + pdbqt_filename
        if o in ('-o', '--o'):
            outputfilename = a 
            if verbose: 
                print 'set output outputfilename to ', a
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print 'set verbose to ', True
        if o in ('-h', '--'):
            usage()
            sys.exit()


    if not pdbqt_filename:
        print 'repair_ligand4: pdbqt_filename must be specified.'
        usage()
        sys.exit()

    mol = Read(pdbqt_filename)[0]
    if verbose: print 'read ', pdbqt_filename
    mol.buildBondsByDistance()
    mol.LPO = AD4LigandPreparation(mol, mode='interactive', repairs='', charges_to_add=None,
                root=0, outputfilename=outputfilename, cleanup='')
    #rebuild torTree
    allAts = mol.allAtoms
    for b in allAts.bonds[0]:
        if b.activeTors: b.activeTors = 0
    torscount = 0
    tM = mol.torTree.torsionMap
    for i in range(len(tM)):
        bnum0, bnum1 = tM[i].bond
        a0 = allAts.get(lambda x: x.number==bnum0 + 1)[0]
        a0.tt_ind = bnum0
        a1 = allAts.get(lambda x: x.number==bnum1 + 1)[0]
        a1.tt_ind = bnum1
        b = AtomSet([a0,a1]).bonds[0]
        if hasattr(b, 'possibleTors'):
            assert b.possibleTors
        else:
            b.possibleTors = 1
        b.activeTors = 1
        torscount = torscount + 1
    mol.torscount = torscount
    mol.ROOT = mol.allAtoms[0]
    mol.ROOT.rnum0 = 0
    mol.LPO.write(outputfilename)


# To execute this command type:
# repair_ligand4.py -f pdbqt_filename [-o outputfilename] -v




