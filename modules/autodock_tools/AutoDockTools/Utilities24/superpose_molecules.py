#!/usr/bin/env python
#
# 
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/superpose_molecules.py,v 1.1 2009/01/07 18:52:50 rhuey Exp $
#

import os, string
from MolKit import Read




if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        print "Usage: superpose_molecules.py -d datafile"
        print "   Each line in the datafile consists of the filename of a molecule and a list of "
        print "   comma-separated atom numbers. Here is an example: "
        print "       reference.pdb 6110,2448,1938,2470,2493,4619  "
        print "       molecule1.pdb 31,15,1,16,18,20  "
        print "       ...........................     "
        print "   The molecule from first line is used as the reference."
        print "   The molecule listed in each subsequent line is superimposed onto the reference"
        print "   and written to a new file with the new coordinates."
        print "   The number of entries in the comma-separated lists must be the same for each entry (and more than 2) "
        print "        -d    data_filename "
        print "   Optional parameters:"
        print "        [-p]    prefix for new filenames (default 'transposed_')"
        print "        [-v]    verbose output"

    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'd:p:vh')
    except getopt.GetoptError, msg:
        print 'superpose_molecules.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    # optional parameters
    verbose = None
    #-d data_filename
    data_filename = None
    #-p prefix
    prefix = "transposed_"

    #'d:p:vh'
    for o, a in opt_list:
        if o in ('-d', '--d'):
            data_filename = a
            if verbose: print 'set data_filename to ', a
        if o in ('-p', '--p'):
            prefix = a
            if verbose: print 'set prefix to ', a
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print 'set verbose to ', True
        if o in ('-h', '--'):
            usage()
            sys.exit()
    if not data_filename:
        print "superpose_molecules.py: data_filename must be specified."
        usage()
        sys.exit()

    #file specified by data_filename contains filenames + atom indices
    #refilename 0,1,2,3,4,5,6
    #filename1 10,31,22,63,44,51,16
    #filename2 0,1,17,23,24,35,56
    #...
    # NOTE each line has filename and list of indices
    # NOTE the first entry is the reference
    # NOTE THE NUMBER OF INDICES must match,indices must be comma-separated
    #                   (with NO intervening spaces)

    #read the data file
    fptr = open(data_filename)
    lines = fptr.readlines()
    fptr.close()

    #setup the reference molecule
    from MolKit import Read
    from MolKit.molecule import AtomSet
    ref_name, ref_ll = lines[0].split()
    ref_mol = Read(ref_name)[0]
    ref_ats = AtomSet()
    ref_indices = map(int, ref_ll.split(','))
    num_ref_indices = len(ref_indices)
    if num_ref_indices<3: 
        print " At least 3 indices are required! only %d indices found!!!"%(num_ref_indices)
        exit()
    for i in ref_indices:
        #ref_ats.append(ref_mol.allAtoms[i-1])
        ref_ats.append(ref_mol.allAtoms.get(lambda x: x.number==i)[0])

    num_ref_ats = len(ref_ats)

    #setup the tool for the rigid fit
    from mglutil.math.rigidFit import RigidfitBodyAligner
    RFA = RigidfitBodyAligner()
    RFA.setRefCoords(ref_ats.coords)

    #transform the rest of the molecules
    for l in lines[1:]:
        name, ll = l.split() 
        ll.strip()
        moving_mol = Read(name)[0]
        indices = map(int, ll.split(','))
        moving_ats = AtomSet()
        for i in indices: 
            new_atoms = moving_mol.allAtoms.get(lambda x: x.number==i)
            if not len(new_atoms):
                print "There is no atom in %s with number %d" %(name, i)
                break
            else:
                moving_ats.append(new_atoms[0])
            #moving_ats.append(moving_mol.allAtoms.get(lambda x: x.number==i)[0])
        num_moving_ats = len(moving_ats)
        if num_moving_ats!=num_ref_ats:
            print "Number of moving_ats (%d) does not match number of reference atoms (%d). Skipping %s" %(num_moving_ats,num_ref_ats,name)
            continue
        #for i in indices: moving_ats.append(moving_mol.allAtoms[i-1])
        RFA.rigidFit(moving_ats.coords)
        newCoords = RFA.transformCoords(moving_mol.allAtoms.coords)
        if not newCoords is None and newCoords.shape[1]!=3:
            newCoords=newCoords[:,:3]
        outputfilename = prefix + name
        moving_mol.parser.write_with_new_coords(newCoords,outputfilename) 
        if verbose: print 'wrote %s' %outputfilename


# To execute this command type:
# superpose_molecules.py  -d data_filename [-p prefix] -v




