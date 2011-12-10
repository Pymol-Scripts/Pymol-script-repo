#!/usr/bin/env python
# $Id: superimpose_based_on_subset.py,v 1.2 2007/12/04 15:03:24 rhuey Exp $

#---------------------------------------------------------------------------#
#   Authors: Ruth Huey, Sandro Cosconati                                    #
#   Copyright: Ruth Huey, Sandro Cosconati  TSRI 2007                       #
#   usage:                                                                  #
#      superimpose_based_on_subsetfile.py -r reffile -m mobfile             #
#               -f subsetfile [-o outputfile -g mobsubsetfile -v verbose]   #
#                                                                           # 
#   subsetfile syntax: each line a single range of residues (zero-based)    #
#       15-17  #use residues at position 15 + 16                            #
#       67-68  #use residue at position  67                                 #
#---------------------------------------------------------------------------#


import os
from MolKit import Read
from MolKit.molecule import AtomSet
from mglutil.math.rigidFit import RigidfitBodyAligner


def get_atoms(mol, list_of_indicies, names_to_use=['N','CA','C'], verbose=False):
    if verbose: 
        print "in get_atoms with list of indicies:" 
        print list_of_indicies
    if not len(list_of_indicies):
        raise 'invalid input: list of indicies is empty!'
    atoms = AtomSet()
    num_res = 0
    for item in list_of_indicies:
        first, second = item
        #check for valid index and for end of chain
        max_index = len(mol.chains.residues)-1
        assert first<=max_index, 'invalid start of residue range'
        assert second<=max_index, 'invalid end of residue range'
        assert second>=first, 'second index cannot be smaller than first'
        if second==max_index:
            #ie mol.chains.residues[second]==mol.chains.residues[-1]:
            these_res = mol.chains.residues[first:]
        else:
            these_res = mol.chains.residues[first:second+1]
        if verbose: print "Adding  %3d residues " %(len(these_res)),
        num_res+=len(these_res)
        if verbose: print "Now there are %d residues total" %(num_res)
        for r in these_res:
            for n in names_to_use:
                atoms.append( r.atoms.get(n)[0])
    assert len(atoms), 'invalid input: lists of indicies did not correspond to any residues!'
    if verbose: print 'returning %d atoms' %(len(atoms))
    return atoms



if __name__ == "__main__":
    import sys
    import getopt

    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: superimpose_based_on_subset_file.py -r reffile -m mobfile"
        print " Description of  command..."
        print "       -r reference filename"
        print "       -m moving filename"
        print "       -f name of file listing reference residue ranges: "
        print "            one range per line such as:\n26-26\n35-39\n"
        print " Optional parameters:"
        print "      [-o] output filename (default is stem of moving file + '_superimposed.pdb')"
        print "      [-g] name of file listing mobile residue ranges, one range per line eg 32-36"
        print "            use this option if the numbering of the residues differs  "
        print "      [-x] string_for_atoms_to_use listing which atoms in each residue in ranges of "
        print "           residues read from file to use for superimposition"
        print "               (default is N_CA_C)"
        print "      [-v] verbose output"


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'r:m:f:o:x:g:vh')
    except getopt.GetoptError, msg:
        print "superimpose_based_on_subsetfile.py: %s" %msg
        usage()
        sys.exit(2)

    # initialize the required parameters
    #  -r:    filename for reference molecule
    reference_filename = None
    #  -m:    filename for moving molecule
    moving_filename = None
    #  -f:    filename for reference residue ranges
    refrange_filename = None
    #  -g:    filename for mobile residue ranges
    mobrange_filename = None
    #  -o:  outputfilename for moving molecule with new coords
    outputfilename = None

    #initialize the optional parameters
    verbose = None
    string_for_atoms_to_use = None
    mobrange_filename = None

    #r:m:f:o:x:vh
    for o, a in opt_list:
       if o in ('-r', '-r'):
           reference_filename = a
       if o in ('-m','-m'):
           moving_filename = a
       if o in ('-f','-f'):
           refrange_filename = a
       if o in ('-g', '-g'):
           mobrange_filename = a
       if o in ('-x', '-x'):
           string_for_atoms_to_use = a
       if o in ('-o', '-o'):
           outputfilename = a
       if o in ('-v', '-v'):
           verbose = True


    if not reference_filename:
        print "superimpose_based_on_subsetfile: reference filename must be specified"
        usage()
        sys.exit()
    if not moving_filename:
        print "superimpose_based_on_subsetfile: moving filename must be specified"
        usage()
        sys.exit()
    if not refrange_filename:
        print "superimpose_based_on_subsetfile: refrange_filename filename must be specified"
        usage()
        sys.exit()
     
    from MolKit import Read       #   Read always returns a MoleculeSet 
    refMol = Read(reference_filename)[0]  # the first Molecule in the set created 
    mobMol = Read(moving_filename)[0] # the first Molecule in the set created from 1hsg.pdb
    fptr = open(refrange_filename)
    #12-20\n36-40\n155-162\n175-178\n233-238\n338-345\n
    reflist = []
    lines = fptr.readlines()
    fptr.close()
    for l in lines:
        #sample line "12-20\n"
        #handle blank lines gracefully
        l = l.strip()
        if len(l):
            first, second  = map(int, l.split('-'))
            reflist.append((first,second))
    if verbose: print "reflist=", reflist
    moblist = []
    if mobrange_filename is not None:
        fptr = open(mobrange_filename)
        lines = fptr.readlines()
        fptr.close()
        for l in lines:
            #handle blank lines gracefully
            l = l.strip()
            if len(l):
                first,second  = map(int, l.split('-'))
                moblist.append((first,second))
        if verbose: print "moblist=", moblist
    else:
        moblist = reflist[:]
    #refMol.chains.residues.atoms       # all levels of hierarchical data structure
    #refMol.allAtoms                    # short cut to lowest level....
    #mobMol.allAtoms                    #each atom has 'coords' attribute
    from mglutil.math.rigidFit import RigidfitBodyAligner
    rigidfitAligner = RigidfitBodyAligner()
    if  string_for_atoms_to_use is None:
        string_for_atoms_to_use = "N_CA_C"
    names_to_use = string_for_atoms_to_use.split('_')
    refAtoms = get_atoms(refMol, reflist, names_to_use, verbose)
    refCoords = refAtoms.coords #or refMol.chains.residues.atoms.get('backbone').coords
    #check that there are some
    if len(refCoords)==0:
        print "NO REF COORDS!!!"
        raise 'Unable to superimpose because no reference coordinates were specified!!!' 
    mobAtoms = get_atoms(mobMol, moblist, names_to_use, verbose)
    mobCoords = mobAtoms.coords
    if len(mobCoords)==0:
        print "NO MOB COORDS!!!"
        raise 'Unable to superimpose because no mobile coordinates were specified!!!' 
    #check that moblist and reflist are of the same length    
    assert len(refCoords)==len(mobCoords)

    rigidfitAligner.setRefCoords(refCoords)  
    rigidfitAligner.rigidFit(mobCoords)
    #now rigidfitAligner has attributes 'rotationMatrix' and 'translationMatrix'

    newCoords = rigidfitAligner.transformCoords(mobMol.allAtoms.coords)
    #make sure the coordinates have the correct shape, ie (num_atoms, 3)
    if not newCoords is None and newCoords.shape[1]!=3:
        newCoords=newCoords[:,:3]   
    #newCoords is 2D Numeric array with dimensions  num_atoms, and length of each coordinate
    #if verbose: print "newCoords[0]=", newCoords[0]
    # end by writing a file with the transformed coordinates:
    #the parser of the moving molecule knows about all the input lines
    # so use it to write a file with the new coords:
    if outputfilename is None:
        EXT = os.path.splitext(refMol.parser.filename)[1]
        outputfilename = refMol.name + "_superimposed"+ EXT
    mobMol.parser.write_with_new_coords(newCoords, outputfilename)
    if verbose: print "wrote ", outputfilename

# To execute this command type:
# superimpose_based_on_subset_file.py -r reffile -m mobfile -o outfile [-x  string_for_atoms_to_use] [-v]
  

