#!/usr/bin/env python
#
# 
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/prepare_gpf.py,v 1.6 2007/10/09 17:30:07 annao Exp $
#

import string
import os.path
from math import ceil
from MolKit import Read
from AutoDockTools.GridParameters import GridParameters, grid_parameter_list
from AutoDockTools.GridParameters import GridParameterFileMaker


def usage():
    print "Usage: prepare_gpf.py -l pdbq_file -r pdbqs_file "
    print "     -l ligand_filename"
    print "     -r receptor_filename"
    print
    print "Optional parameters:"
    print "    [-i reference_gpf_filename]"
    print "    [-o output_gpf_filename]"
    print "    [-p parameter=newvalue]"
    print "    [-v]"
    print
    print "Prepare a grid parameter file (GPF) for AutoDock."
    print
    print "   The GPF will by default be <receptor>.gpf. This"
    print "may be overridden using the -o flag."

    
if __name__ == '__main__':
    import getopt
    import sys

    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'vl:r:i:o:p:')
    except getopt.GetoptError, msg:
        print 'prepare_gpf.py: %s' % msg
        usage()
        sys.exit(2)

    receptor_filename = ligand_filename = None
    list_filename = gpf_filename = gpf_filename = None
    output_gpf_filename = None
    parameters = []
    verbose = None
    for o, a in opt_list:
        if o in ('-v', '--v'):
            verbose = 1
        if o in ('-l', '--l'):
            ligand_filename = a
            if verbose: print 'ligand_filename=', ligand_filename
        if o in ('-r', '--r'):
            receptor_filename = a
            if verbose: print 'receptor_filename=', receptor_filename
        if o in ('-i', '--i'):
            gpf_filename = a
            if verbose: print 'reference_gpf_filename=', gpf_filename
        if o in ('-o', '--o'):
            output_gpf_filename = a
            if verbose: print 'output_gpf_filename=', output_gpf_filename
        if o in ('-p', '--p'):
            parameters.append(a)
            if verbose: print 'parameters=', parameters
        if o in ('-h', '--'):
            usage()
            sys.exit()


    if (not receptor_filename) or (not ligand_filename):
        print "prepare_gpf.py: ligand and receptor filenames"
        print "                    must be specified."
        usage()
        sys.exit()

    gpfm = GridParameterFileMaker(verbose=verbose)
    gpfm.set_ligand(ligand_filename)
    gpfm.set_receptor(receptor_filename)
    if gpf_filename is not None:
        gpfm.read_reference(gpf_filename)
    for p in parameters:
        key,newvalue = string.split(p, '=')
        kw = {key:newvalue}
        apply(gpfm.set_grid_parameters, (), kw)
    #gpfm.set_grid_parameters(spacing=1.0)
    gpfm.write_gpf(output_gpf_filename)


#prepare_gpf.py -l indinavir.pdbq -r 1hsg.pdbqs -p spacing=0.4 -p npts=[60,60,60] -i ref.gpf -o testing.dpf 

