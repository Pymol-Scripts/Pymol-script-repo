#!/usr/bin/env python
#
# 
##############################################################################
#
# Authors: Stefano FORLI , Ruth HUEY
#
# Copyright: A. Olson TSRI 2011
#
#############################################################################

# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/process_VSResultsFox.py,v 1.1.2.1 2011/05/05 18:29:51 rhuey Exp $
#
# $Id: process_VSResultsFox.py,v 1.1.2.1 2011/05/05 18:29:51 rhuey Exp $ 
#

"""
This script creates a pdbqt+ file containing docked coordinates, docked
energy as well as  other important summary information about the docking. 
Input for the script is a 'directory' containing 1 or more docking log files
(*.dlg) from AutoDock calculations and the 'receptor filename'. The script
reclusters the results and includes a summary of the clustering in the pdbqt+.
"""
import os, glob, time
from AutoDockTools.Docking import FoxResultProcessor


if __name__ == '__main__':
    import sys
    import getopt

    #'d:t:f:r:BLx:Dnc vh'

    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: process_VSResultsFox.py -d directory"
        print
        print "    Description of command..."
        print "         -r     receptor filename "
        print "         -d     directory"
        print "    Optional parameters:"
        print "        [-t]    rmsd tolerance (default is 2.0)"
        print "        [-v]    verbose output"


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'r:d:t:vh')
    except getopt.GetoptError, msg:
        print 'process_VSResultsFox.py: %s' %msg
        usage()
        sys.exit(2)
    
    # initialize required parameters
    #-d: directory
    directory =  None
    #-r receptor_filename 
    receptor_filename = None

    # optional parameters
    #-t: rms_tolerance
    rms_tolerance =  2.0

    # -v verbose output
    verbose = None

    # r:d:t:vh
    for o, a in opt_list:
        if o in ('-r', '--r'):
            receptor_filename = a
            if verbose: print 'set receptor_filename to ', receptor_filename
        if o in ('-d', '--d'):
            directory = a
            if verbose: print 'set directory to ', a
        if o in ('-t', '--t'):
            rms_tolerance = float(a)
            if verbose: print 'set rms_tolerance to ', rms_tolerance
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print 'set verbose to ', True
        if o in ('-h', '--'):
            usage()
            sys.exit()

    if not directory:
        print 'process_VSResultsFox: directory must be specified.'
        usage()
        sys.exit()

    if receptor_filename is None:
        print 'process_VSResultsFox: receptor_filename must be specified.'
        usage()
        sys.exit()

    frp = FoxResultProcessor(receptor_filename)
    frp.process(directory=directory, rms_tolerance=rms_tolerance) 

