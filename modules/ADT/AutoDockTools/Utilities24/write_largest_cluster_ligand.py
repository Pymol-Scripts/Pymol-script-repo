#!/usr/bin/env python
#
# 
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/write_largest_cluster_ligand.py,v 1.3 2008/05/28 15:40:15 rhuey Exp $
#
import os, glob

from MolKit import Read
from AutoDockTools.Docking import Docking
from mglutil.math.rmsd import RMSDCalculator




if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        print "Usage: write_largest_cluster_ligand.py "
        print "    This script does the following: "
        print "         (1) read all the files with extension '.dlg' into one Docking"
        print "         (2) compute a clustering at the specified rms tolerance "
        print "         (3) write the ligand with the coordinates of the "
        print "             lowest-energy conformation in the largest cluster to a file" 
        print "    "
        print "    Optional parameters:"
        print "        [-t]    rms_tolerance (default 1.5)"
        print "        [-o pdbqt_filename] (default ligandstem_BC.pdbqt)"
        print "        [-v]    verbose output"

    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 't:o:vh')
    except getopt.GetoptError, msg:
        print 'write_largest_cluster_ligand.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    # optional parameters
    verbose = None
    #-o outputfilename
    outputfilename = None
    #-t rms_tolerance
    rms_tolerance = 1.5

    #'t:o:vh'
    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('-o', '--o'):
            outputfilename = a
            if verbose: print 'set outputfilename to ', a
        if o in ('-t', '--t'):
            rms_tolerance = float(a)
            if verbose: print 'set rms_tolerance to ', a
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print 'set verbose to ', True
        if o in ('-h', '--'):
            usage()
            sys.exit()

    dlg_list = glob.glob('*.dlg')
    d = Docking()
    for dlg in dlg_list:
        d.readDlg(dlg)
    d.clusterer.rmsTool = RMSDCalculator(d.ligMol.allAtoms.coords[:])
    d.clusterer.make_clustering(rms_tolerance)
    clustering = d.clusterer.clustering_dict[rms_tolerance]
    largest = clustering[0]
    for clust in clustering:
        if len(clust)>len(largest):
            largest = clust
    #update the coordinates with those of conf with lowest energy in this largest cluster
    d.ch.set_conformation(largest[0])
    parser = d.ligMol.parser
    lines = []
    #have to add newline character to lines read from dlg
    for l in parser.allLines:
        l+= '\n'
        lines.append(l)
    parser.allLines = lines
    coords = d.ligMol.allAtoms.coords
    if outputfilename is None:
        outputfilename = d.ligMol.name  + '_BC.pdbqt'
    parser.write_with_new_coords(coords, outputfilename) 
    if verbose:
        print 'wrote %s' %outputfilename


# To execute this command type:
# write_largest_cluster_ligand.py  [-t rms_tolerance, -o outputfilename] -v




