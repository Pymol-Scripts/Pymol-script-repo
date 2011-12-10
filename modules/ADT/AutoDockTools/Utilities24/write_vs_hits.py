#!/usr/bin/env python
#
# $Id: write_vs_hits.py,v 1.1 2009/01/09 20:40:32 rhuey Exp $
#

import os

from AutoDockTools.Docking import Docking

if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        """Print helpful, accurate usage statement to stdout."""
        print "Usage: write_vs_hits.py -s summary_filename -d path_to_dockings"
        print
        print "    Write result file of top conformations from VS."
        print
        print "    Required parameters:"
        print "         -s name of summary file, eg 'summary_2.0.sort', produced in exercise10 of VS"
        print "         -p path_to_dockings containing files listed in summary file"
        print "    Optional parameters:"
        print "        [-f pdbqt_filename] file to create. Default is 'vs_results.pdbqt'"
        print "        [-A]             output ATOM and HETATM records only "
        print "        [-v]             verbose output"
        print "        [-D]             debugging output"


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 's:p:f:Avh')
    except getopt.GetoptError, msg:
        print 'write_vs_hits.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    #-s: summary_filename
    summary_filename =  None
    #-p: path_to_dockings
    path_to_dockings =  None

    # optional parameters
    #-f: pdbqtfilename
    pdbqt_filename =  "vs_results.pdbqt"
    verbose = None
    atom_records_only = False

    for o, a in opt_list:
        if o in ('-s', '--s'):
            summary_filename = a
            if verbose: print "set summary_filename  %s" % summary_filename
        if o in ('-p', '--p'):
            path_to_dockings = a
            if verbose: print "set path_to_dockings  %s" % path_to_dockings
        if o in ('-f', '--f'):
            pdbqt_filename = a
            if verbose: print "set pdbqt_filename to  %s" % pdbqt_filename
        if o in ('-A', '--A'):
            atom_records_only = True
            if verbose: print "set print ATOM and HETATM records only to %s" % atom_records_only
        if o in ('-v', '--v'):
            if verbose: print "set verbose to  %s" % verbose
        if o in ('-h', '--'):
            usage()
            sys.exit(2)

    # make sure summary_filename and path_to_dockings were specified
    if not summary_filename:
        print 'write_vs_hits: summary_filename must be specified.'
        usage()
    if not path_to_dockings:
        print 'write_vs_hits: path_to_dockings must be specified.'
        usage()
        sys.exit(2)


    # make sure summary_filename exists
    try:
        os.stat(summary_filename)
    except OSError, msg:
        print "write_vs_hits: %s: %s",  msg.strerror, msg.filename
        sys.exit()

    fileptr = open(summary_filename, 'r')
    all_lines = fileptr.readlines()
    fileptr.close()
    if not len(all_lines): 
        print "write_vs_hits: no useable lines in %s", summary_filename
        sys.exit()
    #ctr will be recorded in the residue field (18-20, 1-based)
    ctr = 1
    optr = open(pdbqt_filename, 'w')
    for line in all_lines:
        #sample line to use:
        #ZINC00057384_x1hpv/ZINC00057384_x1hpv,  20,  4, 10, -9.0900....
        #skip the header line
        if line.find('lowestEnergy')>-1: 
            #lowestEnergy dlgfn     #runs #cl #LEC    LE      rmsd     largestCl_dlgfn ...
            continue
        #write line from summary_filename
        ostr = "REMARK %s"%line
        optr.write(ostr)
        #find the docking log filename
        ll = line.split()
        fn = ll[0].strip()[:-1]
        dlgFN = path_to_dockings + "/" + fn +".dlg"
        #create a docking
        d = Docking()
        d.readDlg(dlgFN)
        #set ligand to lowest energy conformation LE
        cluD = d.clusterer.clustering_dict
        rms_tol = cluD.keys()[0]
        LE = cluD[rms_tol][0][0]
        d.ch.set_conformation(LE)
        #get coordinates of LE conformation
        crds = d.ligMol.allAtoms.coords[:]
        #use the parser to get newlines to output
        lines_to_print = d.ligMol.parser.write_with_new_coords(crds)
        #edit these lines, adding new line + changing the RES field to current count
        for nl in lines_to_print: 
            if nl[-1]!="\n":
                nl += "\n"
            if atom_records_only:
                if nl.find("ATOM")!=0 and nl.find("HETA")!=0: 
                    #skip lines of ligand keywords, eg ROOT
                    continue
            #replace the RES field with current count
            nl = nl[:17] + "%3s"%ctr + nl[20:]
            optr.write(nl)
        if verbose: print "onto next vs hit line"
        ctr+=1
    optr.close()
# To execute this command type:
# write_vs_hits.py -s summary_filename -o pdbqt_filename -v

