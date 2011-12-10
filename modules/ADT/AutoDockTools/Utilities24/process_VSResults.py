#!/usr/bin/env python
#
# 
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/process_VSResults.py,v 1.16.2.1 2011/03/29 17:39:49 rhuey Exp $
#
# $Id: process_VSResults.py,v 1.16.2.1 2011/03/29 17:39:49 rhuey Exp $
#
import os, glob, time

from MolKit import Read
from AutoDockTools.Docking import Docking, DockingResultProcessor
from mglutil.math.rmsd import RMSDCalculator
from AutoDockTools.InteractionDetector import InteractionDetector
from MolKit.hydrogenBondBuilder import HydrogenBondBuilder


def get_filename(d, conf):
    #find the filename of dlg containing this conformation
    for k in d.dlo_list:
        if conf in k.conformations:
            dlo_filename = os.path.splitext(k.parser.filename)[0]
            break
    return dlo_filename




if __name__ == '__main__':
    import sys
    import getopt

    #'d:t:f:r:BLx:Dnc vh'

    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: process_VSResults.py -d directory"
        print
        print "    Description of command..."
        print "         -d     directory"
        print "    Optional parameters:"
        print "        [-t]    rmsd tolerance (default is 2.0)"
        print "        [-f]    rmsd reference filename "
        print "        (default is to use input ligand coordinates from a docking log)"
        print "        [-r]    receptor filename (default is set from gridmap names)"
        print "        [-B]    create best docking pdbqt only (default is create both best energy and largest cluster )"
        print "        [-L]    create largest cluster docking pdbqt only (default is create both best energy and largest cluster )"
        print "        [-l]    stem for largest cluster docking pdbqt file (default is 'ligandname_lc')"
        print "        [-x]    maximum number of clusters to report (default is all)"
        print "        [-Z]    do not include interactions in output pdbqt file(default is to include interactions)"
        print "        [-n]    do not build hydrogen bonds (default is to build hydrogen bonds + report)"
        print "        [-c]    do not detect atoms in close contact (default is to detect + report)"
        print "        [-p]    include detection of pi-pi interactions in close contact (default is not to detect pi-pi and pi-cation)"
        print "        [-H]    custom HydrogenBondBuilder parameters: comma-separated-list eg distCutoff=3.00,distCutoff2=3.25"
        print "        [-v]    verbose output"


    #HydrogenBondBuilder default parameter values:
    #   distCutoff=2.25     #H-Acc Distance cutoff
    #   distCutoff2=3.00    #Donor-Acc Distance cutoff
    #   d2min=120           #sp2-donor-hydrogen-acceptor angle 
    #   d2max=180           #
    #   d3min=120           #sp3-donor-hydrogen-acceptor angle
    #   d3max=180           #
    #   a2min=110           #sp2-donor-acceptor-acceptorN angle, 
    #   a2max=150           #    where acceptorN is 'neighbor' atom bonded to acceptor
    #   a3min=100           #sp3-donor-acceptor-acceptorN angle
    #   a3max=150           #    where acceptorN is 'neighbor' atom bonded to acceptor

    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'd:t:f:r:BLl:x:DZncpH:vh')
    except getopt.GetoptError, msg:
        print 'process_VSResults.py: %s' %msg
        usage()
        sys.exit(2)
    
    # initialize required parameters
    #-d: directory
    directory =  None

    # optional parameters
    #-o: outputfilename
    outputfilename =  None
    #-t: rms_tolerance
    rms_tolerance =  2.0
    #-f: ligand rms reference
    rms_reference = None
    #-r receptor_filename 
    #   if None, set from gridmap names
    receptor_filename = None
    # default is create 2 pdbqt files
    write_both = True
    #-B best only
    best_only = False
    #-L largestCl only
    largestCl_only = False
    #-l stem for largest cluster pdbqtfile 
    lc_stem = ""
    #-x max_cl_to_write 
    #?SHOULD default be all? if so initialize to -1
    #max_cl_to_write = 10
    max_cl_to_write = -1

    #-D include_interactions in output pdbqt file
    #-Z do not include_interactions in output pdbqt file
    include_interactions = True
    #-p include detection of pi-pi and pi-cation interactions 
    detect_pi = False

    #-n do NOT build hydrogen bonds
    build_hbonds = True
    #-c do NOT detect_close_contscts
    detect_close_contacts = True

    #-H custom hydrogen bond parameters
    custom_hb_parameters = ""

    # -v verbose output
    verbose = None

    # d:t:f:r:BLl:x:DZncpH:vh

    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('-d', '--d'):
            directory = a
            if verbose: print 'set directory to ', a
        if o in ('-t', '--t'):
            rms_tolerance = float(a)
            if verbose: print 'set rms_tolerance to ', rms_tolerance
        if o in ('-f', '--f'):
            rms_reference = a
            if verbose: print 'set rms_reference filename to ', rms_reference
        if o in ('-r', '--r'):
            receptor_filename = a
            if verbose: print 'set receptor_filename to ', receptor_filename
        if o in ('-B', '--B'): #create pdbqt file for bestenergy result only
            best_only = True
            write_both = False
            if verbose: print 'set best_only to ', best_only
        if o in ('-L', '--L'): #create pdbqt file for largestcluster result only
            largestCl_only = True
            write_both = False
            if verbose: print 'set largestCL_only to ', largestCL_only
        if o in ('-l', '--l'): #stem for pdbqt file for be in largest cluster
            lc_stem = a
            if verbose: print 'set largest_cluster stem to ', lc_stem

        if o in ('-x', '--x'):
            max_cl_to_write = int(a)
            if verbose: print 'set maximum number of clusters to write to ', max_cl_to_write

        if o in ('-Z', '--Z'):
            include_interactions = False
            if verbose: print 'set include_interactions to ', include_interactions

        if o in ('-n', '--n'):
            build_hbonds = False
            if verbose: print 'set build_hbonds to ', build_hbonds
        if o in ('-c', '--c'):
            detect_close_contacts = False
            if verbose: print 'set detect_close_contacts to ', detect_close_contacts
        if o in ('-p', '--p'):
            detect_pi = True
            if verbose: print 'set detect_pi interactions to ', detect_pi

        if o in ('-H', '--H'):
            custom_hb_parameters = a
            if verbose: print 'custom_hb_parameters to ', custom_hb_parameters
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print 'set verbose to ', True
        if o in ('-h', '--'):
            usage()
            sys.exit()

    if not  directory:
        print 'process_VSResults: directory must be specified.'
        usage()
        sys.exit()

    intF = None
    if len(custom_hb_parameters):
        if verbose: print "using custom hb parameters!"
        #custom_hb_parameters "distCutoff=3.00,distCutoff2=3.25"
        hb_list = custom_hb_parameters.split(',')
        kw = {}
        for item in hb_list:
            param, value = item.split('=')
            kw[param] = float(value)
            if verbose: print param, '=', value
        hbB = apply(HydrogenBondBuilder, (), kw)
        intF = InteractionDetector(receptor_filename, detect_pi=detect_pi, hydrogen_bond_builder=hbB)

    drp = DockingResultProcessor(rms_tolerance=rms_tolerance,
                                 rms_reference=rms_reference,
                                 receptor_filename=receptor_filename,
                                 write_both=write_both,
                                 best_only=best_only,
                                 largestCl_only=largestCl_only,
                                 lc_stem=lc_stem, 
                                 max_cl_to_write=max_cl_to_write, 
                                 include_interactions=include_interactions, 
                                 detect_pi=detect_pi, 
                                 build_hbonds=build_hbonds, 
                                 detect_close_contacts=detect_close_contacts, 
                                 intF=intF) 
    drp.process(directory, verbose=verbose)

