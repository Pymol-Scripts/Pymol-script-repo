#!/usr/bin/env python
#
# 
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/summarize_results.py,v 1.3.4.1 2009/05/26 14:47:22 rhuey Exp $
#
import os, glob

from MolKit import Read
from MolKit.hydrogenBondBuilder import HydrogenBondBuilder
from AutoDockTools.Docking import Docking
from mglutil.math.rmsd import RMSDCalculator
from PyAutoDock.InternalEnergy import InternalEnergy
from PyAutoDock.MolecularSystem import MolecularSystem
from PyAutoDock.AutoDockScorer import AutoDock305Scorer




if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: summarize_results.py -d directory"
        print
        print "    Description of command..."
        print "         -d     directory"
        print "    Optional parameters:"
        print "        [-t]    rmsd tolerance (default is 1.0)"
        print "        [-b]    print best docking info only (default is print all)"
        print "        [-o]    output filename"
        print "                      (default is 'summary_of_results')"
        print "        [-a]    append to  output filename"
        print "                      (default is to open output filename 'w')"
        print "        [-k]    build hydrogen bonds"
        print "        [-r]    receptor filename"
        print "        [-e]    report energy breakdown"
        print "        [-v]    verbose output"
        print "                      (default is leave active)"


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'd:o:t:r:bkaevh')
    except getopt.GetoptError, msg:
        print 'summarize_results.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    #-d: directory
    directory =  None
    #-t: rms_tolerance
    rms_tolerance =  1.0
    # optional parameters
    verbose = None
    #-o outputfilename
    outputfilename = "summary_of_results"
    #-r receptor_filename
    receptor_filename = ""
    #-b print best only
    print_best_only = False
    #-k build hydrogen bonds 
    build_hydrogen_bonds = False
    #-a append to outputfile
    append_to_outputfile = False
    #-e report_energy_breakdown 
    report_energy_breakdown = False


    #'d:o:t:r:bkaevh'
    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('-d', '--d'):
            directory = a
            if verbose: print 'set directory to ', a
        if o in ('-o', '--o'):
            outputfilename = a
            if verbose: print 'set outputfilename to ', a
        if o in ('-t', '--t'):
            rms_tolerance = float(a)
            if verbose: print 'set rms_tolerance to ', a
        if o in ('-r', '--r'):
            receptor_filename = a
            if verbose: print 'set receptor_filename to ', a
        if o in ('-b', '--b'):
            print_best_only = True
            if verbose: print 'set print_best_only to ', True
        if o in ('-k', '--k'):
            build_hydrogen_bonds = True
            if verbose: print 'set build_hydrogen_bonds to ', True
        if o in ('-a', '--a'):
            append_to_outputfile = True
            if verbose: print 'set append_to_outputfile to ', True
        if o in ('-e', '--e'):
            report_energy_breakdown = True
            if verbose: print 'set report_energy_breakdown to ', True
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print 'set verbose to ', True
        if o in ('-h', '--'):
            usage()
            sys.exit()

    if outputfilename=='summary_of_results':
        outputfilename = outputfilename + '_' + str(rms_tolerance)
    if build_hydrogen_bonds is True and receptor_filename is None:
        print 'to build hydrogen bonds, receptor filename must be specified'
        sys.exit()
    if report_energy_breakdown is True and receptor_filename is None:
        print 'to report energy breakdown, receptor filename must be specified'
        sys.exit()
    if not  directory:
        print 'summarize_results: directory must be specified.'
        usage()
        sys.exit()

    #read all the docking logs in as one Docking
    dlg_list = glob.glob(directory + '/*.dlg')
    d = Docking()
    for dlg in dlg_list:
        d.readDlg(dlg)

    #setup rmsd tool
    coords = d.ligMol.allAtoms.coords[:]
    d.clusterer.rmsTool = RMSDCalculator(coords)

    d.clusterer.make_clustering(rms_tolerance) 

    # for building hydrogen bonds or reporting energy breakdown
    # setup receptor:
    if build_hydrogen_bonds or report_energy_breakdown:
        d.ligMol.buildBondsByDistance()
        receptor_filename = directory + '/' + receptor_filename
        receptor = Read(receptor_filename)[0]
        receptor.buildBondsByDistance()
    if build_hydrogen_bonds:
        hbondBuilder = HydrogenBondBuilder()
        #do next two lines for each conf
        #d.ligMol.set_conformation(conf)
        #atDict = hbondBuilder.build(receptor.allAtoms, d.ligMol.allAtoms)
        #num_hbonds= len(atDict)
    if report_energy_breakdown:
        ms = MolecularSystem()
        ms.add_entities(receptor.allAtoms)
        ms.add_entities(d.ligMol.allAtoms)
        ad305scorer = AutoDock305Scorer()
        ad305scorer.set_molecular_system(ms)
        #whenever the ligMol changes conf, 
        #need to call ms.set_coords(1, new_coords)

    mode = 'w'
    if append_to_outputfile:
        mode = 'a'
    fptr = open(outputfilename, mode)
    if verbose: print "print_best_only=", print_best_only
    if print_best_only:
        clust0 = d.clusterer.clustering_dict[rms_tolerance][0]
        c = clust0[0]
        #find the filename of the best result
        for k in d.dlo_list:
            if c in k.conformations:
                dlo_filename = os.path.splitext(k.parser.filename)[0]
                if verbose: print "set dlo_filename to ", dlo_filename
                ostr = dlo_filename + ", "
                fptr.write(ostr)
        score_list = []
        #update the coords if you need to
        if build_hydrogen_bonds or report_energy_breakdown:
            d.ch.set_conformation(c)
        if report_energy_breakdown:
            ms.set_coords(1, c.getCoords()[:])
            score_list = ad305scorer.get_score_per_term()
        if not build_hydrogen_bonds:
            ostr = "% 4d, % 4d, % 8.4f, % 8.4f, % 8.4f" %(len(d.clusterer.data), len(clust0), c.docking_energy, c.binding_energy, c.getRMSD())
        else:
            atDict = hbondBuilder.build(receptor.allAtoms, d.ligMol.allAtoms)
            num_hbonds = len(atDict)
            ostr = "% 4d, % 4d, % 8.4f, % 8.4f, % 8.4f, % 2d" %(len(d.clusterer.data), len(clust0), c.docking_energy, c.binding_energy, c.getRMSD(), num_hbonds)
        for score in score_list:
            ostr = ostr + ", % 8.4f" %score
        ostr = ostr + "\n"
        fptr.write(ostr)
    else:
        for clust in d.clusterer.clustering_dict[rms_tolerance]:
            if verbose:
                print "len(d.clusterer.clustering_dict)=", len(d.clusterer.clustering_dict[rms_tolerance])
                print "len(clust)=", len(clust)
            #for each cluster in the clustering at this rms_tolerance:
            #output length of cluster and energies of the lowest energy result
            c = clust[0]
            if verbose: print "setting c to clust[0]", id(c)
            #write the dlg filename for the best in this cluster
            for k in d.dlo_list:
                if c in k.conformations:
                    dlo_filename = os.path.splitext(k.parser.filename)[0]
                    ostr = dlo_filename + ", " 
                    fptr.write(ostr)
            score_list = []
            if build_hydrogen_bonds or report_energy_breakdown:
                if verbose: print "setting conf to ", id(c)
                d.ch.set_conformation(c)
            if report_energy_breakdown:
                ms.set_coords(1, c.getCoords()[:])
                score_list = ad305scorer.get_score_per_term()
            if not build_hydrogen_bonds:
                ostr = "% 4d, % 8.4f, % 8.4f, % 8.4f" %(len(clust), c.docking_energy, c.binding_energy, c.getRMSD())
            else:
                atDict = hbondBuilder.build(receptor.allAtoms, d.ligMol.allAtoms)
                num_hbonds = len(atDict)
                ostr = "% 4d, % 8.4f, % 8.4f, % 8.4f, % 2d" %(len(clust), c.docking_energy, c.binding_energy, c.getRMSD(), num_hbonds)
            for score in score_list:
                ostr = ostr + ", % 8.4f" %score
            ostr = ostr + "\n"
            fptr.write(ostr)
    fptr.close()
    

# To execute this command type:
# summarize_results.py -d directory -t rmsd tolerance -r receptor filename  
#                         -b report best_docking_only -k build_hydrogen_bonds 
#                         -e report_energy_breakdown 
#                         -v verbose
# NOTE: -d, -t and -r require arguments whereas the others set booleans

