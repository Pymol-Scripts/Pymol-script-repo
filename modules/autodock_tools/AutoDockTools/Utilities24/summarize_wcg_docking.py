#!/usr/bin/env python
#
# 
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/summarize_wcg_docking.py,v 1.3 2007/10/08 18:14:32 rhuey Exp $
#
import os, glob, numpy.oldnumeric as Numeric
from string import split, join

from MolKit import Read
from MolKit.hydrogenBondBuilder import HydrogenBondBuilder
from AutoDockTools.XMLParser import XMLParser
from AutoDockTools.Docking import Docking
from AutoDockTools.cluster import Clusterer
from mglutil.math.rmsd import RMSDCalculator
from PyAutoDock.MolecularSystem import MolecularSystem
from PyAutoDock.AutoDockScorer import AutoDock4Scorer



if __name__ == '__main__':
    import sys
    import getopt

    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: summarize_wcg_docking.py -d directory"
        print
        print "    Description of command..."
        print "         -d     directory"
        print "    Optional parameters:"
        print "        [-t]    rmsd tolerance (default is 1.0)"
        print "        [-b]    print best docking info only (default is print all)"
        print "        [-o]    output filename"
        print "                      (default is directory.txt)"
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
        print 'summarize_wcg_docking.py: %s' %msg
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
    outputfilename = 'results.txt'
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
            directory = a   #faah001....
            if verbose: print 'set directory to ', a
            #CAUTION: THIS IS NOT GENERAL!!!
            #faah001_diversity0001_masterpr_ref...faah0001_xtpv_masterpr_ref
            ll = directory.split('_')
            #ll = ['faah0001', 'diversity0001', 'masterpr', 'ref']
            #OR
            #ll = ['faah0001', 'xtpv', 'masterpr', 'ref'']
            #eg faah0001
            key = ll[0]
            #diversity0001.pdbqt
            lig_fn = ll[1]+".pdbqt"
            #masterpr_ref.pdbqt
            receptor_filename = join(ll[2:],'_') + ".pdbqt"
            #dpf_fn 
            dpf_fn = directory  + ".dpf"
            print "dpf_fn =", dpf_fn
            #dpf_fn = directory + "/" + directory + ".dpf"
            #SET the default outputfilename to directory + .txt here
            if outputfilename=='results.txt':
                outputfilename = directory + '.txt'
        if o in ('-o', '--o'):
            outputfilename = a
            if verbose: print 'set outputfilename to ', a
        if o in ('-t', '--t'):
            rms_tolerance = float(a)
            if verbose: print 'set rms_tolerance to ', a
        if o in ('-r', '--r'):
            #THIS IS SUPERFLUOUS
            receptor_filename = a
            print "warning! setting receptor_filename to ", a
            #if verbose: print 'set receptor_filename to ', a
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

    if not  directory:
        print 'summarize_wcg_docking: directory must be specified.'
        usage()
        sys.exit()

    #read all the xml docking logs in as one Docking
    
    xml_list = glob.glob('*.xml')
    print "xml_list =", xml_list
    p = XMLParser()
    d = Docking(parser=p)
    for xml_file in xml_list:
        print "calling readXMLResults with", xml_file
        d.readXMLResults(xml_file, dpf = dpf_fn)
    ligMol = d.ligMol
    ligAts = ligMol.allAtoms
    #setup rmsd tool
    coords = ligAts.coords[:]
    atom_ct = len(ligAts)
    torsion_ct = len(ligMol.torTree.torsionMap)
    tors_penalty = torsion_ct * 0.2744

    cl = Clusterer(d.ch.conformations)
    d.clusterer = cl
    cl.make_clustering(rms_tolerance)
    ref_coords = cl.clustering_dict[rms_tolerance][0][0].getCoords()[:]

    # for building hydrogen bonds or reporting energy breakdown
    # setup receptor:
    if build_hydrogen_bonds or report_energy_breakdown:
        ligMol.buildBondsByDistance()
        receptor = Read(receptor_filename)[0]
        receptor.buildBondsByDistance()
    if build_hydrogen_bonds:
        hbondBuilder = HydrogenBondBuilder()
    if report_energy_breakdown:
        ms = MolecularSystem()
        ms.add_entities(receptor.allAtoms)
        ms.add_entities(ligMol.allAtoms)
        adscorer = AutoDock4Scorer()
        adscorer.set_molecular_system(ms)

    mode = 'w'
    if append_to_outputfile:
        mode = 'a'
    fptr = open(outputfilename, mode)
    if mode=='w':
        if build_hydrogen_bonds and report_energy_breakdown:
            titles = "#clust  binding_nrg rmsd_v_0          Ki    int_nrg  hbnds sum_py_nrg estat_nrg     hbnd_nrg       vdw_nrg dsolv_nrg#ats #t\n"
        elif build_hydrogen_bonds:
            titles = "#clust  binding_nrg rmsd_v_0          Ki    int_nrg  hbnds #ats #t\n"
        elif report_energy_breakdown:
            titles = "#clust  binding_nrg rmsd_v_0          Ki    int_nrg  sum_py_nrg estat_nrg     hbnd_nrg       vdw_nrg dsolv_nrg#ats #t\n"
        else:
            titles = "#clust  binding_nrg rmsd_v_0          Ki    int_nrg #ats #t\n"
        fptr.write(titles)
    if print_best_only:
        clust0 = d.clusterer.clustering_dict[rms_tolerance][0]
        c = clust0[0]
        #find the filename of the best result
        #for k in d.dlo_list:
        #    if c in k.conformations:
        #        dlo_filename = os.path.splitext(k.parser.filename)[0]
        #        print "set dlo_filename to ", dlo_filename
        #        ostr = dlo_filename + ", "
        #        fptr.write(ostr)
        score_list = []
        #update the coords if you need to
        if build_hydrogen_bonds or report_energy_breakdown:
            d.ch.set_conformation(c)
        if report_energy_breakdown:
            ms.set_coords(1, c.getCoords()[:])
            score_list = adscorer.get_score_per_term()
            total_score = Numeric.add.reduce(score_list) + c.internal_energy + tors_penalty
            #total_score = Numeric.add.reduce(score_list)
        ostr = "% 4d, % 4d, % 12.4f, % 6.4f, % 8.6e,% 9.6f" %(len(d.clusterer.data), len(clust0),  c.binding_energy, c.getRMSD(ref_coords), c.Ki, c.internal_energy)
        if build_hydrogen_bonds:
            atDict = hbondBuilder.build(receptor.allAtoms, ligMol.allAtoms)
            num_hbonds = len(atDict)
            ostr = ostr + ", % 2d" %(num_hbonds)
        if report_energy_breakdown:
            ostr = ostr + ", % 12.4f" %total_score
            #for score in score_list:
            #    ostr = ostr + ", % 12.4f" %score
            ostr = ostr + ", % 6.4f, % 12.4f, % 12.4f,% 6.4f," %(score_list[0],
                        score_list[1], score_list[2], score_list[3])
        #add the number of atoms and the number of active torsions in the ligand
        ostr = ostr + "%4d%2d\n" %(atom_ct, torsion_ct)
        try:
            fptr.write(ostr)
        except:
            fptr.write("\n")
    else:
        for clust in d.clusterer.clustering_dict[rms_tolerance]:
            #for each cluster in the clustering at this rms_tolerance:
            #output length of cluster and energies of the lowest energy result
            c = clust[0]
            #write the dlg filename for the best in this cluster
            #for k in d.dlo_list:
            #    if c in k.conformations:
            #        dlo_filename = os.path.splitext(k.parser.filename)[0]
            #        print "dlo_filename=", dlo_filename
            #        #ostr = dlo_filename + "\n " 
            #        #fptr.write(ostr)
            score_list = []
            if build_hydrogen_bonds or report_energy_breakdown:
                d.ch.set_conformation(c)
            if report_energy_breakdown:
                ms.set_coords(1, c.getCoords()[:])
                score_list = adscorer.get_score_per_term()
                total_score = Numeric.add.reduce(score_list) + c.internal_energy + tors_penalty
            ostr = "% 4d, % 12.4f, % 8.4f, % 8.3e,% 9.6f" %(len(clust),  c.binding_energy, c.getRMSD(ref_coords), c.Ki, c.internal_energy)
            if build_hydrogen_bonds:
                atDict = hbondBuilder.build(receptor.allAtoms, ligMol.allAtoms)
                num_hbonds = len(atDict)
                ostr = ostr + ", % 2d" %(num_hbonds)
            if report_energy_breakdown:
                ostr = ostr + ", % 12.4f" %total_score
                ostr = ostr + ", % 6.4f, % 12.4f, % 12.4f,% 6.4f," %(score_list[0],
                            score_list[1], score_list[2], score_list[3])
                #for score in score_list:
                    #ostr = ostr + ",% 12.4f" %score
            ostr = ostr + "%3d %2d\n" %(atom_ct, torsion_ct)
            fptr.write(ostr)
    fptr.close()

    

# To execute this command type:
# summarize_wcg_docking.py -d directory -t rmsd tolerance -b report best docking
# only -r receptor filename  -k build hydrogen bonds -v

