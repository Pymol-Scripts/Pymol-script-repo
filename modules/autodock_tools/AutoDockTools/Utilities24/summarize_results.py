#!/usr/bin/env python
#
# 
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/summarize_results.py,v 1.5 2009/06/01 17:35:00 rhuey Exp $
#
import os, glob

from MolKit import Read
from MolKit.hydrogenBondBuilder import HydrogenBondBuilder
from AutoDockTools.Docking import Docking
from mglutil.math.rmsd import RMSDCalculator
from PyAutoDock.InternalEnergy import InternalEnergy
from PyAutoDock.MolecularSystem import MolecularSystem
from PyAutoDock.AutoDockScorer import AutoDock305Scorer

def get_filename(d, conf):
    #find the filename of dlg containing this conformation
    for k in d.dlo_list:
        if conf in k.conformations:
            dlo_filename = os.path.splitext(k.parser.filename)[0]
            break
    return dlo_filename


def get_energy_breakdown(receptor_atoms, ligand_atoms, coords_to_use):
    ms = MolecularSystem()
    ms.add_entities(receptor_atoms)
    ms.add_entities(ligand_atoms)
    ms.set_coords(1, coords_to_use) 
    #ms.set_coords(1, c.getCoords()[:])
    ad_scorer = AutoDock305Scorer()
    ad_scorer.set_molecular_system(ms)
    score_list = ad_scorer.get_score_per_term()
    ostr = ""
    for score in score_list:
        ostr = ostr + "% 8.4f," %score
    return ostr


def construct_hydrogen_bonds(receptor_atoms, ligand_atoms):
    atDict = HydrogenBondBuilder().build(receptor_atoms, ligand_atoms)
    num_hbonds = len(atDict)
    ostr = "%2d," % num_hbonds
    return ostr


def get_best_energy_info(d, rms_tolerance, build_hydrogen_bonds=False, \
                                report_energy_breakdown=False,
                                receptor_filename=None, refCoords=None):
    ostr = "" 
    clust0 = d.clusterer.clustering_dict[rms_tolerance][0]
    c = clust0[0]
    #find the filename of the best result
    dlo_filename = get_filename(d, c)
    if verbose: print "set dlo_filename to ", dlo_filename
    ostr = dlo_filename + ", "
    de = c.docking_energy
    if refCoords:
        ostr += "%3d,%3d,%3d,% 8.4f,% 8.4f," %(len(d.clusterer.data), num_clusters, len(clust0),  de, c.getRMSD(refCoords=refCoords))
    else:
        ostr += "%3d,%3d,%3d,% 8.4f,% 8.4f," %(len(d.clusterer.data), num_clusters, len(clust0),  de, c.getRMSD())
    d.ch.set_conformation(c)
    receptor = None
    if build_hydrogen_bonds:
        if not receptor_filename:
            print "receptor_filename must be specified in order to build_hydrogen_bonds"
            return 
        receptor = Read(receptor_filename)[0] 
        receptor.buildBondsByDistance()
        d.ligMol.buildBondsByDistance()
        ostr += construct_hydrogen_bonds(receptor.allAtoms, d.ligMol.allAtoms)
    if report_energy_breakdown:
        if not receptor_filename:
            print "receptor_filename must be specified in order to build_hydrogen_bonds"
        if not receptor:
            receptor = Read(receptor_filename)[0] 
            receptor.buildBondsByDistance()
        #build bonds if necessary 
        d.ligMol.buildBondsByDistance()
        #ostr += get_energy_breakdown(receptor.allAtoms, d.ligMol.allAtoms, c)
        ostr += get_energy_breakdown(receptor.allAtoms, d.ligMol.allAtoms, c.getCoords()[:])

    return ostr


def get_largest_cluster_info(d, rms_tolerance, build_hydrogen_bonds=False, \
                                    report_energy_breakdown=False,
                                    receptor_filename=None, refCoords=None):
    largest = d.clusterer.clustering_dict[rms_tolerance][0]
    if verbose: print "set largest to ", len(largest)
    for clust in d.clusterer.clustering_dict[rms_tolerance]:
        if verbose: print "current largest cluster len= ", len(clust)
        if len(clust)>len(largest): 
            if verbose: print "resetting largest clust: now len=", len(clust)
            largest = clust
    c = largest[0]   #print info about the lowest energy member of this cluster
    ostr = " "
    dlo_filename = get_filename(d, c)
    ostr += dlo_filename + ","
    de = c.docking_energy
    if refCoords:
        ostr += "%3d,%3d,%3d,% 8.4f,% 8.4f," %(len(d.clusterer.data), num_clusters, len(largest),  de, c.getRMSD(refCoords=refCoords))
    else:
        ostr += "%3d,%3d,%3d,% 8.4f,% 8.4f," %(len(d.clusterer.data), num_clusters, len(largest),  de, c.getRMSD())
    if verbose: print "set dlo_filename to ", dlo_filename
    #update the coords only if you need to?
    d.ch.set_conformation(c)
    receptor = None
    if build_hydrogen_bonds:
        if not receptor_filename:
            print "receptor_filename must be specified in order to build_hydrogen_bonds"
            return 
        receptor = Read(receptor_filename)[0] 
        receptor.buildBondsByDistance()
        ostr += construct_hydrogen_bonds(receptor.allAtoms, d.ligMol.allAtoms)
    if report_energy_breakdown:
        if not receptor_filename:
            print "receptor_filename must be specified in order to build_hydrogen_bonds"
            return 
        if not receptor:
            receptor = Read(receptor_filename)[0] 
            receptor.buildBondsByDistance()
        #ostr += get_energy_breakdown(receptor.allAtoms, d.ligMol.allAtoms, c)
        ostr += get_energy_breakdown(receptor.allAtoms, d.ligMol.allAtoms, c.getCoords()[:]) 
    return ostr


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
        print "        [-f]    rmsd reference filename "
        print "        (default is to use input ligand coordinates from docking log)"
        print "        [-b]    print best docking info only (default is print all)"
        print "        [-L]    print largest cluster info only (default is print all)"
        print "        [-B]    print best docking and largest cluster info only (default is print all)"
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
        opt_list, args = getopt.getopt(sys.argv[1:], 'd:o:t:f:r:bLBkaevh')
    except getopt.GetoptError, msg:
        print 'summarize_results.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    #-d: directory
    directory =  None
    #-t: rms_tolerance
    rms_tolerance =  1.0
    #-f: rms reference
    rms_reference = None
    # optional parameters
    verbose = None
    #-o outputfilename
    outputfilename = "summary_of_results"
    #-r receptor_filename
    receptor_filename = None
    #-b print best only
    print_best_only = False
    #-L print Largest cluster only
    print_Largest_only = False
    #-B print best and Largest cluster only
    print_best_and_Largest_only = False
    #-k build hydrogen bonds 
    build_hydrogen_bonds = False
    #-a append to outputfile
    append_to_outputfile = False
    #-e report_energy_breakdown 
    report_energy_breakdown = False


    #'d:o:t:r:f:bLBkaevh'
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
        if o in ('-f', '--f'):
            rms_reference = a
            if verbose: print 'set rms_reference to ', a
        if o in ('-r', '--r'):
            receptor_filename = a
            if verbose: print 'set receptor_filename to ', a
        if o in ('-b', '--b'):
            print_best_only = True
            if verbose: print 'set print_best_only to ', True
        if o in ('-L', '--L'):
            print_Largest_only = True
            if verbose: print 'set print_Largest_only to ', True
        if o in ('-B', '--B'):
            print_best_and_Largest_only = True
            print_best_only = True
            print_Largest_only = True
            if verbose: print 'set print_best_and_Largest_only to ', True
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
    refCoords = None
    if rms_reference is not None:
        file = directory + '/' + rms_reference
        ref = Read(file)
        if not len(ref): 
            print "unable to read ", rms_reference
        else:
            ref = ref[0]
            coords = ref.allAtoms.coords
            refCoords = coords[:]
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
    first = True
    if append_to_outputfile:
        mode = 'a'
        first = not os.path.exists(outputfilename)
    fptr = open(outputfilename, mode)
    num_clusters = len(d.clusterer.clustering_dict[rms_tolerance])
    #output a header if file opened with mode 'w' or if this is first use of outputfile...
    if verbose: print "first is ", first
    if first:
        tstr = ""
        if print_best_only:
            tstr += "lowestEnergy_dlgfn     #runs #cl #LEC    LE      rmsd_LE "
            if build_hydrogen_bonds:
                tstr+= " #hb "
            if report_energy_breakdown:
                tstr+= " #ESTAT      #HB     #VDW   #DSOLV "
        if print_Largest_only:
            if print_best_only:
                tstr += "   "
            tstr += " largestCl_dlgfn    #runs #cl #LC    LE_LC     rmsd_LC "
            if build_hydrogen_bonds:
                tstr+= " #hb "
            if report_energy_breakdown:
                tstr+= " #ESTAT      #HB     #VDW   #DSOLV "
        if not print_best_only and not print_Largest_only:
            tstr = "#dlgfn                      #in cluster #LE   #rmsd "
        tstr+= "#ats #tors"
        tstr+="\n"
        fptr.write(tstr)
        first = False
    ostr = ""
    if print_best_only or print_Largest_only:
        if print_best_only:
            if verbose: print 'print_best_only is True'
            ostr += get_best_energy_info(d, rms_tolerance, build_hydrogen_bonds, report_energy_breakdown,  receptor_filename, refCoords=refCoords)
        if print_Largest_only:
            if verbose: print 'print_Largest_only is True'
            ostr += get_largest_cluster_info(d, rms_tolerance, build_hydrogen_bonds, report_energy_breakdown, receptor_filename, refCoords=refCoords)
        #add number of atoms and number of torsions to output string
        ostr += "%4d,%2d" %(len(d.ligMol.allAtoms), d.ligMol.parser.keys.count('BRANCH'))
        ostr +='\n'
        try:
            fptr.write(ostr)
        except:
            print 'except'
            fptr.write("\n")
    else:
        receptor = None
        if receptor_filename is not None:
            receptor = Read(receptor_filename)
            receptor = receptor[0]
            receptor.buildBondsByDistance()
        for clust in d.clusterer.clustering_dict[rms_tolerance]:
            if verbose:
                print "len(d.clusterer.clustering_dict)=", len(d.clusterer.clustering_dict[rms_tolerance])
                print "len(clust)=", len(clust)
            #for each cluster in the clustering at this rms_tolerance:
            #output length of cluster and energies of the lowest energy result
            c = clust[0]
            if verbose: print "setting c to clust[0]", id(c)
            #write the dlg filename for the best in this cluster
            dlo_filename = get_filename(d, c)
            ostr = dlo_filename + ", " 
            fptr.write(ostr)
            if build_hydrogen_bonds or report_energy_breakdown:
                if verbose: print "setting conf to ", id(c)
                d.ch.set_conformation(c)
            de = c.docking_energy
            ostr = "%3d,% 8.4f,% 8.4f" %(len(clust), de, c.getRMSD(refCoords=refCoords))
            if build_hydrogen_bonds and receptor:
                ostr += construct_hydrogen_bonds(receptor.allAtoms, d.ligMol.allAtoms)
            if report_energy_breakdown and receptor:
                #ostr += get_energy_breakdown(receptor.allAtoms, d.ligMol.allAtoms, c)
                ostr += get_energy_breakdown(receptor.allAtoms, d.ligMol.allAtoms, c.getCoords()[:])
            ostr += "%4d,%3d" %(len(d.ligMol.allAtoms), d.ligMol.parser.keys.count('BRANCH'))
            ostr = ostr + "\n"
            fptr.write(ostr)
    fptr.close()
    

# To execute this command type:
# summarize_results.py -d directory -t rmsd tolerance -r receptor filename  
#                         -b report best_docking_only 
#                         -k build_hydrogen_bonds 
#                         -f rms reference filename
#                         -L report largest cluster information 
#                         -B print both best docking info  and largest cluster info
#                         -e report_energy_breakdown 
#                         -v verbose
# NOTE: -d, -t, -r and -f require arguments whereas the others set booleans

