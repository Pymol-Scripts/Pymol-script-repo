#!/usr/bin/env python
#
# $Id: summarize_docking_directory.py,v 1.5 2007/10/09 17:30:07 annao Exp $
#

import logging as log
import os, glob
import numpy.oldnumeric as Numeric
import string

from MolKit import Read
from AutoDockTools.Docking import Docking

if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        """Print helpful, accurate usage statement to stdout."""
        print "Usage: summarize_docking_directory.py -d directory"
        print
        print "    Write docking result summary to stdout."
        print
        print "    Required parameters:"
        print "         -d directory    docking directory"
        print "    Optional parameters:"
        print "        [-l filename]    reference_ligand_filename"
        print "        [-t rmsd]        rmsd tolerance (default is 1.0)"
        print "        [-f]             write to file (name derived from docking)"
        print "        [-o filename]    write to filename"
        print "        [-v]             verbose output"
        print "        [-D]             debugging output"


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'd:l:fo:t:vhD')
    except getopt.GetoptError, msg:
        print 'summarize_docking_directory.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    #-d: directory
    directory =  None

    # optional parameters

    #-l: ligand_filename
    ligand_filename =  None
    #-t: rms_tolerance
    rms_tolerance =  1.0
    verbose = None
    #-f
    write_to_file = False
    #-o output_filename
    output_filename = None

    for o, a in opt_list:
        if o in ('-d', '--d'):
            directory = a
            log.info("Docking directory: %s" % directory)
        if o in ('-l', '--l'):
            ligand_filename = a
            log.info('set reference ligand_filename to ', ligand_filename)
        if o in ('-o', '--o'):
            output_filename = a
            log.info('set output_filename to ', output_filename)
        if o in ('-t', '--t'):
            rms_tolerance = float(a)
            log.info("Clustering rmsd tolerance %5.2f" % rms_tolerance)
        if o in ('-f', '--f'):
            write_to_file = True
        if o in ('-v', '--v'):
            log.basicConfig(level=log.INFO)
            log.info("verbose output")
        if o in ('-D', '--DEBUG'):
            log.basicConfig(level=log.DEBUG)
            log.debug("debugging output")
        if o in ('-h', '--'):
            usage()
            sys.exit()


    # make sure directory was specified
    if not directory:
        log.error('summarize_docking_directory: docking directory must be specified.')
        usage()
        sys.exit()


    # make sure directory exists
    try:
        os.stat(directory)
    except OSError, msg:
        log.error("summarize_docking_directory: %s: %s",  msg.strerror, msg.filename)
        sys.exit()


    # get absolute, normalized directory_path
    directory_pathname = os.path.abspath(directory) + "/"
    docking_pathname = os.path.dirname(directory_pathname)
    docking_name = string.split(docking_pathname, '/')[-1]

    log.debug("directory: %s" % directory)
    log.debug("directory_pathname: %s" % directory_pathname)
    log.debug("docking_pathname: %s" % docking_pathname)
    log.debug("docking_name: %s" % docking_name)


    output_fileptr = sys.stdout
    if output_filename:
        output_filepath = os.path.join(docking_pathname, output_filename)
        output_fileptr = open(output_filepath, 'w')
        log.info("Writing to %s" % output_filepath)
    else:
        if write_to_file:
            # create output filename from directory and rms_tolerance
            output_filepath = os.path.join(docking_pathname, docking_name) + \
                              '-' + str(rms_tolerance)+ ".csv"
            output_fileptr = open(output_filepath, 'w')
            log.info("Writing to %s" % output_filepath)
            
    # use the directory name to deduce dpf and gpf names
    # get all ligand and receptors filename from dpf and gpf
    
    # @@ ASSUMPTIONS: WCG conventions for dpf and gpf names @@
    # @@ from docking directory name                        @@
    dpf_pathname = os.path.join(docking_pathname, docking_name) + ".dpf"
    dpf_filename = os.path.basename(dpf_pathname)
    # remove "faahNNNN_ from docking name
    gpf_pathname = os.path.join(docking_pathname, docking_name[9:]) + ".gpf"
    gpf_filename = os.path.basename(gpf_pathname)

    #
    # get filenames from the dpf and gpf
    #
    from AutoDockTools.DockingParameters import DockingParameters
    from AutoDockTools.GridParameters import GridParameters

    dpo = DockingParameters()
    log.info("reading dpf: %s" % dpf_pathname)
    dpo.read(dpf_pathname)
    ligand_filename = dpo["move"]["value"]
    ligand_pathname = os.path.join(docking_pathname, ligand_filename)
    log.info ("ligand_filename: %s" % ligand_filename)
    log.debug("ligand_pathname: %s" % ligand_pathname)

    gpo = GridParameters()
    log.info("reading gpf: %s" % gpf_pathname)
    gpo.read(gpf_pathname)
    receptor_filename = gpo["receptor"]["value"]
    receptor_pathname = os.path.join(docking_pathname, receptor_filename)
    log.info ("receptor_filename: %s" % receptor_filename)
    log.debug("receptor_pathname: %s" % receptor_pathname)

    # docking log pathnames
    dlg_pathname_list = glob.glob(os.path.join( docking_pathname, '*.dlg'))

    # At this point, we have pathnames and filenames for all the
    # input and output files, start reading them
    ######################################################################
   
    # the ligand from the dpf
    dpf_ligand = Read(ligand_pathname)[0]

    #read all the docking logs in as one Docking
    d = Docking()
    for dlg in dlg_pathname_list:
        log.debug("reading dlg: %s" % dlg)
##         d.readDlg(dlg, ligand=dpf_ligand)
        d.readDlg(dlg)

    # set up the clusterer to use the HIV PR c2 specific getRMSD_custom method
    d.clusterer.set_get_distance( d.clusterer._get_distance_custom)

    log.info("clustering %s at rms_tol %f" % (docking_name, rms_tolerance))
    d.clusterer.make_clustering(rms_tolerance)
##     d.clusterer.clustering_dict[rms_tolerance].do_stats()
##     d.clusterer.write_summary()


    # create an output dictionary from the clustering

    ol = [] # output list
    seed = d.clusterer.clustering_dict[rms_tolerance][0][0]
    for cx, cluster in enumerate(d.clusterer.clustering_dict[rms_tolerance]):
        od = {} # output dictionary for this cluster
        conf = cluster[0]

        # general info about docking
        od['num_conf_docking'] = "%7d" % (len(d.clusterer.data))
        od['ntors'] = "%7d" % (len(conf.torsions))
        od['natoms'] = "%7d" % (len(conf.coords))

        # general info about this clustering
        od['num_clusters'] = "%7d" % (len(d.clusterer.clustering_dict[rms_tolerance]))
        od['rmstol'] = "%7.3f" % (rms_tolerance)
        
        # specific info about this cluster
        od['num_conf_cluster'] = " %7d" % (len(cluster))
        od['energy'] = "%8.3f" % (conf.binding_energy)
        od['clust_rmsd'] = "%7.3f" % (d.clusterer.get_distance(seed, conf))
        od['cluster_ix'] = "%7d" % (cx+1)

        # search for the seed conformation in the docking logs
        od['dlg'] = "%32s" % (docking_name) # default to this if search (below) fails
        for dlo in d.dlo_list:
            if conf in dlo.conformations:
                od['dlg'] = "%32s" % \
                            (os.path.splitext(os.path.basename(dlo.parser.filename))[0])

        log.debug(od)
        ol.append(od)

    log.debug(od['energy'])
    log.debug(od['clust_rmsd'])
    log.debug(type(od['energy']))
    log.debug(type(od['clust_rmsd']))


    # the ol list of od's should already be sorted by cluster energy
    pass
    
    # initialize output header dictionary (same keys as od)
    oh = {}
    for key in od.keys():
        oh[key] = "<%s header>" % key # init with generic values
    # now put in less-generic? values
    oh['num_conf_docking'] = "%7s" % ( "#conf")
    oh['ntors']            = "%7s" % ( "#tors")
    oh['natoms']           = "%7s" % ( "#ats")
    oh['num_clusters']     = "%7s" % ( "#clust")
    oh['rmstol']           = "%7s" % ( "rmstol")
    oh['num_conf_cluster'] = "%8s" % ( "Csize")
    oh['energy']           = "%8s" % ( "LE")
    oh['clust_rmsd']       = "%7s" % ( "rmsd")
    oh['cluster_ix']       = "%7s" % ( "Crank")
    oh['dlg']              = "%32s" % ( "docking_filename")

    # write the csv
    import csv
    fields = ['dlg',
              'energy',
              'clust_rmsd',
              'num_conf_cluster',
              'cluster_ix',
              'num_clusters',
              'rmstol',
              'num_conf_docking',
              'natoms',
              'ntors']

    writer = csv.DictWriter( output_fileptr, fieldnames=fields, restval="@@")
    writer.writerow(oh)
    writer.writerows(ol)

    log.debug("exiting...")
    sys.exit()


    
    

# To execute this command type:
# summarize_with_rotated_rmsd.py -d directory -t rmsd tolerance -b report best docking
# only -r receptor filename  -k build hydrogen bonds -v

