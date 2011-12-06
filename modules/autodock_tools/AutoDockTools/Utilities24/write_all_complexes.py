#!/usr/bin/env python
#
# 
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/write_all_complexes.py,v 1.2 2007/10/08 18:14:45 rhuey Exp $
#
# $Id: write_all_complexes.py,v 1.2 2007/10/08 18:14:45 rhuey Exp $
#
import os

from MolKit import Read
from AutoDockTools.Docking import Docking


if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: write_all_complexes.py -d directory"
        print
        print "    Description of command..."
        print "         -d     docking filename"
        print "         -r     receptor filename"
        print "    Optional parameters:"
        print "        [-b]    write complex of best result and receptor only (default is write all)"
        print "        [-o]    output file stem"
        print "        (default is receptorname_ligandname plus '_confnum.pdbqt')"
        print "        [-v]    verbose output"


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'd:r:bo:vh')
    except getopt.GetoptError, msg:
        print 'write_all_complexes.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    #-d: docking log filename
    docking_filename =  None
    #-r receptor_filename
    receptor_filename = None

    # initialize optional parameter
    #-b print best only
    write_best_only = False
    #-o output_stem
    output_stem = None
    #-v verbose best only
    verbose = False

    #'d:r:bvh'
    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('-d', '--d'):
            docking_filename = a
            if verbose: print 'set docking_filename to ', a
        if o in ('-r', '--r'):
            receptor_filename = a
            if verbose: print 'set receptor_filename to ', a
        if o in ('-b', '--b'):
            write_best_only = True
            if verbose: print 'set write_best_only to ', True
        if o in ('-o', '--o'):
            output_stem = a
            if verbose: print 'set output_stem to ', a
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print 'set verbose to ', True
        if o in ('-h', '--'):
            usage()
            sys.exit()

    if not docking_filename:
        print 'write_all_complexes: docking_filename must be specified.'
        usage()
        sys.exit()

    if not receptor_filename:
        print 'write_all_complexes: receptor_filename must be specified.'
        usage()
        sys.exit()


    d = Docking()
    d.readDlg(docking_filename)
    receptor_stem = os.path.splitext(os.path.basename(receptor_filename))[0]
    print "receptor_stem =", receptor_stem

    rptr = open(receptor_filename)
    receptor_lines = rptr.readlines()
    rptr.close()
    
    if write_best_only:
        if not hasattr(d, 'clusterer'):
            print docking_filename, ' has no clustering... unable to write complex of best result'
            sys.exit()
        rms = d.clusterer.clustering_dict.keys()[0]
        confs = [d.clusterer.clustering_dict[rms][0][0]]
    else:
        confs = d.ch.conformations
        
    if output_stem is None:
        output_stem = receptor_stem + '_' + d.ligMol.name
    for ix, conf in enumerate(confs):
        d.ch.set_conformation(conf)
        if write_best_only:
            ind = 0
        else:
            ind = d.ch.conformations.index(conf)
        ext = '.pdbq'
        if d.version==4.0:
            ext = '.pdbqt'
        outputfilename = output_stem + '_' + str(ix) + ext
        fptr = open(outputfilename, 'w')
        for l in receptor_lines:
            fptr.write(l)
        ctr = 0
        for l in d.ligMol.parser.allLines:
            if l.find("ATOM")==0 or l.find("HETATM")==0:
                crds = d.ligMol.allAtoms[ctr].coords
                rec = "%s%8.3f%8.3f%8.3f%s\n"%(l[:30],crds[0], crds[1], crds[2],l[54:] ) 
                fptr.write(rec)
                ctr += 1
        fptr.close()
        if verbose:
            print "wrote ", outputfilename


# To execute this command type:
# write_all_complexes.py -d docking -r receptor filename 
# optional arguments
# -b write_complex for best docking result only 
# -o output_stem (default is receptorname_ligandname)
