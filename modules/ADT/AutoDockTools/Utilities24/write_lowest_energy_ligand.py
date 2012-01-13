#!/usr/bin/env python
#
# 
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/write_lowest_energy_ligand.py,v 1.2 2007/10/08 18:15:37 rhuey Exp $
#
import os 

from MolKit import Read
from AutoDockTools.Docking import Docking




if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: write_lowest_energy_ligand.py -f dlgfilename"
        print
        print "    Description of command..."
        print "         -f     dlgfilename"
        print "    Optional parameters:"
        print "        [-v]    verbose output"
        print "        [-o pdbqt_filename] (output filename)"

    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'f:vo:h')
    except getopt.GetoptError, msg:
        print 'write_lowest_energy_ligand.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    #-f: dlgfilename
    dlgfilename =  None
    # optional parameters
    verbose = None
    #-o outputfilename
    outputfilename = None

    #'f:vo:h'
    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('-f', '--f'):
            dlgfilename = a
            if verbose: print 'set dlgfilename to ', a
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print 'set verbose to ', True
        if o in ('-o', '--o'):
            outputfilename = a
            if verbose: print 'set outputfilename to ', a
        if o in ('-h', '--'):
            usage()
            sys.exit()


    if not  dlgfilename:
        print 'write_lowest_energy_ligand: dlgfilename must be specified.'
        usage()
        sys.exit()

    #what about nucleic acids???

    d = Docking()
    d.readDlg(dlgfilename)

    if verbose: print 'read ', dlgfilename
    key0 = d.clusterer.clustering_dict.keys()[0]
    conf0 = d.clusterer.clustering_dict[key0][0][0]
    d.ch.set_conformation(conf0)
    parser = d.ligMol.parser
    lines = []
    #have to add newline character to lines read from dlg
    for l in parser.allLines:
        l+= '\n'
        lines.append(l)
    parser.allLines = lines
    coords = d.ligMol.allAtoms.coords
    if outputfilename is None:
        outputfilename = d.ligMol.name  + '_BE.pdbqt'
    parser.write_with_new_coords(coords, outputfilename) 
    if verbose:
        print 'wrote %s' %outputfilename


# To execute this command type:
# write_lowest_energy_ligand.py -f dlgfilename [-o outputfilename] -v




