#!/usr/bin/env python
#
# 
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/summarize_epdb_results4.py,v 1.2 2010/02/03 00:08:25 rhuey Exp $
#
# $Id: summarize_epdb_results4.py,v 1.2 2010/02/03 00:08:25 rhuey Exp $
#
import os, glob
from numpy import oldnumeric as Numeric
from AutoDockTools.EpdbParser import EpdbParser



if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: summarize_epdb_results4.py -f epdb_dlgfilename"
        print
        print "    Description of command..."
        print "         -f     epdb_dlgfilename"
        print "    Optional parameters:"
        print "        [-o]    output filename"
        print "                      (default is 'summary_of_epdb_results')"
        print "        [-a]    append to  output filename"
        print "                      (default is to open output filename 'w')"
        print "        [-v]    verbose output"


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'f:o:avh')
    except getopt.GetoptError, msg:
        print 'summarize_epdb_results4.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    #-f: epdb_dlgfilename
    epdb_dlgfilename =  None
    # optional parameters
    #-o outputfilename
    outputfilename = "summary_of_epdb_results"
    #-a append to outputfile
    append_to_outputfile = False
    #-v verbose
    verbose = False

    #'f:o:avh'
    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('-f', '--f'):
            epdb_dlgfilename = a
            if verbose: print 'set epdb_dlgfilename to ', a
        if o in ('-o', '--o'):
            outputfilename = a
            if verbose: print 'set outputfilename to ', a
        if o in ('-a', '--a'):
            append_to_outputfile = True
            if verbose: print 'set append_to_outputfile to ', True
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print 'set verbose to ', True
        if o in ('-h', '--'):
            usage()
            sys.exit()

    if not  epdb_dlgfilename:
        print 'summarize_epdb_results4: epdb_dlgfilename must be specified.'
        usage()
        sys.exit()

    ep = EpdbParser(epdb_dlgfilename)
    mode = 'w'
    first = True
    if append_to_outputfile:
        mode = 'a'
        first = not os.path.exists(outputfilename)
    fptr = open(outputfilename, mode)

    if verbose: print "first is ", first
    if first:
        tstr = "      dlgfn      estFE  IntermolE   InternalE     torsFE   sumTotal   sumEstat     sumVDW  #at #tors\n "
        fptr.write(tstr)
    ostr = "%10s,%10.2f,%10.2f,%11.2f,%10.2f,%10.2f,%10.2f,%10.2f,%4d,%5d\n"%(os.path.basename(epdb_dlgfilename), \
                                                                          ep.estFreeEnergy, ep.finalIntermolEnergy, \
                                                                          ep.finalTotalInternalEnergy, ep.torsionalFreeEnergy, \
                                                                          Numeric.add.reduce(ep.total_energies), \
                                                                          Numeric.add.reduce(ep.estat_energies), \
                                                                          Numeric.add.reduce(ep.vdw_energies), ep.atmCtr, ep.ntors)
    fptr.write(ostr)
    fptr.close()
    
    

# To execute this command type:
# summarize_epdb_results4.py -d directory -t rmsd tolerance 
# -f epdb_dlgfilename
# -o outputfilename
# -a append to outputfilename
# -v verbose output
