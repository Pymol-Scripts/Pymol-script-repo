#!/usr/bin/env python
#
# 
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/summarize_time.py,v 1.3 2007/10/08 18:14:20 rhuey Exp $
#
import os, glob, numpy.oldnumeric as Numeric
from AutoDockTools.DlgParser import DlgParser


if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: summarize_time.py -d directory"
        print
        print "    Description of command..."
        print "         -d     directory"
        print "    Optional parameters:"
        print "        [-o]    output filename"
        print "                      (default is 'summary_of_time')"
        print "        [-v]    verbose output"


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'd:o:vh')
    except getopt.GetoptError, msg:
        print 'summarize_time.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    #-d: directory
    directory =  None

    # optional parameters
    #-o outputfilename
    outputfilename = "summary_of_time.txt"
    #-v: verbose
    verbose = None
    #-h: help
    help = None


    #'d:o:vh'
    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('-d', '--d'):
            directory = a
            if verbose: print 'set directory to ', a
        if o in ('-o', '--o'):
            outputfilename = a
            if verbose: print 'set outputfilename to ', a
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print 'set verbose to ', True
        if o in ('-h', '--'):
            usage()
            sys.exit()

    if not  directory:
        print 'summarize_time: directory must be specified.'
        usage()
        sys.exit()

    p = DlgParser()
    dlg_list = glob.glob(directory + '/*.dlg')
    for dlg in dlg_list:
        p.parse(dlg)
    stem = dlg_list[0].split('.')[0]
    ostr = '%s: % 12.4f m \n' %(directory, p.total_time) 
    fptr = open(outputfilename, 'a')
    fptr.write(ostr)
    fptr.close()

# To execute this command type:
# summarize_time.py -d directory -o outputfilename -v

