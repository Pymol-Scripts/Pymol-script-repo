#!/usr/bin/env python
#
# 
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/write_rigid_ligand.py,v 1.1 2008/05/28 16:48:41 rhuey Exp $
#
import os, string
from MolKit import Read




if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        print "Usage: write_rigid_ligand.py "
        print "     This writes a rigid ligandfile with no rotatable bonds"
        print "        -l ligand_filename"
        print "    Optional parameters:"
        print "        [-o pdbqt_filename] (default ligandstem_rigid.pdbqt)"
        print "        [-v]    verbose output"

    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'l:o:vh')
    except getopt.GetoptError, msg:
        print 'write_rigid_ligand.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    # optional parameters
    verbose = None
    #-l ligandfilename
    ligandfilename = None
    #-o outputfilename
    outputfilename = None

    #'l:o:vh'
    for o, a in opt_list:
        if o in ('-l', '--l'):
            ligandfilename = a
            if verbose: print 'set ligandfilename to ', a
        if o in ('-o', '--o'):
            outputfilename = a
            if verbose: print 'set outputfilename to ', a
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print 'set verbose to ', True
        if o in ('-h', '--'):
            usage()
            sys.exit()
    if not ligandfilename:
        print "write_rigid_ligand.py: ligandfilename must be specified."
        usage()
        sys.exit()

    lig = Read(ligandfilename)[0]
    ext = os.path.splitext(ligandfilename)[1]
    if outputfilename is None:
        outputfilename = lig.name + "_rigid" + ext
    if verbose: print 'set outputfilename to ', outputfilename
    parser = lig.parser
    optr = open(outputfilename, 'w')
    #have to add newline character to lines read from dlg
    for l in parser.allLines:
        if l.find("active torsions:")>-1:
            nl = "REMARK 0 active torsions:\n"
            optr.write(nl)
        elif l.find("A    between atoms:")>-1:
            ll = l.split()
            #REMARK,1,A,between,atoms:,N1_1517,and,C31_1555"
            nl = ll[0] + "       I    " + string.join(ll[3:]) + "\n"
            optr.write(nl)
        elif l.find("TORSDOF")==0:
            nl = "ENDROOT\n"
            optr.write(nl)
            optr.write(l)
        elif l.find("BRANCH")==-1 and l.find("ENDBRANCH")==-1 and l.find("ENDROOT")==-1:
            optr.write(l)
    if verbose:
        print 'wrote %s' %outputfilename


# To execute this command type:
# write_rigid_ligand.py  -l ligandfilename [-o outputfilename] -v




