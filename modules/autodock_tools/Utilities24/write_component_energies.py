#!/usr/bin/env python
#
# 
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/write_component_energies.py,v 1.2 2007/10/08 18:15:14 rhuey Exp $
#
import os 

from MolKit import Read
from PyAutoDock.InternalEnergy import InternalEnergy
from PyAutoDock.MolecularSystem import MolecularSystem
from PyAutoDock.AutoDockScorer import AutoDock4Scorer




if __name__ == '__main__':
    import sys
    import getopt


# write_component_energies.py -r receptorfilename -l ligandfilename [-o outputfilename -a append_to_file] -v
    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: write_component_energies.py -r receptorfilename -l ligandfilename"
        print
        print "    Description of command..."
        print "         -r     receptorfilename"
        print "         -l     ligandfilename"
        print "    Optional parameters:"
        print "        [-o]    outputfilename (default is 'component_energies.txt')"
        print "        [-a]    append_to_file (default mode is 'w')"
        print "        [-v]    verbose output"

    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'r:l:avo:h')
    except getopt.GetoptError, msg:
        print 'write_component_energies.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    #-r: receptorfilename
    receptorfilename =  None
    #-l: ligandfilename
    ligandfilename =  None
    #-a append to outputfile
    append_to_outputfile = False
    # optional parameters
    verbose = None
    #-o outputfilename
    outputfilename = 'component_energies.txt'

    #'r:l:avo:h'
    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('-r', '--r'):
            receptorfilename = a
            if verbose: print 'set receptorfilename to ', a
        if o in ('-l', '--l'):
            ligandfilename = a
            if verbose: print 'set ligandfilename to ', a
        if o in ('-a', '--a'):
            append_to_outputfile = True
            if verbose: print 'set append_to_outputfile to ', True
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print 'set verbose to ', True
        if o in ('-o', '--o'):
            outputfilename = a
            if verbose: print 'set outputfilename to ', a
        if o in ('-h', '--'):
            usage()
            sys.exit()


    if not receptorfilename:
        print 'write_component_energies: receptorfilename must be specified.'
        usage()
        sys.exit()

    if not ligandfilename:
        print 'write_component_energies: ligandfilename must be specified.'
        usage()
        sys.exit()

    receptor = Read(receptorfilename)[0]
    if verbose: print 'read ', receptorfilename
    receptor.buildBondsByDistance()

    ligand = Read(ligandfilename)[0]
    if verbose: print 'read ', ligandfilename
    ligand.buildBondsByDistance()
    mode = 'w'
    first = True
    if append_to_outputfile:
        mode = 'a'
        first = not os.path.exists(outputfilename)
    if verbose: print 'first is ', first
    optr = open(outputfilename, mode)
    if first:
        tstr = "RECEPTOR_LIGAND         #ESTAT       #HB      #VDW    #DSOLV \n"
        optr.write(tstr)
    #get the scores
    ms = MolecularSystem()
    ms.add_entities(receptor.allAtoms)
    ms.add_entities(ligand.allAtoms)
    ad_scorer = AutoDock4Scorer()
    ad_scorer.set_molecular_system(ms)
    score_list = ad_scorer.get_score_per_term()
    ostr = "%8s_%6s" %(receptor.name,ligand.name)
    for score in score_list:
        ostr = ostr + " % 8.4f," %score
    ostr = ostr + "\n"   
    optr.write(ostr)
    if verbose:
        print 'wrote %s' %outputfilename


# To execute this command type:
# write_component_energies.py -r receptorfilename -l ligandfilename [-o outputfilename -a append_to_file] -v




