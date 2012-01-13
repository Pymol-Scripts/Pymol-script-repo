#!/usr/bin/env python
#
# 
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/process_VinaResult.py,v 1.5 2010/11/09 00:09:24 rhuey Exp $
#
# $Id: process_VinaResult.py,v 1.5 2010/11/09 00:09:24 rhuey Exp $
#
import os, glob, time

from MolKit import Read
from AutoDockTools.Docking import Docking, VinaResultProcessor
from mglutil.math.rmsd import RMSDCalculator
from AutoDockTools.InteractionDetector import InteractionDetector




if __name__ == '__main__':
    import sys
    import getopt

    #'d:o:t:f:r:BLl:x:Dnc vh'

    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: process_VinaResult.py -f ADVina_resultfile -r receptor_filename"
        print
        print "    Description of command..."
        print "         -f     vina_resultfile"
        print "         -r     receptor filename "
        print "    Optional parameters:"
        print "        [-o]    outputfilestem (default is None yields resultfilenames 'ind_model1.pdbqt')"
        print "        [-m]    remove_model_str (remove 'model' from resultfilenames,eg 'ind_1.pdbqt', 'ind_2.pdbqt'...)"
        print "        [-B]    process best result only 'ind_model1.pdbqt' (default is to output a file for each MODEL) "
        print "        [-x]    maximum number of models to process (default is to process all models)"
        print "        [-c]    do not detect close contacts (default is to detect them)"
        print "        [-p]    report pi-pi interactions (default is NOT to detect pi-pi or pi-cation)"
        print "        [-v]    verbose output"


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'f:r:o:Bx:mcpvh')
    except getopt.GetoptError, msg:
        print 'process_VinaResult.py: %s' %msg
        usage()
        sys.exit(2)
    
    # initialize required parameters
    #-f: vina_resultfile
    vina_resultfile =  None
    #-r receptor_filename 
    receptor_filename = None

    # optional parameters
    #-o: outputfilestem
    outputfilestem =  None
    #-m: remove_model_str
    remove_model_str =  False
    #-B create outputfile for best model only
    best_only = False
    #-c do NOT detect_close_contacts
    detect_close_contacts = True
    #-p include detection of pi-pi and pi-cation interactions 
    detect_pi = False
    #-x maximum number of models to process
    max_models = -1  #use ?None? instead  

    # -v verbose output
    verbose = None


    #'f:r:o:Bcpx vh'
    for o, a in opt_list:
        #required parameters
        if o in ('-f', '--f'):
            vina_resultfile = a
            if verbose: print 'set vina_resultfile to ', a
        if o in ('-r', '--r'):
            receptor_filename = a
            if verbose: print 'set receptor_filename to ', receptor_filename
        #optional parameters
        if o in ('-o', '--o'):
            outputfilestem = a
            if verbose: print 'set outputfilestem to ', outputfilestem
        if o in ('-m', '--m'):
            remove_model_str = True
            if verbose: print 'set remove_model_str to ', remove_model_str
        if o in ('-B', '--B'):
            best_only = True
            if verbose: print 'best_only is ', best_only
        if o in ('-x', '--x'):
            max_models = int(a) 
            if verbose: print 'At most '+ max_models + ' will be created'
        if o in ('-c', '--c'):
            detect_close_contacts = False
            if verbose: print 'set detect_close_contacts to ', detect_close_contacts
        if o in ('-p', '--p'):
            detect_pi = True
            if verbose: print 'set detect_pi interactions to ', detect_pi
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print 'set verbose to ', True
        if o in ('-h', '--'):
            usage()
            sys.exit()

    if not  vina_resultfile:
        print 'process_VinaResult: vina_resultfile must be specified.'
        usage()
        sys.exit()

    if not  receptor_filename:
        print 'process_VinaResult: receptor_filename must be specified.'
        usage()
        sys.exit()


    #read the vina result
    ligMols = Read(vina_resultfile)
    if not len(ligMols):
        print "unable to read vina resultfile ", vina_resultfile
        sys.exit()

    #check that receptor_filename exists
    try:
        assert os.path.exists(receptor_filename)
    except:
        print "Sorry unable to find receptor ", receptor_filename
        sys.exit()


    vrp = VinaResultProcessor( receptor_filename=receptor_filename,
                                 detect_pi=detect_pi, 
                                 detect_close_contacts=detect_close_contacts) 
    vrp.process( vina_resultfile, outputfilestem=outputfilestem, 
                                 remove_model_str=remove_model_str, 
                                 best_only=best_only) 


# process_VinaResult.py 
#       -f vina_resultfile 
#       -r receptor_filename
#optional:
#-o    outputfilestem (default is None resulting in filenames 'ind_model1.pdbqt','ind_model2.pdbqt' ...)"
#-m    remove_model_str (default is False resulting in filenames 'ind_model1.pdbqt','ind_model2.pdbqt' ...)"
#                       (if True results in filenames 'ind_1.pdbqt','ind_2.pdbqt' ...)"
#-B    process best result only (default is create a pdbqt+ file for each MODEL)"
#-x    maximum number of MODELS to process (default is to process all)"
#-c    do not detect atoms in close contact (default is TO detect them)"
#-p    report pi-pi interactions (default is NOT to detect pi-pi or pi-cation)"
#-v    verbose output
#-h    print usage statement

