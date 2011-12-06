#!/usr/bin/env python

import os
from string import find, strip

from MolKit import Read
from MolKit.pdbWriter import PdbqtWriter
from AutoDockTools.atomTypeTools import AutoDock4_AtomTyper
from AutoDockTools.DockingParameters import DockingParameters
from AutoDockTools.DockingParameters import genetic_algorithm_list4
from AutoDockTools.DockingParameters import genetic_algorithm_local_search_list4
from AutoDockTools.DockingParameters import local_search_list4
from AutoDockTools.DockingParameters import simulated_annealing_list4




if __name__ == '__main__':
    import sys
    import getopt

    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: dpf3_to_dpf4.py -s dpf3_stem -r receptor_filename -l ligand_filename"
        print
        print "    Description of command..."
        print "        [-s]    dpf3_stem"
        print "        [-r]    receptor_filename (pdbqt)"
        print "        [-l]    ligand_filename (pdbqt)"
        print "    Optional parameters:"
        print "        [-o]    outputfilename"
        print "        [-p]    parameter_list "
        print "                (default is 'genetic_algorithm_local_search_list4')"
        print "        [-v]    verbose output"


    # process command arguments
    try:
        #opt_list, args = getopt.getopt(sys.argv[1:], 'r:o:bs:ampvHG')
        opt_list, args = getopt.getopt(sys.argv[1:], 's:r:l:o:p:v')

    except getopt.GetoptError, msg:
        print 'dpf3_to_dpf4.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    #-s: pdbq_stem
    dpf3_filename =  None

    #-r: receptor_filename
    receptor_filename =  None

    #-l: ligand_filename
    ligand_filename =  None

    # optional parameters

    #-p: parameter_list
    parameter_list =  'genetic_algorithm_local_search_list4'
    verbose = None
    outputfilename = None

    #'s:r:l:o:p:v'
    for o, a in opt_list:
        if o in ('-s', '--s'):
            dpf3_stem = a
            dpf3_filename = a + '.dpf'
            if verbose: print 'set dpf3_filename to ', dpf3_filename
        if o in ('-r', '--r'):
            receptor_filename = a
            if verbose: print 'set receptor_filename to ', receptor_filename
            extension = os.path.splitext(receptor_filename)[1]
            if extension!='.pdbqt':
                print receptor_filename, " not in pdbqt format"
                usage()
                sys.exit()
        if o in ('-l', '--l'):
            ligand_filename = a
            if verbose: print 'set ligand_filename to ', ligand_filename
            #check that the ligand_filename type is pdbqt
            extension = os.path.splitext(ligand_filename)[1]
            if extension!='.pdbqt':
                print ligand_filename, " not in pdbqt format"
                usage()
                sys.exit()
        if o in ('-p', '--p'):
            parameter_list = a
            if verbose: print 'set parameter_list to ', parameter_list
        if o in ('-o', '--o'):
            outputfilename = a
            if verbose: print 'set outputfilename to ', outputfilename
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print 'set verbose to ', True
        if o in ('-h', '--'):
            usage()
            sys.exit()


    if not ligand_filename:
        print 'dpf3_to_dpf4: ligand_filename must be specified.'
        usage()
        sys.exit()
    if not receptor_filename:
        receptor_filename = dpf3_stem + '.pdbqt'
        print "dpf3_to_dpf4: using default receptor_filename:", receptor_filename
    if not dpf3_filename:
        print 'dpf3_to_dpf4: no dpf3_filename specified using GridParameter defaults'

    parms= {}
    parms['simulated_annealing_list4'] =  simulated_annealing_list4
    parms['genetic_algorithm_list4'] = genetic_algorithm_list4
    parms['genetic_algorithm_local_search_list4'] = genetic_algorithm_local_search_list4
    parms['local_search_list4'] = local_search_list4 

    #what about nucleic acids???
    DPO4 = DockingParameters()
    print "default: ligand_types=", DPO4['ligand_types']['value']
    orig_compute_unbound_extended_flag = DPO4['compute_unbound_extended_flag']['value']
    if dpf3_filename:
        DPO4.read(dpf3_filename)
        if verbose: print 'read ', dpf3_filename
        DPO4['compute_unbound_extended_flag']['value'] = orig_compute_unbound_extended_flag
    else:
        if verbose: print 'using dpo defaults'
    DPO4.set_ligand(ligand_filename)
    DPO4.set_version("4.1")
    DPO4['set_sw1']['value'] = 0
    DPO4['set_sw1_flag']['value'] = 0
    DPO4['set_psw1']['value'] = 1
    DPO4['set_psw1_flag']['value'] = 1
    res = DPO4.set_ligand_types_from_filename(ligand_filename)
    if res=="ERROR": raise 'unreadable file', ligand_filename
    if verbose: 
        print "ligand_types=", DPO4['ligand_types']['value']
        print 'set ligand to  ', ligand_filename
    DPO4.set_receptor(receptor_filename)
    if verbose: print 'set receptor to  ', receptor_filename
    if outputfilename is None:
        outputfilename = dpf3_stem + '_4.dpf'

    DPO4.write4(outputfilename, parms[parameter_list])
    if verbose: 
        print "wrote ", outputfilename, ' using:'
        for p in parameter_list:
            print p
    

# To execute this command type:
# dpf3_to_dpf4.py -s dpf3_stem -r receptor_filename -l ligand_filename -v

