#!/usr/bin/env python
#
# 
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/compute_AutoDock41_score.py,v 1.10 2007/10/10 22:59:26 rhuey Exp $
#
import os 

from MolKit import Read
import MolKit.molecule
import MolKit.protein


from PyAutoDock.InternalEnergy import InternalEnergy
from PyAutoDock.MolecularSystem import MolecularSystem
from PyAutoDock.AutoDockScorer import AutoDock41Scorer

def check_types(molecule,std_types):
    d = {}
    for a in molecule.allAtoms:
        d[a.autodock_element] = 1
    mol_types = d.keys()
    non_std = []
    for t in mol_types:
        if t not in std_types:
            non_std.append(t)
    return non_std            


if __name__ == '__main__':
    import sys
    import getopt


# compute_AutoDock41_score.py -r receptorfilename -l ligandfilename [-o outputfilename -x exclude_torsFreeEnergy -w write_file] -v
    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: compute_AutoDock41_score.py -r receptorfilename -l ligandfilename"
        print
        print "    Description of command..."
        print "         -r     receptorfilename"
        print "         -l     ligandfilename"
        print "    Optional parameters:"
        print "        [-o]    outputfilename (default is 'AutoDock41_scores.txt')"
        print "        [-p]    parameter_library_file (default built-in values from 'AD4.1_bound.dat')"
        print "        [-x]    exclude_torsFreeEnergy (default is to report sum of intermolecular energies plus loss of torsional entropy)"
        print "        [-w]    write_file_mode (default mode is 'a')"
        print "        [-v]    verbose output"

    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'r:l:wvxo:p:h')
    except getopt.GetoptError, msg:
        print 'compute_AutoDock41_score.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    #-r: receptorfilename
    receptorfilename =  None
    #-l: ligandfilename
    ligandfilename =  None
    #-w write new outputfile
    write_file_mode = False
    # optional parameters
    #-o outputfilename
    outputfilename = 'AutoDock41_scores.txt'
    #-p parameter_library_filename
    parameter_library_filename = None
    #-x exclude_torsFreeEnergy, default is to include torsFreeEnergy
    exclude_torsFreeEnergy = False
    #-v verbose
    verbose = None

    #'r:l:wvxo:p:h'
    for o, a in opt_list:
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print 'set verbose to ', True
        if o in ('-r', '--r'):
            receptorfilename = a
            if verbose: print 'set receptorfilename to ', a
        if o in ('-l', '--l'):
            ligandfilename = a
            if verbose: print 'set ligandfilename to ', a
        if o in ('-w', '--w'):
            write_file_mode = True
            if verbose: print 'set write_file_mode to ', True
        if o in ('-x', '--x'):
            exclude_torsFreeEnergy = True
            if verbose: print 'set exclude_torsFreeEnergy to ', exclude_torsFreeEnergy
        if o in ('-o', '--o'):
            outputfilename = a
            if verbose: print 'set outputfilename to ', a
        if o in ('-p', '--p'):
            parameter_library_filename = a
            if verbose: print 'set parameter_library_filename to ', a
        if o in ('-h', '--'):
            usage()
            sys.exit()


    if not receptorfilename:
        print 'compute_AutoDock41_score: receptorfilename must be specified.'
        usage()
        sys.exit()

    if not ligandfilename:
        print 'compute_AutoDock41_score: ligandfilename must be specified.'
        usage()
        sys.exit()

    ad_scorer = AutoDock41Scorer(exclude_torsFreeEnergy=exclude_torsFreeEnergy)
    if parameter_library_filename is not None:
        ad_scorer.read_parameter_file(parameter_library_filename)
    supported_types = ad_scorer.supported_types
    receptor = Read(receptorfilename)[0]
    assert os.path.splitext(receptor.parser.filename)[-1]=='.pdbqt',"receptor file not in required '.pdbqt' format"
    if verbose: print 'read ', receptorfilename
    receptor.buildBondsByDistance()
    rec_non_std = ""
    non_std_types = check_types(receptor, supported_types)
    if len(non_std_types):
        rec_non_std = non_std_types[0]
        if len(non_std_types)>1:
            for t in non_std_types[1:]:
                rec_non_std = rec_non_std + '_' + t
    #check format of receptor
    ligand = Read(ligandfilename)[0]
    assert os.path.splitext(ligand.parser.filename)[-1]=='.pdbqt',"ligand file not in required '.pdbqt' format"
    if verbose: print 'read ', ligandfilename
    ligand.buildBondsByDistance()
    lig_non_std = ""
    non_std_types = check_types(ligand, supported_types)
    if len(non_std_types):
        lig_non_std = non_std_types[0]
        if len(non_std_types)>1:
            for t in non_std_types[1:]:
                lig_non_std = lig_non_std + '_' + t
    mode = 'a'
    first = not os.path.exists(outputfilename)
    if write_file_mode:
        mode = 'w'
        first = True
    if verbose: print 'first is ', first
    optr = open(outputfilename, mode)
    if first:
        tstr = "    Receptor      Ligand   AutoDock4.1Score     estat      hb      vdw   dsolv    tors\n"
        optr.write(tstr)
    #setup the molecular system
    ostr = ""
    if len(lig_non_std):
        ostr = 'ERROR: unable to score ligand "%s" due to presence of non-standard atom type(s): %s\n' %(ligand.name, lig_non_std)
        optr.write(ostr)
    elif len(rec_non_std):
        ostr = 'ERROR: unable to score receptor "%s" due to non-standard atom type(s): %s\n' %(receptor.name, rec_non_std)
        optr.write(ostr)
    else: 
        ms = MolecularSystem()
        ms.add_entities(receptor.allAtoms)
        ms.add_entities(ligand.allAtoms)
        ad_scorer.set_molecular_system(ms)
        #get the scores
        #score per term
        estat, hb,vdw,dsolv = ad_scorer.get_score_per_term()
        torsEnrg = ligand.TORSDOF * ad_scorer.tors_weight
        score = estat + hb + vdw + dsolv + torsEnrg
        ostr = "%12s%12s      % 8.4f        % 6.4f % 6.4f % 6.4f % 6.4f % 6.4f\n" %(receptor.name,ligand.name, score, estat, hb, vdw, dsolv, torsEnrg)
        optr.write(ostr)
        if verbose:
            print 'wrote %s' %outputfilename
    optr.close()

# To execute this command type:
# compute_AutoDock41_score.py -r receptorfilename -l ligandfilename [-o outputfilename -x exclude_torsFreeEnergy -w write_file_mode] -v




