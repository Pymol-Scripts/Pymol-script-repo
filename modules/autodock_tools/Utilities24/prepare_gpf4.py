#!/usr/bin/env python
#
# 
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/prepare_gpf4.py,v 1.10.4.3 2009/03/23 21:54:28 rhuey Exp $
#

import string
import os.path
import glob
from MolKit import Read
from AutoDockTools.GridParameters import GridParameters, grid_parameter_list4
from AutoDockTools.GridParameters import GridParameter4FileMaker
from AutoDockTools.atomTypeTools import AutoDock4_AtomTyper


def usage():
    print "Usage: prepare_gpf4.py -l pdbqt_file -r pdbqt_file "
    print "     -l ligand_filename"
    print "     -r receptor_filename"
    print
    print "Optional parameters:"
    print "    [-i reference_gpf_filename]"
    print "    [-o output_gpf_filename]"
    print "    [-x flexres_filename]"
    print "    [-p parameter=newvalue. For example: -p ligand_types='HD,Br,A,C,OA' ]"
    print "    [-d directory of ligands to use to set types]"
    print "    [-y boolean to center grids on center of ligand]"
    print "    [-n boolean to NOT size_box_to_include_ligand]"
    print "    [-v]"
    print
    print "Prepare a grid parameter file (GPF) for AutoDock4."
    print
    print "   The GPF will by default be <receptor>.gpf. This"
    print "may be overridden using the -o flag."

    
if __name__ == '__main__':
    import getopt
    import sys

    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'vl:r:i:x:o:p:d:yn')
    except getopt.GetoptError, msg:
        print 'prepare_gpf4.py: %s' % msg
        usage()
        sys.exit(2)

    receptor_filename = ligand_filename = None
    list_filename = gpf_filename = gpf_filename = None
    output_gpf_filename = None
    flexres_filename = None
    directory = None
    parameters = []
    verbose = None
    center_on_ligand = False
    size_box_to_include_ligand = True
    for o, a in opt_list:
        if o in ('-v', '--v'):
            verbose = 1
        if o in ('-l', '--l'):
            ligand_filename = a
            if verbose: print 'ligand_filename=', ligand_filename
        if o in ('-r', '--r'):
            receptor_filename = a
            if verbose: print 'receptor_filename=', receptor_filename
        if o in ('-i', '--i'):
            gpf_filename = a
            if verbose: print 'reference_gpf_filename=', gpf_filename
        if o in ('-x', '--x'):
            flexres_filename = a
            if verbose: print 'flexres_filename=', flexres_filename
        if o in ('-o', '--o'):
            output_gpf_filename = a
            if verbose: print 'output_gpf_filename=', output_gpf_filename
        if o in ('-p', '--p'):
            parameters.append(a)
            if verbose: print 'parameters=', parameters
        if o in ('-d', '--d'):
            directory = a
            if verbose: print 'directory=', directory
        if o in ('-y', '--y'):
            center_on_ligand = True
            if verbose: print 'set center_on_ligand to ', center_on_ligand
        if o in ('-n', '--n'):
            size_box_to_include_ligand = False
            if verbose: print 'set size_box_to_include_ligand to ', size_box_to_include_ligand
        if o in ('-h', '--'):
            usage()
            sys.exit()


    if (not receptor_filename) or (ligand_filename is None and directory is None):
        print "prepare_gpf4.py: ligand and receptor filenames"
        print "                    must be specified."
        usage()
        sys.exit()

    gpfm = GridParameter4FileMaker(size_box_to_include_ligand=size_box_to_include_ligand,verbose=verbose)
    if gpf_filename is not None:
        gpfm.read_reference(gpf_filename)
    if ligand_filename is not None:
        gpfm.set_ligand(ligand_filename)
    gpfm.set_receptor(receptor_filename)
    if directory is not None:
        gpfm.set_types_from_directory(directory)
    if flexres_filename is not None:
        flexmol = Read(flexres_filename)[0]
        flexres_types = flexmol.allAtoms.autodock_element
        lig_types = gpfm.gpo['ligand_types']['value'].split()
        all_types = lig_types
        for t in flexres_types:
            if t not in all_types: 
                all_types.append(t)
        all_types_string = all_types[0]
        if len(all_types)>1:
            for t in all_types[1:]:
                all_types_string = all_types_string + " " + t
        gpfm.gpo['ligand_types']['value'] = all_types_string 
    for p in parameters:
        key,newvalue = string.split(p, '=')
        kw = {key:newvalue}
        apply(gpfm.set_grid_parameters, (), kw)
    #gpfm.set_grid_parameters(spacing=1.0)
    if center_on_ligand is True:
        gpfm.gpo['gridcenterAuto']['value'] = 0
        cenx,ceny,cenz = gpfm.ligand.getCenter()
        gpfm.gpo['gridcenter']['value'] = "%.3f %.3f %.3f" %(cenx,ceny,cenz)
    gpfm.write_gpf(output_gpf_filename)


#prepare_gpf4.py -l 1ebg_lig.pdbqt -r 1ebg_rec.pdbqt -p spacing=0.4 -p npts=[60,60,60] -i ref.gpf -o testing.gpf 

