#!/usr/bin/env python

import os
from string import find, strip

from MolKit import Read
from MolKit.molecule import AtomSet
from MolKit.pdbWriter import PdbqtWriter
from AutoDockTools.atomTypeTools import AutoDock4_AtomTyper




if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: pdbq_to_pdbqt.py -s ligand_stem"
        print
        print "    Description of command..."
        print "        [-s]    stem of ligand.pdbq file"
        print "    Optional parameters:"
        print "        [-o]    alternative pdbqt_filename"
        print "        [-v]    verbose output"


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 's:o:v')

    except getopt.GetoptError, msg:
        print 'pdbq_to_pdbqt.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    #-s: pdbq_stem
    pdbq_filename =  None

    # optional parameters
    verbose = None
    pdbqt_filename = None

    #'s:v'
    for o, a in opt_list:
        if o in ('-s', '--s'):
            pdbq_filename = a + '.pdbq'
            if pdbqt_filename is None:
                pdbqt_filename = a + '.pdbqt'
            if verbose: 
                print 'set pdbq_filename to ', a, 
                print " and pdbqt_filename to ", pdbqt_filename
        if o in ('-o', '--o'):
            pdbqt_filename = a 
            if verbose: 
                print 'set pdbqt_filename to ', a
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print 'set verbose to ', True
        if o in ('-h', '--'):
            usage()
            sys.exit()



    if not pdbq_filename:
        print 'pdbq_to_pdbqt: stem of pdbq_filename must be specified.'
        usage()
        sys.exit()

    #what about nucleic acids???

    mols = Read(pdbq_filename)
    if verbose: print 'read ', pdbq_filename
    mol = mols[0]
    mol.buildBondsByDistance()

    #possible clean-up???
    #need to type atoms + assign babel_types
    AD4_typer = AutoDock4_AtomTyper()
    AD4_typer.setAutoDockElements(mol)
    if verbose:
        print "set autodock4 autodock_element for ", mol.name

    writer = PdbqtWriter()
    fptr = open(pdbqt_filename, 'w')
    ctr = 0
    for line in mol.parser.allLines:
        #for at in mol.allAtoms:
        if find(line, 'ATOM')<0 and find(line, "HETA")<0:
            fptr.write(line)
        else:
            name = strip(line[12:16])
            #ats = mol.allAtoms.get(lambda x: x.name==name)
            ##there could be (?) two atoms with the same name
            #if len(ats)==1:
            #    this_atom = ats[0]
            #else:
            x_coords = float(line[30:38])
            y_coords = float(line[38:46])
            z_coords = float(line[46:54])
            #this_atom = mol.allAtoms.get(lambda x:x.coords[0]==x_coords)[0]
            this_atom = mol.allAtoms.get(lambda x:x.coords[0]==x_coords)
            if this_atom is None or len(this_atom)==0:
                print name, x_coords,' do not correspond to an atom in ', pdbq_filename
                raise RuntimeError
            elif len(this_atom)>1:
                if verbose: print "using y coord comparison also for ", this_atom.name[0], this_atom.number
                this_atom = this_atom.get(lambda x:x.coords[1]==y_coords)
                if this_atom is None or len(this_atom)==0:
                    print name, x_coords, y_coords, ' do not correspond to an atom in ', pdbq_filename
                    raise RuntimeError
                elif len(this_atom)>1:
                    if verbose: print "using z coord comparison also for ", this_atom.name[0], this_atom.number
                    this_atom = this_atom.get(lambda x:x.coords[2]==z_coords)
                    if this_atom is None or len(this_atom)==0:
                        print name, x_coords, y_coords, z_coords, ' do not correspond to an atom in ', pdbq_filename
                        raise RuntimeError
                    else:
                        #x and y and z coords were required
                        this_atom = this_atom[0]
                else:
                    #x and y coords were enough
                    this_atom = this_atom[0]

            else:
                #x coords were enough
                this_atom = this_atom[0]
                
            if this_atom.name=='c': 
                this_atom.name = 'Cl'
                this_atom.element = 'Cl'
                this_atom.autodock_element ='CL'
            elif this_atom.name=='b': 
                print "processing ", this_atom.name
                this_atom.name = 'Br'
                this_atom.element = 'Br'
                this_atom.autodock_element ='BR'
            elif this_atom.name=='f': 
                this_atom.name = 'Fe'
                this_atom.element = 'Fe'
                this_atom.autodock_element ='FE'
            writer.write_atom(fptr, this_atom)
        ctr = ctr + 1
    fptr.close()
    if verbose:
        print "wrote ", ctr, " atoms to", pdbqt_filename
    

# To execute this command type:
# pdbq_to_pdbqt.py -s pdbq_file_stem -v




