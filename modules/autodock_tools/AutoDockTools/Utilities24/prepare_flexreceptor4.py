#!/usr/bin/env python
#
# 
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/prepare_flexreceptor4.py,v 1.5.2.1 2009/07/17 16:55:53 rhuey Exp $
#
import os 

from MolKit import Read
from MolKit.protein import ProteinSet, ResidueSet
from MolKit.molecule import BondSet
from MolKit.stringSelector import CompoundStringSelector

from AutoDockTools.MoleculePreparation import AD4FlexibleReceptorPreparation



if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: prepare_flexreceptor4.py -r receptor_filename -s list_of_names_of_residues_to_move"
        print "    Description of command..."
        print "         -r     receptor_filename (.pdbqt)"
        print "         -s     specification for flex residues" 
        print "                Use underscores to separate residue names:"
        print "                  ARG8_ILE84  "
        print "                Use commas to separate 'full names' which uniquely identify residues:"
        print "                  hsg1:A:ARG8_ILE84,hsg1:B:THR4 "
        print "                [syntax is molname:chainid:resname]"
        print "    Optional parameters:"
        print "        [-v]    verbose output"
        print "        [-N]    type(s) of bonds to disallow: "
        print "        [-P]    pairs of atom names bonds between which to disallow: "
        print "        [-g pdbqt_filename] (rigid output filename)"
        print "        [-x pdbqt_filename] (flexible output filename)"

    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'r:vs:N:P:g:x:h')
    except getopt.GetoptError, msg:
        print 'prepare_flexreceptor4.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    #-r: ligand
    receptor_filename =  None
    #-s: residues_to_move
    residues_to_move =  None
    # optional parameters
    verbose = None
    #-N: type of bonds to  disallow
    disallow = ""
    #-P: pairs of atom names bonds between which to  disallow
    disallowed_pairs = ""
    #-g  : rigid output filename
    rigid_filename = None
    #-x  : flexible output filename
    flexres_filename = None

    #'r:vs:N:g:x:h'
    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('-r', '--r'):
            receptor_filename = a
            if verbose: print 'set receptor_filename to ', a
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print 'set verbose to ', True
        if o in ('-s', '--s'):
            residues_to_move = a
            if verbose: print 'set residues_to_move to ', a
        if o in ('-N', '--N'):
            disallow = a
            if verbose: print 'set disallow to ', a
        if o in ('-P', '--P'):
            disallowed_pairs = a
            if verbose: print 'set disallowed_pairs to ', a
        if o in ('-g', '--g'):
            rigid_filename = a
            if verbose: print 'set rigid_filename to ', a
        if o in ('-x', '--'):
            flexres_filename = a
            if verbose: print 'set flexres_filename to ', a
        if o in ('-h', '--'):
            usage()
            sys.exit()


    if not  receptor_filename:
        print 'prepare_flexreceptor4: receptor filename must be specified!\n'
        usage()
        sys.exit()

    if not  residues_to_move:
        print 'prepare_flexreceptor4: residues to move must be specified!\n'
        usage()
        sys.exit()

    extension = os.path.splitext(receptor_filename)[1]
    if extension!=".pdbqt":
        print 'prepare_flexreceptor4: receptor file must be in .pdbqt format\n'
        usage()
        sys.exit()

    r = Read(receptor_filename)[0]
    r.buildBondsByDistance()
    if verbose: print 'read ', receptor_filename

    all_res = ResidueSet()
    # hsg1:A:ARG8_ILE82;hsg1:B:THR4 
    # ARG8_ILE84 
    names_by_chain = residues_to_move.split(',')
    if verbose: print "names_by_chain=", names_by_chain
    #1. ['hsg1:A:ARG8_ILE82','hsg1:B:THR4'] 
    # OR 
    #2. ARG8_ILE84
    for res_names_by_chain in names_by_chain:
        #check for ':'
        if verbose: print "processing res_names_by_chain=", res_names_by_chain
        if res_names_by_chain.find(':')==-1:
            # no ':' so treat simple case
            # ARG8_ILE84 
            res_names.split('_')
            for n in res_names:
                res = r.chains.residues.get(lambda x: x.name==n)
                all_res += res
                if verbose: print "case no ':' get: adding ", res.name, " to ", all_res
        else:
            # 'hsg1:A:ARG8_ILE82'
            n = res_names_by_chain.replace("_", ",")
            if verbose: print "after comma replaced _ n is ", n
            # 'hsg1:A:ARG8,ILE82'
            res, msg = CompoundStringSelector().select(ProteinSet([r]), n)
            if verbose: print "n selection =", res, " msg=", msg
            if len(res):
                all_res += res
            else:
                print "no residue found using string ", n
    if verbose:
        print "all_res=", all_res.full_name(), 'all_res.__class__=', all_res.__class__
    #?check for duplicates
    d = {}
    for res in all_res: d[res] = 1
    all_res = d.keys()
    all_res = ResidueSet(all_res).uniq()
    all_res.sort()
    if verbose: print "after all_res=", all_res.full_name()

    #inactivate specified bonds
    #disallowed_Pairs "CA_CB:CB_CG:C_CA"
    all_bnds = BondSet()
    bnd_pairs = disallowed_pairs.split(':')
    for pair in bnd_pairs:
        names = pair.split('_')
        bnds = all_res.atoms.bonds[0].get(lambda x: x.atom1.name in names and x.atom2.name in names)
        all_bnds += bnds
   
    fdp = AD4FlexibleReceptorPreparation(r, residues=all_res, rigid_filename=rigid_filename, 
                                            flexres_filename=flexres_filename,
                                            non_rotatable_bonds=all_bnds)


# To execute this command type:
# prepare_flexreceptor4.py -l filename -r receptor_filename -s list_of_names_of_residues_to_move"




