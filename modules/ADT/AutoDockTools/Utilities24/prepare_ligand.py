#!/usr/bin/env python
#
# 
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/prepare_ligand.py,v 1.4 2007/10/08 18:12:27 rhuey Exp $
#
import os 

from MolKit import Read

from AutoDockTools.MoleculePreparation import LigandPreparation



if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: prepare_ligand.py -l filename"
        print
        print "    Description of command..."
        print "         -l     ligand_filename"
        print "    Optional parameters:"
        print "        [-v]    verbose output"
        print "        [-o pdbq_filename] (output filename)"
        print "        [-d]    dictionary to write types list and number of active torsions "

        print "        [-A]    type(s) of repairs to make:\n\t\t bonds_hydrogens, bonds, hydrogens, ''"
        print "        [-C]    preserve initial charges (do not add new charges)"
        print "        [-p]    preserve input charges on atom type, eg -p Zn"
        print "        [-K]    add Kollman charges"
        print "        [-U]    cleanup type:\n\t\t nphs_lps, nphs, lps, '' "
        print "        [-B]    type(s) of bonds to allow to rotate"
        print "        [-R]    index for root"
        print "        [-F]    check for and use largest non-bonded fragment (False)"
        print "        [-M]    interactive (default is automatic)"
        print "        [-I]    string of bonds to inactivate composed of "
        print "                   of zero-based atom indices eg 5_13_2_10  "
        print "                   will inactivate atoms[5]-atoms[13] bond "
        print "                               and atoms[2]-atoms[10] bond "
        print "                      (default is '')"
        print "        [-Z]    inactivate all active torsions     "
        print "                      (default is leave active)"


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'l:vo:d:A:Cp:KU:B:R:MFI:Zh')
    except getopt.GetoptError, msg:
        print 'prepare_ligand.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    #-l: ligand
    ligand_filename =  None
    # optional parameters
    verbose = None
    add_bonds = False
    #-b: repairs to make: add bonds and/or hydrogens
    repairs = ""
    #-C  default: add gasteiger charges 
    charges_to_add = 'gasteiger'
    #-p preserve charges on specific atom types
    preserve_charge_types=''
    #-K add_Kollman #no opt
    #-U: cleanup by merging nphs_lps, nphs, lps
    cleanup  = "nphs_lps"
    #-R rotatable bond type(s) to fix
    allowed_bonds = ""
    #-r  root
    root = 'auto'
    #-o outputfilename
    outputfilename = None
    #-F check_for_fragments
    check_for_fragments = False
    #-I bonds_to_inactivate
    bonds_to_inactivate = ""
    #-Z inactivate_all_torsions
    inactivate_all_torsions = False
    #-m mode 
    mode = 'automatic'

    dict = None

    #'l:vo:d:A:Cp:KU:B:R:MFI:Z'
    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('-l', '--l'):
            ligand_filename = a
            if verbose: print 'set ligand_filename to ', a
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print 'set verbose to ', True
        if o in ('-o', '--o'):
            outputfilename = a
            if verbose: print 'set outputfilename to ', a
        if o in ('-d', '--d'):
            dict = a
            if verbose: print 'set dict to ', a
        if o in ('-A', '--A'):
            repairs = a
            if verbose: print 'set repairs to ', a
        if o in ('-C', '--C'):
            charges_to_add = None
            if verbose: print 'do not add charges'
        if o in ('-p', '--p'):
            preserve_charge_types+=a
            preserve_charge_types+=','
            if verbose: print 'preserve initial charges on ', preserve_charge_types
        if o in ('-K', '--K'):
            charges_to_add = 'Kollman'
            if verbose: print 'add Kollman charges'
        if o in ('-p', '--p'):
            preserve_charge_types+=a
            preserve_charge_types+=','
            if verbose: print 'preserve initial charges on ', preserve_charge_types
        if o in ('-U', '--U'):
            cleanup  = a
            if verbose: print 'set cleanup to merge ', a
        if o in ('-B', '--B'):
            allowed_bonds = a
            if verbose: print 'allow ', a, 'bonds set to rotate'
        if o in ('-R', '--R'):
            root = a
            if verbose: print 'set root to ', root
        if o in ('-F', '--F'):
            check_for_fragments = True
            if verbose: print 'set check_for_fragments to True'
        if o in ('-M', '--M'):
            mode = a
            if verbose: print 'set mode to ', a
        if o in ('-I', '--I'):
            bonds_to_inactivate = a
            if verbose: print 'set bonds_to_inactivate to ', a
        if o in ('-Z', '--Z'):
            inactivate_all_torsions = True
            if verbose: print 'set inactivate_all_torsions to ', inactivate_all_torsions
        if o in ('-h', '--'):
            usage()
            sys.exit()


    if not  ligand_filename:
        print 'prepare_ligand: ligand filename must be specified.'
        usage()
        sys.exit()

    #what about nucleic acids???

    mols = Read(ligand_filename)
    if verbose: print 'read ', ligand_filename
    mol = mols[0]
    if len(mols)>1:
        if verbose: 
            print "more than one molecule in file"
        #use the one molecule with the most atoms
        ctr = 1
        for m in mols[1:]:
            ctr += 1
            if len(m.allAtoms)>len(mol.allAtoms):
                mol = m
                if verbose:
                    print "mol set to ", ctr, "th molecule with", len(mol.allAtoms), "atoms"

    #possible clean-up???
    mol.buildBondsByDistance()

    if charges_to_add is not None:
        preserved = {}
        preserved_types = preserve_charge_types.split(',') 
        for t in preserved_types:
            if not len(t): continue
            ats = mol.allAtoms.get(lambda x: x.autodock_element==t)
            for a in ats:
                if a.chargeSet is not None:
                    preserved[a] = [a.chargeSet, a.charge]


    if verbose:
        print "setting up LPO with mode=", mode,
        print "and outputfilename= ", outputfilename
        print "and check_for_fragments=", check_for_fragments
        print "and bonds_to_inactivate=", bonds_to_inactivate
    LPO = LigandPreparation(mol, mode, repairs, charges_to_add, 
                            cleanup, allowed_bonds, root, 
                            outputfilename=outputfilename,
                            dict=dict, check_for_fragments=check_for_fragments,
                            bonds_to_inactivate=bonds_to_inactivate, 
                            inactivate_all_torsions=inactivate_all_torsions)
    #do something about atoms with too many bonds (?)
    #FIX THIS: could be peptide ligand (???)
    #          ??use isPeptide to decide chargeSet??
    if charges_to_add is not None:
        #restore any previous charges
        for atom, chargeList in preserved.items():
            atom._charges[chargeList[0]] = chargeList[1]
            atom.chargeSet = chargeList[0]


# To execute this command type:
# prepare_ligand.py -l pdb_file -v




