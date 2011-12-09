#!/usr/bin/env python
#
# 
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/prepare_ligand_dict.py,v 1.2 2007/10/08 18:12:45 rhuey Exp $
#

import os
from MolKit import Read


class LigandDictionaryMaker:
    """Accept a <ligand>.pdbq and <receptor>.pdbqs and create
    <ligand>_<receptor>.dpf
    """

    def __init__(self, filename="ligand_dict.py", verbose = None):
        self.verbose = verbose
        self.dict = {}
        self.dict['atom_types'] = {}
        self.dict['rbonds'] = {}
        self.dict['zero_charge'] = {}
        self.filename = filename


    def write_ligand_filename(self, ligandfilename):
        ligand = Read(ligandfilename)[0]
        self.write_ligand(ligand)


    def write_ligand(self, ligand):
        if not os.path.exists(self.filename):
            fptr = open(self.filename, 'a')
            ostr = "summary = d = {}\n"
            fptr.write(ostr)
        else:
            fptr = open(self.filename, 'a')
        type_dict = {}
        for a in ligand.allAtoms:
            type_dict[a.autodock_element] = 1
        atom_types = type_dict.keys()
        atom_types.sort()
        ostr = "d['" + ligand.name +"'] = {"
        ostr = ostr + "'atom_types': [" 
        for t in atom_types[:-1]:
            ostr = ostr + "'%s', "%t
        ostr = ostr + "'%s' "%atom_types[-1]
        ostr = ostr + "],\n\t\t\t\t'rbonds':" + str(len(ligand.torTree.torsionMap))
        ostr = ostr + ",\n\t\t\t\t'zero_charge' : ["
        zc_ats = ligand.allAtoms.get(lambda x: x.charge==0)
        for at in zc_ats:
            ostr = ostr + "'%s', "%at.name
        ostr = ostr + "],\n\t\t\t\t}\n"
        fptr.write(ostr)
        fptr.close()

        #to write the dictionary inside out
        #if not os.path.exists(self.filename):
        #    fptr = open(self.filename, 'a')
        #    ostr = "summary = d = {'atom_types':{}, 'rbonds':{}, 'zero_charge':{}}\n"
        #ostr = "d['atom_types']['" + ligand.name +"'] = ["
        #for t in atom_types[:-1]:
        #    ostr = ostr + "'%s', "%t
        #ostr = ostr + "'%s' "%atom_types[-1]
        #ostr = ostr + "]\nd['rbonds']['" + ligand.name + "'] = " + str(len(ligand.torTree.torsionMap))
        #ostr = ostr + "\nd['zero_charge']['" + ligand.name + "'] =  ["
        #zc_ats = ligand.allAtoms.get(lambda x: x.charge==0)
        #for at in zc_ats:
        #    ostr = ostr + "'%s',"%at.name
        #ostr = ostr + "]\n"
        #fptr.write(ostr)
        #fptr.close()

 

def usage():
    print "Usage: prepare_ligand_dict.py -l ligand_filename -d ligand_dict_filename"
    print
    print "Optional parameters:"
    print "    [-d ligand_dict_filename]"
    print "    [-v] verbose output"
    print
    print "Write a python dictionary 'summary' comprised of dictionaries of atom_types, rbonds and zero_charge atoms with an entry for each ligand filename."
    print
    print "   The filename will be ligand_dict.py.  "
    print "This may be overridden using the -d flag."

    
if __name__ == '__main__':
    import getopt
    import sys

    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'hvl:d:')
    except getopt.GetoptError, msg:
        print 'prepare_ligand_dict.py: %s' % msg
        usage()
        sys.exit(2)

    ligand_filename = None
    ligand_dict_filename = None
    verbose = None
    for o, a in opt_list:
        if verbose: print "o=", o, ' a=', a
        if o in ('-v', '--v'):
            verbose = 1
            if verbose: print 'verbose output'
        if o in ('-l', '--l'):
            ligand_filename = a
            if verbose: print 'ligand_filename =', ligand_filename
        if o in ('-d', '--d'):
            ligand_dict_filename = a
            if verbose: print 'ligand_dict_filename =', ligand_dict_filename
        if o in ('-h', '--'):
            usage()
            sys.exit()


    if not ligand_filename:
        print "prepare_ligand_dict.py: ligand filename must be specified."
        usage()
        sys.exit()

    ldm = LigandDictionaryMaker(filename=ligand_dict_filename, verbose=verbose)
    ldm.write_ligand_filename( ligand_filename)
    
#prepare_ligand_dict.py -l indinavir.pdbq -d hiv_ligands_dict.py"

