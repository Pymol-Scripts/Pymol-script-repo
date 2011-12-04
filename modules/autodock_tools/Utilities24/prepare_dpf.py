#!/usr/bin/env python
#
# 
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/prepare_dpf.py,v 1.3.4.1 2009/04/06 18:49:44 rhuey Exp $
#

import string
import os.path
from MolKit import Read
from AutoDockTools.DockingParameters import DockingParameters, genetic_algorithm_list, \
                genetic_algorithm_local_search_list, local_search_list,\
                simulated_annealing_list


class DockingParameterFileMaker:
    """Accept a <ligand>.pdbq and <receptor>.pdbqs and create
    <ligand>_<receptor>.dpf
    """

    def __init__(self, verbose = None):
        self.verbose = verbose
        self.dpo = DockingParameters()


    def set_ligand(self, ligand_filename): 
        self.ligand_filename = os.path.basename(ligand_filename)
        if verbose:
            print "set ligand_filename to", self.ligand_filename
        self.dpo.set_ligand(ligand_filename)
        #expect a filename like ind.out.pdbq: get 'ind' from it
        self.ligand_stem = string.split(self.ligand_filename,'.')[0]
        if verbose: print "set ligand_stem to", self.ligand_stem
        self.ligand = Read(ligand_filename)[0]
        if self.ligand==None:
            print 'ERROR reading: ', ligand_filename
            return 
        if verbose: 
            print "read ", self.ligand.name
        #set dpo:
        #move
        self.dpo['move']['value'] = self.ligand_filename
        if verbose: print "set move to ", self.dpo['move']['value']
        #ndihe
        #assumes ligand has torTree
        self.dpo['ndihe']['value'] = self.ligand.parser.keys.count("BRANCH")
        #self.dpo['ndihe']['value'] = len(self.ligand.torTree.torsionMap)
        if verbose: print "set ndihe to ", self.dpo['ndihe']['value']
        #torsdof
        #caution dpo['torsdof']['value'] is a list [ndihe, 0.3113]
        self.dpo['torsdof']['value'][0] = self.ligand.TORSDOF
        if verbose: print "set torsdof to ", self.dpo['torsdof']['value']
        #types
        d = {}
        for a in self.ligand.allAtoms:
            d[a.autodock_element] = 1
        sortKeyList =  ['C','A','N','O','S','H','P','n','f','F','c','b','I','M']
        lig_types = ""
        for t in sortKeyList:
            if t in d.keys():
                lig_types = lig_types + t
        self.ligand.types = lig_types
        self.dpo['types']['value'] = self.ligand.types
        if verbose: print "set types to ", self.dpo['types']['value']
        #about
        self.ligand.getCenter()
        cen = self.ligand.center
        self.dpo['about']['value'] =  [round(cen[0],4), round(cen[1],4),\
                                        round(cen[2],4)]
        if verbose: print "set about to ", self.dpo['about']['value']
        

    def set_receptor(self, receptor_filename):
        self.receptor_filename = os.path.basename(receptor_filename)
        self.receptor_stem = string.split(self.receptor_filename, '.')[0]
        self.dpo.set_receptor(receptor_filename)


    #def set_docking_parameters(self, newdict={}):
    def set_docking_parameters(self, **kw):
        """Any docking paramters should be set here
        """
        # like this: 
        # newdict = {'ga_num_evals':1750000, 'ga_pop_size':150,
        #            'ga_run':20, 'rmstol':2.0}
        # self.dpo['<parameter>']['value'] = <new value>
        # eg self.dpo['rmstol']['value'] = 2.0
        for parm, newvalue in kw.items():
            #print "parm=", parm, ' newvalue=', newvalue
            self.dpo[parm]['value'] = newvalue
            if parm=='set_sw1':
                self.dpo['set_psw1']['value'] = not newvalue
            if parm=='set_psw1':
                self.dpo['set_sw1']['value'] = not newvalue


    def write_dpf(self, dpf_filename,
              parm_list = genetic_algorithm_local_search_list):
        if not dpf_filename:
            dpf_filename = "%s%s%s%s" % \
                           (self.ligand_stem, "_",
                            self.receptor_stem, ".dpf")
        # now that we have a filename...
        if self.verbose:
            print "writing ", dpf_filename
        self.dpo.write(dpf_filename, parm_list)

 

def usage():
    print "Usage: prepare_dpf.py -l pdbq_file -r pdbqs_file"
    print "    -l ligand_filename"
    print "    -r receptor_filename"
    print
    print "Optional parameters:"
    print "    [-o output dpf_filename]"
    print "    [-i template dpf_filename]"
    print "    [-p parameter_name=new_value]"
    print "    [-k list of parameters to write]"
    print "    [-L] use local search parameters"
    print "    [-S] use simulated annealing search parameters"
    print "    [-v] verbose output"
    print
    print "Prepare a docking parameter file (DPF) for AutoDock."
    print
    print "   The DPF will by default be <ligand>_<receptor>.dpf. This"
    print "may be overridden using the -o flag."

    
if __name__ == '__main__':
    import getopt
    import sys

    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'LShvl:r:i:o:p:k:')
    except getopt.GetoptError, msg:
        print 'prepare_dpf.py: %s' % msg
        usage()
        sys.exit(2)

    receptor_filename = ligand_filename = None
    dpf_filename = None
    template_filename = None
    parameters = []
    parameter_list = genetic_algorithm_local_search_list
    verbose = None
    for o, a in opt_list:
        if verbose: print "o=", o, ' a=', a
        if o in ('-v', '--v'):
            verbose = 1
            if verbose: print 'verbose output'
        if o in ('-l', '--l'):   #ligand filename
            ligand_filename = a
            if verbose: print 'ligand_filename =', ligand_filename
        if o in ('-r', '--r'):   #receptor filename
            receptor_filename = a
            if verbose: print 'receptor_filename =', receptor_filename
        if o in ('-i', '--i'):   #input reference
            template_filename = a
            if verbose: print 'template_filename =', template_filename
        if o in ('-o', '--o'):   #output filename
            dpf_filename = a
            if verbose: print 'output dpf_filename =', dpf_filename
        if o in ('-p', '--p'):   #parameter
            parameters.append(a)
            if verbose: print 'parameters =', parameters
        if o in ('-k', '--k'):   #parameter_list_to_write
            parameter_list = a
            if verbose: print 'parameter_list =', parameter_list
        if o in ('-L', '--L'):   #parameter_list_to_write
            parameter_list = local_search_list
            if verbose: print 'parameter_list =', parameter_list
        if o in ('-S', '--S'):   #parameter_list_to_write
            parameter_list = simulated_annealing_list
            if verbose: print 'parameter_list =', parameter_list
        if o in ('-h', '--'):
            usage()
            sys.exit()


    if (not receptor_filename) or (not ligand_filename):
        print "prepare_dpf.py: ligand and receptor filenames"
        print "                    must be specified."
        usage()
        sys.exit()

    dm = DockingParameterFileMaker(verbose=verbose)
    if template_filename is not None:  #setup values by reading dpf
        dm.dpo.read(template_filename)
    dm.set_ligand(ligand_filename)
    dm.set_receptor(receptor_filename)
    #dm.set_docking_parameters( ga_num_evals=1750000,ga_pop_size=150, ga_run=20, rmstol=2.0)
    for p in parameters:
        key,newvalue = string.split(p, '=')
        #detect string reps of lists: eg "[1.,1.,1.]"
        if newvalue[0]=='[':
            nv = []
            for item in newvalue[1:-1].split(','):
                nv.append(float(item))
            #print "nv=", nv
            newvalue = nv
        kw = {key:newvalue}
        apply(dm.set_docking_parameters, (), kw)
    dm.write_dpf(dpf_filename, parameter_list)
    
#prepare_dpf.py -l indinavir.pdbq -r 1hsg.pdbqs -p ga_num_evals=20000000 -p ga_pop_size=150 -p ga_run=17 -i ref.dpf -o testing.dpf 

