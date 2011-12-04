#!/usr/bin/env python
#
# 
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/prepare_dpf42.py,v 1.1.2.2 2009/05/19 21:03:55 rhuey Exp $
#

import string
import os.path
from MolKit import Read
from AutoDockTools.DockingParameters import DockingParameters, genetic_algorithm_list4_2, \
                genetic_algorithm_local_search_list4_2, local_search_list4_2,\
                simulated_annealing_list4_2

from AutoDockTools.atomTypeTools import AutoDock4_AtomTyper
                


class DockingParameter42FileMaker:
    """Accept a <ligand>.pdbqt and <receptor>.pdbqt and create
    <ligand>_<receptor>42.dpf
    """

    def __init__(self, verbose = None):
        self.verbose = verbose
        self.dpo = DockingParameters()


    def getTypes(self, molecule):
        if not len(molecule.allAtoms.bonds[0]):
            molecule.buildBondsByDistance()
        ad4_typer = AutoDock4_AtomTyper(verbose=self.verbose)
        ad4_typer.setAutoDockElements(molecule)
        dict = {}
        for a in molecule.allAtoms:
            dict[a.autodock_element] = 1
        d_types = dict.keys()
        d_types.sort()
        mol_types = d_types[0]
        for t in d_types[1:]:
            mol_types = mol_types + " " + t
        if self.verbose: print "end of getTypes: types=", mol_types, ' class=', mol_types.__class__
        return mol_types


    def set_write_all(self, value):
        verbose = self.verbose
        self.dpo['write_all_flag']['value'] = True
        if verbose: print "set write_all_flag to", self.dpo['write_all_flag']['value']


    def set_ligand(self, ligand_filename): 
        verbose = self.verbose
        self.ligand_filename = os.path.basename(ligand_filename)
        if verbose: print "set ligand_filename to", self.ligand_filename
        self.dpo.set_ligand(ligand_filename)
        #expect a filename like ind.out.pdbq: get 'ind' from it
        self.ligand_stem = string.split(self.ligand_filename,'.')[0]
        if verbose: print "set ligand_stem to", self.ligand_stem
        self.ligand = Read(ligand_filename)[0]
        if self.ligand==None:
            print 'ERROR reading: ', ligand_filename
            return 
        if verbose: print "read ", self.ligand.name
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
        #caution dpo['torsdof4']['value'] is a list [ndihe, 0.274]
        try:
            self.dpo['torsdof4']['value'][0] = self.ligand.TORSDOF
        except:
            print 'setting torsdof to ligand.ndihe=', self.ligand.ndihe
            self.dpo['torsdof4']['value'][0] = self.ligand.ndihe
        if verbose: print "set torsdof4 to ", self.dpo['torsdof4']['value']
        #types
        self.ligand.types = self.getTypes(self.ligand)
        self.dpo['ligand_types']['value'] = self.ligand.types
        if verbose: print "set types to ", self.dpo['ligand_types']['value']
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


    def set_flexres(self, flexres_filename):
        flexmol = Read(flexres_filename)[0]
        flexres_filename = os.path.basename(flexres_filename)
        self.dpo['flexres_flag']['value'] = True
        self.dpo['flexres']['value'] = flexres_filename
        #make sure each atom type in flexres molecule is in ligand_types
        d = {}
        current_types = self.dpo['ligand_types']['value'].split()
        for t in current_types:
            d[t] = 1
        for a in flexmol.allAtoms:
            d[a.autodock_element] = 1
        self.dpo['ligand_types']['value'] = string.join(d.keys())


    def set_docking_parameters(self, **kw):
        """Any docking parameters should be set here
        """
        # like this: 
        # newdict = {'ga_num_evals':1750000, 'ga_pop_size':150,
        #            'ga_run':20, 'rmstol':2.0}
        # self.mv.dpo['<parameter>']['value'] = <new value>
        for parm, newvalue in kw.items():
            self.dpo[parm]['value'] = newvalue
            if parm=='set_sw1':
                self.dpo['set_psw1']['value'] = not newvalue
            if parm=='set_psw1':
                self.dpo['set_sw1']['value'] = not newvalue
            if parm=='flexres':
                self.set_flexres(newvalue) 
            if parm=='write_all':
                self.set_write_all(newvalue) 


    def write_dpf(self, dpf_filename,
              parm_list = genetic_algorithm_local_search_list4_2, 
              pop_seed = False):
        if not dpf_filename:
            dpf_filename = "%s%s%s%s" % \
                           (self.ligand_stem, "_",
                            self.receptor_stem, ".dpf")
        # now that we have a filename...
        # set initial conformation
        if pop_seed:
            self.dpo['tran0']['value'] = self.dpo['about']['value']
            self.dpo['quat0']['value'] = '1.0 0. 0. 0.'
            dihe0 = '0. '*self.dpo['ndihe']['value']
            dihe0.rstrip()
            self.dpo['dihe0']['value'] = dihe0 
        if self.verbose:
            print "writing ", dpf_filename
        self.dpo.write42(dpf_filename, parm_list)

 

def usage():
    print "Usage: prepare_dpf42.py -l pdbqt_file -r pdbqt_file"
    print "    -l ligand_filename"
    print "    -r receptor_filename"
    print
    print "Optional parameters:"
    print "    [-o output dpf_filename]"
    print "    [-i template dpf_filename]"
    print "    [-x flexres_filename]"
    print "    [-p parameter_name=new_value]"
    print "    [-k list of parameters to write]"
    print "    [-v] verbose output"
    print "    [-L] use local search parameters"
    print "    [-S] use simulated annealing search parameters"
    print "    [-s] seed population using ligand's present conformation"
    print
    print "Prepare a docking parameter file (DPF) for AutoDock42."
    print
    print "   The DPF will by default be <ligand>_<receptor>.dpf. This"
    print "may be overridden using the -o flag."

    
if __name__ == '__main__':
    import getopt
    import sys

    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'sLShvl:r:i:o:x:p:k:')
    except getopt.GetoptError, msg:
        print 'prepare_dpf42.py: %s' % msg
        usage()
        sys.exit(2)

    receptor_filename = ligand_filename = None
    dpf_filename = None
    template_filename = None
    flexres_filename = None
    parameters = []
    parameter_list = genetic_algorithm_local_search_list4_2
    pop_seed = False
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
        if o in ('-x', '--x'):   #flexres_filename 
            flexres_filename = a
            if verbose: print 'flexres_filename =', flexres_filename
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
            parameter_list = local_search_list4_2
            if verbose: print 'parameter_list =', parameter_list
        if o in ('-S', '--S'):   #parameter_list_to_write
            parameter_list = simulated_annealing_list4_2
            if verbose: print 'parameter_list =', parameter_list
        if o in ('-h', '--'):
            usage()
            sys.exit()
        if o in ('-s'):
            pop_seed = True


    if (not receptor_filename) or (not ligand_filename):
        print "prepare_dpf42.py: ligand and receptor filenames"
        print "                    must be specified."
        usage()
        sys.exit()


    dm = DockingParameter42FileMaker(verbose=verbose)
    if template_filename is not None:  #setup values by reading dpf
        dm.dpo.read(template_filename)
    dm.set_ligand(ligand_filename)
    dm.set_receptor(receptor_filename)
    if flexres_filename is not None:
        flexmol = Read(flexres_filename)[0]
        flexres_types = flexmol.allAtoms.autodock_element
        lig_types = dm.dpo['ligand_types']['value'].split()
        all_types = lig_types
        for t in flexres_types:
            if t not in all_types: 
                all_types.append(t)
        all_types_string = all_types[0]
        if len(all_types)>1:
            for t in all_types[1:]:
                all_types_string = all_types_string + " " + t
                if verbose: print "adding ", t, " to all_types->", all_types_string
        dm.dpo['ligand_types']['value'] = all_types_string 
        dm.dpo['flexres']['value'] = flexres_filename
        dm.dpo['flexres_flag']['value'] = True
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
        elif 'flag' in key:
            if newvalue in ['1','0']:
                newvalue = int(newvalue)
            if newvalue =='False':
                newvalue = False
            if newvalue =='True':
                newvalue = True
        kw = {key:newvalue}
        apply(dm.set_docking_parameters, (), kw)
        if key not in parameter_list:
            #special hack for output_pop_file
            if key=='output_pop_file':
                parameter_list.insert(parameter_list.index('set_ga'), key)
            else:
                parameter_list.append(key) 
    dm.write_dpf(dpf_filename, parameter_list, pop_seed)
    
#prepare_dpf42.py -l indinavir.pdbqt -r 1hsg.pdbqt -p ga_num_evals=25000000 -p ga_pop_size=150 -p ga_run=17 -i ref.dpf -o testing.dpf 

