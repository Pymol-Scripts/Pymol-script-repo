## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 


#############################################################################
#
# Author: William LINDSTROM
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/DockingParameters.py,v 1.79.2.11 2009/06/03 17:48:02 rhuey Exp $
#
#
# $Id: DockingParameters.py,v 1.79.2.11 2009/06/03 17:48:02 rhuey Exp $
#
#
#

from energyConstants import Rij, epsij
import UserDict
import string
import os.path
import sys
import types
from MolKit import Read
from AutoDockTools.atomTypeTools import AutoDock4_AtomTyper

class DockingParameters(UserDict.UserDict):
    def __init__(self, receptor_filename='', ligand_filename='', flexres_filename=''):
        UserDict.UserDict.__init__(self)

        # if the docking parameters have been read from or written
        # to a file otherthan stdout,
        # then the following instance variables will be set:
        self.dpf_filename = ''
        self.dpf_written_filename = ''
        self.file_params = []

        # begin dictionary
        self[ 'about' ] = {
            'keyword' : 'about' ,
            'default' : [0., 0., 0.],
            'comment' : "small molecule center",
            'value'   : [0., 0., 0.]
        }
        self[ 'accs' ] = {
            'keyword' : 'accs' ,
            'default' : 100,
            'comment' : "maximum number of accepted steps per cycle",
            'value'   : 100
        }
        self[ 'analysis' ] = {
            'keyword' : 'analysis' ,
            'default' : 1,
            'comment' : 'perform a ranked cluster analysis',
            'value'   : 1, # true or false
        }
        self[ 'axisangle0' ] = {
            'keyword' : 'axisangle0' ,
            'default' : 'random',
            'comment' : "initial orientation",
            'value'   : 'random'
        }
        self[ 'cluster' ] = {
            'keyword' : 'cluster' ,
            'default' : '',
            'comment' : 'structure binning',
            'value'   : ''
        }
        self[ 'compute_unbound_extended' ] = {
            'keyword' : 'compute_unbound_extended' ,
            'default' : '',
            'comment' : "compute extended ligand energy",
            'value'   : '' 
        }
        self[ 'compute_unbound_extended_flag' ] = {
            'keyword' : 'compute_unbound_extended_flag' ,
            'default' : 1,
            'comment' : "whether to compute unbound ligand energy",
            'value'   : 1  # True/False
        }
        self[ 'custom_parameter_file' ] = {
            'keyword' : 'custom_parameter_file' ,
            'default' : 0,
            'comment' : "use custom parameter library",
            'value'   : 0,
        }
        self[ 'cycles' ] = {
            'keyword' : 'cycles' ,
            'default' : 50,
            'comment' : "number of temperature reduction cycles",
            'value'   : 50
        }
        self[ 'desolvmap' ] = {
            'keyword' : 'desolvmap' ,
            'default' : '',
            'comment' : "desolvation map",
            'value'   : ''
        }
        self[ 'dihe0' ] = {
            'keyword' : 'dihe0' ,
            'default' : 'random',
            'comment' : "initial dihedrals (relative) or random",
            'value'   : 'random'
        }
        self[ 'dihrf' ] = {
            'keyword' : 'dihrf' ,
            'default' : 1.0,
            'comment' : "per cycle reduction factor for dihedrals",
            'value'   : 1.0
        }
        self[ 'do_global_only' ] = {
            'keyword' : 'do_global_only' ,
            'default' : 50,
            'comment' : "do this many GA runs",
            'value'   : 50
        }
        self[ 'do_local_only' ] = {
            'keyword' : 'do_local_only' ,
            'default' : 50,
            'comment' : "do this many LS runs",
            'value'   : 50
        }
        self[ 'dstep' ] = {
            'keyword' : 'dstep' ,
            'default' : 50.0,
            'comment' : "torsion step/deg",
            'value'   : 50.0
        }
        self[ 'elecmap' ] = {
            'keyword' : 'elecmap' ,
            'default' : '',
            'comment' : "electrostatics map",
            'value'   : ''
        }
        self[ 'epdb' ] = {
            'keyword' : 'epdb' ,
            'default' : "",
            'comment' : "small molecule to be evaluated",
            'value'   : ""
        }
        self[ 'epdb_flag' ] = {
            'keyword' : 'epdb_flag' ,
            'default' : 0,
            'comment' : "whether to include epdb keyword",
            'value'   : 0  # true/false
        }
        self[ 'e0max' ] = {
            'keyword' : 'e0max' ,
            'default' : [0.0,  10000],
            'comment' : "max initial energy; max number of retries",
            'value'   : [0.0, 10000]
        }
        self[ 'extnrg' ] = {
            'keyword' : 'extnrg' ,
            'default' : 1000.0,
            'comment' : "external grid energy",
            'value'   : 1000.0
        }
        self[ 'fld' ] = {
            'keyword' : 'fld' ,
            'default' : None,
            'comment' : "grid_data_file",
            'value'   : None
        }
        self[ 'flexres_flag' ] = {
            'keyword' : 'flexres_flag' ,
            'default' : 0,
            'comment' : "whether to include flexres file",
            'value'   : 0
        }
        self[ 'flexres' ] = {
            'keyword' : 'flexres' ,
            'default' : flexres_filename,
            'comment' : "file containing flexible residues",
            'value'   : flexres_filename
        }
        self[ 'flexible_residues' ] = {
            'keyword' : 'flexible_residues' ,
            'default' : flexres_filename,
            'comment' : "file containing flexible residues",
            'value'   : flexres_filename
        }
        self[ 'fmap' ] = {
            'keyword' : 'fmap' ,
            'default' : '',
            'comment' : "floating map",
            'value'   : ''
        }
        self[ 'ga_cauchy_alpha' ] = {
            'keyword' : 'ga_cauchy_alpha' ,
            'default' : 0.0,
            'comment' : "Alpha parameter of Cauchy distribution",
            'value'   : 0.0
        }
        self[ 'ga_cauchy_beta' ] = {
            'keyword' : 'ga_cauchy_beta' ,
            'default' : 1.0,
            'comment' : "Beta parameter Cauchy distribution",
            'value'   : 1.0
        }
        self[ 'ga_crossover_rate' ] = {
            'keyword' : 'ga_crossover_rate' ,
            'default' : 0.80,
            'comment' : "rate of crossover",
            'value'   : 0.80
        }
        self[ 'ga_crossover_mode_flag' ] = {
            'keyword' : 'ga_crossover_mode_flag' ,
            'default' : 0,
            'comment' : "mode of crossover",
            'value'   : 0  #False/True
        }
        self[ 'ga_crossover_mode' ] = {
            'keyword' : 'ga_crossover_mode' ,
            'default' : 'twopt',
            'comment' : 'mode of crossover',
            'value'   : 'twopt' #onept,twopt,uniform,arithmetic,branch
        }
        self[ 'ga_elitism' ] = {
            'keyword' : 'ga_elitism' ,
            'default' : 1,
            'comment' : "number of top individuals to survive to next generation",
            'value'   : 1
        }
        self[ 'ga_mutation_rate' ] = {
            'keyword' : 'ga_mutation_rate' ,
            'default' : 0.02,
            'comment' : "rate of gene mutation",
            'value'   : 0.02
        }
        self[ 'ga_num_evals' ] = {
            'keyword' : 'ga_num_evals' ,
            'default' : 2500000,
            'comment' : "maximum number of energy evaluations",
            'value'   : 2500000
        }
        self[ 'ga_num_generations' ] = {
            'keyword' : 'ga_num_generations' ,
            'default' : 27000,
            'comment' : "maximum number of generations",
            'value'   : 27000
        }
        self[ 'ga_pop_size' ] = {
            'keyword' : 'ga_pop_size' ,
            'default' : 150,
            'comment' : "number of individuals in population",
            'value'   : 150
        }
        self[ 'ga_run' ] = {
            'keyword' : 'ga_run' ,
            'default' : 10,
            'comment' : "do this many hybrid GA-LS runs",
            'value'   : 10
        }
        self[ 'ga_window_size' ] = {
            'keyword' : 'ga_window_size' ,
            'default' : 10,
            'comment' : '',
            'value'   : 10
        }
        self[ 'include_1_4_interactions' ] = {
            'keyword' : 'include_1_4_interactions' ,
            'default' : 1.0,
            'comment' : "include internal 1-4 interactions",
            'value'   : 1.0  # weight
        }
        self[ 'include_1_4_interactions_flag' ] = {
            'keyword' : 'include_1_4_interactions_flag' ,
            'default' : 0,
            'comment' : "whether to include internal 1-4 interactions",
            'value'   : 0  # true/false
        }
        self[ 'intelec' ] = {
            'keyword' : 'intelec' ,
            'default' : 0,
            'comment' : "calculate internal electrostatics",
            'value'   : 0  # true/false
        }
        self[ 'intelec4' ] = {
            'keyword' : 'intelec4' ,
            'default' : 0,
            'comment' : "calculate internal electrostatics",
            'value'   : '0.1465'  # from latest force field
        }
        self[ 'intnbp_r_eps' ] = {
            'keyword' : 'intnbp_r_eps' ,
            'default' : None,
            'comment' : '', # to be set at write-time
            'value'   : None,
        }
        self[ 'ligand_types' ] = {
            'keyword' : 'ligand_types' ,
            'default' : 'C A HD OA NA',
            'comment' : "atoms types in ligand",
            'value' : 'C A HD OA NA',
        }
        self[ 'linear_schedule' ] = {
            'keyword' : 'linear_schedule' ,
            'default' : 1,
            'comment' : "use linear, arithmetic temperature reduction",
            'value'   : 1  # true/false
        }
        self[ 'ls_search_freq' ] = {
            'keyword' : 'ls_search_freq' ,
            'default' : 0.06,
            'comment' : "probability of performing local search on individual",
            'value'   : 0.06
        }
        self[ 'map' ] = {
            'keyword' : 'map' ,
            'default' : '',
            'comment' : "atom-specific affinity map",
            'value'   : ''
        }
        self[ 'move' ] = {
            'keyword' : 'move' ,
            'default' : ligand_filename,
            'comment' : "small molecule",
            'value'   : ligand_filename
        }
        self[ 'ndihe' ] = {
            'keyword' : 'ndihe' ,
            'default' : 0,
            'comment' : "number of active torsions",
            'value'   : 0
        }
        self[ 'outlev' ] = {
            'keyword' : 'outlev' ,
            'default' : 1,
            'comment' : "diagnostic output level",
            'value'   : 1
        }
        self[ 'output_pop_file' ] = {
            'keyword' : 'output_pop_file' ,
            'default' : "gen_state.txt",
            'comment' : "generation state file",
            'value'   : "gen_state.txt"
        }
        self[ 'parameter_file' ] = {
            'keyword' : 'parameter_file' ,
            'default' : 'AD4.1_bound.dat',
            'comment' : "parameter library filename ",
            'value'   : 'AD4.1_bound.dat',
        }
        self[ 'psw_trans_scale' ] = {
            'keyword' : 'psw_trans_scale' ,
            'default' : 1,
            'comment' : "pseudo sw translation rho scale",
            'value'   : 1,
        }
        self[ 'psw_rot_scale' ] = {
            'keyword' : 'psw_rot_scale' ,
            'default' : 0.05,
            'comment' : "pseudo sw rotation rho scale",
            'value'   : 0.05,
        }
        self[ 'psw_tors_scale' ] = {
            'keyword' : 'psw_tors_scale' ,
            'default' : 0.1,
            'comment' : "pseudo sw torsion rho scale",
            'value'   : 0.1,
        }
        self[ 'qstep' ] = {
            'keyword' : 'qstep' ,
            'default' : 50.0,
            'comment' : "quaternion step/deg",
            'value'   : 50.0
        }
        self[ 'quarf' ] = {
            'keyword' : 'quarf' ,
            'default' : 1.0,
            'comment' : "per cycle reduction factor for quaternions",
            'value'   : 1.0
        }
        self[ 'quat0' ] = {
            'keyword' : 'quat0' ,
            'default' : 'random',
            'comment' : "initial quaternion",
            'value'   : 'random'
        }
        self[ 'quaternion0' ] = {
            'keyword' : 'quaternion0' ,
            'default' : 'random',
            'comment' : "initial quaternion",
            'value'   : 'random'
        }
        self[ 'rejs' ] = {
            'keyword' : 'rejs' ,
            'default' : 100,
            'comment' : "maximum number of rejected steps per cycle",
            'value'   : 100
        }
        self[ 'reorient_flag' ] = {
            'keyword' : 'reorient_flag' ,
            'default' : 0, #do not routinely include reorient
            'comment' : "whether to reorient ligand at start of each run and how",
            'value'   : 0,
        }
        self[ 'reorient' ] = {
            'keyword' : 'reorient' ,
            'default' : 'random',
            'comment' : "initial orientation of ligand",
            'value'   : 'random'
        }
        self[ 'rmsatoms' ] = {
            'keyword' : 'rmsatoms' ,
            'default' : 'ligand_only', #other choice is 'all'
            'comment' : "cluster reference pdbqt file",
            'value'   : 'ligand_only',
        }
        self[ 'rmsatoms_flag' ] = {
            'keyword' : 'rmsatoms_flag' ,
            'default' : 0, #do not routinely include this keyword
            'comment' : "what to use for cluster reference pdbqt file",
            'value'   : 0,
        }
        self[ 'rmsref' ] = {
            'keyword' : 'rmsref' ,
            'default' : os.path.basename(ligand_filename),
            'comment' : "cluster reference file",
            'value'   : os.path.basename(ligand_filename)
        }
        self[ 'rmsref_flag' ] = {
            'keyword' : 'rmsref_flag' ,
            'default' : 0,
            'comment' : "whether rmsref is present",
            'value'   : 0,
        }
        self[ 'rmstol' ] = {
            'keyword' : 'rmstol' ,
            'default' : 2.0,
            'comment' : "cluster_tolerance/A",
            'value'   : 2.0
        }
        self[ 'rt0' ] = {
            'keyword' : 'rt0' ,
            'default' : 1000.0,
            'comment' : "initial annealing temperature (times gas constant)",
            'value'   : 1000.0
        }
        self[ 'rtrf' ] = {
            'keyword' : 'rtrf' ,
            'default' : 0.95,
            'comment' : "annealing temperature reduction factor",
            'value'   : 0.95
        }
        self[ 'runs' ] = {
            'keyword' : 'runs' ,
            'default' : 10,
            'comment' : "",
            'value'   : 10
        }
        self[ 'seed' ] = {
            'keyword' : 'seed' ,
            'default' : ['pid',  'time'],
            'comment' : "seeds for random generator",
            'value'   : ['pid', 'time']
        }
        self[ 'select' ] = {
            'keyword' : 'select' ,
            'default' : 'm',
            'comment' : "state selection flag: (m)inimum or (l)ast state",
            'value'   : 'm'
        }
        self[ 'set_ga' ] = {
            'keyword' : 'set_ga' ,
            'default' : 1,
            'comment' : "set the above parameters for GA or LGA",
            'value'   : 1  # true/false
        }
        self[ 'set_sw1_flag' ] = {
            'keyword' : 'set_sw1_flag' ,
            'default' : 0,
            'comment' : "set ls to sw",
            'value'   : 0  # true/false
        }
        self[ 'set_sw1' ] = {
            'keyword' : 'set_sw1' ,
            'default' : 0,
            'comment' : "set the above Solis & Wets parameters",
            'value'   : 0  # true/false
        }
        self[ 'set_psw1_flag' ] = {
            'keyword' : 'set_psw1_flag' ,
            'default' : 1,
            'comment' : "set ls to psw",
            'value'   : 1  # true/false
        }
        self[ 'set_psw1' ] = {
            'keyword' : 'set_psw1' ,
            'default' : 1,
            'comment' : "set the above pseudo-Solis & Wets parameters",
            'value'   : 1  # true/false
        }
        self[ 'simanneal' ] = {
            'keyword' : 'simanneal' ,
            'default' : 1,
            'comment' : "do as many SA runs as set by runs keyword above",
            'value'   : 1 # true/false
        }
        self[ 'sw_lb_rho' ] = {
            'keyword' : 'sw_lb_rho' ,
            'default' : 0.01,
            'comment' : "lower bound on rho",
            'value'   : 0.01
        }
        self[ 'sw_max_fail' ] = {
            'keyword' : 'sw_max_fail' ,
            'default' : 4,
            'comment' : "consecutive failures before changing rho",
            'value'   : 4
        }
        self[ 'sw_max_its' ] = {
            'keyword' : 'sw_max_its' ,
            'default' : 300,
            'comment' : "iterations of Solis & Wets local search",
            'value'   : 300
        }
        self[ 'sw_max_succ' ] = {
            'keyword' : 'sw_max_succ' ,
            'default' : 4,
            'comment' : "consecutive successes before changing rho",
            'value'   : 4
        }
        self[ 'sw_rho' ] = {
            'keyword' : 'sw_rho' ,
            'default' : 1.0,
            'comment' : "size of local search space to sample",
            'value'   : 1.0
        }
        self[ 'torsdof' ] = {
            'keyword' : 'torsdof' ,
            'default' : [0,  0.3113],
            'comment' : "torsional degrees of freedom and coeffiecent",
            'value'   : [0, 0.3113]
        }
        self[ 'torsdof4' ] = {
            'keyword' : 'torsdof4' ,
            'default' : [0],
            'comment' : "torsional degrees of freedom",
            'value'   : [0]
        }
        self[ 'tran0' ] = {
            'keyword' : 'tran0' ,
            'default' : 'random',
            'comment' : "initial coordinates/A or random",
            'value'   : 'random'
        }
        self[ 'trnrf' ] = {
            'keyword' : 'trnrf' ,
            'default' : 1.0,
            'comment' : "per cycle reduction factor for translation",
            'value'   : 1.0
        }
        self[ 'tstep' ] = {
            'keyword' : 'tstep' ,
            'default' : [2.0],
            'comment' : "translation step/A",
            'value'   : [2.0]
        }
        self[ 'types' ] = {
            'keyword' : 'types' ,
            'default' : 'ACONSH',
            'comment' : "atom type names",
            'value'   : 'ACONSH'
        }
        self[ 'unbound' ] = {            #4.2
            'keyword' : 'unbound' ,
            'default' : 0.0,
            'comment' : "free energy of ligand's unbound state ",
            'value'   : 0.0
        }
        self[ 'unbound_flag' ] = {
            'keyword' : 'unbound_flag' ,
            'default' : 0,
            'comment' : "whether to include a custom unbound value",
            'value'   : 0
        }
        self[ 'unbound_energy' ] = {      #4.2
            'keyword' : 'unbound_energy' ,
            'default' : 0.0,
            'comment' : "free energy of ligand's unbound state ",
            'value'   : 0.0
        }
        self[ 'unbound_energy_flag' ] = {
            'keyword' : 'unbound_energy_flag' ,
            'default' : 0,
            'comment' : "whether to include a custom unbound_energy value",
            'value'   : 0
        }
        self[ 'unbound_intnbp_coeffs' ] = {
            'keyword' : 'unbound_intnbp_coeffs' ,
            'default' : "20. 0. 1 2",
            'comment' : "intnbp ff coeffs",
            'value'   : "20. 0. 1 2"  # two floats and two integers: 
        }
        self[ 'unbound_intnbp_coeffs_flag' ] = {
            'keyword' : 'unbound_intnbp_coeffs_flag' ,
            'default' : 0,
            'comment' : "whether to include custom intnbp values",
            'value'   : 0
        }

        self[ 'unbound_model' ] = {
            'keyword' : 'unbound_model' ,
            'default' : "bound", # 4.1 default
            'comment' : "state of unbound ligand", 
            'value'   : "bound", #possible values: bound, extended (AD4.0), compact(n/a)
        }
        self[ 'unbound_model_flag' ] = {
            'keyword' : 'unbound_model_flag' ,
            'default' : 0, #
            'comment' : "whether to include unbound_model keyword",
            'value'   : 0, #
        }

        self[ 'autodock_parameter_version' ] = {
            'keyword' : 'autodock_parameter_version' ,
            'default' : "4.2",
            'comment' : "used by autodock to validate parameter set",
            'value'   : "4.2"
        }

        self[ 'write_all' ] = {
            'keyword' : 'write_all' ,
            'default' : "",
            'comment' : "write all conformations in a cluster",
            'value'   : "" # true/false
        }
        self[ 'write_all_flag' ] = {
            'keyword' : 'write_all_flag' ,
            'default' : 0,
            'comment' : "whether to include the write all keyword",
            'value'   : 0 # true/false
        }
        # end dictionary

        self.set_receptor(receptor_filename) # also sets self.receptor_stem
        self.set_ligand(ligand_filename)
        self.boolean_param_list = [
            'analysis' ,
            'epdb_flag',
            'flexres_flag',
            'compute_unbound_extended_flag',
            'include_1_4_interactions_flag',
            'intelec',
            'linear_schedule' ,
            'reorient_flag',
            'rmsatoms_flag',
            'set_ga' ,
            'set_sw1' ,
            'set_psw1' ,
            'simanneal' ,
            'unbound_flag',
            'unbound_energy_flag',
            'unbound_model_flag',
            'unbound_intnbp_coeffs_flag',
            'write_all',
            'write_all_flag'
            ]
        # end __init__


    def set_version(self, version):
        if version in ['3.05', '4.0', '4.1', '4.2']:
            self['autodock_parameter_version']['value'] = version
        else:
            print version, " is not valid. Valid autodock versions are '3.05', '4.1', '4.2']"


    def set_ligand(self, ligand_filename):
        self.ligand_filename = os.path.basename(ligand_filename)
        #self.ligand_filename = ligand_filename
        basename = os.path.basename(ligand_filename)
        self['rmsref']['value'] = basename


    def set_receptor(self, receptor_filename):
        self.receptor_filename = os.path.basename(receptor_filename)
        #self.receptor_filename = receptor_filename
        basename = os.path.basename(receptor_filename)
        self.receptor_stem = basename[:string.rfind(basename, '.')]
        if receptor_filename!='':
            self['fld']['value'] = self.receptor_stem + '.maps.fld'


    def set_ligand_types_from_filename(self, ligand_filename):
        ligand = Read(ligand_filename)[0]
        if ligand is None:
            print "unable to read ", ligand_filename
            return "ERROR"
        d = {}
        for a in ligand.allAtoms:
            d[a.autodock_element] = 1
        keys = d.keys()
        keys.sort()
        type_str = keys[0]
        for t in keys[1:]:
            type_str = type_str + ' ' + t + ' '
        #eg "C A NA N OA SA HD"
        self['ligand_types']['value'] = type_str
        self[ 'torsdof4' ]['value'][0] = ligand.TORSDOF
        #self[ 'torsdof4' ]['value'][0] = ligand.ndihe


    def set_ligand_types_from_filename_v3(self, ligand_filename):
        ligand = Read(ligand_filename)
        if ligand is None:
            print "unable to read ", ligand_filename
            return "ERROR"
        d = {}
        for a in ligand.allAtoms:
            d[a.autodock_element] = 1
        keys = d.keys()
        keys.sort()
        type_str = keys[0]
        for t in keys[1:]:
            type_str = type_str + t 
        #eg "CANOSH"
        self['types']['value'] = type_str
        self[ 'torsdof' ]['value'][0] = ligand.TORSDOF


    def set_ligand_types3_from_ligand_types(self, ligand_types):
        d = {}
        for t in ligand_types:
            if len(t)==1:
                d[t] = 1
            elif t[1] in ['A','D']: #AD4 special cases: NA,SA,OA,HD
                d[t[0]] = 1
            elif t in ['Cl','CL','cl']:  #AD3 special case: chlorine
                d['c'] = 1
            elif t in ['Br','BR','br']:  #AD3 special case: bromine
                d['b'] = 1
            elif t in ['Fe','FE','fe']:  #AD3 special case: iron
                d['f'] = 1
            else:
                print "unrecognized ligand_atom_type:", t
        all_types = d.keys()
        all_types.sort()
        type_str = all_types[0]
        for t in all_types[1:]:
            type_str = type_str + t
        self['types']['value'] = type_str


    #
    # read methods
    #
    def read(self, filename):
        """Read lines from the file and call _parse to set current state.
        """
        self.dpf_filename = os.path.basename(filename)
        self.dpf_written_filename = filename
        dpf_ptr = open(filename)
        lines = dpf_ptr.readlines()
        dpf_ptr.close()
        self._parse(lines)
        

    def _parse(self, lines):
        """set the current state according to lines.
        """
        self.file_params = []
        keys = self.keys()
        found_compute_unbound_extended = False
        found_include_1_4_interactions = False
        found_reorient = False
        for line in lines:
            p = ''
            words = string.split(string.replace(line, '\t', ' '))
            if words!=[] and words[0][0]!='#':
                p = words[0]
                if p not in keys:
                    continue
                # maintain a list of the parameters read from the file
                if self.file_params==[] or p!=self.file_params[-1]:
                    self.file_params.append(p)
                # parse the line
                l = len(words)
                for i in range(l):
                    if words[i][0]=='#':
                        l = i
                values = words[1:l]
                if p=='include_1_4_interactions':
                    self['include_1_4_interactions']['value'] = self._get_val(values[0])
                    self['include_1_4_interactions_flag']['value'] = True
                    found_include_1_4_interactions = True
                elif p=='set_sw1':
                    self[p]['value'] = 1
                    self['set_sw1_flag']['value'] = 1
                    self['set_psw1']['value'] = 0
                    self['set_psw1_flag']['value'] = 0
                elif p=='set_psw1':
                    self[p]['value'] = 1
                    self['set_psw1_flag']['value'] = 1
                    self['set_sw1_flag']['value'] = 0
                    self['set_sw1']['value'] = 0
                elif p=='write_all':
                    self['write_all_flag']['value'] = True
                elif p in self.boolean_param_list:
                    self[p]['value'] = 1
                elif p=='ligand_types':
                    self['ligand_types']['value'] = string.join(words[1:l])
                elif p=='types':
                    self['types']['value'] = string.join(words[1:l])
                    self['autodock_parameter_version']['value'] = "3.05"
                elif p=='rmsref':
                    self['rmsref']['value'] = self._get_val(values[0])
                    self['rmsref_flag']['value'] = True
                elif p=='reorient':
                    self['reorient']['value'] = self._get_val(values[0])
                    self['reorient_flag']['value'] = True
                    found_reorient = True
                elif p=='epdb':
                    if len(values)>0:
                        self['epdb']['value'] = self._get_val(values[0])
                    self['epdb_flag']['value'] = True
                elif p=='compute_unbound_extended':
                    self['compute_unbound_extended_flag']['value'] = True
                    found_compute_unbound_extended = True
                elif p=='ga_crossover_mode':
                    self['ga_crossover_mode']['value'] = self._get_val(values[0])
                    self['ga_crossover_mode_flag']['value'] = True
                elif p=='unbound':
                    self['unbound_flag']['value'] = True
                    self['unbound']['value'] = self._get_val(values[0])
                elif p=='unbound_intnbp_coeffs':
                    self['unbound_intnbp_coeffs_flag']['value'] = True
                    self['unbound_intnbp_coeffs']['value'] = string.join(words[1:l])
                elif p=='rmsatoms' and values[0]=='all':
                    self['rmsatoms_flag']['value'] = True
                    self['rmsatoms']['value'] = values[0]  # this should be 'ligand_only' or'all'
                elif p=='flexres':
                    self['flexres']['value'] = self._get_val(values[0])
                    self['flexres_flag']['value'] = True
                elif p=='parameter_file':
                    self['custom_parameter_file']['value'] = 1
                    self['parameter_file']['value'] = values[0]
                elif ((len(values)==1) and
                      (type(self[p]['default'])!=types.ListType)):
                    self[p]['value'] = self._get_val(values[0])
                elif self.has_key(p):
                    self[p]['value'] = []
                    for v in values:
                        self[p]['value'].append( self._get_val(v))
                    if p=='quat0':  
                        self['axisangle0']['value'] = self[p]['value']
                else:
                    print 'WARNING: unknown keyword=', p
            if p=='fld':
                # so words[1] ends .maps.fld
                ind = string.index(words[1], 'maps') - 1
                self.receptor_stem = words[1][:ind]
            elif p=='map':
                #problem in this case is that elements may be 1 or 2 char
                # so name could end .C.map or .CL.map
                name = words[1]
                if name[-6]=='.':
                    self.receptor_stem = name[:-6]
                elif name[-7]=='.':
                    self.receptor_stem = name[:-7]
        if not found_compute_unbound_extended:
            self['compute_unbound_extended_flag']['value'] = 0
        if not found_include_1_4_interactions:
            self['include_1_4_interactions_flag']['value'] = 0
        if not found_reorient:
            self['reorient_flag']['value'] = 0


    def _get_val(self, val_str):
        try:
            return int(val_str)
        except ValueError:
            pass
        try:
            return float(val_str)
        except ValueError:
            pass
        if type(val_str)==types.StringType:
            return val_str
        else:
            raise NotImplementedError, "value: %s of unsupport type %s" % (val_str, type(val_str).__name__)

    #
    # write methods
    #
    def write(self, filename, param_list):
        """Write the current state to a file

        file is a pointer to a writeable file
        param_list is a list of parameter strings.
        For best results use the parameter_lists supplied by this class.
        """
        if filename=='':
            dpf_ptr = sys.stdout
            self.dpf_filename = ''
            self.dpf_written_filename = None
        else:
            dpf_ptr = open(filename, 'w')
            self.dpf_filename = os.path.basename(filename)
            self.dpf_written_filename = filename
            self.file_params = param_list

        for p in param_list:
            # maps are a special case
            if p=='map':
                for a in self['types']['value']:
                    dpf_ptr.write( self.make_map_string(p, a) )
                # write the electrostatics map; kluge comment
                tmp = self[p]['comment']
                self[p]['comment'] = "electrostatics map"
                dpf_ptr.write( self.make_map_string(p, 'e') )
                self[p]['comment'] = tmp
            elif p=='fmap':
                dpf_ptr.write( self.make_map_string(p, 'f'))
            # intnbp_r_eps is a special case
            elif p=='intnbp_r_eps':
                atoms = self['types']['value']
                #strings don't have index method in python1.5.2
                #for a1 in atoms:
                    #for a2 in atoms[atoms.index(a1):]:
                        #dpf_ptr.write( self.make_intnbp_r_eps_string(a1, a2))
                lenTypes = len(atoms)
                for i in range(lenTypes):
                    for j in range(i, lenTypes):
                        a1 = atoms[i]
                        a2 = atoms[j]
                        dpf_ptr.write( self.make_intnbp_r_eps_string(a1, a2))

            elif p=='intelec' and self[p]['value']:
                dpf_ptr.write('intelec 0.1146                       # calculate internal electrostatics\n')
            elif p=='rmsref_flag':
                if self['rmsref_flag']['value']:
                    dpf_ptr.write(self.make_param_string('rmsref'))
            elif p=='rmsref':
                pass
            elif p=='set_sw1': 
                #print "in write with set_sw1:", self[p]['value']
                if int(self[p]['value'])==0:
                    self['set_psw1']['value']=1
                    dpf_ptr.write( self.make_param_string('set_psw1'))
                else:
                    self[p]['value']=1
                    dpf_ptr.write( self.make_param_string(p))
            # all the other parameters handle themselves
            elif p=='set_psw1': 
                #print "in write with set_psw1:", self[p]['value']
                if int(self[p]['value'])==0:
                    self['set_sw1']['value']=1
                    dpf_ptr.write( self.make_param_string('set_sw1'))
                else:
                    self[p]['value']=1
                    dpf_ptr.write( self.make_param_string(p))
            # all the other parameters handle themselves
            else:
                dpf_ptr.write( self.make_param_string(p))
        if dpf_ptr!=sys.stdout:
            dpf_ptr.close()


    def make_param_string(self, param):
        """return the output string for the given param using the value
           and comment entries in it's dictionary.
        """
        p = self[param]
        vt = type(p['value'])
        if param in self.boolean_param_list:
            if not p['value']:
                return "#\n"
            else:
                val_str = ""
        elif ((vt==types.IntType) or
              (vt==types.LongType) or
              (vt==types.FloatType) or
              (vt==types.StringType)):
            val_str = str(p['value'])
        elif ((vt==types.ListType) or
              (vt==types.TupleType)):
            val_str = ""
            for v in p['value']:
                val_str = val_str + str(v) + " "
        else:
            raise NotImplementedError, "type (%s) of parameter %s unsupported" % (vt.__name__, param)
        return self._make_string(p, val_str)


    def make_intnbp_r_eps_string(self, atom1, atom2):
        p = self[ 'intnbp_r_eps' ]
        index = "lj" + atom1 + atom2
        
        val_str = "%5.2f %9.7f 12 6" % (Rij[index], epsij[index])
        p['comment'] = "%s-%s lj" % (atom1, atom2)
        return self._make_string(p, val_str)


    def make_map_string(self, param, type):
        p = self[param]
        val_str = self.receptor_stem + ".%s.map" % (type)
        return self._make_string(p, val_str)
    

    def _make_string(self, p, val_str):
        return "%s %s%s# %s\n" % (p['keyword'],
                                  val_str,
                                  " "*(36 -(len(p['keyword'])+len(val_str))),
                                  p['comment'])

    def write4(self, filename, param_list):
        """Write the current state to an AutoDock4 dpf file

        file is a pointer to a writeable file
        param_list is a list of parameter strings.
        For best results use the parameter_lists supplied by this class.
        """
        if filename=='':
            dpf_ptr = sys.stdout
            self.dpf_filename = ''
            self.dpf_written_filename = None
        else:
            dpf_ptr = open(filename, 'w')
            self.dpf_filename = os.path.basename(filename)
            self.dpf_written_filename = filename
            self.file_params = param_list

        #to write dpf4, set unbound_model to extended OR unbound_model extended float 
        for p in param_list:
            if p=='autodock_parameter_version':
                oldval = self['autodock_parameter_version']['value']
                self['autodock_parameter_version']['value']=4.2
                dpf_ptr.write( self.make_param_string('autodock_parameter_version'))
                #NEW 4/1/2009
                #also set unbound_model to extended
                self['unbound_model']['value'] = "extended"
                self['unbound_model_flag']['value'] = 1
                self['autodock_parameter_version']['value']=oldval
            elif p=='custom_parameter_file':
                if self['custom_parameter_file']['value']:
                    dpf_ptr.write( self.make_param_string('parameter_file'))
            elif p=='map':
                # maps are a special case
                for a in string.split(self['ligand_types']['value']):
                    dpf_ptr.write( self.make_map_string(p, a) )
                # write the electrostatics map
                dpf_ptr.write( self.make_map_string('elecmap','e') )
                # write the desolvation map
                dpf_ptr.write( self.make_map_string('desolvmap', 'd') )
            elif p=='reorient_flag':
                if self['reorient_flag']['value']:
                    dpf_ptr.write(self.make_param_string('reorient'))
            elif p=='reorient':
                pass
            elif p=='set_psw1' or p=='set_sw1':
                if self['set_psw1_flag']['value']:
                    dpf_ptr.write( self.make_param_string(p) )
                elif self['set_sw1_flag']['value']:
                    dpf_ptr.write( self.make_param_string('set_sw1') )
                else:
                    pass
            elif p=='fmap':
                self[p]['comment'] = "floating point map"
                dpf_ptr.write( self.make_map_string(p, 'f'))
            elif p=='include_1_4_interactions_flag':
                if self['include_1_4_interactions_flag']['value']:
                    dpf_ptr.write(self.make_param_string('include_1_4_interactions'))
            elif p=='include_1_4_interactions':
                pass
            #not support 4/1/2009->
            #elif p=='compute_unbound_extended_flag':
            #    if self['compute_unbound_extended_flag']['value']:
            #        dpf_ptr.write(self.make_param_string('compute_unbound_extended'))
            #elif p=='compute_unbound_extended':
            #    pass
            elif p=='unbound_model_flag':
                self['unbound_model']['value'] = "extended"
                dpf_ptr.write(self.make_param_string('unbound_model'))
            elif p=='unbound_model':
                pass
            #elif p=='unbound_model_flag':
            #    if self['unbound_model_flag']['value']:
            #        dpf_ptr.write(self.make_param_string('unbound model'))
            #elif p=='unbound_model':
            #    pass
            elif p=='unbound_flag':
                if self['unbound_flag']['value']:
                    dpf_ptr.write(self.make_param_string('unbound'))
            elif p=='unbound':
                pass
            elif p=='unbound_intnbp_coeffs_flag':
                if self['unbound_intnbp_coeffs_flag']['value']:
                    dpf_ptr.write(self.make_param_string('unbound_intnbp_coeffs'))
            elif p=='ga_crossover_mode_flag':
                if self['ga_crossover_mode_flag']['value']:
                    dpf_ptr.write(self.make_param_string('ga_crossover_mode'))
            elif p=='ga_crossover_mode':
                pass
            elif p=='unbound_intnbp_coeffs':
                pass
            elif p=='rmsatoms_flag':
                if self['rmsatoms_flag']['value'] and self['rmsatoms']['value']=='all':
                    dpf_ptr.write(self.make_param_string('rmsatoms'))
            elif p=='rmsatoms':
                pass
            elif p=='flexres_flag':
                if self['flexres_flag']['value']:
                    dpf_ptr.write("flexres %s                  # file containing flexible residues\n" %self['flexres']['value'])
            elif p=='flexres':
                pass
            elif p=='unbound_model_flag':
                if self['unbound_model_flag']['value']>0:
                    dpf_ptr.write(self.make_param_string('unbound_model'))
            elif p=='write_all_flag':
                if self['write_all_flag']['value']:
                    dpf_ptr.write("write_all                  # write all conformations in a cluster\n" )
            elif p=='write_all':
                pass
            elif p=='epdb_flag':
                if self['epdb_flag']['value']:
                    if self['epdb']['value']=="":
                        dpf_ptr.write("epdb %s                  # small molecule to be evaluated\n" %self.ligand_filename)
                    else:
                        dpf_ptr.write("epdb %s                  # small molecule to be evaluated\n" %self['epdb']['value'])
            elif p=='epdb':
                pass
            elif p=='rmsref_flag':
                flag = self['rmsref_flag']['value']
                if type(flag)==type(''):
                    flag=eval(flag)
                    self['rmsref_flag']['value']=flag
                if self['rmsref_flag']['value']:
                    if self['rmsref']['value']=="":
                        dpf_ptr.write("rmsref %s                  # reference ligand conformation\n" %self.ligand_filename)
                    else:
                        dpf_ptr.write("rmsref %s                  # reference ligand conformation\n" %self['rmsref']['value'])
            elif p=='rmsref':
                pass
            elif p=='torsdof4':
                dpf_ptr.write('torsdof %d                            # torsional degrees of freedom\n' \
                            %(self['torsdof4']['value'][0]))
            elif p=='intelec':  #always include internal electrostatics
                dpf_ptr.write('intelec                              # calculate internal electrostatics\n')
            # all the other parameters handle themselves
            else:
                dpf_ptr.write( self.make_param_string(p))
        if dpf_ptr!=sys.stdout:
            dpf_ptr.close()


    def write41(self, filename, param_list):
        """Write the current state to an AutoDock4.1 dpf file
        file is a pointer to a writeable file
        param_list is a list of parameter strings.
        For best results use the parameter_lists supplied by this class.
        """
        if filename=='':
            dpf_ptr = sys.stdout
            self.dpf_filename = ''
            self.dpf_written_filename = None
        else:
            dpf_ptr = open(filename, 'w')
            self.dpf_filename = os.path.basename(filename)
            self.dpf_written_filename = filename
            self.file_params = param_list

        
        for p in param_list:
            if p=='custom_parameter_file':
                self['custom_parameter_file']['value'] = 1
                #self['parameter_file']['value'] = 'AD4.1_bound.dat'
                if 'parameter_file' not in param_list:
                    oldval = self['parameter_file']['value']
                    self['parameter_file']['value'] = 'AD4.1_bound.dat'
                    dpf_ptr.write( self.make_param_string('parameter_file'))
                    self['parameter_file']['value'] = oldval
            elif p=='map':
                # maps are a special case
                for a in string.split(self['ligand_types']['value']):
                    dpf_ptr.write( self.make_map_string(p, a) )
                # write the electrostatics map
                dpf_ptr.write( self.make_map_string('elecmap','e') )
                # write the desolvation map
                dpf_ptr.write( self.make_map_string('desolvmap', 'd') )
            elif p=='reorient_flag':
                if self['reorient_flag']['value']:
                    dpf_ptr.write(self.make_param_string('reorient'))
            elif p=='reorient':
                pass
            elif p=='set_psw1' or p=='set_sw1':
                if self['set_psw1_flag']['value']:
                    dpf_ptr.write( self.make_param_string(p) )
                elif self['set_sw1_flag']['value']:
                    dpf_ptr.write( self.make_param_string('set_sw1') )
                else:
                    pass
            elif p=='fmap':
                self[p]['comment'] = "floating point map"
                dpf_ptr.write( self.make_map_string(p, 'f'))
            elif p=='include_1_4_interactions_flag':
                if self['include_1_4_interactions_flag']['value']:
                    dpf_ptr.write(self.make_param_string('include_1_4_interactions'))
            elif p=='include_1_4_interactions':
                pass
            elif p=='unbound_energy_flag':
                #IF user specifies an unbound_energy, write it
                if self['unbound_energy_flag']['value']:
                    dpf_ptr.write(self.make_param_string('unbound_energy'))
            elif p=='unbound_energy':
                pass
            elif p=='unbound_model_flag':
                self['unbound_model']['value'] = "bound"
                dpf_ptr.write(self.make_param_string('unbound_model'))
            elif p=='unbound_model':
                pass
            elif p=='compute_unbound_extended_flag':
                self['compute_unbound_extended_flag']['value']=0
            elif p=='compute_unbound_extended':
                pass
            elif p=='unbound_flag':
                self['unbound_flag']['value']=0
            elif p=='unbound':
                pass
            elif p=='unbound_intnbp_coeffs_flag':
                self['unbound_intnbp_coeffs_flag']['value']=0
            elif p=='ga_crossover_mode_flag':
                if self['ga_crossover_mode_flag']['value']:
                    dpf_ptr.write(self.make_param_string('ga_crossover_mode'))
            elif p=='ga_crossover_mode':
                pass
            elif p=='unbound_intnbp_coeffs':
                pass
            elif p=='rmsatoms_flag':
                if self['rmsatoms_flag']['value'] and self['rmsatoms']['value']=='all':
                    dpf_ptr.write(self.make_param_string('rmsatoms'))
            elif p=='rmsatoms':
                pass
            elif p=='flexres_flag':
                if self['flexres_flag']['value']:
                    dpf_ptr.write("flexres %s                  # file containing flexible residues\n" %self['flexres']['value'])
            elif p=='flexres':
                pass
            elif p=='unbound_model_flag':
                if self['unbound_model_flag']['value']>0:
                    dpf_ptr.write(self.make_param_string('unbound_model'))
            elif p=='unbound_model':
                pass
            elif p=='write_all_flag':
                if self['write_all_flag']['value']:
                    dpf_ptr.write("write_all                  # write all conformations in a cluster\n" )
            elif p=='write_all':
                pass
            elif p=='epdb_flag':
                if self['epdb_flag']['value']:
                    if self['epdb']['value']=="":
                        dpf_ptr.write("epdb %s                  # small molecule to be evaluated\n" %self.ligand_filename)
                    else:
                        dpf_ptr.write("epdb %s                  # small molecule to be evaluated\n" %self['epdb']['value'])
            elif p=='epdb':
                pass
            elif p=='rmsref_flag':
                flag = self['rmsref_flag']['value']
                if type(flag)==type(''):
                    flag=eval(flag)
                    self['rmsref_flag']['value']=flag
                if self['rmsref_flag']['value']:
                    if self['rmsref']['value']=="":
                        dpf_ptr.write("rmsref %s                  # reference ligand conformation\n" %self.ligand_filename)
                    else:
                        dpf_ptr.write("rmsref %s                  # reference ligand conformation\n" %self['rmsref']['value'])
            elif p=='rmsref':
                pass
            elif p=='torsdof4':
                dpf_ptr.write('torsdof %d                            # torsional degrees of freedom\n'  %(self['torsdof4']['value'][0]))
            elif p=='intelec':  #always include internal electrostatics
                dpf_ptr.write('intelec                              # calculate internal electrostatics\n')
            # all the other parameters handle themselves
            else:
                dpf_ptr.write( self.make_param_string(p))
        if dpf_ptr!=sys.stdout:
            dpf_ptr.close()


    def write42(self, filename, param_list):
        """Write the current state to an AutoDock4.2 dpf file
        file is a pointer to a writeable file
        param_list is a list of parameter strings.
        For best results use the parameter_lists supplied by this class.
        """
        if filename=='':
            dpf_ptr = sys.stdout
            self.dpf_filename = ''
            self.dpf_written_filename = None
        else:
            dpf_ptr = open(filename, 'w')
            self.dpf_filename = os.path.basename(filename)
            self.dpf_written_filename = filename
            self.file_params = param_list

        for p in param_list:
            if p=='custom_parameter_file':
                self['custom_parameter_file']['value'] = 1
                #self['parameter_file']['value'] = 'AD4.1_bound.dat'
                if 'parameter_file' not in param_list:
                    oldval = self['parameter_file']['value']
                    self['parameter_file']['value'] = 'AD4.1_bound.dat'
                    dpf_ptr.write( self.make_param_string('parameter_file'))
                    self['parameter_file']['value'] = oldval
            elif p=='map':
                # maps are a special case
                for a in string.split(self['ligand_types']['value']):
                    dpf_ptr.write( self.make_map_string(p, a) )
                # write the electrostatics map
                dpf_ptr.write( self.make_map_string('elecmap','e') )
                # write the desolvation map
                dpf_ptr.write( self.make_map_string('desolvmap', 'd') )
            elif p=='reorient_flag':
                if self['reorient_flag']['value']:
                    dpf_ptr.write(self.make_param_string('reorient'))
            elif p=='reorient':
                pass
            elif p=='set_psw1' or p=='set_sw1':
                if self['set_psw1_flag']['value']:
                    dpf_ptr.write( self.make_param_string(p) )
                elif self['set_sw1_flag']['value']:
                    dpf_ptr.write( self.make_param_string('set_sw1') )
                else:
                    pass
            elif p=='fmap':
                self[p]['comment'] = "floating point map"
                dpf_ptr.write( self.make_map_string(p, 'f'))
            elif p=='include_1_4_interactions_flag':
                if self['include_1_4_interactions_flag']['value']:
                    dpf_ptr.write(self.make_param_string('include_1_4_interactions'))
            elif p=='include_1_4_interactions':
                pass
            elif p=='unbound_energy_flag':
                #IF user specifies an unbound_energy, write it
                if self['unbound_energy_flag']['value']:
                    dpf_ptr.write(self.make_param_string('unbound_energy'))
            elif p=='unbound_energy':
                pass
            elif p=='unbound_model_flag':
                self['unbound_model']['value'] = "bound"
                dpf_ptr.write(self.make_param_string('unbound_model'))
            elif p=='unbound_model':
                pass
            elif p=='compute_unbound_extended_flag':
                self['compute_unbound_extended_flag']['value']=0
            elif p=='compute_unbound_extended':
                pass
            elif p=='unbound_flag':
                self['unbound_flag']['value']=0
            elif p=='unbound':
                pass
            elif p=='unbound_intnbp_coeffs_flag':
                self['unbound_intnbp_coeffs_flag']['value']=0
            elif p=='ga_crossover_mode_flag':
                if self['ga_crossover_mode_flag']['value']:
                    dpf_ptr.write(self.make_param_string('ga_crossover_mode'))
            elif p=='ga_crossover_mode':
                pass
            elif p=='unbound_intnbp_coeffs':
                pass
            elif p=='rmsatoms_flag':
                if self['rmsatoms_flag']['value'] and self['rmsatoms']['value']=='all':
                    dpf_ptr.write(self.make_param_string('rmsatoms'))
            elif p=='rmsatoms':
                pass
            elif p=='flexres_flag':
                if self['flexres_flag']['value']:
                    dpf_ptr.write("flexres %s                  # file containing flexible residues\n" %self['flexres']['value'])
            elif p=='flexres':
                pass
            elif p=='unbound_model_flag':
                if self['unbound_model_flag']['value']>0:
                    dpf_ptr.write(self.make_param_string('unbound_model'))
            elif p=='unbound_model':
                pass
            elif p=='write_all_flag':
                if self['write_all_flag']['value']:
                    dpf_ptr.write("write_all                  # write all conformations in a cluster\n" )
            elif p=='write_all':
                pass
            elif p=='epdb_flag':
                if self['epdb_flag']['value']:
                    if self['epdb']['value']=="":
                        dpf_ptr.write("epdb %s                  # small molecule to be evaluated\n" %self.ligand_filename)
                    else:
                        dpf_ptr.write("epdb %s                  # small molecule to be evaluated\n" %self['epdb']['value'])
            elif p=='epdb':
                pass
            elif p=='rmsref_flag':
                flag = self['rmsref_flag']['value']
                if type(flag)==type(''):
                    flag=eval(flag)
                    self['rmsref_flag']['value']=flag
                if self['rmsref_flag']['value']:
                    if self['rmsref']['value']=="":
                        dpf_ptr.write("rmsref %s                  # reference ligand conformation\n" %self.ligand_filename)
                    else:
                        dpf_ptr.write("rmsref %s                  # reference ligand conformation\n" %self['rmsref']['value'])
            elif p=='rmsref':
                pass
            elif p=='torsdof4':
                dpf_ptr.write('torsdof %d                            # torsional degrees of freedom\n'  %(self['torsdof4']['value'][0]))
            elif p=='intelec':  #always include internal electrostatics
                dpf_ptr.write('intelec                              # calculate internal electrostatics\n')
            # all the other parameters handle themselves
            else:
                dpf_ptr.write( self.make_param_string(p))
        if dpf_ptr!=sys.stdout:
            dpf_ptr.close()

        
# the following lists are class variables describing the keywords
# that must be output for a file of the specified type.
 
cluster_list = ['types', 'rmstol', 'cluster', 'analysis' ] #3.05

genetic_algorithm_list = [ #3.05
    'outlev', 
    'seed',
    'types',
    'fld',
    'map',
    'move',
    'about', 
    'tran0', 
    'quat0', 
    'ndihe', 
    'dihe0', 
    'tstep', 
    'qstep', 
    'dstep', 
    'torsdof', 
    'intnbp_r_eps', 
    'intelec',
    'rmstol', 
    'rmsref',
    'rmsref_flag',
    'extnrg', 
    'e0max', 
    'ga_pop_size', 
    'ga_num_evals', 
    'ga_num_generations', 
    'ga_elitism', 
    'ga_mutation_rate', 
    'ga_crossover_rate', 
    'ga_window_size', 
    'ga_cauchy_alpha', 
    'ga_cauchy_beta', 
    'set_ga', 
    'do_global_only', 
    'analysis' ]

genetic_algorithm_local_search_list = [ #3.05
    'outlev', 
    'seed', 
    'types', 
    'fld', 
    'map', 
    'move', 
    'about', 
    'tran0', 
    'quat0', 
    'ndihe', 
    'dihe0', 
    'tstep', 
    'qstep', 
    'dstep', 
    'torsdof', 
    'intnbp_r_eps', 
    'intelec',
    'rmstol', 
    'rmsref',
    'rmsref_flag',
    'extnrg', 
    'e0max', 
    'ga_pop_size', 
    'ga_num_evals', 
    'ga_num_generations', 
    'ga_elitism', 
    'ga_mutation_rate', 
    'ga_crossover_rate', 
    'ga_window_size', 
    'ga_cauchy_alpha', 
    'ga_cauchy_beta', 
    'set_ga', 
    'sw_max_its', 
    'sw_max_succ', 
    'sw_max_fail', 
    'sw_rho', 
    'sw_lb_rho', 
    'ls_search_freq', 
    'set_sw1', 
    'ga_run', 
    'analysis'
    ]

local_search_list = [ #3.05
    'outlev', 
    'seed', 
    'types', 
    'fld', 
    'map', 
    'move', 
    'about', 
    'tran0', 
    'quat0', 
    'ndihe', 
    'dihe0', 
    'tstep', 
    'qstep', 
    'dstep', 
    'torsdof', 
    'intnbp_r_eps', 
    'intelec',
    'rmstol', 
    'rmsref',
    'rmsref_flag',
    'extnrg', 
    'e0max', 
    'sw_max_its', 
    'sw_max_succ', 
    'sw_max_fail', 
    'sw_rho', 
    'sw_lb_rho', 
    'ls_search_freq', 
    'set_sw1', 
    'do_local_only', 
    'analysis'
    ]

simulated_annealing_list = [ #3.05
    'outlev', 
    'seed', 
    'types', 
    'fld', 
    'map', 
    'move', 
    'about', 
    'tran0', 
    'quat0', 
    'ndihe', 
    'dihe0', 
    'tstep', 
    'qstep', 
    'dstep', 
    'torsdof', 
    'intnbp_r_eps', 
    'intelec',
    'rmstol', 
    'rmsref',
    'rmsref_flag',
    'extnrg', 
    'e0max', 
    'rt0', 
    'linear_schedule', 
    'rtrf', 
    'trnrf', 
    'quarf', 
    'dihrf', 
    'runs', 
    'cycles', 
    'accs', 
    'rejs', 
    'select', 
    'simanneal', 
    'analysis'
    ]

docking_parameter_list = [ #3.05
    'about',
    'accs',
    'analysis',
    'cluster',
    'cycles',
    'dihe0',
    'dihrf',
    'do_global_only',
    'do_local_only',
    'dstep',
    'e0max',
    'extnrg',
    'fld',
    'fmap',
    'ga_cauchy_alpha',
    'ga_cauchy_beta',
    'ga_crossover_rate',
    'ga_elitism',
    'ga_mutation_rate',
    'ga_num_evals',
    'ga_num_generations',
    'ga_pop_size',
    'ga_run',
    'ga_window_size',
    'intnbp_r_eps',
    'intelec',
    'linear_schedule',
    'ls_search_freq',
    'map',
    'move',
    'ndihe',
    'outlev',
    'qstep',
    'quarf',
    'quat0',
    'rejs',
    'rmsref',
    'rmstol',
    'rmsref',
    'rmsref_flag',
    'rt0',
    'rtrf',
    'runs',
    'seed',
    'select',
    'set_ga',
    'set_sw1',
    'set_psw1',
    'simanneal',
    'sw_lb_rho',
    'sw_max_fail',
    'sw_max_its',
    'sw_max_succ',
    'sw_rho',
    'torsdof',
    'tran0',
    'trnrf',
    'tstep',
    'types',
    'write_all',
    'write_all_flag']


#####AUTODOCK4####
 
###???parameter_file???
cluster_list4 = ['autodock_parameter_version', 'custom_parameter_file', 'ligand_types', 'rmstol', 'cluster', 'analysis' ]

genetic_algorithm_list4 = [
    'autodock_parameter_version',
    'outlev', 
    'custom_parameter_file',        #NEW
    'include_1_4_interactions',     #NEW
    'include_1_4_interactions_flag',#NEW
    'intelec',
    'seed',
    'unbound_intnbp_coeffs',
    'unbound_intnbp_coeffs_flag',
    'ligand_types',                 #NEW
    'fld',
    'map',
    'move',
    'flexres_flag',
    'flexres',
    'about', 
    'reorient_flag',
    'tran0', 
    'axisangle0', 
    'dihe0', 
    'tstep', 
    'qstep', 
    'dstep', 
    'torsdof4', 
    'rmstol', 
    'rmsref',
    'rmsref_flag',
    'rmsatoms',
    'rmsatoms_flag',
    'extnrg', 
    'e0max', 
    'ga_pop_size', 
    'ga_num_evals', 
    'ga_num_generations', 
    'ga_elitism', 
    'ga_mutation_rate', 
    'ga_crossover_rate', 
    'ga_crossover_mode_flag', 
    'ga_window_size', 
    'ga_cauchy_alpha', 
    'ga_cauchy_beta', 
    'set_ga', 
    'unbound',                      # set to user value
    'unbound_flag',                 # 0 by default
    'unbound_model',                # set to "extended"
    'unbound_model_flag',           # 1 by default
    'epdb',                         #NEW ???
    'epdb_flag',                    #NEW ???
    #'compute_unbound_extended',
    #'compute_unbound_extended_flag',
    'do_global_only',               #why isn't ga_run here???
    'write_all',
    'write_all_flag',
    'analysis' ]


genetic_algorithm_local_search_list4 = [
    'autodock_parameter_version',
    'outlev', 
    'custom_parameter_file',        #NEW
    'include_1_4_interactions',     #NEW
    'include_1_4_interactions_flag',#NEW
    'intelec',
    'seed',
    'unbound_intnbp_coeffs',
    'unbound_intnbp_coeffs_flag',
    'ligand_types',                 #NEW
    'fld', 
    'map', 
    'move', 
    'flexres_flag',
    'flexres',
    'about', 
    'reorient_flag',
    'tran0', 
    'axisangle0', 
    'dihe0', 
    'tstep', 
    'qstep', 
    'dstep', 
    'torsdof4', 
    'rmstol', 
    'rmsref',
    'rmsref_flag',
    'rmsatoms',
    'rmsatoms_flag',
    'extnrg', 
    'e0max', 
    'ga_pop_size', 
    'ga_num_evals', 
    'ga_num_generations', 
    'ga_elitism', 
    'ga_mutation_rate', 
    'ga_crossover_rate', 
    'ga_crossover_mode_flag', 
    'ga_window_size', 
    'ga_cauchy_alpha', 
    'ga_cauchy_beta', 
    'set_ga', 
    'sw_max_its', 
    'sw_max_succ', 
    'sw_max_fail', 
    'sw_rho', 
    'sw_lb_rho', 
    'ls_search_freq', 
    'set_psw1', 
    'unbound',                      # set to user value
    'unbound_flag',                 # 0 by default
    'unbound_model',                # set to "extended"
    'unbound_model_flag',           # 1 by default
    'epdb',                         #NEW ???
    'epdb_flag',                    #NEW ???
    #'compute_unbound_extended',     
    #'compute_unbound_extended_flag',
    'ga_run', 
    'write_all',
    'write_all_flag',
    'analysis'
    ]

genetic_algorithm_local_search_list4_with_parameter_file = [
    'autodock_parameter_version',
    'outlev', 
    'parameter_file',               #NEW
    'include_1_4_interactions',     #NEW
    'include_1_4_interactions_flag',#NEW
    'intelec',
    'seed', 
    'unbound_intnbp_coeffs',
    'unbound_intnbp_coeffs_flag',
    'ligand_types',                 #NEW
    'fld', 
    'map', 
    'move', 
    'flexres_flag',
    'flexres',
    'about', 
    'reorient_flag',
    'tran0', 
    'axisangle0', 
    'dihe0', 
    'tstep', 
    'qstep', 
    'dstep', 
    'torsdof4', 
    'rmstol', 
    'rmsref',
    'rmsref_flag',
    'rmsatoms',
    'rmsatoms_flag',
    'extnrg', 
    'e0max', 
    'ga_pop_size', 
    'ga_num_evals', 
    'ga_num_generations', 
    'ga_elitism', 
    'ga_mutation_rate', 
    'ga_crossover_rate', 
    'ga_crossover_mode_flag', 
    'ga_window_size', 
    'ga_cauchy_alpha', 
    'ga_cauchy_beta', 
    'set_ga', 
    'sw_max_its', 
    'sw_max_succ', 
    'sw_max_fail', 
    'sw_rho', 
    'sw_lb_rho', 
    'ls_search_freq', 
    'set_psw1', 
    'unbound',                      # set to user value
    'unbound_flag',                 # 0 by default
    'unbound_model',                # set to "extended"
    'unbound_model_flag',           # 1 by default
    'epdb',                         #NEW ???
    'epdb_flag',                    #NEW ???
    #'compute_unbound_extended',     
    #'compute_unbound_extended_flag',
    'ga_run', 
    'write_all',
    'write_all_flag',
    'analysis'
    ]


local_search_list4 = [
    'autodock_parameter_version',
    'outlev', 
    'custom_parameter_file',        #NEW
    'include_1_4_interactions',     #NEW
    'include_1_4_interactions_flag',#NEW
    'intelec',
    'seed', 
    'unbound_intnbp_coeffs',
    'unbound_intnbp_coeffs_flag',
    'ligand_types',                 #NEW
    'fld', 
    'map', 
    'move', 
    'flexres_flag',
    'flexres',
    'about', 
    'reorient_flag',
    'tran0', 
    'axisangle0', 
    'dihe0', 
    'tstep', 
    'qstep', 
    'dstep', 
    'torsdof4', 
    'rmstol', 
    'rmsref',
    'rmsref_flag',
    'rmsatoms',
    'rmsatoms_flag',
    'extnrg', 
    'e0max', 
    'ga_pop_size', 
    'sw_max_its', 
    'sw_max_succ', 
    'sw_max_fail', 
    'sw_rho', 
    'sw_lb_rho', 
    'ls_search_freq', 
    'set_psw1', 
    'unbound',                      # set to user value
    'unbound_flag',                 # 0 by default
    'unbound_model',                # set to "extended"
    'unbound_model_flag',           # 1 by default
    'epdb',                         #NEW
    'epdb_flag',                    #NEW ???
    'do_local_only', 
    'write_all',
    'write_all_flag',
    'analysis'
    ]

simulated_annealing_list4 = [
    'autodock_parameter_version',
    'outlev', 
    'custom_parameter_file',        #NEW
    'include_1_4_interactions',     #NEW
    'include_1_4_interactions_flag',#NEW
    'intelec',
    'seed', 
    'unbound_intnbp_coeffs',
    'unbound_intnbp_coeffs_flag',
    'ligand_types',                 #NEW
    'fld', 
    'map', 
    'move', 
    'flexres_flag',
    'flexres',
    'about', 
    'reorient_flag',
    'tran0', 
    'axisangle0', 
    'dihe0', 
    'tstep', 
    'qstep', 
    'dstep', 
    'torsdof4', 
    'rmsatoms',
    'rmsatoms_flag',
    'rmstol', 
    'rmsref',
    'rmsref_flag',
    'extnrg', 
    'e0max', 
    'rt0', 
    'linear_schedule', 
    'rtrf', 
    'trnrf', 
    'quarf', 
    'dihrf', 
    'runs', 
    'cycles', 
    'accs', 
    'rejs', 
    'select', 
    'unbound',                      # set to user value
    'unbound_flag',                 # 0 by default
    'unbound_model',                # set to "extended"
    'unbound_model_flag',           # 1 by default
    'epdb',                         #NEW
    'epdb_flag',                    #NEW ???
    'simanneal', 
    'write_all',
    'write_all_flag',
    'analysis'
    ]

docking_parameter_list4 = [
    'autodock_parameter_version',
    'about',
    'accs',
    'analysis',
    'cluster',
    'compute_unbound_extended',
    'compute_unbound_extended_flag',#0 by default
    'cycles',
    'dihe0',
    'dihrf',
    'do_global_only',
    'do_local_only',
    'dstep',
    'e0max',
    'epdb',                         #NEW
    'epdb_flag',                    #NEW ???
    'extnrg',
    'fld',
    'fmap',
    'ga_cauchy_alpha',
    'ga_cauchy_beta',
    'ga_crossover_rate',
    'ga_crossover_mode_flag', 
    'ga_elitism',
    'ga_mutation_rate',
    'ga_num_evals',
    'ga_num_generations',
    'ga_pop_size',
    'ga_run',
    'ga_window_size',
    'include_1_4_interactions',     #NEW
    'include_1_4_interactions_flag',#NEW
    'intelec',
    'ligand_types',                 #NEW
    'linear_schedule',
    'ls_search_freq',
    'map',
    'move',
    'reorient_flag',
    'flexres_flag',
    'flexres',
    'outlev',
    'parameter_library',            #NEW
    'qstep',
    'quarf',
    'axisangle0',
    'rejs',
    'rmsatoms',
    'rmsatoms_flag',
    'rmsref',
    'rmstol',
    'rmsref',
    'rmsref_flag',
    'rt0',
    'rtrf',
    'runs',
    'seed',
    'select',
    'set_ga',
    'set_sw1',
    'set_psw1',
    'simanneal',
    'sw_lb_rho',
    'sw_max_fail',
    'sw_max_its',
    'sw_max_succ',
    'sw_rho',
    'torsdof4',
    'unbound',                      # set to user value
    'unbound_flag',                 # 0 by default
    'unbound_model',                # set to "extended"
    'unbound_model_flag',           # 1 by default
    'tran0',
    'trnrf',
    'tstep',
    'write_all',
    'write_all_flag']


#####AUTODOCK4_1####
##### and
#####AUTODOCK4_2####

 
###???parameter_file???
cluster_list4 = ['autodock_parameter_version','custom_parameter_file', 'ligand_types', 'rmstol', 'cluster', 'analysis' ]

genetic_algorithm_list4_1 = [
    'autodock_parameter_version',
    'outlev', 
    #'custom_parameter_file',       #REMOVE WHEN  autodock4.1 is released
    #'parameter_file',              #REMOVE WHEN autodock4.1 is released
    'include_1_4_interactions',
    'include_1_4_interactions_flag',
    'intelec',
    'seed',
    'unbound_intnbp_coeffs',
    'unbound_intnbp_coeffs_flag',
    'ligand_types',               
    'fld',
    'map',
    'move',
    'flexres_flag',
    'flexres',
    'about', 
    'reorient_flag',
    'tran0', 
    'axisangle0', 
    'dihe0', 
    'tstep', 
    'qstep', 
    'dstep', 
    'torsdof4', 
    'rmstol', 
    'rmsref',
    'rmsref_flag',
    'rmsatoms',
    'rmsatoms_flag',
    'extnrg', 
    'e0max', 
    'ga_pop_size', 
    'ga_num_evals', 
    'ga_num_generations', 
    'ga_elitism', 
    'ga_mutation_rate', 
    'ga_crossover_rate', 
    'ga_crossover_mode',
    'ga_crossover_mode_flag',
    'ga_window_size', 
    'ga_cauchy_alpha', 
    'ga_cauchy_beta', 
    'set_ga', 
    'unbound',                   
    'unbound_flag',                 #0 by default
    'unbound_energy',                   
    'unbound_energy_flag',          #0 by default
    'unbound_model',                #'bound' by default
    'unbound_model_flag',
    'epdb',                    
    'epdb_flag',              
    'do_global_only',               #'do_global_only' sets 'DPF_GS' here
    'write_all',
    'write_all_flag',
    'analysis' ]

genetic_algorithm_list4_2 = genetic_algorithm_list4_1


genetic_algorithm_local_search_list4_1 = [
    'autodock_parameter_version',
    'outlev', 
    #'custom_parameter_file',   #REMOVE WHEN  autodock4.1 is released
    #'parameter_file',          #REMOVE WHEN  autodock4.1 is released
    'include_1_4_interactions',
    'include_1_4_interactions_flag',
    'intelec',
    'seed',
    'unbound_intnbp_coeffs',
    'unbound_intnbp_coeffs_flag',
    'ligand_types',     
    'fld', 
    'map', 
    'move', 
    'flexres_flag',
    'flexres',
    'about', 
    'reorient_flag',
    'tran0', 
    'axisangle0', 
    'dihe0', 
    'tstep', 
    'qstep', 
    'dstep', 
    'torsdof4', 
    'rmstol', 
    'rmsref',
    'rmsref_flag',
    'rmsatoms',
    'rmsatoms_flag',
    'extnrg', 
    'e0max', 
    'ga_pop_size', 
    'ga_num_evals', 
    'ga_num_generations', 
    'ga_elitism', 
    'ga_mutation_rate', 
    'ga_crossover_rate', 
    'ga_crossover_mode',
    'ga_crossover_mode_flag',
    'ga_window_size', 
    'ga_cauchy_alpha', 
    'ga_cauchy_beta', 
    'set_ga', 
    'sw_max_its', 
    'sw_max_succ', 
    'sw_max_fail', 
    'sw_rho', 
    'sw_lb_rho', 
    'ls_search_freq', 
    'set_psw1', 
    'epdb',          
    'epdb_flag',    
    'unbound',                   
    'unbound_flag',                 #0 by default
    'unbound_energy',                   
    'unbound_energy_flag',          #0 by default
    'unbound_model',                #'bound' by default
    'unbound_model_flag',
    'ga_run', 
    'write_all',
    'write_all_flag',
    'analysis'
    ]

genetic_algorithm_local_search_list4_2 = genetic_algorithm_local_search_list4_1

genetic_algorithm_local_search_list4_1_with_parameter_file = [
    'autodock_parameter_version',
    'outlev', 
    'custom_parameter_file',  
    'parameter_file',        
    'include_1_4_interactions', 
    'include_1_4_interactions_flag', 
    'intelec',
    'seed', 
    'unbound_intnbp_coeffs',
    'unbound_intnbp_coeffs_flag',
    'ligand_types',                 
    'fld', 
    'map', 
    'move', 
    'flexres_flag',
    'flexres',
    'about', 
    'reorient_flag',
    'tran0', 
    'axisangle0', 
    'dihe0', 
    'tstep', 
    'qstep', 
    'dstep', 
    'torsdof4', 
    'rmstol', 
    'rmsref',
    'rmsref_flag',
    'rmsatoms',
    'rmsatoms_flag',
    'extnrg', 
    'e0max', 
    'ga_pop_size', 
    'ga_num_evals', 
    'ga_num_generations', 
    'ga_elitism', 
    'ga_mutation_rate', 
    'ga_crossover_rate', 
    'ga_crossover_mode',
    'ga_crossover_mode_flag',
    'ga_window_size', 
    'ga_cauchy_alpha', 
    'ga_cauchy_beta', 
    'set_ga', 
    'sw_max_its', 
    'sw_max_succ', 
    'sw_max_fail', 
    'sw_rho', 
    'sw_lb_rho', 
    'ls_search_freq', 
    'set_psw1', 
    'unbound',                   
    'unbound_flag',                 #0 by default
    'unbound_energy',                   
    'unbound_energy_flag',          #0 by default
    'unbound_model',                #'bound' by default
    'unbound_model_flag',
    'epdb',                      
    'epdb_flag',                
    'ga_run', 
    'write_all',
    'write_all_flag',
    'analysis'
    ]

genetic_algorithm_local_search_list4_2_with_parameter_file = genetic_algorithm_local_search_list4_1_with_parameter_file

local_search_list4_1 = [
    'autodock_parameter_version',
    'outlev', 
    #'custom_parameter_file',   #REMOVE WHEN  autodock4.1 is released
    #'parameter_file',          #REMOVE WHEN  autodock4.1 is released
    'include_1_4_interactions',
    'include_1_4_interactions_flag',
    'intelec',
    'seed', 
    'unbound_intnbp_coeffs',
    'unbound_intnbp_coeffs_flag',
    'ligand_types',                
    'fld', 
    'map', 
    'move', 
    'flexres_flag',
    'flexres',
    'about', 
    'reorient_flag',
    'tran0', 
    'axisangle0', 
    'dihe0', 
    'tstep', 
    'qstep', 
    'dstep', 
    'torsdof4', 
    'unbound',                    
    'unbound_flag',              
    'rmstol', 
    'rmsref',
    'rmsref_flag',
    'rmsatoms',
    'rmsatoms_flag',
    'extnrg', 
    'e0max', 
    'ga_pop_size', 
    'sw_max_its', 
    'sw_max_succ', 
    'sw_max_fail', 
    'sw_rho', 
    'sw_lb_rho', 
    'ls_search_freq', 
    'set_psw1', 
    'unbound',                   
    'unbound_flag',                 #0 by default
    'unbound_energy',                   
    'unbound_energy_flag',          #0 by default
    'unbound_model',                #'bound' by default
    'unbound_model_flag',
    'epdb',                     
    'epdb_flag',               
    'do_local_only', 
    'write_all',
    'write_all_flag',
    'analysis'
    ]

local_search_list4_2 = local_search_list4_1

simulated_annealing_list4_1 = [
    'autodock_parameter_version',
    'outlev', 
    #'custom_parameter_file',   #REMOVE WHEN  autodock4.1 is released
    #'parameter_file',          #REMOVE WHEN  autodock4.1 is released
    'include_1_4_interactions',
    'include_1_4_interactions_flag', 
    'intelec',
    'seed', 
    'unbound_intnbp_coeffs',
    'unbound_intnbp_coeffs_flag',
    'ligand_types',                 
    'fld', 
    'map', 
    'move', 
    'flexres_flag',
    'flexres',
    'about', 
    'reorient_flag',
    'tran0', 
    'axisangle0', 
    'dihe0', 
    'tstep', 
    'qstep', 
    'dstep', 
    'torsdof4', 
    'rmsatoms',
    'rmsatoms_flag',
    'rmstol', 
    'rmsref',
    'rmsref_flag',
    'extnrg', 
    'e0max', 
    'rt0', 
    'linear_schedule', 
    'rtrf', 
    'trnrf', 
    'quarf', 
    'dihrf', 
    'runs', 
    'cycles', 
    'accs', 
    'rejs', 
    'select', 
    'unbound',                   
    'unbound_flag',                 #0 by default
    'unbound_energy',                   
    'unbound_energy_flag',          #0 by default
    'unbound_model',                #'bound' by default
    'unbound_model_flag',
    'epdb',                         
    'epdb_flag',                    
    'simanneal', 
    'write_all',
    'write_all_flag',
    'analysis'
    ]
simulated_annealing_list4_2 = simulated_annealing_list4_1

docking_parameter_list4_1 = [
    'autodock_parameter_version',
    'about',
    'accs',
    'analysis',
    'cluster',
    'compute_unbound_extended',
    'compute_unbound_extended_flag',
    'cycles',
    'dihe0',
    'dihrf',
    'do_global_only',
    'do_local_only',
    'dstep',
    'e0max',
    'epdb',                         
    'epdb_flag',                    
    'extnrg',
    'fld',
    'fmap',
    'ga_cauchy_alpha',
    'ga_cauchy_beta',
    'ga_crossover_rate',
    'ga_crossover_mode',
    'ga_crossover_mode_flag',
    'ga_elitism',
    'ga_mutation_rate',
    'ga_num_evals',
    'ga_num_generations',
    'ga_pop_size',
    'ga_run',
    'ga_window_size',
    'include_1_4_interactions',     
    'include_1_4_interactions_flag',     
    'intelec',
    'ligand_types',                 
    'linear_schedule',
    'ls_search_freq',
    'map',
    'move',
    'reorient_flag',
    'flexres_flag',
    'flexres',
    'outlev',
    'parameter_library',            
    'qstep',
    'quarf',
    'axisangle0',
    'rejs',
    'rmsatoms',
    'rmsatoms_flag',
    'rmsref',
    'rmstol',
    'rmsref',
    'rmsref_flag',
    'rt0',
    'rtrf',
    'runs',
    'seed',
    'select',
    'set_ga',
    'set_sw1',
    'set_psw1',
    'simanneal',
    'sw_lb_rho',
    'sw_max_fail',
    'sw_max_its',
    'sw_max_succ',
    'sw_rho',
    'torsdof4',
    'tran0',
    'trnrf',
    'tstep',
    'unbound',                      
    'unbound_flag',                      
    'unbound_energy',                   
    'unbound_energy_flag',          #0 by default
    'unbound_model',
    'unbound_model_flag',
    'write_all',
    'write_all_flag']

docking_parameter_list4_2 = docking_parameter_list4_1


implemented = [
    'accs',
    'about',
    'analysis',
    'axisangle0',
    'cluster',
    'compute_unbound_extended',
    'compute_unbound_extended_flag',
    'custom_parameter_file',
    'cycles',
    'desolvmap',
    'dihe0',
    'do_global_only',
    'do_local_only',
    'dstep',
    'e0max',
    'elecmap',
    'epdb',
    'epdb_flag',
    'extnrg',
    'fld',
    'flexres',
    'flexres_flag',
    'flexible_residues',
    'ga_cauchy_alpha',
    'ga_cauchy_beta',
    'ga_crossover_rate',
    'ga_crossover_mode',
    'ga_crossover_mode_flag',
    'ga_elitism',
    'ga_mutation_rate',
    'ga_num_evals',
    'ga_num_generations',
    'ga_pop_size',
    'ga_run',
    'ga_window_size',
    'include_1_4_interactions',
    'include_1_4_interactions_flag',
    'intelec',
    'intelec4',
    'intnbp_r_eps',
    'ligand_types',
    'linear_schedule',
    'ls_search_freq',
    'map',
    'fmap',
    'move',
    'ndihe',
    'outlev',
    'output_pop_file',
    'parameter_file',
    'psw_trans_scale',
    'psw_rot_scale',
    'psw_tors_scale',
    'qstep',
    'quat0',
    'rejs', 
    'reorient',
    'reorient_flag',
    'rmsatoms',
    'rmsatoms_flag',
    'rmsref',
    'rmsref_flag',
    'rmstol',
    'rt0',
    'rtrf',
    'runs',
    'trnrf', 
    'quarf',
    'dihrf',
    'seed',
    'select',
    'set_ga',
    'set_sw1',
    'set_psw1',
    'simanneal',
    'sw_max_its', 
    'sw_max_succ', 
    'sw_max_fail', 
    'sw_rho', 
    'sw_lb_rho', 
    'torsdof',
    'torsdof4',
    'tran0',
    'tstep',
    'types',
    'unbound',
    'unbound_flag',
    'unbound_energy',                   
    'unbound_energy_flag',          #0 by default
    'unbound_intnbp_coeffs',
    'unbound_intnbp_coeffs_flag',
    'unbound_model',
    'unbound_model_flag',
    'autodock_parameter_version',
    'write_all',
    'write_all_flag']
                

class DockingParameterFileMaker:
    """Accept a <ligand>.pdbq and <receptor>.pdbqs and create
    <ligand>_<receptor>.dpf
    """

    def __init__(self, verbose = None):
        self.verbose = verbose
        self.dpo = DockingParameters()

    def set_ligand(self, ligand_filename): 
        verbose = self.verbose
        self.ligand_filename = os.path.basename(ligand_filename)
        if verbose: print "set ligand_filename to", self.ligand_filename
        self.dpo.set_ligand(ligand_filename)
        #expect a filename like ind.out.pdbq: get 'ind' from it
        self.ligand_stem = string.split(self.ligand_filename,'.')[0]
        if verbose: print "set ligand_stem to", self.ligand_stem
        self.ligand = Read(ligand_filename)[0]
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
        self.dpo['types']['value'] = self.ligand.types


    def set_docking_parameters(self, **kw):
        """Any docking parameters should be set here
        """
        # like this: 
        # newdict = {'ga_num_evals':1750000, 'ga_pop_size':150,
        #            'ga_run':20, 'rmstol':2.0}
        # self.mv.dpo['<parameter>']['value'] = <new value>
        for parm, newvalue in kw.items():
            if self.verbose: print 'set dpo for ', parm, ' to ', newvalue, ' check=', self.dpo[parm]['value']
            self.dpo[parm]['value'] = newvalue
        #self.dpo['ga_num_evals']['value'] = 1750000
        #self.dpo['ga_run']['value'] = 20
        #self.dpo['ga_pop_size']['value'] = 150
        #self.dpo['rmstol']['value'] = 2.0


    def write_dpf(self, dpf_filename,
              parm_list = genetic_algorithm_local_search_list):
        if not dpf_filename:
            dpf_filename = "%s%s%s%s" % \
                           (self.ligand_stem, "_",
                            self.receptor_stem, ".dpf")
        # now that we have a filename...
        if self.verbose: print "writing ", dpf_filename
        self.dpo.write(dpf_filename, parm_list, autodock_parameter_version=self.autodock_parameter_version)


                
class DockingParameter4FileMaker:
    """Accept a <ligand>.pdbqt and <receptor>.pdbqt and create
    <ligand>_<receptor>4.dpf
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
        if value=='True':
            value=1
        if value==True:
            value=1
        if value=='False':
            value=0
        if value==False:
            value=0
        verbose = self.verbose
        self.dpo['write_all']['value'] = value
        if verbose: print "set write_all to", self.dpo['write_all']['value']


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
            print '!unable to use ligand.TORSDOF! set torsdof to ligand.ndihe=', self.ligand.ndihe
            self.dpo['torsdof4']['value'][0] = self.ligand.ndihe
        if verbose: print "set torsdof4 to ", self.dpo['torsdof4']['value']
        #types
        self.ligand.types = self.getTypes(self.ligand)
        self.dpo['ligand_types']['value'] = self.ligand.types
        if verbose: print "set types to ", self.dpo['ligand_types']['value']
        #about
        cen = self.ligand.getCenter()
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
              parm_list = genetic_algorithm_local_search_list4,
              pop_seed=False):
        if not dpf_filename:
            dpf_filename = "%s%s%s%s" % \
                           (self.ligand_stem, "_",
                            self.receptor_stem, ".dpf")
        # set initial conformation
        if pop_seed:
            self.dpo['tran0']['value'] = self.dpo['about']['value']
            self.dpo['axisangle0']['value'] = '1.0 0. 0. 0.'
            dihe0 = '0. '*self.dpo['ndihe']['value']
            dihe0.rstrip()
            self.dpo['dihe0']['value'] = dihe0 
        # now that we have a filename...
        if self.verbose:
            print "writing ", dpf_filename
        self.dpo.write4(dpf_filename, parm_list)


                
class DockingParameter42FileMaker:
    """Accept a <ligand>.pdbqt and <receptor>.pdbqt and create
    <ligand>_<receptor>4.dpf
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
        if value=='True':
            value=1
        if value==True:
            value=1
        if value=='False':
            value=0
        if value==False:
            value=0
        verbose = self.verbose
        self.dpo['write_all']['value'] = value
        if verbose: print "set write_all to", self.dpo['write_all']['value']


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
            print '!unable to use ligand.TORSDOF! set torsdof to ligand.ndihe=', self.ligand.ndihe
            self.dpo['torsdof4']['value'][0] = self.ligand.ndihe
        if verbose: print "set torsdof4 to ", self.dpo['torsdof4']['value']
        #types
        self.ligand.types = self.getTypes(self.ligand)
        self.dpo['ligand_types']['value'] = self.ligand.types
        if verbose: print "set types to ", self.dpo['ligand_types']['value']
        #about
        cen = self.ligand.getCenter()
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
              parm_list = genetic_algorithm_local_search_list4,
              pop_seed=False):
        if not dpf_filename:
            dpf_filename = "%s%s%s%s" % \
                           (self.ligand_stem, "_",
                            self.receptor_stem, ".dpf")
        # set initial conformation
        if pop_seed:
            self.dpo['tran0']['value'] = self.dpo['about']['value']
            self.dpo['axisangle0']['value'] = '1.0 0. 0. 0.'
            dihe0 = '0. '*self.dpo['ndihe']['value']
            dihe0.rstrip()
            self.dpo['dihe0']['value'] = dihe0 
        # now that we have a filename...
        if self.verbose:
            print "writing ", dpf_filename
        self.dpo.write42(dpf_filename, parm_list)

        
if __name__=='__main__':
    import os

    if len(sys.argv) > 1:
        dpo = DockingParameters()
        try:
            dpo.read(sys.argv[1])
        except IOError, msg:
            print "IOError: %s" % (msg)
            exit(-1)

        tmp_filename = '/tmp/%d.dpf' % (os.getpid())
        dpo.write( tmp_filename, dpo.file_params)

        # instantiate a new DP object and read from file just written
        new_dpo = DockingParameters()
        new_dpo.read(tmp_filename)

        if dpo.file_params!=new_dpo.file_params:
            for p in dpo.file_params:
                if not p in new_dpo.file_params:
                    print "%s missing from written dpf!" % p['keyword']
            for p in new_dpo.file_params:
                if not p in dpo.file_params:
                    print "Unexpected %s in written dpf!" % p['keyword']
        for p in new_dpo.file_params:
            if not dpo[p]['value']==new_dpo[p]['value']:
                print "Parameter Discrepancy: %s" % (p)
                print "    I wrote this: %s" % (dpo[p]['value'])
                print "But, I read this: %s" % (new_dpo[p]['value'])
    else:
        dpo = DockingParameters('receptor.test.pdbqs', 'ligand.pdbq')

        # compare docking_parameter_list with keys of DockingParameter
        # dictionary and report discrepancies

        if docking_parameter_list!=dpo.keys():
            print "Comparing docking paramter dictionary and list..."
            for p in dpo.keys():
                if not p in docking_parameter_list:
                    print ">> %s in dictionary but not docking_parameter_list" % (p)
            for p in docking_parameter_list:
                if not p in dpo.keys():
                    print ">> %s in docking_parameter_list but not dictionary" % (p)
            print "done with comparison."

        # This should also run through all the parameters, write default
        # values to a file and then read the file compare with the defaults
        print "Here are the parameters implemented so far (write):"
        dpo.write('', implemented)
        print
        print "Here are the parameters not implemented yet (write):"
        for p in docking_parameter_list:
            if not p in implemented:
                print p

        tmp_filename = '/tmp/%d.dpf' % (os.getpid())
        dpo.write( tmp_filename, implemented)

        new_dpo = DockingParameters()
        new_dpo.read(tmp_filename)

        # write the new dpo to /dev/null to update the value fields
        for p in implemented:
            if not dpo[p]['value']==new_dpo[p]['value']:
                print "%s: %s(%s) -> %s(%s)" % (dpo[p]['keyword'],
                                                dpo[p]['value'],
                                                type(dpo[p]['value']).__name__,
                                                new_dpo[p]['value'],
                                                type(new_dpo[p]['value']).__name__)


