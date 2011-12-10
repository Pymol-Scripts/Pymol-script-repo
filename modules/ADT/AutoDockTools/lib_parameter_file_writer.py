#
# 
#
# $Id: lib_parameter_file_writer.py,v 1.4 2004/08/25 20:44:33 rhuey Exp $
#

import os, string, glob
from MolKit import Read
from AutoDockTools.GridParameters import GridParameters, grid_parameter_list
from AutoDockTools.DockingParameters import DockingParameters
from AutoDockTools.DockingParameters import genetic_algorithm_list, \
                genetic_algorithm_local_search_list, local_search_list,\
                simulated_annealing_list



class LibraryGpfWriter:

    def __init__(self, gpffile, extension='pdbq', err_file="write_lib_gpf_errors", 
                           prefix='',alltypes=None, extratypes_dict={}, verbose=0):
        self.verbose = verbose
        if verbose: print "self.verbose = ", verbose
        self.extratypes_dict = extratypes_dict
        self.setup(gpffile)
        self.alltypes = alltypes
        if self.alltypes is None:
            self.alltypes = self.gettypes(extension)
        self.prefix = prefix
        #real element, Rii, epsii
        #eg {X:['Cl', 4.09, .02345]}
        #WHAT errors should be captured
        #self.error_file = open(err_file, 'w')
    

    def setup(self,gpffile):
        self.gpffile = gpffile
        self.gpo = GridParameters() 
        self.gpo.read(gpffile)
        if len(self.extratypes_dict):
            print "update each type in extratypes here"
        #????RENAME extratype in receptor and write new output file here?


    def set_receptor(self, receptor_filename):
        self.gpo.set_receptor(receptor_filename)
        ty


    def gettypes(self, extension, name=''):
        typeDict = {}
        key = name + '*.' + extension
        ligandfilenames = glob.glob(key)
        if self.verbose:
            print "processing ", len(ligandfilenames), " files matching ", key
        for ligand_file in ligandfilenames:
            m = Read(ligand_file)[0]
            for a in m.allAtoms:
                typeDict[a.autodock_element] = 1
            del(m)
        alltypes = typeDict.keys()
        return alltypes


    def write(self, gpo):
        #set types in groups of 6 + write
        outputlist = []
        number_of_files = len(self.alltypes)/6 + 1
        ctr = 0
        last = 6
        for i in range(number_of_files):
            if i == number_of_files-1:
                current_types = self.alltypes[ctr*6:]
            elif last>len(self.alltypes):
                current_types = self.alltypes[ctr*6:]
            else:
                current_types = self.alltypes[ctr*6:last]
            gpo['types']['value'] =  string.join(current_types,'')
            if not gpo.receptor_stem:
                gpo.receptor_stem = os.path.splitext(gpo['receptor']['value'])[0]
            outputfilename = self.prefix + gpo.receptor_stem + '_'+ str(ctr) + '_library.gpf'
            gpo.write(outputfilename, grid_parameter_list)
            ctr = ctr + 1
            last = last + 6
            outputlist.append(outputfilename)
        return outputlist

#what about checking that the box is big enough for each ligand??



class LibraryDpfWriter:


    def __init__(self, dpffile, extension='pdbq', \
                    err_file="write_lib_dpf_errors", 
                    prefix='', verbose=0):
        self.verbose = verbose
        if verbose: print "self.verbose = ", verbose
        self.prefix = prefix
        self.setup(dpffile)
        self.parm_list_dict = {}
        self.parm_list_dict['GALS'] = genetic_algorithm_local_search_list
        self.parm_list_dict['GA'] = genetic_algorithm_list
        self.parm_list_dict['LS'] = local_search_list
        self.parm_list_dict['SA'] = simulated_annealing_list
        self.ligandfilenames = []
        self.extension = extension
        self.error_file = open(err_file, 'w')


    def setup(self, dpffile):
        self.dpffile = dpffile
        self.dpo = DockingParameters() 
        self.dpo.read(dpffile)


    def getligands(self, extension, name=''):
        typeDict = {}
        key = name + '*.' + extension
        self.ligandfilenames = glob.glob(key)
        return self.ligandfilenames


    def write(self, dpo, style='GALS'):
        #one file per ligand
        dpffiles = []
        receptor_stem = self.dpo.receptor_stem
        parm_list = self.parm_list_dict[style]
        if not len(self.ligandfilenames):
            self.getligands(self.extension)
        #print "len(self.ligandfilenames)=", len(self.ligandfilenames)
        for ligand_file in self.ligandfilenames:
            #get the types 
            typeDict = {}
            m = Read(ligand_file)[0]
            for a in m.allAtoms:
                typeDict[a.autodock_element] = 1
            type_list = typeDict.keys()
            #check that there is a map for each type
            ok = 1
            for t in type_list:
                mapfile = receptor_stem + '.' + t + '.map'
                try:
                    assert os.path.exists(mapfile)
                except:
                    ok = 0
                    ostr = ligand_file + " missing map " + t + "\n"
                    self.error_file.write(ostr)
                    break  #try to continue 'for' statement with next filename
            #update dpo with ligand specific values
            if not ok:
                if self.verbose:
                    print "problem with ", ligand_file
                continue
            types = string.join(type_list,'')
            self.dpo['move']['value'] = ligand_file
            if self.verbose: print "set types to ", types
            self.dpo['types']['value'] = types
            #CAUTION: dpo['torsdof']['value'] is [0,0.3113]
            self.dpo['torsdof']['value'][0] = m.TORSDOF
            if hasattr(m, 'ndihe'):
                ndihe = m.ndihe
            elif hasattr(m, 'torscount'):
                ndihe = m.torscount
            elif hasattr(m, 'torTree'):
                ndihe = len(m.torTree.torsionMap)
            else:
                msg = ligand_file + " is not properly formatted: no torsions!"
                raise AttributeError, msg
            self.dpo['ndihe']['value'] = ndihe
            outputfilename = self.prefix + m.name + "_" + style + ".dpf"
            if self.verbose: print "writing ", outputfilename
            self.dpo.write(outputfilename, parm_list)
            dpffiles.append(outputfilename)
        return dpffiles

            

