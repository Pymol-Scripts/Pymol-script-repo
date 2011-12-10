#########################################################################
#
# Date: Nov 2001 Authors: Michel Sanner
#
#    sanner@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Michel Sanner and TSRI
#
#########################################################################

from NetworkEditor.items import NetworkNode
from AutoDockTools.VisionInterface.Adt.dpf_template import dpf_template
from AutoDockTools.VisionInterface.Adt.gpf_template import gpf_template

class MakeGpfDpfCopies(NetworkNode):
    """
    A node that takes a DPF, a GPF and a directory of pdb/pqr/pdbqt files and 
    makes copies of the DPF/GPF.  The DPF/GPF copies will have the basenames of the
    pdb/pqr/pdbqt files.  
  
    Input 1: GPF template file with the extension .gpf
    Input 2: DPF template file with the extension .dpf
    Input 3: A directory path to the PDB/PQR/PDBQT files
    Output 1:  List of all GPFs in the directory
    Output 2:  List of all DPFs in the directory
    Output 3:  Same as input 3
    """
    
    def __init__(self, name='MakeDPFCopies', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
        ip = self.inputPortsDescr
        ip.append(datatype='gpf_template', name='gpf_file', required=False)
        ip.append(datatype='dpf_template', name='dpf_file', required=False)
        ip.append(datatype='string', name='struct_dir', required=True)

        op = self.outputPortsDescr
        op.append(datatype='list', name='gpf_list')
        op.append(datatype='list', name='dpf_list')
        op.append(datatype='string', name='struct_dir')

        code = """def doit(self, gpf_file, dpf_file, struct_dir):
    import os
    import shutil
    from AutoDockTools.VisionInterface.Adt.dpf_template import dpf_template
    from AutoDockTools.VisionInterface.Adt.gpf_template import gpf_template

    if dpf_file == None and gpf_file == None:
        print "ERROR: DPF and GPF input missing"
        return 'stop'


    if dpf_file != None:
        dpf_file_path = dpf_file.fullpath

        if not(os.path.exists(dpf_file_path)):
            print "ERROR: DPF template file " + dpf_file_path + " does not exist!"
            return 'stop'

    if gpf_file != None:
        gpf_file_path = gpf_file.fullpath

        if not(os.path.exists(gpf_file_path)):
            print "ERROR: GPF template file " + gpf_file_path + " does not exist!"
            return 'stop'

    name_list = set()
    d = os.path.abspath(struct_dir)

    for i in os.listdir(struct_dir):
        if i.endswith(".pdbqt") or i.endswith(".pdb") or i.endswith(".pqr"):
            name_list.add(i.split('.')[0])

    dpf_list = []
    gpf_list = []

    if dpf_file != None:
        for i in name_list:
            d_dst = os.path.join(d, i + '.dpf')
            dpf_list.append(d_dst)       

            if not(os.path.exists(d_dst)):
                shutil.copyfile(dpf_file_path, d_dst)

    if gpf_file != None:
        for i in name_list:
            g_dst = os.path.join(d, i + '.gpf')
            gpf_list.append(g_dst)       

            if not(os.path.exists(g_dst)):
                shutil.copyfile(gpf_file_path, g_dst)
  

    self.outputData(dpf_list=dpf_list, gpf_list=gpf_list, struct_dir=struct_dir)
"""
        self.setFunction(code)
