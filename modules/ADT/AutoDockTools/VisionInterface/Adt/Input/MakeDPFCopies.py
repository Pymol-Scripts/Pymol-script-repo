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

class MakeDPFCopies(NetworkNode):
    """
    A node that takes a DPF and a directory of pdb/pqr/pdbqt files and 
    makes copies of the DPF.  The DPF copies will have the basenames of the
    pdb/pqr/pdbqt files.  
  
    Input 1: DPF template file with the extension .dpf
    Input 2: A directory path to the PDB/PQR/PDBQT files
    Output:  List of all DPFs in the directory
    """
    
    def __init__(self, name='MakeDPFCopies', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
        ip = self.inputPortsDescr
        ip.append(datatype='dpf_template', name='dpf_file')
        ip.append(datatype='string', name='struct_dir')

        op = self.outputPortsDescr
        op.append(datatype='list', name='dpf_list')

        code = """def doit(self, dpf_file, struct_dir):
    import os
    import shutil
    from AutoDockTools.VisionInterface.Adt.dpf_template import dpf_template

    file_path = dpf_file.fullpath

    if not(os.path.exists(file_path)):
        print "ERROR: DPF template file " + file_path + " does not exist!"
        return 'stop'

    name_list = set()
    dpf_list = []

    d = os.path.abspath(struct_dir)

    for i in os.listdir(struct_dir):
        if i.endswith(".pdbqt") or i.endswith(".pdb") or i.endswith(".pqr"):
            name_list.add(i.split('.')[0])

    for i in name_list:
        dst = os.path.join(d, i + '.dpf')
        dpf_list.append(dst)       

        if not(os.path.exists(dst)):
            shutil.copyfile(file_path, dst)
  

    self.outputData(dpf_list=dpf_list)
"""
        self.setFunction(code)
