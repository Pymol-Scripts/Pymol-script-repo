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
from AutoDockTools.VisionInterface.Adt.gpf_template import gpf_template

class MakeGPFCopies(NetworkNode):
    """
    A node that takes a GPF and a directory of pdb/pqr/pdbqt files and 
    makes copies of the GPF.  The GPF copies will have the basenames of the
    pdb/pqr/pdbqt files.  
  
    Input 1: GPF template file with the extension .gpf
    Input 2: A directory path to the PDB/PQR/PDBQT files
    Output:  List of all GPFs in the directory
    """
    
    def __init__(self, name='MakeGPFCopies', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
        ip = self.inputPortsDescr
        ip.append(datatype='gpf_template', name='gpf_file')
        ip.append(datatype='string', name='struct_dir')

        op = self.outputPortsDescr
        op.append(datatype='list', name='gpf_list')

        code = """def doit(self, gpf_file, struct_dir):
    import os
    import shutil
    from AutoDockTools.VisionInterface.Adt.gpf_template import gpf_template

    file_path = gpf_file.fullpath

    if not(os.path.exists(file_path)):
        print "ERROR: GPF template file " + file_path + " does not exist!"
        return 'stop'

    name_list = set()
    gpf_list = []

    d = os.path.abspath(struct_dir)

    for i in os.listdir(struct_dir):
        if i.endswith(".pdbqt") or i.endswith(".pdb") or i.endswith(".pqr"):
            name_list.add(i.split('.')[0])

    for i in name_list:
        dst = os.path.join(d, i + '.gpf')
        gpf_list.append(dst)       

        if not(os.path.exists(dst)):
            shutil.copyfile(file_path, dst)
  

    self.outputData(gpf_list=gpf_list)
"""
        self.setFunction(code)
