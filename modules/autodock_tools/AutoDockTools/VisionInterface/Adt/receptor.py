#############################################################################
#
# Author: J. Ren
#
#############################################################################

"""
Receptor - can be a pdb, pqr, or pdbqt
"""

import os

class receptor:

    def __init__(self, filename=None):
        self.ext_loc = {'pdb': [None, None], 'pqr': [None, None], 'pdbqt': [None, None]}
        self.workdir = None
        self.id = None

        if filename != None:
            self.set_ext_loc(filename)
            self.set_workdir(filename)
 
    def set_workdir(self, d):
        if not(d.startswith('http://')):
            d = os.path.abspath(d)

            if d.endswith('.pdb') or d.endswith('.pqr') or d.endswith('.pdbqt'):
                d = os.path.dirname(d)

            self.workdir = d

    def set_id(self, id):
        bn = os.path.basename(id)
        self.id = bn.split(os.extsep)[0]

    def set_ext_loc(self, filename):
        ext = os.path.splitext(filename)[1].strip(os.extsep)
        self.set_id(filename)        

        if filename.startswith('http://'):
            self.ext_loc[ext] = [filename, 'url']
        else:
            full_path = os.path.abspath(filename)
            dir_name = os.path.dirname(full_path)
            
            ext_list = ['pdb', 'pqr', 'pdbqt']

            for i in ext_list:
                f = dir_name + os.sep + self.id + os.extsep + i
                if os.path.exists(f):
                    self.ext_loc[i] = [f, 'local']
                    self.set_workdir(full_path)

    def get_workdir(self):
        if self.workdir == None:
            exts = ['pdb', 'pqr', 'pdbqt']
            
            for i in exts:
                if self.get_ext_type(i) == "local":
                    self.set_workdir(self.get_ext_loc(i))
                    break

        return self.workdir

    def get_id(self):
        return self.id

    def get_ext_loc(self, ext):
        return self.ext_loc[ext][0]

    def get_ext_type(self, ext):
        return self.ext_loc[ext][1] 

