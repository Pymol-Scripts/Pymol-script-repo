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

from AutoDockTools.VisionInterface.Adt.receptor import receptor

import os

class StructureBrowser(NetworkNode):
    """
    Allows the user to browser a structure file in pdb, pqr, or pdbqt format 
    and returns a receptor object.
    """
    
    def __init__(self, name='StructureBrowser', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
        ip = self.inputPortsDescr
        ip.append(datatype='string', name='receptor_file')

        filetypes = [('All supported files', '*.pdb *.pqr *.pdbqt')]

        self.widgetDescr['receptor_file'] = {
            'class':'NEEntryWithFileBrowser', 'master':'node', 'width':20,
            'filetypes':filetypes,
            'initialValue':'', 'labelCfg':{'text':'Structure file: '}
            }

        op = self.outputPortsDescr
        op.append(datatype='receptor', name='receptor_obj')

        code = """def doit(self, receptor_file):
    import os

    receptor_file = os.path.abspath(receptor_file)

    if not(os.path.exists(receptor_file)):
        print "ERROR: structure file " + receptor_file + " does not exist!"
        return 'stop'

    receptor_obj = receptor(receptor_file)

    self.outputData(receptor_obj=receptor_obj)
"""
        self.setFunction(code)
