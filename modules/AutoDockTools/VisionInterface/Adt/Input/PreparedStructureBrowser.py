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

class PreparedStructureBrowser(NetworkNode):
    """
    Allows the user to browse a structure prepared pdbqt structure.
    """
    
    def __init__(self, name='StructureBrowser', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
        ip = self.inputPortsDescr
        ip.append(datatype='string', name='prepared_structure_file')

        filetypes = [('All supported files', '*.pdbqt')]

        self.widgetDescr['prepared_structure_file'] = {
            'class':'NEEntryWithFileBrowser', 'master':'node', 'width':20,
            'filetypes':filetypes,
            'initialValue':'', 'labelCfg':{'text':'Prepared Structure file: '}
            }

        op = self.outputPortsDescr
        op.append(datatype='receptor_prepared', name='receptor_obj')

        code = """def doit(self, prepared_structure_file):
    import os

    receptor_file = os.path.abspath(prepared_structure_file)

    if not(os.path.exists(prepared_structure_file)):
        print "ERROR: prepared structure file " + prepared_structure_file + " does not exist!"
        return 'stop'

    receptor_obj = receptor(prepared_structure_file)

    self.outputData(receptor_obj=receptor_obj)
"""
        self.setFunction(code)
