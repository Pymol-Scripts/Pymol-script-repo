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

from AutoDockTools.VisionInterface.Adt.autogrid_results import autogrid_results

import os

class AutogridResDirectoryBrowser(NetworkNode):
    """
    Browse local directory that contains Autogrid results.
    The directory will be compressed into a zip file.

    Input: Directory that contains autogrid results
    
    Output:
    port 1: autogrid_result object that contains info about the autogrid result directory
    port 2: path to compressed autogrid directory.  
    """
    
    def __init__(self, name='AutogridResDirectoryBrowser', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )


        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
        ip = self.inputPortsDescr
        ip.append(datatype='string', name='directory')

        self.widgetDescr['directory'] = {
            'class':'NEEntryWithDirectoryBrowser', 'master':'node', 'width':16,
            'initialValue':'', 'labelCfg':{'text':'Autogrid result directory: '}
            }

        op = self.outputPortsDescr
        op.append(datatype='autogrid_results', name='autogrid_result')
        op = self.outputPortsDescr
        op.append(datatype='string', name='zip_file_path')

        code = """def doit(self, directory):
    autogrid_result = autogrid_results(directory, "local")
          
    self.outputData(autogrid_result=autogrid_result, zip_file_path=autogrid_result.path)
"""
        self.setFunction(code)
