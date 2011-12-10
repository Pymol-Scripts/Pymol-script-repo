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

class GPFTemplateBrowser(NetworkNode):
    """
    A node that allows the user to browse for a GPF template with the extension .gpf
  
    Input: GPF template file with the extension .gpf
    Output: gpf_template object that contains info about the GPF template file
    """
    
    def __init__(self, name='GPFTemplateBrowser', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
        ip = self.inputPortsDescr
        ip.append(datatype='string', name='gpf_template_file')

        filetypes = [('All supported files', '*.gpf')]

        self.widgetDescr['gpf_template_file'] = {
            'class':'NEEntryWithFileBrowser', 'master':'node', 'width':20,
            'filetypes':filetypes,
            'initialValue':'', 'labelCfg':{'text':'GPF template file: '}
            }

        op = self.outputPortsDescr
        op.append(datatype='gpf_template', name='gpf_template')

        code = """def doit(self, gpf_template_file):
    import os
    from AutoDockTools.VisionInterface.Adt.gpf_template import gpf_template

    gpf_template_file = os.path.abspath(gpf_template_file)

    if not(os.path.exists(gpf_template_file)):
        print "ERROR: GPF template file " + gpf_template_file + " does not exist!" 
        return 'stop'

    gpf_template = gpf_template(gpf_template_file)

    self.outputData(gpf_template=gpf_template)
"""
        self.setFunction(code)
