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

class DPFTemplateBrowser(NetworkNode):
    """
    A node that allows the user to browse for a DPF template with the extension .dpf
  
    Input: DPF template file with the extension .dpf
    Output: dpf_template object that contains info about the DPF template file
    """
    
    def __init__(self, name='DPFTemplateBrowser', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
        ip = self.inputPortsDescr
        ip.append(datatype='string', name='dpf_template_file')

        filetypes = [('All supported files', '*.dpf')]

        self.widgetDescr['dpf_template_file'] = {
            'class':'NEEntryWithFileBrowser', 'master':'node', 'width':20,
            'filetypes':filetypes,
            'initialValue':'', 'labelCfg':{'text':'DPF template file: '}
            }

        op = self.outputPortsDescr
        op.append(datatype='dpf_template', name='dpf_template')

        code = """def doit(self, dpf_template_file):
    import os
    from AutoDockTools.VisionInterface.Adt.dpf_template import dpf_template

    dpf_template_file = os.path.abspath(dpf_template_file)

    if not(os.path.exists(dpf_template_file)):
        print "ERROR: DPF template file " + dpf_template_file + " does not exist!"
        return 'stop'

    dpf_template = dpf_template(dpf_template_file)

    self.outputData(dpf_template=dpf_template)
"""
        self.setFunction(code)
