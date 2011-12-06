from NetworkEditor.items import NetworkNode

from AutoDockTools.VisionInterface.Adt.autogrid_results import autogrid_results

class AutogridResURL(NetworkNode):
    """
    Creates an autogrid_results object form a URL to the autogrid results
   
    Input: URL to autogrid results
    Output: autogrid_results object that contains info about the URL to autogrid results
    """
    
    def __init__(self, name='AutogridResURL', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )


        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
        ip = self.inputPortsDescr
        ip.append(datatype='string', name='url')

        self.widgetDescr['url'] = {
            'class':'NEEntry', 'master':'node', 'width':40,
            'labelCfg':{'text':'URL to Autogrid Results: '}
            }

        op = self.outputPortsDescr
        op.append(datatype='autogrid_results', name='autogrid_res_obj')

        code = """def doit(self, url):
    autogrid_res_obj=autogrid_results(url, "url")

    self.outputData(autogrid_res_obj=autogrid_res_obj)
"""
        self.setFunction(code)
