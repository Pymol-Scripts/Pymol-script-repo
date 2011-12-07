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

from AutoDockTools.VisionInterface.Adt.LigandDB import LigandDB

import sys
import os

class UrlLigandDB(NetworkNode):
    """
    Input: URL to a ligand DB on the virtual screening server.
    Output: LigandDB object containing info about the input
    """
    
    def __init__(self, name='UrlLigandDB', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )


        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
        ip = self.inputPortsDescr
        ip.append(datatype='string', name='url')
#        ip.append(datatype='string', name='local_dir')

        self.widgetDescr['url'] = {
            'class':'NEEntry', 'master':'node', 'width':40,
            'labelCfg':{'text':'URL to ligands: '}
            }
#        self.widgetDescr['local_dir'] = {
#            'class':'NEEntryWithDirectoryBrowser', 'master':'node', 'width':30,
#            'labelCfg':{'text':'Local directory to save ligands: '}
#            }

        op = self.outputPortsDescr
        op.append(datatype='LigandDB', name='ligDB')

#        code = """def doit(self, url, local_dir):
        code = """def doit(self, url):
#    ligDB = LigandDB(url_compressed_file=url)
    ligDB = LigandDB(url_lib=url)

    bname = os.path.basename(url)

#    if bname.endwith('.tar.gz'):
#        print "TEST tar.gz"
#    elif bname.endswith('.zip'):
#        print "TEST ZIP"

    if not(url.startswith('http://kryptonite.nbcr.net/app')):
         print "ERROR: the URL must be a directory on kryptonite.nbcr.net"
         sys.exit()  
 
    

    self.outputData(ligDB=ligDB)
"""
        self.setFunction(code)
