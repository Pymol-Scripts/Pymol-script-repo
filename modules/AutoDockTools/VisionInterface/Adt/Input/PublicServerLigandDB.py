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
from mglutil.util.packageFilePath import getResourceFolderWithVersion

import os
import time
import urllib2

class PublicServerLigandDB(NetworkNode):
    """
    List of available public libraries on the virtual screening server.
    A description of the ligand libraries can be found on
    http://nbcr.sdsc.edu/pub/wiki/index.php?title=Virtual_Screening_Libraries

    Input: a public ligand library name
    Output: LigandDB object containing info about the info
    """
    
    def __init__(self, name='PublicServerLigandDB', **kw):
        import urllib

        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
        ip = self.inputPortsDescr
        ip.append(datatype='string', name='server_lib', required=True, )

        fqdn = "kryptonite.nbcr.net"
        url = "http://" + fqdn + "/pub_ligand_libs.txt"

        publibdir = os.path.join(getResourceFolderWithVersion(), 'ws')

        if not (os.path.exists(publibdir)):
            os.mkdir(publibdir)
        
        publiblocal = os.path.join(publibdir, 'publibs.txt')

        lock = publiblocal + '.lock'

        if os.path.exists(lock) and time.time() - os.path.getmtime(lock) > 15:
            os.remove(lock)
            
        try:
            if not(os.path.exists(lock)):
                open(lock, 'w').close()
                publibweb = urllib2.urlopen(url)
                outfile = open(publiblocal, 'w')
                outfile.write(publibweb.read())
                outfile.close()
                os.remove(lock)
        except:
            print "[INFO]: Getting list of public server libs from cache"
            pass
            

        try:
            f = open(publiblocal, 'r')
            self.choices = f.read().split()
            f.close()
        except:
            self.choices = []
            print "[ERROR]: Unable to public server libs from the web and from cache"

        self.widgetDescr['server_lib'] = {
            'class':'NEComboBox', 'master':'node',
            'choices':self.choices,
            'fixedChoices':True,
            'entryfield_entry_width':18,
            'labelGridCfg':{'sticky':'w'},
            'widgetGridCfg':{'sticky':'w'},
            'labelCfg':{'text':'Server Libraries:'}}

        op = self.outputPortsDescr
        op.append(datatype='LigandDB', name='ligDB')

        code = """def doit(self, server_lib):
    ligDB = LigandDB(server_lib=server_lib)

    self.outputData(ligDB=ligDB)
"""
        self.setFunction(code)
