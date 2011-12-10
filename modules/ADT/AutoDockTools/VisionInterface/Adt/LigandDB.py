#############################################################################
#
# Author: J. Ren
#
#############################################################################

"""
Ligand Library Object
"""

import os
import tempfile
import zipfile
import string, types, re
from AutoDockTools.filterLigands import CalculateProperties, PropertyTable, LigandList, FilterLigands

class LigandDB:
    """Ligand Library Object
    """
    
#    def __init__(self, server_lib=None, url_lib=None, filter_file=None, accepted=None):
    def __init__(self, server_lib=None, url_lib=None):

        self.server_lib = server_lib
        self.url_lib = url_lib
#        self.accepted = accepted
        self.accepted = None
        self.property_file = None
        self.propertyTable = None

#        if filter_file == None:
#            self.filter_file = os.path.abspath('filtered_ligands.txt')
#        else:
#            self.filter_file = os.path.abspath(filter_file)
        self.filter_file = None

        if self.server_lib != None:
            import urllib
            slu = "http://kryptonite.nbcr.net/ligand_props/" + server_lib + ".prop"
            pt = PropertyTable(url=slu)
            self.propertyTable = pt
            self.loc = self.server_lib
            self.loc_type = "server_lib"
        elif self.url_lib != None:
            import urllib
            slu = url_lib + "/lib.prop"
            #print slu
            pt = PropertyTable(url=slu)
            self.propertyTable = pt  
            self.loc = self.url_lib
            self.loc_type = "url_lib"
        
#        if accepted == None:
#            lfilter = FilterLigands()
#            accepted, rejected = lfilter.filterTable(self.propertyTable) 
#            self.SetAcceptedLigands(accepted.filenames)
#        else:
#            f = open(self.filter_file, 'w')
#            for i in self.accepted:
#                f.write(i + '''
#''')
#            f.close()

    def SetAcceptedLigands(self, accept_list, filter_file=None):
        if filter_file == None:
            workingDir = tempfile.mkdtemp()
            filter_file = workingDir + os.sep + 'filtered_ligands.txt'

        self.accepted = accept_list
        self.filter_file = filter_file

        f = open(self.filter_file, 'w')
        for i in self.accepted:
                f.write(i + '''
''')
        f.close()

    def GetPropertyFile(self):
        return self.property_file


