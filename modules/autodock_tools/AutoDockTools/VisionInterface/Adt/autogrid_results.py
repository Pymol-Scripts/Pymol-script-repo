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

import os
from NetworkEditor.items import NetworkNode
from Vision.make_zip_file import make_zip_file

class autogrid_results:
    """
    Autogrid Result Object
    """
    
    def __init__(self, path=None, type=None):
        self.path = path
        self.type = type
    
        if self.type == "local" and os.path.isdir(path):
            zip_name = os.path.basename(path) + '.zip'
            self.path = make_zip_file(input_directory=path, output_name=zip_name)

         
