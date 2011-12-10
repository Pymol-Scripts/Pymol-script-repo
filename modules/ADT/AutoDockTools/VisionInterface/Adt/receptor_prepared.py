#############################################################################
#
# Author: J. Ren
#
#############################################################################

"""
TBD
"""

import os

class receptor_prepared:

    def __init__(self, receptor_path=None):

        if receptor_path.startswith('http://'):
            self.type = 'url'
            self.path = receptor_path
        else:
            self.type = 'local'
            self.path = os.path.abspath(receptor_path)

        self.basename = os.path.basename(receptor_path)
        self.id = self.basename.rstrip('.pdbqt')
