#############################################################################
#
# Author: J. Ren
#
#############################################################################

"""
TBD
"""

import os

class gpf_template:

    def __init__(self, gpf_template_file=None):

        self.fullpath = os.path.abspath(gpf_template_file)
        self.basename = os.path.basename(gpf_template_file)
        self.id = self.basename.split('.')[0]
