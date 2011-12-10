#############################################################################
#
# Author: J. Ren
#
#############################################################################

"""
TBD
"""

import os

class dpf_template:

    def __init__(self, dpf_template_file=None):

        self.fullpath = os.path.abspath(dpf_template_file)
        self.basename = os.path.basename(dpf_template_file)
        self.id = self.basename.split('.')[0]
