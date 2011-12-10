#############################################################################
#
# Author: Sophie I. COON, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################
#
# $Header: /opt/cvs/python/packages/share1.5/MolKit/Tests/test_molecule.py,v 1.2 2003/08/29 17:50:16 sophiec Exp $
#
# $Id: test_molecule.py,v 1.2 2003/08/29 17:50:16 sophiec Exp $
#
import sys
from mglutil.regression import testplus
"""
This module implements a set of function to test the functionalities
implemented in the molecule module.
"""

harness = testplus.TestHarness( __name__,
                                funs = testplus.testcollect( globals()),
                                )

if __name__ == '__main__':
    print harness
    sys.exit( len( harness))

