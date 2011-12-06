import sys
from mglutil.regression import testplus
import test_htmlparser

harness = testplus.TestHarness( __name__,
                                funs = [],
                                dependents = [test_htmlparser.harness,
                                              ],
                                )

if __name__ == '__main__':
    print harness
    sys.exit( len( harness))
