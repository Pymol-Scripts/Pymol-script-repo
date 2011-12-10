#
#################################################################
#       Author: Sowjanya Karnati
#################################################################
#
#Purpose:To update dependencies list
#
# $Id: test_dependencies.py,v 1.5 2009/06/17 00:03:21 vareille Exp $
from mglutil.TestUtil.Tests.dependenciestest import DependencyTester
import unittest,sys
d = DependencyTester()
result_expected =[]
class test_dep(unittest.TestCase):
    
    def test_dep_1(self):
      if os.name != 'nt': #sys.platform != 'win32':
        result = d.rundeptester('mglutil')    
        if result !=[]:
            print "\nThe Following Packages are not present in CRITICAL or NONCRITICAL DEPENDENCIES of mglutil :\n  %s" %result
            self.assertEqual(result,result_expected) 
        else:
            self.assertEqual(result,result_expected)
    

if __name__ == '__main__':
    unittest.main()


