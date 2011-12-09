#
#
# $Id: test_tester.py,v 1.3 2003/08/29 18:09:15 sophiec Exp $
#

#########################################################################
#
# Date: July 2003  Author: Sophie Coon
#
#       sophiec@scripps.edu
#
# Copyright: TSRI, Sophie Coon.
#
#########################################################################

import unittest, os


class TestSuiteTest(unittest.TestCase):
    pass

class TestLoaderTest(unittest.TestCase):
    """
    TestCase class implementing methods to test Tester functionalities
    """
    def setUp(self):
        from mglutil.TestUtil.tester import TestLoader
        self.tl = TestLoader()

    def tearDown(self):
        self.tl = None

    def test_loadTestsFromTestCase1(self):
        from Data.test_displayCommands import DisplayLinesTest
        suite = self.tl.loadTestsFromTestCase(DisplayLinesTest)
        self.failUnless(isinstance(suite, self.tl.suiteClass))
        self.assertEqual(len(suite._tests), 5)

    
##     def test_loadTestsFromFunctions(self):
##         from Data import test_secondarystructure
##         mod = test_secondarystructure
##         import types
##         funcs = []
##         for name in dir(test_secondarystructure):
##             if name[:4]=='test':
##                 f = getattr(mod,name)
##                 if type(f) is types.FunctionType:
##                     funcs.append(f)
                    
##         setUp = getattr(mod,"startMoleculeViewer")
##         tearDown = getattr(mod,"quitMoleculeViewer")
##         path,testSuite = self.tl.loadTestsFromFunctions(funcs,setUp=setUp,
##                                                    tearDown=tearDown)
##         self.failUnless(isinstance(testSuite, self.tl.suiteClass))
##         self.assertEqual(len(testSuite._tests),19)
        
##         path, testSuite = self.tl.loadTestsFromFunctions(funcs[0])
##         self.failUnless(isinstance(testSuite, self.tl.suiteClass))
##         self.assertEqual(len(testSuite._tests),1)

    def test_loadTestsFromModule_1(self):
        from Data import test_secondarystructure
        mod = test_secondarystructure
        path,testSuite = self.tl.loadTestsFromModule(mod)
        self.failUnless(isinstance(testSuite, self.tl.suiteClass))
        self.assertEqual(len(testSuite._tests), 19)

    def test_loadTestsFromModule_2(self):
        from Data import test_displayCommands
        mod = test_displayCommands
        path,testSuite = self.tl.loadTestsFromModule(mod)
        
    def test_loadTestsFromName_1(self):
        name = 'mglutil.TestUtil.Tests.test_tester'
        path,ts = self.tl.loadTestsFromName(name)
        self.failUnless(isinstance(ts, self.tl.suiteClass))
        self.failUnless(len(ts._tests)!=0)
        
    def test_loadTestsFromName_2(self):
        name = 'mglutil.TestUtil'
        path,ts = self.tl.loadTestsFromName(name)
        self.failUnless(isinstance(ts, self.tl.suiteClass))
        self.failUnless(len(ts._tests)!=0)

    def test_loadTestsFromPackage_1(self):
        from mglutil import TestUtil
        print os.path.abspath(TestUtil.__path__[0])
        r = self.tl.loadTestsFromPackage(TestUtil)
        print r
        #self.failUnless(isinstance(ts, self.tl.suiteClass))
        #self.failUnless(len(ts._tests)!=0)

        
class TesterTest(unittest.TestCase):
    pass


if __name__ == '__main__':
    unittest.main()
