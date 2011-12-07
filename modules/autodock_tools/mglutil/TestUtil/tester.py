# 
#
# $Id: tester.py,v 1.9 2004/08/27 17:48:00 sophiec Exp $
#

#########################################################################
#
# Date: July 2003  Author: Sophie Coon, William Lindstrom
#
#       sophiec@scripps.edu
#       lindy@scripps.edu
#
# Copyright: Michel Sanner, Sophie Coon, William Lindstrom and TSRI
#
#########################################################################


import unittest, sys
import types, os, glob, string

class TestSuite(unittest.TestSuite):
    def __init__(self, tests=(), setUpSuite=None, tearDownSuite=None):
        # Need to know what my tests contain.
        self.setUpSuite = setUpSuite
        self.tearDownSuite = tearDownSuite
        unittest.TestSuite.__init__(self, tests=tests)

    def __call__(self, result=None):
        if not self.setUpSuite is None:
##             if type(self.setUpSuite) is types.MethodType and len(self._tests):
##                 self._tests[1].setUpSuite()
##             else:
            self.setUpSuite()
        for test in self._tests:
            if result.shouldStop:
                break
            test(result)
        if not self.tearDownSuite is None:
##             if type(self.tearDownSuite) is types.MethodType and len(self._tests):
##                 self._tests[0].tearDownSuite()
##             else:
            self.tearDownSuite()

        return result
        
    def __exc_info(self):
        """Return a version of sys.exc_info() with the traceback frame
           minimised; usually the top level of the traceback frame is not
           needed.
        """
        exctype, excvalue, tb = sys.exc_info()
        if sys.platform[:4] == 'java': ## tracebacks look different in Jython
            return (exctype, excvalue, tb)
        newtb = tb.tb_next
        if newtb is None:
            return (exctype, excvalue, tb)
        return (exctype, excvalue, newtb)

class TestLoader(unittest.TestLoader):
    """
    
    """
    testMethodPrefix = 'test_'
    ignore = {}
    ignoredChecked = False
    suiteClass = TestSuite
    def loadTestsFromFunctions(self, functions, setUp=None, tearDown=None):
        """
        The functions needs to be from the same module.
        creates a FunctionTestCase for each function in the sequence and
        returns a TestSuite.
        """
        
        ftc = []
        if not type(functions) in [types.TupleType, types.ListType] and \
               type(functions) is types.FunctionType:
            functions = [functions,]
            
        m = functions[0].__module__
        modName = m.split('.')[-1]
        
        parts = m.split(".")[:-1]
        import string
        p = string.join(parts, '/')
        modPath = os.path.abspath(p)
        for func in functions:
            if not type(func) is types.FunctionType:continue
            ftc.append(unittest.FunctionTestCase(func, setUp=setUp,
                                                 tearDown=tearDown))
        return (modPath, self.suiteClass(ftc))

    def loadTestsFromModule(self, module, funcPrefix=None):
        modName = module.__name__.split('.')[-1]
        tests = []
        modPath = os.path.split(module.__file__)[0]
        modPath = os.path.abspath(modPath)
        ignoring = []
        if self.ignore.has_key(modName):
            # ignore the whole testModule
            ignoring = self.ignore[modName]
            if len(ignoring)==0:
                return (modPath, self.suiteClass(tests))
        if not funcPrefix is None:
            self.testMethodPrefix = funcPrefix

        # Look in the module if a setUp or setUpSuite is defined and a tearDown
        # and tearDownSuite 
        inModule = dir(module)
        if 'setUp' in inModule:
            setUp = getattr(module, 'setUp')
            if not type(setUp) is types.FunctionType:
                setUp=None
            inModule.remove('setUp')
        else:
            setUp = None

        if 'tearDown' in inModule:
            tearDown = getattr(module, 'tearDown')
            if not type(tearDown) is types.FunctionType:
                tearDown=None
            inModule.remove('tearDown')
        else:
            tearDown = None

        if 'setUpSuite' in inModule:
            setUpSuite = getattr(module, 'setUpSuite')
            if not type(setUpSuite) is types.FunctionType:
                setUpSuite=None
            inModule.remove('setUpSuite')
        else:
            setUpSuite = None

        if 'tearDownSuite' in inModule:
            tearDownSuite = getattr(module, 'tearDownSuite')
            if not type(tearDownSuite) is types.FunctionType:
                tearDownSuite=None
            inModule.remove('tearDownSuite')
        else:
            tearDownSuite = None

        testsFunc = []
        for name in dir(module):
            if name in ignoring: continue
            obj = getattr(module, name)
            if (isinstance(obj, (type, types.ClassType)) and
                issubclass(obj, unittest.TestCase)):
##                 inClass = dir(obj)
                # Look if a setUpSuite and a tearDownSuite have been implemented
                # for a testCase.else if one was implemented for the whole
                # module it will be used.
##                 if 'setUpSuite' in inClass:
##                     print 'in setUpSuite'
##                     setUpSuite = getattr(obj, 'setUpSuite')
##                     if not type(setUpSuite) is types.MethodType:
##                         setUpSuite=None

##                 if 'tearDownSuite' in inClass:
##                     tearDownSuite = getattr(obj, 'tearDownSuite')
##                     if not type(tearDownSuite) is types.MethodType:
##                         tearDownSuite=None
                ts = self.suiteClass(tests = map(obj,
                                                 self.getTestCaseNames(obj)),
                                     setUpSuite=setUpSuite, tearDownSuite=tearDownSuite)
                tests.append(ts)
                               
            elif type(obj) is types.FunctionType :
                p = len(self.testMethodPrefix)
                if name[:p]==self.testMethodPrefix:
                    testsFunc.append(unittest.FunctionTestCase(obj, setUp=setUp,
                                                               tearDown=tearDown))

        if len(testsFunc):
            ts = self.suiteClass(tests = testsFunc, setUpSuite=setUpSuite,
                                 tearDownSuite=tearDownSuite)
            tests.append(ts)
        
        return (modPath, self.suiteClass(tests=tests))
        
    def loadTestsFromPackageRecur(self, package, testDirName='Tests',
                                  modPrefix=None, funcPrefix=None):
        # Need to make sure that package is the proper thing
        pathPack = package.__path__[0]
        tests = []
        pName = package.__name__
        for root, dirs, files in os.walk(pathPack):
            if testDirName in dirs:
                #packName = root.replace("/", ".")
                packName = ""
                dir, name = os.path.split(root)
                while name != pName:
                    if packName: packName = name+"."+packName
                    else: packName = name
                    dir, name = os.path.split(dir)
                if packName:
                    packName = pName + "." + packName
                else:
                    packName = pName
                tests.append(self.loadTestsFromName(packName,
                                                    testDirName=testDirName,
                                                    modPrefix=modPrefix,
                                                    funcPrefix=funcPrefix))
        return tests

    def loadTestsFromPackage(self, package, testDirName="Tests",
                             modPrefix=None, funcPrefix=None):
        """
        import all the module from a the given package test directory,
        parse the __init__.py of the Tests directory to get the information
        on which tests to not run. __init__.py empty then takes all the
        pythonmodule with testName.py
        """
        # package : package
        # testDirName : string representing the name of the test directory
        # 1- Needs to get the test directory
        pathPack = package.__path__[0]
        pathTests = os.path.join(pathPack, testDirName)
        if not os.path.exists(pathTests):
            return
        testPackName = package.__name__ + "." + testDirName
        testPack = __import__(testPackName)
        components = testPackName.split('.')
        for comp in components[1:]:
            testPack = getattr(testPack, comp)

        if testPack.__dict__.has_key('ignore') :
            ignore = getattr(testPack, 'ignore')

        if modPrefix is None and testPack.__dict__.has_key('modPrefix'):
            modPrefix = getattr(testPack, "modPrefix")

        if funcPrefix is None and testPack.__dict__.has_key('funcPrefix'):
            funcPrefix = getattr(testPack,'funcPrefix')
            
        # Then need to go in the given directory and get all the python files
        # starting with the proper testMethodPrefix.
        # 2- get the __init__.py and parse the file.
        # Either use glob or walk.
        if modPrefix is None:
            modName = "/*.py"
        else:
            modName = "/%s*.py"%modPrefix
            
        testModules = glob.glob(pathTests+modName)
        ts = []
        for testMod in testModules:
            dir, file = os.path.split(testMod)
            if file in ["__init__.py", "mvAll.log.py"]: continue
            modName = os.path.splitext(file)[0]
            ts.append(self.loadTestsFromName(testPackName+"."+modName,
                                             funcPrefix=funcPrefix)[1])
        # 3- Create a suite of all the  tests in this module.
        packSuite = self.suiteClass(ts)
        return (pathTests, packSuite)

    def getTestsModulesFromPackRecur(self, pack, testDirName='Tests', modPrefix=None):
        pathPack = pack.__path__[0]
        testModules = []
        pName = pack.__name__
        for root, dirs, files in os.walk(pathPack):
            if testDirName in dirs:
                packName = ""
                dir, name = os.path.split(root)
                while name != pName:
                    if packName: packName = name+"."+packName
                    else: packName = name
                    dir, name = os.path.split(dir)
                if packName != "":
                    packName = pName + "." + packName
                else:
                    packName = pName
                pack = self.getObjFromName(packName, testDirName=testDirName)
                testModules = testModules + self.getTestsModulesFromPack(pack, testDirName=testDirName, modPrefix=modPrefix)
        return testModules

    def getTestsModulesFromPack(self, pack, testDirName='Tests', modPrefix=None):
        pathPack = pack.__path__[0]
        pathTests = os.path.join(pathPack, testDirName)
        
        if modPrefix is None:
            modName = "/*.py"
        else:
            modName = "/%s*.py"%modPrefix
        tm = glob.glob(pathTests+modName)
        testModules = []
        for testMod in tm:
            dir, file = os.path.split(testMod)
            if file in ["__init__.py", "mvAll.log.py"]: continue
            modName = os.path.splitext(file)[0]
            pName = pack.__name__+"."+testDirName+"."+modName
            testModules.append(pName)
        return testModules

    def getTestsModules(self, name, recursive=False, module=None, testDirName='Tests',
                        modPrefix=None, funcPrefix=None):
        if funcPrefix:
            self.testMethodPrefix=funcPrefix
        obj = self.getObjFromName(name, module=module, testDirName=testDirName)
        import unittest
        if type(obj) == types.ModuleType:
            # Can either  be a python module or a python package.
            if hasattr(obj,'__path__') and os.path.isdir(obj.__path__[0]):
                if recursive:
                    testModules = self.getTestsModulesFromPackRecur(obj,
                                                                    testDirName=testDirName,
                                                                    modPrefix=modPrefix) 
                    return testModules
                else:
                    testModules = self.getTestsModulesFromPack(obj,
                                                               testDirName=testDirName,
                                                               modPrefix=modPrefix)
                    return testModules
            else:
                return [obj.__name__,]

        elif (isinstance(obj, (type, types.ClassType)) and
              issubclass(obj, unittest.TestCase)):
            return [obj.__module__+'.'+obj.__name__,]

        elif type(obj) == types.FunctionType:
            return [obj.__module__+'.'+obj.__name__,]

    def getObjFromName(self, name, module=None, testDirName='Tests'):
        if name[-3:] == '.py':
            name = name[:-3]
        if '/' in name:
            parts = name.split('/')
            parts = filter(lambda x: x, parts)
        else:
            parts = name.split('.')
            parts = filter(lambda x: x, parts)
        if module is None:
            if not parts:
                raise ValueError, "incomplete test name: %s" % name
            else:
                parts_copy = parts[:]
                while parts_copy:
                    try:
                        module = __import__(string.join(parts_copy,'.'))
                        break
                    except ImportError:
                        del parts_copy[-1]
                        if not parts_copy: raise
                parts = parts[1:]
        obj = module

        for part in parts:
            obj = getattr(obj, part)
            if part==testDirName:
                if obj.__dict__.has_key('ignore'):
                    self.ignore = getattr(obj,'ignore')
        return obj
    
    def loadTestsFromName(self, name, recursive=False, module=None,
                          testDirName='Tests', modPrefix=None,
                          funcPrefix=None):
        """
        Returns a suite of all tests cases given a string specifier.
        The name may resolve either a package, a module, a test case class,
        a test method within a test case class, a test function or a callable
        object which returns a TestCase or TestSuite instance.

        The metod optionally resolves the names relative to a given module.
        """
        if funcPrefix:
            self.testMethodPrefix=funcPrefix

##         if name[-3:] == '.py':
##             name = name[:-3]
##         if '/' in name:
##             parts = name.split('/')
##         else:
##             parts = name.split('.')
        obj = self.getObjFromName(name, module=module, testDirName=testDirName)
                
        import unittest
        if type(obj) == types.ModuleType:
            # Can either  be a python module or a python package.
            if hasattr(obj,'__path__') and os.path.isdir(obj.__path__[0]):
                if recursive:
                    
                    return self.loadTestsFromPackageRecur(obj,
                                                          testDirName=testDirName,
                                                          modPrefix=modPrefix
                                                          )
                else:
                    return self.loadTestsFromPackage(obj,
                                                     testDirName=testDirName,
                                                     modPrefix=modPrefix)
            else:
                return self.loadTestsFromModule(obj)

        elif (isinstance(obj, (type, types.ClassType)) and
              issubclass(obj, unittest.TestCase)):
            m = obj.__module__
            parts = m.split(".")[:-1]
            p = string.join(parts, "/")
            return (p, self.loadTestsFromTestCase(obj))

        elif type(obj) == types.FunctionType:
            # need to get the setUp and tearDown method.
            m = obj.__module__
            module = __import__(m)
            parts = m.split('.')
            p = string.join(parts[:-1], '/')
            modPath = os.path.abspath(p)
            for part in parts[1:]:
                module = getattr(module , part)

            setUp = None
            tearDown = None
            setUpSuite = None
            tearDownSuite = None

            if module .__dict__.has_key('setUp'):
                setUp = getattr(module , 'setUp')
            if module .__dict__.has_key('tearDown'):
                tearDown = getattr(module , 'tearDown')

            if module .__dict__.has_key('setUpSuite'):
                setUpSuite = getattr(module, 'setUpSuite')
                if not type(setUpSuite) is types.FunctionType:
                    setUpSuite=None

            if  module .__dict__.has_key('tearDownSuite'):
                tearDownSuite = getattr(module, 'tearDownSuite')
                if not type(tearDownSuite) is types.FunctionType:
                    tearDownSuite=None
            
            tfc = unittest.FunctionTestCase(obj, setUp=setUp,
                                            tearDown=tearDown)
            ts = self.suiteClass(tests = [tfc,], setUpSuite=setUpSuite,
                                 tearDownSuite=tearDownSuite)
            return (modPath, ts)
            
        elif type(obj) == types.UnboundMethodType:
            newobj = obj.im_class(obj.__name__)
            m = newobj.__module__
            parts = m.split(".")[:-1]
            p = string.join(parts, "/")
           
            return (p, obj.im_class(obj.__name__))

        elif callable(obj):
            test = obj()
            if not isinstance(test, unittest.TestCase) and \
               not isinstance(test, unittest.TestSuite):
                raise ValueError, \
                      "calling %s returned %s, not a test" % (obj,test)
            return (None,test)
        else:
            raise ValueError, "don't know how to make test from: %s" % obj
        

