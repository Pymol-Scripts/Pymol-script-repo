#!/usr/bin/env python
############################################################################
# Joshua R. Boverhof, LBNL
# See Copyright for copyright notice!
###########################################################################
import unittest, warnings, os
from ZSI import version
from ZSI.wstools.logging import gridLog
from ServiceTest import main, CONFIG_PARSER, DOCUMENT, LITERAL, BROKE, TESTS

os.environ['GRIDLOG_ON'] = '1'
os.environ['GRIDLOG_DEST'] = "gridlog-udp://portnoy.lbl.gov:15100"

# General targets
def dispatch():
    """Run all dispatch tests"""
    return _dispatchTestSuite(broke=False)

def local():
    """Run all local tests"""
    return _localTestSuite(broke=False)
    
def net():
    """Run all network tests"""
    return _netTestSuite(broke=False)
    
def all():
    """Run all tests"""
    return _allTestSuite(broke=False)


# Specialized binding targets
def docLitTestSuite():
    """Run all doc/lit network tests"""
    return _netTestSuite(broke=False, document=True, literal=True)

def rpcLitTestSuite():
    """Run all rpc/lit network tests"""
    return _netTestSuite(broke=False, document=False, literal=True)

def rpcEncTestSuite():
    """Run all rpc/enc network tests"""
    return _netTestSuite(broke=False, document=False, literal=False)
    
    
# Low level functions
def _allTestSuite(document=None, literal=None, broke=None):
    return _makeTestSuite('all', document, literal, broke)

def _netTestSuite(document=None, literal=None, broke=None):
    return _makeTestSuite('net', document, literal, broke)

def _localTestSuite(document=None, literal=None, broke=None):
    return _makeTestSuite('local', document, literal, broke)
    
def _dispatchTestSuite(document=None, literal=None, broke=None):
    return _makeTestSuite('dispatch', document, literal, broke)

    
def _makeTestSuite(test, document=None, literal=None, broke=None):
    """Return a test suite containing all test cases that satisfy 
    the parameters. None means don't check.
    
    Parameters:
       test -- "net" run network tests, "local" run local tests,
           "dispatch" run dispatch tests, "all" run all tests.
       document -- None, True, False
       literal -- None, True, False
       broke -- None, True, False
    """
    assert test in ['net', 'local', 'dispatch', 'all'],(
        'test must be net, local, dispatch, or all')
    
    cp = CONFIG_PARSER
    testSections = []
    sections = [\
        'rpc_encoded' , 'rpc_encoded_broke',
        'rpc_literal', 'rpc_literal_broke', 'rpc_literal_broke_interop',
        'doc_literal', 'doc_literal_broke', 'doc_literal_broke_interop',
    ]
    boo = cp.getboolean
    for s,d,l,b in map(\
        lambda sec: \
            (sec, (None,boo(sec,DOCUMENT)), (None,boo(sec,LITERAL)), (None,boo(sec,BROKE))), sections):
        if document in d and literal in l and broke in b:
            testSections.append(s)
        
    suite = unittest.TestSuite()
    for section in testSections:
        moduleList = cp.get(section, TESTS).split()
        for module in  map(__import__, moduleList):
            def _warn_empty():
                warnings.warn('"%s" has no test "%s"' %(module, test))
                return unittest.TestSuite()
                
            s = getattr(module, test, _warn_empty)()
            suite.addTest(s)
    return suite

   
if __name__ == "__main__": 
    gridLog(prog="runTests.py", zsi="v%d.%d.%d" % version.Version, event="zsi.test.wsdl2py.runTests.ping")
    main()
    

