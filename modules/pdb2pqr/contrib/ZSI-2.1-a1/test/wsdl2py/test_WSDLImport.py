#!/usr/bin/env python
############################################################################
# Joshua R. Boverhof, LBNL
# See LBNLCopyright for copyright notice!
###########################################################################
import os, sys, unittest, time
from ServiceTest import main, ServiceTestCase, ServiceTestSuite
from ZSI import FaultException
from ZSI.TC import _get_global_element_declaration as GED
from ZSI.writer import SoapWriter
from ZSI.parse import ParsedSoap

"""
Unittest for Bug Report 
[ ] 

test_WSDLImport.wsdl
test_WSDLImport2.wsdl

"""

# General targets
def dispatch():
    """Run all dispatch tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(WSDLImportTestCase, 'test_dispatch'))
    return suite

def local():
    """Run all local tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(WSDLImportTestCase, 'test_local'))
    return suite

def net():
    """Run all network tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(WSDLImportTestCase, 'test_net'))
    return suite
    
def all():
    """Run all tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(WSDLImportTestCase, 'test_'))
    return suite


class WSDLImportTestCase(ServiceTestCase):
    name = "test_WSDLImport"
    types_file_name = "OutSchemaTest_client.py"

    def __init__(self, methodName):
        ServiceTestCase.__init__(self, methodName)
        self.wsdl2py_args.append('-b')

    def test_local_attribute1(self):
        """
        """
        return

if __name__ == "__main__" :
    main()

