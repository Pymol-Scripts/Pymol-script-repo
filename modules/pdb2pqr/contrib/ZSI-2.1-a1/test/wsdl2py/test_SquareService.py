#!/usr/bin/env python

############################################################################
# Joshua R. Boverhof, LBNL
# See LBNLCopyright for copyright notice!
# Tests for Holger's Square Service
#
###########################################################################
import sys, unittest
from ServiceTest import main, ServiceTestCase, ServiceTestSuite, TestException

"""
Unittest for contacting the SquareService rpc/literal tests.

From the paper "Interoperable WSDL/SOAP web services introduction: 
Python ZSI, Excel XP, gSOAP C/C++ & Applix SS", Holger Joukl

WSDL: SquareService.wsdl

"""
# General targets
def dispatch():
    """Run all dispatch tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(Test, 'test_dispatch'))
    return suite

def local():
    """Run all local tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(Test, 'test_local'))
    return suite

def net():
    """Run all network tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(Test, 'test_net'))
    return suite
    
def all():
    """Run all tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(Test, 'test_'))
    return suite


class Test(ServiceTestCase):
    """Test case for Holger's SquareService
    """
    name = "test_SquareService"
    client_file_name = "SquareService_client.py"
    types_file_name = "SquareService_types.py"
    server_file_name = "SquareService_server.py"

    def __init__(self, methodName):
        ServiceTestCase.__init__(self, methodName)
        self.wsdl2py_args.append('-b')
    
    def test_local_getSquare(self):
        from ZSI.writer import SoapWriter
        
    def test_dispatch_getSquare(self):
        loc = self.client_module.SquareServiceLocator()
        port = loc.getSquarePort(**self.getPortKWArgs())

        msg = self.client_module.getSquareRequest()
        msg.X = 4.0
        rsp = port.getSquare(msg)
        
        self.failUnless(rsp.Return == msg.X**2,
                        "Square Failed:  got %d, expecting %d" %(rsp.Return,msg.X**2))
        

if __name__ == "__main__" :
    main()
