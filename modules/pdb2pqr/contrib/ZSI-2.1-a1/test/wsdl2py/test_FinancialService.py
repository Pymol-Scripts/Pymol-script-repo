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
    name = "test_FinancialService"
    client_file_name = "FinancialService_client.py"
    types_file_name = "FinancialService_types.py"
    server_file_name = "FinancialService_server.py"

    def __init__(self, methodName):
        ServiceTestCase.__init__(self, methodName)
        self.wsdl2py_args.append('-b')
    
    def test_dispatch_getPV(self):
        loc = self.client_module.FinancialServiceLocator()
        port = loc.getFinancialService_Port(**self.getPortKWArgs())

        msg = self.client_module.getPVRequest()
        msg.Irate = 4
        msg.CFSequence = cfs = msg.new_CFSequence()
        cfs.CF = [100.0,5.0,5.0,105.0]

        rsp = port.getPV(msg)
        self.failUnless(rsp == 202.775091, "Received %d" %rsp)
        

if __name__ == "__main__" :
    main()
