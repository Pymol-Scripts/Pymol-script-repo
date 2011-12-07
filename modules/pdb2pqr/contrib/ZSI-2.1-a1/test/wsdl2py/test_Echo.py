#!/usr/bin/env python
############################################################################
# Joshua R. Boverhof, LBNL
# See LBNLCopyright for copyright notice!
###########################################################################
import os, sys, unittest
from ServiceTest import main, ServiceTestCase, ServiceTestSuite
from ZSI import FaultException
"""
Unittest 

WSDL:  ../../samples/Echo/Echo.wsdl
"""

# General targets
def dispatch():
    """Run all dispatch tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(EchoTestCase, 'test_dispatch'))
    return suite

def local():
    """Run all local tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(EchoTestCase, 'test_local'))
    return suite

def net():
    """Run all network tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(EchoTestCase, 'test_net'))
    return suite
    
def all():
    """Run all tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(EchoTestCase, 'test_'))
    return suite


class EchoTestCase(ServiceTestCase):
    name = "test_Echo"
    client_file_name = "EchoServer_client.py"
    types_file_name  = "EchoServer_types.py"
    server_file_name = "EchoServer_server.py"

    def __init__(self, methodName):
        ServiceTestCase.__init__(self, methodName)
        self.wsdl2py_args.append('-b')

    def test_local_Echo(self):
        msg = self.client_module.EchoRequest()
        rsp = self.client_module.EchoResponse()

    def test_dispatch_Echo(self):
        loc = self.client_module.EchoServerLocator()
        port = loc.getEchoServer(**self.getPortKWArgs())
        
        msg = self.client_module.EchoRequest()
        msg.EchoIn = 'bla bla bla'
        rsp = port.Echo(msg)
        self.failUnless(rsp.EchoResult == msg.EchoIn, "Bad Echo")


if __name__ == "__main__" :
    main()

