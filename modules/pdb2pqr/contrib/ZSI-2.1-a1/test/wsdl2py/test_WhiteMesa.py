#!/usr/bin/env python

############################################################################
# Joshua R. Boverhof, LBNL
# See LBNLCopyright for copyright notice!
###########################################################################
import sys, unittest
from ServiceTest import main, ServiceTestCase, ServiceTestSuite, TestException

"""
Unittest for contacting the WhiteMesa web service for rpc/literal tests.

WSDL: http://www.whitemesa.net/wsdl/test-rpc-lit.wsdl

"""
# General targets
def dispatch():
    """Run all dispatch tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(WhiteMesaTest, 'test_dispatch'))
    return suite

def local():
    """Run all local tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(WhiteMesaTest, 'test_local'))
    return suite

def net():
    """Run all network tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(WhiteMesaTest, 'test_net'))
    return suite
    
def all():
    """Run all tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(WhiteMesaTest, 'test_'))
    return suite


class WhiteMesaTest(ServiceTestCase):
    """Test case for ZipCodeResolver Web service
    """
    name = "test_WhiteMesa"
    client_file_name = "RPC_Literal_TestDefinitions_client.py"
    types_file_name = "RPC_Literal_TestDefinitions_types.py"
    server_file_name = "RPC_Literal_TestDefinitions_server.py"

    def __init__(self, methodName):
        ServiceTestCase.__init__(self, methodName)
        self.wsdl2py_args.append('-b')
    
    def test_local_EchoBoolean(self):
        from ZSI.writer import SoapWriter
        msg = self.client_module.echoBooleanRequest()
        msg._inputBoolean = True
        sw = SoapWriter()
        sw.serialize(msg)
        
    def test_net_EchoBoolean(self):
        msg = self.client_module.echoBooleanRequest()
        msg._inputBoolean = True
        
        loc = self.client_module.WhiteMesaSoapRpcLitTestSvcLocator()
        port = loc.getSoap11TestRpcLitPort(**self.getPortKWArgs())
        rsp = port.echoBoolean(msg)
        
        self.failUnless(msg._inputBoolean == rsp._return,
                        "EchoBoolean Failed")
        
    def test_dispatch_EchoBoolean(self):
        self.test_net_EchoBoolean()
    

if __name__ == "__main__" :
    main()
