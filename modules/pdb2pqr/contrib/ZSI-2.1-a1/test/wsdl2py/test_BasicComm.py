#!/usr/bin/env python

import os, sys, unittest
from ServiceTest import main, ServiceTestCase, ServiceTestSuite
from ZSI import FaultException, Fault
from ConfigParser import ConfigParser, NoSectionError, NoOptionError
"""
Unittest 

WSDL:  BasicComm.wsdl
"""
CONFIG_FILE = 'config.txt'
CONFIG_PARSER = ConfigParser()
SECTION_DISPATCH = 'dispatch'

CONFIG_PARSER.read(CONFIG_FILE)

# General targets
def dispatch():
    """Run all dispatch tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(BasicCommTestCase, 'test_dispatch'))
    return suite

def local():
    """Run all local tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(BasicCommTestCase, 'test_local'))
    return suite

def net():
    """Run all network tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(BasicCommTestCase, 'test_net'))
    return suite
    
def all():
    """Run all tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(BasicCommTestCase, 'test_'))
    return suite


class BasicCommTestCase(ServiceTestCase):
    name = "test_BasicComm"
    client_file_name = "BasicServer_client.py"
    types_file_name  = "BasicServer_types.py"
    server_file_name = "BasicServer_server.py"

    def __init__(self, methodName):
        ServiceTestCase.__init__(self, methodName)
        self.wsdl2py_args.append('-b')

    def test_local_Basic(self):
        msg = self.client_module.BasicRequest()
        rsp = self.client_module.BasicResponse()

    def test_dispatch_Basic(self):
        loc = self.client_module.BasicServerLocator()
        port = loc.getBasicServer(**self.getPortKWArgs())
        
        msg = self.client_module.BasicRequest()
        msg._BasicIn = 'bla bla bla'
        rsp = port.Basic(msg)
        self.failUnless(rsp._BasicResult == msg._BasicIn, "Bad Echo")

        # test whether we get an HTTP response on a message with
        # no soap response.
        import httplib
        msg = u"""
            <SOAP-ENV:Envelope 
               xmlns:SOAP-ENC="http://schemas.xmlsoap.org/soap/encoding/" 
               xmlns:SOAP-ENV="http://schemas.xmlsoap.org/soap/envelope/" 
               xmlns:ZSI="http://www.zolera.com/schemas/ZSI/" 
               xmlns:xsd="http://www.w3.org/2001/XMLSchema" 
               xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
                 <SOAP-ENV:Header></SOAP-ENV:Header>
                 <SOAP-ENV:Body xmlns:ns1="urn:ZSI:examples">
                   <ns1:BasicOneWay>
                     <ns1:BasicIn>bla bla bla</ns1:BasicIn>
                   </ns1:BasicOneWay>
                 </SOAP-ENV:Body>
            </SOAP-ENV:Envelope>""".encode('utf-8')
        headers = {"Content-type": 'text/xml; charset="utf-8"', 'Content-Length': str(len(msg))}

        host = CONFIG_PARSER.get(SECTION_DISPATCH, 'host')
        port = CONFIG_PARSER.get(SECTION_DISPATCH, 'port')
        path = CONFIG_PARSER.get(SECTION_DISPATCH, 'path')

        conn = httplib.HTTPConnection("%s:%s" % (host, port))
        conn.request('POST', '/' + path, msg, headers)
        try:
            response = conn.getresponse()
        except httplib.BadStatusLine:
            conn.close()
            self.fail('No HTTP Response')

        conn.close()
        self.failUnless(response.status == 200, 'Wrong HTTP Result')

    def test_dispatch_BasicOneWay(self):
        loc = self.client_module.BasicServerLocator()
        port = loc.getBasicServer(**self.getPortKWArgs())
        
        msg = self.client_module.BasicOneWayRequest()
        msg.BasicIn = 'bla bla bla'
        rsp = port.BasicOneWay(msg)
        self.failUnless(rsp == None, "Bad One-Way")

    def test_dispatch_BasicOneWay_fault(self):
        """server will send back a soap:fault
        """
        loc = self.client_module.BasicServerLocator()
        port = loc.getBasicServer(**self.getPortKWArgs())
        
        msg = self.client_module.BasicOneWayRequest()
        msg.BasicIn = 'fault'
        self.failUnlessRaises(FaultException, port.BasicOneWay, msg)


if __name__ == "__main__" :
    main()

