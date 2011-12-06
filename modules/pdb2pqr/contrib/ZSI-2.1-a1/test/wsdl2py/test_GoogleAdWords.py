#!/usr/bin/env python
############################################################################
# Joshua R. Boverhof, LBNL
# See LBNLCopyright for copyright notice!
###########################################################################
import os, sys, unittest
from ServiceTest import main, ServiceTestCase, ServiceTestSuite
from ZSI import FaultException
"""
Unittest for contacting google adwords

WSDL:  
"""

# General targets
def dispatch():
    """Run all dispatch tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(ServiceTest, 'test_dispatch'))
    return suite

def local():
    """Run all local tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(ServiceTest, 'test_local'))
    return suite

def net():
    """Run all network tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(ServiceTest, 'test_net'))
    return suite
    
def all():
    """Run all tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(ServiceTest, 'test_'))
    return suite


class TrafficEstimatorServiceTest(ServiceTestCase):
    """Test case for Google AdWords, sandbox v8
    Reads header information from a file "adwords.properties", need to format this for ConfigParser

[test_GoogleAdWords]
email = 
password = 
useragent = 
applicationtoken = 

    """
    name = "test_GoogleAdWords"
    client_file_name = "TrafficEstimatorService_client.py"
    types_file_name = "TrafficEstimatorService_types.py"
    server_file_name = "TrafficEstimatorService_server.py"

    header_info = os.path.join(os.getenv('HOME'), 'adwords.properties')

    def __init__(self, methodName):
        ServiceTestCase.__init__(self, methodName)
        self.wsdl2py_args.append('-b')

    def _get_soap_headers(self):
        from ConfigParser import ConfigParser
        cp = ConfigParser(); cp.read(self.header_info)
        p,e,a,u = map(lambda var: cp.get(self.__class__.name, var), 'password email applicationtoken useragent'.split())
        tns,GED = "https://adwords.google.com/api/adwords/v10", self.client_module.GED

        password = GED(tns, "password").pyclass(p)
        email = GED(tns, "email").pyclass(e)
        atoken = GED(tns, "applicationToken").pyclass(a)
        useragent = GED(tns, "useragent").pyclass(u)

        # google sandbox uses these conventions...
        dtoken = GED(tns, "developerToken").pyclass('%s++USD' %e) ## v8 sandbox syntax isnt working for v10
        cemail = GED(tns, "clientEmail").pyclass('client_1+'+e)

        return (email, password, useragent, dtoken, atoken, cemail)

    def test_net_KeywordEstimate(self):
        loc = self.client_module.TrafficEstimatorServiceLocator()
        port = loc.getTrafficEstimatorService(**self.getPortKWArgs())
        msg =  self.client_module.estimateKeywordListRequest()

        kwd =  msg.new_keywordRequests()
        kwd.Text = "flowers"
        kwd.MaxCpc = 50000L
        kwd.Type = "Broad"
        msg.KeywordRequests = [ kwd ]

        rsp = port.estimateKeywordList(msg, soapheaders=self._get_soap_headers())

ServiceTest = TrafficEstimatorServiceTest

if __name__ == "__main__" :
    main()

