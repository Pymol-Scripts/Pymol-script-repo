#!/usr/bin/env python
############################################################################
# Joshua R. Boverhof, LBNL
# See LBNLCopyright for copyright notice!
###########################################################################
import os, sys, unittest
from ServiceTest import main, ServiceTestCase, ServiceTestSuite
from ZSI import FaultException

"""
Unittest for contacting Clearspace blog webservice

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


# NEED TO CREATE WSSE typecodes
from ZSI.generate.commands import wsdl2py
if not os.path.isdir('stubs'): os.makedirs('stubs')
wsdl2py(['--complexType', '--schema','--output-dir=stubs', 'http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-secext-1.0.xsd'])


class BlogServiceTest(ServiceTestCase):
    """Test case for Clearspace sandbox, example how to use client WSSE:Security UsernameToken Profile
    
<wsdl:Envelope xmlns:soap="..." xmlns:wsse="..." >
   <wsdl:Header>
      <wsse:Security>
         <wsse:UsernameToken>
            <wsse:Username>admin</wsse:Username>
            <wsse:Password>password</wsse:Password>
         </wsse:UsernameToken>
      </wsse:Security>
   </wsdl:Header>
</wsdl:Envelope>

    """
    name = "test_Clearspace"
    client_file_name = "BlogService_client.py"
    types_file_name = "BlogService_types.py"
    server_file_name = "BlogService_server.py"

    def __init__(self, methodName):
        ServiceTestCase.__init__(self, methodName)
        self.wsdl2py_args.append('-b')

    def _get_soap_headers(self):
        import oasis_200401_wss_wssecurity_secext_1_0_xsd_types
        from ZSI.schema import GED
        security = GED("http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-secext-1.0.xsd", "Security").pyclass()
        token = GED("http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-secext-1.0.xsd", "UsernameToken").pyclass()
        security.Any = [token]
        token.Username = 'billy'
        klass = GED("http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-secext-1.0.xsd", "Password").pyclass
        token.Any = [klass('guest'),]

        return (security,)

    def test_net_Blogcount(self):
        loc = self.client_module.BlogServiceLocator()
        msg = self.client_module.getBlogCountRequest()
        port = loc.getBlogServiceHttpPort(**self.getPortKWArgs())
        rsp = port.getBlogCount(msg, soapheaders=self._get_soap_headers(),)

    def test_local_(self):
        import oasis_200401_wss_wssecurity_secext_1_0_xsd_types
        return

ServiceTest = BlogServiceTest

if __name__ == "__main__" :
    main()

