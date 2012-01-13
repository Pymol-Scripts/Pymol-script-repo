#!/usr/bin/env python
############################################################################
# Joshua R. Boverhof, LBNL
# See LBNLCopyright for copyright notice!
###########################################################################
import sys, unittest
from ServiceTest import main, ServiceTestCase, ServiceTestSuite
from ZSI.auth import AUTH
"""
Unittest for contacting the Map Point Service.  

WSDL:  
"""
# General targets
def dispatch():
    """Run all dispatch tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(MapPointTest, 'test_dispatch'))
    return suite

def local():
    """Run all local tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(MapPointTest, 'test_local'))
    return suite

def net():
    """Run all network tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(MapPointTest, 'test_net'))
    return suite
    
def all():
    """Run all tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(MapPointTest, 'test_'))
    return suite



class MapPointTest(ServiceTestCase):
    """Test case for OPCService Web service
    
    """
    name = "test_MapPoint"
    client_file_name = "CommonService_client.py"
    types_file_name = "CommonService_types.py"
    server_file_name = "CommonService_server.py"

    def __init__(self, methodName):
        ServiceTestCase.__init__(self, methodName)
        self.wsdl2py_args.append('-b')

    def test_net_GetVersionInfo(self):
        """expect this to fail cause i'm not doing http authentication.
        """
        loc = self.client_module.CommonServiceLocator()
        kw = self.getPortKWArgs()
        #port = loc.getCommonServiceSoap(auth=(AUTH.httpdigest, "USERNAME", "PASSWORD"), **kw)
        port = loc.getCommonServiceSoap(**kw)
        
        msg = self.client_module.GetVersionInfoSoapIn()
        try:
            rsp = port.GetVersionInfo(msg)
        except RuntimeError:
            # RuntimeError: HTTP Digest Authorization Failed
            pass

        port.binding.SetAuth(AUTH.httpdigest, user="USERNAME", password="PASSWORD")
        print ">> DIGEST AUTH"
        try:
            rsp = port.GetVersionInfo(msg)
        except RuntimeError:
            pass
        

if __name__ == "__main__" :
    main()

