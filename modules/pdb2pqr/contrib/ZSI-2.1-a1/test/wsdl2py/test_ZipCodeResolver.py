#!/usr/bin/env python

############################################################################
# David W. Robertson, LBNL
# See LBNLCopyright for copyright notice!
###########################################################################
import sys, unittest
from ServiceTest import main, ServiceTestCase, ServiceTestSuite

"""
Unittest for contacting the ZipCodeResolver Web service.

WSDL: http://webservices.eraserver.net/zipcoderesolver/zipcoderesolver.asmx?WSDL

"""

# General targets
def dispatch():
    """Run all dispatch tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(ZipCodeResolverTest, 'test_dispatch'))
    return suite

def local():
    """Run all local tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(ZipCodeResolverTest, 'test_local'))
    return suite

def net():
    """Run all network tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(ZipCodeResolverTest, 'test_net'))
    return suite
    
def all():
    """Run all tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(ZipCodeResolverTest, 'test_'))
    return suite


class ZipCodeResolverTest(ServiceTestCase):
    """Test case for ZipCodeResolver Web service
    """
    name = "test_ZipCodeResolver"
    client_file_name = "ZipCodeResolver.py"
    types_file_name = "ZipCodeResolver_types.py"
    #server_file_name = "ZipCodeResolver_server.py"
    
    def test_net_CorrectedAddressHtml(self):
        loc = self.client_module.ZipCodeResolverLocator()
        port = loc.getZipCodeResolverSoap(**self.getPortKWArgs())
        
        msg = self.client_module.CorrectedAddressHtmlSoapIn()
        msg._address = '636 Colusa Avenue'
        msg._city = 'Berkeley'
        msg._state = 'California'
        rsp = port.CorrectedAddressHtml(msg)

    def test_net_CorrectedAddressXml(self):
        loc = self.client_module.ZipCodeResolverLocator()
        port = loc.getZipCodeResolverSoap(**self.getPortKWArgs())
        
        msg = self.client_module.CorrectedAddressXmlSoapIn()
        msg._address = '636 Colusa Avenue'
        msg._city = 'Berkeley'
        msg._state = 'California'
        rsp = port.CorrectedAddressXml(msg)
         
    def test_net_FullZipCode(self):
        loc = self.client_module.ZipCodeResolverLocator()
        port = loc.getZipCodeResolverSoap(**self.getPortKWArgs())
        
        msg = self.client_module.FullZipCodeSoapIn()
        msg._address = '636 Colusa Avenue'
        msg._city = 'Berkeley'
        msg._state = 'California'
        rsp = port.FullZipCode(msg)
    
    def test_net_ShortZipCode(self):
        loc = self.client_module.ZipCodeResolverLocator()
        port = loc.getZipCodeResolverSoap(**self.getPortKWArgs())
        
        msg = self.client_module.ShortZipCodeSoapIn()
        msg._address = '636 Colusa Avenue'
        msg._city = 'Berkeley'
        msg._state = 'California'
        rsp = port.ShortZipCode(msg)
    
    def test_net_VersionInfo(self):
        loc = self.client_module.ZipCodeResolverLocator()
        port = loc.getZipCodeResolverSoap(**self.getPortKWArgs())
        
        msg = self.client_module.VersionInfoSoapIn()
        rsp = port.VersionInfo(msg)


if __name__ == "__main__" :
    main()

