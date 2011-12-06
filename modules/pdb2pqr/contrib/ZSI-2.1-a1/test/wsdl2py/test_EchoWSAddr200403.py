#!/usr/bin/env python
############################################################################
# Joshua R. Boverhof, LBNL
# See LBNLCopyright for copyright notice!
###########################################################################
import os, sys, unittest
from ServiceTest import main, ServiceTestCase, ServiceTestSuite
from ZSI import FaultException, TC
from ZSI.schema import GED, GTD

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
    name = "test_EchoWSAddr200403"
    client_file_name = "EchoWSAddr200403Server_client.py"
    types_file_name  = "EchoWSAddr200403Server_types.py"
    server_file_name = "EchoWSAddr200403Server_server.py"

    def __init__(self, methodName):
        ServiceTestCase.__init__(self, methodName)
        self.wsdl2py_args.append('-ab')

    def getPortKWArgs(self):
        kw = ServiceTestCase.getPortKWArgs(self)
        kw['wsAddressURI'] = 'http://schemas.xmlsoap.org/ws/2004/03/addressing'
        return kw

    def test_local_Echo(self):
        msg = self.client_module.EchoRequest()
        rsp = self.client_module.EchoResponse()

    def test_dispatch_Echo(self):
        loc = self.client_module.EchoWSAddr200403ServerLocator()
        port = loc.getport(**self.getPortKWArgs())
        
        msg = self.client_module.EchoRequest()
        msg.EchoIn = 'bla bla bla'
        rsp = port.Echo(msg)
        self.failUnless(rsp.EchoResult == msg.EchoIn, "Bad Echo")

    def test_dispatch_Echo_MIH_EPR(self):
        epr = GED('http://schemas.xmlsoap.org/ws/2004/03/addressing','EndpointReference').pyclass()
        epr.Address = 'urn:whatever'

        loc = self.client_module.EchoWSAddr200403ServerLocator()
        port = loc.getport(endPointReference=epr, **self.getPortKWArgs())

        msg = self.client_module.EchoRequest()
        msg.EchoIn = 1
        rsp = port.Echo(msg)
        self.failUnless(rsp.EchoResult == msg.EchoIn, "Bad Echo")

    def test_dispatch_Echo_MIH_EPR2(self):
        epr = GED('http://schemas.xmlsoap.org/ws/2004/03/addressing','EndpointReference').pyclass()
        epr.Address = 'urn:whatever'
        epr.ReferenceProperties = epr.new_ReferenceProperties()

        loc = self.client_module.EchoWSAddr200403ServerLocator()
        port = loc.getport(endPointReference=epr, **self.getPortKWArgs())

        msg = self.client_module.EchoRequest()
        msg.EchoIn = 1
        rsp = port.Echo(msg)
        self.failUnless(rsp.EchoResult == msg.EchoIn, "Bad Echo")

    def test_dispatch_Echo_MIH_EPR3_BadHeader(self):
        """Unqualified element "mystr" in Header
        """
        epr = GED('http://schemas.xmlsoap.org/ws/2004/03/addressing','EndpointReference').pyclass()
        epr.Address = 'urn:whatever'
        epr.ReferenceProperties = epr.new_ReferenceProperties()
        class Xstr(str): 
            typecode = TC.String('mystr')

        epr.ReferenceProperties.Any = [Xstr('whatever'),]

        loc = self.client_module.EchoWSAddr200403ServerLocator()
        self._setUpDispatch()
        port = loc.getport(endPointReference=epr, **self.getPortKWArgs())

        msg = self.client_module.EchoRequest()
        self.failUnlessRaises(FaultException, port.Echo,msg)

    def test_dispatch_Echo_MIH_EPR3(self):
        epr = GED('http://schemas.xmlsoap.org/ws/2004/03/addressing','EndpointReference').pyclass()
        epr.Address = 'urn:whatever'
        epr.ReferenceProperties = epr.new_ReferenceProperties()
        class Xstr(str): 
            typecode = TC.String(('urn:josh','mystr'))

        epr.ReferenceProperties.Any = [Xstr('whatever'),]

        loc = self.client_module.EchoWSAddr200403ServerLocator()
        self._setUpDispatch()
        port = loc.getport(endPointReference=epr, **self.getPortKWArgs())

        msg = self.client_module.EchoRequest()
        epr2 = GTD('http://schemas.xmlsoap.org/ws/2004/03/addressing','EndpointReferenceType')(None).pyclass()
        epr2.Address = epr.Address
        epr2.ReferenceProperties = epr.ReferenceProperties

        msg.EchoIn = epr2
        rsp = port.Echo(msg)
        self.failUnless(rsp.EchoResult.Address == msg.EchoIn.Address, "Bad Echo")
        self.failUnless(rsp.EchoResult.ReferenceProperties.Any == msg.EchoIn.ReferenceProperties.Any, "Bad Echo")


if __name__ == "__main__" :
    main()

