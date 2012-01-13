#!/usr/bin/env python

############################################################################
# Joshua R. Boverhof
# See LBNLCopyright for copyright notice!
###########################################################################
import sys, unittest
from ServiceTest import main, ServiceTestCase, ServiceTestSuite


"""
Unittest for contacting the threatService Web service.

WSDL:  http://www.boyzoid.com/threat.cfc?wsdl
"""

# General targets
def dispatch():
    """Run all dispatch tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(HomelandTestCase, 'test_dispatch'))
    return suite

def local():
    """Run all local tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(HomelandTestCase, 'test_local'))
    return suite

def net():
    """Run all network tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(HomelandTestCase, 'test_net'))
    return suite
    
def all():
    """Run all tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(HomelandTestCase, 'test_'))
    return suite


class HomelandTestCase(ServiceTestCase):
    """Test case for ZipCodeResolver Web service
    """
    name = "test_ThreatService"
    client_file_name = "Current_Homeland_Security_Threat_Level_client.py"
    types_file_name = "Current_Homeland_Security_Threat_Level_types.py"
    server_file_name = None

    def __init__(self, methodName):
        ServiceTestCase.__init__(self, methodName)
        self.wsdl2py_args.append('-b')
    
    def test_net_threatLevel(self):
        loc = self.client_module.Current_Homeland_Security_Threat_LevelLocator()
        port = loc.getthreat_cfc(**self.getPortKWArgs())

        msg = self.client_module.threatLevelRequest()
        rsp = port.threatLevel(msg)
        for item in rsp.ThreatLevelReturn.Item:
            item.Key
            item.Value
        
    

if __name__ == "__main__" :
    main()
