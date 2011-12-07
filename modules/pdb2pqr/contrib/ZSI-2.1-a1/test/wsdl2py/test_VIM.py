#!/usr/bin/env python
############################################################################
# Joshua R. Boverhof, LBNL
# See LBNLCopyright for copyright notice!
###########################################################################
import os, sys, unittest
from ServiceTest import main, ServiceTestCase, ServiceTestSuite
from ZSI import FaultException, ParsedSoap, SoapWriter
"""
Unittest 

WSDL:  wsdl/vim.wsdl
"""

# General targets
def dispatch():
    """Run all dispatch tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(VIMTestCase, 'test_dispatch'))
    return suite

def local():
    """Run all local tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(VIMTestCase, 'test_local'))
    return suite

def net():
    """Run all network tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(VIMTestCase, 'test_net'))
    return suite
    
def all():
    """Run all tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(VIMTestCase, 'test_'))
    return suite


class VIMTestCase(ServiceTestCase):
    name = "test_VIM"
    client_file_name = "VIM_client.py"
    types_file_name  = "VIM_types.py"
    server_file_name = "VIM_server.py"

    def __init__(self, methodName):
        ServiceTestCase.__init__(self, methodName)
        self.wsdl2py_args.append('--lazy')
        self.wsdl2py_args.append('-b')

    def test_local_substitute_SessionManager(self):
        # BUG [ 1755740 ] Multiple calls to the same method
        MSG = """<?xml version="1.0" encoding="UTF-8"?>
<soapenv:Envelope xmlns:soapenc="http://schemas.xmlsoap.org/soap/encoding/"
 xmlns:soapenv="http://schemas.xmlsoap.org/soap/envelope/"
 xmlns:xsd="http://www.w3.org/2001/XMLSchema"
 xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
<soapenv:Body>
<RetrieveServiceContentResponse xmlns="urn:vim2">
  <returnval>
    <rootFolder type="Folder">group-d1</rootFolder>
    <propertyCollector type="PropertyCollector">propertyCollector</propertyCollector>
    <about>
      <name>VMware VirtualCenter</name>
      <fullName>VMware VirtualCenter 2.0.1 build-32042</fullName>
      <vendor>VMware, Inc.</vendor>
      <version>2.0.1</version>
      <build>32042</build>
      <localeVersion>INTL</localeVersion>
      <localeBuild>000</localeBuild>
      <osType>win32-x86</osType>
      <productLineId>vpx</productLineId>
      <apiType>VirtualCenter</apiType>
      <apiVersion>2.0.0</apiVersion>
    </about>
    <setting type="OptionManager">VpxSettings</setting>
    <userDirectory type="UserDirectory">UserDirectory</userDirectory>
    <sessionManager type="SessionManager">SessionManager</sessionManager>
    <authorizationManager type="AuthorizationManager">AuthorizationManager</authorizationManager>
    <perfManager type="PerformanceManager">PerfMgr</perfManager>
    <scheduledTaskManager type="ScheduledTaskManager">ScheduledTaskManager</scheduledTaskManager>
    <alarmManager type="AlarmManager">AlarmManager</alarmManager>
    <eventManager type="EventManager">EventManager</eventManager>
    <taskManager type="TaskManager">TaskManager</taskManager>
    <customizationSpecManager type="CustomizationSpecManager">CustomizationSpecManager</customizationSpecManager>
    <customFieldsManager type="CustomFieldsManager">CustomFieldsManager</customFieldsManager>
    <diagnosticManager type="DiagnosticManager">DiagMgr</diagnosticManager>
    <licenseManager type="LicenseManager">LicenseManager</licenseManager>
    <searchIndex type="SearchIndex">SearchIndex</searchIndex>
  </returnval>
</RetrieveServiceContentResponse>
</soapenv:Body>
</soapenv:Envelope>"""

        # Parse it out 
        ps = ParsedSoap(MSG)
        pyobj = ps.Parse( self.client_module.RetrieveServiceContentResponseMsg.typecode )
        sessionMgr = pyobj.Returnval.SessionManager

        # Serialize SessionManager in different context
        msg = self.client_module.LogoutRequestMsg()
        msg._this = sessionMgr
        SoapWriter().serialize(msg)

        # Parse it out: was failing
        # ZSI.EvaluateException: Element "__this" missing from complexType
        # [Element trace: /soapenv:Envelope/soapenv:Body/RetrieveServiceContentResponse/returnval]
        ps = ParsedSoap(MSG)
        pyobj = ps.Parse( self.client_module.RetrieveServiceContentResponseMsg.typecode )


if __name__ == "__main__" :
    main()

