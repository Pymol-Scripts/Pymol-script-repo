#!/usr/bin/env python
############################################################################
# Joshua R. Boverhof, LBNL
# See LBNLCopyright for copyright notice!
###########################################################################
import sys, unittest
from ServiceTest import main, ServiceTestCase, ServiceTestSuite
import ZSI
from ZSI import FaultException
"""
Unittest for contacting the OPC XML-DA Service.

WSDL:  http://tswinc.us/XMLDADemo/ts_sim/OpcDaGateway.asmx?WSDL
"""

# General targets
def dispatch():
    """Run all dispatch tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(OPCServiceTest, 'test_dispatch'))
    return suite

def local():
    """Run all local tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(OPCServiceTest, 'test_local'))
    return suite

def net():
    """Run all network tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(OPCServiceTest, 'test_net'))
    return suite
    
def all():
    """Run all tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(OPCServiceTest, 'test_'))
    return suite


class OPCServiceTest(ServiceTestCase):
    """Test case for OPCService Web service
    
    def GetProperties(self, request):
    def Subscribe(self, request):
    def SubscriptionPolledRefresh(self, request):
    def SubscriptionCancel(self, request):
    def GetStatus(self, request):
    def Browse(self, request):
    def Read(self, request):
    def Write(self, request):
    """
    name = "test_OpcDaGateway"
    client_file_name = "OpcXmlDaSrv_client.py"
    types_file_name = "OpcXmlDaSrv_types.py"
    server_file_name = "OpcXmlDaSrv_server.py"

    def __init__(self, methodName):
        ServiceTestCase.__init__(self, methodName)
        self.wsdl2py_args.append('-b')

    def test_local_anyType(self):
        """serialize an int via anyType, then parse it back.
        """
        import time
        pyobj = self.client_module.ReadSoapOut()
        pyobj.RItemList = pyobj.new_RItemList()
        item = pyobj.RItemList.new_Items()
        pyobj.RItemList.Items = [item,]
        item.typecode.ofwhat[1].processContents = 'lax'
        item.Value = 123
        s = str(ZSI.SoapWriter().serialize(pyobj))
        
        ps = ZSI.ParsedSoap(s)
        pyobj = ps.Parse(pyobj.typecode)
        for item in pyobj.RItemList.Items:
            item.Value
        
    def test_net_Browse(self):
        """FaultException: The item path is not known to the server.
        """
        loc = self.client_module.OpcXmlDaSrvLocator()
        port = loc.getOpcXmlDaSrvSoap(**self.getPortKWArgs())
        
        msg = self.client_module.BrowseSoapIn()
        msg._PropertyNames=['Static']
        msg._attrs = {'ItemPath':'Static'}
        
        self.failUnless(\
            getattr(msg.typecode, 'attribute_typecode_dict', None) is not None,
             )

        self.failUnlessRaises(FaultException, port.Browse, msg)


#    def hangs_test_net_Read(self):
#        """FaultException: The item path is not known to the server.
#        """
#        msg = self.getInputMessageInstance('Read')
#        
#        #msg = ReadSoapIn()
#        op = msg.new_Options()
#        msg.Options = op
#        op.set_attribute_ReturnItemTime(True)
#        op.set_attribute_ReturnItemName(True)
#        op.set_attribute_ClientRequestHandle("")
#        op.set_attribute_LocaleID('en-us')
#
#        item_list = msg.new_ItemList()
#        msg.ItemList = item_list
#        item_list.set_attribute_MaxAge(1000)
#
#        item = item_list.new_Items()
#        item_list.Items = item
#        item.set_attribute_ItemPath("")
#        item.set_attribute_ItemName("Staic.Analog Types.Int")
#        item.set_attribute_ClientItemHandle("")
#        
#        self.failUnless(\
#            getattr(msg.typecode, 'attribute_typecode_dict', None) is not None,
#             )
#
#        self._ports[0].Read(msg)
#        

#    def hangs_test_net_GetProperties(self):
#        """sending an empty GetProperties request, 
#        receiving empty response.
#        
#        <GetPropertiesResult RcvTime="2005-12-16T13:23:55.8593750-05:00" 
#            ReplyTime="2005-12-16T13:24:00.1093750-05:00" 
#            RevisedLocaleID="en-us" ServerState="running" />
#        """
#        operationName = 'GetProperties'
#        msg = self.getInputMessageInstance(operationName)
#        msg._ItemIDs
#        msg._PropertyNames
#        response = self._ports[0].GetProperties(msg)
#        result = response._GetPropertiesResult
#        
#        # not sure these attributes are required but check for them.
#        self.failUnless(isinstance(getattr(result, '_attrs', None), dict))
#        for k in ['RcvTime','ReplyTime','RevisedLocaleID','ServerState']:
#            self.failUnless(result._attrs.has_key(k))
#            
#        self.failUnless(len(response._PropertyLists) == 0)
#        self.failUnless(len(response._Errors) == 0)
        



if __name__ == "__main__" :
    main()

