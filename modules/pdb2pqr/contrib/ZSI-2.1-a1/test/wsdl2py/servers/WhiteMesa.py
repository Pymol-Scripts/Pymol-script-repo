#!/usr/bin/env python
############################################################################
# Joshua R. Boverhof, LBNL
# See LBNLCopyright for copyright notice!
###########################################################################
import sys
from ZSI.ServiceContainer import AsServer
from RPC_Literal_TestDefinitions_server import WhiteMesaSoapRpcLitTestSvc as WhiteMesa
"""
WhiteMesa web service for rpc/literal tests.

WSDL: http://www.whitemesa.net/wsdl/test-rpc-lit.wsdl

"""

class Service(WhiteMesa):
    def soap_echoStruct(self, ps):
        request,response = WhiteMesa.soap_echoStruct(self, ps)
        return request,response
    def soap_echoStructArray(self, ps):
        request,response = WhiteMesa.soap_echoStructArray(self, ps)
        return request,response
    def soap_echoStructAsSimpleTypes(self, ps):
        request,response = WhiteMesa.soap_echoStructAsSimpleTypes(self, ps)
        return request,response
    def soap_echoSimpleTypesAsStruct(self, ps):
        request,response = WhiteMesa.soap_echoSimpleTypesAsStruct(self, ps)
        return request,response
    def soap_echoNestedStruct(self, ps):
        request,response = WhiteMesa.soap_echoNestedStruct(self, ps)
        return request,response
    def soap_echoNestedArray(self, ps):
        request,response = WhiteMesa.soap_echoNestedArray(self, ps)
        return request,response
    def soap_echoStringArray(self, ps):
        request,response = WhiteMesa.soap_echoStringArray(self, ps)
        return request,response
    def soap_echoIntegerArray(self, ps):
        request,response = WhiteMesa.soap_echoIntegerArray(self, ps)
        return request,response
    def soap_echoBoolean(self, ps):
        request,response = WhiteMesa.soap_echoBoolean(self, ps)
        response._return = request._inputBoolean
        return request,response
    def soap_echoString(self, ps):
        request,response = WhiteMesa.soap_echoString(self, ps)
        return request,response

if __name__ == "__main__" :
    port = int(sys.argv[1])
    AsServer(port, (Service('test'),))


