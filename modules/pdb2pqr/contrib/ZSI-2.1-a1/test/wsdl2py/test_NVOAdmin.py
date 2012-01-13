#!/usr/bin/env python
############################################################################
# Joshua R. Boverhof, LBNL
# See LBNLCopyright for copyright notice!
###########################################################################
import os, sys, unittest
from ServiceTest import main, ServiceTestCase, ServiceTestSuite
from ZSI import SoapWriter, ParsedSoap, TC,  FaultException
from ZSI import _child_elements, _get_element_nsuri_name, _is_element
from ZSI.schema import ElementDeclaration
from ZSI.wstools.Namespaces import SCHEMA

from xml.dom.ext.reader import PyExpat
from xml.dom import Node
"""
Unittest for NVO Admin.

WSDL:  
"""
#from ZSI.wstools.logging import setBasicLoggerDEBUG; setBasicLoggerDEBUG()

class schema(ElementDeclaration, TC.XML):
    """Create an element for dealing with <xsd:schema>
    """
    schema = SCHEMA.XSD3
    literal = "schema"

    def __init__(self, *args, **kw):
        # minOccurs=1, maxOccurs=1, nillable=False, encoded=kw.get("encoded")
        TC.XML.__init__(self, pname=(SCHEMA.XSD3, "schema"),  wrapped=False, **kw)


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


class NVOAdmin(ServiceTestCase):
    """Test case for NVO Admin

    """
    name = "test_NVOAdmin"
    client_file_name = "RegistryAdmin_client.py"
    types_file_name = "RegistryAdmin_types.py"
    #server_file_name = "RegistryAdmin_server.py"

    def __init__(self, methodName):
        ServiceTestCase.__init__(self, methodName)
        self.wsdl2py_args.append('-b')

    def test_local_serialize_schema(self):
        from ZSI import SoapWriter    
        from ZSI import _child_elements
        from xml.dom.ext.reader import PyExpat
        msg = self.client_module.DSQueryRegistrySoapOut()
        msg.DSQueryRegistryResult = msg.new_DSQueryRegistryResult()
        msg.DSQueryRegistryResult.Any = 'hi'

        input = open('wsdl/nvo-admin.wsdl').read()
        reader = PyExpat.Reader()
        dom = reader.fromString(input)

        dnode =  _child_elements(dom)[0]
        tnode =  _child_elements(dnode)[0]
        snode =  _child_elements(tnode)[0]

        msg.DSQueryRegistryResult.Schema = snode

        sw = SoapWriter()
        sw.serialize(msg)
        soap = str(sw)
        print soap        

        ps = ParsedSoap(soap)
        pyobj = ps.Parse(msg.typecode)
        self.failUnlessEqual(pyobj.DSQueryRegistryResult.Any, msg.DSQueryRegistryResult.Any)
        self.failUnless(_is_element(pyobj.DSQueryRegistryResult.Schema))
        print _get_element_nsuri_name(pyobj.DSQueryRegistryResult.Schema)
        self.failUnlessEqual(_get_element_nsuri_name(pyobj.DSQueryRegistryResult.Schema), (u'http://www.w3.org/2001/XMLSchema', u'schema'))
        

ServiceTest = NVOAdmin

if __name__ == "__main__" :
    main()

