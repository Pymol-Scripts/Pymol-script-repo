#!/usr/bin/env python
############################################################################
# Joshua R. Boverhof, LBNL
# See LBNLCopyright for copyright notice!
###########################################################################
import os, sys, unittest
from ServiceTest import main, ServiceTestCase, ServiceTestSuite
from ZSI import EvaluateException, FaultException
from ZSI.writer import SoapWriter
from ZSI.parse import ParsedSoap
from ZSI.TC import _get_type_definition as GTD
from ZSI.TC import _get_global_element_declaration as GED

"""
Unittest 

WSDL:  derivedTypes.
"""

# General targets
def dispatch():
    """Run all dispatch tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(DTTestCase, 'test_dispatch'))
    return suite

def local():
    """Run all local tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(DTTestCase, 'test_local'))
    return suite

def net():
    """Run all network tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(DTTestCase, 'test_net'))
    return suite
    
def all():
    """Run all tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(DTTestCase, 'test_'))
    return suite


class DTTestCase(ServiceTestCase):
    name = "test_DerivedTypes"
    client_file_name = None
    types_file_name  = "test_DerivedTypes_xsd_types.py"
    server_file_name = None

    def __init__(self, methodName):
        ServiceTestCase.__init__(self, methodName)
        self.wsdl2py_args.append('-x')
        self.wsdl2py_args.append('-b')

    def test_local_ged_substitution(self):
        """This test is designed to fail, trying to dump
        a GED in via type substitution.
        """
        self.types_module
        pyobj = GED('urn:test', 'test').pyclass()
        
        # use GED of a derived type
        pyobj.Actor = sub = GED('urn:test', 'MiddleActor').pyclass()
        sub.Element1 = 'foo'
        sub.Element2 = 'bar'
        
        sw = SoapWriter()
        self.failUnlessRaises(TypeError, sw.serialize, pyobj)
        
    def test_local_type_substitution_test2(self):
        """test extension of extension"""

        attr1 = 'aone'
        attr2 = 'atwo'
        attr3 = 'athree'
        self.types_module
        pyobj = GED('urn:test', 'test2').pyclass()

        # Test maxOccurs>1 for substitution 
        # 
        pyobj.Actor = [GTD('urn:test', 'TopActor')(None).pyclass()]
        sub1 = pyobj.Actor[0]
        sub1.Element1 = 'one'
        sub1.Element2 = 'two'
        sub1.Element3 = 'three'
        sub1.set_attribute_attr1(attr1)
        sub1.set_attribute_attr2(attr2)
        sub1.set_attribute_attr3(attr3)
        
        sw = SoapWriter()
        sw.serialize(pyobj)
        xml = str(sw)
        ps = ParsedSoap(xml)
        pyobj2 = ps.Parse(pyobj.typecode)
        sub2 = pyobj2.Actor[0]

        self.failUnless(sub2.get_attribute_attr1() == attr1, 'bad attribute 1')
        self.failUnless(sub2.get_attribute_attr2() == attr2, 'bad attribute 2')
        self.failUnless(sub2.get_attribute_attr3() == attr3, 'bad attribute 3')

        self.failUnless(sub2.Element1 == sub1.Element1, 'bad element 1')
        self.failUnless(sub2.Element2 == sub1.Element2, 'bad element 2')
        self.failUnless(sub2.Element3 == sub1.Element3, 'bad element 3')
                
        # check parsed out correct type
        self.failUnless(isinstance(sub2.typecode, sub1.typecode.__class__), 
            'local element actor "%s" must be an instance of "%s"'%
                (sub2.typecode, sub1.typecode.__class__))
        
        # check local element is derived from base
        base = GTD('urn:test', 'BaseActor')
        self.failUnless(isinstance(sub2.typecode, base), 
            'local element actor must be a derived type of "%s"'%
                base)

        
    def test_local_type_substitution2(self):
        """test extension of extension"""

        attr1 = 'aone'
        attr2 = 'atwo'
        attr3 = 'athree'
        self.types_module
        pyobj = GED('urn:test', 'test').pyclass()

        # [ 1489129 ] Unexpected subsitution error message
        #  try to parse before type ever initialized
        """
        ps = ParsedSoap(MSG1)
        pyobj0 = ps.Parse(pyobj.typecode)
        sub0 = pyobj0.Actor
        self.failUnless(sub0.get_attribute_attr1() == attr1, 'bad attribute1')
        self.failUnless(sub0.get_attribute_attr2() == attr2, 'bad attribute2')
        """

        # [ 1489090 ] Derived type attributes don't populate the attr dictionary
        # [ 1489677 ] Derivation from derived type missing derived element
        # 
        pyobj.Actor = sub1 = GTD('urn:test', 'TopActor')(None).pyclass()
        sub1.Element1 = 'one'
        sub1.Element2 = 'two'
        sub1.Element3 = 'three'
        sub1.set_attribute_attr1(attr1)
        sub1.set_attribute_attr2(attr2)
        sub1.set_attribute_attr3(attr3)
        
        sw = SoapWriter()
        sw.serialize(pyobj)
        xml = str(sw)
        ps = ParsedSoap(xml)
        pyobj2 = ps.Parse(pyobj.typecode)
        sub2 = pyobj2.Actor

        self.failUnless(sub2.get_attribute_attr1() == attr1, 'bad attribute 1')
        self.failUnless(sub2.get_attribute_attr2() == attr2, 'bad attribute 2')
        self.failUnless(sub2.get_attribute_attr3() == attr3, 'bad attribute 3')

        self.failUnless(sub2.Element1 == sub1.Element1, 'bad element 1')
        self.failUnless(sub2.Element2 == sub1.Element2, 'bad element 2')
        self.failUnless(sub2.Element3 == sub1.Element3, 'bad element 3')
                
        # check parsed out correct type
        self.failUnless(isinstance(sub2.typecode, sub1.typecode.__class__), 
            'local element actor "%s" must be an instance of "%s"'%
                (sub2.typecode, sub1.typecode.__class__))
        
        # check local element is derived from base
        base = GTD('urn:test', 'BaseActor')
        self.failUnless(isinstance(sub2.typecode, base), 
            'local element actor must be a derived type of "%s"'%
                base)

    def test_local_parse_missing_type_substitution(self):
        """attempt to substitute an unregistered/unknown type """
        attr1 = 'myclass'
        attr2 = 'whatever'
        self.types_module
        pyobj = GED('urn:test', 'test').pyclass()

        ps = ParsedSoap(NO_SUB_MSG)
        self.failUnlessRaises(EvaluateException, ps.Parse, pyobj.typecode)

    def test_local_type_substitution1(self):
        """test extension.   Parse known instance, serialize an equivalent, Parse it back. """
        attr1 = 'myclass'
        attr2 = 'whatever'
        self.types_module
        pyobj = GED('urn:test', 'test').pyclass()

        # [ 1489129 ] Unexpected subsitution error message
        #  try to parse before type ever initialized
        ps = ParsedSoap(MSG1)
        pyobj0 = ps.Parse(pyobj.typecode)
        sub0 = pyobj0.Actor
        self.failUnless(sub0.get_attribute_attr1() == attr1, 'bad attribute1')
        self.failUnless(sub0.get_attribute_attr2() == attr2, 'bad attribute2')

        # [ 1489090 ] Derived type attributes don't populate the attr dictionary
        # 
        pyobj.Actor = sub1 = GTD('urn:test', 'MiddleActor')(None).pyclass()
        sub1.Element1 = 'foo'
        sub1.Element2 = 'bar'
        sub1.set_attribute_attr1(attr1)
        sub1.set_attribute_attr2(attr2)
        
        sw = SoapWriter()
        sw.serialize(pyobj)
        xml = str(sw)
        ps = ParsedSoap(xml)
        pyobj2 = ps.Parse(pyobj.typecode)
        sub2 = pyobj2.Actor

        self.failUnless(sub2.get_attribute_attr1() == attr1, 'bad attribute class')
        self.failUnless(sub2.get_attribute_attr2() == attr2, 'bad attribute name')
                
        # check parsed out correct type
        self.failUnless(isinstance(sub2.typecode, sub1.typecode.__class__), 
            'local element actor "%s" must be an instance of "%s"'%
                (sub2.typecode, sub1.typecode.__class__))
        
        # check local element is derived from base
        base = GTD('urn:test', 'BaseActor')
        self.failUnless(isinstance(sub2.typecode, base), 
            'local element actor must be a derived type of "%s"'%
                base)
        

MSG1 = """<SOAP-ENV:Envelope xmlns:SOAP-ENC="http://schemas.xmlsoap.org/soap/encoding/" xmlns:SOAP-ENV="http://schemas.xmlsoap.org/soap/envelope/" xmlns:ZSI="http://www.zolera.com/schemas/ZSI/" xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"><SOAP-ENV:Header></SOAP-ENV:Header><SOAP-ENV:Body xmlns:ns1="urn:test"><ns1:test><actor attr1="myclass" attr2="whatever" xsi:type="ns1:MiddleActor"><element1 xsi:type="xsd:string">foo</element1><element2 xsi:type="xsd:string">bar</element2></actor></ns1:test></SOAP-ENV:Body></SOAP-ENV:Envelope>"""

NO_SUB_MSG = """<SOAP-ENV:Envelope xmlns:SOAP-ENC="http://schemas.xmlsoap.org/soap/encoding/" xmlns:SOAP-ENV="http://schemas.xmlsoap.org/soap/envelope/" xmlns:ZSI="http://www.zolera.com/schemas/ZSI/" xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"><SOAP-ENV:Header></SOAP-ENV:Header><SOAP-ENV:Body xmlns:ns1="urn:test"><ns1:test><actor attr1="myclass" attr2="whatever" xsi:type="ns1:Bogus"><element1 xsi:type="xsd:string">foo</element1><element2 xsi:type="xsd:string">bar</element2></actor></ns1:test></SOAP-ENV:Body></SOAP-ENV:Envelope>"""

if __name__ == "__main__" :
    main()

