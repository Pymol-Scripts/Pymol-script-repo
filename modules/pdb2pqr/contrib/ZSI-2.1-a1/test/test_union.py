#!/usr/bin/env python
import unittest, sys, sha, base64
import ZSI 
from ZSI import _get_element_nsuri_name
from ZSI.schema import GED, TypeDefinition, ElementDeclaration
from ZSI.parse import ParsedSoap
from ZSI.wstools.c14n import Canonicalize
from ZSI.wstools.Namespaces import WSA200403, SOAP
from cStringIO import StringIO


# 
# Generated code
class ns3:
    class localPAssertionId_Dec(ElementDeclaration):
        literal = "localPAssertionId"
        schema = "http://www.pasoa.org/schemas/version024/PStruct.xsd"
        def __init__(self, **kw):
            kw["pname"] = ("http://www.pasoa.org/schemas/version024/PStruct.xsd","localPAssertionId")
            kw["aname"] = "_localPAssertionId"
            if ns3.LocalPAssertionId_Def not in ns3.localPAssertionId_Dec.__bases__:
                bases = list(ns3.localPAssertionId_Dec.__bases__)
                bases.insert(0, ns3.LocalPAssertionId_Def)
                ns3.localPAssertionId_Dec.__bases__ = tuple(bases)

            ns3.LocalPAssertionId_Def.__init__(self, **kw)
            if self.pyclass is not None: self.pyclass.__name__ = "localPAssertionId_Dec_Holder"


    class LocalPAssertionId_Def(ZSI.TC.Union, TypeDefinition):
        memberTypes = [(u'http://www.w3.org/2001/XMLSchema', u'long'), (u'http://www.w3.org/2001/XMLSchema', u'string'), (u'http://www.w3.org/2001/XMLSchema', u'anyURI')]
        schema = "http://www.pasoa.org/schemas/version024/PStruct.xsd"
        type = (schema, "LocalPAssertionId")
        def __init__(self, pname, **kw):
            ZSI.TC.Union.__init__(self, pname, **kw)



class UnionTestCase(unittest.TestCase):
    "test Union TypeCode"

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def check_union_long(self):
        import time
        typecode = GED("http://www.pasoa.org/schemas/version024/PStruct.xsd", "localPAssertionId")
        for value in (1234455, "whatever", "urn:whatever"):
            sw = ZSI.SoapWriter()
            sw.serialize(value, typecode)

            xml = str(sw)
            ps = ParsedSoap(xml)
            pyobj = ps.Parse(typecode)

            # Union Limitation:  
            #     currently it tries to parse it sequentially via memberTypes,
            #     so string is going to parse the URI when we want anyURI
            self.failUnless(value == pyobj, 'Expected equivalent')



def makeTestSuite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(UnionTestCase, "check"))
    return suite

def main():
    unittest.main(defaultTest="makeTestSuite")

if __name__ == '__main__': 
    main()

