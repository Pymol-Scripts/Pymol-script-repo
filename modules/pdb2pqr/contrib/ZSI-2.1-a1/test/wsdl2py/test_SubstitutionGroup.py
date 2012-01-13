#!/usr/bin/env python
############################################################################
# Joshua R. Boverhof, LBNL
# See LBNLCopyright for copyright notice!
###########################################################################
import os, sys, unittest, time
from ServiceTest import main, ServiceTestCase, ServiceTestSuite
from ZSI import FaultException
from ZSI.TC import _get_global_element_declaration as GED
from ZSI.writer import SoapWriter
from ZSI.parse import ParsedSoap

"""
Unittest for substitutionGroup
[ ] 

XSD: 
"""

# General targets
def dispatch():
    """Run all dispatch tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(SubstitutionGroupTestCase, 'test_dispatch'))
    return suite

def local():
    """Run all local tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(SubstitutionGroupTestCase, 'test_local'))
    return suite

def net():
    """Run all network tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(SubstitutionGroupTestCase, 'test_net'))
    return suite
    
def all():
    """Run all tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(SubstitutionGroupTestCase, 'test_'))
    return suite


class SubstitutionGroupTestCase(ServiceTestCase):
    name = "test_SubstitutionGroup"
    types_file_name = "test_SubstitutionGroup_xsd_types.py"

    def __init__(self, methodName):
        ServiceTestCase.__init__(self, methodName)
        self.wsdl2py_args.append('-b')
        self.wsdl2py_args.append('-x')

    def test_local_attribute1(self):
        """
        """
        self.types_module

        xml = """<?xml version="1.0" encoding="UTF-8"?>
        <holder xmlns='urn:subGroup:types'>
        <baseElt><base>from base</base></baseElt>
        <childElt><base>from base</base><child>from child</child></childElt></holder>"""

        ps = ParsedSoap(xml, envelope=False)
        p1 = ps.Parse(GED("urn:subGroup:types", "holder"))

        b1 = p1.BaseElt[0]
        c1 = p1.BaseElt[1]

        sw = SoapWriter(envelope=False)
        sw.serialize(p1)

        ps = ParsedSoap(str(sw), envelope=False)
        p2 = ps.Parse(GED("urn:subGroup:types", "holder"))
        b2 = p2.BaseElt[0]
        c2 = p2.BaseElt[1]
 
        self.failUnlessEqual(b1.Base, b2.Base)
        self.failUnlessEqual(c1.Base, c2.Base)
        self.failUnlessEqual(c1.Child, c2.Child)



if __name__ == "__main__" :
    main()

