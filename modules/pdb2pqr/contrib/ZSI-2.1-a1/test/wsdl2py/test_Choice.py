#!/usr/bin/env python
############################################################################
# Joshua R. Boverhof, LBNL
# See LBNLCopyright for copyright notice!
###########################################################################
import os, sys, unittest
from ServiceTest import main, ServiceTestCase, ServiceTestSuite
from ZSI import FaultException
from ZSI.TC import _get_global_element_declaration as GED
from ZSI.writer import SoapWriter

"""
Unittest for Bug Report 
[ 1441574 ] ZSI assumes minOccurs(1) for all parts

WSDL:  
"""

# General targets
def dispatch():
    """Run all dispatch tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(ChoiceTestCase, 'test_dispatch'))
    return suite

def local():
    """Run all local tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(ChoiceTestCase, 'test_local'))
    return suite

def net():
    """Run all network tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(ChoiceTestCase, 'test_net'))
    return suite
    
def all():
    """Run all tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(ChoiceTestCase, 'test_'))
    return suite


class ChoiceTestCase(ServiceTestCase):
    name = "test_Choice"
    types_file_name = "test_Choice_xsd_types.py"
 
    def __init__(self, methodName):
        ServiceTestCase.__init__(self, methodName)
        self.wsdl2py_args.append('-b')
        self.wsdl2py_args.append('-x')
 
    def test_local_choice_default_facets_legal1(self):
        """<choice minOccurs=1 maxOccurs=1>
        """
        pyobj = GED("urn:example", "Easy").pyclass()
        pyobj.Rank = 1
        sw = SoapWriter()
        sw.serialize(pyobj)
        print str(sw)

    def test_local_choice_maxOccurs_unbounded(self):
        """<choice minOccurs=1 maxOccurs=unbounded>
        """
        pyobj = GED("urn:example", "Hard").pyclass()
        pyobj.Name = ["steve", "mark"]
        pyobj.Any = ["whatever"]
        pyobj.Rank = [2,3,4]
        sw = SoapWriter()
        sw.serialize(pyobj)
        print str(sw)


if __name__ == "__main__" :
    main()

