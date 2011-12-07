#!/usr/bin/env python
############################################################################
# Joshua R. Boverhof, LBNL
# See LBNLCopyright for copyright notice!
# Tests for Holger's DateService
###########################################################################
import sys, unittest, time
from ServiceTest import main, ServiceTestCase, ServiceTestSuite, TestException

"""
Unittest for contacting the DateService rpc/literal tests.

From the paper "Interoperable WSDL/SOAP web services introduction: 
Python ZSI, Excel XP, gSOAP C/C++ & Applix SS", Holger Joukl

WSDL: DateService.wsdl

"""
# General targets
def dispatch():
    """Run all dispatch tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(Test, 'test_dispatch'))
    return suite

def local():
    """Run all local tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(Test, 'test_local'))
    return suite

def net():
    """Run all network tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(Test, 'test_net'))
    return suite
    
def all():
    """Run all tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(Test, 'test_'))
    return suite


class Test(ServiceTestCase):
    """Test case for Holger's DateService
    """
    name = "test_DateService"
    client_file_name = "DateService_client.py"
    types_file_name = "DateService_types.py"
    server_file_name = "DateService_server.py"

    def __init__(self, methodName):
        ServiceTestCase.__init__(self, methodName)
        self.wsdl2py_args.append('-b')
    
    #def test_local_getDate(self):
    #    from ZSI.writer import SoapWriter
        
    def test_dispatch_getCurrentDate_getDate(self):
        offset = 9
        loc = self.client_module.simple_Date_ServiceLocator()
        port = loc.getDateService_Port(**self.getPortKWArgs())
        print "START"
        msg = self.client_module.getCurrentDateRequest()
        msg.Input = "Test"
        rsp = port.getCurrentDate(msg)
        
        today = rsp.Today
        today.Month
        today.Day
        today.Hour
        today.Minute
        today.Second
        today.Weekday
        today.DayOfYear
        today.Dst

        dateRequest = self.client_module.getDateRequest()
        # We use the current date as input to getDate
        dateRequest.Someday = today
        dateRequest.Offset = offset
        date = port.getDate(dateRequest)

        print '\n\nRESULT'
        print '%10s = %s' % ('today', _make_asctime(today))
        print '%6s + %d = %s' % ('today', dateRequest.Offset, _make_asctime(date.Day))


def _make_asctime(date_object):
    timeTuple = (date_object._year, date_object._month, date_object._day,
                 date_object._hour, date_object._minute, date_object._second,
                 date_object._weekday, date_object._dayOfYear, date_object._dst
                 )
    return time.asctime(timeTuple)

        

if __name__ == "__main__" :
    main()
