#!/usr/bin/env python

############################################################################
# Joshua R. Boverhof, LBNL
# See LBNLCopyright for copyright notice!
###########################################################################
import sys, unittest, time
from ServiceTest import main, ServiceTestCase, ServiceTestSuite, TestException
from ZSI.schema import ElementDeclaration, GED
from ZSI import ParsedSoap

"""
Unittest for contacting the Amazon ECommerce Service

WSDL: 

"""
# General targets
def dispatch():
    """Run all dispatch tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(AmazonTestCase, 'test_dispatch'))
    return suite

def local():
    """Run all local tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(AmazonTestCase, 'test_local'))
    return suite

def net():
    """Run all network tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(AmazonTestCase, 'test_net'))
    return suite
    
def all():
    """Run all tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(AmazonTestCase, 'test_'))
    return suite


class AmazonTestCase(ServiceTestCase):
    """Test case for AmazonS3 web service
    """
    name = "test_AmazonS3"
    client_file_name = "AmazonS3_client.py"
    types_file_name  = "AmazonS3_types.py"
    server_file_name = "AmazonS3_server.py"

    def __init__(self, methodName):
        ServiceTestCase.__init__(self, methodName)
        self.wsdl2py_args.append('-b')
        self.wsdl2py_args.append('--lazy')

    def test_local_import(self):
        pass
    
    def test_net_CreateBucket(self):
        loc = self.client_module.AmazonS3Locator()
        port = loc.getAmazonS3(**self.getPortKWArgs())

        msg = self.client_module.CreateBucketRequest()
        #msg.SubscriptionId = '0HP1WHME000749APYWR2'
        msg.Bucket = "HoneyPot"
        acl = msg.AccessControlList = msg.new_AccessControlList()
        grant = acl.new_Grant()
        acl.Grant = [grant]
        grant.Grantee = grant.new_Grantee()
        grant.Permission = grant.new_Permission("YES")

        msg.AWSAccessKeyId = '0HP1WHME000749APYWR2'
        msg.Timestamp = time.gmtime()
        msg.Signature = 'whatever'

        rsp = port.CreateBucket(msg)



if __name__ == '__main__':
    main()
