#!/usr/bin/env python
import unittest, sys
from ZSI import *


class t7TestCase(unittest.TestCase):
    "Test case wrapper for old ZSI t7 test case"

    def checkt7(self):
        ps = ParsedSoap(text)

        tcdict = TC.Apache.Map('c-gensym1')
        tclist = TC.Apache.Map('c-gensym1', aslist=1)

        d = tcdict.parse(ps.body_root, ps)
        self.assertEqual(d, { u'a':123, '\x00\x01':456 })
        print 'as dictionary\n', d

        l = tclist.parse(ps.body_root, ps)
        self.assertEqual(l, [('\x00\x01', 456), (u'a', 123)])
        print '\n', '=' * 30
        print 'as list\n', l

        print '\n', '=' * 30
        sw = SoapWriter()
        sw.serialize(d, tcdict)
        print >>sys.stdout, sw

        print '\n', '=' * 30
        sw = SoapWriter()
        sw.serialize(l, tclist)
        print >>sys.stdout, sw

def makeTestSuite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(t7TestCase, "check"))
    return suite

def main():
    unittest.main(defaultTest="makeTestSuite")

text = '''
<SOAP-ENV:Envelope xmlns:SOAP-ENC="http://schemas.xmlsoap.org/soap/encoding/"
    xmlns:SOAP-ENV="http://schemas.xmlsoap.org/soap/envelope/"
    xmlns:xsi="http://www.w3.org/1999/XMLSchema-instance"
    xmlns:xsd="http://www.w3.org/1999/XMLSchema"
    xmlns:xmlsoap="http://xml.apache.org/xml-soap">
<SOAP-ENV:Body>
<c-gensym1 xsi:type="xmlsoap:Map">
  <item>
    <key xsi:type="SOAP-ENC:base64">AAE=</key>
    <value xsi:type="xsd:int">456</value>
  </item>
  <item>
    <key xsi:type="xsd:string">a</key>
    <value xsi:type="xsd:int">123</value>
  </item>
</c-gensym1>
</SOAP-ENV:Body>
</SOAP-ENV:Envelope>
'''


if __name__ == "__main__" : main()


