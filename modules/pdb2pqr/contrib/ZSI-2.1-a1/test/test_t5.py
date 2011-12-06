#!/usr/bin/env python
import unittest, multifile, mimetools
from ZSI import *
from ZSI import resolvers
from xml.dom import Node
#from xml.dom.ext.reader import PyExpat
from ZSI.parse import DefaultReader as Reader


try:
    import cStringIO as StringIO
except ImportError:
    import StringIO

class t5TestCase(unittest.TestCase):
    "Test case wrapper for old ZSI t5 test case"

    def checkt5(self):
        istr = StringIO.StringIO(intext)
        m = mimetools.Message(istr)
        if  m.gettype()[0:10] == "multipart/":
            cid = resolvers.MIMEResolver(m['content-type'], istr)
            xml = cid.GetSOAPPart()
            print 'xml=', xml.getvalue()
            for h,b in cid.parts:
                print h, b.read()
            dom = Reader.fromStream(xml)
            print dom

def makeTestSuite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(t5TestCase, "check"))
    return suite

def main():
    unittest.main(defaultTest="makeTestSuite")

intext = '''Content-Type: multipart/mixed; boundary="sep"
Subject: testing

skip this

--sep
Content-type: text/xml

<foo xmlns='www.zolera.com'>hello world</foo>
--sep
Content-Type: text/plain
content-id: <part111@example.zolera.com>

this is some plain text
--sep
content-type: application/skipme

do not see this
okay?
--sep
Content-Type: text/xml
Content-ID: <part2@example.zolera.com>

<xml>spoken</xml>
--sep--
'''

if __name__ == "__main__" : main()


