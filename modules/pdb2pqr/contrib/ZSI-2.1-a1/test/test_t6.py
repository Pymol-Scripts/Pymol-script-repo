#!/usr/bin/env python
import unittest, sys, multifile, mimetools, base64
from ZSI import *
from ZSI import resolvers
try:
    import cStringIO as StringIO
except ImportError:
    import StringIO

class t6TestCase(unittest.TestCase):
    "Test case wrapper for old ZSI t6 test case"

    def checkt6(self):
        try:
            istr = StringIO.StringIO(intext)
            m = mimetools.Message(istr)
            cid = resolvers.MIMEResolver(m['content-type'], istr)
            xml = cid.GetSOAPPart()
            ps = ParsedSoap(xml, resolver=cid.Resolve)
        except ParseException, e:
            print >>OUT, FaultFromZSIException(e).AsSOAP()
            self.fail()
        except Exception, e:
            # Faulted while processing; assume it's in the header.
            print >>OUT, FaultFromException(e, 1, sys.exc_info()[2]).AsSOAP() 
            self.fail()

        try:
            dict = ps.Parse(typecode)
        except Exception, e:
            # Faulted while processing; now it's the body
            print >>OUT, FaultFromException(e, 0, sys.exc_info()[2]).AsSOAP()
            self.fail()

        self.failUnlessEqual(dict['stringtest'], strExtTest, 
                            "Failed to extract stringtest correctly")
        print base64.encodestring(cid['partii@zolera.com'].read()) 
        v = dict['b64']
        print type(v), 'is type(v)' 
        self.failUnlessEqual(cid['partii@zolera.com'].getvalue(), v,
                                    "mismatch")
        print base64.encodestring(v)             
        from ZSI.wstools.c14n import Canonicalize 
        z = dict['xmltest'] 
        print type(z), z 
        print Canonicalize(z)

def makeTestSuite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(t6TestCase, "check"))
    return suite

def main():
    unittest.main(defaultTest="makeTestSuite")

OUT = sys.stdout
typecode = TC.Struct(None, [
                TC.String('b64'),
                TC.String('stringtest'),
                TC.XML('xmltest'),
            ])
                                                            
intext='''Return-Path: <rsalz@zolera.com>
Received: from zolera.com (os390.zolera.com [10.0.1.9])
        by zolera.com (8.11.0/8.11.0) with ESMTP id f57I2sf00832
        for <rsalz@zolera.com>; Thu, 7 Jun 2001 14:02:54 -0400
Sender: rsalz@zolera.com
Message-ID: <3B1FC1D1.FF6B21B4@zolera.com>
Date: Thu, 07 Jun 2001 14:02:57 -0400
From: Rich Salz <rsalz@zolera.com>
X-Mailer: Mozilla 4.72 [en] (X11; U; Linux 2.2.14-5.0 i686)
X-Accept-Language: en
MIME-Version: 1.0
To: rsalz@zolera.com
Subject: mime with attachments
Content-Type: multipart/mixed;
 boundary="------------68E4BAC5B266315E42428C64"
Status: R

This is a multi-part message in MIME format.
--------------68E4BAC5B266315E42428C64
Content-Type: text/plain; charset=us-ascii
Content-Transfer-Encoding: 7bit

<SOAP-ENV:Envelope
 xmlns="http://www.example.com/schemas/TEST"
 xmlns:SOAP-ENV="http://schemas.xmlsoap.org/soap/envelope/"
 xmlns:SOAP-ENC="http://schemas.xmlsoap.org/soap/encoding/"
 xmlns:ZSI="http://www.zolera.com/schemas/ZSI/">
<SOAP-ENV:Body>
<hreftest>
    <stringtest href="cid:part1@zolera.com"/>
    <b64 href="cid:partii@zolera.com"/>
    <xmltest href="cid:12@zolera.com"/>
</hreftest>
</SOAP-ENV:Body>
</SOAP-ENV:Envelope>
--------------68E4BAC5B266315E42428C64
Content-Type: text/plain; charset=us-ascii;
 name="abs.txt"
Content-Transfer-Encoding: 7bit
Content-ID: <part1@zolera.com>
Content-Disposition: inline;
 filename="abs.txt"


Digitial Signatures in a Web Services World

An influential Forrestor report created the term inverted security: it's
not about who you keep out, it's about who you let in.  Customer portals,
without a costly PKI deployment or application integration issues.

--------------68E4BAC5B266315E42428C64
Content-Type: application/pdf;
 name="gmpharma.pdf"
Content-Transfer-Encoding: base64
Content-ID: <partii@zolera.com>
Content-Disposition: inline;
 filename="gmpharma.pdf"

JVBERi0xLjINJeLjz9MNCjQzIDAgb2JqDTw8IA0vTGluZWFyaXplZCAxIA0vTyA0NSANL0gg
WyAxMTQ0IDM5NiBdIA0vTCA2NjkwMiANL0UgMTAyODIgDS9OIDkgDS9UIDY1OTI0IA0+PiAN
RB3nwVOQH9JpmFv6Ri2Zq7mlddSS2B5WcZwvAP+gy9QtuYlfqj1rsi9WqJOszzHXmXZ8fXxK
XBBztIpgbkRrd+SGtY4QXo0fX0VN86uKXwtrkd7h1qiq2FUtXl6uNfnCoyX1Dve1O3RPRyhG
sKn6fLMb+uSSIHPQkClRBwu5gechz/1PBUBSB34jXbPdMTIb+/wRP+pauSAhLBzFELDOgk5b
PaIPAnIudFovQTc7Df2Ws9Atz4Bua+oINphIOojogG5LP3Tb3oNu8bsmuK+wFXEdbfgFIx+G
gKULYx5A2WnaDXB5JeoRQg90S0HcX2dCPmRCqDXB/aX34KujsPwJ/UpRdxXPeAftDkQS6hag
bh/yTOiUyqBz9CzxnyMYQGDO0jrUZ47kkWfmYvVg
--------------68E4BAC5B266315E42428C64
Content-ID: <12@zolera.com>

<foo xmlns="example.com" xmlns:Z="zolera">
    this is a foo
    <b xmlns:Z="zolera">redundnant ns decl</b>
    <b Z:x="this was first" Z:a="2nd-orig">b test</b>
</foo>

--------------68E4BAC5B266315E42428C64--
'''

strExtTest = '''
Digitial Signatures in a Web Services World

An influential Forrestor report created the term inverted security: it's
not about who you keep out, it's about who you let in.  Customer portals,
without a costly PKI deployment or application integration issues.

'''

if __name__ == "__main__" : main()


