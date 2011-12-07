#!/usr/bin/env python
import unittest, sys
from ZSI import *
from ZSI import resolvers

OUT = sys.stdout
IN='''<SOAP-ENV:Envelope
 xmlns="http://www.example.com/schemas/TEST"
 xmlns:SOAP-ENV="http://schemas.xmlsoap.org/soap/envelope/"
 xmlns:SOAP-ENC="http://schemas.xmlsoap.org/soap/encoding/"
 xmlns:ZSI="http://www.zolera.com/schemas/ZSI/">
 <SOAP-ENV:Body>
 <hreftest>
    <xmltest href="http://www-itg.lbl.gov/~kjackson/zsi.xsd"/>
    <stringtest href="http://www.microsoft.com"/>
 </hreftest>
 </SOAP-ENV:Body>
 </SOAP-ENV:Envelope>
 '''

class t4TestCase(unittest.TestCase):
    "Test case wrapper for old ZSI t4 test case"

    def checkt4(self):
        try:
            r = resolvers.NetworkResolver(['http:'])
            ps = ParsedSoap(IN, resolver=r.Resolve)
        except ParseException, e:
            print >>OUT, FaultFromZSIException(e).AsSOAP()
            self.fail() 
        except Exception, e: 
            # Faulted while processing; assume it's in the header.  
            print >>OUT, FaultFromException(e, 1, sys.exc_info()[2]).AsSOAP() 
            self.fail() 
        print 'resolving' 
        typecode = TC.Struct(None, [ TC.XML('xmltest'), 
                           TC.String('stringtest', resolver=r.Opaque), ]) 
        try: 
            dict = ps.Parse(typecode) 
        except EvaluateException, e: 
            print >>OUT, FaultFromZSIException(e).AsSOAP() 
            self.fail() 
        except Exception, e: 
            # Faulted while processing; now it's the body 
            print >>OUT, FaultFromException(e, 0, sys.exc_info()[2]).AsSOAP() 
            self.fail() 
        ##PrettyPrint(dict['xmltest']) 
        print '**', dict['stringtest'], '**'

def makeTestSuite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(t4TestCase, "check"))
    return suite

def main():
    unittest.main(defaultTest="makeTestSuite")


if __name__ == "__main__" : main()


