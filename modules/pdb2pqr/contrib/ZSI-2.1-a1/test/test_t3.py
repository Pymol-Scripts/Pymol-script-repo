#!/usr/bin/env python
import unittest
from ZSI import *
from ZSI.wstools.logging import setBasicLoggerDEBUG
setBasicLoggerDEBUG()

class t3TestCase(unittest.TestCase):
    "Test case wrapper for old ZSI t3 test case"

    def checkt3(self):
        a = None
        try: 
            3 / 0
        except Exception, e:
            a = e
        f = FaultFromException(a, 0)
        text = f.AsSOAP()
        i = 0
        for l in text.split('\n'):
            print i, l
            i += 1
        ps = ParsedSoap(text)
        if ps.IsAFault():
            f = FaultFromFaultMessage(ps)
            print f.AsSOAP()
            self.failUnless(f.AsSOAP().find(str(a)) > 0)
        print '--'*20


def makeTestSuite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(t3TestCase, "check"))
    return suite

def main():
    unittest.main(defaultTest="makeTestSuite")


if __name__ == "__main__" : main()


