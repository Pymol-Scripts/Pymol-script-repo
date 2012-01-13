#!/usr/bin/env python
import unittest, sys, tests_good, tests_bad, time
from ZSI import *
try:
    import cStringIO as StringIO
except ImportError:
    import StringIO


"""Bug [ 1520092 ] URI Bug: urllib.quote escaping reserved chars

From rfc2396:

"If the data for a URI component would conflict with the reserved 
purpose, then the conflicting data must be escaped before forming the 
URI."

reserved = ";" | "/" | "?" | ":" | "@" | "&" | "=" | "+" |
"$" | ","



This implies that if ":" is used for a reserved purpose,

if scheme is defined then
append scheme to result
append ":" to result

, then it should not be escaped.
"""


class TestCase(unittest.TestCase):
    def check_uri_quoting(self):
        """ all reserved characters used for reserved purpose.
        """ 
        sw1 = SoapWriter(envelope=False)
        tc1= TC.URI('sourceforge')
        orig = 'https://sourceforge.net/tracker/index.php?func=detail&aid=1520092&group_id=26590&atid=387667'
        sw1.serialize(orig, typecode=tc1, typed=False)
        s1 = str(sw1)

        sw2 = SoapWriter(envelope=False)
        tc2= TC.String('sourceforge')
        sw2.serialize(orig, typecode=tc2, typed=False)
        s2 = str(sw2)

        print s1
        print s2
        self.failUnless(s1 == s2,
            'reserved characters used for reserved purpose should not be escaped.')

        ps = ParsedSoap(s2, envelope=False)
        pyobj = ps.Parse(tc2)

        self.failUnless(pyobj == orig, 'parsed object should be equivalent to original')



#
# Creates permutation of test options: "check", "check_any", etc
#
_SEP = '_'
for t in [i[0].split(_SEP) for i in filter(lambda i: callable(i[1]), TestCase.__dict__.items())]:
    test = ''
    for f in t:
        test += f
        if globals().has_key(test): test += _SEP; continue
        def _closure():
            name = test
            def _makeTestSuite():
                suite = unittest.TestSuite()
                suite.addTest(unittest.makeSuite(TestCase, name))
                return suite
            return _makeTestSuite

        globals()[test] = _closure()
        test += _SEP


makeTestSuite = check
def main():
    unittest.main(defaultTest="makeTestSuite")
if __name__ == "__main__" : main()


