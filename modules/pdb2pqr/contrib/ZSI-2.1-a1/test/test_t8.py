#!/usr/bin/env python
import unittest, sys, types, time
from ZSI import TC, SoapWriter, ParsedSoap, EvaluateException
from ZSI.wstools.Namespaces import SCHEMA, SOAP

NSDICT = {'tns':'xmlns:tns="urn:a"',
    'xsi':'xmlns:xsi="%s"' %SCHEMA.XSI3,
    'xsd':'xmlns:xsd="%s"' %SCHEMA.XSD3,
    'soap':'xmlns:SOAP-ENC="%s"' %SOAP.ENC,
}

class AnyTestCase(unittest.TestCase):
    "Test Any serialize and parse"

    def check_empty_array(self):
        """Empty Array returned as list()
        """
        data = []
        s = str(SoapWriter().serialize(data,TC.Any(aslist=True)))
        p = ParsedSoap(s).Parse(TC.Any())
        self.failUnless(data==p, 'expecting "%s", got "%s"' %(data,p))

    def check_empty_struct(self):
        """Empty Struct is None, maybe dict() makes more sense, but this
        is fairly hard to determine if not typed (which is the norm).
        """
        data = {}
        s = str(SoapWriter().serialize(data,TC.Any()))
        p = ParsedSoap(s).Parse(TC.Any())
        self.failUnless(p==None, 'expecting "%s", got "%s"' %(None,p))

    def check_parse_empty_all(self):
        # None
        skip = [TC.FPEnumeration, TC.Enumeration, TC.IEnumeration, TC.List, TC.Integer]
        for typeclass in filter(lambda c: type(c) in [types.ClassType,type] and not issubclass(c, TC.String) and issubclass(c, TC.SimpleType), TC.__dict__.values()):
            if typeclass in skip: continue
            tc = typeclass()
            sw = SoapWriter()
            sw.serialize(None, typecode=tc, typed=True)
            soap = str(sw)
            ps = ParsedSoap(soap)
            parsed = ps.Parse(TC.Any())
            self.assertEqual(None, parsed)

    def check_parse_empty_string(self):
        # Empty String
        typecodes = TC.Any.parsemap.values()
        for tc in filter(lambda c: isinstance(c, TC.String), TC.Any.parsemap.values()):
            sw = SoapWriter()
            sw.serialize("", typecode=tc, typed=True)
            soap = str(sw)
            ps = ParsedSoap(soap)
            parsed = ps.Parse(TC.Any())
            self.assertEqual("", parsed)

    def check_builtins(self):
        myInt,myLong,myStr,myDate,myFloat = 123,2147483648,\
            u"hello", time.gmtime(), 1.0001
        orig = [myInt,myLong,myStr,myDate,myFloat]

        sw = SoapWriter()
        sw.serialize(orig, typecode=TC.Any(pname="builtins", aslist=True))
        
        ps = ParsedSoap(str(sw)) 
        parsed = ps.Parse(TC.Any())
        self.assertEqual(len(orig), len(parsed))

        self.assertEqual(myInt, parsed[0])
        self.assertEqual(myLong, parsed[1])
        self.assertEqual(myStr, parsed[2])
        self.assertEqual(myDate[0:6], parsed[3][0:6])
        self.assertEqual(myFloat, parsed[4])
        
        self.assertEqual(type(myInt), type(parsed[0]))
        self.assertEqual(type(myLong), type(parsed[1]))
        self.assertEqual(str, type(parsed[2]))
        self.assertEqual(tuple, type(parsed[3]))
        self.assertEqual(type(myFloat), type(parsed[4]))

    def check_any_nill(self):
        result = ['23', {'a' : None, 'b': 5}]
        soap = str(SoapWriter().serialize(result, TC.Any(pname="NilRequest", nillable=True, aslist=True)))

        ps = ParsedSoap(soap)
        tc = TC.Any(nillable=True)
        pyobj = ps.Parse(tc)

    def check_any_compound(self):
        # from zsi developer's guide
        xml = """
<tns:foo %(tns)s %(xsi)s %(soap)s>
    <tns:i xsi:type="SOAP-ENC:integer">12</tns:i>
    <tns:name xsi:type="SOAP-ENC:string">Hello world</tns:name>
</tns:foo>""" %NSDICT

        ps = ParsedSoap(xml, envelope=False)
        self.failUnless(ps.Parse(TC.Any()) == {'i': 12, 'name': 'Hello world'})
        self.failUnless(ps.Parse(TC.Any(aslist=True)) == [12, 'Hello world'])

    def check_any_typed_soap_integer(self):
        # from zsi developer's guide
        value = 12
        d = dict(value=value)
        d.update(NSDICT)
        xml = """<tns:i xsi:type="SOAP-ENC:integer" %(xsi)s %(soap)s %(tns)s>%(value)d</tns:i>""" %d
        ps = ParsedSoap(xml, envelope=False)
        self.failUnless(ps.Parse(TC.Any()) == value)

    def check_any_typed_xsd_int(self):
        # from zsi developer's guide
        value = 12
        d = dict(value=value)
        d.update(NSDICT)
        xml = """<tns:i xsi:type="xsd:int" %(xsi)s %(soap)s %(tns)s %(xsd)s>%(value)d</tns:i>""" %d
        ps = ParsedSoap(xml, envelope=False)
        self.failUnless(ps.Parse(TC.Any()) == value)

    def check_any_typed_nonNegativeInteger(self):
        # from zsi developer's guide
        value = 12
        d = dict(value=value)
        d.update(NSDICT)
        xml = """<tns:i xsi:type="xsd:nonNegativeInteger" %(xsi)s %(soap)s %(tns)s %(xsd)s>%(value)d</tns:i>""" %d
        ps = ParsedSoap(xml, envelope=False)
        self.failUnless(ps.Parse(TC.Any()) == value)

    def check_any_untyped_int(self):
        # from zsi developer's guide
        d = dict(value=12)
        d.update(NSDICT)
        xml = """<tns:i %(tns)s>12</tns:i>""" %NSDICT
        ps = ParsedSoap(xml, envelope=False)
        self.failUnless(int(ps.Parse(TC.Any())) == 12)

    def check_any_dict_list_rpcenc(self):
        sw = SoapWriter()
        testObj = [{"a":1,"b":2}, {"d":4,"e":5}, {"f":{"x":9}, "g":[6,7.0]}]
        typecode = TC.Any(aslist=True)
        sw.serialize(testObj, typecode=typecode)
        xml = str(sw)
        ps = ParsedSoap(xml)
        result = TC.Any().parse(ps.body_root, ps)
        self.failUnless(result == testObj)

#
# Creates permutation of test options: "check", "check_any", etc
#
_SEP = '_'
for t in [i[0].split(_SEP) for i in filter(lambda i: callable(i[1]), AnyTestCase.__dict__.items())]:
    test = ''
    for f in t:
        test += f
        if globals().has_key(test): test += _SEP; continue
        def _closure():
            name = test
            def _makeTestSuite():
                suite = unittest.TestSuite()
                suite.addTest(unittest.makeSuite(AnyTestCase, name))
                return suite
            return _makeTestSuite

        globals()[test] = _closure()
        test += _SEP

makeTestSuite = check

def main():
    unittest.main(defaultTest="makeTestSuite")

if __name__ == "__main__" : main()


