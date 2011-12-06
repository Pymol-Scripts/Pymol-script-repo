#!/usr/bin/env python
import unittest, sys, tests_good, tests_bad, time
from ZSI import *
try:
    import cStringIO as StringIO
except ImportError:
    import StringIO

class t1TestCase(unittest.TestCase):
    "Test case wrapper for old ZSI t1 test case"

    def setUp(self):
        self.goodTests = []
        self.badTests = []
        for key,val in tests_good.__dict__.items():
            try:
                if key[0:4] == "test" and int(key[4:]) > 0:
                    self.goodTests.append((key,val))
            except:
                pass
        for key,val in tests_bad.__dict__.items():
            try:
                if key[0:4] == "test" and int(key[4:]) > 0:
                    self.badTests.append((key,val))                    
            except:
                pass
        self.goodTests.sort(lambda a,b: cmp(a[0], b[0]))
        self.badTests.sort(lambda a,b: cmp(a[0], b[0]))

    def check_bug1724481(self):
        # [ 1724481 ] Error handling of enum serialization broken"
        enum = TC.Enumeration(['Red', 'Blue', 'Green'], 'color')
        sw = SoapWriter()
        self.assertRaises(EvaluateException, sw.serialize,'ddd', enum)

    def checkt1(self):
        for key,val in self.badTests:
            print "\n", "." * 60, key
            self.failUnlessRaises(ParseException, ParsedSoap, val)
        for key,val in self.goodTests:
            print "\n", "." * 60, key
            ps = ParsedSoap(val)

        ps = ParsedSoap(datatest)
        elts = ps.data_elements

        self.failUnlessEqual(TC.Integer(None, nillable=True).parse(elts[10], ps),
                                                                        None)
        self.failUnlessEqual(TC.Ibyte(None, nillable=True).parse(elts[10], ps),
                                                                        None)
        B = [ TC.Integer('Price'), TC.Integer('p2'), TC.String(unique=1) ]
        self.failUnlessEqual(TC.Integer(('test-uri', 'Price')).parse(elts[0], ps), 
                                                        34)
        self.failUnlessEqual(B[0].parse(elts[0], ps), 34)
        self.failUnlessEqual(B[1].parse(elts[1], ps), 44)
        self.failUnlessEqual(B[2].parse(elts[2], ps), u"This is the name")
        self.failUnlessEqual(TC.HexBinaryString().parse(elts[9], ps), "? A")
        self.failUnlessEqual(TC.String('Name').parse(elts[2], ps), 
                                                    u"This is the name")
        self.failUnlessEqual(TC.Any('Price').parse(elts[0], ps), 34)
        self.failUnlessEqual(TC.Any('n3').parse(elts[4], ps), 
                                                    u"The value of n3")
        TC.XML('n2').parse(elts[3], ps)
        nodelist = TC.XML('a2').parse(elts[7], ps)
        self.failUnlessEqual(TC.String('n3').parse(elts[4], ps), 
                                                    u"The value of n3")
        self.failUnlessEqual(TC.Base64String('n64').parse(elts[5], ps),
                                                    u"hello")
        self.failUnlessEqual(TC.String('n64').parse(elts[5], ps),
                                                    u"a GVsbG8=")
        enum = TC.Enumeration(['Red', 'Blue', 'Green'], 'color')
        self.failUnlessEqual(enum.parse(elts[6], ps), u'Red')
        self.failUnlessEqual(TC.IEnumeration([44,46,47]).parse(elts[1],ps),
                                                        44)
        S = TC.Struct(None, [TC.String('t'), TC.Integer('i')], inorder=0)
        pyobj = S.parse(elts[8], ps)
        S2 = TC.Struct(myclass, [TC.IunsignedShort('i'), TC.String('q:z',
        minOccurs=0), TC.String('t')], 'a2', typed=0)
        pyobj2 = S2.parse(elts[8], ps)
        self.failUnlessEqual(TC.URI().parse(elts[12], ps), 
                                            u'"http://foo.com/~salz"')
        self.failUnlessEqual(pyobj["i"], pyobj2.i)
        self.failUnlessEqual(pyobj["t"], pyobj2.t)

        tcary = TC.Array('SOAP-ENC:int', TC.Integer())
        nsa = tcary.parse(elts[14],ps)
        self.failUnlessEqual(nsa, [None, None, None, 12, 13, 14, 15, 16, 17])
        tcary.sparse = 1
        sa = tcary.parse(elts[14],ps)
        self.failUnlessEqual(sa, 
                    [(3, 12), (4, 13), (5, 14), (6, 15), (7, 16), (8, 17)])

        """
        mychoice = TC.Choice([
            TC.String('n3'),
            TC.URI('uri'),
            TC.Integer('Price'),
        ])

        b = mychoice.parse(elts[0], ps)
        self.failUnlessEqual(b[0], 'Price')
        self.failUnlessEqual(b[1], 34)
        b = mychoice.parse(elts[12], ps)
        self.failUnlessEqual(b[0], 'uri')
        self.failUnlessEqual(b[1], u'"http://foo.com/~salz"')
        b = mychoice.parse(elts[4], ps)
        self.failUnlessEqual(b[0], 'n3')
        self.failUnlessEqual(b[1], u'The value of n3')
        """

        self.failUnlessEqual(TC.Array(('test-uri','x'), TC.Any()).parse(elts[15], ps),
                                            [u'The value of n3', u'rich salz', 13])
        self.failUnlessEqual(TC.Struct(None,(TC.FPfloat('a'), TC.Decimal('b'),
                                            TC.FPdouble('c'))).parse(elts[13],ps),
                                            {'a': 6.9000000000000004, 'c':
                                                TC._make_inf(), 'b': 0.0})
        nsdict = ps.GetElementNSdict(ps.header)
        nsdict[''] = "http://www.zolera.com/ns/"
        nsdict['q'] = 'q-namespace-uri' 
        sio = StringIO.StringIO()
        z = SoapWriter(sio, header=ps.header_elements, nsdict=nsdict) 
        z.serialize(pyobj2, S2) 
        S2.inline = 1 
        S2.typed = 0 
        tc = TC.gDateTime('dt') 
        z.serialize(pyobj2, S2) 
        z.serialize(pyobj, S) 
        #z.serialize(('n3', '******this is the value of a union'), mychoice) 
        z.serialize('uri:some/place/special', TC.XML('foo', nsdict=nsdict)) 
        tcary.sparse = False
        z.serialize(nsa, tcary, childnames='tt') 
        tcary.sparse = True
        z.serialize(sa, tcary, name='MYSPARSEARRAY') 
        z.serialize(time.time(), tc) 
        z.serialize(time.time(), TC.gTime('monthday')) 
        z.serialize('$$$$$foo<', TC.String(textprotect=0)) 
        self.failUnlessEqual(TC.Any().parse(elts[11], ps),
                                        {'urt-i': 12, 'urt-t': u'rich salz'})

        try: 
            a = bar() 
        except Exception, e: 
            f = FaultFromException(e, 0, sys.exc_info()[2]) 
            print f.AsSOAP() 
        print
        print
        print FaultFromNotUnderstood('myuri', 'dalocalname', actor='cher').AsSOAP() 
        print
        print
        print FaultFromActor('actor:i:dont:understand').AsSOAP()


def makeTestSuite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(t1TestCase, "check"))
    return suite

##  exceptions
def foo():
    '''dummy'''
    return 3 / 0

def bar():
    return foo() + 2

class zParseException: pass

class myclass:
    def __init__(self, name=None):
        self.name = name or id(self)
        self.z = 'z value'
    def __str__(self):
        return 'myclass-%s-(%d,"%s")' % (self.name, self.i, self.t) + \
                                                        str(self.z)

datatest = '''<SOAP-ENV:Envelope
  xmlns:SOAP-ENV="http://schemas.xmlsoap.org/soap/envelope/"
  xmlns:SOAP-ENC="http://schemas.xmlsoap.org/soap/encoding/"
  xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
  xmlns:xsd="http://www.w3.org/2001/XMLSchema"
  SOAP-ENV:encodingStyle="http://schemas.xmlsoap.org/soap/encoding/">
  <SOAP-ENV:Header xmlns:t="http://www.zolera.com/ns/" xmlns:q='"'>
  <t:sometag SOAP-ENV:mustUnderstand="1">you must grok sometag</t:sometag>
  </SOAP-ENV:Header>
    <SOAP-ENV:Body xmlns='test-uri'>
        <root SOAP-ENC:root='1'/>
        <Price xsi:type='xsd:integer'>34</Price> <!-- 0 -->
        <SOAP-ENC:byte>44</SOAP-ENC:byte> <!-- 1 -->
        <Name>This is the name</Name>   <!-- 2 -->
        <n2><xmldoc><![CDATA[<greeting>Hello</greeting>]]></xmldoc></n2> <!-- 3 -->
        <n3 href='#zzz' xsi:type='SOAP-ENC:string'/> <!-- 4 -->
        <n64>a GVsbG8=</n64> <!-- 5 -->
        <SOAP-ENC:string>Red</SOAP-ENC:string> <!-- 6 -->
        <a2 href='#tri2'/> <!-- 7 --> 
        <a2><i>12</i><t>rich salz</t></a2> <!-- 8 --> 
        <xsd:hexBinary>3F2041</xsd:hexBinary> <!-- 9 --> 
        <nullint xsi:nil='1'/> <!-- 10 --> 
        <Anytest><urt-i xsi:type='SOAP-ENC:byte'>12</urt-i> 
        <urt-t id="urtid" 
            xsi:type="xsd:string">rich salz</urt-t></Anytest> <!-- 11 --> 
        <uri>"http://foo.com/%7Esalz"</uri> <!-- 12 --> 
        <floattest> <!-- 13 --> 
            <a>6.9</a> <b>-0</b> <c>INF</c> 
        </floattest> 
        <atest SOAP-ENC:offset='[3]' SOAP-ENC:arrayType="x"> <!-- 14 --> 
            <i>12</i> 
            <SOAP-ENC:integer id='n13'>13</SOAP-ENC:integer> 
            <i>14</i> 
            <i>15</i> 
            <i>16</i> 
            <i>17</i> 
        </atest> 
        <sarray SOAP-ENC:arrayType="struct"> <!-- 15 --> 
            <i href="#zzz" xsi:type='xsd:string'/> 
            <i href="#urtid"/> 
            <thing href="#n13"/> 
        </sarray> 
        <xpath>//sarray</xpath> <!-- 16 --> 
        <z xmlns='myns' xsi:type='SOAP-ENC:string' id='zzz'>The value of n3</z> 
        <zz xmlns='myns2' id='tri2'>
            <inner xmlns='myns2' >
                <f1>content</f1>
                <sec xmlns='myns2' >ond</sec >
            </inner>
        </zz> 
    </SOAP-ENV:Body> 
</SOAP-ENV:Envelope>'''

def main():
    unittest.main(defaultTest="makeTestSuite")


if __name__ == "__main__" : main()


