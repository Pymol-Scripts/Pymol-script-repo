#!/usr/bin/env python
import unittest, sys, sha, base64
from ZSI import _get_element_nsuri_name
from ZSI.parse import ParsedSoap
from ZSI.wstools.c14n import Canonicalize
from ZSI.wstools.Namespaces import WSA200403, SOAP
from cStringIO import StringIO


class CanonicalizeFromTestCase(unittest.TestCase):
    "c14n tests, this has nothing to do with ws-addressing."

    def setUp(self):
        self.ps = ParsedSoap(XML_INST1)
        self.el = filter(lambda el: _get_element_nsuri_name(el) == (WSA200403.ADDRESS, "From"), 
                      self.ps.header_elements)[0]

    def tearDown(self):
        del self.ps
        del self.el

    def check_c14n(self):
        """http://www.w3.org/TR/xml-c14n
        """
        s = StringIO()
        Canonicalize(self.el, s, unsuppressedPrefixes=None)
        cxml = s.getvalue()
        d1 = base64.encodestring(sha.sha(C14N_INC1).digest()).strip()
        d2 = base64.encodestring(sha.sha(cxml).digest()).strip()
        self.assertEqual(d1, d2)
        self.assertEqual(d1, C14N_INC1_DIGEST)

    def check_c14n_exc(self):
        """http://www.w3.org/TR/xml-exc-c14n/
        """
        s = StringIO()
        Canonicalize(self.el, s, unsuppressedPrefixes=[])
        cxml = s.getvalue()
        d1 = base64.encodestring(sha.sha(C14N_EXCL1).digest()).strip()
        d2 = base64.encodestring(sha.sha(cxml).digest()).strip()
        self.assertEqual(d1, C14N_EXCL1_DIGEST)
        self.assertEqual(d1, d2)

    def check_c14n_exc2_unsuppressed(self):
        """http://www.w3.org/TR/xml-exc-c14n/
        The method of canonicalization described in this specification receives 
        an InclusiveNamespaces PrefixList parameter, which lists namespace prefixes 
        that are handled in the manner described by the Canonical XML Recommendation 
        """
        s = StringIO()
        Canonicalize(self.el, s, unsuppressedPrefixes=['xsi', 'xsd'])
        cxml = s.getvalue()
        d1 = base64.encodestring(sha.sha(C14N_EXCL2).digest()).strip()
        d2 = base64.encodestring(sha.sha(cxml).digest()).strip()
        self.assertEqual(d1, C14N_EXCL2_DIGEST)
        self.assertEqual(d1, d2)

    def check_c14n_exc3(self):
        """http://www.w3.org/TR/xml-exc-c14n/
        tests if a namespace defined in a parent node to the top node 
        to be canonicalized is added when discovered that this namespace
        is used.
        """
        self.ps = ParsedSoap(XML_INST2)
        self.el = self.ps.body

        s = StringIO()
        Canonicalize(self.el, s, unsuppressedPrefixes=[])
        cxml = s.getvalue()
        print cxml
        d1 = base64.encodestring(sha.sha(C14N_EXCL3).digest()).strip()
        d2 = base64.encodestring(sha.sha(cxml).digest()).strip()
        self.assertEqual(d1, C14N_EXCL3_DIGEST)
        self.assertEqual(d1, d2)

    def xcheck_c14n_exc4(self):
        RCVDIGEST = "jhTbi7gWlY9oLqsRr+EZ0bokRFA="
        CALDIGEST = "IkMyI4zCDlK41qE7sZxvkFHJioU="

        d1 = base64.encodestring(sha.sha(WRONG).digest()).strip()
        d2 = base64.encodestring(sha.sha(CORRECT).digest()).strip()

        ps = ParsedSoap(XML_INST4)
        el = filter(lambda el: _get_element_nsuri_name(el) == (WSA200403.ADDRESS, "MessageID"), 
                      ps.header_elements)[0]

        s = StringIO()
        Canonicalize(el, s, unsuppressedPrefixes=[])
        cxml = s.getvalue()
        print "-- "*20
        print cxml
        print "-- "*20
        d3 = base64.encodestring(sha.sha(cxml).digest()).strip()

        self.assertEqual(d1, RCVDIGEST)
        self.assertEqual(d2, CALDIGEST)
        self.assertEqual(d3, CALDIGEST)


def makeTestSuite():
    suite = unittest.TestSuite()
    #suite.addTest(unittest.makeSuite(CanonicalizeFromTestCase, "check"))
    suite.addTest(unittest.makeSuite(CanonicalizeFromTestCase, "xcheck"))
    return suite


C14N_EXCL1_DIGEST =  "xSOXT+dlQwo5uT9PbK08of6W9PM="
C14N_EXCL1 = """<wsa:From xmlns:ns3="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-utility-1.0.xsd" xmlns:soapenv="http://schemas.xmlsoap.org/soap/envelope/" xmlns:wsa="http://schemas.xmlsoap.org/ws/2004/03/addressing" ns3:Id="id-7680063" soapenv:mustUnderstand="0"><wsa:Address>http://bosshog.lbl.gov:9999/wsrf/services/SecureCounterService</wsa:Address><wsa:ReferenceProperties><ns1:CounterKey xmlns:ns1="http://counter.com" ns3:Id="10112">10577413</ns1:CounterKey></wsa:ReferenceProperties></wsa:From>"""

C14N_INC1_DIGEST =  "qdU4f7/+BeHV/JlVGIPM90fNeV8="
C14N_INC1 = """<wsa:From xmlns:ns1="http://counter.com" xmlns:ns3="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-utility-1.0.xsd" xmlns:soapenv="http://schemas.xmlsoap.org/soap/envelope/" xmlns:wsa="http://schemas.xmlsoap.org/ws/2004/03/addressing" xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" ns3:Id="id-7680063" soapenv:mustUnderstand="0"><wsa:Address>http://bosshog.lbl.gov:9999/wsrf/services/SecureCounterService</wsa:Address><wsa:ReferenceProperties><ns1:CounterKey ns3:Id="10112">10577413</ns1:CounterKey></wsa:ReferenceProperties></wsa:From>"""

C14N_EXCL2_DIGEST = "+IEqF6DRo36Bh93A06S7C4Cmcuo="
C14N_EXCL2 = """<wsa:From xmlns:ns3="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-utility-1.0.xsd" xmlns:soapenv="http://schemas.xmlsoap.org/soap/envelope/" xmlns:wsa="http://schemas.xmlsoap.org/ws/2004/03/addressing" xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" ns3:Id="id-7680063" soapenv:mustUnderstand="0"><wsa:Address>http://bosshog.lbl.gov:9999/wsrf/services/SecureCounterService</wsa:Address><wsa:ReferenceProperties><ns1:CounterKey xmlns:ns1="http://counter.com" ns3:Id="10112">10577413</ns1:CounterKey></wsa:ReferenceProperties></wsa:From>"""


C14N_EXCL3_DIGEST = "VJvTr+Mx3TeWsQY6iwGbhAJ9/eA="
C14N_EXCL3 = """<soapenv:Body xmlns:soapenv="http://schemas.xmlsoap.org/soap/envelope/" xmlns:wsu="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-utility-1.0.xsd" wsu:Id="id-28219008"><RequestSecurityTokenResponse xmlns="http://schemas.xmlsoap.org/ws/2004/04/trust"><wsa:EndpointReference xmlns:wsa="http://schemas.xmlsoap.org/ws/2004/03/addressing"><wsa:Address>http://131.243.2.147:8888/wsrf/services/DelegationService</wsa:Address><wsa:ReferenceProperties><ns1:DelegationKey xmlns:ns1="http://www.globus.org/08/2004/delegationService">8adaa710-ba01-11da-bc99-cbed73daa755</ns1:DelegationKey></wsa:ReferenceProperties></wsa:EndpointReference></RequestSecurityTokenResponse></soapenv:Body>"""


XML_INST1 = """<soapenv:Envelope xmlns:soapenv="http://schemas.xmlsoap.org/soap/envelope/" xmlns:wsa="http://schemas.xmlsoap.org/ws/2004/03/addressing" xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"><soapenv:Header>
	<wsse:Security soapenv:mustUnderstand="1" xmlns:wsse="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-secext-1.0.xsd"><ds:Signature xmlns:ds="http://www.w3.org/2000/09/xmldsig#">
	<ds:SignedInfo>
	<ds:CanonicalizationMethod Algorithm="http://www.w3.org/2001/10/xml-exc-c14n#"><ec:InclusiveNamespaces PrefixList="soapenv wsa xsd xsi" xmlns:ec="http://www.w3.org/2001/10/xml-exc-c14n#"/></ds:CanonicalizationMethod>
	<ds:SignatureMethod Algorithm="http://www.globus.org/2002/04/xmlenc#gssapi-sign"/>
	<ds:Reference URI="#id-8409752">
	<ds:Transforms>
	<ds:Transform Algorithm="http://www.w3.org/2001/10/xml-exc-c14n#"><ec:InclusiveNamespaces PrefixList="wsa xsd xsi" xmlns:ec="http://www.w3.org/2001/10/xml-exc-c14n#"/></ds:Transform>
	</ds:Transforms>
	<ds:DigestMethod Algorithm="http://www.w3.org/2000/09/xmldsig#sha1"/>
	<ds:DigestValue>m9pihAqIBdcdk7ytDvccj89eWi8=</ds:DigestValue>
	</ds:Reference>
	<ds:Reference URI="#id-11434871">
	<ds:Transforms>
	<ds:Transform Algorithm="http://www.w3.org/2001/10/xml-exc-c14n#"><ec:InclusiveNamespaces PrefixList="xsd xsi" xmlns:ec="http://www.w3.org/2001/10/xml-exc-c14n#"/></ds:Transform>
	</ds:Transforms>
	<ds:DigestMethod Algorithm="http://www.w3.org/2000/09/xmldsig#sha1"/>
	<ds:DigestValue>ofD+Ket5kzR2u/5jWKbFTMtmigk=</ds:DigestValue>
	</ds:Reference>
	<ds:Reference URI="#id-19645447">
	<ds:Transforms>
	<ds:Transform Algorithm="http://www.w3.org/2001/10/xml-exc-c14n#"><ec:InclusiveNamespaces PrefixList="xsd xsi" xmlns:ec="http://www.w3.org/2001/10/xml-exc-c14n#"/></ds:Transform>
	</ds:Transforms>
	<ds:DigestMethod Algorithm="http://www.w3.org/2000/09/xmldsig#sha1"/>
	<ds:DigestValue>SoQ7RlJa3r94weDWBuWAg/BvydQ=</ds:DigestValue>
	</ds:Reference>
	<ds:Reference URI="#id-5428820">
	<ds:Transforms>
	<ds:Transform Algorithm="http://www.w3.org/2001/10/xml-exc-c14n#"><ec:InclusiveNamespaces PrefixList="xsd xsi" xmlns:ec="http://www.w3.org/2001/10/xml-exc-c14n#"/></ds:Transform>
	</ds:Transforms>
	<ds:DigestMethod Algorithm="http://www.w3.org/2000/09/xmldsig#sha1"/>
	<ds:DigestValue>z6sCEkkRJrCuY/C0S5b+46WfyMs=</ds:DigestValue>
	</ds:Reference>
	<ds:Reference URI="#id-7680063">
	<ds:Transforms>
	<ds:Transform Algorithm="http://www.w3.org/2001/10/xml-exc-c14n#"><ec:InclusiveNamespaces PrefixList="xsd xsi" xmlns:ec="http://www.w3.org/2001/10/xml-exc-c14n#"/></ds:Transform>
	</ds:Transforms>
	<ds:DigestMethod Algorithm="http://www.w3.org/2000/09/xmldsig#sha1"/>
	<ds:DigestValue>+IEqF6DRo36Bh93A06S7C4Cmcuo=</ds:DigestValue>
	</ds:Reference>
	<ds:Reference URI="#id-28476580">
	<ds:Transforms>
	<ds:Transform Algorithm="http://www.w3.org/2001/10/xml-exc-c14n#"><ec:InclusiveNamespaces PrefixList="xsd xsi" xmlns:ec="http://www.w3.org/2001/10/xml-exc-c14n#"/></ds:Transform>
	</ds:Transforms>
	<ds:DigestMethod Algorithm="http://www.w3.org/2000/09/xmldsig#sha1"/>
	<ds:DigestValue>NFltkKAJpmMkPbJQj5MW1qVceto=</ds:DigestValue>
	</ds:Reference>
	</ds:SignedInfo>
	<ds:SignatureValue>AAAAAAAAAAMAAAvZTrXlZjRSO7tP12tId+lehprEKgk=</ds:SignatureValue>
	<ds:KeyInfo>
	<wsse:SecurityTokenReference><wsse:Reference URI="#SecurityContextToken-32970611" ValueType="http://www.globus.org/ws/2004/09/security/sc#GSSAPI_CONTEXT_TOKEN"/></wsse:SecurityTokenReference>
	</ds:KeyInfo>
	</ds:Signature><wsc:SecurityContextToken wsu:Id="SecurityContextToken-32970611" xmlns:wsc="http://schemas.xmlsoap.org/ws/2004/04/sc" xmlns:wsu="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-utility-1.0.xsd"><wsc:Identifier>3b1ef410-ab3d-11da-9436-88b687faed94</wsc:Identifier></wsc:SecurityContextToken></wsse:Security><wsa:MessageID wsu:Id="id-11434871" soapenv:mustUnderstand="0" xmlns:wsu="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-utility-1.0.xsd">uuid:3d592ca0-ab3d-11da-9436-88b687faed94</wsa:MessageID><wsa:To wsu:Id="id-19645447" soapenv:mustUnderstand="0" xmlns:wsu="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-utility-1.0.xsd">http://schemas.xmlsoap.org/ws/2004/03/addressing/role/anonymous</wsa:To><wsa:Action wsu:Id="id-5428820" soapenv:mustUnderstand="0" xmlns:wsu="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-utility-1.0.xsd">http://counter.com/CounterPortType/addResponse</wsa:Action><wsa:From ns3:Id="id-7680063" soapenv:mustUnderstand="0" xmlns:ns1="http://counter.com" xmlns:ns3="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-utility-1.0.xsd"><wsa:Address>http://bosshog.lbl.gov:9999/wsrf/services/SecureCounterService</wsa:Address><wsa:ReferenceProperties><ns1:CounterKey ns3:Id="10112">10577413</ns1:CounterKey></wsa:ReferenceProperties></wsa:From><wsa:RelatesTo RelationshipType="wsa:Reply" wsu:Id="id-28476580" soapenv:mustUnderstand="0" xmlns:wsu="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-utility-1.0.xsd">uuid:1141449047.05</wsa:RelatesTo></soapenv:Header><soapenv:Body wsu:Id="id-8409752" xmlns:wsu="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-utility-1.0.xsd"><addResponse xmlns="http://counter.com">13</addResponse></soapenv:Body></soapenv:Envelope>"""


XML_INST2 = """<?xml version="1.0" encoding="UTF-8"?>
        <soapenv:Envelope 
xmlns:soapenv="http://schemas.xmlsoap.org/soap/envelope/" 
xmlns:wsa="http://schemas.xmlsoap.org/ws/2004/03/addressing" 
xmlns:xsd="http://www.w3.org/2001/XMLSchema" 
xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"><soapenv:Header>
        <wsse:Security soapenv:mustUnderstand="1" 
xmlns:wsse="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-secext-1.0.xsd"><wsse:BinarySecurityToken 
EncodingType="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-soap-message-security-1.0#Base64Binary" 
ValueType="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-x509-token-profile-1.0#X509PKIPathv1" 
wsu:Id="CertId-1851922" 
xmlns:wsu="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-utility-1.0.xsd">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</wsse:BinarySecurityToken><ds:Signature 
xmlns:ds="http://www.w3.org/2000/09/xmldsig#">
        <ds:SignedInfo>
        <ds:CanonicalizationMethod 
Algorithm="http://www.w3.org/2001/10/xml-exc-c14n#"/>
        <ds:SignatureMethod 
Algorithm="http://www.w3.org/2000/09/xmldsig#rsa-sha1"/>
        <ds:Reference URI="#id-25484440">
        <ds:Transforms>
        <ds:Transform Algorithm="http://www.w3.org/2001/10/xml-exc-c14n#"/>
        </ds:Transforms>
        <ds:DigestMethod 
Algorithm="http://www.w3.org/2000/09/xmldsig#sha1"/>
        <ds:DigestValue>xqoPUGjk97yY+StAheOFmeaHgbw=</ds:DigestValue>
        </ds:Reference>
        <ds:Reference URI="#id-28219008">
        <ds:Transforms>
        <ds:Transform Algorithm="http://www.w3.org/2001/10/xml-exc-c14n#"/>
        </ds:Transforms>
        <ds:DigestMethod 
Algorithm="http://www.w3.org/2000/09/xmldsig#sha1"/>
        <ds:DigestValue>VJvTr+Mx3TeWsQY6iwGbhAJ9/eA=</ds:DigestValue>
        </ds:Reference>
        <ds:Reference URI="#id-18539969">
        <ds:Transforms>
        <ds:Transform Algorithm="http://www.w3.org/2001/10/xml-exc-c14n#"/>
        </ds:Transforms>
        <ds:DigestMethod 
Algorithm="http://www.w3.org/2000/09/xmldsig#sha1"/>
        <ds:DigestValue>W1PrEK32GMCbF6FTEmlYiYwqAeQ=</ds:DigestValue>
        </ds:Reference>
        <ds:Reference URI="#id-14816181">
        <ds:Transforms>
        <ds:Transform Algorithm="http://www.w3.org/2001/10/xml-exc-c14n#"/>
        </ds:Transforms>
        <ds:DigestMethod 
Algorithm="http://www.w3.org/2000/09/xmldsig#sha1"/>
        <ds:DigestValue>lWhqYlKqBnB4LwkRWyXMwHy18hc=</ds:DigestValue>
        </ds:Reference>
        <ds:Reference URI="#id-8120088">
        <ds:Transforms>
        <ds:Transform Algorithm="http://www.w3.org/2001/10/xml-exc-c14n#"/>
        </ds:Transforms>
        <ds:DigestMethod 
Algorithm="http://www.w3.org/2000/09/xmldsig#sha1"/>
        <ds:DigestValue>2Zjgz4McHaxMLfpBbqelAqWvRsU=</ds:DigestValue>
        </ds:Reference>
        <ds:Reference URI="#id-8450175">
        <ds:Transforms>
        <ds:Transform Algorithm="http://www.w3.org/2001/10/xml-exc-c14n#"/>
        </ds:Transforms>
        <ds:DigestMethod 
Algorithm="http://www.w3.org/2000/09/xmldsig#sha1"/>
        <ds:DigestValue>knsi7QmfOvjrn5mWClmsbCpZ32A=</ds:DigestValue>
        </ds:Reference>
        <ds:Reference URI="#id-19744521">
        <ds:Transforms>
        <ds:Transform Algorithm="http://www.w3.org/2001/10/xml-exc-c14n#"/>
        </ds:Transforms>
        <ds:DigestMethod 
Algorithm="http://www.w3.org/2000/09/xmldsig#sha1"/>
        <ds:DigestValue>/yBXQ0yxPqpgSGYym/DA08k0dXM=</ds:DigestValue>
        </ds:Reference>
        </ds:SignedInfo>
        <ds:SignatureValue>
        
LuGNfoVBzUIoF0AU0lzJkH9kAOi+PQVG8hMrCIEjWh1lifSG/bquhu/qZVq78x3UR+tGK411hWuQ
        nGle1GvY0A==
        </ds:SignatureValue>
        <ds:KeyInfo Id="KeyId-17834932">
        <wsse:SecurityTokenReference wsu:Id="STRId-9973812" 
xmlns:wsu="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-utility-1.0.xsd"><wsse:Reference 
URI="#CertId-1851922" 
ValueType="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-x509-token-profile-1.0#X509PKIPathv1"/></wsse:SecurityTokenReference>
        </ds:KeyInfo>
        </ds:Signature>
        <wsu:Timestamp wsu:Id="id-25484440" 
xmlns:wsu="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-utility-1.0.xsd"><wsu:Created>2006-03-23T00:11:14Z</wsu:Created><wsu:Expires>2006-03-23T00:16:14Z</wsu:Expires></wsu:Timestamp></wsse:Security><wsa:MessageID 
wsu:Id="id-8120088" soapenv:mustUnderstand="0" 
xmlns:wsu="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-utility-1.0.xsd">uuid:8aec8160-ba01-11da-bc99-cbed73daa755</wsa:MessageID><wsa:To 
wsu:Id="id-18539969" soapenv:mustUnderstand="0" 
xmlns:wsu="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-utility-1.0.xsd">http://schemas.xmlsoap.org/ws/2004/03/addressing/role/anonymous</wsa:To><wsa:Action 
wsu:Id="id-19744521" soapenv:mustUnderstand="0" 
xmlns:wsu="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-utility-1.0.xsd">http://www.globus.org/08/2004/delegationService/DelegationFactoryPortType/RequestSecurityTokenResponse</wsa:Action><wsa:From 
wsu:Id="id-14816181" soapenv:mustUnderstand="0" 
xmlns:wsu="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-utility-1.0.xsd"><wsa:Address>http://bosshog.lbl.gov:8888/wsrf/services/DelegationFactoryService</wsa:Address></wsa:From><wsa:RelatesTo 
RelationshipType="wsa:Reply" wsu:Id="id-8450175" 
soapenv:mustUnderstand="0" 
xmlns:wsu="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-utility-1.0.xsd">uuid:1143072675.25</wsa:RelatesTo></soapenv:Header><soapenv:Body 
wsu:Id="id-28219008" 
xmlns:wsu="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-utility-1.0.xsd"><RequestSecurityTokenResponse 
xmlns="http://schemas.xmlsoap.org/ws/2004/04/trust"><wsa:EndpointReference 
xmlns:ns1="http://www.globus.org/08/2004/delegationService"><wsa:Address>http://131.243.2.147:8888/wsrf/services/DelegationService</wsa:Address><wsa:ReferenceProperties><ns1:DelegationKey>8adaa710-ba01-11da-bc99-cbed73daa755</ns1:DelegationKey></wsa:ReferenceProperties></wsa:EndpointReference></RequestSecurityTokenResponse></soapenv:Body></soapenv:Envelope>"""



CORRECT = """<ns1:MessageID xmlns:ns1="http://schemas.xmlsoap.org/ws/2004/03/addressing" xmlns:ns2="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-utility-1.0.xsd" ns2:Id="10102">uuid:1143760705.98</ns1:MessageID>"""

WRONG = """<ns1:MessageID xmlns:ns1="http://schemas.xmlsoap.org/ws/2004/03/addressing" xmlns:ns2="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-secext-1.0.xsd" ns2:Id="10102">uuid:1143760705.98</ns1:MessageID>"""

XML_INST4 = """<?xml version="1.0" encoding="UTF-8"?>
<SOAP-ENV:Envelope xmlns:SOAP-ENV="http://schemas.xmlsoap.org/soap/envelope/" xmlns:SOAP-ENC="http://schemas.xmlsoap.org/soap/encoding/" xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ZSI="http://www.zolera.com/schemas/ZSI/"><SOAP-ENV:Header xmlns:ns2="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-secext-1.0.xsd" xmlns:ns1="http://schemas.xmlsoap.org/ws/2004/03/addressing"><ns1:MessageID xmlns:ns2="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-utility-1.0.xsd" ns2:Id="10102">uuid:1143760705.98</ns1:MessageID><ns1:Action xmlns:ns2="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-utility-1.0.xsd" ns2:Id="10103">http://counter.com/CounterPortType/createCounterRequest</ns1:Action><ns1:To xmlns:ns2="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-utility-1.0.xsd" ns2:Id="10104">http://131.243.2.159:9080/wsrf/services/SecureCounterService</ns1:To><ns1:From xmlns:ns2="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-utility-1.0.xsd" ns2:Id="10105"><ns1:Address>http://schemas.xmlsoap.org/ws/2004/03/addressing/role/anonymous</ns1:Address></ns1:From><ns2:Security xmlns:ns3="http://www.w3.org/2000/09/xmldsig#" xmlns:ns4="http://schemas.xmlsoap.org/ws/2004/04/sc"><ns3:Signature><ns3:SignedInfo xsi:type="ns3:SignedInfoType"><ns3:CanonicalizationMethod xsi:type="ns3:CanonicalizationMethodType" Algorithm="http://www.w3.org/2001/10/xml-exc-c14n#"></ns3:CanonicalizationMethod><ns3:SignatureMethod xsi:type="ns3:SignatureMethodType" Algorithm="http://www.globus.org/2002/04/xmlenc#gssapi-sign"></ns3:SignatureMethod><ns3:Reference xsi:type="ns3:ReferenceType" URI="#10102"><ns3:DigestMethod xsi:type="ns3:DigestMethodType" Algorithm="http://www.w3.org/2000/09/xmldsig#sha1"></ns3:DigestMethod><ns3:DigestValue xsi:type="ns3:DigestValueType">
        IkMyI4zCDlK41qE7sZxvkFHJioU=
        </ns3:DigestValue></ns3:Reference><ns3:Reference xsi:type="ns3:ReferenceType" URI="#10103"><ns3:DigestMethod xsi:type="ns3:DigestMethodType" Algorithm="http://www.w3.org/2000/09/xmldsig#sha1"></ns3:DigestMethod><ns3:DigestValue xsi:type="ns3:DigestValueType">
        DyEF6Pa7w3SSEVJ98LIoX2LW85k=
        </ns3:DigestValue></ns3:Reference><ns3:Reference xsi:type="ns3:ReferenceType" URI="#10104"><ns3:DigestMethod xsi:type="ns3:DigestMethodType" Algorithm="http://www.w3.org/2000/09/xmldsig#sha1"></ns3:DigestMethod><ns3:DigestValue xsi:type="ns3:DigestValueType">
        p/2PhmYP+/1UPcpwsRcdlvLmOAg=
        </ns3:DigestValue></ns3:Reference><ns3:Reference xsi:type="ns3:ReferenceType" URI="#10105"><ns3:DigestMethod xsi:type="ns3:DigestMethodType" Algorithm="http://www.w3.org/2000/09/xmldsig#sha1"></ns3:DigestMethod><ns3:DigestValue xsi:type="ns3:DigestValueType">
        KFLeYjf5ohGUIoPoZV/oew9SuUM=
        </ns3:DigestValue></ns3:Reference><ns3:Reference xsi:type="ns3:ReferenceType" URI="#10106"><ns3:DigestMethod xsi:type="ns3:DigestMethodType" Algorithm="http://www.w3.org/2000/09/xmldsig#sha1"></ns3:DigestMethod><ns3:DigestValue xsi:type="ns3:DigestValueType">
        7Gg0SC1wltHVAwiOfdgZsGM9W5g=
        </ns3:DigestValue></ns3:Reference></ns3:SignedInfo><ns3:SignatureValue xsi:type="ns3:SignatureValueType">
        AAAAAAAAAAEAAAdrBxzrHLZG4NglRglL9F3rKQu0658=
        </ns3:SignatureValue><ns3:KeyInfo xsi:type="ns3:KeyInfoType"><ns2:SecurityTokenReference><ns2:Reference URI="#CertId-10107" ValueType="http://www.globus.org/ws/2004/09/security/sc#GSSAPI_CONTEXT_TOKEN"></ns2:Reference></ns2:SecurityTokenReference></ns3:KeyInfo></ns3:Signature><ns4:SecurityContextToken xmlns:ns5="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-utility-1.0.xsd" ns5:Id="CertId-10107"><ns4:Identifier xsi:type="xsd:anyURI">1000</ns4:Identifier></ns4:SecurityContextToken></ns2:Security></SOAP-ENV:Header><SOAP-ENV:Body xmlns:ns1="http://counter.com"><ns1:createCounter xmlns:ns2="http://docs.oasis-open.org/wss/2004/01/oasis-200401-wss-wssecurity-utility-1.0.xsd" ns2:Id="10106"></ns1:createCounter></SOAP-ENV:Body></SOAP-ENV:Envelope>"""

def main():
    unittest.main(defaultTest="makeTestSuite")

if __name__ == '__main__': 
    main()

