###########################################################################
# Joshua R. Boverhof, LBNL
# See Copyright for copyright notice!
# $Id: WSsecurity.py 1134 2006-02-24 00:23:06Z boverhof $
###########################################################################

import sys, time, warnings
import sha, base64

# twisted & related imports
from zope.interface import classProvides, implements, Interface
from twisted.python import log, failure
from twisted.web.error import NoResource
from twisted.web.server import NOT_DONE_YET
from twisted.internet import reactor
import twisted.web.http
import twisted.web.resource

# ZSI imports
from ZSI import _get_element_nsuri_name, EvaluateException, ParseException
from ZSI.parse import ParsedSoap
from ZSI.writer import SoapWriter
from ZSI.TC import _get_global_element_declaration as GED
from ZSI import fault
from ZSI.wstools.Namespaces import OASIS, DSIG
from WSresource import DefaultHandlerChain, HandlerChainInterface,\
    WSAddressCallbackHandler, DataHandler, WSAddressHandler


#
# Global Element Declarations
# 
UsernameTokenDec = GED(OASIS.WSSE, "UsernameToken")
SecurityDec = GED(OASIS.WSSE, "Security")
SignatureDec = GED(DSIG.BASE, "Signature")
PasswordDec = GED(OASIS.WSSE, "Password")
NonceDec = GED(OASIS.WSSE, "Nonce")
CreatedDec = GED(OASIS.UTILITY, "Created")

if None in [UsernameTokenDec,SecurityDec,SignatureDec,PasswordDec,NonceDec,CreatedDec]:
    raise ImportError, 'required global element(s) unavailable: %s ' %({
        (OASIS.WSSE, "UsernameToken"):UsernameTokenDec,
        (OASIS.WSSE, "Security"):SecurityDec,
        (DSIG.BASE, "Signature"):SignatureDec,
        (OASIS.WSSE, "Password"):PasswordDec,
        (OASIS.WSSE, "Nonce"):NonceDec,
        (OASIS.UTILITY, "Created"):CreatedDec,
        })
    
    
# 
# Stability: Unstable, Untested, Not Finished.
# 

class WSSecurityHandler:
    """Web Services Security: SOAP Message Security 1.0
    
    Class Variables:
        debug -- If True provide more detailed SOAP:Fault information to clients.
    """
    classProvides(HandlerChainInterface)
    debug = True
    
    @classmethod
    def processRequest(cls, ps, **kw):
        if type(ps) is not ParsedSoap:
            raise TypeError,'Expecting ParsedSoap instance'
        
        security = ps.ParseHeaderElements([cls.securityDec])
        
        # Assume all security headers are supposed to be processed here.
        for pyobj in security or []:
            for any in pyobj.Any or []:
                
                if any.typecode is UsernameTokenDec:
                    try:
                        ps = cls.UsernameTokenProfileHandler.processRequest(ps, any)
                    except Exception, ex:
                        if cls.debug: raise
                        raise RuntimeError, 'Unauthorized Username/passphrase combination'
                    continue
                
                if any.typecode is SignatureDec:
                    try:
                        ps = cls.SignatureHandler.processRequest(ps, any)
                    except Exception, ex:
                        if cls.debug: raise
                        raise RuntimeError, 'Invalid Security Header'
                    continue
                
                raise RuntimeError, 'WS-Security, Unsupported token %s' %str(any)
            
        return ps

    @classmethod
    def processResponse(cls, output, **kw):
        return output


    class UsernameTokenProfileHandler:
        """Web Services Security UsernameToken Profile 1.0
        
        Class Variables:
            targetNamespace --
        """
        classProvides(HandlerChainInterface)
        
        # Class Variables
        targetNamespace = OASIS.WSSE
        sweepInterval = 60*5
        nonces = None
            
        # Set to None to disable
        PasswordText = targetNamespace + "#PasswordText"
        PasswordDigest = targetNamespace + "#PasswordDigest"
            
        # Override passwordCallback 
        passwordCallback = lambda cls,username: None
        
        @classmethod
        def sweep(cls, index):
            """remove nonces every sweepInterval.
            Parameters:
                index -- remove all nonces up to this index.
            """
            if cls.nonces is None: 
                cls.nonces = []
            
            seconds = cls.sweepInterval
            cls.nonces = cls.nonces[index:]
            reactor.callLater(seconds, cls.sweep, len(cls.nonces))
        
        @classmethod
        def processRequest(cls, ps, token, **kw):
            """
            Parameters:
                ps -- ParsedSoap instance
                token -- UsernameToken pyclass instance
            """
            if token.typecode is not UsernameTokenDec:
                raise TypeError, 'expecting GED (%s,%s) representation.' %(
                    UsernameTokenDec.nspname, UsernameTokenDec.pname)
                    
            username = token.Username
            
            # expecting only one password
            # may have a nonce and a created
            password = nonce = timestamp = None
            for any in token.Any or []:
                if any.typecode is PasswordDec:
                    password = any
                    continue
                
                if any.typecode is NonceTypeDec:
                    nonce = any
                    continue
                
                if any.typecode is CreatedTypeDec:
                    timestamp = any
                    continue
                
                raise TypeError, 'UsernameTokenProfileHander unexpected %s' %str(any)

            if password is None:
                raise RuntimeError, 'Unauthorized, no password'
            
            # TODO: not yet supporting complexType simpleContent in pyclass_type
            attrs = getattr(password, password.typecode.attrs_aname, {})
            pwtype = attrs.get('Type', cls.PasswordText)
            
            # Clear Text Passwords
            if cls.PasswordText is not None and pwtype == cls.PasswordText:
                if password == cls.passwordCallback(username):
                    return ps
                
                raise RuntimeError, 'Unauthorized, clear text password failed'
            
            if cls.nonces is None: cls.sweep(0)
            if nonce is not None:
                if nonce in cls.nonces:
                    raise RuntimeError, 'Invalid Nonce'
                
                # created was 10 seconds ago or sooner
                if created is not None and created < time.gmtime(time.time()-10):
                    raise RuntimeError, 'UsernameToken created is expired' 
                
                cls.nonces.append(nonce)
            
            # PasswordDigest, recommended that implemenations
            # require a Nonce and Created
            if cls.PasswordDigest is not None and pwtype == cls.PasswordDigest:
                digest = sha.sha()
                for i in (nonce, created, cls.passwordCallback(username)):
                    if i is None: continue
                    digest.update(i)

                if password == base64.encodestring(digest.digest()).strip():
                    return ps
                
                raise RuntimeError, 'Unauthorized, digest failed'
            
            raise RuntimeError, 'Unauthorized, contents of UsernameToken unknown'
            
        @classmethod
        def processResponse(cls, output, **kw):
            return output
        
    @staticmethod
    def hmac_sha1(xml):
        return 
    
    class SignatureHandler:
        """Web Services Security UsernameToken Profile 1.0
        """
        digestMethods = {
            DSIG.BASE+"#sha1":sha.sha,
            }
        signingMethods = {
            DSIG.BASE+"#hmac-sha1":hmac_sha1,
            }
        canonicalizationMethods = {
            DSIG.C14N_EXCL:lambda node: Canonicalize(node, unsuppressedPrefixes=[]),
            DSIG.C14N:lambda node: Canonicalize(node),
            }
            
        @classmethod
        def processRequest(cls, ps, signature, **kw):
            """
            Parameters:
                ps -- ParsedSoap instance
                signature -- Signature pyclass instance
            """
            if token.typecode is not SignatureDec:
                raise TypeError, 'expecting GED (%s,%s) representation.' %(
                    SignatureDec.nspname, SignatureDec.pname)
                    
            si = signature.SignedInfo
            si.CanonicalizationMethod
            calgo = si.CanonicalizationMethod.get_attribute_Algorithm()
            for any in si.CanonicalizationMethod.Any:
                pass
            
            # Check Digest
            si.Reference
            context = XPath.Context.Context(ps.dom, processContents={'wsu':OASIS.UTILITY})
            exp = XPath.Compile('//*[@wsu:Id="%s"]' %si.Reference.get_attribute_URI())
            nodes = exp.evaluate(context)
            if len(nodes) != 1:
                raise RuntimeError, 'A SignedInfo Reference must refer to one node %s.' %(
                    si.Reference.get_attribute_URI())
                    
            try:
                xml = cls.canonicalizeMethods[calgo](nodes[0])
            except IndexError:
                raise RuntimeError, 'Unsupported canonicalization algorithm'
            
            try:
                digest = cls.digestMethods[salgo]
            except IndexError:
                raise RuntimeError, 'unknown digestMethods Algorithm'
            
            digestValue = base64.encodestring(digest(xml).digest()).strip()
            if si.Reference.DigestValue != digestValue:
                raise RuntimeError, 'digest does not match'
            
            if si.Reference.Transforms:
                pass
            
            signature.KeyInfo
            signature.KeyInfo.KeyName
            signature.KeyInfo.KeyValue
            signature.KeyInfo.RetrievalMethod
            signature.KeyInfo.X509Data
            signature.KeyInfo.PGPData
            signature.KeyInfo.SPKIData
            signature.KeyInfo.MgmtData
            signature.KeyInfo.Any 
            
            signature.Object
            
            # TODO: Check Signature
            signature.SignatureValue
            si.SignatureMethod
            salgo = si.SignatureMethod.get_attribute_Algorithm()
            if si.SignatureMethod.HMACOutputLength:
                pass
            for any in si.SignatureMethod.Any:
                pass
            
            # <SignedInfo><Reference URI="">
            exp = XPath.Compile('//child::*[attribute::URI = "%s"]/..' %(
                                 si.Reference.get_attribute_URI()))
            nodes = exp.evaluate(context)
            if len(nodes) != 1:
                raise RuntimeError, 'A SignedInfo Reference must refer to one node %s.' %(
                    si.Reference.get_attribute_URI())
                    
            try:
                xml = cls.canonicalizeMethods[calgo](nodes[0])
            except IndexError:
                raise RuntimeError, 'Unsupported canonicalization algorithm'
            
            # TODO: Check SignatureValue
            
        @classmethod
        def processResponse(cls, output, **kw):
            return output
        

    class X509TokenProfileHandler:
        """Web Services Security UsernameToken Profile 1.0
        """
        targetNamespace = DSIG.BASE
        
        # Token Types
        singleCertificate = targetNamespace + "#X509v3"
        certificatePath = targetNamespace + "#X509PKIPathv1"
        setCerticatesCRLs = targetNamespace + "#PKCS7"
        
        @classmethod
        def processRequest(cls, ps, signature, **kw):
            return ps



"""
<element name="KeyInfo" type="ds:KeyInfoType"/>
<complexType name="KeyInfoType" mixed="true">
  <choice maxOccurs="unbounded">
    <element ref="ds:KeyName"/>
    <element ref="ds:KeyValue"/>
    <element ref="ds:RetrievalMethod"/>
    <element ref="ds:X509Data"/>
    <element ref="ds:PGPData"/>
    <element ref="ds:SPKIData"/>
    <element ref="ds:MgmtData"/>
    <any processContents="lax" namespace="##other"/>
    <!-- (1,1) elements from (0,unbounded) namespaces -->
  </choice>
  <attribute name="Id" type="ID" use="optional"/>
</complexType>



<element name="Signature" type="ds:SignatureType"/>
<complexType name="SignatureType">
  <sequence>
    <element ref="ds:SignedInfo"/>
    <element ref="ds:SignatureValue"/>
    <element ref="ds:KeyInfo" minOccurs="0"/>
    <element ref="ds:Object" minOccurs="0" maxOccurs="unbounded"/>
  </sequence>
  <attribute name="Id" type="ID" use="optional"/>
</complexType>

  <element name="SignatureValue" type="ds:SignatureValueType"/>
  <complexType name="SignatureValueType">
    <simpleContent>
      <extension base="base64Binary">
        <attribute name="Id" type="ID" use="optional"/>
      </extension>
    </simpleContent>
  </complexType>

<!-- Start SignedInfo -->

<element name="SignedInfo" type="ds:SignedInfoType"/>
<complexType name="SignedInfoType">
  <sequence>
    <element ref="ds:CanonicalizationMethod"/>
    <element ref="ds:SignatureMethod"/>
    <element ref="ds:Reference" maxOccurs="unbounded"/>
  </sequence> 
  <attribute name="Id" type="ID" use="optional"/>
</complexType>
"""


class WSSecurityHandlerChainFactory:
    protocol = DefaultHandlerChain
    
    @classmethod
    def newInstance(cls):
            
        return cls.protocol(WSAddressCallbackHandler, DataHandler, 
            WSSecurityHandler, WSAddressHandler())



