#! /usr/bin/env python
# $Header$
'''Faults.
'''

from ZSI import _copyright, _children, _child_elements, \
        _get_idstr, _stringtypes, _seqtypes, _Node, SoapWriter, ZSIException

from ZSI.TCcompound import Struct
from ZSI.TC import QName, URI, String, XMLString, AnyElement, UNBOUNDED

from ZSI.wstools.Namespaces import SOAP, ZSI_SCHEMA_URI
from ZSI.wstools.c14n import Canonicalize
from ZSI.TC import ElementDeclaration

import traceback, cStringIO as StringIO


class Detail:
    def __init__(self, any=None):
        self.any = any

Detail.typecode = Struct(Detail, [AnyElement(aname='any',minOccurs=0, maxOccurs="unbounded",processContents="lax")], pname='detail', minOccurs=0)

class FaultType:
    def __init__(self, faultcode=None, faultstring=None, faultactor=None, detail=None):
        self.faultcode = faultcode
        self.faultstring= faultstring
        self.faultactor = faultactor
        self.detail = detail
        
FaultType.typecode = \
    Struct(FaultType,
        [QName(pname='faultcode'), 
         String(pname='faultstring'),
         URI(pname=(SOAP.ENV,'faultactor'), minOccurs=0),
         Detail.typecode,
         AnyElement(aname='any',minOccurs=0, maxOccurs=UNBOUNDED, processContents="lax"),
        ], 
        pname=(SOAP.ENV,'Fault'), 
        inline=True,
        hasextras=0, 
    )

class ZSIHeaderDetail:
    def __init__(self, detail):
        self.any = detail

ZSIHeaderDetail.typecode =\
    Struct(ZSIHeaderDetail, 
           [AnyElement(aname='any', minOccurs=0, maxOccurs=UNBOUNDED, processContents="lax")], 
           pname=(ZSI_SCHEMA_URI, 'detail'))


class ZSIFaultDetailTypeCode(ElementDeclaration, Struct):
    '''<ZSI:FaultDetail>
           <ZSI:string>%s</ZSI:string>
           <ZSI:trace>%s</ZSI:trace>
       </ZSI:FaultDetail>
    '''
    schema = ZSI_SCHEMA_URI
    literal = 'FaultDetail'
    def __init__(self, **kw):
        Struct.__init__(self, ZSIFaultDetail, [String(pname=(ZSI_SCHEMA_URI, 'string')), 
            String(pname=(ZSI_SCHEMA_URI, 'trace'),minOccurs=0),], 
            pname=(ZSI_SCHEMA_URI, 'FaultDetail'), **kw
        )

class ZSIFaultDetail:
    def __init__(self, string=None, trace=None):
        self.string = string
        self.trace = trace

    def __str__(self):
        if self.trace:
            return self.string + '\n[trace: ' + self.trace + ']'
        return self.string

    def __repr__(self):
        return "<%s.ZSIFaultDetail %s>" % (__name__, _get_idstr(self))
ZSIFaultDetail.typecode = ZSIFaultDetailTypeCode()


class URIFaultDetailTypeCode(ElementDeclaration, Struct):
    '''
    <ZSI:URIFaultDetail>
        <ZSI:URI>uri</ZSI:URI>
        <ZSI:localname>localname</ZSI:localname>
    </ZSI:URIFaultDetail>
    '''
    schema = ZSI_SCHEMA_URI
    literal = 'URIFaultDetail'
    def __init__(self, **kw):
        Struct.__init__(self, URIFaultDetail, 
           [String(pname=(ZSI_SCHEMA_URI, 'URI')), String(pname=(ZSI_SCHEMA_URI, 'localname')),], 
            pname=(ZSI_SCHEMA_URI, 'URIFaultDetail'), **kw
        )

class URIFaultDetail:
    def __init__(self, uri=None, localname=None):
        self.URI = uri
        self.localname = localname
URIFaultDetail.typecode = URIFaultDetailTypeCode()


class ActorFaultDetailTypeCode(ElementDeclaration, Struct):
    '''
    <ZSI:ActorFaultDetail>
        <ZSI:URI>%s</ZSI:URI>
    </ZSI:ActorFaultDetail>
    '''
    schema = ZSI_SCHEMA_URI
    literal = 'ActorFaultDetail'
    def __init__(self, **kw):
        Struct.__init__(self, ActorFaultDetail, [String(pname=(ZSI_SCHEMA_URI, 'URI')),],
            pname=(ZSI_SCHEMA_URI, 'ActorFaultDetail'), **kw
        )

class ActorFaultDetail:
    def __init__(self, uri=None):
        self.URI = uri
ActorFaultDetail.typecode = ActorFaultDetailTypeCode()


class Fault(ZSIException):
    '''SOAP Faults.
    '''

    Client = "SOAP-ENV:Client"
    Server = "SOAP-ENV:Server"
    MU     = "SOAP-ENV:MustUnderstand"

    def __init__(self, code, string,
                actor=None, detail=None, headerdetail=None):
        if detail is not None and type(detail) not in _seqtypes:
            detail = (detail,)
        if headerdetail is not None and type(headerdetail) not in _seqtypes:
            headerdetail = (headerdetail,)
        self.code, self.string, self.actor, self.detail, self.headerdetail = \
                code, string, actor, detail, headerdetail
        ZSIException.__init__(self, code, string, actor, detail, headerdetail)

    def DataForSOAPHeader(self):
        if not self.headerdetail: return None
        # SOAP spec doesn't say how to encode header fault data.
        return ZSIHeaderDetail(self.headerdetail)

    def serialize(self, sw):
        '''Serialize the object.'''
        detail = None
        if self.detail is not None: 
            detail = Detail()
            detail.any = self.detail

        pyobj = FaultType(self.code, self.string, self.actor, detail)
        sw.serialize(pyobj, typed=False)

    def AsSOAP(self, **kw):

        header = self.DataForSOAPHeader() 
        sw = SoapWriter(**kw)
        self.serialize(sw)
        if header is not None:
            sw.serialize_header(header, header.typecode, typed=False)
        return str(sw)

    def __str__(self):
        strng = str(self.string) + "\n"
        if hasattr(self, 'detail'):
            if hasattr(self.detail, '__len__'):
                for d in self.detail:
                    strng += str(d)
            else:
                strng += str(self.detail)
        return strng

    def __repr__(self):
        return "<%s.Fault at %s>" % (__name__, _get_idstr(self))

    AsSoap = AsSOAP


def FaultFromNotUnderstood(uri, localname, actor=None):
    detail, headerdetail = None, URIFaultDetail(uri, localname)
    return Fault(Fault.MU, 'SOAP mustUnderstand not understood',
                actor, detail, headerdetail)


def FaultFromActor(uri, actor=None):
    detail, headerdetail = None, ActorFaultDetail(uri)
    return Fault(Fault.Client, 'Cannot process specified actor',
                actor, detail, headerdetail)


def FaultFromZSIException(ex, actor=None):
    '''Return a Fault object created from a ZSI exception object.
    '''
    mystr = getattr(ex, 'str', None) or str(ex)
    mytrace = getattr(ex, 'trace', '')
    elt = '''<ZSI:ParseFaultDetail>
<ZSI:string>%s</ZSI:string>
<ZSI:trace>%s</ZSI:trace>
</ZSI:ParseFaultDetail>
''' % (mystr, mytrace)
    if getattr(ex, 'inheader', 0):
        detail, headerdetail = None, elt
    else:
        detail, headerdetail = elt, None
    return Fault(Fault.Client, 'Unparseable message',
                actor, detail, headerdetail)


def FaultFromException(ex, inheader, tb=None, actor=None):
    '''Return a Fault object created from a Python exception.

    <SOAP-ENV:Fault>
      <faultcode>SOAP-ENV:Server</faultcode>
      <faultstring>Processing Failure</faultstring>
      <detail>
        <ZSI:FaultDetail>
          <ZSI:string></ZSI:string>
          <ZSI:trace></ZSI:trace>
        </ZSI:FaultDetail>
      </detail>
    </SOAP-ENV:Fault>
    '''
    tracetext = None
    if tb:
        try:
            lines = '\n'.join(['%s:%d:%s' % (name, line, func)
                        for name, line, func, text in traceback.extract_tb(tb)])
        except:
            pass
        else:
            tracetext = lines
  
    exceptionName = ""
    try:
        exceptionName = ":".join([ex.__module__, ex.__class__.__name__])
    except: pass
    elt = ZSIFaultDetail(string=exceptionName + "\n" + str(ex), trace=tracetext)
    if inheader:
        detail, headerdetail = None, elt
    else:
        detail, headerdetail = elt, None
    return Fault(Fault.Server, 'Processing Failure',
                actor, detail, headerdetail)


def FaultFromFaultMessage(ps):
    '''Parse the message as a fault.
    '''
    pyobj = ps.Parse(FaultType.typecode)

    if pyobj.detail == None:  detailany = None
    else:  detailany = pyobj.detail.any

    return Fault(pyobj.faultcode, pyobj.faultstring,
                pyobj.faultactor, detailany)


if __name__ == '__main__': print _copyright
