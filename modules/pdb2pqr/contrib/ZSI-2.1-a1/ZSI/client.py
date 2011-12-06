#! /usr/bin/env python
# $Header$
#
# Copyright (c) 2001 Zolera Systems.  All rights reserved.

from ZSI import _copyright, _seqtypes, ParsedSoap, SoapWriter, TC, ZSI_SCHEMA_URI,\
    EvaluateException, FaultFromFaultMessage, _child_elements, _attrs, _find_arraytype,\
    _find_type, _get_idstr, _get_postvalue_from_absoluteURI, FaultException, WSActionException,\
    UNICODE_ENCODING
from ZSI.auth import AUTH
from ZSI.TC import AnyElement, AnyType, String, TypeCode, _get_global_element_declaration,\
    _get_type_definition
from ZSI.TCcompound import Struct
import base64, httplib, Cookie, types, time, urlparse
from ZSI.address import Address
from ZSI.wstools.logging import getLogger as _GetLogger
_b64_encode = base64.encodestring

class _AuthHeader:
    """<BasicAuth xmlns="ZSI_SCHEMA_URI">
           <Name>%s</Name><Password>%s</Password>
       </BasicAuth>
    """
    def __init__(self, name=None, password=None):
        self.Name = name
        self.Password = password
_AuthHeader.typecode = Struct(_AuthHeader, ofwhat=(String((ZSI_SCHEMA_URI,'Name'), typed=False), 
        String((ZSI_SCHEMA_URI,'Password'), typed=False)), pname=(ZSI_SCHEMA_URI,'BasicAuth'), 
        typed=False)
  

class _Caller:
    '''Internal class used to give the user a callable object
    that calls back to the Binding object to make an RPC call.
    '''

    def __init__(self, binding, name, namespace=None):
        self.binding = binding
        self.name = name
        self.namespace = namespace

    def __call__(self, *args):
        nsuri = self.namespace
        if nsuri is None:
            return self.binding.RPC(None, self.name, args, 
                            encodingStyle="http://schemas.xmlsoap.org/soap/encoding/",
                            replytype=TC.Any(self.name+"Response"))
            
        return self.binding.RPC(None, (nsuri,self.name), args, 
                   encodingStyle="http://schemas.xmlsoap.org/soap/encoding/",
                   replytype=TC.Any((nsuri,self.name+"Response")))
    

class _NamedParamCaller:
    '''Similar to _Caller, expect that there are named parameters
    not positional.
    '''
        
    def __init__(self, binding, name, namespace=None):
        self.binding = binding
        self.name = name
        self.namespace = namespace
        
    def __call__(self, **params):
        # Pull out arguments that Send() uses
        kw = {}
        for key in [ 'auth_header', 'nsdict', 'requesttypecode', 'soapaction' ]:
            if params.has_key(key):
                kw[key] = params[key]
                del params[key]
        
        nsuri = self.namespace
        if nsuri is None:
            return self.binding.RPC(None, self.name, None, 
                        encodingStyle="http://schemas.xmlsoap.org/soap/encoding/",
                        _args=params, 
                        replytype=TC.Any(self.name+"Response", aslist=False),
                        **kw)      
                
        return self.binding.RPC(None, (nsuri,self.name), None, 
                   encodingStyle="http://schemas.xmlsoap.org/soap/encoding/",
                   _args=params, 
                   replytype=TC.Any((nsuri,self.name+"Response"), aslist=False),
                   **kw)


class _Binding:
    '''Object that represents a binding (connection) to a SOAP server.
    Once the binding is created, various ways of sending and
    receiving SOAP messages are available.
    '''
    defaultHttpTransport = httplib.HTTPConnection
    defaultHttpsTransport = httplib.HTTPSConnection
    logger = _GetLogger('ZSI.client.Binding')

    def __init__(self, nsdict=None, transport=None, url=None, tracefile=None,
                 readerclass=None, writerclass=None, soapaction='', 
                 wsAddressURI=None, sig_handler=None, transdict=None, **kw):
        '''Initialize.
        Keyword arguments include:
            transport -- default use HTTPConnection. 
            transdict -- dict of values to pass to transport.
            url -- URL of resource, POST is path 
            soapaction -- value of SOAPAction header
            auth -- (type, name, password) triplet; default is unauth
            nsdict -- namespace entries to add
            tracefile -- file to dump packet traces
            cert_file, key_file -- SSL data (q.v.)
            readerclass -- DOM reader class
            writerclass -- DOM writer class, implements MessageInterface
            wsAddressURI -- namespaceURI of WS-Address to use.  By default 
            it's not used.
            sig_handler -- XML Signature handler, must sign and verify.
            endPointReference -- optional Endpoint Reference.
        '''
        self.data = None
        self.ps = None
        self.user_headers = []
        self.nsdict = nsdict or {}
        self.transport = transport
        self.transdict = transdict or {}
        self.url = url
        self.trace = tracefile
        self.readerclass = readerclass
        self.writerclass = writerclass
        self.soapaction = soapaction
        self.wsAddressURI = wsAddressURI
        self.sig_handler = sig_handler
        self.address = None
        self.endPointReference = kw.get('endPointReference', None)
        self.cookies = Cookie.SimpleCookie()
        self.http_callbacks = {}

        if kw.has_key('auth'):
            self.SetAuth(*kw['auth'])
        else:
            self.SetAuth(AUTH.none)

    def SetAuth(self, style, user=None, password=None):
        '''Change auth style, return object to user.
        '''
        self.auth_style, self.auth_user, self.auth_pass = \
            style, user, password
        return self

    def SetURL(self, url):
        '''Set the URL we post to.
        '''
        self.url = url
        return self

    def ResetHeaders(self):
        '''Empty the list of additional headers.
        '''
        self.user_headers = []
        return self

    def ResetCookies(self):
        '''Empty the list of cookies.
        '''
        self.cookies = Cookie.SimpleCookie()

    def AddHeader(self, header, value):
        '''Add a header to send.
        '''
        self.user_headers.append((header, value))
        return self

    def __addcookies(self):
        '''Add cookies from self.cookies to request in self.h
        '''
        for cname, morsel in self.cookies.items():
            attrs = []
            value = morsel.get('version', '')
            if value != '' and value != '0':
                attrs.append('$Version=%s' % value)
            attrs.append('%s=%s' % (cname, morsel.coded_value))
            value = morsel.get('path')
            if value:
                attrs.append('$Path=%s' % value)
            value = morsel.get('domain')
            if value:
                attrs.append('$Domain=%s' % value)
            self.h.putheader('Cookie', "; ".join(attrs))

    def RPC(self, url, opname, obj, replytype=None, **kw):
        '''Send a request, return the reply.  See Send() and Recieve()
        docstrings for details.
        '''
        self.Send(url, opname, obj, **kw)
        return self.Receive(replytype, **kw)

    def Send(self, url, opname, obj, nsdict={}, soapaction=None, wsaction=None, 
             endPointReference=None, soapheaders=(), **kw):
        '''Send a message.  If url is None, use the value from the
        constructor (else error). obj is the object (data) to send.
        Data may be described with a requesttypecode keyword, the default 
        is the class's typecode (if there is one), else Any.

        Try to serialize as a Struct, if this is not possible serialize an Array.  If 
        data is a sequence of built-in python data types, it will be serialized as an
        Array, unless requesttypecode is specified.

        arguments:
            url -- 
            opname -- struct wrapper
            obj -- python instance

        key word arguments:
            nsdict -- 
            soapaction --
            wsaction -- WS-Address Action, goes in SOAP Header.
            endPointReference --  set by calling party, must be an 
                EndPointReference type instance.
            soapheaders -- list of pyobj, typically w/typecode attribute.
                serialized in the SOAP:Header.
            requesttypecode -- 

        '''
        url = url or self.url
        endPointReference = endPointReference or self.endPointReference

        # Serialize the object.
        d = {}
        d.update(self.nsdict)
        d.update(nsdict)

        sw = SoapWriter(nsdict=d, header=True, outputclass=self.writerclass, 
                 encodingStyle=kw.get('encodingStyle'),)
        
        requesttypecode = kw.get('requesttypecode')
        if kw.has_key('_args'): #NamedParamBinding
            tc = requesttypecode or TC.Any(pname=opname, aslist=False)
            sw.serialize(kw['_args'], tc)
        elif not requesttypecode:
            tc = getattr(obj, 'typecode', None) or TC.Any(pname=opname, aslist=False)
            try:
                if type(obj) in _seqtypes:
                    obj = dict(map(lambda i: (i.typecode.pname,i), obj))
            except AttributeError:
                # can't do anything but serialize this in a SOAP:Array
                tc = TC.Any(pname=opname, aslist=True)
            else:
                tc = TC.Any(pname=opname, aslist=False)

            sw.serialize(obj, tc)
        else:
            sw.serialize(obj, requesttypecode)
            
        for i in soapheaders:
           sw.serialize_header(i) 
            
        # 
        # Determine the SOAP auth element.  SOAP:Header element
        if self.auth_style & AUTH.zsibasic:
            sw.serialize_header(_AuthHeader(self.auth_user, self.auth_pass),
                _AuthHeader.typecode)

        # 
        # Serialize WS-Address
        if self.wsAddressURI is not None:
            if self.soapaction and wsaction.strip('\'"') != self.soapaction:
                raise WSActionException, 'soapAction(%s) and WS-Action(%s) must match'\
                    %(self.soapaction,wsaction)

            self.address = Address(url, self.wsAddressURI)
            self.address.setRequest(endPointReference, wsaction)
            self.address.serialize(sw)

        # 
        # WS-Security Signature Handler
        if self.sig_handler is not None:
            self.sig_handler.sign(sw)

        scheme,netloc,path,nil,nil,nil = urlparse.urlparse(url)
        transport = self.transport
        if transport is None and url is not None:
            if scheme == 'https':
                transport = self.defaultHttpsTransport
            elif scheme == 'http':
                transport = self.defaultHttpTransport
            else:
                raise RuntimeError, 'must specify transport or url startswith https/http'

        # Send the request.
        if issubclass(transport, httplib.HTTPConnection) is False:
            raise TypeError, 'transport must be a HTTPConnection'

        soapdata = str(sw)
        self.h = transport(netloc, None, **self.transdict)
        self.h.connect()
        self.SendSOAPData(soapdata, url, soapaction, **kw)

    def SendSOAPData(self, soapdata, url, soapaction, headers={}, **kw):
        # Tracing?
        if self.trace:
            print >>self.trace, "_" * 33, time.ctime(time.time()), "REQUEST:"
            print >>self.trace, soapdata

        url = url or self.url
        request_uri = _get_postvalue_from_absoluteURI(url)
        self.h.putrequest("POST", request_uri)
        self.h.putheader("Content-Length", "%d" % len(soapdata))
        self.h.putheader("Content-Type", 'text/xml; charset="%s"' %UNICODE_ENCODING)
        self.__addcookies()

        for header,value in headers.items():
            self.h.putheader(header, value)

        SOAPActionValue = '"%s"' % (soapaction or self.soapaction)
        self.h.putheader("SOAPAction", SOAPActionValue)
        if self.auth_style & AUTH.httpbasic:
            val = _b64_encode(self.auth_user + ':' + self.auth_pass) \
                        .replace("\012", "")
            self.h.putheader('Authorization', 'Basic ' + val)
        elif self.auth_style == AUTH.httpdigest and not headers.has_key('Authorization') \
            and not headers.has_key('Expect'):
            def digest_auth_cb(response):
                self.SendSOAPDataHTTPDigestAuth(response, soapdata, url, request_uri, soapaction, **kw)
                self.http_callbacks[401] = None
            self.http_callbacks[401] = digest_auth_cb

        for header,value in self.user_headers:
            self.h.putheader(header, value)
        self.h.endheaders()
        self.h.send(soapdata)

        # Clear prior receive state.
        self.data, self.ps = None, None

    def SendSOAPDataHTTPDigestAuth(self, response, soapdata, url, request_uri, soapaction, **kw):
        '''Resend the initial request w/http digest authorization headers.
        The SOAP server has requested authorization.  Fetch the challenge, 
        generate the authdict for building a response.
        '''
        if self.trace:
            print >>self.trace, "------ Digest Auth Header"
        url = url or self.url
        if response.status != 401:
            raise RuntimeError, 'Expecting HTTP 401 response.'
        if self.auth_style != AUTH.httpdigest:
            raise RuntimeError,\
                'Auth style(%d) does not support requested digest authorization.' %self.auth_style

        from ZSI.digest_auth import fetch_challenge,\
            generate_response,\
            build_authorization_arg,\
            dict_fetch

        chaldict = fetch_challenge( response.getheader('www-authenticate') )
        if dict_fetch(chaldict,'challenge','').lower() == 'digest' and \
            dict_fetch(chaldict,'nonce',None) and \
            dict_fetch(chaldict,'realm',None) and \
            dict_fetch(chaldict,'qop',None):
            authdict = generate_response(chaldict,
                request_uri, self.auth_user, self.auth_pass, method='POST')
            headers = {\
                'Authorization':build_authorization_arg(authdict),
                'Expect':'100-continue',
            }
            self.SendSOAPData(soapdata, url, soapaction, headers, **kw)
            return

        raise RuntimeError,\
            'Client expecting digest authorization challenge.'

    def ReceiveRaw(self, **kw):
        '''Read a server reply, unconverted to any format and return it.
        '''
        if self.data: return self.data
        trace = self.trace
        while 1:
            response = self.h.getresponse()
            self.reply_code, self.reply_msg, self.reply_headers, self.data = \
                response.status, response.reason, response.msg, response.read()
            if trace:
                print >>trace, "_" * 33, time.ctime(time.time()), "RESPONSE:"
                for i in (self.reply_code, self.reply_msg,):
                    print >>trace, str(i)
                print >>trace, "-------"
                print >>trace, str(self.reply_headers)
                print >>trace, self.data
            saved = None
            for d in response.msg.getallmatchingheaders('set-cookie'):
                if d[0] in [ ' ', '\t' ]:
                    saved += d.strip()
                else:
                    if saved: self.cookies.load(saved)
                    saved = d.strip()
            if saved: self.cookies.load(saved)
            if response.status == 401:
                if not callable(self.http_callbacks.get(response.status,None)):
                    raise RuntimeError, 'HTTP Digest Authorization Failed'
                self.http_callbacks[response.status](response)
                continue
            if response.status != 100: break

            # The httplib doesn't understand the HTTP continuation header.
            # Horrible internals hack to patch things up.
            self.h._HTTPConnection__state = httplib._CS_REQ_SENT
            self.h._HTTPConnection__response = None
        return self.data

    def IsSOAP(self):
        if self.ps: return 1
        self.ReceiveRaw()
        mimetype = self.reply_headers.type
        return mimetype == 'text/xml'

    def ReceiveSOAP(self, readerclass=None, **kw):
        '''Get back a SOAP message.
        '''
        if self.ps: return self.ps
        if not self.IsSOAP():
            raise TypeError(
                'Response is "%s", not "text/xml"' % self.reply_headers.type)
        if len(self.data) == 0:
            raise TypeError('Received empty response')

        self.ps = ParsedSoap(self.data, 
                        readerclass=readerclass or self.readerclass, 
                        encodingStyle=kw.get('encodingStyle'))

        if self.sig_handler is not None:
            self.sig_handler.verify(self.ps)

        return self.ps

    def IsAFault(self):
        '''Get a SOAP message, see if it has a fault.
        '''
        self.ReceiveSOAP()
        return self.ps.IsAFault()

    def ReceiveFault(self, **kw):
        '''Parse incoming message as a fault. Raise TypeError if no
        fault found.
        '''
        self.ReceiveSOAP(**kw)
        if not self.ps.IsAFault():
            raise TypeError("Expected SOAP Fault not found")
        return FaultFromFaultMessage(self.ps)

    def Receive(self, replytype, **kw):
        '''Parse message, create Python object.

        KeyWord data:
            faults   -- list of WSDL operation.fault typecodes
            wsaction -- If using WS-Address, must specify Action value we expect to
                receive.
        '''
        self.ReceiveSOAP(**kw)
        if self.ps.IsAFault():
            msg = FaultFromFaultMessage(self.ps)
            raise FaultException(msg)

        tc = replytype
        if hasattr(replytype, 'typecode'):
            tc = replytype.typecode

        reply = self.ps.Parse(tc)
        if self.address is not None:
            self.address.checkResponse(self.ps, kw.get('wsaction'))
        return reply

    def __repr__(self):
        return "<%s instance %s>" % (self.__class__.__name__, _get_idstr(self))


class Binding(_Binding):
    '''Object that represents a binding (connection) to a SOAP server.  
    Can be used in the "name overloading" style.
    
    class attr:
        gettypecode -- funcion that returns typecode from typesmodule,
            can be set so can use whatever mapping you desire.
    '''
    gettypecode = staticmethod(lambda mod,e: getattr(mod, str(e.localName)).typecode)
    logger = _GetLogger('ZSI.client.Binding')

    def __init__(self, url, namespace=None, typesmodule=None, **kw):
        """
        Parameters:
            url -- location of service
            namespace -- optional root element namespace
            typesmodule -- optional response only. dict(name=typecode), 
                lookup for all children of root element.  
        """
        self.typesmodule = typesmodule
        self.namespace = namespace
        
        _Binding.__init__(self, url=url, **kw)

    def __getattr__(self, name):
        '''Return a callable object that will invoke the RPC method
        named by the attribute.
        '''
        if name[:2] == '__' and len(name) > 5 and name[-2:] == '__':
            if hasattr(self, name): return getattr(self, name)
            return getattr(self.__class__, name)
        return _Caller(self, name, self.namespace)

    def __parse_child(self, node):
        '''for rpc-style map each message part to a class in typesmodule
        '''
        try:
            tc = self.gettypecode(self.typesmodule, node)
        except:
            self.logger.debug('didnt find typecode for "%s" in typesmodule: %s', 
                node.localName, self.typesmodule)
            tc = TC.Any(aslist=1)
            return tc.parse(node, self.ps)

        self.logger.debug('parse child with typecode : %s', tc)
        try:
            return tc.parse(node, self.ps)
        except Exception:
            self.logger.debug('parse failed try Any : %s', tc)

        tc = TC.Any(aslist=1)
        return tc.parse(node, self.ps)

    def Receive(self, replytype, **kw):
        '''Parse message, create Python object.

        KeyWord data:
            faults   -- list of WSDL operation.fault typecodes
            wsaction -- If using WS-Address, must specify Action value we expect to
                receive.
        ''' 
        self.ReceiveSOAP(**kw)
        ps = self.ps
        tp = _find_type(ps.body_root)
        isarray = ((type(tp) in (tuple,list) and tp[1] == 'Array') or _find_arraytype(ps.body_root))
        if self.typesmodule is None or isarray:
            return _Binding.Receive(self, replytype, **kw)

        if ps.IsAFault():
            msg = FaultFromFaultMessage(ps)
            raise FaultException(msg)

        tc = replytype
        if hasattr(replytype, 'typecode'):
            tc = replytype.typecode

        #Ignore response wrapper
        reply = {}
        for elt in _child_elements(ps.body_root):
            name = str(elt.localName)
            reply[name] = self.__parse_child(elt)

        if self.address is not None:
            self.address.checkResponse(ps, kw.get('wsaction'))

        return reply


class NamedParamBinding(Binding):
    '''Like Binding, except the argument list for invocation is
    named parameters.
    '''
    logger = _GetLogger('ZSI.client.Binding')

    def __getattr__(self, name):
        '''Return a callable object that will invoke the RPC method
        named by the attribute.
        '''
        if name[:2] == '__' and len(name) > 5 and name[-2:] == '__':
            if hasattr(self, name): return getattr(self, name)
            return getattr(self.__class__, name)
        return _NamedParamCaller(self, name, self.namespace)


if __name__ == '__main__': print _copyright
