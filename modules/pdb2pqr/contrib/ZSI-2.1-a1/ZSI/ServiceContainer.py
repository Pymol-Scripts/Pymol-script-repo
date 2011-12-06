#! /usr/bin/env python
'''Simple Service Container
   -- use with wsdl2py generated modules.
'''

import urlparse, types, os, sys, cStringIO as StringIO, thread,re
from BaseHTTPServer import BaseHTTPRequestHandler, HTTPServer
from ZSI import ParseException, FaultFromException, FaultFromZSIException, Fault
from ZSI import _copyright, _seqtypes, _get_element_nsuri_name, resolvers
from ZSI import _get_idstr
from ZSI.address import Address
from ZSI.parse import ParsedSoap
from ZSI.writer import SoapWriter
from ZSI.dispatch import _ModPythonSendXML, _ModPythonSendFault, _CGISendXML, _CGISendFault
from ZSI.dispatch import SOAPRequestHandler as BaseSOAPRequestHandler

"""
Functions:
    _Dispatch
    AsServer
    GetSOAPContext

Classes:
    SOAPContext
    NoSuchService
    PostNotSpecified
    SOAPActionNotSpecified
    ServiceSOAPBinding
    WSAResource
    SimpleWSResource
    SOAPRequestHandler
    ServiceContainer
"""
class NoSuchService(Exception): pass
class UnknownRequestException(Exception): pass
class PostNotSpecified(Exception): pass
class SOAPActionNotSpecified(Exception): pass
class WSActionException(Exception): pass
class WSActionNotSpecified(WSActionException): pass
class NotAuthorized(Exception): pass
class ServiceAlreadyPresent(Exception): pass


class SOAPContext:
    def __init__(self, container, xmldata, ps, connection, httpheaders,
                 soapaction):

        self.container = container
        self.xmldata    = xmldata
        self.parsedsoap = ps
        self.connection = connection
        self.httpheaders= httpheaders
        self.soapaction = soapaction

_contexts = dict()
def GetSOAPContext():
    global _contexts
    return _contexts[thread.get_ident()]

def _Dispatch(ps, server, SendResponse, SendFault, post, action, nsdict={}, **kw):
    '''Send ParsedSoap instance to ServiceContainer, which dispatches to
    appropriate service via post, and method via action.  Response is a
    self-describing pyobj, which is passed to a SoapWriter.

    Call SendResponse or SendFault to send the reply back, appropriately.
        server -- ServiceContainer instance

    '''
    localURL = 'http://%s:%d%s' %(server.server_name,server.server_port,post)
    address = action
    service = server.getNode(post)
    isWSResource = False
    if isinstance(service, WSAResource):
        isWSResource = True
        service.setServiceURL(localURL)
        address = Address()
        try:
            address.parse(ps)
        except Exception, e:
            return SendFault(FaultFromException(e, 0, sys.exc_info()[2]), **kw)
        if action and action != address.getAction():
            e = WSActionException('SOAP Action("%s") must match WS-Action("%s") if specified.' \
                %(action,address.getAction()))
            return SendFault(FaultFromException(e, 0, None), **kw)
        action = address.getAction()

    if isinstance(service, ServiceInterface) is False:
        e = NoSuchService('no service at POST(%s) in container: %s' %(post,server))
        return SendFault(FaultFromException(e, 0, sys.exc_info()[2]), **kw)

    if not service.authorize(None, post, action):
        return SendFault(Fault(Fault.Server, "Not authorized"), code=401)
        #try:
        #    raise NotAuthorized()
        #except Exception, e:
            #return SendFault(FaultFromException(e, 0, None), code=401, **kw)
            ##return SendFault(FaultFromException(NotAuthorized(), 0, None), code=401, **kw)

    try:
        method = service.getOperation(ps, address)
    except Exception, e:
        return SendFault(FaultFromException(e, 0, sys.exc_info()[2]), **kw)

    try:
        if isWSResource is True: 
            request,result = method(ps, address)
        else: 
            request,result = method(ps)
    except Exception, e:
        return SendFault(FaultFromException(e, 0, sys.exc_info()[2]), **kw)

    # Verify if Signed
    service.verify(ps)

    # If No response just return.
    if result is None:
        return SendResponse('', **kw)

    sw = SoapWriter(nsdict=nsdict)
    try:
        sw.serialize(result)
    except Exception, e:
        return SendFault(FaultFromException(e, 0, sys.exc_info()[2]), **kw)

    if isWSResource is True:
        action = service.getResponseAction(ps, action)
        addressRsp = Address(action=action)
        try:
            addressRsp.setResponseFromWSAddress(address, localURL)
            addressRsp.serialize(sw)
        except Exception, e:
            return SendFault(FaultFromException(e, 0, sys.exc_info()[2]), **kw)

    # Create Signatures
    service.sign(sw)

    try:
        soapdata = str(sw)
        return SendResponse(soapdata, **kw)
    except Exception, e:
        return SendFault(FaultFromException(e, 0, sys.exc_info()[2]), **kw)


def AsServer(port=80, services=()):
    '''port --
       services -- list of service instances
    '''
    address = ('', port)
    sc = ServiceContainer(address, services)
    sc.serve_forever()


class ServiceInterface:
    '''Defines the interface for use with ServiceContainer Handlers.

    class variables:
        soapAction -- dictionary of soapAction keys, and operation name values.
           These are specified in the WSDL soap bindings. There must be a 
           class method matching the operation name value.  If WS-Action is
           used the keys are WS-Action request values, according to the spec
           if soapAction and WS-Action is specified they must be equal.
           
        wsAction -- dictionary of operation name keys and WS-Action 
           response values.  These values are specified by the portType.

        root -- dictionary of root element keys, and operation name values.

    '''
    soapAction = {}
    wsAction = {}
    root = {}

    def __init__(self, post):
        self.post = post

    def authorize(self, auth_info, post, action):
        return 1

    def __str__(self):
        return '%s(%s) POST(%s)' %(self.__class__.__name__, _get_idstr(self), self.post)

    def sign(self, sw):
        return

    def verify(self, ps):
        return

    def getPost(self):
        return self.post

    def getOperation(self, ps, action):
        '''Returns a method of class.
           action -- soapAction value
        '''
        opName = self.getOperationName(ps, action)
        return getattr(self, opName)

    def getOperationName(self, ps, action):
        '''Returns operation name.
           action -- soapAction value
        '''
        method = self.root.get(_get_element_nsuri_name(ps.body_root)) or \
            self.soapAction.get(action)
        if method is None:
            raise UnknownRequestException, \
                'failed to map request to a method: action(%s), root%s' %(action,_get_element_nsuri_name(ps.body_root))
        return method


class ServiceSOAPBinding(ServiceInterface):
    '''Binding defines the set of wsdl:binding operations, it takes as input a
    ParsedSoap instance and parses it into a pyobj.  It returns a response pyobj.
    '''
    def __init__(self, post):
        ServiceInterface.__init__(self, post)

    def __call___(self, action, ps):
        return self.getOperation(ps, action)(ps)


class WSAResource(ServiceSOAPBinding):
    '''Simple WSRF service, performs method resolutions based
    on WS-Action values rather than SOAP Action.

    class variables:
        encoding  
        wsAction -- Must override to set output Action values.
        soapAction -- Must override to set input Action values.
    '''
    encoding = "UTF-8"

    def __init__(self, post):
        '''
        post -- POST value
        '''
        assert isinstance(self.soapAction, dict), "soapAction must be a dict"
        assert isinstance(self.wsAction, dict), "wsAction must be a dict"
        ServiceSOAPBinding.__init__(self, post)

    def __call___(self, action, ps, address):
        return self.getOperation(ps, action)(ps, address)

    def getServiceURL(self):
        return self._url

    def setServiceURL(self, url):
        self._url = url

    def getOperation(self, ps, address):
        '''Returns a method of class.
        address -- ws-address 
        '''
        action = address.getAction()
        opName = self.getOperationName(ps, action)
        return getattr(self, opName)

    def getResponseAction(self, ps, action):
        '''Returns response WS-Action if available
           action -- request WS-Action value.
        '''
        opName = self.getOperationName(ps, action)
        if self.wsAction.has_key(opName) is False:
            raise WSActionNotSpecified, 'wsAction dictionary missing key(%s)' %opName
        return self.wsAction[opName]

    def do_POST(self):
        '''The POST command.  This is called by HTTPServer, not twisted.
        action -- SOAPAction(HTTP header) or wsa:Action(SOAP:Header)
        '''
        global _contexts

        soapAction = self.headers.getheader('SOAPAction')
        post = self.path
        if not post:
            raise PostNotSpecified, 'HTTP POST not specified in request'
        if soapAction:
            soapAction = soapAction.strip('\'"')
        post = post.strip('\'"')
        try:
            ct = self.headers['content-type']
            if ct.startswith('multipart/'):
                cid = resolvers.MIMEResolver(ct, self.rfile)
                xml = cid.GetSOAPPart()
                ps = ParsedSoap(xml, resolver=cid.Resolve, readerclass=DomletteReader)
            else:
                length = int(self.headers['content-length'])
                ps = ParsedSoap(self.rfile.read(length), readerclass=DomletteReader)
        except ParseException, e:
            self.send_fault(FaultFromZSIException(e))
        except Exception, e:
            # Faulted while processing; assume it's in the header.
            self.send_fault(FaultFromException(e, 1, sys.exc_info()[2]))
        else:
            # Keep track of calls
            thread_id = thread.get_ident()
            _contexts[thread_id] = SOAPContext(self.server, xml, ps,
                                               self.connection,
                                               self.headers, soapAction)

            try:
                _Dispatch(ps, self.server, self.send_xml, self.send_fault,
                    post=post, action=soapAction)
            except Exception, e:
                self.send_fault(FaultFromException(e, 0, sys.exc_info()[2]))

            # Clean up after the call
            if _contexts.has_key(thread_id):
                del _contexts[thread_id]


class SOAPRequestHandler(BaseSOAPRequestHandler):
    '''SOAP handler.
    '''
    def do_POST(self):
        '''The POST command.
        action -- SOAPAction(HTTP header) or wsa:Action(SOAP:Header)
        '''
        soapAction = self.headers.getheader('SOAPAction')
        post = self.path
        if not post:
            raise PostNotSpecified, 'HTTP POST not specified in request'
        if soapAction:
            soapAction = soapAction.strip('\'"')
        post = post.strip('\'"')
        try:
            ct = self.headers['content-type']
            if ct.startswith('multipart/'):
                cid = resolvers.MIMEResolver(ct, self.rfile)
                xml = cid.GetSOAPPart()
                ps = ParsedSoap(xml, resolver=cid.Resolve)
            else:
                length = int(self.headers['content-length'])
                xml = self.rfile.read(length)
                ps = ParsedSoap(xml)
        except ParseException, e:
            self.send_fault(FaultFromZSIException(e))
        except Exception, e:
            # Faulted while processing; assume it's in the header.
            self.send_fault(FaultFromException(e, 1, sys.exc_info()[2]))
        else:
            # Keep track of calls
            thread_id = thread.get_ident()
            _contexts[thread_id] = SOAPContext(self.server, xml, ps,
                                               self.connection,
                                               self.headers, soapAction)

            try:
                _Dispatch(ps, self.server, self.send_xml, self.send_fault, 
                    post=post, action=soapAction)
            except Exception, e:
                self.send_fault(FaultFromException(e, 0, sys.exc_info()[2]))

            # Clean up after the call
            if _contexts.has_key(thread_id):
                del _contexts[thread_id]

    def do_GET(self):
        '''The GET command.
	'''
        if self.path.lower().endswith("?wsdl"):
            service_path = self.path[:-5]
            service = self.server.getNode(service_path)
            if hasattr(service, "_wsdl"):
                wsdl = service._wsdl
                # update the soap:location tag in the wsdl to the actual server
                #   location
                # - default to 'http' as protocol, or use server-specified protocol
                proto = 'http'
                if hasattr(self.server,'proto'):
                    proto = self.server.proto
                serviceUrl = '%s://%s:%d%s' % (proto,
                                                self.server.server_name,
                                                self.server.server_port,
                                                service_path)
                soapAddress = '<soap:address location="%s"/>' % serviceUrl
                wsdlre = re.compile('\<soap:address[^\>]*>',re.IGNORECASE)
                wsdl = re.sub(wsdlre,soapAddress,wsdl)
                self.send_xml(wsdl)
            else:
                self.send_error(404, "WSDL not available for that service [%s]." % self.path)
        else:
            self.send_error(404, "Service not found [%s]." % self.path)

class ServiceContainer(HTTPServer):
    '''HTTPServer that stores service instances according 
    to POST values.  An action value is instance specific,
    and specifies an operation (function) of an instance.
    '''
    class NodeTree:
        '''Simple dictionary implementation of a node tree
        '''
        def __init__(self):
            self.__dict = {}

        def __str__(self):
            return str(self.__dict)

	def listNodes(self):
	    print self.__dict.keys()

        def getNode(self, url):
            path = urlparse.urlsplit(url)[2]
            if path.startswith("/"):
                path = path[1:]

            if self.__dict.has_key(path):
                return self.__dict[path]
            else:
                raise NoSuchService, 'No service(%s) in ServiceContainer' %path

        def setNode(self, service, url):
            path = urlparse.urlsplit(url)[2]
            if path.startswith("/"):
                path = path[1:]

            if not isinstance(service, ServiceSOAPBinding):
               raise TypeError, 'A Service must implement class ServiceSOAPBinding'
            if self.__dict.has_key(path):
                raise ServiceAlreadyPresent, 'Service(%s) already in ServiceContainer' % path
            else:
                self.__dict[path] = service

        def removeNode(self, url):
            path = urlparse.urlsplit(url)[2]
            if path.startswith("/"):
                path = path[1:]

            if self.__dict.has_key(path):
                node = self.__dict[path]
                del self.__dict[path]
                return node
            else:
                raise NoSuchService, 'No service(%s) in ServiceContainer' %path
            
    def __init__(self, server_address, services=[], RequestHandlerClass=SOAPRequestHandler):
        '''server_address -- 
           RequestHandlerClass -- 
        '''
        HTTPServer.__init__(self, server_address, RequestHandlerClass)
        self._nodes = self.NodeTree()
        map(lambda s: self.setNode(s), services)

    def __str__(self):
        return '%s(%s) nodes( %s )' %(self.__class__, _get_idstr(self), str(self._nodes))

    def __call__(self, ps, post, action, address=None):
        '''ps -- ParsedSoap representing the request
           post -- HTTP POST --> instance
           action -- Soap Action header --> method
           address -- Address instance representing WS-Address 
        '''
        method = self.getCallBack(ps, post, action)
        if (isinstance(method.im_self, WSAResource) or 
            isinstance(method.im_self, SimpleWSResource)):
            return method(ps, address)
        return method(ps)


    def setNode(self, service, url=None):
        if url is None: 
            url = service.getPost()
        self._nodes.setNode(service, url)

    def getNode(self, url):
        return self._nodes.getNode(url)

    def removeNode(self, url):
        self._nodes.removeNode(url)


class SimpleWSResource(ServiceSOAPBinding):

    def getNode(self, post):
        '''post -- POST HTTP value
        '''
        return self._nodes.getNode(post)

    def setNode(self, service, post):
        '''service -- service instance
           post -- POST HTTP value
        '''
        self._nodes.setNode(service, post)

    def getCallBack(self, ps, post, action):
        '''post -- POST HTTP value
           action -- SOAP Action value
        '''
        node = self.getNode(post)
        if node is None:
            raise NoSuchFunction
        if node.authorize(None, post, action):
            return node.getOperation(ps, action)
        else:
            raise NotAuthorized, "Authorization failed for method %s" % action


if __name__ == '__main__': print _copyright
