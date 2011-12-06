############################################################################
# Joshua R. Boverhof, LBNL
# See Copyright for copyright notice!
# $Id: $
###########################################################################
import os, sys, types, inspect
from StringIO import StringIO

# twisted & related imports
from zope.interface import classProvides, implements, Interface

# ZSI imports
from ZSI import _get_element_nsuri_name, EvaluateException, ParseException,\
    fault, ParsedSoap, SoapWriter
from ZSI.twisted.reverse import DataHandler, ReverseHandlerChain,\
    HandlerChainInterface

"""
EXAMPLES:

     See zsi/samples/WSGI


"""

def soapmethod(requesttypecode, responsetypecode, soapaction='', 
               operation=None, **kw):
    """@soapmethod
    decorator function for soap methods
    """
    def _closure(func_cb):
        func_cb.root = (requesttypecode.nspname,requesttypecode.pname)
        func_cb.action = soapaction
        func_cb.requesttypecode = requesttypecode
        func_cb.responsetypecode = responsetypecode
        func_cb.soapmethod = True
        func_cb.operation = None
        return func_cb

    return _closure


class SOAPCallbackHandler:
    """ ps --> pyobj, pyobj --> sw
    class variables:
        writerClass -- ElementProxy implementation to use for SoapWriter instances.
    """
    classProvides(HandlerChainInterface)
    writerClass = None

    @classmethod
    def processRequest(cls, ps, **kw):
        """invokes callback that should return a (request,response) tuple.
        representing the SOAP request and response respectively.
        ps -- ParsedSoap instance representing HTTP Body.
        request -- twisted.web.server.Request
        """
        resource = kw['resource']
        request = kw['request']

        root = _get_element_nsuri_name(ps.body_root)
        for key,method in inspect.getmembers(resource, inspect.ismethod):
            if (getattr(method, 'soapmethod', False) and method.root == root):
                break
        else:
            raise RuntimeError, 'Missing soap callback method for root "%s"' %root

        try:
            req = ps.Parse(method.requesttypecode)
        except Exception, ex:
            raise
        try:
            rsp = method.responsetypecode.pyclass()
        except Exception, ex:
            raise
        
        try:
            req,rsp = method(req, rsp)
        except Exception, ex:
            raise

        return rsp

    @classmethod
    def processResponse(cls, output, **kw):
        sw = SoapWriter(outputclass=cls.writerClass)
        sw.serialize(output)
        return sw
    
    
class SOAPHandlerChainFactory:
    protocol = ReverseHandlerChain

    @classmethod
    def newInstance(cls):
        return cls.protocol(DataHandler, SOAPCallbackHandler)


class WSGIApplication(dict):
    encoding = "UTF-8"
    
    def __call__(self, env, start_response):
        """do dispatching, else process
        """
        script = env['SCRIPT_NAME'] # consumed
        ipath = os.path.split(env['PATH_INFO'])[1:]
        for i in range(1, len(ipath)+1):
            path = os.path.join(*ipath[:i])
            print "PATH: ", path
            application = self.get(path)
            if application is not None:
                env['SCRIPT_NAME'] = script + path
                env['PATH_INFO'] =  ''
                print "SCRIPT: ", env['SCRIPT_NAME']
                return application(env, start_response)
            
        return self._request_cb(env, start_response)

    def _request_cb(self, env, start_response):
        """callback method, override
        """
        start_response("404 ERROR", [('Content-Type','text/plain')])
        return ['Move along people, there is nothing to see to hear']
    
    def putChild(self, path, resource):
        """
        """
        path = path.split('/')
        lp = len(path)
        if lp == 0:
            raise RuntimeError, 'bad path "%s"' %path
        
        if lp == 1:
            self[path[0]] = resource
        
        for i in range(len(path)):
            if not path[i]: continue
            break
        
        next = self.get(path[i], None)
        if next is None:
            next = self[path[i]] = WSGIApplication()
            
        next.putChild('/'.join(path[-1:]), resource)
        
        


class SOAPApplication(WSGIApplication):
    """
    """
    factory = SOAPHandlerChainFactory
    
    def __init__(self, **kw):
        dict.__init__(self, **kw)
        self.delegate = None
        
    def _request_cb(self, env, start_response):
        """process request, 
        """
        if env['REQUEST_METHOD'] == 'GET':
            return self._handle_GET(env, start_response)

        if env['REQUEST_METHOD'] == 'POST':
            return self._handle_POST(env, start_response)
            
        start_response("500 ERROR", [('Content-Type','text/plain')])
        s = StringIO()
        h = env.items(); h.sort()
        for k,v in h:
            print >>s, k,'=',`v`
        return [s.getvalue()]

    def _handle_GET(self, env, start_response):
        if env['QUERY_STRING'].lower() == 'wsdl':
            start_response("200 OK", [('Content-Type','text/plain')])
            r = self.delegate or self
            return _resourceToWSDL(r)

        start_response("404 ERROR", [('Content-Type','text/plain')])
        return ['NO RESOURCE FOR GET']
    
    def _handle_POST(self, env, start_response):
        """Dispatch Method called by twisted render, creates a
        request/response handler chain.
        request -- twisted.web.server.Request
        """
        input = env['wsgi.input']
        data = input.read( int(env['CONTENT_LENGTH']) )
        mimeType = "text/xml"
        if self.encoding is not None:
            mimeType = 'text/xml; charset="%s"' % self.encoding

        request = None
        resource = self.delegate or self
        chain = self.factory.newInstance()
        try:
            pyobj = chain.processRequest(data, request=request, resource=resource)
        except Exception, ex:
            start_response("500 ERROR", [('Content-Type',mimeType)])
            return [fault.FaultFromException(ex, False, sys.exc_info()[2]).AsSOAP()]

        try:
            soap = chain.processResponse(pyobj, request=request, resource=resource)
        except Exception, ex:
            start_response("500 ERROR", [('Content-Type',mimeType)])
            return [fault.FaultFromException(ex, False, sys.exc_info()[2]).AsSOAP()]
        
        start_response("200 OK", [('Content-Type',mimeType)])
        return [soap]
    
    
def test(app, port=8080, host="localhost"):
    """
    """
    from twisted.internet import reactor
    from twisted.python import log
    from twisted.web2.channel import HTTPFactory
    from twisted.web2.server import Site
    from twisted.web2.wsgi import WSGIResource
    
    log.startLogging(sys.stdout)
    reactor.listenTCP(port, 
        HTTPFactory( Site(WSGIResource(app)) ),
        interface=host,
    )
    reactor.run()


def _issoapmethod(f):
    return type(f) is types.MethodType and getattr(f, 'soapmethod', False)

def _resourceToWSDL(resource):
    from xml.etree import ElementTree
    from xml.etree.ElementTree import Element, QName
    from ZSI.wstools.Namespaces import WSDL
    
    r = resource
    methods = filter(_issoapmethod, map(lambda i: getattr(r, i), dir(r)))
    tns = ''
    
    #tree = ElementTree()
    defs = Element("{%s}definitions" %WSDL.BASE)
    defs.attrib['name'] = 'SampleDefs'
    defs.attrib['targetNamespace'] = tns
    #tree.append(defs)
    
    porttype = Element("{%s}portType" %WSDL)
    porttype.attrib['name'] = QName("{%s}SamplePortType" %tns)
    
    binding = Element("{%s}binding" %WSDL)
    defs.append(binding)
    binding.attrib['name'] = QName("{%s}SampleBinding" %tns)
    binding.attrib['type'] = porttype.get('name')
    
    for m in methods:
        m.action
        
    service = Element("{%s}service" %WSDL.BASE)
    defs.append(service)
    service.attrib['name'] = 'SampleService'
    
    port = Element("{%s}port" %WSDL.BASE)
    service.append(port)
    port.attrib['name'] = "SamplePort"
    port.attrib['binding'] = binding.get('name')
    
    soapaddress = Element("{%s}address" %WSDL.BIND_SOAP)
    soapaddress.attrib['location'] = 'http://localhost/bla'
    port.append(soapaddress)
    
    return [ElementTree.tostring(defs)]
    
   
    
"""
<?xml version="1.0" encoding="UTF-8"?>
<wsdl:definitions name="Counter" targetNamespace="http://counter.com/bindings" xmlns:wsdl="http://schemas.xmlsoap.org/wsdl/" xmlns:porttype="http://counter.com" xmlns:soap="http://schemas.xmlsoap.org/wsdl/soap/">
  <wsdl:import namespace="http://counter.com" location="counter_flattened.wsdl"/>
  <wsdl:binding name="CounterPortTypeSOAPBinding" type="porttype:CounterPortType">
    <soap:binding style="document" transport="http://schemas.xmlsoap.org/soap/http"/>
    <wsdl:operation name="createCounter">
      <soap:operation soapAction="http://counter.com/CounterPortType/createCounterRequest"/>
      <wsdl:input>
        <soap:body use="literal"/>
      </wsdl:input>
      <wsdl:output>
        <soap:body use="literal"/>
      </wsdl:output>
    </wsdl:operation>


<wsdl:definitions name="Counter" targetNamespace="http://counter.com/service" 
xmlns:wsdl="http://schemas.xmlsoap.org/wsdl/" xmlns:soap="http://schemas.xmlsoap.org/wsdl/soap/" xmlns:binding="http://counter.com/bindings">
  <wsdl:import namespace="http://counter.com/bindings" location="counter_bindings.wsdl"/>
  <wsdl:service name="CounterService">
    <wsdl:port name="CounterPortTypePort" binding="binding:CounterPortTypeSOAPBinding">
      <soap:address location="http://localhost:8080/wsrf/services/"/>
    </wsdl:port>
  </wsdl:service>
</wsdl:definitions>
"""
    
