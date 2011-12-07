############################################################################
# Joshua R. Boverhof, LBNL
# See LBNLCopyright for copyright notice!
###########################################################################
import time

# twisted & related imports
from zope.interface import classProvides, implements, Interface
from twisted.web import client
from twisted.internet import defer
from twisted.internet import reactor
from twisted.python import log
from twisted.python.failure import Failure

from ZSI.parse import ParsedSoap
from ZSI.writer import SoapWriter
from ZSI.fault import FaultFromFaultMessage
from ZSI.wstools.Namespaces import WSA

from WSresource import HandlerChainInterface, CheckInputArgs


# 
# Stability: Unstable
# 

class HTTPPageGetter(client.HTTPPageGetter):
    def handleStatus_500(self):
        """potentially a SOAP:Fault.
        """
        log.err('HTTP Error 500')
    def handleStatus_404(self):
        """client error, not found
        """
        log.err('HTTP Error 404')


client.HTTPClientFactory.protocol = HTTPPageGetter


def getPage(url, contextFactory=None, *args, **kwargs):
    """Download a web page as a string.

    Download a page. Return a deferred, which will callback with a
    page (as a string) or errback with a description of the error.

    See HTTPClientFactory to see what extra args can be passed.
    """
    scheme, host, port, path = client._parse(url)
    factory = client.HTTPClientFactory(url, *args, **kwargs)
    if scheme == 'https':
        if contextFactory is None:
            raise RuntimeError, 'must provide a contextFactory'
        conn = reactor.connectSSL(host, port, factory, contextFactory)
    else:
        conn = reactor.connectTCP(host, port, factory)

    return factory


class ClientDataHandler:
    """
    class variables:
        readerClass -- factory class to create reader for ParsedSoap instances.
        writerClass -- ElementProxy implementation to use for SoapWriter 
            instances.
    """
    classProvides(HandlerChainInterface)
    readerClass = None
    writerClass = None

    @classmethod
    def processResponse(cls, soapdata, **kw):
        """called by deferred, returns pyobj representing reply.
        Parameters and Key Words:
          soapdata -- SOAP Data
          replytype -- reply type of response
        """
        if len(soapdata) == 0:
            raise TypeError('Received empty response')

#        log.msg("_" * 33, time.ctime(time.time()), 
#                  "RESPONSE: \n%s" %soapdata, debug=True)

        ps = ParsedSoap(soapdata, readerclass=cls.readerClass)
        if ps.IsAFault() is True:
            log.msg('Received SOAP:Fault', debug=True)
            raise FaultFromFaultMessage(ps)

        return ps

    @classmethod
    def processRequest(cls, obj, nsdict={}, header=True, 
                       **kw):
        tc = None
        if kw.has_key('requesttypecode'):
            tc = kw['requesttypecode']
        elif kw.has_key('requestclass'):
            tc = kw['requestclass'].typecode
        else:
            tc = getattr(obj.__class__, 'typecode', None)

        sw = SoapWriter(nsdict=nsdict, header=header,
                        outputclass=cls.writerClass)
        sw.serialize(obj, tc)
        return sw
    
    
class WSAddressHandler:
    """Minimal WS-Address handler.  Most of the logic is in
    the ZSI.address.Address class.
    
    class variables:
        uri -- default WSA Addressing URI
    """
    implements(HandlerChainInterface)
    uri = WSA.ADDRESS
    
    def processResponse(self, ps, wsaction=None, soapaction=None, **kw):
        addr = self.address 
        addr.parse(ps)
        action = addr.getAction()
        if not action:
            raise WSActionException('No WS-Action specified in Request')

        if not soapaction:
            return ps
    
        soapaction = soapaction.strip('\'"')
        if soapaction and soapaction != wsaction:
            raise WSActionException(\
                'SOAP Action("%s") must match WS-Action("%s") if specified.'%(
                    soapaction, wsaction)
            )

        return ps

    def processRequest(self, sw, wsaction=None, url=None, endPointReference=None, **kw):
        from ZSI.address import Address
        if sw is None:
            self.address = None
            return
        
        if not sw.header:
            raise RuntimeError, 'expecting SOAP:Header'
        
        self.address = addr = Address(url, wsAddressURI=self.uri)
        addr.setRequest(endPointReference, wsaction)
        addr.serialize(sw, typed=False)
        
        return sw


class DefaultClientHandlerChain:

    @CheckInputArgs(HandlerChainInterface)
    def __init__(self, *handlers):
        self.handlers = handlers
        self.debug = len(log.theLogPublisher.observers) > 0
        self.flow = None
        
    @staticmethod
    def parseResponse(ps, replytype):
        return ps.Parse(replytype)
        
    def processResponse(self, arg, replytype, **kw):
        """
        Parameters:
            arg -- deferred 
            replytype -- typecode
        """
        if self.debug:
            log.msg('--->PROCESS REQUEST\n%s' %arg, debug=1)

        for h in self.handlers:
            arg.addCallback(h.processResponse, **kw)
            
        arg.addCallback(self.parseResponse, replytype)
            
    def processRequest(self, arg, **kw):
        """
        Parameters:
            arg -- XML Soap data string
        """
        if self.debug:
            log.msg('===>PROCESS RESPONSE: %s' %str(arg), debug=1)

        if arg is None:
            return

        for h in self.handlers:
            arg = h.processRequest(arg, **kw)
            
        s = str(arg)
        if self.debug:
            log.msg(s, debug=1)

        return s


class DefaultClientHandlerChainFactory:
    protocol = DefaultClientHandlerChain
    
    @classmethod
    def newInstance(cls):
        return cls.protocol(ClientDataHandler)
        

class WSAddressClientHandlerChainFactory:
    protocol = DefaultClientHandlerChain
    
    @classmethod
    def newInstance(cls):
        return cls.protocol(ClientDataHandler, 
            WSAddressHandler())


class Binding:
    """Object that represents a binding (connection) to a SOAP server.
    """
    agent='ZSI.twisted client'
    factory = DefaultClientHandlerChainFactory
    defer = False

    def __init__(self, url=None, nsdict=None, contextFactory=None, 
                 tracefile=None, **kw):
        """Initialize.
        Keyword arguments include:
            url -- URL of resource, POST is path 
            nsdict -- namespace entries to add
            contextFactory -- security contexts
            tracefile -- file to dump packet traces
        """
        self.url = url
        self.nsdict = nsdict or {}
        self.contextFactory = contextFactory
        self.http_headers  = {'content-type': 'text/xml',}
        self.trace = tracefile

    def addHTTPHeader(self, key, value):
        self.http_headers[key] = value
   
    def getHTTPHeaders(self):
        return self.http_headers
        
    def Send(self, url, opname, pyobj, nsdict={}, soapaction=None, chain=None, 
             **kw):
        """Returns a ProcessingChain which needs to be passed to Receive if 
        Send is being called consecutively.
        """
        url = url or self.url
        cookies = None
        if chain is not None:
            cookies = chain.flow.cookies
        
        d = {}
        d.update(self.nsdict)
        d.update(nsdict)
         
        if soapaction is not None:
            self.addHTTPHeader('SOAPAction', soapaction)
        
        chain = self.factory.newInstance()
        soapdata = chain.processRequest(pyobj, nsdict=nsdict, 
                                        soapaction=soapaction, **kw)
            
        if self.trace:
            print >>self.trace, "_" * 33, time.ctime(time.time()), "REQUEST:"
            print >>self.trace, soapdata

        f = getPage(str(url), contextFactory=self.contextFactory, 
                    postdata=soapdata, agent=self.agent, 
                    method='POST', headers=self.getHTTPHeaders(), 
                    cookies=cookies)
        
        if isinstance(f, Failure):
            return f
        
        chain.flow = f
        self.chain = chain
        return chain
            
    def Receive(self, replytype, chain=None, **kw):
        """This method allows code to act in a synchronous manner, it waits to 
        return until the deferred fires but it doesn't prevent other queued 
        calls from being executed.  Send must be called first, which sets up 
        the chain/factory.  
        
        WARNING: If defer is set to True, must either call Receive
        immediately after Send (ie. no intervening Sends) or pass
        chain in as a paramter.
        
        Parameters:
            replytype -- TypeCode
        KeyWord Parameters:
            chain -- processing chain, optional
            
        """        
        chain = chain or self.chain
        d = chain.flow.deferred
        if self.trace:
            def trace(soapdata):
                print >>self.trace, "_" * 33, time.ctime(time.time()), "RESPONSE:"
                print >>self.trace, soapdata
                return soapdata
            
            d.addCallback(trace)
            
        chain.processResponse(d, replytype, **kw)
        if self.defer:
            return d

        failure = []
        append = failure.append
        def errback(result):
            """Used with Response method to suppress 'Unhandled error in 
            Deferred' messages by adding an errback.
            """
            append(result)
            return None
        
        d.addErrback(errback)

        # spin reactor
        while not d.called:
            reactor.runUntilCurrent()
            t2 = reactor.timeout()
            t = reactor.running and t2
            reactor.doIteration(t)

        pyobj = d.result
        if len(failure):
            failure[0].raiseException()
        
        return pyobj

def trace():
        if trace:
            print >>trace, "_" * 33, time.ctime(time.time()), "RESPONSE:"
            for i in (self.reply_code, self.reply_msg,):
                print >>trace, str(i)
            print >>trace, "-------"
            print >>trace, str(self.reply_headers)
            print >>trace, self.data
