###########################################################################
# Joshua R. Boverhof, LBNL
# See Copyright for copyright notice!
# $Id: $
###########################################################################
from ZSI import _get_element_nsuri_name, SoapWriter, ParsedSoap
from interfaces import HandlerChainInterface
from zope.interface import classProvides, implements, Interface


class DataHandler:
    """ str --> ps, sw --> str
    class variables:
        readerClass -- factory class to create reader for ParsedSoap instances.
    """
    classProvides(HandlerChainInterface)
    readerClass = None

    @classmethod
    def processRequest(cls, input, **kw):
        return ParsedSoap(input, readerclass=cls.readerClass)

    @classmethod
    def processResponse(cls, sw, **kw):
        return str(sw)


class CallbackHandler:
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
        method =  getattr(resource, 'soap_%s' %
                           _get_element_nsuri_name(ps.body_root)[-1])
                                              
        try:
            req,rsp = method(ps, request=request)
        except Exception, ex:
            raise
        
        return rsp
    
    @classmethod
    def processResponse(cls, output, **kw):
        sw = SoapWriter(outputclass=cls.writerClass)
        sw.serialize(output)
        return sw
    

class ReverseHandlerChain:

    def __init__(self, *handlers):
        self.in_handlers = handlers
        handlers = list(handlers); handlers.reverse()
        self.out_handlers = tuple(handlers)
        
    def processRequest(self, arg, **kw):
        for h in self.in_handlers:
            arg = h.processRequest(arg, **kw)
            
        return arg
            
    def processResponse(self, arg, **kw):
        if arg is None: 
            return

        for h in self.out_handlers:
            arg = h.processResponse(arg, **kw)
            
        return arg


class ReverseHandlerChainFactory:
    protocol = ReverseHandlerChain
    
    @classmethod
    def newInstance(cls):
        return cls.protocol(DataHandler, CallbackHandler)

