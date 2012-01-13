###########################################################################
# Joshua R. Boverhof, LBNL
# See Copyright for copyright notice!
# $Id: $
###########################################################################

import sys, warnings

# twisted & related imports
from zope.interface import classProvides, implements, Interface

# ZSI imports
from ZSI import EvaluateException, ParseException, ParsedSoap, SoapWriter


# 
# Stability: Unstable
# 

def CheckInputArgs(*interfaces):
    """Must provide at least one interface, the last one may be repeated.
    """
    l = len(interfaces)
    def wrapper(func):
        def check_args(self, *args, **kw):
            for i in range(len(args)):
                if (l > i and interfaces[i].providedBy(args[i])) or interfaces[-1].providedBy(args[i]):
                    continue
                if l > i: raise TypeError, 'arg %s does not implement %s' %(args[i], interfaces[i])
                raise TypeError, 'arg %s does not implement %s' %(args[i], interfaces[-1])
            func(self, *args, **kw)
        return check_args
    return wrapper


class HandlerChainInterface(Interface):
    
    def processRequest(self, input, **kw):
        """returns a representation of the request, the 
        last link in the chain must return a response
        pyobj with a typecode attribute.
        Parameters:
            input --
        Keyword Parameters:
            request -- HTTPRequest instance
            resource  -- Resource instance
        """
    def processResponse(self, output, **kw):
        """returns a string representing the soap response.
        Parameters
            output --
        Keyword Parameters:
            request -- HTTPRequest instance
            resource  -- Resource instance
        """

class CallbackChainInterface(Interface):
    
    def processRequest(self, input, **kw):
        """returns a response pyobj with a typecode 
        attribute.
        Parameters:
            input --
        Keyword Parameters:
            request -- HTTPRequest instance
            resource  -- Resource instance
        """

class DataHandler:
    """
    class variables:
        readerClass -- factory class to create reader for ParsedSoap instances.
        writerClass -- ElementProxy implementation to use for SoapWriter instances.
    """
    classProvides(HandlerChainInterface)
    readerClass = None
    writerClass = None

    @classmethod
    def processRequest(cls, input, **kw):
        return ParsedSoap(input, readerclass=cls.readerClass)

    @classmethod
    def processResponse(cls, output, **kw):
        sw = SoapWriter(outputclass=cls.writerClass)
        sw.serialize(output)
        return sw


class DefaultHandlerChain:

    @CheckInputArgs(CallbackChainInterface, HandlerChainInterface)
    def __init__(self, cb, *handlers):
        self.handlercb = cb
        self.handlers = handlers
        
    def processRequest(self, arg, **kw):
        
        for h in self.handlers:
            arg = h.processRequest(arg, **kw)
            
        return self.handlercb.processRequest(arg, **kw)
            
    def processResponse(self, arg, **kw):

        if arg is None: 
            return

        for h in self.handlers:
            arg = h.processResponse(arg, **kw)
            
        s = str(arg)
        
        return s

