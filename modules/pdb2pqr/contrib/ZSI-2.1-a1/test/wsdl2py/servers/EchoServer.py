#!/usr/bin/env python
############################################################################
# Joshua R. Boverhof, LBNL
# See LBNLCopyright for copyright notice!
###########################################################################
import sys
from ZSI import ServiceContainer, Fault
from ZSI.ServiceContainer import AsServer, ServiceSOAPBinding
from EchoServer_server import EchoServer

from ZSI.fault import Fault, ZSIFaultDetail
def WSDLFaultFromException(ex, inheader, tb=None, actor=None):
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

    if isinstance(ex, Fault):
        return ex

    elt = ZSIFaultDetail(string=exceptionName + "\n" + str(ex), trace=tracetext)
    if inheader:
        detail, headerdetail = None, elt
    else:
        detail, headerdetail = elt, None
    return Fault(Fault.Server, 'Processing Failure',
                actor, detail, headerdetail)

ServiceContainer.FaultFromException = WSDLFaultFromException

class Service(EchoServer):
    def soap_Echo(self, ps, **kw):
        request,response = EchoServer.soap_Echo(self, ps, **kw)
        response.EchoResult = request.EchoIn
        #return request,response
        #raise RuntimeError, 'hi'
        #raise Fault(911, "EMERGENCY", None, response)
        return request,response


def twisted_main(port):
    from twisted.internet import reactor
    from twisted.application import service, internet
    from twisted.web.resource import Resource
    from twisted.web.server import Site

    root = Resource()
    root.putChild('test', Service())
    reactor.listenTCP(port, Site(root))

def main():
    port = int(sys.argv[1])
    if issubclass(EchoServer, ServiceSOAPBinding):
        AsServer(port, (Service('test'),))
        return

    #from ZSI.twisted.WSresource import WSResource
    #if issubclass(EchoServer, WSResource):
    from twisted.internet import reactor
    reactor.callWhenRunning(twisted_main, port)
    reactor.run()


if __name__ == "__main__" :
    main()
