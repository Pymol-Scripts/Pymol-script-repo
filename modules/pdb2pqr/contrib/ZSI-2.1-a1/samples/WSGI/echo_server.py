############################################################################
# Joshua R. Boverhof, LBNL
# See Copyright for copyright notice!
# $Id: __init__.py 1132 2006-02-17 01:55:41Z boverhof $
###########################################################################
import sys
from EchoServer_client import *
from ZSI.twisted.wsgi import SOAPApplication, soapmethod, SOAPHandlerChainFactory

class EchoService(SOAPApplication):
    factory = SOAPHandlerChainFactory
    wsdl_content = dict(name='Echo', targetNamespace='urn:echo', imports=(), portType='')

    @soapmethod(EchoRequest.typecode, EchoResponse.typecode, operation='Echo', soapaction='Echo')
    def soap_Echo(self, request, response, **kw):
        response = request
        return request,response


def main():
    from wsgiref.simple_server import make_server, demo_app
    from ZSI.twisted.wsgi import WSGIApplication
    application = WSGIApplication()
    httpd = make_server('', 8000, application)
    application['echo'] = EchoService()
    httpd.serve_forever()

def main_twisted():
    from ZSI.twisted.wsgi import test, WSGIApplication
    app = WSGIApplication()
    app['echo'] = EchoService()
    test(app, port=8000)

if __name__ == '__main__':
    if len(sys.argv) == 1:
        main()
    else: 
        var = sys.argv[1] 
        try:
            getattr(sys.modules[__name__], 'main_%s' %var)(*sys.argv[2:])
        except Exception, ex:
            print>>sys.stderr, ex
            sys.exit(1)
