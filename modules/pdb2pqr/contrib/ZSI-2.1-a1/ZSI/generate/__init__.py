# $Id: __init__.py 1261 2006-10-10 23:45:43Z boverhof $

__all__ = ['wsdl2python', 'utility', 'containers', 'commands']

class WSDLFormatError(Exception):
   '''Detected errors in the WSDL document.
   '''

class WsdlGeneratorError(Exception):
    pass

class Wsdl2PythonError(Exception):
    pass

class WSInteropError(Exception):
   '''Conformance to WS-I Basic-Profile 1.0 specification
   '''

class WSISpec:
    R2203 = 'An rpc-literal binding in a DESCRIPTION MUST refer, in its soapbind:body element(s), only to wsdl:part element(s) that have been defined using the type attribute.'
    R2710 = 'The operations in a wsdl:binding in a DESCRIPTION MUST result in wire signatures that are different from one another.'
    R2717 = 'An rpc-literal binding in a DESCRIPTION MUST have the namespace attribute specified, the value of which MUST be an absolute URI, on contained  soapbind:body elements.'
    R2729 = 'A MESSAGE described with an rpc-literal binding that is a response message MUST have a wrapper element whose name is the corresponding wsdl:operation name suffixed with the string "Response"'
    
