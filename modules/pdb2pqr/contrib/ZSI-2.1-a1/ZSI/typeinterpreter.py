###########################################################################
# Joshua R. Boverhof, LBNL
# See LBNLCopyright for copyright notice!
###########################################################################

import ZSI
from ZSI import TC, TCtimes, TCcompound
from ZSI.TC import TypeCode
from ZSI import _copyright, EvaluateException 
from ZSI.wstools.Utility import SplitQName
from ZSI.wstools.Namespaces import SOAP, SCHEMA

###########################################################################
# Module Classes: BaseTypeInterpreter
###########################################################################

class NamespaceException(Exception): pass
class BaseTypeInterpreter:
    """Example mapping of xsd/soapenc types to zsi python types.
    Checks against all available classes in ZSI.TC.  Used in 
    wsdl2python, wsdlInterpreter, and ServiceProxy.
    """

    def __init__(self):
        self._type_list = [TC.Iinteger, TC.IunsignedShort, TC.gYearMonth, \
                           TC.InonNegativeInteger, TC.Iint, TC.String, \
                           TC.gDateTime, TC.IunsignedInt, TC.Duration,\
                           TC.IpositiveInteger, TC.FPfloat, TC.gDay, TC.gMonth, \
                           TC.InegativeInteger, TC.gDate, TC.URI, \
                           TC.HexBinaryString, TC.IunsignedByte, \
                           TC.gMonthDay, TC.InonPositiveInteger, \
                           TC.Ibyte, TC.FPdouble, TC.gTime, TC.gYear, \
                           TC.Ilong, TC.IunsignedLong, TC.Ishort, \
                           TC.Token, TC.QName]

        self._tc_to_int = [
            ZSI.TCnumbers.IEnumeration,
            ZSI.TCnumbers.Iint,
            ZSI.TCnumbers.Iinteger,
            ZSI.TCnumbers.Ilong,
            ZSI.TCnumbers.InegativeInteger,
            ZSI.TCnumbers.InonNegativeInteger,
            ZSI.TCnumbers.InonPositiveInteger,
            ZSI.TC.Integer,
            ZSI.TCnumbers.IpositiveInteger,
            ZSI.TCnumbers.Ishort]
  
        self._tc_to_float = [
            ZSI.TC.Decimal,
            ZSI.TCnumbers.FPEnumeration,
            ZSI.TCnumbers.FPdouble,
            ZSI.TCnumbers.FPfloat]
        
        self._tc_to_string = [
            ZSI.TC.Base64String,
            ZSI.TC.Enumeration,
            ZSI.TC.HexBinaryString,
            ZSI.TCnumbers.Ibyte,
            ZSI.TCnumbers.IunsignedByte,
            ZSI.TCnumbers.IunsignedInt,
            ZSI.TCnumbers.IunsignedLong,
            ZSI.TCnumbers.IunsignedShort,
            ZSI.TC.String,
            ZSI.TC.URI,
            ZSI.TC.XMLString,
            ZSI.TC.Token]

        self._tc_to_tuple = [
            ZSI.TC.Duration,
            ZSI.TC.QName,
            ZSI.TCtimes.gDate,
            ZSI.TCtimes.gDateTime,
            ZSI.TCtimes.gDay,
            ZSI.TCtimes.gMonthDay,
            ZSI.TCtimes.gTime,
            ZSI.TCtimes.gYear,
            ZSI.TCtimes.gMonth,
            ZSI.TCtimes.gYearMonth]
        
        return
    
    def _get_xsd_typecode(self, msg_type):
        untaged_xsd_types = {'boolean':TC.Boolean, 
            'decimal':TC.Decimal, 
            'base64Binary':TC.Base64String}
        if untaged_xsd_types.has_key(msg_type):
            return untaged_xsd_types[msg_type]
        for tc in self._type_list:
            if tc.type == (SCHEMA.XSD3,msg_type):
                break
        else:
            tc = TC.AnyType
        return tc

    def _get_soapenc_typecode(self, msg_type):
        if msg_type == 'Array':
            return TCcompound.Array
        if msg_type == 'Struct':
            return TCcompound.Struct

        return self._get_xsd_typecode(msg_type)

    def get_typeclass(self, msg_type, targetNamespace):
        prefix, name = SplitQName(msg_type)
        if targetNamespace in SCHEMA.XSD_LIST:
            return self._get_xsd_typecode(name)
        elif targetNamespace in [SOAP.ENC]:
            return self._get_soapenc_typecode(name)
        return None

    def get_pythontype(self, msg_type, targetNamespace, typeclass=None):
        if not typeclass:
            tc = self.get_typeclass(msg_type, targetNamespace)
        else:
            tc = typeclass
        if tc in self._tc_to_int:
            return 'int'
        elif tc in self._tc_to_float:
            return 'float'
        elif tc in self._tc_to_string:
            return 'str'
        elif tc in self._tc_to_tuple:
            return 'tuple'
        elif tc in [TCcompound.Array]:
            return 'list'
        elif tc in [TC.Boolean]:
            return 'bool'
        elif isinstance(tc, TypeCode):
            raise EvaluateException,\
               'failed to map zsi typecode to a python type'
        return None


