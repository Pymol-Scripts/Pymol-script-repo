#! /usr/bin/env python
# $Header$
'''ZSI:  Zolera Soap Infrastructure.

Copyright 2001, Zolera Systems, Inc.  All Rights Reserved.
'''

_copyright = """ZSI:  Zolera Soap Infrastructure.

Copyright 2001, Zolera Systems, Inc.  All Rights Reserved.
Copyright 2002-2003, Rich Salz. All Rights Reserved.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, and/or
sell copies of the Software, and to permit persons to whom the Software
is furnished to do so, provided that the above copyright notice(s) and
this permission notice appear in all copies of the Software and that
both the above copyright notice(s) and this permission notice appear in
supporting documentation.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT
OF THIRD PARTY RIGHTS. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR HOLDERS
INCLUDED IN THIS NOTICE BE LIABLE FOR ANY CLAIM, OR ANY SPECIAL INDIRECT
OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
OR PERFORMANCE OF THIS SOFTWARE.

Except as contained in this notice, the name of a copyright holder
shall not be used in advertising or otherwise to promote the sale, use
or other dealings in this Software without prior written authorization
of the copyright holder.


Portions are also:

Copyright (c) 2003, The Regents of the University of California,
through Lawrence Berkeley National Laboratory (subject to receipt of
any required approvals from the U.S. Dept. of Energy). All rights
reserved. Redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following
conditions are met:

(1) Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
(2) Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.
(3) Neither the name of the University of California, Lawrence Berkeley
National Laboratory, U.S. Dept. of Energy nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes,
patches, or upgrades to the features, functionality or performance of
the source code ("Enhancements") to anyone; however, if you choose to
make your Enhancements available either publicly, or directly to
Lawrence Berkeley National Laboratory, without imposing a separate
written license agreement for such Enhancements, then you hereby grant
the following license: a non-exclusive, royalty-free perpetual license
to install, use, modify, prepare derivative works, incorporate into
other computer software, distribute, and sublicense such Enhancements
or derivative works thereof, in binary and source code form.


For wstools also:

Zope Public License (ZPL) Version 2.0
-----------------------------------------------

This software is Copyright (c) Zope Corporation (tm) and
Contributors. All rights reserved.

This license has been certified as open source. It has also
been designated as GPL compatible by the Free Software
Foundation (FSF).

Redistribution and use in source and binary forms, with or
without modification, are permitted provided that the
following conditions are met:

1. Redistributions in source code must retain the above
   copyright notice, this list of conditions, and the following
   disclaimer.

2. Redistributions in binary form must reproduce the above
   copyright notice, this list of conditions, and the following
   disclaimer in the documentation and/or other materials
   provided with the distribution.

3. The name Zope Corporation (tm) must not be used to
   endorse or promote products derived from this software
   without prior written permission from Zope Corporation.

4. The right to distribute this software or to use it for
   any purpose does not give you the right to use Servicemarks
   (sm) or Trademarks (tm) of Zope Corporation. Use of them is
   covered in a separate agreement (see
   http://www.zope.com/Marks).

5. If any files are modified, you must cause the modified
   files to carry prominent notices stating that you changed
   the files and the date of any change.

Disclaimer

  THIS SOFTWARE IS PROVIDED BY ZOPE CORPORATION ``AS IS''
  AND ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT
  NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
  AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN
  NO EVENT SHALL ZOPE CORPORATION OR ITS CONTRIBUTORS BE
  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
  OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
  DAMAGE.


This software consists of contributions made by Zope
Corporation and many individuals on behalf of Zope
Corporation.  Specific attributions are listed in the
accompanying credits file.
"""

##
##  Stuff imported from elsewhere.
from xml.dom import Node as _Node
import types as _types

##
##  Public constants.
from ZSI.wstools.Namespaces import ZSI_SCHEMA_URI
UNICODE_ENCODING = 'utf-8'

##
##  Not public constants.
_inttypes = [ _types.IntType, _types.LongType ]
_floattypes = [ _types.FloatType ]
_seqtypes = [ _types.TupleType, _types.ListType ]
_stringtypes = [ _types.StringType, _types.UnicodeType ]

##
##  Low-level DOM oriented utilities; useful for typecode implementors.
_attrs = lambda E: (E.attributes and E.attributes.values()) or []
_children = lambda E: E.childNodes or []
_child_elements = lambda E: [ n for n in (E.childNodes or [])
                        if n.nodeType == _Node.ELEMENT_NODE ]

##
##  Stuff imported from elsewhere.
from ZSI.wstools.Namespaces import SOAP as _SOAP, SCHEMA as _SCHEMA, XMLNS as _XMLNS

##
##  Low-level DOM oriented utilities; useful for typecode implementors.
_find_arraytype = lambda E: E.getAttributeNS(_SOAP.ENC, "arrayType")
_find_encstyle = lambda E: E.getAttributeNS(_SOAP.ENV, "encodingStyle")
try:
    from xml.dom import EMPTY_NAMESPACE
    _empty_nsuri_list = [ EMPTY_NAMESPACE ]
    #if '' not in _empty_nsuri_list: __empty_nsuri_list.append('')
    #if None not in _empty_nsuri_list: __empty_nsuri_list.append(None)
except:
    _empty_nsuri_list = [ None, '' ]
def _find_attr(E, attr):
    for nsuri in _empty_nsuri_list:
        try:
            v = E.getAttributeNS(nsuri, attr)
            if v: return v
        except: pass
    return None

def _find_attrNS(E, namespaceURI, localName):
    '''namespaceURI
       localName
    '''
    try:
        v = E.getAttributeNS(namespaceURI, localName)
        if v: return v
    except: pass
    return None


def _find_attrNodeNS(E, namespaceURI, localName):
    '''Must grab the attribute Node to distinquish between
    an unspecified attribute(None) and one set to empty string("").
       namespaceURI
       localName
    '''
    attr = E.getAttributeNodeNS(namespaceURI, localName)
    if attr is None: return None
    try:
        return attr.value
    except: pass
    return E.getAttributeNS(namespaceURI, localName)




_find_href = lambda E: _find_attr(E, "href")
_find_xsi_attr = lambda E, attr: \
                E.getAttributeNS(_SCHEMA.XSI3, attr) \
                or E.getAttributeNS(_SCHEMA.XSI1, attr) \
                or E.getAttributeNS(_SCHEMA.XSI2, attr)
_find_type = lambda E: _find_xsi_attr(E, "type")

_find_xmlns_prefix = lambda E, attr: E.getAttributeNS(_XMLNS.BASE, attr)
_find_default_namespace = lambda E: E.getAttributeNS(_XMLNS.BASE, None)

#_textprotect = lambda s: s.replace('&', '&amp;').replace('<', '&lt;')

_get_element_nsuri_name = lambda E: (E.namespaceURI, E.localName)

_is_element = lambda E: E.nodeType == _Node.ELEMENT_NODE

def _resolve_prefix(celt, prefix):
    '''resolve prefix to a namespaceURI.  If None or 
    empty str, return default namespace or None.

    Parameters:
      celt -- element node
      prefix -- xmlns:prefix, or empty str or None
    '''
    namespace = None
    while _is_element(celt):
        if prefix:
            namespaceURI = _find_xmlns_prefix(celt, prefix)
        else:
            namespaceURI = _find_default_namespace(celt)
        if namespaceURI: break
        celt = celt.parentNode
    else:
        if prefix:  
            raise EvaluateException, 'cant resolve xmlns:%s' %prefix
    return namespaceURI

def _valid_encoding(elt):
    '''Does this node have a valid encoding?
    '''
    enc = _find_encstyle(elt)
    if not enc or enc == _SOAP.ENC: return 1
    for e in enc.split():
        if e.startswith(_SOAP.ENC):
            # XXX Is this correct?  Once we find a Sec5 compatible
            # XXX encoding, should we check that all the rest are from
            # XXX that same base?  Perhaps.  But since the if test above
            # XXX will surely get 99% of the cases, leave it for now.
            return 1
    return 0

def _backtrace(elt, dom):
    '''Return a "backtrace" from the given element to the DOM root,
    in XPath syntax.
    '''
    s = ''
    while elt != dom:
        name, parent = elt.nodeName, elt.parentNode
        if parent is None: break
        matches = [ c for c in _child_elements(parent) 
                        if c.nodeName == name ]
        if len(matches) == 1:
            s = '/' + name + s
        else:
            i = matches.index(elt) + 1
            s = ('/%s[%d]' % (name, i)) + s
        elt = parent
    return s

def _get_idstr(pyobj):
    '''Python 2.3.x generates a FutureWarning for negative IDs, so
    we use a different prefix character to ensure uniqueness, and
    call abs() to avoid the warning.'''
    x = id(pyobj)
    if x < 0:
        return 'x%x' % abs(x)
    return 'o%x' % x


def _get_postvalue_from_absoluteURI(url):
    """Bug [ 1513000 ] POST Request-URI not limited to "abs_path"
    Request-URI = "*" | absoluteURI | abs_path | authority
    
    Not a complete solution, but it seems to work with all known
    implementations.  ValueError thrown if bad uri.
    """
    cache = _get_postvalue_from_absoluteURI.cache
    path = cache.get(url, '')
    if not path:
        scheme,authpath = url.split('://')
        s = authpath.split('/', 1)
        if len(s) == 2: path = '/%s' %s[1]
        if len(cache) > _get_postvalue_from_absoluteURI.MAXLEN:cache.clear()
        cache[url] = path
    return path
_get_postvalue_from_absoluteURI.cache = {}
_get_postvalue_from_absoluteURI.MAXLEN = 20


##
##  Exception classes.
class ZSIException(Exception):
    '''Base class for all ZSI exceptions.
    '''
    pass

class ParseException(ZSIException):
    '''Exception raised during parsing.
    '''

    def __init__(self, str, inheader, elt=None, dom=None):
        Exception.__init__(self)
        self.str, self.inheader, self.trace = str, inheader, None
        if elt and dom:
            self.trace = _backtrace(elt, dom)

    def __str__(self):
        if self.trace:
            return self.str + '\n[Element trace: ' + self.trace + ']'
        return self.str

    def __repr__(self):
        return "<%s.ParseException %s>" % (__name__, _get_idstr(self))


class EvaluateException(ZSIException):
    '''Exception raised during data evaluation (serialization).
    '''

    def __init__(self, str, trace=None):
        Exception.__init__(self)
        self.str, self.trace = str, trace

    def __str__(self):
        if self.trace:
            return self.str + '\n[Element trace: ' + self.trace + ']'
        return self.str

    def __repr__(self):
        return "<%s.EvaluateException %s>" % (__name__, _get_idstr(self))

class FaultException(ZSIException):
    '''Exception raised when a fault is received.
    '''

    def __init__(self, fault):
        self.fault = fault

    def __str__(self):
        return str(self.fault)

    def __repr__(self):
        return "<%s.FaultException %s>" % (__name__, _get_idstr(self))

class WSActionException(ZSIException):
    '''Exception raised when WS-Address Action Header is incorrectly
    specified when received by client or server.
    '''
    pass

##
##  Importing the rest of ZSI.
import version
def Version():
    return version.Version

from writer import SoapWriter
from parse import ParsedSoap
from fault import Fault, \
    FaultFromActor, FaultFromException, FaultFromFaultMessage, \
    FaultFromNotUnderstood, FaultFromZSIException
import TC
TC.RegisterType(TC.String, minOccurs=0, nillable=False)
TC.RegisterType(TC.URI, minOccurs=0, nillable=False)
TC.RegisterType(TC.Base64String, minOccurs=0, nillable=False)
TC.RegisterType(TC.HexBinaryString, minOccurs=0, nillable=False)

#TC.RegisterType(TC.Integer)
#TC.RegisterType(TC.Decimal)
for pyclass in (TC.IunsignedByte, TC.IunsignedShort, TC.IunsignedInt, TC.IunsignedLong,
    TC.Ibyte, TC.Ishort, TC.Iint, TC.Ilong, TC.InegativeInteger,
    TC.InonPositiveInteger, TC.InonNegativeInteger, TC.IpositiveInteger,
    TC.Iinteger, TC.FPfloat, TC.FPdouble, ):

    TC.RegisterType(pyclass, minOccurs=0, nillable=False)

TC.RegisterType(TC.Boolean, minOccurs=0, nillable=False)
TC.RegisterType(TC.Duration, minOccurs=0, nillable=False)
TC.RegisterType(TC.gDateTime, minOccurs=0, nillable=False)
TC.RegisterType(TC.gDate, minOccurs=0, nillable=False)
TC.RegisterType(TC.gYearMonth, minOccurs=0, nillable=False)
TC.RegisterType(TC.gYear, minOccurs=0, nillable=False)
TC.RegisterType(TC.gMonthDay, minOccurs=0, nillable=False)
TC.RegisterType(TC.gDay, minOccurs=0, nillable=False)
TC.RegisterType(TC.gTime, minOccurs=0, nillable=False)
TC.RegisterType(TC.Apache.Map, minOccurs=0, nillable=False)

##
## Register Wrappers for builtin types.
## TC.AnyElement wraps builtins so element name information can be saved
##
import schema
for i in [int,float,str,tuple,list,unicode]:
    schema._GetPyobjWrapper.RegisterBuiltin(i)

## Load up Wrappers for builtin types
schema.RegisterAnyElement()


#try:
#    from ServiceProxy import *
#except:
#    pass

if __name__ == '__main__': print _copyright
