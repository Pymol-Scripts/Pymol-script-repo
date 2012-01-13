#! /usr/bin/env python
# $Header$
'''General typecodes.
'''

from ZSI import _copyright, _children, _child_elements, \
    _floattypes, _stringtypes, _seqtypes, _find_attr, _find_attrNS, _find_attrNodeNS, \
    _find_arraytype, _find_default_namespace, _find_href, _find_encstyle, \
    _resolve_prefix, _find_xsi_attr, _find_type, \
    _find_xmlns_prefix, _get_element_nsuri_name, _get_idstr, \
    _Node, EvaluateException, UNICODE_ENCODING, \
    _valid_encoding, ParseException
    
from ZSI.wstools.Namespaces import SCHEMA, SOAP
from ZSI.wstools.Utility import SplitQName
from ZSI.wstools.c14n import Canonicalize
from ZSI.wstools.logging import getLogger as _GetLogger

import re, types, time, copy

from base64 import decodestring as b64decode, encodestring as b64encode
from urllib import unquote as urldecode, quote as urlencode
from binascii import unhexlify as hexdecode, hexlify as hexencode
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO


_is_xsd_or_soap_ns = lambda ns: ns in [
                        SCHEMA.XSD3, SOAP.ENC, SCHEMA.XSD1, SCHEMA.XSD2, ]
_find_nil = lambda E: _find_xsi_attr(E, "null") or _find_xsi_attr(E, "nil")

def _get_xsitype(pyclass):
    '''returns the xsi:type as a tuple, coupled with ZSI.schema
    '''
    if hasattr(pyclass,'type') and type(pyclass.type) in _seqtypes:
        return pyclass.type
    elif hasattr(pyclass,'type') and hasattr(pyclass, 'schema'):
        return (pyclass.schema, pyclass.type)

    return (None,None)


# value returned when xsi:nil="true"
Nilled = None
UNBOUNDED = 'unbounded'


class TypeCode:
    '''The parent class for all parseable SOAP types.
    Class data:
        typechecks -- do init-time type checking if non-zero
    Class data subclasses may define:
        tag -- global element declaration
        type -- global type definition
        parselist -- list of valid SOAP types for this class, as
            (uri,name) tuples, where a uri of None means "all the XML
            Schema namespaces"
        errorlist -- parselist in a human-readable form; will be
            generated if/when needed
        seriallist -- list of Python types or user-defined classes
            that this typecode can serialize.
        logger -- logger instance for this class.
    '''
    tag = None
    type = (None,None)
    typechecks = True
    attribute_typecode_dict = None
    logger = _GetLogger('ZSI.TC.TypeCode')

    def __init__(self, pname=None, aname=None, minOccurs=1,
         maxOccurs=1, nillable=False, typed=True, unique=True, 
         pyclass=None, attrs_aname='_attrs', **kw):
        '''Baseclass initialization.
        Instance data (and usually keyword arg)
            pname -- the parameter name (localname).
            nspname -- the namespace for the parameter;
                None to ignore the namespace
            typed -- output xsi:type attribute
            unique -- data item is not aliased (no href/id needed)
            minOccurs -- min occurances
            maxOccurs -- max occurances
            nillable -- is item nillable?
            attrs_aname -- This is variable name to dictionary of attributes
            encoded -- encoded namespaceURI (specify if use is encoded)
        '''
        if type(pname) in _seqtypes:
            self.nspname, self.pname = pname
        else:
            self.nspname, self.pname = None, pname

        if self.pname:
            self.pname = str(self.pname).split(':')[-1]

        self.aname = aname or self.pname
        self.minOccurs = minOccurs
        self.maxOccurs = maxOccurs
        self.nillable = nillable
        self.typed = typed
        self.unique = unique
        self.attrs_aname = attrs_aname
        self.pyclass = pyclass

        # Need this stuff for rpc/encoded.
        encoded = kw.get('encoded')
        if encoded is not None:
            self.nspname = kw['encoded']

    def parse(self, elt, ps):
        '''
        Parameters:
            elt -- the DOM element being parsed
            ps -- the ParsedSoap object.
        '''
        raise EvaluateException("Unimplemented evaluation", ps.Backtrace(elt))

    def serialize(self, elt, sw, pyobj, name=None, orig=None, **kw):
        '''
        Parameters:
           elt -- the current DOMWrapper element 
           sw -- soapWriter object
           pyobj -- python object to serialize

        '''
        raise EvaluateException("Unimplemented evaluation", sw.Backtrace(elt))

    def text_to_data(self, text, elt, ps):
        '''convert text into typecode specific data.
        Parameters:
            text -- text content
            elt -- the DOM element being parsed
            ps -- the ParsedSoap object.
        '''
        raise EvaluateException("Unimplemented evaluation", ps.Backtrace(elt))

    def serialize_as_nil(self, elt):
        '''
        Parameters:
            elt -- the current DOMWrapper element 
        '''
        elt.setAttributeNS(SCHEMA.XSI3, 'nil', '1')

    def SimpleHREF(self, elt, ps, tag):
        '''Simple HREF for non-string and non-struct and non-array.
        Parameters:
            elt -- the DOM element being parsed
            ps -- the ParsedSoap object.
            tag -- 
        '''
        if len(_children(elt)): return elt
        href = _find_href(elt)
        if not href:
            if self.minOccurs is 0: return None
            raise EvaluateException('Required' + tag + ' missing',
                    ps.Backtrace(elt))
        return ps.FindLocalHREF(href, elt, 0)

    def get_parse_and_errorlist(self):
        """Get the parselist and human-readable version, errorlist is returned,
        because it is used in error messages.
        """
        d = self.__class__.__dict__
        parselist = d.get('parselist')
        errorlist = d.get('errorlist')
        if parselist and not errorlist:
            errorlist = []
            for t in parselist:
                if t[1] not in errorlist: errorlist.append(t[1])
            errorlist = ' or '.join(errorlist)
            d['errorlist'] = errorlist
        return (parselist, errorlist)

    def checkname(self, elt, ps):
        '''See if the name and type of the "elt" element is what we're
        looking for.   Return the element's type.
        Parameters:
            elt -- the DOM element being parsed
            ps -- the ParsedSoap object.
        '''

        parselist,errorlist = self.get_parse_and_errorlist()
        ns, name = _get_element_nsuri_name(elt)
        if ns == SOAP.ENC:
            # Element is in SOAP namespace, so the name is a type.
            if parselist and \
            (None, name) not in parselist and (ns, name) not in parselist:
                raise EvaluateException(
                'Element mismatch (got %s wanted %s) (SOAP encoding namespace)' % \
                        (name, errorlist), ps.Backtrace(elt))
            return (ns, name)

        # Not a type, check name matches.
        if self.nspname and ns != self.nspname:
            raise EvaluateException('Element NS mismatch (got %s wanted %s)' % \
                (ns, self.nspname), ps.Backtrace(elt))

        if self.pname and name != self.pname:
            raise EvaluateException('Element Name mismatch (got %s wanted %s)' % \
                (name, self.pname), ps.Backtrace(elt))
        return self.checktype(elt, ps)

    def checktype(self, elt, ps):
        '''See if the type of the "elt" element is what we're looking for.
        Return the element's type.
        Parameters:
            elt -- the DOM element being parsed
            ps -- the ParsedSoap object.
        '''
        typeName = _find_type(elt)
        if typeName is None or typeName == "":
            return (None,None)

        # Parse the QNAME.
        prefix,typeName = SplitQName(typeName)
        uri = ps.GetElementNSdict(elt).get(prefix)
        if uri is None:
            raise EvaluateException('Malformed type attribute (bad NS)',
                    ps.Backtrace(elt))

        #typeName = list[1]
        parselist,errorlist = self.get_parse_and_errorlist()
        if not parselist or \
        (uri,typeName) in parselist or \
        (_is_xsd_or_soap_ns(uri) and (None,typeName) in parselist):
            return (uri,typeName)
        raise EvaluateException(
                'Type mismatch (%s namespace) (got %s wanted %s)' % \
                (uri, typeName, errorlist), ps.Backtrace(elt))

    def name_match(self, elt):
        '''Simple boolean test to see if we match the element name.
        Parameters:
            elt -- the DOM element being parsed
        '''
        return self.pname == elt.localName and \
                    self.nspname in [None, '', elt.namespaceURI]

    def nilled(self, elt, ps):
        '''Is the element NIL, and is that okay?
        Parameters:
            elt -- the DOM element being parsed
            ps -- the ParsedSoap object.
        '''
        if _find_nil(elt) not in [ "true",  "1"]: return False
        if self.nillable is False:
            raise EvaluateException('Non-nillable element is NIL',
                    ps.Backtrace(elt))
        return True

    def simple_value(self, elt, ps, mixed=False):
        '''Get the value of the simple content of this element.
        Parameters:
            elt -- the DOM element being parsed
            ps -- the ParsedSoap object.
            mixed -- ignore element content, optional text node
        '''
        if not _valid_encoding(elt):
            raise EvaluateException('Invalid encoding', ps.Backtrace(elt))
        c = _children(elt)
        if mixed is False:
            if len(c) == 0:
                raise EvaluateException('Value missing', ps.Backtrace(elt))
            for c_elt in c:
                if c_elt.nodeType == _Node.ELEMENT_NODE:
                    raise EvaluateException('Sub-elements in value',
                        ps.Backtrace(c_elt))

        # It *seems* to be consensus that ignoring comments and
        # concatenating the text nodes is the right thing to do.
        return ''.join([E.nodeValue for E in c
                if E.nodeType 
                in [ _Node.TEXT_NODE, _Node.CDATA_SECTION_NODE ]])

    def parse_attributes(self, elt, ps):
        '''find all attributes specified in the attribute_typecode_dict in
        current element tag, if an attribute is found set it in the 
        self.attributes dictionary.  Default to putting in String.
        Parameters:
            elt -- the DOM element being parsed
            ps -- the ParsedSoap object.
        '''
        if self.attribute_typecode_dict is None: 
            return
        
        attributes = {}
        for attr,what in self.attribute_typecode_dict.items():
            namespaceURI,localName = None,attr
            if type(attr) in _seqtypes: 
                namespaceURI,localName = attr
            value = _find_attrNodeNS(elt, namespaceURI, localName)
            self.logger.debug("Parsed Attribute (%s,%s) -- %s", 
                               namespaceURI, localName, value)

            # For Now just set it w/o any type interpretation.
            if value is None: continue
            attributes[attr] = what.text_to_data(value, elt, ps)

        return attributes

    def set_attributes(self, el, pyobj):
        '''Instance data attributes contains a dictionary 
        of keys (namespaceURI,localName) and attribute values.
        These values can be self-describing (typecode), or use
        attribute_typecode_dict to determine serialization.
        Paramters:
            el -- MessageInterface representing the element
            pyobj -- 
        '''
        if not hasattr(pyobj, self.attrs_aname):
            return

        if not isinstance(getattr(pyobj, self.attrs_aname), dict):
            raise TypeError,\
                'pyobj.%s must be a dictionary of names and values'\
                % self.attrs_aname

        for attr, value in getattr(pyobj, self.attrs_aname).items():
            namespaceURI,localName = None, attr
            if type(attr) in _seqtypes:
                namespaceURI, localName = attr

            what = None
            if getattr(self, 'attribute_typecode_dict', None) is not None:
                what = self.attribute_typecode_dict.get(attr)
                if what is None and namespaceURI is None:
                    what = self.attribute_typecode_dict.get(localName)

            # allow derived type
            if hasattr(value, 'typecode') and not isinstance(what, AnyType):
                if what is not None and not isinstance(value.typecode, what):
                    raise EvaluateException, \
                        'self-describing attribute must subclass %s'\
                        %what.__class__

                what = value.typecode
                
            self.logger.debug("attribute create -- %s", value)
            if isinstance(what, QName):
                what.set_prefix(el, value)
            
            #format the data
            if what is None:
                value = str(value)
            else:
                value = what.get_formatted_content(value)

            el.setAttributeNS(namespaceURI, localName, value)

    def set_attribute_xsi_type(self, el, **kw):
        '''if typed, set the xsi:type attribute 
        Paramters:
            el -- MessageInterface representing the element
        '''
        if kw.get('typed', self.typed):
            namespaceURI,typeName = kw.get('type', _get_xsitype(self))
            if namespaceURI and typeName:
                self.logger.debug("attribute: (%s, %s)", namespaceURI, typeName)
                el.setAttributeType(namespaceURI, typeName)

    def set_attribute_href(self, el, objid):
        '''set href attribute
        Paramters:
            el -- MessageInterface representing the element
            objid -- ID type, unique id
        '''
        el.setAttributeNS(None, 'href', "#%s" %objid)

    def set_attribute_id(self, el, objid):
        '''set id attribute
        Paramters:
            el -- MessageInterface representing the element
            objid -- ID type, unique id
        '''
        if self.unique is False:
            el.setAttributeNS(None, 'id', "%s" %objid)

    def get_name(self, name, objid):
        '''
        Paramters:
            name -- element tag
            objid -- ID type, unique id
        '''
        if type(name) is tuple:
            return name

        ns = self.nspname
        n = name or self.pname or ('E' + objid)
        return ns,n

    def has_attributes(self):
        '''Return True if Attributes are declared outside
        the scope of SOAP ('root', 'id', 'href'), and some
        attributes automatically handled (xmlns, xsi:type).
        '''
        if self.attribute_typecode_dict is None: return False
        return len(self.attribute_typecode_dict) > 0


class SimpleType(TypeCode):
    '''SimpleType -- consist exclusively of a tag, attributes, and a value
    class attributes:
        empty_content -- value representing an empty element.
    '''
    empty_content = None
    logger = _GetLogger('ZSI.TC.SimpleType')
    
    def parse(self, elt, ps):
        self.checkname(elt, ps)
        if len(_children(elt)) == 0:
            href = _find_href(elt)
            if not href:
                if self.nilled(elt, ps) is False:
                    # No content, no HREF, not NIL:  empty string
                    return self.text_to_data(self.empty_content, elt, ps)
                    
                # No content, no HREF, and is NIL...
                if self.nillable is True: 
                    return Nilled
                raise EvaluateException('Requiredstring missing',
                        ps.Backtrace(elt))
                        
            if href[0] != '#':
                return ps.ResolveHREF(href, self)
            
            elt = ps.FindLocalHREF(href, elt)
            self.checktype(elt, ps)
            if self.nilled(elt, ps): return Nilled
            if len(_children(elt)) == 0: 
                v = self.empty_content
            else:
                v = self.simple_value(elt, ps)
        else:
            v = self.simple_value(elt, ps)
            
        pyobj = self.text_to_data(v, elt, ps)
        
        # parse all attributes contained in attribute_typecode_dict 
        # (user-defined attributes), the values (if not None) will 
        # be keyed in self.attributes dictionary.
        if self.attribute_typecode_dict is not None:
            attributes = self.parse_attributes(elt, ps)
            if attributes:
                setattr(pyobj, self.attrs_aname, attributes)
        
        return pyobj

    def get_formatted_content(self, pyobj):
        raise NotImplementedError, 'method get_formatted_content is not implemented'

    def serialize_text_node(self, elt, sw, pyobj):
        '''Serialize without an element node.
        '''
        textNode = None
        if pyobj is not None:
            text = self.get_formatted_content(pyobj)
            if type(text) not in _stringtypes:
                raise TypeError, 'pyobj must be a formatted string'

            textNode = elt.createAppendTextNode(text)

        return textNode

    def serialize(self, elt, sw, pyobj, name=None, orig=None, **kw):
        '''Handles the start and end tags, and attributes.  callout
        to get_formatted_content to get the textNode value.
        Parameters:
            elt -- ElementProxy/DOM element 
            sw -- SoapWriter instance
            pyobj -- processed content
            
        KeyWord Parameters:
            name -- substitute name, (nspname,name) or name
            orig --
            
        '''
        objid = _get_idstr(pyobj)
        ns,n = self.get_name(name, objid)

        # nillable
        el = elt.createAppendElement(ns, n)
        if self.nillable is True and pyobj is Nilled:
            self.serialize_as_nil(el)
            return None

        # other attributes
        self.set_attributes(el, pyobj)

        # soap href attribute
        unique = self.unique or kw.get('unique', False)
        if unique is False and sw.Known(orig or pyobj):
            self.set_attribute_href(el, objid)
            return None

        # xsi:type attribute 
        if kw.get('typed', self.typed) is True:
            self.set_attribute_xsi_type(el, **kw)

        # soap id attribute
        if self.unique is False:
            self.set_attribute_id(el, objid)

        #Content, <empty tag/>c
        self.serialize_text_node(el, sw, pyobj)

        return el


class Any(TypeCode):
    '''When the type isn't defined in the schema, but must be specified
    in the incoming operation.
        parsemap -- a type to class mapping (updated by descendants), for
                parsing
        serialmap -- same, for (outgoing) serialization
    '''
    logger = _GetLogger('ZSI.TC.Any')
    parsemap, serialmap = {}, {}

    def __init__(self, pname=None, aslist=False, minOccurs=0, unique=False, **kw):
        TypeCode.__init__(self, pname, minOccurs=minOccurs, unique=unique, **kw)
        self.aslist = aslist
        self.kwargs = dict(aslist=aslist, unique=unique)
        self.kwargs.update(kw)

    # input arg v should be a list of tuples (name, value).
    def listify(self, v):
        if self.aslist: return [ k for j,k in v ]
        return dict(v)

    def parse_into_dict_or_list(self, elt, ps):
        c = _child_elements(elt)
        count = len(c)
        v = []
        if count == 0:
            href = _find_href(elt)
            if not href: return v
            elt = ps.FindLocalHREF(href, elt)
            self.checktype(elt, ps)
            c = _child_elements(elt)
            count = len(c)
            if count == 0: return self.listify(v)
        if self.nilled(elt, ps): return Nilled

        for c_elt in c:
            v.append((str(c_elt.localName), self.__class__(**self.kwargs).parse(c_elt, ps)))

        return self.listify(v)

    def parse(self, elt, ps):
        (ns,type) = self.checkname(elt, ps)
        if not type and self.nilled(elt, ps): return Nilled
        if len(_children(elt)) == 0:
            href = _find_href(elt)
            if not href:
                if self.minOccurs < 1:
                    if _is_xsd_or_soap_ns(ns):
                        parser = Any.parsemap.get((None,type))
                        if parser: return parser.parse(elt, ps)
                    if ((ns,type) == (SOAP.ENC,'Array') or 
                        (_find_arraytype(elt) or '').endswith('[0]')):
                        return []
                    return None
                raise EvaluateException('Required Any missing',
                        ps.Backtrace(elt))
            elt = ps.FindLocalHREF(href, elt)
            (ns,type) = self.checktype(elt, ps)
        if not type and elt.namespaceURI == SOAP.ENC:
            ns,type = SOAP.ENC, elt.localName
        if not type or (ns,type) == (SOAP.ENC,'Array'):
            if self.aslist or _find_arraytype(elt):
                return [ self.__class__(**self.kwargs).parse(e, ps)
                            for e in _child_elements(elt) ]
            if len(_child_elements(elt)) == 0:
                #raise EvaluateException("Any cannot parse untyped element",
                #        ps.Backtrace(elt))
                return self.simple_value(elt, ps)
            return self.parse_into_dict_or_list(elt, ps)
        parser = Any.parsemap.get((ns,type))
        if not parser and _is_xsd_or_soap_ns(ns):
            parser = Any.parsemap.get((None,type))
        if not parser:
            raise EvaluateException('''Any can't parse element''',
                    ps.Backtrace(elt))
        return parser.parse(elt, ps)

    def get_formatted_content(self, pyobj):
        tc = type(pyobj)
        if tc == types.InstanceType:
            tc = pyobj.__class__
            if hasattr(pyobj, 'typecode'):
                #serializer = pyobj.typecode.serialmap.get(tc)
                serializer = pyobj.typecode
            else:
                serializer = Any.serialmap.get(tc)
            if not serializer:
                tc = (types.ClassType, pyobj.__class__.__name__)
                serializer = Any.serialmap.get(tc)
        else:
            serializer = Any.serialmap.get(tc)
            if not serializer and isinstance(pyobj, time.struct_time):
                from ZSI.TCtimes import gDateTime
                serializer = gDateTime()
        if serializer:
            return serializer.get_formatted_content(pyobj)
        raise EvaluateException, 'Failed to find serializer for pyobj %s' %pyobj

    def serialize(self, elt, sw, pyobj, name=None, **kw):
        if hasattr(pyobj, 'typecode') and pyobj.typecode is not self:
            pyobj.typecode.serialize(elt, sw, pyobj, **kw)
            return

        objid = _get_idstr(pyobj)
        ns,n = self.get_name(name, objid)
        kw.setdefault('typed', self.typed)
        tc = type(pyobj)
        self.logger.debug('Any serialize -- %s', tc)
        if tc in _seqtypes:
            if self.aslist:
                array = elt.createAppendElement(ns, n)
                array.setAttributeType(SOAP.ENC, "Array")
                array.setAttributeNS(self.nspname, 'SOAP-ENC:arrayType', 
                    "xsd:anyType[" + str(len(pyobj)) + "]" )
                for o in pyobj:
                    #TODO maybe this should take **self.kwargs...
                    serializer = getattr(o, 'typecode', Any(**self.kwargs))
                    serializer.serialize(array, sw, o, name='element', **kw)
            else:
                struct = elt.createAppendElement(ns, n)
                for o in pyobj:
                    #TODO maybe this should take **self.kwargs...
                    serializer = getattr(o, 'typecode', Any(**self.kwargs))
                    serializer.serialize(struct, sw, o, **kw)
            return

        kw['name'] = (ns,n)
        if tc == types.DictType:
            el = elt.createAppendElement(ns, n)
            parentNspname = self.nspname # temporarily clear nspname for dict elements
            self.nspname = None
            for o,m in pyobj.items():
                if type(o) != types.StringType and type(o) != types.UnicodeType:
                    raise Exception, 'Dictionary implementation requires keys to be of type string (or unicode).' %pyobj
                kw['name'] = o
                kw.setdefault('typed', True)
                self.serialize(el, sw, m, **kw)
            # restore nspname
            self.nspname = parentNspname
            return
                
        if tc == types.InstanceType:
            tc = pyobj.__class__
            if hasattr(pyobj, 'typecode'):
                #serializer = pyobj.typecode.serialmap.get(tc)
                serializer = pyobj.typecode
            else:
                serializer = Any.serialmap.get(tc)
            if not serializer:
                tc = (types.ClassType, pyobj.__class__.__name__)
                serializer = Any.serialmap.get(tc)
        else:
            serializer = Any.serialmap.get(tc)
            if not serializer and isinstance(pyobj, time.struct_time):
                from ZSI.TCtimes import gDateTime
                serializer = gDateTime()

        if not serializer:
            # Last-chance; serialize instances as dictionary
            if pyobj is None:
                self.serialize_as_nil(elt.createAppendElement(ns, n))
            elif type(pyobj) != types.InstanceType:
                raise EvaluateException('''Any can't serialize ''' + \
                        repr(pyobj))
            else:
                self.serialize(elt, sw, pyobj.__dict__, **kw)
        else:
            # Try to make the element name self-describing
            tag = getattr(serializer, 'tag', None)
            if self.pname is not None:
                #serializer.nspname = self.nspname
                #serializer.pname = self.pname
                if "typed" not in kw:
                    kw['typed'] = False
            elif tag:
                if tag.find(':') == -1: tag = 'SOAP-ENC:' + tag
                kw['name'] = tag
                kw['typed'] = False

            serializer.unique = self.unique
            serializer.serialize(elt, sw, pyobj, **kw)
            # Reset TypeCode
            #serializer.nspname = None
            #serializer.pname = None


class String(SimpleType):
    '''A string type.
    '''
    empty_content = ''
    parselist = [ (None,'string') ]
    seriallist = [ types.StringType, types.UnicodeType ]
    type = (SCHEMA.XSD3, 'string')
    logger = _GetLogger('ZSI.TC.String')

    def __init__(self, pname=None, strip=True, **kw):
        TypeCode.__init__(self, pname, **kw)
        if kw.has_key('resolver'): self.resolver = kw['resolver']
        self.strip = strip

    def text_to_data(self, text, elt, ps):
        '''convert text into typecode specific data.
        Encode all strings as UTF-8, which will be type 'str'
        not 'unicode'
        '''
        if self.strip: text = text.strip()
        if self.pyclass is not None:
            return self.pyclass(text.encode(UNICODE_ENCODING))
        return text.encode(UNICODE_ENCODING)

    def get_formatted_content(self, pyobj):
        if type(pyobj) not in _stringtypes:
            pyobj = str(pyobj)
        if type(pyobj) == unicode: 
            return pyobj.encode(UNICODE_ENCODING)
        return pyobj


class URI(String):
    '''A URI.
    Class data:
        reserved -- urllib.quote will escape all reserved characters
             regardless of whether they are used for the reserved purpose.

    '''
    parselist = [ (None,'anyURI'),(SCHEMA.XSD3, 'anyURI')]
    type = (SCHEMA.XSD3, 'anyURI')
    logger = _GetLogger('ZSI.TC.URI')
    reserved = ";/?:@&=+$,"

    def text_to_data(self, text, elt, ps):
        '''text --> typecode specific data.
        '''
        return String.text_to_data(self, urldecode(text), elt, ps)   

    def get_formatted_content(self, pyobj):
        '''typecode data --> text
        '''
        u = urlencode(pyobj, self.reserved)
        return String.get_formatted_content(self, 
            u)


class QName(String):
    '''A QName type
    '''
    parselist = [ (None,'QName') ]
    type = (SCHEMA.XSD3, 'QName')
    logger = _GetLogger('ZSI.TC.QName')

    def __init__(self, pname=None, strip=1, **kw):
        String.__init__(self, pname, strip, **kw)
        self.prefix = None

    def get_formatted_content(self, pyobj):
        value = pyobj
        if isinstance(pyobj, tuple):
            namespaceURI,localName = pyobj
            if self.prefix is not None:
                value = "%s:%s" %(self.prefix,localName)
        return String.get_formatted_content(self, value)

    def set_prefix(self, elt, pyobj):
        '''use this method to set the prefix of the QName,
        method looks in DOM to find prefix or set new prefix.
        This method must be called before get_formatted_content.
        '''
        if isinstance(pyobj, tuple):
            namespaceURI,localName = pyobj
            self.prefix = elt.getPrefix(namespaceURI)

    def text_to_data(self, text, elt, ps):
        '''convert text into typecode specific data.
        '''
        prefix,localName = SplitQName(text)
        nsdict = ps.GetElementNSdict(elt)
        prefix = prefix or ''
        try:
            namespaceURI = nsdict[prefix]
        except KeyError, ex:
            raise EvaluateException('cannot resolve prefix(%s)'%prefix,
                ps.Backtrace(elt))
                
        v = (namespaceURI,localName)
        if self.pyclass is not None:
            return self.pyclass(v)    
        return v

    def serialize_text_node(self, elt, sw, pyobj):
        '''Serialize without an element node.
        '''
        self.set_prefix(elt, pyobj)
        return String.serialize_text_node(self, elt, sw, pyobj)


class Token(String):
    '''an xsd:token type
    '''
    parselist = [ (None, 'token') ]
    type = (SCHEMA.XSD3, 'token')
    logger = _GetLogger('ZSI.TC.Token')


class Base64String(String):
    '''A Base64 encoded string.
    '''
    parselist = [ (None,'base64Binary'), (SOAP.ENC, 'base64') ]
    type = (SOAP.ENC, 'base64')
    logger = _GetLogger('ZSI.TC.Base64String')

    def text_to_data(self, text, elt, ps):
        '''convert text into typecode specific data.
        '''
        val = b64decode(text.replace(' ', '').replace('\n','').replace('\r',''))
        if self.pyclass is not None:
            return self.pyclass(val)
        return val

    def get_formatted_content(self, pyobj):
        pyobj = '\n' + b64encode(pyobj)
        return String.get_formatted_content(self, pyobj)


class Base64Binary(String):
    parselist = [ (None,'base64Binary'), ]
    type = (SCHEMA.XSD3, 'base64Binary')
    logger = _GetLogger('ZSI.TC.Base64Binary')

    def text_to_data(self, text, elt, ps):
        '''convert text into typecode specific data.
        '''
        val = b64decode(text)
        if self.pyclass is not None:
            return self.pyclass(val) 
        return val

    def get_formatted_content(self, pyobj):
        pyobj = b64encode(pyobj).strip()
        return pyobj


class HexBinaryString(String):
    '''Hex-encoded binary (yuk).
    '''
    parselist = [ (None,'hexBinary') ]
    type = (SCHEMA.XSD3, 'hexBinary')
    logger = _GetLogger('ZSI.TC.HexBinaryString')

    def text_to_data(self, text, elt, ps):
        '''convert text into typecode specific data.
        '''
        val = hexdecode(text)
        if self.pyclass is not None:
            return self.pyclass(val) 
        return val

    def get_formatted_content(self, pyobj):
        pyobj = hexencode(pyobj).upper()
        return String.get_formatted_content(self, pyobj)


class XMLString(String):
    '''A string that represents an XML document
    '''
    logger = _GetLogger('ZSI.TC.XMLString')
    
    def __init__(self, pname=None, readerclass=None, **kw):
        String.__init__(self, pname, **kw)
        self.readerclass = readerclass

    def parse(self, elt, ps):
        if not self.readerclass:
            from xml.dom.ext.reader import PyExpat
            self.readerclass = PyExpat.Reader
        v = String.parse(self, elt, ps)
        return self.readerclass().fromString(v)

    def get_formatted_content(self, pyobj):
        #pyobj = Canonicalize(pyobj)
        return String.get_formatted_content(self, pyobj)


class Enumeration(String):
    '''A string type, limited to a set of choices.
    '''
    logger = _GetLogger('ZSI.TC.Enumeration')
    
    def __init__(self, choices, pname=None, **kw):
        String.__init__(self, pname, **kw)
        t = type(choices)
        if t in _seqtypes:
            self.choices = tuple(choices)
        elif TypeCode.typechecks:
            raise TypeError(
                'Enumeration choices must be list or sequence, not ' + str(t))
        if TypeCode.typechecks:
            for c in self.choices:
                if type(c) not in _stringtypes:
                    raise TypeError(
                        'Enumeration choice ' + str(c) + ' is not a string')

    def parse(self, elt, ps):
        val = String.parse(self, elt, ps)
        if val not in self.choices:
            raise EvaluateException('Value not in enumeration list',
                    ps.Backtrace(elt))
        return val

    def serialize(self, elt, sw, pyobj, name=None, orig=None, **kw):
        if pyobj not in self.choices:
            raise EvaluateException('Value not in enumeration list',
                    sw.Backtrace(elt))
        String.serialize(self, elt, sw, pyobj, name=name, orig=orig, **kw)



# This is outside the Integer class purely for code esthetics.
_ignored = []

class Integer(SimpleType):
    '''Common handling for all integers.
    '''

    ranges = {
        'unsignedByte':         (0, 255),
        'unsignedShort':        (0, 65535),
        'unsignedInt':          (0, 4294967295L),
        'unsignedLong':         (0, 18446744073709551615L),

        'byte':                 (-128, 127),
        'short':                (-32768, 32767),
        'int':                  (-2147483648L, 2147483647),
        'long':                 (-9223372036854775808L, 9223372036854775807L),

        'negativeInteger':      (_ignored, -1),
        'nonPositiveInteger':   (_ignored, 0),
        'nonNegativeInteger':   (0, _ignored),
        'positiveInteger':      (1, _ignored),

        'integer':              (_ignored, _ignored)
    }
    parselist = [ (None,k) for k in ranges.keys() ]
    seriallist = [ types.IntType, types.LongType ]
    logger = _GetLogger('ZSI.TC.Integer')

    def __init__(self, pname=None, format='%d', **kw):
        TypeCode.__init__(self, pname, **kw)
        self.format = format

    def text_to_data(self, text, elt, ps):
        '''convert text into typecode specific data.
        '''
        if self.pyclass is not None:
            v = self.pyclass(text) 
        else:
            try:
                v = int(text)
            except:
                try:
                    v = long(text)
                except:
                    raise EvaluateException('Unparseable integer', 
                        ps.Backtrace(elt))
        return v

    def parse(self, elt, ps):
        (ns,type) = self.checkname(elt, ps)
        if self.nilled(elt, ps): return Nilled
        elt = self.SimpleHREF(elt, ps, 'integer')
        if not elt: return None

        if type is None:
           type = self.type[1] 
        elif self.type[1] is not None and type != self.type[1]:
            raise EvaluateException('Integer type mismatch; ' \
                'got %s wanted %s' % (type,self.type[1]), ps.Backtrace(elt))
        
        v = self.simple_value(elt, ps)
        v = self.text_to_data(v, elt, ps)

        (rmin, rmax) = Integer.ranges.get(type, (_ignored, _ignored))
        if rmin != _ignored and v < rmin:
            raise EvaluateException('Underflow, less than ' + repr(rmin),
                    ps.Backtrace(elt))
        if rmax != _ignored and v > rmax:
            raise EvaluateException('Overflow, greater than ' + repr(rmax),
                    ps.Backtrace(elt))
        return v

    def get_formatted_content(self, pyobj):
        return self.format %pyobj



# See credits, below.
def _make_inf():
    x = 2.0
    x2 = x * x
    i = 0
    while i < 100 and x != x2:
        x = x2
        x2 = x * x
        i = i + 1
    if x != x2:
        raise ValueError("This machine's floats go on forever!")
    return x

# This is outside the Decimal class purely for code esthetics.
_magicnums = { }
try:
    _magicnums['INF'] = float('INF')
    _magicnums['-INF'] = float('-INF')
except:
    _magicnums['INF'] = _make_inf()
    _magicnums['-INF'] = -_magicnums['INF']

# The following comment and code was written by Tim Peters in
# article <001401be92d2$09dcb800$5fa02299@tim> in comp.lang.python,
# also available at the following URL:
# http://groups.google.com/groups?selm=001401be92d2%2409dcb800%245fa02299%40tim
# Thanks, Tim!

# NaN-testing.
#
# The usual method (x != x) doesn't work.
# Python forces all comparisons thru a 3-outcome cmp protocol; unordered
# isn't a possible outcome.  The float cmp outcome is essentially defined
# by this C expression (combining some cross-module implementation
# details, and where px and py are pointers to C double):
#   px == py ? 0 : *px < *py ? -1 : *px > *py ? 1 : 0
# Comparing x to itself thus always yields 0 by the first clause, and so
# x != x is never true.
# If px and py point to distinct NaN objects, a strange thing happens:
# 1. On scrupulous 754 implementations, *px < *py returns false, and so
#    does *px > *py.  Python therefore returns 0, i.e. "equal"!
# 2. On Pentium HW, an unordered outcome sets an otherwise-impossible
#    combination of condition codes, including both the "less than" and
#    "equal to" flags.  Microsoft C generates naive code that accepts
#    the "less than" flag at face value, and so the *px < *py clause
#    returns true, and Python returns -1, i.e. "not equal".
# So with a proper C 754 implementation Python returns the wrong result,
# and under MS's improper 754 implementation Python yields the right
# result -- both by accident.  It's unclear who should be shot <wink>.
#
# Anyway, the point of all that was to convince you it's tricky getting
# the right answer in a portable way!
def isnan(x):
    """x -> true iff x is a NaN."""
    # multiply by 1.0 to create a distinct object (x < x *always*
    # false in Python, due to object identity forcing equality)
    if x * 1.0 < x:
        # it's a NaN and this is MS C on a Pentium
        return 1
    # Else it's non-NaN, or NaN on a non-MS+Pentium combo.
    # If it's non-NaN, then x == 1.0 and x == 2.0 can't both be true,
    # so we return false.  If it is NaN, then assuming a good 754 C
    # implementation Python maps both unordered outcomes to true.
    return 1.0 == x and x == 2.0


class Decimal(SimpleType):
    '''Parent class for floating-point numbers.
    '''

    parselist = [ (None,'decimal'), (None,'float'), (None,'double') ]
    seriallist = _floattypes
    type = None
    ranges =  {
        'float': ( 7.0064923216240861E-46,
                        -3.4028234663852886E+38, 3.4028234663852886E+38 ),
        'double': ( 2.4703282292062327E-324,
                        -1.7976931348623158E+308, 1.7976931348623157E+308),
    }
    zeropat = re.compile('[1-9]')
    logger = _GetLogger('ZSI.TC.Decimal')

    def __init__(self, pname=None, format='%f', **kw):
        TypeCode.__init__(self, pname, **kw)
        self.format = format


    def text_to_data(self, text, elt, ps):
        '''convert text into typecode specific data.
        '''
        v = text
        if self.pyclass is not None:
            return self.pyclass(v)

        m = _magicnums.get(v)
        if m: return m

        try:
            return float(v)
        except:
            raise EvaluateException('Unparseable floating point number',
                    ps.Backtrace(elt))

    def parse(self, elt, ps):
        (ns,type) = self.checkname(elt, ps)
        elt = self.SimpleHREF(elt, ps, 'floating-point')
        if not elt: return None
        tag = getattr(self.__class__, 'type')
        if tag:
            if type is None:
                type = tag
            elif tag != (ns,type):
                raise EvaluateException('Floating point type mismatch; ' \
                        'got (%s,%s) wanted %s' % (ns,type,tag), ps.Backtrace(elt))
        # Special value?
        if self.nilled(elt, ps): return Nilled
        v = self.simple_value(elt, ps)
        try:
            fp = self.text_to_data(v, elt, ps)
        except EvaluateException, ex:
            ex.args.append(ps.Backtrace(elt))
            raise ex
   
        m = _magicnums.get(v)
        if m: 
            return m

        if str(fp).lower() in [ 'inf', '-inf', 'nan', '-nan' ]:
            raise EvaluateException('Floating point number parsed as "' + \
                    str(fp) + '"', ps.Backtrace(elt))
        if fp == 0 and Decimal.zeropat.search(v):
            raise EvaluateException('Floating point number parsed as zero',
                    ps.Backtrace(elt))
        (rtiny, rneg, rpos) = Decimal.ranges.get(type, (None, None, None))
        if rneg and fp < 0 and fp < rneg:
            raise EvaluateException('Negative underflow', ps.Backtrace(elt))
        if rtiny and fp > 0 and fp < rtiny:
            raise EvaluateException('Positive underflow', ps.Backtrace(elt))
        if rpos and fp > 0 and fp > rpos:
            raise EvaluateException('Overflow', ps.Backtrace(elt))
        return fp

    def get_formatted_content(self, pyobj):
        if pyobj == _magicnums['INF']:
            return 'INF'
        elif pyobj == _magicnums['-INF']:
            return '-INF'
        elif isnan(pyobj):
            return 'NaN'
        else:
            return self.format %pyobj


class Boolean(SimpleType):
    '''A boolean.
    '''

    parselist = [ (None,'boolean') ]
    seriallist = [ bool ]
    type = (SCHEMA.XSD3, 'boolean')
    logger = _GetLogger('ZSI.TC.Boolean')
    
    def text_to_data(self, text, elt, ps):
        '''convert text into typecode specific data.
        '''
        v = text
        if v == 'false': 
            if self.pyclass is None:
                return False
            return self.pyclass(False)

        if v == 'true': 
            if self.pyclass is None:
                return True
            return self.pyclass(True)

        try:
            v = int(v)
        except:
            try:
                v = long(v)
            except:
                raise EvaluateException('Unparseable boolean', 
                        ps.Backtrace(elt))

        if v:
            if self.pyclass is None:
                return True
            return self.pyclass(True)

        if self.pyclass is None:
             return False
        return self.pyclass(False)

    def parse(self, elt, ps):
        self.checkname(elt, ps)
        elt = self.SimpleHREF(elt, ps, 'boolean')
        if not elt: return None
        if self.nilled(elt, ps): return Nilled

        v = self.simple_value(elt, ps).lower()
        return self.text_to_data(v, elt, ps)

    def get_formatted_content(self, pyobj):
        if pyobj: return 'true'
        return 'false'


#XXX NOT FIXED YET
class XML(TypeCode):
    '''Opaque XML which shouldn't be parsed.
        comments -- preserve comments
        inline -- don't href/id when serializing
        resolver -- object to resolve href's
        wrapped -- put a wrapper element around it
    '''

    # Clone returned data?
    copyit = 0
    logger = _GetLogger('ZSI.TC.XML')
    
    def __init__(self, pname=None, comments=0, inline=0, wrapped=True, **kw):
        TypeCode.__init__(self, pname, **kw)
        self.comments = comments
        self.inline = inline
        if kw.has_key('resolver'): self.resolver = kw['resolver']
        self.wrapped = wrapped
        self.copyit = kw.get('copyit', XML.copyit)

    def parse(self, elt, ps):
        if self.wrapped is False:
            return elt
        c = _child_elements(elt)
        if not c:
            href = _find_href(elt)
            if not href:
                if self.minOccurs == 0: return None
                raise EvaluateException('Embedded XML document missing',
                        ps.Backtrace(elt))
            if href[0] != '#':
                return ps.ResolveHREF(href, self)
            elt = ps.FindLocalHREF(href, elt)
            c = _child_elements(elt)
        if _find_encstyle(elt) != "":
            #raise EvaluateException('Embedded XML has unknown encodingStyle',
            #       ps.Backtrace(elt)
            pass
        if len(c) != 1:
            raise EvaluateException('Embedded XML has more than one child',
                    ps.Backtrace(elt))
        if self.copyit: return c[0].cloneNode(1)
        return c[0]

    def serialize(self, elt, sw, pyobj, name=None, unsuppressedPrefixes=[], **kw):
        objid = _get_idstr(pyobj)
        ns,n = self.get_name(name, objid)

        xmlelt = elt
        if self.wrapped:
            xmlelt = elt.createAppendElement(ns, n)

        #if type(pyobj) in _stringtypes:
        #    self.set_attributes(xmlelt, pyobj)
        #    self.set_attribute_href(xmlelt, objid)
        #elif kw.get('inline', self.inline):
        #    self.cb(xmlelt, sw, pyobj, unsuppressedPrefixes)
        #else:
        #    self.set_attributes(xmlelt, pyobj)
        #    self.set_attribute_href(xmlelt, objid)
        #    sw.AddCallback(self.cb, elt, sw, pyobj, unsuppressedPrefixes)

        self.cb(xmlelt, sw, pyobj, unsuppressedPrefixes)

    def cb(self, elt, sw, pyobj, unsuppressedPrefixes=[]):
        """pyobj -- xml.dom.Node.ELEMENT_NODE
        """
        #if sw.Known(pyobj): 
        #    return

        if type(pyobj) in _stringtypes:
            elt.createAppendTextNode(pyobj)
            return

        ## grab document and import node, and append it
        doc = elt.getDocument()
        node = doc.importNode(pyobj, deep=1)
        child = elt.node.appendChild(node)

        ## copy xmlns: attributes into appended node
        parent = pyobj.parentNode
        while parent.nodeType == _Node.ELEMENT_NODE:
            for attr in filter(lambda a: a.name.startswith('xmlns:') and a.name not in child.attributes.keys(), parent.attributes): 
                child.setAttributeNode(attr.cloneNode(1))

            parent = parent.parentNode


class AnyType(TypeCode):
    """XML Schema xsi:anyType type definition wildCard.
       class variables: 
          all -- specifies use of all namespaces.
          other -- specifies use of other namespaces
          type --
    """
    all = '#all'
    other = '#other'
    type = (SCHEMA.XSD3, 'anyType')
    logger = _GetLogger('ZSI.TC.AnyType')
    
    def __init__(self, pname=None, namespaces=['#all'],
    minOccurs=1, maxOccurs=1, strip=1, **kw):
        TypeCode.__init__(self, pname=pname, minOccurs=minOccurs, 
              maxOccurs=maxOccurs, **kw)
        self.namespaces = namespaces

    def get_formatted_content(self, pyobj):
        # TODO: not sure this makes sense, 
        # parse side will be clueless, but oh well..
        what = getattr(pyobj, 'typecode', Any())
        return what.get_formatted_content(pyobj)

    def text_to_data(self, text, elt, ps):
        '''convert text into typecode specific data.  Used only with
        attributes so will not know anything about this content so
        why guess?
        Parameters:
            text -- text content
            elt -- the DOM element being parsed
            ps -- the ParsedSoap object.
        '''
        return text

    def serialize(self, elt, sw, pyobj, **kw):
        nsuri,typeName = _get_xsitype(pyobj)
        if self.all not in self.namespaces and nsuri not in self.namespaces:
            raise EvaluateException(
                '<anyType> unsupported use of namespaces "%s"' %self.namespaces)
        
        what = getattr(pyobj, 'typecode', None)
        if what is None:
            # TODO: resolve this, "strict" processing but no 
            # concrete schema makes little sense.
            #what = _AnyStrict(pname=(self.nspname,self.pname))
            what = Any(pname=(self.nspname,self.pname), unique=True, 
                       aslist=False)
            kw['typed'] = True
            what.serialize(elt, sw, pyobj, **kw)
            return

        # Namespace if element AnyType was namespaced.
        what.serialize(elt, sw, pyobj, 
           name=(self.nspname or what.nspname, self.pname or what.pname), **kw)

    def parse(self, elt, ps):
        #element name must be declared ..
        nspname,pname = _get_element_nsuri_name(elt)
        if nspname != self.nspname or pname != self.pname:
            raise EvaluateException('<anyType> instance is (%s,%s) found (%s,%s)' %(
                    self.nspname,self.pname,nspname,pname), ps.Backtrace(elt))

        #locate xsi:type
        prefix, typeName = SplitQName(_find_type(elt))
        namespaceURI = _resolve_prefix(elt, prefix)
        pyclass = GTD(namespaceURI, typeName)
        if not pyclass:
            if _is_xsd_or_soap_ns(namespaceURI):
                pyclass = Any
            elif (str(namespaceURI).lower()==str(Apache.Map.type[0]).lower())\
                and (str(typeName).lower() ==str(Apache.Map.type[1]).lower()):
                pyclass = Apache.Map
            else:
                # Unknown type, so parse into a dictionary
                pyobj = Any().parse_into_dict_or_list(elt, ps)
                return pyobj
                    
        what = pyclass(pname=(self.nspname,self.pname))
        pyobj = what.parse(elt, ps)
        return pyobj


class AnyElement(AnyType):
    """XML Schema xsi:any element declaration wildCard.
       class variables: 
            tag -- global element declaration
    """
    tag = (SCHEMA.XSD3, 'any')
    logger = _GetLogger('ZSI.TC.AnyElement')
    
    def __init__(self, namespaces=['#all'],pname=None, 
        minOccurs=1, maxOccurs=1, strip=1, processContents='strict',
        **kw):
        
        if processContents not in ('lax', 'skip', 'strict'):
            raise ValueError('processContents(%s) must be lax, skip, or strict')
            
        self.processContents = processContents
        AnyType.__init__(self, namespaces=namespaces,pname=pname,
            minOccurs=minOccurs, maxOccurs=maxOccurs, strip=strip, **kw)
       
    def serialize(self, elt, sw, pyobj, **kw):
        '''Must provice typecode to AnyElement for serialization, else
        try to use TC.Any to serialize instance which will serialize 
        based on the data type of pyobj w/o reference to XML schema 
        instance.
        '''
        if isinstance(pyobj, TypeCode):
            raise TypeError, 'pyobj is a typecode instance.'
        
        what = getattr(pyobj, 'typecode', None)
        if what is not None and type(pyobj) is types.InstanceType:
            tc = pyobj.__class__
            what = Any.serialmap.get(tc)
            if not what:
                tc = (types.ClassType, pyobj.__class__.__name__)
                what = Any.serialmap.get(tc)
        
        self.logger.debug('processContents: %s', self.processContents)

        # failed to find a registered type for class
        if what is None:
            #TODO: seems incomplete.  what about facets.
            #if self.processContents == 'strict':
            what = Any(pname=(self.nspname,self.pname))
                
        self.logger.debug('serialize with %s', what.__class__.__name__)
        what.serialize(elt, sw, pyobj, **kw)

    def parse(self, elt, ps):
        '''
        processContents -- 'lax' | 'skip' | 'strict', 'strict'
        1) if 'skip' check namespaces, and return the DOM node.
        2) if 'lax' look for declaration, or definition.  If
           not found return DOM node.
        3) if 'strict' get declaration, or raise.
        '''
        skip = self.processContents == 'skip'
        nspname,pname = _get_element_nsuri_name(elt)
        what = GED(nspname, pname)
        if not skip and what is not None:
            pyobj = what.parse(elt, ps)
            try:
                pyobj.typecode = what
            except AttributeError, ex:
                # Assume this means builtin type.
                pyobj = WrapImmutable(pyobj, what)
            return pyobj
        
        # Allow use of "<any>" element declarations w/ local
        # element declarations
        prefix, typeName = SplitQName(_find_type(elt))
        if not skip and typeName:
            namespaceURI = _resolve_prefix(elt, prefix or 'xmlns')
            # First look thru user defined namespaces, if don't find
            # look for 'primitives'.
            pyclass = GTD(namespaceURI, typeName) or Any
            what = pyclass(pname=(nspname,pname))
            pyobj = what.parse(elt, ps)
            try:
                pyobj.typecode = what
            except AttributeError, ex:
                # Assume this means builtin type.
                pyobj = WrapImmutable(pyobj, what)
                
            what.typed = True
            return pyobj

        if skip:
            what = XML(pname=(nspname,pname), wrapped=False)
        elif self.processContents == 'lax':
            what = Any(pname=(nspname,pname), unique=True)
        else:
            what = Any(pname=(nspname,pname), unique=True)

        try:
            pyobj = what.parse(elt, ps)
        except EvaluateException, ex:
            self.logger.debug("error parsing:  %s" %str(ex))

            if len(_children(elt)) != 0:
                self.logger.debug('parse <any>, return as dict')
                return Any(aslist=False).parse_into_dict_or_list(elt, ps)

            self.logger.debug("Give up, parse (%s,%s) as a String", 
                  what.nspname, what.pname)
            what = String(pname=(nspname,pname), typed=False)
            return WrapImmutable(what.parse(elt, ps), what)

        if pyobj is None: 
            return

        # dict is elementName:value pairs
        if type(pyobj) is dict:
            return pyobj

        try:
            pyobj.typecode = what
        except AttributeError:
            pyobj = WrapImmutable(pyobj, what)

        return pyobj  



class Union(SimpleType):
    '''simpleType Union

    class variables:
        memberTypes -- list [(namespace,name),] tuples, each representing a type defintion.
    '''
    memberTypes = None
    logger = _GetLogger('ZSI.TC.Union')
    
    def __init__(self, pname=None, minOccurs=1, maxOccurs=1, **kw):
        SimpleType.__init__(self, pname=pname, minOccurs=minOccurs, maxOccurs=maxOccurs, **kw)
        self.memberTypeCodes = []

    def setMemberTypeCodes(self):
        if len(self.memberTypeCodes) > 0: 
            return
        if self.__class__.memberTypes is None:
            raise EvaluateException, 'uninitialized class variable memberTypes [(namespace,name),]'
        for nsuri,name in self.__class__.memberTypes:
            tcclass = GTD(nsuri,name)
            if tcclass is None:
                tc = Any.parsemap.get((nsuri,name)) or Any.parsemap.get((None, name))
                typecode = tc.__class__(pname=(self.nspname,self.pname))
            else:
                typecode = tcclass(pname=(self.nspname,self.pname))

            if typecode is None:
                raise EvaluateException, \
                    'Typecode class for Union memberType (%s,%s) is missing' %(nsuri,name)
            if isinstance(typecode, Struct):
                raise EvaluateException, \
                    'Illegal: Union memberType (%s,%s) is complexType' %(nsuri,name)
            self.memberTypeCodes.append(typecode)

    def parse(self, elt, ps, **kw):
        '''attempt to parse sequentially.  No way to know ahead of time
        what this instance represents.  Must be simple type so it can
        not have attributes nor children, so this isn't too bad.
        '''
        self.setMemberTypeCodes()
        (nsuri,typeName) = self.checkname(elt, ps)

        #if (nsuri,typeName) not in self.memberTypes:
        #    raise EvaluateException(
        #            'Union Type mismatch got (%s,%s) not in %s' % \
        #            (nsuri, typeName, self.memberTypes), ps.Backtrace(elt))

        for indx in range(len(self.memberTypeCodes)):
            typecode = self.memberTypeCodes[indx]
            try:
                pyobj = typecode.parse(elt, ps)
            except ParseException, ex:
                continue
            except Exception, ex:
                continue

            if indx > 0:
                self.memberTypeCodes.remove(typecode)
                self.memberTypeCodes.insert(0, typecode)
            break

        else:
            raise

        return pyobj

    def get_formatted_content(self, pyobj, **kw): 
        self.setMemberTypeCodes()
        for indx in range(len(self.memberTypeCodes)):
            typecode = self.memberTypeCodes[indx]
            try:
                content = typecode.get_formatted_content(copy.copy(pyobj))
                break
            except (ParseException, TypeError):
                pass

            if indx > 0:
                self.memberTypeCodes.remove(typecode)
                self.memberTypeCodes.insert(0, typecode)

        else:
            raise

        return content


class List(SimpleType):
    '''simpleType List
    Class data:
        itemType -- sequence (namespaceURI,name) or a TypeCode instance
            representing the type definition
    '''
    itemType = None
    logger = _GetLogger('ZSI.TC.List')
    
    def __init__(self, pname=None, itemType=None, **kw):
        '''Currently need to require maxOccurs=1, so list
        is interpreted as a single unit of data.
        '''
        assert kw.get('maxOccurs',1) == 1, \
            'Currently only supporting SimpleType Lists with  maxOccurs=1'

        SimpleType.__init__(self, pname=pname, **kw)
        self.itemType = itemType or self.itemType
        self.itemTypeCode = self.itemType

        itemTypeCode = None
        if type(self.itemTypeCode) in _seqtypes:
            namespaceURI,name = self.itemTypeCode
            try:
                itemTypeCode = GTD(*self.itemType)(None)
            except:
                if _is_xsd_or_soap_ns(namespaceURI) is False:
                    raise
                for pyclass in TYPES:
                    if pyclass.type == self.itemTypeCode:
                        itemTypeCode =  pyclass(None)
                        break
                    elif pyclass.type[1] == name:
                        itemTypeCode =  pyclass(None)

                if itemTypeCode is None:
                    raise EvaluateException('Failed to locate %s' %str(self.itemTypeCode))

            if hasattr(itemTypeCode, 'text_to_data') is False:
                raise EvaluateException('TypeCode class %s missing text_to_data method' %itemTypeCode)

            self.itemTypeCode = itemTypeCode


    def text_to_data(self, text, elt, ps):
        '''convert text into typecode specific data.  items in
        list are space separated.
        '''
        v = []
        items = text.split()
        for item in items:
            v.append(self.itemTypeCode.text_to_data(item, elt, ps))

        if self.pyclass is not None:
            return self.pyclass(v)
        return v

    def parse(self, elt, ps):
        '''elt -- the DOM element being parsed
        ps -- the ParsedSoap object.
        '''
        self.checkname(elt, ps)
        if len(_children(elt)) == 0:
            href = _find_href(elt)
            if not href:
                if self.nilled(elt, ps) is False:
                    return []
                if self.nillable is True: 
                    return Nilled
                raise EvaluateException('Required string missing',
                        ps.Backtrace(elt))
            if href[0] != '#':
                return ps.ResolveHREF(href, self)
            elt = ps.FindLocalHREF(href, elt)
            self.checktype(elt, ps)

        if self.nilled(elt, ps): return Nilled
        if len(_children(elt)) == 0: return []

        v = self.simple_value(elt, ps)
        return self.text_to_data(v, elt, ps)


    def serialize(self, elt, sw, pyobj, name=None, orig=None, **kw):
        '''elt -- the current DOMWrapper element 
           sw -- soapWriter object
           pyobj -- python object to serialize
        '''
        if pyobj is not None and type(pyobj) not in _seqtypes:
            raise EvaluateException, 'expecting a list or None'

        objid = _get_idstr(pyobj)
        ns,n = self.get_name(name, objid)
        el = elt.createAppendElement(ns, n)
        if self.nillable is True and pyobj is None:
            self.serialize_as_nil(el)
            return None
        
        tc = self.itemTypeCode
        s = StringIO(); sep = ' '
        for item in pyobj:
            s.write(tc.get_formatted_content(item))
            s.write(sep)

        el.createAppendTextNode(s.getvalue())


def RegisterType(C, clobber=0, *args, **keywords):
    instance = apply(C, args, keywords)
    for t in C.__dict__.get('parselist', []):
        prev = Any.parsemap.get(t)
        if prev:
            if prev.__class__ == C: continue
            if not clobber:
                raise TypeError(
                    str(C) + ' duplicating parse registration for ' + str(t))
        Any.parsemap[t] = instance
    for t in C.__dict__.get('seriallist', []):
        ti = type(t)
        if ti in [ types.TypeType, types.ClassType]:
            key = t
        elif ti in _stringtypes:
            key = (types.ClassType, t)
        else:
            raise TypeError(str(t) + ' is not a class name')
        prev = Any.serialmap.get(key)
        if prev:
            if prev.__class__ == C: continue
            if not clobber:
                raise TypeError(
                    str(C) + ' duplicating serial registration for ' + str(t))
        Any.serialmap[key] = instance


#def _DynamicImport(moduleName, className):
#    '''
#    Utility function for RegisterTypeWithSchemaAndClass
#    '''
#    mod = __import__(moduleName)
#    components = moduleName.split('.')
#    for comp in components[1:]:
#        mod = getattr(mod, comp)
#    return getattr(mod, className)
#
#def _RegisterTypeWithSchemaAndClass(importedSchemaTypes, schemaTypeName, classModuleName, className, generatedClassSuffix="_"):
#    '''
#    Used by RegisterGeneratedTypesWithMapping.
#    Helps register classes so they can be serialized and parsed as "any".
#    Register a type by providing its schema and class.  This allows
#       Any and AnyType to reconstruct objects made up of your own classes.
#       Note: The class module should be able to be imported (by being in your
#       pythonpath).  Your classes __init__ functions shoud have default
#       arguments for all extra parameters.
#    Example of use:
#        import SchemaToPyTypeMap # Mapping written by you.  Also used with wsdl2py -m
#             # mapping = {"SomeDescription":("Descriptions", "SomeDescription"),
#             #             schemaTypeName  :  moduleName   ,  className 
#        # The module on the next line is generated by wsdl2py
#        from EchoServer_services_types import urn_ZSI_examples as ExampleTypes
#
#        for key,value in SchemaToPyTypeMap.mapping.items():
#        ZSI.TC.RegisterTypeWithSchemaAndClass(importedSchemaTypes = ExampleTypes, schemaTypeName=key, classModuleName=value[0], className=value[1])
#
#    '''
#    # Doing this: (schemaTypeName="ExampleTypes", classModuleName="Description",
#    #               className="SomeDescription")
#    # sd_instance = ExampleTypes.SomeDescription_(pname="SomeDescription")
#    # Any.serialmap[Descriptions.SomeDescription] = sd_instance
#    # Any.parsemap[(None,'SomeDescription')] = sd_instance
#    classDef = _DynamicImport(classModuleName, className)
#    interfaceDef = getattr(importedSchemaTypes, schemaTypeName + generatedClassSuffix)
#
#    instance = interfaceDef(pname=className)
#    Any.serialmap[classDef] = instance
#    Any.parsemap[(None,schemaTypeName)] = instance
#
#def RegisterGeneratedTypesWithMapping(generatedTypes, mapping, generatedClassSuffix="_"):
#    '''
#    Registers python classes so they can be serialized and parsed as "any".
#        generatedTypes is a class containing typecode classes generated by zsi.
#        mapping is a dictionary that maps
#        {schemaTypeName : moduleName, className}
#        and is also used with wsdl2py -m
#
#    Example of use:
#        import SchemaToPyTypeMap      # See RegisterTypeWithSchemaAndClass for description
#        # The module on the next line is generated by wsdl2py and
#        #    contains generated typecodes.
#        from EchoServer_services_types import urn_ZSI_examples as ExampleTypes
#        RegisterGeneratedTypesWithMapping(generatedTypes = ExampleTypes, mapping=SchemaToPyTypeMap.mapping)
#    '''
#    for key,value in mapping.items():
#        _RegisterTypeWithSchemaAndClass(importedSchemaTypes = generatedTypes, schemaTypeName=key, classModuleName=value[0], className=value[1], generatedClassSuffix=generatedClassSuffix)


from TCnumbers import *
from TCtimes import *
from schema import GTD, GED, WrapImmutable
from TCcompound import *
from TCapache import *

# aliases backwards compatiblity
_get_type_definition, _get_global_element_declaration, Wrap  = GTD, GED, WrapImmutable

f = lambda x: type(x) == types.ClassType and issubclass(x, TypeCode) and getattr(x, 'type', None) is not None
TYPES = filter(f, map(lambda y:eval(y),dir()))


if __name__ == '__main__': print _copyright

