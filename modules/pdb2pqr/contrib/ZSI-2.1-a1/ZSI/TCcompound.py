#! /usr/bin/env python
# $Header$
'''Compound typecodes.
'''

from ZSI import _copyright, _children, _child_elements, \
    _inttypes, _stringtypes, _seqtypes, _find_arraytype, _find_href, \
    _find_type, _find_xmlns_prefix, _get_idstr, EvaluateException, \
    ParseException
    
from TC import _get_element_nsuri_name, \
     _get_xsitype, TypeCode, Any, AnyElement, AnyType, \
     Nilled, UNBOUNDED
    
from schema import GED, ElementDeclaration, TypeDefinition, \
    _get_substitute_element, _get_type_definition, _is_substitute_element

from ZSI.wstools.Namespaces import SCHEMA, SOAP
from ZSI.wstools.Utility import SplitQName
from ZSI.wstools.logging import getLogger as _GetLogger
import re, types
from copy import copy as _copy

_find_arrayoffset = lambda E: E.getAttributeNS(SOAP.ENC, "offset")
_find_arrayposition = lambda E: E.getAttributeNS(SOAP.ENC, "position")

_offset_pat = re.compile(r'\[[0-9]+\]')
_position_pat = _offset_pat

def _check_typecode_list(ofwhat, tcname):
    '''Check a list of typecodes for compliance with Struct
    requirements.'''
    for o in ofwhat:
        if callable(o): #skip if _Mirage
            continue
        if not isinstance(o, TypeCode):
            raise TypeError(
                tcname + ' ofwhat outside the TypeCode hierarchy, ' +
                str(o.__class__))
        if o.pname is None and not isinstance(o, AnyElement):
            raise TypeError(tcname + ' element ' + str(o) + ' has no name')


def _get_type_or_substitute(typecode, pyobj, sw, elt):
    '''return typecode or substitute type for wildcard or
    derived type.  For serialization only.
    '''
    sub = getattr(pyobj, 'typecode', typecode)
    if sub is typecode or sub is None:
        return typecode

    # Element WildCard
    if isinstance(typecode, AnyElement):
        return sub
 
    # Global Element Declaration
    if isinstance(sub, ElementDeclaration):
        if (typecode.nspname,typecode.pname) == (sub.nspname,sub.pname):
            raise TypeError(\
                'bad usage, failed to serialize element reference (%s, %s), in: %s' %
                 (typecode.nspname, typecode.pname, sw.Backtrace(elt),))

        # check substitutionGroup 
        if _is_substitute_element(typecode, sub):
            return sub

        raise TypeError(\
            'failed to serialize (%s, %s) illegal sub GED (%s,%s): %s' %
             (typecode.nspname, typecode.pname, sub.nspname, sub.pname,
              sw.Backtrace(elt),))

    # Local Element
    if not isinstance(typecode, AnyType) and not isinstance(sub, typecode.__class__):
        raise TypeError(\
            'failed to serialize substitute %s for %s,  not derivation: %s' %
             (sub, typecode, sw.Backtrace(elt),))

    # Make our substitution type match the elements facets,
    # since typecode is created for a single existing pyobj
    # some facets are irrelevant.
    sub = _copy(sub)
    sub.nspname = typecode.nspname
    sub.pname = typecode.pname
    sub.aname = typecode.aname
    sub.minOccurs = sub.maxOccurs = 1
    return sub


        
class ComplexType(TypeCode):
    '''Represents an element of complexType, potentially containing other 
    elements.
    '''
    logger = _GetLogger('ZSI.TCcompound.ComplexType')
    
    def __init__(self, pyclass, ofwhat, pname=None, inorder=False, inline=False,
    mutable=True, mixed=False, mixed_aname='_text', **kw):
        '''pyclass -- the Python class to hold the fields
        ofwhat -- a list of fields to be in the complexType
        inorder -- fields must be in exact order or not
        inline -- don't href/id when serializing
        mutable -- object could change between multiple serializations
        type -- the (URI,localname) of the datatype
        mixed -- mixed content model? True/False
        mixed_aname -- if mixed is True, specify text content here. Default _text
        '''
        TypeCode.__init__(self, pname, pyclass=pyclass, **kw)
        self.inorder = inorder
        self.inline = inline
        self.mutable = mutable
        self.mixed = mixed
        self.mixed_aname = None
        if mixed is True:
            self.mixed_aname = mixed_aname

        if self.mutable is True: self.inline = True
        self.type = kw.get('type') or _get_xsitype(self)
        t = type(ofwhat)
        if t not in _seqtypes:
            raise TypeError(
                'Struct ofwhat must be list or sequence, not ' + str(t))
        self.ofwhat = tuple(ofwhat)
        if TypeCode.typechecks:
            # XXX Not sure how to determine if new-style class..
            if self.pyclass is not None and \
                type(self.pyclass) is not types.ClassType and not isinstance(self.pyclass, object):
                raise TypeError('pyclass must be None or an old-style/new-style class, not ' +
                        str(type(self.pyclass)))
            _check_typecode_list(self.ofwhat, 'ComplexType')

    def parse(self, elt, ps):
        debug = self.logger.debugOn()
        debug and self.logger.debug('parse')
        
        xtype = self.checkname(elt, ps)
        if self.type and xtype not in [ self.type, (None,None) ]:
            if not isinstance(self, TypeDefinition):
                raise EvaluateException(\
                    'ComplexType for %s has wrong type(%s), looking for %s' %
                        (self.pname, self.checktype(elt,ps), self.type), 
                                        ps.Backtrace(elt))
            else:
                #TODO: mabye change MRO to handle this 
                debug and self.logger.debug('delegate to substitute type')
                what = TypeDefinition.getSubstituteType(self, elt, ps)
                return what.parse(elt, ps)
            
        href = _find_href(elt)
        if href:
            if _children(elt):
                raise EvaluateException('Struct has content and HREF',
                        ps.Backtrace(elt))
            elt = ps.FindLocalHREF(href, elt)
        c = _child_elements(elt)
        count = len(c)
        if self.nilled(elt, ps): return Nilled

        # Create the object.
        v = {}

        # parse all attributes contained in attribute_typecode_dict (user-defined attributes),
        # the values (if not None) will be keyed in self.attributes dictionary.
        attributes = self.parse_attributes(elt, ps)
        if attributes:
            v[self.attrs_aname] = attributes

        #MIXED
        if self.mixed is True:
            v[self.mixed_aname] = self.simple_value(elt,ps, mixed=True)

        # Clone list of kids (we null it out as we process)
        c, crange = c[:], range(len(c))
        # Loop over all items we're expecting
        
        if debug:
            self.logger.debug("ofwhat: %s",str(self.ofwhat))
            
        any = None
        for i,what in [ (i, self.ofwhat[i]) for i in range(len(self.ofwhat)) ]:
            
            # retrieve typecode if it is hidden
            if callable(what): what = what()
            
            # Loop over all available kids
            if debug: 
                self.logger.debug("what: (%s,%s)", what.nspname, what.pname)
                
            for j,c_elt in [ (j, c[j]) for j in crange if c[j] ]:
                # Parse value, and mark this one done. 
                if debug:
                    self.logger.debug("child node: (%s,%s)", c_elt.namespaceURI, c_elt.tagName)

                match = False
                if what.name_match(c_elt):
                    match = True
                    value = what.parse(c_elt, ps)
                else:
                    # substitutionGroup head must be a global element declaration
                    # if successful delegate to matching GED
                    subwhat = _get_substitute_element(what, c_elt, ps)
                    if subwhat:
                        match = True
                        value = subwhat.parse(c_elt, ps)

                    if debug: 
                        self.logger.debug("substitutionGroup: %s", subwhat)

                if match:
                    if what.maxOccurs > 1:
                        if v.has_key(what.aname):
                            v[what.aname].append(value)
                        else:
                            v[what.aname] = [value]
                        c[j] = None
                        continue
                    else:
                        v[what.aname] = value
                    c[j] = None
                    break

                if debug:
                    self.logger.debug("no element (%s,%s)", what.nspname, what.pname)

                # No match; if it was supposed to be here, that's an error.
                if self.inorder is True and i == j:
                    raise EvaluateException('Out of order complexType',
                            ps.Backtrace(c_elt))
            else:
                # only supporting 1 <any> declaration in content.
                if isinstance(what,AnyElement):
                    any = what
                elif hasattr(what, 'default'):
                    v[what.aname] = what.default
                elif what.minOccurs > 0 and not v.has_key(what.aname):
                    raise EvaluateException('Element "' + what.aname + \
                        '" missing from complexType', ps.Backtrace(elt))

        # Look for wildcards and unprocessed children
        # XXX Stick all this stuff in "any", hope for no collisions
        if any is not None:
            occurs = 0
            v[any.aname] = []
            for j,c_elt in [ (j, c[j]) for j in crange if c[j] ]:
                value = any.parse(c_elt, ps)
                if any.maxOccurs == UNBOUNDED or any.maxOccurs > 1:
                    v[any.aname].append(value)
                else:
                    v[any.aname] = value

                occurs += 1

            # No such thing as nillable <any>
            if any.maxOccurs == 1 and occurs == 0:
                v[any.aname] = None
            elif occurs < any.minOccurs or (any.maxOccurs!=UNBOUNDED and any.maxOccurs<occurs):
                raise EvaluateException('occurances of <any> elements(#%d) bound by (%d,%s)' %(
                    occurs, any.minOccurs,str(any.maxOccurs)), ps.Backtrace(elt))

        if not self.pyclass: 
            return v

        # type definition must be informed of element tag (nspname,pname),
        # element declaration is initialized with a tag.
        try:
            pyobj = self.pyclass()
        except Exception, e:
            raise TypeError("Constructing element (%s,%s) with pyclass(%s), %s" \
                %(self.nspname, self.pname, self.pyclass.__name__, str(e)))
        for key in v.keys():
            setattr(pyobj, key, v[key])
        return pyobj

    def serialize(self, elt, sw, pyobj, inline=False, name=None, **kw):
        if inline or self.inline:
            self.cb(elt, sw, pyobj, name=name, **kw)
        else:
            objid = _get_idstr(pyobj)
            ns,n = self.get_name(name, objid)
            el = elt.createAppendElement(ns, n)
            el.setAttributeNS(None, 'href', "#%s" %objid)
            sw.AddCallback(self.cb, elt, sw, pyobj)

    def cb(self, elt, sw, pyobj, name=None, **kw):
        debug = self.logger.debugOn()
        if debug:
            self.logger.debug("cb: %s" %str(self.ofwhat))

        objid = _get_idstr(pyobj)
        ns,n = self.get_name(name, objid)
        if pyobj is None:
            if self.nillable is True:
                elem = elt.createAppendElement(ns, n)
                self.serialize_as_nil(elem)
                return
            raise EvaluateException, 'element(%s,%s) is not nillable(%s)' %(
                self.nspname,self.pname,self.nillable)

        if self.mutable is False and sw.Known(pyobj): 
            return
        
        if debug:
            self.logger.debug("element: (%s, %s)", str(ns), n)
            
        if n is not None:
            elem = elt.createAppendElement(ns, n)
            self.set_attributes(elem, pyobj)
            if kw.get('typed', self.typed) is True:
                self.set_attribute_xsi_type(elem)

            #MIXED For now just stick it in front.
            if self.mixed is True and self.mixed_aname is not None:
               if hasattr(pyobj, self.mixed_aname):
                   textContent = getattr(pyobj, self.mixed_aname)
                   if hasattr(textContent, 'typecode'):
                       textContent.typecode.serialize_text_node(elem, sw, textContent)
                   elif type(textContent) in _stringtypes:
                       if debug:
                           self.logger.debug("mixed text content:\n\t%s", 
                                             textContent)
                       elem.createAppendTextNode(textContent)
                   else:
                       raise EvaluateException('mixed test content in element (%s,%s) must be a string type' %(
                           self.nspname,self.pname), sw.Backtrace(elt))
               else:
                   if debug:
                       self.logger.debug("mixed NO text content in %s", 
                                         self.mixed_aname)
        else:
            #For information items w/o tagNames 
            #  ie. model groups,SOAP-ENC:Header
            elem = elt

        if self.inline:
            pass
        elif not self.inline and self.unique:
            raise EvaluateException('Not inline, but unique makes no sense. No href/id.',
                sw.Backtrace(elt))
        elif n is not None:
            self.set_attribute_id(elem, objid)

        if self.pyclass and type(self.pyclass) is type:
            f = lambda attr: getattr(pyobj, attr, None)
        elif self.pyclass:
            d = pyobj.__dict__
            f = lambda attr: d.get(attr)
        else:
            d = pyobj
            f = lambda attr: pyobj.get(attr)
            if TypeCode.typechecks and type(d) != types.DictType:
                raise TypeError("Classless complexType didn't get dictionary")

        indx, lenofwhat = 0, len(self.ofwhat)
        if debug:
            self.logger.debug('element declaration (%s,%s)', self.nspname, 
                              self.pname)
            if self.type:
                self.logger.debug('xsi:type definition (%s,%s)', self.type[0], 
                                  self.type[1])
            else:
                self.logger.warning('NO xsi:type')

        while indx < lenofwhat:
            occurs = 0
            what = self.ofwhat[indx]
            
            # retrieve typecode if hidden
            if callable(what): what = what()
            
            if debug:
                self.logger.debug('serialize what -- %s', 
                                  what.__class__.__name__)

            # No way to order <any> instances, so just grab any unmatched
            # anames and serialize them.  Only support one <any> in all content.
            # Must be self-describing instances

            # Regular handling of declared elements
            aname = what.aname
            v = f(aname)
            indx += 1
            if what.minOccurs == 0 and v is None: 
                continue

            # Default to typecode, if self-describing instance, and check 
            # to make sure it is derived from what.
            whatTC = what
            if whatTC.maxOccurs > 1 and v is not None:
                if type(v) not in _seqtypes:
                    raise EvaluateException('pyobj (%s,%s), aname "%s": maxOccurs %s, expecting a %s' %(
                         self.nspname,self.pname,what.aname,whatTC.maxOccurs,_seqtypes), 
                         sw.Backtrace(elt))

                for v2 in v: 
                    occurs += 1
                    if occurs > whatTC.maxOccurs:
                        raise EvaluateException('occurances (%d) exceeded maxOccurs(%d) for <%s>' %(
                                occurs, whatTC.maxOccurs, what.pname), 
                                sw.Backtrace(elt))
                        
                    what = _get_type_or_substitute(whatTC, v2, sw, elt)
                    if debug and what is not whatTC:
                        self.logger.debug('substitute derived type: %s' %
                                          what.__class__)
                        
                    what.serialize(elem, sw, v2, **kw)
#                    try:
#                        what.serialize(elem, sw, v2, **kw)
#                    except Exception, e:
#                        raise EvaluateException('Serializing %s.%s, %s %s' %
#                            (n, whatTC.aname or '?', e.__class__.__name__, str(e)))

                if occurs < whatTC.minOccurs:
                    raise EvaluateException(\
                        'occurances(%d) less than minOccurs(%d) for <%s>' %
                        (occurs, whatTC.minOccurs, what.pname), sw.Backtrace(elt))
                    
                continue

            if v is not None or what.nillable is True:
                what = _get_type_or_substitute(whatTC, v, sw, elt)
                if debug and what is not whatTC:
                    self.logger.debug('substitute derived type: %s' %
                                      what.__class__)
                what.serialize(elem, sw, v, **kw)
#                try:
#                    what.serialize(elem, sw, v, **kw)
#                except (ParseException, EvaluateException), e:
#                    raise
#                except Exception, e:
#                    raise EvaluateException('Serializing %s.%s, %s %s' %
#                        (n, whatTC.aname or '?', e.__class__.__name__, str(e)),
#                        sw.Backtrace(elt))
                continue

            raise EvaluateException('Got None for nillable(%s), minOccurs(%d) element (%s,%s), %s' %
                    (what.nillable, what.minOccurs, what.nspname, what.pname, elem),
                    sw.Backtrace(elt))


    def setDerivedTypeContents(self, extensions=None, restrictions=None):
        """For derived types set appropriate parameter and 
        """
        if extensions:
            ofwhat = list(self.ofwhat)
            if type(extensions) in _seqtypes:
                ofwhat += list(extensions)
            else:
                ofwhat.append(extensions)
        elif restrictions:
            if type(restrictions) in _seqtypes:
                ofwhat = restrictions
            else:
                ofwhat = (restrictions,)
        else:
            return
        self.ofwhat = tuple(ofwhat)
        self.lenofwhat = len(self.ofwhat)


class Struct(ComplexType):
    '''Struct is a complex type for accessors identified by name. 
       Constraint: No element may have the same name as any other,
       nor may any element have a maxOccurs > 1.
       
      <xs:group name="Struct" >
        <xs:sequence>
          <xs:any namespace="##any" minOccurs="0" maxOccurs="unbounded" processContents="lax" />
        </xs:sequence>
      </xs:group>

      <xs:complexType name="Struct" >
        <xs:group ref="tns:Struct" minOccurs="0" />
        <xs:attributeGroup ref="tns:commonAttributes"/>
      </xs:complexType> 
    '''
    logger = _GetLogger('ZSI.TCcompound.Struct')
    
    def __init__(self, pyclass, ofwhat, pname=None, inorder=False, inline=False,
        mutable=True, **kw):
        '''pyclass -- the Python class to hold the fields
        ofwhat -- a list of fields to be in the struct
        inorder -- fields must be in exact order or not
        inline -- don't href/id when serializing
        mutable -- object could change between multiple serializations
        '''
        ComplexType.__init__(self, pyclass, ofwhat, pname=pname, 
            inorder=inorder, inline=inline, mutable=mutable, 
            **kw
            )
        
        # Check Constraints
        whats = map(lambda what: (what.nspname,what.pname), self.ofwhat)
        for idx in range(len(self.ofwhat)):
            what = self.ofwhat[idx]
            key = (what.nspname,what.pname)
            if not isinstance(what, AnyElement) and what.maxOccurs > 1:
                raise TypeError,\
                    'Constraint: no element can have a maxOccurs>1'
            if key in whats[idx+1:]:
                raise TypeError,\
                    'Constraint: No element may have the same name as any other'


class Array(TypeCode):
    '''An array.
        atype -- arrayType, (namespace,ncname) 
        mutable -- object could change between multiple serializations
        undeclared -- do not serialize/parse arrayType attribute.
    '''
    logger = _GetLogger('ZSI.TCcompound.Array')
    
    def __init__(self, atype, ofwhat, pname=None, dimensions=1, fill=None,
    sparse=False, mutable=False, size=None, nooffset=0, undeclared=False,
    childnames=None, **kw):
        TypeCode.__init__(self, pname, **kw)
        self.dimensions = dimensions
        self.atype = atype
        if undeclared is False and self.atype[1].endswith(']') is False:
            self.atype = (self.atype[0], '%s[]' %self.atype[1])
        # Support multiple dimensions
        if self.dimensions != 1:
            raise TypeError("Only single-dimensioned arrays supported")
        self.fill = fill
        self.sparse = sparse
        #if self.sparse: ofwhat.minOccurs = 0
        self.mutable = mutable
        self.size = size
        self.nooffset = nooffset
        self.undeclared = undeclared
        self.childnames = childnames
        if self.size:
            t = type(self.size)
            if t in _inttypes:
                self.size = (self.size,)
            elif t in _seqtypes:
                self.size = tuple(self.size)
            elif TypeCode.typechecks:
                raise TypeError('Size must be integer or list, not ' + str(t))

        # by default use Any
        ofwhat = ofwhat or Any()

        if TypeCode.typechecks:
            if self.undeclared is False and type(atype) not in _seqtypes and len(atype) == 2:
                raise TypeError("Array type must be a sequence of len 2.")
            t = type(ofwhat)
            if not isinstance(ofwhat, TypeCode):
                raise TypeError(
                    'Array ofwhat outside the TypeCode hierarchy, ' +
                    str(ofwhat.__class__))
            if self.size:
                if len(self.size) != self.dimensions:
                    raise TypeError('Array dimension/size mismatch')
                for s in self.size:
                    if type(s) not in _inttypes:
                        raise TypeError('Array size "' + str(s) +
                                '" is not an integer.')
        self.ofwhat = ofwhat

    def parse_offset(self, elt, ps):
        o = _find_arrayoffset(elt)
        if not o: return 0
        if not _offset_pat.match(o):
            raise EvaluateException('Bad offset "' + o + '"',
                        ps.Backtrace(elt))
        return int(o[1:-1])

    def parse_position(self, elt, ps):
        o = _find_arrayposition(elt)
        if not o: return None
        if o.find(',') > -1:
            raise EvaluateException('Sorry, no multi-dimensional arrays',
                    ps.Backtrace(elt))
        if not _position_pat.match(o):
            raise EvaluateException('Bad array position "' + o + '"',
                    ps.Backtrace(elt))
        return int(o[1:-1])

    def parse(self, elt, ps):
        href = _find_href(elt)
        if href:
            if _children(elt):
                raise EvaluateException('Array has content and HREF',
                        ps.Backtrace(elt))
            elt = ps.FindLocalHREF(href, elt)
        if self.nilled(elt, ps): return Nilled
        if not _find_arraytype(elt) and self.undeclared is False:
            raise EvaluateException('Array expected', ps.Backtrace(elt))
        t = _find_type(elt)
        if t:
            pass # XXX should check the type, but parsing that is hairy.
        offset = self.parse_offset(elt, ps)
        v, vlen = [], 0
        if offset and not self.sparse:
            while vlen < offset:
                vlen += 1
                v.append(self.fill)
        for c in _child_elements(elt):
            item = self.ofwhat.parse(c, ps)
            position = self.parse_position(c, ps) or offset
            if self.sparse:
                v.append((position, item))
            else:
                while offset < position:
                    offset += 1
                    v.append(self.fill)
                v.append(item)
            offset += 1
        return v

    def serialize(self, elt, sw, pyobj, name=None, childnames=None, **kw):
        debug = self.logger.debugOn()
        if debug:
            self.logger.debug("serialize: %r" %pyobj)
        
        if self.mutable is False and sw.Known(pyobj): return
        objid = _get_idstr(pyobj)
        ns,n = self.get_name(name, objid)
        el = elt.createAppendElement(ns, n)

        # nillable
        if self.nillable is True and pyobj is None:
            self.serialize_as_nil(el)
            return None

        # other attributes
        self.set_attributes(el, pyobj)

        # soap href attribute
        unique = self.unique or kw.get('unique', False)
        if unique is False and sw.Known(pyobj):
            self.set_attribute_href(el, objid)
            return None

        # xsi:type attribute 
        if kw.get('typed', self.typed) is True:
            self.set_attribute_xsi_type(el, **kw)

        # soap id attribute
        if self.unique is False:
            self.set_attribute_id(el, objid)

        offset = 0
        if self.sparse is False and self.nooffset is False:
            offset, end = 0, len(pyobj)
            while offset < end and pyobj[offset] == self.fill:
                offset += 1
            if offset: 
                el.setAttributeNS(SOAP.ENC, 'offset', '[%d]' %offset)

        if self.undeclared is False:
            el.setAttributeNS(SOAP.ENC, 'arrayType', 
                '%s:%s' %(el.getPrefix(self.atype[0]), self.atype[1])
            )

        if debug:
            self.logger.debug("ofwhat: %r" %self.ofwhat)

        d = {}
        kn = childnames or self.childnames
        if kn:
            d['name'] = kn
        elif not self.ofwhat.aname:
            d['name'] = 'element'
            
        if self.sparse is False:
            for e in pyobj[offset:]: self.ofwhat.serialize(el, sw, e, **d)
        else:
            position = 0
            for pos, v in pyobj:
                if pos != position:
                    el.setAttributeNS(SOAP.ENC, 'position', '[%d]' %pos)
                    position = pos

                self.ofwhat.serialize(el, sw, v, **d)
                position += 1


if __name__ == '__main__': print _copyright
