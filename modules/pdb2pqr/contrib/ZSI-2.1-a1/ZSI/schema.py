#! /usr/bin/env python
# $Header$
'''XML Schema support
'''

from ZSI import _copyright, _seqtypes, _find_type, _get_element_nsuri_name, EvaluateException
from ZSI.wstools.Namespaces import SCHEMA, SOAP
from ZSI.wstools.Utility import SplitQName


def _get_type_definition(namespaceURI, name, **kw):
    return SchemaInstanceType.getTypeDefinition(namespaceURI, name, **kw)

def _get_global_element_declaration(namespaceURI, name, **kw):
    return SchemaInstanceType.getElementDeclaration(namespaceURI, name, **kw)

def _get_substitute_element(head, elt, ps):
    '''if elt matches a member of the head substitutionGroup, return 
    the GED typecode.

    head -- ElementDeclaration typecode, 
    elt -- the DOM element being parsed
    ps -- ParsedSoap Instance
    '''
    if not isinstance(head, ElementDeclaration):
        return None

    return ElementDeclaration.getSubstitutionElement(head, elt, ps)

def _has_type_definition(namespaceURI, name):
    return SchemaInstanceType.getTypeDefinition(namespaceURI, name) is not None

def _is_substitute_element(head, sub):
    '''if head and sub are both GEDs, and sub declares 
    head as its substitutionGroup then return True.

    head -- Typecode instance
    sub  -- Typecode instance
    '''
    if not isinstance(head, ElementDeclaration) or not isinstance(sub, ElementDeclaration):
        return False

    try:
        group = sub.substitutionGroup 
    except (AttributeError, TypeError):
        return False

    ged = GED(*group)

    # TODO: better way of representing element references.  Wrap them with
    # facets, and dereference when needed and delegate to..
    print (head.nspname == ged.nspname and head.pname == ged.pname)
    if head is ged or not (head.nspname == ged.nspname and head.pname == ged.pname):
        return False

    return True

#
# functions for retrieving schema items from 
# the global schema instance.
#
GED = _get_global_element_declaration
GTD = _get_type_definition


def WrapImmutable(pyobj, what):
    '''Wrap immutable instance so a typecode can be
    set, making it self-describing ie. serializable.
    '''
    return _GetPyobjWrapper.WrapImmutable(pyobj, what)

def RegisterBuiltin(arg):
    '''Add a builtin to be registered, and register it
    with the Any typecode.
    '''
    _GetPyobjWrapper.RegisterBuiltin(arg)
    _GetPyobjWrapper.RegisterAnyElement()

def RegisterAnyElement():
    '''register all Wrapper classes with the Any typecode.
    This allows instances returned by Any to be self-describing.
    ie. serializable.  AnyElement falls back on Any to parse
    anything it doesn't understand.
    '''
    return _GetPyobjWrapper.RegisterAnyElement()


class SchemaInstanceType(type):
    '''Register all types/elements, when hit already defined 
    class dont create a new one just give back reference.  Thus 
    import order determines which class is loaded.

    class variables:
        types -- dict of typecode classes definitions 
            representing global type definitions.
        elements -- dict of typecode classes representing 
            global element declarations.
        element_typecode_cache -- dict of typecode instances 
            representing global element declarations.
    '''
    types = {}
    elements = {}
    element_typecode_cache = {}
    #substitution_registry = {}
    
    def __new__(cls,classname,bases,classdict):
        '''If classdict has literal and schema register it as a
        element declaration, else if has type and schema register
        it as a type definition.
        '''
        if classname in ['ElementDeclaration', 'TypeDefinition', 'LocalElementDeclaration',]:
            return type.__new__(cls,classname,bases,classdict)

        if ElementDeclaration in bases:
            if classdict.has_key('schema') is False  or classdict.has_key('literal') is False: 
                raise AttributeError, 'ElementDeclaration must define schema and literal attributes'

            key = (classdict['schema'],classdict['literal'])
            if SchemaInstanceType.elements.has_key(key):
                return SchemaInstanceType.elements[key]

            # create global element declaration
            ged = SchemaInstanceType.elements[key] = type.__new__(cls,classname,bases,classdict)

            # TODO: Maybe I want access to all registrants??
            # 
            #if classdict.has_key('substitutionGroup'):
            #    sub = classdict.has_key('substitutionGroup')
            #    if not SchemaInstanceType.substitution_registry.has_key(sub):
            #        SchemaInstanceType.substitution_registry[sub] = [ged]
            #    else:
            #        SchemaInstanceType.substitution_registry[sub].append(ged)

            return ged

        if TypeDefinition in bases:
            if classdict.has_key('type') is None:
                raise AttributeError, 'TypeDefinition must define type attribute'

            key = classdict['type']
            if SchemaInstanceType.types.has_key(key) is False:
                SchemaInstanceType.types[key] = type.__new__(cls,classname,bases,classdict)
            return SchemaInstanceType.types[key]

        if LocalElementDeclaration in bases:
                return type.__new__(cls,classname,bases,classdict)

        raise TypeError, 'SchemaInstanceType must be an ElementDeclaration or TypeDefinition '

    def getTypeDefinition(cls, namespaceURI, name, lazy=False):
        '''Grab a type definition, returns a typecode class definition
        because the facets (name, minOccurs, maxOccurs) must be provided.
 
        Parameters:
           namespaceURI -- 
           name -- 
        '''
        klass = cls.types.get((namespaceURI, name), None)
        if lazy and klass is not None:
            return _Mirage(klass)
        return klass
    getTypeDefinition = classmethod(getTypeDefinition)

    def getElementDeclaration(cls, namespaceURI, name, isref=False, lazy=False):
        '''Grab an element declaration, returns a typecode instance
        representation or a typecode class definition.  An element 
        reference has its own facets, and is local so it will not be
        cached.

        Parameters:
            namespaceURI -- 
            name -- 
            isref -- if element reference, return class definition.
        '''
        key = (namespaceURI, name)
        if isref:
            klass = cls.elements.get(key,None)
            if klass is not None and lazy is True:
                return _Mirage(klass)
            return klass
 
        typecode = cls.element_typecode_cache.get(key, None)
        if typecode is None:
            tcls = cls.elements.get(key,None)
            if tcls is not None:
                typecode = cls.element_typecode_cache[key] = tcls()
                typecode.typed = False
            
        return typecode
    getElementDeclaration = classmethod(getElementDeclaration)


class ElementDeclaration:
    '''Typecodes subclass to represent a Global Element Declaration by
    setting class variables schema and literal.

    schema = namespaceURI
    literal = NCName
    substitutionGroup -- GED reference of form, (namespaceURI,NCName)
    '''
    __metaclass__ = SchemaInstanceType

    def checkSubstitute(self, typecode):
        '''If this is True, allow typecode to be substituted
        for "self" typecode.
        '''
        if not isinstance(typecode, ElementDeclaration): 
            return False

        try:
            nsuri,ncname = typecode.substitutionGroup
        except (AttributeError, TypeError):
            return False

        if (nsuri,ncname) != (self.schema,self.literal):
            # allow slop with the empty namespace 
            if not nsuri and not self.schema and ncname == self.literal:
                 return True

            return False

        sub = GED(self.schema, self.literal)
        if sub is None or sub is not typecode:
            return False

        return True

    def getSubstitutionElement(self, elt, ps):
        '''if elt matches a member of the head substitutionGroup, return 
        the GED typecode representation of the member.

        head -- ElementDeclaration typecode, 
        elt -- the DOM element being parsed
        ps -- ParsedSoap instance
        '''
        nsuri,ncname = _get_element_nsuri_name(elt)
        typecode = GED(nsuri,ncname)
        if typecode is None:
            return

        try:
            nsuri,ncname = typecode.substitutionGroup
        except (AttributeError, TypeError):
            return

        if (ncname == self.pname) and (nsuri == self.nspname or 
           (not nsuri and not self.nspname)):
             return typecode
       
        return 
 
 
class LocalElementDeclaration:
    '''Typecodes subclass to represent a Local Element Declaration.
    '''
    __metaclass__ = SchemaInstanceType
    

class TypeDefinition:
    '''Typecodes subclass to represent a Global Type Definition by
    setting class variable type.

    type = (namespaceURI, NCName)
    '''
    __metaclass__ = SchemaInstanceType
    
    def getSubstituteType(self, elt, ps):
        '''if xsi:type does not match the instance type attr,
        check to see if it is a derived type substitution.
        
        DONT Return the element's type.
        
        Parameters:
            elt -- the DOM element being parsed
            ps -- the ParsedSoap object.
        '''
        pyclass = SchemaInstanceType.getTypeDefinition(*self.type)
        if pyclass is None:
            raise EvaluateException(
                    'No Type registed for xsi:type=(%s, %s)' %
                    (self.type[0], self.type[1]), ps.Backtrace(elt))
            
        typeName = _find_type(elt)
        prefix,typeName = SplitQName(typeName)
        uri = ps.GetElementNSdict(elt).get(prefix)
        subclass = SchemaInstanceType.getTypeDefinition(uri, typeName)
        if subclass is None:
            raise EvaluateException(
                    'No registered xsi:type=(%s, %s), substitute for xsi:type=(%s, %s)' %
                    (uri, typeName, self.type[0], self.type[1]), ps.Backtrace(elt))
                    
        if not issubclass(subclass, pyclass) and subclass(None) and not issubclass(subclass, pyclass):
            raise TypeError(
                    'Substitute Type (%s, %s) is not derived from %s' %
                    (self.type[0], self.type[1], pyclass), ps.Backtrace(elt))

        return subclass((self.nspname, self.pname))
    
    

class _Mirage:
    '''Used with SchemaInstanceType for lazy evaluation, eval during serialize or 
    parse as needed.  Mirage is callable, TypeCodes are not.  When called it returns the
    typecode.  Tightly coupled with generated code.
    
    NOTE: **Must Use ClassType** for intended MRO of __call__ since setting it in
    an instance attribute rather than a class attribute (will not work for object).
    '''
    def __init__(self, klass):
        self.klass = klass
        self.__reveal = False
        self.__cache = None
        if issubclass(klass, ElementDeclaration):
            self.__call__ = self._hide_element
            
    def __str__(self):
        msg = "<Mirage id=%s, Local Element %s>"
        if issubclass(self.klass, ElementDeclaration):
            msg = "<Mirage id=%s, GED %s>"
        return  msg %(id(self), self.klass)
        
    def _hide_type(self, pname, aname, minOccurs=0, maxOccurs=1, nillable=False, 
                   **kw):
        self.__call__ = self._reveal_type
        self.__reveal = True
        
        # store all attributes, make some visable for pyclass_type
        self.__kw = kw
        self.minOccurs,self.maxOccurs,self.nillable = minOccurs,maxOccurs,nillable
        self.nspname,self.pname,self.aname = None,pname,aname
        if type(self.pname) in (tuple,list):
            self.nspname,self.pname = pname
        
        return self
        
    def _hide_element(self, minOccurs=0, maxOccurs=1, nillable=False, **kw):
        self.__call__ = self._reveal_element
        self.__reveal = True
        
        # store all attributes, make some visable for pyclass_type
        self.__kw = kw
        self.nspname = self.klass.schema
        self.pname = self.klass.literal
        #TODO: Fix hack
        #self.aname = '_%s' %self.pname
        self.minOccurs,self.maxOccurs,self.nillable = minOccurs,maxOccurs,nillable
        
        return self
    
    def _reveal_type(self):
        if self.__cache is None:
            self.__cache = self.klass(pname=self.pname, 
                            aname=self.aname, minOccurs=self.minOccurs, 
                            maxOccurs=self.maxOccurs, nillable=self.nillable, 
                            **self.__kw)
        return self.__cache
        
    def _reveal_element(self):
        if self.__cache is None:
            self.__cache = self.klass(minOccurs=self.minOccurs, 
                            maxOccurs=self.maxOccurs, nillable=self.nillable, 
                            **self.__kw)
        return self.__cache
    
    __call__ = _hide_type


class _GetPyobjWrapper:
    '''Get a python object that wraps data and typecode.  Used by
    <any> parse routine, so that typecode information discovered
    during parsing is retained in the pyobj representation
    and thus can be serialized.
    '''
    types_dict = dict()

    def RegisterBuiltin(cls, arg):
        '''register a builtin, create a new wrapper.
        '''
        if arg in cls.types_dict:
            raise RuntimeError, '%s already registered' %arg
        class _Wrapper(arg):
            'Wrapper for builtin %s\n%s' %(arg, cls.__doc__)
        _Wrapper.__name__ = '_%sWrapper' %arg.__name__
        cls.types_dict[arg] = _Wrapper
    RegisterBuiltin = classmethod(RegisterBuiltin)
        
    def RegisterAnyElement(cls):
        '''If find registered TypeCode instance, add Wrapper class 
        to TypeCode class serialmap and Re-RegisterType.  Provides
        Any serialzation of any instances of the Wrapper.
        '''
        for k,v in cls.types_dict.items():
            what = Any.serialmap.get(k)
            if what is None: continue
            if v in what.__class__.seriallist: continue
            what.__class__.seriallist.append(v)
            RegisterType(what.__class__, clobber=1, **what.__dict__)
    RegisterAnyElement = classmethod(RegisterAnyElement)

    def WrapImmutable(cls, pyobj, what):
        '''return a wrapper for pyobj, with typecode attribute set.
        Parameters:
            pyobj -- instance of builtin type (immutable)
            what -- typecode describing the data
        '''
        d = cls.types_dict
        if type(pyobj) is bool:  
            pyclass = d[int]
        elif d.has_key(type(pyobj)) is True:
            pyclass = d[type(pyobj)]
        else:
            raise TypeError,\
               'Expecting a built-in type in %s (got %s).' %(
                d.keys(),type(pyobj))

        newobj = pyclass(pyobj)
        newobj.typecode = what
        return newobj
    WrapImmutable = classmethod(WrapImmutable)
    

from TC import Any, RegisterType

if __name__ == '__main__': print _copyright

