# Copyright (c) 2001 Zope Corporation and Contributors. All Rights Reserved.
#
# This software is subject to the provisions of the Zope Public License,
# Version 2.0 (ZPL).  A copy of the ZPL should accompany this distribution.
# THIS SOFTWARE IS PROVIDED "AS IS" AND ANY AND ALL EXPRESS OR IMPLIED
# WARRANTIES ARE DISCLAIMED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF TITLE, MERCHANTABILITY, AGAINST INFRINGEMENT, AND FITNESS
# FOR A PARTICULAR PURPOSE.

import weakref, re, os, sys
from ConfigParser import SafeConfigParser as ConfigParser,\
    NoSectionError, NoOptionError
from urlparse import urlparse

from ZSI import TC
from ZSI.client import _Binding
from ZSI.generate import commands,containers
from ZSI.schema import GED, GTD

import wstools


#url_to_mod = re.compile(r'<([^ \t\n\r\f\v:]+:)?include\s+location\s*=\s*"(\S+)"')
def _urn_to_module(urn): return '%s_types' %re.sub(_urn_to_module.regex, '_', urn)
_urn_to_module.regex = re.compile(r'[\W]')
    

class ServiceProxy:
    """A ServiceProxy provides a convenient way to call a remote web
       service that is described with WSDL. The proxy exposes methods
       that reflect the methods of the remote web service."""

    def __init__(self, wsdl, url=None, service=None, port=None,
                 cachedir=os.path.join(os.path.expanduser('~'), '.zsi_service_proxy_dir'), 
                 asdict=True, lazy=False, pyclass=False, force=False, **kw):
        """
        Parameters:
           wsdl -- URL of WSDL.
           url -- override WSDL SOAP address location
           service -- service name or index
           port -- port name or index
           cachedir -- where to store generated files
           asdict -- use dicts, else use generated pyclass
           lazy -- use lazy typecode evaluation
           pyclass -- use pyclass_type metaclass adds properties, "new_", "set_, 
               "get_" methods for schema element and attribute declarations.
           force -- regenerate all WSDL code, write over cache.
           
        NOTE: all other **kw will be passed to the underlying 
        ZSI.client._Binding constructor.
           
        """
        self._asdict = asdict
        
        # client._Binding
        self._url = url
        self._kw = kw
        
        # WSDL
        self._wsdl = wstools.WSDLTools.WSDLReader().loadFromURL(wsdl)
        self._service = self._wsdl.services[service or 0]
        self.__doc__ = self._service.documentation
        self._port = self._service.ports[port or 0]
        self._name = self._service.name
        self._methods = {}
        self._cachedir = cachedir
        self._lazy = lazy
        self._pyclass = pyclass
        self._force = force
        
        # Set up rpc methods for service/port
        port = self._port
        binding = port.getBinding()
        portType = binding.getPortType()
        for port in self._service.ports:
            for item in port.getPortType().operations:
                try:
                    callinfo = wstools.WSDLTools.callInfoFromWSDL(port, item.name)
                except:
                    # ignore non soap-1.1  bindings
                    continue
                
                method = MethodProxy(self, callinfo)
                setattr(self, item.name, method)
                self._methods.setdefault(item.name, []).append(method)
       
        self._mod = self._load(wsdl)
        
    def _load(self, location):
        """
        location -- URL or file location
        isxsd -- is this a xsd file?
        """
        cachedir = self._cachedir
        # wsdl2py: deal with XML Schema
        if not os.path.isdir(cachedir): os.mkdir(cachedir)
    
        file = os.path.join(cachedir, '.cache')
        section = 'TYPES'
        cp = ConfigParser()
        try:
            cp.readfp(open(file, 'r'))
        except IOError:
            del cp;  cp = None
            
        option = location.replace(':', '-') # colons seem to screw up option
        if (not self._force and cp is not None and cp.has_section(section) and 
            cp.has_option(section, option)):
            types = cp.get(section, option)
        else:
            # dont do anything to anames
            if not self._pyclass:
                containers.ContainerBase.func_aname = lambda instnc,n: str(n)
                
            args = ['-o', cachedir, location]
            if self._lazy: args.insert(0, '-l')
            if self._pyclass: args.insert(0, '-b')
            files = commands.wsdl2py(args)

            if cp is None: cp = ConfigParser()
            if not cp.has_section(section): cp.add_section(section)
            types = filter(lambda f: f.endswith('_types.py'), files)[0]
            cp.set(section, option, types)
            cp.write(open(file, 'w'))
            
        if os.path.abspath(cachedir) not in sys.path:
            sys.path.append(os.path.abspath(cachedir))
            
        mod = os.path.split(types)[-1].rstrip('.py')
        return __import__(mod)
    
    def _load_schema(self, location, xml=None):
        """
        location -- location of schema, also used as a key
        xml -- optional string representation of schema
        """
        cachedir = self._cachedir
        # wsdl2py: deal with XML Schema
        if not os.path.isdir(cachedir): os.mkdir(cachedir)
    
        file = os.path.join(cachedir, '.cache')
        section = 'TYPES'
        cp = ConfigParser()
        try:
            cp.readfp(open(file, 'r'))
        except IOError:
            del cp;  cp = None
            
        option = location.replace(':', '-') # colons seem to screw up option
        if (cp is not None and cp.has_section(section) and 
            cp.has_option(section, option)):
            types = cp.get(section, option)
        else:
            # dont do anything to anames
            if not self._pyclass:
                containers.ContainerBase.func_aname = lambda instnc,n: str(n)
                
            from ZSI.wstools import XMLSchema
            reader = XMLSchema.SchemaReader(base_url=location)
            if xml is not None and isinstance(xml, basestring):
                schema = reader.loadFromString(xml)
            elif xml is not None:
                raise RuntimeError, 'Unsupported: XML must be string'
            elif not os.path.isfile(location):
                schema = reader.loadFromURL(location)
            else:
                schema = reader.reader.loadFromFile(location)
                
            # TODO: change this to keyword list
            class options:
                output_dir = cachedir
                schema = True
                simple_naming = False
                address = False
                lazy = self._lazy
                complexType = self._pyclass

            schema.location = location
            files = commands._wsdl2py(options, schema)
            if cp is None: cp = ConfigParser()
            if not cp.has_section(section): cp.add_section(section)
            types = filter(lambda f: f.endswith('_types.py'), files)[0]
            cp.set(section, option, types)
            cp.write(open(file, 'w'))
            
        if os.path.abspath(cachedir) not in sys.path:
            sys.path.append(os.path.abspath(cachedir))
            
        mod = os.path.split(types)[-1].rstrip('.py')
        return __import__(mod)
            
    def _call(self, name, soapheaders):
        """return the Call to the named remote web service method.
        closure used to prevent multiple values for name and soapheaders 
        parameters 
        """
        
        def call_closure(*args, **kwargs):
            """Call the named remote web service method."""
            if len(args) and len(kwargs):
                raise TypeError, 'Use positional or keyword argument only.'
                
            if len(args) > 0:
                raise TypeError, 'Not supporting SOAPENC:Arrays or XSD:List'
            
            if len(kwargs): 
                args = kwargs

            callinfo = getattr(self, name).callinfo
    
            # go through the list of defined methods, and look for the one with
            # the same number of arguments as what was passed.  this is a weak
            # check that should probably be improved in the future to check the
            # types of the arguments to allow for polymorphism
            for method in self._methods[name]:
                if len(method.callinfo.inparams) == len(kwargs):
                    callinfo = method.callinfo
    
            binding = _Binding(url=self._url or callinfo.location,
                              soapaction=callinfo.soapAction,
                              **self._kw)
    
            kw = dict(unique=True)
            if callinfo.use == 'encoded':
                kw['unique'] = False
    
            if callinfo.style == 'rpc':
                request = TC.Struct(None, ofwhat=[], 
                                 pname=(callinfo.namespace, name), **kw)
                
                response = TC.Struct(None, ofwhat=[], 
                                 pname=(callinfo.namespace, name+"Response"), **kw)
                
                if len(callinfo.getInParameters()) != len(args):
                    raise RuntimeError('expecting "%s" parts, got %s' %(
                           str(callinfo.getInParameters(), str(args))))
                
                for msg,pms in ((request,callinfo.getInParameters()), 
                                (response,callinfo.getOutParameters())):
                    msg.ofwhat = []
                    for part in pms:
                        klass = GTD(*part.type)
                        if klass is None:
                            if part.type:
                                klass = filter(lambda gt: part.type==gt.type,TC.TYPES)
                                if len(klass) == 0:
                                    klass = filter(lambda gt: part.type[1]==gt.type[1],TC.TYPES)
                                    if not len(klass):klass = [TC.Any]
                                if len(klass) > 1: #Enumerations, XMLString, etc
                                    klass = filter(lambda i: i.__dict__.has_key('type'), klass)
                                klass = klass[0]
                            else:
                                klass = TC.Any
                    
                        msg.ofwhat.append(klass(part.name))
                        
                    msg.ofwhat = tuple(msg.ofwhat)
                if not args: args = {}
            else:
                # Grab <part element> attribute
                ipart,opart = callinfo.getInParameters(),callinfo.getOutParameters()
                if ( len(ipart) != 1 or not ipart[0].element_type or 
                    ipart[0].type is None ):
                    raise RuntimeError, 'Bad Input Message "%s"' %callinfo.name
        
                if ( len(opart) not in (0,1) or not opart[0].element_type or 
                    opart[0].type is None ):
                    raise RuntimeError, 'Bad Output Message "%s"' %callinfo.name
                
#                if ( len(args) > 1 ):
#                    raise RuntimeError, 'Message has only one part:  %s' %str(args)
                
                ipart = ipart[0]
                request,response = GED(*ipart.type),None
                if opart: response = GED(*opart[0].type)
    
            msg = args
            if self._asdict: 
                if not msg: msg = dict()
                self._nullpyclass(request)
            elif request.pyclass is not None:
                if type(args) is dict:
                    msg = request.pyclass()
                    msg.__dict__.update(args)
                elif type(args) is list and len(args) == 1: 
                    msg = request.pyclass(args[0])
                else: 
                    msg = request.pyclass()
                    
            binding.Send(None, None, msg,
                         requesttypecode=request,
                         soapheaders=soapheaders,
                         encodingStyle=callinfo.encodingStyle)
            
            if response is None: 
                return None
            
            if self._asdict: self._nullpyclass(response)
            return binding.Receive(replytype=response,
                         encodingStyle=callinfo.encodingStyle)
            
        return call_closure

    def _nullpyclass(cls, typecode):
        typecode.pyclass = None
        if not hasattr(typecode, 'ofwhat'): return
        if type(typecode.ofwhat) not in (list,tuple): #Array
            cls._nullpyclass(typecode.ofwhat)
        else: #Struct/ComplexType
            for i in typecode.ofwhat: cls._nullpyclass(i)    
    _nullpyclass = classmethod(_nullpyclass)


class MethodProxy:
    """ """
    def __init__(self, parent, callinfo):
        self.__name__ = callinfo.methodName
        self.__doc__ = callinfo.documentation
        self.callinfo = callinfo
        self.parent = weakref.ref(parent)
        self.soapheaders = []

    def __call__(self, *args, **kwargs):
        return self.parent()._call(self.__name__, self.soapheaders)(*args, **kwargs)

    def add_headers(self, **headers):
        """packing dicts into typecode pyclass, may fail if typecodes are
        used in the body (when asdict=True)
        """
        class _holder: pass
        def _remap(pyobj, **d):
            pyobj.__dict__ = d
            for k,v in pyobj.__dict__.items():
                if type(v) is not dict: continue
                pyobj.__dict__[k] = p = _holder()
                _remap(p, **v)
    
        for k,v in headers.items():
            h = filter(lambda i: k in i.type, self.callinfo.inheaders)[0]
            if h.element_type != 1: 
                raise RuntimeError, 'not implemented'

            typecode = GED(*h.type)
            if typecode is None: 
                raise RuntimeError, 'no matching element for %s' %str(h.type)

            pyclass = typecode.pyclass
            if pyclass is None: 
                raise RuntimeError, 'no pyclass for typecode %s' %str(h.type)

            if type(v) is not dict:
                pyobj = pyclass(v)
            else:
                pyobj = pyclass()
                _remap(pyobj, **v)

            self.soapheaders.append(pyobj)
