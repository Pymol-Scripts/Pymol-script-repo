############################################################################
# Monte M. Goode, LBNL
# See LBNLCopyright for copyright notice!
###########################################################################

# main generator engine for new generation generator

# $Id: wsdl2python.py 1402 2007-07-06 22:51:32Z boverhof $

import os, sys, warnings
from ZSI import _get_idstr
from ZSI.wstools.logging import getLogger as _GetLogger
from ZSI.wstools import WSDLTools
from ZSI.wstools.WSDLTools import SoapAddressBinding,\
    SoapBodyBinding, SoapBinding,MimeContentBinding,\
    HttpUrlEncodedBinding
from ZSI.wstools.XMLSchema import SchemaReader, ElementDeclaration, SchemaError
from ZSI.typeinterpreter import BaseTypeInterpreter
from ZSI.generate import WsdlGeneratorError, Wsdl2PythonError
from containers import *
from ZSI.generate import utility
from ZSI.generate.utility import NamespaceAliasDict as NAD
from ZSI.generate.utility import GetModuleBaseNameFromWSDL

"""
classes:
    WriteServiceModule 
    -- composes/writes out client stubs and types module.

    ServiceDescription
    -- represents a single WSDL service.

    MessageWriter
    -- represents a single WSDL Message and associated bindings
    of the port/binding.

    SchemaDescription 
    -- generates classes for defs and decs in the schema instance.

    TypeWriter
    -- represents a type definition.

    ElementWriter
    -- represents a element declaration.

"""

class WriteServiceModule:
    """top level driver class invoked by wsd2py
    class variables:
        client_module_suffix -- suffix of client module.
        types_module_suffix -- suffix of types module.
    """
    client_module_suffix = '_client'
    messages_module_suffix = '_messages'
    types_module_suffix = '_types'
    logger = _GetLogger("WriteServiceModule")
    
    def __init__(self, wsdl, addressing=False, notification=False,
                 do_extended=False, extPyClasses=None, configParser = None):
        self._wsdl = wsdl
        self._addressing = addressing
        self._notification = notification
        self._configParser = configParser
        self.usedNamespaces = None
        self.services = []
        self.client_module_path = None
        self.types_module_name = None
        self.types_module_path = None
        self.messages_module_path = None # used in extended generation
        self.do_extended = do_extended
        self.extPyClasses = extPyClasses
       
    def getClientModuleName(self):
        """client module name.
        """
        name = GetModuleBaseNameFromWSDL(self._wsdl)
        if not name:
            raise WsdlGeneratorError, 'could not determine a service name'
        
        if self.client_module_suffix is None:
            return name

        return '%s%s' %(name, self.client_module_suffix)

#    def getMessagesModuleName(self):
#        name = GetModuleBaseNameFromWSDL(self._wsdl)
#        if not name:
#            raise WsdlGeneratorError, 'could not determine a service name'
#        
#        if self.messages_module_suffix is None:
#            return name
#
#        if len(self.messages_module_suffix) == 0:
#            return self.getClientModuleName()
#
#        return '%s%s' %(name, self.messages_module_suffix)

    def setTypesModuleName(self, name):
        self.types_module_name = name

    def getTypesModuleName(self):
        """types module name.
        """
        if self.types_module_name is not None:
            return self.types_module_name

        name = GetModuleBaseNameFromWSDL(self._wsdl)
        if not name:
            raise WsdlGeneratorError, 'could not determine a service name'
        
        if self.types_module_suffix is None:
            return name

        return '%s%s' %(name, self.types_module_suffix)

    def setClientModulePath(self, path):
        """setup module path to where client module before calling fromWsdl.
        module path to types module eg. MyApp.client
        """
        self.client_module_path = path

    def getTypesModulePath(self):
        """module path to types module eg. MyApp.types
        """
        return self.types_module_path 

#    def getMessagesModulePath(self):
#        '''module path to messages module
#           same as types path
#        '''
#        return self.messages_module_path 

    def setTypesModulePath(self, path):
        """setup module path to where service module before calling fromWsdl.
        module path to types module eg. MyApp.types
        """
        self.types_module_path = path

#    def setMessagesModulePath(self, path):
#        """setup module path to where message module before calling fromWsdl.
#        module path to types module eg. MyApp.types
#        """
#        self.messages_module_path = path

    def gatherNamespaces(self):
        '''This method must execute once..  Grab all schemas
        representing each targetNamespace.
        '''
        if self.usedNamespaces is not None:
            return

        self.logger.debug('gatherNamespaces')
        self.usedNamespaces = {}
            
        # Add all schemas defined in wsdl
        # to used namespace and to the Alias dict
        for schema in self._wsdl.types.values():
            tns = schema.getTargetNamespace()
            self.logger.debug('Register schema(%s) -- TNS(%s)'\
                %(_get_idstr(schema), tns),)
            if self.usedNamespaces.has_key(tns) is False:
                self.usedNamespaces[tns] = []
            self.usedNamespaces[tns].append(schema)
            NAD.add(tns)
            
        # Add all xsd:import schema instances
        # to used namespace and to the Alias dict
        for k,v in SchemaReader.namespaceToSchema.items():
            self.logger.debug('Register schema(%s) -- TNS(%s)'\
                %(_get_idstr(v), k),)
            if self.usedNamespaces.has_key(k) is False:
                self.usedNamespaces[k] = []
            self.usedNamespaces[k].append(v)
            NAD.add(k)
            
    def writeClient(self, fd, sdClass=None, **kw):
        """write out client module to file descriptor.
        Parameters and Keywords arguments:
            fd -- file descriptor
            sdClass -- service description class name
            imports -- list of imports
            readerclass -- class name of ParsedSoap reader
            writerclass -- class name of SoapWriter writer
        """
        sdClass = sdClass or ServiceDescription
        assert issubclass(sdClass, ServiceDescription), \
            'parameter sdClass must subclass ServiceDescription'

#        header = '%s \n# %s.py \n# generated by %s\n%s\n'\
#                  %('#'*50, self.getClientModuleName(), self.__module__, '#'*50)
        print >>fd, '#'*50
        print >>fd, '# file: %s.py' %self.getClientModuleName()
        print >>fd, '# '
        print >>fd, '# client stubs generated by "%s"' %self.__class__
        print >>fd, '#     %s' %' '.join(sys.argv)
        print >>fd, '# '
        print >>fd, '#'*50

        self.services = []
        for service in self._wsdl.services:
            sd = sdClass(self._addressing, do_extended=self.do_extended, 
                         wsdl=self._wsdl)
            if len(self._wsdl.types) > 0:
                sd.setTypesModuleName(self.getTypesModuleName(), 
                                      self.getTypesModulePath())
#                sd.setMessagesModuleName(self.getMessagesModuleName(), 
#                                         self.getMessagesModulePath())

            self.gatherNamespaces()
            sd.fromWsdl(service, **kw)
            sd.write(fd)
            self.services.append(sd)

    def writeTypes(self, fd):
        """write out types module to file descriptor.
        """
        print >>fd, '#'*50
        print >>fd, '# file: %s.py' %self.getTypesModuleName()
        print >>fd, '#'
        print >>fd, '# schema types generated by "%s"' %self.__class__
        print >>fd, '#    %s' %' '.join(sys.argv)
        print >>fd, '#'
        print >>fd, '#'*50
                  
        print >>fd, TypesHeaderContainer()
        self.gatherNamespaces()
        for l in self.usedNamespaces.values():
            sd = SchemaDescription(do_extended=self.do_extended, 
                                   extPyClasses=self.extPyClasses)
            for schema in l:
                sd.fromSchema(schema)
            sd.write(fd)

            
class ServiceDescription:
    """client interface - locator, port, etc classes"""
    separate_messages = False
    logger = _GetLogger("ServiceDescription")

    def __init__(self, addressing=False, do_extended=False, wsdl=None):
        self.typesModuleName = None
        self.messagesModuleName = None
        self.wsAddressing = addressing
        self.imports   = ServiceHeaderContainer()
        self.messagesImports   = ServiceHeaderContainer()
        self.locator   = ServiceLocatorContainer()
        self.bindings   = []
        self.messages  = []
        self.do_extended=do_extended
        self._wsdl = wsdl # None unless do_extended == True

    def setTypesModuleName(self, name, modulePath=None):
        """The types module to be imported.
        Parameters
        name -- name of types module
        modulePath -- optional path where module is located.
        """
        self.typesModuleName = '%s' %name
        if modulePath is not None:
            self.typesModuleName = '%s.%s' %(modulePath,name)

#    def setMessagesModuleName(self, name, modulePath=None):
#        '''The types module to be imported.
#        Parameters
#        name -- name of types module
#        modulePath -- optional path where module is located.
#        '''
#        self.messagesModuleName = '%s' %name
#        if modulePath is not None:
#            self.messagesModuleName = '%s.%s' %(modulePath,name)

    def fromWsdl(self, service, **kw):
        self.imports.setTypesModuleName(self.typesModuleName)
#        if self.separate_messages:
#            self.messagesImports.setMessagesModuleName(self.messagesModuleName)
        self.imports.appendImport(kw.get('imports', []))

        self.locator.setUp(service)

        try:
            bindings =  map(lambda p: p.binding, service.ports)
        except:
            warnings.warn('not all ports have binding declared,')
            bindings = ()

        for port in service.ports:
            if port.binding not in bindings:
                continue
            while port.binding in bindings: 
                bindings.remove(port.binding)

            desc = BindingDescription(useWSA=self.wsAddressing, 
                                      do_extended=self.do_extended, 
                                      wsdl=self._wsdl)
            try:
                desc.setUp(port.getBinding())
            except Wsdl2PythonError, ex:
                self.logger.warning('Skipping port(%s)' %port.name)
                if len(ex.args): 
                    self.logger.warning(ex.args[0])
                continue
   
            desc.setReaderClass(kw.get('readerclass'))
            desc.setWriterClass(kw.get('writerclass'))
            for soc in desc.operations:
                if soc.hasInput() is True:
                    mw = MessageWriter(do_extended=self.do_extended)
                    mw.setUp(soc, port, input=True)
                    self.messages.append(mw)

                    if soc.hasOutput() is True:
                        mw = MessageWriter(do_extended=self.do_extended)
                        mw.setUp(soc, port, input=False)
                        self.messages.append(mw)

            self.bindings.append(desc)

 
    def write(self, fd, msg_fd=None):
        """write out module to file descriptor.
        fd -- file descriptor to write out service description.
        msg_fd -- optional file descriptor for messages module.
        """
#        if msg_fd != None:
#            print >>fd, self.messagesImports
#            print >>msg_fd, self.imports
#        else:
        print >>fd, self.imports
            
        print >>fd, self.locator
        for m in self.bindings:
            print >>fd, m

#        if msg_fd != None:
#            for m in self.messages:
#                print >>msg_fd, m
#        else:
        for m in self.messages:
            print >>fd, m


class MessageWriter:
    logger = _GetLogger("MessageWriter")

    def __init__(self, do_extended=False):
        """Representation of a WSDL Message and associated WSDL Binding.
        operation --
        boperation --
        input --
        rpc --
        literal --
        simple --
        """
        self.content = None
        self.do_extended = do_extended
       
    def __str__(self):
        if not self.content:
            raise Wsdl2PythonError, 'Must call setUp.'
        return self.content.getvalue()
        
    def setUp(self, soc, port, input=False):
        assert isinstance(soc, ServiceOperationContainer),\
            'expecting a ServiceOperationContainer instance'
        assert isinstance(port, WSDLTools.Port),\
            'expecting a WSDL.Port instance'

        rpc,literal = soc.isRPC(), soc.isLiteral(input)
        kw,klass = {}, None
        
        if rpc and literal:
            klass = ServiceRPCLiteralMessageContainer
        elif not rpc and literal:
            kw['do_extended'] = self.do_extended
            klass = ServiceDocumentLiteralMessageContainer
        elif rpc and not literal:
            klass = ServiceRPCEncodedMessageContainer
        else:
            raise WsdlGeneratorError, 'doc/enc not supported.'
                                
        self.content = klass(**kw)
        self.content.setUp(port, soc, input)


class SchemaDescription:
    """generates classes for defs and decs in the schema instance.
    """
    logger = _GetLogger("SchemaDescription")

    def __init__(self, do_extended=False, extPyClasses=None):
        self.classHead = NamespaceClassHeaderContainer()
        self.classFoot = NamespaceClassFooterContainer()
        self.items = []
        self.__types = []
        self.__elements = []
        self.targetNamespace = None
        self.do_extended=do_extended
        self.extPyClasses = extPyClasses

    def fromSchema(self, schema):
        ''' Can be called multiple times, but will not redefine a
        previously defined type definition or element declaration.
        '''  
        ns = schema.getTargetNamespace()
        assert self.targetNamespace is None or self.targetNamespace == ns,\
            'SchemaDescription instance represents %s, not %s'\
            %(self.targetNamespace, ns)

        if self.targetNamespace is None:
            self.targetNamespace = ns
 
        self.classHead.ns = self.classFoot.ns = ns
        for item in [t for t in schema.types if t.getAttributeName() not in self.__types]:
            self.__types.append(item.getAttributeName())
            self.items.append(TypeWriter(do_extended=self.do_extended, extPyClasses=self.extPyClasses))
            self.items[-1].fromSchemaItem(item)

        for item in [e for e in schema.elements if e.getAttributeName() not in self.__elements]:
            self.__elements.append(item.getAttributeName())
            self.items.append(ElementWriter(do_extended=self.do_extended))
            self.items[-1].fromSchemaItem(item)

    def getTypes(self):
        return self.__types

    def getElements(self):
        return self.__elements

    def write(self, fd):
        """write out to file descriptor.
        """
        print >>fd, self.classHead
        for t in self.items:
            print >>fd, t
        print >>fd, self.classFoot

class SchemaItemWriter:
    """contains/generates a single declaration"""
    logger = _GetLogger("SchemaItemWriter")

    def __init__(self, do_extended=False, extPyClasses=None):
        self.content = None
        self.do_extended=do_extended
        self.extPyClasses=extPyClasses
        
    def __str__(self):
        '''this appears to set up whatever is in self.content.localElements,
        local elements simpleType|complexType.
        '''
        assert self.content is not None, 'Must call fromSchemaItem to setup.'
        return str(self.content)

    def fromSchemaItem(self, item):
        raise NotImplementedError, ''


class ElementWriter(SchemaItemWriter):
    """contains/generates a single declaration"""
    logger = _GetLogger("ElementWriter")

    def fromSchemaItem(self, item):
        """set up global elements.
        """
        if item.isElement() is False or item.isLocal() is True:
            raise TypeError, 'expecting global element declaration: %s' %item.getItemTrace()

        local = False
        qName = item.getAttribute('type')
        if not qName:
            etp = item.content
            local = True
        else:
            etp = item.getTypeDefinition('type')
            
        if etp is None:
            if local is True: 
                self.content = ElementLocalComplexTypeContainer(do_extended=self.do_extended)
            else: 
                self.content = ElementSimpleTypeContainer()
        elif etp.isLocal() is False:
            self.content = ElementGlobalDefContainer()
        elif etp.isSimple() is True: 
            self.content = ElementLocalSimpleTypeContainer()
        elif etp.isComplex():
            self.content = ElementLocalComplexTypeContainer(do_extended=self.do_extended)
        else:
            raise Wsdl2PythonError, "Unknown element declaration: %s" %item.getItemTrace()

        self.logger.debug('ElementWriter setUp container "%r", Schema Item "%s"' %(
            self.content, item.getItemTrace()))
        
        self.content.setUp(item)


class TypeWriter(SchemaItemWriter):
    """contains/generates a single definition"""
    logger = _GetLogger("TypeWriter")
        
    def fromSchemaItem(self, item):
        if item.isDefinition() is False or item.isLocal() is True:
            raise TypeError, \
                'expecting global type definition not: %s' %item.getItemTrace()

        self.content = None
        if item.isSimple():
            if item.content.isRestriction():
                self.content = RestrictionContainer()
            elif item.content.isUnion():
                self.content = UnionContainer()
            elif item.content.isList():
                self.content = ListContainer()
            else:
                raise Wsdl2PythonError,\
                    'unknown simple type definition: %s' %item.getItemTrace()
                    
            self.content.setUp(item)
            return
        
        if item.isComplex():
            kw = {}
            if item.content is None or item.content.isModelGroup():
                self.content = \
                    ComplexTypeContainer(\
                        do_extended=self.do_extended, 
                        extPyClasses=self.extPyClasses
                        )
                kw['empty'] = item.content is None
            elif item.content.isSimple():
                    self.content = ComplexTypeSimpleContentContainer()
            elif item.content.isComplex():
                    self.content = \
                        ComplexTypeComplexContentContainer(\
                            do_extended=self.do_extended
                            )
            else:
                raise Wsdl2PythonError,\
                    'unknown complex type definition: %s' %item.getItemTrace()

            self.logger.debug('TypeWriter setUp container "%r", Schema Item "%s"' %(
                self.content, item.getItemTrace()))
            
            try:
                self.content.setUp(item, **kw)
            except Exception, ex:
                args = ['Failure in setUp: %s' %item.getItemTrace()]
                args += ex.args
                ex.args = tuple(args)
                raise
            
            return

        raise TypeError,\
            'expecting SimpleType or ComplexType: %s' %item.getItemTrace()

        

