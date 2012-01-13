#!/usr/bin/env python
import inspect
from cStringIO import StringIO
import ZSI, string, sys, getopt, urlparse, types, warnings
from ZSI.wstools import WSDLTools
from ZSI.ServiceContainer import ServiceSOAPBinding, SimpleWSResource, WSAResource

from ZSI.generate import WsdlGeneratorError, Wsdl2PythonError
from utility import TextProtect, GetModuleBaseNameFromWSDL, \
    NCName_to_ClassName, GetPartsSubNames, TextProtectAttributeName
from containers import BindingDescription
from wsdl2python import MessageWriter, WriteServiceModule,\
    MessageTypecodeContainer, SchemaDescription

# Split last token
rsplit = lambda x,sep,: (x[:x.rfind(sep)], x[x.rfind(sep)+1:],)
if sys.version_info[0:2] == (2, 4, 0, 'final', 0)[0:2]:
    rsplit = lambda x,sep,: x.rsplit(sep, 1)


class SOAPService:
    def __init__(self, service):
        self.classdef = StringIO()
        self.initdef  = StringIO()
        self.location = ''
        self.methods  = []

    def newMethod(self):
        '''name -- operation name
        '''
        self.methods.append(StringIO())
        return self.methods[-1]


class ServiceModuleWriter:
    '''Creates a skeleton for a SOAP service instance.
    '''
    indent = ' '*4
    server_module_suffix = '_server'
    func_aname = TextProtectAttributeName
    func_aname = staticmethod(func_aname)
    separate_messages = False 

    def __init__(self, base=ServiceSOAPBinding, prefix='soap', 
                 service_class=SOAPService):
        '''
        parameters:
            base -- either a class definition, or a str representing a qualified 
                class name (eg. module.name.classname)
            prefix -- method prefix.
        '''
        if inspect.isclass(base):
            self.base_class_name = base.__name__
            self.base_module_name = inspect.getmodule(base).__name__
        else:
            self.base_module_name, self.base_class_name  = base.rsplit('.', 1)

        self.wsdl = None
        self.method_prefix = prefix
        self._service_class = SOAPService

        self.header  = None
        self.imports  = None
        self.messages = []
        self._services = None
        self.types_module_path = None
        self.types_module_name = None
        self.messages_module_name = None

    def reset(self):
        self.header  = StringIO()
        self.imports  = StringIO()
        self.message = []
        self._services = {}

    def getIndent(self, level=1):
        '''return indent.
        '''
        assert 0 < level < 10, 'bad indent level %d' %level
        return self.indent*level

    def getMethodName(self, method):
        '''return method name.
        '''
        return '%s_%s' %(self.method_prefix, TextProtect(method))

    def getClassName(self, name):
        '''return class name.
        '''
        return NCName_to_ClassName(name)

    def setTypesModuleName(self, name):
        self.types_module_name = name

    # Backwards compatibility
    setClientModuleName = setTypesModuleName
    
    def getTypesModuleName(self):
        '''return module name.
        '''
        assert self.wsdl is not None, 'initialize, call fromWSDL'
        if self.types_module_name is not None:
            return self.types_module_name

        wsm = WriteServiceModule(self.wsdl)
        return wsm.getTypesModuleName()

    def getServiceModuleName(self):
        '''return module name.
        '''
        name = GetModuleBaseNameFromWSDL(self.wsdl)
        if not name:
            raise WsdlGeneratorError, 'could not determine a service name'
        
        if self.server_module_suffix is None:
            return name
        return '%s%s' %(name, self.server_module_suffix)

    def getTypesModulePath(self):
        return self.types_module_path
    getClientModulePath = getTypesModulePath
    
    def setTypesModulePath(self, path):
        '''setup module path to where client module before calling fromWSDL.
        '''
        self.types_module_path = path
    setClientModulePath = setTypesModulePath
    
    def setUpClassDef(self, service):
        '''set class definition and class variables.
        service -- ServiceDescription instance
        '''
        assert isinstance(service, WSDLTools.Service) is True,\
            'expecting WSDLTools.Service instance.'

        s = self._services[service.name].classdef

        print >>s, 'class %s(%s):' %(self.getClassName(service.name), self.base_class_name)

        print >>s, '%ssoapAction = {}' % self.getIndent(level=1)
        print >>s, '%sroot = {}' % self.getIndent(level=1)
        
    def setUpImports(self):
        '''set import statements
        '''
        i = self.imports
        print >>i, 'from ZSI.schema import GED, GTD'
        print >>i, 'from ZSI.TCcompound import ComplexType, Struct'

        module = self.getTypesModuleName()
        package = self.getTypesModulePath()
        if package:
            module = '%s.%s' %(package, module)
            
        print >>i, 'from %s import *' %(module)
            
        print >>i, 'from %s import %s' %(self.base_module_name, self.base_class_name)

    def setUpInitDef(self, service):
        '''set __init__ function
        '''
        assert isinstance(service, WSDLTools.Service), \
            'expecting WSDLTools.Service instance.'
            
        sd = self._services[service.name]
        d = sd.initdef
 
        if sd.location is not None:
            scheme,netloc,path,params,query,fragment = urlparse.urlparse(sd.location)
            print >>d, '%sdef __init__(self, post=\'%s\', **kw):' %(self.getIndent(level=1), path)
        else:
            print >>d, '%sdef __init__(self, post, **kw):' %self.getIndent(level=1)

        # Require POST initialization value for test implementation
        if self.base_module_name == inspect.getmodule(ServiceSOAPBinding).__name__:
            print >>d, '%s%s.__init__(self, post)' %(self.getIndent(level=2), self.base_class_name)
            return 

        # No POST initialization value, obtained from HTTP Request in twisted or wsgi
        print >>d, '%s%s.__init__(self)' %(self.getIndent(level=2), self.base_class_name)

    def mangle(self, name):
        return TextProtect(name)

    def getAttributeName(self, name):
        return self.func_aname(name)

    def setUpMethods(self, port):
        '''set up all methods representing the port operations.
        Parameters:
            port -- Port that defines the operations.
        '''
        assert isinstance(port, WSDLTools.Port), \
            'expecting WSDLTools.Port not: ' %type(port)

        sd = self._services.get(port.getService().name)
        assert sd is not None, 'failed to initialize.'

        binding = port.getBinding()
        portType = port.getPortType()
        action_in = ''
        for bop in binding.operations:
            try:
                op = portType.operations[bop.name]
            except KeyError, ex:
                raise WsdlGeneratorError,\
                    'Port(%s) PortType(%s) missing operation(%s) defined in Binding(%s)' \
                    %(port.name,portType.name,bop.name,binding.name)

            for ext in bop.extensions:
                 if isinstance(ext, WSDLTools.SoapOperationBinding):
                     action_in = ext.soapAction
                     break
            else:
                warnings.warn('Port(%s) operation(%s) defined in Binding(%s) missing soapAction' \
                    %(port.name,op.name,binding.name)
                )

            msgin = op.getInputMessage()
            msgin_name = TextProtect(msgin.name)
            method_name = self.getMethodName(op.name)

            m = sd.newMethod()
            print >>m, '%sdef %s(self, ps, **kw):' %(self.getIndent(level=1), method_name)
            if msgin is not None:
                print >>m, '%srequest = ps.Parse(%s.typecode)' %(self.getIndent(level=2), msgin_name)
            else:
                print >>m, '%s# NO input' %self.getIndent(level=2)

            msgout = op.getOutputMessage()
            if msgout is not None:
                msgout_name = TextProtect(msgout.name)
                print >>m, '%sreturn request,%s()' %(self.getIndent(level=2), msgout_name)
            else:
                print >>m, '%s# NO output' % self.getIndent(level=2)
                print >>m, '%sreturn request,None' % self.getIndent(level=2)

            print >>m, ''
            print >>m, '%ssoapAction[\'%s\'] = \'%s\'' %(self.getIndent(level=1), action_in, method_name)
            print >>m, '%sroot[(%s.typecode.nspname,%s.typecode.pname)] = \'%s\'' \
                     %(self.getIndent(level=1), msgin_name, msgin_name, method_name)

        return

    def setUpHeader(self):
        print >>self.header, '#'*50
        print >>self.header, '# file: %s.py' %self.getServiceModuleName()
        print >>self.header, '#'
        print >>self.header, '# skeleton generated by "%s"' %self.__class__
        print >>self.header, '#      %s' %' '.join(sys.argv)
        print >>self.header, '#'
        print >>self.header, '#'*50

    def write(self, fd=sys.stdout):
        '''write out to file descriptor, 
        should not need to override.
        '''
        print >>fd, self.header.getvalue()
        print >>fd, self.imports.getvalue()
        
        print >>fd, '# Messages ',
        for m in self.messages:
            print >>fd, m
        
        print >>fd, ''
        print >>fd, ''
        print >>fd, '# Service Skeletons'
        for k,v in self._services.items():
            print >>fd, v.classdef.getvalue()
            print >>fd, v.initdef.getvalue()
            for s in v.methods:
                print >>fd, s.getvalue()

    def fromWSDL(self, wsdl):
        '''setup the service description from WSDL,
        should not need to override.
        '''
        assert isinstance(wsdl, WSDLTools.WSDL), 'expecting WSDL instance'

        if len(wsdl.services) == 0:
            raise WsdlGeneratorError, 'No service defined'
            
        self.reset() 
        self.wsdl = wsdl
        self.setUpHeader()
        self.setUpImports()
                
        for service in wsdl.services:
            sd = self._service_class(service.name)
            self._services[service.name] = sd

            for port in service.ports:
                desc = BindingDescription(wsdl=wsdl)
                try:
                    desc.setUp(port.getBinding())
                except Wsdl2PythonError, ex:
                    continue
                
                for soc in desc.operations:
                    if not soc.hasInput(): continue
                    
                    self.messages.append(MessageWriter())
                    self.messages[-1].setUp(soc, port, input=True)
                    if soc.hasOutput():
                        self.messages.append(MessageWriter())
                        self.messages[-1].setUp(soc, port, input=False)
                
                for e in port.extensions:
                    if isinstance(e, WSDLTools.SoapAddressBinding):
                        sd.location = e.location

                self.setUpMethods(port)

            self.setUpClassDef(service)
            self.setUpInitDef(service)


class WSAServiceModuleWriter(ServiceModuleWriter):
    '''Creates a skeleton for a WS-Address service instance.
    '''
    def __init__(self, base=WSAResource, prefix='wsa', service_class=SOAPService, 
                 strict=True):
        '''
        Parameters:
            strict -- check that soapAction and input ws-action do not collide.
        '''
        ServiceModuleWriter.__init__(self, base, prefix, service_class)
        self.strict = strict

    def createMethodBody(msgInName, msgOutName, **kw):
        '''return a tuple of strings containing the body of a method.
        msgInName -- None or a str
        msgOutName --  None or a str
        '''
        body = []
        if msgInName is not None:
            body.append('request = ps.Parse(%s.typecode)' %msgInName)
            
        if msgOutName is not None:
            body.append('return request,%s()' %msgOutName)
        else: 
            body.append('return request,None')
            
        return tuple(body)
    createMethodBody = staticmethod(createMethodBody)

    def setUpClassDef(self, service):
        '''use soapAction dict for WS-Action input, setup wsAction
        dict for grabbing WS-Action output values.
        '''
        assert isinstance(service, WSDLTools.Service), \
            'expecting WSDLTools.Service instance'

        s = self._services[service.name].classdef
        print >>s, 'class %s(%s):' %(self.getClassName(service.name), self.base_class_name)
        print >>s, '%ssoapAction = {}' % self.getIndent(level=1)
        print >>s, '%swsAction = {}' % self.getIndent(level=1)
        print >>s, '%sroot = {}' % self.getIndent(level=1)

    def setUpMethods(self, port):
        '''set up all methods representing the port operations.
        Parameters:
            port -- Port that defines the operations.
        '''
        assert isinstance(port, WSDLTools.Port), \
            'expecting WSDLTools.Port not: ' %type(port)

        binding = port.getBinding()
        portType = port.getPortType()
        service = port.getService()
        s = self._services[service.name]
        for bop in binding.operations:
            try:
                op = portType.operations[bop.name]
            except KeyError, ex:
                raise WsdlGeneratorError,\
                    'Port(%s) PortType(%s) missing operation(%s) defined in Binding(%s)' \
                    %(port.name, portType.name, op.name, binding.name)

            soap_action = wsaction_in = wsaction_out = None
            if op.input is not None:
                wsaction_in = op.getInputAction()
            if op.output is not None:
                wsaction_out = op.getOutputAction()

            for ext in bop.extensions:
                if isinstance(ext, WSDLTools.SoapOperationBinding) is False: continue
                soap_action = ext.soapAction
                if not soap_action: break
                if wsaction_in is None: break
                if wsaction_in == soap_action: break
                if self.strict is False:
                    warnings.warn(\
                        'Port(%s) operation(%s) in Binding(%s) soapAction(%s) != WS-Action(%s)' \
                         %(port.name, op.name, binding.name, soap_action, wsaction_in),
                    )
                    break
                raise WsdlGeneratorError,\
                    'Port(%s) operation(%s) in Binding(%s) soapAction(%s) MUST match WS-Action(%s)' \
                     %(port.name, op.name, binding.name, soap_action, wsaction_in)

            method_name = self.getMethodName(op.name)

            m = s.newMethod()
            print >>m, '%sdef %s(self, ps, address):' %(self.getIndent(level=1), method_name)
            
            msgin_name = msgout_name = None
            msgin,msgout = op.getInputMessage(),op.getOutputMessage()
            if msgin is not None: 
                msgin_name = TextProtect(msgin.name)
            if msgout is not None: 
                msgout_name = TextProtect(msgout.name)
        
            indent = self.getIndent(level=2)
            for l in self.createMethodBody(msgin_name, msgout_name):
                print >>m, indent + l

            print >>m, ''
            print >>m, '%ssoapAction[\'%s\'] = \'%s\'' %(self.getIndent(level=1), wsaction_in, method_name)
            print >>m, '%swsAction[\'%s\'] = \'%s\'' %(self.getIndent(level=1), method_name, wsaction_out)
            print >>m, '%sroot[(%s.typecode.nspname,%s.typecode.pname)] = \'%s\'' \
                     %(self.getIndent(level=1), msgin_name, msgin_name, method_name)
 
