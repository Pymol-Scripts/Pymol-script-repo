############################################################################
# Monte M. Goode, LBNL
# See LBNLCopyright for copyright notice!
###########################################################################

# utility classes used by new generator - mostly 'sugar' classes
# that are actually imported by the generated code.  also includes
# utilities used by wsdl2python itself.

# $Id: utility.py 1226 2006-05-26 18:11:19Z boverhof $

import re
from ZSI import EvaluateException
from ZSI.TCcompound import Struct
from ZSI.generate import WsdlGeneratorError, Wsdl2PythonError
from ZSI.wstools.Utility import SplitQName
from ZSI.wstools.Namespaces import SCHEMA

NCName_to_ModuleName = lambda name: re.sub('\.', '_', TextProtect(name))
NCName_to_ClassName = lambda name: re.sub('\.', '_', TextProtect(name))
TextProtect = lambda s: re.sub('[-./:# ]', '_', s)
TextProtectAttributeName = lambda name: TextProtect('_%s' %name)
Namespace2ModuleName = lambda ns: TextProtect(ns.lstrip('http://')).rstrip('_')


def GetModuleBaseNameFromWSDL(wsdl):
    """By default try to construct a reasonable base name for all
    generated modules.  Otherwise return None.
    """
    base_name = wsdl.name or wsdl.services[0].name
    base_name = SplitQName(base_name)[1]
    if base_name is None: 
        return None
    return NCName_to_ModuleName(base_name)

namespace_name = lambda cls, ns: 'ns%s' % len(cls.alias_list)

class NamespaceAliasDict:
    """a lookup table to store relevant namespaces and their aliases"""
    alias_dict = {}
    alias_list = []

    def add(cls, ns):
        if cls.alias_dict.has_key(ns):
            return
        cls.alias_dict[ns] = (Namespace2ModuleName(ns), '%s' % namespace_name(cls,ns))
        cls.alias_list.append(ns)
    add = classmethod(add)
            
    def getModuleName(cls, ns):
        if cls.alias_dict.has_key(ns):
            return cls.alias_dict[ns][0]
                                 
        msg = 'failed to find import for schema "%s"'%ns +\
        'possibly missing @schemaLocation attribute.'
        if ns in SCHEMA.XSD_LIST:
            msg = 'missing built-in typecode for schema "%s"' %ns
            
        raise WsdlGeneratorError, msg
                                 
    getModuleName = classmethod(getModuleName)
        
    def getAlias(cls, ns):
        if cls.alias_dict.has_key(ns):
            return cls.alias_dict[ns][1]
                                 
        msg = 'failed to find import for schema "%s"'%ns +\
        'possibly missing @schemaLocation attribute.'
        if ns in SCHEMA.XSD_LIST:
            msg = 'missing built-in typecode for schema "%s"' %ns
            
        raise WsdlGeneratorError, msg
    
    getAlias = classmethod(getAlias)

    def getNSList(cls):
        return tuple(cls.alias_list)
    getNSList = classmethod(getNSList)


class StringWriter:
    """generator util"""
    def __init__(self, val=None):
        self.data = []
        if val:
            self.data.append(val)

    def set(self, val):
        if self.data:
            # in some cases the empty list reassignment fails, so....
            self.data = None
            self.data = []

        self.data.append(val)

    def write(self, val):
        self.data.append(val)

    def getvalue(self):
        if self.data:
            return ''.join(self.data)
        else:
            return ''

    def __iadd__(self, val):
        self.data.append(val)
        return self

    def __str__(self):
        return self.getvalue()


# ---- generated code utils

class MessageContainer:
    """generator util - used by address.py"""
    pass

# Extract sub names from message parts so they can be used when mapping
#    a message's contents to a function's arguments.
# Args is a list of Message Parts.  i.e.: op.getInputMessage().parts.values()
def GetPartsSubNames(args, wsdl):
    do_extended = True
    from wsdl2python import WriteServiceModule, SchemaDescription
    wsm = WriteServiceModule(wsdl, do_extended=do_extended)
    wsm.gatherNamespaces()
    toReturn = []
    for arg in args:
        argSubnames = []
        for l in wsm.usedNamespaces.values():
            for schema in l:
                sd = SchemaDescription(do_extended=do_extended)
                sd.fromSchema(schema)
                argNamespace = arg.element[0]
                if (sd.targetNamespace == argNamespace):
                    for i in sd.items:
                        # arg.name is the part name, but we want it's type
                        argElementType = arg.element[1]
                        if str(argElementType) == str(i.content.name):
                            argSubnames = []
			    # I'm not sure when the name attribute was dropped
			    # but at some point, or in some circumstance it's not
			    # there, but instead a ref attribute is there which is
		     	    # tuple of (namespace, name). This hack fixes things, 
			    # but I'm not sure why this happens or has happened.
			    # IRJ - 2005-05-25
                            if i.content.mgContent != None:
                                for c in i.content.mgContent:
                                    nValue = "None"
                                    if c.isWildCard():
                                        nValue="any"
                                    elif c.attributes.has_key("name"):
                                        nValue = c.attributes["name"]
                                    elif c.attributes.has_key("ref"):
                                        nValue = c.attributes["ref"][1]
                                    argSubnames.append(nValue)

        toReturn.append(argSubnames)
    return toReturn
