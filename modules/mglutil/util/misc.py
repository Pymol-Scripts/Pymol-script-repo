## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/mglutil/util/misc.py,v 1.18 2009/11/06 21:08:50 annao Exp $
#
# $Id: misc.py,v 1.18 2009/11/06 21:08:50 annao Exp $
#

import types
import sys
import numpy.oldnumeric as Numeric

def issequence(a):
    return type(a) is types.TupleType or \
           type(a) is types.ListType or \
           isinstance(a, Numeric.ArrayType)

def isnumericstring(a):
    try:
        float(a)
        return 1
    except:
        return 0

def uniq(objectSequence):
    """Remove the duplicates from a list while keeping the original
    list order """
    l = []
    d = {}
    for o in objectSequence:
        if not d.has_key(o):
            d[o] = None
            l.append(o)
    return l


def deepCopySeq(sequence):
    """ returns the deep copy of the given sequence """
    import numpy.oldnumeric as Numeric
    from types import TupleType, ListType
    assert type(sequence) in (TupleType, ListType, type(Numeric.array([1,2,3])))
    if hasattr(sequence, 'copy'):
        dcSeq = sequence.copy()
    else:
        dcSeq = sequence[:]

    return dcSeq


def ensureFontCase(font):
    return font
#    from Tkinter import TkVersion
#    lFont = font[0].upper() + font[1:].lower()
#    if TkVersion == '8.4' and sys.platform != "win32":
#        lFont = font.lower()
#    return lFont


def isInstance(lObject):

    import types
    if sys.version.startswith('2.5'): #detect python25
        if type(lObject) == types.InstanceType:
            return True
        else:
            return False
    else:
            import inspect
            ltype = type(lObject)
            if ltype == types.InstanceType:
                return True
            elif inspect.isclass(lObject) is False \
              and isinstance(lObject, ltype) is True:
                from abc import ABCMeta
                if ltype == types.ClassType is True:
                    return True
                elif type(ltype) == ABCMeta:
                    return True
                else:
                    return False
            else:
                return False


def importMainOrIPythonMain():
    try:
        from IPython import ipapi
        mainDict = ipapi.get().user_ns
    except:
        mainDict = __import__('__main__').__dict__
    return mainDict


def suppressMultipleQuotes(aString):
    lStringToSimplify = aString
    lSimplifiedString = lStringToSimplify
    while type(lStringToSimplify) == types.StringType:
        lSimplifiedString = lStringToSimplify
        try:
           lStringToSimplify = eval(lSimplifiedString)
        except:
           break
    return lSimplifiedString


class IntVar:
    def __init__(self, val=0):
        self.set(val)

    def get(self):
        return self.val

    def set(self,val):
        self.val = int(val)


class StringVar:
    def __init__(self, val=""):
        self.set(val)

    def get(self):
        return self.val

    def set(self,val):
        self.val = str(val)


class BooleanVar:
    def __init__(self, val=False):
        self.set(val)

    def get(self):
        return self.val

    def set(self,val):
        self.val = (val==True)
