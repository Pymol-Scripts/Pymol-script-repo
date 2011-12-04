## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/mglutil/util/misc.py,v 1.5.4.3 2009/04/10 22:16:14 vareille Exp $
#
# $Id: misc.py,v 1.5.4.3 2009/04/10 22:16:14 vareille Exp $
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
    if sys.platform == "win32":
        lFont = font[0].upper() + font[1:].lower()
    else:
        lFont = font.lower()
    return lFont


def isInstance(lObject):
    import types
    import inspect
    ltype = type(lObject)
    if ltype == types.InstanceType:
        return True
    elif (ltype == types.ClassType) \
      and inspect.isclass(lObject) is False \
      and isinstance(lObject, ltype) is True:
        return True # this solves the UserList case in python2.6
    else:
        return False


def importMainOrIPythonMain():
    try:
        from IPython import ipapi
        mainDict = ipapi.get().user_ns
    except:
        mainDict = __import__('__main__').__dict__
    return mainDict
