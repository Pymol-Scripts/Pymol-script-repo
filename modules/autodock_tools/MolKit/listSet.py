#############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/MolKit/listSet.py,v 1.32 2010/10/08 22:05:37 sanner Exp $
#
# $Id: listSet.py,v 1.32 2010/10/08 22:05:37 sanner Exp $
#

"""
This module implements a Set class that uses a List to store the objects
in the set.
"""

import types
from UserList import UserList
from mglutil.util import misc

verbose = False

class ListSet(UserList):  # would (list) work???
    """Class to represent Sets of objects stored in a list. There is an
    implicit order amongst the objects and there can be duplicate objects.

    __getattr__, __setattr__ and __delattr__ have been modified to operate on
    the list of objects rather than the TreeNodeSet itself, i.e. if atm is an
    instance of a ListSet a.xxx will not return the member xxx of the object
    atm but rather a list of the members xxx from each object in the set atm.
    xxx can be a member of a function that requires no argument.

    Example:

      if atm is an instance of a ListSet:
      atm.name            return the name attribute of each Atom in atm
      atm.newprop = 7.2   creates a newprop attribute for each Atom in atm
                          with an initial value of 7.2
      atm.newIndex = range(len(atm)) create a newIndex attribute for each Atom
                          in atm with values 0 for the first atom, 1 for the
                          second, etc...
      del atm.newIndex
      
    This class also implement  boolean operations on ListSets. These operation
    overload some operators.
    
    A uniq() method returns a list with the double removed.
    A makeUnique() method removes duplicates from list (in place).
    """
    
    def __init__(self, data=None, elementType=None, stringRepr=None, 
                 comments="", keywords=[]):
        if data is not None and len(data)!=0:
            assert( hasattr(data[0], '__class__') )
            if elementType is not None:
                assert( isinstance(data[0], elementType) )
                
        UserList.__init__(self, data)
        # save call to __setattr__
        self.__dict__['elementType'] = elementType
        #self.elementType = elementType # class of elements in that set
        self.stringRepr = stringRepr # this will hold a concise string repr of
                                     # the set of objects in that list
        self.__dict__['comments'] = comments
        self.__dict__['keywords'] = keywords
        self.selector = None

        
    def buildRepr(self):
        # walk up the tree and add repr strings when no all children
        # are selected
        names = ""
        selDict = {}.fromkeys(self.data) # dict of objects in Set
        parents = self.parent
        if parents[0] != None:
            parents = parents.uniq() # uniq parents
        else:
            for obj in self.data:
                names += obj.name+';'
            return names[:-1]
        level = 0
        while True:
            pselDict = {}
            for p in parents:
                selChildren = [x for x in p.children if selDict.has_key(x)]
                if len(selChildren) == len(p.children):
                    pselDict[p] = True # all children are in so it might
                                       # be added at next level
                else: # not all children are in so we put the full name
                      # of each child and add ':' to get to the level of self
                    name = p.full_name()+':'
                    for c in selChildren:
                        name += c.name + ','
                    names += name[:-1] + ':'*level + ';'
            if len(pselDict)==0:
                break
            level = level+1
            nparents = {}.fromkeys([x.parent for x in pselDict.keys()]).keys()
            #if len(pselDict)<10:
            #    print 'AAAAAAAAAA', pselDict, nparents, parents
            if nparents[0]==None:
                for p in parents:
                    selChildren = [x for x in p.children if selDict.has_key(x)]
                    if len(selChildren) == len(p.children):
                        names += p.full_name()+':'*level+';'
                break
            selDict = pselDict
            parents = nparents

        return names[:-1]
##         # try to identify residue ranges
##         newnames = ''
##         # at this point parents holds chains

##         # split each selection string using ';'
##         for selector in names[:-1].split(';'):
##             # find out how many levels
##             levels = selector.split(':')
##             if len(levels)==3 or len(levels)==4 and len(levels[3])==0:
##                 # look at residue level
##                 chainId = levels[1]
##                 chain = [c for c in parents if c.id==chainId][0]
##                 allResidues = chain.residues
##                 residues = allResidues.get(levels[2])
##                 indices = [allResidues.index(x) for x in residues]
##                 resNames = levels[2].split(',')
##                 nrn = resNames[0]
##                 ci = indices[0]
##                 inRange = False
##                 for ind, rname in zip(indices[1:], resNames[1:]):
##                     if ind!=ci+1: # not consecutive
##                         if inRange:
##                             nrn += prn+','
##                         else:
##                             nrn += ','+rname
##                         inRange = False
##                     else:
##                         if not inRange:
##                             nrn += '-'
##                         inRange = True
##                     ci = ind
##                     prn = rname
##                 nrn += rname
##                 newnames += levels[0]+':'+levels[1]+':'+nrn
##                 if len(levels)==4:
##                     newnames += ':'
##                 newnames += ';'
##             else: # no residues
##                 newnames += selector+';'
##         return newnames[:-1]


    def setStringRepr(self, string):
        """set the string representation of this set"""
        assert type(string) in types.StringTypes
        self.stringRepr = string


    def getStringRepr(self):
        """return the string representation of this set"""
        return self.stringRepr


    def copy(self):
        """return a copy of the set"""
        copy = self.__class__(self.data, stringRepr = self.stringRepr)
        return copy

 
    def __str__(self):
        """add here because __str__ is missing in UserList which creates a pb
        in jpython"""
        return str(self.data)

        
    def __delattr__(self, member):
        if member[:2]=='__' or  member in ['data', 'elementType']:
            return
        func = 'if hasattr(o,"%s"): del o.%s' % (member,member)
        for o in self.data:
            exec (func)
        

    def __iter__(self, *cfg, **kw):
        return iter(self.data)
    
        
    def __getattr__(self, member):
        """Extract the specified member from each objects in the set and
        returns them as a list
        """

        if member[:2]=='__':
            if self.__dict__.has_key(member):
                return self.__dict__[member]
            else:
                raise AttributeError('member %s not found'%member)
        elif member in ['data', 'elementType', 'stringRepr','comments', 'keywords', 'selector']:
            return self.__dict__[member]
        else:

            result = []
#            if len(self.data) and callable( eval('self.data[0].%s' % member )):
            if len(self.data) and callable( getattr(self.data[0], member) ):
                m = self.data[0].__class__.__dict__[member]
                for o in self.data:
                    result.append( m(o) )
#                result = map( eval ('lambda x: x.%s()' % member), self.data )
            else:
                for o in self.data:
                    result.append( o.__dict__[member] )
#                result = map( eval ('lambda x: x.%s' % member), self.data )

            return result

    def getAll(self, member):
        return self.__getattr__(member)

    def setAll(self, member, value):
        return self.__setattr__(member, value)


    def setSetAttribute(self, name, value):
        """
        set an attribute for the Set, rather than for the objects in the set
        """
        self.__dict__[name] = value

        
    def __setattr__(self, member, value):
        """Set or create member in each object in this set.
        If value is a sequence it has to be of the same length as the set.
        else the new member in each object in the set is set to 'value'
        """

        if member[:2]=='__':
            self.__dict__[member] = value
        elif member in ['data', 'elementType', 'stringRepr', 'comments', 'keywords', 'selector']:
            self.__dict__[member] = value
        else:
            l = len(self.data)
            if not misc.issequence(value): # simple number of string
                #map( setThem, self.data, (value,)*l, (member,)*l )
                #for o in self.data: exec ( 'o.%s = value'%member )
                for o in self.data: o.__dict__[member] = value

            else: # value is a sequence

                if len(value) == 0: # empty sequence .. assign sequence to attribute
                     #map( setThem, self.data, value*l, (member,)*l )
                     #for o in self.data: exec ( 'o.%s = value'%member )
                     self.__dict__[member] = value

                elif len(value) == 1: # sequence of 1 .. treat like single number
                     #map( setThem, self.data, value*l, (member,)*l )
                     #for o in self.data: exec ( 'o.%s = value'%member )
                     for o in self.data: o.__dict__[member] = value[0]

                else: # sequence of many values
                    assert len(self.data) == len(value)
                    for o,v in map(None, self.data, value):
                        setattr(o, member, v)
#                    for i in xrange(len(self.data)):
#                        self.data[i].__dict__[member] = value[i]


    def append(self, item):
        if self.elementType is not None:
            assert isinstance(item, self.elementType)
        if len(self.data)>0 and self.stringRepr and hasattr(item, 'full_name'):
            self.stringRepr = self.stringRepr+'/+/'+item.full_name()
        elif hasattr(item, 'full_name'):
            self.stringRepr = item.full_name()
        self.data.append(item)
    
    def insert(self, i, item):
        if self.elementType is not None:
            assert isinstance(item, self.elementType)
        if len(self.data)>0 and self.stringRepr and hasattr(item, 'full_name'):
            self.stringRepr = self.stringRepr+'/+/'+item.full_name()
        elif hasattr(item, 'full_name'):
            self.stringRepr = item.full_name()
        self.data.insert(i, item)
    
    def pop(self, i=-1):
        if self.elementType is not None:
            assert isinstance(item, self.elementType)
        item = self.data.pop(i)
        if len(self.data)>0 and self.stringRepr and hasattr(item, 'full_name'):
            self.stringRepr = self.stringRepr+'/-/'+item.full_name()
        else:
            self.stringRepr = None
        return item
    
    def remove(self, item):
        if self.elementType is not None:
            assert isinstance(item, self.elementType)
        self.data.remove(item)
        if len(self.data)>0 and self.stringRepr and hasattr(item, 'full_name'):
            self.stringRepr = self.stringRepr+'/-/'+item.full_name()
        else:
            self.stringRepr = None

    def __getslice__(self, i, j):
        to_end = False
        if j>len(self.data)-1:
            j=len(self.data)
            if i==0: to_end = True
        #else:
        #    j += 1  # %d-%d in selection language includes start and end
                    # while [i:j] excludes j
        if self.stringRepr:
            # build a string repr which is the and of current set and
            # range [i,j] on surrent level
            stringRepr = self.stringRepr+'/&/'
            ind = self.stringRepr.rfind(':')
            if to_end:
                stringRepr = self.stringRepr
            elif ind > 0: # we found ':'
                if self.stringRepr[-1]==':': #nothing beyond last ':'
                    # everything at this level is in the set
                    # we only have to add the range
                    if i<j-1:
                        stringRepr = self.stringRepr[:ind+1]+'%d-%d'%(i,j-1)
                    else:  #could be only the last item
                        stringRepr = self.stringRepr[:ind+1]+'%d'%(i)
                else:
                    # only a subset of the current level is in self
                    # we have to do a sub-selection
                    if i<j-1:
                        stringRepr = self.stringRepr+'\\s\\'+'%d-%d'%(i,j-1)
                    else:  #could be only the last item
                        stringRepr = self.stringRepr+'\\s\\'+'%d'%(i)
            else: # we are the root level
                stringRepr = self.data[i].name
                for m in self.data[i+1:j]:
                    stringRepr += ','+m.name
        else:
            if verbose: 
                print 'WARNING long stringRepr due to getslice'
            stringRepr = ''
            for obj in self.data[i:j]:
                stringRepr += obj.full_name()+';'
            stringRepr = stringRepr[:-1] # remove last semi colon
        return self.__class__(self.data[i:j], stringRepr=stringRepr)
        

    def __delslice__(self, i, j):
        if verbose:
            print 'WARNING long stringRepr due to delslice'
        del self.data[i:j]
        stringRepr = ''
        for obj in self.data:
            stringRepr += obj.full_name()+';'
        self.stringRepr = stringRepr[:-1] # remove last semi colon


    def __mul__(self, n):
        # returns a new set
        if len(self.data)==0:
            return self.__class__([])
        origStringRepr = self.stringRepr
        for i in range(i-1):
            stringRepr += '/+/'+origStringRepr
        return self.__class__(self.data*n, stringRepr=stringRepr)


    def __imul__(self, n):
        # multiplies the set in place
        if len(self.data)==0:
            return self
        self.data *= n
        origStringRepr = self.stringRepr
        for i in range(i-1):
            stringRepr += '/+/'+origStringRepr
        self.stringRepr = stringRepr
        return self
        

    def extend(self, right):
        assert isinstance(right, self.__class__)
        if len(right.data)==0: return
        self.data.extend(right.data)
        if self.stringRepr and right.stringRepr:
            self.stringRepr = self.stringRepr+'/+/'+right.stringRepr
        elif verbose:
            import traceback
            traceback.print_stack()
            print 'extending sets with no stringRepr:', repr(self), repr(right)
        
        
    def __iadd__(self, right):
        """See add: overloads += operator"""
        self.extend(right)
        return self
    

    def __add__(self, right):
        """See add: overloads + operator"""
        assert isinstance(right, self.__class__)
        if len(right.data)==0: return self.copy()
        if len(self.data)==0: return right.copy()
        stringRepr = None
        if self.stringRepr and right.stringRepr:
            stringRepr = self.stringRepr+'/+/'+right.stringRepr
        elif verbose:
            import traceback
            traceback.print_stack()
            print 'adding sets with no stringRepr:', repr(self), repr(right)
            stringRepr = None
        return self.__class__( self.data + right.data, stringRepr=stringRepr)


    def union(self, right):
        """Union: returns a Set holding objects appearing in either list"""

        assert isinstance(right, self.__class__)
        stringRepr = None
        if len(right.data)==0: return self.copy()
        if len(self.data)==0: return right.copy()
        if self.stringRepr and right.stringRepr:
            if self.stringRepr == right.stringRepr:
                stringRepr = self.stringRepr
            else:
                stringRepr = self.stringRepr+'/|/'+right.stringRepr
        elif verbose:
            import traceback
            traceback.print_stack()
            print 'union of sets with no stringRepr:', repr(self), repr(right)
            stringRepr = None
        return self.__class__( misc.uniq(self.data + right.data),
                              stringRepr=stringRepr )


    def __or__(self, right):
        """See union: overloads | operator"""
        return self.union(right)


    def xor(self, right):
        """XOR operation: Returns a set made of the elements appearing in first
        or second set but not in both"""

        assert isinstance(right, self.__class__)
        if len(right.data)==0: return self.copy()
        if len(self.data)==0: return right.copy()
        stringRepr = None
        l1 = ListSet.__sub__(self, right)
        l2 = ListSet.__sub__(right, self)
        if self.stringRepr and right.stringRepr:
            stringRepr = self.stringRepr+'/^/'+right.stringRepr
        elif verbose:
            import traceback
            traceback.print_stack()
            print 'xoring sets with no stringRepr:', repr(self), repr(right)
            stringRepr = None
        return self.__class__( l1.data + l2.data, stringRepr=stringRepr )


    def __xor__(self, right):
        """See union: overloads ^ operator"""
        return self.xor(right)
        

    def inter(self, right):
        """Intersection: returns a Set holding objects appearing in both sets
        """

        assert isinstance(right, self.__class__)
        if len(right.data)==0: return self.copy()
        if len(self.data)==0: return right.copy()
        l1 = self
        l2 = right
        if len(l1.data) > len(right.data):
            l1 = right
            l2 = self
        # l1 is always shorter list
        for o in l2.data: o._setFlag = 0
        for o in l1.data: o._setFlag = 1
        newlist = filter( lambda x: x._setFlag==1, l2.data )
        for o in l2.data:
            if hasattr(o, '_setFlag'):
                del o._setFlag
        for o in l1.data:
            if hasattr(o, '_setFlag'):
                del o._setFlag
        stringRepr = None
        if self.stringRepr and right.stringRepr:
            stringRepr = self.stringRepr+'/&/'+right.stringRepr
        elif verbose:
            import traceback
            traceback.print_stack()
            print 'intersecting sets with no stringRepr:', repr(self), repr(right)
            stringRepr = None
        return self.__class__(misc.uniq(newlist), stringRepr=stringRepr)


    def __and__(self, right):
        """See inter: overloads & operator"""
        return self.inter(right)
        

    def subtract(self, right):
        """Returns a set made of the elements of the first set not appearing
        in the second set"""

        stringRepr = None
        assert isinstance(right, self.__class__)
        if len(right.data)==0: return self.copy()
        if len(self.data)==0: return self.copy()
        for o in self.data: o._setFlag = 1
        for o in right.data: o._setFlag = 0
        newlist = filter( lambda x: x._setFlag==1, self.data )
        for o in self.data:
            if hasattr(o, '_setFlag'):
                del o._setFlag
        for o in right.data:
            if hasattr(o, '_setFlag'):
                del o._setFlag
        if self.stringRepr and right.stringRepr:
            stringRepr = self.stringRepr+'/-/'+right.stringRepr
        elif verbose:
            import traceback
            traceback.print_stack()
            print 'subtracting sets with no stringRepr:', repr(self), repr(right)
            #stringRepr = None
        return self.__class__(newlist, stringRepr=stringRepr)


    def __sub__(self, right):
        """See subtract: overloads - operator"""
        return self.subtract(right)
 
        
    def makeUniq(self):
        """removes duplicates from set (in place)"""
        l = []
        d = {}
        for value in self.data:
            if not d.has_key(id(value)):
                d[id(value)]=value
                l.append(value)
        self.__dict__['data'] = l


    def uniq(self):    # Fastest order preserving
        set = {}
        
        return self.__class__([set.setdefault(e,e) for e in self.data if e not in set])

##     def uniq(self):
##         l = []
##         d = {}
##         for value in self.data:
##             # Here we use the id(value) as a key rather than value itself
##             # because if self.data is a TreeNodeSet each time
##             # we test if d.has_key(value) it calls the __hash__ method
##             # of tree.py which has been overwritten
##             # Function calls in python are costly.
##             if not d.has_key(id(value)):
##                 d[id(value)]=value
##                 l.append(value)
##         return self.__class__(l)

