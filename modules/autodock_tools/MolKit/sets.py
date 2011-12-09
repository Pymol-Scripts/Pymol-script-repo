#############################################################################
#
# Authors: Michel F. SANNER, Ruth Huey
#
# Copyright: M. Sanner TSRI 2005
#
#############################################################################
#
# $Header: /opt/cvs/python/packages/share1.5/MolKit/sets.py,v 1.2 2006/06/12 18:28:35 sargis Exp $
#
# 
#
#  
#

from tree import TreeNodeSet
import types

class Sets(dict):
    """
Object used to manage a collection of explicit sets of TreeNodes
"""

    def add(self, name, set, overwrite=True):
        assert isinstance(set, TreeNodeSet)
        assert type(name) in types.StringTypes
               
        if not overwrite:
            assert name not in self.keys()
        self[name] = set


    def remove(self, name):
        # remove a set by name. Silently ignore non existing sets but returns
        # true when a set gets deleted, else returns False
        if name in self.keys():
            del self[name]
            return True
        return False


    def removeByInstance(self, set):
        # remove a set that is specified by a TreeNodeSet.
        # Silently ignore non existing sets but returns
        # true when a set gets deleted, else returns False
        for n,s in self.items():
            if s==set:
                del self[n]
                return True
        return False


    def get(self, stype=None):
        # return a dict of sets optionally restricted to a user specified type
        # if stype is specified it has to be a subclass of TreeNodeSet
        if stype is None:
            return self
        else:  # select the sets of a given type
            assert issubclass(stype, TreeNodeSet)
            result = {}
            for name,set in self.items():
                if isinstance(set, stype):
                    result[name] = set
            return result
        
