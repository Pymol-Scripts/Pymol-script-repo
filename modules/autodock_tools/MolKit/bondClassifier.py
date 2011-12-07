#############################################################################
#
# Author: Ruth Huey, William M. Lindstrom
#
# Copyright: R. Huey, W. M. Lindstrom TRSI 2003
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/MolKit/bondClassifier.py,v 1.3 2004/04/06 23:29:50 rhuey Exp $
#
# $Id: bondClassifier.py,v 1.3 2004/04/06 23:29:50 rhuey Exp $
#
#
#

"""
This module implements a classifier which select bonds based on a 
dictionary of key, bondSelector.
It returns  a dictionary with keys the specified bond types and 
values the bonds which have been classified.
"""

from MolKit.molecule import BondSet



class BondClassifier:
    """ Base class that sorts bonds based on an input dictionary with keys
    and bondSelector values
    """


    def __init__(self, d={}):
        self.dict = d


    def classify(self, bonds=None):
        """ 
        select using each bondselector (the values of the dict); store result
        in resultDict and return the result dict when finished...
        """
        #make sure that classify is called with some bonds
        assert isinstance(bonds, BondSet)
        resultDict = {}
        for k, v in self.dict.items():
            resultDict[k] = v.select(bonds)
        return resultDict

