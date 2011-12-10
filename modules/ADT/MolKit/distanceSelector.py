## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#############################################################################
#
# Author: Ruth Huey
#
# Copyright: M. Sanner TSRI 2001
#
#############################################################################


# $Header: /opt/cvs/python/packages/share1.5/MolKit/distanceSelector.py,v 1.10 2007/07/24 17:30:40 vareille Exp $
#
# $Id: distanceSelector.py,v 1.10 2007/07/24 17:30:40 vareille Exp $
#

"""
This module implements classes which compare coordinates of two sets of atoms.
They return two dictionaries which have keys from the first set of atoms. The 
first of these dictionaries is 'pairDict' which has values of lists of atoms 
in the second set and the second 'distDict' which has values of distances 
corresponding to key-value pairs in pairDict.
"""

import numpy.oldnumeric as Numeric, types, string, math
from molecule import AtomSet
from protein import ResidueSet



class DistanceSelector:
    """ Class that detects atoms closer than cutoff distance.
    """


    def __init__(self, return_dist=1, **kw):
        self.func = Numeric.less
        self.return_dist = return_dist


    def setupCutoff(self, bAts, sAts, cutoff):
        # when called len(bAts)>len(sAts)
        llen = len(bAts)
        slen = len(sAts)
        #test whether cutoff is 1 value or an array
        if not hasattr(cutoff, 'shape') and type(cutoff)==types.FloatType:
            #build array with shape (llen, slen)  of value cutoff
            return Numeric.array(Numeric.ones((llen, slen))*cutoff).astype('f')
        else:
            #might need to reshape cutoff ?
            assert cutoff.shape==(llen, slen)
            return cutoff


    def mul(self, coords, mat):
        shape = coords.shape
        assert shape[-1]==3
        #coords = Numeric.reshape(coords, (-1, shape[-1]))
        one = Numeric.ones((shape[0],1),'f')
        c = Numeric.concatenate((coords, one),1)
        coords = Numeric.array(Numeric.reshape(coords, shape))
        return Numeric.array(Numeric.dot(c, Numeric.transpose(mat))[:,:3])
        

    def select(self, keyAts, checkAts, cutoff=3.0, percentCutoff=1.0, 
                keyMat=None, checkMat=None):
        """ keyAts, checkAts, cutoff, percentCutoff
            keyAts: first set of atoms
            checkAts: a second set of atoms which is checked vs. keyAts
            cutoff: 
                either a single float by default 3.0 
                or a matrix with shape:
            (max(len(keyAts),len(checkAts)), min(len(keyAts),len(checkAts)))
            percentCutoff: by default 1.0 (cutoff is multiplied by this value)
            keyMat: transformation of keyAts
            checkMat: transformation of checkAts

        returns 'pairDict' whose keys are atoms used as reference points and whose 
        values are atoms within cutoff distance of corresponding key. 

        If 'return_dist' flag is set, 
        'distDict' is returned also, whose keys are 
        the same atoms which are used as reference points and whose values are 
        lists of distances to atoms  within cutoff distance of corresponding key

        """
        lenK = len(keyAts)
        lenC = len(checkAts)

        #data arrays are used to find atoms with given indices quickly
        atar = Numeric.array(checkAts.data)
        keyAtar = Numeric.array(keyAts.data)

        #basic arrays of coords used to build others
        c = Numeric.array(checkAts.coords, 'f')
        if checkMat:
            c = self.mul(c, checkMat)
        k = Numeric.array(keyAts.coords, 'f')
        if keyMat:
            k = self.mul(k, keyMat)

        # first build matrix of distances between all pairs of ats
        # rows correspond to ats in larger set, columns to those in smaller
        # first build square matrix
        if lenC >= lenK:
            bigC = Numeric.resize(c, (lenC, lenC, 3))
            k.shape = (lenK,1,3)
            bigM = bigC[:lenK]
            smallM = k

            cutoff = self.setupCutoff(checkAts, keyAts, cutoff)
            #print "0a:cutoff[0][0]=", cutoff[0][0]
            cutoff.shape = (lenK, -1)

        else:
            bigK = Numeric.resize(k, (lenK, lenK, 3))
            c.shape = (lenC,1,3)
            bigM = bigK[:lenC]
            smallM = c
            cutoff = self.setupCutoff(keyAts, checkAts, cutoff)
            #print "0b:cutoff[0][0]=", cutoff[0][0]
            cutoff.shape = (lenC, -1)

        # distance matrix
        d = bigM - smallM
        # distance squared matrix
        dSQ = d * d
        # next step sums deltaX**2, deltaY**2, deltaZ**2
        dSQMAT = Numeric.sum(dSQ,2)

        #percentCutoff lets user relax sum of radii
        #the smaller the percentCutoff the smaller the key 
        #dSQ has to be less than
        cutoff = cutoff * percentCutoff
        cutoffSQMAT = cutoff * cutoff
        #cutoffSQMAT = cutoffSQMAT * percentCutoff

        # ansMat has 1 where sq dist. is smaller than cutoff
        ansMat = Numeric.logical_and(self.func(dSQMAT, cutoffSQMAT) , \
                    Numeric.not_equal(dSQMAT, 0.))
        if lenK > lenC:
            # in this case need to rearrange matrix
            # which got shuffled in if-else above
            ansMat = Numeric.swapaxes(ansMat, 0, 1)
            dSQMAT = Numeric.swapaxes(dSQMAT, 0, 1)

        # finally, build result dictionaries which have atom keys:
        #   pairDict has values which are lists of close atoms
        #   distDict has values which are lists of distances
        pairDict = {}
        distDict = {}
        # get a list of rows which have non-zero entries 
        # to loop over in next section
        rowIndices = Numeric.nonzero(Numeric.sum(ansMat,1))
        # rows correspond to ats in keyAts
        # columns correspond to ats in checkAts
        for i in rowIndices:
            # atindex is a list [7 8 9] indexing into checkAts 
            atindex = Numeric.nonzero(ansMat[i])
            # keyAtar[i] is ith atom in keyAts
            keyAt = keyAtar[i]
            pairDict[keyAt] = Numeric.take(atar, atindex)
            if self.return_dist:
                distDict[keyAt] = []
                for ind in atindex:
                    distDict[keyAt].append(math.sqrt(dSQMAT[i][ind]))
        
        #getting distDict back is optional
        if self.return_dist: return pairDict, distDict
        else: return pairDict


    def dictToResidues(self, pairDict):
        """ 
        The input, 'pairDict',  is a dictionary with atoms as keys and lists of atoms close to each key atom as values.  This method converts these keys and values to residue sets. The first returned value is a unique ResidueSet of the 'control' or 'key' atoms used for the distance selection and the second a unique ResidueSet of parents of atoms close to the keys
        """ 
        if not len(pairDict):
            return ResidueSet(), ResidueSet()
        from MolKit.molecule import AtomSet
        #parents of keys
        key_parents = AtomSet(pairDict.keys()).parent.uniq()
        key_parents.sort()
        #parents of values
        #build unique list of close atoms
        ats = {}
        for k in pairDict.keys():
            for rec_at in pairDict[k]:
                ats[rec_at] = 1
        close_ats = AtomSet(ats.keys())
        close_parents = close_ats.parent.uniq()
        close_parents.sort()
        return key_parents, close_parents
        


class FartherThanSelector(DistanceSelector):
    """ Class that detects atoms farther apart than cutoff distance.
    """

    def __init__(self, return_dist=1):
        DistanceSelector.__init__(self)
        self.func = Numeric.greater 
        self.return_dist = return_dist



class CloserThanVDWSelector(DistanceSelector):
    """ Class that detects atoms closer than sum of van der waals radii
        distance.
    """


    def setupCutoff(self, bAts, sAts, cutoff):
        #len(bAts)>len(sAts)
        lenC = len(bAts)
        lenK = len(sAts)
        #set up sum of radii cutoff matrix:
        checkRadii = Numeric.array(bAts.vdwRadius, 'f')
        bigRC = Numeric.resize(checkRadii, (lenC, lenC))
        bigR = bigRC[:lenK]

        keyRadii = Numeric.array(sAts.vdwRadius, 'f')
        keyRadii.shape = (lenK, 1)
        smallR = keyRadii
        cutoff = bigR + smallR
        return cutoff

 
class CloserThanVDWPlusConstantSelector(CloserThanVDWSelector):
    """ Class that detects atoms closer than sum of van der waals radii
        distance plus a constant.
    """


    def __init__(self, return_dist=1, constant=1):
        CloserThanVDWSelector.__init__(self)
        self.constant = constant



    def setupCutoff(self, bAts, sAts, cutoff):
        #len(bAts)>len(sAts)
        lenC = len(bAts)
        lenK = len(sAts)
        #set up sum of radii cutoff matrix:
        checkRadii = Numeric.array(bAts.vdwRadius, 'f')
        bigRC = Numeric.resize(checkRadii, (lenC, lenC))
        bigR = bigRC[:lenK]

        keyRadii = Numeric.array(sAts.vdwRadius, 'f')
        keyRadii.shape = (lenK, 1)
        smallR = keyRadii
        cutoff = bigR + smallR
        #create Numeric.ones of appropriate shape * self.const
        #add it to cutoff
        constFactor = Numeric.ones(len(bigR.ravel()))*self.constant
        #print "len(constFactor)=", len(constFactor)
        constFactor.shape = cutoff.shape
        #print "constFactor.shape=" , constFactor.shape
        #print "1:cutoff[0][0]=", cutoff[0][0]
        cutoff = cutoff+constFactor
        #print "2:cutoff[0][0]=", cutoff[0][0]
        return cutoff

