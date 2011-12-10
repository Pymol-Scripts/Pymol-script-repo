## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#############################################################################
#
# Author: Sophie I. COON, William LINDSTROM, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################
#
# Last modified on Thu Aug  9 20:00:31 PDT 2001 by lindy
#
# $Header: /opt/cvs/python/packages/share1.5/mglutil/math/rmsd.py,v 1.6 2007/07/24 17:30:40 vareille Exp $
#

import numpy.oldnumeric as Numeric
import math


class RMSDCalculator:
    """
    This class implements method to compute RMSD and distance vector
    between two given lists of coordinates.
    """
    def __init__(self, refCoords = None):
        self.refCoords = refCoords

    def setRefCoords(self, refCoords):
        self.refCoords = refCoords
        
    def computeRMSD(self, listCoords):
        """rmsd <- computRMSD(listCoords)
        rmsd returns the overall root mean square distance (rmsd) and
        also sets self.distVect as the vector of distances between each
        pair of points.
        """
        if self.refCoords is None:
            raise ValueError("no reference coordinates set")
        if len(self.refCoords) != len(listCoords):
            raise ValueError("input vector length mismatch")

        deltaVect = Numeric.array(self.refCoords) - Numeric.array(listCoords)
        distSquaredVect = Numeric.sum(Numeric.transpose(deltaVect*deltaVect))
        self.distVect = Numeric.sqrt(distSquaredVect)
        self.rmsd = math.sqrt(Numeric.sum(distSquaredVect)/len(self.refCoords))
        return self.rmsd
