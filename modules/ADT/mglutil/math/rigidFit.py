## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#############################################################################
#
# Author: Sophie I. COON, William LINDSTROM, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

import numpy.oldnumeric as Numeric, math
import numpy.oldnumeric.linear_algebra as LinearAlgebra
from mglutil.math.rmsd import RMSDCalculator

class RigidfitBodyAligner:
    """
    This class implements a set of method to compute transformation matrices
    to superimpose a set of mobile coordinates onto a set of reference
    coordinates, apply the resulting transformation matrices to a given set of
    coordinates and finally compute the RMSD and the distance vectors between
    the set of reference coordinates and a given set of coordinates transformed by
    the computed transformation matrices.
    """
    def __init__(self, refCoords=None):
        """The constructor of the rigidfitBodyAligner takes one required argument:
        refCoords: cartesian coordinates of reference structure (3,n) (input)
        This method creates:
        self.superimposed: flag when set to 1 the transformation matrices to superimpose
                           a mobile set of coordinates onto the reference coords have been
                           computed.
        self.rmsdCalculator: instance of the RMSDCalculator class that computes the rmsd
                             between two lists of coordinates. This object will be used by
                             the rmsdAfterSuperimposition method.
        """
        self.refCoords = refCoords
        self.superimposed = 0
        self.rmsdCalculator = RMSDCalculator()


    def setRefCoords(self, refCoords):
        """ The setRefCoords method allows you to lists the reference coordinates."""
        self.refCoords = refCoords
        # New set of ref coords so set the flag back to 0:
        self.superimposed = 0
        

    def rigidFit(self, mobileCoords):
        """
        the rigidFit method computes the necessary
        transformation matrices to superimpose the list of mobileCoords
        onto the list of referenceCoords, and stores the resulting matrices

        (rot, trans) <- rigidFit(mobileCoords)
        Rigid body fit. Finds transformation (rot,trans) such that
        r.m.s dist(x,rot*y+trans) --> min !
        mobileCoords: cartesian coordinates of mobile structure (3,n) (input)
        rot   : rotation matrix (3,3) (output)
        trans : translation vector (3) (output)
        status: 0 if OK, 1 if singular problem (n<3 ...)

        Method: W.Kabsch, Acta Cryst. (1976). A32,922-923
        W.Kabsch, Acta Cryst. (1978). A34,827-828
        """
        if self.refCoords is None:
            raise RuntimeError(" no reference coordinates specified")
        
        refCoords = self.refCoords
        if len(refCoords) != len(mobileCoords):
            raise RuntimeError("input vector length mismatch")

        refCoords = Numeric.array(refCoords)
        mobileCoords = Numeric.array(mobileCoords)

        #
        # 1. Compute centroids:
        refCentroid = Numeric.sum(refCoords)/ len(refCoords)
        mobileCentroid = Numeric.sum(mobileCoords)/ len(mobileCoords)

        #
        # 2. Wolfgang Kabsch's method for rotation matrix rot:
        rot = Numeric.identity(3).astype('f')
        # LOOK how to modify that code.
        for i in xrange(3):
            for j in xrange(3):
                rot[j][i] = Numeric.sum((refCoords[:,i]-refCentroid[i])*
                                        (mobileCoords[:,j]-mobileCentroid[j]))

        rotTransposed = Numeric.transpose(rot)
        e = Numeric.dot(rot, rotTransposed)

        evals, evecs = LinearAlgebra.eigenvectors(e)

        ev = Numeric.identity(3).astype('d')
        # set ev[0] to be the evec or the largest eigenvalue
        # and ev[1] to be the evec or the second largest eigenvalue
        eigenValues = list(evals)
        discard = eigenValues.index(min(eigenValues))
        i = j =0
        while i<3:
            if i==discard:
                i = i + 1
                continue
            ev[j] = evecs[i]
            j = j + 1
            i = i + 1
        evecs = ev

        evecs[2][0] = evecs[0][1]*evecs[1][2] - evecs[0][2]*evecs[1][1]
        evecs[2][1] = evecs[0][2]*evecs[1][0] - evecs[0][0]*evecs[1][2]
        evecs[2][2] = evecs[0][0]*evecs[1][1] - evecs[0][1]*evecs[1][0]

        b = Numeric.dot(evecs, rot)

        norm = math.sqrt(b[0][0]*b[0][0] + b[0][1]*b[0][1] + b[0][2]*b[0][2])
        if math.fabs(norm)<1.0e-20: return -1, -1
        b[0] = b[0]/norm

        norm = math.sqrt(b[1][0]*b[1][0] + b[1][1]*b[1][1] + b[1][2]*b[1][2])
        if math.fabs(norm)<1.0e-20: return -1, -1
        b[1] = b[1]/norm

        # vvmult(b[0],b[1],b[2])
        b[2][0] = b[0][1]*b[1][2] - b[0][2]*b[1][1]
        b[2][1] = b[0][2]*b[1][0] - b[0][0]*b[1][2]
        b[2][2] = b[0][0]*b[1][1] - b[0][1]*b[1][0]
        # mtrans3(e)
        e = evecs
        tempo=e[0][1]; e[0][1]=e[1][0]; e[1][0]=tempo
        tempo=e[0][2]; e[0][2]=e[2][0]; e[2][0]=tempo
        tempo=e[1][2]; e[1][2]=e[2][1]; e[2][1]=tempo
        # mmmult3(b,e,rot)
        rot[0][0] = b[0][0]*e[0][0] + b[1][0]*e[0][1] + b[2][0]*e[0][2]
        rot[0][1] = b[0][1]*e[0][0] + b[1][1]*e[0][1] + b[2][1]*e[0][2]
        rot[0][2] = b[0][2]*e[0][0] + b[1][2]*e[0][1] + b[2][2]*e[0][2]

        rot[1][0] = b[0][0]*e[1][0] + b[1][0]*e[1][1] + b[2][0]*e[1][2]
        rot[1][1] = b[0][1]*e[1][0] + b[1][1]*e[1][1] + b[2][1]*e[1][2]
        rot[1][2] = b[0][2]*e[1][0] + b[1][2]*e[1][1] + b[2][2]*e[1][2]

        rot[2][0] = b[0][0]*e[2][0] + b[1][0]*e[2][1] + b[2][0]*e[2][2]
        rot[2][1] = b[0][1]*e[2][0] + b[1][1]*e[2][1] + b[2][1]*e[2][2]
        rot[2][2] = b[0][2]*e[2][0] + b[1][2]*e[2][1] + b[2][2]*e[2][2]

        #
        # Compute translation vector trans:
        # mvmult3(rot,cy,cy);
        trans3 = [0, 0, 0]
        for i in range(3):
            trans3[i] = mobileCentroid[0]*rot[0][i] + mobileCentroid[1]*rot[1][i] + \
                        mobileCentroid[2]*rot[2][i]

        #bcopy(t3,cy,sizeof(t3));
        #vvdiff(cx,cy,trans);
        trans = ( refCentroid[0]-trans3[0], refCentroid[1]-trans3[1], refCentroid[2]-trans3[2] )
        #
        #   That's it...

        self.rotationMatrix = rot
        self.translationMatrix = trans
        self.superimposed = 1
            
    def transformCoords(self, setCoords):
        """ The transformCoords method applies the transformation matrices
        computed by rigidFit method to the given list of coordinates.
        """
        # This can only be done if the transformation matrices have been computed
        # by the rigidFit method.
        if not self.superimposed: return

        # 1- apply the rotation and the translation matrix to the given set of coords.
        transfoMatrix = Numeric.identity(4, 'd')
        transfoMatrix[:3,:3] = self.rotationMatrix
        transfoMatrix[3][0] = self.translationMatrix[0]
        transfoMatrix[3][1] = self.translationMatrix[1]
        transfoMatrix[3][2] = self.translationMatrix[2]
        
        # 2- now apply the transformation to the list of given coordinates list:
        # make homogeneous coords
        homoCoords = Numeric.concatenate((setCoords,
                                          Numeric.ones( (len(setCoords), 1), 'd')), 1)
        # 3- apply the transformation matrix to the homogeneous coords.
        transformedCoords = Numeric.dot(homoCoords, transfoMatrix)

        return transformedCoords
        
    def rmsdAfterSuperimposition(self, setCoords):
        """ The computeRMSD method computes the overall root mean square
        distance (rmsd) and also the distance between each pair of points
        (distVect) between the refCoords and the given list of Coords.
        The transformation matrices will be applied to the given list of coords."""
        if not self.superimposed: return

        transformedCoords = self.transformCoords(setCoords)
        self.rmsdCalculator.setRefCoords(self.refCoords)
        #return self.rmsdCalculator.computeRMSD(setCoords)
        return self.rmsdCalculator.computeRMSD(transformedCoords[:,:3].tolist())

    
if __name__ == '__main__':
    import numpy.oldnumeric as Numeric

    a = [ [ 2.92666,  -5.8288,  4.7963 ],
          [ 3.4784 ,  -7.6197,  4.9499 ],
          [-10.6508,  -4.7754,  1.1695 ],
          [ -9.4381,  -4.5545,  1.18   ] ]

    b = [ [-10.6508,  -4.7754, 1.1695 ],
          [ -9.4381,  -4.5545,   1.18 ],
          [  2.9266,  -5.8288,  4.7963],
          [  3.4784,  -7.6197,  4.9499] ] 

    a = [ [ 0.,  0.,  0. ],
          [ 10.,  0.,  0. ],
          [ 0.,  10.,  0. ],
          [ 0.,  0.,  10. ] ]

    b = [ [ 0.,  0.,  0. ],
          [ 10.,  0.,  0. ],
          [ 0.,  0.,  -10. ],
          [ 0.,  10.,  0. ] ]

    #o = RigidfitBodyAligner(a)
    #o.rigidFit(b)
    #print o.rotationMatrix
    #print o.translationMatrix

#[[-0.91369569 -0.13739982 -0.38240376]
# [-0.13729294 -0.7812196   0.60899734]
# [-0.38244215  0.60897326  0.69491529]]
#(-6.17204211132, -12.4548639806, 3.08235584151)
