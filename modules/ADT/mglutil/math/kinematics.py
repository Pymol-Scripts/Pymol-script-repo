## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#
# Last modified on Mon Oct 15 15:33:49 PDT 2001 by lindy
#
# $Header: /opt/cvs/python/packages/share1.5/mglutil/math/kinematics.py,v 1.16 2007/07/24 17:30:40 vareille Exp $
#

"""kinematics.py - kinematic manipulation of chains of points

All transformations happen in the local coordinate space.
The refCoords supplied to the constructor and returned by the object
are local to the object. Clients should handle putting the points into
world coordinates (using translation, orientation, and origin).
"""


#from mglutil.math.ncoords import Ncoords
from mglutil.math.rotax import rotax
import numpy.oldnumeric as Numeric, math


class Kinematics:
    rads_per_degree = Numeric.pi/180.

    def __init__(self, allAtomsCoords, torTree, tolist=1):
        """refCoords is an nx3 list of n points
        
        resultCoords is set up and maintained as homogeneous coords
        """
        self.allAtomsCoords = allAtomsCoords
        self.torTree = torTree

    def __applyTorsion(self, node, parent_mtx):
        """Transform the subtree rooted at node.

        The new torsion angle must be pre-set.
        Children of the node are transformed recursively.
        """
        # get rotation matrix for node
        # my_mtx = self.rotax(node)
        mtx = rotax( Numeric.array(node.a.coords),
                     Numeric.array(node.b.coords),
                     node.angle * self.rads_per_degree, transpose=1)

        # node_mtx = Numeric.dot(parent_mtx, mtx)
        node_mtx = self.mult4_3Mat(parent_mtx, mtx)

        # set-up for the transformation
        mm11 = node_mtx[0][0]; mm12 = node_mtx[0][1]; mm13 = node_mtx[0][2]
        mm21 = node_mtx[1][0]; mm22 = node_mtx[1][1]; mm23 = node_mtx[1][2]
        mm31 = node_mtx[2][0]; mm32 = node_mtx[2][1]; mm33 = node_mtx[2][2]
        mm41 = node_mtx[3][0]; mm42 = node_mtx[3][1]; mm43 = node_mtx[3][2]
        atomSet = node.atomSet        
        # transform the coordinates for the node
        for i in node.atomRange:
            x,y,z = node.coords[i][:3] # get origin-subtracted originals
            c = atomSet[i].coords
            c[0] = x*mm11 + y*mm21 + z*mm31 + mm41
            c[1] = x*mm12 + y*mm22 + z*mm32 + mm42
            c[2] = x*mm13 + y*mm23 + z*mm33 + mm43

        # recurse through children
        for child in node.children:
            self.__applyTorsion(child, node_mtx)


    def applyAngList(self, angList, mtx):
        """
        """
        # pre-set the torsion angles
        self.torTree.setTorsionAngles(angList)

        # set-up for the transformation
        mm11 = mtx[0][0]; mm12 = mtx[0][1]; mm13 = mtx[0][2]
        mm21 = mtx[1][0]; mm22 = mtx[1][1]; mm23 = mtx[1][2]
        mm31 = mtx[2][0]; mm32 = mtx[2][1]; mm33 = mtx[2][2]
        mm41 = mtx[3][0]; mm42 = mtx[3][1]; mm43 = mtx[3][2]
        root = self.torTree.rootNode
        atomSet = root.atomSet        
        # transform the coordinates for the node
        for i in root.atomRange:
            x,y,z = root.coords[i][:3]
            c = atomSet[i].coords
            c[0] = x*mm11 + y*mm21 + z*mm31 + mm41
            c[1] = x*mm12 + y*mm22 + z*mm32 + mm42
            c[2] = x*mm13 + y*mm23 + z*mm33 + mm43

        # traverse children of rootNode
        for child in root.children:
            self.__applyTorsion(child, mtx)


    def mult4_3Mat(self, m1, m2):
        ma11 = m1[0][0]
        ma12 = m1[0][1]
        ma13 = m1[0][2]
        ma21 = m1[1][0]
        ma22 = m1[1][1]
        ma23 = m1[1][2]
        ma31 = m1[2][0]
        ma32 = m1[2][1]
        ma33 = m1[2][2]
        ma41 = m1[3][0]
        ma42 = m1[3][1]
        ma43 = m1[3][2]

        mb11 = m2[0][0]
        mb12 = m2[0][1]
        mb13 = m2[0][2]
        mb21 = m2[1][0]
        mb22 = m2[1][1]
        mb23 = m2[1][2]
        mb31 = m2[2][0]
        mb32 = m2[2][1]
        mb33 = m2[2][2]
        mb41 = m2[3][0]
        mb42 = m2[3][1]
        mb43 = m2[3][2]

        # first line of resulting matrix
        val1 = ma11*mb11 + ma12*mb21 + ma13*mb31
        val2 = ma11*mb12 + ma12*mb22 + ma13*mb32
        val3 = ma11*mb13 + ma12*mb23 + ma13*mb33
        result = [[val1, val2, val3, 0.0]]

        # second line of resulting matrix
        val1 = ma21*mb11 + ma22*mb21 + ma23*mb31
        val2 = ma21*mb12 + ma22*mb22 + ma23*mb32
        val3 = ma21*mb13 + ma22*mb23 + ma23*mb33
        result.append([val1, val2, val3, 0.0])

        # third line of resulting matrix
        val1 = ma31*mb11 + ma32*mb21 + ma33*mb31
        val2 = ma31*mb12 + ma32*mb22 + ma33*mb32
        val3 = ma31*mb13 + ma32*mb23 + ma33*mb33
        result.append([val1, val2, val3, 0.0])

        # fourth line of resulting matrix
        val1 = ma41*mb11 + ma42*mb21 + ma43*mb31 + mb41
        val2 = ma41*mb12 + ma42*mb22 + ma43*mb32 + mb42
        val3 = ma41*mb13 + ma42*mb23 + ma43*mb33 + mb43
        result.append([val1, val2, val3, 1.0])

        return result


    def rotax(self, node):
        """
        Build 4x4 matrix of clockwise rotation about axis a-->b
        by angle tau (radians).
        a and b are numeric arrys of floats of shape (3,)
        Result is a homogenous 4x4 transformation matrix.

        NOTE: This has been changed by Brian, 8/30/01: rotax now returns
        the rotation matrix, _not_ the transpose. This is to get
        consistency across rotax, mat_to_quat and the classes in
        transformation.py
        """
        tau = node.angle * self.rads_per_degree
        ct = math.cos(tau)
        ct1 = 1.0 - ct
        st = math.sin(tau)
        v = node.torUnitVector

        rot = Numeric.zeros( (4,4), 'f' )
        # Compute 3x3 rotation matrix

        v2 = v*v
        v3 = (1.0-v2)*ct
        rot[0][0]=v2[0]+v3[0]
        rot[1][1]=v2[1]+v3[1]
        rot[2][2]=v2[2]+v3[2]
        rot[3][3] = 1.0;

        v2 = v*st
        rot[1][0]=v[0]*v[1] * ct1-v2[2]
        rot[2][1]=v[1]*v[2] * ct1-v2[0]
        rot[0][2]=v[2]*v[0] * ct1-v2[1]
        rot[0][1]=v[0]*v[1] * ct1+v2[2]
        rot[1][2]=v[1]*v[2] * ct1+v2[0]
        rot[2][0]=v[2]*v[0] * ct1+v2[1]

        # add translation
        a = node.torBase.coords
        print "    torBase (%2d) %4f, %4f, %4f:" % (node.bond[0], a[0], a[1], a[2])
        for i in (0,1,2):
            rot[3][i] = a[i]
            for j in (0,1,2): rot[3][i] = rot[3][i]-rot[j][i]*a[j]
            rot[i][3]=0.0

        return rot


