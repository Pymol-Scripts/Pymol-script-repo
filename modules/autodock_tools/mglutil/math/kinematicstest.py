## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#
# Last modified on Thu Jan 24 16:14:31 PST 2002 by lindy
#
# $Header: /opt/cvs/python/packages/share1.5/mglutil/math/kinematicstest.py,v 1.8 2007/07/24 17:30:40 vareille Exp $
#

"""Unit test for kinematics.py

Requirements for kinematics.py:
    A. Kinematics.__init__:
        1. make Numeric, homogenious coordinates out of refCoords
        2. raise ValueError if refCoords is bad
    B. Kinematics.reset():
        3. nx3 slice of resultCoords must be equal to refCoords
    C. Kinematic.getResultCoords:
        4. return nx3 (not nx4) coordinates
        5. return as Numeric.array or ListType accorinding to self.tolist

    D. Kinematics.setTorTree:
        0. should ...
    E. Kinematics.applyState:
        0. return the coordinates correctly
        0. raise ValueError for bad State parameters
    F. Kinematics.applyAngList:
        0. return coorectly transformed coordinates
        0. rasie ValueError if len(angList) != len(self.torsions)
    G. Kinematis.applyTorsion:
        0. return the coordinates correctly
    H. Kinematics.applyOrientation
        0. return the coordinates correctly
    I. Kinematics.applyQuaternion
        0. return the coordinates correctly
    J. Kinematics.applyTranslation
        0. return the coordinates correctly
        
"""

from mglutil.math.kinematics import Kinematics
from mglutil.math.transformation import Transformation
import unittest, math
import numpy.oldnumeric as Numeric, numpy.oldnumeric.random_array as RandomArray
from UserList import UserList



class TestTorTree(UserList):
    def __init__(self, data=None):
        UserList.__init__(self, data)
        self.tList = self.data

    def setTorsionAngles(self, angList):
        pass



class TestTorsion:
    def __init__(self, pivot1, pivot2, points):
        self.atm1_ix = pivot1
        self.atm2_ix = pivot2
        self.atms_to_move = points



class KinematicsTest(unittest.TestCase):
    def setUp(self):
        """Called for every test."""
        self.decimals = 4 # for Numeric.around; 7 is SciPy default.
        self.idmtx = Transformation().getMatrix(transpose=1)

        self.known_points = [ [1., 0., 2.],
                              [1., 0., 1.],
                              [1., 1., 1.],
                              [0., 0., 1.],
                              [0., 0., 0.],
                              [0., 1., 0.],
                              [0., 1., -1.],
                              [1., 1., -1.],
                              [1., 2., -1.],
                              [1., 1., -2.]]
        npts = len(self.known_points)
        dim = 3
        self.max = 9999999.
        self.min = -self.max
        self.random_points = RandomArray.uniform(self.min,
                                                 self.max, (npts,dim)).tolist()

        # create a simple torsion system for both point lists
        torTree = TestTorTree()
        torTree.append(TestTorsion(4, 3, [0,1,2]))
        torTree.append(TestTorsion(3, 1, [0,2]))
        torTree.append(TestTorsion(6, 7, [8,9]))
        self.torTree = torTree


    def assertArrayEqual(self, a1, a2, decimals=None):
        """Round the arrays according to decimals and compare

        Use self.decimals if decimals is not supplied"""
        if decimals:
            d = decimals
        else:
            d = self.decimals

        for point in xrange(len(a1)):
            for axis in [0, 1, 2]:
                self.assertEqual(round(a1[point][axis],d),
                                 round(a2[point][axis],d))


    def tearDown(self):
        pass



class InputOutputValues(KinematicsTest):


    def test_setTorTree(self):
        """setTorTree       -- I/O validity untested"""
        self.assertEqual(0, 1)



class ComputedValues(KinematicsTest):


    def test_applyAngList00(self):
        """applyAngList     -- return correct coords (0., 0., 0.)"""
        ko = Kinematics(self.known_points, self.torTree, tolist=1)
        self.assertEqual(self.known_points, ko.applyAngList([0,0,0], self.idmtx))


    def test_applyAngList01(self):
        """applyAngList     -- return correct coords (360., 360., 360.)"""
        ko = Kinematics(self.known_points, self.torTree, tolist=0)
        result = ko.applyAngList([360.,360.,360.], self.idmtx)[:,:3]
        self.assertArrayEqual(self.known_points, result)


    def test_applyAngList020(self):
        """applyAngList     -- return correct coords (-90, 0, 0)"""
        ko = Kinematics(self.known_points, self.torTree, tolist=0)
        result = ko.applyAngList([-90.,0.,0.], self.idmtx)[:,:3]

        self.known_points[0] = [0., -1., 2.]
        self.known_points[1] = [0., -1., 1.]
        self.known_points[2] = [1., -1., 1.]
        self.assertArrayEqual(self.known_points, result)


    def test_applyAngList021(self):
        """applyAngList     -- return correct coords (90, 0, 0)"""
        ko = Kinematics(self.known_points, self.torTree, tolist=0)
        result = ko.applyAngList([90.,0.,0.], self.idmtx)[:,:3]

        self.known_points[0] = [0.,  1., 2.]
        self.known_points[1] = [0.,  1., 1.]
        self.known_points[2] = [-1., 1., 1.]
        self.assertArrayEqual(self.known_points, result)


    def test_applyAngList022(self):
        """applyAngList     -- return correct coords ( 0,-90, 0)"""
        ko = Kinematics(self.known_points, self.torTree, tolist=0)
        result = ko.applyAngList([0.,-90.,0.], self.idmtx)[:,:3]

        self.known_points[0] = [1., 1., 1.]
        self.known_points[2] = [1., 0., 0.]
        self.assertArrayEqual(self.known_points, result)


    def test_applyAngList03(self):
        """applyAngList     -- return correct coords (90, -90, 0)"""
        ko = Kinematics(self.known_points, self.torTree, tolist=0)
        result = ko.applyAngList([90.,-90.,0.], self.idmtx)[:,:3]

        self.known_points[0] = [-1.,  1., 1.]
        self.known_points[1] = [ 0.,  1., 1.]
        self.known_points[2] = [ 0.,  1., 0.]
        self.assertArrayEqual(self.known_points, result)


    def test_applyAngList04(self):
        """applyAngList     -- one-at-time equals all-at-once"""
        ko1 = Kinematics(self.known_points, self.torTree, tolist=0)
        ko2 = Kinematics(self.known_points, self.torTree, tolist=0)        
        result1 = ko1.applyAngList([90.,-90.,0.], self.idmtx)[:,:3]

        # do the rotations separately
        ko2.applyAngList([90.,0.,0.])
        result2 = ko2.applyAngList([0.,-90.,0.], self.idmtx)[:,:3]
        self.assertArrayEqual(result1, result2)


    def test_applyAngList05(self):
        """applyAngList     -- ValueError if len(angList) != len(torsions)"""
        ko = Kinematics(self.known_points, self.torTree, tolist=1)
        self.assertRaises(ValueError, ko.applyAngList, ([0,0], self.idmtx))


##      def test_applyTorsion00(self):
##          """applyAngList, applyTorsion give same results"""
##          ko = Kinematics(self.known_points, self.torTree, tolist=0)
##          result = Numeric.around(ko.applyAngList([90.,-90.,0.])[:,:3],
##                                  self.decimals).tolist()
##          jo = Kinematics(self.known_points, self.torTree, tolist=0)
##          jo.applyTorsion(jo.torsions[0], 90.)
##          jo.applyTorsion(jo.torsions[1], -90.)
##          #
##          self.known_points[0] = [-1.,  1., 1.]
##          self.known_points[1] = [ 0.,  1., 1.]
##          self.known_points[2] = [ 0.,  1., 0.]
##          self.assertEqual( self.known_points, result)


    def test_applyTorsion01(self):
        """applyTorsion     -- return correct coords torsion=0, angle=0."""
        ko = Kinematics(self.known_points, self.torTree, tolist=1)
        self.assertEqual(self.known_points, ko.applyTorsion( ko.torsions[0], 0.) )


    def test_applyTorsion02(self):
        """applyTorsion     -- return correct coords torsion=1, angle=90."""
        ko = Kinematics(self.known_points, self.torTree, tolist=0)
        result = ko.applyTorsion(ko.torsions[1], 90.)[:,:3]

        self.known_points[0] = [1., -1., 1.]
        self.known_points[2] = [1., 0., 2.]
        self.assertArrayEqual(self.known_points, result)


    def test_applyTorsion02(self):
        """applyTorsion     -- apply torsions in any order"""
        ko1 = Kinematics(self.known_points, self.torTree, tolist=0)
        ko1.applyTorsion(ko1.torsions[0], 90.)
        result1 = ko1.applyTorsion(ko1.torsions[1], -90.)

        ko2 = Kinematics(self.known_points, self.torTree, tolist=0)
        ko2.applyTorsion(ko2.torsions[1], -90.)
        result2 = ko2.applyTorsion(ko2.torsions[0], 90.)
        self.assertArrayEqual(result1, result2)


if __name__ == '__main__':
    unittest.main()   

# for example:
#     py mglutil/math/kinematicstest.py -v
# or, to redirect output to a file:
#     py kinematicstest.py -v > & ! /tmp/kt.out






