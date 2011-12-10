## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#
# Last modified on Thu Sep  6 12:14:22 PDT 2001 by lindy
#
# $Header: /opt/cvs/python/packages/share1.5/mglutil/math/statetocoordstest.py,v 1.3 2007/07/24 17:30:40 vareille Exp $
#

"""Unit test for state.py
"""

from mglutil.math.statetocoords import StateToCoords
import unittest, math
import numpy.oldnumeric as Numeric, numpy.oldnumeric.random_array as RandomArray
from UserList import UserList



class TestState:
    def __init__(self, q, t, o, tList):
        self.quaternion = q
        self.translation = t
        self.origin = o
        self.torsions = tList



class TestTorTree(UserList):
    def __init__(self, data=None):
        UserList.__init__(self, data)
        self.tList = self.data



class TestTorsion:
    def __init__(self, pivot1, pivot2, points):
        self.atm1_ix = pivot1
        self.atm2_ix = pivot2
        self.atms_to_move = points


        
class StateToCoordsTest(unittest.TestCase):
    def setUp(self):
        """Called for every test."""
        self.decimals = 2 # for rounding; 7 is SciPy default.

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
        # create a simple torsion system for known_points
        torTree = TestTorTree()
        torTree.append(TestTorsion(4, 3, [0,1,2]))
        torTree.append(TestTorsion(3, 1, [0,2]))
        torTree.append(TestTorsion(6, 7, [8,9]))
        self.torTree = torTree

        npts = 5
        dim = 3
        self.max = 9999.
        self.min = -self.max
        self.random_points = RandomArray.uniform(self.min,
                                                 self.max, (npts,dim)).tolist()


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



class ComputedValues(StateToCoordsTest):


    def test_applyState(self):
        """applyState         -- computed values untested"""
        value = TestState(q=(1.0, 0.0, 0.0, 90.0),
                          t=(10., 10., 10.),
                          o=(1., 0., 1.),
                          tList=[0., 0., 270.])
        state = StateToCoords(self.known_points, self.torTree, tolist=1)
        result = state.applyState(value)
        expected = [[11.,  9., 11.],
                    [11., 10., 11.],
                    [11., 10., 12.],
                    [10., 10., 11.],
                    [10., 11., 11.],
                    [10., 11., 12.],
                    [10., 12., 12.],
                    [11., 12., 12.],
                    [11., 13., 12.],
                    [11., 12., 11.]]
        self.assertArrayEqual(expected, result)


    def test_applyOrientation00(self):
        """applyOrientation00 -- random pts with defaults"""
        state = StateToCoords(self.random_points, tolist=1)
        self.assertEqual(self.random_points, state.applyOrientation())


    def test_applyOrientation01(self):
        """applyOrientation01 -- Q, T, O combined"""
        state = StateToCoords(self.known_points, tolist=1)
        result = state.applyOrientation( quat=(1.0, 0.0, 0.0, 90.0),
                                         trans= (10., 10., 10.),
                                         origin=(1., 0., 1.))
        expected = [[11.,  9., 11.],
                    [11., 10., 11.],
                    [11., 10., 12.],
                    [10., 10., 11.],
                    [10., 11., 11.],
                    [10., 11., 12.],
                    [10., 12., 12.],
                    [11., 12., 12.],
                    [11., 12., 13.],
                    [11., 13., 12.]]
        self.assertArrayEqual(expected, result)


    def test_applyQuaternion01(self):
        """applyQuaternion01  -- known pts 360 about x-axis"""
        state = StateToCoords(self.known_points, tolist=1)
        result = state.applyQuaternion((1.0, 0.0, 0.0, 360.0))
        self.assertArrayEqual(self.known_points, result)


    def test_applyQuaternion02(self):
        """applyQuaternion02  -- known pts  90 about x-axis"""
        state = StateToCoords(self.known_points, tolist=1)
        result = state.applyQuaternion((1.0, 0.0, 0.0, 90.0))
        expected = [[1., -2., 0.],
                    [1., -1., 0.],
                    [1., -1., 1.],
                    [0., -1., 0.],
                    [0.,  0., 0.],
                    [0.,  0., 1.],
                    [0.,  1., 1.],
                    [1.,  1., 1.],
                    [1.,  1., 2.],
                    [1.,  2., 1.]]
        self.assertArrayEqual(expected, result)


    def test_applyQuaternion021(self):
        """applyQuaternion021 -- known pts  90 about x-axis; origin(1., 0., 1.)"""
        state = StateToCoords(self.known_points, tolist=1)
        result = state.applyQuaternion((1.0, 0.0, 0.0, 90.0), (1., 0., 1.))
        expected = [[1., -1., 1.],
                    [1.,  0., 1.],
                    [1.,  0., 2.],
                    [0.,  0., 1.],
                    [0.,  1., 1.],
                    [0.,  1., 2.],
                    [0.,  2., 2.],
                    [1.,  2., 2.],
                    [1.,  2., 3.],
                    [1.,  3., 2.]]
        self.assertArrayEqual(expected, result)


    def test_applyQuaternion03(self):
        """applyQuaternion02  -- known pts  90 about z-axis"""
        state = StateToCoords(self.known_points, tolist=1)
        result = state.applyQuaternion((0.0, 0.0, 1.0, 90.0))
        expected = [[ 0.,  1., 2.],
                    [ 0.,  1., 1.],
                    [-1.,  1., 1.],
                    [ 0.,  0., 1.],
                    [ 0.,  0., 0.],
                    [-1.,  0., 0.],
                    [-1.,  0., -1.],
                    [-1.,  1., -1.],
                    [-2.,  1., -1.],
                    [-1.,  1., -2.]]
        self.assertArrayEqual(expected, result)


    def test_applyQuaternion04(self):
        """applyQuaternion04  -- random pts 360 about random-axis"""
        state = StateToCoords(self.random_points, tolist=1)
        q = RandomArray.uniform(self.min, self.max, (4,))
        q[3] = 360.0
        result = state.applyQuaternion(q)
        self.assertArrayEqual(self.random_points, result)


    def test_applyQuaternion05(self):
        """applyQuaternion05  -- random pts 3*120 about x-axis"""
        state = StateToCoords(self.random_points, tolist=1)
        q = (0.0, 0.0, 1.0, 120.0)
        for n in xrange(3):
            result = state.applyQuaternion(q)
        self.assertArrayEqual(self.random_points, result,0)


    def test_applyQuaternion06(self):
        """applyQuaternion06  -- random pts 2*180 about random-axis"""
        state = StateToCoords(self.random_points, tolist=1)
        q = RandomArray.uniform(self.min, self.max, (4,))
        q[3] = 180.0
        for n in xrange(2):
            result = state.applyQuaternion(q)
        self.assertArrayEqual(self.random_points, result)


    def test_applyTranslation00(self):
        """applyTranslation00 -- random pts x (0., 0., 0.)"""
        state = StateToCoords(self.random_points, tolist=1)
        zzz = Numeric.zeros((3,), 'f')
        self.assertEqual(self.random_points, state.applyTranslation(zzz))


    def test_applyTranslation01(self):
        """applyTranslation01 -- random pts x (random translation)"""
        state = StateToCoords(self.random_points, tolist=0)
        trn = RandomArray.uniform(self.min, self.max, (3,))
        expected = (Numeric.array(self.random_points) + trn)
        self.assertArrayEqual(expected, state.applyTranslation(trn))


    def test_applyTranslation02(self):
        """applyTranslation02 -- known pts x (1., 1., 1.)"""
        state = StateToCoords(self.known_points, tolist=1)
        ones = Numeric.ones((3,), 'f')
        expected = (Numeric.array(self.known_points) + ones).tolist()
        self.assertEqual(expected, state.applyTranslation(ones))


    def test_applyTranslation03(self):
        """applyTranslation03 -- random pts x (random there and back)"""
        state = StateToCoords(self.random_points, tolist=1)
        trn = RandomArray.uniform(self.min, self.max, (3,))
        state.applyTranslation(trn)
        trn = -1*trn
        self.assertArrayEqual(self.random_points, state.applyTranslation(trn))


if __name__ == '__main__':
    unittest.main()   

# for example:
#     py mglutil/math/statetest.py -v
# or, to redirect output to a file:
#     py statetest.py -v > & ! /tmp/st.out






