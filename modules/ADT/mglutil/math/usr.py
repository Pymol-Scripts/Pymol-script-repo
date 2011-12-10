## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#!/usr/bin/env python
#
# Last modified on Wed Mar 21 13:32:05 PDT 2007 by lindy
#
# $Id: usr.py,v 1.3 2007/07/24 17:30:40 vareille Exp $
#

"""usr.py - ultrafast shape recognition

Implements method described in:
  Ballester, PJ & Richards, WG (2007) Proc.R.Soc.A
  doi:10.1098/rspa.2007.1823
"""

import numpy.oldnumeric as N
import math

    
def centroid(points):
    """reuturn the centroid of points
    """
    points = N.array(points)
    num, dim = points.shape
    return N.add.reduce(points)/float(num)


def dist_array(point, points):
    """return array of distances from point to each of points
    """
    point  = N.array(point)
    points = N.array(points)
    return N.sqrt(N.add.reduce((points - point)**2, axis=1))


def mean_var_skew(data):
    """return the mean, variance, and skewness of data
    """
    data = N.array(data)
    num = data.shape[0]
    mean = N.add.reduce(data)/float(num)
    var = N.add.reduce((data - mean)**2)/float(num-1)
    #
    std = math.sqrt(var)
    skew = N.add.reduce( (data - mean)**3)/(float(num-1)* std**3)
    return mean, var, skew

        
def usr_descriptors(points):
    """return 12-tuple of geoemtric descriptors for points

Reference for method:
Ballester, PJ & Richards, WG (2007) Proc.R.Soc.A
doi:10.1098/rspa.2007.1823
    """
    # centroid
    ctr = centroid(points)
    ctr_da = dist_array(ctr, points)
    ctr_m, ctr_v, ctr_s = mean_var_skew(ctr_da)

    # closest to centroid
    cst = points[ N.argmin(ctr_da) ]
    cst_da = dist_array(cst, points)
    cst_m, cst_v, cst_s = mean_var_skew(cst_da)
    
    # farthest from centroid
    fct = points [ N.argmax(ctr_da) ]
    fct_da = dist_array(fct, points)
    fct_m, fct_v, fct_s = mean_var_skew(fct_da)
    
    # farthest from fct
    ftf = points [ N.argmax(fct_da) ]
    ftf_da = dist_array(ftf, points)
    ftf_m, ftf_v, ftf_s = mean_var_skew(ftf_da)

    return (ctr_m, ctr_v, ctr_s,
            cst_m, cst_v, cst_s,
            fct_m, fct_v, fct_s,
            ftf_m, ftf_v, ftf_s)


def usr_similarity(x, y):
    """return similarity of two usr_descriptor vectors, x and y
    """
    x = N.array(x)
    num = float(x.shape[0])
    y = N.array(y)
    # normalized and montonically inverted Manhattan distance
    return num/(num + N.add.reduce(N.absolute(x-y)))


def print_usr_desc(val):
    print "%0.3f, %0.3f, %0.3f, %0.3f, %0.3f, %0.3f, %0.3f, %0.3f, %0.3f, %0.3f, %0.3f, %0.3f" %val


#
# everything down here is just scaffold and
# technically not needed for the USR code above.
#

import unittest


def distance_matrix(points):
    """return NxN point to point distance matrix
    this works but is not needed for USR
    """
    points = N.array(points)
    num, dim = points.shape
    delta = N.zeros((num,num), 'd')
    for d in xrange(dim):
        data = points[:,d]
        delta += (data - data[:,N.NewAxis])**2
    return N.sqrt(delta)


def closest( point, points):
    """works, but @@DEPRECATED!!
    return index into points of closest-to-point
    """
    da = dist_array(point, points)
    return N.argmin(da)

def farthest( point, points):
    """works, but @@DEPRECATED!!
    return index into points of farthest-from-point
    """
    da = dist_array(point, points)
    return N.argmax(da)


class USRTest(unittest.TestCase):
    def setUp(self):
        self.points_2D = [[0, 0], [1, 1], [4, 5]]
        self.points_3D = [[0, 0, 0], [1.0, 1, 1], [4, 5, 6], [10,10,10]]

        # from Ballester & Richards, Figure 2. 
        self.mq = (4.44, 2.98, 1.04,
                   4.55, 4.70, 0.23,
                   8.30, 16.69, -22.97,
                   7.37, 15.64, 0.51)
        self.mi = (4.39, 3.11, 1.36,
                   4.50, 4.44, 0.09,
                   8.34, 16.78, -23.20,
                   7.15, 16.52, 0.13)
        self.sqi = 0.811359026369 # but in Figure 2, Sqi = 0.812 !!

        
class DistanceMatrixTest(USRTest):
    def test_2D(self):
        dm = distance_matrix(self.points_2D)
        self.assertAlmostEqual(1.414213562373095049, dm[0][1])
        self.assertAlmostEqual(6.403124237432848686, dm[0][2])
        self.assertAlmostEqual(5, dm[1][2])
        self._testSymmetry(dm)
        
    def test_3D(self):
        dm = distance_matrix(self.points_3D)
        self.assertAlmostEqual(1.732050807568877294, dm[0][1])
        self.assertAlmostEqual(8.77496438739212206, dm[0][2])
        self.assertAlmostEqual(17.32050807568877294, dm[0][3])
        self.assertAlmostEqual(7.071067811865475244, dm[1][2])
        self.assertAlmostEqual(15.58845726811989564, dm[1][3])
        self.assertAlmostEqual(8.77496438739212206, dm[2][3])
        self._testSymmetry(dm)

    def _testSymmetry(self, dm):
        for i in range(len(dm)):
            for j in range(len(dm)):
                self.assertEqual(dm[i][j], dm[j][i])

                
class CentroidTest(USRTest):
    def test_2D(self):
        ctr = centroid(self.points_2D)
        self.assertAlmostEqual((5./3.),  ctr[0])
        self.assertAlmostEqual((6./3.),  ctr[1])


    def test_3D(self):
        ctr = centroid(self.points_3D)
        self.assertAlmostEqual((15./4.),   ctr[0])
        self.assertAlmostEqual((16./4.),   ctr[1])
        self.assertAlmostEqual((17./4.),  ctr[2])

    
class ClosestTest(USRTest):
    def test_2D(self):
        point = [0.5, 0.]
        ix = closest(point, self.points_2D)
        self.assertEqual( 0, ix)
        point = [1.5, 1.]
        ix = closest(point, self.points_2D)
        self.assertEqual( 1, ix)
        point = [10., 10.]
        ix = closest(point, self.points_2D)
        self.assertEqual( 2, ix)

    def test_3D(self):
        point = [0.5, 0., 0.]
        ix = closest(point, self.points_3D)
        self.assertEqual( 0, ix)
        point = [1.5, 1., 1.]
        ix = closest(point, self.points_3D)
        self.assertEqual( 1, ix)
        point = [4.5, 5., 6.8]
        ix = closest(point, self.points_3D)
        self.assertEqual( 2, ix)
        point = [10., 10., 11.]
        ix = closest(point, self.points_3D)
        self.assertEqual( 3, ix)


class FarthestTest(USRTest):
    def test_2D(self):
        point = [0.5, 0.]
        ix = farthest(point, self.points_2D)
        self.assertEqual( 2, ix)
        point = [1.5, 1.]
        ix = farthest(point, self.points_2D)
        self.assertEqual( 2, ix)
        point = [10., 10.]
        ix = farthest(point, self.points_2D)
        self.assertEqual( 0, ix)

    def test_3D(self):
        point = [0.5, 0., 0.]
        ix = farthest(point, self.points_3D)
        self.assertEqual( 3, ix)
        point = [1.5, 1., 1.]
        ix = farthest(point, self.points_3D)
        self.assertEqual( 3, ix)
        point = [4.5, 5., 6.8]
        ix = farthest(point, self.points_3D)
        self.assertEqual( 0, ix)
        point = [10., 10., 11.]
        ix = farthest(point, self.points_3D)
        self.assertEqual( 0, ix)


class SimilarityTest(USRTest):
    def test_fig2(self):
        s = usr_similarity(self.mq, self.mi)
        self.assertAlmostEqual( self.sqi, s)


if __name__ == '__main__':
    unittest.main()
