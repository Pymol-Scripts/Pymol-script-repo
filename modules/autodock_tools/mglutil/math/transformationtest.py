## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#
# Last modified on Wed Aug 22 14:22:48 PDT 2001 by lindy
#
# $Header: /opt/cvs/python/packages/share1.5/mglutil/math/transformationtest.py,v 1.3 2007/07/24 17:30:40 vareille Exp $
#

"""Unit test for transformation.py

Transformation.py defines three classes:
    Quaternion
    UnitQuaternion(Quaternion)
    Transformation(UnitQuaternion)

"""

from mglutil.math.transformation import UnitQuaternion, Quaternion
from mglutil.math.transformation import Transformation
import unittest
#import numpy.oldnumeric as N
import numpy.oldnumeric as Numeric
N = Numeric
import numpy.oldnumeric.random_array as RA
import math

#
# Unit tests for the Quaternion class
#

class QuaternionTest(unittest.TestCase):
    def setUp(self):
        """The Quaternion class is tested through the UnitQuaternion class"""
        pass



class UnitQuaternionTest(unittest.TestCase):
    def setUp(self):
        self.max = 999999999.
        self.min = -self.max

        

class UnitQuaternionKnownValues(UnitQuaternionTest):
##          theta = 360.0*RA.random()
##          knownValues = (
##              # identity
##              ((0., 0., 0., 0.), ( (1.0, 0.0, 0.0, 0.0),
##                                   (0.0, 1.0, 0.0, 0.0),
##                                   (0.0, 0.0, 1.0, 0.0),
##                                   (0.0, 0.0, 0.0, 1.0))),
##              # rotation about x
##              ((0., 0., 0., 0.), ( (1.0, 0.0, 0.0, 0.0),
##                                   (0.0, 1.0, 0.0, 0.0),
##                                   (0.0, 0.0, 1.0, 0.0),
##                                   (0.0, 0.0, 0.0, 1.0))),
##              # rotation about y
##              ((0., 0., 0., 0.), ( (1.0, 0.0, 0.0, 0.0),
##                                   (0.0, 1.0, 0.0, 0.0),
##                                   (0.0, 0.0, 1.0, 0.0),
##                                   (0.0, 0.0, 0.0, 1.0))),
##              # rotation about z
##              ((0., 0., 0., 0.), ( (1.0, 0.0, 0.0, 0.0),
##                                   (0.0, 1.0, 0.0, 0.0),
##                                   (0.0, 0.0, 1.0, 0.0),
##                                   (0.0, 0.0, 0.0, 1.0))),
##              ((0., 0., 0., 0.), ( (1.0, 0.0, 0.0, 0.0),
##                                   (0.0, 1.0, 0.0, 0.0),
##                                   (0.0, 0.0, 1.0, 0.0),
##                                   (0.0, 0.0, 0.0, 1.0))),
##              ((0., 0., 0., 0.), ( (1.0, 0.0, 0.0, 0.0),
##                                   (0.0, 1.0, 0.0, 0.0),
##                                   (0.0, 0.0, 1.0, 0.0),
##                                   (0.0, 0.0, 0.0, 1.0))),
##              ((0., 0., 0., 0.), ( (1.0, 0.0, 0.0, 0.0),
##                                   (0.0, 1.0, 0.0, 0.0),
##                                   (0.0, 0.0, 1.0, 0.0),
##                                   (0.0, 0.0, 0.0, 1.0))),
##              ((0., 0., 0., 0.), ( (1.0, 0.0, 0.0, 0.0),
##                                   (0.0, 1.0, 0.0, 0.0),
##                                   (0.0, 0.0, 1.0, 0.0),
##                                   (0.0, 0.0, 0.0, 1.0))),
##              ((0., 0., 0., 0.), ( (1.0, 0.0, 0.0, 0.0),
##                                   (0.0, 1.0, 0.0, 0.0),
##                                   (0.0, 0.0, 1.0, 0.0),
##                                   (0.0, 0.0, 0.0, 1.0))))

    def tearDown(self):
        pass


    def testKnown00(self):
        """<describe what's being tested here>"""
        pass



class UnitQuaternionProperties(UnitQuaternionTest):

    def testProperties00(self):
        """The product of the conjugate is the conjucate of the product"""
        q1 = Quaternion(tuple(RA.uniform(self.min,
                                         self.max, (1, 4)).tolist()[0]))
        q2 = Quaternion(tuple(RA.uniform(self.min,
                                         self.max, (1, 4)).tolist()[0]))
##          pc = q1.conjugate()*q2.conjugate()
##          qp = q1*q2
##          cp = qp.conjugate()
##          self.assertEqual( pc, cp)
        # the commented code fails with the same error as this line...
        self.assertEqual( q1.conjugate()*q2.conjugate(), (q2*q1).conjugate())
        

    def testProperties01(self):
        """The magnitudue of the product is the product of the magnitudes"""
        q1 = Quaternion(tuple(RA.uniform(self.min,
                                         self.max, (1, 4)).tolist()[0]))
        q2 = Quaternion(tuple(RA.uniform(self.min,
                                         self.max, (1, 4)).tolist()[0]))
        self.assertEqual((q1*q2).magnitude(), q1.magnitude()*q2.magnitude())


    def testProperties02(self):
        """The conjugate of a unit quaternion is it's inverse"""
        q1 = Quaternion(tuple(RA.uniform(self.min,
                                         self.max, (1, 4)).tolist()[0]))
        self.assertEqual( q1.conjugate(), q1.inverse())

    


class TransformationTest(unittest.TestCase):
    def setUp(self):
        pass


    def tearDown(self):
        pass

    def test(self):
        pass



if __name__ == '__main__':
    unittest.main()   

# for example: py mglutil/math/transformationtest.py -v


