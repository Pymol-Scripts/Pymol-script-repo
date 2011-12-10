## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#
#
#
#$Id: test_interpolate.py,v 1.4 2007/07/24 17:30:40 vareille Exp $

from mglutil.math.rotax import interpolate3DTransform, rotax
from math import pi, sin, cos, sqrt
import numpy.oldnumeric as N
import unittest

degtorad = pi/180.
                    
class Interpolate3DBaseTest(unittest.TestCase):
    def diff(self,res, expect):
        return res-expect < 1.0e-6  # close enough -> true
        
    def test_interpolate3D(self):
        mat1=rotax([0,0,0], [0,0,1],30.0*degtorad)
        mat2=rotax([0,0,0], [0,0,1],60.0*degtorad)
        mat3=rotax([0,0,0], [0,0,1],90.0*degtorad)
        # add translation (0,1,0) to mat2 
        mat2 = N.array([
               [ 0.5       ,  0.86602539,  0.        ,  0.        ],
               [-0.86602539,  0.5       ,  0.        ,  0.        ],
               [ 0.        ,  0.        ,  1.        ,  0.        ],
               [ 0.        ,  1.        ,  0.        ,  1.        ]],'f')
    
        matList=[mat1, mat2, mat3]
        indexList = [0.33333333, 0.66666666667, 1.0]
        data = [[0.,0.,0.,1.],[1.0, 0.0, 0.0,1.0],[2.,0.,0.,1.]]
        p=0.5
        M = interpolate3DTransform(matList, indexList, p)
        
        res=N.dot(data, M)[1]
        
        self.assertEqual( self.diff(res[0], 0.70710677 ),True)
        self.assertEqual(self.diff(res[1], 1.20710677 ),True) # 50% translation along Y axis
        self.assertEqual(self.diff(res[2], 0.0),True)
    
        p=1.5
        M = interpolate3DTransform(matList, indexList, p)
        res=N.dot(data, M)[1]
        self.assertEqual(self.diff(res[0], -0.70710677 ),True)
        self.assertEqual(self.diff(res[1],  0.70710677 ),True)
        self.assertEqual(self.diff(res[2],  0.0),True)


if __name__ == '__main__':
    unittest.main()
