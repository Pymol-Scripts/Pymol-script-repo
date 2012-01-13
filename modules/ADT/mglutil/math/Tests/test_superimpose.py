## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

from mglutil.regression import testplus
from mglutil.math.rotax import interpolate3DTransform, rotax
#from math import pi, sin, cos, sqrt
import numpy.oldnumeric as N
#degtorad = pi/180.
                    
def diff(res, expect):
    return res-expect < 1.0e-6  # close enough -> true




def test_superimpose():
    from mglutil.math.rigidFit import RigidfitBodyAligner
    rigidfitAligner = RigidfitBodyAligner()
    refCoords=[[0,0,0] , [1,0,0], [0,1,0], [0,0,1]]
    mobCoords=[[10,0,0] , [11,0,0], [10,1,0], [10,0,1]]
    rigidfitAligner.setRefCoords(refCoords)                
    rigidfitAligner.rigidFit(mobCoords)
    rmsd=rigidfitAligner.rmsdAfterSuperimposition(mobCoords)
    assert diff(rmsd, 0.0 )
