## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

from mglutil.regression import testplus
from FlexTree.FTMotions import FTMotion_RotationAboutAxis, FTMotion_Hinge
import numpy.oldnumeric as Numeric

def diff(res, expect):
    return res-expect < 1.0e-6  # close enough -> true


def rotateObject():
    from mglutil.math.rotax import rotax
    from math import sin, cos, pi, sqrt, fabs
    import numpy.oldnumeric as N
    import random

    degtorad = pi/180.
    point1 = [random.random()*8-4, random.random()*8-4.,random.random()*10-5]
    point2 = [random.random()*8-4, random.random()*8-4.,random.random()*10-5.]
    angle  = random.random()*360.
    transf = rotax(point1, point2, angle*degtorad, transpose=1)

    from mglutil.math.rotax import mat_to_axis_angle
    vector , pp, angle = mat_to_axis_angle(transf)

    from FlexTree.FTMotions import FTMotion_Hinge
    motion = FTMotion_Hinge(axis={'vector':vector,'point':pp})
    motion.configure(angle=angle)

    m1=Numeric.array(motion.getMatrix()).ravel()
    m2=transf.ravel()
    bSame = True
    for i in range(len(m1)):
        if fabs(m1[i]-m2[i]) > 1e-4:
            bSame = False
    assert (bSame, True)


def test_rotateObject():
    for i in range(10):
        rotateObject()
    

harness = testplus.TestHarness( __name__,
                              funs = testplus.testcollect( globals()),
                              )

if __name__ == '__main__':
    print harness
    sys.exit( len( harness))
