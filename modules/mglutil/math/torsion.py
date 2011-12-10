## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#taken from Pmv/measureCommands.py

def torsion( x1, x2, x3, x4):
    """
    Compute the torsion angle between x1, x2, x3, x4.
    All coordinates are cartesian; result is in degrees.
    Raises a ValueError if angle is not defined.
    """
    from math import sqrt, acos
    import numpy.oldnumeric as Numeric
    N=Numeric
    
    tang=0.0
    x1 = N.array(x1, 'f')
    x2 = N.array(x2, 'f')
    x3 = N.array(x3, 'f')
    x4 = N.array(x4, 'f')
    
    assert x1.shape == (3, )
    assert x2.shape == (3, )
    assert x3.shape == (3, )
    assert x4.shape == (3, )

    a = x1-x2
    b = x3-x2
    c = vvmult(a, b)

    a = x2-x3
    b = x4-x3
    d = vvmult(a, b)

    dd=sqrt(Numeric.sum(c*c))
    de=sqrt(Numeric.sum(d*d))

    if dd<0.001 or de<0.001:
        raise ValueError ( 'Torsion angle undefined, degenerate points')

    vv = Numeric.dot(c, d) / (dd*de);
    if vv<1.0: tang=vv
    else: tang= 1.0
    if tang<-1.0: tang=-1.0
    tang = acos(tang)
    tang = tang*57.296

    b = vvmult(c, d)
    if Numeric.dot(a, b) > 0.0: tang = -tang
    return tang

def vvmult( a, b):
    """
    Compute a vector product for 3D vectors
    """
    import numpy.oldnumeric as Numeric
    res = Numeric.zeros(3, 'f')
    res[0] = a[1]*b[2] - a[2]*b[1]
    res[1] = a[2]*b[0] - a[0]*b[2]
    res[2] = a[0]*b[1] - a[1]*b[0]
    return res
