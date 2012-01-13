#############################################################################
#
# Author: Michel F. SANNER
# Reimplemented from Babel v1.6 from Pat Walters and Math Stahl
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

"""
This file defines constants and a few functions used throughout the package


reimplmentation of Babel1.6 in Python by Michel Sanner April 2000
Original code by W. Patrick Walters and Matthew T. Stahl 
"""

#
# $Header: /opt/cvs/python/packages/share1.5/PyBabel/util.py,v 1.3 2010/10/29 17:05:26 sargis Exp $
#
# $Id: util.py,v 1.3 2010/10/29 17:05:26 sargis Exp $
#


import math

RAD_TO_DEG = 180.0/math.pi


def distance(a, b):
    """float <- distance(a,b) returns the distance between 3D points a and b"""

    dx = b[0]-a[0]
    dy = b[1]-a[1]
    dz = b[2]-a[2]
    return math.sqrt( dx*dx + dy*dy +dz*dz)


def bond_angle(a, b, c):
    """
    float <- bond_angle(a, b, c)
    returns the angle in degrees between 3D pts a,b,c
    """
    
    dist = distance(a,b) * distance(b,c)
    if dist == 0: # SD 2010 - http://mgl.scripps.edu/forum/viewtopic.php?f=11&t=245&p=1882
        raise ZeroDivisionError("Input used:", a, b, c)

    cos_theta = ( (a[0] - b[0]) * (c[0] - b[0]) +
                  (a[1] - b[1]) * (c[1] - b[1]) +
                  (a[2] - b[2]) * (c[2] - b[2]) ) / dist
    if cos_theta  + 1.0 < 0.0001: angle = 180.0
    else: angle = (math.acos(cos_theta)) * RAD_TO_DEG
    return angle


def torsion_angle(c1, c2, c3, c4):
    """
    float <- torsion_angle(a, b, c, d)
    returns the torsion angle in degrees between 3D pts a,b,c,d
    """

    v1 = (c1[0]-c2[0], c1[1]-c2[1], c1[2]-c2[2])
    v2 = (c2[0]-c3[0], c2[1]-c3[1], c2[2]-c3[2])
    v3 = (c3[0]-c4[0], c3[1]-c4[1], c3[2]-c4[2])

    p = (v2[1]*v1[2] - v1[1]*v2[2],
         v1[0]*v2[2] - v2[0]*v1[2],
         v2[0]*v1[1] - v1[0]*v2[1])

    q = (v3[1]*v2[2] - v2[1]*v3[2],
         v2[0]*v3[2] - v3[0]*v2[2],
         v3[0]*v2[1] - v2[0]*v3[1])

    n = 1.0 / math.sqrt( p[0]*p[0] + p[1]*p[1] + p[2]*p[2] )
    p = (p[0]*n, p[1]*n, p[2]*n )
    
    n = 1.0 / math.sqrt( q[0]*q[0] + q[1]*q[1] + q[2]*q[2] )
    q = (q[0]*n, q[1]*n, q[2]*n )

    xtheta = p[0]*q[0] + p[1]*q[1] + p[2]*q[2]
  
    if xtheta > 1.0: xtheta = 1.0
    if xtheta < -1.0: xtheta = -1.0
    theta = math.acos(xtheta) * 57.29578
    absth = math.fabs(theta)

    if absth < 0.001:
        return 0.0
    elif math.fabs(absth - 180.0) < 0.001:
        return 180.0

    s = v1[0]*q[0] + v1[1]*q[1] + v1[2]*q[2]
    if s < 0.0:
        theta = 360.0 - theta
  
    if theta > 180.0:
        theta = theta - 360.0

    return theta


def vec3(a, b, norm=1.0):
    """
    x,y,z <- vec3(a, b, norm=1.0)
    returns the vector a, b scale to norm
    """
    dx = b[0]-a[0]
    dy = b[1]-a[1]
    dz = b[2]-a[2]
    l = norm / math.sqrt( dx*dx + dy*dy +dz*dz)
    return [dx*l, dy*l, dz*l] 


def determinant_3x3(m):
    """
    float <- determinant_3x3(m)
    returns the determinant of the 3x3 matrix m
    """
    x = m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2])
    y = m[1][0] * (m[2][1] * m[0][2] - m[0][1] * m[2][2])
    z = m[2][0] * (m[0][1] * m[1][2] - m[1][1] * m[0][2])
    return (x + y + z)


def invert_3x3(m):
    """
    matrix3x3 <- invert_3x3(m)
    returns the inverse of a 3x3 matrix
    """

    det = determinant_3x3(m)
    if (det != 0.0):
        t = [ [0,0,0], [0,0,0], [0,0,0] ]
        det = 1.0/det
      
        t[0][0] = m[1][1]*m[2][2] - m[2][1]*m[1][2]
        t[0][1] = m[2][1]*m[0][2] - m[0][1]*m[2][2]
        t[0][2] = m[0][1]*m[1][2] - m[1][1]*m[0][2]
        t[1][0] = m[1][2]*m[2][0] - m[2][2]*m[1][0]
        t[1][1] = m[2][2]*m[0][0] - m[0][2]*m[2][0]
        t[1][2] = m[0][2]*m[1][0] - m[1][2]*m[0][0]
        t[2][0] = m[1][0]*m[2][1] - m[2][0]*m[1][1]
        t[2][1] = m[2][0]*m[0][1] - m[0][0]*m[2][1]
        t[2][2] = m[0][0]*m[1][1] - m[1][0]*m[0][1]
        
        m[0][0] = t[0][0]*det
        m[0][1] = t[0][1]*det
        m[0][2] = t[0][2]*det
        m[1][0] = t[1][0]*det
        m[1][1] = t[1][1]*det
        m[1][2] = t[1][2]*det
        m[2][0] = t[2][0]*det
        m[2][1] = t[2][1]*det
        m[2][2] = t[2][2]*det
        return m
