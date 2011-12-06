"""
    Quatfit routines for PDB2PQR

    This module is used to find the coordinates of a new
    atom based on a reference set of
    coordinates and a definition set of coordinates.

    Original Code by David J. Heisterberg
    The Ohio Supercomputer Center
    1224 Kinnear Rd.
    Columbus, OH  43212-1163
    (614)292-6036
    djh@osc.edu    djh@ohstpy.bitnet    ohstpy::djh
    
    Translated to C from fitest.f program and interfaced with
    Xmol program by Jan Labanowski,  jkl@osc.edu   jkl@ohstpy.bitnet
    ohstpy::jkl

    ----------------------------
   
    PDB2PQR -- An automated pipeline for the setup, execution, and analysis of
    Poisson-Boltzmann electrostatics calculations

    Copyright (c) 2002-2011, Jens Erik Nielsen, University College Dublin; 
    Nathan A. Baker, Battelle Memorial Institute, Developed at the Pacific 
    Northwest National Laboratory, operated by Battelle Memorial Institute, 
    Pacific Northwest Division for the U.S. Department Energy.; 
    Paul Czodrowski & Gerhard Klebe, University of Marburg.

	All rights reserved.

	Redistribution and use in source and binary forms, with or without modification, 
	are permitted provided that the following conditions are met:

		* Redistributions of source code must retain the above copyright notice, 
		  this list of conditions and the following disclaimer.
		* Redistributions in binary form must reproduce the above copyright notice, 
		  this list of conditions and the following disclaimer in the documentation 
		  and/or other materials provided with the distribution.
        * Neither the names of University College Dublin, Battelle Memorial Institute,
          Pacific Northwest National Laboratory, US Department of Energy, or University
          of Marburg nor the names of its contributors may be used to endorse or promote
          products derived from this software without specific prior written permission.

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
	ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
	WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
	IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
	INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
	BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
	DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
	LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
	OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED 
	OF THE POSSIBILITY OF SUCH DAMAGE.

    ----------------------------

"""

__date__ = "28 February 2006"
__author__ = "David Heisterberg, Jan Labanowski, Jens Erik Nielsen, Todd Dolinsky"

import math
from utilities import *

def findCoordinates(numpoints, refcoords, defcoords, defatomcoords):
    """
        Driver for the quaternion file.  Provide the coordinates as inputs
        and obtain the coordinates for the new atom as output.

        Parameters
            numpoints:     The number of points in each list (int)
            refcoords:     The reference coordinates, a list of lists of form
                           [x,y,z] (list)
            defcoords:     The definition coordinates, a list of lists of form
                           [x,y,z] (list)
            defatomcoords: The definition coordinates for the atom to be
                           placed in the reference frame (list)
        Returns
            newcoords:     The coordinates of the new atom in the
                           reference frame (list)
    """
    refcenter, fitcenter, rotation = qfit(numpoints, refcoords, defcoords)
    newcoords = qtransform(1, defatomcoords, refcenter, fitcenter, rotation)

    # Only return the first coordinates
    return newcoords[0]
    
def qtransform(numpoints, defcoords, refcenter, fitcenter, rotation):
    """
        Transform the set of defcoords using the reference center, the fit
        center, and a rotation matrix.

        Parameters
            numpoints: The number of points in each list (int)
            defcoords: Definition coordinates (list)
            refcenter: The reference center (list)
            defcenter: The definition center (list)
            rotation:  The rotation matrix (list)
        Returns
            newcoords: The coordinates of the new point (list)
    """
    
    if numpoints == 1:
        defcoords = [defcoords]

    fitcoords = translate(numpoints, defcoords, fitcenter, 1)
    rotated = rotmol(numpoints, fitcoords, rotation)
    newcoords = translate(numpoints, rotated, refcenter, 2)

    return newcoords


def qfit(numpoints, refcoords, defcoords):
    """
        Method for getting new atom coordinates from sets of reference
        and definition coordinates.

        Parameters
            numpoints: The number of points in each list (int)
            refcoords: List of reference coordinates, with each set
                       a list of form [x,y,z] (list)
            defcoords: List of definition coordinates, with each set
                       a list of form [x,y,z] (list)
    """
    nrot = 30
    
    refcenter, refcoords = center(numpoints, refcoords)
    defcenter, defcoords = center(numpoints, defcoords)

    q, u = qtrfit(numpoints, defcoords, refcoords, nrot)
    rotated = rotmol(numpoints, defcoords, u)    
    newcoords = translate(numpoints, rotated, refcenter, 2)

    return refcenter, defcenter, u

def qchichange(initcoords, refcoords, angle):
    """
        Change the chiangle of the reference coordinate using the
        initcoords and the given angle

        Parameters
            initcoords: Coordinates based on the point and basis atoms
                        (one dimensional list)
            difchi    : The angle to use (float)
            refcoords : The atoms to analyze (list of many coordinates)
        Returns
            newcoords : The new coordinates of the atoms (list of many coords)
    """
    # Initialize

    L,R = [],[]
    for i in range(3):
        L.append(0.0)
        R.append([0.0,0.0,0.0])

    # Convert to radians and normalize
    
    radangle = math.pi * angle/180.0
    normalized = normalize(initcoords)

    L[0] = normalized[0]
    L[1] = normalized[1]
    L[2] = normalized[2]
 
    # Construct the rotation matrix

    R[0][0] = math.cos(radangle) + L[0]*L[0] * (1.0 - math.cos(radangle))
    R[1][1] = math.cos(radangle) + L[1]*L[1] * (1.0 - math.cos(radangle))
    R[2][2] = math.cos(radangle) + L[2]*L[2] * (1.0 - math.cos(radangle))
    R[1][0] = L[0]*L[1]*(1.0 - math.cos(radangle)) - L[2] * math.sin(radangle)
    R[2][0] = L[0]*L[2]*(1.0 - math.cos(radangle)) + L[1] * math.sin(radangle)
    R[0][1] = L[1]*L[0]*(1.0 - math.cos(radangle)) + L[2] * math.sin(radangle)
    R[2][1] = L[1]*L[2]*(1.0 - math.cos(radangle)) - L[0] * math.sin(radangle)
    R[0][2] = L[2]*L[0]*(1.0 - math.cos(radangle)) - L[1] * math.sin(radangle)
    R[1][2] = L[2]*L[1]*(1.0 - math.cos(radangle)) + L[0] * math.sin(radangle)
    

    numpoints = len(refcoords)
    newcoords = rotmol(numpoints, refcoords, R)

    return newcoords

def rotmol(numpoints, x, u):
    """
        Rotate a molecule

        Parameters
            numpoints:  The number of points in the list (int)
            x:          The input coordinates (list)
            u:          The left rotation matrix (list)
        Returns
            out:        The rotated coordinates out=u * x (list)
    """
    out = []
    for i in range(numpoints):
        out.append([])
        out[i].append(u[0][0] *x[i][0] + u[1][0] * x[i][1] + u[2][0] * x[i][2])
        out[i].append(u[0][1] *x[i][0] + u[1][1] * x[i][1] + u[2][1] * x[i][2])
        out[i].append(u[0][2] *x[i][0] + u[1][2] * x[i][1] + u[2][2] * x[i][2])

    return out

def qtrfit(numpoints, defcoords, refcoords, nrot):
    """
        Find the quaternion, q, [and left rotation matrix, u] that minimizes
            | qTXq - Y | ^ 2 [|uX - Y| ^ 2]
        This is equivalent to maximizing Re (qTXTqY)
        The left rotation matrix, u, is obtained from q by
            u = qT1q

        Parameters
            numpoints: The number of points in each list (int)
            defcoords: List of definition coordinates, with each set
                       a list of form [x,y,z] (list)
            refcoords: List of fitted coordinates, with each set
                       a list of form [x,y,z] (list)
            nrot     : The maximum number of jacobi sweeps
        Returns
            q        : The best-fit quaternion
            u        : The best-fit left rotation matrix
    """
    xxyx = 0.0
    xxyy = 0.0
    xxyz = 0.0
    xyyx = 0.0
    xyyy = 0.0
    xyyz = 0.0
    xzyx = 0.0
    xzyy = 0.0
    xzyz = 0.0

    q = []
    c = []

    for i in range(numpoints):
         xxyx = xxyx + defcoords[i][0] * refcoords[i][0]
         xxyy = xxyy + defcoords[i][0] * refcoords[i][1]
         xxyz = xxyz + defcoords[i][0] * refcoords[i][2]
         xyyx = xyyx + defcoords[i][1] * refcoords[i][0]
         xyyy = xyyy + defcoords[i][1] * refcoords[i][1]
         xyyz = xyyz + defcoords[i][1] * refcoords[i][2]
         xzyx = xzyx + defcoords[i][2] * refcoords[i][0]
         xzyy = xzyy + defcoords[i][2] * refcoords[i][1]
         xzyz = xzyz + defcoords[i][2] * refcoords[i][2]

    for i in range(4):
        c.append([])
        for j in range(4): 
            c[i].append(0.0)

    c[0][0] = xxyx + xyyy + xzyz
    
    c[0][1] = xzyy - xyyz
    c[1][1] = xxyx - xyyy - xzyz
    
    c[0][2] = xxyz - xzyx
    c[1][2] = xxyy + xyyx
    c[2][2] = xyyy - xzyz - xxyx
    
    c[0][3] = xyyx - xxyy
    c[1][3] = xzyx + xxyz
    c[2][3] = xyyz + xzyy
    c[3][3] = xzyz - xxyx - xyyy

    d,v = jacobi(c, nrot) # diagonalize c

    for i in range(4):
        q.append(v[i][3])

    u = q2mat(q)

    return q,u

def jacobi(a, nrot):
    """
        Jacobi diagonalizer with sorted output, only good for 4x4 matrices

        Parameters
            a:    Matrix to diagonalize (4x4 list)
            nrot: Maximum number of sweeps
        Returns
            d:    Eigenvalues
            v:    Eigenvectors
    """
    v = []
    d = []
    for j in range(4):
        d.append(0)
        v.append([])
        for i in range(4):
            v[j].append(0.0)
        v[j][j] = 1.0
        d[j] = a[j][j]

    for l in range(nrot):
        dnorm = 0.0
        onorm = 0.0
        for j in range(4):
            dnorm = dnorm + abs(d[j])
            for i in range(j):
                onorm = onorm + abs(a[i][j])

        if dnorm != 0:
            if onorm/dnorm <= 1e-12: break
            
        for j in range(1,4):
            for i in range(j):
                b = a[i][j]
                if abs(b) > 0.0:
                    dma = d[j] - d[i]
                    if abs(dma) + abs(b) <= abs(dma):
                        t = b / dma
                    else:
                        q = 0.5 * dma/b
                        t = 1.0/(abs(q) + math.sqrt(1 + q*q))
                        if q < 0:
                            t = t * -1
                    c = 1.0/math.sqrt(t*t + 1)
                    s = t*c
                    a[i][j] = 0.0
                    for k in range(i):
                        atemp = c * a[k][i] - s * a[k][j]
                        a[k][j] = s * a[k][i] + c * a[k][j]
                        a[k][i] = atemp
                    for k in range(i+1 ,j):
                        atemp = c * a[i][k] - s * a[k][j]
                        a[k][j] = s * a[i][k] + c * a[k][j]
                        a[i][k] = atemp
                    for k in range(j+1, 4):
                        atemp = c * a[i][k] - s * a[j][k]
                        a[j][k] = s * a[i][k] + c * a[j][k]
                        a[i][k] = atemp
                    for k in range(4):
                        vtemp = c * v[k][i] - s * v[k][j]
                        v[k][j] = s * v[k][i] + c * v[k][j]
                        v[k][i] = vtemp
                            
                    dtemp = c*c*d[i] + s*s*d[j] - 2.0*c*s*b
                    d[j] = s*s*d[i] + c*c*d[j] +  2.0*c*s*b
                    d[i] = dtemp

    nrot = l
    for j in range(3):
        k = j
        dtemp = d[k]
        for i in range(j+1,4):
            if d[i] < dtemp:
                k = i
                dtemp = d[k]
        if k > j:
            d[k] = d[j]
            d[j] = dtemp
            for i in range(4):
                dtemp = v[i][k]
                v[i][k] = v[i][j]
                v[i][j] = dtemp
                        
    return d,v

def q2mat(q):
    """
        Generate a left rotation matrix from a normalized quaternion

        Parameters
            q:  The normalized quaternion (list)
        Returns
            u:  The rotation matrix (2-dimensional list)
    """
    u = []
    for i in range(3):
        u.append([])
        for j in range(3):
            u[i].append(0.0)

    u[0][0] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3]
    u[0][1] = 2.0 * (q[1] * q[2] - q[0] * q[3])
    u[0][2] = 2.0 * (q[1] * q[3] + q[0] * q[2])
    
    u[1][0] = 2.0 * (q[2] * q[1] + q[0] * q[3])
    u[1][1] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3]
    u[1][2] = 2.0 * (q[2] * q[3] - q[0] * q[1])
    
    u[2][0] = 2.0 *(q[3] * q[1] - q[0] * q[2])
    u[2][1] = 2.0 * (q[3] * q[2] + q[0] * q[1])
    u[2][2] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3]

    return u
    
def center(numpoints, refcoords):
    """
        Center a molecule using equally weighted points

        Parameters
            numpoints: Number of points
            refcoords: List of reference coordinates, with each set
                       a list of form [x,y,z] (list)
        Returns
            refcenter: Center of the set of points (list)
            relcoords: Moved refcoords relative to refcenter (list)
    """
    refcenter = []
    relcoords = []
    
    for i in range(3):
        refcenter.append(0.0)

    for i in range(numpoints):
        refcenter[0] += refcoords[i][0]
        refcenter[1] += refcoords[i][1]
        refcenter[2] += refcoords[i][2]

    for i in range(3):
        refcenter[i] = refcenter[i] / numpoints
        
    for i in range(numpoints):
        relcoords.append([])
        relcoords[i].append(refcoords[i][0] - refcenter[0])
        relcoords[i].append(refcoords[i][1] - refcenter[1])
        relcoords[i].append(refcoords[i][2] - refcenter[2])
        
    return refcenter, relcoords


def translate(numpoints, refcoords, center, mode):
    """
        Translate a molecule using equally weighted points

        Parameters
            numpoints: Number of points
            refcoords: List of reference coordinates, with each set
                       a list of form [x,y,z] (list)
            center:    Center of the system(list)
            mode:      If 1, center will be subtracted from refcoords
                       If 2, center will be added to refcoords
        Returns
            relcoords: Moved refcoords relative to refcenter (list)
    """
    relcoords = []
    
    if mode == 1:
        modif = -1
    elif mode == 2:
        modif = 1

    for i in range(numpoints):
        relcoords.append([])
        relcoords[i].append(refcoords[i][0] + modif * center[0])
        relcoords[i].append(refcoords[i][1] + modif * center[1])
        relcoords[i].append(refcoords[i][2] + modif * center[2])

    return relcoords
