## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

""" This file contains the following functions,
rotax
mat_to_quat
mat_to_axis_angle
inverse4X4
rotVectToVect
interpolate3DTransform
"""


from math import pi, sin, cos, sqrt
import numpy.oldnumeric as N
degtorad = pi/180.


def rotax( a, b, tau, transpose=1 ):
    """
    Build 4x4 matrix of clockwise rotation about axis a-->b
    by angle tau (radians).
    a and b are sequences of 3 floats each
    Result is a homogenous 4x4 transformation matrix.
    NOTE: This has been changed by Brian, 8/30/01: rotax now returns
    the rotation matrix, _not_ the transpose. This is to get
    consistency across rotax, mat_to_quat and the classes in
    transformation.py
    when transpose is 1 (default) a C-style rotation matrix is returned
    i.e. to be used is the following way Mx (opposite of OpenGL style which
    is using the FORTRAN style)
    """

    assert len(a) == 3
    assert len(b) == 3
    if tau <= -2*pi or tau >= 2*pi:
        tau = tau%(2*pi)

    ct = cos(tau)
    ct1 = 1.0 - ct
    st = sin(tau)

    # Compute unit vector v in the direction of a-->b. If a-->b has length
    # zero, assume v = (1,1,1)/sqrt(3).

    v = [b[0]-a[0], b[1]-a[1], b[2]-a[2]]
    s = v[0]*v[0]+v[1]*v[1]+v[2]*v[2]
    if s > 0.0:
        s = sqrt(s)
	v = [v[0]/s, v[1]/s, v[2]/s]
    else:
        val = sqrt(1.0/3.0)
        v = (val, val, val)

    rot = N.zeros( (4,4), 'f' )
    # Compute 3x3 rotation matrix

    v2 = [v[0]*v[0], v[1]*v[1], v[2]*v[2]]
    v3 = [(1.0-v2[0])*ct, (1.0-v2[1])*ct, (1.0-v2[2])*ct]
    rot[0][0]=v2[0]+v3[0]
    rot[1][1]=v2[1]+v3[1]
    rot[2][2]=v2[2]+v3[2]
    rot[3][3] = 1.0;

    v2 = [v[0]*st, v[1]*st, v[2]*st]
    rot[1][0]=v[0]*v[1] * ct1-v2[2]
    rot[2][1]=v[1]*v[2] * ct1-v2[0]
    rot[0][2]=v[2]*v[0] * ct1-v2[1]
    rot[0][1]=v[0]*v[1] * ct1+v2[2]
    rot[1][2]=v[1]*v[2] * ct1+v2[0]
    rot[2][0]=v[2]*v[0] * ct1+v2[1]

    # add translation
    for i in (0,1,2):
        rot[3][i] = a[i]
	for j in (0,1,2):
            rot[3][i] = rot[3][i]-rot[j][i]*a[j]
	rot[i][3]=0.0

    if transpose:
        return rot
    else:
        return N.transpose(rot)


from math import asin

def rotVectToVect(vect1, vect2, i=None):
    """returns a 4x4 transformation that will align vect1 with vect2
vect1 and vect2 can be any vector (non-normalized)
"""
    v1x, v1y, v1z = vect1
    v2x, v2y, v2z = vect2
    
    # normalize input vectors
    norm = 1.0/sqrt(v1x*v1x + v1y*v1y + v1z*v1z )
    v1x *= norm
    v1y *= norm
    v1z *= norm    
    norm = 1.0/sqrt(v2x*v2x + v2y*v2y + v2z*v2z )
    v2x *= norm
    v2y *= norm
    v2z *= norm
    
    # compute cross product and rotation axis
    cx = v1y*v2z - v1z*v2y
    cy = v1z*v2x - v1x*v2z
    cz = v1x*v2y - v1y*v2x

    # normalize
    nc = sqrt(cx*cx + cy*cy + cz*cz)
    if nc==0.0:
        return [ [1., 0., 0., 0.],
                 [0., 1., 0., 0.],
                 [0., 0., 1., 0.],
                 [0., 0., 0., 1.] ]

    cx /= nc
    cy /= nc
    cz /= nc
    
    # compute angle of rotation
    if nc<0.0:
        if i is not None:
            print 'truncating nc on step:', i, nc
        nc=0.0
    elif nc>1.0:
        if i is not None:
            print 'truncating nc on step:', i, nc
        nc=1.0
        
    alpha = asin(nc)
    if (v1x*v2x + v1y*v2y + v1z*v2z) < 0.0:
        alpha = pi - alpha

    # rotate about nc by alpha
    # Compute 3x3 rotation matrix

    ct = cos(alpha)
    ct1 = 1.0 - ct
    st = sin(alpha)
    
    rot = [ [0., 0., 0., 0.],
            [0., 0., 0., 0.],
            [0., 0., 0., 0.],
            [0., 0., 0., 0.] ]


    rv2x, rv2y, rv2z = cx*cx, cy*cy, cz*cz
    rv3x, rv3y, rv3z = (1.0-rv2x)*ct, (1.0-rv2y)*ct, (1.0-rv2z)*ct
    rot[0][0] = rv2x + rv3x
    rot[1][1] = rv2y + rv3y
    rot[2][2] = rv2z + rv3z
    rot[3][3] = 1.0;

    rv4x, rv4y, rv4z = cx*st, cy*st, cz*st
    rot[0][1] = cx * cy * ct1 - rv4z
    rot[1][2] = cy * cz * ct1 - rv4x
    rot[2][0] = cz * cx * ct1 - rv4y
    rot[1][0] = cx * cy * ct1 + rv4z
    rot[2][1] = cy * cz * ct1 + rv4x
    rot[0][2] = cz * cx * ct1 + rv4y

    return rot


def mat_to_quat(matrix,transpose=1):
    """ takes a four by four matrix (optionally with shape (16,) and
    converts it into the axis of rotation and angle to rotate by
    (x,y,z,theta). It does not expect an OpenGL style, transposed
    matrix, so is consistent with rotax
    """
    
    if N.shape(matrix) not in ((16,),(4,4)):
        raise ValueError("Argument must Numeric array of shape (4,4) or (16,)")

    if N.shape(matrix) == (4, 4):
        matrix = N.reshape(matrix,(16,))

    cofactor1 = matrix[5]*matrix[10] - matrix[9]*matrix[6]
    cofactor2 = matrix[8]*matrix[6] - matrix[4]*matrix[10] 
    cofactor3 = matrix[4]*matrix[9] - matrix[8]*matrix[5]

    det = matrix[0]*cofactor1 + matrix[1]*cofactor2 + matrix[2]*cofactor3
    if not (0.999 < det < 1.001):
        print "Not a unit matrix: so not a pure rotation"
        print 'Value of Determinant is: ',det
    trace = matrix[0] + matrix[5] + matrix[10] + matrix[15]       
    if trace > 0.0000001:# rotation other than 180deg
        S = 0.5/sqrt(trace)
        Qw = 0.25/S
        Qx = (matrix[9]-matrix[6])*S
        Qy = (matrix[2]-matrix[8])*S
        Qz = (matrix[4]-matrix[1])*S
    else: #180deg rotation, just need to figure out the axis
        Qw = 0.
        diagonal = ((matrix[0],0),
                    (matrix[5],5),
                    (matrix[10],10))
        idx = max(diagonal)[1]
        if idx == 0:
            S  = sqrt(1.0 + matrix[0] - matrix[5] - matrix[10])*2
            Qy = (matrix[1] + matrix[4] ) / S
            Qz = (matrix[2] + matrix[8] ) / S
            Qx = N.sqrt(1-Qy*Qy-Qz*Qz)
        elif idx==5:
            S  = sqrt( 1.0 + matrix[5] - matrix[0] - matrix[10] )*2
            Qx = (matrix[1] + matrix[4] ) / S
            Qz = (matrix[6] + matrix[9] ) / S
            Qy = N.sqrt(1-Qx*Qx-Qz*Qz)
        elif idx==10:
            S  = sqrt( 1.0 + matrix[10] - matrix[0] - matrix[5] )*2
            Qx = (matrix[2] + matrix[8] ) / S
            Qy = (matrix[6] + matrix[9] ) / S
            Qz = N.sqrt(1-Qx*Qx-Qy*Qy)
    # check if identity or not
    if Qw != 1.:
        angle = N.arccos(Qw)
        theta = angle*360./N.pi
        Z = sqrt(Qx*Qx + Qy*Qy + Qz*Qz)
        if transpose:
            Qx = -Qx/Z
            Qy = -Qy/Z
            Qz = -Qz/Z
        else:
            Qx = Qx/Z
            Qy = Qy/Z
            Qz = Qz/Z                
        Qw = theta
        return [Qx,Qy,Qz,Qw]
    else:
        return [0.,0.,0.,0.]



def inverse4X4(matrix):
    """ returns the inverse of the given 4x4 transformation matrix
t_1: the negetive of Translation vector
r_1: the inverse of rotation matrix

inversed transformation is
1) t_1 applied first
2) then r_1 is applied

to validate the result, N.dot(matrix, mat_inverse)==N.identity(4,'f')
"""
    # check the shape
    if matrix.shape !=(4,4) and matrix.shape !=(16,) :
        raise ValueError("Argument must Numeric array of shape (4,4) or (16,)")
        return None
    if matrix.shape ==(16,):            
        matrix=N.array(matrix,'f')
        matrix=N.reshape(matrix,(4,4))  # force the matrix to be (4,4)
    t_1=N.identity(4,'f')
    t_1[:2,3]= - matrix[:2, 3]
    r_1=N.identity(4,'f')
    r_1[:3,:3] = N.transpose(matrix[:3,:3])
    mat_inverse=N.dot(r_1, t_1)
    #asert N.dot(matrix, mat_inverse) is N.identity(4,'f')
    return mat_inverse


def mat_to_axis_angle( matrix ):
    """
    NOTE: This function is added by Yong 2/01/04: given a 4x4 transformation
matrix of hinge motion, now returns the rotation angle and axis (defined by
vector and a point) Please be noticed that if the motion is not hinge, the
function will complain and return none
"""
    if matrix.shape != (16,) and matrix.shape != (4,4):
        raise ValueError("matrix should be of shape (4,4) or (16,)")
        return None

    if matrix.shape == (16,):
        matrix = N.reshape(matrix, (4,4))

    from math import sin, cos, pi, sqrt, fabs, acos
    degtorad = pi/180.

    transf = matrix
    from mglutil.math.rotax import mat_to_quat
    rotMat =  N.identity(4, 'f')
    rotMat[:3,:3] = transf[:3,:3]
    qB = mat_to_quat(matrix=N.array(rotMat).ravel())
    angle = qB[3]
    sa=sin(angle*degtorad/2.0)
    if(fabs(angle) < 0.0005):
        sa = 1
    if angle > 180.:
        vector=[-qB[0]/sa, -qB[1]/sa, -qB[2]/sa]
    else:
        vector=[qB[0]/sa, qB[1]/sa, qB[2]/sa]
    tranMat = transf[3,:3]

    # check if the transformation is a hinge motion
    a=vector
    b=tranMat
    c =[a[0]-b[0], a[1]-b[1], a[2]-b[2]]
    a2= a[0]*a[0] + a[1]*a[1] + a[2]*a[2]
    b2= b[0]*b[0] + b[1]*b[1] + b[2]*b[2]
    c2= c[0]*c[0] + c[1]*c[1] + c[2]*c[2]
    theta = acos((c2-a2-b2)/(2* sqrt(a2*b2))) / pi * 180
    if fabs(theta -90) > 1e-4:
        raise ValueError("The given transformation is not a hinge motion")
        return None, None, None
    
    ratio = sqrt( 1. / (2. * (1.- cos(degtorad * angle ))))
    p = [tranMat[0]*ratio, tranMat[1]*ratio, tranMat[2]*ratio]

    ang = 90. - angle/2.0
    rot = rotax([0.,0.,0.], vector, ang*degtorad, transpose=1)
    rot = rot[:3,:3]
    point = N.dot(p, rot)
    
    return vector, point, angle


def interpolate3DTransform(matrixList, indexList, percent):
    """ This function gets input of two list and a percent value.
Return value is a 4x4 matrix corresponding to percent% of the transformation.

matrixList: a list of 4x4 transformation matrix
indexList : a list of sorted index (positive float number)
percent   : a positive float number.
if only one matrix in the matrix list:
percent =   0.0  means no transformation (identity)
            1.0  means 100% of the transformation (returns mat)
            0.58 means 58% of translation and rotatetion 58% of rotation angle
            along the same rotation axis
percent can go above 1.0

If matrixList has more than one matrix:
matrixList=[M1,  M2,  M3]     #Attention: All M uses the same reference frame
indexList =[0.2, 0.5, 1.0]    #Attention: assume the list sorted ascendingly
p = 0.5 means apply M2
p = 0.8 means apply M3
p = 0.9 means apply M2 first, then apply 50% of M'.  M' is the transformation
                    from M2 to M3.   50% = (0.9-0.8) / (1.0-0.8)
                    M2 x M' = M3
                    -->  M2.inverse x M2 x M'= M2.inverse x M3 
                    -->  M'= M2.inverse x M
"""
    listLen = len(matrixList)
    if listLen != len(indexList):
        raise ValueError("matrix list should have same length of index list")
    if listLen == 0:
        raise ValueError("no matrix found in the matrix list")

    offset = -1
    for i in range(listLen):
        if indexList[i] >= percent:
            offset = i
            break

    prevMat = nextMat = N.identity(4,'f')
    if offset == -1:
        prevMat = matrixList[-1]
        p = percent/indexList[-1]
        return _interpolateMat(matrixList[-1], p)
    elif offset == 0:
        nextMat = matrixList[0]
        p = percent/indexList[0]
        return _interpolateMat(N.array(matrixList[0]), p)
    else:
        prevMat = matrixList[offset-1]
        nextMat = matrixList[offset]
        p = (percent-indexList[offset-1])/(
                                    indexList[offset]-indexList[offset-1])
        from numpy.oldnumeric.linear_algebra import inverse
        M = N.dot(inverse(prevMat), nextMat)
        Mat = _interpolateMat(M, p)
        return N.dot(prevMat, Mat)

    
def interpolate3DTransform1(matrixList, indexList, percent):
    # MS version that does not assume identity as fist matrix and does
    # not wrap around

    if percent <= indexList[0]:
        return matrixList[0]

    if percent >=indexList[-1]:
        return matrixList[-1]
    
    listLen = len(indexList)
    for i in range(listLen):
        if indexList[i] > percent:
            break

    prevMat = matrixList[i-1]
    nextMat = matrixList[i]
    from numpy.oldnumeric.linear_algebra import inverse
    M = N.dot(inverse(prevMat), nextMat)
    p = (percent-indexList[i-1]) / (indexList[i]-indexList[i-1])
    Mat = _interpolateMat(M, p)
    return N.dot(prevMat, Mat)



def _interpolateMat(mat, percent):
    """ called only by interpolate3DTransform()
    """
    if mat.shape != (16,) and mat.shape != (4,4):
        raise ValueError("matrix should be of shape (4,4) or (16,)")
        return None

    if mat.shape == (16,):
        mat = N.reshape(mat, (4,4))

  #  if percent <0.0:
  #      raise ValueError('The parameter percent should be a positive float"')
  #      return

    p = percent
    transf = mat[:,:]
    rotMat =  N.identity(4, 'f')
    rotMat[:3,:3]=transf.astype(N.Float32)[:3,:3]

    from mglutil.math.rotax import mat_to_quat
    quat = mat_to_quat(matrix=N.array(rotMat).ravel())
    angle = quat[3] * p

    newRotMat = rotax([0.,0.,0.], quat[:3], angle*degtorad, transpose = 1)
    newTranMat =  N.identity(4, 'f')
    newTranMat[3][0] = transf[3][0]*p
    newTranMat[3][1] = transf[3][1]*p
    newTranMat[3][2] = transf[3][2]*p

    transform = newRotMat
    transform[3][0] = newTranMat[3][0]
    transform[3][1] = newTranMat[3][1]
    transform[3][2] = newTranMat[3][2]
    # That's it.. the self.transform is now updated.    
    
    return transform
