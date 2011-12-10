## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

import numpy.oldnumeric as Numeric, types

def EulerAnglesToMat( angles ):
    """Builds a rotation matrix from Euler angles given in degrees"""

    from math import pi, cos, sin
    
    t1 = angles[0]* (pi/180.0)
    t2 = angles[1]* (pi/180.0)
    t3 = angles[2]* (pi/180.0)
    # from Lattman, Meth. Enzymology V. 115 p. 63
    rot = Numeric.identity(3).astype('d')
    rot[0][0]= -sin(t1)*cos(t2)*sin(t3) + cos(t1)*cos(t3)
    rot[0][1]= cos(t1)*cos(t2)*sin(t3) + sin(t1)*cos(t3)
    rot[0][2]= sin(t2)*sin(t3)
    rot[1][0]= -sin(t1)*cos(t2)*cos(t3) - cos(t1)*sin(t3)
    rot[1][1]= cos(t1)*cos(t2)*cos(t3) - sin(t1)*sin(t3)
    rot[1][2]= sin(t2)*cos(t3)
    rot[2][0]= sin(t1)*sin(t2)
    rot[2][1]= -cos(t1)*sin(t2)
    rot[2][2]= cos(t2)
    return rot



class Crystal:
    """Class to provide functionalities related to crystal structures"""


    def __init__(self, length, angles):
        """Constructor for Crystal"""

        self.ctof = None # transpose of Xformation to go from cart. to fract.
        self.ftoc = None # transpose of Xformation to go from fract. to cart.
        assert len(length)==3 and type(length[0])==types.FloatType
        self.length = length
        assert len(angles)==3 and type(angles[0])==types.FloatType
        self.angles = angles
        self.ctof = self._orthog( self.length, self.angles) # build ctof
        self.ftoc = self._uinv( Numeric.transpose(self.ctof) ) # build ftoc


    def _orthog(self, length, angles):
        """
           let u be a 3x3 transformation matrix which relates a new cell
           with orthogonal axes of unit dimensions to the original unit
           cell (crystal system).  
            Aug 5, 1991 changed to XPLOR convention was:
            ut[0][0]=as;
            ut[0][1]=0.;
            ut[0][2]=0.;
            ut[1][0]=bs*cosgas;
            ut[1][1]=1./b;
            ut[1][2]= -(1./tan(al))/b;
            ut[2][0]=cs*cosbes;
            ut[2][1]=0.;
            ut[2][2]=1./(c*sin(al));

            June 1, 1996:  Corrected LFT
            The correct orthogonalization matrix for the default
            PDB convention (Cartesion x along A, y in A-B plane,
            and z along c*) is

            1/a   -cos(ga)/(a sin(ga))   as cosbes
             0          1/(b sin(ga))    bs cosals
             0               0              cs

            and the deorthogonalization matrix is

            a   b cos(ga)   c cos(be)
            0   b sin(ga)   -c sin(be) cosals
            0       0         1/cs

            usage:
            xf = x * ut[0][0] + y * ut[0][1] + z * ut[0][2]
            yf = x * ut[1][0] + y * ut[1][1] + z * ut[1][2]
            zf = x * ut[2][0] + y * ut[2][1] + z * ut[2][2]
            """

        from math import pi, cos, sin, sqrt
        degtor=pi/180.
        al=angles[0]*degtor
        be=angles[1]*degtor
        ga=angles[2]*degtor
        vol= length[0]* length[1]* length[2] *sqrt(1.-cos(al)*cos(al)-cos(be)*cos(be)
                          -cos(ga)*cos(ga)+2.*(cos(al)*cos(be)*cos(ga)))
        lAs=(length[1]* length[2]*sin(al))/vol
        bs=(length[0]* length[2]*sin(be))/vol
        cs=(length[0]* length[1]*sin(ga))/vol
        cosals=(cos(be)*cos(ga)-cos(al))/(sin(be)*sin(ga))
        cosbes=(cos(ga)*cos(al)-cos(be))/(sin(ga)*sin(al))
        cosgas=(cos(al)*cos(be)-cos(ga))/(sin(al)*sin(be))

        ut = Numeric.identity(3).astype('f')
        # sabgs1 = sqrt(1.0-cos(al)*cos(al))
        ut[0][0]=1./length[0]
        ut[0][1]= -(cos(ga))/(sin(ga)*length[0])
        # ut[0][2]= -(cos(ga)*sin(be)*cosals + cos(be)*sin(ga))/ */
        #	(sin(be)*sabgs1*sin(ga)*a) */
        ut[0][2]= lAs * cosbes
        ut[1][0]=0.0
        ut[1][1]=1./(sin(ga)*length[1])
        # ut[1][2]=cosals/(sabgs1*sin(ga)*b)
        ut[1][2]= bs * cosals
        ut[2][0]=0.0
        ut[2][1]=0.0
        # ut[2][2]=1.0/(sin(be)*sabgs1*c)
        ut[2][2]= cs
        return Numeric.transpose(ut)


    def _uinv(self, mat):
        """Inverts a transformation matrix used to go from cartesian to cell
        coordinates"""

        dmat = Numeric.identity(3).astype('d')
        imat = Numeric.identity(3).astype('d')

        dmat[0][0]=mat[1][1]*mat[2][2]-mat[2][1]*mat[1][2]
        dmat[1][0]=mat[1][2]*mat[2][0]-mat[1][0]*mat[2][2]
        dmat[2][0]=mat[1][0]*mat[2][1]-mat[1][1]*mat[2][0]
        dmat[0][1]=mat[0][2]*mat[2][1]-mat[0][1]*mat[2][2]
        dmat[1][1]=mat[0][0]*mat[2][2]-mat[0][2]*mat[2][0]
        dmat[2][1]=mat[0][1]*mat[2][0]-mat[0][0]*mat[2][1]
        dmat[0][2]=mat[0][1]*mat[1][2]-mat[0][2]*mat[1][1]
        dmat[1][2]=mat[0][2]*mat[1][0]-mat[0][0]*mat[1][2]
        dmat[2][2]=mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0]
        det = mat[0][0]*mat[1][1]*mat[2][2]+mat[1][0]*mat[2][1]*mat[0][2]+ \
              mat[2][0]*mat[0][1]*mat[1][2]-mat[2][2]*mat[0][1]*mat[1][0]- \
              mat[0][0]*mat[2][1]*mat[1][2]-mat[2][0]*mat[1][1]*mat[0][2]
        for i in (0,1,2):
            for j in (0,1,2):
                imat[j][i]=dmat[j][i]/det

        return Numeric.transpose(imat)


    def toCartesian(self, coords):
        """Transform coords from fractional coordinates to cartesian space"""

        return Numeric.dot(coords, self.ftoc)
        

    def toFractional(self, coords):
        """Transform coords from cartesian coordinates to fractional space"""

        return Numeric.dot(coords, self.ctof)


    def rotate(self, rot, coords):
        """Apply rotation to the coordinates. Rotation can be a 3x3 matrix or
        3 Euler angles"""

        if isinstance(rot[0], types.FloatType):
            rot = EulerAnglesToMat(rot)
        return Numeric.dot(coords, Numeric.transpose(rot)) 


    def translate(self, trans, coords):
        """Apply a translation to the coordinates"""
        
        return coords + trans


if __name__ == '__main__':
    cryst = Crystal( (120.3, 120.3, 78.4), (90., 90., 120.) )

    from MolKit import Read
    mol = Read('./1ctt_CA.pdb')[0]
    c = Numeric.array(mol.chains.residues.atoms.coords)

    fc = cryst.toFractional( c )
    fc1 = cryst.rotate( ((0., 1., 0.), (1., 0., 0.), (0., 0., -1.)), fc )
    fc2 = cryst.toCartesian( fc1 )
