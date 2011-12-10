## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

##  Copyright 1997-1999 by Konrad Hinsen, except as noted below.

##  Permission to use, copy, modify, and distribute this software and its
##  documentation for any purpose and without fee is hereby granted,
##  provided that the above copyright notice appear in all copies and that
##  both that copyright notice and this permission notice appear in
##  supporting documentation.

##  THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
##  INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO
##  EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, INDIRECT OR
##  CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
##  USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
##  OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
##  PERFORMANCE OF THIS SOFTWARE.

# This module defines 3d geometrical tensors with the standard
# operations on them. The elements are stored in an array.
#
# Written by: Konrad Hinsen <hinsen@cnrs-orleans.fr>
# Last revision: 1999-7-23
#

_undocumented = 1

import numpy.oldnumeric as Numeric, types
import VectorModule

class Tensor:

    """Tensor in 3D space

    Constructor: Tensor([[xx, xy, xz], [yx, yy, yz], [zx, zy, zz]])

    Tensors support the usual arithmetic operations
    ('t1', 't2': tensors, 'v': vector, 's': scalar): 

    -  't1+t2'        (addition)
    -  't1-t2'        (subtraction)
    -  't1*t2'        (tensorial (outer) product)
    -  't1*v'         (contraction with a vector, same as t1.dot(v.asTensor()))
    -  's*t1', 't1*s' (multiplication with a scalar)
    -  't1/s'         (division by a scalar)

    The coordinates can be extracted by indexing; a tensor of rank N
    can be indexed like an array of dimension N.

    Tensors are *immutable*, i.e. their elements cannot be changed.

    Tensor elements can be any objects on which the standard
    arithmetic operations are defined. However, eigenvalue calculation
    is supported only for float elements.
    """

    is_tensor = 1

    def __init__(self, elements, nocheck = None):
	self.array = Numeric.array(elements)
	if nocheck is None:
	    if not Numeric.logical_and.reduce(
		Numeric.equal(Numeric.array(self.array.shape), 3)):
		raise ValueError, 'Tensor must have length 3 along any axis'
	self.rank = len(self.array.shape)

    def __repr__(self):
	return 'Tensor(' + str(self) + ')'

    def __str__(self):
	return str(self.array)

    def __add__(self, other):
	return Tensor(self.array+other.array, 1)
    __radd__ = __add__

    def __neg__(self):
	return Tensor(-self.array, 1)

    def __sub__(self, other):
	return Tensor(self.array-other.array, 1)

    def __rsub__(self, other):
	return Tensor(other.array-self.array, 1)

    def __mul__(self, other):
	if isTensor(other):
	    a = self.array[self.rank*(slice(None),)+(Numeric.NewAxis,)]
	    b = other.array[other.rank*(slice(None),)+(Numeric.NewAxis,)]
	    return Tensor(Numeric.innerproduct(a, b), 1)
	elif VectorModule.isVector(other):
	    return other.__rmul__(self)
	else:
	    return Tensor(self.array*other, 1)

    def __rmul__(self, other):
	return Tensor(self.array*other, 1)

    def __div__(self, other):
	if isTensor(other):
	    raise TypeError, "Can't divide by a tensor"
	else:
	    return Tensor(self.array/(1.*other), 1)

    def __rdiv__(self, other):
        raise TypeError, "Can't divide by a tensor"

    def __cmp__(self, other):
	if self.rank != other.rank:
	    return 1
	else:
	    return not Numeric.logical_and.reduce(
		Numeric.equal(self.array, other.array).ravel())

    def __len__(self):
	return 3

    def __getitem__(self, index):
	elements = self.array[index]
	if type(elements) == type(self.array):
	    return Tensor(elements)
	else:
	    return elements

    def asVector(self):
        "Returns an equivalent vector object (only for rank 1)."
	if self.rank == 1:
	    return VectorModule.Vector(self.array)
	else:
	    raise ValueError, 'rank > 1'

    def dot(self, other):
        "Returns the contraction with |other|."
	if isTensor(other):
	    a = self.array
	    b =  Numeric.transpose(other.array, range(1, other.rank)+[0])
	    return Tensor(Numeric.innerproduct(a, b), 1)
	else:
	    return Tensor(self.array*other, 1)

    def diagonal(self, axis1=0, axis2=1):
	if self.rank == 2:
	    return Tensor([self.array[0,0], self.array[1,1], self.array[2,2]])
	else:
	    if axis2 < axis1: axis1, axis2 = axis2, axis1
	    raise ValueError, 'Not yet implemented'

    def trace(self, axis1=0, axis2=1):
        "Returns the trace of a rank-2 tensor."
	if self.rank == 2:
	    return self.array[0,0]+self.array[1,1]+self.array[2,2]
	else:
	    raise ValueError, 'Not yet implemented'

    def transpose(self):
        "Returns the transposed (index reversed) tensor."
	return Tensor(Numeric.transpose(self.array))

    def symmetricalPart(self):
        "Returns the symmetrical part of a rank-2 tensor."
	if self.rank == 2:
	    return Tensor(0.5*(self.array + \
			       Numeric.transpose(self.array,
						    Numeric.array([1,0]))),
			  1)
	else:
	    raise ValueError, 'Not yet implemented'

    def asymmetricalPart(self):
        "Returns the asymmetrical part of a rank-2 tensor."
	if self.rank == 2:
	    return Tensor(0.5*(self.array - \
			       Numeric.transpose(self.array,
						    Numeric.array([1,0]))),
			  1)
	else:
	    raise ValueError, 'Not yet implemented'

    def eigenvalues(self):
        "Returns the eigenvalues of a rank-2 tensor in an array."
	if self.rank == 2:
	    from numpy.oldnumeric.linear_algebra import eigenvalues
	    return eigenvalues(self.array)
	else:
	    raise ValueError, 'Undefined operation'

    def diagonalization(self):
        """Returns the eigenvalues of a rank-2 tensor and a tensor
        representing the rotation matrix to the diagonalized form."""
	if self.rank == 2:
	    from numpy.oldnumeric.linear_algebra import eigenvectors
            ev, vectors = eigenvectors(self.array)
	    return ev, Tensor(vectors)
	else:
	    raise ValueError, 'Undefined operation'

    def inverse(self):
        "Returns the inverse of a rank-2 tensor."
	if self.rank == 2:
	    from numpy.oldnumeric.linear_algebra import inverse
	    return Tensor(inverse(self.array))
	else:
	    raise ValueError, 'Undefined operation'

# Type check

def isTensor(x):
    "Return 1 if |x| is a tensor."
    return hasattr(x,'is_tensor')


# Some useful constant tensors

delta = Tensor([[1,0,0], [0,1,0], [0,0,1]])
_epsilon = Numeric.zeros((3,3,3))
_epsilon[0,1,2] = 1
_epsilon[1,2,0] = 1
_epsilon[2,0,1] = 1
_epsilon[0,2,1] = -1
_epsilon[2,1,0] = -1
_epsilon[1,0,2] = -1
epsilon = Tensor(_epsilon)
del _epsilon
