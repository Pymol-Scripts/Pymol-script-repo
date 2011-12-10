## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#
# Last modified on Tue Sep  4 17:02:59 PDT 2001 by lindy
#
# $Header: /opt/cvs/python/packages/share1.5/mglutil/math/ncoordstest.py,v 1.2 2007/07/24 17:30:40 vareille Exp $
#

"""Unit test for ncoords.py

Requirements for ncoords.py:
    A. __init__:
        1. make Numeric, homogenious coordinates out of refCoords
        2. raise ValueError if refCoords is bad
    B. reset():
        3. nx3 slice of resultCoords must be equal to refCoords
    C. getResultCoords:
        4. return nx3 (not nx4) coordinates
        5. return as Numeric.array or ListType accorinding to self.tolist
"""

from mglutil.math.ncoords import Ncoords
import unittest, math
import numpy.oldnumeric as Numeric, numpy.oldnumeric.random_array as RandomArray




class NcoordsTest(unittest.TestCase):
    def setUp(self):
        """Called for every test."""

        npts = 500
        dim = 3
        self.max = 9999999.
        self.min = -self.max
        self.random_points = RandomArray.uniform(self.min,
                                                 self.max, (npts,dim)).tolist()


    def tearDown(self):
        pass



class InputOutputValues(NcoordsTest):


    def test_constructor_shape(self):
        """__init__         -- make refCoords and resultCoords homogeneous"""
        n = len(self.random_points)
        ncoords = Ncoords( self.random_points) ### tested call ###
        # confirm shape to be nx4
        self.assertEqual( (n, 4), Numeric.shape(ncoords.resultCoords))
        self.assertEqual( (n, 4), Numeric.shape(ncoords.refCoords))
        # cofirm that the last column is all ones
        self.assertEqual(Numeric.ones(n).tolist(),
                         ncoords.resultCoords[:,3].tolist())
        self.assertEqual(Numeric.ones(n).tolist(),
                         ncoords.refCoords[:,3].tolist())


    def test_input_error(self):
        """__init__         -- ValueError on bad input"""
        self.assertRaises(ValueError, Ncoords, range(10))
        self.assertRaises(ValueError, Ncoords, [(1,1,1),(1,1)] )


    def test_reset_values(self):
        """reset            -- points equal input values after reset"""
        nc = Ncoords( self.random_points, tolist=1)
        nc.reset() ### tested call ###
        result = nc.getResultCoords()
        # compare input and output point lists
        self.assertEqual( self.random_points, result)


    def test_getResultCoords_shape(self):
        """getResultCoords  -- if tolist: return nx3 ListType"""
        n = len(self.random_points)
        nc = Ncoords(self.random_points, tolist=0)
        nc.tolist=1
        result = nc.getResultCoords() ### tested call ###
        # confirm shape
        self.assertEqual((n, 3), Numeric.shape(result))
        # confirm type
        self.assertEqual(type([]), type(result))


    def test_getResultCoords_type(self):
        """getResultCoords  -- if not tolist: return nx4 Numeric.array"""
        n = len(self.random_points)
        nc = Ncoords(self.random_points, tolist=1)
        nc.tolist=0
        result = nc.getResultCoords() ### tested call ###
        # confirm shape
        self.assertEqual((n, 4), Numeric.shape(result))
        # confirm type
        self.assertEqual(type(Numeric.array([])), type(result))



if __name__ == '__main__':
    unittest.main()   

# for example:
#     py mglutil/math/ncoordstest.py -v
# or, to redirect output to a file:
#     py ncoordstest.py -v > & ! /tmp/nct.out






