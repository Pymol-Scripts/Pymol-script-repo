'''
See more here: http://www.pymolwiki.org/index.php/perp_maker
####################################################################################
#
# perp_maker.py:  Creates perpendicular planes.
# =============
#
# Nothing to do with cops.  Given a simple PyMol scene, attempts to
# create a CGO background triangle perpendicular to the vector created - which is
# parallel to the line segment drawn through the camera point and current center of
# mass - as obtained by "get_position," or "get_view."
#
# @COPYRIGHT: Jason Vertrees (C), 2005-2007
# @LICENSE: Released under GPL:
# This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
# Street, Fifth Floor, Boston, MA 02110-1301, USA
#
#
#

#
# To use: Load your scene.  Orient the scene as you wish.  Run the script.
# Could it be any simpler?!
#

# The TTT Matrix has to be the identity, do achieve this result.  So,
# run the following:
#   -- 'reset'
#   -- then orient your molecule as desired using the EDITING features!
#     -- before running this script, make sure 'get_view' shows the identity
#     -- matrix for the first 9 elements.
#   -- then run the script
#
'''

from __future__ import print_function

import pymol
import random
from pymol.cgo import *

############################################################
#
# Methods
#
############################################################

#
# Given the viewVector and center, creates a random sized plane
# perpendicular to the viewVector through the origin.  It is then
# the next step's responsibility to move the plane back some so
# it dosen't cut the molecule/scene in half.
#


def getPPlane(viewVector, center, side_length=100):
    """Returns a 3-tuple of 3D points representing the perp. plane."""

    # for reproduceable testing
    # random.seed(10)
    #
    # The formula for a plane with our chacteristics is defined by
    #
    # A(x - x') + B(y - y') + C(y - y') + D = 0, where
    # A, B and C are not all zero coefficients in the vector
    # Ai + Bj + Ck such that the plane is perpendicular to this
    # vector; x, y, and z are points on the plane; x', y', and z'
    # are the coordiates through which the plane shall run.
    #

    # This is fool-ass.  Gotta' be a better way to do this.
    # Declaring that rVal is a 3-Tuple.
    rVal = [[], [], [], [], [], []]

    # Compose two triangles into a square.
    # Never learned any GFX coding, so I'm sure there's something
    # better than this; but, this works.
    for i in range(0, 6):
        if (i == 0) or (i == 5):
            x = -side_length + center[0]
            y = -side_length + center[1]
        elif (i == 1):
            x = -side_length + center[0]
            y = side_length + center[1]
        elif (i == 2) or (i == 3):
            x = side_length + center[0]
            y = side_length + center[1]
        elif (i == 4):
            x = side_length + center[0]
            y = -side_length + center[1]

        if (viewVector[2] != 0):
            z = -(((viewVector[0] * (x - center[0])) - (viewVector[1] * (y - center[1]))) /
                  viewVector[2]) + center[2]

        else:
            print("Z-component of viewVector is zero.  Now, I need a nonzero value here \
so I'm just making one up. :)")
            z = random.randint(-200, 200)

        rVal[i] = [x, y, z]

    return rVal

############################################################
#
# End methods
#
############################################################


def perp_maker(name='pPlane', quiet=1):
    '''
DESCRIPTION

    Creates perpendicular planes
    '''
    quiet = int(quiet)

# First, get the center and camera locations
    view = cmd.get_view()
    camera = [view[9], view[10], view[11]]
    center = [view[12], view[13], view[14]]

# Sanity check
    if not quiet:
        print("Camera is: " + str(camera))
        print("Center is: " + str(center))

# Create the vector through the two points directed
# from the camera to the center - the viewVector
    viewVector = [center[0] - camera[0],
                   center[1] - camera[1],
                   center[2] - camera[2]]

    if not quiet:
        print("ViewVector is: " + str(viewVector))

# Create the plane perpendicular to the viewVector
# running through the origin
    pPlane = getPPlane(viewVector, center, side_length=100)
    if not quiet:
        print("Plane points calculated as: " + str(pPlane))

# now create the CGO and load from the points
    obj = [
        BEGIN, TRIANGLES,
        COLOR, 0.2, 0.4, 1,

        VERTEX, pPlane[0][0], pPlane[0][1], pPlane[0][2],
        VERTEX, pPlane[1][0], pPlane[1][1], pPlane[1][2],
        VERTEX, pPlane[2][0], pPlane[2][1], pPlane[2][2],

        VERTEX, pPlane[3][0], pPlane[3][1], pPlane[3][2],
        VERTEX, pPlane[4][0], pPlane[4][1], pPlane[4][2],
        VERTEX, pPlane[5][0], pPlane[5][1], pPlane[5][2],

        END
    ]

    cmd.load_cgo(obj, name)
    cmd.set_view(view)

if __name__ in ['pymol', '__main__']:
    print('__name__ =', __name__)
    perp_maker(quiet=0)

cmd.extend('perp_maker', perp_maker)

# vi:expandtab:smarttab
