'''
See more here: http://www.pymolwiki.org/index.php/com

DESCRIPTION
        get the center of mass of selection or move selection to the origin.
 
ARGUMENTS
        selection = string: a valid PyMOL selection {default: all}
        center = 0 or 1: if center=1 center the selection {default: 0}
        returns: center of mass: [ xCOM, yCOM, zCOM ]
 
SEE ALSO
        get_extent, get_position, http://pymolwiki.org/index.php/Center_Of_Mass

# @AUTHOR: Jason Vertrees
# Copyright (c) 2008, Jason Vertrees
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following
# conditions are met:
#
#     * Redistributions of source code must retain the above copyright notice, this list of conditions and the following
#     * disclaimer.
#     * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following
#     * disclaimer in the documentation and/or other materials provided with the distribution.
#     * Neither the name of the <ORGANIZATION> nor the names of its contributors may be used to endorse or promote products derived
#     * from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
# THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
# OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# DATE  : 2008-09-26
# REV   : 1

'''
from pymol import cmd
from pymol import stored
from chempy import cpv


def COM(selection='all', center=0, quiet=1):

    model = cmd.get_model(selection)
    nAtom = len(model.atom)

    COM = cpv.get_null()

    for a in model.atom:
        COM = cpv.add(COM, a.coord)
    COM = cpv.scale(COM, 1. / nAtom)

    if not int(quiet):
        print ' COM: [%8.3f,%8.3f,%8.3f]' % tuple(COM)

    if int(center):
        cmd.alter_state(1, selection, "(x,y,z)=sub((x,y,z), COM)",
                        space={'COM': COM, 'sub': cpv.sub})

    return COM

cmd.extend("COM", COM)
