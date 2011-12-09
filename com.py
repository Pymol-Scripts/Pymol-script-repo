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
'''
from pymol import cmd
from chempy import cpv
def COM(selection='all', center=0, quiet=1):

        model = cmd.get_model(selection)
        nAtom = len(model.atom)
 
        COM = cpv.get_null()
 
        for a in model.atom:
                COM = cpv.add(COM, a.coord)
        COM = cpv.scale(COM, 1./nAtom)
 
        if not int(quiet):
                print ' COM: [%8.3f,%8.3f,%8.3f]' % tuple(COM)
 
        if int(center):
                cmd.alter_state(1, selection, "(x,y,z)=sub((x,y,z), COM)",
                        space={'COM': COM, 'sub': cpv.sub})
 
        return COM
 
cmd.extend("COM", COM)
