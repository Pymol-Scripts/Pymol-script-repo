'''
See more here: http://www.pymolwiki.org/index.php/color_by_conservation

    PARAMETERS
    aln
        (string) the alignment object name
    names
        (list) a list of object names that are in the alignment;
               if (), then PyMOL will attempt to glean the names
               from the alignment object
    color
        (string) valid PyMOL spectrum name

    as_putty
        (0 or 1) if 0 display is not changed, else participating objects are shown
                 as cartoon putty, colored by the 'color' field
'''

from __future__ import print_function
from pymol import cmd


def color_by_conservation(aln, names=(), color="rainbow", as_putty=0, _self=cmd):
    # PyMOL doesn't yet know about object:alignment
    # but we need to check that this exists or we might crash
    if _self.get_type(aln) not in ("object:", "object:alignment"):
        print("Error: Bad or incorrectly specified alignment object.")
        return None

    r = cmd.get_raw_alignment(aln)

    if names == ():
        known_objs = []
        list(map(known_objs.extend, [[y[0] for y in x] for x in r]))
        known_objs = set(known_objs)

        # highest number of matches seen
        M = max(list(map(len, r))) + 1
    else:
        known_objs = set(names)
        M = len(known_objs) + 1

    for obj in known_objs:
        _self.alter(obj, "b=0.0")

    for af in r:
        c = float(1.0 + len(af)) / float(M)
        for y in af:
            _self.alter("%s and index %s" % (y[0], y[1]), "b=c", space={'c': c})

    if as_putty != 0:
        for obj in known_objs:
            _self.show_as("cartoon", "%s" % obj)
            _self.cartoon("putty", "%s" % obj)
            _self.spectrum('b', color, obj)
            _self.sort()
            _self.rebuild()
    return None
cmd.extend("color_by_conservation", color_by_conservation)
