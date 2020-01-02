'''
More information at: http://www.pymolwiki.org/index.php/elbow_angle

Calculate the elbow angle of an antibody Fab complex and optionally draw a
graphical representation of the vectors used to determine the angle.

NOTE: There is no automatic checking of the validity of limit_l and limit_h
values or of the assignment of light and heavy chain IDs. If these are entered
incorrectly or omitted, the reported angle will likely be incorrect.

As always with these things, your mileage may vary.  Use at your own risk!

REQUIREMENTS

    numpy, version 1.6
        http://numpy.scipy.org

    transformations.py, version 2012.01.01
        by Christoph Gohlke
        www.lfd.uci.edu/~gohlke/code

        May also require an edit to transformations.py:
        Changes `1e-8` to `1e-7` in lines 357 & 363 to avoid a numerical error.

    com.py
        by Jason Vertrees
        http://www.pymolwiki.org/index.php/com


CHANGELOG
=========
* 0.1.1
  - fixed: Print an error message instead of failing silently when a limit
    residue doesn't exist.

'''

from __future__ import print_function

__author__ = 'Jared Sampson'
__version__ = '0.1.1'


from pymol import cmd
import pymol
import numpy

try:
    import transformations
except ImportError:
    from . import transformations

try:
    COM = cmd.centerofmass
except AttributeError:
    from com import COM


################################################################################
def calc_super_matrix(mobile, static):
    '''

DESCRIPTION

    Aligns two objects (or selections), returns the transformation matrix,
    and resets the matrix of the mobile object.

    Uses CEAlign PyMOL function for alignment.

ARGUMENTS

    mobile = string: selection describing the mobile object whose rotation
    matrix will be reported

    static = string: selection describing the static object onto which the
    mobile object will be aligned

REQUIRES: numpy
    '''

    cmd.cealign(static, mobile)
#    cmd.super(mobile,static)
    T = cmd.get_object_matrix(mobile)

    R = numpy.identity(4)
    k = 0
    for i in range(0, 4):
        for j in range(0, 4):
            R[i][j] = T[k]
            k += 1

    return R


################################################################################
def elbow_angle(obj, light='L', heavy='H', limit_l=107, limit_h=113, draw=0):
    """

DESCRIPTION

    Calculates the integer elbow angle of an antibody Fab complex and
    optionally draws a graphical representation of the vectors used to
    determine the angle.

ARGUMENTS

    obj = string: object

    light/heavy = strings: chain ID of light and heavy chains, respectively

    limit_l/limit_h = integers: residue numbers of the last residue in the
    light and heavy chain variable domains, respectively

    draw = boolean: Choose whether or not to draw the angle visualization

REQUIRES: com.py, transformations.py, numpy (see above)


    """

    # store current view
    orig_view = cmd.get_view()

    limit_l = int(limit_l)
    limit_h = int(limit_h)
    draw = int(draw)

    # for temp object names
    tmp_prefix = "tmp_elbow_"

    prefix = tmp_prefix + obj + '_'

    # names
    vl = prefix + 'VL'
    vh = prefix + 'VH'
    cl = prefix + 'CL'
    ch = prefix + 'CH'

    # selections
    vl_sel = 'polymer and %s and chain %s and resi 1-%i' % (obj, light, limit_l)
    vh_sel = 'polymer and %s and chain %s and resi 1-%i' % (obj, heavy, limit_h)
    cl_sel = 'polymer and %s and chain %s and not resi 1-%i' % (obj, light, limit_l)
    ch_sel = 'polymer and %s and chain %s and not resi 1-%i' % (obj, heavy, limit_h)
    v_sel = '((' + vl_sel + ') or (' + vh_sel + '))'
    c_sel = '((' + cl_sel + ') or (' + ch_sel + '))'

    # create temp objects
    cmd.create(vl, vl_sel)
    cmd.create(vh, vh_sel)
    cmd.create(cl, cl_sel)
    cmd.create(ch, ch_sel)

    # superimpose vl onto vh, calculate axis and angle
    Rv = calc_super_matrix(vl, vh)
    angle_v, direction_v, point_v = transformations.rotation_from_matrix(Rv)

    # superimpose cl onto ch, calculate axis and angle
    Rc = calc_super_matrix(cl, ch)
    angle_c, direction_c, point_c = transformations.rotation_from_matrix(Rc)

    # delete temporary objects
    cmd.delete(vl)
    cmd.delete(vh)
    cmd.delete(cl)
    cmd.delete(ch)

    # if dot product is positive, angle is acute
    if (numpy.dot(direction_v, direction_c) > 0):
        direction_c = direction_c * -1   # ensure angle is > 90 (need to standardize this)

        # TODO: make both directions point away from the elbow axis.

    elbow = int(numpy.degrees(numpy.arccos(numpy.dot(direction_v, direction_c))))
    # while (elbow < 90):
    # elbow = 180 - elbow   # limit to physically reasonable range

    # compare the direction_v and direction_c axes to the vector defined by
    # the C-alpha atoms of limit_l and limit_h of the original fab
    hinge_l_sel = "%s//%s/%s/CA" % (obj, light, limit_l)
    hinge_h_sel = "%s//%s/%s/CA" % (obj, heavy, limit_h)

    try:
        hinge_l = cmd.get_atom_coords(hinge_l_sel)
        hinge_h = cmd.get_atom_coords(hinge_h_sel)
    except pymol.CmdException:
        # Either hinge_l_sel or hinge_h_sel atom did not exist.
        raise pymol.CmdException(
            'Unable to calculate elbow angle. Please check '
            'your limit and chain selections and try again.'
        )
    hinge_vec = numpy.array(hinge_h) - numpy.array(hinge_l)

    test = numpy.dot(hinge_vec, numpy.cross(direction_v, direction_c))
    if (test > 0):
        elbow = 360 - elbow

    print("    Elbow angle: %i degrees" % elbow)

    if (draw == 1):
        # there is probably a more elegant way to do this, but
        # it works so I'm not going to mess with it for now

        pre = obj + '_elbow_'

        # draw hinge vector
        cmd.pseudoatom(pre + "hinge_l", pos=hinge_l)
        cmd.pseudoatom(pre + "hinge_h", pos=hinge_h)
        cmd.distance(pre + "hinge_vec", pre + "hinge_l", pre + "hinge_h")
        cmd.set("dash_gap", 0)

        # draw the variable domain axis
        com_v = COM(v_sel)
        start_v = [a - 10 * b for a, b in zip(com_v, direction_v)]
        end_v = [a + 10 * b for a, b in zip(com_v, direction_v)]
        cmd.pseudoatom(pre + "start_v", pos=start_v)
        cmd.pseudoatom(pre + "end_v", pos=end_v)
        cmd.distance(pre + "v_vec", pre + "start_v", pre + "end_v")

        # draw the constant domain axis
        com_c = COM(c_sel)
        start_c = [a - 10 * b for a, b in zip(com_c, direction_c)]
        end_c = [a + 10 * b for a, b in zip(com_c, direction_c)]
        cmd.pseudoatom(pre + "start_c", pos=start_c)
        cmd.pseudoatom(pre + "end_c", pos=end_c)
        cmd.distance(pre + "c_vec", pre + "start_c", pre + "end_c")

        # customize appearance
        cmd.hide("labels", pre + "hinge_vec")
        cmd.hide("labels", pre + "v_vec")
        cmd.hide("labels", pre + "c_vec")
        cmd.color("green", pre + "hinge_l")
        cmd.color("red", pre + "hinge_h")
        cmd.color("black", pre + "hinge_vec")
        cmd.color("black", pre + "start_v")
        cmd.color("black", pre + "end_v")
        cmd.color("black", pre + "v_vec")
        cmd.color("black", pre + "start_c")
        cmd.color("black", pre + "end_c")
        cmd.color("black", pre + "c_vec")
        # draw spheres
        cmd.show("spheres", pre + "hinge_l or " + pre + "hinge_h")
        cmd.show("spheres", pre + "start_v or " + pre + "start_c")
        cmd.show("spheres", pre + "end_v or " + pre + "end_c")
        cmd.set("sphere_scale", 2)
        cmd.set("dash_gap", 0, pre + "hinge_vec")
        cmd.set("dash_width", 5)
        cmd.set("dash_radius", 0.3)

        # group drawing objects
        cmd.group(pre, pre + "*")

    # restore original view
    cmd.set_view(orig_view)

    return 0

cmd.extend("elbow_angle", elbow_angle)
