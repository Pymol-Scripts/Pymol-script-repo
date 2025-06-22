'''
Described at PyMOL wiki:
https://pymolwiki.org/index.php/Plane_Wizard

Authors : Troels Schwarz-Linnet
Date    : Dec 2016
Modified: From previous contributors.
'''

import pymol
from pymol import cmd
from pymol.wizard import Wizard
from chempy import cpv
from pymol.cgo import COLOR, SPHERE, CYLINDER, BEGIN, TRIANGLE_STRIP, NORMAL, VERTEX, END, ALPHA

def makePrimitive(cgo, name):
    az = cmd.get('auto_zoom', quiet=1)
    cmd.set('auto_zoom', 0, quiet=1)
    cmd.load_cgo(cgo, name)
    cmd.set('auto_zoom', az, quiet=1)

def point(p):
    x, y, z = p
    return [COLOR, 1, 1, 1, SPHERE, float(x), float(y), float(z), 0.5]


def line(p1, p2):
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    return [CYLINDER, float(x1), float(y1), float(z1), float(x2), float(y2), float(z2), 0.25, 1, 1, 1, 1, 1, 1]


def plane(corner1, corner2, corner3, corner4, normal, settings):
    planeObj = []
    planeObj.extend(point(corner1))
    planeObj.extend(point(corner2))
    planeObj.extend(point(corner3))
    planeObj.extend(point(corner4))
    planeObj.extend(line(corner1, corner2))
    planeObj.extend(line(corner2, corner3))
    planeObj.extend(line(corner3, corner4))
    planeObj.extend(line(corner4, corner1))

    # Make settings
    if 'ALPHA' in settings:
        planeObj.extend([ALPHA, settings['ALPHA']])

    if 'COLOR' in settings:
        planeObj.extend([COLOR, settings['COLOR'][0], settings['COLOR'][1], settings['COLOR'][2]])
    else:
        planeObj.extend([COLOR, 0.8, 0.8, 0.8]) # greyish

    planeObj.extend([BEGIN, TRIANGLE_STRIP])
    planeObj.append(NORMAL)

    if 'INVERT' in settings:
        if settings['INVERT']==True:
            planeObj.extend(cpv.negate(normal))
        else:
            planeObj.extend(normal)
    else:
        planeObj.extend(normal)


    for corner in [corner1, corner2, corner3, corner4, corner1]:
        planeObj.append(VERTEX)
        planeObj.extend(corner)
    planeObj.append(END)
    return planeObj


def planeFromPoints(p1, p2, p3, vm1=None, vm2=None, center=True, settings={}):
    v1 = cpv.sub(p1, p2)
    v2 = cpv.sub(p3, p2)
    normal = cpv.cross_product(v1, v2)

    if 'translate' in settings:
        vtran = cpv.scale(cpv.normalize(normal), settings['translate'])
        p1_t = cpv.sub(p1, vtran)
        p2_t = cpv.sub(p2, vtran)
        p3_t = cpv.sub(p3, vtran)
        print("New coordinates are:")
        print_info("New", p1_t, p2_t, p3_t)
        print("New coordinates are for normalized plane:")
        v1_t = cpv.normalize(cpv.sub(p1_t, p2_t))
        v2_t = cpv.normalize(cpv.sub(p3_t, p2_t))
        normal_t = cpv.normalize(cpv.cross_product(v1_t, v2_t))
        v2_t = cpv.normalize(cpv.cross_product(normal_t, v1_t))
        p1_t2 = cpv.add(v1_t, p2_t)
        p3_t2 = cpv.add(v2_t, p2_t)
        print_info("Newnormal", p1_t2, p2_t, p3_t2)

    if vm1!=None:
        v1 = cpv.scale(cpv.normalize(v1), vm1)
    if vm2!=None:
        v2 = cpv.scale(cpv.normalize(v2), vm2)

    centrum = p2
    if center:
        corner1 = cpv.add(cpv.add(centrum, v1), v2)
        corner2 = cpv.sub(cpv.add(centrum, v1), v2)
        corner3 = cpv.sub(cpv.sub(centrum, v1), v2)
        corner4 = cpv.add(cpv.sub(centrum, v1), v2)
    else:
        corner1 = cpv.add(cpv.add(centrum, v1), v2)
        corner2 = cpv.add(centrum, v1)
        corner3 = centrum
        corner4 = cpv.add(centrum, v2)

    return plane(corner1, corner2, corner3, corner4, normal, settings)


def print_info(name, coor1, coor2, coor3):
    cs1 = (map(float, [ '%.2f' % elem for elem in coor1 ]) )
    cs2 = (map(float, [ '%.2f' % elem for elem in coor2 ]) )
    cs3 = (map(float, [ '%.2f' % elem for elem in coor3 ]) )
    print("You can also use the function calls with these coordinates")
    print("plane.make_plane_points(name='%s', l1=%s, l2=%s, l3=%s)"%(name, cs1, cs2, cs3))


def make_plane(name,a1='(pk1)',a2='(pk2)',a3='(pk3)', vm1=None, vm2=None, center=True, makepseudo=True, settings={}):
    """
    DESCRIPTION
    Create a CGO plane from three atomic coordinates

    USAGE
    make_plane name, a1, a2, a3

    where each atom is a standard PyMOL selection (defaults to pk1,pk2 and pk3)
    """
    # get coordinates from atom selections
    coor1 = cmd.get_model(a1).get_coord_list()[0]
    coor2 = cmd.get_model(a2).get_coord_list()[0]
    coor3 = cmd.get_model(a3).get_coord_list()[0]

    # Help with alternative
    print_info(name, coor1, coor2, coor3)

    # Get the plane
    plane = planeFromPoints(p1=coor1, p2=coor2, p3=coor3, vm1=vm1, vm2=vm2, center=center, settings=settings)
    makePrimitive(plane, name)
    #cmd.show("cgo", "plane*")

    if makepseudo:
        cmd.pseudoatom("%s_%s"%(name, "l1"), color="tv_blue", pos=coor1)
        cmd.pseudoatom("%s_%s"%(name, "l2"), color="tv_green", pos=coor2)
        cmd.pseudoatom("%s_%s"%(name, "l3"), color="red", pos=coor3)

# Extend function to be called inside pymol
cmd.extend("make_plane", make_plane)

def make_plane_points(name,l1=None,l2=None,l3=None, vm1=None, vm2=None, center=True, makepseudo=True, settings={}):
    """
    DESCRIPTION
    Create a CGO plane from three atomic coordinates

    USAGE
    make_plane name, l1, l2, l3

    where each xys is a list with floats of x,y,z coordinates
    """
    if l1==None or l2==None or l3==None:
        print("Please provide a list of xyz floats for each 3 positions")
        return
    if type(l1) is not list or type(l2) is not list or type(l3) is not list:
        print(type(l1),type(l2),type(l3))
        print("Please provide 3 list of xyz floats for each 3 positions")
        return

    plane = planeFromPoints(p1=l1, p2=l2, p3=l3, vm1=vm1, vm2=vm2, center=center, settings=settings)
    makePrimitive(plane, name)

    if makepseudo:
        cmd.pseudoatom("%s_%s"%(name, "l1"), color="tv_blue", pos=l1)
        cmd.pseudoatom("%s_%s"%(name, "l2"), color="tv_green", pos=l2)
        cmd.pseudoatom("%s_%s"%(name, "l3"), color="red", pos=l3)

# Extend function to be called inside pymol
cmd.extend("make_plane_points", make_plane_points)

class PlaneWizard(Wizard):
    def __init__(self):
        Wizard.__init__(self)

        # some attributes to do with picking
        self.pick_count = 0
        self.object_count = 0
        self.object_prefix = "pw"

        self.selection_mode = cmd.get_setting_legacy("mouse_selection_mode")
        cmd.set("mouse_selection_mode",0) # set selection mode to atomic
        cmd.deselect()

    def reset(self):
        cmd.delete(self.object_prefix + "*")
        cmd.delete("sele*")
        cmd.delete("_indicate*")
        cmd.unpick()
        cmd.refresh_wizard()

    def delete_all(self):
        cmd.delete("plane*")

    def cleanup(self):
        cmd.set("mouse_selection_mode",self.selection_mode) # restore selection mode
        self.reset()
        self.delete_all()

    def get_prompt(self):
        self.prompt = None
        if self.pick_count == 0:
            self.prompt = [ 'Please click on the first atom...']
        elif self.pick_count == 1:
            self.prompt = [ 'Please click on the second atom...' ]
        elif self.pick_count == 2:
            self.prompt = [ 'Please click on the third atom...' ]
        return self.prompt

    def do_select(self, name):
        # "edit" only this atom, and not others with the object prefix
        try:
            cmd.edit("%s and not %s*" % (name, self.object_prefix))
            self.do_pick(0)
        except pymol.CmdException as pmce:
            print(pmce)

    def pickNextAtom(self, atom_name):
        # transfer the click selection to a named selection
        cmd.select(atom_name, "(pk1)")

        # delete the click selection
        cmd.unpick()

        # using the magic of indicate, highlight stuff
        indicate_selection = "_indicate" + self.object_prefix
        cmd.select(indicate_selection, atom_name)
        cmd.enable(indicate_selection)

        self.pick_count += 1
        self.error = None

        # necessary to force update of the prompt
        cmd.refresh_wizard()

    def do_pick(self, picked_bond):

        # this shouldn't actually happen if going through the "do_select"
        if picked_bond:
            self.error = "Error: please select bonds, not atoms"
            print(self.error)
            return

        atom_name = self.object_prefix + str(self.pick_count)
        if self.pick_count < 2:
            self.pickNextAtom(atom_name)
        else:
            self.pickNextAtom(atom_name)

            point1 = cmd.get_atom_coords("(%s%s)" % (self.object_prefix, "0"))
            point2 = cmd.get_atom_coords("(%s%s)" % (self.object_prefix, "1"))
            point3 = cmd.get_atom_coords("(%s%s)" % (self.object_prefix, "2"))
            plane = planeFromPoints(point1, point2, point3)

            planeName = "plane-%02d" % self.object_count

            print_info(planeName, point1, point2, point3)

            self.object_count += 1
            makePrimitive(plane, planeName)
            cmd.show("cgo", "plane*")

            self.pick_count = 0
            self.reset()

    def get_panel(self):
        return [
            [ 1, 'Plane Wizard',''],
            [ 2, 'Reset','cmd.get_wizard().reset()'],
            [ 2, 'Delete All Planes' , 'cmd.get_wizard().delete_all()'],
            [ 2, 'Done','cmd.set_wizard()'],
        ]
