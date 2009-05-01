import pymol
from pymol import cmd
from pymol.wizard import Wizard
from chempy import cpv
from cgo import *
 
def makePrimitive(cgo, name):
    cmd.set('auto_zoom', 0, quiet=1)
    cmd.load_cgo(cgo, name)
    cmd.set('auto_zoom', 1, quiet=1)
 
def point(p):
    x, y, z = p
    return [COLOR, 1, 1, 1, SPHERE, float(x), float(y), float(z), 0.5]
 
def line(p1, p2):
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    return [CYLINDER, float(x1), float(y1), float(z1), float(x2), float(y2), float(z2), 0.25, 1, 1, 1, 1, 1, 1]
 
def plane(corner1, corner2, corner3, corner4, normal):
    planeObj = []
    planeObj.extend(point(corner1))
    planeObj.extend(point(corner2))
    planeObj.extend(point(corner3))
    planeObj.extend(point(corner4))
    planeObj.extend(line(corner1, corner2))
    planeObj.extend(line(corner2, corner3))
    planeObj.extend(line(corner3, corner4))
    planeObj.extend(line(corner4, corner1))
 
    planeObj.extend([COLOR, 0.8, 0.8, 0.8])
    planeObj.extend([BEGIN, TRIANGLE_STRIP])
    planeObj.append(NORMAL)
    planeObj.extend(normal)
    for corner in [corner1, corner2, corner3, corner4, corner1]:
        planeObj.append(VERTEX)
        planeObj.extend(corner)
    planeObj.append(END)
    return planeObj
 
def planeFromPoints(point1, point2, point3, facetSize):
    v1 = cpv.normalize(cpv.sub(point2, point1))
    v2 = cpv.normalize(cpv.sub(point3, point1))
    normal = cpv.cross_product(v1, v2)
    v2 = cpv.cross_product(normal, v1)
    x = cpv.scale(v1, facetSize)
    y = cpv.scale(v2, facetSize)
    center = point2
    corner1 = cpv.add(cpv.add(center, x), y)
    corner2 = cpv.sub(cpv.add(center, x), y)
    corner3 = cpv.sub(cpv.sub(center, x), y)
    corner4 = cpv.add(cpv.sub(center, x), y)
    return plane(corner1, corner2, corner3, corner4, normal)
 
 
class PlaneWizard(Wizard):
 
    def __init__(self):
        Wizard.__init__(self)
 
        # some attributes to do with picking
        self.pick_count = 0
        self.object_count = 0
        self.object_prefix = "pw"
 
        # the plane facet size (the 'radius' of the section of plane we show)
        self.facetSize = 5
 
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
        except pymol.CmdException, pmce:
            print pmce
 
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
            print self.error
            return
 
        atom_name = self.object_prefix + str(self.pick_count)
        if self.pick_count < 2:
            self.pickNextAtom(atom_name)
        else:
            self.pickNextAtom(atom_name)
 
            point1 = cmd.get_atom_coords("(%s%s)" % (self.object_prefix, "0"))
            point2 = cmd.get_atom_coords("(%s%s)" % (self.object_prefix, "1"))
            point3 = cmd.get_atom_coords("(%s%s)" % (self.object_prefix, "2"))
            plane = planeFromPoints(point1, point2, point3, self.facetSize)
 
            planeName = "plane-%02d" % self.object_count
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
 
# create an instance
 
wiz = PlaneWizard()
 
# make this the active wizard
 
cmd.set_wizard(wiz)
