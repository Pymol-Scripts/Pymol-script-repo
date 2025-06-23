# Author: Takanori Nakane
# License: BSD-2 Clause
# Version: 0.3.20120830

'''
Dynamic Mesh

This script was tested on PyMOL 1.2 and 1.5.

Example:

run dynamic_mesh.py
fetch 1hwk, async=0
fetch 1hwk, 1hwk_map, type=2fofc, async=0
dynamic_mesh 1hwk_map, sym_source=1hwk
show sticks, resn 117
show ribbon
zoom chain A and resn 117

Note: On PyMOL <= 1.4, you have to download the electron density
map from the Uppsala Electron Density Server manually.
'''

from __future__ import print_function
from pymol.callback import Callback
from pymol import cmd
from chempy import cpv


class DynamicMesh(Callback):

    def __init__(self, map_name, level, radius, name, sym_source):
        self.level = level
        self.radius = radius
        self.map_name = map_name
        self.name = name
        self.center_name = cmd.get_unused_name('_center')
        self.callback_name = cmd.get_unused_name('_cb')

        cmd.set("auto_zoom", 0)
        cmd.pseudoatom(self.center_name)
        cmd.hide("everything", self.center_name)

        symmetry = cmd.get_symmetry(sym_source or map_name)
        if symmetry:
            cmd.set("map_auto_expand_sym", 1)
            cmd.set_symmetry(self.center_name, *symmetry)

        cmd.set_key("pgup", self.contour_plus)
        cmd.set_key("pgdn", self.contour_minus)

        self.update()

    def load(self):
        cmd.load_callback(self, self.callback_name)

    def contour_plus(self, d=0.1):
        self.level += d
        print("Map level: " + str(self.level))
        self.update()

    def contour_minus(self):
        if self.level < 0.15:
            return
        self.contour_plus(-0.1)

    def update(self):
        self.center = cmd.get_position()
        cmd.alter_state(0, self.center_name, "(x, y, z) = p", space={'p': self.center})
        cmd.isomesh(self.name, self.map_name, self.level, self.center_name, carve=self.radius)

    def __call__(self):
        if self.name not in cmd.get_names('objects'):
            cmd.delete(self.callback_name)
            cmd.set_key("pgup", lambda: None)
            cmd.set_key("pgdn", lambda: None)
            return

        tmp = cmd.get_position()
        r = cpv.distance_sq(self.center, tmp)
        if (r > 0.3):  # increase this number if it is too slow
            self.update()

    def get_extent(self):
        tmp = cmd.get_position()
        return [[i - self.radius for i in tmp], [i + self.radius for i in tmp]]


def dynamic_mesh(map_name, level=1.0, radius=8, name='dynamic_mesh', sym_source=None):
    '''
DESCRIPTION

    Make 'dynamic' mesh from volumetric data such as electron density map.
    The mesh will dynamically follow the center of the view.
    Contour level of isomesh can be changed by PageDown and PageUp keys.

    NOTE: Crystallographic operations can be applied to the map.

USAGE

    dynamic_mesh map_name [, level [, radius [, name [, sym_source ]]]]

ARGUMENTS

    map_name = string: name of volumetric object(map) to display

    level = float: contour level of isomesh {default: 1.0}

    radius = float: radius of isomesh around the center of the view {default: 8}

    name = string: name of mesh object {default: dynamic_mesh}

    sym_source = string: name of object from which symmetry
                         information is derived {default: map_name}

EXAMPLE

    fetch 1hwk, async=0
    fetch 1hwk, 1hwk_map, type=2fofc, async=0
    dynamic_mesh 1hwk_map

SEE ALSO

    isomesh
    '''
    cmd.delete(name)
    callback = DynamicMesh(map_name, float(level), float(radius), name, sym_source)
    callback.load()

cmd.extend('dynamic_mesh', dynamic_mesh)
cmd.auto_arg[0]['dynamic_mesh'] = cmd.auto_arg[1]['isomesh']
