'''
Square and Tetrahedra representations

(c) 2013 Thomas Holder

License: BSD-2-Clause
'''

from __future__ import print_function
from pymol import cmd, cgo
from chempy import cpv


def cgo_cube(x, y, z, r):
    r *= 3 ** -.5
    return [
        cgo.BEGIN, cgo.TRIANGLE_STRIP,
        cgo.NORMAL, 0., 0., 1.,
        cgo.VERTEX, x + r, y + r, z + r,
        cgo.VERTEX, x + r, y - r, z + r,
        cgo.VERTEX, x - r, y + r, z + r,
        cgo.VERTEX, x - r, y - r, z + r,
        cgo.END,
        cgo.BEGIN, cgo.TRIANGLE_STRIP,
        cgo.NORMAL, 1., 0., 0.,
        cgo.VERTEX, x + r, y - r, z - r,
        cgo.VERTEX, x + r, y + r, z - r,
        cgo.VERTEX, x + r, y - r, z + r,
        cgo.VERTEX, x + r, y + r, z + r,
        cgo.END,
        cgo.BEGIN, cgo.TRIANGLE_STRIP,
        cgo.NORMAL, 0., 1., 0.,
        cgo.VERTEX, x + r, y + r, z - r,
        cgo.VERTEX, x - r, y + r, z - r,
        cgo.VERTEX, x + r, y + r, z + r,
        cgo.VERTEX, x - r, y + r, z + r,
        cgo.END,
        cgo.BEGIN, cgo.TRIANGLE_STRIP,
        cgo.NORMAL, 0., 0., -1.,
        cgo.VERTEX, x - r, y - r, z - r,
        cgo.VERTEX, x - r, y + r, z - r,
        cgo.VERTEX, x + r, y - r, z - r,
        cgo.VERTEX, x + r, y + r, z - r,
        cgo.END,
        cgo.BEGIN, cgo.TRIANGLE_STRIP,
        cgo.NORMAL, -1., 0., 0.,
        cgo.VERTEX, x - r, y + r, z + r,
        cgo.VERTEX, x - r, y - r, z + r,
        cgo.VERTEX, x - r, y + r, z - r,
        cgo.VERTEX, x - r, y - r, z - r,
        cgo.END,
        cgo.BEGIN, cgo.TRIANGLE_STRIP,
        cgo.NORMAL, 0., -1., 0.,
        cgo.VERTEX, x - r, y - r, z + r,
        cgo.VERTEX, x + r, y - r, z + r,
        cgo.VERTEX, x - r, y - r, z - r,
        cgo.VERTEX, x + r, y - r, z - r,
        cgo.END,
    ]


def cgo_tetrahedron(x, y, z, r):
    vertices = [cpv.add((x, y, z), cpv.scale(v, r)) for v in [
        [0., 1., 0.],
        [0.0, -0.3338068592337708, 0.9426414910921784],
        [0.8163514779470693, -0.3338068592337708, -0.471320745546089],
        [-0.816351477947069, -0.3338068592337708, -0.4713207455460897]
    ]]
    return [
        cgo.BEGIN, cgo.TRIANGLES,
        cgo.NORMAL, 0.8165448970931916, 0.33317549135767066, 0.4714324161421696,
        cgo.VERTEX] + vertices[0] + [
        cgo.VERTEX] + vertices[1] + [
        cgo.VERTEX] + vertices[2] + [
        cgo.NORMAL, 0., 0.3331754913576707, -0.9428648322843389,
        cgo.VERTEX] + vertices[0] + [
        cgo.VERTEX] + vertices[2] + [
        cgo.VERTEX] + vertices[3] + [
        cgo.NORMAL, -0.8165448970931919, 0.3331754913576705, 0.4714324161421693,
        cgo.VERTEX] + vertices[0] + [
        cgo.VERTEX] + vertices[3] + [
        cgo.VERTEX] + vertices[1] + [
        cgo.NORMAL, 0., -1., 0.,
        cgo.VERTEX] + vertices[1] + [
        cgo.VERTEX] + vertices[2] + [
        cgo.VERTEX] + vertices[3] + [
        cgo.END,
    ]


def cubes(selection='all', name='', state=0, scale=0.5, atomcolors=1, _func=cgo_cube):
    '''
DESCRIPTION

    Create a cube representation CGO for all atoms in selection.

ARGUMENTS

    selection = string: atom selection {default: all}

    name = string: name of CGO object to create

    state = int: object state {default: 0 = all states}

    scale = float: scaling factor. If scale=1.0, the corners of the cube will
    be on the VDW surface of the atom {default: 0.5}

    atomcolors = 0/1: use atom colors (cannot be changed), otherwise
    apply one color to the object (can be changed with color command)
    {default: 1}

SEE ALSO

    tetrahedra
    '''
    if not name:
        name = cmd.get_unused_name('cubes')
    state, scale, atomcolors = int(state), float(scale), int(atomcolors)
    if state < 0:
        state = cmd.get_setting_int('state')
    states = [state] if state else list(range(1,
                                         cmd.count_states(selection) + 1))

    def callback(x, y, z, vdw, color):
        if atomcolors:
            obj.append(cgo.COLOR)
            obj.extend(cmd.get_color_tuple(color))
        obj.extend(_func(x, y, z, vdw * scale))
    space = {'xcb': callback}
    for state in states:
        obj = []
        cmd.iterate_state(state, selection,
                          'xcb(x, y, z, vdw, color)', space=space)
        cmd.load_cgo(obj, name, state)
    if not atomcolors:
        cmd.color('auto', name)


def tetrahedra(selection='all', name='', state=0, scale=0.5, atomcolors=1):
    '''
DESCRIPTION

    Create a tetrahedra representation CGO for all atoms in selection.

SEE ALSO

    cubes
    '''
    if not name:
        name = cmd.get_unused_name('tetrahedra')
    return cubes(selection, name, state, scale, atomcolors, cgo_tetrahedron)

cmd.extend('cubes', cubes)
cmd.extend('tetrahedra', tetrahedra)
