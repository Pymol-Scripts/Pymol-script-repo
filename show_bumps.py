'''
http://pymolwiki.org/index.php/show_bumps

(c) 2011 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

from pymol import cmd


def show_bumps(selection='(all)', name='bump_check', quiet=1):
    '''
DESCRIPTION

    Visualize VDW clashes

ARGUMENTS

    selection = string: atom selection {default: all}

    name = string: name of CGO object to create {default: bump_check}
    '''
    cmd.delete(name)
    cmd.create(name, selection, zoom=0)
    cmd.set('sculpt_vdw_vis_mode', 1, name)
    cmd.set('sculpt_field_mask', 0x020)  # cSculptVDW
    for state in range(1, 1 + cmd.count_states('%' + name)):
        cmd.sculpt_activate(name, state)
        strain = cmd.sculpt_iterate(name, state, cycles=0)
        if not int(quiet):
            print('VDW Strain in state %d: %f' % (state, strain))
    cmd.show_as('cgo', name)

cmd.extend('show_bumps', show_bumps)
