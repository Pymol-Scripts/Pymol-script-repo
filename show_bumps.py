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
    cmd.sculpt_activate(name)
    cmd.show_as('cgo', name)
    cmd.set('sculpt_vdw_vis_mode', 1, name)
    cmd.set('sculpt_field_mask', 0x020) # cSculptVDW
    strain = cmd.sculpt_iterate(name, cycles=1)
    if not int(quiet):
        print 'VDW Strain:', strain
    return strain

cmd.extend('show_bumps', show_bumps)
