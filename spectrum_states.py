'''
http://pymolwiki.org/index.php/spectrum_states

(c) 2011 Takanori Nakane and Thomas Holder

License: BSD-2-Clause
'''

from pymol import cmd, CmdException

def spectrum_states(selection='all', representations='cartoon ribbon',
        color_list='blue cyan green yellow orange red',
        first=1, last=0, quiet=1):
    '''
DESCRIPTION

    Color each state in a multi-state object different.

USAGE

    spectrum_states [ selection [, representations [, color_list [, first [, last ]]]]]

ARGUMENTS

    selection = string: object names (works with complete objects only)
    {default: all}

    representations = string: space separated list of representations
    {default: cartoon ribbon}

    color_list = string: space separated list of colors {default: blue cyan
    green yellow orange red}

SEE ALSO

    spectrum, spectrumany
    '''
    from math import floor, ceil

    first, last, quiet = int(first), int(last), int(quiet)
    colors = color_list.split()
    if len(colors) < 2:
        print ' Error: please provide at least 2 colors'
        raise CmdException

    colvec = [cmd.get_color_tuple(i) for i in colors]

    # filter for valid <repr>_color settings
    settings = []
    for r in representations.split():
        if r[-1] == 's':
            r = r[:-1]
        s = r + '_color'
        if s in cmd.setting.name_list:
            settings.append(s)
        elif not quiet:
            print ' Warning: no such setting:', s

    # object names only
    selection = ' '.join(cmd.get_object_list('(' + selection + ')'))
    if cmd.count_atoms(selection) == 0:
        print ' Error: empty selection'
        raise CmdException

    if last < 1:
        last = cmd.count_states(selection)

    val_range = int(last - first + 1)
    if val_range < 2:
        print ' Error: no spectrum possible, need more than 1 state'
        raise CmdException

    for i in range(val_range):
        p = float(i) / (val_range - 1) * (len(colvec) - 1)
        p0, p1 = int(floor(p)), int(ceil(p))
        ii = (p - p0)
        col_list = [colvec[p1][j] * ii + colvec[p0][j] * (1.0 - ii) for j in range(3)]
        col_name = '0x%02x%02x%02x' % (col_list[0] * 255, col_list[1] * 255, col_list[2] * 255)
        for s in settings:
            cmd.set(s, col_name, selection, state=i+first)

cmd.extend('spectrum_states', spectrum_states)

# tab-completion of arguments
cmd.auto_arg[0]['spectrum_states'] = cmd.auto_arg[0]['disable']
cmd.auto_arg[1]['spectrum_states'] = [ cmd.auto_arg[0]['show'][0], 'representation', ' ' ]
cmd.auto_arg[2]['spectrum_states'] = [ cmd.auto_arg[0]['color'][0], 'color', ' ' ]

# vi:expandtab:smarttab
