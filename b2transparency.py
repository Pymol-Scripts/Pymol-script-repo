'''
http://pymolwiki.org/index.php/B2transparency

(c) 2012 Thomas Holder

License: BSD-2-Clause
'''

from pymol import cmd, CmdException


def b2transparency(selection='all', setting='transparency', minimum=None,
                   maximum=None, var='b', quiet=1):
    '''
DESCRIPTION

    Set surface (or other) transparency for each atom scaled by b-factor.

    Does not work for all, but for some transparency settings (for example
    transparency, sphere_transparency)

ARGUMENTS

    selection = string: atom selection {default: all}

    setting = string: setting name {default: transparency}

    minimum = float: b-factor range minimum {default: automatic}

    maximum = float: b-factor range maximum {default: automatic}

    var = string: numeric atomic property like b or q {default: b}

SEE ALSO

    spectrum, cartoon putty
    '''
    quiet = int(quiet)

    if minimum is None or maximum is None:
        b_list = []
        if not cmd.iterate(selection, 'b_list.append(%s)' % var, space=locals()):
            if not quiet:
                print(' b2transparency: empty selection')
            return

        if minimum is None:
            minimum = min(b_list)
        if maximum is None:
            maximum = max(b_list)

    minimum, maximum = float(minimum), float(maximum)
    if not quiet:
        print(' b2transparency: range (%.5f to %.5f)' % (minimum, maximum))

    cmd.iterate(selection, """cmd.set('%s', min(max((%s - %f) / %f, 0), 0.9),
            '(%%s`%%d)' %% (model, index))
            """ % (setting, var, minimum, maximum - minimum), quiet=quiet)

cmd.extend('b2transparency', b2transparency)

cmd.auto_arg[0]['b2transparency'] = cmd.auto_arg[0]['zoom']
cmd.auto_arg[1]['b2transparency'] = cmd.auto_arg[0]['set']
