'''
http://pymolwiki.org/index.php/Extra_fit

(c) 2011 Thomas Holder

License: BSD-2-Clause
'''

__version__ = '1.0'

from pymol import cmd, CmdException


def extra_fit(selection='(all)', reference=None, method='align', zoom=1,
              quiet=0, _self=cmd, **kwargs):
    '''
DESCRIPTION

    Like "intra_fit", but for multiple objects instead of
    multiple states.

ARGUMENTS

    selection = string: atom selection of multiple objects {default: all}

    reference = string: reference object name {default: first object in selection}

    method = string: alignment method (command that takes "mobile" and "target"
    arguments, like "align", "super", "cealign" {default: align}

    ... extra arguments are passed to "method"

SEE ALSO

    alignto, cmd.util.mass_align, align_all.py from Robert Campbell
    '''
    zoom, quiet = int(zoom), int(quiet)
    models = cmd.get_object_list('(%s)' % selection)
    if reference is None:
        reference = models[0]
        models = models[1:]
    elif reference in models:
        models.remove(reference)
    if cmd.is_string(method):
        if method in cmd.keyword:
            method = cmd.keyword[method][0]
        else:
            print('Unknown method ' + str(method))
            raise CmdException
    for model in models:
        x = method(mobile='(%s) and model %s' % (selection, model),
                   target='(%s) and model %s' % (selection, reference), **kwargs)
        if not quiet:
            if cmd.is_sequence(x):
                print('%-20s RMS = %8.3f (%d atoms)' % (model, x[0], x[1]))
            elif isinstance(x, float):
                print('%-20s RMS = %8.3f' % (model, x))
            elif _self.is_dict(x) and 'RMSD' in x:
                natoms = x.get('alignment_length', 0)
                suffix = (' (%s atoms)' % natoms) if natoms else ''
                print('%-20s RMS = %8.3f' % (model, x['RMSD']) + suffix)
            else:
                print('%-20s' % (model,))

    if zoom:
        cmd.zoom(selection)

cmd.extend('extra_fit', extra_fit)

# vi:expandtab:smarttab
