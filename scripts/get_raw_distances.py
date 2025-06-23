'''
http://pymolwiki.org/index.php/get_raw_distances

(c) 2012 Takanori Nakane and Thomas Holder

License: BSD-2-Clause
'''

from pymol import cmd, CmdException


def get_raw_distances(names='', state=1, selection='all', quiet=1):
    '''
DESCRIPTION

    Get the list of pair items from distance objects. Each list item is a
    tuple of (index1, index2, distance).

    Based on a script from Takanori Nakane, posted on pymol-users mailing list.
    http://www.mail-archive.com/pymol-users@lists.sourceforge.net/msg10143.html

ARGUMENTS

    names = string: names of distance objects (no wildcards!) {default: all
    measurement objects}

    state = integer: object state {default: 1}

    selection = string: atom selection {default: all}

SEE ALSO

    select_distances, cmd.find_pairs, cmd.get_raw_alignment
    '''
    from chempy import cpv

    state, quiet = int(state), int(quiet)
    if state < 1:
        state = cmd.get_state()

    valid_names = cmd.get_names_of_type('object:measurement')
    if names == '':
        names = ' '.join(valid_names)
    else:
        for name in names.split():
            if name not in valid_names:
                print(' Error: no such distance object: ' + name)
                raise CmdException

    raw_objects = cmd.get_session(names, 1, 1, 0, 0)['names']

    xyz2idx = {}
    cmd.iterate_state(state, selection, 'xyz2idx[x,y,z] = (model,index)',
                      space=locals())

    r = []
    for obj in raw_objects:
        try:
            points = obj[5][2][state - 1][1]
            if points is None:
                raise ValueError
        except (KeyError, ValueError):
            continue
        for i in range(0, len(points), 6):
            xyz1 = tuple(points[i:i + 3])
            xyz2 = tuple(points[i + 3:i + 6])
            try:
                r.append((xyz2idx[xyz1], xyz2idx[xyz2], cpv.distance(xyz1, xyz2)))
                if not quiet:
                    print(' get_raw_distances: ' + str(r[-1]))
            except KeyError:
                if quiet < 0:
                    print(' Debug: no index for %s %s' % (xyz1, xyz2))
    return r


def select_distances(names='', name='sele', state=1, selection='all', cutoff=-1, quiet=1):
    '''
DESCRIPTION

    Turns a distance object into a named atom selection.

ARGUMENTS

    names = string: names of distance objects (no wildcards!) {default: all
    measurement objects}

    name = a unique name for the selection {default: sele}

SEE ALSO

    get_raw_distances
    '''
    state, cutoff, quiet = int(state), float(cutoff), int(quiet)

    sele_dict = {}
    distances = get_raw_distances(names, state, selection)
    for idx1, idx2, dist in distances:
        if cutoff <= 0.0 or dist <= cutoff:
            sele_dict.setdefault(idx1[0], set()).add(idx1[1])
            sele_dict.setdefault(idx2[0], set()).add(idx2[1])

    cmd.select(name, 'none')
    tmp_name = cmd.get_unused_name('_')

    r = 0
    for model in sele_dict:
        cmd.select_list(tmp_name, model, list(sele_dict[model]), mode='index')
        r = cmd.select(name, tmp_name, merge=1)
        cmd.delete(tmp_name)

    if not quiet:
        print(' Selector: selection "%s" defined with %d atoms.' % (name, r))
    return r

cmd.extend('get_raw_distances', get_raw_distances)
cmd.extend('select_distances', select_distances)

_auto_arg0_distances = [
    lambda: cmd.Shortcut(cmd.get_names_of_type('object:measurement')),
    'distance object', '']

cmd.auto_arg[0].update([
    ('get_raw_distances', _auto_arg0_distances),
    ('select_distances', _auto_arg0_distances),
])

# vi: ts=4:sw=4:smarttab:expandtab
