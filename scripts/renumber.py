'''
http://pymolwiki.org/index.php/Renumber

(c) 2012 Thomas Holder

License: BSD-2-Clause
'''

from pymol import cmd, CmdException


def renumber(selection='all', start=1, startsele=None, quiet=1):
    '''
DESCRIPTION

    Set residue numbering (resi) based on connectivity.

ARGUMENTS

    selection = string: atom selection to renumber {default: all}

    start = integer: counting start {default: 1}

    startsele = string: residue to start counting from {default: first in
    selection}
    '''
    start, quiet = int(start), int(quiet)
    model = cmd.get_model(selection)
    cmd.iterate(selection, 'next(atom_it).model = model',
            space={'atom_it': iter(model.atom), 'next': next})
    if startsele is not None:
        startidx = cmd.index('first (' + startsele + ')')[0]
        for atom in model.atom:
            if (atom.model, atom.index) == startidx:
                startatom = atom
                break
        else:
            print(' Error: startsele not in selection')
            raise CmdException
    else:
        startatom = model.atom[0]
    for atom in model.atom:
        atom.adjacent = []
        atom.visited = False
    for bond in model.bond:
        atoms = [model.atom[i] for i in bond.index]
        atoms[0].adjacent.append(atoms[1])
        atoms[1].adjacent.append(atoms[0])
    minmax = [start, start]

    def traverse(atom, resi):
        atom.resi = resi
        atom.visited = True
        for other in atom.adjacent:
            if other.visited:
                continue
            if (atom.name, other.name) in [('C', 'N'), ("O3'", 'P')]:
                minmax[1] = resi + 1
                traverse(other, resi + 1)
            elif (atom.name, other.name) in [('N', 'C'), ('P', "O3'")]:
                minmax[0] = resi - 1
                traverse(other, resi - 1)
            elif (atom.name, other.name) not in [('SG', 'SG')]:
                traverse(other, resi)
    traverse(startatom, start)
    cmd.alter(selection, 'resi = next(atom_it).resi',
            space={'atom_it': iter(model.atom), 'next': next})
    if not quiet:
        print(' Renumber: range (%d to %d)' % tuple(minmax))

cmd.extend('renumber', renumber)

# vi:expandtab:smarttab
