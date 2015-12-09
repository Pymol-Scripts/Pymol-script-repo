'''
http://pymolwiki.org/index.php/Save_pdb_with_anisou

(c) 2012 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

from __future__ import print_function

from pymol import cmd, CmdException


def save_pdb_with_anisou(filename, selection='(all)', state=1, quiet=1):
    '''
DESCRIPTION

     Save in PDB format including ANISOU records.

SEE ALSO

    save
    '''
    state, quiet = int(state), int(quiet)

    pdbstr = cmd.get_pdbstr(selection, state)
    atom_it = iter(cmd.get_model(selection, state).atom)

    def mergeaniso():
        for line in pdbstr.splitlines(True):
            yield line
            if line[:6] in ['ATOM  ', 'HETATM']:
                u_str = ''.join('%7.0f' % (u * 1e4) for u in atom_it.next().u_aniso)
                yield 'ANISOU' + line[6:28] + u_str + line[70:]

    pdbstr = ''.join(mergeaniso())

    filename = cmd.exp_path(filename)
    f = open(filename, 'w')
    f.write(pdbstr)
    f.close()

    if not quiet:
        print(' Save with ANISOU: wrote "%s"' % (filename))


cmd.extend('save_pdb_with_anisou', save_pdb_with_anisou)

cmd.auto_arg[1]['save_pdb_with_anisou'] = cmd.auto_arg[1]['save']

# vi:expandtab:smarttab
