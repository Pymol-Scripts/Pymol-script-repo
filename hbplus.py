'''
http://pymolwiki.org/index.php/hbplus

(c) 2012 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

from __future__ import print_function

from pymol import cmd, CmdException


def hbplus(selection='all', exe='hbplus', prefix='hb_', state=-1, quiet=1):
    '''
DESCRIPTION

    HBPLUS wrapper
    '''
    import subprocess
    import tempfile
    import os
    import shutil

    state, quiet = int(state), int(quiet)

    tempdir = tempfile.mkdtemp()
    pdbfile = os.path.join(tempdir, 'foo.pdb')
    hb2file = os.path.join(tempdir, 'foo.hb2')

    try:
        for model in cmd.get_object_list('(' + selection + ')'):
            cmd.save(pdbfile, 'model %s and (%s)' % (model, selection), state)
            process = subprocess.Popen([exe, pdbfile], cwd=tempdir,
                                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            stdout, _ = process.communicate()

            if not quiet:
                print(stdout)

            for line in open(hb2file):
                if line[5] != '-' or line[19] != '-':
                    continue
                chain1, resi1, name1 = line[0], line[1:5], line[10:14]
                chain2, resi2, name2 = line[14], line[15:19], line[24:28]
                sele1 = '/%s//%s/%s/%s' % (model, chain1, resi1.lstrip('0'), name1.rstrip())
                sele2 = '/%s//%s/%s/%s' % (model, chain2, resi2.lstrip('0'), name2.rstrip())
                cat = line[33:35]
                cmd.distance(prefix + cat, sele1, sele2)

    except OSError:
        print('Error: Cannot execute exe=' + exe)
        raise CmdException
    finally:
        shutil.rmtree(tempdir)

cmd.extend('hbplus', hbplus)
