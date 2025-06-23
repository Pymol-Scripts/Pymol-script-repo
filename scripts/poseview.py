'''
http://pymolwiki.org/index.php/PoseView

(c) 2012 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

from __future__ import print_function

from pymol import cmd, CmdException


def poseview(ligand='organic inorganic', protein='polymer', width=0, height=0,
             filename='', exe='poseview', state=-1, quiet=1):
    '''
DESCRIPTION

    PoseView wrapper

    http://www.biosolveit.de/poseview/

USAGE

    poseview [ ligand [, protein [, width [, height [, exe [, state ]]]]]]

ARGUMENTS

    ligand = string: atom selection {default: organic inorganic}

    protein = string: atom selection {default: polymer}

    width = int: image width {default: viewport width}

    height = int: image height {default: viewport height}

    filename = string: output PNG file name {default: temporary}

    exe = string: path to executable {default: poseview}

SETUP

    1) Put poseview executable to PATH (e.g. /usr/bin/poseview)
    2) Set environment variable BIOSOLVE_LICENSE_FILE="/path/to/poseview.lic"
    '''
    import tempfile
    import subprocess
    import os
    import shutil

    width, height = int(width), int(height)
    state, quiet = int(state), int(quiet)

    if width == 0 or height == 0:
        viewport = cmd.get_viewport()
        if width != 0:
            height = int(viewport[0] / float(width) * viewport[1])
        elif height != 0:
            width = int(viewport[1] / float(height) * viewport[0])
        else:
            width, height = viewport

    exe = cmd.exp_path(exe)
    tempdir = tempfile.mkdtemp()

    try:
        ligand_filename = os.path.join(tempdir, 'ligand.sdf')
        protein_filename = os.path.join(tempdir, 'protein.pdb')

        if not filename:
            filename = os.path.join(tempdir, 'image.png')

        cmd.save(ligand_filename, ligand, state)
        cmd.save(protein_filename, protein, state)

        args = [exe, '-l', ligand_filename, '-p', protein_filename, '-t', '',
                '-o', filename, '-s', str(width), str(height)]

        if not quiet:
            print(' poseview: running...')

        process = subprocess.Popen(args,
                                   universal_newlines=True,
                                   stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
        stdout, _ = process.communicate()

        if not quiet:
            print(stdout)
            print(' poseview: done')

        if filename.endswith('.png'):
            cmd.load(filename)
        elif not quiet:
            print(' Warning: cannot load "%s" into PyMOL, unsupported file type' % (filename))

    except OSError:
        print(' Error: Cannot execute "%s", please provide full path to poseview executable' % (exe))
        raise CmdException
    finally:
        shutil.rmtree(tempdir)

cmd.extend('poseview', poseview)
