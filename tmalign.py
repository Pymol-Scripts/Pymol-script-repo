'''
(c) 2011-2012 Thomas Holder, MPI for Developmental Biology
'''

from __future__ import print_function

__author__ = 'Thomas Holder'
__version__ = '1.1'
__license__ = 'BSD-2-Clause'

from pymol import cmd, CmdException


def save_pdb_without_ter(filename, selection, **kwargs):
    '''
DESCRIPTION

    Save PDB file without TER records. External applications like TMalign and
    DynDom stop reading PDB files at TER records, which might be undesired in
    case of missing loops.
    '''
    v = cmd.get_setting_boolean('pdb_use_ter_records')
    if v:
        cmd.set('pdb_use_ter_records', 0)
    cmd.save(filename, selection, **kwargs)
    if v:
        cmd.set('pdb_use_ter_records')


def alignwithanymethod(mobile, target, methods='align super cealign tmalign',
                       async_=1, quiet=1, **kwargs):
    '''
DESCRIPTION

    Align copies of mobile to target with several alignment methods

ARGUMENTS

    mobile = string: atom selection

    target = string: atom selection

    methods = string: space separated list of PyMOL commands which take
    arguments "mobile" and "target" (in any order) {default: align super
    cealign tmalign}
    '''
    import threading
    import time
    methods = methods.split()
    async_, quiet = int(kwargs.pop('async', async_)), int(quiet)
    mobile_obj = cmd.get_object_list('first (' + mobile + ')')[0]

    def myalign(method):
        newmobile = cmd.get_unused_name(mobile_obj + '_' + method)
        cmd.create(newmobile, mobile_obj)
        start = time.time()
        cmd.do('%s mobile=%s in %s, target=%s' % (method, newmobile, mobile, target))
        if not quiet:
            print('Finished: %s (%.2f sec)' % (method, time.time() - start))

    for method in methods:
        if async_:
            t = threading.Thread(target=myalign, args=(method,))
            t.setDaemon(1)
            t.start()
        else:
            myalign(method)


def tmalign(mobile, target, args='', exe='TMalign', ter=0, transform=1, object=None, quiet=0):
    '''
DESCRIPTION

    TMalign wrapper

    Reference: Y. Zhang and J. Skolnick, Nucl. Acids Res. 2005 33, 2302-9
    http://zhanglab.ccmb.med.umich.edu/TM-align/

USAGE

    tmalign mobile, target [, args [, exe ]]

ARGUMENTS

    mobile, target = string: atom selections

    args = string: Extra arguments like -d0 5 -L 100

    exe = string: Path to TMalign executable {default: TMalign}

    ter = 0/1: If ter=0, then ignore chain breaks because TMalign will stop
    at first TER record {default: 0}

SEE ALSO

    tmscore, mmalign
    '''
    import subprocess
    import tempfile
    import os
    import re

    ter, quiet = int(ter), int(quiet)

    mobile_filename = tempfile.mktemp('.pdb', 'mobile')
    target_filename = tempfile.mktemp('.pdb', 'target')
    matrix_filename = tempfile.mktemp('.txt', 'matrix')
    mobile_ca_sele = '(%s) and (not hetatm) and name CA and alt +A' % (mobile)
    target_ca_sele = '(%s) and (not hetatm) and name CA and alt +A' % (target)

    if ter:
        save = cmd.save
    else:
        save = save_pdb_without_ter
    save(mobile_filename, mobile_ca_sele)
    save(target_filename, target_ca_sele)

    exe = cmd.exp_path(exe)
    args = [exe, mobile_filename, target_filename, '-m', matrix_filename] + args.split()

    try:
        process = subprocess.Popen(args, stdout=subprocess.PIPE,
                universal_newlines=True)
        lines = process.stdout.readlines()
    except OSError:
        print('Cannot execute "%s", please provide full path to TMscore or TMalign executable' % (exe))
        raise CmdException
    finally:
        os.remove(mobile_filename)
        os.remove(target_filename)

    # TMalign >= 2012/04/17
    if os.path.exists(matrix_filename):
        lines += open(matrix_filename).readlines()
        os.remove(matrix_filename)

    r = None
    re_score = re.compile(r'TM-score\s*=\s*(\d*\.\d*)')
    rowcount = 0
    matrix = []
    line_it = iter(lines)
    headercheck = False
    alignment = []
    for line in line_it:
        if 4 >= rowcount > 0:
            if rowcount >= 2:
                a = list(map(float, line.split()))
                matrix.extend(a[2:5])
                matrix.append(a[1])
            rowcount += 1
        elif not headercheck and line.startswith(' * '):
            a = line.split(None, 2)
            if len(a) == 3:
                headercheck = a[1]
        elif (' rotation matrix') in line.lower():
            rowcount = 1
        elif line.startswith('(":" denotes'):
            alignment = [next(line_it).rstrip() for i in range(3)]
        else:
            match = re_score.search(line)
            if match is not None:
                r = float(match.group(1))
        if not quiet:
            print(line.rstrip())

    if not quiet:
        for i in range(0, len(alignment[0]) - 1, 78):
            for line in alignment:
                print(line[i:i + 78])
            print('')

    assert len(matrix) == 3 * 4
    matrix.extend([0.0, 0.0, 0.0, 1.0])

    if int(transform):
        cmd.transform_selection('byobject (%s)' % (mobile), matrix, homogenous=1)

    # alignment object
    if object is not None:
        mobile_idx, target_idx = [], []
        space = {'mobile_idx': mobile_idx, 'target_idx': target_idx}
        cmd.iterate(mobile_ca_sele, 'mobile_idx.append("%s`%d" % (model, index))', space=space)
        cmd.iterate(target_ca_sele, 'target_idx.append("%s`%d" % (model, index))', space=space)
        for i, aa in enumerate(alignment[0]):
            if aa == '-':
                mobile_idx.insert(i, None)
        for i, aa in enumerate(alignment[2]):
            if aa == '-':
                target_idx.insert(i, None)
        if (len(mobile_idx) == len(target_idx) == len(alignment[2])):
            cmd.rms_cur(
                ' '.join(idx for (idx, m) in zip(mobile_idx, alignment[1]) if m in ':.'),
                ' '.join(idx for (idx, m) in zip(target_idx, alignment[1]) if m in ':.'),
                cycles=0, matchmaker=4, object=object)
        else:
            print('Could not load alignment object')

    if not quiet and r is not None:
        print('Found in output TM-score = %.4f' % (r))

    return r


def tmscore(mobile, target, args='', exe='TMscore', quiet=0, **kwargs):
    '''
DESCRIPTION

    TMscore wrapper

    Reference: Yang Zhang and Jeffrey Skolnick, Proteins 2004 57: 702-710
    http://zhanglab.ccmb.med.umich.edu/TM-score/

ARGUMENTS

    mobile, target = string: atom selections

    args = string: Extra arguments like -d 5

    exe = string: Path to TMscore executable {default: TMscore}

    ter = 0/1: If ter=0, then ignore chain breaks because TMscore will stop
    at first TER record {default: 0}

SEE ALSO

    tmalign, mmalign
    '''
    kwargs.pop('_self', None)
    return tmalign(mobile, target, args, exe, quiet=quiet, **kwargs)


def mmalign(mobile, target, args='', exe='MMalign', ter=0, transform=1, quiet=0):
    '''
DESCRIPTION

    MMalign wrapper

    Reference: S. Mukherjee and Y. Zhang, Nucleic Acids Research 2009; 37: e83
    http://zhanglab.ccmb.med.umich.edu/MM-align/

SEE ALSO

    tmalign, tmscore
    '''
    return tmalign(mobile, target, args, exe, ter, transform, quiet=quiet)

# pymol commands
cmd.extend('alignwithanymethod', alignwithanymethod)
cmd.extend('tmalign', tmalign)
cmd.extend('tmscore', tmscore)
cmd.extend('mmalign', tmalign)

# autocompletion
cmd.auto_arg[0].update({
    'tmalign': cmd.auto_arg[0]['align'],
    'tmscore': cmd.auto_arg[0]['align'],
    'mmalign': cmd.auto_arg[0]['align'],
})
cmd.auto_arg[1].update({
    'tmalign': cmd.auto_arg[1]['align'],
    'tmscore': cmd.auto_arg[1]['align'],
    'mmalign': cmd.auto_arg[1]['align'],
})
