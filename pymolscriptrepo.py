'''
Support module for Pymol-script-repo

Adds repo binary directories to PATH if PYMOL_GIT_MOD is defined.
'''

import os

if 'PYMOL_GIT_MOD' in os.environ:
    import sys
    import platform

    repopath = [os.environ.get('PATH', '.')]

    ##
    # this is ugly, we need proper paths for each platform, but not for
    # every app an extra path!
    ##
    if sys.platform.startswith('linux'):
        if platform.architecture()[0] == '32bit':
            repopath.append(os.path.join(os.environ['PYMOL_GIT_MOD'], 'autodock_423', 'i86Linux2'))
        else:
            repopath.append(os.path.join(os.environ['PYMOL_GIT_MOD'], 'autodock_423', 'ia64Linux2'))
        repopath.append(os.path.join(os.environ['PYMOL_GIT_MOD'], 'autodock_vina', 'autodock_vina_1_1_2_linux_x86'))
    elif sys.platform.startswith('darwin'):
        version = platform.release().split('.')[0]
        if version in ('10', '9', '8'):
            repopath.append(os.path.join(os.environ['PYMOL_GIT_MOD'], 'autodock_423', 'universalDarwin' + version))
        repopath.append(os.path.join(os.environ['PYMOL_GIT_MOD'], 'autodock_vina', 'autodock_vina_1_1_2_mac'))
    elif sys.platform.startswith('win'):
        repopath.append(os.path.join(os.environ['PYMOL_GIT_MOD'], 'autodock_423', 'win32'))
        repopath.append(os.path.join(os.environ['PYMOL_GIT_MOD'], 'autodock_vina', 'autodock_vina_1_1_2_win32'))

    if len(repopath) > 1:
        os.environ['PATH'] = os.pathsep.join(repopath)


def which(*names):
    '''
    Return full path to executable or empty string if not found in PATH.
    '''
    if 'PATHEXT' in os.environ:
        pathext = [''] + os.environ['PATHEXT'].split(';')
        names = [n + ext for n in names for ext in pathext]
    path = os.environ['PATH'].split(os.pathsep)
    for n in names:
        for p in path:
            full = os.path.join(p, n)
            if os.path.isfile(full):
                return full
    return ''
