'''
See more here: http://www.pymolwiki.org/index.php/ex
###############################################
#  File:          ex.py
#  Author:        Troels E. Linnet
#  Creation Date: 11/12/11
#
#  Notes: Easy opening of example files
#
###############################################
import ex
ex rotkit_1.pml
'''

from __future__ import print_function
from pymol import cmd
import os
if 'PYMOL_GIT_MOD' in os.environ:
    os.environ['PYMOL_GIT_EX'] = os.path.join(os.path.split(os.environ['PYMOL_GIT_MOD'])[0], "examples")


def ex(filename=''):
    if os.path.splitext(filename)[1] == '':
        filename = os.path.splitext(filename)[0] + '.pml'
        print(("filename is: %s" % filename))
    if 'PYMOL_GIT_EX' in os.environ:
        expath = os.environ['PYMOL_GIT_EX']
        cmdrun = os.path.join(expath, filename)
        cmd.do("@%s" % cmdrun)
    else:
        cmdrun = filename
        cmd.do("@%s" % cmdrun)
    return(cmdrun)
cmd.extend("ex", ex)
