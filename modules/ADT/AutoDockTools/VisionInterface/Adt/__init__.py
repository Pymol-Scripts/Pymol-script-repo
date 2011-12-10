########################################################################
#
# Date: Jan 2006 Authors: Guillaume Vareille, Michel Sanner
#
#    vareille@scripps.edu
#    sanner@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Guillaume Vareille, Michel Sanner and TSRI
#
#    Vision Library Loader
#
#########################################################################
#
# /home/vareille/.mgltools/1.6.0/Vision/UserLibs/MyDefaultLib/__init__.py
# Vision will generate this file automatically if it can't find it
#

from os import sep, path
from Vision.VPE import NodeLibrary
from Vision.UserLibBuild import userLibBuild, addDirToSysPath, addTypes

dependents = {} # {'scipy':'0.6.0',} the numbers indicate the highest tested version of the needed packages
libraryColor = '#4444FF'

addDirToSysPath(path.dirname(__file__)+sep+'..')
fileSplit = __file__.split(sep)
if fileSplit[-1] == '__init__.pyc' or fileSplit[-1] == '__init__.py':
    libInstanceName = fileSplit[-2]
else:
    libInstanceName = path.splitext(fileSplit[-1])[0]
try:
    from Vision import ed
except:
    ed = None
if ed is not None and ed.libraries.has_key(libInstanceName):
    locals()[libInstanceName] = ed.libraries[libInstanceName]
else:
    locals()[libInstanceName] = NodeLibrary(libInstanceName, libraryColor, mode='readWrite')
success = userLibBuild(eval(libInstanceName), __file__, dependents=dependents)
if success is False:
    locals().pop(libInstanceName)
elif path.isfile(path.dirname(__file__)+sep+'libTypes.py'):
    addTypes(locals()[libInstanceName], libInstanceName + '.libTypes')

