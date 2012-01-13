########################################################################
#
# Date: 2000 Authors: Guillaume Vareille, Michel Sanner
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
#########################################################################
#
# $Header: /opt/cvs/python/packages/share1.5/mglutil/gui/__init__.py,v 1.13 2010/04/14 22:39:42 sanner Exp $
#
# $Id: __init__.py,v 1.13 2010/04/14 22:39:42 sanner Exp $
#

import os
import sys
import warnings

mglutilrcText = """########################################################################
#
# Date: Decembre 2006 Authors: Guillaume Vareille, Michel Sanner
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
#    mglutil Resource File
#
#########################################################################
# To customize mglutil, you can modify the _mglutilrc file:
# unix: ~/.mgltools/[version number]/mglutil/_mglutilrc
# windows: \Documents and Settings\(user name)\.mgltools\(version numer)\mglutil\_mglutilrc
# DejaVu will generate this file automatically if it can't find it
# Do not modify the original source file
##################################################################

import os

widgetsOnBackWindowsCanGrabFocus = (os.name != 'nt') # True, False

"""
        
def ensureMglutilResourceFile():
    """verify or generate _mglutilrc file
"""
    #print "ensureMglutilResourceFile"

    from mglutil.util.packageFilePath import getResourceFolderWithVersion
    rcFolder = getResourceFolderWithVersion()
    if rcFolder is None:
        return None
    rcFolder += os.sep + 'mglutil'
    if not os.path.isdir(rcFolder):
        try:
            os.mkdir(rcFolder)
        except:
            txt = "Cannot create the Resource Folder %s" %rcFolder
            warnings.warn(txt)
            return None

    rcFile = rcFolder + os.sep + '_mglutilrc'
    if os.path.isfile(rcFile) is False:
        try:
            f = open(rcFile, "w")
            global mglutilrcText
            map( lambda x, f=f: f.write(x), mglutilrcText )
            f.close()
        except:
            txt = "can not create _mglutilrc"
            warnings.warn(txt)
            return None

    return rcFile


# after this we can access variables in _mglutilrc with
# from mglutil.gui import widgetsOnBackWindowsCanGrabFocus
rcFile = ensureMglutilResourceFile()
if rcFile is None:
    exec( mglutilrcText )
else:
    execfile( rcFile )

    if os.name != 'nt': #sys.platform != 'win32':
        # we create a symbolic link to the shell script that launch python
        # this is link is used by the self-running saved vision network to
        # find mgltools
        # (we create this each time we run mglutils,
        # this way the last running vision will be used)
        mgltoolsDir = os.path.split(os.path.split(os.path.split(os.path.split(
            __file__)[0])[0])[0])[0]
        pythonshFile = os.path.join(mgltoolsDir, 'bin', 'pythonsh')
        if os.path.isfile(pythonshFile) is False:
            pythonshFile = os.path.realpath(sys.executable)+'sh'
        try:
            os.system('ln -f -s ' + pythonshFile + ' ~' + os.sep +'.mgltools'+ os.sep +'pythonsh')
        except:
            pass
    

