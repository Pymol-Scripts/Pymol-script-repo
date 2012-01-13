#############################################################################
#
# Author: Ruth HUEY, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Objects/__init__.py,v 1.2 2009/06/17 00:03:24 vareille Exp $
#
# $Id: __init__.py,v 1.2 2009/06/17 00:03:24 vareille Exp $

# create hostDict with hostMacros accessible by anyone
from AutoDockTools.adthosts import hostMacros
from AutoDockTools.autodockHosts import AutoDockHosts

hostDict = AutoDockHosts(hostMacros)

import socket
h= socket.gethostname()
hostDict[h]=hostDict['localhost']
hostDict[h]['host']=h
del hostDict['localhost']

# try to extend that dictionary with user specific hostMacros

# first try to find a adthost file in current directory
import os,sys
if os.path.isfile('./adthosts.py'):
    execfile('./adthosts.py')
    if globals().has_key('hostMacros'):
        hostDict.update(hostMacros)
elif os.name != 'nt':#sys.platform!='win32':
    # try to find the user's home directory
    import posix
    if 'HOME' in posix.environ.keys():
	try:
            execfile(os.path.join(posix.environ['HOME'],'adthosts.py'))
            if globals().has_key('hostMacros'):
                hostDict.update(hostMacros)
	except:
	    pass
    
                 
