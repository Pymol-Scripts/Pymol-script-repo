#############################################################################
#
# Author: Ruth HUEY, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2008
#
#############################################################################


# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/autostart4Commands.py,v 1.3 2009/06/17 00:03:24 vareille Exp $
#
# $Id: autostart4Commands.py,v 1.3 2009/06/17 00:03:24 vareille Exp $ 
#
#
#
#
#
#

"""
This Module facilitates starting autogrid and autodock jobs and managing them

"""
from ViewerFramework.VFCommand import CommandGUI
from AutoDockTools.autostartCommands import ADKill, ADProcessManager,\
AutoStarter, AutoGridStarter, AutoDockStarter, AddAutoDockHost, removePCs,\
entropiaPresent, menuText



if entropiaPresent:
    from Entropia.EntropiaDef import entropia_job_dir
    from Entropia.EntropiaUI import EntropiaUI
    from Entropia.EntropiaEx import EntropiaError
    import ftplib

ADProcessManagerGUI=CommandGUI()
ADProcessManagerGUI.addMenuCommand('AutoTools4Bar', menuText['StartMB'], menuText['processManagerMB'])

AutoGridStarterGUI=CommandGUI()
AutoGridStarterGUI.addMenuCommand('AutoTools4Bar', menuText['StartMB'], menuText['startGridMB'])


AutoDockStarterGUI=CommandGUI()
AutoDockStarterGUI.addMenuCommand('AutoTools4Bar', menuText['StartMB'], menuText['startDockMB'])

        
AddAutoDockHostGUI=CommandGUI()
AddAutoDockHostGUI.addMenuCommand('AutoTools4Bar', menuText['StartMB'], menuText['editHostsMB'])

commandList = [
    {'name':'AD4start_autogrid','cmd':AutoGridStarter(),'gui':AutoGridStarterGUI},
    {'name':'AD4start_autodock','cmd':AutoDockStarter(),'gui':AutoDockStarterGUI},
    {'name':'AD4start_editHostMacros','cmd':AddAutoDockHost(),'gui':AddAutoDockHostGUI},
    ]

import sys, os
if os.name != 'nt': #not sys.platform == 'win32':
    commandList.insert(2,
    {'name':'ADstart_manage','cmd':ADProcessManager(),'gui':ADProcessManagerGUI})
else:
    import binaries
    os.environ['PATH'] = binaries.__path__[0]+";"+os.environ['PATH']
        
def initModule(vf):

    for dict in commandList:
        vf.addCommand(dict['cmd'],dict['name'],dict['gui'])

    if hasattr(vf,'GUI'):
        for item in vf.GUI.menuBars['AutoTools4Bar'].menubuttons.values():
            item.configure(background = 'tan')
            item.configure(underline = '-1')
    else:
        vf.addCommand(ADProcessManager(),'AD4start_manage')
        vf.addCommand(AutoGridStarter(), 'AD4start_autogrid')
        vf.addCommand(AutoDockStarter(), 'AD4start_autodock')
