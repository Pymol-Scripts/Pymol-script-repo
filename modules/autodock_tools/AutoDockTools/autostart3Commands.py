#############################################################################
#
# Author: Ruth HUEY, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################


# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/autostart3Commands.py,v 1.3 2009/06/17 00:03:24 vareille Exp $
#
# $Id: autostart3Commands.py,v 1.3 2009/06/17 00:03:24 vareille Exp $
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
ADProcessManagerGUI.addMenuCommand('AutoTools3Bar', menuText['StartMB'], menuText['processManagerMB'])


AutoGridStarterGUI=CommandGUI()
AutoGridStarterGUI.addMenuCommand('AutoTools3Bar', menuText['StartMB'], menuText['startGridMB'])


AutoDockStarterGUI=CommandGUI()
AutoDockStarterGUI.addMenuCommand('AutoTools3Bar', menuText['StartMB'], menuText['startDockMB'])

        
AddAutoDockHostGUI=CommandGUI()
AddAutoDockHostGUI.addMenuCommand('AutoTools3Bar', menuText['StartMB'], menuText['editHostsMB'])

commandList = [
    {'name':'AD3start_autogrid','cmd':AutoGridStarter(),'gui':AutoGridStarterGUI},
    {'name':'AD3start_autodock','cmd':AutoDockStarter(),'gui':AutoDockStarterGUI},
    {'name':'AD3start_editHostMacros','cmd':AddAutoDockHost(),'gui':AddAutoDockHostGUI},
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
        for item in vf.GUI.menuBars['AutoTools3Bar'].menubuttons.values():
            item.configure(background = 'tan')
            item.configure(underline = '-1')
    else:
        vf.addCommand(ADProcessManager(),'AD3start_manage')
        vf.addCommand(AutoGridStarter(), 'AD3start_autogrid')
        vf.addCommand(AutoDockStarter(), 'AD3start_autodock')
