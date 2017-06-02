'''
Described at: http://www.pymolwiki.org/index.php/resicolor_plugin
plugin contributed by Philippe Garteiser garteiserp@omrf.org
'''
from pymol import cmd

import sys
if sys.version_info[0] < 3:
    import tkSimpleDialog
else:
    from tkinter import simpledialog as tkSimpleDialog


def __init__(self):
    # Simply add the menu entry and callback
    self.menuBar.addmenuitem('Plugin', 'command',
                             'resicolor',
                             label='resicolor',
                             command=lambda s=self: getselection(s))


def resicolor(selection):
    if selection:   # None is returned for user cancel
        cmd.select('calcium', 'resn ca or resn cal')
        cmd.select('acid', 'resn asp or resn glu or resn cgu')
        cmd.select('basic', 'resn arg or resn lys or resn his')
        cmd.select('nonpolar', 'resn met or resn phe or resn pro or resn trp or resn val or resn leu or resn ile or resn ala')
        cmd.select('polar', 'resn ser or resn thr or resn asn or resn gln or resn tyr')
        cmd.select('cys', 'resn cys or resn cyx')
        cmd.select('backbone', 'name ca or name n or name c or name o')
        cmd.select('none')
        code = {'acid': 'red',
              'basic': 'blue',
              'nonpolar': 'orange',
              'polar': 'green',
              'cys': 'yellow'}
        cmd.select('none')
        for elem in code:
            line = 'color ' + code[elem] + ',' + elem + '&' + selection
            print(line)
            cmd.do(line)
        word = 'color white,backbone &' + selection
        print(word)
        cmd.do(word)


def getselection(app):
    selection = tkSimpleDialog.askstring('resicolor',
                                       'Please enter a selection',
                                       parent=app.root)
    resicolor(selection)
