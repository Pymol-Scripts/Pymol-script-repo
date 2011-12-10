#########################################################################
#
# Date: Jul 2002 Author: Daniel Stoffler
#
#    stoffler@scripps.edu
#
# Copyright:  Daniel Stoffler, and TSRI
#
#########################################################################

import sys,unittest
from time import sleep
from mglutil.gui.BasicWidgets.Tk.toolbarbutton import ToolBarButton
import Tkinter,Pmw


buttonFuncs = {}
widget = None
pushedButton = 0

def pause(sleepTime=0.2):
    widget.master.update()
    sleep(sleepTime)


def selectFunc(name):
    global buttonFuncs
    curFunc = buttonFuncs[name]
    if curFunc:
        curFunc()


def pushButton():
    global pushedButton
    pushedButton = 1

class ToolBarButtonBaseTest(unittest.TestCase):
    def test_01_constructor(self):
        # test if we can display a button
        global widget
        widget = Tkinter.Frame()
        widget.balloons = Pmw.Balloon(widget)
    
        ToolBarButton(balloonmaster=widget, master=widget, name='smiley',
                  icon1='smiley.gif') 

        widget.pack()
        pause(0.4)
        widget.master.destroy()


    def test_02_multipleButtons(self):
        # test if we can have multiple buttons
        global widget
        global buttonFuncs
        widget = Tkinter.Frame()
        widget.balloons = Pmw.Balloon(widget)
        # add separator, 5 smilies, separator
        for key, func, balloon in [
        ('sep1', None, None),
        ('smiley1', None, 'This is icon1'),
        ('smiley2', None, 'This is icon2'),
        ('smiley3', None, 'This is icon3'),
        ('smiley4', None, 'This is icon4'),
        ('smiley5', None, 'This is icon5'),
        ('sep2', None, None)
        ]:
    
            ToolBarButton(widget, widget, name=key, icon1='%s.gif' % key[:-1],
                      command=selectFunc, balloonhelp=balloon)
            buttonFuncs[key] = func
        widget.pack()
        pause(0.6)
        widget.master.destroy()


    def test_03_buttonCommand(self):
    # test if we can activate button function
        global widget
        global buttonFuncs
        global pushedButton

        widget = Tkinter.Frame()
        widget.balloons = Pmw.Balloon(widget)


        tb=ToolBarButton(balloonmaster=widget, master=widget, name='smiley',
                  icon1='smiley.gif', command=selectFunc, state='normal')
        buttonFuncs['smiley'] = pushButton

        widget.pack()

        # call the button command which will call the method pushButton
        # which will set the value pushedButton from 0 to 1
        tb.buttonFocus = 1 # simulate mouse pointer entering button, else the
                       # callback is not called
        tb.buttonUp(None)
        self.assertEqual(pushedButton == 1,True)
        pause(0.4)
        widget.master.destroy()
    
    
if __name__ == '__main__':
    unittest.main()
