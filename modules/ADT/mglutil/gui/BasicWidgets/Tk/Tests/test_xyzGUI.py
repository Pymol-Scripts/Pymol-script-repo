#########################################################################
#
# Date: Jun 2002 Authors: Daniel Stoffler, Michel Sanner
#
#    stoffler@scripps.edu
#    sanner@scripps.edu
#
# Copyright:  Daniel Stoffler, Michel Sanner, and TSRI
#
#########################################################################

import sys,unittest
from time import sleep
from mglutil.gui.BasicWidgets.Tk.xyzGUI import xyzGUI
widget = None
wasCalled = 0

def pause(sleepTime=0.2):
    widget.master.update()
    sleep(sleepTime)

class  XYZGUIBaseTest(unittest.TestCase):
    def test_constructor(self):
        # test if we can display the xyzGUI (consists of 3 thumbwheels)
        global widget
        widget = xyzGUI()
        pause(0.6)
        widget.master.destroy()

if __name__ == '__main__':
    unittest.main()
