## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

########################################################################
#
# Date: Nov 2001 Author: Daniel Stoffler
#
#    stoffler@scripps.edu
#
# Copyright: Daniel Stoffler, TSRI
#
#########################################################################

import Tkinter, numpy.oldnumeric as Numeric, string, math

from mglutil.util.callback import CallbackManager
from thumbwheel import ThumbWheel


class xyzGUI(Tkinter.Frame):

    """ 
    """
    def __init__(self, master=None, name='XYZ',callback=None,
                 callbackX=None, widthX=100, heightX=26, wheelPadX=4,
                 labcfgX={'text':None}, widthY=100, heightY=26, wheelPadY=4,
                 labcfgY={'text':None}, callbackY=None, widthZ=100,
                 heightZ=26, wheelPadZ=4, callbackZ=None,
                 labcfgZ={'text':None}, **kw):


        self.callback = callback   # user specified callback
        self.name=name               # title inside canvas

        self.widthX=widthX
        self.heightX=heightX
        self.wheelPadX=wheelPadX
        self.labcfgX=labcfgX
        self.widthY=widthY
        self.heightY=heightY
        self.wheelPadY=wheelPadY
        self.labcfgY=labcfgY
        self.widthZ=widthZ
        self.heightZ=heightZ
        self.wheelPadZ=wheelPadZ
        self.labcfgZ=labcfgZ

	Tkinter.Frame.__init__(self, master)
        Tkinter.Pack.config(self)

        self.callbacks = CallbackManager() # object to manage callback
                                        # functions. They get called with the
                                        # current value as an argument

        self.frame = Tkinter.Frame(self, relief = 'sunken', borderwidth=5)
        self.frame.pack(expand=1, fill='x')
        self.createEntries(self.frame)

        if self.callback:
            self.callbacks.AddCallback(self.callback)

   
    def createEntries(self, master):
        self.f = Tkinter.Frame(master)
	self.f.grid(column=1, rowspan=3)

        self.thumbx = ThumbWheel(master=self.f, width=self.widthX,
                                 height=self.heightX, labcfg=self.labcfgX,
                                 wheelPad=self.wheelPadX)
        self.thumbx.callbacks.AddCallback(self.thumbx_cb)
        self.thumbx.grid(row=0, column=0)

        self.thumby = ThumbWheel(master=self.f, width=self.widthY,
                                 height=self.heightY, labcfg=self.labcfgY,
                                 wheelPad=self.wheelPadY)
        self.thumby.callbacks.AddCallback(self.thumby_cb)
        self.thumby.grid(row=1, column=0)

        self.thumbz = ThumbWheel(master=self.f, width=self.widthZ,
                                 height=self.heightZ, labcfg=self.labcfgZ,
                                 wheelPad=self.wheelPadZ)
        self.thumbz.callbacks.AddCallback(self.thumbz_cb)
        self.thumbz.grid(row=2, column=0)

        self.f.pack(side='top', expand=1)


    def thumbx_cb(self, events=None):
        self.callbacks.CallCallbacks(self.thumbx.value)


    def thumby_cb(self, events=None):
        self.callbacks.CallCallbacks(self.thumby.value)


    def thumbz_cb(self, events=None):
        self.callbacks.CallCallbacks(self.thumbz.value)


    def set(self, x, y, z):
        # called from outside
        self.thumbx.setValue(x)
        self.thumby.setValue(y)
        self.thumbz.setValue(z)


 #####################################################################
 # the 'configure' methods:
 #####################################################################

    def configure(self, **kw):
        for key,value in kw.items():
            # the 'set parameter' callbacks
            if key=='labcfgX': self.setLabel(value,'x')
            elif key=='labcfgY': self.setLabel(value,'y')
            elif key=='labcfgZ': self.setLabel(value,'z')

            elif key=='continuousX': self.setContinuous(value,'x')
            elif key=='continuousY': self.setContinuous(value,'y')
            elif key=='continuousZ': self.setContinuous(value,'z')

            elif key=='precisionX': self.setPrecision(value,'x')
            elif key=='precisionY': self.setPrecision(value,'y')
            elif key=='precisionZ': self.setPrecision(value,'z')

            elif key=='typeX': self.setType(value,'x')
            elif key=='typeY': self.setType(value,'y')
            elif key=='typeZ': self.setType(value,'z')

            elif key=='minX': self.setMin(value,'x')
            elif key=='minY': self.setMin(value,'y')
            elif key=='minZ': self.setMin(value,'z')

            elif key=='maxX': self.setMax(value,'x')
            elif key=='maxY': self.setMax(value,'y')
            elif key=='maxZ': self.setMax(value,'z')

            elif key=='oneTurnX': self.setOneTurn(value,'x')
            elif key=='oneTurnY': self.setOneTurn(value,'y')
            elif key=='oneTurnZ': self.setOneTurn(value,'z')

            elif key=='showLabelX': self.setShowLabel(value,'x')
            elif key=='showLabelY': self.setShowLabel(value,'y')
            elif key=='showLabelZ': self.setShowLabel(value,'z')

            elif key=='incrementX': self.setIncrement(value,'x')
            elif key=='incrementY': self.setIncrement(value,'y')
            elif key=='incrementZ': self.setIncrement(value,'z')

            ####################################################

            elif key=='lockTypeX': self.lockType(value,'x')
            elif key=='lockTypeY': self.lockType(value,'y')
            elif key=='lockTypeZ': self.lockType(value,'z')

            elif key=='lockMinX': self.lockMin(value,'x')
            elif key=='lockMinY': self.lockMin(value,'y')
            elif key=='lockMinZ': self.lockMin(value,'z')

            elif key=='lockBMinX': self.lockBMin(value,'x')
            elif key=='lockBMinY': self.lockBMin(value,'y')
            elif key=='lockBMinZ': self.lockBMin(value,'z')

            elif key=='lockMaxX': self.lockMax(value,'x')
            elif key=='lockMaxY': self.lockMax(value,'y')
            elif key=='lockMaxZ': self.lockMax(value,'z')

            elif key=='lockBMaxX': self.lockBMax(value,'x')
            elif key=='lockBMaxY': self.lockBMax(value,'y')
            elif key=='lockBMaxZ': self.lockBMax(value,'z')

            elif key=='lockIncrementX': self.lockIncrement(value,'x')
            elif key=='lockIncrementY': self.lockIncrement(value,'y')
            elif key=='lockIncrementZ': self.lockIncrement(value,'z')

            elif key=='lockBIncrementX': self.lockBIncrement(value,'x')
            elif key=='lockBIncrementY': self.lockBIncrement(value,'y')
            elif key=='lockBIncrementZ': self.lockBIncrement(value,'z')

            elif key=='lockPrecisionX': self.lockPrecision(value,'x')
            elif key=='lockPrecisionY': self.lockPrecision(value,'y')
            elif key=='lockPrecisionZ': self.lockPrecision(value,'z')

            elif key=='lockShowLabelX': self.lockShowLabel(value,'x')
            elif key=='lockShowLabelY': self.lockShowLabel(value,'y')
            elif key=='lockShowLabelZ': self.lockShowLabel(value,'z')

            elif key=='lockValueX': self.lockValue(value,'x')
            elif key=='lockValueY': self.lockValue(value,'y')
            elif key=='lockValueZ': self.lockValue(value,'z')

            elif key=='lockContinuousX': self.lockContinuous(value,'x')
            elif key=='lockContinuousY': self.lockContinuous(value,'y')
            elif key=='lockContinuousZ': self.lockContinuous(value,'z')

            elif key=='lockOneTurnX': self.lockOneTurn(value,'x')
            elif key=='lockOneTurnY': self.lockOneTurn(value,'y')
            elif key=='lockOneTurnZ': self.lockOneTurn(value,'z')


    def setLabel(self, label, mode):
        if mode == 'x':
            self.thumbx.setLabel(label)
        elif mode == 'y':
            self.thumby.setLabel(label)
        elif mode == 'z':
            self.thumbz.setLabel(label)


    def setContinuous(self, value, mode):
        if mode == 'x':
            self.thumbx.setContinuous(value)
        elif mode == 'y':
            self.thumby.setContinuous(value)
        elif mode == 'z':
            self.thumbz.setContinuous(value)


    def setPrecision(self, value, mode):
        if mode == 'x':
            self.thumbx.setPrecision(value)
        elif mode == 'y':
            self.thumby.setPrecision(value)
        elif mode == 'z':
            self.thumbz.setPrecision(value)


    def setType(self, type, mode):
        if type == 'int': type = int
        else: type = float

        if mode == 'x':
            self.thumbx.setType(type)
        elif mode == 'y':
            self.thumby.setType(type)
        elif mode == 'z':
            self.thumbz.setType(type)


    def setMin(self, value, mode):
        if mode == 'x':
            self.thumbx.setMin(value)
        elif mode == 'y':
            self.thumby.setMin(value)
        elif mode == 'z':
            self.thumbz.setMin(value)


    def setMax(self, value, mode):
        if mode == 'x':
            self.thumbx.setMax(value)
        elif mode == 'y':
            self.thumby.setMax(value)
        elif mode == 'z':
            self.thumbz.setMax(value)


    def setOneTurn(self, value, mode):
        if mode == 'x':
            self.thumbx.setOneTurn(value)
        elif mode == 'y':
            self.thumby.setOneTurn(value)
        elif mode == 'z':
            self.thumbz.setOneTurn(value)


    def setShowLabel(self, value, mode):
        if mode == 'x':
            self.thumbx.setShowLabel(value)
        if mode == 'y':
            self.thumby.setShowLabel(value)
        if mode == 'z':
            self.thumbz.setShowLabel(value)
            

    def setIncrement(self, value, mode):
        if mode == 'x':
            self.thumbx.setIncrement(value)
        if mode == 'y':
            self.thumby.setIncrement(value)
        if mode == 'z':
            self.thumbz.setIncrement(value)



 #####################################################################
 # the 'lock' methods:
 #####################################################################


    def lockType(self, value, mode):
        if mode == 'x':
            self.thumbx.lockTypeCB(value)
        if mode == 'y':
            self.thumby.lockTypeCB(value)
        if mode == 'z':
            self.thumbz.lockTypeCB(value)


    def lockMin(self, value, mode):
        if mode == 'x':
            self.thumbx.lockMinCB(value)
        if mode == 'y':
            self.thumby.lockMinCB(value)
        if mode == 'z':
            self.thumbz.lockMinCB(value)


    def lockBMin(self, value, mode):
        if mode == 'x':
            self.thumbx.lockBMinCB(value)
        if mode == 'y':
            self.thumby.lockBMinCB(value)
        if mode == 'z':
            self.thumbz.lockBMinCB(value)


    def lockMax(self, value, mode):
        if mode == 'x':
            self.thumbx.lockMaxCB(value)
        if mode == 'y':
            self.thumby.lockMaxCB(value)
        if mode == 'z':
            self.thumbz.lockMaxCB(value)


    def lockBMax(self, value, mode):
        if mode == 'x':
            self.thumbx.lockBMaxCB(value)
        if mode == 'y':
            self.thumby.lockBMaxCB(value)
        if mode == 'z':
            self.thumbz.lockBMaxCB(value)


    def lockIncrement(self, value, mode):
        if mode == 'x':
            self.thumbx.lockIncrementCB(value)
        if mode == 'y':
            self.thumby.lockIncrementCB(value)
        if mode == 'z':
            self.thumbz.lockIncrementCB(value)


    def lockBIncrement(self, value, mode):
        if mode == 'x':
            self.thumbx.lockBIncrementCB(value)
        if mode == 'y':
            self.thumby.lockBIncrementCB(value)
        if mode == 'z':
            self.thumbz.lockBIncrementCB(value)


    def lockPrecision(self, value, mode):
        if mode == 'x':
            self.thumbx.lockPrecisionCB(value)
        if mode == 'y':
            self.thumby.lockPrecisionCB(value)
        if mode == 'z':
            self.thumbz.lockPrecisionCB(value)


    def lockShowLabel(self, value, mode):
        if mode == 'x':
            self.thumbx.lockShowLabelCB(value)
        if mode == 'y':
            self.thumby.lockShowLabelCB(value)
        if mode == 'z':
            self.thumbz.lockShowLabelCB(value)


    def lockValue(self, value, mode):
        if mode == 'x':
            self.thumbx.lockValueCB(value)
        if mode == 'y':
            self.thumby.lockValueCB(value)
        if mode == 'z':
            self.thumbz.lockValueCB(value)


    def lockContinuous(self, value, mode):
        if mode == 'x':
            self.thumbx.lockContinuousCB(value)
        if mode == 'y':
            self.thumby.lockContinuousCB(value)
        if mode == 'z':
            self.thumbz.lockContinuousCB(value)


    def lockOneTurn(self, value, mode):
        if mode == 'x':
            self.thumbx.lockOneTurnCB(value)
        if mode == 'y':
            self.thumby.lockOneTurnCB(value)
        if mode == 'z':
            self.thumbz.lockOneTurnCB(value)


if __name__ == '__main__':
    test = xyzGUI()
    def foo(val):
        print val
    test.callbacks.AddCallback(foo)
    test.configure(lockOneTurnX=1)
