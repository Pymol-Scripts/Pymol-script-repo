## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#########################################################################
#
# Date: Mai 2001 Authors: Michel Sanner, Daniel Stoffler
#
#    sanner@scripps.edu
#    stoffler@scripps.edu
#
# Copyright: Michel Sanner, Daniel Stoffler and TSRI
#
#########################################################################

import Tkinter
import math
import types
import sys
import os

from mglutil.util.callback import CallbackManager
from mglutil.util.misc import ensureFontCase
from optionsPanel import OptionsPanel

from KeyboardEntry import KeyboardEntry

class Dial(Tkinter.Frame, KeyboardEntry):
    """This class implements a Dial widget.
The widget has a pointer that can be moved around a circle.
The range corresponding to one full turn can be specified as well as the min
and max values that are allowed. By defaults these are set to None meaning that
there is no min and no max. One turn corresponds to 360 units by default.
A dial can also operate in discrete mode (if self.increment is set to x). In
this mode the values will be restrained to be multiples of self.increment.

The Widget has a Callback manager. Callback functions get called at every value
change if self.contiguous is set to 1, else they get called when the mouse
button is released. They always get called with the current value as an
argument.

An optional label can be displayed at the center of the Dial widget.
The size of the dial has to be specified at instanciation. Other parameters
can be set after the widget has been created.
The widget tried to adjust automatically the size of the arrow according to
the size of the dial.

The widget has a configure() method: type, min, max, increment, precision,
showLabel, value, continuous, oneTurn can be set this way.

master, labCfg and size can be passed only to the constructor.

a lock() method is used to disable the various gui components of the
options panel. Usage: <instance>.lock(<component>=<value>)
components see configure(). value is 0 or 1. 1 disables,
0 enables.

Setting values with increment enabled:
if using the method set(), the actual value will 'snap' to the next increment.
i.e., if the value is set to 3, and the increment is set to 2, setting the
value to 6 will actually result in 7  (3,5,7,9,.....)
To still be able to set the value, disregarding the current active increment,
the set method understands the optional keyword force=True, i.e.
dial.set(<value>, force=True)), which will set the value to <value>. The
increment will now be added to this new <value>
 
"""
    def __init__(self, master=None, type='float',
                 labCfg={'fg':'black','side':'left', 'text':None},
                 min=None, max=None, increment=.0, precision=2,
                 showLabel=1, value=0.0, continuous=1, oneTurn=360.,
                 size=50, callback=None,
                 lockMin=0, lockBMin=0, lockMax=0, lockBMax=0,
                 lockIncrement=0, lockBIncrement=0,
                 lockPrecision=0, lockShowLabel=0, lockValue=0,
                 lockType=0, lockContinuous=0, lockOneTurn=0, **kw):

	Tkinter.Frame.__init__(self, master)
        Tkinter.Pack.config(self)

        self.callbacks = CallbackManager() # object to manage callback
                                        # functions. They get called with the
                                        # current value as an argument

        # initialize various attributes with default values
        self.precision = 2              # decimal places
        self.min = None                 # minimum value
        self.max = None                 # maximum value
        self.increment = increment            # value increment
        self.minOld = 0.                # used to store old values 
        self.maxOld = 0.
        self.incrementOld = increment
        self.size = 50                  # defines widget size
        self.offsetValue = 0.           # used to set increment correctly
        self.lab = None                 # label
        self.callback = None            # user specified callback
        self.opPanel = None             # option panel widget
        self.oneTurn = 360.             # value increment for 1 full turn
        self.value = 0.0                # current value of widget
        self.oldValue = 0.0             # old value of widget
        self.showLabel = 1              # turn on to display label on
        self.continuous = 1             # set to 1 to call callbacks at 
                                        # each value change, else gets called 
                                        # on button release event
        self.angle = 0.                 # angle corresponding to value

        self.labCfg = labCfg            # Tkinter Label options
        self.labelFont = (
            ensureFontCase('helvetica'), 14, 'bold')    # label font
        self.labelColor = 'yellow'      # label color
        self.canvas = None              # the canvas to create the widget in
        self.usedArcColor = '#aaaaaa'   # filled arc color of used portion
        self.unusedArcColor = '#cccccc' # filled arc color of unused portion
        self.pyOver180 = math.pi/180.0  # constants used in various places
        self.threeSixtyOver1turn = 1
        self.piOver1turn = math.pi/360.

        self.lockMin = lockMin          # lock<X> vars are used in self.lock()
        self.lockMax = lockMax          # to lock/unlock entries in optionpanel
        self.lockIncrement = lockIncrement
        self.lockBMin = lockBMin
        self.lockBMax = lockBMax
        self.lockBIncrement = lockBIncrement
        self.lockPrecision = lockPrecision
        self.lockShowLabel = lockShowLabel
        self.lockValue = lockValue
        self.lockType = lockType
        self.lockContinuous = lockContinuous
        self.lockOneTurn = lockOneTurn

        self.setArrow()

        # configure with user-defined values
        self.setSize(size)
        self.setCallback(callback)     
        self.setContinuous(continuous)

        self.setType(type)              
        self.setPrecision(precision)
        self.setOneTurn(oneTurn)        
        self.setMin(min)
        self.setMax(max)
        self.setIncrement(increment)
        
        self.setShowLabel(showLabel)
        
        self.setValue(value)
        
        self.setLabel(self.labCfg)

        self.createCanvas(master)

        canvas = self.canvas
	canvas.bind("<ButtonPress-1>", self.mouseDown)
	canvas.bind("<ButtonRelease-1>", self.mouseUp)
	canvas.bind("<B1-Motion>", self.mouseMove)
        canvas.bind("<Button-3>", self.toggleOptPanel)

        if os.name == 'nt': #sys.platform == 'win32':
            canvas.bind("<MouseWheel>", self.mouseWheel)
        else:
            canvas.bind("<Button-4>", self.mouseWheel)
            canvas.bind("<Button-5>", self.mouseWheel)

        KeyboardEntry.__init__(self, (canvas,), self.setFromEntry)

        self.opPanel = OptionsPanel(master = self, title="Dial Options")
        
##         if self.callback:
##             self.callbacks.AddCallback(self.callback)


    def setFromEntry(self, valueString):
        try:
            self.set(self.type(valueString))
        except ValueError:
            # fixme we would like to pop this up in a window maybe
            import traceback
            traceback.print_stack()
            traceback.print_exc()


    def handleKeyStroke(self, event):
        # handle key strokes for numbers only in widget keyboard entry label
        key = event.keysym

        if key.isdigit() or key=='period' or key=='minus' or key=='plus':
            if key == 'period':
                key = '.'
            elif key == 'minus':
                key = '-'
            elif key == 'plus':
                key = '+'
            self.typedValue += key
            self.typedValueTK.configure(text=self.typedValue)
            
        else:
            KeyboardEntry.handleKeyStroke(self, event)


    def setSize(self, size):
        """Set widget size. Size must be of type int and greater than 0"""

        assert isinstance(size, types.IntType),\
               "Illegal size: expected type %s, got %s"%(type(1), type(size) )
        assert size > 0, "Illegal size: must be > 0, got %s"%size
        self.size = size


    def setCallback(self, cb):
        """Set widget callback. Must be callable function. Callback is called
every time the widget value is set/modified"""

        assert cb is None or callable(cb) or type(cb) is types.ListType,\
               "Illegal callback: must be either None or callable, or list. Got %s"%cb
        if cb is None: return
        elif type(cb) is types.ListType:
            for func in cb:
                assert callable(func), "Illegal callback must be callable. Got %s"%func
                self.callbacks.AddCallback(func)
        else:
            self.callbacks.AddCallback(cb)
        self.callback = cb
        

    def toggleOptPanel(self, event=None):
        if self.opPanel.flag:
           self.opPanel.Dismiss_cb()
        else:
            if not hasattr(self.opPanel, 'optionsForm'):
                self.opPanel.displayPanel(create=1)
            else:
                self.opPanel.displayPanel(create=0)


    def setArrow(self, size=None):
        if size is not None:
            self.setSize(size)
        aS = self.size/40
        self.arrowLength = max(3, 3*aS) # arrow head length
        self.arrowWidth = max(2, aS)    # half the arrow body width
        self.arrowBorderwidth = max(1, self.arrowWidth/2)  # width of arrow
                                                           # shadow lines
        self.arrowHeadWidth = 2*self.arrowWidth # width of arrow head base


    def mouseDown(self, event):
	# remember where the mouse went down
	self.lastx = event.x
	self.lasty = event.y


    def mouseUp(self, event):
	# call callbacks if not in continuous mode
        if not self.continuous:
            self.callbacks.CallCallbacks(self.opPanel.valInput.get())
        if self.showLabel == 2:
            # no widget labels on mouse release
            self.canvas.itemconfigure(self.labelId2, text='')
            self.canvas.itemconfigure(self.labelId, text='')


    def mouseMove(self, event):
        dx = event.x-self.xm
        dy = self.ym-event.y
    
        n = math.sqrt(dx*dx+dy*dy)
        if n == 0.0: v = [0.0, 0.0]
        else: v = [dx/n, dy/n]

        # find the cosine of the angle between new hand position and previous
        # hand position
        ma = v[0]*self.vector[0] + v[1]*self.vector[1]
        
        # assure no rounding errors
        if ma > 1.0: ma = 1.0
        elif ma < -1.0: ma = -1.0

        # compute angle increment compared to current vector
        ang = math.acos(ma)

        # find the sign of the rotation, sign of z component of vector prod.
        oldv = self.vector
        normz = oldv[0]*v[1] - oldv[1]*v[0]
        if normz>0: ang = -1. * ang
        
        # compute the new value
        val = self.value + ang*self.oneTurnOver2pi

        self.set(val)

        self.lastx = event.x
        self.lasty = event.y


    def mouseWheel(self, event):
        #print "mouseWheel", event, event.num
        if os.name == 'nt': #sys.platform == 'win32':
            if event.delta > 0:
                lEventNum = 4
            else:
                lEventNum = 5
        else:
            lEventNum = event.num
        if lEventNum == 4:
            self.set(self.value+self.oneTurn)
        else:
            self.set(self.value-self.oneTurn)


    def get(self):
        return self.type(self.value)


    def printLabel(self):
        if self.canvas is None:
            return
        self.canvas.itemconfigure(self.labelId2,
                                  text=self.labelFormat%self.value)#newVal)
        self.canvas.itemconfigure(self.labelId,
                                  text=self.labelFormat%self.value)#newVal)
        
    
    def set(self, val, update=1, force=0):
        # if force is set to 1, we call this method regardless of the
        # widget configuration. This is for example the case if the dial
        # is set to continuous=0, but the value is set in the options panel

        # snap to closest increment
        if self.increment is not None and self.increment != 0. and not force:
            offset = self.offsetValue%self.increment
            dval = round(val/self.increment) * self.increment

            if val < dval:
                dval = dval + offset - self.increment
            else:
                dval = dval + offset
        
            if self.min is not None and dval < self.min:
                dval = self.min
            elif self.max is not None and dval > self.max:
                dval = self.max

            # recompute vector and angle corresponding to val
            self.angle = (dval%self.oneTurn)*self.threeSixtyOver1turn
            if dval <0.0:
                self.angle = self.angle - 360.0
            a = self.angle*self.pyOver180
            self.vector = [math.sin(a), math.cos(a)]
            self.value = dval
            self.offsetValue = dval

        else:
            # 'regular' mode, i.e. no step-wise increment
            if self.min is not None and val < self.min: val = self.min
            elif self.max is not None and val > self.max: val = self.max
            
            # recompute vector and angle corresponding to val
            self.angle = (val%self.oneTurn)*self.threeSixtyOver1turn
            if val <0.0: self.angle = self.angle - 360.0
            a = self.angle*self.pyOver180
            self.vector = [math.sin(a), math.cos(a)]
            self.value = val
            self.offsetValue = val

        #update arrow in display
        self.drawArrow()
        newVal = self.get()
        
        if self.continuous or force:
            if update and self.oldValue != newVal or force:
                self.oldValue = newVal
                self.callbacks.CallCallbacks(newVal)
            if self.showLabel==2:
                self.printLabel()
        else:
            if self.showLabel==2:
                self.printLabel()
                
        if self.showLabel==1:
            self.printLabel()
        if self.opPanel:
            self.opPanel.valInput.set(self.labelFormat%newVal)


    def drawArrow(self):
        if self.canvas is None:
            return
        # end point
        x1 = self.xm + self.vector[0]*self.rad
        y1 = self.ym - self.vector[1]*self.rad
        # point at arrow head base
        xb = self.xm + self.vector[0]*self.radNoArrow
        yb = self.xm - self.vector[1]*self.radNoArrow
        # vector orthogonal to arrow
        n = [-self.vector[1], -self.vector[0]]
        pts1 = [ self.xm+n[0]*self.arrowWidth, self.ym+n[1]*self.arrowWidth,
                 xb+n[0]*self.arrowWidth, yb+n[1]*self.arrowWidth,
                 xb+n[0]*self.arrowHeadWidth, yb+n[1]*self.arrowHeadWidth,
                 x1, y1 ]
        pts2 = [ x1, y1,
                 xb-n[0]*self.arrowHeadWidth, yb-n[1]*self.arrowHeadWidth,
                 xb-n[0]*self.arrowWidth, yb-n[1]*self.arrowWidth,
                 self.xm-n[0]*self.arrowWidth, self.ym-n[1]*self.arrowWidth ]
        canvas = self.canvas
        if self.vector[0] > 0.0:
            col1 = '#DDDDDD'
            col2 = 'black'
        else:
            col1 = 'black'
            col2 = '#DDDDDD'
        apply( canvas.coords, (self.arrowPolId,) + tuple(pts1+pts2) )
        apply( canvas.coords, (self.arrowPolborder1,) + tuple(pts1) )
        canvas.itemconfigure( self.arrowPolborder1, fill=col1 )
        apply( canvas.coords, (self.arrowPolborder2,) + tuple(pts2) )
        canvas.itemconfigure( self.arrowPolborder2, fill=col2 )
        canvas.itemconfigure(self.arcId, extent = 0.0-self.angle)

       
    def createCanvas(self, master):
        size = self.size
        self.frame = Tkinter.Frame(self, borderwidth=3, relief='sunken')

	self.canvas = Tkinter.Canvas(self.frame, width=size+2, height=size+2)

        self.xm = self.ym = size/2+2
        self.rad = size/2
        self.radNoArrow = self.rad-self.arrowLength
        self.vector = [0, 1]
        x1 = self.xm + self.vector[0]*self.rad
        y1 = self.ym + self.vector[1]*self.rad
        canvas = self.canvas
        self.circleId = canvas.create_oval(2,2,size,size, width=1,
                                           fill=self.unusedArcColor)
        self.arcId = canvas.create_arc(2,2,size,size, start=90.,
                                            extent=0, fill=self.usedArcColor)
        canvas.create_line(2, self.ym, size+2, self.ym)
        canvas.create_line(self.xm, 2, self.ym, size+2)

        self.arrowPolId = canvas.create_polygon( 0,0,0,0,0,0,0,0,
                                                 0,0,0,0,0,0,0,0,
                                                 fill='gray75' )
        self.arrowPolborder1 = canvas.create_line( 0,0,0,0,0,0,0,0,
                                                   fill='black',
                                             width = self.arrowBorderwidth)
        self.arrowPolborder2 = canvas.create_line( 0,0,0,0,0,0,0,0,
                                                   fill='white',
                                             width = self.arrowBorderwidth )

        r = size/20
        off = self.arrowBorderwidth
        canvas.create_oval(self.xm-r,self.ym-r-off/2,self.xm+r,self.ym+r-off/2,
                           fill='#DDDDDD', outline='white')

        canvas.create_oval(self.xm-r,self.ym-r+off,self.xm+r,self.ym+r+off,
                           fill='black', outline='black')
        
        canvas.create_oval(self.xm-r,self.ym-r,self.xm+r,self.ym+r,
                           fill='gray70', outline='#DDDDDD')
       
        
        self.labelId2 = canvas.create_text(self.xm+2, self.ym+2,
                                           fill='black',
                                           justify='center', text='',
                                           font = self.labelFont)
        self.labelId = canvas.create_text(self.xm, self.ym,
                                          fill=self.labelColor,
                                          justify='center', text='',
                                          font = self.labelFont)

        self.drawArrow()

        self.opPanel = OptionsPanel(master = self, title="Dial Options")

        # pack em up
        self.canvas.pack(side=Tkinter.TOP)
        self.frame.pack(expand=1, fill='x')
        self.toggleWidgetLabel(self.showLabel)


    def toggleWidgetLabel(self, val):
        if val == 0:
            # no widget labels
            self.showLabel=0
            self.canvas.itemconfigure(self.labelId2,
                                      text='')
            self.canvas.itemconfigure(self.labelId,
                                      text='')

        if val == 1:
            # show always widget labels
            self.showLabel=1
            self.printLabel()

        if val == 2:
            # show widget labels only when mouse moves
            self.showLabel=2
            self.canvas.itemconfigure(self.labelId2,
                                      text='')
            self.canvas.itemconfigure(self.labelId,
                                      text='')


    def setValue(self, val):

        if type(val) == types.StringType:
            val = float(val)

        assert type(val) in [types.IntType, types.FloatType],\
               "Illegal type for value: expected %s or %s, got %s"%(
                   type(1), type(1.0), type(val) )

        # setValue does NOT call a callback!
        if self.min is not None and val < self.min: val = self.min
        if self.max is not None and val > self.max: val = self.max
        self.value = self.type(val)
        self.offsetValue=self.value
        self.oldValue = self.value

        #update arrow in display
        self.angle = (self.value%self.oneTurn)*self.threeSixtyOver1turn
        if self.value <0.0: self.angle = self.angle - 360.0
        a = self.angle*self.pyOver180
        self.vector = [math.sin(a), math.cos(a)]
        self.drawArrow()

        if self.showLabel == 1:
            self.printLabel()
        if self.opPanel:
            self.opPanel.valInput.set(self.labelFormat%self.value)


    def setLabel(self, labCfg):
        self.labCfg = labCfg

        text = labCfg.get('text', None)
        if text is None or text=='':
            return

        d={}
        for k, w in self.labCfg.items():
            if k == 'side': continue
            else: d[k] = w
        if not 'side' in self.labCfg.keys():
            self.labCfg['side'] = 'left'

        if not self.lab:
            self.lab = Tkinter.Label(self, d)
            self.lab.pack(side=self.labCfg['side'])
            self.lab.bind("<Button-3>", self.toggleOptPanel)
        else:
            self.lab.configure(text)

            
 #####################################################################
 # the 'configure' methods:
 #####################################################################

    def configure(self, **kw):
        for key,value in kw.items():
            # the 'set' parameter callbacks
            if key=='labCfg': self.setLabel(value)
            elif key=='type': self.setType(value)
            elif key=='min': self.setMin(value)
            elif key=='max': self.setMax(value)
            elif key=='increment': self.setIncrement(value)
            elif key=='precision': self.setPrecision(value)
            elif key=='showLabel': self.setShowLabel(value)
            elif key=='continuous': self.setContinuous(value)
            elif key=='oneTurn': self.setOneTurn(value)

            # the 'lock' entries callbacks
            elif key=='lockType': self.lockTypeCB(value)
            elif key=='lockMin': self.lockMinCB(value)
            elif key=='lockBMin': self.lockBMinCB(value)
            elif key=='lockMax': self.lockMaxCB(value)
            elif key=='lockBMax': self.lockBMaxCB(value)
            elif key=='lockIncrement': self.lockIncrementCB(value)
            elif key=='lockBIncrement': self.lockBIncrementCB(value)
            elif key=='lockPrecision': self.lockPrecisionCB(value)
            elif key=='lockShowLabel': self.lockShowLabelCB(value)
            elif key=='lockValue': self.lockValueCB(value)
            elif key=='lockContinuous': self.lockContinuousCB(value)
            elif key=='lockOneTurn': self.lockOneTurnCB(value)

            
    def setType(self, Type):

        assert type(Type) in [types.StringType, types.TypeType],\
               "Illegal type for datatype. Expected %s or %s, got %s"%(
                   type('a'), type(type), type(Type) )
        
        if type(Type) == type(""): # type str
            assert Type in ('int','float'),\
            "Illegal type descriptor. Expected 'int' or 'float', got '%s'"%Type
            self.type = eval(Type)
        else:
            self.type = Type
        
        if self.type == int:
            self.labelFormat = "%d"
            self.int_value = self.value
        else:
            self.labelFormat = "%."+str(self.precision)+"f"

        if hasattr(self.opPanel, 'optionsForm'):
            w = self.opPanel.idf.entryByName['togIntFloat']['widget']
            if self.type == int:
                w.setvalue('int')
            elif self.type == 'float':
                w.setvalue('float')
            
        if self.opPanel:
            self.opPanel.updateDisplay()

        # and update the printed label
        if self.canvas and self.showLabel == 1:
            self.printLabel()


    def setMin(self, min):
        if min is not None:

            assert type(min) in [types.IntType, types.FloatType],\
                 "Illegal type for minimum. Expected type %s or %s, got %s"%(
                     type(0), type(0.0), type(min) )
            
            if self.max and min > self.max:
                min = self.max

            self.min = self.type(min)

            if self.showLabel == 1:
                self.printLabel()

            if self.value < self.min:
                self.set(self.min)

            if hasattr(self.opPanel, 'optionsForm'):
                self.opPanel.minInput.set(self.labelFormat%self.min) 

                self.opPanel.toggleMin.set(1)
                self.opPanel.min_entry.configure(state='normal', fg='gray0')
            self.minOld = self.min    
                
        else:
            self.min = None
            if hasattr(self.opPanel, 'optionsForm'):
                self.opPanel.toggleMin.set(0)
                self.opPanel.min_entry.configure(state='disabled',
                                                 fg='gray40')


    def setMax(self, max):
        if max is not None:

            assert type(max) in [types.IntType, types.FloatType],\
                 "Illegal type for maximum. Expected type %s or %s, got %s"%(
                     type(0), type(0.0), type(max) )
             
            if self.min and max < self.min:
                max = self.min

            self.max = self.type(max)

            if self.showLabel == 1:
                self.printLabel()

            if self.value > self.max:
                self.set(self.max)

            if hasattr(self.opPanel, 'optionsForm'):
                self.opPanel.maxInput.set(self.labelFormat%self.max) 

                self.opPanel.toggleMax.set(1)
                self.opPanel.max_entry.configure(state='normal', fg='gray0')
            self.maxOld = self.max
                
        else:
            self.max = None
            if hasattr(self.opPanel, 'optionsForm'):
                self.opPanel.toggleMax.set(0)
                self.opPanel.max_entry.configure(state='disabled', fg='gray40')


    def setIncrement(self, incr):
        if incr is not None:

            assert type(incr) in [types.IntType, types.FloatType],\
                "Illegal type for increment. Expected type %s or %s, got %s"%(
                    type(0), type(0.0), type(incr) )

            self.increment = self.type(incr)
            self.offsetValue = self.value
            self.incrementOld = self.increment
            if hasattr(self.opPanel, 'optionsForm'):
                self.opPanel.incrInput.set(self.labelFormat%self.increment) 
                self.opPanel.toggleIncr.set(1)
                self.opPanel.incr_entry.configure(state='normal', fg='gray0')
                
        else:
            self.increment = self.type(0)
            if hasattr(self.opPanel, 'optionsForm'):
                self.opPanel.toggleIncr.set(0)
                self.opPanel.incrInput.set(self.labelFormat%0)
                self.opPanel.incr_entry.configure(state='disabled',
                                                  fg='gray40')


    def setPrecision(self, val):
        assert type(val) in [types.IntType, types.FloatType],\
               "Illegal type for precision. Expected type %s or %s, got %s"%(
                     type(0), type(0.0), type(val) )
        val = int(val)
        
        if val > 10:
            val = 10
        if val < 1:
            val = 1
        self.precision = val

        if self.type == float:
            self.labelFormat = "%."+str(self.precision)+"f"
        else:
            self.labelFormat = "%d"

        if hasattr(self.opPanel, 'optionsForm'):
            w = self.opPanel.idf.entryByName['selPrec']['widget']
            w.setvalue(val)
            
        if self.opPanel:
            self.opPanel.updateDisplay()

        # and update the printed label
        if self.canvas and self.showLabel == 1:
            self.printLabel()
            
    def setContinuous(self, cont):
        """ cont can be None, 0 or 1 """

        assert cont in [None, 0, 1],\
             "Illegal value for continuous: expected None, 0 or 1, got %s"%cont

        if cont != 1:
            cont = None
        self.continuous = cont
        if hasattr(self.opPanel, 'optionsForm'):
            w = self.opPanel.idf.entryByName['togCont']['widget']
            if cont:
                w.setvalue('on')#i=1
            else:
                w.setvalue('off')#i=0

        if self.opPanel:
            self.opPanel.updateDisplay()


    def setShowLabel(self, val):
        """Show label can be 0, 1 or 2
0: no label
1: label is always shown
2: show label only when value changes"""
        
        assert val in [0,1,2],\
               "Illegal value for showLabel. Expected 0, 1 or 2, got %s"%val
        
        if val != 0 and val != 1 and val != 2:
            print "Illegal value. Must be 0, 1 or 2"
            return
        self.showLabel = val
        self.toggleWidgetLabel(val)

        if hasattr(self.opPanel, 'optionsForm'):
            w = self.opPanel.idf.entryByName['togLabel']['widget']
            if self.showLabel == 0:
                label = 'never'
            elif self.showLabel == 1:
                label = 'always'
            elif self.showLabel == 2:
                label = 'move'
            w.setvalue(label)

        if self.opPanel:
            self.opPanel.updateDisplay()


    def setOneTurn(self, oneTurn):
        assert type(oneTurn) in [types.IntType, types.FloatType],\
               "Illegal type for oneTurn. Expected %s or %s, got %s"%(
                   type(0), type(0.0), type(oneTurn) )
        
        self.oneTurn = oneTurn
        self.threeSixtyOver1turn = 360./oneTurn
        self.piOver1turn = math.pi/oneTurn
        self.oneTurnOver2pi = oneTurn / (2*math.pi)

        if self.opPanel:
            self.opPanel.updateDisplay()


 #####################################################################
 # the 'lock' methods:
 #####################################################################


    def lockTypeCB(self, mode):
        if mode != 0: mode = 1
        self.lockType = mode
        if hasattr(self.opPanel, 'optionsForm'):
            self.opPanel.lockUnlockDisplay()
        

    def lockMinCB(self, mode): #min entry field
        if mode != 0: mode = 1
        self.lockMin = mode
        if hasattr(self.opPanel, 'optionsForm'):
            self.opPanel.lockUnlockDisplay()
        

    def lockBMinCB(self, mode): # min checkbutton
        if mode != 0: mode = 1
        self.lockBMin = mode
        if hasattr(self.opPanel, 'optionsForm'):
            self.opPanel.lockUnlockDisplay()


    def lockMaxCB(self, mode): # max entry field
        if mode != 0: mode = 1
        self.lockMax = mode
        if hasattr(self.opPanel, 'optionsForm'):
            self.opPanel.lockUnlockDisplay()


    def lockBMaxCB(self, mode): # max checkbutton
        if mode != 0: mode = 1
        self.lockBMax = mode
        if hasattr(self.opPanel, 'optionsForm'):
            self.opPanel.lockUnlockDisplay()


    def lockIncrementCB(self, mode): # increment entry field
        if mode != 0: mode = 1
        self.lockIncrement = mode
        if hasattr(self.opPanel, 'optionsForm'):
            self.opPanel.lockUnlockDisplay()


    def lockBIncrementCB(self, mode): # increment checkbutton
        if mode != 0: mode = 1
        self.lockBIncrement = mode
        if hasattr(self.opPanel, 'optionsForm'):
            self.opPanel.lockUnlockDisplay()


    def lockPrecisionCB(self, mode):
        if mode != 0: mode = 1
        self.lockPrecision = mode
        if hasattr(self.opPanel, 'optionsForm'):
            self.opPanel.lockUnlockDisplay()
            

    def lockShowLabelCB(self, mode):
        if mode != 0: mode = 1
        self.lockShowLabel = mode
        if hasattr(self.opPanel, 'optionsForm'):
            self.opPanel.lockUnlockDisplay()


    def lockValueCB(self, mode):
        if mode != 0: mode = 1
        self.lockValue = mode
        if hasattr(self.opPanel, 'optionsForm'):
            self.opPanel.lockUnlockDisplay()


    def lockContinuousCB(self, mode):
        if mode != 0: mode = 1
        self.lockContinuous = mode
        if hasattr(self.opPanel, 'optionsForm'):
            self.opPanel.lockUnlockDisplay()


    def lockOneTurnCB(self, mode):
        if mode != 0: mode = 1
        self.lockOneTurn = mode
        if hasattr(self.opPanel, 'optionsForm'):
            self.opPanel.lockUnlockDisplay()


if __name__ == '__main__':
    def foo(val):
        print val
    d = Dial(size=50)
    d.configure(showLabel=1)
    d.callbacks.AddCallback(foo)




