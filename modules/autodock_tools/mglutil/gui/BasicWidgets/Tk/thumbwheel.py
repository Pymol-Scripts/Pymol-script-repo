## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#########################################################################
#
# Date: May 2001 Authors: Michel Sanner, Daniel Stoffler
#
#    sanner@scripps.edu
#    stoffler@scripps.edu
#
# Copyright: Michel Sanner, Daniel Stoffler and TSRI
#
#########################################################################

import Tkinter
import numpy.oldnumeric as Numeric
import math
import types
import sys
import os

from mglutil.util.callback import CallbackManager
from mglutil.util.misc import ensureFontCase
from optionsPanel import OptionsPanel

from KeyboardEntry import KeyboardEntry

class ThumbWheel(Tkinter.Frame, KeyboardEntry):
    """ This class builds a thumbwheel put onto a wheelpad.
constructor options:
- master: the master into the thumwheel can be packed. If one is specified,
  the widget gets packed in a toplevel() window
- height, width, orient specify the size and orientation of the widget.
- wheelpad specifies the pad onto which the thumbwheel gets packed
  optional is a canvascfg to configure the wheelpad canvas Tkinter options.
  for example if the wheelpad needs a blue background: canvascfg={'bg':'blue'}
- continuous: boolean. When set to True, continuous is 'on' and the callback functions
  will be called each time the value changes. Otherwise continuous will be 'off' and
  callback functions will be called on mouse release.
- callback: (None, callable function or list of callable functions).to specify function
  to be called when the value of the thumbwheel is modified.
- type: ('float', 'int' ...) string describing the type of the thumbwheel
- min, max, increment, precision, oneTurn specify the parameters of the thumbwheel.
- labCfg describes the label of the thumbwheel which will be packed to the left of the widget by default unless 'side' is specified. Possible keywords are text and side.

- wheelLabCfg describes the label located on the wheel where the value of the thumbwheel
  will be displayed
- showLabel Flag to specify whether to display the wheelLabel always 1, never 0,
  only on motion 2.
- canvasCfg describe the canvas containing the thumbwheel.
An option panel is available to the user to modify the thumbwheel settings by right clicking
on the widget
- lockMin, lockBMin, lockMax, lockBMax, lockIncrement, lockBIncrement,
  lockPrecision, lockShowLabel, lockValue,  lockType, lockContinuous, lockOneTurn:
  These flags specifies whether (when set to 1) or not (when set to 0) the user will be
  allowed to change the setting of the thumbwheel.
- reportDelta is flag to specify whether value differences should be reported
  rathe than absolute values
  
The widget has a configure() method: type, min, max, increment, precision,
showLabel, value, continuous, oneTurn, orient, reportDelta  can be set this
way.

a lock() method is used to disable the various gui components of the
options panel. Usage: <instance>.configure(<component>=<value>)
components see configure(). value is 0 or 1. 1 disables,0 enables.    
    """

#FIXME: this would be the mechanism to remove all arguments from the
# init (see DejaVu)
#    keywords = [
#        'master',
#        'oneTurn',
#        ]
    def __init__(self, master=None, labCfg={'fg':'black','side':'left',
                                            'text':None},
                 wheelLabCfg={}, canvasCfg={}, callback=None,
                 type='float', min=None, max=None, increment=0.0, precision=2,
                 showLabel=1, value=0.0, continuous=True, width=200,
                 height=40, oneTurn=10., wheelPad=6,
                 lockMin=0, lockBMin=0, lockMax=0,lockBMax=0,
                 lockIncrement=0, lockBIncrement=0,
                 lockPrecision=0, lockShowLabel=0, lockValue=0,
                 lockType=0, lockContinuous=0, lockOneTurn=0,
                 orient=None, reportDelta=0, **kw):

# See FIXME init
#        if __debug__:
#            checkkeywords(kw)
            
	Tkinter.Frame.__init__(self, master)
	Tkinter.Pack.config(self, side='left', anchor='w')

        #FIXME: nblines are not dynamically computed
        self.nblines = 30
        self.callback = None
        self.callbacks = CallbackManager()  # object to manage callback
                                        # functions. They get called with the
                                        # current value as an argument

        self.width = 200
        self.height = 40
        self.setWidth(width)
        self.setHeight(height)
        # set widget orientation: either horizontal or vertical
        self.setOrient(orient)


        # initialize various attributes with default values
        self.precision = 2              # decimal places
        self.min = None                 # minimum value
        self.max = None                 # maximum value
        self.increment = increment      # value increment
        self.minOld = 0.                # used to store old values 
        self.maxOld = 0.
        self.incrementOld = increment
        self.size = 50                  # defines widget size
        self.offsetValue = 0.           # used to set increment correctly
        self.lab = None                 # label
        self.opPanel = None             # option panel widget
        self.oneTurn = 360.             # value increment for 1 full turn
        self.value = 0.0                # current value of widget
        self.oldValue = 0.0             # old value of widget
        self.showLabel = 1              # turn on to display label on
        self.continuous = True          # set to 1 to call callbacks at 
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

        self.wheelLabCfg = wheelLabCfg  # Tkinter wheel label options
        self.canvasCfg = canvasCfg      # Tkinter Canvas options
        self.wheelPad = wheelPad        # pad between wheel and widget borders
        self.deltaVal = 0.              # delta with value at last callback

        self.valueLabel = None          # Tkinter Label
        self.setType(type)              # can be float or int

        self.discreteValue = 0.         # value in discrete space

        self.lockMin = lockMin          # lock<X> variables in configure()
        self.lockMax = lockMax          # to lock/unlock entries in options
                                        # panel
        

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
        self.reportDelta = reportDelta

        self.setLabel(self.labCfg)

        self.createCanvas(master, wheelPad=wheelPad)        
	self.canvas.bind("<ButtonPress-1>", self.mouseDown)
	self.canvas.bind("<ButtonRelease-1>", self.mouseUp)
	self.canvas.bind("<B1-Motion>", self.mouseMove)
        self.canvas.bind("<Button-3>", self.toggleOptPanel)

 	self.valueLabel.bind("<ButtonPress-1>", self.mouseDown)
	self.valueLabel.bind("<ButtonRelease-1>", self.mouseUp)
	self.valueLabel.bind("<B1-Motion>", self.mouseMove)
        self.valueLabel.bind("<Button-3>", self.toggleOptPanel)

        if os.name == 'nt': #sys.platform == 'win32':
            self.canvas.bind("<MouseWheel>", self.mouseWheel)
            self.valueLabel.bind("<MouseWheel>", self.mouseWheel)
        else:
            self.canvas.bind("<Button-4>", self.mouseWheel)
            self.valueLabel.bind("<Button-4>", self.mouseWheel)
            self.canvas.bind("<Button-5>", self.mouseWheel)
            self.valueLabel.bind("<Button-5>", self.mouseWheel)


        self.bind("<Button-3>", self.toggleOptPanel)

        KeyboardEntry.__init__(self, (self, self.canvas, self.valueLabel),
                               self.setFromEntry)
        
        self.opPanel = OptionsPanel(master = self, title="Thumbwheel Options")

        # now set the constructor options correctly using the configure method
        self.setValue(value)
        apply( self.configure, (),
               {'type':type, 'min':min, 'max':max, 'increment':increment,
                'precision':precision, 'showLabel':showLabel,
                'continuous':continuous, 'oneTurn':oneTurn,
                'lockType':lockType, 'lockMin':lockMin, 'lockBMin':lockBMin,
                'lockMax':lockMax, 'lockBMax':lockBMax,
                'lockIncrement':lockIncrement, 'lockBIncrement':lockBIncrement,
                'lockPrecision':lockPrecision, 'lockShowLabel':lockShowLabel,
                'lockValue':lockValue, 'lockContinuous':lockContinuous,
                'lockOneTurn':lockOneTurn, 'orient':orient,
                'reportDelta':reportDelta,
                'callback':callback
                })
        

    def setFromEntry(self, valueString):
        try:
            self.set(self.type(valueString), force=1)
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


    def setWidth(self, width):
        assert isinstance(width, types.IntType),\
               "Illegal width: expected %s, got %s"%(type(1),type(width))
        assert width > 0,"Illegal width: must be > 0, got %s"%width
        self.width = width


    def setHeight(self, height):
        assert isinstance(height, types.IntType),\
               "Illegal height: expected %s, got %s"%(type(1),type(height))
        assert height > 0,"Illegal height: must be > 0, got %s"%height
        self.height = height


    def setCallbacks(self, cb):
        """Set widget callback. Must be callable function. Callback is called
every time the widget value is set/modified"""

        assert cb is None or callable(cb) or type(cb) is types.ListType,\
               "Illegal callback: must be either None or callable. Got %s"%cb
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

 
    def mouseDown(self, event):
	# remember where the mouse went down
        if self.orient=='horizontal':
            self.lastx = event.x
        else:
            self.lastx = event.y


    def mouseUp(self, event):
	# call callbacks if not in continuous mode
        newVal = self.get()
        if self.oldValue != newVal:
            if not self.continuous:
                self.oldValue = newVal
                self.callbacks.CallCallbacks(newVal)
            if self.showLabel == 2:
                # no widget labels on mouse release
                self.valueLabel.place_forget()


    def mouseMove(self, event):
        if self.orient=='horizontal':
            dx = event.x - self.lastx
            self.lastx = event.x
        else:
            dx = event.y - self.lastx
            self.lasty = event.y
        dang = dx * math.pi / 180.
        val = dang*self.oneTurnOver2pi
        self.computeNewValue(val)


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
            self.computeNewValue(self.oneTurn/10)
        else:
            self.computeNewValue(-self.oneTurn/10)


    def get(self):
        if self.reportDelta:
            return self.type(self.deltaVal)
        else:
            return self.type(self.value)


    def printLabel(self):
        if not self.showLabel in [1,2]:
            return
        hl = self.valueLabel.winfo_reqheight()
        wl = self.valueLabel.winfo_reqwidth()
        h = self.canvas.winfo_reqheight()
        w = self.canvas.winfo_reqwidth()
        self.valueLabel.place(in_=self.canvas, x=(w-wl)*.5, y=((h-hl)*0.5)-
                              self.height/4)
        self.valueLabel.configure(text=self.labelFormat%self.value)


    def set(self, val, update=1, force=0):
        # if force is set to 1, we call this method regardless of the
        # widget configuration. This is for example the case if the thumwheel
        # is set to continuous=0, but the value is set in the options panel
        if self.min is not None and val <= self.min:
            val = self.min
        elif self.max is not None and val >= self.max:
            val = self.max
      
        oldval = self.value
        self.value = val
        self.deltaVal = self.value - oldval
        
        newVal = self.get()
        if update and self.continuous or force:
            if self.oldValue != self.value or force:
                self.oldValue = newVal
                self.callbacks.CallCallbacks(newVal)
            if self.showLabel==2:
                self.printLabel()
 
        if self.showLabel==1:
            self.printLabel()
            
        if self.opPanel:
            self.opPanel.valInput.set(self.labelFormat%self.value)
        
    #FIXME: this could be renamed increment (to parallel the set method)
    # some code in this method is duplicated in set method
    def computeNewValue(self, val):
        # snap to closest increment
        oldval = self.value
        if self.increment is not None and self.increment != 0.:
            self.discreteValue = self.discreteValue + val
            if self.discreteValue>=self.increment:
                self.value = self.value+self.increment
                self.discreteValue=0.
            elif -self.discreteValue>=self.increment:
                self.value = self.value-self.increment
                self.discreteValue=0.
        else:
            self.value = self.value + val

        if self.min is not None and self.value <= self.min:
            self.value = self.min
            val = 0
        elif self.max is not None and self.value >= self.max:
            self.value = self.max
            val = 0
        
        self.deltaVal = self.value - oldval

        # self.angle is used in drawLines()
        self.angle = val/self.oneTurnOver2pi
        self.drawLines()

        if self.showLabel>0: # ALWAYS
            self.printLabel()
            
        if self.opPanel:
            self.opPanel.valInput.set(self.labelFormat%self.value)

        newVal = self.get()
        if self.oldValue != newVal and not self.reportDelta:
            # this part is called either when the value is float OR when
            # continuous is off
            if self.continuous:
                self.oldValue = newVal
                self.callbacks.CallCallbacks(newVal)
            else:
                pass
        else:
            # this part is usually called when the datatype is INT AND
            # continuous is on
            if self.oldValue != newVal:
                # only call this part when we "reach" a new value
                self.oldValue = newVal
                self.callbacks.CallCallbacks(newVal)
            else:
                # else we need nothing to do
                pass


    def createCanvas(self, master, wheelPad=6):
        bw = self.borderWidth = wheelPad # distance between wheel and raise
        
        cd={'width':self.width+bw, 'height':self.height+bw, 'relief':'raised',
                'borderwidth':3}
        
        for k, w in self.canvasCfg.items():
            cd[k] = w

        self.canvas = Tkinter.Canvas(self, cd)
        cbdw = int(self.canvas.cget('borderwidth'))
        bd = self.borderWidth + cbdw + 1 # +1 for pixel0 that is not drawn

        height = self.height-bw
        width = self.width-bw
        cp = self.canvas.create_polygon
	self.outline1 = cp( bd, bd, bd, height+cbdw, width+cbdw, height+cbdw,
                            width+cbdw, bd, bd, bd,
                            width=1, outline='gray60', fill='gray85')

        ul = bd+1 # upper left pixel
        l = (width+cbdw-1) - (bd+1) # length of the inner box
        cl25 = 2.*l/25.
        cl = self.canvas.create_line
	self.outline2 = cl( ul+int(cl25), ul, ul, ul, ul,
                            height+cbdw-1, ul+int(cl25),
                            height+cbdw-1,
                            width=1, fill='gray20')
 	self.outline3 = cl( ul+int(cl25), ul, ul+int(3*cl25), ul,
                            width=1, fill='gray60')
 	self.outline4 = cl( ul+int(cl25), height+cbdw-1,
                            ul+int(3*cl25), height+cbdw-1,
                            width=1, fill='gray60')
 	self.outline5 = cl( ul+int(5*cl25), ul, ul+int(7.5*cl25), ul,
                            width=1, fill='white')
 	self.outline4 = cl(  ul+int(5*cl25), height+cbdw-1, ul+int(7.5*cl25),
                            height+cbdw-1,
                            width=1, fill='white')
 	self.outline6 = cl( ul+int(9.5*cl25), ul, ul+int(11.5*cl25), ul,
                            width=1, fill='gray60')
 	self.outline7 = cl( ul+int(9.5*cl25), height+cbdw-1,
                            ul+int(11.5*cl25), height+cbdw-1,
                            width=1, fill='gray60')

        re = ul+l
	self.outline8 = cl( ul+int(11.5*cl25), ul, re, ul, re,
                            height+cbdw-1, ul+int(11.5*cl25),
                            height+cbdw-1,
                            width=1, fill='gray20')

        # corners of the box where the lines have to be drawn
        self.innerbox = (ul+1, ul+1, re-1, height+cbdw-1)

        self.circlePtsAngles = []
        inc = 2*math.pi/float(self.nblines)
        for i in range(self.nblines):
            self.circlePtsAngles.append( i*inc )
        self.circlePtsAngles = Numeric.array(self.circlePtsAngles, 'f')
        self.linesIds = []
        self.shLinesIds = []
        # this table should depend on the number of lines
        # currently it is of length 15 (self.nblines/2)
        # It should be resized automatically to self.nblines/2
        self.shadowLinesOptions = [
            ( 0, 'black', 1), # offset, color, width
            ( 2, 'gray30', 1),
            ( 2, 'gray30', 1),
            ( 0, 'gray30', 1),
            (-1, 'white', 1),
            (-1, 'white', 1),
            (-1, 'white', 2),
            (-1, 'white', 2),
            (-1, 'white', 2),
            (-1, 'white', 2),
            (-1, 'white', 1),
            (-1, 'white', 1),
            ( 0, 'gray30', 1),
            (-2, 'gray30', 1),
            (-2, 'gray30', 1),
            ( 0, 'black', 1), # offset, color, width
            ( 0, 'black', 1), # offset, color, width
            ]

        for i in range(self.nblines):
            self.linesIds.append(cl(0,0,0,0,width=1, fill='black'))
            self.shLinesIds.append(cl(0,0,0,0))
        
        wlabCfg = {'padx':0,'pady':0}
        wlabCfg.update(self.wheelLabCfg)
        self.valueLabel = apply( Tkinter.Label, (self.master,), wlabCfg)

        self.drawLines()
	self.canvas.pack(side=Tkinter.LEFT)
        self.toggleWidgetLabel(self.showLabel)

        
    def drawLines(self):
        # angle has to be provided in radians
        angle = self.angle
            
        self.circlePtsAngles = self.circlePtsAngles+angle
        self.circlePtsAngles = Numeric.remainder(self.circlePtsAngles,
                                                 2*math.pi)
        xcoords = Numeric.cos(self.circlePtsAngles)
        xsin = Numeric.sin(self.circlePtsAngles)
        if self.orient=='horizontal':
            w = self.innerbox[2] - self.innerbox[0]
            r = w/2
            c = self.innerbox[0] + r
            y1 = self.innerbox[1]
            y2 = self.innerbox[3]
        else:
            w = self.innerbox[3] - self.innerbox[1]
            r = w/2
            c = self.innerbox[1] + r
            y1 = self.innerbox[0]
            y2 = self.innerbox[2]
            
        cl = self.canvas.create_line
        setCoords = self.canvas.coords
        setOpt = self.canvas.itemconfigure
        pi2 = math.pi/2.
        drawLines = 0
        co = Numeric.take(xcoords,
                          Numeric.nonzero(Numeric.greater_equal(xsin, 0.0)))
        co = Numeric.sort(co)
        co = [-1]+list(co)
        for i in range(len(co)):
            x = c - int(co[i]*r)
            if self.orient=='horizontal':
                setCoords(self.linesIds[i], x, y1, x, y2)
            else:
                setCoords(self.linesIds[i], y1, x, y2, x)
            shopt = self.shadowLinesOptions[i]
            x = x + shopt[0]
            if self.orient=='horizontal':
                setCoords(self.shLinesIds[i], x, y1, x, y2)
            else:
                setCoords(self.shLinesIds[i], y1, x, y2, x)
            setOpt(self.shLinesIds[i], fill = shopt[1], width=shopt[2])
            

    def toggleWidgetLabel(self, val):
        if val == 0:
            # no widget labels
            self.showLabel=0
            self.valueLabel.place_forget()
            
        if val == 1:
            # show always widget labels
            self.showLabel=1
            self.printLabel()

        if val == 2:
            # show widget labels only when mouse moves
            self.showLabel=2
            self.valueLabel.place_forget()

    def setValue(self, val):
        assert type(val) in [types.IntType, types.FloatType],\
               "Illegal type for value: expected %s or %s, got %s"%(
                   type(1), type(1.0), type(val) )

        # setValue does NOT call a callback!
        if self.min is not None and val < self.min: val = self.min
        if self.max is not None and val > self.max: val = self.max
        self.value = self.type(val)
        self.oldValue = self.value
        self.offsetValue=self.value
        if self.showLabel == 1:
            self.printLabel()
            

 #####################################################################
 # the 'configure' methods:
 #####################################################################

    def configure(self, **kw):
        if 'type' in kw.keys(): # make sure type is set first
            self.setType(kw['type'])
            del kw['type']
            
        for key,value in kw.items():
            # the 'set' parameter callbacks
            if key=='labCfg': self.setLabel(value)
            elif key=='callback': self.setCallbacks(value)
            elif key=='wheelLabCfg':
                self.wheelLablCfg.update(value)
                self.printLabel()
            elif key=='min': self.setMin(value)
            elif key=='max': self.setMax(value)
            elif key=='increment': self.setIncrement(value)
            elif key=='precision': self.setPrecision(value)
            elif key=='showLabel': self.setShowLabel(value)
            elif key=='continuous':
                self.setContinuous(value)
            elif key=='oneTurn': self.setOneTurn(value)
            elif key=='orient': self.setOrient(value)
            elif key=='reportDelta': self.setReportDelta(value)

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


        if self.valueLabel and self.showLabel == 1:
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
            self.increment = None
            if hasattr(self.opPanel, 'optionsForm'):
                self.opPanel.toggleIncr.set(0)
                self.opPanel.incr_entry.configure(state='disabled', fg='gray40')


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

        assert cont in [True, False, 0, 1],\
             "Illegal value for continuous: expecting a boolean True, False, got %s"%cont

        self.continuous = cont
        if hasattr(self.opPanel, 'optionsForm'):
            w=self.opPanel.idf.entryByName['togCont']['widget']
            if cont:
                w.setvalue('on')
            else:
                w.setvalue('off')


    def setShowLabel(self, val):
        """Show label can be 0, 1 or 2
0: no label
1: show label only when value changes
2: label is always shown"""
        
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
        

    def setOrient(self, orient):
        if orient is None:
            if self.width > self.height:
                orient='horizontal'
            else:
                orient = 'vertical'
        assert orient in ['horizontal', 'vertical'],\
               "Expected 'vertical' or 'horizontal', got '%s'"%orient
        self.orient = orient
        

    def setReportDelta(self, rD):
        assert rD in [None, 0, 1],\
               "Expected None, 0 or 1, got %s"%rD
        if rD is None or rD == 0:
            self.reportDelta = 0
        else:
            self.reportDelta = 1


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
    tw = ThumbWheel(width=100, height=20, labCfg={'text':'X:'},
#    tw = ThumbWheel(width=20, height=100, labCfg={'text':'X:'},
#                    orient='vertical',
                    wheelPad=2, type=float,
                    min=10)
    tw.configure(showLabel=1)
    tw.callbacks.AddCallback(foo)
