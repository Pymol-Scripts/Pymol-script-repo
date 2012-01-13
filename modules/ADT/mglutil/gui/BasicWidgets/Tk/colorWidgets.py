#############################################################################
#
# Author: Sophie COON, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

#$Header: /opt/cvs/python/packages/share1.5/mglutil/gui/BasicWidgets/Tk/colorWidgets.py,v 1.27.2.7 2011/04/26 20:34:41 sargis Exp $
#
#$Id: colorWidgets.py,v 1.27.2.7 2011/04/26 20:34:41 sargis Exp $
#

import Tkinter, Pmw, os
import numpy.oldnumeric as Numeric
from types import ListType, TupleType
from mglutil.util.callback import CallbackFunction, CallbackManager
from mglutil.util.colorUtil import ToRGB, ToHSV, ToHEX
from mglutil.gui.InputForm.Tk.gui import InputFormDescr,InputForm,evalString
from mglutil.gui.BasicWidgets.Tk.customizedWidgets import ExtendedSliderWidget
from mglutil.gui.BasicWidgets.Tk.customizedWidgets import ListChooser
from mglutil.gui.BasicWidgets.Tk.fileBrowsers import fileOpenAsk, fileSaveAsk
import tkMessageBox
import os, math
from mglutil.util.packageFilePath import getResourceFolder

class ColorWheel:        
    def __init__(self, master, title=None, callback=None, immediate=1):
        if not master:
            master = Tkinter.Toplevel()

        if title is not None:
            master.title(title)

        f = self.frame = Tkinter.Frame(master)
        path = __import__('mglutil').__path__
        iconfile = os.path.join(path[0],'gui/BasicWidgets/Tk/cw.ppm')
        self.iconfile = iconfile
        self.cwim = Tkinter.PhotoImage(file=iconfile, master=master)
        self.width = self.cwim.width()
        self.height = self.cwim.height()
        self.cwcanvas = Tkinter.Canvas(f, width=self.width,
                                       height=self.height, ###relief='sunken',
                                       borderwidth=3 )
        self.cwcanvas.create_image(3, 3, anchor=Tkinter.NW, image=self.cwim)
        self.cwcanvas.pack()
        self.frame.pack()

        #self.callback = None
        self.cbManager = CallbackManager()
        if callback:
            if type(callback) in [ListType, TupleType]:
                map(self.cbManager.AddCallback, callback)
            else:
                self.cbManager.AddCallback(callback)
        self.afterID = None
        self.immediate = immediate
        self.x = 0
        self.y = 0
        self.radius = 55
        cx = self.cx = self.width/2 + 3
        cy = self.cy = self.height/2 + 3

        self.cursor = self.cwcanvas.create_line(
            
            cx-3, cy-3, cx-3, cy+3, cx+3,cy+3, cx+3, cy-3, cx-3, cy-3 )

        self.hsvColor = [1.,1.,1.]
        
        self.cwcanvas.bind('<ButtonPress-1>', self.mouse1Down)
        #self.cursor = self.cwcanvas.create_line(cx, cy, cx-55, cy)


    def _MoveCursor(self, x, y, trigger=1):
	# find the saturation based on distance
	s = math.sqrt(x*x + y*y) / self.radius
	if s > 1.0:
	    x = x / s
	    y = y / s
	    s = 1.0

	# now find the hue based on the angle 
	if x or y:
	    angle = math.atan2(y, x)
	    if angle < 0.0:
                angle = angle + (2.*math.pi)
	    h = 1. - angle / (2.0 * math.pi)
	else:
	    h = 0
	# check if redraw and callback are needed
        if self.hsvColor[0] != h or self.hsvColor[1] != s:
            if trigger==1:
                self.hsvColor[0] = h
                self.hsvColor[1] = s
            cx = self.cx+x
            cy = self.cy+y
            self.cwcanvas.coords( self.cursor, cx-3, cy-3, cx-3, cy+3,
                                  cx+3,cy+3, cx+3, cy-3, cx-3, cy-3 )

	    if self.cbManager.callbacks:
                self.cbManager.CallCallbacks(self.get('RGB'))

    def mouse1Down(self, event=None):
        self.cwcanvas.bind('<B1-Motion>', self.mouse1Move)
        self.cwcanvas.bind('<ButtonRelease-1>', self.mouse1Up)
        self._MoveCursor(event.x - self.cx, event.y - self.cy)

        
    def mouse1Move(self, event=None):
        if self.immediate:
            if self.afterID is not None:
                self.cwcanvas.after_cancel(self.afterID)
                self.afterID = None
            else:
                self.afterID = self.cwcanvas.after(15,self._MoveCursor ,
                                                   event.x - self.cx,
                                                   event.y - self.cy)
            
        else:
            self._MoveCursor(event.x - self.cx, event.y - self.cy,
                             trigger=0)

        
    def mouse1Up(self, event=None):
        self._MoveCursor(event.x - self.cx, event.y - self.cy)
        self.cwcanvas.unbind('<B1-Motion>')
        self.cwcanvas.unbind('<ButtonRelease-1>')


    def get(self, mode='HSV'):
	"""Get the current color"""
	if mode == 'RGB':
	    rgb = ToRGB(self.hsvColor)
	    #return OneColor(rgb)
            return rgb

	elif mode == 'HSV':
	    #return OneColor(self.hsvColor)
            return self.hsvColor

	elif mode == 'HEX':
	    col = Numeric.array(ToRGB(self.hsvColor[:]), 'f') * 255
	    return ToHEX(col)


    def set(self, color, mode='HSV', trigger=1):
	"""Set the current color"""
	assert len(color)==3
	#color = OneColor(color)
	if mode=='RGB': color = ToHSV(color[:])
	self.hsvColor = list(color[:])

        # update cursor
	rad = self.hsvColor[1] * self.radius
	angle = 2.0 * math.pi * (1. - self.hsvColor[0])
	cx = self.cx + int(rad * math.cos(angle))
	cy = self.cy + int(rad * math.sin(angle))
        self.cwcanvas.coords( self.cursor, cx-3, cy-3, cx-3, cy+3,
                              cx+3,cy+3, cx+3, cy-3, cx-3, cy-3 )

        if trigger==1 and self.immediate and self.cbManager.callbacks:
            self.cbManager.CallCallbacks(self.get('RGB'))


class ColorEditor:
    """
    The ColorEditor is a widget providing a colorwheel, a value scale,
    HSV entries, RGB entries and HEX(HexTriplet) entries
    """
    def __init__(self, master=None, currentColor=(1.0,1.0,1.0), mode='RGB',
                 commands = None, immediate=1):
        assert mode in ['RGB', 'HSV','HEX']
        self.mode = mode
        if not master:
            self.master = Tkinter.Toplevel()
        else:
            self.master = master
        self.afterID = None
        self.immediate = immediate
        # The editFrame is the main Frame of the widget
        self.editFrame = Tkinter.Frame(self.master, borderwidth=2,
                                      relief='ridge')
        self.cbManager = CallbackManager()
        if commands:
            if type(commands) in [ListType, TupleType]:
                map(self.cbManager.AddCallback, commands)
            else:
                self.cbManager.AddCallback(commands)

        if mode == 'HSV':
            self.currentHSV = list(currentColor)
            self.currentRGB = list(ToRGB(currentColor))
            self.currentHEX = ToHEX(currentColor, mode = 'HSV')

        elif mode == 'RGB':
            self.currentRGB = list(currentColor)
            self.currentHSV = list(ToHSV(currentColor))
            self.currentHEX = ToHEX(currentColor, mode='RGB')

        elif mode == 'HEX':
            self.currentRGB = ToHEX(currentColor, mode='RGB')
            self.currentHSV = ToHEX(currentColor, mode='HSV')
            self.currentHEX = currentColor
        else:
            print 'mode not recognized mode set to RGB'
        self.createColorEditor()


    def createColorEditor(self):
        # The chooserFrame containinsg the colorwheel and the value scale
        chooserFrame = Tkinter.Frame(self.editFrame)
        # Create a Tkinter Scale widget bound a callback : self.scale_cb
        self.vScale = Tkinter.Scale(chooserFrame,
                                    from_ = 1.0, to_ = 0.0,
                                    orient='vertical', resolution=0.01,)
##                                     command = self.scale_cb)
        if self.immediate:
            self.vScale.configure(command=self.scale_cb)
        else:
            self.vScale.bind('<ButtonRelease-1>', self.scaleUp_cb)

        # pack the scale to be on the left side of the colorwheel
        self.vScale.pack(side='right', padx=5, pady=5, expand=1,
                         fill='both')
        self.vScale.set(1.0)

        # Create the colorWheel Frame
        wheelFrame = Tkinter.Frame(chooserFrame, relief='ridge',
                                   borderwidth=1)
        # Pack the ColorWheel
        wheelFrame.pack(side='left',pady=5, padx=10,fill='both',
                        expand = 1)
        # Create the ColorWheel
        
        
        # to silent an error report from pychecker
        self.cw = ColorWheel(wheelFrame,None,self.colorWidget_cb,self.immediate)
        #self.cw = ColorWheel(wheelFrame,title=None, 
        #                     callback=self.colorWidget_cb, 
        #                     immediate=self.immediate)
                             
                             
        self.cw.set(self.currentRGB, mode='RGB', trigger=0)
        # Bind the colorwidget to a callback function self.colorWidget_cb
        #self.cw.callback = self.colorWidget_cb
        # Pack the chooserFrame
        chooserFrame.pack(expand=1,fill='both')

        bottomFrame = Tkinter.Frame(self.editFrame)
        #The preview frame will contain the frame to show the choosen color
        previewFrame = Tkinter.Frame(bottomFrame)
        previewFrame.pack(side='left', fill='both', expand=1)
        preview = Tkinter.Frame(previewFrame,)
        bg = self.currentHEX
        self.chip = Tkinter.Frame(previewFrame, 
                                  borderwidth=3, width=50,
                                  height=30, bg=bg, relief='ridge')
        # Pack the chipFrame
        self.chip.pack(fill='both', expand = 1)
        

        #The entriesFrame will contain all the entryFields
        entriesFrame = Tkinter.Frame(bottomFrame)
        entriesOption = {'labelpos':'w',
                         'validate':{'validator':'real',
                                     'min':0.0, 'max':1.0},
                         'entry_width':4,
                         }
                         
        # the hsvFrame will contain the H,S,V entryFields
        hsvFrame = Tkinter.Frame(entriesFrame)
        
        entriesOption['label_text'] = 'H'
        entriesOption['value'] = "%4.2f"%self.currentHSV[0]
        entriesOption['command'] = self.hVal_cb
        self.hVal = apply(Pmw.EntryField, (hsvFrame,), entriesOption)
        self.hVal.pack(side = 'left')

        entriesOption['label_text'] = 'S'
        entriesOption['value'] = "%4.2f"%self.currentHSV[1]
        entriesOption['command'] = self.sVal_cb
        self.sVal = apply(Pmw.EntryField, (hsvFrame,), entriesOption)
        self.sVal.pack(side = 'left')

        entriesOption['label_text'] = 'V'
        entriesOption['value'] = "%4.2f"%self.currentHSV[2]
        entriesOption['command'] = self.vVal_cb
        self.vVal= apply(Pmw.EntryField, (hsvFrame,), entriesOption)
        self.vVal.pack(side = 'left')

        hsvFrame.pack(padx=4, pady=4,fill='both',expand=1)

        rgbFrame = Tkinter.Frame(entriesFrame)
        # RGB entries
        entriesOption['label_text'] = 'R'
        entriesOption['value'] = "%4.2f"%self.currentRGB[0]
        entriesOption['command'] = self.rVal_cb
        self.rVal = apply(Pmw.EntryField, (rgbFrame,), entriesOption)
        self.rVal.pack(side = 'left')

        entriesOption['label_text'] = 'G'
        entriesOption['value'] = "%4.2f"%self.currentRGB[1]
        entriesOption['command'] = self.gVal_cb
        self.gVal = apply(Pmw.EntryField, (rgbFrame,), entriesOption)
        self.gVal.pack(side = 'left')

        entriesOption['label_text'] = 'B'
        entriesOption['value'] = "%4.2f"%self.currentRGB[2]
        entriesOption['command'] = self.bVal_cb
        self.bVal = apply(Pmw.EntryField, (rgbFrame,), entriesOption)
        self.bVal.pack(side = 'left')
        rgbFrame.pack(padx=4, pady=4,fill='both',expand=1)

        hexFrame = Tkinter.Frame(entriesFrame)
        entriesOption['label_text'] = 'Hex triplet'
        entriesOption['value'] = self.currentHEX
        entriesOption['command'] = self.hexVal_cb
        del entriesOption['validate']
        entriesOption['entry_width']=8
        #entriesOption['validate']='alphanumeric'
        self.hexVal = apply(Pmw.EntryField, (hexFrame,), entriesOption)
        self.hexVal.pack(padx=4, pady=4,side = 'left')
        hexFrame.pack(fill='both',expand=1)

        entriesFrame.pack(side = 'right', fill='both',expand=1)
        bottomFrame.pack(side='bottom',fill='both', expand=1)


    ###############################################################
    ####           COLOR CHOOSER UTILITY FUNCTIONs             ####
    ###############################################################
    #def caluculate(self):
        
    def set(self, color, mode='RGB', trigger=1):
	"""Set the current color"""
	#assert len(color)==3
        assert mode in ['HSV', 'RGB', 'HEX']
        self.mode = mode
        if mode == 'HSV':
            newRGB = map(lambda x: float("%4.2f"%x), ToRGB(color))
        elif mode == 'HEX':
            newRGB = ToRGB(color, mode='HEX')
        else: newRGB = color
        if newRGB != self.currentRGB:
            self.updateWidgetsColor(newRGB)
            if trigger==1 and self.immediate and self.cbManager.callbacks:
                self.cbManager.CallCallbacks(self.currentRGB)

    def get(self, mode = 'RGB'):
        assert mode in ['RGB','HSV', 'HEX']
        self.mode = mode
        if mode == 'RGB':
            return self.currentRGB
        elif mode == 'HSV':
            return self.currentHSV
        elif mode == 'HEX':
            col = ToHEX(self.currentRGB, mode='HEX')
            return col

    def pack(self,*args, **kw):
        apply(self.editFrame.pack, args, kw)

    def pack_forget(self,*args, **kw):
        apply(self.editFrame.pack_forget, args, kw)
                         
    def grid(self,*args, **kw):
        apply(self.editFrame.grid, args, kw)

    def grid_forget(self,*args, **kw):
        apply(self.editFrame.grid_forget, args, kw)

    ###############################################################
    ####               WIDGETS CALLBACK FUNCTIONS              ####
    ###############################################################
        
    def colorWidget_cb(self, rgbcolor):
        # Do this test because updateCurrent is called after the set
        # cw.
        # if color is different from current color then update chip
        color=list(ToHSV(rgbcolor))[:]
        color[2] = float(self.vScale.get())
        newrgb = list(ToRGB(color))
        if newrgb != self.currentRGB:
            self.updateWidgetsColor(newrgb, who='cw')

    def scale_cb(self, val):
        if self.afterID is not None:
            self.vScale.after_cancel(self.afterID)
            self.afterID = None
        else:
            self.afterID = self.vScale.after(17, self.scaleImm_cb, val)
            
        
    def scaleImm_cb(self, val):
        newHSV = [float(self.hVal.get()),float(self.sVal.get()), float(val)]
        if newHSV != self.currentHSV:
            self.updateWidgetsColor(ToRGB(newHSV), who='scale')
        
    def scaleUp_cb(self, event=None):
        val = float(self.vScale.get())
        newHSV = [float(self.hVal.get()),float(self.sVal.get()), float(val)]
        if newHSV != self.currentHSV:
            self.updateWidgetsColor(ToRGB(newHSV), who='scale')

    def hVal_cb(self):
        val = float(self.hVal.get())
        newColor = self.currentHSV
        newColor[0] = val
        newHSV = map(lambda x: float("%4.2f"%x), newColor)
        if (not (float(self.vVal.get())==0.00 or \
                (float(self.sVal.get())==0 and float(self.vVal.get())==1.0))):
            self.updateWidgetsColor(ToRGB(newHSV), who='h')

    def sVal_cb(self):
        val = float(self.sVal.get())
        newColor = list(ToHSV(self.currentRGB[:]))
        newColor[1] = val
        newHSV = map(lambda x: float("%4.2f"%x), newColor)
        if (not float(self.vVal.get())==0) and newHSV != self.currentHSV:
            self.updateWidgetsColor(ToRGB(newHSV), who='s')

    def vVal_cb(self):
        newColor = [float(self.hVal.get()),
                    float(self.sVal.get()),
                    float(self.vVal.get())]
        newHSV = map(lambda x: float("%4.2f"%x), newColor)
        if newHSV != self.currentHSV:
            self.updateWidgetsColor(ToRGB(newHSV), who='v')

    def rVal_cb(self):
        val = float(self.rVal.get())
        newColor = self.currentRGB[:]
        newColor[0] = val
        newRGB = map(lambda x: float("%4.2f"%x), newColor)
        if newRGB != self.currentRGB:
            self.updateWidgetsColor(newRGB, who='r')


    def gVal_cb(self):
        val = float(self.gVal.get())
        newColor = self.currentRGB[:]
        newColor[1] = val
        newRGB = map(lambda x: float("%4.2f"%x), newColor)
        if newRGB != self.currentRGB:
            self.updateWidgetsColor(newRGB, who='g')

    def bVal_cb(self):
        val = float(self.bVal.get())
        newColor = self.currentRGB[:]
        newColor[2] = val
        newRGB = map(lambda x: float("%4.2f"%x), newColor)
        if newRGB != self.currentRGB:
            self.updateWidgetsColor(newRGB, who='b')


    def hexVal_cb(self):
        val = self.hexVal.get()
        if val[0] !='#' or len(val)!=7:
            val = self.currentHEX
        newRGB = ToRGB(val, 'HEX')
        if newRGB != self.currentRGB:
            self.updateWidgetsColor(newRGB, who='hex')
        
    ###############################################################
    ####               WIDGETS UPDATE FUNCTIONS                ####
    ###############################################################
    def updateWidgetsColor(self, rgbcolor, who = 'set', trigger=1):
        oldRGB = list(self.currentRGB)
        self.currentRGB = map(lambda x: float("%4.2f"%x), rgbcolor)
        # If newcolor is the same than old color nothing to update.
        if oldRGB == self.currentRGB : return
        hsvcolor = ToHSV(rgbcolor[:])
        self.currentHSV = map(lambda x: float("%4.2f"%x), hsvcolor)
        self.currentHEX = ToHEX(self.currentRGB)
        # Update the preview chip
        self.chip.configure( bg = self.currentHEX )

        # ColorWidget:
        cwColor = self.cw.get(mode='RGB')
        newRGB = map(lambda x: float("%4.2f"%x),cwColor)
        if newRGB != self.currentRGB and not who in ['v', 'scale']:
            self.cw.set(self.currentRGB, mode='RGB',trigger=0)

        # Value Scale:
        scaleCol = self.vScale.get()
        if scaleCol != self.currentHSV[2]:
            self.vScale.set(self.currentHSV[2])

        # H Entry:
        h = float(self.hVal.get())
        hCol = float("%4.2f"%h)
        if hCol != self.currentHSV[0] and not who in ['v', 'scale']:
            self.hVal.setentry(self.currentHSV[0])

        # S Entry:
        s = float(self.sVal.get())
        sCol = float("%4.2f"%s)
        if sCol != self.currentHSV[1] and self.currentHSV[2] !=0 \
           and not who in ['v', 'scale']:
            self.sVal.setentry(self.currentHSV[1])

        # V Entry:
        v = float(self.vVal.get())
        vCol = float("%4.2f"%v)
        if vCol != self.currentHSV[2]and self.currentHSV[2] !=0:
            self.vVal.setentry(self.currentHSV[2])

        # R Entry:
        r = float(self.rVal.get())
        rCol = float("%4.2f"%r)
        if rCol != self.currentRGB[0]:
            self.rVal.setentry(self.currentRGB[0])
        # G Entry:
        g = float(self.gVal.get())
        gCol=float("%4.2f"%g)
        if gCol != self.currentRGB[1]:
            self.gVal.setentry(self.currentRGB[1])
        
        # B Entry:
        b = float(self.bVal.get())
        bCol = float("%4.2f"%b)
        if bCol != self.currentRGB[2]:
            self.bVal.setentry(self.currentRGB[2])
        
        # Hex Entry:
        hexCol = self.hexVal.get()
        if hexCol != self.currentHEX:
            self.hexVal.setentry(self.currentHEX)

        # This might depend of the mode. ?
        if trigger==1 and self.immediate and self.cbManager.callbacks:
            self.cbManager.CallCallbacks(self.currentRGB)



class Chooser(object):


    def __init__(self, master=None, title = 'Chooser', commands = None,
                 immediate=0, exitFunction=None):
        if master is None:
            self.master = Tkinter.Toplevel()
            self.ownmaster=1
            self.master.title(title)
            #self.master.protocol('WM_DELETE_WINDOW', self.dismiss)
        else:
            self.ownmaster=0
            self.master = master
        # The editFrame is the main Frame of the widget
        self.masterFrame = Tkinter.Frame(self.master,
                                         borderwidth=2, relief='ridge')
        self.immediate=immediate
        # Create a cbManager
        self.cbManager = CallbackManager()
        self.createCommon()
        self.createChooser()
        if commands:
            if type(commands) in [ListType, TupleType]:
                map(self.cbManager.AddCallback, commands)
                map(self.ce.cbManager.AddCallback, commands)
            else:
                self.cbManager.AddCallback(commands)
                self.ce.cbManager.AddCallback(commands)
    

    def createCommon(self):
        # Create the Menu Bar
        self.menuBar = Pmw.MenuBar(self.masterFrame,
                                   hull_relief = 'raised',
                                   hull_borderwidth = 1)
        self.menuBar.addmenu('File', 'Close this window or exit')
        
        self.mainFrame = Tkinter.Frame(self.masterFrame,
                                       borderwidth=2, relief='ridge',
                                        width=150, height=200)
        
        self.menuBar.pack(fill = 'x')

        ## Create the ColorEditor
        self.ce = ColorEditor(self.mainFrame, immediate=self.immediate)
        
        self.hidden=1

        self.mainFrame.pack(fill='both', expand=1)

    ######################################################################
    ####         UTILITY FUNCTIONS                                    ####
    ######################################################################


    def createChooser(self):
        pass

    def pack(self, *args, **kw):
        apply(self.masterFrame.pack, args, kw)

    def pack_forget(self, *args, **kw):
        apply(self.masterFrame.pack_forget, args, kw)

    def grid(self ,*args, **kw):
        apply(self.masterFrame.grid, args, kw)

    def grid_forget(self, *args, **kw):
        apply(self.masterFrame.grid_forget, args, kw)
        
colorChoosersDict = {} #colorChoosersDict stores previous ColorChoosers

class ColorChooser(Chooser):

    colors = [
        '#FF8284', '#ffff84', '#84ff84', '#00ff84', '#84ffff', '#0082ff', '#ff82c6', '#ff82ff',
        '#ff0000', '#ffff00', '#84ff00', '#00ff42', '#00ffff', '#0082c6', '#8482c6', '#ff00ff',
        '#844142', '#ff8242', '#00ff00', '#008284', '#004184', '#8482ff', '#840042', '#ff0084',
        '#840000', '#ff8200', '#008200', '#008242', '#0000ff', '#0000a5', '#840084', '#8400ff',
        '#420000', '#844100', '#004100', '#004142', '#000084', '#000042', '#420042', '#420084',
        '#000000', '#848200', '#848242', '#848284', '#428284', '#c6c3c6', '#6b0c94', '#ffffff',
        '#FFFFFF', '#FFFFFF', '#FFFFFF', '#FFFFFF', '#FFFFFF', '#FFFFFF', '#FFFFFF', '#FFFFFF',
        '#FFFFFF', '#FFFFFF', '#FFFFFF', '#FFFFFF', '#FFFFFF', '#FFFFFF', '#FFFFFF', '#FFFFFF',
        '#FFFFFF', '#FFFFFF', '#FFFFFF', '#FFFFFF', '#FFFFFF', '#FFFFFF', '#FFFFFF', '#FFFFFF'
    ]

    def __new__(cls, master=None, title = 'Chooser', commands = None,
                 immediate=0, exitFunction=None):

        #the following code was added to avoid having multiple copies of ColorChoosers
        global colorChoosersDict
        if  colorChoosersDict.has_key(title):
            try:
                colorChoosersDict[title].master.destroy()
            except:
                pass

        self= object.__new__(cls)
        self.exitFunc = exitFunction
        self.mapping = {}
        self.currentTag = None
        self.rcPath = getResourceFolder()
        self.customColorsPath = os.path.join(self.rcPath,'customColors')
        if os.path.exists(self.customColorsPath):
            customColors = open(self.customColorsPath).read().strip()
            customCols = customColors.split()
            self.colors[48:48+len(customCols)] = customCols
        
        Chooser.__init__(self, master=master, title=title, commands = commands,
                         immediate=immediate, exitFunction=exitFunction)

        if self.ownmaster:
            ## create dismiss button
            if exitFunction:
                cb = exitFunction
            else:
                cb = self.master.destroy
            self.dismissb = Tkinter.Button(self.master, text='DISMISS',
                                           command=cb)
            self.dismissb.pack(side='bottom', expand=1, fill='x')
        
        try:
            self.master.protocol('WM_DELETE_WINDOW', self.exit)
        except:
            pass 
        
        colorChoosersDict[title] = self
        return self
            

    def __init__(self, master=None, title = 'Chooser', commands = None,
                 immediate=0, exitFunction=None):
        pass
                
            
    def createChooser(self):
        self.ccFrame = Tkinter.Frame(self.mainFrame)
        self.menuBar.forget()

##         self.menuBar.addmenuitem('File', 'command', 'Load custom',
##                                  command = self.load,
##                                  label='Load')

##         self.menuBar.addmenuitem('File', 'command', 'Save Colors',
##                                  command = self.save_cb,
##                                  label='Save')

##         self.menuBar.addmenuitem('File', 'separator')

##         self.menuBar.addmenuitem('File', 'command', 'Exit',
##                                  command=self.exit,
##                                  label='Exit')
        
##         self.menuBar.addmenu('Edit', 'editing commands')

##         self.menuBar.addmenuitem('Edit', 'command','Add New Color',
##                                  command = self.addColor,
##                                  label='Add New Color')
##         self.add = 0

##         self.menuBar.addmenuitem('Edit', 'command','Edit Custom Color',
##                                  command = self.editColor,
##                                  label='Edit Selected Color')
        self.edit = 0

##         self.menuBar.addmenuitem('Edit', 'command','Hide Color Editor',
##                                  command = self.hideCE,
##                                  label='Hide Color Editor')
        # Here we are creating the part of the widget that will contain the
        # chips.
        # ccFrame is the left part of the widget
        self.ccFrame.pack(side='left',expand=1, fill='both')
        self.addButton = Tkinter.Button(self.ccFrame, text='Add to custom',
                                        command=self.addToCustom_cb)
        self.addHidden = 1

        # The chips frame will contain the color chips which are RadioSelect
        # and it is a scrolledFrame.
        chipsSFrame = Pmw.Group(self.ccFrame, tag_text='Basic colors:')
        #chipsSFrame = Pmw.ScrolledFrame(self.ccFrame, usehullsize=1,
        #                                hull_width=130,
        #                                hull_height=200,
        #                                hscrollmode = 'none')
        chipsSFrame.pack(padx = 5, pady = 3)#, fill = 'both', expand = 1)
        
        # Create the RadioSelect empty
        self.chipsFrame = chipsSFrame.interior()

        self.colorVar = colorVar = Tkinter.IntVar(0)
        colors = self.colors
        for i in range(6):
            for j in range(8):
                val = i*8+j
                col = colors[val]
                b = Tkinter.Radiobutton(
                    self.chipsFrame, text="", variable=colorVar, value=val,
                    bg = col, activebackground=col,
                    fg = col, activeforeground = col,
                    indicatoron=0, selectcolor=col,
                    width=2, height=1, command=self.selectColor)
                b.grid(row=i, column=j)

        customFrame = Pmw.Group(self.ccFrame, tag_text='Custom colors')
        self.customchipsFrame = customFrame.interior()
        self.custombuttons = []
        for i in range(3):
            for j in range(8):
                val = 48 + i*8+j
                col = colors[val]
                b = Tkinter.Radiobutton(
                    self.customchipsFrame, text="", variable=colorVar, value=val,
                    bg=col, activebackground=col,
                    fg=col, activeforeground=col,
                    indicatoron=0, selectcolor=col,
                    width=2, height=1, command=self.selectColor)
                b.grid(row=i, column=j)
                self.custombuttons.append(b)
                
        customFrame.pack(padx = 5, pady = 3)#, fill = 'both', expand = 1)

        #self.colorChips=Pmw.RadioSelect(self.chipsFrame,
        #                                label_text='Default Colors',
        #                                labelpos='nw', orient='vertical',
        #                                buttontype='radiobutton',
        #                                command = self.colButton_cb)

        #self.mod = {}
        #execfile(self.customFilename, self.mod)
        #self.cFlag=0
        #self.addCustomCol(paletteName = self.colorsName)
        #self.doubleClick = False
        #self.editColor()
        #self.ce.pack(side = 'right', fill='both', expand=1)
        self.ce.immediate=1
        self.ce.cw.immediate=1
        self.currentEditingCB = None

        
    def selectColor(self, event=None):
        colNum = self.colorVar.get()

        hcol = self.colors[colNum]
        rgb = int(hcol[1:3], 16), int(hcol[3:5], 16), int(hcol[5:7], 16)
        col = [x/255. for x in rgb]

        if colNum > 47: # custom color
            self.ce.pack(side = 'right', fill='both', expand=1)
            if self.currentEditingCB:
                self.ce.cbManager.RemoveCallback(self.currentEditingCB)
            cb = CallbackFunction(self.editCustom, custColNum=colNum-48)
            self.ce.cbManager.AddCallback(cb)
            self.ce.set(col, trigger=0)
            self.currentEditingCB = cb
        else:
            self.ce.pack_forget()
            
        self.cbManager.CallCallbacks(col)
        

    def editCustom(self, col, custColNum=0):
        hexcol = ToHEX(col)
        self.custombuttons[custColNum].configure(
            bg=hexcol, activebackground=hexcol,
            fg=hexcol, activeforeground=hexcol, selectcolor=hexcol)
        self.colors[48+custColNum] = hexcol
        self.cbManager.CallCallbacks(col)
        try:
            outFile = open(self.customColorsPath, 'w')
        except Exception, inst: 
            print inst
            print "Can't save custom colors in ", self.customColorsPath
            return
        outFile.write(' '.join(self.colors[48:]))

    ## def load(self):
    ##     ftypes = [ ('Python files', '*.py') ]
    ##     filename = fileOpenAsk(self.master, types=ftypes,
    ##                            title='Load custom colors' )
    ##     # Open the module
    ##     if filename is None: return
    ##     self.customFilename = filename
    ##     self.mod = {}
    ##     execfile( self.customFilename, self.mod)
    ##     self.mod.keys()
    ##     colName = filter(lambda x: x[:2]!='__',self.mod.keys())
    ##     entries = map(lambda x: (x, None), colName)
    ##     # From the module display the available colorPalette.
    ##     self.showChooser(entries)
                                 
    def showChooser(self, entries):
        if self.cFlag == 1:
            self.palChooser.clear()
            map(self.palChooser.add, entries)
            self.root.deiconify()
        else:
            self.root = Tkinter.Toplevel()
            self.chooserFrame = Tkinter.Frame(self.root)
            self.palChooser = ListChooser(self.chooserFrame, mode = 'extended',
                                         title='Customized colors groups',
                                          entries = entries,
                                          command=self.addCustomCol,)
            self.cFlag=1
            dismissChooser = Tkinter.Button(self.chooserFrame,
                                            text='Dismiss',
                                            command=self.root.withdraw)
            self.palChooser.pack()
            dismissChooser.pack()
            self.chooserFrame.pack()
            
    ## def hideCE(self):
    ##     if self.hidden == 0:
    ##         self.ce.pack_forget()
    ##         self.hidden=1
    ##     if self.addHidden == 0:
    ##         self.addButton.pack_forget()
    ##         self.addHidden=1
    ##     self.add =0
    ##     self.edit=0
        

    ## def addCustomCol(self, event=None, paletteName=None):
    ##     if paletteName is None:
    ##         paletteName = self.palChooser.get()[0]
    ##         self.colorsName = paletteName
    ##     # first clean up what is there:
    ##     if not self.mod.has_key(paletteName):
    ##         self.colDict={}
    ##         return
    ##     else:
    ##         self.colDict = self.mod[paletteName]
    ##     if len(self.colorChips._buttonList)!=0:
    ##         self.colorChips.deleteall()

    ##     self.colorChips.configure(label_text=paletteName)
    ##     items = self.colDict.items()
    ##     items.sort()
    ##     for name, value in items:
    ##         fg = col = ToHEX(value)
    ##         try:
    ##             int(name)
    ##         except:
    ##             fg = 'black'
    ##         self.colorChips.add(name, bg = col,
    ##                             activebackground=col,
    ##                             fg = fg, activeforeground = col,
    ##                             indicatoron=0,selectcolor=col,
    ##                             width=10,height=1,
    ##                             value = name)
    ##     self.colorChips.pack(fill='x', expand=1)
    ##     self.ce.cbManager.AddCallback(self.editChip)
        
    ##     if hasattr(self, 'chooserFrame'):
    ##         self.chooserFrame.master.withdraw()

    ## def save_cb(self, fileName = None, paletteName=None):
    ##     """Save the color palette """
    ##     if paletteName is None or fileName is None:
    ##         if hasattr(self, 'saveHidden'):
    ##             if self.saveHidden == 1:
    ##                 self.root.deiconify()
    ##                 self.saveHidden = 0
    ##         else:
    ##             self.root = Tkinter.Toplevel()
    ##             self.saveFrame = Tkinter.Frame(self.root, )
    ##             groupFrame = Tkinter.Frame(self.saveFrame)
    ##             self.groupEntry = Pmw.EntryField(groupFrame,
    ##                                              label_text='Group name:',
    ##                                              labelpos='w',
    ##                                              value=self.colorsName)
    ##             label = Tkinter.Label(groupFrame, text="\t\t",
    ##                                   )
    ##             self.groupEntry.pack(side='left')
    ##             label.pack(side='right', fill='x', expand=1)
    ##             groupFrame.pack(fill='x', expand=1)

    ##             fileFrame = Tkinter.Frame(self.saveFrame)

    ##             self.fileEntry = Pmw.EntryField(fileFrame, 
    ##                                             label_text='Python File name:',
    ##                                             labelpos='w',
    ##                                             value=self.customFilename)

    ##             browseBut = Tkinter.Button(fileFrame, text='Browse',
    ##                                        command=self.browse)
    ##             self.fileEntry.pack(side='left')
    ##             browseBut.pack(side='right', fill='x', expand=1)
    ##             fileFrame.pack(fill='x', expand=1)

    ##             buttonFrame = Tkinter.Frame(self.saveFrame)
    ##             ok = Tkinter.Button(buttonFrame, text='OK', command=self.ok)
    ##             cancel = Tkinter.Button(buttonFrame, text='Cancel',
    ##                                     command=self.cancel)
    ##             ok.pack(side='left', fill='x',expand=1)
    ##             cancel.pack(side='right', fill='x', expand=1)
    ##             buttonFrame.pack(side='bottom', fill='both', expand=1)
    ##             self.saveFrame.pack(fill='both', expand=1)
                
    ##             self.saveHidden = 0

    def ok(self):
        filename= self.fileEntry.get()
        groupname = self.groupEntry.get()
        if not groupname:
            print 'ERROR'
        self.save(filename, groupname)
        self.root.withdraw()
        
    def save(self, filename, groupname):
        if filename is None: return
        if not os.path.isfile(filename):
            f = open(filename, 'w')
            s = groupname+'='+repr(self.colDict)
            f.write(s)
            f.write('\n')
            f.close()
        else:
            mod = {}
            execfile(filename, mod)
            colName = filter(lambda x: x[:2]!='__',dir(mod.keys()))
            if groupname in colName:
                f = open(filename, 'w')
                for name in colName:
                    if name == groupname:
                        s = groupname+'='+repr(self.colDict)
                    else:
                        dict = getattr(self.mod,name)
                        s = name+'='+repr(dict)
            else:
                f = open(filename, 'w')
                f.write('\n')
                s = groupname+'='+repr(self.colDict)
            f.write(s)
            f.write('\n')
            f.close()
                
    def cancel(self):
        self.root.withdraw()
        self.saveHidden=1

    def browse(self):
        ftypes = [ ('Python files', '*.py') ]

        filename = fileSaveAsk(self.master,  types=ftypes,
                               title='Save custom colors')
        if filename:
            self.fileEntry.setentry(filename)

    def hide(self):
        if hasattr(self.masterFrame.master,'withdraw'):
            try:
                self.masterFrame.master.withdraw()
            except:
                pass
        
    ## def exit(self):
    ##     self.hideCE()
    ##     if self.exitFunc is None:
    ##         if hasattr(self.masterFrame.master,'withdraw'):
    ##             self.masterFrame.master.withdraw()
    ##     else:
    ##         self.exitFunc()

    ## def editColor(self):
    ##     if self.edit == 1:
    ##         return
    ##     self.edit = 1
    ##     if self.hidden == 1:
    ##         self.ce.pack(side = 'right', fill='both', expand=1)
    ##         self.hidden = 0
    ##     else:
    ##         if self.add == 0:
    ##             self.ce.pack_forget()
    ##             self.hidden=1
    ##     if self.addHidden == 0:
    ##         self.addButton.pack_forget()
    ##         self.addHidden=1
    ##     self.ce.immediate=1
    ##     self.add = 0

    def addColor(self):
        if self.add ==1: return
        if self.hidden:
            self.ce.pack(side = 'right', fill='both', expand=1)
            self.hidden=0
        else:
            if self.edit == 0:
                self.ce.pack_forget()
                self.hidden=1
        if self.addHidden:
            self.addButton.pack(side = 'bottom', fill='x', expand=1)
            self.addHidden=0
        self.ce.immediate=0
        self.add = 1
        self.edit = 0

        
    #####################################################################
    #####     CALLBACKS
    ####################################################################
    def editChip(self, col):
        hexcol = ToHEX(col)
        chipName = self.colorChips.getcurselection()
        if chipName is None: return
        chip = self.colorChips.button(chipName)
        chip.configure(bg=hexcol, fg=hexcol, activebackground=hexcol,
                       activeforeground=hexcol, selectcolor=hexcol)
        self.colDict[chipName]=col
            
    def colButton_cb(self, tag):
        #this is needed to self.colorChips._buttonList working
        if tag in self.mapping.keys():
            tag = self.mapping[tag]
        
        col = self.colDict[tag]
        self.ce.set(col)
        if self.edit == 0 or self.add==0:
            self.cbManager.CallCallbacks(col)
            
        if self.doubleClick:
            if self.currentTag != tag:
                self.currentTag = tag 
                return
            self.doubleClick = False
            dialog = Pmw.PromptDialog(self.master, 
            title = 'Color Name Dialog',
            label_text = 'Enter text to label this color:',
            entryfield_labelpos = 'n',
            buttons = ('OK', 'Cancel'))
            result = dialog.activate()      
            if result == 'OK':
                txt = dialog.get() 
                fg='black' 
                if not txt: 
                    txt = str(len(self.colDict)+1)
                    fg = col
                self.colorChips.button(self.colorChips.index(
                                                             self.colorChips.selection
                                                             )).configure(text=txt,
                                                                          fg=fg,
                                                                          value=txt)
                color = self.colDict.pop(tag)
                self.colDict[txt] = color
                self.mapping[tag] = txt
                self.save(self.customFilename, self.colorsName)
                
        else:
            self.doubleClick = True
            self.master.after(1000, self.setDoubleClick)
    
    def setDoubleClick(self):
        self.doubleClick = False
        
    def addToCustom_cb(self):
        rgbcol = self.ce.get()
        newKey = str(len(self.colDict)+1)
        while newKey in self.colDict:
            newKey = str(int(newKey)+1)
        self.colDict[newKey]=rgbcol
        col = ToHEX(rgbcol)
        self.colorChips.add(newKey, bg = col,
                            activebackground=col,
                            fg = col, activeforeground = col,
                            indicatoron=0,selectcolor=col,
                            width=10,height=1, value = newKey)
        try:
            self.save(self.customFilename, self.colorsName)
        except Exception, inst:
            print inst

class BackgroundColorChooser(ColorChooser):
    def __init__(self, master=None, title = 'Chooser', commands = None,
                 immediate=0, exitFunction=None):
        ColorChooser.__init__(self, master, title, commands, immediate, exitFunction)
        tmpFrame = Tkinter.Frame(self.master)
        tmpFrame.pack(side='bottom',expand=0, fill='x')
        button = Tkinter.Button(tmpFrame, text='Make Default',  command=self.makeDefault)
        button.pack(side='left', expand=1, fill='x')
        button = Tkinter.Button(tmpFrame, text='Restore Default',  command=self.restoreDefault)
        button.pack(side='right', expand=1, fill='x')
        
        
    def makeDefault(self):
        colNum = self.colorVar.get()
        hcol = self.colors[colNum]
        path = os.path.join(self.rcPath,'backgroundColor')
        open(path,'w').write(hcol)
        
    def restoreDefault(self):
        path = os.path.join(self.rcPath, 'backgroundColor')
        if os.path.exists(path):
            os.remove(path)
        self.cbManager.CallCallbacks((.0,.0,.0))
            
            
class PaletteChooser():    
    def __init__(self, master=None, title = 'Palette Chooser', makeDefault_cb=None,
                 restoreDefault_cb=None, apply_cb=None, ramp=None, labels=None):
        self.ramp = ramp
        self.labels = labels
        self.makeDefault_cb = makeDefault_cb
        self.restoreDefault_cb = restoreDefault_cb
        self.apply_cb = apply_cb
        if master is None:
            self.master = Tkinter.Toplevel()
            self.ownmaster=1
            self.master.title(title)
            #self.master.protocol('WM_DELETE_WINDOW', self.dismiss)
        else:
            self.ownmaster=0
            self.master = master
        self.balloon = Pmw.Balloon(self.master)
        self.masterFrame = Tkinter.Frame(self.master,
                                         borderwidth=2, relief='ridge')
        self.masterFrame.pack()
        self.colorVar = Tkinter.IntVar()
        self.colorVar.set(-1)
        self.buttons = []
        for index, item in enumerate(self.labels):            
            col = ToHEX(self.ramp[index])
            b = Tkinter.Radiobutton(
                self.masterFrame, text=item, variable=self.colorVar, value=index,
                bg = col, activebackground=col,
                font=("Times",-16,'bold'), 
                indicatoron=0, selectcolor=col,
                width=5, height=1, command=self.selectColor)
            row = index%4
            b.grid(row=index/5, column=index%5)
            self.buttons.append(b)
                    
        if self.ownmaster:
            ## create dismiss button
            self.cb = self.master.destroy
            okCancelFrame = Tkinter.Frame(self.master,  borderwidth=2, relief='ridge')

            makeDefault = Tkinter.Button(okCancelFrame, text='Make Default',  command=self.makeDefault)
            makeDefault.grid(row=0, column=0)
            self.balloon.bind(makeDefault, 'Saves this palette as a default')
            restoreDefault = Tkinter.Button(okCancelFrame, text='Restore Default',   command=self.restoreDefault)
            restoreDefault.grid(row=0, column=1)
            self.balloon.bind(restoreDefault, 'Restores original palette as a default')
            Tkinter.Button(okCancelFrame, text='   Dismiss   ',  
                           command=self.cb).grid(row=1, column=0,columnspan=2)
            okCancelFrame.pack(fill='x')
        try:
            self.master.protocol('WM_DELETE_WINDOW', self.exit)
        except:
            pass 

    def selectColor(self, event=None):
        colNum = self.colorVar.get()
        self.chooser = ColorChooser(None, commands=self.color_cb, exitFunction=self.onColorChooserDismiss)
        self.chooser.pack(expand=1, fill='both')
    
    def restoreDefault(self):
        "Must be implemented for each instance"
        self.restoreDefault_cb()
    
    def onColorChooserDismiss(self):
        self.chooser.master.destroy()
        index = self.colorVar.get()
        self.buttons[index].deselect()
        
    def color_cb(self, color):
        col = ToHEX(color)
        index = self.colorVar.get()
        self.buttons[index].configure(bg = col, activebackground=col, selectcolor=col)
        self.ramp[index] = [color[0], color[1], color[2], 1] 
        self.apply_cb()
        
    def makeDefault(self):
        if self.makeDefault_cb:
            self.makeDefault_cb()
        self.cb()