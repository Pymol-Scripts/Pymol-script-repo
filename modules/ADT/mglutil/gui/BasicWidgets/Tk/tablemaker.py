## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#############################################################################
#
# Author: Anna Omelchenko
#
# Copyright: M. Sanner TSRI 2003
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/mglutil/gui/BasicWidgets/Tk/tablemaker.py,v 1.8 2008/07/14 21:58:25 vareille Exp $
#
# $Id: tablemaker.py,v 1.8 2008/07/14 21:58:25 vareille Exp $
#

"""tablemaker.py module provides widgets (transfer function editor) used
to edit color and opacity lookup tables.
Example:
    root = Tkinter.Tk()
    def AlphaCallback(value):
        # value is a 2 element list. value[0]-starting index of alpha
        # entries in the lookup table, value[1] - an integer array containing
        # new alpha values
        pass

    def RGBCallback(value):
        # value is a 2 element list. value[0]-starting index of color
        # entries in the lookup table, value[1] - a 2D array containing
        # new RGB values
        
    t = TableManager(master=root, xmaxval=255, xminval = 0, ymaxval=4095,
                     alphaCallback=AlphaCallback, colorCallback=ColorCallback)
"""

SelectPointTxt="Clicking the left mouse button on a dot selects it.\n"
SetOpacityTxt = "Dragging a dot with the left mouse button modifies\nthe  opacity.\n"
SetColorTxt = "Left clicking or moving the circle in the area labeled\n Hue/Saturation changes the color of the selected dot.\n The red-green-blue-opacity values can be directly typed\n in the text entries.\n"
AddDotShapeTxt = "Double left click on a line connecting two dots \nadds a dot to the shape.\nShift-left button click in the transfer function window\nadds a new set of four dots (a shape) to the function.\n"
MoveShapeTxt= "Middle clicking and draggind a dot of a shape moves the\nshape left or right in the window.\n"
DeleteDotTxt ="Dragging a dot up out of the window deletes it.\nDelete button of Edit menu can be used to delete a selected dot.\n"
SplitMenuTxt ="The 'Split function' button of Edit menu allows to\nzoom in on a particular interval of values.\nThe user must enter starting and ending value points of a\nnew interval in the pop up form. This will result in\nsplitting the initial interval of values into two or three\nintervals that will appear in separate transfer function boxes.\n"
ResetMenuTxt =  "The 'Reset' button of Edit menu sets the lookup table to\ngrayscale-ramp color values and linear opacity values.\n"
SaveTFTxt = "The 'Save TF' button of File menu allows to save current color\nand opacity values to a file.\n"
SaveLUTTxt = "The 'Save LUT (.lut)' button of File menu allows to save the\nlookup table to a .lut file, which can be used later to restore\ncurrent transfer function.\n"
LoadLUTTxt = "The 'Load LUT (.lut)' button loads a .lut file.\n"
    
AddISOTxt = "An ISO value can be added by pressing 'Add new isovalue'\nbutton of Edit menu.\nIn the pop up dialog the user can specify the scalar value,\nopacity and RGB values for the new ISO dot.\nThe RGBA values of selected ISOval can be typed in the\ncorresponding entry fields.\n"

DescriptionText = "Widget description:\nThe function is drawn at the top of the interface\nas a series of circles (dots) connected by line segments.\nThe X location of each dot indicates the data value,\nthe Y location indicates the opacity of that data value.\n" + SelectPointTxt + SetOpacityTxt + SetColorTxt + AddDotShapeTxt + MoveShapeTxt + DeleteDotTxt + LoadLUTTxt+ SaveLUTTxt+ SaveTFTxt + SplitMenuTxt+ ResetMenuTxt + AddISOTxt

import Tkinter
import numpy.oldnumeric as Numeric
import Pmw
import math

#from DejaVu.EventHandler import CallbackFunctions
from mglutil.gui.BasicWidgets.Tk.eventHandler import CallbackFunctions
#from DejaVu import colorTool
from mglutil.util import colorUtil
import os
import tkFileDialog
import cPickle
from mglutil.util.packageFilePath import findFilePath
from mglutil.util.colorUtil import ToRGB, ToHSV

#SELECTED_COLOR = "green"
SELECTED_COLOR = "gray68"
UNSELECTED_COLOR = '#CCCCCC'
_PI = math.pi


class TableManager:
    """Class for manipulating canvas widgets(lookup table), colormap widget
    and entry fields on the input form"""
    
    def __init__(self, viewer=None, master=None, xminval=0, xmaxval=255,
                 ymaxval=4095, alphaCallback = None, colorCallback = None,
                 iconpath=None):
        self.viewer = viewer
        self.xminval = int(xminval)
        self.xmaxval = int(xmaxval)
        self.ymaxval = int(ymaxval)
        self.ymaxval_limit = int(ymaxval)
        self.master = master
        self.lutcallbacks = {}

        balloon = Pmw.Balloon(self.master)
        #place menu buttons

        menuBar = Tkinter.Frame(self.master, relief = Tkinter.RAISED,
                                borderwidth=2)
        menuBar.pack(side = Tkinter.TOP, fill=Tkinter.X)
        fileBtn = Tkinter.Menubutton(menuBar, text="File", underline=0)
        fileBtn.pack(side=Tkinter.LEFT, padx ="2m")
        fileBtn.menu = Tkinter.Menu(fileBtn)
        fileBtn.menu.add_command(label="Load LUT (.lut)",
                                 command=self.ask_open_file_cb)
        fileBtn.menu.add_command(label="Save LUT (.lut)",
                                 command=self.saveLUT)
        fileBtn.menu.add_command(label='Save TF (.clu)', command=self.saveTF)
        fileBtn['menu'] = fileBtn.menu

        editBtn = Tkinter.Menubutton(menuBar, text="Edit", underline=0)
        editBtn.pack(side=Tkinter.LEFT, padx ="2m")
        editBtn.menu = Tkinter.Menu(editBtn)
        editBtn.menu.add_command(label="Reset", command=self.reset)
        editBtn.menu.add_command(label="Split function",
                                 command=self.show_splitDialog)
        editBtn.menu.add_command(label="Merge function",
                                 command=self.merge_function,
                                 state="disabled")

        editBtn.menu.drawchoices= Tkinter.Menu(editBtn.menu)
        editBtn.menu.add_cascade(label="Draw function...",
                                 menu=editBtn.menu.drawchoices)
        editBtn.menu.drawchoices.add_command(label ="ramp",
                      command=(lambda self=self, num=0: self.draw_ramp(num)))
        editBtn.menu.drawchoices.add_command(label ="reversed ramp",
                      command=(lambda self=self, num=1: self.draw_ramp(num)))
        editBtn.menu.add_command(label="Flip function",
                                 command = self.flip_function_cb)
        editBtn.menu.sizechoices= Tkinter.Menu(editBtn.menu)
        editBtn.menu.add_command(label="Set new opacity range",
                                 command=self.change_yaxis)
        editBtn.menu.add_command(label="Add new ISO value",
                                 command=self.addIsoVal_cb)
        editBtn.menu.add_cascade(label="Set size of ...",
                                 menu=editBtn.menu.sizechoices)
        editBtn.menu.sizechoices.add_command(label ="dots",
                                          command=(lambda self=self,
                                          st="dot": self.set_dotsize_cb(st)))
        editBtn.menu.sizechoices.add_command(label ="isovalue dots",
                                       command=(lambda self=self,
                                       st="isodot": self.set_dotsize_cb(st)))
        editBtn.menu.add_command(label="Delete selected dot",
                                 command=self.del_selected_dot)

        #editBtn.menu.options = Tkinter.Menu(editBtn.menu)
        #editBtn.menu.add_cascade(label="Options...",
        #                         menu=editBtn.menu.options)
        self.continVar = Tkinter.IntVar()
        self.continVar.set(1)
        editBtn.menu.add_checkbutton(label="Continuous update",
                                     var = self.continVar,
                                     command=self.set_continuous) 
        editBtn['menu'] = editBtn.menu
        self.editBtn = editBtn
        self.applylutBtn = Tkinter.Button(menuBar, text = "Apply LUT",
                                          command=self.applylut_cb,
                                          relief=Tkinter.FLAT, borderwidth=1,
                                          state="disabled")
        self.applylutBtn.pack(side=Tkinter.LEFT, padx='2m')
        balloon.bind(self.applylutBtn, "Applies LUT (non continuous mode)")
        
        helpBtn = Tkinter.Menubutton(menuBar, text="Help", underline=0)
        helpBtn.pack(side=Tkinter.RIGHT, padx ="2m")
        helpBtn.menu = Tkinter.Menu(helpBtn) 
        helpBtn.menu.add_command(label="Widget description",
                                 command = (lambda self=self,
                                            st="description": self.show_help_txt(st)))
        helpBtn.menu.choices = Tkinter.Menu(helpBtn.menu)
        helpBtn.menu.add_cascade(label="How to ...",
                                 menu=helpBtn.menu.choices)
        helpBtn.menu.choices.add_command(label="Add dot/shape",
                                        command=(lambda self=self,
                                                 st="addDotShape":self.show_help_txt(st)))
        helpBtn.menu.choices.add_command(label="Delete dot",
                                         command=(lambda self=self,
                                                  st="deleteDot": self.show_help_txt(st)))
        helpBtn.menu.choices.add_command(label="Move shape",
                                         command=(lambda self=self,
                                      st="moveShape": self.show_help_txt(st)))
        helpBtn.menu.choices.add_command(label="Set color",
                                         command=(lambda self=self,
                                      st="setColor": self.show_help_txt(st)))
        
        
        helpBtn['menu'] = helpBtn.menu
        
        font = self.font = '-*-Helvetica-Bold-R-Normal-*-*-160-*-*-*-*-*-*'
        self.set_font(menuBar, font)

        self.intervals_list = [(self.xminval, self.xmaxval),]
        self.parent_interval = (self.xminval, self.xmaxval)
        
        self.create_splitDialog()

        #place Lookup Table editor on the form
        self.f1 = Tkinter.Frame(self.master)
        self.f1.pack(side=Tkinter.TOP)
        lut = LUT(self.f1, xmaxval=self.xmaxval, xminval=xminval, grid_row=1,
                  num_of_entries = self.xmaxval+1, ymaxval=ymaxval)
        lut.canvas.bind('<Button-1>', self.selected)
        self.canvas_list = [lut.canvas,]
        lut.canvas.itemconfigure('outline', outline='blue')                    
        self.with_focus = lut
        self.lut_list = [lut,]

        #place entries for ISO value info
        
        isoFrame = Tkinter.Frame(self.master)
        isoFrame.pack(side=Tkinter.TOP)
        self.place_ISOentries(isoFrame)
        self.isoVal = Tkinter.StringVar()
        self.isoA = Tkinter.StringVar()
        self.isoR = Tkinter.StringVar()
        self.isoG = Tkinter.StringVar()
        self.isoB = Tkinter.StringVar()
        self.isoValDialog = None
        
        #place color editor on the form
        f = Tkinter.Frame(self.master)
        f.pack(side=Tkinter.TOP)
        colormapFrame = Tkinter.Frame(f)        
        colormapFrame.pack(side=Tkinter.LEFT,anchor=Tkinter.NW)#,padx=20)
        if not iconpath:
            iconpath = 'mglutil.gui.BasicWidgets.Tk'
            
        file = findFilePath('icons', iconpath)
        file = os.path.join(file,'colors.gif')
        self.colormap = Colormap(colormapFrame, file=file)
        self.colormap.callbacks = [self.set_hsv]
                           
        #place a group of entry fields on the form that display
        #alpha and scalar values
                           
        entrFrame = Tkinter.Frame(f, takefocus=0)        
        entrFrame.pack(side=Tkinter.LEFT, anchor=Tkinter.NW)#,padx=10)
        self.place_TFentries(entrFrame)
        self.lutcallbacks['entries'] = self.entries_update
        self.lutcallbacks['rgbmap'] = self.colormap.show_color
        if alphaCallback:
            self.lutcallbacks['setalpha'] = alphaCallback
        else:
            self.lutcallbacks['setalpha'] = self.lut_alpha_cb
        if colorCallback :
             self.lutcallbacks['setcolor'] = colorCallback
        else :
            self.lutcallbacks['setcolor'] = self.lut_color_cb
        lut.setCallbacks(self.lutcallbacks)
        self.data = [] #data saved to Undo splitting intervals

    def set_font(self, wid, newfont):
        try:
            wid.config(font=newfont)
        except :
            pass
        if len(wid.children)==0: return
        for item in wid.children.values():
            self.set_font(item, newfont)

    def create_splitDialog(self):
        """Creates 'Split interval' dialog used for obtaining
        interval split point values from the user."""
        
        self.splitDialog = Pmw.Dialog(self.master,
	    buttons=('OK','Cancel'), defaultbutton='OK',
	    title='Split interval dialog', command=self.execute)
	self.splitDialog.withdraw()

        Tkinter.Label(self.splitDialog.interior(),
              text = 'Enter min and max values of new interval\n'
                     '(length of the resulting intervals \n'
                     'must be greater or equal to 50)',
              font=('verdana',12),
              anchor=Tkinter.W).pack(expand = 1,
                                     fill = 'both', padx = 4, pady = 4)
        
        self.entry_maxval = int(self.xmaxval)
        self.entry_minval = int(self.xminval)
        self.min_val = Pmw.EntryField(self.splitDialog.interior(),
		labelpos = 'w',	value = '%d'%(self.xminval),
		label_text = 'Minval-int:',
                validate = self.validate_minval,
                modifiedcommand = self.min_val_changed)
        self.max_val = Pmw.EntryField(self.splitDialog.interior(),
		labelpos = 'w',	value = '%d'%(self.xmaxval),
                label_text = 'Maxval-int:',
                validate = self.validate_maxval)
        self.min_val.pack(fill='x', expand=1, padx=10, pady=8)
        self.max_val.pack(fill='x', expand=1, padx=10, pady=8)
        Pmw.alignlabels((self.min_val, self.max_val))


    def validate_minval(self, text):
        """Validate minval entry of the 'Split interval' dialog."""
        
        if text in ('', '-', '+'):
            return Pmw.PARTIAL
        try:
            value = int(text)
            for i in self.intervals_list:
                if value >= i[0]:
                    self.parent_interval = i
            if value < 0:
                return Pmw.ERROR
            if value != self.parent_interval[0]:
                if value < self.parent_interval[0]+49:
                    return Pmw.PARTIAL
            if value > self.parent_interval[1]-49:
                return Pmw.PARTIAL
            self.entry_minval=value
            return Pmw.OK
        except ValueError:
            return Pmw.ERROR

    def min_val_changed(self):
        """Called when minval entry of the 'Split interval' dialog is
        changed. """
        
        if self.min_val.valid():
            self.max_val.setentry(str(self.parent_interval[1]))


    def validate_maxval(self, text):
        """Validate maxval entry of the 'Split interval' dialog."""
        
        for i in self.intervals_list:
                if self.entry_minval >= i[0]:
                    self.parent_interval = i
        if text in ('', '-', '+'):
            return Pmw.PARTIAL
        try:
            value = int(text)
            if value > self.parent_interval[1]:
                return Pmw.ERROR
            if value < self.entry_minval+49:
                return Pmw.PARTIAL
            if value != self.parent_interval[1]:
                if value > self.parent_interval[1]-49:
                    return Pmw.ERROR
            self.entry_maxval=value
            return Pmw.OK
        except ValueError:
            return Pmw.ERROR

    def execute(self, result):
        """Called when 'OK' or 'Cancel' button of the 'Split dialog'
        is pressed. """
        
        if result == 'OK':
            entry_minval = int(self.min_val.get())
            entry_maxval = int(self.max_val.get())
            d1 = entry_minval-self.parent_interval[0]
            d2 = entry_maxval - entry_minval
            d3 = self.parent_interval[1] - entry_maxval
            if entry_minval == self.parent_interval[0] and \
               entry_maxval == self.parent_interval[1]:
                print "interval (%d,%d) exists" % (entry_minval, entry_maxval)
                self.splitDialog.deactivate(result)
                return
            elif d2< 49:
                print "Error: in 'Split dialog' (maxval-minval) must be greater or equal to 50 ; %d, %d" % (entry_maxval, entry_minval)
                return
            elif d1 != 0 and d1< 49:
                print "Error: in 'Split dialog' length of the resulting intervals must be greater or equal to 50(%d-%d < 50)"% (entry_minval, self.parent_interval[0])
                return
            elif d3 != 0 and d3 < 49:
                print "error: in 'Split dialog' length of the resulting intervals must be greater or equal to 50(%d-%d < 50)"% (self.parent_interval[1] ,entry_maxval)
                return
            else:
                #print "splitting"
                self.execute_split(entry_minval, entry_maxval)
                #self.split_cb(entry_minval, entry_maxval, split=1)
        self.splitDialog.deactivate(result) 

    def execute_split(self, entry_minval, entry_maxval,split=1):
        
        self.split_cb(entry_minval, entry_maxval, split=split)
    
    def show_splitDialog(self):
        """Activates the 'Split interval' dialog"""

        self.min_val.setentry(self.intervals_list[0][0])
        self.max_val.setentry(self.intervals_list[0][1])
        self.splitDialog.activate()
        
    def lut_color_cb(self, value):
        """Modifies the color table entries for given interval of data
        values. value is a list: value[0] - starting index,
        value[1] - an  array of RGB values.
        Each entry is in range 0.0 ... 1.0"""        
        pass

    def lut_alpha_cb(self, value):
        """Modifies the alpha(opacity) table entries for given interval of data
        values. value is a list: value[0] - starting index,
        value[1] - an integer array containing new alpha values.
        Each entry is in range 0...4095."""
        
        pass

    def change_yaxis(self):
        dialog = Pmw.CounterDialog(self.master,
                                   label_text = "Enter new Y axis max value",
                                   counter_labelpos='n',
                                   counter_datatype='numeric',
                                   entryfield_value = self.ymaxval,
                                 entryfield_validate={'validator':'numeric',
                                         'min':0, 'max':self.ymaxval_limit},
                                   buttons = ("OK", "Cancel"),
                                   defaultbutton = "OK",
                                   title="New Opacity Range(y_axis)")
        dialog. tkraise()
        result = dialog.activate()
        if result == "OK":
            new_y = int(dialog.get())
            #print  "new_y:", new_y
            if new_y != self.ymaxval:
                for lut in self.lut_list:
                    lut.update_yaxis(new_y)
                self.entry_y.configure(validate=
                                     {'validator':'integer',
                                      'min':0,'max': new_y})
                self.ymaxval = new_y

    def place_ISOentries(self, f):
        """Place text entries for opacity, color and scalar values
        of ISO points on the input form."""
        
        font = '-*-Helvetica-Bold-R-Normal-*-*-160-*-*-*-*-*-*'
        ew = 4
        xminval=self.with_focus.xminval
        xmaxval=self.with_focus.xmaxval
        self.entry_isoVal = Pmw.EntryField(f, labelpos = 'w',
                                           entry_font = font,
                                           entry_width = ew+1,
                                           label_text = 'ISO value: (None)',
                                           label_font = font,
                                           value = str(xminval),
                                           validate={'validator':'integer',
                                              'min': xminval,
                                              'max': xmaxval},
                                           command=self.iso_value_set)
        self.entry_isoY = Pmw.EntryField(f, labelpos = 'w', label_text = 'A',
                                         entry_width = ew+1, value = '0',
                                         entry_font = font, label_font = font,
                                         validate={'validator':'integer',
                                                   'min':0, 'max':self.ymaxval},
                                         command=self.iso_alpha_set)
        self.entry_isoR = Pmw.EntryField(f, labelpos = 'w', label_text = 'R',
                                         entry_width = ew, value = '1.0',
                                         entry_font = font, label_font = font,
                                         validate={'validator':'real',
                                                   'min':0, 'max':1.0},
                                         command=self.iso_color_set)
        self.entry_isoG = Pmw.EntryField(f, labelpos = 'w', label_text = 'G',
                                         entry_width = ew, value = '1.0',
                                         entry_font = font, label_font = font,
                                         validate={'validator':'real',
                                                   'min':0, 'max':1.0},
                                         command=self.iso_color_set)
        self.entry_isoB = Pmw.EntryField(f, labelpos = 'w', label_text = 'B',
                                         entry_width = ew, value = '1.0',
                                         entry_font = font, label_font = font,
                                         validate={'validator':'real',
                                                   'min':0, 'max':1.0},
                                         command=self.iso_color_set)
        
        entries = [self.entry_isoVal, self.entry_isoY,
                   self.entry_isoR, self.entry_isoG, self.entry_isoB]
        i = 0
        for entry in entries:
            entry.grid(row=0, column=i,sticky='w')
            i = i+1
##          Tkinter.Button(f, text='Add new', command=self.addIsoVal_cb,
##                         takefocus=0, font=font).grid(row=0, column=i,
##                                                      sticky='w') 

    def iso_color_set(self):
        """Callback for ISO color text entries. Called when
        the entry text is chaneged. """
        
        rgb = (float(self.entry_isoR.get()),
               float(self.entry_isoG.get()),
               float(self.entry_isoB.get()))
        self.with_focus.set_ISOcolor(rgb)
        hsv = ToHSV(rgb)
        self.colormap.show_color(hsv)

    def iso_alpha_set(self):
        """Callback for ISO alpha(opacity) text entry. Called when
        the entry text is chaneged. """
        
        alpha_val = int(self.entry_isoY.get())
        self.with_focus.moveIsoVal_alpha(alpha_val)

    def iso_value_set(self):
        """Callback for ISO scalar value text entry. Called when
        the entry text is chaneged. """
        
        val = int(self.entry_isoVal.get())
        self.with_focus.moveIsoVal_value(val)

    def addIsoVal_cb(self):
        """Callback function of 'Add new' button.
        Creates a dialog for entering RGBA values of a new
        ISO point. """
        
        font = '-*-Helvetica-Bold-R-Normal-*-*-160-*-*-*-*-*-*'
        dialog = Pmw.Dialog(self.master,
                            buttons=('OK','Cancel'), defaultbutton='OK',
                            title='New ISOval dialog', 
                            command=self.addIsoVal)
	#dialog.withdraw()
        ew = 7
        root =dialog.interior()
        Pmw.EntryField(root, labelpos = 'w', 
                       entry_width = ew, entry_textvariable = self.isoVal,
                       label_text = 'ISO val:', 
                       #label_font = font, entry_font = font,
                       validate={'validator':'integer',
                                 'min':0, 'max':self.xmaxval},
                       value = 0).pack(side=Tkinter.TOP)
        Pmw.EntryField(root, labelpos = 'w', 
                       entry_width = ew, value = '0',
                       #entry_font = font, label_font = font,
                       entry_textvariable = self.isoA,
                       label_text = 'Alpha:',
                       validate={'validator':'integer',
                         'min':0, 'max':self.ymaxval}).pack(side=Tkinter.TOP)
        Pmw.EntryField(root, labelpos = 'w', 
                       entry_width = ew, value = '1.0',
                       #entry_font = font, label_font = font,
                       entry_textvariable = self.isoR,  label_text = 'Red:',
                       validate={'validator':'real',
                                 'min':0, 'max':1.0}).pack(side=Tkinter.TOP)
        Pmw.EntryField(root, labelpos = 'w', 
                       entry_width = ew, value = '1.0',
                       #entry_font = font, label_font = font,
                       entry_textvariable = self.isoG,  label_text = 'Green:',
                       validate={'validator':'real',
                                 'min':0, 'max':1.0}).pack(side=Tkinter.TOP)
        Pmw.EntryField(root, labelpos = 'w',
                       entry_width = ew, value = '1.0',
                       #entry_font = font, label_font = font,
                       entry_textvariable = self.isoB,   label_text = 'Blue:',
                       validate={'validator':'real',
                                 'min':0, 'max':1.0}).pack(side=Tkinter.TOP)
        self.isoValDialog = dialog
        dialog.activate()

    def addIsoVal(self, val):
        """Adds an ISO value to the transfer function editor
        (based on the users input in isoValDialog)"""
        
        if val == "OK":
            isoVal = int(self.isoVal.get())
            isoA = int(self.isoA.get())
            isoRGB = ( float(self.isoR.get()),
                       float(self.isoG.get()),
                       float(self.isoB.get()) )
            i = 0
            ind = 0
            num_isovals = 0
            for lut in self.lut_list:
                if isoVal in range(lut.xminval, lut.xmaxval+1):
                    #self.with_focus = lut
                    ind = i 
                i = i + 1
                num_isovals = num_isovals + len(lut.isoVals)
            lut = self.lut_list[ind]
            self.make_selected(lut.canvas)
            lut.drawIsoVal(isoVal, isoA, isoRGB)
            if num_isovals == 0:
                self.entry_isoVal.configure(label_text = 'ISO value:')

        self.isoValDialog.deactivate()
        
    
    def rmIsoVal_cb(self):
        self.with_focus.removeIsoVal()
        num_isovals = 0
        for lut in self.lut_list:
            num_isovals = num_isovals + len(lut.isoVals)
        if num_isovals == 0:
                self.entry_isoVal.configure(label_text = 'ISO value (None):')
                
    def place_TFentries(self, f):
        """Place RGBA textentries of the transfer function points
        on the input form."""
        
        font = '-*-Helvetica-Bold-R-Normal-*-*-160-*-*-*-*-*-*'
        ew = 10
        self.entry_y = Pmw.EntryField(f, labelpos = 'w', label_text = 'Alpha',
                                      entry_width = ew, value = '0',
                                      entry_font = font, label_font = font,
                                      validate={'validator':'integer',
                                                'min':0, 'max':self.ymaxval},
                                      command=self.entry_alpha_set)
        self.entry_R = Pmw.EntryField(f, labelpos = 'w', label_text = 'Red',
                                      entry_width = ew, value = '1.0',
                                      entry_font = font, label_font = font,
                                      validate={'validator':'real',
                                                'min':0, 'max':1.0},
                                      command=self.entry_color_set)
        self.entry_G = Pmw.EntryField(f, labelpos = 'w', label_text = 'Green',
                                      entry_width = ew, value = '1.0',
                                      entry_font = font, label_font = font,
                                      validate={'validator':'real',
                                                'min':0, 'max':1.0},
                                      command=self.entry_color_set)
        self.entry_B = Pmw.EntryField(f, labelpos = 'w', label_text = 'Blue',
                                      entry_width = ew, value = '1.0',
                                      entry_font = font, label_font = font,
                                      validate={'validator':'real',
                                                'min':0, 'max':1.0},
                                      command=self.entry_color_set)
        self.entry_x1 = Pmw.EntryField(f, labelpos = 'w',
                                       entry_width = ew, entry_font = font,
                                       label_text = 'Scalar',
                                       label_font = font, value = '0',
                                       validate={'validator':'integer',
                                                'min':0, 'max':self.xmaxval},
                                       command = self.entry_value_set)
        self.entry_x2=Pmw.EntryField(f, labelpos = 'w',
                                     entry_width = ew, entry_font = font,
                                     label_text = 'Scalar',
                                     label_font = font, value = '0')
        entries = [self.entry_y, self.entry_R, self.entry_G, self.entry_B,
                   self.entry_x1, self.entry_x2]
        i = 1                              
        for entry in entries:
            entry.grid(column=0, row = i, sticky = 'NSEW', pady=4, padx=20)
            i = i+1
        Pmw.alignlabels(entries)

    def entry_alpha_set(self):
        """Callback function of 'Alpha' text entry. Sets typed in alpha
        value and moves the dot."""

        alpha = int(self.entry_y.get())
        self.with_focus.move_dot_alpha(alpha)

    def entry_value_set(self):
        """Callback function of 'Scalar' text entry. Sets typed in scalar
        value and moves the dot."""
        try:
            val = int(self.entry_x1.get())
        except ValueError:
            self.entry_x1.cget("invalidcommand")() 
            self.entry_x1.clear()
            return
        res = self.with_focus.check_scalarval(val)
        if not res:
            self.entry_x1.cget("invalidcommand")() 
            self.entry_x1.clear()
        else:
            self.with_focus.move_dot_value(val)

    def entry_color_set(self):
        """Callback function of 'Red', 'Green' and 'Blue' text entries.
        Sets typed in color value to selected scalar value."""
        
        rgb = (float(self.entry_R.get()),
               float(self.entry_G.get()),
               float(self.entry_B.get()))
        hsv = ToHSV(rgb)
        self.with_focus.set_red_green_blue(rgb)
        self.colormap.show_color(hsv)

    def selected(self, event):
        """Activates selected canvas widget."""
        curr_widget = event.widget
        self.make_selected(curr_widget)

    def make_selected(self, curr_widget):
        if curr_widget == self.with_focus.canvas:
            return
        curr_widget.itemconfigure('outline', outline='blue')
        self.with_focus.canvas.itemconfigure('outline', outline='')
        ind = self.canvas_list.index(curr_widget)
        i = 0
        for lut in self.lut_list:
            if i == ind:
                lut.bind_tags()
            else:
                lut.unbind_tags()
            i=i+1
        self.with_focus = self.lut_list[ind]


    def split_cb(self, entry_minval, entry_maxval, split=1):
        """Executing command of the 'Split Interval Dialog'.
        Determines new intervals of values. Destroys canvas widget
        representing the original interval. Calls a function creating
        canvas widgets for the new intervals."""
        
        parent_interval = self.parent_interval
        try:
            ind = self.intervals_list.index(parent_interval)
        except ValueError:
            print 'ValueError: interval',parent_interval,'is not in self.intervals_list'
            return
        if entry_minval == parent_interval[0]:
            intervals = [(entry_minval, entry_maxval),
                         (entry_maxval+1, parent_interval[1])]
            curr_interval = ind
        elif entry_maxval == parent_interval[1]:
            intervals = [(parent_interval[0], entry_minval-1),
                         (entry_minval, entry_maxval)]
            curr_interval = ind+1
        else:
            intervals = [(parent_interval[0], entry_minval-1),
                         (entry_minval, entry_maxval),
                         (entry_maxval+1, parent_interval[1])]
            curr_interval = ind+1 

        self.intervals_list.pop(ind)
        self.with_focus.canvas.itemconfigure('outline', outline='')
        self.canvas_list[ind].destroy()
        self.canvas_list.pop(ind)
        old_lut = self.lut_list.pop(ind)
        #print "ind :", ind
        i = ind
        for interval in intervals:
            self.intervals_list.insert(i, interval)
            lut = LUT(self.f1, xminval=interval[0],
                      xmaxval=interval[1], ymaxval=self.ymaxval,
                      grid_row=i+1,
                      num_of_entries=self.xmaxval+1, initdraw = 0)
            lut.canvas.bind('<Button-1>', self.selected)
            self.lut_list.insert(i, lut)
            self.canvas_list.insert(i, lut.canvas)
            i = i + 1
            lut.setCallbacks(self.lutcallbacks)
        no_intervals = len(self.intervals_list)
        if i < no_intervals:
            for n in range(i, no_intervals):
                self.canvas_list[n].grid(row=n, sticky=W)
        self.canvas_list[curr_interval].itemconfigure('outline',
                                                     outline='blue')
        self.with_focus = self.lut_list[curr_interval]
        i = 0
        for interval in intervals:
            self.split_function(old_lut, interval, ind+i)
            i = i + 1
        old_lut = None
        self.editBtn.menu.entryconfigure("Merge function", state='normal')
        

    def entries_update(self,**values):
        """Updates RGBA entries on the input form."""
        
        for k in values.keys():
            if k == 'val_x1':
                self.entry_x1.setentry(str(values[k]))
            elif k == 'val_x2':
                self.entry_x2.setentry(str(values[k]))
            elif k == 'val_y':
                self.entry_y.setentry(str(values[k]))
            elif k == 'val_r':
                self.entry_R.setentry(str(values[k]))
            elif k == 'val_g':
                self.entry_G.setentry(str(values[k]))
            elif k == 'val_b':
                self.entry_B.setentry(str(values[k]))
            elif k == 'iso_val':
                self.entry_isoVal.setentry(str(values[k]))
            elif k == 'iso_y':
                self.entry_isoY.setentry(str(values[k]))
            elif k == 'iso_rgb':
                self.entry_isoR.setentry(str(values[k][0]))
                self.entry_isoG.setentry(str(values[k][1]))
                self.entry_isoB.setentry(str(values[k][2]))
    
    def set_hsv(self, hsv):
        """Callback function of the Colormap widget.
        Sets hue, saturation and value to selected scalar value."""
        
        self.with_focus.set_hue_satur_value(hsv)


    def split_function(self, old_lut, interval, new_ind):
        """Find values, points, alphas of a resulting interval
        after interval splitting."""
        
        lut = self.lut_list[new_ind]
        max_val = interval[1]
        min_val = interval[0]
##          print 'split_function: min_val =',min_val,'max_val =', max_val 
	#indexes of min_val and max_val in alpha and color arrays
	min_arrind = min_val - old_lut.xminval
	max_arrind = max_val - old_lut.xminval
        new_values = []
        for val in old_lut.values:
            if val >  min_val and val <  max_val:
                new_values.append(val)
            elif val >= max_val:
                next_max_val = val
                break
##          print "new_values1 =", new_values
        sclx = lut.sclx
        left = lut.left
        right = lut.right
        bott = lut.bott            
        if len(new_values) != 0:
            min_ind = old_lut.values.index(new_values[0])
            max_ind = old_lut.values.index(new_values[-1])
        
            new_shapes = []
            for s in old_lut.shapes:
                if s >= min_ind and s <= max_ind:
                    new_shapes.append(s)
##          print "new_shapes1=", new_shapes
            new_points = []
            i = min_ind
            for val in new_values:
                new_points.append(((val-min_val)*sclx+left,
                                  old_lut.points[i][1]))
                i = i + 1
            if not(min_ind in old_lut.shapes and min_ind-1 in old_lut.shapes):
                Y = bott-old_lut.scly*old_lut.alpha_arr[min_arrind]
                new_points.insert(0, (left,Y))
                new_values.insert(0, min_val)
                min_ind = min_ind - 1
                new_shapes.insert(0, min_ind)
            
            if not(max_ind in old_lut.shapes and max_ind+1 in old_lut.shapes):
                Y = bott-old_lut.scly*old_lut.alpha_arr[max_arrind]
                new_points.append((right, Y))
                new_values.append(max_val)
                max_ind = max_ind + 1
                new_shapes.append(max_ind)
            
            new_values.insert(0, min_val)
            new_values.append(max_val)
            new_points.insert(0, (left, bott))
            new_points.append((right, bott))
            min_ind = min_ind - 1
            max_ind = max_ind + 1
            new_shapes.insert(0, min_ind)
            new_shapes.append(max_ind)
##          print "new_values2 =", new_values
##          print "min_ind=", min_ind
##          print "new_shapes2 =", new_shapes
            new_shapes = map(lambda x, y=min_ind: x-y, new_shapes)
##          print "new_shapes3 =", new_shapes
            
        else:
          max_ind = old_lut.values.index(next_max_val)
          min_ind = max_ind-1
          if min_ind in old_lut.shapes and max_ind in old_lut.shapes:
              new_values = [min_val, min_val, max_val, max_val]
              new_points = [(left,bott), (left,bott), (right,bott), (right,bott)]
              new_shapes = [0,1,2,3]
##                return
          else:
              new_values = [min_val, min_val, max_val, max_val]
              Y1 = bott-old_lut.scly*old_lut.alpha_arr[min_arrind]
              Y2 = bott-old_lut.scly*old_lut.alpha_arr[max_arrind]
              new_points = [(left,bott), (left,Y1), (right,Y2), (right,bott)]
              new_shapes = [0,1,2,3]
        
        v_len = len(new_values)
        if v_len == 4:
            if new_values[2] - new_values[1] <= 1:
                new_values = [min_val, min_val, max_val, max_val]
                new_points = [(left,bott), (left,bott), (right,bott), (right,bott)]
                new_shapes = [0,1,2,3]
##                  return
            else:
                new_points.insert(2,
                       (new_points[1][0]+(new_points[2][0]-new_points[1][0])/2,
                       new_points[1][1]+(new_points[2][1]-new_points[1][1])/2))
                new_values.insert(2,
                       int(round((new_points[2][0]-left)/sclx))+min_val )
                new_shapes = [0,1,3,4]
        else:
            if new_shapes[-2]-new_shapes[-3] == 1:
                if new_values[-2]-new_values[-3] <= 1:
                    new_points.pop(-2)
                    new_points.pop(-2)
                    new_values.pop(-2)
                    new_values.pop(-2)
                    new_shapes.pop(-2)
                    new_shapes.pop(-2)
                else:
                    new_points.insert(v_len-2,
                                      (new_points[-3][0]+
                                       (new_points[-2][0]-new_points[-3][0])/2,
                                       new_points[-3][1]+
                                      (new_points[-2][1]-new_points[-3][1])/2))
                    new_values.insert(v_len-2,
                             int(round((new_points[-3][0]-left)/sclx))+min_val)
                    new_shapes[-2] = new_shapes[-2]+1
                    new_shapes[-1] = new_shapes[-1]+1
            if new_shapes[2]-new_shapes[1] == 1:
                if new_values[2]-new_values[1] <= 1:
                    new_points.pop(1)
                    new_points.pop(1)
                    new_values.pop(1)
                    new_values.pop(1)
                    new_shapes.pop(1)
                    new_shapes.pop(1)
                    new_shapes = map(lambda x: x-2, new_shapes)
                else:
                    new_points.insert(2,
                       (new_points[1][0]+(new_points[2][0]-new_points[1][0])/2,
                       new_points[1][1]+(new_points[2][1]-new_points[1][1])/2))
                    new_values.insert(2,
                             int(round((new_points[2][0]-left)/sclx))+min_val )
                    for i in range(2, len(new_shapes)):
                        new_shapes[i] = new_shapes[i]+1
        lut.shapes = new_shapes
        lut.values = new_values
        lut.points = new_points
        lut.alpha_arr[:] = \
                         old_lut.alpha_arr[min_arrind:max_arrind+1]
        lut.color_arr[:] = \
                         old_lut.color_arr[min_arrind:max_arrind+1]
        lut.alpha_arr[:] = \
                         old_lut.alpha_arr[min_arrind:max_arrind+1]
##          print 'new_interval:', interval
##          print 'old_shapes:', old_lut.shapes
##          print 'old_points:', old_lut.points
##          print 'old_values:', old_lut.values
##          print 'old_alphas:', old_lut.alpha_arr
##          print '**********'
##          print 'new_shapes:', lut.shapes
##          print 'new_points:', lut.points
##          print 'new_values:', lut.values
##          print 'new_alphas:', lut.alpha_arr
        lut.isoVals= []
        if len(old_lut.isoVals):
            for iso_val in old_lut.isoVals:
                v = iso_val['val']
                if v >= lut.xminval and v <= lut.xmaxval:
                    lut.isoVals.append(iso_val)
        lut.redraw()

        
    def ask_save_file(self, filetypes, title):
        """pop up tkFileDialog"""
        file = tkFileDialog.asksaveasfilename(filetypes=filetypes,
                                                      title=title)
        return file

    def saveLUT_old(self, file = None):
        """saves lookup table in .lut file that can be used to restore
        current transfer function """

        if not file:
            file = self.ask_save_file([("LUT Files","*.lut")], 'Save LUT')
        if file:
            of = open(file, 'w')
            cPickle.dump(self.intervals_list, of)
            for lut in self.lut_list:
                data = [lut.points, lut.shapes, lut.color_arr]
                cPickle.dump(data, of)
            of.close()

    def saveLUT(self, file=None):
        """saves lookup table in .lut file that can be used to restore
        current transfer function """

        if not file:
            file = self.ask_save_file([("LUT Files","*.lut")], 'Save LUT')
        of = open(file, 'w')
        of.write("Transfer function\n")
        of.write("Maxalpha %d\n" % self.ymaxval)
        intervals = self.intervals_list
        of.write("NumIntervals %d\n" % len(intervals) )
        intervals_str = "Intervals"
        for i in intervals:
            intervals_str = intervals_str+" %d %d"%(i[0], i[1])
        intervals_str = intervals_str+"\n"
        of.write(intervals_str)
        numbers = "DataSizes"
        for lut in self.lut_list:
            numbers = numbers+" %d %d %d"%( len(lut.shapes),
                                               len(lut.values),
                                               len(lut.color_arr)*3)
            #print "shapes:", lut.shapes
        numbers = numbers+"\n"
        #print "numbers: ", numbers
        of.write(numbers)
        of.write("End header\n")
        import struct
        for lut in self.lut_list:
            fmt_shapes = ">%di"%len(lut.shapes)
            fmt_colors = ">%df"% (len(lut.color_arr)*3 ,)
            fmt_values = ">%di"%len(lut.values)
            alpha_inds = map(lambda x: x-lut.xminval, lut.values)
            alphas = Numeric.take(lut.alpha_arr, alpha_inds)
            #print fmt_values, fmt_shapes, fmt_colors
            of.write(apply(struct.pack, (fmt_shapes,)+tuple(lut.shapes)))
            of.write(apply(struct.pack, (fmt_values,)+tuple(lut.values)))
            of.write(apply(struct.pack, (fmt_values,)+tuple(alphas)))
            of.write( apply(struct.pack, (fmt_colors,)+tuple(lut.color_arr.ravel())) )
        of.close()
            

    def load_file(self, file):
        from string import split, atoi
        of = open(file, 'r')
        line = split(of.readline())
        warning = "Warning: wrong file format. Could not read file %s" %(file,)
        if len(line):
            if line[0] != "Transfer":
                of.close()
                if self.load_file_old(file):
                    return
                else:
                    print warning
                return
        else:
            of.close()
            print warning
            return
        ymaxval = 0
        while(1):
            line = of.readline()
            line = split(line)
            if len(line):
                if line[0] == "End": break
                else:
                    if line[0] == "NumIntervals":
                        nintervals = atoi(line[1])
                    elif line[0] == "Intervals":
                        intervals_str= line[1:]
                    elif line[0] == "DataSizes":
                        sizes_str = line[1:]
                    elif line[0] == "Maxalpha":
                        ymaxval = atoi(line[1])
        
        intervals = []
        sizes = []
        for n in range(nintervals):
            intervals.append( (atoi(intervals_str[n*2]),
                                    atoi(intervals_str[n*2+1]) ))
            sizes.append( (atoi(sizes_str[n*3]),
                           atoi(sizes_str[n*3+1]),
                           atoi(sizes_str[n*3+2])) )
        #print "nintervals: ", nintervals
        #print "intervals: ", intervals
        #print "sizes: ", sizes
        
        data = []
        xmaxval = intervals[nintervals-1][1]
        if xmaxval != self.xmaxval:
            if self.viewer.hasGui:
                text = "WARNING: number of LUT entries in\n%s\n is %d,\n current number of LUT entries is %d.\nLoad new LUT?" %(file,xmaxval+1,self.xmaxval+1)
                dialog = Pmw.MessageDialog(self.master,
                                           buttons=('OK','Cancel'),
                                           defaultbutton='OK',
                                           title='Load New LUT',
                                           message_text=text)
                result=dialog.activate()
                if result=='Cancel':
                    of.close()
                    return
                else:
                    self.xmaxval = xmaxval
            else:
                print "WARNING: number of LUT entries in %s is %d, current number of LUT entries is %d." % (file, xmaxval+1, self.xmaxval+1)

        shapes = []
        values = []
        alphas = []
        colors = []
        from struct import unpack, calcsize
        for n in range(nintervals):
            fmt_shapes = ">%di"%sizes[n][0]
            fmt_values = ">%di"%sizes[n][1]
            fmt_colors = ">%df"%sizes[n][2]
            l_shapes = of.read(calcsize(fmt_shapes))
            #values and alphas have the same format
            l_values= of.read(calcsize(fmt_values))
            l_alphas = of.read(calcsize(fmt_values))
            l_colors = of.read(calcsize(fmt_colors))
            shapes.append(list(unpack(fmt_shapes, l_shapes)))
            values.append(unpack(fmt_values, l_values))
            alphas.append(unpack(fmt_values, l_alphas))
            colors.append(unpack(fmt_colors, l_colors))
        #print 'shapes:', shapes    
        #print 'values: ', values
        #print 'alphas: ', alphas
        
        of.close()
        d = 1
        if ymaxval:
##              d = (ymaxval*1.0)/self.ymaxval
            self.ymaxval = ymaxval
        #print "ymaxval: ", ymaxval, "self.ymaxval: ", self.ymaxval
        self.xmaxval = xmaxval
        for canvas in self.canvas_list:
            canvas.destroy()
        self.canvas_list = []
        self.lut_list = []
        self.intervals_list = intervals
        i = 0
        for interval in intervals:
##              print 'in load_file: interval =', interval
            lut = LUT(self.f1, xminval=interval[0],
                      xmaxval=interval[1], ymaxval = self.ymaxval,
                      grid_row=i,
                      num_of_entries=xmaxval+1, initdraw = 0)
            lut.canvas.bind('<Button-1>', self.selected)
            self.lut_list.append(lut)
            self.canvas_list.append(lut.canvas)
            if d != 1 :
                lut.calculate_points(values[i],
                                     map(lambda x: x/d, alphas[i]))
            else:
                lut.calculate_points(values[i], alphas[i])
            lut.shapes = shapes[i]
            colors_size = (len(colors[i])/3, 3)
            #print "colors_size: ", len(colors[i]), colors_size
            lut.color_arr = Numeric.reshape(Numeric.array(colors[i], 'f'),
                                            colors_size)
            lut.calculate_alphas()
            lut.redraw()
            lut.setCallbacks(self.lutcallbacks)
            lut.callbacks['setalpha']([interval[0], lut.alpha_arr])
            lut.callbacks['setcolor']([interval[0], lut.color_arr])
            i = i + 1
        self.lut_list[0].canvas.itemconfigure('outline', outline='blue')
        self.with_focus = self.lut_list[0]
        if len(intervals)>1:
            self.editBtn.menu.entryconfigure("Merge function", state='normal')
 
    def saveData(self):
        intervals = self.intervals_list[:]
        data ={'intervals':intervals}
        i = 0
        for interval in self.intervals_list:
            lut = self.lut_list[i]
            points = lut.points[:]
            shapes = lut.shapes[:]
            color_arr = (lut.color_arr*1).astype('f')
            data[interval] = [points, shapes, color_arr]
            i = i+1
        self.data.append(data)

    def saveTF(self):
        """saves transfer function (opacity and color values)in .clu file
        ('Clut1999a'  Mitsubishi Electric format).
        Range of color values: 0 ... 255, range of alpha values: 0.0 ... 0.1"""
        
        filename = self.ask_save_file([("TF Files","*.clu")],
                                  'Save Transfer Function')
        if filename:
            data = self.getLookupTable()
            if data:
                file = open(filename, 'w')
                ymaxval = float(self.ymaxval)
                if file:
                    file.write("Clut1999a\n")
                    file.write("%d\n"%len(data))
                    file.write("RGBRange %d\n"%(self.xmaxval,)) 
                    file.write("AlphaRange 1\n")
                    file.write("Title CLUT\n")
                    file.write("Copyright 1998 Mitsubishi Electric I.T.A.\n")
                    file.write("##\n\n")
                    for i in range(len(data)):
                        file.write("%d %d %d %d %.3f\n"%(i,data[i][0],
                                                      data[i][1],data[i][2],
                                                      data[i][3]/ymaxval) )
                    file.close()
        else:
            return

    def getLookupTable(self):
        """build an array of RGBA for the current lookup table.
        """
        rgba=Numeric.zeros((self.xmaxval+1, 4), 'f')
        alpha = Numeric.zeros(self.xmaxval+1,'i')
        for lut in self.lut_list:
            xminval = lut.xminval
            xmaxval = lut.xmaxval
            rgba[xminval:xmaxval+1, :3] = lut.color_arr[:]
            alpha[xminval:xmaxval+1] = lut.alpha_arr[:]
        rgba=(rgba*255).astype('i')
        rgba[:,-1]=alpha
        return rgba
        
    def ask_open_file_cb(self):
        file = tkFileDialog.askopenfilename\
               (filetypes=[("LUT Files", "*.lut")], title='Load LUT')
        if file:
            if os.path.splitext(file)[1] == '.lut':
                self.load_lutfile(file)
            else:
                print "Wrong file name", file
                return
            
    def load_lutfile(self, file):
        self.load_file(file)

    def load_file_old(self, file):
        of = open(file, 'r')
        try:
            intervals = cPickle.load(of)
        except:
            print "Load LUT ERROR: could not read file: ", file
            return
        data = []
        nintervals = len(intervals)
        xmaxval = intervals[nintervals-1][1]
        if xmaxval != self.xmaxval:
            if self.viewer:
                if self.viewer.hasGui:
                    text = "WARNING: number of LUT entries in\n%s\n is %d,\n current number of LUT entries is %d.\nLoad new LUT?" %(file,xmaxval+1,self.xmaxval+1)
                    dialog = Pmw.MessageDialog(self.master,
                                               buttons=('OK','Cancel'),
                                               defaultbutton='OK',
                                               title='Load New LUT',
                                               message_text=text)
                    result=dialog.activate()
                    if result=='Cancel':
                        of.close()
                        return
                    else:
                        self.xmaxval = xmaxval
            else:
                print "WARNING: number of LUT entries in %s is %d, current number of LUT entries is %d." % (file, xmaxval+1, self.xmaxval+1)
        for n in range(nintervals):
            try:
                data.append(cPickle.load(of))
            except:
                print "Load LUT Error: could not read the data."
                of.close()
                return 0
##              print 'lut No %s'%n
##              print 'points:', data[n][0]
##              print 'shapes:', data[n][1]
        of.close()
        self.xmaxval = xmaxval
        for canvas in self.canvas_list:
            canvas.destroy()
        self.canvas_list = []
        self.lut_list = []
        self.intervals_list = intervals
        i = 0
        for interval in intervals:
##              print 'in load_file: interval =', interval
            lut = LUT(self.f1, xminval=interval[0],
                      xmaxval=interval[1], ymaxval = self.ymaxval,
                      grid_row=i,
                      num_of_entries=xmaxval+1, initdraw = 0)
            lut.canvas.bind('<Button-1>', self.selected)
            self.lut_list.append(lut)
            self.canvas_list.append(lut.canvas)
            lut.points = data[i][0]
            lut.shapes = data[i][1]
            lut.color_arr = data[i][2]
            lut.calculate_alphas()
            lut.redraw()
            lut.setCallbacks(self.lutcallbacks)
            lut.callbacks['setalpha']([interval[0], lut.alpha_arr])
            lut.callbacks['setcolor']([interval[0], lut.color_arr])
            i = i + 1
        self.lut_list[0].canvas.itemconfigure('outline', outline='blue')
        self.with_focus = self.lut_list[0]
        return 1

    def undo_split(self, data = None):
        calculatealpha = 0
        if not data:
            data = self.data
            calculatealpha = 1
        if len(data) == 0:
            return
        for canvas in self.canvas_list:
            canvas.destroy()
        self.canvas_list = []
        self.lut_list = []
        intervals = data[-1]['intervals']
        #print "in undo_split, intervals_list:", self.intervals_list
        #print "data['intervals']=", self.data[-1]['intervals']
        self.intervals_list = intervals
        i=0
        for interval in intervals:
##              print 'in load_file: interval =', interval
            lut = LUT(self.f1, xminval=interval[0], ymaxval=self.ymaxval,
                      xmaxval=interval[1], grid_row=i,
                      num_of_entries=self.xmaxval+1, initdraw = 0)
            lut.canvas.bind('<Button-1>', self.selected)
            self.lut_list.append(lut)
            self.canvas_list.append(lut.canvas)
            lut.points = data[-1][interval][0]
            lut.shapes = data[-1][interval][1]
            lut.color_arr = data[-1][interval][2]
            if calculatealpha:
                lut.calculate_alphas()
            else:
                lut.alpha_arr = data[-1][interval][3]
            lut.redraw()
            lut.setCallbacks(self.lutcallbacks)
            if calculatealpha:
                lut.callbacks['setalpha']([interval[0], lut.alpha_arr])
                lut.callbacks['setcolor']([interval[0], lut.color_arr])
            i=i+1
        self.lut_list[0].canvas.itemconfigure('outline', outline='blue')
        self.with_focus = self.lut_list[0]
        data.pop()
        
    def reset(self):
        for lut in self.lut_list:
            lut.canvas.delete('line')
            lut.canvas.delete('dot')
            lut.canvas.delete('isoline')
            lut.canvas.delete('isodot')
            lut.isoVals = []
            lut.draw_initshape()
            lut.callbacks['setalpha']([lut.xminval, lut.alpha_arr])
            lut.callbacks['setcolor']([lut.xminval, lut.color_arr])

    def getTransferFunc(self):
        """ returns a list containing LUT values(scalar values of points on
        the widget transfer function), corresponding alphas and colors"""
        values = []
        alphas = []
        colors = []
        for lut in self.lut_list:
            values.extend(lut.values[1:-1])
            xminval = lut.xminval
            for v in lut.values[1:-1]:
                alphas.append(lut.alpha_arr[v-xminval])
                colors.append(lut.color_arr[v-xminval].tolist())
        return values, alphas, colors

    def getIsoVals(self):
       isoVals = []
       for lut in self.lut_list:
           if len(lut.isoVals):
               isoVals.extend(lut.isoVals)
       return isoVals

    def draw_ramp(self, reverse = 0):
        """Draw a ramp that goes from left lower corner to right upper corner
        (or if 'reverse' is specified - from left upper corner to right lower
        corner) of currenly selected LUT frame."""

        curr_widget =self.with_focus.canvas
        ind = self.canvas_list.index(curr_widget)
        lut = self.lut_list[ind]
        xmaxval=lut.xmaxval-lut.xminval
        lut.canvas.delete('isoline')
        lut.canvas.delete('isodot')
        lut.isoVals = []
        lut.draw_initshape(xminval = 0, xmaxval=xmaxval,
                           num_of_entries=xmaxval+1, gray=0, reverse = reverse)
        lut.values = [lut.xminval, lut.xminval,
                      lut.xminval+ xmaxval/2,
                       lut.xmaxval, lut.xmaxval]
        lut.callbacks['setalpha']([lut.xminval, lut.alpha_arr])
        lut.callbacks['setcolor']([lut.xminval, lut.color_arr])

    def set_dotsize_cb(self, type):
        """Set the size of the LUT widget's dots (type='dot') or
        isovalue dots (type='isodot')."""
        if type == "dot":
            text = "Set dot size"
            comm = lambda v, s=self: s.set_dotsize(v)
            self.currDotrad = self.lut_list[0].dotrad
        elif type == "isodot":
            text = "Set isovalue dot size"
            comm = lambda v, s=self: s.set_isodotsize(v)
            self.currDotrad = self.lut_list[0].isodotrad
        dotsizeDialog = Pmw.Dialog(self.master,
                                   buttons=('OK',), defaultbutton='OK',
                                   title=text,)# command=self.ok_dotsize)
        interior = dotsizeDialog.interior()
        
        sizescale = Tkinter.Scale(interior,
                                  label = "Set new dot radius",
                                  orient = 'horizontal',
                                  length = 200,
                                  from_= 2., to=10., tickinterval=0,
                                  font=self.font, resolution=0.5,
                                  command=comm)
        sizescale.pack(side=Tkinter.LEFT)
        sizescale.set(self.currDotrad)
        dotsizeDialog.activate()

    def set_dotsize(self, val):
        val = float(val)
        for lut in self.lut_list:
            lut.set_dotsize(val)

    def set_isodotsize(self, val):
        val = float(val)
        for lut in self.lut_list:
            lut.set_dotsize(val, "isodot")

    def ok_dotsize(self, val):
        print val

    def flip_function_cb(self):
        self.with_focus.flip_function()

    def del_selected_dot(self):
        lut = self.with_focus
        c = lut.canvas
        tags = c.gettags('selected')
        if len(tags) == 0: return
        if tags[0] == 'isodot':
            lut.removeIsoVal()
            num_isovals = 0
            for lut in self.lut_list:
                num_isovals = num_isovals + len(lut.isoVals)
            if num_isovals == 0:
                self.entry_isoVal.configure(label_text = 'ISO value (None):')
        elif tags[0] == 'dot':
            if lut.curr_ind-1 in lut.shapes and \
               lut.curr_ind+1 in lut.shapes:
                lut.deleteShape()
            else:
                lut.deleteDot()
                
    def show_help_txt(self, txt):
        text_opts = {"addDotShape":AddDotShapeTxt, "deleteDot":DeleteDotTxt ,
                     "moveShape":MoveShapeTxt, "setColor": SetColorTxt,
                     "description":DescriptionText}
        opt = text_opts[txt]
        helpDialog = Pmw.TextDialog(self.master, scrolledtext_labelpos='n',
                                    title = "Widget description",
                                    defaultbutton=0,
                                    scrolledtext_usehullsize=1,
                                    scrolledtext_hull_width=600,
                                    scrolledtext_hull_height=200,
                                    label_text="Widget description")
        helpDialog.insert(Tkinter.END, opt)
        helpDialog.configure(text_state='disabled')
        #helpDialog.activate()
        helpDialog.show()
        
    def set_continuous(self):
        val = self.continVar.get()
        LUT.continuous = val
        if val:
            self.applylutBtn.configure(state="disabled")
        else:
            self.applylutBtn.configure(state="normal")

    def applylut_cb(self):
        for lut in self.lut_list:
            lut.apply_lut()
        alpha_arr = self.lut_list[0].alpha_arr
        color_arr = self.lut_list[0].color_arr
        if len(self.lut_list)>1:
            for lut in self.lut_list[1:]:
                alpha_arr = Numeric.concatenate((alpha_arr, lut.alpha_arr))
                color_arr = Numeric.concatenate((color_arr, lut.color_arr))
        self.lutcallbacks['setalpha'] ([0, alpha_arr])
        self.lutcallbacks['setcolor'] ([0, color_arr])

    def merge_function (self):
        if len(self.intervals_list)<2:
            return
        interval = (self.xminval, self.xmaxval)
        #data ={'intervals':[interval,]}
        points = self.lut_list[0].points[1:-1]
        values = self.lut_list[0].values[1:-1]
        shapes = self.lut_list[0].shapes[1:-1]
        color_arr = self.lut_list[0].color_arr
        alpha_arr = self.lut_list[0].alpha_arr
        
        for lut in self.lut_list[1:]:
            points.extend(lut.points[1:-1])
            
            alpha_arr = Numeric.concatenate((alpha_arr, lut.alpha_arr))
            color_arr = Numeric.concatenate((color_arr, lut.color_arr))
            old_shapes =lut.shapes[1:-1]
            old_vals = lut.values[1:-1]
            val1 = values[-1]
            val2 = old_vals[0]
            if alpha_arr[val1] == 0 and alpha_arr[val2] == 0:
                d = shapes[-1]
                for s in old_shapes:
                    shapes.append(s+d)
            else:
                d = shapes[-1]
                shapes[-1] = d+old_shapes[1]
                if len(old_shapes) > 2:
                    for s in old_shapes[2:]:
                        shapes.append(s+d)
            values.extend(old_vals)    
        shapes.insert(0,0)
        shapes.append(shapes[-1]+1)
        ########
        for canvas in self.canvas_list:
            canvas.destroy()
        self.canvas_list = []
        self.lut_list = []
        self.intervals_list = [interval,]
        i=0
        lut = LUT(self.f1, xminval=interval[0], ymaxval=self.ymaxval,
                  xmaxval=interval[1], grid_row=i,
                  num_of_entries=self.xmaxval+1, initdraw = 0)
        lut.canvas.bind('<Button-1>', self.selected)
        self.lut_list.append(lut)
        self.canvas_list.append(lut.canvas)
        right = lut.right
        left = lut.left
        bott = lut.bott
        sclx = lut.sclx
        new_points = [(left, bott),]
        for i in range(len(points)):
            new_points.append(((values[i]-self.xminval)*sclx+left,
                                  points[i][1]))
        new_points.append((right, bott))
        lut.points = new_points
        lut.shapes = shapes
        lut.color_arr = color_arr
        lut.alpha_arr = alpha_arr
        values.insert(0, self.xminval)
        values.append(self.xmaxval)
        lut.values = values
        lut.redraw()
        lut.setCallbacks(self.lutcallbacks)
        self.lut_list[0].canvas.itemconfigure('outline', outline='blue')
        self.with_focus = self.lut_list[0]
        #######
        self.editBtn.menu.entryconfigure("Merge function", state='disabled')
        #data[interval] = [new_points, shapes, color_arr, alpha_arr]

        
class LUT:
    #callbacks = {'entries':None, 'rgbmap':None, 'setalpha':None,
    #             'setcolor':None}
                       #'entries' - TableManager.entries_update; 
                       #'rgbmap' - Colormap.show_color;
                       #'setalpha'-  should be reserved for setting
                       #alpha values, as in transferCommands.py(setAlpha):
                       #    setAlpha(value)  ,
                       #   where value is a list: value[0] - starting index,
                       #   value[1] - an integer array containing new
                       #   alpha values. Each entry is in range 0...4095.
                       #'setcolors' - for setting colors :
                       #    setColors(value) ,
                       #   where value is a list: value[0] - starting index,
                       #   value[1] - an  array of RGB values.
                       #   Each entry is in range 0.0 ... 1.0 
    continuous = 1
    
    def __init__(self, master=None, xminval=0,xmaxval=255, ymaxval=4095,
                 width = '13c', height = '3c', grid_row=0,
                 num_of_entries = 256, initdraw = 1):
        """create canvas widget """
        
        self.callbacks = {'entries':None, 'rgbmap':None, 'setalpha':None,
                 'setcolor':None}
        self.master = master
        self.xmaxval = xmaxval
        self.xminval = xminval
        self.ymaxval = ymaxval
        self.num_of_entries = num_of_entries
        assert self.num_of_entries > 0
        self.font = '-*-Helvetica-Bold-R-Normal-*-*-160-*-*-*-*-*-*'

        self.canvas = Tkinter.Canvas(master,width = width, height=height,
                             relief = Tkinter.FLAT, borderwidth = 2)
        c = self.canvas
        c.grid(row=grid_row, sticky=Tkinter.W)
        
        width = float(c.cget('width')) # 461.0
        height = float(c.cget('height'))#106.0
        diff = c.winfo_fpixels('0.5c')
        self.right = width-diff
        self.top = diff
        self.left = diff
        self.bott = height-diff
        
        c.create_rectangle(3,3, width+3,height+3,
                           outline ='', width=4,
                           fill='', tags='outline')
        c.create_rectangle(self.left,self.top, self.right,self.bott,
                           outline ='black', width=1,
                           fill='white', tags='box')
        c.create_text(self.right-15, self.bott+10, text=str(xmaxval),
                      anchor=Tkinter.W, font=self.font)
        c.create_text(10, self.bott+10, text=str(xminval), anchor=Tkinter.W,
                      font=self.font)
        c.create_text(10, self.top-7, text=str(ymaxval), anchor=Tkinter.W,
                      font=self.font, tags = 'ytext')
        
        # scale x and y axes 
        self.sclx = (self.right-self.left)/(xmaxval-xminval)
        self.scly = (self.bott-self.top)/ymaxval
        self.delta = self.sclx+1e-6
##          c.create_line(self.left,self.bott,self.right,
##                        self.bott, width = 2, fill='black', tags=('line','line0'))
##          self.line_count = 1
##          self.dot_count = 0

        c.tag_bind('dot', '<Button-1>', self.mouse1Down)
        c.tag_bind('dot', '<ButtonRelease-1>', self.mouse1Up)
        c.tag_bind('isodot', '<Button-1>', self.pickIsoVal)
        c.tag_bind('isodot', '<B1-Motion>', self.moveIsoVal)
        #c.tag_bind('line', '<Button-1>', self.pick)
        self.bind_tags()
        arr_len = xmaxval-xminval+1
        self.color_arr = Numeric.zeros((arr_len,3),'f')
        self.alpha_arr = Numeric.zeros(arr_len).astype(Numeric.Int16)
        self.points = []
        self.values = []
        self.shapes = []
        self.isoVals = []
        self.dotrad = 5
        self.isodotrad = 5
        self.isoval_inrange = []
        if initdraw:
            self.draw_initshape()
        self.last_event = ''
        self.curr_ind = 1
        #self.continuous = 1


    def setCallbacks(self, func_dict):
        alpha = func_dict.get('setalpha')
        color = func_dict.get('setcolor')
        entr = func_dict.get('entries')
        rgbmap = func_dict.get('rgbmap')
        if alpha:
            self.callbacks['setalpha']=alpha
        if color:
            self.callbacks['setcolor']=color
        if entr:
            self.callbacks['entries']=entr
        if rgbmap:
            self.callbacks['rgbmap']=rgbmap

    def compute_alpha_ramp(self, xminval, num_of_entries, arr_len,
                           reverse = 0):
        #print "xminval=", xminval
        #print "numentr:", num_of_entries
        #print "arr_len:", arr_len
        assert num_of_entries > 1
        #d = self.ymaxval/(num_of_entries-1.0)
        #print "d:", d
        ymaxval = self.ymaxval
        n = num_of_entries-1.0
        if not reverse:
            for i in range(arr_len):
                #self.alpha_arr[i] = (xminval+i)*d
                self.alpha_arr[i] = ymaxval * (xminval+i)/ n
        else:
            j = arr_len-1
            for i in range(arr_len):
                #self.alpha_arr[i] = (xminval+j)*d
                self.alpha_arr[i] = ymaxval * (xminval+j)/n
                j=j-1
##          print self.alpha_arr[0], self.alpha_arr[-1]

    def draw_initshape(self, xminval = None, xmaxval = None,
                       num_of_entries=None, gray=1, reverse = 0):
        if xmaxval is None:
            xmaxval = self.xmaxval
        if xminval is None:
            xminval = self.xminval
        if num_of_entries is None:
            num_of_entries = self.num_of_entries
        arr_len = xmaxval-xminval+1
        self.compute_alpha_ramp(xminval, num_of_entries, arr_len, reverse)
        first_dot = (self.left,self.bott-self.scly*self.alpha_arr[0] )
        last_dot = (self.right,self.bott-self.scly*self.alpha_arr[-1])
        mid_dot = (self.left+(self.right-self.left)/2,
                   first_dot[1]-(first_dot[1]-last_dot[1])/2)
        self.points = [(self.left,self.bott),
                       first_dot, mid_dot, last_dot,
                       (self.right, self.bott)]
        #print 'points=', self.points
        self.values = [xminval, xminval, xminval+(xmaxval-xminval)/2,
                       xmaxval, xmaxval]
        self.shapes = [0,1,3,4]
        if gray:
            self.grayRamp(xminval, xmaxval)
        self.redraw(xminval, xmaxval)
        

    def grayRamp(self, xminval=None, xmaxval=None):
        if xminval is None:
            xminval = self.xminval
        if xmaxval is None:
            xmaxval = self.xmaxval
        n = self.num_of_entries/256.
##          print 'n=',n
        j = int(xminval/n)
        k=0
        for i in range(xminval, xmaxval+1):
            if i%n == 0 and i != xminval : j = j+1 
            c = j/255.
            self.color_arr[k] = [c,c,c]
            k = k+1
##          print self.color_arr

    
    def check_bounds(self, x, l_end, r_end):
        diff = r_end - l_end
        diff1 = r_end-x
        if diff < 10: return 0
        elif diff <= 30:
            x1=l_end+3
            x2=r_end-3
        elif diff1 < 8: return 0
        elif diff1 <= 30:
            x1=x
            x2=r_end-3
        else:
            x1=x
            x2=x+30
        return x1,x2

    def update_linetags(self, ind1, ind2):
        """update tags of the lines when a dot or a shape added to/
        (deleted from) the canvas"""
        line_ind = range(ind1, self.line_count)
        if ind2 > 0:
            line_ind.reverse()
        for i in line_ind:
            self.canvas.itemconfigure('line%d'%(i,),
                                      tags=('line', 'line%d'%(i+ind2,)))

    def update_dottags(self, ind1, ind2):
        """update tags of the dots when a dot or a shape added to/
        (deleted from) the canvas"""
        dot_ind =range(ind1, self.line_count)
        if len(dot_ind) == 0: return
        r_edge = None
        if ind2 > 0:
            dot_ind.reverse()
            if 'edge' in self.canvas.gettags('dot%d'%(dot_ind[0],)):
              self.canvas.itemconfigure('dot%d'%(dot_ind[0],),
                           tags=('dot', 'dot%d'%(dot_ind[0]+ind2,),'edge'))
              dot_ind.pop(0)
        else:
            if 'edge' in self.canvas.gettags('dot%d'%(dot_ind[-1],)):
              r_edge = dot_ind.pop(-1)
        for i in dot_ind:
            self.canvas.itemconfigure('dot%d'%(i,),
                                      tags=('dot', 'dot%d'%(i+ind2,)))
        if r_edge:
            self.canvas.itemconfigure('dot%d'%(r_edge,),
                           tags=('dot', 'dot%d'%(r_edge+ind2,),'edge'))
            
    def drawShape(self,event):
        """ Draw a shape that consists of four dots
        and five lines, find coordinates and values of the points
        and update the entry fields on the form"""
        self.last_event = 'Shift-Button-1'
        c = self.canvas
        # memorize the x and y coordinates of the cursor
        x = c.canvasx(event.x)
        y = c.canvasy(event.y)
        points = self.points
        shapes = self.shapes
        # find point where shape will be placed
        i=0
        for point in points:
            if x < point[0]:
                r_point = point
                break
            i=i+1
        ind = i-1
        l_point = points[ind]
        # if pointer is over existing shape - do nothing
        if ind in shapes:
            newind=shapes.index(ind)
            rind=shapes[newind+1]
            if (rind-ind)>1: return
        else: return
        bounds = self.check_bounds(x, l_point[0], r_point[0])
        if bounds == 0: return
        else : x1, x2 = bounds
        c.delete('line%d'%(ind,))
##          c.itemconfig('selected', fill=UNSELECTED_COLOR,
##                       outline=UNSELECTED_COLOR)
        #c.itemconfig('selected', outline='')
        c.itemconfig('selected', outline=UNSELECTED_COLOR)
        c.dtag('selected')
        if self.line_count > 1: self.update_linetags(ind+1,4)
        if self.dot_count > 1: self.update_dottags(ind+1,4)
        bott = self.bott
        
        # draw lines of the shape
        eps =1e-6
        x1_1 = x1+self.sclx+eps
        x2_1 = x2-self.sclx+eps

        c.create_line(l_point[0], bott, x1, bott, width = 2,
                      activefill='gray68',
                      fill='black', tags= ('line', 'line%d'%(ind,)))
        c.create_line(x1, bott, x1_1, y, width = 2, fill='black',
                      activefill='gray68',
                      tags=('line', 'line%d'%(ind+1,)))
        c.create_line(x1_1, y, x2_1, y,
                      width = 2, fill='black',
                      activefill='gray68',
                      tags=('line','line%d'%(ind+2,)))
        c.create_line(x2_1, y, x2, bott, width = 2, fill='black',
                      activefill='gray68',
                      tags=('line','line%d'%(ind+3,)))
        c.create_line(x2, bott, r_point[0], bott, width = 2,
                      activefill='gray68',
                      fill='black', tags=('line','line%d'%(ind+4,)))

        # update list of values and list of points
        xminval =self.xminval
        tkcolors = []
        val_x1 = int(round((x1-self.left)/self.sclx))+xminval
        val_x2 = val_x1+int(round((x2-x1)/self.sclx))        
        for j, point, val in ([i,(x1,bott),val_x1], [i+1,(x1_1,y),val_x1+1],\
                         [i+2,(x2_1,y),val_x2-1], [i+3,(x2,bott),val_x2]):
            tkcolors.append(colorUtil.TkColor(self.color_arr[val-xminval]))
            points.insert(j,point)
            self.values.insert(j,val)

        # draw dots
        rad = self.dotrad
        c.create_oval(x1-rad, bott-rad, x1+rad, bott+rad,
                      outline=UNSELECTED_COLOR,
                      width=1, fill=tkcolors[0],
                      activefill='gray68',
                      tags=('dot','dot%d'%(ind+1,)))
        c.create_oval(x1_1-rad, y-rad,
                      x1_1+rad, y+rad,
                      outline=UNSELECTED_COLOR, width=1,
                      fill=tkcolors[1], tags=('dot', 'dot%d'%(ind+2,)))
        c.create_oval(x2_1-rad, y-rad,
                      x2_1+rad, y+rad,
                      activefill='gray68',
                      outline=UNSELECTED_COLOR, width=1,
                      fill=tkcolors[2], tags=('dot','dot%d'%(ind+3,)))
        c.create_oval(x2-rad, bott-rad, x2+rad, bott+rad,
                      outline=UNSELECTED_COLOR,
                      activefill='gray68',
                      width=1,fill=tkcolors[3],
                      tags=('dot','dot%d'%(ind+4,)))
        
        # update list of shapes
        self.update_shapes(i,4)
        self.line_count = self.line_count+4
        self.dot_count = self.dot_count+4
        if ind > 0:
            c.lift('dot%d'%(ind,))
        if ind+5 < self.line_count:
            c.lift('dot%d'%(ind+5,))
            
        # calculate alphas
        val_y = round((bott-y)/self.scly)
        for j in range(val_x1+1-xminval, val_x2-xminval):
            self.alpha_arr[j] = val_y
        
        # update entry fields on the form
        self.callbacks['entries'](val_x1=val_x1, val_x2=val_x2,
                            val_y=int(val_y)) # , val_r='', val_g='', val_b='')
        self.callbacks['setalpha']([val_x1+1,
                           self.alpha_arr[val_x1+1-xminval:val_x2-xminval]])

    def drawIsoVal(self, isoVal, isoA, isoRGB):
        """Draws a line with a circle representing an Iso Value"""
        
        y = self.bott-self.scly*isoA
        x = (isoVal - self.xminval)*self.sclx + self.left
        tkcolor = colorUtil.TkColor(isoRGB)
        n = len(self.isoVals)
        rad = self.isodotrad
        
        self.canvas.itemconfig('selected', outline=UNSELECTED_COLOR)
        self.canvas.dtag('selected')
        self.canvas.create_line(x, self.bott, x, y, fill="black",
                                activefill='gray68',
                                tags=('isoline', 'isoline%d'%(n,)) )
        self.canvas.create_oval( x-rad, y-rad, x+rad, y+rad,
                                 fill=tkcolor, outline=SELECTED_COLOR,
                                 activefill='gray68',
                                 tags=('isodot','isodot%d'%(n,),'selected'))
        self.isoVals.append({'val':isoVal, 'alpha':isoA,
                             'rgb':isoRGB})
        if self.continuous:
            self.callbacks['setalpha']([isoVal,
                            Numeric.array([isoA]).astype(Numeric.Int16)])
            self.callbacks['setcolor']([isoVal, Numeric.array([isoRGB],'f')])
        self.callbacks['entries'](iso_val=isoVal, iso_y=isoA, iso_rgb = isoRGB)
        

    def removeIsoVal(self):
        c = self.canvas
        tags = c.gettags('selected')
        if len(tags) == 0: return
        if tags[0] != 'isodot': return
        iso_ind = int(tags[1][6:])
        c.delete('selected')
        c.delete('isoline%d'%(iso_ind,))
        isoVal = self.isoVals.pop(iso_ind)
        dots_id = c.find_withtag('isodot')
        lines_id = c.find_withtag('isoline')
        for i , j in map(None, dots_id, lines_id):
            d_ind = int(c.gettags(i)[1][6:])
            l_ind = int(c.gettags(j)[1][7:])
            if d_ind > iso_ind:
                c.itemconfigure(i, tags = ('isodot','isodot%d'%(d_ind-1,)))
            if l_ind > iso_ind:
                c.itemconfigure(j, tags = ('isoline','isoline%d'%(l_ind-1,)))
        val = isoVal['val']
        alpha = self.alpha_arr[val-self.xmaxval]
        rgb = self.color_arr[val-self.xmaxval]
        if self.continuous:
            self.callbacks['setalpha']([val,
                                Numeric.array([alpha]).astype(Numeric.Int16)])
            self.callbacks['setcolor']([val, Numeric.array([rgb],'f')])
        self.callbacks['entries'](iso_val=0, iso_y=0, iso_rgb = (0.0, 0.0, 0.0))
        
    def flip_function(self):
        if self.values[0] !=self.values[1]:
            val0 = self.values[0]
            self.values.insert(0, val0)
            point0 = self.points[0]
            self.points.insert(0, point0)
        if self.values[-1] !=self.values[-2]:
            val1 = self.values[-1]
            self.values.append(val1)
            point1 = self.points[-1]
            self.points.append(point1)

        y_points = []
        for point in self.points:
            y_points.append(point[1])
        y_min = min(y_points)
        d = self.bott+y_min
        for i in range(1,len(self.points)-1):
            px, py = self.points[i]
            self.points[i] = (px, d-py)
        alpha_max= Numeric.maximum.reduce(self.alpha_arr)
        self.alpha_arr = (self.alpha_arr*(-1)+alpha_max).astype(Numeric.Int16)
        shapes = [0,1]
        num_points = len(self.points)
        if len(self.shapes) > 4:
            if num_points >=8:
                i = 3
                while i < num_points-4:
                    val = self.values[i]-self.xminval
                    if self.alpha_arr[val] == 0:
                        val1 = self.values[i+1]-self.xminval
                        if self.alpha_arr[val1] == 0:
                            shapes.append(i)
                            shapes.append(i+1)

                            self.points[i] = (self.points[i][0], self.bott)
                            self.points[i+1] = (self.points[i+1][0], self.bott)
                            i = i+3
                        else:
                            i=i+1
                    else:
                        i=i+1
                if num_points - shapes[-1] > 3:
                    shapes.append(num_points-2)
            else:
                shapes.append(num_points-2)
        else:
            shapes.append(num_points-2)
        shapes.append(num_points-1)
        self.shapes = shapes
        self.redraw(redrawIsoVals=0)
        self.callbacks['setalpha'] ([self.xminval, self.alpha_arr])
        
        
    def redraw(self, xminval = None, xmaxval = None, redrawIsoVals = 1):
        """redraw canvas widget after interval splitting"""
        c = self.canvas
        if xminval is None:
            xminval =self.xminval
        if xmaxval is None:
            xmaxval = self.xmaxval
##          c.delete('line0')
        c.delete('line')
        c.delete('dot')
        points = self.points
        p_len = len(points)
        if p_len == 2:
            c.create_line(points[0][0], points[0][1],
                          points[1][0], points[1][1],
                          width = 2, fill='black',
                          activefill='gray68',
                          tags= ('line', 'line0'))
##              print 'points=', self.points
##              print 'shapes=', self.shapes
##              print 'values=', self.values
            self.line_count = p_len-1
            self.dot_count = p_len-2
            return 
        rad = self.dotrad
        for i in range(p_len-1):
            c.create_line(points[i][0], points[i][1],
                          points[i+1][0], points[i+1][1],
                          activefill='gray68',
                          width = 2, fill='black',
                          tags= ('line', 'line%d'%(i,)))
        i = 1
        for point in points[1:-1]:
            tkcolor = colorUtil.TkColor(self.color_arr[self.values[i]-xminval])
            c.create_oval(point[0]-rad, point[1]-rad,
                          point[0]+rad, point[1]+rad,
                          fill=tkcolor, outline=UNSELECTED_COLOR,
                          activefill='gray68',
                          tags=('dot','dot%d'%(i,)))
            i = i+1
        #if points[-2][0] == self.right:
        if self.values[-2] == xmaxval:
            c.addtag_withtag('edge', 'dot%d'%(p_len-2,))
        #if points[1][0] == self.left:
        if self.values[1] == xminval:
            c.addtag_withtag('edge', 'dot1')
        self.line_count = p_len-1
        self.dot_count = p_len-2
        if redrawIsoVals:
            if len(self.isoVals):
                rad = self.isodotrad
                i = 0
                for iso_val in self.isoVals:
                    y = self.bott-self.scly*iso_val['alpha']
                    x = (iso_val['val'] - xminval)*self.sclx + self.left
                    tkcolor = colorTool.TkColor(iso_val['rgb'])
                    self.canvas.create_line(x, self.bott, x, y, fill="black",
                                            activefill='gray68',
                                            tags=('isoline', 'isoline%d'%(i,)) )
                    self.canvas.create_oval( x-rad, y-rad, x+rad, y+rad,
                                      fill=tkcolor, outline=UNSELECTED_COLOR,
                                      activefill='gray68',
                                             tags=('isodot','isodot%d'%(i,)) )
                    i = i + 1

    def update_yaxis(self, ymaxval):
        self.scly = (self.bott-self.top)/ymaxval
        self.canvas.itemconfigure('ytext', text = str(ymaxval))
        #self.alpha_arr = (self.alpha_arr*ymaxval/self.ymaxval).astype(Numeric.Int16)
        self.alpha_arr = Numeric.clip(self.alpha_arr, 0, ymaxval).astype(Numeric.Int16)
        for i in range(1, len(self.points)-1):
            arr_ind = self.values[i]-self.xminval
            alpha = self.alpha_arr[arr_ind]
            self.points[i] = (self.points[i][0], 
                              self.bott-self.scly*alpha)
        self.ymaxval = ymaxval
        self.redraw()
        self.callbacks['setalpha']([self.xminval, self.alpha_arr])
        ind = self.values[self.curr_ind]-self.xminval
        self.callbacks['entries'](val_y=int(self.alpha_arr[ind]))

    def update_shapes(self, point_ind, num):
        """Update list of shapes """
        if point_ind not in self.shapes:
            print "point %d is not in self.shapes"%point_ind
            return 0
        shape_ind = self.shapes.index(point_ind)
        if num > 0:
            self.shapes.insert(shape_ind,point_ind)
            self.shapes.insert(shape_ind+1,point_ind+num-1)
            delta = 2
        elif num <0:
            self.shapes.pop(shape_ind)
            self.shapes.pop(shape_ind)
            delta = 0
        shapes_len = len(self.shapes)
        for i in range(shape_ind+delta, shapes_len):
            self.shapes[i]=self.shapes[i]+num

    def mouse1Down(self, event):
        self.last_event = 'Button-1'
        #print 'last_event =', self.last_event
        c= self.canvas
        tag = c.gettags(Tkinter.CURRENT)[1]
        curr_ind = int(tag[3:])
        
        #memorize last position of the mouse-click on a point(dot)
        self.lastx = self.points[curr_ind][0]
        self.lasty = self.points[curr_ind][1]
        
        #find the scalar values and alpha values of the selected point
        #and its adjoining points
        self.next_l = self.points[curr_ind-1] #coords of left neighbor. point
        self.next_r = self.points[curr_ind+1] #coords of right neighbor. point
        scalars = [self.values[curr_ind-1], self.values[curr_ind],
                   self.values[curr_ind+1]]
        xminval = self.xminval
        self.alphas = [self.alpha_arr[scalars[0]-xminval],
                       self.alpha_arr[scalars[1]-xminval],
                       self.alpha_arr[scalars[2]-xminval]]
        
        # change outline of the selected dot
        c.itemconfig('selected', outline=UNSELECTED_COLOR)
        c.dtag('selected')
        c.addtag_withtag('selected', Tkinter.CURRENT)
        c.itemconfig('selected', outline=SELECTED_COLOR)

        # get color of the selected (current) point
        curr_color = self.curr_color = list(self.color_arr[scalars[1]-xminval])
        curr_hsv = ToHSV(curr_color)
        #print "curr_hsv in mouse1Down:", curr_hsv
        self.scalars = scalars
        self.curr_ind = curr_ind
        self.curr_ind_status = None
        if curr_ind in self.shapes:
            if 'edge' in c.gettags('dot%d'%(curr_ind,)):
                self.curr_ind_status = 'edge'
            else : self.curr_ind_status = 'shape'

        # update entry fields on the form
        self.callbacks['entries'](val_x1=scalars[1], val_x2='',
                            val_y=int(self.alphas[1]),
                            val_r=round(curr_color[0],2),
                            val_g=round(curr_color[1],2),
                            val_b=round(curr_color[2],2))
        self.callbacks['rgbmap'](curr_hsv)
        #check if there is an Iso val in the interval scalars[0]...scalars[2]
        self.isoval_inrange = []
        if len(self.isoVals):
            for iso_val in self.isoVals:
                v = iso_val['val']
                if v >= scalars[0] and v <= scalars[2]:
                    self.isoval_inrange.append(iso_val)

    def mouse1Move(self, event):
        self.last_event = 'B1-Motion'
        curr_ind = self.curr_ind
        c = self.canvas
        #current position of the mouse
        x = c.canvasx(event.x)
        y = c.canvasy(event.y)
        
        next_l = self.next_l
        next_r = self.next_r
        delta = self.delta
        top = self.top
        curr_tag = 'dot%d'%(curr_ind,)

        # the dot must not move beyond the canvas outline
        if y >= self.bott: y = self.bott
        elif y > top-15 and y <= top: y = top
        elif y <= top-15:
            if curr_ind in self.shapes: return
            #draggind the dot beyond the upper boundary delets either the shape
            #(if there is only three points in the shape) or the dot (if
            #there is more than three points in the shape).
            elif curr_ind-1 in self.shapes and curr_ind+1 in self.shapes:
                self.deleteShape()
            else: self.deleteDot()
            return
        
        if x < next_l[0]+delta:
            if curr_ind-1 == 0: x = self.left
            else: x = next_l[0]+delta
        elif x > next_r[0]-delta:
            if curr_ind+1 == len(self.points)-1: x = self.right
            else: x = next_r[0]-delta
            
        if self.curr_ind_status != None:
            if self.curr_ind_status == 'edge':
                lasty = self.lasty
                lastx = self.lastx
                if lasty < self.bott:
                    c.move(curr_tag, 0, y-lasty)
                else :
                    if abs(y-lasty)>abs(x-lastx):
                        c.move(curr_tag, 0, y-lasty)
                    else:
                        if next_r[0] == self.right: self.alphas[2] = 0
                        elif next_l[0] == self.left :
                            self.alphas[0] = 0
                        c.move(curr_tag, x-lastx, 0)
                        self.curr_ind_status = 'shape'
                        c.dtag(curr_tag, 'edge')
            elif self.curr_ind_status == 'shape':
                c.move(curr_tag, x-self.lastx, 0)
        else:
            c.move(curr_tag, x-self.lastx, y-self.lasty)

        dot = c.coords('dot%d'%(curr_ind,))
        rad = self.dotrad
        X = dot[0]+rad
        Y = dot[1]+rad
        c.coords('line%d'%(curr_ind-1,), next_l[0], next_l[1], X, Y)
        c.coords('line%d'%(curr_ind,), X, Y, next_r[0], next_r[1])
        if x == self.right:
            self.lastx = x
        elif x == self.left and dot[0]<self.left:
            # the case when the center of the current dot is on
            #the left canvas border
            self.lastx = self.left
        else: self.lastx = X

        self.lasty = y
        c.update_idletasks()

        # calculate alphas and colors in the interval between two neighboring
        # dots of the current dot (scalars[0] ... scalars[2])
        xminval = self.xminval
        scalars = self.scalars
        alphas = self.alphas
        scalars[1] = int(round((X-self.left)/self.sclx))+xminval
        if scalars[1] >= scalars[2]:
            if scalars[2] == self.xmaxval:
                scalars[1] = self.xmaxval
            else :
                scalars[1] = scalars[2]-1
        elif scalars[1] <= scalars[0]:
            if scalars[0] == xminval:
                scalars[1] = xminval
            else:
                scalars[1] = scalars[0]+1
        self.alphas[1] = round((self.bott-Y)/self.scly)
        self.values[self.curr_ind] = scalars[1]
        self.compute_rgba(scalars, alphas)
        # update entry fields on the form
        self.callbacks['entries'](val_x1=scalars[1], val_y=int(alphas[1]))

    def compute_rgba(self, scalars, alphas):
        xminval = self.xminval
        xmaxval = self.xmaxval
        sc0 = scalars[0]-xminval
        sc1 = scalars[1]-xminval
        sc2 = scalars[2]-xminval
        red1, green1, blue1 = self.color_arr[sc0]
        red2, green2, blue2 = self.curr_color
        red3, green3, blue3 = self.color_arr[sc2]
        d = scalars[0]-scalars[1]
        if d != 0 :
            K = (alphas[0]-alphas[1])/d
            C = alphas[0]-K*scalars[0]
            Kr = (red1-red2)/d
            Kg = (green1-green2)/d
            Kb = (blue1-blue2)/d
            Cr = red1-Kr*scalars[0]
            Cg = green1-Kg *scalars[0]
            Cb = blue1-Kb*scalars[0]
            s = Numeric.arange(scalars[0], scalars[1])
            self.alpha_arr[sc0:sc1] =\
                (K*s+C).astype(Numeric.Int16)
            if scalars[1] > xminval or scalars[1] < xmaxval:
                self.color_arr[sc0:sc1] = \
                         Numeric.transpose(Numeric.array((Kr*s+Cr, Kg*s+Cg,
                                                    Kb*s+Cb)).astype('f'))
        d = scalars[1]-scalars[2]
        if d != 0 :
            K = (alphas[1]-alphas[2])/d
            C = alphas[1]-K*scalars[1]
            Kr = (red2-red3)/d
            Kg = (green2-green3)/d
            Kb = (blue2-blue3)/d
            Cr = red2-Kr*scalars[1]
            Cg = green2-Kg *scalars[1]
            Cb = blue2-Kb*scalars[1]
            s = Numeric.arange(scalars[1], scalars[2])
            self.alpha_arr[sc1:sc2] =\
                (K*s+C).astype(Numeric.Int16)
            if scalars[1] > xminval or scalars[1] < xmaxval:
                self.color_arr[sc1:sc2] = \
                   Numeric.transpose(Numeric.array((Kr*s+Cr, Kg*s+Cg,
                                                    Kb*s+Cb)).astype('f'))
        else:
            s = sc1
            last_ind = xmaxval-xminval
            if s > last_ind: s=last_ind
            self.alpha_arr[s]=alphas[1]
              
        #if there is any Iso vals in the range:
        # memorize current rgba values , make changes in self.alpha_arr
        # and self.color_arr, do callbacks with new alphas and colors
        # then assign old values to self.alpha_arr and self.color_arr
        old_rgba = None
        if len(self.isoval_inrange):
            old_rgba = []
            for iso_val in self.isoval_inrange:
                v = iso_val['val']
                old_rgba.append((v, self.alpha_arr[v-xminval],
                                 self.color_arr[v-xminval].tolist()))
                self.alpha_arr[v-xminval] = iso_val['alpha']
                self.color_arr[v-xminval] = iso_val['rgb']
        if self.continuous:
            self.callbacks['setalpha']([scalars[0], self.alpha_arr[sc0:sc2+1] ])
        if scalars[1] > xminval or scalars[1] < xmaxval:
            self.color_arr[sc1] = [red2, green2, blue2]
##              curr_color = self.color_arr[sc1]
##              tkcolor = colorUtil.TkColor(curr_color)
##              l_color = colorUtil.TkColor(self.color_arr[sc0])
##              r_color = colorUtil.TkColor(self.color_arr[sc2])
##              self.canvas.itemconfig('dot%d'%(self.curr_ind+1,), fill=r_color)
##              self.canvas.itemconfig('dot%d'%(self.curr_ind-1,), fill=l_color)
##              self.canvas.itemconfig('selected', fill=tkcolor)
            if self.continuous:
                self.callbacks['setcolor']([scalars[0], self.color_arr[sc0:sc2+1] ])
        if old_rgba:
            for rgba in old_rgba:
                v = rgba[0]
                self.alpha_arr[v-xminval] = rgba[1]
                self.color_arr[v-xminval] = rgba[2]
                
        #self.values[self.curr_ind] = scalars[1]
         # update entry fields on the form
        #self.callbacks['entries'](val_x1=scalars[1], val_y=int(alphas[1]))
        #self.alphas = alphas    

    def mouse1Up(self, event):
        self.last_event = 'ButtonRelease-1'
        curr_ind = self.curr_ind
        if self.lastx == self.right or self.lastx == self.left:
            if 'edge' not in self.canvas.gettags('dot%d'%(curr_ind,)):
                self.canvas.addtag_withtag('edge', 'dot%d'%(curr_ind,))
        if self.lastx != self.points[curr_ind][0] or \
           self.lasty != self.points[curr_ind][1]:
            dot = self.canvas.coords('dot%d'%(curr_ind,))
            X = dot[0]+self.dotrad
            Y = dot[1]+self.dotrad
            self.points[curr_ind] = (X,Y)
        

    def pickIsoVal(self, event):
        c= self.canvas
        tag = c.gettags(Tkinter.CURRENT)[1]
        #print "tag:", tag
        self.iso_ind = int(tag[6:])
        #c.itemconfig('selected', outline='')
        c.itemconfig('selected', outline=UNSELECTED_COLOR)
        c.dtag('selected')
        c.addtag_withtag('selected', Tkinter.CURRENT)
        c.itemconfig('selected', outline=SELECTED_COLOR)
        rgb = self.isoVals[self.iso_ind]['rgb']
        self.callbacks['entries'](iso_val=self.isoVals[self.iso_ind]['val'],
                             iso_y=self.isoVals[self.iso_ind]['alpha'],
                             iso_rgb=rgb)
        hsv = ToHSV(rgb)
        self.callbacks['rgbmap'](hsv)

    def moveIsoVal(self, event):
        """moves Iso dot when it is draged with mouse button 1"""
        x = self.canvas.canvasx(event.x)
        y = self.canvas.canvasy(event.y)
        if y >= self.bott: y = self.bott
        elif y > self.top-25 and y <= self.top: y = self.top
        elif y <= self.top-25:
            self.removeIsoVal()
            return
        if x <= self.left: x = self.left
        elif x >= self.right: x = self.right
        # calculate new scalar value and opacity of moved point 
        val = int(round((x-self.left)/self.sclx))+ self.xminval
        alpha = int(round((self.bott-y)/self.scly))
        self.doMoveIsoVal(x, y ,val, alpha)
        self.callbacks['entries'](iso_val=val, iso_y=alpha)

    def doMoveIsoVal(self, x, y ,val, alpha):
        rad = self.isodotrad
        iso_ind = self.iso_ind
        xminval = self.xminval
        #find values of all Iso dots
        #vals = map(lambda i: i['val'], self.isoVals)
        isoVal = self.isoVals[iso_ind]
        # move the dot in canvas widget
        c = self.canvas
        c.coords('isoline%d'%(iso_ind,), x, self.bott, x, y)
        c.coords('isodot%d'%(iso_ind,), x-rad, y+rad, x+rad, y-rad)
        c.update_idletasks()
        # get last values of the dot
        oldval = isoVal['val']
        oldalpha = self.alpha_arr[oldval-xminval]
        oldrgb = self.color_arr[oldval-xminval]
        
        # update values for the dot
        isoVal['val'] = val
        isoVal['alpha'] = alpha
        rgb = isoVal['rgb']
        #find values of all Iso dots
        vals = map(lambda i: i['val'], self.isoVals)
        #update LUT
        if oldval < val: #dot moved to the right
            d = val-oldval
            alphas = Numeric.zeros(d+1).astype(Numeric.Int16)
            rgbs = Numeric.zeros((d+1,3),'f')
            alphas[0] = oldalpha
            alphas[-1] = alpha
            rgbs[0] = oldrgb
            rgbs[-1] = rgb
            if d > 1:
                alphas[1:-1] = self.alpha_arr[oldval-xminval+1:val-xminval]
                rgbs[1:-1] = self.color_arr[oldval-xminval+1:val-xminval]
            # check if the dot moved across another iso dot
            if len(vals)>1:
                for v in vals:
                    if v in range(oldval, val):
                        other_ind = vals.index(v)
                        alphas[v-oldval] = self.isoVals[other_ind]['alpha']
                        rgbs[v-oldval] = self.isoVals[other_ind]['rgb']
            if self.continuous:
                self.callbacks['setalpha']([oldval, alphas])
                self.callbacks['setcolor']([oldval, rgbs])
        elif oldval > val: #dot moved to the left
            d = oldval - val
            alphas = Numeric.zeros(d+1).astype(Numeric.Int16)
            rgbs = Numeric.zeros((d+1,3),'f')
            alphas[0] = alpha
            alphas[-1] = oldalpha
            rgbs[0] = rgb
            rgbs[-1] = oldrgb
            if d > 1:
                alphas[1:-1] = self.alpha_arr[val-xminval+1:oldval-xminval]
                rgbs[1:-1] = self.color_arr[val-xminval+1:oldval-xminval]
            # check if the dot moved across another iso dot
            if len(vals)>1:
                for v in vals:
                    if v in range(val+1, oldval+1):
                        other_ind = vals.index(v)
                        alphas[v-val] = self.isoVals[other_ind]['alpha']
                        rgbs[v-val] = self.isoVals[other_ind]['rgb']
            if self.continuous:
                self.callbacks['setalpha']([val, alphas])
                self.callbacks['setcolor']([val, rgbs])
        #self.callbacks['entries'](iso_val=val, iso_y=alpha)

    def moveIsoVal_value(self, val):
        """Move selected iso dot when vol entry is changed"""
        c = self.canvas
        tags = c.gettags('selected')
        if len(tags) == 0: #no dot is selected
            return
        if tags[0] != 'isodot': return
        x = (val - self.xminval)*self.sclx + self.left
        y = c.coords('isoline%d'%self.iso_ind)[3]
        alpha = self.isoVals[self.iso_ind]['alpha']
        self.doMoveIsoVal(x, y, val, alpha)
        
        
    def move_dot_alpha(self,alpha):
        """Move selected dot when alpha entry is changed"""
        c = self.canvas
        tags = c.gettags('selected')
        if len(tags) == 0: #no dot is selected
            self.callbacks['entries'](val_y='0')
            return
        if tags[0] != 'dot': return        
        curr_ind = self.curr_ind
        curr_tag = 'dot%d'%(curr_ind,)
        if curr_ind in self.shapes:
            if 'edge' not in c.gettags(curr_tag):
                self.callbacks['entries'](val_y='0')
                return
            
        next_l = self.points[curr_ind-1]
        next_r = self.points[curr_ind+1]        
        lasty = self.points[curr_ind][1]
        y = self.bott-self.scly*alpha
        c.move(curr_tag, 0, y-lasty)

        dot = c.coords(curr_tag)
        X = dot[0]+self.dotrad
        Y = dot[1]+self.dotrad
        c.coords('line%d'%(curr_ind-1,), next_l[0], next_l[1], X, Y)
        c.coords('line%d'%(curr_ind,), X, Y, next_r[0], next_r[1])
        self.lasty = y
        c.update_idletasks()

        # calculate alphas in the interval between two neighboring
        # dots of the current (selected)dot 
        xminval = self.xminval
        scalars = [self.values[curr_ind-1], self.values[curr_ind],
                   self.values[curr_ind+1]]
        alphas = [self.alpha_arr[scalars[0]-xminval],
                  round((self.bott-Y)/self.scly),
                  self.alpha_arr[scalars[2]-xminval]]
        d = scalars[0]-scalars[1]
        if d != 0 :
            K = (alphas[0]-alphas[1])/d
            C = alphas[0]-K*scalars[0]
            s = Numeric.arange(scalars[0], scalars[1])
            self.alpha_arr[scalars[0]-xminval:scalars[1]-xminval] =\
                    (K*s+C).astype(Numeric.Int16)
        d = scalars[1]-scalars[2]
        if d != 0 :
            K = (alphas[1]-alphas[2])/d
            C = alphas[1]-K*scalars[1]
            s = Numeric.arange(scalars[1], scalars[2])
            self.alpha_arr[scalars[1]-xminval:scalars[2]-xminval] =\
                    (K*s+C).astype(Numeric.Int16)
        else:
            s = scalars[1]-xminval
            last_ind = self.xmaxval-xminval
            if s > last_ind: s=last_ind
            self.alpha_arr[s]=alphas[1]
##          f=s+1-xminval
##          try:
##              self.alpha_arr[f] = K*(s+1)+C
##          except:
##              print 'index out of bounds:self.alpha_arr[%d]' % f
##              print 'scalar[1]=%d, scalar[2]=%d'%(scalars[1],scalars[2])
            
        #if there is any Iso vals in the range:
        # memorize current alpha values , make changes in self.alpha_arr,
        # do callbacks with new alphas ,
        # then assign old values to self.alpha_arr
        old_alpha = None
        if len(self.isoval_inrange):
            old_alpha = []
            for iso_val in self.isoval_inrange:
                v = iso_val['val']
                old_alpha.append((v, self.alpha_arr[v-xminval]))
                self.alpha_arr[v-xminval] = iso_val['alpha']
        if self.continuous:
            self.callbacks['setalpha']([scalars[0],
                   self.alpha_arr[scalars[0]-xminval:scalars[2]+1-xminval]])
        if old_alpha:
            for a  in old_alpha:
                self.alpha_arr[a[0]-xminval] = a[1]

        self.last_event = 'ButtonRelease-1'
        self.points[curr_ind] = (X,Y)

    def check_scalarval(self, value):
        """Check if a new  value for the selected dot lies in the range of
        its neighbors' values."""

        val_l = self.values[self.curr_ind - 1]
        val_r = self.values[self.curr_ind + 1]
        
        if val_r == self.xmaxval and self.curr_ind == self.shapes[-2]:
            val_r = val_r + 1
        if val_l == self.xminval and self.curr_ind == 1:
            val_l = val_l - 1
        if value > val_l and value < val_r:
            x_ind = self.values[self.curr_ind]-self.xminval
            if 'edge' in self.canvas.gettags("dot%d"%self.curr_ind) and \
               self.alpha_arr[x_ind] != 0:
                return 0
            else:
                return 1
        else : return 0

    def move_dot_value(self, value):
        curr_ind = self.curr_ind
        old_val = self.values[curr_ind]
        if value == old_val:
            return
        c = self.canvas
        tags = c.gettags('selected')
        if len(tags) == 0: #no dot is selected
            self.callbacks['entries'](val_x='0')
            return
        if tags[0] != 'dot': return        
        curr_tag = 'dot%d'%(curr_ind,)
        next_l = self.points[curr_ind-1]
        next_r = self.points[curr_ind+1]        
        lastx = self.points[curr_ind][0]
        x = (value-self.xminval)*self.sclx+self.left
        c.move(curr_tag, x-lastx, 0)
        y = self.points[curr_ind][1]
        c.coords('line%d'%(curr_ind-1,), next_l[0], next_l[1], x, y)
        c.coords('line%d'%(curr_ind,), x, y, next_r[0], next_r[1])
        if "edge" in tags:
            c.dtag(curr_tag, 'edge')
        self.lastx = x
        c.update_idletasks()
        self.points[curr_ind] = (x,y)
        self.values[curr_ind] = value
        scalars = [self.values[curr_ind-1], value,
                   self.values[curr_ind+1]]
        alphas = [self.alpha_arr[scalars[0]-self.xminval],
                  round((self.bott-y)/self.scly),
                  #self.alpha_arr[old_val-self.xminval],
                  self.alpha_arr[scalars[2]-self.xminval]]
        self.curr_color = list(self.color_arr[old_val-self.xminval])
        self.compute_rgba(scalars, alphas)
        self.alphas = alphas 
        self.last_event = 'ButtonRelease-1'
        

    def moveIsoVal_alpha(self, alpha):
        """Move selected Iso dot when alpha entry is changed"""
        c = self.canvas
        tags = c.gettags('selected')
        if len(tags) == 0: #no dot is selected
            return
        if tags[0] != 'isodot': return
        iso_ind = int(tags[1][6:])
        y = self.bott-self.scly*alpha
        rad = self.dotrad
        x = c.coords('isoline%d'%(iso_ind,))[0]
        c.coords('isoline%d'%(iso_ind,), x, self.bott, x, y)
        c.coords('isodot%d'%(iso_ind,), x-rad, y+rad, x+rad, y-rad)
        c.update_idletasks()
        val = self.isoVals[iso_ind]['val']
        self.isoVals[iso_ind]['alpha'] = alpha
        if self.continuous:
            self.callbacks['setalpha']([val,
                              Numeric.array([alpha]).astype(Numeric.Int16)])
        
    def mouse2Down(self, event):
        
        c= self.canvas
        points = self.points
        self.lastx = c.canvasx(event.x)
        self.lasty = c.canvasy(event.y)
        tag = c.gettags(Tkinter.CURRENT)[1]
        curr_ind = int(tag[3:])
        shapes = self.shapes
        
        # find all points of selected shape & store them in select_shape
        i=0
        for shape in shapes:
            if curr_ind <= shape: break
            i=i+1
        l_shape = shapes[i-1]
        r_shape = shapes[i+1]
        if (shape-l_shape) > 1:
            select_shape = range(l_shape,shape+1)
        elif (r_shape-shape) >1:
            select_shape = range(shape, r_shape+1)
        # find coords of adjoining points of the shapes next to selected one
        next_l = select_shape[0]-1
        next_r = select_shape[-1]+1
        self.next_l = points[next_l][0]
        self.next_r = points[next_r][0]
        first_point = select_shape[0]
        
        # find distance between event.x and first and last points of shape
        rad = self.dotrad
        self.diffr = points[select_shape[-1]][0]+rad-self.lastx
        self.diffl = self.lastx-(points[first_point][0]-rad)
        self.firstx = self.lastx
                
        # calculate alphas, colors
        first_scal = self.values[select_shape[0]]
        last_scal = self.values[select_shape[-1]]

        # get alphas(colors) for every entry(scalar) in LUT in interval between
        # first and last points of the selected shape
        xminval = self.xminval
        if not(self.last_event == 'ButtonRelease-2'\
           and select_shape == self.select_shape):
            #list of alphas of the shape
            self.alpha_list = list(self.alpha_arr[first_scal-xminval:
                                                 last_scal-xminval+1])
            #color arrays of the shape
##              self.rgb = (self.color_arr[first_scal-xminval:
##                                         last_scal+1-xminval]*1).astype('f')
            self.rgb = Numeric.array(self.color_arr[first_scal-xminval:
                                       last_scal+1-xminval])
            
                
        # memorize initial value(scalar) of the first point 
        self.old_valx1 = first_scal
        #c.itemconfig('selected', outline='')
        c.itemconfig('selected', outline=UNSELECTED_COLOR)
        ##  c.itemconfig('selected', fill=UNSELECTED_COLOR,
##                       outline=UNSELECTED_COLOR)
        c.dtag('selected')
        c.addtag_withtag('selected', Tkinter.CURRENT)
        c.itemconfig('selected', outline=SELECTED_COLOR)
        ##  c.itemconfig('selected', fill=SELECTED_COLOR, outline=SELECTED_COLOR)
        # get color of current poit
        curr_scal = self.values[curr_ind]-xminval
        curr_color = list(self.color_arr[curr_scal])
        curr_hsv = ToHSV(curr_color)
        #print "curr_hsv in mouse2Down:", curr_hsv
        self.select_shape = select_shape
        self.counter = 0
        self.curr_ind = curr_ind
        self.last_event = 'Button-2'
        
        #check if there is an Iso val in the interval
        self.isoval_inrange = []
        if len(self.isoVals):
            for iso_val in self.isoVals:
                v = iso_val['val']
                if v >= self.values[next_l] and v <= self.values[next_r]:
                    self.isoval_inrange.append(iso_val)
                    
        self.callbacks['entries'](val_x1=first_scal,
                            val_x2=first_scal+len(self.alpha_list)-1,
                            val_y = int(self.alpha_arr[curr_scal]),
                            val_r=round(curr_color[0],2),
                            val_g=round(curr_color[1],2),
                            val_b=round(curr_color[2],2))
        self.callbacks['rgbmap'](curr_hsv)

    def mouse2Move(self, event):
        self.last_event = 'B2-Motion'
        select_shape = self.select_shape
        c = self.canvas
        if 'edge' in c.gettags('dot%d'%(select_shape[-1],)) or \
           'edge' in c.gettags('dot%d'%(select_shape[0],)):
            return
        self.counter = self.counter+1
        x = c.canvasx(event.x)
        eps = 1e-6
        rad = self.dotrad
        if x+self.diffr >= self.next_r:
            if self.next_r == self.right: x=self.next_r-self.diffr+rad
            else: x=self.next_r-self.diffr-self.sclx+eps+rad
        if x-self.diffl <= self.next_l:
            if self.next_l == self.left: x=self.next_l+self.diffl-rad
            else: x=self.next_l+self.diffl+self.sclx+eps-rad
        dist = x-self.lastx
        # move current shape
        for i in select_shape:
            c.move('dot%d'%(i,), dist, 0)
            if i == select_shape[-1]:
                break
            else:
                c.move('line%d'%(i,), dist, 0)
        j=select_shape[0]
        # update side lines
        l_line = c.coords('line%d'%(j-1,))
        r_line = c.coords('line%d'%(i,))
        
        x1 = l_line[2]+dist
        x2 = r_line[0]+dist
        c.coords('line%d'%(j-1,), l_line[0], l_line[1],
                x1 , l_line[1])
        c.coords('line%d'%(i,), x2, r_line[1],
                 r_line[2], r_line[3])
        self.lastx = x
        # update val --- alpha values
        if self.counter%20 == 0:
            old_valx1 = self.old_valx1
            len_list = len(self.alpha_list)
            # get new value(scalar) of the first point of the shape 
            new_valx1 = int(round((x1-self.left)/self.sclx))+self.xminval
            delta = old_valx1-new_valx1
            end = new_valx1+len_list-1 
            xminval = self.xminval
            if delta < 0:
                r_val = self.values[select_shape[-1]+1]
                if end >= r_val:
                    if r_val == self.xmaxval:  end = r_val
                    else:  end = r_val-1
                    new_valx1 = end-len_list+1
                start_p = old_valx1-xminval
                mid_p = new_valx1-xminval
                stop_p = mid_p+len_list
                self.alpha_arr[start_p:mid_p] = 0
                self.alpha_arr[mid_p:stop_p] = self.alpha_list
            elif delta > 0:
                l_val = self.values[select_shape[0]-1]
                if new_valx1 <= l_val:
                    if l_val == xminval:  new_valx1 = xminval
                    else: new_valx1 = l_val+1
                    end = new_valx1+len_list-1
                    delta = old_valx1-new_valx1
                start_p = new_valx1-xminval
                mid_p = start_p+len_list
                stop_p = mid_p+delta
                self.alpha_arr[start_p:mid_p] = self.alpha_list
                self.alpha_arr[mid_p:stop_p] = 0
            elif delta == 0: return
            self.color_arr[new_valx1-xminval:end-xminval+1] = self.rgb[:][:]
            self.old_valx1 = new_valx1
            #if there is any Iso vals in the range:
            # memorize current rgba values , make changes in self.alpha_arr
            # and self.color_arr, do callbacks with new alphas and colors
            # then assign old values to self.alpha_arr and self.color_arr
            old_rgba = None
            if len(self.isoval_inrange):
                old_rgba = []
                for iso_val in self.isoval_inrange:
                    v = iso_val['val']
                    old_rgba.append((v, self.alpha_arr[v-xminval],
                                     self.color_arr[v-xminval].tolist()))
                    self.alpha_arr[v-xminval] = iso_val['alpha']
                    self.color_arr[v-xminval] = iso_val['rgb']
        
            if self.continuous:
                self.callbacks['setcolor']([new_valx1,
                             self.color_arr[new_valx1-xminval:end-xminval+1]])
                self.callbacks['setalpha']([start_p+xminval,
                                  self.alpha_arr[start_p:stop_p]])
            if old_rgba:
                for rgba in old_rgba:
                    v = rgba[0]
                    self.alpha_arr[v-xminval] = rgba[1]
                    self.color_arr[v-xminval] = rgba[2]
            self.callbacks['entries'](val_x1=new_valx1, val_x2=new_valx1+len_list-1)
            
            
    def mouse2Up(self, event):
        self.last_event = 'ButtonRelease-2'
        # update list of points
        dist = self.lastx - self.firstx
        sclx = self.sclx
        if dist == 0 : return
        if self.counter%20 == 0:
            diff = self.old_valx1 - self.values[self.select_shape[0]]
        else:
            diff = int(round(dist/sclx)) 
        for k in self.select_shape:
            self.points[k] = (self.points[k][0]+dist, self.points[k][1])
            self.values[k] = self.values[k]+diff
        # get alphas and colors
        left = self.left
        # get new value(scalar) of the first and last points of the shape 
        new_valx1 = self.values[self.select_shape[0]] #first point
        len_list = len(self.alpha_list)
        end = self.values[self.select_shape[-1]] #last point
        if self.counter%20 != 0:
            delta = self.old_valx1-new_valx1
            if delta < 0: #the shape moved to the right
                r_val = self.values[self.select_shape[-1]+1]
                if end >= r_val:
                    if r_val == self.xmaxval:
                        d = end - self.xmaxval
                        end = r_val
                    else:
                        d = end-r_val+1
                        end = r_val-1
                    new_valx1 = end-len_list+1
                    if d>0: 
                        self.values[self.select_shape[0]:
                                    self.select_shape[-1]+1] = \
                                    map(lambda val, d=d: val-d,
                                        self.values[self.select_shape[0]:
                                                    self.select_shape[-1]+1])
                start_p = self.old_valx1-self.xminval
                mid_p = new_valx1-self.xminval
                stop_p = mid_p+len_list
                self.alpha_arr[start_p:mid_p] = 0
                self.alpha_arr[mid_p:stop_p] = self.alpha_list
            elif delta > 0: #the shape moved to the left
                l_val = self.values[self.select_shape[0]-1]
                if new_valx1 <= l_val:
                    if l_val == self.xminval:
                        d = self.xminval - new_valx1
                        new_valx1 = self.xminval
                    else:
                        d = l_val - new_valx1-1
                        new_valx1 = l_val+1
                    end = new_valx1+len_list-1
                    delta = self.old_valx1-new_valx1
                    if d < 0:
                        self.values[self.select_shape[0]:
                                    self.select_shape[-1]+1] = \
                                    map(lambda val, d=d: val+d,
                                        self.values[self.select_shape[0]:
                                                    self.select_shape[-1]+1])
                start_p = new_valx1-self.xminval
                mid_p = start_p+len_list
                stop_p = mid_p+delta
                self.alpha_arr[start_p:mid_p] = self.alpha_list
                self.alpha_arr[mid_p:stop_p] = 0
            elif delta == 0: return
            self.color_arr[new_valx1-self.xminval:end-self.xminval+1] = \
                                                          self.rgb[:][:]
            
        # find scalar values of neighboring points of current shape
        next_l = self.values[self.select_shape[0]-1]
        next_r = self.values[self.select_shape[-1]+1]
        #print 'call to update_colors(%d,%d)'%(next_l-xminval,new_valx1-xminval)
        #update colors in the intervals outside of the shape
        
        self.update_colors(next_l-self.xminval, new_valx1-self.xminval)
        #print 'call to update_colors(%d,%d)'%(end-1-xminval, next_r-xminval)
        self.update_colors(end-self.xminval, next_r-self.xminval)
        
        #if there is any Iso vals in the range:
        # memorize current rgba values , make changes in self.alpha_arr
        # and self.color_arr, do callbacks with new alphas and colors
        # then assign old values to self.alpha_arr and self.color_arr
        old_rgba = None
        if len(self.isoval_inrange):
            old_rgba = []
            for iso_val in self.isoval_inrange:
                v = iso_val['val']
                old_rgba.append((v, self.alpha_arr[v-self.xminval],
                                 self.color_arr[v-self.xminval].tolist()))
                self.alpha_arr[v-self.xminval] = iso_val['alpha']
                self.color_arr[v-self.xminval] = iso_val['rgb']
        if self.counter%20 != 0:
            self.callbacks['entries'](val_x1=new_valx1, val_x2=new_valx1+len_list-1)
            if self.continuous:
                self.callbacks['setalpha']([start_p+self.xminval,
                               self.alpha_arr[start_p:stop_p]])
        if self.continuous:
            self.callbacks['setcolor']([next_l,
                 self.color_arr[next_l-self.xminval:next_r+1-self.xminval]])
        if old_rgba:
            for rgba in old_rgba:
                v = rgba[0]
                self.alpha_arr[v-self.xminval] = rgba[1]
                self.color_arr[v-self.xminval] = rgba[2]

    def apply_lut(self):
        """Called in non continuous mode (when the 'Apply LUT' button
                                       is pressed). """
        if len(self.isoVals):
            for iso_val in self.isoVals:
                v = iso_val['val']
                self.alpha_arr[v-self.xminval] = iso_val['alpha']
                self.color_arr[v-self.xminval] = iso_val['rgb']
        #self.callbacks['setalpha'] ([self.xminval, self.alpha_arr])
        #self.callbacks['setcolor'] ([self.xminval, self.color_arr])


    def addDot(self, event):
##          self.last_event = 'Double-Button-1'
        self.last_event = 'ButtonRelease-1' 
        c= self.canvas
        x = c.canvasx(event.x)
        y = c.canvasy(event.y)
        tag = c.gettags(Tkinter.CURRENT)[1]
        curr_ind = int(tag[4:])
        # if line between two dots which values differ less than 2 is
        # picked - do nothing 
        if self.points[curr_ind+1][0]-self.points[curr_ind][0]<self.sclx*2:
            return
        # find coords of selected line
        curr_coords = c.coords('current')
        # if side horizontal line is picked --- do nothing
        if curr_ind in self.shapes and curr_ind+1 in self.shapes:
            return
        if curr_coords[0] == curr_coords[2]:
            x = curr_coords[0]
        curr_val = int(round((x-self.left)/self.sclx))+self.xminval
        if curr_val == self.values[curr_ind+1]:
            curr_val = self.values[curr_ind+1]-1
        elif curr_val == self.values[curr_ind]:
            curr_val = self.values[curr_ind]+1
        #print "curr. scalar:", curr_val-self.xminval
        # else --- remove the line, update lines and dots tags, 
        # draw a dot and two side lines
        c.delete('current')
        #c.itemconfig('selected', outline='')
        c.itemconfig('selected', outline=UNSELECTED_COLOR)
        ##  c.itemconfig('selected', fill=UNSELECTED_COLOR,
##                       outline=UNSELECTED_COLOR)
        c.dtag('selected')
        self.update_linetags(curr_ind+1, 1)
        self.update_dottags(curr_ind+1, 1)

        rad = self.dotrad
        c.create_line(curr_coords[0], curr_coords[1], x, y, width = 2,
                      activefill='gray68',
                      fill='black', tags= ('line', 'line%d'%(curr_ind,)))
        c.create_line(x, y, curr_coords[2], curr_coords[3], width = 2,
                      activefill='gray68',
                      fill='black', tags= ('line', 'line%d'%(curr_ind+1,)))
##          c.create_oval(x-rad, y-rad, x+rad, y+rad, outline=SELECTED_COLOR,
##                        fill=SELECTED_COLOR,
##                        tags=('dot','dot%s'%str(curr_ind+1),'selected'))
        
        curr_color = list(self.color_arr[curr_val-self.xminval])
        #print "curr_color:", curr_color
        tkcolor = colorUtil.TkColor(curr_color)
        c.create_oval(x-rad, y-rad, x+rad, y+rad, outline=SELECTED_COLOR,
                      fill=tkcolor, activefill='gray68',
                      tags=('dot','dot%d'%(curr_ind+1,),'selected'))        
        self.line_count = self.line_count+1
        self.dot_count = self.dot_count+1
        c.lift('dot%d'%(curr_ind+2,))
        c.lift('dot%d'%(curr_ind,))
        # update list of points
        self.points.insert(curr_ind+1, (x,y))
        # update list of values(scalars)
        self.values.insert(curr_ind+1, curr_val)
        # update list of shapes
        i=0
        for shape in self.shapes:
            if curr_ind < shape:
                break
            i=i+1
        shapes_len = len(self.shapes)
        for j in range(i, shapes_len):
            self.shapes[j]=self.shapes[j]+1
        
        curr_hsv = ToHSV(curr_color)
        self.callbacks['entries'](val_x1=curr_val,
                             val_y=int(self.alpha_arr[curr_val-self.xminval]),
                             val_r=round(curr_color[0],2),
                             val_g=round(curr_color[1],2),
                             val_b=round(curr_color[2],2))
        self.callbacks['rgbmap'](curr_hsv)
        self.curr_ind = curr_ind + 1
        
    def deleteDot(self):
        curr_ind = self.curr_ind
        if curr_ind in self.shapes:
            return
        c = self.canvas
        c.delete('dot%d'%(curr_ind,))
        c.delete('line%d'%(curr_ind,))
        c.coords('line%d'%(curr_ind-1,), self.points[curr_ind-1][0],
                                           self.points[curr_ind-1][1],
                                           self.points[curr_ind+1][0],
                                           self.points[curr_ind+1][1])
        self.update_linetags(curr_ind+1, -1)
        self.update_dottags(curr_ind+1, -1)
        if curr_ind < self.line_count:
            c.lift('dot%d'%(curr_ind,))
        self.line_count = self.line_count - 1
        self.dot_count = self.dot_count - 1
        self.points.pop(curr_ind)
        self.values.pop(curr_ind)
        # update list of shapes
        i=0
        for shape in self.shapes:
            if curr_ind < shape:
                break
            i=i+1
        shapes_len = len(self.shapes)
        for j in range(i, shapes_len):
            self.shapes[j]=self.shapes[j]-1
        # set alphas 
        xminval = self.xminval
        d = self.scalars[0]-self.scalars[2]
        if d == 0: return
        K = (self.alphas[0]-self.alphas[2])/d
        C = self.alphas[0]-K*self.scalars[0]
        for s in range(self.scalars[0], self.scalars[2]):
            self.alpha_arr[s-xminval] = K*s+C
        self.update_colors(self.scalars[0]-xminval, self.scalars[2]-xminval)
        #if there is any Iso vals in the range:
        # memorize current rgba values , make changes in self.alpha_arr
        # and self.color_arr, do callbacks with new alphas and colors
        # then assign old values to self.alpha_arr and self.color_arr
        old_rgba = None
        if len(self.isoval_inrange):
            old_rgba = []
            for iso_val in self.isoval_inrange:
                v = iso_val['val']
                old_rgba.append((v, self.alpha_arr[v-xminval],
                                 self.color_arr[v-xminval].tolist()))
                self.alpha_arr[v-xminval] = iso_val['alpha']
                self.color_arr[v-xminval] = iso_val['rgb']
        if self.continuous:
            self.callbacks['setalpha']([self.scalars[0],
                  self.alpha_arr[self.scalars[0]-xminval:self.scalars[2]-1-xminval]])
            self.callbacks['setcolor']([self.scalars[0],
                   self.color_arr[self.scalars[0]-xminval:self.scalars[2]-xminval+1]])
        if old_rgba:
            for rgba in old_rgba:
                v = rgba[0]
                self.alpha_arr[v-xminval] = rgba[1]
                self.color_arr[v-xminval] = rgba[2]
        
        self.callbacks['entries'](val_y='0', val_x1='0', val_x2='',
                            val_r='1.0', val_g='1.0', val_b='1.0')

    def deleteShape(self):
        c = self.canvas
        curr_ind = self.curr_ind
        tags = ['dot%d'%(curr_ind,), 'line%d'%(curr_ind,),
                'line%d'%(curr_ind-1,), 'dot%d'%(curr_ind-1,),
                'dot%d'%(curr_ind+1,), 'line%d'%(curr_ind+1,)]
        for tag in tags:
            c.delete(tag)
        c.coords('line%d'%(curr_ind-2,), self.points[curr_ind-2][0],
                                           self.points[curr_ind-2][1],
                                           self.points[curr_ind+2][0],
                                           self.points[curr_ind+2][1])
        self.update_linetags(curr_ind+1, -3)
        self.update_dottags(curr_ind+1, -3)
        self.line_count = self.line_count - 3
        self.dot_count = self.dot_count - 3
        if curr_ind-1 < self.line_count:
            c.lift('dot%d'%(curr_ind-1,))
            
        self.points.pop(curr_ind)
        self.points.pop(curr_ind-1)
        self.points.pop(curr_ind-1)
        
        self.values.pop(curr_ind)
        self.values.pop(curr_ind-1)
        self.values.pop(curr_ind-1)
        
        # update list of shapes
        if self.update_shapes(curr_ind-1, -3) == 0: return
        
        # update alphas and colors
        scalars = self.scalars
        xminval = self.xminval
        self.alpha_arr[scalars[0]-xminval:scalars[2]+1-xminval] = 0
        #values of the two neighbouring dots of the deleted shape:
        l_scalar = self.values[curr_ind-2]
        r_scalar = self.values[curr_ind-1]
        self.update_colors(l_scalar-xminval, r_scalar-xminval)
            
        #if there is any Iso vals in the range(l_scalar, r_scalar):
        # memorize current rgba values , make changes in self.alpha_arr
        # and self.color_arr, do callbacks with new alphas and colors
        # then assign old values to self.alpha_arr and self.color_arr
        old_rgba = None
        if len(self.isoVals):
            for iso_val in self.isoVals:
                v = iso_val['val']
                if v >= l_scalar and v <= r_scalar:
                    old_rgba.append((v, self.alpha_arr[v-xminval],
                                 self.color_arr[v-xminval].tolist()))
                    self.alpha_arr[v-xminval] = iso_val['alpha']
                    self.color_arr[v-xminval] = iso_val['rgb']
        if self.continuous:
            self.callbacks['setalpha']([scalars[0],
                     self.alpha_arr[scalars[0]-xminval:scalars[2]+1-xminval]])
            self.callbacks['setcolor']([l_scalar,
                      self.color_arr[l_scalar-xminval: r_scalar-xminval+1]] )
        if old_rgba:
            for rgba in old_rgba:
                v = rgba[0]
                self.alpha_arr[v-xminval] = rgba[1]
                self.color_arr[v-xminval] = rgba[2]
        self.callbacks['entries'](val_y='0', val_x1='0', val_x2='',
                            val_r='1.0', val_g='1.0', val_b='1.0')



    def pick(self, event):
        current = self.canvas.gettags(Tkinter.CURRENT)
        print 'current' , current

    def set_colors(self, rgb, scalar1, scalar2, scalar3):
        """find range of data values(scalars) where colors are to be changed"""
        xminval = self.xminval
        if self.last_event == 'ButtonRelease-1' or \
           self.last_event == 'ButtonRelease-2' or \
           self.last_event == 'B1-Motion':
##              print 'scalars(in set_colors)=[%d,%d,%d]'%(scalar1,scalar2,scalar3)
            red1, green1, blue1 = self.color_arr[scalar1]
            red2, green2, blue2 = rgb
            red3, green3, blue3 = self.color_arr[scalar3]
            d = scalar1-scalar2
            if d != 0:
                Kr = (red1-red2)/d
                Kg = (green1-green2)/d
                Kb = (blue1-blue2)/d
                Cr = red1-Kr*scalar1
                Cg = green1-Kg *scalar1
                Cb = blue1-Kb*scalar1
                s = Numeric.arange(scalar1, scalar2)
                self.color_arr[scalar1:scalar2] = \
                   Numeric.transpose(Numeric.array((Kr*s+Cr, Kg*s+Cg,
                                                    Kb*s+Cb)).astype('f'))
            d = scalar2-scalar3
            if d != 0:
                Kr = (red2-red3)/d
                Kg = (green2-green3)/d
                Kb = (blue2-blue3)/d
                Cr = red2-Kr*scalar2
                Cg = green2-Kg *scalar2
                Cb = blue2-Kb*scalar2
                s = Numeric.arange(scalar2, scalar3)
                self.color_arr[scalar2:scalar3] = \
                   Numeric.transpose(Numeric.array((Kr*s+Cr, Kg*s+Cg,
                                                    Kb*s+Cb)).astype('f'))

            self.color_arr[scalar2] = rgb
            curr_color = self.color_arr[scalar2]
            tkcolor = colorUtil.TkColor(curr_color)
            l_color = colorUtil.TkColor(self.color_arr[scalar1])
            r_color = colorUtil.TkColor(self.color_arr[scalar3])
            self.canvas.itemconfig('dot%d'%(self.curr_ind+1,), fill=r_color)
            self.canvas.itemconfig('dot%d'%(self.curr_ind-1,), fill=l_color)
            self.canvas.itemconfig('selected', fill=tkcolor)
        else : return

    def set_hue_satur_value(self, hsv):
        """find range of data values(scalars) where hue and saturation
        are to be changed. (Colormap object callback)."""
        
        tags = self.canvas.gettags('selected')
        if len(tags) == 0: # nothing is selected
            return
        if tags[0] == 'isodot':
            rgb = ToRGB(hsv)
            iso_ind = int(tags[1][6:])
            self.setIsoColor(iso_ind, rgb)
            return
        xminval = self.xminval
        if self.last_event == 'ButtonRelease-1' or \
           self.last_event == 'ButtonRelease-2':
            curr_ind = self.curr_ind
            scalar1 = self.values[curr_ind-1]-xminval
            scalar2 = self.values[curr_ind]-xminval
            scalar3 = self.values[curr_ind+1]-xminval
            self.last_event = 'ButtonRelease-1' 
##              print 'scalars(in set_colors)=[%d,%d,%d]'%(scalar1,scalar2,scalar3)
        else : return
        rgb = ToRGB(hsv)
        #print "hsv in set_hue_satur.. ", hsv, "rgb", rgb
        #print "scalar1: %d, scalar2: %d, scalar3: %d"%(scalar1,scalar2,scalar3)
        self.set_colors(rgb, scalar1, scalar2, scalar3)
        #if there is any Iso vals in the range:
        # memorize current rgb values , make changes in self.color_arr,
        #do callback with new colors, then assign old values to self.color_arr
        old_rgb = None
        if len(self.isoval_inrange):
            old_rgb = []
            for iso_val in self.isoval_inrange:
                v = iso_val['val']
                old_rgb.append((v, self.color_arr[v-xminval].tolist()))
                self.color_arr[v-xminval] = iso_val['rgb']
        if self.continuous:
            self.callbacks['setcolor']([scalar1+xminval,
                              self.color_arr[scalar1:scalar3+1]])
        if old_rgb:
            for c in old_rgb:
                self.color_arr[c[0]-xminval] = c[1]
                
        self.callbacks['entries'](val_r=round(rgb[0],2),
                             val_g=round(rgb[1],2),
                             val_b=round(rgb[2],2))

    def set_red_green_blue(self, rgb):
        """Sets new rgb (typed in entries) to selected value(dot),
        and updates(interpolates) colors between this dot and two
        neighboring dots"""
        
        tags = self.canvas.gettags('selected')
        if len(tags) == 0: # nothing is selected
            return
        if tags[0] != 'dot': return
        xminval = self.xminval
        if self.last_event == 'ButtonRelease-1' or \
           self.last_event == 'ButtonRelease-2':
            curr_ind = self.curr_ind
            scalar1 = self.values[curr_ind-1]-xminval
            scalar2 = self.values[curr_ind]-xminval
            scalar3 = self.values[curr_ind+1]-xminval
            self.last_event = 'ButtonRelease-1' 
##              print 'scalars(in set_colors)=[%d,%d,%d]'%(scalar1,scalar2,scalar3)
        else : return
        self.set_colors(rgb, scalar1, scalar2, scalar3)

        #if there is any Iso vals in the range:
        # memorize current rgb values , make changes in self.color_arr,
        #do callback with new colors, then assign old values to self.color_arr
        old_rgb = None
        if len(self.isoval_inrange):
            old_rgb = []
            for iso_val in self.isoval_inrange:
                v = iso_val['val']
                old_rgb.append((v, self.color_arr[v-xminval].tolist()))
                self.color_arr[v-xminval] = iso_val['rgb']
        if self.continuous:
            self.callbacks['setcolor']([scalar1+xminval,
                              self.color_arr[scalar1:scalar3+1]])
        if old_rgb:
            for c in old_rgb:
                self.color_arr[c[0]-xminval] = c[1]
        

    def update_colors(self,scalar1, scalar2):
        red1 = self.color_arr[scalar1][0]
        red2 = self.color_arr[scalar2][0]
        green1 = self.color_arr[scalar1][1]
        green2 = self.color_arr[scalar2][1]
        blue1 = self.color_arr[scalar1][2]
        blue2 = self.color_arr[scalar2][2]
        d = scalar2-scalar1
        if d == 0: return [scalar1, scalar2, Numeric.zeros((0,3),'f')]
        r = (red2-red1)/d
        g = (green2-green1)/d
        b = (blue2-blue1)/d
        s = Numeric.arange(d)
        self.color_arr[scalar1:scalar2] = \
          Numeric.transpose(Numeric.array((r*s+red1,
                                           g*s+green1, b*s+blue1)).astype('f'))

    def set_ISOcolor(self, rgb):
        """change the color of selected Iso value when new rgb values
        are entered in ISO entries """
        
        tags = self.canvas.gettags('selected')
        if len(tags) == 0: # nothing is selected
            return
        if tags[0] != 'isodot': return
        iso_ind = int(tags[1][6:])
        self.setIsoColor(iso_ind, rgb)


    def setIsoColor(self, iso_ind, rgb):
        tkcolor = colorUtil.TkColor(rgb)
        self.canvas.itemconfig('isodot%d'%(iso_ind,), fill=tkcolor)
        val = self.isoVals[iso_ind]['val']
        iso_rgb = (round(rgb[0],2), round(rgb[1],2), round(rgb[2],2))
        self.isoVals[iso_ind]['rgb'] = iso_rgb
        if self.continuous:
            self.callbacks['setcolor']([val, Numeric.array([rgb],'f')])
        self.callbacks['entries'](iso_rgb = iso_rgb)
        
    
    def bind_tags(self):
        self.canvas.tag_bind('box', '<Shift-Button-1>', self.drawShape)
        self.canvas.tag_bind('line', '<Double-Button-1>', self.addDot)
        self.canvas.tag_bind('dot', '<Button-2>', self.mouse2Down)
        self.canvas.tag_bind('dot','<B1-Motion>', self.mouse1Move)
        self.canvas.tag_bind('dot','<B2-Motion>', self.mouse2Move)
        self.canvas.tag_bind('dot', '<ButtonRelease-2>', self.mouse2Up)

    def unbind_tags(self):
        self.canvas.tag_unbind('box', '<Shift-Button-1>')
        self.canvas.tag_unbind('line', '<Double-Button-1>')
        self.canvas.tag_unbind('dot', '<Button-2>')
        self.canvas.tag_unbind('dot','<B1-Motion>')
        self.canvas.tag_unbind('dot','<B2-Motion>')
        self.canvas.tag_unbind('dot', '<ButtonRelease-2>')

    def calculate_points(self, values, alphas):
        """Calculates points from given alpha values and scalar values at LUT points"""
        self.points = [(self.left, self.bott)] 
        for i in range(1, len(values)-1):
            self.points.append(((values[i]-self.xminval)*self.sclx+self.left,
                                self.bott-self.scly*alphas[i]))
        self.points.append((self.right, self.bott))
        #print "points: ", self.points


    def calculate_alphas(self, values = None, alphas_points=None):
        """Calculate alpha array and values of the points from self.points list"""
        scly = self.scly
        bott = self.bott
        sclx = self.sclx
        left = self.left
        xminval = self.xminval
        if not values:
            values = []
            for point in self.points:
                s = int(round((point[0]-left)/self.sclx))+xminval
                if s > self.xmaxval: s = self.xmaxval
                values.append(s)
            
            values[0] = xminval
            values[-1] = self.xmaxval
        if not alphas_points:
            alphas_points = []
            for point in self.points:
                alphas_points.append(round((bott-point[1])/scly))
##          print 'values=', values
        k=0
        for j in range(len(alphas_points)-1):
            d = (values[j]-values[j+1])
            if d == 0:
##                  print 'values[%s],[%s]:%s,%s'%(j,j+1,values[j],values[j+1]), 'd=0'
                continue
            K = (alphas_points[j]-alphas_points[j+1])/d
            C = alphas_points[j]-K*values[j]
            for s in range(values[j], values[j+1]):
                self.alpha_arr[k] = K*s+C
                k=k+1
##          print 'k=',k, 's=',s, 'j=',j
        self.alpha_arr[k] = K*values[j+1]+C
        self.values = values

    def set_dotsize(self, rad, type = "dot"):
        """Set a new radius of the widget's dot (type = 'dot') or
        isodot (type = 'isodot')."""
        c = self.canvas
        if type == "dot":
            if rad == self.dotrad: return
            d = self.dotrad - rad
            self.dotrad = rad
            points=self.points
            num_dots = self.shapes[-1]
            for ind in range(1, num_dots):
                x0 = points[ind][0]-rad
                x1 = points[ind][0]+rad
                y0 = points[ind][1]-rad
                y1 = points[ind][1]+rad
                c.coords("dot%d"%ind , x0, y0, x1, y1)
            
        elif type == "isodot":
            if rad == self.isodotrad: return
            d = self.isodotrad - rad
            self.isodotrad = rad
            num_isodots = len(self.isoVals)
            for ind in range(num_isodots):
                coords = c.coords("isoline%d"%ind)[2:]
                x0 = coords[0]-rad
                x1=coords[0]+rad
                y0=coords[1]-rad
                y1=coords[1]+rad
                c.coords("isodot%d"%ind , x0, y0, x1, y1)
                

            
class Colormap(CallbackFunctions):
    def __init__(self, master, file=None):
        self.callbacks = []
        self.img = Tkinter.PhotoImage(file=file, master=master)
        width = self.img.width()
        height = self.img.height()
        self.axis_len = max(width, height)/2
        self.X0 = width/2
        self.Y0 = height/2
        self.lastx = self.X0
        self.lasty = self.Y0
        self.hsv =[0.0, 0.0, 0.0]
        font = '-*-Helvetica-Bold-R-Normal-*-*-160-*-*-*-*-*-*'
        Tkinter.Label(master, text='Hue/Saturation',
              font=font,width=20).grid(row=0,column=0,sticky='NW')
        Tkinter.Label(master, text='Value',
              font=font).grid(row=0,column=1,sticky='NW')
        self.canvas = Tkinter.Canvas(master, width=width,
                             height=height)
        self.canvas.grid(row=1, column=0, sticky='W')        
        
        self.canvas.create_image(0,0, anchor='nw',image=self.img)
        rad=4
        self.canvas.create_oval(self.X0-rad, self.Y0-rad, self.X0+rad,
                                self.Y0+rad, fill='white', outline='black',
                                activefill='gray68',
                                tags='oval')

        Tkinter.Widget.bind(self.canvas, "<Button-1>", self.mouseDown)
        Tkinter.Widget.bind(self.canvas,'<Button1-Motion>', self.mouseMove)
        Tkinter.Widget.bind(self.canvas,'<Button1-ButtonRelease>', self.mouseUp)
        
        self.scale = Tkinter.Scale(master, orient=Tkinter.VERTICAL,
                                   length=height-5,
                                   from_=0.0, to=1.0, tickinterval=0,
                                   font=font, resolution=0.05,
                                   command=lambda v,s=self: s.get_value(v))
        self.scale.grid(row=1, column=1, sticky='W')
        self.scale.set(1.0)
        self.updateval = 1

    def get_value(self, val):
        hsv = [self.hsv[0], self.hsv[1], float(val)]
        if self.updateval:
            self.callbacks[0](hsv)
        else: self.updateval = 1
        self.hsv = hsv

    def mouseDown(self, event):
	x = self.canvas.canvasx(event.x)
	y = self.canvas.canvasy(event.y)
        self.canvas.delete('oval')
        rad = 4
        self.canvas.create_oval(x-rad, y-rad, x+rad, y+rad,
                                activefill='gray68',
                                fill='white', outline='black', tags='oval')
        self.lastx = x
        self.lasty = y
        
    def get_hsv(self, X, Y):
        if X!= 0 or Y!=0:
            hue = math.acos(X/math.sqrt(X*X+Y*Y))
            if Y < 0:
                hue = 2*_PI-hue
        else: hue = 0
        hue = hue/(2*_PI)
        dist = math.sqrt(X*X+Y*Y)
        if dist > self.axis_len: dist = self.axis_len
        sat = dist/self.axis_len
        val = self.scale.get()
        if hue > 1.0: hue = 1.0
        elif hue < 0.0 : hue = 0.0
        if sat > 1.0: sat = 1.0
        elif sat < 0.0 : sat = 0.0
        if val > 1.0: val = 1.0
        elif val < 0.0 : val = 0.0
        self.hsv = [hue, sat, val]
        return [hue, sat, val]

    def mouseMove(self, event):
        c = self.canvas
        x = c.canvasx(event.x)
	y = c.canvasy(event.y)
        if y < 0: y = 0
        elif y > self.Y0*2: y=self.Y0*2
        if x < 0: x = 0
        elif x > self.X0*2: x = self.X0 * 2
        
        c.move('oval', x-self.lastx, y-self.lasty)
        self.lastx = x
        self.lasty = y

    def mouseUp(self, event):
        X = self.lastx-self.X0
        Y = self.Y0-self.lasty
        hsv = self.get_hsv(X,Y)
        self.callbacks[0](hsv)

    def show_color(self, hsv):
        h = hsv[0]*2*_PI
        s = hsv[1]*self.axis_len
        v = round(hsv[2],2)
        x = math.cos(h)*s+self.X0
        y = self.Y0-math.sin(h)*s
        self.canvas.move('oval', x-self.lastx, y-self.lasty)
        self.updateval = 0
        self.scale.set(v)
        self.lastx = x
        self.lasty = y
        self.hsv = hsv


if __name__ == '__main__':
    def MyCallback1(value):
        print 'callback with alpha argument'
        #print 'alpha = ', value

    def MyCallback2(value):
        #print 'colors =', value
        print 'callback with color argument'
            
    root = Tkinter.Tk()
    root.title('Transfer Function')
    #t = TableManager(root,xmaxval=50)
    #t = TableManager(root,xmaxval=4095)
    t = TableManager(master=root, ymaxval=255, xminval = 0, xmaxval=255,
                     alphaCallback=MyCallback1,
                     colorCallback=MyCallback2)
#    mainloop()
