#########################################################################
#
# Date: May 2001 Author: Daniel Stoffler
#
#    stoffler@scripps.edu
#
# Copyright: Daniel Stoffler, TSRI
#
#########################################################################


import Tkinter, Pmw
from mglutil.gui.InputForm.Tk.gui import InputFormDescr, InputForm
from mglutil.util.misc import ensureFontCase

class OptionsPanel:

    """ This class builds the options panel used in various GUIs.
master is the widget for which the panel is created.
title is the title of the 
    """

    def __init__(self, master=None, title=None):
        self.master = master
        self.root = None

        self.title = title
        if self.title is None:
            self.title = 'Options Panel'

        if self.master is None: # this is to test this widget without a widget
            self.master = Tkinter.Tk()
            self.master.continuous = 1

            def configure(**kw):
                pass

            self.master.configure = configure
            self.master.min = None        # min value allowed for val
            self.master.max = None        # max value allowed for val
            self.master.increment = None  # discrete value allowed for val
            self.master.incrementOld = 0.0
            self.master.oneTurn = 1.0
            self.master.value = 0.0
            self.master.type = float
            self.master.precision = 3
            self.master.showLabel = 0
            self.master.labelFormat = "%d"
            self.master.minOld = 0
            self.master.maxOld = 0
            self.master.lockContinuous = 0
            self.master.lockMin = 0
            self.master.lockMax = 0
            self.master.lockBMin = 0
            self.master.lockBMax = 0
            self.master.lockIncrement = 1
            self.master.lockBIncrement = 0
            self.master.lockValue = 0
            self.master.lockOneTurn = 0
            self.master.lockShowLabel = 0
            self.master.lockType = 0
            self.master.lockPrecision = 0
            self.master.opPanel = self
            #self.master.withdraw()
            
        self.flag = 0     # this flag toggles display/undisplay
                          # used in the various GUIs
                          
        self.type = float # can be float or int
        
        self.toggleMin  =  Tkinter.IntVar()
        self.minInput   =  Tkinter.StringVar()
        self.min_entry  =  None

        self.toggleMax  =  Tkinter.IntVar()
        self.maxInput   =  Tkinter.StringVar()
        self.max_entry  =  None

        self.toggleIncr =  Tkinter.IntVar()
        self.incrInput  =  Tkinter.StringVar()
        self.incr_entry =  None

        self.valInput   =  Tkinter.StringVar()
        self.sensInput  =  Tkinter.StringVar()

        self.toggleMin.set(0)
        self.toggleMax.set(0)
        self.toggleIncr.set(0)
       
        self.idf = InputFormDescr(title=self.title)
        
        self.idf.append({'widgetType':Tkinter.Label,
                         'wcfg':{'text': '\n'},
                         'gridcfg':{'columnspan':3, 'row':0, 'column':0},
                         })


        self.idf.append({'name':'togCont',
                         'widgetType':Pmw.OptionMenu,
                         'wcfg':{'labelpos':'w',
                                 'label_text':'Continuous    ',
                                 'menubutton_width':3,
                                 'items':('on', 'off'),
                                 'command': self.toggleCont_cb},
                         'gridcfg':{'sticky':'we',
                                    'columnspan':2, 'row':1, 'column':0},
                         })


        self.idf.append({'name':'togMin',
                         'widgetType':Tkinter.Checkbutton,
                         'wcfg':{'text':'Minimum',
                                 'variable':self.toggleMin,
                                 'command':self.toggleMin_cb},
                         'gridcfg':{'sticky':'w',
                                    'columnspan':2, 'row':2, 'column':0}
                         })


        self.idf.append({'name':'inpMin',
                         'widgetType':Tkinter.Entry,
                         'defaultValue':'0.0',
                         'wcfg':{'font':(ensureFontCase('Courier'),10),
                                 'width':5, 'textvariable':self.minInput,
                                 'command':self.inputMin_cb,
                                 'eventType':'<Return>'},
                         'gridcfg':{'sticky':'we',
                                    'columnspan':2, 'row':2, 'column':1 }
                         })


        self.idf.append({'name':'togMax',
                         'widgetType':Tkinter.Checkbutton,
                         'wcfg':{'text':'Maximum',
                                 'variable':self.toggleMax,
                                 'command':self.toggleMax_cb,
                                 },
                         'gridcfg':{'sticky':'w',
                                    'columnspan':2, 'row':3, 'column':0}
                         })


        self.idf.append({'name':'inpMax',
                         'widgetType':Tkinter.Entry,
                         'defaultValue':'0.0',
                         'wcfg':{'font':(ensureFontCase('Courier'),10),
                                 'width':5,'textvariable':self.maxInput,
                                 'command':self.inputMax_cb,},
                         'gridcfg':{'sticky':'we',
                                    'columnspan':2, 'row':3, 'column':1 }
                         })

        
        self.idf.append({'name':'togIncr',
                         'widgetType':Tkinter.Checkbutton,
                         'wcfg':{'text':'Increment',
                                 'variable':self.toggleIncr,
                                 'command':self.toggleIncr_cb},
                         'gridcfg':{'sticky':'w',
                                    'columnspan':2, 'row':4, 'column':0}
                         })


        self.idf.append({'name':'inpIncr',
                         'widgetType':Tkinter.Entry,
                         'defaultValue':'0.1',
                         'wcfg':{'font':(ensureFontCase('Courier'),10),
                                 'width':5,'textvariable':self.incrInput,
                                 'command':self.inputIncr_cb,
                                 },
                         'gridcfg':{'sticky':'we',
                                    'columnspan':2, 'row':4, 'column':1 }
                         })


        self.idf.append({'widgetType':Tkinter.Label,
                         'wcfg':{'text': 'Value'},
                         'gridcfg':{'sticky':'w',
                                    'columnspan':2, 'row':5, 'column':0},
                         })


        self.idf.append({'name':'inpVal',
                         'widgetType':Tkinter.Entry,
                         'defaultValue':'0.0',
                         'wcfg':{'font':(ensureFontCase('Courier'),10),
                                 'width':5,'textvariable':self.valInput,
                                 'command':self.inputVal_cb,
                                 },
                         'gridcfg':{'sticky':'we',
                                    'columnspan':2, 'row':5, 'column':1 }
                         })


        self.idf.append({'widgetType':Tkinter.Label,
                         'wcfg':{'text': 'Sensitivity'},
                         'gridcfg':{'sticky':'w',
                                    'columnspan':2, 'row':6, 'column':0},
                         })


        self.idf.append({'name':'inpSens',
                         'widgetType':Tkinter.Entry,
                         'defaultValue':self.master.oneTurn,
                         'wcfg':{'font':(ensureFontCase('Courier'),10),
                                 'width':5,'textvariable':self.sensInput,
                                 'command':self.inputSens_cb,
                                 },
                         'gridcfg':{'sticky':'we',
                                    'columnspan':2, 'row':6, 'column':1 }
                         })


        self.idf.append({'name':'togLabel',
                         'widgetType':Pmw.OptionMenu,
                         'wcfg':{'labelpos':'w',
                                 'label_text':'Show label',
                                 'menubutton_width':5,
                                 'items':('never', 'always', 'move'),
                                 'command': self.toggleLabel_cb},
                         'gridcfg':{'sticky':'we',
                                    'columnspan':2},
                         })


        self.idf.append({'name':'togIntFloat',
                         'widgetType':Pmw.OptionMenu,
                         'wcfg':{'labelpos':'w',
                                 'label_text':'Type',
                                 'menubutton_width':3,
                                 'items':('float', 'int'),
                                 'command': self.toggleIntFloat_cb},
                         'gridcfg':{'sticky':'we',
                                    'columnspan':2},
                         })


        self.idf.append({'name':'selPrec',
                         'widgetType':Pmw.OptionMenu,
                         'wcfg':{'labelpos':'w',
                                 'label_text':'Precision',
                                 'menubutton_width':3,
                                 'items':("1", '2', '3', '4', '5', '6', '7', '8',
                                          '9', '10' ),
                                 'command': self.selPrec_cb},
                         'gridcfg':{'sticky':'we',
                                    'columnspan':2},
                         })


#########Buttons##########

        self.idf.append({'name':'OKButton',
                         'widgetType':Tkinter.Button,
                         'text':'  OK  ',
                         'wcfg':{'bd':3},
                         'gridcfg':{'sticky':'we', 'row':10,'column':0},
                         'command': self.OK_cb})
        
	self.idf.append({'name':'ApplyButton',
                         'widgetType':Tkinter.Button,
                         'text':'Apply',
                         'wcfg':{'bd':3},
                         'gridcfg':{'sticky':'we', 'row':10,'column':1},
                         'command': self.Apply_cb})

	self.idf.append({'name':'CancelButton',
                         'widgetType':Tkinter.Button,
                         'text':'Dismiss',
                         'wcfg':{'bd':3},
                         'gridcfg':{'sticky':'we', 'row':11, 'column':0,
                                    'columnspan':2},
                         'command': self.Dismiss_cb})



    def toggleCont_cb(self, val):
        if val == 'off':
            self.master.configure(continuous=0)
        else:
            self.master.configure(continuous=1)


    def toggleMin_cb(self):
        if self.toggleMin.get() != 1:
            self.master.min = None
            self.min_entry.configure(state='disabled', fg='gray40')
        else:
            self.min_entry.configure(state='normal', fg='gray0')
            self.inputMin_cb()


    def inputMin_cb(self, event=None):
        val = self.minInput.get()
        
        if len(val) == 0: val = self.master.min
        try:
            val = float(val)
            self.master.configure(min=val)
        except ValueError:
            # put back original value if someone types garbage
            self.minInput.set(self.master.labelFormat%self.master.minOld)


    def toggleMax_cb(self):
        if self.toggleMax.get() != 1:
            self.master.max = None
            self.max_entry.configure(state='disabled', fg='gray40')
        else:
            self.max_entry.configure(state='normal', fg='gray0')
            self.inputMax_cb()

            
    def inputMax_cb(self, event=None):
        val = self.maxInput.get()
        if len(val) == 0: val = self.master.max
        try:
            val = float(val)
            self.master.configure(max=val)
        except ValueError:
            # put back original value if someone types garbage
            self.maxInput.set(self.master.labelFormat%self.master.maxOld)


    def toggleIncr_cb(self):
        if self.toggleIncr.get() != 1:
            self.master.increment = None#self.master.type(0)
            self.incr_entry.configure(state='disabled', fg='gray40')
            self.toggleIncr.set(0)
        else:
            self.incr_entry.configure(state='normal', fg='gray0')
            self.master.offsetValue = self.master.value
            self.inputIncr_cb()


    def inputIncr_cb(self, event=None):
        val = self.incrInput.get()
        if len(val) == 0: val = self.master.increment
        try:
            val = float(val)
            self.master.configure(increment=val)
        except ValueError:
            # put back original value if someone types garbage
            self.incrInput.set(self.master.labelFormat%self.master.incrementOld)


    def inputVal_cb(self, event=None):
        val = self.valInput.get()
        if len(val) == 0: val = self.master.value
        try:
            val = float(val)
            self.master.set(val, force=1) # force execution, even if mode
                                          # 'continuous' is OFF
        except ValueError:
            self.valInput.set(self.master.labelFormat%self.master.value)


    def inputSens_cb(self, event=None):
        val = self.sensInput.get()
        if len(val) == 0: val = self.master.oneTurn
        try:
            val = float(val)
            self.master.configure(oneTurn=val)
        except ValueError:
            # put back original value if someone types garbage
            self.sensInput.set(self.master.labelFormat%self.master.oneTurn)


    def toggleLabel_cb(self, val):
        if val == 'never':
            self.master.configure(showLabel=0)
        if val == 'always':
            self.master.configure(showLabel=1)
        if val == 'move':
            self.master.configure(showLabel=2)


    def toggleIntFloat_cb(self, val):
        if val == 'float':
            type=float
        else:
            type=int

        self.master.configure(type=type)
        self.showHidePrec(val)
        

    def selPrec_cb(self, val):
        val = int(val)
        self.master.configure(precision=val)
        

    def showHidePrec(self, val):
        if val == 'int':
            self.idf.entryByName['selPrec']['widget'].grid_forget()
        else:
            apply( self.idf.entryByName['selPrec']['widget'].grid,(),
                   self.idf.entryByName['selPrec']['gridcfg'])


    def setTitle(self, title):
        self.title = title
        self.idf.title = title
        if hasattr(self, 'optionsForm'):
            self.optionsForm.root.title(self.title)
        

    def OK_cb(self, event=None):
        self.Apply_cb()
        self.Dismiss_cb()


    def Apply_cb(self, event=None):
        # read all the Tkinter.Entry widgets and configure the widget
        # note: combobox widgets get applied automatically upon user selection

        if self.master.min is not None:
            self.inputMin_cb()
        if self.master.max is not None:
            self.inputMax_cb()
        if self.master.increment is not None:
            self.inputIncr_cb()
        self.inputVal_cb()
        self.inputSens_cb()
        

    def Dismiss_cb(self, event=None):
        self.flag = 0
	self.optionsForm.withdraw()


    def displayPanel(self, create):
        self.flag = 1
        if create == 0:
            self.optionsForm.deiconify()
        else:
            self.optionsForm = InputForm(self.master, self.root,
                                         descr = self.idf,
                                         modal = 0, blocking = 0)

            self.cont_entry  = self.idf.entryByName['togCont']['widget']
            self.min_entry   = self.idf.entryByName['inpMin']['widget']
            self.bmin_entry  = self.idf.entryByName['togMin']['widget']
            self.max_entry   = self.idf.entryByName['inpMax']['widget']
            self.bmax_entry  = self.idf.entryByName['togMax']['widget']
            self.incr_entry  = self.idf.entryByName['inpIncr']['widget']
            self.bincr_entry = self.idf.entryByName['togIncr']['widget']
            self.val_entry   = self.idf.entryByName['inpVal']['widget']
            self.sens_entry  = self.idf.entryByName['inpSens']['widget']
            self.lab_entry   = self.idf.entryByName['togLabel']['widget']
            self.if_entry    = self.idf.entryByName['togIntFloat']['widget']
            self.prec_entry  = self.idf.entryByName['selPrec']['widget']

            if self.master.min is not None:
                val = self.master.type(self.master.min)
                self.minInput.set(self.master.labelFormat%val)
                self.toggleMin.set(1)
                self.min_entry.configure(state='normal', fg='gray0')
            else:
                self.min_entry.configure(state='disabled', fg='gray40')
            self.bmin_entry.configure(state='normal', fg='gray0')

            if self.master.max is not None:
                val = self.master.type(self.master.max)
                self.maxInput.set(self.master.labelFormat%val)
                self.toggleMax.set(1)
                self.max_entry.configure(state='normal', fg='gray0')
            self.bmax_entry.configure(state='normal', fg='gray0')

            if self.master.increment is not None:
                val = self.master.type(self.master.increment)
                self.incrInput.set(self.master.labelFormat%val)
                if self.master.increment != 0:
                    self.toggleIncr.set(1)
                self.incr_entry.configure(state='normal', fg='gray0')
            self.bincr_entry.configure(state='normal', fg='gray0')

                               
            menus = (self.cont_entry, self.lab_entry, self.if_entry,
                     self.prec_entry)
            Pmw.alignlabels(menus)

            if self.master.continuous == None or self.master.continuous == 0:\
               i=1
            else: i=0
            self.cont_entry.setitems(
                items=('on','off'),index=i)

            if self.master.type == int:
                self.if_entry.setitems(
                    items=('float', 'int'), index = 1)

            prc = int(self.master.precision) - 1
            if prc > 9: prc = 9
            self.prec_entry.setitems(
                items=('1', '2', '3', '4', '5', '6', '7', '8', '9', '10'),
                index = prc)

            self.lab_entry.setitems(
                items=('never', 'always', 'move'),
                index=self.master.showLabel)


            if self.if_entry.getcurselection() == 'int':
                self.showHidePrec('int')

            self.updateDisplay()
            self.lockUnlockDisplay()


    def lockUnlockDisplay(self):
        if self.master.lockContinuous == 1:
            self.cont_entry.component('menubutton').configure(state='disabled')
        else:
            self.cont_entry.component('menubutton').configure(state='normal')
            
        if self.master.lockMin == 0:
            self.min_entry.configure(state='disabled', fg='gray40')
        else:
            self.min_entry.configure(state='normal', fg='gray0')

        if self.master.lockBMin == 1:
            self.bmin_entry.configure(state='disabled')
        else:
            self.bmin_entry.configure(state='normal') 
            
        if self.master.lockMax == 0:
            self.max_entry.configure(state='disabled', fg='gray40')
        else: 
            self.max_entry.configure(state='normal', fg='gray0')
            
        if self.master.lockBMax == 1:
            self.bmax_entry.configure(state='disabled')
        else:
            self.bmax_entry.configure(state='normal') 

        if self.master.lockIncrement == 0:
            self.incr_entry.configure(state='disabled', fg='gray40')
        else:
            self.incr_entry.configure(state='normal', fg='gray0')

        if self.master.lockBIncrement == 1:
            self.bincr_entry.configure(state='disabled')
        else:
            self.bincr_entry.configure(state='normal')

        if self.master.lockValue == 1:
            self.val_entry.configure(state='disabled', fg='gray40')
        else:
            self.val_entry.configure(state='normal', fg='gray0')

        if self.master.lockOneTurn == 1:
            self.sens_entry.configure(state='disabled', fg='gray40')
        else:
            self.sens_entry.configure(state='normal', fg='gray0')

        if self.master.lockShowLabel == 1:
            self.lab_entry.component('menubutton').configure(state='disabled')
        else:
            self.lab_entry.component('menubutton').configure(state='normal')

        if self.master.lockType == 1:
            self.if_entry.component('menubutton').configure(state='disabled')
        else:
            self.if_entry.component('menubutton').configure(state='normal')

        if self.master.lockPrecision == 1:
            self.prec_entry.component('menubutton').configure(state='disabled')
        else:
            self.prec_entry.component('menubutton').configure(state='normal')


    def updateDisplay(self):
        if self.master.opPanel.flag == 1:
            if self.master.min is None:
                self.minInput.set(self.master.labelFormat%self.master.minOld)
                self.min_entry.configure(state='disabled', fg='gray40')
            else:  
                self.minInput.set(self.master.labelFormat%self.master.min)
                self.min_entry.configure(state='normal', fg='gray0')

            if self.master.max is None:
                self.maxInput.set(self.master.labelFormat%self.master.maxOld)
                self.max_entry.configure(state='disabled', fg='gray40')
            else:  
                self.maxInput.set(self.master.labelFormat%self.master.max)
                self.max_entry.configure(state='normal', fg='gray0')

            if self.master.increment is None or self.master.increment == 0.:
                self.incrInput.set(
                    self.master.labelFormat%self.master.incrementOld)
                self.incr_entry.configure(state='disabled', fg='gray40')
            else:
                self.incrInput.set(self.master.labelFormat%self.master.increment)
                self.incr_entry.configure(state='normal', fg='gray0')

            if self.master.value == None:
                self.valInput.set(self.master.labelFormat%0)
            else:  
                self.valInput.set(self.master.labelFormat%self.master.value)

            self.sensInput.set(self.master.labelFormat%self.master.oneTurn)    


class VectorOptionsPanel:

    """ This class builds the options panel used in vectorGUI
    """

    def __init__(self, master=None, title=None):
        self.master = master
        self.root = None
        self.title = title
        if self.title is None:
            self.title = 'Options Panel'

        if self.master is None:
            self.master = Tkinter.Tk()
            self.master.continuous = 1             
            self.master.withdraw()

        self.flag = 0     # this flag toggles display/undisplay

        self.idf = InputFormDescr(title=self.title)
        
        self.idf.append({'widgetType':Tkinter.Label,
                         'wcfg':{'text': '\n'},
                         'gridcfg':{'columnspan':2, 'row':0, 'column':0},
                         })


        self.idf.append({'name':'togCont',
                         'widgetType':Pmw.OptionMenu,
                         'wcfg':{'labelpos':'w',
                                 'label_text':'Continuous       ',
                                 'menubutton_width':3,
                                 'items':('on', 'off'),
                                 'command': self.toggleCont_cb},
                         'gridcfg':{'sticky':'we',
                                    'columnspan':2, 'row':1, 'column':0},
                         })


        self.idf.append({'name':'togAxes',
                         'widgetType':Pmw.OptionMenu,
                         'wcfg':{'labelpos':'w',
                                 'label_text':'Axis Mode       ',
                                 'menubutton_width':3,
                                 'items':('XY', 'X', 'Y', 'Z'),
                                 'command': self.toggleAxes_cb},
                         'gridcfg':{'sticky':'we',
                                    'columnspan':2, 'row':2, 'column':0},
                         })


        self.idf.append({'name':'selPrec',
                         'widgetType':Pmw.OptionMenu,
                         'wcfg':{'labelpos':'w',
                                 'label_text':'Precision',
                                 'menubutton_width':3,
                                 'items':('1', '2', '3', '4', '5', '6', '7',
                                          '8', '9', '10' ),
                                 'command': self.selPrec_cb},
                         'gridcfg':{'sticky':'we',
                                    'columnspan':2, 'row':3, 'column':0},
                         })
        

        self.idf.append({'name':'DismissButton',
                         'widgetType':Tkinter.Button,
                         'text':'Dismiss',
                         'wcfg':{'bd':3},
                         'gridcfg':{'columnspan':2},
                         'command': self.Dismiss_cb})


    def toggleCont_cb(self, val):
        if val == 'off':
            self.master.configure(continuous=0)
        else:
            self.master.configure(continuous=1)


    def toggleAxes_cb(self, val):
        self.master.configure(mode=val)


    def selPrec_cb(self, val):
        self.master.configure(precision=val)


    def Dismiss_cb(self):
        self.flag = 0
	self.optionsForm.withdraw()


    def displayPanel(self, create):
        self.flag = 1
        if create == 0:
            self.optionsForm.deiconify()
        else:
            self.optionsForm = InputForm(self.master, self.root,
                                         descr = self.idf,
                                         modal = 0, blocking = 0)

            self.cont_entry  = self.idf.entryByName['togCont']['widget']
            self.mode_entry  = self.idf.entryByName['togAxes']['widget']
            self.prec_entry  = self.idf.entryByName['selPrec']['widget']

            menus = (self.cont_entry, self.mode_entry, self.prec_entry)
            Pmw.alignlabels(menus)

            w = self.idf.entryByName['togCont']['widget']
            if self.master.continuous == None or self.master.continuous == 0:
                w.setvalue('on')#i=1
            else:
                w.setvalue('off')#i=0
 
            axe = self.master.mode
            if axe == 'XY': axe=0
            elif axe == 'X': axe=1
            elif axe == 'Y': axe=2
            elif axe == 'Z': axe=3
            self.idf.entryByName['togAxes']['widget'].setitems(
                items=('XY', 'X', 'Y', 'Z'),
                index=axe)

            prc = int(self.master.precision) - 1
            if prc > 9: prc = 9
            self.idf.entryByName['selPrec']['widget'].setitems(
                items=('1', '2', '3', '4', '5', '6', '7', '8', '9', '10'),
                index = prc)

            self.updateDisplay()
            self.lockUnlockDisplay()


    def lockUnlockDisplay(self):
        if self.master.lockContinuous == 1:
            self.cont_entry.component('menubutton').configure(state='disabled')
        else:
            self.cont_entry.component('menubutton').configure(state='normal')

        if self.master.lockMode == 1:
            self.mode_entry.component('menubutton').configure(state='disabled')
        else:
            self.mode_entry.component('menubutton').configure(state='normal')
 
        if self.master.lockPrecision == 1:
            self.prec_entry.component('menubutton').configure(state='disabled')
        else:
            self.prec_entry.component('menubutton').configure(state='normal')
        

            
    def updateDisplay(self):
        self.master.entryXTk.set(
            self.master.thumbx.labelFormat%self.master.vector[0])
        self.master.entryYTk.set(
            self.master.thumby.labelFormat%self.master.vector[1])
        self.master.entryZTk.set(
            self.master.thumbz.labelFormat%self.master.vector[2])
        

if __name__ == '__main__':
    test = OptionsPanel(title="My Options Panel")
    test.displayPanel(create=1)
