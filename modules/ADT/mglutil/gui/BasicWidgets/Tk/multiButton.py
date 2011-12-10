#########################################################################
#
# Date: Nov 2002 Author: Daniel Stoffler
#
#    stoffler@scripps.edu
#
# Copyright: Daniel Stoffler and TSRI
#
#########################################################################

import types, re
import Pmw, Tkinter
from mglutil.gui.BasicWidgets.Tk.customizedWidgets import KeySelectableScrolledFrame

class kbScrolledFrame(Pmw.ScrolledFrame, KeySelectableScrolledFrame):
    """Pmw.ScrolledFrame with support to type the name of an item which will
    scroll the canvas to the requested item. Disclaimer: this works only if the items
    in the frame have an attribute name of type string. It is up to the user
    to add this attribute when building the ScrolledFrame widget."""
    
    def __init__(self, *args, **kw):

        # Pmw widgets are very delicate when it comes to subclassing!
        apply( Pmw.ScrolledFrame.__init__, (self,)+args, kw)
        
        # now remove the Pmw initialization only keywords and configure the
        # widget
        if kw.has_key('borderframe'):
            del kw['borderframe']
        if kw.has_key('labelpos'):
            del kw['labelpos']
        if kw.has_key('usehullsize'):
            del kw['usehullsize']

        if len(kw):
            apply( self.configure, (self,), kw)

        myWidget = self
        KeySelectableScrolledFrame.__init__(self, myWidget)


class MultiButtons(Tkinter.Frame):
    """This base class provides everything to build a multi-button panel.
    Buttons are packed in a Pmw.ScrolledFrame widget.
    Further information is found in MultiCheckbutton and MultiRadiobutton."""
    

    def __init__(self, master=None, valueList=None, callback=None,
                 sfcfg=None, immediate=1, **kw):
        
	Tkinter.Frame.__init__(self, master)
	Tkinter.Pack.config(self, fill='both', expand=1, side='left')

        self.master = master

        self.callback = callback # external callback which gets called
                                 # when checkbutton is checked
                                 # or unchecked
                                 # Note that the callback function to be
                                 # written by the user needs an argument
                                 # which is the widget itself. Example:
                                 # def myCallback(self, widget):
                                 #    values = widget.get()
                                 #    print values


        self.immediate = immediate #if set to 0, checking a button does not
                                   #call the callback method
                                 
        if sfcfg is None: # scrolledFrame Tkinter configuration dict
            self.sfcfg = {}
        else:
            self.sfcfg = sfcfg
            
        self.scrolledFrame = None # the Tkinter ScrolledFrame widget

        # the next 2 items are set by the widget and not by the user
        self.buttonList = [] # tuple of the button names and original names
                             
        self.buttonDict = {} # dictionary of the buttons:
                             # keys:
                             #   names of the buttons
                             # values:
                             #   origName: original names of buttons
                             #   index: number of this button 
                             #   value (value is 0 or 1: on or off)
                             #   button (value is the Tkinter.Checkbutton)
                             #   buttoncfg: dictionary that stored the Tkinter
                             #              checkbox configuration
                             

        self.rebuildOK = 1     # used in the rebuild() method
        
        # parse the valueList, build self.buttonList
        self.buttonList, self.buttonDict = self.buildButtonDict(valueList)

        # build the widget
        self.buildPanel()

        # now set the constructor options correctly using the configure method
        apply( self.configure, (),
               {'immediate':immediate,
                })

        
    def callCallback(self, event=None):
        if self.callback is not None:
            if self.immediate:
                self.callback(self)


    def buildButtonDict(self, valueList):
        if valueList is None or len(valueList) == 0:
            return [], {}
        bList = [] 
        bDict = {}

        i = 0 # used to generate unique button names
        for data in valueList:
            if data is None:
                continue

            buttoncfg = {} # dict with tkinter data for checkbutton
            dataDict = {}  # dict that stores everything about a button
 
            if type(data) == types.StringType:
                origName = data
                buttonName = data + str(i)
                buttoncfg['variable'] = Tkinter.IntVar()
                buttoncfg['command'] = self.callCallback
                dataDict['value'] = 0 # default button status: off (0)
                dataDict['name'] = buttonName

            elif type(data) == types.ListType or \
                 type(data) == types.TupleType:

                origName = data[0]
                buttonName = data[0] + str(i)
                if type(data[1]) == types.DictType:
                    if data[1].has_key('buttoncfg'):
                        buttoncfg = data[1]['buttoncfg']
                        if not buttoncfg.has_key('variable'):
                            buttoncfg['variable'] = Tkinter.IntVar()
                        if not buttoncfg.has_key('command'):
                            buttoncfg['command'] = self.callCallback
                    else:
                        buttoncfg['variable'] = Tkinter.IntVar()
                        buttoncfg['command'] = self.callCallback
                        
                    if not data[1].has_key('value'):
                        dataDict['value'] = 0 # default button status: off (0)
                    else:
                        dataDict['value'] = data[1]['value']
                    # and set the value to the Tkinter variable
                    buttoncfg['variable'].set(dataDict['value'])
            
            # now append the data to the list and dict
            dataDict['buttoncfg'] = buttoncfg
            dataDict['index'] = i
            dataDict['origName'] = origName
            bDict[buttonName] = dataDict
            bList.append( (buttonName, origName) )
            i = i + 1
            
        return bList, bDict
               
                
    def buildPanel(self):
        """To be implemented by subclass"""
        pass


    def rebuild(self, valueList):
        """This method rebuilds the panel only if the valueList has
        changed."""
        bList, bDict = self.buildButtonDict(valueList)
        self.rebuildOK = 0 # do not rebuild the panel

        statusList = self.get() # holds the origNames, check/uncheck of button

        # now lets test if we have to rebuild the panel
        if len(bList) != len(self.buttonList):
            self.rebuildOK = 1 # rebuild the panel

        for name, origName in bList:
            if not self.buttonDict.has_key(name): # a new button was found
                self.rebuildOK = 1 # rebuild the panel 
          
            for i in range(len(statusList)):
                if origName == statusList[i][0]:
                    bDict[name]['value'] = statusList[i][1]
                    del statusList[i] # and remove this entry from the list
                    break
                    
        if self.rebuildOK == 1: # we can go ahead and rebuild the panel
            self.buttonList = bList
            self.buttonDict = bDict
            self.buildPanel()


    def get(self, mode='all', event=None):
        blist = []
        if mode == 'ViPEr':
            for tupl in self.buttonList:
                name = tupl[0]
                value = self.buttonDict[name]['button'].var.get()
                blist.append( (name, value ) )
        if mode == 'checked' or mode =='checkedNames':
            for tupl in self.buttonList:
                name = tupl[0]
                value = self.buttonDict[name]['button'].var.get()
                if value == 1:
                    if mode == 'checkedNames':
                        blist.append(self.buttonDict[name]['origName'])
                    else:
                        blist.append( (self.buttonDict[name]['origName'],
                                       value ) )
            
        elif mode == 'unchecked' or mode =='uncheckedNames':
            for tupl in self.buttonList:
                name = tupl[0]
                value = self.buttonDict[name]['button'].var.get()
                if value == 0:
                    if mode == 'uncheckedNames':
                        blist.append(self.buttonDict[name]['origName'])
                    else:
                        blist.append( (self.buttonDict[name]['origName'],
                                       value ) )

        elif mode == 'all' or mode == 'allNames':
            for tupl in self.buttonList:
                name = tupl[0]
                value = self.buttonDict[name]['button'].var.get()
                if mode == 'allNames':
                    blist.append(self.buttonDict[name]['origName'])
                else:
                    blist.append( (self.buttonDict[name]['origName'], value ) )
        return blist
        

    def getName(self, index):
        return self.buttonList[index]
        

    def getIndex(self, name):
        return self.buttonDict[name]['index']
    

    def configure(self, **kw):
        for key,value in kw.items():
            if key=='immediate': self.setImmediate(value)
 

    def setImmediate(self, val):
        assert val in [0,1]
        self.immediate = val

        
    def getButtonStatus(self, buttonName):
        if len(self.buttonDict) == 0:
            return 0
        if not self.buttonDict.has_key(buttonName):
            return 0
        else:
            return self.buttonDict[buttonName]['button'].var.get()



class MultiCheckbuttons(MultiButtons):
    """This class builds a multi-checkbutton panel. Checkbuttons are packed in
a Pmw.ScrolledFrame widget. Buttons are created usinge the valueList argrument:
the valueList can be either:
    1) a list of strings: every string is used as the name of the button,
       Example: valueList=['apple', 'banana', 'cherry']. This would create 3
       checkbuttons (default status is: unchecked)
or:

    2) a list of tuples: the first item has to be a string used as the name of
       the button, the second item has to be a dictionary with the following
       architecture:
       key: 'value' and its value can be 0 or 1: this is used to check/uncheck
                    the checkbutton
       key: 'buttoncfg' and its value is another dictionary where the user
                        can specify any Tkinter Checkbutton options
   Example: valueList=[('apple', {'value':1, buttoncfg:{'command':foo,
                                     'variable':Tkinter.IntVar()} }), 'banana']
   This example creates 2 buttons and the first button is also checked.

A callback can be specified that is bound to every checkbutton. This function
gets called with 1 argument, which is the widget itself. This way, the user
can call the methods of the widget in the callback method.

Optional Tkinter configuration options can be passed into the ScrolledFrame
using the sfcfg argument (has to be a dictionary of key:value pairs).


USAGE: mb = MultiCheckbuttons(valueList=myList, callback=myCallback)
"""
    

    def __init__(self, master=None, valueList=None, callback=None,
                 sfcfg=None, immediate=1, **kw):

        MultiButtons.__init__(self, master, valueList, callback, sfcfg,
                              immediate)
        self.reGUI = None    # a Toplevel window for the regular expression


    def buildPanel(self):
        if self.buttonList is None or len(self.buttonList) == 0:
            return
        
        if self.scrolledFrame:
            self.scrolledFrame.destroy()
        self.scrolledFrame = apply( kbScrolledFrame, (self.master,),
                                    self.sfcfg)
        self.scrolledFrame.pack(padx=3, pady=3, fill='both', expand=1)
        self.frame = self.scrolledFrame.interior()

        row = 0
        col = 0
        
        for i in range(len(self.buttonList)):
            name, origName = self.buttonList[i]
            buttoncfg = self.buttonDict[name]['buttoncfg']
            value = self.buttonDict[name]['value']
            labelName = origName#+' ('+name+')'
            label = apply( Tkinter.Label, (self.frame,), {'text':labelName, } )
            label.grid(sticky='E', row=row, column=col)
            
            checkbutton = apply(Tkinter.Checkbutton, (self.frame,), buttoncfg)
            checkbutton.var = buttoncfg['variable']
            checkbutton.var.set(value)
            checkbutton.name = origName
            checkbutton.grid(row=row, column=col+1)
            self.buttonDict[name]['button'] = checkbutton
            row = row + 1


    def checkAll(self, event=None):
        """Check all buttons"""
        for tupl in self.buttonList:
            name = tupl[0]
            self.buttonDict[name]['button'].var.set(1)


    def uncheckAll(self, event=None):
        """Uncheck all buttons"""
        for tupl in self.buttonList:
            name = tupl[0]
            self.buttonDict[name]['button'].var.set(0)

        
    def invertAll(self, event=None):
        """Toggle the current state of all buttons"""
        for tupl in self.buttonList:
            name = tupl[0]
            value = self.buttonDict[name]['button'].var.get()
            if value == 0:
                value = 1
            else:
                value = 0
            self.buttonDict[name]['button'].var.set(value)
            

    def reSelect(self, pat, mode=None, event=None):
        # FIXME: regexp search should be done in the RegexpGUI class
        #        so that it can be used by others too
        """Select buttons using regular expression syntax"""
        if pat is None or pat == '':
            return
        if type(pat) != types.StringType:
            return
        pattern = re.compile(pat)
        for tupl in self.buttonList:
            name = tupl[0]
            match = pattern.search(self.buttonDict[name]['origName'])
            if match:
                if mode == 'check':
                    self.buttonDict[name]['button'].var.set(1)
                elif mode == 'uncheck':
                    self.buttonDict[name]['button'].var.set(0)


class MultiRadiobuttons(MultiButtons):
    """This class builds a multi-radiobutton panel. For further information,
    look at the doc string of MultiCheckbuttons.
"""

    def __init__(self, master=None, valueList=None, callback=None,
                 sfcfg=None, immediate=1, **kw):

        MultiButtons.__init__(self, master, valueList, callback, sfcfg,
                              immediate)

    def buildPanel(self):
        if self.buttonList is None or len(self.buttonList) == 0:
            return
        
        if self.scrolledFrame:
            self.scrolledFrame.destroy()
        self.scrolledFrame = apply( kbScrolledFrame, (self.master,),
                                    self.sfcfg)
        self.scrolledFrame.pack(padx=3, pady=3, fill='both', expand=1)
        self.frame = self.scrolledFrame.interior()

        row = 0
        col = 0

        var = Tkinter.IntVar()
        
        for i in range(len(self.buttonList)):
            name, origName = self.buttonList[i]
            buttoncfg = self.buttonDict[name]['buttoncfg']

            label = apply( Tkinter.Label, (self.frame,), {'text':origName, } )
            label.grid(sticky='E', row=row, column=col)

            buttoncfg['variable']=var
            buttoncfg['value'] = i

            radiobutton = apply(Tkinter.Radiobutton, (self.frame,), buttoncfg)
            radiobutton.var = buttoncfg['variable']
            radiobutton.name = origName
            radiobutton.grid(row=row, column=col+1)
            self.buttonDict[name]['button'] = radiobutton
            row = row + 1


class RegexpGUI:
    """This class builds a regular expression matching GUI."""
    
    def __init__(self, master=None, callback=None, **kw):

        if master is None:
            master = Tkinter.Toplevel()
        self.master = master

        self.callback = callback # the method called in self.input_cb
        self.visible = 1         # used to toggle show/hide

        self.frame = Tkinter.Frame(self.master)
	self.frame.pack(fill='both', expand=1)
        self.master.protocol('WM_DELETE_WINDOW', self.hide )

        self.inputTk = Tkinter.StringVar()
        self.reEntry = Tkinter.Entry(self.frame, textvariable=self.inputTk)
        self.reEntry.bind('<Return>', self.input_cb)
        self.reEntry.pack()

        self.radioTk = Tkinter.StringVar()
        self.radioTk.set('check')
	self.buttonCheck = Tkinter.Radiobutton(self.frame, text='Check',
                                               variable=self.radioTk,
                                               indicatoron=1,
                                               value='check').pack(side='left')
	self.buttonUncheck = Tkinter.Radiobutton(self.frame, text='Uncheck',
                                             variable=self.radioTk,
                                             indicatoron=1,
                                             value='uncheck').pack(side='left')

    def show(self):
        self.master.deiconify()
        self.visible = 1


    def hide(self, event=None):
        self.master.withdraw()
        self.visible = 0


    def toggleVisibility(self, event=None):
        if self.visible:
            self.hide()
        else:
            self.show()


    def input_cb(self, event=None):
        self.callback(self.inputTk.get(),mode=self.radioTk.get())


if __name__ == '__main__':
    def myCallback(widget, event=None):
        values =  widget.get()
        print '*****myCallback was called'
        print values
        

    import types
    txt = dir(types)

##      mb1 = MultiCheckbuttons(
##          valueList=['apple','banana',('orange',{'value':1}),'cherry'],
##          callback=myCallback)

    mb2 = MultiCheckbuttons(valueList=txt, callback=myCallback)
    # lets also open the optional regexp GUI
    panel = RegexpGUI(master=None, callback=mb2.reSelect) 

    #mb3 = MultiRadiobuttons(valueList=txt, callback=myCallback)



