#########################################################################
#
# Date: Oct 2006 Authors: Michel Sanner
#
#    sanner@scripps.edu
#
# Copyright: Michel Sanner and TSRI
#
#########################################################################

import sys
from Tkinter import Toplevel, Label

from mglutil.gui import widgetsOnBackWindowsCanGrabFocus

class KeyboardEntry:
    """Mixin class that can be used to add the ability to type values for
widget.  When the cursor enters the widget the focus is set to the widget
and the widget that had the focus previousy is saved in .lastFocus.
When the cursor leaves the widget focus is restored to self.lastFocus.

key_cb(event) is called upon KeyRelease events. If this is the first keytroke
an undecorated window is created next to the cursor and handleKeyStroke(event)
is called.  This method should be subclassed to modify the way the string is
typed in (see thumbwheel.py for an example of handling numbers only).
When the Return key is pressed, the entry window is destroyed and if the
typed striong is not empty self.callSet(self.typedValue) is called.

If the cursor moves out of the widget before Rrturn is pressed the value
in the entry window is discarded and the entry window is destroyed.
"""
    def __init__(self, widgets, callSet):
        # callSet is a function taking one argument (string) which will be
        # when the Returnm key is hit (if the typed string is not empty)
        # widgets is a list of Tkinter wigets for which <Enter>, <Leave>
        # <Return> and <KeyRelease> are handled. These widgets has to exist
        # before thsi constructor can be called
        assert callable(callSet)
        self.callSet = callSet
        for w in widgets:
            w.bind("<Enter>", self.enter_cb)
            w.bind("<Leave>", self.leave_cb)
            w.bind("<Return>", self.return_cb)
            w.bind("<KeyRelease>", self.key_cb)

        self.typedValue = '' # string accumulating valuescharacters typed
        self.lastFocus = None # widget to which the focus will be restored
                              # when the mouse leaves the widget
        self.topEntry = None # top level for typing values
        self.typedValueTK = None #Tk label used as entry for value


    def key_cb(self, event=None):
        # call back function for keyboard events (except for Return)
        # handle numbers to set value
        key = event.keysym
        if key=='Return': return

        #print 'Key:', key
        if len(self.typedValue)==0 and self.topEntry is None:
            # create Tk label showing typed value
            x = event.x
            y = event.y

            #print 'create entry'
            self.topEntry = Toplevel()
            self.topEntry.overrideredirect(1)
            w = event.widget
            if event.x >= 0:
                self.topEntry.geometry('+%d+%d'%(w.winfo_rootx() + event.x+10,
                                                 w.winfo_rooty() + event.y))
            else:
                self.topEntry.geometry('+%d+%d'%(w.winfo_rootx() +10,
                                                 w.winfo_rooty() ))
            self.typedValueTK = Label(
                master=self.topEntry, text='', relief='sunken', bg='yellow')
            self.typedValueTK.pack()

        self.handleKeyStroke(event)
        

    def leave_cb(self, event=None):
        # make sure widget gets keyboard events
        # print 'leave', event.widget

        if self.topEntry:
            self.typedValueTK.destroy()
            self.topEntry.destroy()
            self.topEntry = None
        self.typedValue = ''
        if self.lastFocus:
            #print 'restoring focus'
            if widgetsOnBackWindowsCanGrabFocus is False:
                lActiveWindow = self.lastFocus.focus_get()
                if    lActiveWindow is not None \
                  and ( lActiveWindow.winfo_toplevel() != self.lastFocus.winfo_toplevel() ):
                    return

            self.lastFocus.focus_set()
            self.lastFocus = None
            event.widget.config(cursor='')


    def enter_cb(self, event=None):
        # make sure widget gets keyboard events
        #print 'enter', event.widget

        if widgetsOnBackWindowsCanGrabFocus is False:
            lActiveWindow = event.widget.focus_get()
            if    lActiveWindow is not None \
              and ( lActiveWindow.winfo_toplevel() != event.widget.winfo_toplevel() ):
                return

        if self.lastFocus is None:
            #print 'setting focus'
            self.lastFocus = self.focus_lastfor()
        event.widget.focus_set()
        event.widget.config(cursor='xterm')


    def return_cb(self, event):
        # return should destroy the topEntry
        #print "return_cb"
        if self.typedValueTK is not None:
            self.typedValueTK.destroy()
        if self.topEntry is not None:
            self.topEntry.destroy()
        self.topEntry = None
        if len(self.typedValue):
            #print 'setting to', self.type(self.typedValue)
            self.callSet(self.typedValue)
        self.typedValue = ''


    ## TO BE SUBCLASSED

    def handleKeyStroke(self, event):
        # by default we handle delete character keys
        key = event.keysym
        if key=='BackSpace' or key=='Delete':
            self.typedValue = self.typedValue[:-1]
            self.typedValueTK.configure(text=self.typedValue)
