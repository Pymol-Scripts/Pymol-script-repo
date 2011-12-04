## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#############################################################################
#
# Author: Sophie COON, Michel F. SANNER, Daniel Stoffler
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

#$Header: /opt/cvs/python/packages/share1.5/mglutil/gui/BasicWidgets/Tk/customizedWidgets.py,v 1.48 2008/03/21 22:30:54 annao Exp $
#
#$Id: customizedWidgets.py,v 1.48 2008/03/21 22:30:54 annao Exp $
#

import sys
import Tkinter, Pmw
import types
import tkMessageBox
from mglutil.util.callback import CallBackFunction
from mglutil.gui.InputForm.Tk.gui import InputFormDescr,InputForm,evalString
from mglutil.gui import widgetsOnBackWindowsCanGrabFocus


class KeySelectable:
    """Adds the ability to use keystrokes to quickly select items in a list.
    The widget has to be a widget supporting .bind .after, .focus_set and
    .focus_get"""
    
    def __init__(self, widget, command=None):
        self.matchCharIndex = 1 # number of char in matchString to be matched
                                # by entry
        self.afterID = None
        self.matchString = ''
        self.widget = widget
        self.command = command
        self.oldfocus = None
        widget.bind('<ButtonPress-1>', self.setFocus_cb)
	#widget.bind('<Enter>', self.Enter_cb)
	#widget.bind('<Leave>', self.Leave_cb)
	widget.bind('<KeyPress>', self.key_cb)

        if command:
            widget.bind('<Return>', self.enterKeyCb)

        self.matchCharIndex = 1 # number of char in matchString to be matched
                                # by entry
        self.afterID = None
        self.matchString = ''


    def enterKeyCb(self, event):
        if self.command:
            #self.command(self._entryfield.get())
            self.command(self.get())


    def timeOut(self, event=None):
        """resets self.matchCharIndex to 0, called after a short period of
        time if no new character has been typed"""
        self.matchString = ''
        self.matchCharIndex = 1
        self.afterID = None

        
    def key_cb(self, event=None):
        # use key strokes to select entry in listbox
        # strokes placed within 500 miliseconds are concatenated
        #print self.matchCharIndex, '|', self.matchString, '|', event.keysym
        try:
            l = self.getAllValues()
            i = 0
            str = self.matchString + event.keysym
            for e in l:
                if e[:self.matchCharIndex]==str:
                    #print e, i
                    self.clearSelection()
                    self.selectItem(i)
                    self.matchCharIndex = self.matchCharIndex + 1
                    if self.afterID is not None:
                        self.widget.after_cancel(self.afterID)
                    self.afterID = self.widget.after(500, self.timeOut)
                    self.matchString = str
                    break
                i = i + 1
        except:
            import traceback
            traceback.print_exc()

        
    def setFocus_cb(self, event=None):
        if widgetsOnBackWindowsCanGrabFocus is False:
            lActiveWindow = self.widget.focus_get()
            if    lActiveWindow is not None \
              and ( lActiveWindow.winfo_toplevel() != self.widget.winfo_toplevel() ):
                return

        self.oldfocus = self.widget.focus_get()
 	self.widget.focus_set()


    def releaseFocus_cb(self, event=None):
        if self.oldfocus is not None:
            if widgetsOnBackWindowsCanGrabFocus is False:
                lActiveWindow = self.oldfocus.focus_get()
                if    lActiveWindow is not None \
                  and ( lActiveWindow.winfo_toplevel() != self.oldfocus.winfo_toplevel() ):
                    return

            self.oldfocus.focus_set()


    def Enter_cb(self, event=None):
        if widgetsOnBackWindowsCanGrabFocus is False:
            lActiveWindow = self.widget.focus_get()
            if    lActiveWindow is not None \
              and ( lActiveWindow.winfo_toplevel() != self.widget.winfo_toplevel() ):
                return

        self.oldfocus = self.widget.focus_get()
 	self.widget.focus_set()


    def Leave_cb(self, event=None):
        if self.oldfocus is not None:
            if widgetsOnBackWindowsCanGrabFocus is False:
                lActiveWindow = self.oldfocus.focus_get()
                if    lActiveWindow is not None \
                  and ( lActiveWindow.winfo_toplevel() != self.oldfocus.winfo_toplevel() ):
                    return

            self.oldfocus.focus_set()
       
    ## These functions have to be sub-classed:
    
    def getAllValues(self):
        """returns all values in the chooser"""
        pass


    def selectItem(self, i):
        """do what has to be done to show what matches the typed string"""
        pass


    def clearSelection(self):
        """clear the current selection in the widget"""
        pass


class KeySelectableListBox(KeySelectable):
    """Adds the ability to use keystrokes to quickly select items in a
    Tkinter.ListBox."""

    def __init__(self, widget, command=None):
        KeySelectable.__init__(self, widget, command)
        

    def getAllValues(self):
        """returns all values in the chooser"""
        return self.widget.get(0, 'end')


    def selectItem(self, i):
        """do what has to be done to show what matches the typed string"""
        self.widget.see(i)
        self.widget.selection_set(i)


    def clearSelection(self):
        """clear the current selection in the widget"""
        sel = self.widget.curselection()
        if len(sel):
            for j in sel: self.widget.selection_clear(j)

    
class KeySelectableScrolledCanvas(KeySelectable):
    """Adds the ability to use keystrokes to quickly select items in a
    Pmw.ScrolledCanvas. The names of the items on the canvas have to be
    passed as a separate argument"""

    def __init__(self, widget, itemNames):
        command = None
        KeySelectable.__init__(self, widget, command)
        self.itemNames = itemNames


    def getAllValues(self):
        """returns all values in the chooser"""
        return self.itemNames


    def selectItem(self, i):
        """do what has to be done to show what matches the typed string"""
        scrollTo = float(i) / float(len(self.itemNames))
        # compensate for the first pixel:
        # 2.1 seems to work better than 2.0
        scrollTo = scrollTo - ( (scrollTo /
                                 float(len(self.itemNames)) ) /2.1 )
        self.widget.yview_moveto(scrollTo)


    def clearSelection(self):
        """clear the current selection in the widget"""
        pass


class KeySelectableScrolledFrame(KeySelectable):
    """Adds the ability to use keystrokes to quickly select items in a
    Pmw.ScrolledFrame. Disclaimer: this works only if the items
    in the frame have an attribute name of type string. It is up to the user
    to add this attribute when building the ScrolledFrame widget."""

    def __init__(self, widget):
        command = None
        KeySelectable.__init__(self, widget, command)
        self.itemNames = []


    def getAllValues(self):
        """returns all the names of all children of this frame. Please note
        that this only works if the children have an attribute name of
        type string."""
        
        itemNames = []
        itemList = self.widget.interior().children.values()
        for item in itemList:
            if hasattr(item, 'name'):
                itemNames.append(item.name)
        itemNames.sort()
        self.itemNames = itemNames
        return itemNames


    def selectItem(self, i):
        """do what has to be done to show what matches the typed string"""
        scrollTo = float(i) / float(len(self.itemNames))
        # there is some drift that has to be fixed:
        # this factor 2.2 is empiric...
        scrollTo = scrollTo - ( (scrollTo /
                                 float(len(self.itemNames)) ) /2.2 ) #FIXME
        # I am upset!! In the official Pmw documentation the method is called
        # 'yview()' but looking at the source code of Pmw widget I found out
        # they call it '_yview()'!!! Grrr!! And it doesnt work as described
        # in the manual. Even the official demos dont work.
        # this is the part of the scrolling method I am interested:
        frameHeight = self.widget._frame.winfo_reqheight()
        self.widget.startY = scrollTo * float(frameHeight)
        self.widget.reposition()


    def clearSelection(self):
        """clear the current selection in the widget"""
        pass


class ListChooser(KeySelectableListBox):
    """class to present a list of objects in a scrolled listbox.
    entries can be given as a list of items associated with a description,
    (items,description) or (item,None) if no description.
    double clicking on an entry display the associated comment in a
    textfield."""
    
    def __init__(self, root, mode='single', title='Choose', withComment=0,
                 entries=[], cwcfg=None,
                 command=None, commandEvent="<ButtonRelease-1>",
                 lbpackcfg={'fill':'both', 'expand':1}, lbwcfg={}):

        assert mode in ['single', 'browse', 'multiple', 'extended' ]

        cwcfg1={'usehullsize':1, 'hull_width':100,'hull_height':80,
                'text_wrap':'none'}
        if cwcfg is not None:
            cwcfg1.update( cwcfg1 )
        cwcfg = cwcfg1
            
        # put everything inside a Frame
        self.top = Tkinter.Frame(root, bd=4)

        # Create a scrolled text widget to display comments if not None
        if withComment:
            self.comments = apply(Pmw.ScrolledText,(self.top,),cwcfg)
            self.comments.pack(fill='both', expand=1)
        self.withComment = withComment
        
        # add a title
        self.title = Tkinter.Label(self.top,text=title)
        self.title.pack(side='top', anchor='n')

        # build the scroll list
        f = Tkinter.Frame(self.top)
        self.mode = mode
        scrollbarY = Tkinter.Scrollbar(f, orient='vertical')
        scrollbarX = Tkinter.Scrollbar(f, orient='horizontal')
        lbwcfg['selectmode'] = mode
        lbwcfg['yscrollcommand'] = scrollbarY.set
        lbwcfg['xscrollcommand'] = scrollbarX.set
        self.lb = apply( Tkinter.Listbox, (f, ), lbwcfg)

        # assign command to event
        # command can either be a function name, in that case will get
        # the event from
        # commandEvent.command can also be a list of tuple (command,
        # commandEvent) use to assign more than one command.
        if command:
            # list of tuple [(command,commandEvent),]
            if type(command) == types.ListType:
                for c in command:
                    if type(c) == types.TupleType:
                        self.lb.bind(c[1],c[0],'+')
            else:
                self.lb.bind(commandEvent,command,'+')
            
        scrollbarY.config(command=self.lb.yview)
        scrollbarY.pack(side='right', fill='y')
        scrollbarX.config(command=self.lb.xview)
        scrollbarX.pack(side='bottom', fill='x')
        apply( self.lb.pack, (), lbpackcfg)
        f.pack(fill='both', expand=1)
        self.hasInfo = 0
        KeySelectableListBox.__init__(self, self.lb, command)

        self.entries = []
        # add entries
        if entries:
            for e in entries:
                self.add(e)

        # If no comments are associated to the item the scrolledText
        # is not shown.
        #if self.hasInfo==0 and withComment==0:
        #    self.comments.forget()

        
    def add(self, entry):
        """add an entry to the list"""

        assert(type(entry)==types.TupleType)
        name = entry[0]
        if not isinstance (entry[0], str):
            name = repr(entry[0])
        self.entries.append((name,)+tuple(entry[1:]))
        self.lb.insert('end', name)
        #No comments or comments==None.
        if self.withComment and len(entry) > 1 and entry[1] != None:
            self.lb.bind("<ButtonRelease-1>", self.info, '+')
            self.hasInfo=1
            #self.comments.pack()

            
    def insert(self, pos, entry, comment=None):
        """insert an entry to the list"""
        #assert pos < len(self.entries)
        if pos == 'end':
            self.entries.append( (entry, comment) )
        else:
            self.entries.insert(pos, (entry, comment))
        if not type(entry) in types.StringTypes:
            entry = repr(entry)
        self.lb.insert(pos, entry)


    def info(self, event=None):
        """dispaly information about the curent selection"""
        if not self.withComment:
            return
        self.comments.delete(0.0, 'end')
        if len(self.lb.curselection())!=0:
            s=self.entries[map( int, self.lb.curselection())[0]][1]
            self.comments.insert('end', s+'\n')


    def clear(self):
        """clear all entries"""
        #self.lb.delete(0.0, 'end')
        self.lb.delete(0, 'end')
        self.entries = []

    def setlist(self, items):
        """ Replace all the items of the listbox component with those
        specified by the items sequence of the following type
        (entryName, comment) or(entryName, None) """
        self.clear()
        map(self.add, items)
        
    def clearComments(self):
        """clear all entries"""
        self.comments.delete(0.0, 'end')

        
    def remove(self, entry):
        """remove an entry"""
        if not len(self.entries): return
        ## ind = map(lambda x: x[0],self.entries).index(entry)
##         self.entries.remove(self.entries[ind])
##         self.lb.delete( ind )
        if type(entry) in types.StringTypes:
            ind = map(lambda x: x[0],self.entries).index(entry)
            #self.entries.remove(self.entries[ind])
        elif isinstance(entry, types.IntType):
            ind = entry

        del self.entries[ind]
        self.lb.delete( ind )

    def getInd(self,event=None):
        """ get index of current selection"""
        res =[]
        for ent in map( int, self.lb.curselection() ):
            res.append( ent)
        return res
    
    def get(self, event=None, index=0):
        """get the current selection"""
        res = []
        for ent in map( int, self.lb.curselection() ):
            res.append( self.entries[ ent ][index] )
        return res

    
    def getAll(self, event=None, index=0):
        """get the current selection"""
        res = []
        for ent in self.entries:
            res.append( ent[index] )
        return res

    
    def select(self, first, last=None):
        """add an entry to the current selection"""
        """select the given entry"""
        # need to get the index of first
        if first in self.entries:
            firstInd = self.entries.index(first)
        elif first in map(lambda x: x[0], self.entries):
            firstInd = self.entries.index((first,None))
        else:
            firstInd = first

        # need to get the index of last
        if last in self.entries:
            lastInd = self.entries.index(last)
        elif last in map(lambda x: x[0], self.entries):
            lastInd = self.entries.index((last,None))
        else:
            lastInd = last

        if lastInd and self.mode!='single':
            try:
                self.lb.selection_set(firstInd, lastInd)
            except:
                print "WARNING: Could not select the given entries %s, %s"%(first, last)
        else:
            try:
                self.lb.select_set(firstInd)
                self.lb.activate(firstInd)
                if self.withComment and self.hasInfo:
                    self.info()
            except:
                print "WARNING: Could not select the given entries %s"%(first)

    def deselect(self, first, last=None):
        """remove an entry to the current selection"""
        # need to get the index of first
        if first in self.entries:
            firstInd = self.entries.index(first)
        elif first in map(lambda x: x[0], self.entries):
            firstInd = self.entries.index((first,None))
        else:
            firstInd = first

        # need to get the index of last
        if last in self.entries:
            lastInd = self.entries.index(last)
        elif last in map(lambda x: x[0], self.entries):
            lastInd = self.entries.index((last,None))
        else:
            lastInd = last

        if lastInd:
            try:
                self.lb.select_clear(firstInd, lastInd)
            except:
                print "WARNING: Could not unselect the given entries %s, %s"%(first, last)
        else:
            try:
                self.lb.select_clear(firstInd)
                if self.withComment and self.hasInfo:
                    self.info()
            except:
                print "WARNING: Could not unselect the given entries %s"%(first)

    def set(self, item):
        """ item must be the entrie showing in """
        assert item in map(lambda x: x[0], self.entries)
        # get the entry of the listBox corresponding to the value.
        defaultEntry = filter(lambda x, df = item:
                              x[0] == df, self.entries)[0]
        
        # Get its index
        indexDefaultValue = self.entries.index(defaultEntry)

        # Set the the selection at this index using the select_set of the
        # listBox widget.
        self.lb.selection_set(indexDefaultValue)
        self.lb.activate(indexDefaultValue)
        # If has a comment show it:
        if self.withComment and len(defaultEntry)>1 and defaultEntry[1]!=None:
            self.info()
        

    def grid(self,**kw):
         #apply(self.top.grid,args,kw)
        self.top.grid(kw)
        if not kw.has_key('row'):
            row = 0
        else:
            row = kw['row']
        if not kw.has_key('column'):
            column=0
        else:
            column=kw['column']
        if kw.has_key('weight'):
            weight = kw['weight']
        else:
            weight = 1
        self.top.rowconfigure(row, weight=weight)
        self.top.columnconfigure(column, weight=weight)
        

    def pack_forget(self):
        self.top.pack_forget()

    def pack(self, **kw):
##          apply(self.top.grid,(),kw)
        self.top.pack(kw)

    def grid_forget(self):
        self.top.grid_forget()


class ObjectChooser(ListChooser):
    """Same as a ListChooser but we pass (name,object) pairs, the names are displayed but the objects handles are returned when selected.
"""
    def __init__(self, root, mode='single', title='Choose',
                 objects=None, cwcfg={'usehullsize':1,
                                      'hull_width':100,'hull_height':80,
                                      'text_wrap':'none'}, 
                 lbpackcfg={'fill':'x', 'expand':1}, lbwcfg={}, command=None):

        self.nameObject = {} # key is name in chooser, values object
                             # used to insure name unicity
        self.entries = [] # order list of (name, object) pairs

        ListChooser.__init__(self, root, mode, title, command=command,
                             cwcfg=cwcfg, lbpackcfg=lbpackcfg, lbwcfg=lbwcfg)
        for o in objects:
            self.add(o)


    def add(self, object):
        assert len(object)>=2
        name, object = object
        self.nameObject[name]=object
        if not self.nameObject.has_key(name):
            self.nameObject[name] = object

        ListChooser.add(self, (name,))

    def insert(self, pos, object):
        name, object  = object[0]
        if not self.nameObject.has_key(name):
            self.nameObject[name]=object
        ListChooser.insert(self,pos,(name,))

    def get(self, event = None):
        res = []
        if not self.lb.curselection():
            return res

        for ent in map( int, self.lb.curselection() ):
            obj= self.nameObject[self.entries[ ent ][0]]
            res.append(obj)
        if self.mode=='single' and res:
            res = res[0]

        return res

class LoadOrSaveText:
    """Class to present a scrolled text widget with a save button allowing
    you to save the content of the scrolled text in a file or display the
    content of a file in the scrolled text."""

    def __init__(self, parent,
                 textwcfg={'usehullsize':1,
                         'hull_width':400,
                         'hull_height':300,
                         'text_wrap':'none'},
                 textpackcfg={},
                 savewcfg = {},
                 loadwcfg = {}):
        # Put everything inside a Frame
        self.top = Tkinter.Frame(parent, bd=4)

        # Create the Scrolled Text Region.
        self.text = apply(Pmw.ScrolledText, (self.top,),textwcfg)
        apply(self.text.pack, (), textpackcfg)

        # Create the Save and LoadButton
        savewcfg['callback'] = self.text.exportfile
        saveButton = apply(SaveButton,(self.top,),savewcfg)
        loadwcfg['callback'] = self.text.importfile
        loadButton = apply(LoadButton,(self.top,),loadwcfg)
        saveButton.pack(side='right',expand = 1,fill = 'both')
        loadButton.pack(side='left',expand = 1,fill='both')
        
    def get(self,first = None, last = None):
        return self.text.get(first,last)

    def set(self, value):
        if not type(value) in types.StringTypes:
            return
        try:
            self.text.importfile(value)
        except:
            self.text.settext(value)

    def grid(self,**kw):
##          apply(self.top.grid,(),kw)
        self.top.grid(kw)
        
    def grid_forget(self):
        self.top.grid_forget()

    def pack(self,**kw):
        apply(self.top.pack,(),kw)

    def pack_forget(self):
        self.top.pack_forget()

class SaveButton:
    def __init__(self, root, buttonType = Tkinter.Button,
                 widgetwcfg={'text':'Save'},
                 title='Save In File ',
                 idir = None, ifile = None, types = [('All types', '*.*')],
                 callback=None):
        widgetwcfg['command']= CallBackFunction(
            self.save_cb, title, idir, ifile, types, callback)
        self.button = apply(buttonType, (root,), widgetwcfg)
        #self.master = master
        self.root = root
        self.filename = None

    def configure(self,title ='', idir=None, ifile = None, types=[('All types', '*.*')], callback=None):
        newcb = CallBackFunction(self.save_cb,title = title,idir=idir,
                                ifile=ifile,types = types,
                                callback=callback)
        self.button.configure(command =newcb)        
    
    def save_cb(self, title = '',idir=None,ifile=None,
                types = [('All types', '*.*')],
                callback = None):
        from ViewerFramework.VFGUI import fileSaveAsk
        newfile = fileSaveAsk(self.root,idir=idir,ifile=ifile,
                              types = types,title = title)
        if callback and newfile:
            callback(newfile)
        self.filename = newfile

    def get(self):
        return self.filename

    def pack(self, **kw):
        apply(self.button.pack,(),kw)

    def grid(self, **kw):
        #apply(self.button.grid,(),kw)
        self.button.grid(kw)

    def grid_forget(self):
        self.button.grid_forget()


class LoadButton:
    def __init__(self, root, buttonType =Tkinter.Button,
                 widgetwcfg={'text':'Load'},
                 title='Load File',
                 idir = None, ifile = None, types = [('All types', '*.*')],
                 callback = None):
        self.root = root
        widgetwcfg['command']= CallBackFunction(self.load_cb,
                                           title,
                                           idir,
                                           ifile,
                                           types,
                                           callback)
        self.button = apply(buttonType, (root,), widgetwcfg)
        #self.master = master
        self.filename = None

    def load_cb(self, title = '',idir=None,ifile=None,
                types = [('All types', '*.*')],
                callback = None):
        from ViewerFramework.VFGUI import fileOpenAsk
        newfile = fileOpenAsk(self.root,idir=idir,ifile=ifile,
                              types = types,title = title)
        if callback and newfile:
            callback(newfile)
        self.filename = newfile

    def get(self):
        return self.filename

    def pack(self, **kw):
        apply(self.button.pack,(),kw)

    def grid(self, **kw):
        self.button.grid(kw)
##          apply(self.button.grid,(),kw)

    def grid_forget(self):
        self.button.grid_forget()

class FunctionButton:

    def __init__(self, parent, buttonType = Tkinter.Checkbutton,
                 callback=None,
                 buttonwcfg={'text':'function'},textwcfg = {}):
        buttonwcfg['command'] = CallBackFunction(self.buildText,parent,
                                                 textwcfg,
                                                 buttonwcfg)
        self.button = apply(buttonType,(parent,),buttonwcfg)
        self.callback = callback

    def buildText(self,parent,textwcfg, buttonwcfg):
        # Build the LoadOrSaveText widget or ScrolledText
        self.idf = InputFormDescr()
        textwcfg['name']='UserFunction'
        self.text = buttonwcfg['text']
        self.idf.append(textwcfg)
        #windows = InputForm(parent, self.idf)
        master = parent
        root = None
        windows = InputForm(master,root,  self.idf)
        self.vals = windows.go()
        # Uncheck the checkbutton if needed
        if isinstance(self.button,Tkinter.Checkbutton) and \
           buttonwcfg.has_key('variable'):
            if isinstance(buttonwcfg['variable'], Tkinter.StringVar):
                buttonwcfg['variable'].set('0')
            elif isinstance(buttonwcfg['variable'], Tkinter.IntVar):
                buttonwcfg['variable'].set(0)

            

    def get(self):
        if not hasattr(self, 'vals'):
            return
        # Apply if Cancel button was checked self.result = None
        elif not self.vals:
            #self.result =(self.text,None,None)
            self.result = {}

        # No function defined is the text field.
        elif self.vals['UserFunction'] in ['\012','']:
            #self.result = (self.text,None,None)
            self.result = {}

        else:
            # Text lets remove comments
            # Function defined
            funcString = self.vals['UserFunction']
            function = evalString(funcString)

            textwcfg = self.idf[0]
            if self.callback:
##                 self.result = (self.text,
##                                self.callback(function), funcString)
                self.result = (self.text,
                               self.callback(function),function)
                               
            else:
##                 self.result = (self.text,None,funcString)
                self.result = (self.text,None,function)
        # Have to send back a tuple (name of the widget, result, function)

        return self.result
    
    def set(self, value):
        pass

    def pack(self, **kw):
        apply(self.button.pack,(),kw)

    def grid(self, **kw):
##          apply(self.button.grid,(),kw)
        self.button.grid(kw)

    def grid_forget(self):
        self.button.grid_forget()


#from DejaVu.extendedSlider import ExtendedSlider
#from DejaVu.Slider import Slider

#from DejaVu.EventHandler import CallbackFunctions
from mglutil.gui.BasicWidgets.Tk.eventHandler import CallbackFunctions
import types


class SliderWidget(CallbackFunctions):
    """Class for a simple slider"""

    options = ['labelfont', 'cursorfont', 'immediate', 'command', 'incr',
               'width', 'minval', 'maxval', 'height', 'label',
               'withValue', 'labelsCursorFormat', 'sliderType', 'lookup']
    

    def __init__(self, master=None, label=None, minval=0.0, maxval=100.0,
		 incr=None, init=None, width=150, height=20, withValue=1,
		 immediate=1, left=10, right=10 ,
                 labelsCursorFormat = '%4.2f', sliderType = 'float',
                 command = None, lookup=None):
        
        # Calls the base class constructor
	CallbackFunctions.__init__(self)

        # Create a toplevel if no master specified
        if master is None:
            self.master = Tkinter.Toplevel()
        else:
            self.master = master

        # Create the master frame
	self.frame = Tkinter.Frame(master)

        # Set some attributes
	self.lastx=0
	self.immediate = immediate
        self.labelsCursorFormat = labelsCursorFormat
        self.sliderType = eval(sliderType)

        # Create the label widget
	if not label: label = 'slider'
	fnt='-*-helvetica-medium-r-narrow-*-*-120-*-*-*-*-*-*'
	self.label = Tkinter.Label(self.frame, text=label, font=fnt)
	self.label.pack(side=Tkinter.LEFT)

	if withValue: height = max(30, height)
	self.withValue = withValue

        # Create the canvas.
	self.draw = Tkinter.Canvas(self.frame, width=width, height=height,
                                   relief=Tkinter.SUNKEN)
	self.width = width
	self.height = height
	self.left = left
	self.right = width-right

        # MS May 1st add support for discrete set of values
        self.lookup = lookup # can be set to a list of values indexed by
                             # int(value)
        if lookup is not None:
            # force min and max to provide proper indices
            # SC minval needs to be an int !
            minval = 0
            maxval = len(lookup)-1
            incr = 1
            
	self.max=maxval
	self.min=minval
	self._ComputeCst()
	self.incr = incr
	self.incrCanvas = None

        # Compute the increment on the canvas based on the given incr value.
	if incr:
	    self.incrCanvas = round(incr*self.cst)

	if withValue: m = self.height / 2
	else: m = int(self.height * 0.75)
	self.middle = m

	self.draw.create_line( self.left, m, self.right, m,
			       width=2, fill='black', tag = 'line')

	y = m-10 # 10 is the cursor's height
	l = self.left
	self.cursor = self.draw.create_polygon( l, m, l+5, y, l-5, y, l, m,
						fill='blue', tag='cursor')
        self.valueLabel = None
        # create the cursor label which needs to be formatted properly.
	if withValue:
	    y = self.middle+10
            lv = (self.labelsCursorFormat)%minval
	    self.valueLabel = self.draw.create_text( l, y,
                                                     text= lv,
						     font=fnt, tag='cursor')

	self.draw.pack(side = Tkinter.RIGHT)

	if init is None:
            if self.lookup:
                # SC init must be a value and not an index
                init = self.lookup[int(round(self.min))]
            else:
                init = self.min

	self.set(init)
	Tkinter.Widget.bind(self.draw, "<1>", self.MoveCursor)
	Tkinter.Widget.bind(self.draw, "<B1-Motion>", self.MoveCursor)
	if not immediate:
	    Tkinter.Widget.bind(self.draw, "<ButtonRelease-1>", self.MouseUp)
        self.command = command
        if command:
            self.AddCallback(command)


    def Callbacks(self):
	"""Implement call to all callbacks"""
        if self.lookup:
            val = self.lookup[int(round(self.val))]
        else:
            val = self.val

	for f in self.callbacks:
	    f(val)


    def MoveCursor(self, event):
	"""Callback function for left mouse button"""

	x = event.x - self.left
	self._MoveCursor(x)


    def DrawCursor(self):
	"""Update graphics representatin of the cursor"""

	# compute position
	x = (self.val-self.min)*self.cst
	deltax = x - self.lastx

	if self.withValue:
            if self.lookup:
                val = self.lookup[int(round(self.val))]
            else:
                val = self.val
	    self.draw.itemconfig(self.valueLabel,
                                 text = (self.labelsCursorFormat)%val )
	self.draw.move('cursor', deltax, 0 )
	self.lastx = self.lastx + deltax


    def _MoveCursor(self, x):
	"""Compute value and move the cursor to the new position"""

	# compute the new value
	val = self.min+self.cst1 * x
	if self.incrCanvas:
	    val = round(val/self.incr)*self.incr
	if val<self.min: val = self.min
	elif val>self.max: val = self.max

	# check if redraw and callback are needed
	if self.val != val:
	    self.val = val
	    self.DrawCursor()
	    if self.immediate: self.Callbacks()


    def set(self, val, update=1):
	"""Set the cursor"""
        # Here needs to get the index of the given val in the
        # lookup table.
        # If the given val is smaller than the min then
        # set to the min if val bigger than the max then
        # set to the max....
        if self.lookup:
            if not val in self.lookup:
                raise ValueError
            val = self.lookup.index(val)
            if val<self.min:
                val = self.min
            elif val>self.max:
                val = self.max
        else:
            if val<self.min: val = self.min
            elif val>self.max: val = self.max
        
	if self.incrCanvas:
	    val = round(val/self.incr)*self.incr
	self.val = val

	self.DrawCursor()
	if update: self.Callbacks()
	#return self.val

    def setMin(self, val):
        """Set the minimum value of the slider"""
        if self.lookup:
            if val > self.lookup[int(round(self.max))]: return
            if val in self.lookup: self.min = self.lookup.index(val)
            else:
                raise ValueError

        else:
            if val>self.max: return
            self.min = val

        self._ComputeCst()
        if self.val < val:
            self.set(val)
        else:
            self.DrawCursor()

    def setMax(self, val):
        """Set the maximum value of the slider"""
        if self.lookup:
            if val < self.lookup[int(round(self.min))]: return
            if val in self.lookup: self.max = self.lookup.index(val)
            else:
                raise ValueError
        else:
            if val<self.min: return
            self.max = val

        self._ComputeCst()
        if self.val > val:
            self.set(val)
        else:
            self.DrawCursor()

    def get(self):
	"""Get the slider's value"""
        if self.lookup:
            val = self.lookup[int(round(self.val))]
        else:
            val = self.val
        
	return val


    def _ComputeCst(self, event=None):
	"""Conversion factor between values and pixels"""
	self.cst = 1.0*(self.right-self.left) / (self.max-self.min)
	self.cst1 = 1.0 / self.cst



    def configure(self, **kw):
        """Method to modify the configuration options of the widget."""
        
        if kw == {} :
            return self.options
        else:
            keys = kw.keys()
            for key in keys:
                if key not in self.options:
                    print "WARNING: option %s is not available for the Slider widget" % key
            lf = kw.get('labelfont')
            if lf:
                self.label.configure(font=lf)

            cf = kw.get('cursorfont')
            if cf:
                if self.withValue:
                    self.draw.itemconfig(self.valueLabel, font=cf)
                    
            imm = kw.get('immediate', None)
            if imm is not None:
                assert imm == 1 or imm == 0
                if self.immediate != imm:
                    self.immediate = imm
                    if imm == 1: #immediate mode
                        Tkinter.Widget.unbind(self.draw, "<ButtonRelease-1>")
                    else: #not immediate mode
                        Tkinter.Widget.bind(self.draw, "<ButtonRelease-1>",
                                            self.MouseUp)

            com = kw.get('command')
            if com:
                if self.command:
                    self.RemoveCallback(self.command)
                self.AddCallback(com)
                self.command = com
            
            incr = kw.get('incr')
            if incr:
                self.incrCanvas = round(incr*self.cst)
                self.incr = incr

            width = kw.get('width')
            if width:
                w_type = type(width)
                assert w_type == types.IntType or w_type == types.FloatType
                if width != self.width:
                    self.draw.configure(width = width)
                    old_width = self.width
                    right = old_width - self.right
                    self.right = width - right
                    self._ComputeCst()
                    self.draw.coords('line', self.left, self.middle,
                                     self.right, self.middle)
                    self.width = width
                    self.DrawCursor()

            height = kw.get('height')
            if height:
                h_type = type(height)
                assert h_type == types.IntType or h_type == types.FloatType
                if self.withValue: height = max(30, height)
                if height != self.height:
                    self.draw.configure(height=height)
                    self.height = height
                    if self.withValue: m = height / 2
                    else: m = int(height * 0.75)
                    self.middle = m
                    l = self.left
                    self.draw.coords('line', l, m, self.right, m)
                    y = m-10 # 10 is the cursor's height
                    x = self.draw.coords(self.cursor)[0]
                    self.draw.coords(self.cursor, x, m, x+5, y,
                                     x-5, y, x, m)
                    if self.withValue:
                        y = m+10
                        self.draw.coords(self.valueLabel, x, y)
                    
            maxval = kw.get('maxval', None)
            if maxval is not None:
                self.setMax(maxval)

            minval = kw.get('minval', None)
            if minval is not None:
                self.setMin(minval)

            label = kw.get('label')
            if label:
                self.label.configure(text = label)

            withval = kw.get('withValue', None)
            if withval is not None:
                if withval != self.withValue:
                    self.withValue = withval
                    x = self.draw.coords(self.cursor)[0]
                    if withval:
                        height = max(30, self.height)
                        if height != self.height:
                            self.draw.configure(height=height)
                            self.height = height
                        m = self.height / 2
                        y = m+10
                        fnt='-*-helvetica-medium-r-narrow-*-*-120-*-*-*-*-*-*'
                        sv = (self.labelsCursorFormat)%withval
                        self.valueLabel = self.draw.create_text(x, y,
                                                     text= sv,
						     font=fnt, tag='cursor')
                    else:
                        m = int(self.height * 0.75)
                        if self.valueLabel:
                            self.draw.delete(self.valueLabel)
                            self.valueLabel = None
                    y = m-10 # 10 is the cursor's height
                    self.draw.coords('line', self.left, m,
                                     self.right, m)
                    self.draw.coords(self.cursor, x, m, x+5, y,
                                     x-5, y, x, m)
                    self.middle = m

            lcf = kw.get('labelsCursorFormat')
            if lcf:
                self.labelsCursorFormat = lcf
                self.set(self.val, update=0)

            st = kw.get('sliderType')
            if st:
                self.sliderType = st
                self.set(self.val, update=0)


    def grid(self, **kw):
        #apply(self.frame.grid,(),kw)
        self.frame.grid(kw)
        
    def grid_forget(self):
        self.frame.grid_forget()
    def pack(self, **kw):
        self.frame.pack(kw)

    def pack_forget(self):
        self.frame.pack_forget()

class ExtendedSliderWidget(SliderWidget):
	#"""Composite slider which includes a simple slider 
	#   and an entry window used to set self.value """
	
    def __init__(self, master=None, label=None,
                 minval=0.0, maxval=100.0,
		 incr=None, init=None, width=150, height=20, withValue=1,
                 immediate=0, left=10, right=10, entrywcfg = None,
                 entrypackcfg = {'side':'top'},
                 labelsCursorFormat = '%4.2f',
                 sliderType = 'float',command = None, lookup=None):

        if not entrywcfg:
            entrywcfg = self.entrywcfg = {'width':5}

        if entrywcfg.has_key('textvariable'):
            self.entryContent = entrywcfg['textvariable']
        else:
            self.entryContent = entrywcfg['textvariable'] = Tkinter.StringVar()
            
	try:
            SliderWidget.__init__(self, master=master, label=label,
                                  minval=minval,maxval=maxval,
                                  incr=incr, init=init, width=width,
                                  height=height, withValue=withValue, 
                                  immediate=immediate,left=left,
                                  right=right,
                                  labelsCursorFormat=labelsCursorFormat,
                                  sliderType=sliderType,
                                  command=command, lookup=lookup)
	except ValueError:
            print "problem w/Slider.__init__"

        # CREATE THE ENTRY
        if entrywcfg.has_key('label'):
            # Label defined for the entry
            entrylabel = Tkinter.Label(self.frame)
            #FIXME!!!!
            #apply(entrylabel.pack, (), entrypackcfg)
            entrylabel.pack(side='top', before = self.draw, anchor = 'w')
            label = entrywcfg['label']
            entrylabel.config(text=label)
            # remove this entry from the configuration dictionary
            del entrywcfg['label']
        else:
            label = None

	self.entry = apply(Tkinter.Entry, (self.frame,), entrywcfg)
        entrypackcfg['before'] = self.draw
        apply(self.entry.pack, (), entrypackcfg)
        # Bind the setval method to the entry and Return is the event.
        self.entry.bind('<Return>', self.setval)

        # Set the entryContent to the init 
        self.entryContent.set(self.labelsCursorFormat%self.get())
        

    def MoveCursor(self, event):
        """Callback function for left mouse button"""
        SliderWidget.MoveCursor(self,event)
##         x = event.x - self.left
##         self._MoveCursor(x)
## 	val =  x
## 	val = self.min+self.cst1 * x
## 	if self.incrCanvas:
## 	    val = round(val/self.incr)*self.incr
## 	if val<self.min: val = self.min
## 	elif val>self.max: val = self.max
##         if self.lookup:
##             val = self.lookup[int(round(self.val))]
##         else:
##             val = self.val
        
        if self.lookup:
            sval = self.lookup[int(round(self.val))]
        else:
            sval = self.val
 	self.entryContent.set(self.labelsCursorFormat%sval)

    def _FixedMoveCursor(self, x):
	"""move cursor to new position which is typed 
	   into the entry window"""
	# compute the new value
        # Problem with the lookup
	val =  x
	if self.incrCanvas:
	    val = round(val/self.incr)*self.incr
	if val<self.min: val = self.min
	elif val>self.max: val = self.max
	# check if redraw and callback are needed
	if self.val != val:
	    self.val = val
	    self.DrawCursor()
	    if self.immediate: self.Callbacks()

    def set(self, val, update=1):
        """Set Both the cursor and entry window"""
        SliderWidget.set(self, val, update=update)
        if self.lookup:
            sval = self.lookup[int(round(self.val))]
        else:
            sval = self.val
        self.entryContent.set(self.labelsCursorFormat%sval)
   
    def setval(self, event):
	"""Bound to Button-1"""

	try:
            newx=self.sliderType(self.entryContent.get())
            if self.lookup:
                if not newx in self.lookup:
                    raise ValueError
                    #return
                val = self.lookup.index(newx)
            else:
                val = newx
            self._FixedMoveCursor(val)
            if not self.immediate:
                self.Callbacks()

	except ValueError:
            if self.lookup:
                val = self.lookup[int(round(self.val))]
            else:
                val = self.val
            self.entryContent.set('')

    def grid(self, **kw):
        #apply(self.frame.grid,(),kw)
        self.frame.grid(kw)

    def grid_forget(self):
        self.frame.grid_forget()


    def pack(self, **kw):
        self.frame.pack(kw)

    def pack_forget(self):
        self.frame.pack_forget()



class kbComboBox(Pmw.ComboBox, KeySelectableListBox):

    def __init__(self, *args, **kw):
        apply( Pmw.ComboBox.__init__, (self,)+args, kw)
        cmd=None
        if kw.has_key('selectioncommand'):
            cmd = kw['selectioncommand']
        lb = self.component('scrolledlist').component('listbox')
        KeySelectableListBox.__init__(self, lb, command = cmd)
        self.master = self.widget.master

    def selectItem(self, i, run=1):
        self.selectitem(i)
        # check if i is an index or an item in the list before calling
        # the command with it.
        val = self.get()
        if self.command and run:
            self.command(val)
        

    def set(self, value, run=1):
        self.selectitem(value)
        if self.command and run:
            self.command(value)


    def _addHistory(self):
        """ overwrite Pmw.ComboBox._addHistory, 
so the case where input == '' is accessible
"""
        input = self._entryWidget.get()
        if input == '': 
            self.command(input)
        else:
            apply( Pmw.ComboBox._addHistory, (self,) , )



class kbScrolledListBox(Pmw.ScrolledListBox, KeySelectableListBox):
    """
    Adds the ability to use keystrokes to quickly select items in a
    Pmw.ScrolledListBox.
    """
    def __init__(self, *args, **kw):
        apply(Pmw.ScrolledListBox.__init__, (self,)+args, kw)
        cmd=None
        if kw.has_key('selectioncommand'):
            cmd = kw['selectioncommand']
        lb = self.component('listbox')
        KeySelectableListBox.__init__(self, lb, command = cmd)
        self.master = self.widget.master
	self.widget.bind('<Enter>', self.Enter_cb)
	self.widget.bind('<Leave>', self.Leave_cb)

    def selectItem(self, i, run=1):
        self.select_set(i)
        self.see(i)
        self.activate(i)
        if self.command and run:
            self.command()
        

    
