## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#############################################################################
#
# author: Yong Zhao
#
# Copyright: Yong Zhao, TSRI 2003
#
#############################################################################
"""
Package: 
Module : 
This module provides 4 classes. TK based GUI (graphical user interface) of
1) tree-view
   The tree view is designed to display hierarchical data structure.
   each level of the tree can be expanded / collapsed to show / hide subtree
   one or more items of the tree could be selected by mouse click

2) list-view
   The list view is designed to display a list of objects
   each item in the list can be locked (to prevent deletion) and unlocked
   the list is mainteined by FIFO principle (first in, first out), unless locked by user
  
"""
import sys
import os
import Tkinter, types
from warnings import warn
from Tkinter import Scrollbar, Listbox, Button, Label, Toplevel
import Pmw
from time import sleep
import types

import numpy.oldnumeric as Numeric

try:
    from PIL import Image, ImageTk
except:
    pass

from mglutil.util.packageFilePath import findFilePath
from os import path
# find path to TreeWidget.icons directory
ICONDIR= path.split( findFilePath('folder.gif',
                  'mglutil.gui.BasicWidgets.Tk.TreeWidget.icons') )[0]

OFFSET= 20

class ListItem:
    """
ListItem Class, defines the items in ListView class.

Key properties include draw (display), toggle (save or unsaved).
the saved / unsaved mark is drawn to the left of the item name.

    """
    def __init__(self, list=None, name=None, highlight='yellow'):
        """
Constructor of ListItem Class.
list: the list that self (ListItem object) belongs to
name: name of the item that will be displayed
"""
        self.name   =name
        self.icon   =None    # canvas ID of the "save" icons
        self.caption=None    # canvas ID of the "save" icons
        self.list   = list   # the list that self (listItem) belongs to
        self.y      = 0
        self.x      = 2  
        self.locked = False
        self.tree   = None
        #self.selected = False
        self.canvasIDs = []
        self.highlight = highlight
        
    
    def Draw(self):
        """
Draw the ListItem Object on canvas.
"""
##         if self.icon:   # if has a representation, remove them
##             canvas.delete(self.icon)
##             canvas.delete(self.caption)
        canvas = self.list.canvas

        if len(self.canvasIDs):
            for id in self.canvasIDs[:]:
                canvas.delete(id)
            #self.canvasIDs=[]
            
        if self.list.selection == self:
            idx = (self.y*OFFSET -2 + OFFSET/2)/OFFSET

            if idx <1 or idx > self.list.length:
                return

            box =self.list.selectionBox
            canvas=self.list.canvas
            if box:
                canvas.delete(box)

            if self.highlight is not None:
                box=canvas.create_rectangle(
                    2+OFFSET, 2+OFFSET * idx - OFFSET/2,
                    400, 2+OFFSET * (idx+1) -OFFSET/2, 
                    fill = self.highlight, outline="")

            self.list.selectionBox = box
            self.canvasIDs.append(box)

            
        if self.locked:
            img = self.list.pinDown_icon
        else:
            img = self.list.pinUp_icon
        
        h= self.y
        self.icon = canvas.create_image(2, 2+OFFSET*h,
                                image=img, anchor='w')
        self.caption = canvas.create_text(2+OFFSET ,2+OFFSET*h,
                                text=self.name, anchor='w') 

        self.canvasIDs.append(self.icon)
        self.canvasIDs.append(self.caption)
        
        canvas.tag_bind(self.icon, "<1>", self.Toggle_cb)
        #canvas.tag_bind(self.caption, "<Double-Button-1>", self.Chosen_cb)
        canvas.tag_bind(self.caption, "<1>", self.PickItem_cb)
        
        lcanvas = canvas.component('canvas')
        balloon = Pmw.Balloon(lcanvas)
        balloon.tagbind(lcanvas, self.icon, "if checked it won't go down in history")

         
    def Toggle_cb(self,event=None):
        """
call back function when icon is clicked.
toggle the status between 'locked' and 'unlocked'
"""
        self.locked = not self.locked
        self.Draw()

        
    def Chosen_cb(self, event=None):
        """ call back function when picking the current history item, 
the node cooresponding to this item should be selected in the TreeView
"""
        # fixme.. will be removed..
        #print 'DOUBLE PICK'
        return
        if self.tree:
            self.tree.Select(self.name)
            
            
    def PickItem_cb(self, event=None):
        """ Pick the current history item, (hilight with yellow box)
the node cooresponding to this item should be selected in the TreeView
"""
##         """callback for left mouse picking event"""
##         print "bar"
##         # check if we clicked inside the list
##         if event.x<OFFSET+2 or event.y < 2+OFFSET/2:
##             return

        #print 'PICK'
        self.list.selection = self
        self.Draw()
        
        if self.tree:
            self.tree.Select(self.name)
            
 
        
            
    
class ListView:
    """
    The ListView Class defines a Microsoft-style view of a list. 
    Items in the list ( ListItem object ) can be:
        hightlighted, interted, deleted, saved / unsaved.
"""
    
    def __init__(self, tree, master=None, name='ListView', multi_choice=False,
                list_length=10, width=200 ):
        """ Constructor, build the GUI of a list
tree               : the TreeView object associated with this list
master=None        : TK master
name='ListView'    : Name of this ListView object
multi_choice=False : allow multiple choice (more than one selection)
list_length=10     : maximum length of the list 
width=200          : whe width of the list, in pixels
    
"""
        self.master     = master
        self.length     = 0     # length of the list
        self.ItemList   = []    # list of all the items (ListItem object)
        self.ObjList    = []    # object list
        self.locked     = []    # whether the item is locked (saved)
        self.current    = -1    # index of current selection
        self.tree       = tree  # the cooresponding TreeView object
        self.number=list_length # how many items can be stored in the list
        self.selectionBox = None
        self.width=width
        self.selection  = None  # the selected list item
        
	self.canvas = Pmw.ScrolledCanvas(
            master, borderframe=1, #labelpos='n', label_text='main',
            usehullsize=0, hull_width=200, hull_height=400,
            vscrollmode='dynamic', hscrollmode='dynamic')
##         self.canvas = Tkinter.Canvas(master, width=200, height=400,
##                                      borderwidth=2)
        
##         canvas=self.canvas
##         canvas.scrollX = Tkinter.Scrollbar(self.master,
##                                                 orient=Tkinter.HORIZONTAL)
##         canvas.scrollY = Tkinter.Scrollbar(self.master,
##                                                 orient=Tkinter.VERTICAL)

##         canvas['xscrollcommand'] = self.canvas.scrollX.set
##         canvas['yscrollcommand'] = self.canvas.scrollY.set
##         canvas.scrollX['command'] = self.canvas.xview
##         canvas.scrollY['command'] = self.canvas.yview
##         canvas.scrollX.pack(side='bottom', fill='x')
##         canvas.scrollY.pack(side='right', fill='y')
        
        #self.canvas.bind("<Button-1>", self.pickNode_cb)
        self.canvas.pack(side='left', expand=1, fill='both')

        join = path.join
        down_icon  = Image.open(join(ICONDIR, "pinDown.gif"))
        self.pinDown_icon = ImageTk.PhotoImage(master=master, image=down_icon)
        up_icon  = Image.open(join(ICONDIR, "pinUp.gif"))
        self.pinUp_icon = ImageTk.PhotoImage(master=master, image=up_icon)
                
        
    def FindListItembyName(self, name):
        """
name :  the full name of the item to be found 
returns the index of a given item in the list
returns -1 if not found
"""
        list=self.ItemList
        for idx in range(self.length):
            if list.name==name:
                return idx
        return -1
    
    def Insert(self, name=None, object=None,locked=False):
        """ Insert a ListItem object in the list.
name  : name of the new ListItem to be displayed
object: the object to be added
locked: lock the object
"""
        if name==None:
            return
        for i in self.ItemList:
            if i.name ==name:
                warn( name +" is already on the list")
                return
        
        # if list is full, delete the first unsaved item
        if self.length == self.number:
            space_available = False
            idx = 0
            for i in self.ItemList:
                if i.locked==False:
                    self.Delete(idx)
                    space_available=True
                    break
                idx += 1
            if not space_available:
                warn( "List is full, please remove some saved items first")
                return
        
        # append to the end of list
        new_item=ListItem(self, name, highlight=None) # create a new obj
        new_item.tree=self.tree
        new_item.locked=locked
        new_item.y=self.length+1
        self.length +=1
        
        self.ItemList.append(new_item)
        self.ObjList.append(object)
        
        new_item.Draw()
        
        bb = self.canvas.bbox(Tkinter.ALL)
        canvas = self.canvas.component('canvas')
        canvas.configure( scrollregion=(0, 0,bb[2]+OFFSET, bb[3]+OFFSET))


##     def pickNode_cb(self, event=None):
##         """callback for left mouse picking event"""
##         print "bar"
##         # check if we clicked inside the list
##         if event.x<OFFSET+2 or event.y < 2+OFFSET/2:
##             return

##         # find out which item in the list was picked
##         idx = (event.y -2 + OFFSET/2)/OFFSET
##         if idx <1 or idx > self.length:
##             return
        
##         box =self.selectionBox
##         canvas=self.canvas
##         if box:
##             canvas.delete(box)
                
##         self.selectionBox=canvas.create_rectangle(
##                                     2+OFFSET, 2+OFFSET * idx - OFFSET/2,
##                                     400, 2+OFFSET * (idx+1) -OFFSET/2, 
##                                     fill = 'yellow', outline="")
##         canvas.lower(self.selectionBox)
##         self.current = idx-1

        
    
    def Delete(self, index):
        """
Delete the index-th ListItem object from the history list 
index : specifies which item to be removed,  range [0, len_of_list-1]
        """
        if index >= self.length or index <0:
            warn( "The list index specified is out of range")
            return
        to_remove = self.ItemList[index]
        if to_remove.locked:
            warn( "Can't delete saved item. Uncheck the save mark")
            return
        # delete the representation from canvas
        self.canvas.delete(to_remove.icon)
        self.canvas.delete(to_remove.caption)
        
        # If the item to be deleted is selected, remove the selection box
        if self.current==index:
            self.canvas.delete(self.selectionBox)
            self.current_selected = None
        
        self.ItemList.remove(to_remove)
        self.length -= 1
        if index <= self.length:
            self.MoveUp(index, self.length)
        return            

    def Delete_multi(self, index=[]):
        """ Delete the items from the list. 
The items to be removed are specified by index
        """
        removeList=[]
        for idx in index:
            if idx >= self.length or idx <0:
                warn( "The list index specified is out of range")
                return
            to_remove = self.ItemList[idx]
            removeList.append(to_remove)
            if to_remove.locked:
                warn( "Can't delete saved item. Uncheck the save mark")
                return
            # delete the representation from canvas
            self.canvas.delete(to_remove.icon)
            self.canvas.delete(to_remove.caption)
            
            # If the item to be deleted is selected, remove the selection box
            if self.current==idx:
                self.canvas.delete(self.selectionBox)
                self.current_selected = None
            
        for r in removeList:
            self.ItemList.remove(r)
            #del r
            
        # Update GUI of the list
        self.length -= len(index)
        i=1
        for item in self.ItemList:
            item.y=i
            item.Draw()
            i+=1


    def MoveUp(self, begin, end):
        """After deleting an item, move the items below up by one unit
a private function: not meant to be called by user"""
        for n in range(begin, end):
            to_move = self.ItemList[n]
            to_move.y -= 1
            to_move.Draw()

        
    def GetCurrentSelection(self):
        """
returns current selected ListItem object, None if nothing is selected
        """
        if self.current != -1:
            return self.ItemList[self.current]
        else:
            return None
        
    

    def SaveItem(self, index = 0):
        """
Save the n-th ListItem object in the history list 
n is specified as index, range [0, len_of_list-1]
After saving, the item cannot be deleted unless unsaved by UnSaveItem function
        """
        if index >= self.length or index <0:
            warn( "The list index specified is out of range")
            return
        item = self.ItemList[index]
        if item.locked ==False:
            item.locked = True
            item.Draw()

    
    def UnSaveItem(self, index = 0):
        """
Unsave the n-th ListItem object in the history list 
n is specified as index, range [0, len_of_list-1]
After unsaving, the item can be deleted unless saved by SaveItem function
        """
        if index >= self.length or index <0:
            warn( "The list index specified is out of range")
            return
        item = self.ItemList[index]
        if item.locked ==True:
            item.locked = False
            item.Draw()

    def LockItem(self, index = 0):
        """ same as calling SaveItem """
        self.SaveItem(index=index)

        
    def UnlockItem(self, index = 0):
        """ same as calling UnSaveItem """
        self.UnSaveItem(index=index)        
            
                
 
class Node:
    """ Tree node object
        Integrated with TreeView. (Each item in TreeView is a Node object)
"""
    def __init__(self, name, object=None, mouseBinding=None,
                 hasChildren=False, firstExpand_cb=None):
        """Constructor
name               # name of the node to be displayed
object = object    # the object associated with this Node        
objectKey= None    # the key to look up object
haschildren =False # any children node ?
firstExpand_cb     # call back function when first expand the node
        """
        self.name   = name
        self.visible = False
        self.object = object    # the object associated with this Node        
        self.objectKey= None    # the key to look up object
        self.children = []
        self.parent = None      # handle to parent node
        self.maxy   = None      # maximum y value used by this node
        self.canvasIDs  = []    #   will store list of Canvas items 
                                #   used to draw this node
        self.uniqueID   = None  # unique number associated with 
                                #   this node in a tree
        self.tag=None           # Tkinter tag
        self.tree=None          # tree object to which this node is added
        self.expanded=False
        self.x=self.y=OFFSET    # the x,y of the nodes
        self.deleted=False

        # FIXME this should really be fullName of this node
        self.parentFullname=None# Full name of parent
        self.height=1           # height of the subtree, 
                                #   including this node,all its children
        self.update=False       # if true, will redraw the node
                                #   when expanding or colapsing
        self.inserted=False     # this node is being inserted/deleted
        self.selected=False     # this node is selected ?
        self.selectboxID=None   # the canvas ID for selection box

        # if node added as 'hasChildren', + and folder will be drawn
        # however, the children will NOT be added to save loading time
        self.hasChildren = hasChildren

        self.firstExpanded = True
        if self.hasChildren:
            self.firstExpanded = False


        # the first time to expand a node that 'hasChildren'=True,
        # the call_back function
        self.firstExpand_cb = None
        # FIXME: We do not pass this through the constructor right now
        # but use the method self.addFirstExpandCallback()
        if firstExpand_cb is not None:
            try:
                assert callable(firstExpand_cb)
                self.firstExpand_cb = firstExpand_cb
            except:
                warn( "ERROR: Callback %s not callable!"%firstExpand_cb)
            
        self.mouseBinding = {'<1>':self.pickNode_cb,
                             '<Double-Button-1>':self.double1PickNode_cb,
                             '<Triple-Button-1>':self.addToHistory_cb,
                             '<Double-Button-2>':self.double2PickNode_cb,
                             '<Double-Button-3>':self.double3PickNode_cb,
                             '<Alt-Button-1>':self.addToHistory_cb,
                             '<Control-Button-1>':self.ctrlPickNode_cb,
                             '<Shift-Button-1>':self.shiftPickNode_cb}
        if mouseBinding:
            self.mouseBinding.update(mouseBinding)
        
        # temp varibles used for deleting a node
        self.x1 = self.x2 = self.y1 = self.y2 = self.move_up =0
        

    def addFirstExpandCallback(self, cb):
        """ """
        try:
            assert callable(cb)
            self.firstExpand_cb = cb
        except:
            warn(
                    "ERROR: Callback %s not callable!"%cb)
            self.firstExpand_cb = None
            
        
    def findChildByName(self, master,name):
        """ returns the child node anmesd 'name'. 
return None if not found
"""
        path = name.split('|')
        for c in master.children:
            tmp = c.name.split('|')
            if tmp[-1]==path[0]:
                newName = name[len(path[0])+1:]
                if len(newName):
                    res = self.findChildByName(c,newName)
                    return res
                else:
                    return c
        return None


    def invoke(self,event=None):
        """
Callback function when the plus/minus icon is clicked
switch between the expanded / collapsed status
expand the tree if not expanded or collapsed the tree if already expanded
        """
        
        if self.expanded:
            self.Collapse()
        else:
            self.Expand()
            
        
    def ctrlPickNode_cb(self, event=None):
        """ used for multiple selection, 
        control + left mouse button down
        Attention: not fully implimented yet, 12-1-2004
        """
        #print 'CTRL PICK'
        tree = self.tree
        tree.selectionHistory.append(tree.list_selected)
        if self.selected:
            tree.list_selected.remove(self)
            self.selected = False
            tree.canvas.delete(self.selectboxID)
        else:
            self.tree.list_selected.append(self)
            self.selected = True
            self.drawSelectionBox()

            
    def triple1PickNode_cb(self, event=None):
        """
        customizable, override
        """
        #print 'TRIPLE'
        #self.tree.double1PickNode_Func(self)
        


    def double1PickNode_cb(self, event=None):
        """
        customizable, override
        """
        #print 'DOUBLE'
        self.tree.double1PickNode_Func(self)
        

    def double2PickNode_cb(self, event=None):
        """
        customizable, override
        """
        self.tree.double2PickNode_Func(self)
        

    def double3PickNode_cb(self, event=None):
        """
        customizable, override
        """
        self.tree.double3PickNode_Func(self)
        

    def shiftPickNode_cb(self, event=None):
        """ used for multiple selection, 
        shift + left mouse button down
        Attention: not fully implimented yet, 12-1-2004
        """
        tree = self.tree
        if len(tree.list_selected)==0:
            return
        last = tree.list_selected[-1]
        all = self.parent.children
        try:
            i1 = all.index(self)
            i2 = all.index(last)
        except ValueError:
            raise ValueError, "range only work over sibblings"

        if i1 > i2:
            tmp=i1; i1=i2; i2=tmp

        tree.selectionHistory.append(tree.list_selected)
        for node in all[i1:i2+1]:
            if node==last: continue
            if node.selected:
                tree.list_selected.remove(node)
                node.selected = False
                tree.canvas.delete(node.selectboxID)
            else:
                tree.list_selected.append(node)
                node.selected = True
                node.drawSelectionBox()


    def pickNode_cb(self,event=None):
        """callback function bound to left mouse picking on tree canvas"""
        #self.selected = not self.selected
        #print 'PICK'
        tree = self.tree
        for node in tree.list_selected:
            node.selected = False
            node.tree.canvas.delete(node.selectboxID)
        tree.selectionHistory.append(tree.list_selected)
        tree.list_selected = [self]
        self.selected = True
        self.drawSelectionBox()

        if self.selected:
            if tree.actions['select'] is not None:
                tree.actions['select'](self)
        else:
            if tree.actions['deselect'] is not None:
                tree.actions['deselect'](self)

    
    def displayChildren(self,master):
        """ Display all the children nodes of given master node
"""
        for c in master.children:
            c.draw()
            
            if c.children and c.expanded:
                c.displayChildren(c)
            
            
    def updateHeight(self,master):
        """ Recursively update the height of parents
"""
        if master is not None:
            master.height =1
            if master.expanded:
                for c in master.children:
                    #if c.expanded:
                    master.height += c.height
                self.updateHeight(master.parent)
            else:
                master.height =1
            
                
    def hideChildren(self):
        """ Recursively remove the representation of all the children
"""
        for c in self.children:
            c.visible = False
            for id in c.canvasIDs[:]:
                self.tree.canvas.delete(id)
                c.canvasIDs.remove(id)            
            c.hideChildren()

    def updateX (self):
        """ Update the coordination x value
"""
        if self.parent==None:
            self.x=OFFSET
        else:
            self.x=self.parent.x+OFFSET
            
        for c in self.children:
            c.x = c.parent.x + OFFSET
            c.updateX()

        return

        

    def updateVLines (self, current, master):
        """ Recursively redraw the vertical lines to the left of the node,
which is being expanded or collapsed
"""        
        if master is not None:
            if current==master.children[-1]:
                self.updateVLines(master, master.parent)
                return
            # if current node is not the last child of parent    
            index = master.children.index(current)
            x1 = master.children[index+1].x
            y1 = master.children[index+1].y
            v_line = self.tree.canvas.create_line(current.x,current.y+5,
                        x1,y1, fill='gray') 
            self.tree.canvas.lower(v_line)
            master.children[index+1].canvasIDs.append(v_line)
            # recursively updateVLines
            self.updateVLines(master, master.parent)
            return
        else:
            # we are now at the roots..
            if current==self.tree.roots[-1]:
                return
            index = self.tree.roots.index(current)
            x1 = self.tree.roots[index+1].x
            y1 = self.tree.roots[index+1].y
            v_line = self.tree.canvas.create_line(current.x,current.y+5,
                        x1,y1, fill='gray') 
            self.tree.canvas.lower(v_line)
            self.tree.roots[index+1].canvasIDs.append(v_line)
            return
        # end of UpdateVLines
        

    def updateVLines_after_Del (self, current, master):
        """private function, not meant to be called by user."""
        if master is not None:
            if current==master.children[-1]:
                # the current was already erased
                # not more vertical line needed under current node
                self.updateVLines_after_Del(master, master.parent)
                return
            # if current node is not the last child of parent 
            index = master.children.index(current)
            x1 = current.x
            # print "from", current.name,
            if index ==0 :
                if current.deleted:
                    y1 = master.y
                else:
                    y1 = master.children[index-1].y
            else:
                if current.deleted :
                    y1 = master.children[index-1].y
                else:
                    y1 = master.children[index+1].y

            v_line = self.tree.canvas.create_line(current.x,current.y,
                        x1,y1, fill='gray') 
            self.tree.canvas.lower(v_line)
            master.children[index+1].canvasIDs.append(v_line)
            # recursively updateVLines
            self.updateVLines_after_Del(master, master.parent)
            return
        
        else: # we are now at the roots..
            if current==self.tree.roots[-1] or current==self.tree.roots[0]:
                return
            index = self.tree.roots.index(current)
            x1 = current.x
            if current.deleted:
                y1 = self.tree.roots[index-1].y
            else:
                y1 = self.tree.roots[index+1].y
            
            v_line = self.tree.canvas.create_line(current.x,current.y,
                        x1,y1, fill='gray') 
            self.tree.canvas.lower(v_line)
            self.tree.roots[index+1].canvasIDs.append(v_line)
            return        
        
        
    def Expand(self):
        """ The Expand() function will expand current node, show the "minus"
        icon, display all the children node.
"""
        if self.expanded:
            return
        canvas = self.tree.canvas
        canvas.configure(cursor='watch')
        canvas.update_idletasks()
        prev_selection = self.tree.current_selected
        self.expanded = not self.expanded
        self.update = True
        self.draw()

        if self.firstExpanded==False:
            if self.firstExpand_cb:
                self.firstExpand_cb(node=self, object=self.object) 
                self.firstExpanded = True
            #else:
            #    warn("Warning: firstExpand_cb not defined!")
        
        # redraw the selection node.. 
        # to make sure the yellow box is in lower layer of background
        if prev_selection and prev_selection != self:
            prev_selection.draw()
        canvas.configure(cursor='')
      
    def Collapse(self):
        """ The Collapse() function will collaps current node, show the "plus"
        icon, hide all the children node.
"""
        if not self.expanded:
            return
        self.expanded = not self.expanded
        prev_selction = self.tree.current_selected
        self.update = True
        self.draw()
         
    def getHeight(self):
        """Returns the Height of node.
        height, stored in self.height, is the layer number of all the children
        and grand-children. 
        e.g. Every leaf node displayed has height of 1. The height of a 
        parent node with 4 leaf nodes is 4+1=5 """
        h=0
        if self.children:
            for c in self.children:
                if c.expanded:
                    h += c.getHeight()
                else:
                    h += 1
        return h + 1 # the node itself has a representation, so,  + 1
    
    
    def UpdateNode(self):
        """self.update is a flag indicating whether expanding or collapsing
is triggered. The UpdateNode function handle the expanding or collapsing
of the tree 
self.update=True means toggle the expanding status. ( expand <-> collapse )
"""
        self.update = False
        all_items = self.tree.canvas.bbox(Tkinter.ALL)
        xx = all_items[2]+100
        yy = all_items[3]+OFFSET

        if self.expanded:                
            # now delete / re-draw the vertical lines
            self.tree.canvas.addtag_overlapping("VLines", 
                        0,self.y + 6, self.x, self.y + 6)
            v_line_removed = self.tree.canvas.find_withtag("VLines")
            for i in v_line_removed:
                self.tree.canvas.dtag(i, "VLines")
                self.tree.canvas.delete(i)

            #self.tree.canvas.addtag_overlapping("tmp",0,self.y, xx, yy)
            self.tree.canvas.addtag_enclosed("tag_"+self.name,0,self.y, xx, yy)

            items_moved = self.tree.canvas.find_withtag("tag_"+self.name)

            self.height = self.getHeight()            
            self.updateHeight(self.parent)
            self.tree.updateY()

            #self.displayChildren(self) # display all children

            # update the vertical lines to the left of the node
            self.updateVLines(self,self.parent)
            # move down the partition 
            # that is right below the expanded part
            deltaY = self.height - 1
            if deltaY > 0:
                self.tree.canvas.move("tag_"+self.name, 0, deltaY*OFFSET)

            self.displayChildren(self) # display all children

            # delete tag "tmp"                 
            for i in items_moved:
                self.tree.canvas.dtag(i, "tag_"+self.name) 

            #resize the canvas' scrolling region
            bb = self.tree.canvas.bbox(Tkinter.ALL)
            self.tree.canvas.configure(
                    scrollregion=(0, 0,bb[2]+OFFSET, bb[3]+OFFSET) )
            # end of Expand
        else:
            deltaY = self.getHeight() - 1
            self.height = 1
            self.updateHeight(self.parent)
            self.tree.updateY()

            # now delete / re-draw the vertical lines
            vTag = 'V_' + self.name
            self.tree.canvas.addtag_overlapping(vTag, 
                        0,self.y - 6, self.x, self.y + 6)
            #self.tree.canvas.create_line(0,self.y + 10, self.x, self.y + 10)

            v_line_removed = self.tree.canvas.find_withtag(vTag)

            for i in v_line_removed:
                self.tree.canvas.dtag(i, vTag)
                self.tree.canvas.delete(i)

            self.hideChildren()

            mTag='M_' + self.name
            self.tree.canvas.addtag_enclosed(mTag, 0, self.y, xx,yy)
            if deltaY > 0:
                self.tree.canvas.move(mTag, 0,-deltaY*OFFSET)

            # update the vertical lines to the left of the node
            self.updateVLines(self,self.parent)

            # delete tag "tmp" 
            items_moved = self.tree.canvas.find_withtag(mTag)
            for i in items_moved:
                self.tree.canvas.dtag(i, mTag) 

            # resize the canvas' scrolling region
            bb = self.tree.canvas.bbox(Tkinter.ALL)
            self.tree.canvas.configure(
                    scrollregion=(0, 0,bb[2]+OFFSET, bb[3]+OFFSET) )
            # end of Collapse

                
    def Erase(self):
        """ erase the representation of the node from canvas
            All the children nodes will also be erased.
        """
        if self.parent:
            if self.parent.expanded == False:
                return
        if len(self.canvasIDs): # the node has a representation            
            for id in self.canvasIDs[:]:
                self.tree.canvas.delete(id)
                self.canvasIDs.remove(id)
                
        self.hideChildren()
        
        self.x1     = 0
        self.y1     = self.y + (self.height -1)*OFFSET
        all_items   = self.tree.canvas.bbox(Tkinter.ALL)
        if all_items:
            self.x2     = all_items[2]+100
            self.y2     = all_items[3]+OFFSET
            self.move_up= -self.height * OFFSET
        else:   # all_items is None, nothing left on canvas
            self.x2 = self.x1
            self.x2 = self.x1
        
        # save the old height for future use (in vertical line drawing)
        self.deleted = True 
        self.height  = 0
        self.updateHeight(self.parent)
        self.tree.updateY()
        
        # Move up the region below the current node
        if self.isShown():
            self._moveUp()
        
        
    def _moveUp(self):
        """find the vertical lines to delete / re-draw
NOTE: private function. should not be called directly from user
        """
        self.tree.canvas.addtag_overlapping("V_lines", 
                0, self.y+10, self.x, self.y+10)
        v_line_removed = self.tree.canvas.find_withtag("V_lines")
        for i in v_line_removed:
            self.tree.canvas.dtag(i, "V_lines")
            self.tree.canvas.delete(i)
        
        # find the region to move upward
        self.tree.canvas.addtag_enclosed("tmp", 
                    self.x1, self.y1, self.x2, self.y2)
        self.tree.canvas.move("tmp", 0, self.move_up)
        items_moved = self.tree.canvas.find_withtag("tmp")
        for i in items_moved:
            self.tree.canvas.dtag(i, "tmp")    
    
        #update the vertical lines to the left of the node
        self.updateVLines_after_Del(self,self.parent)


    def draw_new_insert(self, num=1, mode='single'):
        """this function draws the canvas representation when one new node is
added as the child of 'self'.
Two modes are alowed:
single : a new child node just added
batch   : a set of children nodes added.
        """
        if self.parent:
            if self.parent.expanded == False:
                return
        if len(self.canvasIDs): # the node has a representation
            for id in self.canvasIDs[:]:
                self.tree.canvas.delete(id)
                self.canvasIDs.remove(id)
        
        if self.expanded:
            if mode == 'single':  #
                newNode = self.children[-1]
                #insert to the end of children list
                last_y = newNode.y
                last_x = newNode.x
            elif mode =='batch':
                last_x = self.children[0].x
                last_y = self.children[0].y
            else:
                warn('Unknown mode')                

            all_items = self.tree.canvas.bbox(Tkinter.ALL)
            if all_items:
                xx = all_items[2]+100
                yy = all_items[3]+OFFSET

                self.tree.canvas.addtag_enclosed("tmp", 0, last_y-10, xx, yy)
                self.tree.canvas.move("tmp", 0, OFFSET*num)

                # now delete the vertical lines 
                self.tree.canvas.addtag_overlapping("V_lines", 
                            0,last_y + 10, last_x, last_y + 10)
                v_line_removed = self.tree.canvas.find_withtag("V_lines")
                for i in v_line_removed:
                    self.tree.canvas.dtag(i, "V_lines")
                    self.tree.canvas.delete(i)
                #re-draw the vertical lines to the left of the node
                self.updateVLines(self,self.parent)
                # delete tag "tmp" 
                items_moved = self.tree.canvas.find_withtag("tmp")
                for i in items_moved:
                      self.tree.canvas.dtag(i, "tmp")        

            if mode == 'single':
                newNode.draw()
            else:
                for n in self.children:
                    n.draw()
 
        # draw the line between plus/minus icon and the folder icon        
        self.canvasIDs.append( 
                self.tree.canvas.create_line(self.x,self.y,
                self.x + OFFSET,self.y,fill='gray') )
        # draw vertical line above node
        if self.parent:
            parent=self.parent
            if self == parent.children[0]:
                x1 = parent.children[0].x
                y1 = parent.y
            else:
                index= parent.children.index(self)
                x1 = parent.children[index-1].x
                y1 = parent.children[index-1].y
            
            v_line = self.tree.canvas.create_line(
                        self.x,self.y,x1,y1,fill='gray') 
            self.tree.canvas.lower(v_line)
            self.canvasIDs.append(v_line)    
        else:
            roots = self.tree.roots
            if self != roots[0]:
                index=roots.index(self)
                x1   =roots[index-1].x
                y1   =roots[index-1].y
                v_line = self.tree.canvas.create_line(
                            self.x,self.y,x1,y1,fill='gray') 
                self.tree.canvas.lower(v_line)
                self.canvasIDs.append(v_line)    
            
        self.draw_caption_icon()      


##     def getValueOfObject(self):
##         """ returns a string that represents the value of the object associated with current Node"""
##         obj = self.object
##         t=type(obj)
##         if t is types.StringType: return obj
##         if t is types.NoneType  : return 'None'        
##         if t in [types.BooleanType, types.FloatType, types.IntType,
##                  types.LongType, types.NoneType]:
##             return str(obj)
##         if hasattr(obj,'__class__'):            
##             return '<'+obj.__class__.__name__+' Instance>'
##         else:
##             return "Unknown "+str(t)


    def draw_caption_icon(self):
        """Display the name of the node, as well as the folder icon
"""
        self.visible = True
        fullText = self.name
        # add ... to the end if caption is too long
        if len(fullText) > 60:
                caption = fullText[:57] + '...'
        else:   caption = fullText

        if self.expanded:
            signText='-'
#            img = self.tree.open_folder_icon
            plus_minus = self.tree.minus_icon
        else:
            signText='+'
#            img = self.tree.close_folder_icon
            plus_minus =self.tree.plus_icon
        
        canvas=self.tree.canvas
        if len(self.children) or self.hasChildren:
            # draw folder and +/-
##             folder=canvas.create_image(self.x+OFFSET,self.y,
##                     image=img, anchor=Tkinter.CENTER) 
##             self.canvasIDs.append(folder)
##             canvas.lift(folder)
##             signBox = canvas.create_rectangle(self.x+OFFSET-25,self.y+5,
##                                               self.x+OFFSET-15,self.y-5,
##                                               fill='white')
##             sign=canvas.create_text(self.x+OFFSET-20,self.y, text=signText)
##             canvas.tag_bind(sign, "<1>", self.invoke)            

            sign = canvas.create_image(self.x,self.y,
                image=plus_minus, anchor=Tkinter.CENTER)
            self.canvasIDs.append(sign)
            canvas.lift(sign)

            # bind the click action with function "invoke"
            canvas.tag_bind(sign, "<1>", self.invoke)            
##             canvas.tag_bind(folder,"<Double-Button-1>", self.invoke)
##             txt = canvas.create_text(self.x + 2*OFFSET ,self.y,
##                                      text=caption,anchor=Tkinter.W) 
##             canvas.lift(txt)
##         else:
##             txt = canvas.create_text(self.x + OFFSET ,self.y,
##                                      text=caption,anchor=Tkinter.W)
        txt = canvas.create_text(self.x + OFFSET ,self.y,
                                 text=caption,anchor=Tkinter.W)         
        canvas.tag_bind(txt, "<Button-3>", self.showFullText1_cb)
        self.canvasIDs.append(txt)

        # THIS IS EXPANSIVE
        #balloon = Pmw.Balloon(self.tree.topFrame_1)
        #balloon.tagbind(canvas, txt, repr(self.object).replace('\\n','\n'))

        #'This is help for\nan arc item')#fullText)

        # Set the mouse binding
        for action,func in self.mouseBinding.iteritems():
            if action == "<Control-Button-1>":
                if self.tree.multi_choice:
                    canvas.tag_bind(txt, action,func , '+')
            else:
                canvas.tag_bind(txt, action,func , '+')

##        canvas.tag_bind(txt, "<1>", self.pickNode_cb, '+')
##        canvas.tag_bind(txt, "<Double-Button-1>", self.addToHistory_cb,'+')
       
##         if self.tree.multi_choice:
##             canvas.tag_bind(txt, 
##                     "<Control-Button-1>", self.ctrlPickNode_cb, '+')

        if self.tree.displayValue:
            self.displayValue()

            
    def displayValue(self):
        """Add things to the rigth side of the tree
"""
        canvas = self.tree.canvas
        if type(self.object) == Numeric.arraytype:
            sh = self.object.shape
            dim=1
            for i in range(len(sh)):
                dim *= sh[i]
            #print "Large array found !!", dim
            if dim > 10:
                counter = 0
                text = str(sh) 
                text +=' array('
                text += self._getElement(self.object, \
                                    counter=counter)[0]                    
                text += '...)'
            else:
                text = repr(self.object)
        else:
            text = repr(self.object)
        text= text.replace('\n', '')
        #if len(text)>80:
        #    text = text[:77] + "..."
        valueStr = canvas.create_text(self.x + OFFSET + 150,self.y, \
                                      text=text, anchor=Tkinter.W)
        canvas.tag_bind(valueStr, "<1>", self.pickNode_cb, '+')
        canvas.tag_bind(valueStr, "<Button-3>", self.showFullText2_cb)
        # add canvas ID to list so they will be erased and moved properly
        self.canvasIDs.append(valueStr)            


    def _getElement(self, array,  counter):
        """recursive function for displaying first 20 elements of large,
multi-dimentional array"""
        shape = array.shape
        sh = len(shape)
        element = '['
        if sh ==1:
            c=0
            for i in range(shape[0]):
                e = repr(array[i])
                counter +=1
                element += e +', '
                if counter >= 20:
                    return element, counter
            element += '], '
        else:
            c=counter
            for i in range(shape[0]):
                e,c = self._getElement(array=array[i],  counter=c)
                element += e               
                counter = c
                if c >= 20:
                    element += '], '
                    break
            
        return element, counter
        

    def showFullText1_cb(self, event=None):
        print self.name
        

    def showFullText2_cb(self, event=None):
        print repr(self.object).replace('\\n','\n')

        
    def addToHistory_cb(self, event=None):
        """ Callback function. When double clicked on the tree node, the node 
        will be added to the history list"""
        tree =self.tree
        tree.AddToHistoryList()
        
        # please be noticed that by double-clicking the history item
        # the system would first take action on the click event, then the
        # double click event.
        self.selected=True
        self.draw()
        
        
    def rename(self, newName):
        """this function renames the node and updates the caption on canvas """

        if type(newName) is not types.StringType:
            return
        self.name = newName
        self.draw()


    def draw(self):
        """place items on the canvas to represent that node.
if the node has a current representation it will be deleted first.
put a folder is node has children"""

        parent = self.parent
        canvas = self.tree.canvas
        canvasIDs = self.canvasIDs
        if parent:
            if parent.expanded == False:
                return
        if self.update:
            self.UpdateNode()       
        
        if len(canvasIDs):        # erase the current representation of 'self'
            for id in canvasIDs:
                canvas.delete(id)
            #if self.selectboxID:
            #    canvas.delete(self.selectboxID)

        # draw the '-' line between plus/minus icon and the folder icon        
        canvasIDs.append( 
                self.tree.canvas.create_line(self.x,self.y,
                self.x + OFFSET,self.y,fill='gray') )
        # draw vertical line above node
        if parent:
            if self == parent.children[0]:
                x1 = parent.children[0].x
                y1 = parent.y
            else:
                index= parent.children.index(self)
                x1 = parent.children[index-1].x
                y1 = parent.children[index-1].y
                
            v_line = self.tree.canvas.create_line(
                        self.x,self.y,x1,y1,fill='gray') 
            canvas.lower(v_line)
            canvasIDs.append(v_line)    
        else:
            # if this is a new root
            if self != self.tree.roots[0]:
                index = self.tree.roots.index(self)
                x1 = self.tree.roots[index-1].x
                y1 = self.tree.roots[index-1].y
                h_line = self.tree.canvas.create_line(self.x,self.y,
                        x1,y1, fill='gray') 
                canvas.lower(h_line)
                canvasIDs.append(h_line)    
        
        self.visible = True
        self.draw_caption_icon()
        self.drawSelectionBox()


    def drawSelectionBox(self):
        # Draw the selection box
        canvas = self.tree.canvas
        canvasIDs = self.canvasIDs
        prev_selection = self.tree.current_selected
        if self.selected:
            if prev_selection and not self.tree.multi_choice:
                if prev_selection != self:
                    prev_selection.selected = False
                    if prev_selection.parent:
                        if prev_selection.parent.expanded:
                            prev_selection.draw()
                    else:
                        prev_selection.draw()
                    
            # save current node as current_selected
            self.tree.current_selected = self
            bb = self.tree.canvas.bbox(Tkinter.ALL)
            if len(self.tree.list_selected)>0:
                yoff=1
            else:
                yoff=0

            selectionBox = canvas.create_rectangle(
                1, self.y - OFFSET/2, bb[2]-yoff, self.y+OFFSET/2, 
                tags="selected",
                fill = "yellow")
            canvas.lower(selectionBox)
            self.selectboxID = selectionBox
            canvasIDs.append(selectionBox)
        
    
    def DoSelect(self, parent, path=None):    
        """ Called by TreeView function Select()
            The DoSelect function select a node, if the node's parent 
            is collapsed, the parent will be expanded (recursively)
        """
        result = False
        for c in parent.children:
            if c.name == path[0]:
                path.remove(path[0])
                if len(path)==0:
                    c.selected = True
                    c.draw()
                    return True
                else:
                    if not c.expanded:
                        c.invoke()      
                    return c.DoSelect(c,path)
        
    def GetFullName(self):
        """Returns the full name of current node """
        if self.parent:
            return self.parentFullname +'|'+self.name
        else:
            return self.name

    def _expandAll(self):
        """recursively expand all the subtrees"""
        for c in self.children:
            if not c.expanded:
                c.invoke()
            c._expandAll()
            

    def _collapseAll(self):
        """recursively collapse all the subtrees
"""
        if self.isShown() is False: return 
        for c in self.children:            
            if c.expanded:
                c.Collapse()
                c._collapseAll() # collapse all the children will do the job

    def increaseParentHeight(self, offset=1):
        """ After adding a new node, the height of the parent is increased by 1
        This functions also update the height of all the acestors."""
        self.height += offset
        if self.parent:
            # recursively increase the height of parent by 1
            self.parent.increaseParentHeight(offset=offset)


    def isShown(self):
        """check if self is shown in the tree """
        
        parent=self.parent
        if parent is not None:
            if parent.expanded == False: # collapsed
                return False
            else:
                return self.parent.isShown()
        else:
            # root node is always shown
            return True       


    def getRootNode(self):
        """ returns the root node of current branch where the node liad in"""

        name=self.GetFullName().split('|')[0]
        for root in self.tree.roots:
            if root.name == name:
                return root
        return None
        


    def moveUp(self):
        """ move self up by one position in the children nodes list of self.parent"""
        ## NOTE: This is achieved by collapse parent node and rebuild the children node list.
        if self.parent is None: # self is one of the root nodes
            try:
                nodeList=self.tree.roots
                index=nodeList.index(self)
            except:
                return
        else:
            nodeList=self.parent.children
            index=nodeList.index(self)

        if index==0:
            print "already at the top of children node list"
            return
        
        nodeList.remove(self)
        nodeList.insert(index-1, self)
            
        if self.isShown():
            if self.parent:
                self.parent.invoke()
                self.parent.invoke()
            else: # no parent: one of the roots
                self.redraw_all_roots()
        else:
            pass

    def moveDown(self):
        """ move self down by one position in the children nodes list of self.parent"""
        ## NOTE: This is achieved by collapse parent node and rebuild the children node list.
        if self.parent is None: # self is one of the root nodes
            try:
                nodeList=self.tree.roots
                index=nodeList.index(self)
            except:
                return
        else:
            nodeList=self.parent.children
            index=nodeList.index(self)

        if index==len(nodeList)-1:
            print "already at the end of children node list"
            return
        
        nodeList.remove(self)
        nodeList.append(self)
            
        if self.isShown():
            if self.parent:
                self.parent.invoke()
                self.parent.invoke()
            else: # no parent: one of the roots
                self.redraw_all_roots()
        else:
            pass


    def redraw_all_roots(self):
        """ redraw all the nodes in root.nodes """
        ## FIXME: the selection box is still in wrong position
        
        nodeList=self.tree.roots
        for n in nodeList:
            n.Collapse()
            n.Erase()
            n.height=1
        nodeList[0].y=OFFSET
        self.tree.updateY()
        for n in nodeList:
            n.draw_caption_icon()

    
class TreeView:
    """Tree Widget class
    A TreeWiget contains tree nodes (object of Class Node). 
    Each node can have children. Nodes that do have children can be expanded 
and collapsed using the + - icon placed before de nodes' icon. Each no has an 
icon and a name. It is possible to associate an arbitrary Python object to each 
node. The node in the tree is the graphical representation of the object.
"""
    def __init__(self, master=None, name='Tree', multi_choice=False,
                 width=200, height=200, treeWidth=140, treeHeight=100,
                 historyWidth=100, historyHeight=100, mode='Extended',
                 historyVisible=False,nohistory=False,
                 mouseBinding=None,obj2Node=True, displayValue=False,
                 offx=0, offy=0, canvasHeaderCol=False):
        """Constructor:
Creates the canvas used to draw the tree
"""     
        self.multi_choice = multi_choice
                                # True: ctrl or shift click to select
                                # False: only one node can be selected
        
        self.name=name 	        # name of the window if this widget 
                                # is in its own window
                                
        self.numberOfNodes=0    # high water mark of number of nodes
                                # used to create a unique id for each node
                                # which can be used to create tags as well
        
        self.roots=[]           # define list of root nodes
        self.selectionHistory = [] # only works in multi_choice mode
        self.list_selected=[]   # the list of nodes picked by mouse
                                # in multi_choice mode
        self.current_selected=None  # the seleced node
                                    # in single selection mode
        self.width = width      # width of paned widget (Tree and history)
        self.height = height    # height of paned widget

        self.obj2Node = obj2Node # True or False,false we do not make
                                 # the correspondance between object and a node
        self.objToNode = {}  # lookup to find the TreeNode that
                             # correspond to an object for a given application
        self.objIndex=0      # index for objects in this tree
        self.displayValue =displayValue  # whether to display value of objects

        self.nohistory = nohistory # if True do not create the history frame
        self.mouseBinding = mouseBinding # dictionary of event to be
                                         # associate with a mouse button event
                                         # when a node is selected
                                         # see class Node for default value
        
        # This is the multiple selection that is dsiable currently
        self.mode = mode        # selection mode TO BE IMPLMEENTED

        # offset of tree from upper left corner
        self.offx = offx
        self.offy = offy
        
        # functions to be called upon events
        self.actions = {'select':None, 'deselect':None,
                        'expand':None, 'collapse':None }

        self.ownsMaster = True  # when set to True, destroying the tree will
                                # destroy the Tree's master

        if master is None:
            master = Toplevel()
        else:
            self.ownsMaster = False
            
        self.master = master

        # GUI related
        self.topFrame = Pmw.PanedWidget(
            master, orient='horizontal', hull_relief='sunken',
            hull_width=width, hull_height=height,
            command=self.savePaneWidth_cb,
            )

        if hasattr(self.topFrame.component('hull').master, 'title'):
            self.topFrame.component('hull').master.title(name)

        if not nohistory:
            self.topFrame.bind("<Configure>", self.configureTopFrame_cb, '+')
        
        self.historyWidth=historyWidth
        self.treeWidth=treeWidth
        if not nohistory:
            self.topFrame_2 = self.topFrame.add('history', 
                                                min=0, size=historyWidth)
        self.topFrame_1 = self.topFrame.add('tree', min=20, size=treeWidth)
        
        # create the Canvas
        if canvasHeaderCol:
            self.canvasHeaderCol = Tkinter.Canvas(self.topFrame_1, height=50)
            #self.canvasHeaderCol.grid(column=0, row=0, sticky='ew')
            self.canvasHeaderCol.pack(side='top', expand=0, fill='x')

            self.scrolledCanvas = Pmw.ScrolledCanvas(
                self.topFrame_1, vscrollmode='dynamic', hscrollmode='none',
                horizscrollbar_command=self.scrollAll, canvas_bg='white')
                #canvas_bg='#c3d0a6')
        else:
            self.scrolledCanvas = Pmw.ScrolledCanvas(
                self.topFrame_1, vscrollmode='dynamic', hscrollmode='dynamic',
                horizscrollbar_command=self.scrollAll, canvas_bg='white')
            #canvas_bg='#c3d0a6')
            
            self.canvasHeaderCol = None

        self.scrolledCanvas.pack(side='left', expand=1, fill='both')
        self.canvas = self.scrolledCanvas.interior()
        if os.name == 'nt': #sys.platform == 'win32':
            self.canvas.bind("<MouseWheel>", self.scrollUpDown)
        else:
            self.canvas.bind("<Button-4>", self.scrollUp)
            self.canvas.bind("<Button-5>", self.scrollDown)

        #self.canvas.configure(width=treeWidth, height=treeHeight,
        #    borderwidth=2)

        #if canvasHeaderCol:
            #self.sbx=Scrollbar(self.topFrame_1, orient="horizontal")
            #self.sbx.grid(column=nextCol, row=2, sticky='ew')
            #self.sbx.pack(side='bottom', expand=1, fill='x')
            #self.canvasHeaderCol['xscrollcommand']=self.sb.set

        join = path.join
        plus_icon  = Image.open(join(ICONDIR, "plusnode.gif"))
        self.plus_icon  = ImageTk.PhotoImage( master=master,image=plus_icon)
        minus_icon  = Image.open(join(ICONDIR, "minusnode.gif"))
        self.minus_icon  = ImageTk.PhotoImage( master=master,image=minus_icon)

        # create the frame for "History" Label and the two hide/show labels
        #frame = Tkinter.Frame(master=self.topFrame_1, )
        #self.showHideHistoryIcon=Label(frame, text="Show History")
        #self.showHideHistoryIcon.pack(side="left", anchor='nw')
        #frame.pack(side='top')
        #self.showHideHistoryIcon.bind('<1>', self.toggleHistoryVisibility)

        if not nohistory:
            frame = Tkinter.Frame(master=self.topFrame_2, )
            Label(frame, text = "History").pack(side="left")
            frame.pack(side='top')
        
            # Add HistoryList box
            self.historyList = ListView(self, master=self.topFrame_2, 
                                        width=self.historyWidth)
        
        self.show()
        self.topFrame.update()

        self.currentPanesWidth = [historyWidth, treeWidth]

        # hide history is that was specified:
        if not historyVisible and not nohistory:
            self.currentPanesWidth[0] = 0.
            self.topFrame.configurepane('history', size=0.)
        #self.historyShown=True
        #self.toggleHistoryVisibility()

        #self.canvas.bind_all("<Key>", self.handle_key)
        # class level binding . P105 Grayson TKinter book
        self.canvas.bind_class('Canvas', "<Up>",    self.handle_key)
        self.canvas.bind_class('Canvas', "<Down>",  self.handle_key)
        self.canvas.bind_class('Canvas', "<Left>",  self.handle_key)
        self.canvas.bind_class('Canvas', "<Right>", self.handle_key)
        self.canvas.bind('<Enter>', self.Enter_cb)

        self.double1PickNode_Func_cb = None # call back for left double click on
                                            # tree node label
        self.double2PickNode_Func_cb = None # call back for left double click
        self.double3PickNode_Func_cb = None # call back for left double click
        

    def double1PickNode_Func(self, node):
        if self.double1PickNode_Func_cb:
            self.double1PickNode_Func_cb(node)

    def double2PickNode_Func(self, node):
        if self.double2PickNode_Func_cb:
            self.double2PickNode_Func_cb(node)
        

    def double3PickNode_Func(self, node):
        if self.double3PickNode_Func_cb:
            self.double3PickNode_Func_cb(node)

    def scrollAll(self, *args):        
        #print "scrollAll", args
        apply( self.canvas.xview, args)
        canvas = self.canvasHeaderCol
        if canvas:
            region = self.canvas.configure('scrollregion')[4]
            #print region, self.canvasHeaderCol.configure('scrollregion')[4]
            #print self.canvas.configure('width')
            #print canvas.configure('width')
            canvas.configure(scrollregion=region)
            apply( canvas.xview, args)


    def scrollUp(self, event):        
        self.canvas.yview_scroll(-1, "units")


    def scrollDown(self, event):        
        self.canvas.yview_scroll(1, "units")


    def scrollUpDown(self, event):
        if event.delta < 0:
            self.scrollDown(event)
        else:
            self.scrollUp(event)


    def handle_key(self, event):
        """handle the keyboard operations..
Up / Down / Left / Right arrow keys are supported"""
        # widget-wide key dispatcher

##         atFocus = self.canvas.focus()
##         if not atFocus:
##             return

        # navigation
        if self.multi_choice: return
        if event.keysym == "Up":
            self.moveSelectionUp()
        elif event.keysym == "Down":
            self.moveSelectionDown()
        elif event.keysym == "Right":
            sel = self.GetSelected()
            if sel:
                sel.Expand()
        elif event.keysym == "Left":
            sel = self.GetSelected()
            if sel:
                sel.Collapse()
        else:
            pass # print event.keysym

    def Enter_cb(self, event):
	"""Call back function trigger when the mouse enters the cavas """
        #print 'entering tree'
	self.canvas.focus_set()
        #atFocus = self.canvas.focus()

        
    def moveSelectionUp(self):
        """handle the 'up' key..  """
        selection = self.GetSelected()
        if selection.parent is None :
            idx = self.roots.index(selection)
            if idx ==0 :
                return
            else:
                self.Select( self.lastVisibleNodeOf( self.roots[ idx -1 ]))
                return
        children = selection.parent.children
        idx = children.index(selection)
        if idx is not 0:
            self.Select ( self.lastVisibleNodeOf(  children[ idx -1 ] ) )
            return
        else:
            self.Select(selection.parent)
            return

    def lastVisibleNodeOf(self, node):
        """ return the last visible node of 'node'  """
        if len(node.children ) == 0 or not node.expanded:
            return node
        return self.lastVisibleNodeOf(node.children[-1])
        

    def moveSelectionDown(self):
        """handle the 'down' key..  """
        selection = self.GetSelected()
        if selection.expanded and len(selection.children) is not 0:
            self.Select(selection.children[0])
            return

        self.Select (self.nextVisibleNodeOf( selection ) )
##         if selection.parent is None :
##             idx = self.roots.index(selection)
##             if idx ==  len(self.roots) -1: # last root
##                 pass
##             else:
##                 self.Select( self.roots[ idx +1 ])
##         else:
##             children = selection.parent.children
##             idx = children.index(selection)
##             if idx is len(children)-1:            
##                 self.Select ( self.nextVisibleNodeOf( selection ) )
##             else:
##                 self.Select(children[idx + 1 ])            

        return


######## fix me fix me

        


    def nextVisibleNodeOf(self, node):
        """ return the next visible node of 'node'  """
        if node.parent is None:
            idx = self.roots.index(node)
            if idx ==  len(self.roots) -1: # last root
                return node
            else:
                return self.roots[idx+1]
        else:
            children = node.parent.children
            idx = children.index(node)
            if idx is len(children)-1:
                return  self.nextVisibleNodeOf( node.parent ) 
            else:
                return children[idx + 1 ]

            
        if len(node.children ) == 0 or not node.expanded:
            return node
        return self.lastVisibleNodeOf(node.children[-1])



    def setAction(self, event, function):
        """define function to be called when an event occurs
current events are 'select', 'deselect', 'expand', 'collapse'
function has to be a callable and accept one argument which wil be the node
on which the event occurs.
Function can be None to not unregister the callback
"""
        assert event in ['select', 'deselect', 'expand', 'collapse']
        assert function is None or callable(function)
        self.actions[event] = function

        
    def destroy(self):
        """Destroy the frame, including TreeView and history list"""
        if self.ownsMaster:
            self.master.destroy()
        else:
            self.topFrame.destroy()


    def hide(self):
        """ Hide the frame, including TreeView and history list """
        self.topFrame.forget()

    def forget(self):
        self.hide()
    
    def show(self):
        """ Display the frame, including TreeView and history list """
        self.topFrame.pack(expand=1, fill='both')


    def findNodeFromName(self, name):
        """walks the tree and find the node matching the name
        name is the full name of the node"""
        path = name.split('|')
        for root in self.roots:
            if root.name==path[0]:
                newName = name[len(path[0])+1:]
                if len(newName):
                    n = root.findChildByName(root, newName)
                    if n:  return n
                else:
                    return root
                return None

    def draw_new_root(self, newRoot):
        """ this function draw the canvas representation of a root node 
        when a new root node is added """
        root = newRoot
        canvasIDs= root.canvasIDs

        if len(canvasIDs): # the node has a representation
            for id in canvasIDs:
                self.canvas.delete(id)
        
        # draw the line between plus/minus icon and the folder icon        
        canvasIDs.append( 
                self.canvas.create_line(root.x,root.y,
                root.x + OFFSET, root.y,fill='gray') )
        # draw vertical line above node
        # if this is a new root and not the first root
        if root != self.roots[0]:
            index = self.roots.index(root)
            x1 = self.roots[index-1].x
            y1 = self.roots[index-1].y
            h_line = self.canvas.create_line(root.x, root.y,
                    x1,y1, fill='gray') 
            self.canvas.lower(h_line)
            canvasIDs.append(h_line)    
    
        root.draw_caption_icon()


    def addNodeSet(self, name, object=None, parent=None, mouseBinding={},\
                   hasChildren=False, firstExpand_cb=None, nodeClass=Node):
        """Add a set of children nodes to a node and returns handle of parenet
The node has a name, an optional object that will be associated with the node
in the tree.
parent      : the node to which children nodes are added.  Cane be a string, a node isntance or an object already added to the tree.
name        : list of children names
objects     : list of objects associated with the children nodes
hasChildren : list of True/False, whether the added child has children
mouseBinding: list of mouseBinding
firstExpand_cb : call_back fucntions
                 (not a list, always objectBrower.expandTreeNode_cb)
"""
        
        if (type(object) is not types.ListType) or \
           (type(name) is not types.ListType) or \
           (type(hasChildren) is not types.ListType):
            warn("List of children needed, non-list type found")
            return None
           
        if self.mouseBinding is not None:
            mouseBinding.update(self.mouseBinding)

        num = len(name)
        nodeList=[]
        for i in range(num):
            if self.mouseBinding is not None:
                mouseBinding.update(self.mouseBinding[i])
            node = nodeClass(name[i], object[i],  \
                  hasChildren=hasChildren[i], firstExpand_cb=firstExpand_cb)
            nodeList.append(node)
            node.tree = self
            try:
                hash(object[i])
                node.objectKey = object[i]
            except TypeError:
                node.objectKey = self.objIndex
                self.objIndex +=1

    ##         if type(object) is not types.InstanceType:
    ##             node.objectKey = self.objIndex
    ##             self.objIndex +=1
    ##         else:
    ##             node.objectKey = object

            if self.obj2Node:
                self.objToNode[node.objectKey] = node

            self.numberOfNodes += 1
            node.uniqueID = self.numberOfNodes
            node.tag = [str(node.uniqueID)]
        
            # if parent given as a string, find the Node obj of the parent
            if type(parent) is types.StringType:
                input=parent
                parent = self.findNodeFromName(parent)
                if parent is None:
                    node.parentFullname = None
                    warn( "error in addNode, check name of parent: "+ input)  
                    return
                else:
                    node.parentFullname = input
            elif self.objToNode.has_key(parent):
                parent = self.objToNode[parent]
            elif not isinstance(parent, Node) and parent is not None:
                raise RuntimeError('bad parent')

            # if parent is given as None,we have a new root node
            # The new root is added to the end(bottom) of the tree
            if parent is None:
                node.parentFullname = None
                h = 0
                for r in self.roots:
                    if r.name == name :
                        warn("The node with name"+name + "already exists")
                        return
                    h += r.height
                # calc the Y offset of current node
                node.y += h * OFFSET + self.offy  
                node.x += self.offx
                self.roots.append(node)
            else:
                assert isinstance(parent, Node)
                if parent.parentFullname != None:
                    node.parentFullname = parent.parentFullname + '|' + \
                                          parent.name
                else:
                    node.parentFullname = parent.name

            node.parent = parent
      
        if parent is not None:
            # check duplicated node
            # FIXME ... this is expensive
##             for c in parent.children:
##                 if c.name == node.name:
##                     print "The node with name", name, "already exists"
##                     return           

            for node in nodeList:
                node.x = parent.x + OFFSET
                parent.children.append(node)
            if parent.expanded:
                parent.increaseParentHeight(offset=num)
                parent.inserted = True
            self.updateY()
            if parent.inserted:
                parent.draw_new_insert(num=num, mode = 'batch')
                parent.inserted = False
            parent.draw()
        else:
            for i in range(num):
                self.draw_new_root(nodeList[i])
        
        bb = self.canvas.bbox(Tkinter.ALL)
        self.canvas.configure(scrollregion=(0, 0,bb[2]+OFFSET, bb[3]+OFFSET))
            
        return nodeList


    def addNode(self, name, object=None, parent=None, mouseBinding={},\
                hasChildren=False, firstExpand_cb=None, nodeClass=Node):
        """Add a node to the tree and returns a handle to it
The node has a name, an optional object that will be associated with the node
in the tree.
Parent can be a string specifying the full name of a node or an isntance of
a Node. If no parent is specified this node will become the next root node
"""
        # the '|' is not allowed as name of the node
        if name.find('|')!=-1:
            warn( "No '|' is allowed in node name ")
            return

        if self.mouseBinding is not None:
            mouseBinding.update(self.mouseBinding)

        node = nodeClass(name, object, mouseBinding=mouseBinding, \
                         hasChildren=hasChildren, firstExpand_cb=firstExpand_cb)

        node.tree = self
        try:
            hash(object)
            node.objectKey = object
        except TypeError:
            node.objectKey = self.objIndex
            self.objIndex +=1
            
##         if type(object) is not types.InstanceType:
##             node.objectKey = self.objIndex
##             self.objIndex +=1
##         else:
##             node.objectKey = object
        
        if self.obj2Node:
            self.objToNode[node.objectKey] = node
        
        self.numberOfNodes += 1
        node.uniqueID = self.numberOfNodes
        node.tag = [str(node.uniqueID)]
                
        # if parent is given as None,we have a new root node
        # The new root is added to the end(bottom) of the tree
        if parent is None:
            node.parentFullname = None
            h = 0
            for r in self.roots:
                if r.name == name :
                    warn( "The node with name"+ name + "already exists")
                    return
                h += r.height
            # calc the Y offset of current node
            node.y += h * OFFSET + self.offy  
            node.x += self.offx
            self.roots.append(node)
            self.draw_new_root(node)
            
        else:
            # if parent given as a string, find the Node obj of the parent
            if type(parent) is types.StringType:
                input=parent
                parent = self.findNodeFromName(parent)
                if parent is None:
                    node.parentFullname = None
                    warn( "error in addNode, check name of parent:"+ input)
                    return            
            elif self.objToNode.has_key(parent):
                parent = self.objToNode[parent]
            elif not isinstance(parent, Node):
                raise RuntimeError('bad parent')
            #else:
            #    # only Node type is accepted.
            #    assert isinstance(parent, Node)

            if parent.parentFullname != None:
                node.parentFullname = parent.parentFullname + '|' + parent.name
            else:
                node.parentFullname = parent.name
                
            node.parent = parent      
            # check duplicated node
            # FIXME ... this is expensive
##             for c in parent.children:
##                 if c.name == node.name:
##                     print "The node with name", name, "already exists"
##                     return           
            node.x = parent.x + OFFSET
            parent.children.append(node)
            if parent.expanded:
                parent.increaseParentHeight()
                parent.inserted = True
            self.updateY()
            if parent.inserted:
                parent.draw_new_insert()
                parent.inserted = False
            # FIXME erasing the parent is very expensif, we only need to
            # draw from node to end of children and move everything below
            # parent down
            parent.draw()            
        
        bb = self.canvas.bbox(Tkinter.ALL)
        self.canvas.configure(
                            scrollregion=(0, 0,bb[2]+OFFSET, bb[3]+OFFSET))
            
        return node

    
    def deleteNode_byName(self, name, parent=None):#, object=None):
        """Delete a node from TreeView. using the full name of the node
parent should be given as full name of parent node"""
        if name==None :
            warn( "Node cannot be None")
            return
        
        if type(parent) is types.StringType:
            input=parent
            parent=self.findNodeFromName(parent)
            if parent is None:
                warn( "error in Delete, check name of parent: " +input)
                return
        
        bFound=False
        if parent is None:
            for r in self.roots:
                if r.name==name:
                    r.Erase()
                    self.roots.remove(r)
                    bFound=True
                    break
        else:  
            assert isinstance(parent, Node)
            # check for duplicated node
            for c in parent.children:
                if c.name==name:
                    c.Erase()
                    parent.children.remove(c)
                    bFound=True
                    break
        if bFound==False:
            warn( "The node with name"+ name+ "is not found", stacklevel=2)
            return
        

    def deleteNode(self, node):
        """ Delete a node from the TreeView. 
        Node can be given as the Node object, or full name of the node"""

        if type(node) is types.StringType:
            n = self.findNodeFromName(node)
            if n is None:
                warn( node+" not found")
            else:
                node = n

        if node in self.list_selected:
            self.list_selected.remove(node)
            
        assert isinstance(node, Node)
        if not node.tree:
            warn(node.name +" is not part of tree")
            return

        fullname=node.GetFullName()
        if self.current_selected:
            current_select_name= self.current_selected.GetFullName()
            # if the selected node is (grand) child of the node to be deleted
            # their fullnames should share same prefix
            if current_select_name.find(fullname)==0: 
                self.current_selected=None

        node.Erase()
        self.cleanUpOnNodeDelete(node)
        parent = node.parent
        if parent is None: # root node
            self.roots.remove(node)
        else:
            parent.children.remove(node)
            if len(parent.children) == 0:
                parent.draw()
        
        if self.nohistory: return
        # remove all the related items in the History list
        list=self.historyList.ItemList
        index=0
        to_delete=[] # the list
        for item in list:
            if item.name.find(fullname)==0:   # All items begin with node.name
                item.locked=False              # unlock the item
                to_delete.append(index)
            index+=1

        self.historyList.Delete_multi(to_delete)
        del node #remove node from namespace
        return    


    def cleanUpOnNodeDelete(self, node):
        # recursively traverse sub tree and update tree to reflect deletion
        # of subtree
        for c in node.children:
            self.cleanUpOnNodeDelete(c)

        if node.object is not None and self.obj2Node:
            del self.objToNode[node.objectKey] # remove object from dict
        node.tree = None   # remove handle to tree
        del node.object # remove handle to object
        self.numberOfNodes -= 1


    def do_updateY(self, from_root):
        """Called by updateY(), 
        recursively update the y coordination of all the children"""
        
        h=1
        for node in from_root.children:
            node.y=node.parent.y + OFFSET * h
            if node.height==0: # if the node was deleted.. height=0 
                continue
            if node.expanded and len(node.children) >0:
                h +=node.height
                self.do_updateY(node)                
            else:
                h +=1           
            
        
    def updateY(self):    
        """Recursively update the y coordination all the root nodes"""
        index=0
        for r in self.roots:
            index=self.roots.index(r)
            # if this is the first root node in Canvas            
            if index>0:
                r.y =self.roots[index-1].y+OFFSET*(self.roots[index-1].height)
            index +=1    
            self.do_updateY(r)


    def GetSelected(self):
        """returns the selected Node object """
        return self.current_selected
        

    def ExpandNode(self, node):
        """ Expand a given node in TreeView. used only by commandline 
        node can be full name or a handel to Node object"""
        if Node == None:
            return
        self.Select(node)
        tmp=self.GetSelected()
        if tmp is not None:
            tmp.Expand()
            self.deSelect()
            
        
    def CollapseNode(self, node):
        """ Collapse a given node in TreeView. used only by commandline 
        node can be full name or a handel to Node object"""
        if Node == None:
            return
        self.Select(node)
        tmp=self.GetSelected()
        if tmp is not None:
            tmp.Collapse()
            self.deSelect()
            
    def clearSelection(self):
        if not self.multi_choice:
            return
        canvas = self.canvas
        for node in self.list_selected:
            node.selected = False
            canvas.delete(node.selectboxID)
        self.list_selected = []

    def expandPath(self, node):
        while node and node.parent:
            if not node.parent.expanded:
                node.parent.Expand()
            node = node.parent
            self.expandPath(node.parent)

    def selectNodes(self, nodeList):
        # given a list of objects or node select them if the tree is in
        # multi_choice mode
        if not self.multi_choice:
            return
        firstNode = None
        self.selectionHistory.append(self.list_selected)
        for node in nodeList:
            if isinstance(node, Node):
                if firstNode==None:
                    firstNode = node
                if not node.selected:
                    if len(node.children): #prevent last elevel expansion
                        self.expandPath(node)
                    if node.visible:
                        node.selected= True
                        self.list_selected.append(node)
                        node.draw()
            elif self.objToNode.has_key(node):
                node = self.objToNode[node]
                if firstNode==None:
                    firstNode = node
                if not node.selected:
                    if len(node.children): #prevent last elevel expansion
                        self.expandPath(node)
                    if node.visible:
                        node.selected= True
                        self.list_selected.append(node)
                        node.draw()
            else:
                raise ValueError, 'Node %s not found'%str(node)

        if firstNode is not None:
            self.showNode(firstNode)
        

    def undoSelect(self, event=None):
        canvas = self.canvas
        for n in self.list_selected:
            n.selected = False
            canvas.delete(n.selectboxID)

        if len(self.selectionHistory):
            newSelection = self.selectionHistory.pop()
            self.list_selected = newSelection
            for n in newSelection:
                n.selected = True
                n.drawSelectionBox()


    def Select(self, node=None):
        """ Select a node, given as a string (full name), or as a Node object.
        """
        if self.multi_choice:
            if self not in self.list_selected:
                self.selectionHistory.append(self.list_selected)
                self.list_selected.append(self)
                
        prev_selection=self.GetSelected()
   
        result = False
        if type(node) is types.StringType:
            fullname=node
        else:
            assert isinstance(node, Node)
            if node.parent:
                fullname = node.parentFullname + "|"+node.name
            else:
                fullname = node.name
                
        path = fullname.split('|')        
        for root in self.roots:
            if root.name == path[0]:
                path.remove(path[0])
                if len(path)==0:
                    root.selected = True
                    root.draw()
                    result = True
                else:
                    if not root.expanded:
                        root.invoke()
                    result = root.DoSelect(root,path)
                break

        # make sure the selected is visible in canvas
        n = self.GetSelected()
        if result and n==prev_selection:
            return
            #warn( "The node" + n.name + " was already selected")
        
        self.showNode(n)

        if not result:
            warn( "Selection failed. Check the input")
            return
        else:
            if self.actions['select'] is not None:
                self.actions['select'](self.current_selected)
        

    def showNode(self, n):
        """move the tree to make sure the node n is visible
"""
        y = n.y - OFFSET + 0.0
        bb = self.canvas.bbox(Tkinter.ALL)
        self.canvas.yview_moveto(y/bb[3])
        
        
    def deSelect(self):
        if self.current_selected:
            self.current_selected.selected = False
            self.current_selected.draw()
            if self.actions['deselect'] is not None:
                self.actions['deselect'](self.current_selected)
        
        self.current_selected = None
        return


    def ButtonAdd_cb(self, event = None):
        """Callback function, when the Button is clicked"""
        self.AddToHistoryList()

        
    def AddToHistoryList(self):
        """Add current selelcted node to the history list"""
        if self.nohistory: return
        current=self.current_selected
        if current==None:
            return
        if current.parentFullname is not None:
            fullname=current.parentFullname + "|"+current.name
        else:
            fullname=current.name
        
        self.historyList.Insert(fullname, current)
        return


    def DeleteFromHistoryList(self, index=0):
        self.historyList.Delete(index)

        
    def LockHistoryItem(self, index):
        """Lock the n-th ListItem object in the history list 
           n is specified as index, range [0, len_of_list-1]
           After locking, the item cannot be deleted 
        """
        self.historyList.SaveItem(index)


    def UnLockHistoryItem(self, index):
        """Unlock the n-th ListItem object in the history list 
           n is specified as index, range [0, len_of_list-1]
           After unlocking, the item can be deleted 
        """
        self.historyList.UnSaveItem(index)


    def savePaneWidth_cb(self, sizes):
        # only way I found to keep track of pane sizes
        self.currentPanesWidth = sizes

        
    def configureTopFrame_cb(self, event=None):
        # prevent history from growing when width changes
        w = self.currentPanesWidth[0]
        # width = 0 makes pane full width while width=0. make it disappear
        # savePaneWidth_cb saved 0 when panes is configured with 0.
        if w==0:
            w = 0.
        self.topFrame.configurepane('history', size=w)


    def copy(self,newtree):
        """ copy the tree root  to another specify NewTree.
        Also allow to update the new tree.
        """
        assert isinstance(newtree, TreeView)
        for root in self.roots:
            self.copyNode(newtree,root,root.children)
            
    def copyNode(self,newtree,node,children,parent=None):
        """ add a node from tree to a new tree
        Function use recursively to go through all the children
        of a node
        """
        # test if node already in newtree
        new_node = newtree.findNodeFromName(node.GetFullName())
        if not new_node:
            new_node = newtree.addNode(node.name,object=node.object,
                                       parent=parent)
            
        for child in children:
            self.copyNode(newtree,child,child.children,parent=new_node)
        


    def deleteAllChildren(self,node):
        """ delete all the children of a node """
        # we make a deep copy of the children list
        # so we do not work on a truncated one after
        # a node has been deleted
        children = node.children[:]
        for child in children:
            self.deleteNode(child)

##     def toggleHistoryVisibility(self, event=None):
##         """Callback function. Show or hide the HistoryList
## """
##         if self.historyShown:
##             self.historyShown=False
##             self.showHideHistoryIcon.configure(text="Hide History")
##             self.oldHistSize = self.currentPanesWidth[0]
##             self.topFrame.configurepane('history', size=0.)
##         else:
##             self.historyShown=True
##             self.showHideHistoryIcon.configure(text="Show History")
##             self.topFrame.configurepane('history', size=self.oldHistSize)


    def moveNode(self,nodeFrom, newParent=None):
        """ this function change the parent of node 'nodeFrom' to be newParent. The whole sub-tree of nodeFrom should be appended as the child of newParent node.

nodeFrom : a Node object or fullName of Node
newParent: a Node object or fullName of Node, if None, add as new root

    """
        if nodeFrom == None: return
        # if newParent is None: add as new root node
        
        if type(newParent) is types.StringType:
            input=newParent
            newParent = self.findNodeFromName(newParent)
            if newParent is None:
                warn( "error in moveNode, check name of parent: "+ input  )
                return

        if type(nodeFrom) is types.StringType:
            input=nodeFrom
            nodeFrom = self.findNodeFromName(nodeFrom)
            if nodeFrom is None:
                warn("error in moveNode, check name of node to move: "+input)
                return

        # now the nodeFrom and newParent both are Node object.
        if nodeFrom.parent == newParent: return  # same parent?
        node = nodeFrom
        parent = newParent
        oldParent=node.parent
        
        node._collapseAll()
            
        if parent is not None:
            if parent.expanded and parent.isShown():
                parent.invoke()            
            self.deleteNode_byName(nodeFrom.name, nodeFrom.parentFullname)
            parent.children.append(node)
            node.parent=parent
            node.parentFullname = parent.GetFullName()

            node.height=1
            node.expanded =False
            node.updateHeight(parent)
            node.updateHeight(oldParent)
            
            node.updateX()            
            root = node.getRootNode()
            self.do_updateY(root)

            if parent.isShown():
                parent.draw_new_insert()

            if oldParent is not None and oldParent.isShown():
                oldParent.draw()            
            
        else:
            self.deleteNode_byName(nodeFrom.name, nodeFrom.parentFullname)
            self.roots.append(node)
            node.parent=None
            node.parentFullname = None            

            node.height=1
            node.expanded=False
            node.updateHeight(oldParent)
            node.updateX()
            self.updateY()
            self.draw_new_root(node)

        bb = self.canvas.bbox(Tkinter.ALL)
        self.canvas.configure(
                            scrollregion=(0, 0,bb[2]+OFFSET, bb[3]+OFFSET))


        return
    

    def expandAllNodes(self):
        """ expanding all the tree nodes
        """
        canvas = self.tree.canvas
        canvas.configure(cursor='watch')
        canvas.update_idletasks()
        for r in self.roots:
            if not r.expanded:
                r.invoke()
            r._expandAll()
        canvas.configure(cursor='')


    def collapseAllNodes(self):
        """only showing the root nodes"""
        for r in self.roots:
            if r.expanded:
                r.invoke()
            #r._collapseAll()
            # Notice: this function is never called..
            # collapse the root nodes will do the job

    def displayValueInTree(self, display=False):
        if display == self.displayValue: return
        self.displayValue = display
        for r in self.roots:
            r.invoke()

## end of file
