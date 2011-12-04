########################################################################
##
## Author Michel Sanner  Copyright TSRI (C) 2006 
##
########################################################################
#
# $Header: /opt/cvs/python/packages/share1.5/mglutil/gui/BasicWidgets/Tk/trees/tree.py,v 1.28 2008/10/02 17:10:14 sargis Exp $
#
# $Id: tree.py,v 1.28 2008/10/02 17:10:14 sargis Exp $
#
"""Efficient TreeWidget based on Tkinter and Pmw

The speed comes from only drawing visibles parts of the tree.
"""

import os, sys, types, Pmw
from PIL import Image, ImageTk
from weakref import ref
from Pmw import ScrolledCanvas
from Tkinter import Tk, PhotoImage, Label, ALL, Menu, IntVar
from mglutil.util.packageFilePath import findFilePath
from mglutil.gui import widgetsOnBackWindowsCanGrabFocus

## TODO
## - add size subdirectory to IconManager
## - cannot select 2 ranges, might want to add drag selection

##
## BUGS
##

class IconsManager:
    """The IconsManager object simplifies creating PhotoImage icons from
images in the directory.  The directory can be set by specifying a path,
a path within a Python package, or a Pyton package.
Images are loaded upon requests and a reference to the icon is saved to
avoid its garbage collection."""


    def __init__(self, path=[], packageName=''):
        """Create an IconManager object.
obj <- IconManager(path='', packageName='')
If packageName is specified, the path is assumed to be relative to the location
of the package."""
        if packageName:
            mod = __import__(packageName)
            components = packageName.split('.')
            for comp in components[1:]:
                mod = getattr(mod, comp)
            packagePath = mod.__path__[0]
        else:
            packagePath=''
            
        self.directory = os.path.join( packagePath, os.path.join(*path))
        
        assert os.path.exists(self.directory)
        assert os.path.isdir(self.directory)
        self.icons = {}


    def get(self, iconName, master=None):
        """ return an Tk PhotoImage for a given icon name and saves a reference
"""
        icon = self.icons.get(iconName, None)
        if icon is None:
            filename = os.path.join(self.directory, iconName)
            image = Image.open(filename)
            icon = ImageTk.PhotoImage(master=master, image=image)
            self.icons[iconName] = icon
        return icon



class KeySelectable:
    """Adds the ability to use keystrokes to quickly select items in a list.
root has to be a widget supporting .bind .after
"""
    
    def __init__(self, KeyRootTk):
        self.KeyRootTk = KeyRootTk
        self.afterID = None
        self.matchString = ''
        self.lastMatchString = ''
	KeyRootTk.bind('<KeyPress>', self.key_cb)
	KeyRootTk.bind('<KeyRelease>', self.keyUp_cb)
        self.isControl = False
        self.isShift = False
        self.isAlt = False
        self.ctrlModCallback = None
        self.shiftModCallback = None
        self.altModCallback = None
        

    def timeOut(self, event=None):
        """resets self.matchCharIndex to 0, called after a short period of
        time if no new character has been typed"""
        #print 'timeout'
        self.lastMatchString = self.matchString
        self.matchString = ''
        self.matchCharIndex = 1
        self.afterID = None
        

    def keyUp_cb(self, event=None):
        if event.keysym=='Control_L' or event.keysym=='Control_R':
            self.isControl = False
        elif event.keysym=='Shift_L' or event.keysym=='Shift_R':
            self.isShift = False
        elif event.keysym=='Alt_L' or event.keysym=='Alt_R':
            self.isAlt = False

        
    def key_cb(self, event=None):
        # use key strokes to select entry in listbox
        # strokes placed within 500 miliseconds are concatenated
        #print self.matchCharIndex, '|', self.matchString, '|', event.keysym
        if event.keysym=='Control_L' or event.keysym=='Control_R':
            self.isControl = True
            return
        elif event.keysym=='Shift_L' or event.keysym=='Shift_R':
            self.isShift = True
            return
        elif event.keysym=='Alt_L' or event.keysym=='Alt_R':
            self.isAlt = True
            return

        if self.isControl:
            if self.ctrlModCallback:
                self.ctrlModCallback(event)
            return
        elif self.isShift:
            if self.shiftModCallback:
                self.shiftModCallback(event)
            return
        elif self.isAlt:
            if self.altModCallback:
                self.altModCallback(event)
            return
            
        if event.keysym=='Return':
            str = self.lastMatchString
        else:
            str = self.matchString + event.keysym
        #print str
        item = self.match(str)
        if item:
            self.selectItem(item)
            if self.afterID is not None:
                self.KeyRootTk.after_cancel(self.afterID)
            self.afterID = self.KeyRootTk.after(1000, self.timeOut)
            self.matchString = str

    # SUBCLASS THIS
    def match(self, name):
        """has to return None if no match or an object that matches"""
        print 'jumping to ', name
        return None


    def selectItem(self, item):
        """do what has to be done to show what matches the typed string"""
        print 'selecting item', item



class Tree(ScrolledCanvas, KeySelectable):
    """Tree widget for Tk and Pmw"""


    def __init__(self, master, root, iconsManager=None, selectionMode='single',
                 idleRedraw=True, nodeHeight=20, **kw):
        """Tree( master, root, idleRedraw=True, kw = {}, **opts)
- Master can be a Tk frame in which the tree is displayed
- root has to be a Node object
- iconsManageer has to be an instance of an IconManager object or None.
If None is passed, an IconManager with the default icons directory will be
created.  If an IconManager is passed we expect to find ...
- selection mode cane be 'single' or 'multiple'. In multiple mode, the shift
modifier is use to select/deselect ranges (ranges can only be defined on
sibling nodes) and the Control modifier is used to add/remove to the current
selection.  Node are selected by click withthe left mouse button on the node's
label or icon.
- Set idleRedraw to False to disable background drawing.
- nodeHeight is the heigh of a node int he Tree in pixels
- kw can contain any keywork allowable for a Pmw.ScrolledCanvas
"""
        assert isinstance(root, Node)
        assert selectionMode in ['single', 'multiple']

        self.master = master
        
        self.selectionMode = selectionMode
        self.lastPickedNode = None
        
        self.idleRedraw  = idleRedraw  # when set to true redraw operation
                                       # occur when CPU is idle
        
        if iconsManager is None:
            iconsManager = IconsManager(
                ['gui','BasicWidgets','Tk','trees','Icons']
                , 'mglutil')

        assert isinstance(iconsManager, IconsManager)
        self.iconsManager = iconsManager

        # cache the icons use to expand and collapse since we will
        # use them for each node
        self.collapsedIcon  = self.iconsManager.get("1rightarrow.png")
        self.expandedIcon = self.iconsManager.get("1downarrow.png")
        self.iconHalfWidth = self.collapsedIcon.width()/2 
        self.selectedNodes = []
        self.selectionHistory = [] # used to undo selection operations
        self.nodeHeight = nodeHeight
        self.headerHeight = 30
        
        self.fisrtVisibleNodeLineNum = 0 # height of the first visible node
                                # i.e. how many nodes drown above, 0 for root

        self.firstVisibleNode = root # first node visible in window
        
        self.displayedNodes = [] # list of nodes having a graphical represenation
        self.pending_draw = None
        self.suspendRedraw = False

        self.root = root
        root.currenty = 0
        root.tree = ref(self)

        root.childNumber = 0
        self.nbLines = 1 # height of the tree (i.e. how many lines in canvas)

        # build the GUI
        if not kw.has_key('vscrollmode'):
            kw['vscrollmode'] ='dynamic'
        if not kw.has_key('hscrollmode'):
            kw['hscrollmode'] ='dynamic'
        if not kw.has_key('hull_width'):
            kw['hull_width'] = 550
        if not kw.has_key('hull_height'):
            kw['hull_height'] = 100
        if not kw.has_key('background') and not kw.has_key('bg') and \
               not kw.has_key('canvas_bg'):
            kw['canvas_bg']='white'
        if not kw.has_key('usehullsize'):
            kw['usehullsize'] = 1
            
        self.nbLinesPerPage = (kw['hull_height']-self.headerHeight) / nodeHeight

        if not kw.has_key('yscrollincrement'):
            kw['canvas_yscrollincrement'] = nodeHeight

        kw['horizscrollbar_command'] = self.xview
        kw['vertscrollbar_command'] = self.yview
        kw['borderframe'] = 5
        #kw['canvasmargin'] = 10
        #kw['hscrollmode'] = 'none'
        
        apply( ScrolledCanvas.__init__, (self, master), kw)

        canvas = self.canvas = self.component('canvas')

        KeySelectable.__init__(self, canvas)
        self.ctrlModCallback = self.handleControlKey

        canvas.bind("<Configure>", self.configure_cb)
        canvas.bind("<Key-Prior>", self.pageUp)
        canvas.bind("<Key-Next>", self.pageDown)
        canvas.bind("<Key-Up>", self.lineUp)
        canvas.bind("<Key-Down>", self.lineDown)
        if sys.platform == 'win32':
            canvas.bind("<MouseWheel>", self.lineUpDown)
        else:
            canvas.bind("<Button-4>", self.lineUp)
            canvas.bind("<Button-5>", self.lineDown)
        canvas.bind("<Enter>", self.enter_cb)

        self.isControl = False
        self.isShift = False
        self.isAlt = False

        if widgetsOnBackWindowsCanGrabFocus is False:
            lActiveWindow = canvas.focus_get()
            if    lActiveWindow is not None \
              and ( lActiveWindow.winfo_toplevel() != canvas.winfo_toplevel() ):
                return

        canvas.focus_set()


    def redrawHeader(self, *args):
        pass
    
    def yview(self, *args):
        ## callback for vertscrollbar_command
        # args can be ('scroll', number, type) where type is page or unit
        # or ('moveto', '0.2854982')
        #rint args
        if args[0] == "scroll":
            # ('scroll', '1', 'pages') when click on back of bar
            # ('scroll', '1', 'units') ehwn click on bar arrow
            self.yview_scroll(args[1], args[2])
        else:
            # ('moveto', '0.350939') when dragging bar
            percent = float(args[1])
            line = int(self.nbLines * float(percent))

            if self.fisrtVisibleNodeLineNum != line:
                self.scrollView(line - self.fisrtVisibleNodeLineNum)
                # compute percentage of motion
                total = self.nbLines+(self.headerHeight/float(self.nodeHeight))
                percent = self.fisrtVisibleNodeLineNum/total
                ScrolledCanvas.yview_moveto(self, percent)
                self.redraw()

                self.redrawHeader(args)
        

    def pageUp(self, event):
        self.yview_scroll(-1, "pages")
        return "break"


    def pageDown(self, event):
        self.yview_scroll(1, "pages")
        return "break"


    def lineUp(self, event):
        self.yview_scroll(-1, "units")
        return "break"


    def lineDown(self, event):
        self.yview_scroll(1, "units")
        return "break"


    def lineUpDown(self, event):
        if event.delta < 0:
            return self.lineDown(event)
        else:
            return self.lineUp(event)


    def yview_scroll(self, *args):
        # mouse wheel or arrow keys up and down args = (1, 'units'), (-1, 'units')
        # page up and down keys (1, 'pages')(13, 'units')
        #                       (-1, 'pages')(-13, 'units')
        # end key args = "jumping to  End"
        # home key args = "jumping to  Home"
        height = self.winfo_height() - self.headerHeight
        #height = int(self.canvas['height']) - self.headerHeight
        if self.nbLines * self.nodeHeight <= height:
            return

        if args[1] == "pages":
            # FIXME we scroll the whole thing one page
            linesPerPage = int(args[0]) * (height / self.nodeHeight)
            self.yview_scroll(linesPerPage, "units")
        else: # it is a unit scroll of args[0] units
            self.scrollView(int(args[0]))
            ScrolledCanvas.yview_scroll(self, *args)
            self.redraw()

        self.redrawHeader(args)


    def scrollView(self, deltalines):
        """move the visible part of the tree up of down by deltalines lines
"""
        if deltalines > 0:
            # compute the index of the line for which the bottom of the tree
            # is drawn at the bottom of the visible window
            last = self.nbLines - self.nbLinesPerPage
            # clamp deltalines so we do no go too far down
            deltalines = min(last-self.fisrtVisibleNodeLineNum, deltalines)
            for i in range(deltalines):
                node = self.firstVisibleNode.nextNode()
                if node is None:
                    break
                self.fisrtVisibleNodeLineNum += 1
                self.firstVisibleNode = node
        else:
            for i in range(-deltalines):
                node = self.firstVisibleNode.previousNode()
                if node is None:
                    break
                self.fisrtVisibleNodeLineNum -= 1
                self.firstVisibleNode = node


    def enter_cb(self, event=None):
        if widgetsOnBackWindowsCanGrabFocus is False:
            lActiveWindow = self.focus_get()
            if    lActiveWindow is not None \
              and ( lActiveWindow.winfo_toplevel() != self.winfo_toplevel() ):
                return

        self.canvas.focus_set()


    def configure_cb(self, event=None):
        
        nl = (self.winfo_height()-self.headerHeight) / self.nodeHeight
        self.nbLinesPerPage = nl
        last = self.nbLines - nl
        if self.fisrtVisibleNodeLineNum > last:
            self.scrollView(last - self.fisrtVisibleNodeLineNum)
        self.redraw()


    def destroy(self):
        if self.root:
            self.root.destroy()
        ScrolledCanvas.destroy(self)

    ##
    ##  Drawing
    ##
    def undisplay(self):
        for node in self.displayedNodes:
            node.deleteNodeIcons()

    def reparentNodes(self, tree):
        n = self.root
        while n:
            n.tree = ref(tree)
            n = n.nextNode()

    def redraw(self):
        """post a redraw or actually redraw depending on self.idleRedraaw
"""
        if self.suspendRedraw: return
        if self.root:
            if self.idleRedraw:
                if self.pending_draw:
                    self.after_cancel(self.pending_draw)
                self.pending_draw = self.after(10, self.reallyRedraw)
            else:
                self.reallyRedraw()


    def reallyRedraw(self):
        """actually redraw"""
        self.pending_draw = None

        # destroy representation of visible nodes
        for node in self.displayedNodes:
            #print 'destroying', node.object
            node.deleteNodeIcons()

        # build list of visible nodes
        nodeHeight = self.nodeHeight
        nodes = []
        node = self.firstVisibleNode
        for i in range( self.nbLinesPerPage + 1):
            nodes.append(node)
            node = node.nextNode()
            if node is None:
                break

##  this approach does not destoy some canvas items of node that move
## for instance after an expand

##         # destroy representation of nodes no longer visible
##         for node in self.displayedNodes:
##             if nodesd.has_key(node):
##                 print 'deleteIcon', node.object
##                 node.deleteNodeIcons()

        # Draw visible nodes
        ypos = self.fisrtVisibleNodeLineNum * nodeHeight + self.headerHeight
        for node in nodes:
            if node.currenty is None:
                node.currenty = 0
            #print 'redraw', node.needsRedraw, node.object
            node.redraw(ypos)
            ypos += nodeHeight
        
        self.displayedNodes = nodes
      
        # Update the scroll-bar size
        if self.root:
            canvas = self.canvas
            x1, y1, x2, y2 = canvas.bbox(ALL)
            canvas.configure( scrollregion=(0, 0, x2, self.nbLines*nodeHeight
                                            +self.headerHeight))


    def updateTreeHeight(self):
        """ """
        self.nbLines = self.root.countSubtreeLines()


    def clearSelection(self, history=True):
        """Deselect all selected Nodes."""
        canvas = self.canvas
        if history:
            self.selectionHistory.append(self.selectedNodes[:])
        selectedNodes = self.selectedNodes[:]
        for node in selectedNodes:
            node.deselect(history=False)


    def undoSelect(self):
        """ """
        if len(self.selectionHistory):
            self.clearSelection(history=False)
            self.selectedNodes = self.selectionHistory.pop()
            for n in self.selectedNodes:
                n.select(only=False, history=False)


    def handleControlKey(self, event):
        """ """
        if event.keysym in ['z', 'Z']:
            self.undoSelect()

##     def clearTree(self):
##         if self.pending_draw:
##             self.after_cancel(self.pending_draw)
##             self.pending_draw = None
      
##         if self.root:
##             self.root.destroy()
##             for n in self.displayedNodes:
##                 n.deleteNodeIcons()


##     def refresh(self):
##         """Refresf the tree. Notice that it is speeder to update only some Nodes (with Node.update() or Node.updatetree() and then tree.redraw()), if you know which Nodes have changed."""
##         if self.root:
##             self.root.refresh()
##             self.redraw()



class Node:
    """Base class for Nodes in a Tree"""

    def __init__(self, object, parent):
        """Create a Node.
object is the Python object represented by this node in the tree.
Object is expected to provide an interable sequence in its .children attribute.
Parent can be either the parent's Node, or the tree (for the root Node).
"""
        self.object = object
        self.isExpanded = False
        self.isSelected = False
        self.hasBeenExpanded = False
        self.children = []
        self.childNumber = 0 # index of this node in the list of children of its parent
        
        self.iconWidth = 0

        self.currenty = None     # y value on canvas when node is visible
                                 # if None the node is not drawn
        self.needsRedraw = False
        
        self.parent = parent
        if parent is None:
            self.generation = 0   # how many ancestors
            self.tree = None      # weakref to tree
        else:
            self.generation = parent.generation + 1
            self.tree = parent.tree

        # canvas ids
        self.labelTkid = 0  # canvas id of node name
        self.iconTkid = 0   # canvas id if node's icon
        self.expandIconTkid = 0 # canvas id of expand/collapse icon
        self.selectionBoxId = 0 
        self.canvasIDs = [] # list of other ids
        
    ##
    ## to be overridden
    ##
    def isExpandable(self):
        """Returns true if this node has children"""
        return 0


    def createChildrenNodes(self):
        """Create node for all children of self.object"""
        children = []
        klass = self.__class__
        for child in self.object.children:
            children.append(klass(child, self))
        return children


    ##
    ##  recursive traversals
    ##
    def countSubtreeLines(self):
        """Recusrsively traverse the tree and count the number of lines needed
to draw this subtree on a canvas"""
        nbLines = 1
        if self.isExpanded:
            for child in self.children:
                nbLines = nbLines + child.countSubtreeLines()
        return nbLines


    def lastVisibleChild(self):
        """Return the last visible child of this node"""
        if self.hasBeenExpanded and self.isExpanded:
            return self.children[-1].lastVisibleChild()
        else:
            return self


    def findNextNode(self, node):
        """recursively traverse the tree to find the first node visible
below self"""
        # if this is not the last child, return the next one
        if node.childNumber + 1 < len(self.children):
            return self.children[node.childNumber + 1]

        if self.parent: # we need to walk the tree to find te next node
            return self.parent.findNextNode(self)
        else: # root node
            return None


    ##
    ## tree navigation methods
    ##
    def previousNode(self):
        """return the node right above this one in the tree"""
        if self.parent is None: # root node
            return None
        if self.childNumber==0: # node is first child, return parent
            return self.parent
        else: # return last visible child of sibling before self
            return self.parent.children[self.childNumber-1].lastVisibleChild()


    def nextNode(self):
        """return the node right below this one in the tree"""
        if self.isExpanded and self.children:
            return self.children[0]
        else:
            if self.parent:
                return self.parent.findNextNode(self)
            else:
                return None # Root Node


    ##
    ##  expanding and collapsing the tree
    ##
    def expand(self, event = None):
        """Expand this node and show its children"""
        if self.isExpanded:
            return

        tree = self.tree()
        if not self.hasBeenExpanded:
            self.children = self.createChildrenNodes()
            for i,c in enumerate(self.children):
                c.childNumber = i
                tree.objectToNode[c.object] = c
            self.hasBeenExpanded = True

        if len(self.children):
            self.needsRedraw = True
            self.isExpanded = True
            tree.updateTreeHeight()
            tree.redraw()


    def collapse(self, event = None):
        """Collapse this node, i.e hid its children"""
        if not self.isExpanded:
            return
        
        self.isExpanded = False
        self.needsRedraw = True
        tree = self.tree()
        tree.updateTreeHeight()
        tree.redraw()


    def toggleExpansion(self, event=None):
        """Toggles expanded/collapsed state of this node"""
        if self.isExpanded:
            self.collapse()
        else:
            self.expand()


    ##
    ##  Drawing
    ##
    def getIcon(self):
        """return node's icons"""
        iconsManager = self.tree().iconsManager
        if self.isExpandable():
            if self.isExpanded:
                icon = iconsManager.get("expandedIcon.pgm")
            else:
                icon = iconsManager.get("collapsedIcon.pgm")
        else:
            icon = iconsManager.get("python.pgm")

        if icon:
            self.iconWidth = icon.width()
        else:
            self.iconWidth = 0
        return icon


    def drawExpandCollapseIcon(self, x, y):
        """Draw the icon used to expand and collapse nodes"""
        tree = self.tree()
        canvas = tree.canvas

        # delete old icon
        if self.expandIconTkid:
            canvas.delete(self.expandIconTkid)

        # draw new one
        if self.isExpandable():
            if self.isExpanded:
                icon = tree.expandedIcon
            else:
                icon = tree.collapsedIcon

            tkid = self.expandIconTkid = canvas.create_image(x, y, image=icon)
            canvas.tag_bind(tkid, "<Button-1>", self.toggleExpansion)
            iconWidth = icon.width()
        else:
            iconWidth = 0
            self.expandIconTkid = None

        return iconWidth


    def drawNodeIcon(self, x, y):
        """Draw the node's icon"""
        tree = self.tree()
        canvas = tree.canvas

        if self.iconTkid:
            canvas.delete(self.iconTkid)

        icon = self.getIcon()
        if icon:
            self.iconTkid = canvas.create_image(x, y, image=icon, anchor="nw")
            canvas.tag_bind(self.iconTkid, "<Double-1>", self.toggleExpansion)
            canvas.tag_bind(self.iconTkid, "<Shift-Button-1>",
                            self.selectRange_cb)
            canvas.tag_bind(self.iconTkid, "<Control-Button-1>",
                            self.modifySelection_cb)
            canvas.tag_bind(self.iconTkid, "<Button-1>", self.toggleSelection)
            if sys.platform != 'win32':
                canvas.tag_bind(self.iconTkid, "<Button-4>", self.tree().lineUp)
                canvas.tag_bind(self.iconTkid, "<Button-5>", self.tree().lineDown)
        else:
            self.iconTkid = None

        return self.iconWidth


    def drawNodeLabel(self, x, y):
        """Draw the node's label"""
        tree = self.tree()
        canvas = tree.canvas

        if self.labelTkid:
            canvas.delete(self.labelTkid)
        try:
            label = self.label
        except AttributeError:
            text = repr(self)
            if not text:
                self.needsRedraw = True
                return 0
            balloon = None
            if len(text) > 18:
                fullTxt = text
                text = text[:15]+"..."
                balloon = Pmw.Balloon(canvas)
                
            label = self.label = Label(canvas, text=text,
                                       bd=0, padx=2, pady=2, bg=canvas["bg"])
            if balloon:
                balloon.bind(self.label, fullTxt)
            label.bind("<Button-1>", self.toggleSelection)
            label.bind("<Shift-Button-1>", self.selectRange_cb)
            label.bind("<Control-Button-1>", self.modifySelection_cb)
            label.bind("<Double-1>", self.toggleExpansion)
            label.bind("<Button-3>", self.button3OnLabel)
            if sys.platform == 'win32':
                 label.bind("<MouseWheel>", self.tree().lineUpDown)
            else:
                label.bind("<Button-4>", self.tree().lineUp)
                label.bind("<Button-5>", self.tree().lineDown)
        self.labelTkid = canvas.create_window(x, y, anchor='nw',
                                              window=self.label)

        bb = canvas.bbox(self.labelTkid)
        return bb[2]-bb[0]

    
    def drawNodeCustomization(self, x, y):
        """Draw additional things on the canvas"""
        return 0
  

    def redraw(self, y):
        """Redraw this node at position y on the canvas"""
        if self.currenty is None:
            return

        tree = self.tree()
        if y != self.currenty or self.needsRedraw:
            x = self.generation * tree.nodeHeight
            x = x + tree.iconHalfWidth
            x += self.drawExpandCollapseIcon(x, y+tree.iconHalfWidth)
            x += self.drawNodeIcon(x, y)
            x += self.drawNodeLabel(x, y)
            x += self.drawNodeCustomization(x, y)
            self.currenty = y
            if self.isSelected:
                if self.selectionBoxId:
                    tree.canvas.delete(self.selectionBoxId)
                self.drawSelectionBox()
        self.needsRedraw = False


    def deleteNodeIcons(self):
        """Delect canvas items representing this node"""
        canvas = self.tree().canvas

        if self.iconTkid:
            canvas.delete(self.iconTkid)
            self.iconTkid = None

        if self.expandIconTkid:
            canvas.delete(self.expandIconTkid)
            self.expandIconTkid = None

        if self.labelTkid:
            canvas.delete(self.labelTkid)
            self.labelTkid = None

        try:
            canvas.delete(self.label)
            del self.label
        except AttributeError:
            pass
        
        if self.selectionBoxId:
            canvas.delete(self.selectionBoxId)
            self.selectionBoxId = None
            self.labelTkid = None

        for id in self.canvasIDs:
            canvas.delete(id)
        self.canvasIDs = []

        self.currenty = None


    def destroy(self):
        for child in self.children:
            child.destroy()


    ##
    ## selection
    ##
    def drawSelectionBox(self):
        """draw a yellow box"""
        tree = self.tree()
        canvas = tree.canvas
        w = tree.winfo_width()
        y = self.currenty
        id = canvas.create_rectangle(0, y-2, w, y+20, outline='yellow',
                                     fill='yellow') 
        self.label.configure(bg='yellow')
        canvas.lower(id)
        self.selectionBoxId = id


    def button3OnLabel(self, event=None):
        # override in subclass
        return

    
    def toggleSelection(self, only=True, event=None):
        tree = self.tree()
        tree.lastPickedNode = self
        if self.isSelected:
            self.deselect()
        else:
            self.select(only=only)


    def select(self, only=True, history=True):
        if self.isSelected:
            return

        tree = self.tree()
        if history:
            tree.selectionHistory.append(tree.selectedNodes[:])
        if only:
            tree.clearSelection(history=False)
        self.isSelected = True
        tree.selectedNodes.append(self)
        if self in tree.displayedNodes:
            self.drawSelectionBox()


    def deselect(self, history=True):
        if not self.isSelected:
            return
        tree = self.tree()
        if history:
            tree.selectionHistory.append(tree.selectedNodes[:])
        tree.selectedNodes.remove(self)
        if self in tree.displayedNodes:
            tree.canvas.delete(self.selectionBoxId)
            self.label.configure(bg='white')
        self.isSelected = False
        self.selectionBoxId = None


    def selectRange_cb(self, event=None):
        tree = self.tree()
        if tree.selectionMode=='single':
            return

        if len(tree.selectedNodes)==0:
            return

        last = tree.lastPickedNode
        all = self.parent.children
        try:
            i1 = all.index(self)
            i2 = all.index(last)
        except ValueError:
            raise ValueError, "range only work over siblings"

        if i1 > i2:
            tmp=i1; i1=i2; i2=tmp

        tree.selectionHistory.append(tree.selectedNodes[:])

        #print last, i1, i2+1
        for node in all[i1:i2+1]:
            if node==last:
                continue
            if node.isSelected:
                node.deselect(history=False)
            else:
                node.select(only=False, history=False)

    
    def modifySelection_cb(self, event=None):
        tree = self.tree()
        tree.lastPickedNode = self
        if tree.selectionMode=='single':
            return
        self.toggleSelection(only=False)


    def undoSelect(self, event=None):
        tree.lastPickedNode = None
        self.tree().undoSelect()


    def refresh(self):
        self.needsRedraw = True
        try:
            self.label["text"] = unicode(self)
        except:
            pass


    def refreshSubTree(self):
        self.refresh()
        for child in self.children:
            child.refreshSubTree()
        self.tree().redraw()
        

    def refreshChildren(self, redraw=True):
        ## bizarre
        if self.hasBeenExpanded:
            if self.isExpanded:
                tree = self.tree()
                # save current children
                oldchildren = self.children
                # delete all icons for current children
                for child in self.children:
                    child.deleteNodeIcons()
                # create Nodes for new children while saving existing ones
                oldChildObj = {}
                for c in oldchildren: # dict of obj:existing nodes
                    oldChildObj[c.object]= c
                i = 0
                newchildren = [] # build the list and create new children
                for childobj in self.object.children:
                    try:
                        c = oldChildObj[childobj]
                    except KeyError:
                        c = self.__class__(childobj, self)
                        tree.objectToNode[c.object] = c
                    newchildren.append(c)
                    c.childNumber = i
                    i = i + 1
                self.children = newchildren

                if len(self.children)==0:
                    self.isExpanded = 0 # Cannot be expanded if no children
                if redraw:
                    self.tree().redraw()
            else:
                for child in self.children:
                    child.destroy()
                self.children = []
                self.hasBeenExpanded = False
                if redraw:
                    self.redraw(self.currenty)
        else:
            if redraw:
                self.redraw(self.currenty)


##     def selectTree(self):
##         if not self.isSelected: self.select()
##          for child in self.children:
##              child.selectTree()


##     def deselectTree(self):
##         self.deselect()
##         for child in self.children:
##             child.deselectTree()

## #root = Tk()
## #iconsmanager = IconsManager('Icons/32x32/', 'Pmv')
## #ico = iconmanager.get('ss.png', root)

if __name__=='__main__':

    from mglutil.util.callback import CallbackFunction

    from MolKit.molecule import Atom, Molecule
    from MolKit.protein import Residue, Chain

    class ObjectTree(Tree):
        """Each node in the tree has an object associated in the node's .object
    attribute.  The objects are expected to have a .parent and a .children
    attribute describing the hierarchy."""

        def __init__(self, master, root, iconsManager=None, idleRedraw=True,
                     nodeHeight=20, **kw):
            Tree.__init__(self, master, root, iconsManager=None, idleRedraw=True,
                          nodeHeight=20, **kw)
            self.objectToNode = {}  # key is object, values if Node instance
            self.objectToNode[root.object] = root

            self.nbCol = 10
            self.menuEntries = []
            for i in range(self.nbCol):
                self.menuEntries.append([])
            self.menuEntries[0] = [
                'displayLines', 'display S&B', 'display CPK',
                'display second. struct.', 'display molecular surface',
                'display second and S&B'
                'undisplayLines', 'undisplay S&B', 'undisplay CPK',
                'undisplay second. struct.', 'undisplay molecular surface',
                'undisplay second and S&B'
                ]
            self.menuEntries[1] = [
                'label Atoms', 'label residues', 'label chains', 'label molecules'
                ]
            self.menuEntries[2] = [
                'color by atom types', 'color by molecule', 'color by chains',
                'color by residue (RASMOL)', 'color by residue (SHAPELY)',
                'color by DG', 'color by instance', 'color by second. struct.'
                ]

        def createNodeForObject(self, object):
            """given an object if the corresponding Node is not yet created
    for its creation by expanding all its ancestors"""
            try:
                return self.objectToNode[object]
            except KeyError:
                p = object.parent
                ancestors = [p]
                while p and not self.objectToNode.get(p, None):
                    ancestors.append(p)
                    p = p.parent
                ancestors.append(p)
                ancestors.reverse()
                for p in ancestors:
                    self.objectToNode[p].expand()
                return self.objectToNode[object]


    class ObjectNode(Node):
        """ """

        def __repr__(self):
            return self.object.name

        def getIcon(self):
            """return node's icons"""
            iconsManager = self.tree().iconsManager
            object = self.object
            
            if isinstance(object, Atom):
                icon = iconsManager.get("atom.png", self.tree().master)
            elif isinstance(object, Residue):
                icon = iconsManager.get("residue.png", self.tree().master)
            elif isinstance(object, Chain):
                icon = iconsManager.get("chain.png", self.tree().master)
            elif isinstance(object, Molecule):
                icon = iconsManager.get("ms.png", self.tree().master)
            else:
                icon = None
                
            if icon:
                self.iconWidth = icon.width()
            else:
                self.iconWidth = 0
            return icon



        def isExpandable(self):
            """Returns true if this node has children"""
            return len(self.object.children)


        def createChildrenNodes(self):
            """Create node for all children of self.object"""
            children = []
            for child in self.object.children:
                children.append(ObjectNode(child, self))
            return children


        def drawNodeCustomization(self, x, y):
            """Add things to the rigth side of the tree
    """
            tree = self.tree()
            canvas = tree.canvas
            level = 2
            col = ['gray75', 'red', 'cyan', 'green', 'yellow']

            nbButtons = 10
    ##         if not hasattr(self, 'chkbt'):
    ##             self.chkbtVar = []
    ##             self.chkbt = []
    ##             for i in range(nbButtons):
    ##                 v = IntVar()
    ##                 self.chkbtVar.append(v)
    ##                 cb = CallbackFunction(self.buttonClick, i )
    ##                 button = Checkbutton(
    ##                     canvas, variable=v, command=cb, padx=0, pady=0,
    ##                     background=col[level-1], anchor='nw')
    ##                 self.chkbt.append(button)

            x = 150
            fill = ['white', 'green']
            if not hasattr(self, 'chkbt'):
                nbButtons = tree.nbCol
                self.chkbt = [0]*nbButtons
                self.chkbtid = [None]*nbButtons

                # create a pull down menu and allocate variables for each column
                self.menu = Menu(canvas, title='Choose', tearoff=False)
                self.menuVars = []  # list of [column][menu entries] IntVar
                for i, entries in enumerate(tree.menuEntries):
                    l = []
                    for j in range(len(entries)):
                        l.append(IntVar())
                    self.menuVars.append(l) # a list for each column

            for i, val in enumerate(self.chkbt):
                xo = x+i*35
                cid = canvas.create_oval(xo, y, xo+15,y+15, fill=fill[val])
                cb = CallbackFunction(self.buttonClick, i )
                canvas.tag_bind(cid, "<Button-1>", cb)
                cb = CallbackFunction(self.shiftButtonClick, i )
                canvas.tag_bind(cid, "<Shift-Button-1>", cb)
                self.chkbtid[i] = cid
                self.canvasIDs.append(cid)
    ##             button = self.chkbt[i]
    ##             cid = canvas.create_window( x+i*35, y, window=button,
    ##                                         width=20, height=15)
    ##             self.canvasIDs.append(cid)

                # add secondary structure glyph
                molFrag = self.object
                if isinstance(molFrag, Residue):
                    if hasattr(molFrag, 'secondarystructure'):
                        ssname = molFrag.secondarystructure.name
                        if ssname[:6]=='Strand': color = '#FFF700'
                        elif ssname[:4]=='Coil': color = 'grey45'
                        elif ssname[:5]=='Helix': color = '#FF198C'
                        elif ssname[:4]=='Turn': color = 'blue'

                        cid = canvas.create_rectangle(
                            130, y-10, 140, y+10, outline=color,fill=color)
                        self.canvasIDs.append(cid)

                        cid = canvas.create_text( 152, y, text=ssname, anchor='nw')
                        self.canvasIDs.append(cid)
            return x+i*35
    ##             func = tree.buttonValFunc[i]
    ##             if func:
    ##                 func(self)

        def menu_cb(self, column, what, menuEntryIndex, event=None):
            print 'FFFF', column, what, menuEntryIndex


        def buttonClick(self, column, event=None):
            # get called for each checkbutton
            tree = self.tree()

            if self.chkbt[column]==0:
                tree.canvas.itemconfigure(self.chkbtid[column], fill='green')
                self.chkbt[column] = 1
            else:
                tree.canvas.itemconfigure(self.chkbtid[column], fill='white')
                self.chkbt[column] = 0

            for i,v in enumerate(self.menuVars[column]):
                if v.get():
                    print 'DDD', tree.menuEntries[column][i], self.chkbt[column]


        def shiftButtonClick(self, column, event=None):
            # get called for each checkbutton
            tree = self.tree()
            menu = self.menu

            if menu.index(0)==0: # there is something in the menu, remove it
                menu.delete(0, 100)

            for i, menuEntry in enumerate(tree.menuEntries[column]):
                v = self.menuVars[column][i]
                cb = CallbackFunction(self.menu_cb, column, menuEntry, i)
                menu.add_checkbutton(label=menuEntry, variable=v, command=cb)
            menu.post(event.x_root, event.y_root)

    ##         if self.chkbt[column]==0:
    ##             tree.canvas.itemconfigure(self.chkbtid[column], fill='green')
    ##             self.chkbt[column] = 1
    ##         else:
    ##             tree.canvas.itemconfigure(self.chkbtid[column], fill='white')
    ##             self.chkbt[column] = 0

    ##   def getIcon(self):
    ##       """return node's icons"""
    ##       return None

    from MolKit import Read
    from MolKit.molecule import MolecularSystem
    syst = MolecularSystem ('world')
    mols = Read('../dev23/1crn.pdb')
    syst.adopt(mols[0])
    
    #mols = Read('../dev23/2plv.pdb')
    #syst.adopt(mols[0])
    #mols = Read('../dev23/1gav.pdb')

    root = Tk()
    rootnode = ObjectNode(syst, None)
    tree = ObjectTree(root, rootnode, selectionMode='multiple')
    tree.pack(expand=1, fill="both")
