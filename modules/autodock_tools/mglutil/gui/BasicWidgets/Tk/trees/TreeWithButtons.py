from tree import Tree, Node
from Tkinter import Tk, Label, ALL, Menu, IntVar, StringVar, Toplevel
from Pmw import Balloon
from mglutil.util.callback import CallbackFunction

from MolKit.molecule import Atom, Molecule
from MolKit.protein import Residue, Chain


## BUGS: select a range and expand inside the range
##       select 2 ranges

class ColumnDescriptor:
    """ """
    def __init__(self, name, cmd, btype='checkbutton', buttonSize=(15,15),
                 buttonShape='circle', buttonColors = ['white', 'green'],
                 inherited=True, title=None, color='black',
                 objClassHasNoButton=None):

        self.tree = None
        self.canvas = None
        self.inherited = inherited
        self.title = title
        self.name = name
        self.color = 'black', #color # column color for outline and label
        
        self.cmd = cmd  # should be (cmd, *args, **kw)

        # menu for command selection panel
        self.menu = None

        # menu for command option setting
        self.optMenu = None

        assert btype in ['checkbutton', 'button']
        self.commandType = btype

        self.buttonSize = buttonSize
        self.buttonShape = buttonShape
        self.buttonColors = buttonColors

        # list classes of objects for which this column should have no buttons
        if objClassHasNoButton is None:
            self.objClassHasNoButton = []
        else:
            self.objClassHasNoButton = objClassHasNoButton


    def isOn(self, node):
        # subclass to implement decision for checkbutton value on/off
        # when node is first created
        return 0

    
    def optMenu_cb(self, node, column, event=None):
        cmd, args, kw = self.cmd
        values = cmd.showForm(posx=event.x_root, posy=event.y_root, master=self.tree.master)
        if values and len(values)==0: return # Cancel was pressed
        cmd.lastUsedValues['default'].update(values)
        node.buttonClick(column, val=1) # always call with button on
        

    def execute(self, node, colInd):
        objects = node.getObjects(colInd)
        val = node.chkbtval[colInd]
        cmd, args, kw = self.cmd
        defaultValues = cmd.getLastUsedValues()
        defaultValues.update( kw )
        
        if  self.commandType == 'checkbutton':
            if defaultValues.has_key('negate'):
                defaultValues['negate'] = not val
            elif not val:
                return
        for objs in objects:
            cmd ( *((objs,)+args), **defaultValues)



class TreeWithButtons(Tree):
    """Each node in the tree has an object associated in the node's .object
attribute.  The objects are expected to have a .parent and a .children
attribute describing the hierarchy."""
    
    def __init__(self, master, root, iconsManager=None, idleRedraw=True,
                 nodeHeight=20, **kw):
        Tree.__init__(self, master, root, iconsManager=None, idleRedraw=True,
                      nodeHeight=20, **kw)
        self.objectToNode = {}  # key is object, values if Node instance
        self.objectToNode[root.object] = root
        
        self.nbCol = 0  # number of columns of buttons
        self.columns = [] # list of ColumnDescriptor instance
        self.balloon = Balloon(master)

        self.colLabIds = []

        self.colWidth = 25 # width of each column in pixels
        self.treeWidth = 180 # width of the tree part of the widget
        self.prevX = 0
        self.prevY = 0

        self.crosshairTk = None
        self.canvas.bind("<Leave>", self.leave_cb)


    ## CROSSHAIR
    def enter_cb(self, event=None):
        """add crosshair"""
        #print 'enter', event.widget
        Tree.enter_cb(self)
        canvas = self.canvas
        self.crosshairTk = canvas.create_line(0,0,0,0,0,0,0,0,0,0,
                                              fill='#35A8FF', width=2)
        self.master.bind("<Motion>", self.move_cb)
        canvas.tag_lower(self.crosshairTk)


    def leave_cb(self, event=None):
        """destroy crosshair"""
        #print 'leave', event.widget
        self.master.unbind("<Motion>")
        self.canvas.delete(self.crosshairTk)


    def move_cb(self,event):
        """move crosshair"""
        canvas = self.canvas
        x1, y1, x2, y2 = canvas.bbox(ALL)
        x = canvas.canvasx(event.x)
        y = canvas.canvasy(event.y)
        if isinstance(event.widget, Label):
            x = x + event.widget.winfo_x()
            y = y + event.widget.winfo_y()
        h = self.headerHeight - 10
        x2 += 100
        y2 += 100
        canvas.coords( self.crosshairTk, x1, y, x2, y, x2, y2, x, y2, x, h )


    def redrawHeader(self, *args):
        ypos = (self.fisrtVisibleNodeLineNum * self.nodeHeight) + \
               (self.headerHeight * 0.5)
        x, y = self.canvas.coords('ColHeaders')
        self.canvas.move('ColHeaders', 0, ypos-y)
        self.canvas.move('ColHeadersWidgets', 0, ypos-y)


    def undisplay(self):
        for node in self.displayedNodes:
            node.deleteNodeIcons()
        
        
    def resetDisplayedNodes(self):
        # delete these list to force the nodes to recreate them
        for node in self.displayedNodes:
            node.chkbtid = []
            node.needsRedraw = True


    def addColumnDescriptor(self, columnDescr):
        assert isinstance(columnDescr, ColumnDescriptor)
        columnDescr.tree = self
        columnDescr.canvas = self.canvas

        # make titles use alternate lines
        if self.nbCol%2:
            columnDescr.title = '\n'+columnDescr.title
        else:
            columnDescr.title += '\n'

        self.createColHeader(columnDescr, self.nbCol)

        #columnDescr.setMenuEntries(columnDescr.menuEntries)
        self.columns.append(columnDescr)
        self.root.chkbtval.append(0)
        self.nbCol += 1
        self.resetDisplayedNodes()
        self.redraw()


    def createColHeader(self, columnDescr, number):
        cw = self.colWidth
##         x = self.treeWidth + len(self.columns)*cw - 5
##         cid = self.canvas.create_line(x, 50000, x, self.headerHeight,
##                                       x+cw, self.headerHeight, x+cw, 50000)

        # label is centered on self.treeWidth + self.colWidth/2 *i*self.colWidth
        x = self.treeWidth + (0.5*self.colWidth)
        id = self.canvas.create_text(x+number*cw-4, 18,
                                     text=columnDescr.title, justify='center',
                                     tags=('ColHeaders',),
                                     fill=columnDescr.color,
                                     font="Arial 9")
        cb = CallbackFunction(self.rightButtonClick, columnDescr)
        self.canvas.tag_bind(id, "<Button-3>", cb)
        self.colLabIds.append(id)


    def rightButtonClick(self, columnDescr, event):
        print 'right click on', columnDescr

        
    def deleteColumnDescriptor(self, columnDescr):
        if isinstance(columnDescr, int):
            columnDescr = self.columns[columnDescr]
            ind = columnDescr
        else:
            ind = self.columns.index(columnDescr)

        assert isinstance(columnDescr, ColumnDescriptor)

        columnDescr.tree = None
        columnDescr.canvas = None

        self.columns.remove(columnDescr)
        chkbtval = self.root.chkbtval
        self.root.chkbtval = chkbtval[:ind] + chkbtval[ind+1:]
        self.nbCol -= 1
        self.resetDisplayedNodes()
        self.redraw()


    def insertColumnDescriptor(self, column, columnDescr):

        assert isinstance(columnDescr, ColumnDescriptor)
        columnDescr.tree = self
        columnDescr.canvas = tree.canvas

        self.columns.insert(column, columnDescr)
        self.root.chkbtval.insert(column, 0)
        self.nbCol += 1
        self.resetDisplayedNodes()
        self.redraw()


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


    def manageChildren(self, node, column):
        # sets the checkbutton of all children to value of parent (i.e. node)
        if len(node.children)==0:
            return

        val = node.chkbtval[column]
        for c in node.children:
            c.chkbtval[column] = val
            if c in self.displayedNodes:
                c.needsRedraw = True
            self.manageChildren(c, column)
        self.redraw()



class NodeWithButtons(Node):
    """ """
    def __init__(self, object, parent):

        Node.__init__(self, object, parent)

        # list of values used to emulate buttons and associated Tk ids
        # created in drawNodeCustomization
        if parent:
            chkbtval = []
            pchkbtval = parent.chkbtval
            for i, c in enumerate(parent.tree().columns):
                if c.inherited:
                    chkbtval.append(pchkbtval[i])
                else:
                    chkbtval.append(c.isOn(self))
            self.chkbtval = chkbtval
        else:
            self.chkbtval = []
        self.chkbtid = []


    def __repr__(self):
        if hasattr(self.object,'name'):
            return self.object.name
        else:
            return ""
            


    def isExpandable(self):
        """Returns true if this node has children"""
        if hasattr(self.object,'children'):
            return len(self.object.children)
        else:
            return False


    def drawNodeCustomization(self, x, y):
        """Add things to the rigth side of the tree
"""
        tree = self.tree()
        canvas = tree.canvas
        level = 2
        col = ['gray75', 'red', 'cyan', 'green', 'yellow']

        nbButtons = tree.nbCol
        x = tree.treeWidth
        # create Tk variables if needed
        if len(self.chkbtid)==0:
            nbButtons = tree.nbCol
            self.chkbtid = [None]*nbButtons

        i = 0
        cw = tree.colWidth
        for i, val in enumerate(self.chkbtval):
            col = tree.columns[i]
            if self.object.__class__ in col.objClassHasNoButton:
                continue
            xo = x+i*cw
            shape = col.buttonShape
            sx, sy = col.buttonSize
            fill = col.buttonColors
            color = col.color
            width = 1
            if shape=='circle':
                cid = canvas.create_oval(xo, y, xo+sx,y+sx, outline=color,
                                         width=width, fill=fill[val])
            elif shape=='oval':
                cid = canvas.create_oval(xo, y, xo+sx,y+sy, outline=color,
                                         width=width, fill=fill[val])
            elif shape=='square':
                cid = canvas.create_rectangle(xo, y, xo+sx,y+sx, outline=color,
                                              width=width, fill=fill[val])
            elif shape=='rectangle':
                cid = canvas.create_rectangle(xo, y, xo+sx,y+sy, outline=color,
                                              width=width, fill=fill[val])
            elif shape=='diamond':
                hw = sx*.5
                hh = sy*.5
                cid = canvas.create_polygon( xo,    y+hh, xo+hw, y,
                                             xo+sx, y+hh, xo+hw, y+sy,
                                             xo,    y+hh, fill=fill[val],
                                             width=width, outline=color)

            cb = CallbackFunction(self.buttonClick, i )
            canvas.tag_bind(cid, "<ButtonRelease-1>", cb)
            #cb = CallbackFunction(self.postCmdMenu, i )
            #canvas.tag_bind(cid, "<Button-3>", cb)
            #cb = CallbackFunction(self.postCmdOptMenu, i )
            #cb = CallbackFunction(col.optMenu_cb)
            cb = CallbackFunction(col.optMenu_cb, self, i )
            canvas.tag_bind(cid, "<ButtonRelease-3>", cb)
            balloon = """%s
left-click execute
right click set options"""%col.name
            tree.balloon.tagbind(canvas, cid, balloon)
            self.canvasIDs.append(cid)
            self.chkbtid[i] = cid

            nodeHeight = tree.nodeHeight

            # add secondary structure glyph
            molFrag = self.object
            if isinstance(molFrag, Residue):
                if hasattr(molFrag, 'secondarystructure'):
                    ssname = molFrag.secondarystructure.name
                    if ssname[:6]=='Strand': color = '#FFFC6D'
                    elif ssname[:4]=='Coil': color = 'grey45'
                    elif ssname[:5]=='Helix': color = '#FF198C'
                    elif ssname[:4]=='Turn': color = 'blue'

                    tw = tree.treeWidth
                    cid = canvas.create_rectangle(
                        tw-15, y, tw-5, y+nodeHeight, outline=color, fill=color)

                    self.canvasIDs.append(cid)

                    tree.balloon.tagbind(canvas, cid, ssname)
                    
#                    cid = canvas.create_text( 152, y, text=ssname, anchor='nw')
#                    self.canvasIDs.append(cid)
        return x + (i*cw)
##             func = tree.buttonValFunc[i]
##             if func:
##                 func(self)

    def deleteNodeIcons(self):
        Node.deleteNodeIcons(self)
        self.chkbtid = []


    def set(self, column, value):
        tree = self.tree()
        col = tree.columns[column]
        fill = col.buttonColors
        val = self.chkbtval[column]
        
        if value:
            if val: return
            val = self.chkbtval[column] = 1
        else:
            if not val: return
            val = self.chkbtval[column] = 0

        if len(self.chkbtid): # length==0 when node is not drawn
            if self.chkbtid[column]:
                tree.canvas.itemconfigure(self.chkbtid[column], fill=fill[val])

            if col.commandType=='button' and val==1:
                cb = CallbackFunction(self.resetButtons, column)
                tree.master.after(100, cb)

            if col.inherited:
                tree.manageChildren(self, column)
        

    def toggle(self, column):
        if not self.chkbtval[column]:
            self.set(column, 1)
        else:
            self.set(column, 0)
##         if val==0:
##             val = self.chkbtval[column] = 1
##         else:
##             val = self.chkbtval[column] = 0

##         if len(self.chkbtid): # length==0 when node is not drawn
##             tree.canvas.itemconfigure(self.chkbtid[column], fill=fill[val])

##             if col.commandType=='button' and val==1:
##                 cb = CallbackFunction(self.resetButtons, column)
##                 tree.master.after(100, cb)

##             if col.inherited:
##                 tree.manageChildren(self, column)


    def buttonClick(self, column, event=None, val=None, execute=True):
        # get called for each checkbutton on left or right click
        # on right click the value is forced to 1

        #print 'Click', event.widget
        if val is None:
            self.toggle(column)
        else:
            self.set(column, val)
        tree = self.tree()
        col = tree.columns[column]
        
        if execute:
            col.execute(self, column)


    def resetButtons(self, column):
        if self.isSelected:
            for n in self.tree().selectedNodes:
                n.chkbtval[column] = 1
                n.buttonClick(column, execute=False)
        else:
            self.chkbtval[column] = 1
            self.buttonClick(column, execute=False)
        
        
    def postCmdMenu(self, column, event=None):
        tree = self.tree()
        menu = tree.columns[column].menu
        menu.post(event.x_root, event.y_root)


    def postCmdOptMenu(self, column, event=None):
        tree = self.tree()
        menu = tree.columns[column].optMenu
        menu.post(event.x_root, event.y_root)
