from tree import Tree, Node
from Tkinter import Tk, Label, ALL, Menu, IntVar, StringVar, Toplevel, TclError
from Pmw import Balloon
from mglutil.util.callback import CallbackFunction

from MolKit.molecule import Atom, Molecule
from MolKit.protein import Residue, Chain


## BUGS: select a range and expand inside the range
##       select 2 ranges

class ColumnDescriptor:
    """ """
    def __init__(self, name, cmd, btype='checkbutton', 
                 buttonShape='circle', buttonColors = ['white', 'green'],
                 inherited=True, title=None, color='black',
                 objClassHasNoButton=None, showPercent=False,
                 buttonBalloon=None, onButtonBalloon=None,
                 offButtonBalloon=None):

        self.tree = None
        self.canvas = None
        self.inherited = inherited
        self.title = title
        self.name = name
        self.color = 'black', #color # column color for outline and label
        self.showPercent = showPercent
        
        self.cmd = cmd  # should be (cmd, *args, **kw)

        self.cbOn = None # used to specify a call back for the on button
        self.cbOff = None # used to specify a call back for the on button
        self.onOnly = False # when True we create only 1 button for OnOff

        if buttonBalloon is None:
            buttonBalloon = """left-click execute cmd on %s
right click set options"""
        self.buttonBalloon = buttonBalloon

        if onButtonBalloon is None:
            onButtonBalloon = """left-click execute cmd on %s
right click set options"""
        self.onButtonBalloon = onButtonBalloon

        if offButtonBalloon is None:
            offButtonBalloon = """left-click execute cmd on %s
right click set options"""
        self.offButtonBalloon = offButtonBalloon

        
        # menu for command selection panel
        self.menu = None

        # menu for command option setting
        self.optMenu = None

        assert btype in ['checkbutton', 'button'], btype
        self.commandType = btype

        self.buttonShape = buttonShape
        self.buttonColors = buttonColors

        # list classes of objects for which this column should have no buttons
        if objClassHasNoButton is None:
            self.objClassHasNoButton = [NodeWithoutButtons]
        else:
            self.objClassHasNoButton = objClassHasNoButton


    def isOn(self, node):
        # subclass to implement decision for checkbutton value on/off
        # when node is first created
        return 0

    
    def optMenu_cb(self, node, column, event=None):
        cmd, args, kw = self.cmd
        values = cmd.showForm(posx=event.x_root, posy=event.y_root,
                              master=self.tree.master)
        if values=={}: return # Cancel was pressed
        val = 1
        if values.has_key('display'):
            if values['display']=='undisplay':
                val = 0

        cmd.lastUsedValues['default'].update(values)
        node.buttonClick(column, val=val) # always call with button on
        

    def _getNodes(self, node, colInd):
        # added so that MVColumnDescriptor can override
        #print 'TREEWITHBUTTONS _getNodes', node
        return node.getObjects(colInd)


    def execute(self, node, colInd):
        #print 'ColumnDescriptor.execute', node, colInd
        objects = self._getNodes(node, colInd)
        val = node.chkbtval[colInd]
        cmd, args, kw = self.cmd
        defaultValues = cmd.getLastUsedValues()
        defaultValues.update( kw )
        
        if  self.commandType == 'checkbutton':
            if defaultValues.has_key('negate'):
                defaultValues['negate'] = not val
            elif not val:
                return
        for objs in objects: # loop over list that might contain an atom set
                             # residue set, chain set and/or molecule set
            #print 'ColumnDescriptor execute GGGG', cmd, objs, args, defaultValues
            cmd ( *((objs,)+args), **defaultValues)



class TreeWithButtons(Tree):
    """Each node in the tree has an object associated in the node's .object
attribute.  The objects are expected to have a .parent and a .children
attribute describing the hierarchy."""
    
    def __init__(self, master, root, iconsManager=None, idleRedraw=True,
                 nodeHeight=15, headerHeight=30, treeWidth=140, **kw):
        Tree.__init__(self, master, root, iconsManager=iconsManager,
                      idleRedraw=idleRedraw, nodeHeight=nodeHeight,
                      headerHeight=headerHeight, **kw)

        self.objectToNode = {}  # key is object, values if Node instance
        self.objectToNode[root.object] = root
        
        self.nbCol = 0  # number of columns of buttons
        self.columns = [] # list of ColumnDescriptor instance
        self.balloon = Balloon(master)

        self.colLabIds = []

        self.colWidth = 17 # width of each column in pixels
        w = self.treeWidth = treeWidth # width of the tree part of the widget
        self.newTreeWidth = 0 # used to widen tree whebn labels gets long
        self.prevX = 0
        self.prevY = 0

        # draw the divider between tree and buttons
        canvas = self.canvas

        id_ = canvas.create_line( w-6, 0, w-6, 1000, fill='grey75', width=3)
        self.dividerCanvasId = id_
        canvas.tag_bind(id_,"<Enter>", self.enterDivider_cb)
        canvas.tag_bind(id_,"<Leave>", self.leaveDivider_cb)
        canvas.tag_bind(id_,"<Button-1>", self.dividePress_cb)
        canvas.tag_bind(id_,"<ButtonRelease-1>", self.divideRelease_cb)
        #id_ = canvas.create_line( w+1, 0, w+1, 1000, fill='grey75')
        #self.dividerCanvasIds.append(id_)
        
        self.crosshairTk = None
        self.canvas.bind("<Leave>", self.leave_cb)
        self.circlesForOnOff = True
        self.lastHighligted = None

    def enterDivider_cb(self, event):
        self.canvas.config(cursor='sb_h_double_arrow')


    def leaveDivider_cb(self, event):
        self.canvas.config(cursor='')


    def dividePress_cb(self, event):
        self.treeWidthOff = off = event.x_root-self.treeWidth
        self.canvas.bind("<Motion>", self.moveDivider_cb)


    def divideRelease_cb(self, event):
        self.canvas.unbind("<Motion>")
        #self.redraw()

        
    def moveDivider_cb(self, event):
        tw = event.x_root-self.treeWidthOff
        off = tw - self.treeWidth
        self.setTreeWidth(event.x_root - self.treeWidthOff)


    def setTreeWidth(self, width):
        off = width-self.treeWidth
        self.treeWidth = width
        self.canvas.move(self.dividerCanvasId, off, 0)
        self.canvas.move('ColHeaders', off, 0)
        self.canvas.move('backgroundStripes', off, 0)
        self.canvas.move('button', off, 0)

        
    def reallyRedraw(self):
        Tree.reallyRedraw(self)
        if self.newTreeWidth > self.treeWidth:
            self.setTreeWidth(self.newTreeWidth)
            self.newTreeWidth = 0
            self.redraw()

        
    def redraw(self, event=None):
        Tree.redraw(self)
        if self.crosshairTk and event:
            self.move_cb(event)

    ## CROSSHAIR
    def enter_cb(self, event=None):
        """add crosshair"""
        #print 'enter', event.widget
        Tree.enter_cb(self)
        canvas = self.canvas
        self.crosshairTk = canvas.create_line(0,0,0,0,0,0,0,0,0,0,
                                              fill='#35A8FF', width=2)
        canvas.bind("<Motion>", self.move_cb)
        canvas.tag_lower(self.crosshairTk)
        canvas.tag_lower('backgroundStripes')


    def leave_cb(self, event=None):
        """destroy crosshair"""
        #print 'leave', event.widget
        self.canvas.unbind("<Motion>")
        self.canvas.delete(self.crosshairTk)
        self.crosshairTk = None


    def move_cb(self,event):
        """move crosshair"""
        #print 'move crosshair', event.widget
        canvas = self.canvas
        x1, y1, x2, y2 = canvas.bbox(ALL)
        #print 'move crosshair',  x1, y1, x2, y2
        x = canvas.canvasx(event.x)
        y = canvas.canvasy(event.y)
        if isinstance(event.widget, Label):
            x = x + event.widget.winfo_x()
            y = y + event.widget.winfo_y()
        h = 0#self.headerHeight - 10
        x2 += 100000
        y2 += 100000
        canvas.coords( self.crosshairTk, x1, y, x2, y, x2, y2, x, y2, x, h )

        # take care of higlighting on or off button circle for OnOff buttons
        if self.lastHighligted:
            canvas.itemconfig(self.lastHighligted,
                              outline=self.lastHighligtedColor, width=1)
        tags = canvas.gettags('current')
        if len(tags):
            if tags[0]=='on' or tags[0]=='off':
                cid = int(tags[2])
                self.lastHighligted = cid
                self.lastHighligtedColor = canvas.itemcget(cid, 'outline')
                canvas.itemconfig(cid, outline='blue', width=2)
            else:
                self.lastHighligted = None
            
            
    def redrawHeader(self, *args):
        ypos = (self.firstVisibleNodeLineNum * self.nodeHeight) + \
               (self.headerHeight * 0.5)
        dims = self.canvas.coords('ColHeaders')
        if len(dims):
            x, y = dims
        else:
            y = 0
        self.canvas.move('ColHeaders', 0, ypos-y-5.5)
        #self.canvas.move('backgroundStripes', 0, ypos-y-5.5)
        dims = self.canvas.coords('dashButtons')
        if len(dims):
            x, y = dims
        else:
            y = 0
        #print 'REDRAW1', x, y, ypos, 0, ypos-y
        self.canvas.move('dashButtons', 0, ypos-y)


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
        #if self.nbCol%2:
        #    columnDescr.title = '\n'+columnDescr.title
        #else:
        #    columnDescr.title += '\n'

        self.createColHeader(columnDescr, self.nbCol)

        #columnDescr.setMenuEntries(columnDescr.menuEntries)
        self.columns.append(columnDescr)
        self.root.chkbtval.append(0)
        self.nbCol += 1
        self.resetDisplayedNodes()
        self.redraw()


    def setColumnWidth(self, cw):
        self.colWidth = cw
        self.canvas.delete('backgroundStripes')
        self.canvas.delete( 'ColHeaders')
        for i, col in enumerate(self.columns):
            self.createColHeader(col, i)
            
        self.redraw()


    def createColHeader(self, columnDescr, number):
        cw = self.colWidth
        x = self.treeWidth + (0.5*self.colWidth)

        # add column background
        if number%2:
            _id = self.canvas.create_rectangle(
                x+number*cw-cw/2, 0, x+number*cw+cw/2, 10000,
                tags=('backgroundStripes',), fill='#EEEEEE', outline='#EEEEEE')
            
        # add title and icon
        _id = self.canvas.create_text(
            x+number*cw, 7, text=columnDescr.title, justify='center',
            tags=('ColHeaders',), fill=columnDescr.color, font=self.font)
        cb = CallbackFunction(self.rightButtonClick, columnDescr)
        self.canvas.tag_bind(_id, "<Button-3>", cb)
        self.colLabIds.append(_id)

        # add column header icons
        if columnDescr.iconfile:
            if columnDescr.icon is None:
                columnDescr.getIcon(columnDescr.iconfile)
            _id = self.canvas.create_image(
                x+number*cw, 19,tags=('ColHeaders',),
                    image=columnDescr.icon)
        self.colLabIds.append(_id)


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
    def __init__(self, object, parent, buttonType=None):

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
        self.buttonType = buttonType # can be 'OnOffButtons'


    def __repr__(self):
        if hasattr(self,'name'):
            return self.name
        if hasattr(self.object,'name'):
            return self.object.name
        else:
            return ""
            

    def isExpandable(self):
        """Returns true if this node has children"""
        if hasattr(self.object,'children'):
            return len(self.getChildren())
        else:
            return False


    def drawNodeLabel(self, x, y):
        # sub class to increase tree width of label gets wider
        x = Node.drawNodeLabel(self, x, y)       
        tree = self.tree()
        canvas = tree.canvas
        bb = canvas.bbox(self.labelTkid)
        if bb[2]>tree.treeWidth:
            tree.newTreeWidth = bb[2]
        return x
    

    def drawNodeCustomization(self, x, y):
        """
        Add things to the rigth side of the tree
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
            cw2 = cw/2
            midx = xo+cw2
            midy = y + tree.fontSize/2 + 2
            cid1 = None
            color = col.color
            width = 1

            if self.buttonType=='OnOffButtons' and not col.onOnly:
                # overrides col.commandType and put 2 buttons for on off
                #uniq = '%d_%d'%(midx,midy)
                cid = canvas.create_oval(
                    midx-cw2+1, midy-cw2+1, midx+1, midy+1, outline='#3D8E54',
                    width=width, fill='#44A05E', tags=('on','button'))
                canvas.addtag_withtag(str(cid), cid)

                if col.cbOn is None:
                    cb1 = CallbackFunction(self.onOffButtonClick, i, 1)
                else:
                    cb1 = CallbackFunction(col.cbOn, self, i, 1)
                    
                cb2 = CallbackFunction(col.optMenu_cb, self, i)
                balloon = col.onButtonBalloon%(self.object.name)

                cid1 = canvas.create_oval(
                    midx-1, midy-1, midx+cw2-1, midy+cw2-1, outline='#B52121',
                    width=width, fill='#E2282D', tags=('off','button'))
                canvas.addtag_withtag(str(cid1), cid1)

                if col.cbOff is None:
                    cb1_1 = CallbackFunction(self.onOffButtonClick, i, 0 )
                else:
                    cb1_1 = CallbackFunction(col.cbOff, self, i, 0)
                cb2_1 = None
                balloon1 = col.offButtonBalloon%(self.object.name)

            else:
                shape = col.buttonShape
                sx2 = cw/2-2
                sy2 = tree.fontSize/2
                fill = col.buttonColors
                from MolKit.listSet import ListSet
                if shape=='circle':
                    sy2 = min(sx2, sy2)
                    if col.showPercent: fillCol = fill[0]
                    else: fillCol = fill[val]
                    cid = canvas.create_oval(
                        midx-sx2, midy-sy2, midx+sx2, midy+sy2, outline=color,
                        width=width, fill=fillCol, tags=('button',))
                    if col.showPercent and not isinstance(self.object, ListSet):
                        if hasattr(self.object,col.showPercent):
                            extent = getattr(self.object,col.showPercent)*360
                            if extent>0.0:
                                #print 'YES', self.object.name, col.showPercent, extent
                                if extent==360.:
                                    # special case full disc as sometimes
                                    # extent 360 shows nothing
                                    cid1 = canvas.create_oval(
                                        midx-sx2, midy-sy2, midx+sx2, midy+sy2,
                                        fill=fill[1], outline='',
                                        tags=('button',))
                                    self.chkbtval[i] = 1
                                else:
                                    cid1 = canvas.create_arc(
                                        midx-sx2, midy-sy2, midx+sx2, midy+sy2,
                                        extent=-extent, fill=fill[1],
                                        outline='', tags=('button',))
                                    self.chkbtval[i] = 0
                            else:
                                self.chkbtval[i] = 0

                elif shape=='oval':
                    cid = canvas.create_oval(
                        midx-dx, midy-dx, midx+dx, midy+dy, outline=color,
                        width=width, fill=fill[val], tags=('button',))

                elif shape=='square':
                    dx = min(sx2, sy2)
                    cid = canvas.create_rectangle(
                        midx-dx, midy-dx, midx+dx, midy+dy, outline=color,
                        width=width, tags=('button',))# fill=fill[val])

                elif shape=='downTriangle':
                    cid = canvas.create_polygon(
                        midx-sx2, midy-sy2, midx+sx2, midy-sy2,
                        midx, midy+sy2, midx-sx2, midy-sy2,
                        outline=color, width=width, fill=fill[val],
                        tags=('button',))

                elif shape=='rectangle':
                    if col.showPercent: fillCol = fill[0]
                    else: fillCol = fill[val]
                    cid = canvas.create_rectangle(
                        midx-sx2, midy-sy2, midx+sx2, midy+sy2, outline=color,
                        width=width, fill=fillCol, tags=('button',))
                    if col.showPercent and not isinstance(self.object, ListSet):
                        if hasattr(self.object,col.showPercent):
                            extent = getattr(self.object, col.showPercent)
                            #print 'Redraw 1', self.object, extent
                            if extent>0.0:
                                cid1 = canvas.create_rectangle(
                                    midx-sx2+1, midy+sy2,
                                    midx+sx2-1, midy-sy2+2*sy2*(1.-extent),
                                    fill=fill[1], outline='', tags=('button',))
                                self.chkbtval[i] = int(extent)
                            else:
                                self.chkbtval[i] = 0
                                
                elif shape=='diamond':
                    hw = sx2
                    hh = sy2
                    cid = canvas.create_polygon(
                        midx, midy-sy, midx+sx, midy,
                        midx, midy+sy, midx-sx, midy, midx, midy-sy,
                        fill=fill[val], width=width, outline=color,
                        tags=('button',))

                if cid1:
                    cb1_1 = CallbackFunction(self.buttonClick, i )
                    cb2_1 = CallbackFunction(col.optMenu_cb, self, i )
                    balloon1 = col.buttonBalloon%(self.object.full_name())
 
                if col.cbOn is None:
                    cb1 = CallbackFunction(self.buttonClick, i )
                else:
                    cb1 = CallbackFunction(col.cbOn, self, i, 1)
                    
                if col.cbOff is None:
                    cb2 = CallbackFunction(col.optMenu_cb, self, i )
                else:
                    cb2 = CallbackFunction(col.cbOff, self, i, 0)

                balloon = col.buttonBalloon%(self.object.full_name())

            if cid1:
                canvas.tag_raise(cid1)
                canvas.tag_bind(cid1, "<ButtonRelease-1>", cb1_1)
                if cb2_1: canvas.tag_bind(cid1, "<ButtonRelease-3>", cb2_1)
                tree.balloon.tagbind(canvas, cid1, balloon1)
                self.canvasIDs.append(cid1)
                self.chkbtid[i] = cid1
                
            canvas.tag_bind(cid, "<ButtonRelease-1>", cb1)
            canvas.tag_bind(cid, "<ButtonRelease-3>", cb2)

            tree.balloon.tagbind(canvas, cid, balloon)
            self.canvasIDs.append(cid)
            self.chkbtid[i] = cid

            nodeHeight = tree.nodeHeight

            ## MS 09/10 this is molecular stuff, should not be in mglutil
            # add secondary structure glyph
            molFrag = self.object
            if isinstance(molFrag, Residue):
                if hasattr(molFrag, 'secondarystructure'):
                    ssname = molFrag.secondarystructure.name
                    if ssname[:6]=='Strand': color = '#FFFC6D'
                    elif ssname[:4]=='Coil': color = 'grey45'
                    elif ssname[:5]=='Helix': color = '#FF198C'
                    elif ssname[:4]=='Turn': color = 'blue'

                    #tw = tree.treeWidth
                    tw = self.nodeStartX
                    cid = canvas.create_rectangle(
                        tw-20, y, tw-10, y+nodeHeight, outline=color,
                        fill=color)

                    self.canvasIDs.append(cid)

                    tree.balloon.tagbind(canvas, cid, ssname)
                    cb = CallbackFunction(self.toggleResSelection,
                                          molFrag.secondarystructure.residues)
                    canvas.tag_bind(cid, "<ButtonRelease-1>", cb)
                    #cid = canvas.create_text( 152, y, text=ssname, anchor='nw')
#                    self.canvasIDs.append(cid)
        return x + (i*cw)
##             func = tree.buttonValFunc[i]
##             if func:
##                 func(self)

    def toggleResSelection(self, residues, event=None):
        tree = self.tree()
        for res in residues:
            node = tree.objectToNode[res]
            if node.isSelected:
                node.deselect()
            else:
                node.select(only=False)


    def deleteNodeIcons(self):
        Node.deleteNodeIcons(self)
        self.chkbtid = []


    def set(self, column, value):
        tree = self.tree()
        col = tree.columns[column]
        fill = col.buttonColors
        val = self.chkbtval[column]

        #print 'AAAAAAAA SET', col.name, value, val, len(self.chkbtid)
        if value:
            if val:
                #tree.redraw()
                return
            val = self.chkbtval[column] = 1
        else:
            if not val:
                #tree.redraw()
                return
            val = self.chkbtval[column] = 0

        if len(self.chkbtid): # length==0 when node is not drawn
            if self.chkbtid[column] and not col.showPercent:
                tree.canvas.itemconfigure(self.chkbtid[column], fill=fill[val])

            if col.commandType=='button' and val==1:
                cb = CallbackFunction(self.resetButtons, column)
                tree.master.after(100, cb)

            if col.inherited:
                tree.manageChildren(self, column)
        # MS added the redraw to force buttons with showPercent to update
        # a cheaper version might be possible
        #tree.redraw()
        

    def toggle(self, column):
        if not self.chkbtval[column]:
            self.set(column, 1)
        else:
            self.set(column, 0)
        self.tree().redraw()
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


    def onOffButtonClick(self, column, value, event=None, val=None, execute=True):
        # get called for onOffButtons on left click
        #print 'onOffButtonClick', value, execute, value
        self.set(column, value)
        tree = self.tree()
        col = tree.columns[column]
        
        if execute:
            col.pickEvent = event
            #print 'executing', col, column
            col.execute(self, column)
        if event:
            try: # sometimes this generates a TclError. Not sure why MS 09/2010
                tree.move_cb(event) # to redraw crosshair
            except TclError:
                pass


    def buttonClick(self, column, event=None, val=None, execute=True):
        # get called for each checkbutton on left or right click
        # on right click the value is forced to 1

        tree = self.tree()
        col = tree.columns[column]

        #print 'button click', column, val, self.chkbtval[column]
##         if col.showPercent:
##             if self.chkbtval[column]==1:
##                 self.set(column, 0)
##             else:
##                 self.set(column, 1)
##         else:
##             if val is None:
##                 self.toggle(column)
##             else:
##                 self.set(column, val)
        if val is None:
            self.toggle(column)
        else:
            self.set(column, val)
         
        if execute:
            col.pickEvent = event
            col.execute(self, column)
        tree.redraw(event)


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



class NodeWithoutButtons(NodeWithButtons):
 

    def drawNodeCustomization(self, x, y):
        """
        Add things to the rigth side of the tree
        """
        return 0
