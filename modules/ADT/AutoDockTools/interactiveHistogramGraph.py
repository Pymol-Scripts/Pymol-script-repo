#############################################################################
# 
# Graph class: 
# Written by  John E. Grayson 
#
# InteractiveHistogramBuilder class written by Ruth Huey
#
#
#############################################################################


import sys
import Tkinter
from Canvas import Line, CanvasText, Rectangle
import tkFileDialog
import string, math, os

from mglutil.util.misc import ensureFontCase
from mglutil.gui import widgetsOnBackWindowsCanGrabFocus

def minBound(nodeList):
    x = 10000000
    y = 10000000
    for x1, y1 in  nodeList:
        if x1 < x: x= x1
        if y1 < y: y= y1
    return x, y


def maxBound(nodeList):
    x = -10000000
    y = -10000000
    for x1, y1 in  nodeList:
        if x1 > x: x= x1
        if y1 > y: y= y1
    return x, y




class InteractiveGraphBuilder(Tkinter.Frame):


    def __init__(self, name='NoName', master=None, nodeList=[], **kw):
        """
        """
        self.name = name
        self.nodes = {}
        self.geoms = {}

        if kw.has_key('visibleWidth'):
            self.visibleWidth = visibleWidth
            del kw['visibleWidth']
        else:
            self.visibleWidth = 400

        if kw.has_key('visibleHeight'):
            self.visibleWidth = visibleHeight
            del kw['visibleHeight']
        else:
            self.visibleHeight = 400

        if kw.has_key('totalWidth'):
            self.visibleWidth = totalWidth
            del kw['totalWidth']
        else:
            self.totalWidth = 400

        if kw.has_key('totalHeight'):
            self.totalHeight = totalHeight
            del kw['totalHeight']
        else:
            self.totalHeight = 400

        if kw.has_key('hasScroll'):
            self.hasScroll = hasScroll
            del kw['hasScroll']
        else:
            self.hasScroll = 0

        if kw.has_key('font'):
            self.font = font
            del kw['font']
        else:
            self.font = (ensureFontCase('helvetica'), 14)

        Tkinter.Frame.__init__(self, master)
        Tkinter.Pack.config(self, expand=1, fill=Tkinter.BOTH)
        apply(self.createCanvas, (), kw)
        self.paramPanelObject = None
        self.plotarea_size = [None, None]
        Tkinter.Widget.bind(self.draw, '<Enter>', self.enter)
        # bind mouse button callbacks
        self.draw.bind("<Any-ButtonPress-1>", self.mouse1DownCanvas)
        self.draw.bind("<Any-ButtonPress-3>", self.mouse1DownCanvas)
        self.draw.bind("<Any-ButtonPress-2>", self.mouse1DownCanvas)
        self.draw.bind("<Any-Double-ButtonPress-1>", self.mouse2DownCanvas)
        self.draw.bind("<Any-Double-ButtonPress-2>", self.mouse2DownCanvas)
        self.draw.bind("<Any-Double-ButtonPress-3>", self.mouse2DownCanvas)
        self.draw.old_current = None

    
    def mouse1DownCanvas(self, event=None):
        #print 'in mouse1DC'
        obj = self.draw.find_withtag('current')
        if self.draw.old_current:
            self.draw.itemconfig(self.draw.old_current, fill = 'blue')
        self.draw.itemconfig('current', fill = 'red')
        self.draw.old_current = obj
        #print obj
        if len(obj) and obj[0] in self.geoms.keys():
            self.currentNode = self.geoms[obj[0]]
            #print 'currentNode=', self.currentNode
            if hasattr(self, 'reverseIndex'):
                #print self.reverseIndex[self.currentNode.info]
                self.currentInfo = self.reverseIndex[self.currentNode.info]


    def mouse2DownCanvas(self, event=None):
        #print 'in mouse2DownCanvas'
        pass



    def enter(self, event=None):
        try:
            if widgetsOnBackWindowsCanGrabFocus is False:
                lActiveWindow = self.draw.focus_get()
                if    lActiveWindow is not None \
                  and \
                    ( lActiveWindow.winfo_toplevel() != self.draw.winfo_toplevel() ):
                    return
        except:
            pass
        self.draw.focus_set()

    
    def createCanvas(self):
        """ create the Canvas and Title widgets"""
        self.widgetArea = Tkinter.Frame(self, borderwidth=2, relief='sunken')
        self.createMenus()
        self.canvasFrame = Tkinter.Frame(self)
        #if self.hasScroll:
        self.scrollregion = [0, 0, self.totalWidth, self.totalHeight]
        self.visibleOriginX = 0
        self.visibleOriginY = 0
        self.draw = Tkinter.Canvas(self.canvasFrame, width=self.visibleWidth,
                        height=self.visibleHeight,
                        background='white',
                        scrollregion=tuple(self.scrollregion))
        self.draw.scrollX = Tkinter.Scrollbar(self.canvasFrame,
                                        orient=Tkinter.HORIZONTAL)
        self.draw.scrollY = Tkinter.Scrollbar(self.canvasFrame,
                                        orient=Tkinter.VERTICAL)
        self.draw.bind('<Configure>', self.configure)

        self.draw['xscrollcommand'] = self.draw.scrollX.set
        self.draw['yscrollcommand'] = self.draw.scrollY.set
        self.draw.scrollX['command'] = self.myxview
        self.draw.scrollY['command'] = self.myyview
        self.canvasFrame.pack(side=Tkinter.LEFT,expand=1,fill=Tkinter.BOTH)
        self.draw.scrollX.pack(side=Tkinter.BOTTOM,fill=Tkinter.X)
        self.draw.scrollY.pack(side=Tkinter.LEFT,expand=1,fill=Tkinter.BOTH)
        self.draw.pack(side=Tkinter.LEFT,expand=1,fill=Tkinter.BOTH)
        border_w = self.draw.winfo_reqwidth() - self.visibleWidth
        border_h = self.draw.winfo_reqheight() - self.visibleHeight
        self.border = (border_w, border_h)


    def myxview(self, command, value):
        self.visibleOriginX = int(self.totalWidth*float(value))
        self.draw.xview(command, value)


    def myyview(self, command, value):
        self.visibleOriginY = int(self.totalHeight*float(value))
        self.draw.yview(command, value)


    def configure(self, event=None):
        new_width = event.width - self.border[0]
        new_height = event.height - self.border[1]
        width = int(self.draw.cget('width'))
        height = int(self.draw.cget('height'))
        if new_width==width and new_height==height:
            return
        self.draw.configure(width=new_width, height=new_height)
        self._setsize()
        self.visibleWidth = new_width
        self.visibleHeight = new_height
        self.clear()
        self.replot()

    

    def _setsize(self):
        self.width = int(self.draw.cget('width'))
        self.height = int(self.draw.cget('height'))
        #self.plotarea_size[0] = 0.97 * self.width
        #self.plotarea_size[1] = 0.97 * -self.height      
        #xo = 0.5*(self.width-self.plotarea_size[0])
        #yo = 0.5*(self.height-self.plotarea_size[1])
        #self.plotarea_origin = (xo, yo)
        self.plotarea_size = [0.96 * self.width, -0.90 * self.height]
        self.plotarea_origin = (0.030*self.width, 0.96*self.height)
        #self.plotarea_origin = (0.016*self.width, 0.96*self.height)


    def createMenus(self):
        self.mBar = Tkinter.Frame(self, relief=Tkinter.RAISED,borderwidth=2)
        self.mBar.pack(fill=Tkinter.X)
        self.menuButtons = {}
        self.makeFileMenu()
        self.makeEditMenu()
        apply(self.mBar.tk_menuBar, self.menuButtons.values())
        self.title = Tkinter.Label(self.mBar, text=self.name)
        self.title.pack(side=Tkinter.RIGHT)


    def makeFileMenu(self):
        File_button = Tkinter.Menubutton(self.mBar, text='File',underline=0)
        self.menuButtons['File'] = File_button
        File_button.pack(side = Tkinter.LEFT, padx='1m')
        File_button.menu = Tkinter.Menu(File_button)
        #File_button.menu.add_command(label='Load...', underline=0,
                                    #command = self.loadFile)
        File_button.menu.add_command(label='Exit...', underline=0,
                                    command = self.exit)
        File_button['menu'] = File_button.menu


    def loadFile(self, event=None):
        pass


    def exit(self, event=None):
        self.master.destroy()


    def makeEditMenu(self, event=None):
        Edit_button = Tkinter.Menubutton(self.mBar, text='Edit',underline=0)
        self.menuButtons['Edit'] = Edit_button
        Edit_button.pack(side = Tkinter.LEFT, padx='1m')
        Edit_button.menu = Tkinter.Menu(Edit_button)
        Edit_button.menu.add_command(label='Clear...', underline=0,
                                    command = self.clear)
        Edit_button.menu.add_command(label='Redraw...', underline=0,
                                    command = self.replot)
        #Edit_button.menu.add_command(label='Delete...', underline=0,
                                    #command = self.delete)
        Edit_button.menu.add_command(label='Write...', underline=0,
                                    command = self.printCanvas)
        Edit_button['menu'] = Edit_button.menu


    def clear(self, event=None):
        self.geoms = {}
        self.draw.delete('all')


    def delete(self, event=None):
        pass


    def replot(self):
        if self.last_drawn is not None:
            apply(self.buildIt, (self.last_drawn,), {})


    def printCanvas(self, idir=None, ifile=None, title=None):
        #other options are 'colormap', 'colormode', 'fontmap', 'height'
        #'pageanchor', 'pageheight','pagewidth', 'pagex','pagey','rotate'
        #'width', 'x', 'y'
        types = [('Postscript Files', '*.ps')]
        if idir==None:
            idir = os.curdir
        ifile = ifile
        if title==None:
            title='Postscript File'
        file = tkFileDialog.asksaveasfilename( filetypes=types,
                                           initialdir=idir,
                                           initialfile=ifile,
                                           title=title)
        self.draw.postscript({'file':file,'colormode':'color'})


    def buildIt(self, graphics, xaxis='automatic', yaxis='custom',
            xlabels='float', ylabels='int'):
        #NB: p1 is top left point of bounding box
        p1, p2 = graphics.boundingBox()
        xaxis = self._axisInterval(xaxis, p1[0], p2[0])
        yaxis = self._axisInterval(yaxis, p1[1], p2[1])
        text_width = [0., 0.]
        text_height = [0., 0.]
        if xaxis is not None:
            p1 = xaxis[0], p1[1]
            p2 = xaxis[1], p2[1]
            xticks = self._ticks(xaxis[0], xaxis[1], form=xlabels)
            bb = self._textBoundingBox(xticks[0][1])
            text_height[1] = bb[3]-bb[1]
            text_width[0] = 0.5*(bb[2]-bb[0])
            bb = self._textBoundingBox(xticks[-1][1])
            text_width[1] = 0.5*(bb[2]-bb[0])
        else:
            xticks = None
        if yaxis is not None:
            p1 = p1[0], yaxis[0]
            p2 = p2[0], yaxis[1]
            yticks = self._ticks(yaxis[0], yaxis[1], form=ylabels)
            for y in yticks:
                bb = self._textBoundingBox(y[1])
                w = bb[2]-bb[0]
                text_width[0] = max(text_width[0], w)
            h = 0.5*(bb[3]-bb[1])
            text_height[0] = h
            text_height[1] = max(text_height[1], h)
        else:
            yticks = None
        text1 = [text_width[0], -text_height[1]]
        text2 = [text_width[1], -text_height[0]]
        self.p1 = p1
        self.p2 = p2
        scale = ((self.plotarea_size[0]-text1[0]-text2[0]) / \
                 (p2[0]-p1[0]),
                 (self.plotarea_size[1]-text1[1]-text2[1]) / \
                 (p2[1]-p1[1]))
        shift = ((-p1[0]*scale[0]) + self.plotarea_origin[0] + \
                 text1[0],
                 (-p1[1]*scale[1]) + self.plotarea_origin[1] + \
                 text1[1])
        self._drawAxes(self.draw, xaxis, yaxis, p1, p2,
                        scale, shift, xticks, yticks)
        graphics.fitToScale(scale, shift, p1, p2)
        graphics.buildIcons(self) 
        if self.label_text:
            pt1 = self.plotarea_origin[0] + self.width/2
            pt2 = (self.height + self.plotarea_size[1])/2.0
            self.draw.create_text(pt1,pt2,text=self.label_text,anchor='s',font=self.font)
        if self.xlabel_text:
            pt1 = self.plotarea_origin[0] + self.width/2
            pt2 = self.plotarea_origin[1] - self.spacing
            self.draw.create_text(pt1,pt2,text=self.xlabel_text,anchor='n',font=self.font)
        if self.ylabel_text:
            pt1 = (self.plotarea_origin[0] - self.visibleOriginX)/3.0
            #pt1 = self.plotarea_origin[0]
            #pt1 = self.plotarea_origin[0] + self.spacing/8
            pt2 = self.plotarea_origin[1]/2
            self.draw.create_text(pt1,pt2,text=self.ylabel_text,anchor='w',font=self.font)
        self.last_drawn = graphics


    def _axisInterval(self, spec, lower, upper):
        if spec==None: return None
        if spec=='minimal':
            if lower==upper:
                return lower-0.5, upper+0.5
            else:
                return lower, upper
        if spec=='automatic' or 'custom':
            range = upper-lower
            if range==0:
                return lower-0.5, upper+0.5
            log = math.log10(range)
            power = math.floor(log)
            fraction = log - power
            if fraction<=0.05:
                power = power-1
            grid = 10.**power
            lower = lower - lower% grid
            mod = upper % grid
            if mod!=0:
                upper = upper - mod + grid
            if spec == 'custom':
                lower = 0.0
            return lower, upper
        if type(spec)==type(()):
            lower, upper = spec
            if lower <= upper:
                return lower, upper
            else: return upper, lower
        raise ValueError, str(spec)+': illegal axis specification'


    def _drawAxes(self, canvas, xaxis, yaxis,
                  bb1, bb2, scale, shift, xticks, yticks):
        dict = {'anchor': Tkinter.N, 'fill': 'black'}
        if self.font is not None:
            dict['font'] = self.font
        if xaxis is not None:
            lower, upper = xaxis
            text = 1
            for y, d in [(bb1[1], -3), (bb2[1], 3)]:
                p1 = (scale[0]*lower)+shift[0], (scale[1]*y)+shift[1]
                p2 = (scale[0]*upper)+shift[0], (scale[1]*y)+shift[1]
                Line(self.draw, p1[0], p1[1], p2[0], p2[1],
                     fill = 'black', width = 1)
                if xticks:
                    for x, label in xticks:
                        p = (scale[0]*x)+shift[0], \
                            (scale[1]*y)+shift[1]
                        Line(self.draw, p[0], p[1], p[0], p[1]+d,
                             fill = 'black', width = 1)
                        if text:
                            dict['text'] = label
                            apply(CanvasText, (self.draw, p[0],
                                               p[1]), dict)
                text = 0
        dict['anchor'] = Tkinter.E
        if yaxis is not None:
            lower, upper = yaxis
            text = 1
            for x, d in [(bb1[0], -3), (bb2[0], 3)]:
                p1 = (scale[0]*x)+shift[0], (scale[1]*lower)+shift[1]
                p2 = (scale[0]*x)+shift[0], (scale[1]*upper)+shift[1]
                Line(self.draw, p1[0], p1[1], p2[0], p2[1],
                     fill = 'black', width = 1)
                if yticks:
                    for y, label in yticks:
                        p = (scale[0]*x)+shift[0], \
                            (scale[1]*y)+shift[1]
                        Line(self.draw, p[0], p[1], p[0]-d, p[1],
                             fill = 'black', width = 1)
                        if text:
                            dict['text'] = label
                            apply(CanvasText,(self.draw,
                                              p[0]-2,p[1]), dict)
                text = 0


    def _ticks(self, lower, upper, form='float'):
        ideal = (upper-lower)/7.
        log = math.log10(ideal)
        power = math.floor(log)
        fraction = log-power
        factor = 1.
        error = fraction
        for f, lf in self._multiples:
            e = math.fabs(fraction-lf)
            if e < error:
                error = e
                factor = f
        grid = factor * 10.**power
        if power > 3 or power < -3:
            format = '%+7.0e'
        elif power >= 0:
            digits = max(1, int(power))
            format = '%' + `digits`+'.0f'
        else:
            digits = -int(power)
            format = '%'+`digits+2`+'.'+`digits`+'f'
        ticks = []
        t = -grid*math.floor(-lower/grid)
        while t <= upper and len(ticks) < 200:
            #ticks.append(t, format % (t,))
            if form=='float':
                ticks.append((t,str(t)))
            else:
                ticks.append((int(t),str(int(t))))
            t = t + grid
        return ticks


    def _textBoundingBox(self, text):
        bg = self.draw.cget('background')
        dict = {'anchor': Tkinter.NW, 'text': text, 'fill': bg}
        if self.font is not None:
            dict['font'] = self.font
        item = apply(CanvasText, (self.draw, 0., 0.), dict)
        bb = self.draw.bbox(item)
        self.draw.delete(item)
        return bb



class GraphPoints:


    def __init__(self, points, attr):
        self.points = points
        self.scaled = self.points
        self.attributes = {}
        for name, value in self._attributes.items():
            try:
                value = attr[name]
            except KeyError: pass
            self.attributes[name] = value


    def boundingBox(self):
        return minBound(self.points),  maxBound(self.points)


    def fitToScale(self, scale=(1,1), shift=(0,0)):
        self.scaled = []
        for x,y in self.points:
            self.scaled.append(((scale[0]*x)+shift[0],\
                               (scale[1]*y)+shift[1]))
        self.anchor = scale[1]*self.attributes.get('anchor', 0.0)+\
                      shift[1]



class GraphLine(GraphPoints):


    def __init__(self, points, **attr):
        GraphPoints.__init__(self, points, attr)

    _attributes = {'color':       'black',
                   'width':        1,
                   'smooth':       0,
                   'splinesteps': 12}

    def draw(self, canvas):
        color  = self.attributes['color']
        width  = self.attributes['width']
        smooth = self.attributes['smooth']
        steps  = self.attributes['splinesteps']
        arguments = (canvas,)
        if smooth:
            for i in range(len(self.points)):
                x1, y1 = self.scaled[i]
                arguments = arguments + (x1, y1)
        else:
            for i in range(len(self.points)-1):
                x1, y1 = self.scaled[i]
                x2, y2 = self.scaled[i+1]
                arguments = arguments + (x1, y1, x2, y2)
        apply(Line, arguments, {'fill': color, 'width': width,
                                'smooth': smooth, 'splinesteps':steps})



class InteractiveGraphObject:


    def __init__(self,  **kw):
        self.id = None
        if kw.has_key('name'):
            self.name = kw['name']
            del kw['name']
        else:
            self.name = 'NoName'
        if kw.has_key('editor'):
            self.editor = kw['editor']
            del kw['editor']
        else:
            self.editor = None
        self.selected = 0


    def editNoneMenu(self, type, *args, **kw):
        apply(eval('self.menu'+type), args, kw)


    def __repr__(self):
        return '<%s %s>'%(self.__class__.__name__, self.name)


    def rename(self, name):
        self.name=name


    def select(self):
        self.selected = 1
        self.editor.draw.addtag_withtag('selected', self.iconTag)


    def deselect(self):
        self.selected = 0
        self.editor.draw.dtag(self.iconTag, 'selected')


    def delete(self):
        self.editor.draw.delete(self.iconTag)


    def gettags(self):
        self.editor.draw.gettags(self.iconTag)
        
        
    def instantiate(self, editor):
        """method to be implmented by subclass
        this method should create the geometric primitive for this item
        and tag them all with self.iconTag
        define its uniqueTag
        set its tags
        """
        self.editor = editor
        self.menu = Tkinter.Menu(editor, title = self.name)
        self.menu.add_command(label='info', command= self.showInfo)
        self.buildNodeIcon(editor)


    def buildNodeIcon(self, editor):
        print 'in buildNodeIcon'



class GraphObjects:


    def __init__(self, objects):
        self.objects = objects


    def boundingBox(self):
        c1, c2 = self.objects[0].boundingBox()
        for object in self.objects[1:]:
            c1o, c2o = object.boundingBox()
            c1 = minBound([c1, c1o])

            c2 = maxBound([c2, c2o])
        return c1, c2


    def fitToScale(self, scale=(1,1), shift=(0,0)):
        for object in self.objects:
            object.fitToScale(scale, shift)


    def draw(self, canvas):
        for object in self.objects:
            object.draw(canvas)



class HistogramNode(InteractiveGraphObject):


    def __init__(self, width=10, height=40, name='NoName',info=None,
                        editor=None, **kw):
        kw['name'] = name
        kw['editor'] = editor
        apply(InteractiveGraphObject.__init__,(self,), kw)
        self.width = width
        self.height = height
        #info is a link to state info
        self.info = info
        self.point = (width, height)
        self.instantiate(editor, info)


    def instantiate(self, editor, info=None):
        """method to be implmented by subclass
        this method should create the geometric primitive for this item
        and tag them all with self.iconTag
        define its uniqueTag
        set its tags
        """
        self.editor = editor
        self.info = info
        self.menu = Tkinter.Menu(editor, title = self.name)
        self.menu.add_command(label='info', command= self.showInfo)


    def showInfo(self, event=None):
        print 'in showInfo with ', self.info


    def buildNodeIcon(self, editor, x1, y1, x2, y2, fill='red',
            width=1, outline='black', stipple = ''):
        canvas = editor.draw
        g = canvas.create_rectangle(x1, y1, x2, y2, fill=fill,
                width=width, outline=outline, stipple=stipple)
        editor.geoms[g] = self
        canvas.lift(g)
        self.id = g



class HistogramNodeSet:


    def __init__(self, nodes, **kw):
        points = []
        self.nodes = nodes
        for n in nodes:
            points.append(n.point)
        self.points = points
        self.scaled = points
        self.attributes = {}
        #anchor is y axis of histogram base
        #spread is width of each bar
        self.attributes['anchor'] = 0
        self.attributes['spread'] = 5
        self.attributes['fill'] = 'blue'
        self.attributes['stipple'] = ''
        self.attributes['outline'] = 'black'
        self.attributes['width'] = 1
        for name, value in kw.items():
            self.attributes[name] = value


    def boundingBox(self):
        return minBound(self.points), maxBound(self.points)


    def fitToScale(self, scale=(1,1), shift=(0,0), p1=(1,1), p2=(1,1)):
        self.scaled = []
        for x, y in self.points:
            self.scaled.append(((scale[0]*x) + shift[0],\
                                (scale[1]*y) + shift[1]))
    
        #need bb[0]
        newanchor = scale[1]*p1[1]+ shift[1]
        self.attributes['anchor'] = newanchor


    def buildIcons(self, editor):
        fill = self.attributes['fill']
        width = self.attributes['width']
        stipple = self.attributes['stipple']
        outline = self.attributes['outline']
        spread = self.attributes['spread']
        anchor = self.attributes['anchor']
        canvas = editor.draw
        argument = (canvas,)
        p1, p2 = self.boundingBox()

        for i in range(len(self.nodes)):
            node = self.nodes[i]
            x1, y1 = self.scaled[i]
            node.buildNodeIcon(editor, x1-spread, y1, x1+spread,
                    anchor, fill=fill, width=width,
                    outline=outline, stipple=stipple)


                
class InteractiveHistogramGraph(InteractiveGraphBuilder):


    def __init__(self, name='InteractiveHistogramBuilder', master=None, 
               nodeList=[], yoffset=200, maxwidth=200, 
               reverseIndex=None, label_text=None,
               xlabel_text=None, ylabel_text=None, **kw):
        apply( InteractiveGraphBuilder.__init__,(self, name, master), kw)
        #then instantiate the nodes in nodeList
        #nodeList is list of tuples-> (width, height, info)
        self.spacing = 5
        self._multiples = [(2., math.log10(2.)), (5., math.log10(5.))] 
        self.reverseIndex = reverseIndex

        num = 0
        nodes = []
        #for w, h, info  in nodeList:
        for w, h  in nodeList:
            newname = '%d'%num
            newnode = HistogramNode(width=w, height=h,
                name=newname, info=num)
            self.nodes[newname] = newnode
            num = num + 1
        
        self.nodeSet = HistogramNodeSet(self.nodes.values())
        self._setsize()
        self.label_text = label_text
        self.xlabel_text = xlabel_text
        self.ylabel_text = ylabel_text
        self.buildIt(self.nodeSet, xaxis='automatic', yaxis='custom',
                ylabels='int')


if __name__== '__main__':
    root = Tkinter.Tk()
    #nodeList  = [(0,0,'origin'),(32, 20, 'first'), (15, 35, 'second'),(20, 10,'third')]
    #nodeList  = [(0,0),(32, 20 ), (15, 35 ),(20, 10,)]
    #histNB = InteractiveHistogramGraph('Test', master=root, nodeList=nodeList)
    import energy_histogram
    ehist = energy_histogram.EnergyHistogram()
    ehist.get_input_from_energies('ebind_top.hist')
    nodeList = ehist.histogram.array
    ehist.createReverseIndex()
    reverseIndex = ehist.reverseIndex
    histNB = InteractiveHistogramGraph('test', master=root, nodeList=nodeList,
        reverseIndex = reverseIndex)
    

#from InteractiveGraph import InteractiveHistogramGraph
#import energy_histogram
#ehist = energy_histogram.EnergyHistogram()
#ehist.get_input_from_energies('ebind_top.hist')
#nodeList = ehist.histogram.array
#import Tkinter
#root = Tkinter.Tk()
#histNB = InteractiveHistogramGraph('test', master=root, nodeList=nodeList)



