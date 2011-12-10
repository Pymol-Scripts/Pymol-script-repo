
import tkColorChooser
import Tkinter
T = Tkinter

from mglutil.util.callback import CallBackFunction
from mglutil.gui.BasicWidgets.Tk.thumbwheel import ThumbWheel
##  from ViewerFramework.gui import InputFormDescr, InputForm
from mglutil.gui.InputForm.Tk.gui import InputFormDescr, InputForm

class IntThumbWheel(ThumbWheel):

    def __init__(self, master=None, width=200, height=40,
                 nblines=40):
        ThumbWheel.__init__(self, master, width, height, nblines)
        self.showLabelInWidget = 1

    def get(self):
        return int(self.val)


class MovableWidget:

    def __init__(self, widget, posx=0, posy=0, gridSize=1):
        self.widget = widget
        self.widget.place(x=posx, y=posy)
        self.widget.bind('<Enter>', self.enter)
        self.widget.bind('<Button-3>', self.showMenu)
        #self.widget.bind('<Leave>', self.leave)
        # Creating the frame that will resize when resizing the widget.
        self.resizeDraw = T.Frame(widget.master,
                                  relief='ridge', background='red')
##          self.resizeDraw.bind('<B2-Motion>', self.resizeTop)
##          self.resizeDraw.bind('<ButtonRelease-2>', self.endResize)
        self.posx = posx
        self.posy = posy
        self.posxg = posx
        self.posyg = posy
        self.gs = gridSize
        self.optionWidgets = {
            'activebackground':self.setColor,
            'activeforeground':self.setColor,
            'background':self.setColor,
            'foreground':self.setColor,
            'disabledbackground':self.setColor,
            'disabledforeground':self.setColor,
            'highlightbackground':self.setColor,
            'highlightforeground':self.setColor,
            'highlightcolor':self.setColor,
            'borderwidth':self.getInt,
       
            }
        self.cc = tkColorChooser.Chooser()
        self.buildMenu()
        

    def setColor(self):
        return self.cc.show()[1]


    def getInt(self):
        idf = InputFormDescr(title='Choose a value')
        idf.append({'widgetType':IntThumbWheel,
                    'name':'tw',
                    'wcfg':{'width':125, 'height':30, 'nblines':30}})
        form  = InputForm(master=self.widget.master, root = None, descr=idf)
        values = form.go()
        return values['tw']

    
    def setOption_cb(self, option, event=None):
        value = self.optionWidgets[option]()
##          print option, value
        apply( self.widget.configure, (), {option:value} )


    def buildMenu(self):
        self.menu = T.Menu(self.widget.master)
        for k in self.widget.keys():
            if k in self.optionWidgets.keys():
                cb = CallBackFunction(self.setOption_cb, k)
                self.menu.add_command(label=k, command=cb)

    def showMenu(self, event):
        self.menu.post(event.x_root, event.y_root)
        

    def place(self, widget, x=None, y=None, w=None, h=None, **kw):
        gs = self.gs
        if x: kw['x'] = (x/gs)*gs
        if y: kw['y'] = (y/gs)*gs
        if w: kw['width'] = w
        if h: kw['height'] = h
##          print kw['x'], kw['y'], self.posx, self.posy
        apply( widget.place, (), kw)

        
    def enter(self, event=None):
        self.widget.focus_set()
        self.widget.bind('<ButtonPress-2>', self.moveOnButton)


    def moveOnButton(self, event=None):
        x0 = self.posx
        y0 = self.posy
        w = self.widget.winfo_width()
        h = self.widget.winfo_height()
        self.origx = x = event.x
        self.origy = y = event.y
        if y>0 and y<10:
            self.widget.place_forget()
            self.resizeDraw.configure(cursor='top_side')
            self.place(self.resizeDraw, self.posx, self.posy, w, h)
            self.widget.update_idletasks()
            self.resizeDraw.grab_set()
            self.resizeDraw.bind('<B2-Motion>', self.resizeTop)
            self.resizeDraw.bind('<ButtonRelease-2>', self.endResize)

        elif y<h and y>h-10:
            self.resizeDraw.configure(cursor='bottom_side')
            self.widget.place_forget()
            self.place(self.resizeDraw, self.posx, self.posy, w, h)
            self.widget.update_idletasks()
            self.resizeDraw.grab_set()
            self.resizeDraw.bind('<B2-Motion>', self.resizeBottom)
            self.resizeDraw.bind('<ButtonRelease-2>', self.endResize)

        elif x>0 and x<10:
            self.resizeDraw.configure(cursor='left_side')
            self.widget.place_forget()
            self.place(self.resizeDraw, self.posx, self.posy, w, h)
            self.widget.update_idletasks()
            self.resizeDraw.grab_set()
            self.resizeDraw.bind('<B2-Motion>', self.resizeLeft)
            self.resizeDraw.bind('<ButtonRelease-2>', self.endResize)

        elif x<w and x>w-10:
            self.resizeDraw.configure(cursor='right_side')
            self.widget.place_forget()
            self.place(self.resizeDraw, self.posx, self.posy, w, h)
            self.widget.update_idletasks()
            self.resizeDraw.grab_set()
            self.resizeDraw.bind('<B2-Motion>', self.resizeRight)
            self.resizeDraw.bind('<ButtonRelease-2>', self.endResize)

        else:
            self.widget.configure(cursor='')


    def resizeTop(self, event=None):
        dy = event.y-self.origy
        self.posy = self.posy+dy
        w = self.resizeDraw.winfo_width()
        h = self.resizeDraw.winfo_height()        
        self.place(self.resizeDraw, self.posx, self.posy, w, h-dy)

    def resizeLeft(self, event=None):
        dx = event.x-self.origx
        self.posx = self.posx+dx
        w = self.resizeDraw.winfo_width()
        h = self.resizeDraw.winfo_height()
        self.place(self.resizeDraw, self.posx, self.posy, w-dx, h)

    def resizeRight(self, event=None):
        dx = event.x-self.origx
        w = self.resizeDraw.winfo_width()
        h = self.resizeDraw.winfo_height()        
        self.place(self.resizeDraw, self.posx, self.posy, w+dx, h)
        self.origx = self.origx+dx
        
    def resizeBottom(self, event=None):
        dy = event.y-self.origy
        w = self.resizeDraw.winfo_width()
        h = self.resizeDraw.winfo_height()        
        self.place(self.resizeDraw, self.posx, self.posy, w, h+dy)
        self.origy = self.origy+dy
        
    def endResize(self, event=None):
        self.resizeDraw.grab_release()
        w = self.resizeDraw.winfo_width()
        h = self.resizeDraw.winfo_height()
        self.resizeDraw.place_forget()
        self.widget.place(x=self.posx, y=self.posy, width=w, height=h)
        
    def moveRight(self, event=None):
        event.widget.place_forget()
        self.posx = self.posx+1
        event.widget.place(x=self.posx, y=self.posy)
        
    def moveLeft(self, event=None):
        event.widget.place_forget()
        self.posx = self.posx-1
        event.widget.place(x=self.posx, y=self.posy)

    def moveUp(self, event=None):
        event.widget.place_forget()
        self.posy = self.posy-1
        event.widget.place(x=self.posx, y=self.posy)

    def moveDown(self, event=None):
        event.widget.place_forget()
        self.posy = self.posy+1
        event.widget.place(x=self.posx, y=self.posy)

    def recordPosition(self, event=None):
        self.x0 = event.x
        self.y0 = event.y
##          print 'Origin', self.x0, self.y0
#        self.widget.update_idletasks()
        self.widget.grab_set()

    def moveButton(self, event):
##          print 'event', event.x, event.y
        dx = event.x - self.x0
        dy = event.y - self.y0
        self.widget.place_forget()
        self.posx = self.posx + dx
        self.posy = self.posy + dy
        gs = self.gs


        w = self.widget.winfo_width()
        h = self.widget.winfo_height()
        self.widget.place(x=self.posx, y=self.posy, width = w, height = h)


##      def moveButton(self, event):
##          print 'event', event.x, event.y
##          dx = event.x - self.x0
##          dy = event.y - self.y0
##          self.widget.place_forget()
##          self.posx = self.posx + dx
##          self.posy = self.posy + dy
##          gs = self.gs

##          posxg = round((self.posx/float(gs)))*gs
##          dx = self.posxg - posxg
##          self.x0 = self.x0 + dx
        
##          posyg = round((self.posy/float(gs)))*gs
##          dy = self.posyg - posyg
##          self.y0 = self.y0 + dy

##          w = self.widget.winfo_width()
##          h = self.widget.winfo_height()
##          self.widget.place(x=self.posxg, y=self.posyg, width = w, height = h)

    def endMove(self, event=None):
        self.widget.grab_release()
