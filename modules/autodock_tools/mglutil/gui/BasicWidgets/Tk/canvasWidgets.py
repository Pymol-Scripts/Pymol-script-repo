#############################################################################
#
# Author: Alexandre T. GILLET
#
# Copyright: TSRI 2003
#
#############################################################################

import types
import Tkinter, Pmw, os
from mglutil.util.callback import CallBackFunction, CallbackManager




class PatternCanvas:
    
    def __init__(self,master,title=None,callback=None,
                 size=5,border=1,tkCol=None,immediate=1):
        #tkCol is a list of list where each elemnt is a color
        # tkCol[0][0] = a Tk color
        # size = size of the sqare representing the Tkcolor on the canvas
        # border = border size
        # 
        if not master:
            master = Tkinter.Toplevel()

        if title is not None:
            master.title(title)
        
        f = self.frame = Tkinter.Frame(master)

        
        if tkCol:
            self.Xelemt = len(tkCol)
            self.Yelemt = len(tkCol[0])
            self.size = size
            self.border = border
            self.width =(self.Xelemt*self.size+2*self.border)
            self.height =(self.Yelemt*self.size+2*self.border)
            self.cwcanvas = Tkinter.Canvas(f,width=self.width,
                                           height=self.height, 
                                           borderwidth=3 )
                                           
            xo = self.border
            for i in range(self.Xelemt):
                xe = xo + self.size
                yo = self.border
                for j in range(self.Yelemt):
                    ye = yo + self.size
                    self.cwcanvas.create_rectangle(xo, yo, xe, ye,
                                                   fill=tkCol[i][j])
                    yo = ye
                xo = xe

        self.cwcanvas.pack()
        #self.frame.pack()

        # CallBack Manager
        self.cbManager = CallbackManager()
        if callback:
            if type(callback) in [types.ListType, types.TupleType]:
                map(self.cbManager.AddCallback, callback)
            else:
                self.cbManager.AddCallback(callback)


    def get(self):
        pass
    
    def set(self):
        pass
    
    def pack(self,*args, **kw):
        apply(self.frame.pack, args, kw)
        
    def pack_forget(self,*args, **kw):
        apply(self.frame.pack_forget, args, kw)
                         
    def grid(self,*args, **kw):
        apply(self.frame.grid, args, kw)

    def grid_forget(self,*args, **kw):
        apply(self.frame.grid_forget, args, kw)

    def destroy(self,*args, **kw):
        apply(self.frame.destroy, args, kw)

class CoefCanvas:
    
    def __init__(self,master,title=None,width=150,height=50,
                 callback=None):
        #tkCol is a list of list where each elemnt is a color
        # tkCol[0][0] = a Tk color
        # size = size of the sqare representing the Tkcolor on the canvas
        # border = border size
        # 
        if not master:
            master = Tkinter.Toplevel()

        if title is not None:
            master.title(title)
        
        f = self.frame = Tkinter.Frame(master)

        
        self.coeff = 0.0
        self.value = "%.1f"%self.coeff+'%'
        self.width = width
        self.height = height
        
        self.cwcanvas = Tkinter.Canvas(f,width=self.width,
                                       height=self.height)
        self.wrect = self.cwcanvas.create_rectangle(0,0,100,10,
                                                    outline='black',
                                                    fill='white')
        self.rrect = self.cwcanvas.create_rectangle(0,0,self.coeff,10,
                                                    outline='black',
                                                    fill='red')
        self.labelvalue = self.cwcanvas.create_text(120,10,text=self.value)

        self.cwcanvas.pack(side='top')
        #self.frame.pack()

        # CallBack Manager
        self.cbManager = CallbackManager()
        if callback:
            if type(callback) in [types.ListType, types.TupleType]:
                map(self.cbManager.AddCallback, callback)
            else:
                self.cbManager.AddCallback(callback)


    def get(self):
        pass
    
    def set(self,coeff,color):

        self.cwcanvas.coords(self.rrect,0,0,coeff,10)
        self.cwcanvas.itemconfigure(self.rrect,fill=color)
        self.cwcanvas.itemconfigure(self.labelvalue,text="%.1f"%coeff+'%')
    
    def pack(self,*args, **kw):
        apply(self.frame.pack, args, kw)
        
    def pack_forget(self,*args, **kw):
        apply(self.frame.pack_forget, args, kw)
                         
    def grid(self,*args, **kw):
        apply(self.frame.grid, args, kw)

    def grid_forget(self,*args, **kw):
        apply(self.frame.grid_forget, args, kw)

    def destroy(self,*args, **kw):
        apply(self.frame.destroy, args, kw)
