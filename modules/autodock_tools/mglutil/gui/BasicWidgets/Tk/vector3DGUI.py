## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#########################################################################
#
# Date: Nov 2001 Authors: Michel Sanner, Daniel Stoffler
#
#    sanner@scripps.edu
#    stoffler@scripps.edu
#
# Copyright: Michel Sanner, Daniel Stoffler and TSRI
#
#########################################################################

import Tkinter, numpy.oldnumeric as Numeric, string, math
    
from mglutil.math.rotax import rotax
from mglutil.gui.Misc.Tk.KeybdModMonitor import KeyboardModifierMonitor
from mglutil.util.callback import CallbackManager
from thumbwheel import ThumbWheel
from optionsPanel import VectorOptionsPanel

class vectorGUI(Tkinter.Frame, KeyboardModifierMonitor):

    """ This class implements a vector widget.
    The widget has a vector which can be moved within a sphere to generate
    a 3D vector. Values are normalized and stored in self.vector
    In addition, the vector can be rotated with 3 thumbwheels.
    Values can be entered directly by typing them into the 3 entry forms.
    Then, the 'normalize and set' button has to be pressed in order to
    normalize and set the new vector.

    The widget has a configure() method: vector, mode, precision and
    continuous can be set this way.
    vector is a list of 3 floating values, e.g. [0., 0., 1.]
    mode describes the axis movement (rotation around an axis): is type
    string and can be either 'X', 'Y' or 'Z'. Free movement (standard
    value) is 'XY'.
    continuous can be either 0 (or None) or 1. Default is 0
    precision is type int and ranges from 1 - 10
    master, name and size can be passed only to the constructor.

    a lock() method is used to disable the various gui components of the
    options panel. Usage: <instance>.lock(<component>=<value>)
    component is continuous, precision or mode. value is 0 or 1. 1 disables,
    0 enables.
    """
    def __init__(self, master=None, name='vector', size=200, continuous = 1,
                 vector=[0.0, 0.0, 1.0], mode='XY', precision=5,
                 lockContinuous=0,  lockPrecision=0, lockMode=0,
                 callback=None, labelSide='top'):
        
	KeyboardModifierMonitor.__init__(self)

        self.callback = callback # user specified callback
        self.name=name             # title inside canvas
        self.labelSide=labelSide   # where title gets packed
	self.mode=mode             # axe mode: can be 'XY', 'X', 'Y' or 'Z'
        self.precision=precision   # floating number digits
        self.continuous=continuous # can be 1 or 0
        self.vector=vector         # initial vector value
        self.size=size             # size of vector widget

        self.lockContinuous = lockContinuous  # set to 1 to lock menus in
                                              # option panel
        self.lockPrecision = lockPrecision
        self.lockMode = lockMode

        self.r = self.size/2
        self.r2 = self.r*self.r

        self.drawShadowX = 0
        self.drawShadowY = 1
        self.drawShadowZ = 0
	self.fillShadowPlanes = 1

	Tkinter.Frame.__init__(self, master)
        Tkinter.Pack.config(self)

        self.callbacks = CallbackManager() # object to manage callback
                                        # functions. They get called with the
                                        # current value as an argument
        self.zeros = Numeric.array( (0,0,0), 's')
        self.viewingMatInv = Numeric.array(
            [[  0.96770716, -0.03229283, -0.25      ,  0.        ],
             [  0.03229283, -0.96770716,  0.25      ,  0.        ],
             [  0.25      ,  0.25      ,  0.93541437,  0.        ],
             [  0.        ,  0.        ,  0.        ,  1.        ]],'f')
        self.viewingMat = Numeric.transpose(self.viewingMatInv)
        self.createCanvas(master, size)
        self.createEntries(self.frame)
	Tkinter.Widget.bind(self.canvas, "<ButtonPress-1>", self.mouseDown)
	Tkinter.Widget.bind(self.canvas, "<ButtonRelease-1>", self.mouseUp)
	Tkinter.Widget.bind(self.canvas, "<B1-Motion>", self.mouseMove)

        self.setEntries()

        self.opPanel = VectorOptionsPanel(master = self,
                                          title="Vector GUI Options")
        Tkinter.Widget.bind(self.canvas, "<Button-3>", self.toggleOptPanel)

        if self.callback:
            self.callbacks.AddCallback(self.callback)


    def toggleOptPanel(self, event=None):
        # opens and closes options panel by right clicking on widget
        if self.opPanel.flag:
           self.opPanel.Dismiss_cb()
        else:
            if not hasattr(self.opPanel, 'optionsForm'):
                self.opPanel.displayPanel(create=1)
            else:
                self.opPanel.displayPanel(create=0)
                

    def mouseUp(self, event):
        if not self.continuous:
            self.callbacks.CallCallbacks(self.vector)

	
    def mouseDown(self, event):
	# remember where the mouse went down
        xc = event.x - self.xm
        yc = self.ym - event.y
        # compute the intersection point between
        z2 = self.r2-(xc*xc)-(yc*yc)

        if z2>=0: # we picked inside the sphere. going for a XY rotation
            self.lastPt3D = (xc, yc, math.sqrt(z2))
        else: # going for a Z rotation
            pass


    def mouseMove(self, event):
        # simple trackball, only works inside cirle
        # creates an XY rotation defined by pts intersecting the spheres
        xc = event.x - self.xm
        yc = self.ym - event.y
        # compute the intersection point between
        xc2 = xc*xc
        yc2 = yc*yc
        z2 = self.r2-xc2-yc2

        if z2 < 0:
            lInvMag = 1./math.sqrt(xc2 + yc2)
            xc *= lInvMag * (self.r)
            yc *= lInvMag * (self.r)
            z2 = 0

        # compute rotation angle
	a = self.lastPt3D
	b = (xc, yc, math.sqrt(z2))
        ang = math.acos((a[0]*b[0]+a[1]*b[1]+a[2]*b[2])/self.r2)
        if self.mode=='XY':
	    #compute rotation axis
            rotaxis = Numeric.array( (a[1]*b[2] - a[2]*b[1],
				  a[2]*b[0] - a[0]*b[2],
				  a[0]*b[1] - a[1]*b[0] ), 'f' )
        elif self.mode=='X': rotaxis = Numeric.array( (1.,0.,0.), 'f')
        elif self.mode=='Y': rotaxis = Numeric.array( (0.,1.,0.), 'f')
        elif self.mode=='Z': rotaxis = Numeric.array( (0.,0.,1.), 'f')
        mat = rotax( self.zeros, rotaxis, ang )
        self.lastPt3D = b
        self.updateVector(mat)


    def updateVector(self, mat):
        mat = Numeric.reshape(mat, (4,4))
        newPts = self.vector + [1]
        newPts = Numeric.dot( [newPts], mat )[0]
        self.vector = [newPts[0], newPts[1], newPts[2]]
        self.setEntries()
        self.drawVector()
        if self.continuous:
            self.callbacks.CallCallbacks(self.vector)
            
        
    def drawVector(self):
        coords3D = self.vector + [1]
	# apply viewing transformation to vector
        newPtsWithView = Numeric.dot( [coords3D],
                                                 self.viewingMat)[0]
	# compute 2D projection of vector (broken on 2 segments for
	# depth cueing
        x1 = self.xm+int(newPtsWithView[0]*(self.xm))
        y1 = self.ym+int(newPtsWithView[1]*(self.ym))

	# change vector's segments coordinates
        self.canvas.coords(self.lId1, self.xm, self.ym, x1, y1)

	# update vector shadows
	# Y=0 plane
	if self.drawShadowY:
	    pt = [coords3D[0], 0, coords3D[2], 1.]
	    newPtsWithView = Numeric.dot( [pt], self.viewingMat)[0]
	    xm = self.xm+int(newPtsWithView[0]*(self.xm))
	    ym = self.ym+int(newPtsWithView[1]*(self.ym))
	    if self.fillShadowPlanes:
		self.canvas.coords(self.shadowPY,self.xm,self.ym,xm,ym,x1,y1)
	    self.canvas.coords(self.shadowY,self.xm,self.ym,xm,ym,x1,y1)

	# X=0 plane
	if self.drawShadowX:
	    pt = [0, coords3D[1], coords3D[2], 1.]
	    newPtsWithView = Numeric.dot( [pt], self.viewingMat)[0]
	    xm = self.xm+int(newPtsWithView[0]*(self.xm))
	    ym = self.ym+int(newPtsWithView[1]*(self.ym))
	    if self.fillShadowPlanes:
		self.canvas.coords(self.shadowPX,self.xm,self.ym,xm,ym,x1,y1)
	    self.canvas.coords(self.shadowX, self.xm, self.ym, xm, ym, x1,y1)

	# Z=0 plane
	if self.drawShadowZ:
	    pt = [coords3D[0], coords3D[1], 0, 1.]
	    newPtsWithView = Numeric.dot( [pt], self.viewingMat)[0]
	    xm = self.xm+int(newPtsWithView[0]*(self.xm))
	    ym = self.ym+int(newPtsWithView[1]*(self.ym))
	    if self.fillShadowPlanes:
		self.canvas.coords(self.shadowPZ,self.xm,self.ym,xm,ym,x1,y1)
	    self.canvas.coords(self.shadowZ, self.xm, self.ym, xm, ym, x1,y1)

	if self.vector[0]<0.0:
	    self.canvas.tag_raise('verticalCircle', 'moving')
	else:
	    self.canvas.tag_lower('verticalCircle', 'moving')

	if self.vector[1]<0.0:
	    self.canvas.tag_raise('horizontalCircle', 'moving')
	else:
	    self.canvas.tag_lower('horizontalCircle', 'moving')

	if self.vector[2]<0.0 or self.vector[1]<0.0:
	    self.canvas.tag_raise('axis', 'moving')
	else:
	    self.canvas.tag_lower('axis', 'moving')


    def thumbx_cb(self, events=None):
        val=self.thumbx.value

##          valX=self.thumbx.value
##          valY=self.thumby.value
##          valZ=self.thumbz.value

##          n = math.sqrt(valX*valX+valY*valY+valZ*valZ)
##          if n == 0.0: v = [0.0, 0.0, 1.0]
##          else: v = [valX/n, valY/n, valZ/n]
##          val = v[0]

	rot = Numeric.zeros( (4,4), 'f' )
        rot[0][0] = 1.0
        rot[1][1] = math.cos(val)
        rot[1][2] = -math.sin(val)
        rot[2][1] = math.sin(val)
        rot[2][2] = math.cos(val)
        self.updateVector(rot)

    def thumby_cb(self, events=None):
        val=self.thumby.value
	rot = Numeric.zeros( (4,4), 'f' )
        rot[0][0] = math.cos(val)
        rot[0][2] = -math.sin(val)
        rot[1][1] = 1.0
        rot[2][0] = math.sin(val)
        rot[2][2] = math.cos(val)
        self.updateVector(rot)


    def thumbz_cb(self, events=None):
        val=self.thumbz.value
	rot = Numeric.zeros( (4,4), 'f' )
        rot[0][0] = math.cos(val)
        rot[0][1] = -math.sin(val)
        rot[1][0] = math.sin(val)
        rot[1][1] = math.cos(val)
        rot[2][2] = 1.0
        self.updateVector(rot)


    def entryX_cb(self, event=None):
        val = self.entryXTk.get()
        if len(val) == 0: val = self.vector[0]
        try:
            val = float(val)
            self.entryXTk.set(self.thumbx.labelFormat%val)
        except ValueError:
            # put back original value if someone types garbage
            self.entryXTk.set(self.thumbx.labelFormat%self.vector[0])


    def entryY_cb(self, event=None):
        val = self.entryYTk.get()
        if len(val) == 0: val = self.vector[1]
        try:
            val = float(val)
            self.entryYTk.set(self.thumby.labelFormat%val)
        except ValueError:
            # put back original value if someone types garbage
            self.entryYTk.set(self.thumby.labelFormat%self.vector[1])


    def entryZ_cb(self, event=None):
        val = self.entryZTk.get()
        if len(val) == 0: val = self.vector[2]
        try:
            val = float(val)
            self.entryZTk.set(self.thumbz.labelFormat%val)
        except ValueError:
            # put back original value if someone types garbage
            self.entryZTk.set(self.thumbz.labelFormat%self.vector[2])


    def entryV_cb(self, event=None):
        v = self.entryVTk.get()
        try: val = string.split(v)
        except: 
            self.setEntries()
            return

        if val is None or len(val)!= 3:
            self.setEntries()
            return

        try:
            valX = float(val[0])
            valY = float(val[1])
            valZ = float(val[2])
        except:
            self.setEntries()
            return

        # compute normalized vector
        n = math.sqrt(valX*valX+valY*valY+valZ*valZ)
        if n == 0.0: v = [0.0, 0.0, 1.0]
        else: v = [valX/n, valY/n, valZ/n]
        self.vector = v
        self.setEntries()
        self.drawVector()
        if self.continuous:
            self.callbacks.CallCallbacks(self.vector)

    def setButton_cb(self, event=None):
        valX = float(self.entryXTk.get())
        valY = float(self.entryYTk.get())
        valZ = float(self.entryZTk.get())

        # compute normalized vector
        n = math.sqrt(valX*valX+valY*valY+valZ*valZ)
        if n == 0.0: v = [0.0, 0.0, 1.0]
        else: v = [valX/n, valY/n, valZ/n]
        self.vector = v
        self.setEntries()
        self.drawVector()
        if self.continuous:
            self.callbacks.CallCallbacks(self.vector)
        

    def createEntries(self, master):
        self.f = Tkinter.Frame(master)
	self.f.grid(column=3, rowspan=3)

        def fX(): self.vector = [1.,0.,0.]; self.setEntries(); self.callbacks.CallCallbacks(self.vector)
        def fY(): self.vector = [0.,1.,0.]; self.setEntries(); self.callbacks.CallCallbacks(self.vector)
        def fZ(): self.vector = [0.,0.,1.]; self.setEntries(); self.callbacks.CallCallbacks(self.vector)
        lX = Tkinter.Button(master=self.f, text='x', command=fX)
        lY = Tkinter.Button(master=self.f, text='y', command=fY)
        lZ = Tkinter.Button(master=self.f, text='z', command=fZ)
        lX.grid(row=0, column=0)
        lY.grid(row=1, column=0)
        lZ.grid(row=2, column=0)

        self.thumbx = ThumbWheel(master=self.f, width=50,
                                 height=20, labcfg={'text':'X:','side':'left'},
                                 wheelPad=2, oneTurn=.1, min=-1, max=1,
                                 showLabel=0, precision=5, type=float)
        self.thumbx.callbacks.AddCallback(self.thumbx_cb)
        self.thumbx.unbind("<Button-3>")
        self.thumbx.canvas.unbind("<Button-3>")
        self.thumbx.grid(row=0, column=1)

        self.thumby = ThumbWheel(master=self.f, width=50,
                                 height=20, labcfg={'text':'Y:','side':'left'},
                                 wheelPad=2, oneTurn=.1, min=-1, max=1,
                                 showLabel=0, precision=5, type=float)
        self.thumby.callbacks.AddCallback(self.thumby_cb)
        self.thumby.unbind("<Button-3>")
        self.thumby.canvas.unbind("<Button-3>")
        self.thumby.grid(row=1, column=1)

        self.thumbz = ThumbWheel(master=self.f, width=50,
                                 height=20, labcfg={'text':'Z:','side':'left'},
                                 wheelPad=2, oneTurn=.1, min=-1, max=1,
                                 showLabel=0, precision=5, type=float)
        self.thumbz.callbacks.AddCallback(self.thumbz_cb)
        self.thumbz.unbind("<Button-3>")
        self.thumbz.canvas.unbind("<Button-3>")
        self.thumbz.grid(row=2, column=1)

        self.entryXTk = Tkinter.StringVar()
        self.entryX = Tkinter.Entry(master=self.f, textvariable=self.entryXTk,
                                    width=8)
        self.entryX.bind('<Return>', self.entryX_cb)
        self.entryX.grid(row=0, column=2)

        self.entryYTk = Tkinter.StringVar()
        self.entryY = Tkinter.Entry(master=self.f, textvariable=self.entryYTk,
                                    width=8)
        self.entryY.bind('<Return>', self.entryY_cb)
        self.entryY.grid(row=1, column=2)

        self.entryZTk = Tkinter.StringVar()
        self.entryZ = Tkinter.Entry(master=self.f, textvariable=self.entryZTk,
                                    width=8)
        self.entryZ.bind('<Return>', self.entryZ_cb)
        self.entryZ.grid(row=2, column=2)

        self.entryVTk = Tkinter.StringVar()
        self.entryV = Tkinter.Entry(master, textvariable=self.entryVTk,
                                    width=18)
        
        self.entryV.bind('<Return>', self.entryV_cb)
        
        self.f.pack(side='top', expand=1)

        self.entryV.pack()
        
	self.setButton=Tkinter.Button(master, text='normalize and set',
                                      command = self.setButton_cb)
	self.setButton.pack(side='bottom')


    def setEntries(self):
        self.entryXTk.set(self.thumbx.labelFormat%self.vector[0])
        self.entryYTk.set(self.thumby.labelFormat%self.vector[1])
        self.entryZTk.set(self.thumbz.labelFormat%self.vector[2])

        lf = '%.3f'
        self.entryVTk.set(lf%self.vector[0]+' '+lf%self.vector[1]+' '\
                          +lf%self.vector[2])
        self.drawVector()


    def createCanvas(self, master, size=200):

        self.frame = Tkinter.Frame(self, relief = 'sunken', borderwidth=5)

        if self.name is not None:
            self.title = Tkinter.Label(self.frame, text=self.name)
            self.title.pack(side=self.labelSide)

	self.canvas = Tkinter.Canvas(self.frame, width=size, height=size)

        # set the focus so that we get keyboard events, and add callbacks
        self.canvas.bind('<KeyPress>', self.modifierDown)
        self.canvas.bind("<KeyRelease>", self.modifierUp)

        xm = self.xm = ym = self.ym = self.r
	self.canvas.create_oval(0, 0, size, size)
	self.canvas.create_oval(xm-(xm/4), 0, xm+(xm/4), size,
				tags='verticalCircle')
	self.canvas.create_oval(0, ym-(ym/4), size, ym+(ym/4),
				tags='horizontalCircle')

        # apply viewing transformation to vector
        XaxisWithView = Numeric.dot([(1.,0.,0.,1.)],self.viewingMat)[0]
        x1 = self.xm+int(XaxisWithView[0]*(self.xm))
        y1 = self.ym+int(XaxisWithView[1]*(self.ym))
        self.canvas.create_line(xm, ym, x1, y1, fill='red', tags='axis')

        XaxisWithView = Numeric.dot([(0.,1.,0.,1.)],self.viewingMat)[0]
        x2 = self.xm+int(XaxisWithView[0]*(self.xm))
        y2 = self.ym+int(XaxisWithView[1]*(self.ym))
        self.canvas.create_line(xm, ym, x2, y2, fill='green', tags='axis')

        XaxisWithView = Numeric.dot([(0.,0.,1.,1.)],self.viewingMat)[0]
        x3 = self.xm+int(XaxisWithView[0]*(self.xm))
        y3 = self.ym+int(XaxisWithView[1]*(self.ym))
        self.canvas.create_line(xm, ym, x3, y3, fill='blue', tags='axis')

	self.textId = self.canvas.create_text(0, size, anchor='sw', text="XY")

	# shadow line in X=0 plane
	self.shadowPX = self.canvas.create_polygon(0,0,0,0,0,0, fill='red',
					       tag='moving')
	self.shadowPY = self.canvas.create_polygon(0,0,0,0,0,0, fill='green',
					       tag='moving')
	self.shadowPZ = self.canvas.create_polygon(0,0,0,0,0,0, fill='blue',
					       tag='moving')

	self.shadowX = self.canvas.create_line(0, 0, 0, 0, fill='black',
					       tag='moving')
	self.shadowY = self.canvas.create_line(0, 0, 0, 0, fill='black',
					       tag='moving')
	self.shadowZ = self.canvas.create_line(0, 0, 0, 0, fill='black',
					       tag='moving')

	self.lId1 = self.canvas.create_line(0, 0, 0, 0, fill='black', width=3,
					    arrow='last')
	self.canvas.pack(side='top')
        self.frame.pack(expand=1, fill='x')
        self.xm = self.ym = self.r
        self.drawVector()


    def setVector(self, value):
        #setVector does not call a callback!
        v = value
        # compute normalized vector
        n = math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])
        if n == 0.0: v = [0.0, 0.0, 1.0]
        else: v = [v[0]/n, v[1]/n, v[2]/n]
        self.vector = v
        self.setEntries()
        self.drawVector()


 #####################################################################
 # the 'configure' methods:
 #####################################################################

    def configure(self, **kw):
        for key,value in kw.items():
            # the 'set parameter' callbacks
            if key=='continuous': self.setContinuous(value)
            elif key=='mode': self.setMode(value)
            elif key=='precision': self.setPrecision(value)

            # the 'lock entries' callbacks
            elif key=='lockContinuous': self.lockContinuousCB(value)
            elif key=='lockMode': self.lockModeCB(value)
            elif key=='lockPrecision': self.lockPrecisionCB(value)


    def setContinuous(self, cont):
        """ cont can be None, 0 or 1 """
        if cont != 1:
            cont = None
        self.continuous = cont
        if hasattr(self.opPanel, 'optionsForm'):
            w=self.opPanel.idf.entryByName['togCont']['widget']
            if cont:
                w.setvalue('on')
            else:
                w.setvalue('off')


    def setMode(self, mode):
        if mode!='XY' and mode!='X' and mode!='Y' and mode!='Z': mode = 'XY'
	self.canvas.itemconfigure( self.textId, text=mode)
	self.mode = mode

        if hasattr(self.opPanel, 'optionsForm'):
            w=self.opPanel.idf.entryByName['togAxes']['widget']
            w.setvalue(mode)
          

    def setPrecision(self, val):
        val = int(val)
        if val > 10: val = 10
        if val < 1:  val = 1
        
        self.thumbx.configure(precision=val)
        self.thumby.configure(precision=val)
        self.thumbz.configure(precision=val)

        self.entryXTk.set(self.thumbx.labelFormat%self.vector[0])
        self.entryYTk.set(self.thumby.labelFormat%self.vector[1])
        self.entryZTk.set(self.thumbz.labelFormat%self.vector[2])
        
        if hasattr(self.opPanel, 'optionsForm'):
            w = self.opPanel.idf.entryByName['selPrec']['widget']
            w.setvalue(val)

        if self.opPanel:
            self.opPanel.updateDisplay()


 #####################################################################
 # the 'lock' methods:
 #####################################################################
 

    def lockContinuousCB(self, mode):
        if mode != 0: mode = 1
        self.lockContinuous = mode
        if hasattr(self.opPanel, 'optionsForm'):
            self.opPanel.lockUnlockDisplay()


    def lockPrecisionCB(self, mode):
        if mode != 0: mode = 1
        self.lockPrecision = mode
        if hasattr(self.opPanel, 'optionsForm'):
            self.opPanel.lockUnlockDisplay()


    def lockModeCB(self, mode):
        if mode != 0: mode = 1
        self.lockMode = mode
        if hasattr(self.opPanel, 'optionsForm'):
            self.opPanel.lockUnlockDisplay()

 
if __name__ == '__main__':
    test = vectorGUI(size = 200)
    def foo(val):
        print val
    test.callbacks.AddCallback(foo)
