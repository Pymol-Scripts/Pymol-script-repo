## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

##############################################################################
#
#
#   Authors: Sowjanya Karnati,Michel F Sanner
#
#
###############################################################################
#
#
#
#
#$Id: graphtool.py,v 1.47 2007/12/04 21:28:04 vareille Exp $

#Graph Tool is a widget  with movable graph curve
#app=GraphApp(root)
#app.caluculate_ramp() returns current ramp
#from ViewerFramework.VFCommand import Command, CommandGUI
from Tkinter import *
import Tkinter
import tkFileDialog
import types,os
from mglutil.util.callback import CallBackFunction
from mglutil.util.callback import CallbackManager
from mglutil.util.misc import ensureFontCase
import numpy.oldnumeric as Numeric
from mglutil.gui.BasicWidgets.Tk.thumbwheel import ThumbWheel
import Pmw
from Pmw import *
from mglutil.util.misc import deepCopySeq

class GraphApp:
    
      
    # Initialization
    def __init__(self,master=None,callback=None,continuous=1):
               
        self.master=master
        self.callback = None
        self.callbacks = CallbackManager()
        self.canvas=canvas = Canvas(self.master,width=345,height=320,bg='white') 
        self.toolbar = Frame(master)             # Create Toolbar
        self.toolbar.pack(side='top', expand=1, fill='both') 
        self.menuFrame1 = Tkinter.Frame(self.toolbar, relief='raised', borderwidth=3)
        self.menuFrame1.pack(side='top', expand=1, fill='x')
        self.filebutton = Tkinter.Menubutton(self.menuFrame1, text='File')
        self.filebutton.pack(side='left')
        self.filemenu = Tkinter.Menu(self.filebutton, {})
        self.filemenu.add_command(label='Read', command=self.read_cb)
        self.filemenu.add_command(label='Write', command=self.write_cb)
        self.filebutton['menu'] = self.filemenu
        self.editbutton = Tkinter.Menubutton(self.menuFrame1, text='Edit')
        self.editbutton.pack(side='left', anchor='w')
        self.editmenu = Tkinter.Menu(self.editbutton, {})
        self.editmenu.add_command(label='Reset to first in history', command=self.resetAll_cb)
        self.editmenu.add_command(label='Step back in history loop', command=self.stepBack_cb)
        self.editmenu.add_command(label='Default Curve', command=self.defaultcurve_cb)
        self.editmenu.add_command(label='Invert Curve',command=self.invertGraph)
        self.histvar=IntVar()
        self.histvar.set(1)
        self.editmenu.add_checkbutton(label='Histogram',var=self.histvar,command=self.drawHistogram)
        self.editbutton['menu'] = self.editmenu
        self.optionType = IntVar()
        self.updatebutton = Tkinter.Menubutton(self.menuFrame1, text='Update')
        self.updatebutton.pack(side='left', anchor='w')
        self.updatemenu = Tkinter.Menu(self.updatebutton,{} )
        for v,s in {0:'Continuous',1:'MouseButtonUp',2:'Update'}.items():
            self.updatemenu.add_radiobutton(label=s,
                                var=self.optionType,
                                value = v,command=self.calloption)
        if continuous==1:
            self.optionType.set(0)
        self.updatebutton['menu'] = self.updatemenu
        #Curve Type
        self.CurveType = IntVar()
        self.CurveType.set(0) 
        self.Smooth=1
        self.curvebutton = Tkinter.Menubutton(self.menuFrame1, text='Curve')
        self.curvebutton.pack(side='left', anchor='w')
        self.curvemenu = Tkinter.Menu(self.curvebutton,{} )
        for v,s in {0:'Smooth',1:'Freehand'}.items():
            self.curvemenu.add_radiobutton(label=s,
                                var=self.CurveType,
                                value = v,command=self.curveoption)
        
            
        self.curvebutton['menu'] = self.curvemenu
        f1 = Tkinter.Frame(self.master)
        f1.pack(side='bottom', fill='both', expand=1)
        self.d1scalewheellab=Label(f1,text="Sensitivity")
        self.d1scalewheellab.pack(side="left")
        self.d1scalewheel=ThumbWheel(width=100, height=26,wheelPad=4,master=f1,labcfg={'fg':'black', 'side':'left', 'text':'Test:'},wheelLabcfg1={'font':(ensureFontCase('times'),14,'bold')},wheelLabcfg2={'font':(ensureFontCase('times'),14,'bold')},canvascfg={'bg':'blue'},min = 0.0,max = 1.0,precision =4,showlabel =0,value =0.013,continuous =0,oneTurn =0.01,size = 200)
        self.d1scalewheel.pack(side="left")
        #tooltip
        self.balloon = Pmw.Balloon(f1)
        self.balloon.bind(self.d1scalewheel,"cutoff value for differences in Z xoordinates,small values generate more contours")
        self.Updatebutton=Button(f1,text='  Update  ',command=self.Update)
        self.Updatebutton.pack(side=LEFT)
        
        self.Quitbutton=Button(f1,text=' Dismiss ',command=self.dismiss_cb)
        self.Quitbutton.pack(side=RIGHT)
        self.canvas.bind("<Button-1>", self.OnCanvasClicked)
        self.canvas.bind("<B1-Motion>", self.OnCanvasMouseDrag)
        self.canvas.bind("<ButtonRelease-1>", self.OnCanvasMouseUp)
        self.canvas.config(closeenough=2.0)
        self.canvas.pack(side=BOTTOM, fill=BOTH,expand=1)
        self.startpoint=(px,py)=(50,275)
        self.endpoint=(px1,py1)=(305,20)
        self.newpoints=[(px,py),(px1,py1)]
        self.canvas.create_rectangle([(px-1,py),(px1+1,py1)],fill='white',outline="black",width=1)
        self.canvas.create_text(46,281,text=0,anchor=N)
                  
        #Drawing Graph Sheet
        for i in range(1,6):
            x=50+i*50
            canvas.create_line(x,280,x,275,width=1)
            canvas.create_text(x,281,text='%d' %(50*i),anchor=N)
        for i in range(1,5):
            x=50+i*50
            canvas.create_line(x,275,x,20,width=1,fill="gray80")
        for i in range(1,6):
            y=275-i*50
            canvas.create_line(45,y,50,y,width=1)
            canvas.create_text(44,y,text='%d' %(50*i),anchor=E)
        for i in range(1,5):
            y=275-i*50
            canvas.create_line(50,y,305,y,width=1,fill="gray80")
        (x,y)=self.newpoints[0]
        (x1,y1)=self.newpoints[-1]
        self.curline=canvas.create_line(self.newpoints,fill='black',width=1)
        
        #GRAY SCALE
        grays=[]
        for i in range(0,100,1):
            grays.append("gray"+"%d" %i)
        #grays.reverse()
        #bottom one
        x1=48
        x2=51
        self.canvas.create_rectangle([(50,315),(307,300)],fill='white',outline="black",width=0.5)
        for a in grays:
            if x1>306:
                
               x1=x2=306
            self.canvas.create_rectangle([(x1+2.5,314),(x2+2.5,301)],fill=a,outline=a,width=1)
            x1=x1+2.5
            x2=x2+2.5
        #left one
        y1=274
        y2=271
        self.canvas.create_rectangle([(20,275),(5,20)],fill='black',outline="black",width=0.5)
        for a in grays:
            if y1>275:
                y1=y2=275
            self.canvas.create_rectangle([(19,y1-2.5),(6,y2-2.5)],fill=a,outline=a,width=1)    
            y1=y1-2.5
            y2=y2-2.5
              
        self.oldpoints=[]
        self.canvas.configure(cursor='cross')
        self.curovals=[]
        self.default_points=[(50,275),(88, 238), (101, 150), (154, 78), (75, 271),(305,20)]
        # now set the constructor options correctly using the configure method
        apply( self.configure, (),{'callback':callback,'continuous':continuous})
        self.continuous=continuous
        self.mousebuttonup=0
        self.update=0
        self.range_points=[]
        self.history=[]
        self.bars=[]
        self.default_ramp=[]
        self.histvalues=[]
        
    def calloption(self):
        tag=self.optionType.get()
        self.continuous=0
        self.mousebuttonup=0
        self.update=0
        if tag==0:
            self.continuous=1
        elif tag==1:
            self.mousebuttonup=1
        elif tag==2:    
            self.update=1
        
    def curveoption(self):
        tag=self.CurveType.get()
        self.Smooth=0
        self.Freehand=0
        if tag==0:
            self.Smooth=1
            self.canvas.delete(self.curline)
            self.curline=self.canvas.create_line(self.getControlPoints(),smooth=1)
        elif tag==1:
            self.Freehand=1
            self.canvas.delete(self.curline)
            self.curline=self.canvas.create_line(self.getControlPoints())
            
    def OnCanvasClicked(self,event):
        """Appends last drag point to controlpoint list if not appended by mouseleave func.
        when clicked on any controlpoint removes control point and draws line with remaining 
        control points.""" 
        self.CLICK_NODRAG=1 
        if self.history!=[]:
            if self.history[-1][1]!=self.d1scalewheel.get():
                self.history[-1]=(self.history[-1][0],self.d1scalewheel.get())
        if hasattr(self,"curx"): 
            if (self.curx,self.cury) :#not in [self.startpoint,self.endpoint]:
                (self.ox,self.oy)=(self.curx,self.cury)
                if (self.ox,self.oy) not in self.oldpoints:
                    self.oldpoints.append((self.ox,self.oy))
                    if hasattr(self,"curoval"):
                        self.curovals.append(self.curoval)
        self.OrgX=event.x
        self.OrgY=event.y
        CtlPoints=[]
        xcoords=[]
        ycoords=[]
        #Limiting points not to cross
        self.limit_xcoord=[]
        for i in range(0,10):
            xcoords.append(self.OrgX-i)
            ycoords.append(self.OrgY-i)
            xcoords.append(self.OrgX+i)
            ycoords.append(self.OrgY+i)
            
        if xcoords!=[] and ycoords!=[]:
            for x in xcoords:
                for y in ycoords:
                    CtlPoints.append((x,y))
            self.range_points=self.oldpoints
            self.range_points.sort()
            for c in CtlPoints:
             if  c in self.range_points:
                    index_c=self.range_points.index(c)
                    if index_c<len(self.range_points)-1:
                        self.limit_xcoord.append(self.range_points[index_c+1])
                        
                        if index_c>0:
                            self.limit_xcoord.append(self.range_points[index_c-1])
                            return
                        else:
                            self.limit_xcoord.append(self.startpoint)
                            return
                    elif index_c==len(self.range_points)-1:
                        self.limit_xcoord.append(self.range_points[index_c-1])
                        self.limit_xcoord.append(self.endpoint)
                        return 
        self.newd1ramp= self.caluculate_ramp()
        
    
    def OnCanvasMouseUp(self,event):
             
        CtlPoints=[]
        xcoords=[]
        ycoords=[]
        if hasattr(self,"curx"): 
            (self.ox,self.oy)=(self.curx,self.cury)
            if (self.ox,self.oy) not in self.oldpoints :#not in [self.startpoint,self.endpoint] :
                self.oldpoints.append((self.ox,self.oy))     
                if hasattr(self,"curoval"):
                    if self.curoval not in self.curovals:
                        self.curovals.append(self.curoval)
        
        if self.CLICK_NODRAG==1:
         #finding out points around the selected point
         for i in range(0,10):
            xcoords.append(self.OrgX-i)
            ycoords.append(self.OrgY-i)
            xcoords.append(self.OrgX+i)
            ycoords.append(self.OrgY+i)
            
         if xcoords!=[] and ycoords!=[]:
            for x in xcoords:
                for y in ycoords:
                    CtlPoints.append((x,y))
            
            for c in CtlPoints:
             if  c in self.oldpoints:
                ind=self.oldpoints.index(c)
                op=self.oldpoints[ind]
                if ind>0:
                    prev_oldpoint=self.oldpoints[ind-1]
                else:
                    prev_oldpoint=self.endpoint
                del self.oldpoints[ind]
                
                for co in self.curovals:
                    ov_point1=self.canvas.coords(co)
                    if len(ov_point1)!=0:
                        ov_point=(int(ov_point1[0]+2),int(ov_point1[1]+2))
                                    
                        if ov_point==c  and ov_point not in [self.startpoint,self.endpoint]:
                            self.canvas.delete(co)
                            self.curovals.remove(co) 
                            if hasattr(self,"curx"): 
                             if ov_point==(self.curx,self.cury):
                                (self.curx,self.cury)=prev_oldpoint
                                                            
                            self.draw()
        if  self.mousebuttonup:
            self.newd1ramp=self.caluculate_ramp()
            self.callbacks.CallCallbacks(self.newd1ramp)
        
        self.history.append((deepCopySeq(self.oldpoints),self.d1scalewheel.get()))
        
    
        
    
    def OnCanvasMouseDrag(self,event):
        self.CLICK_NODRAG=0
        CtlPoints=[]
        xcoords=[]
        ycoords=[]
        #making active clickrange to be ten points around clicked point  
        for i in range(0,10):
            xcoords.append(self.OrgX-i)
            ycoords.append(self.OrgY-i)
            xcoords.append(self.OrgX+i)
            ycoords.append(self.OrgY+i)
            
        if xcoords!=[] and ycoords!=[]:
            for x in xcoords:
                for y in ycoords:
                    CtlPoints.append((x,y))
            
            for c in CtlPoints:
             if  c in self.oldpoints:
                ind=self.oldpoints.index(c)
                op=self.oldpoints[ind]
                del self.oldpoints[ind]
                
                for co in self.curovals:
                    ov_point1=self.canvas.coords(co)
                    if len(ov_point1)!=0:
                        ov_point=(int(round(ov_point1[0],3))+2,int(round(ov_point1[1],3))+2)
                                    
                        if ov_point==c :
                            self.canvas.delete(co)
                            self.curovals.remove(co)
                            
        self.curx=dx=event.x
        self.cury=dy=event.y
        self.draw() 
        
    def draw(self):
        """Draws line,ovals with current controlpoints. """
        new1points=[]
        curve_points=[]
        self.smoothened_points=[]
        if self.CLICK_NODRAG==0:
            dx=self.curx
            dy=self.cury
        else:
            (dx,dy)=(self.curx,self.cury)=(self.endpoint)
        ###Limiting  xcoords of the current point not to cross adjacent points
        if hasattr(self,"limit_xcoord"):
         if self.limit_xcoord!=[]:
            self.limit_xcoord.sort()
            if  (self.curx,self.cury) not in [self.startpoint,self.endpoint]:
                if (self.curx,self.cury)< self.limit_xcoord[0]:
                    if self.curx<=self.limit_xcoord[0][0] and self.cury<self.limit_xcoord[0][1]:
                        dx=self.curx=self.limit_xcoord[0][0]+1
                    if  self.curx<=self.limit_xcoord[0][0] and self.cury>self.limit_xcoord[0][1]:
                        dx=self.curx=self.limit_xcoord[0][0]
                if  (self.curx,self.cury)> self.limit_xcoord[1]:
                    if self.curx>=self.limit_xcoord[1][0] and self.cury>self.limit_xcoord[1][1]:
                        dx=self.curx=self.limit_xcoord[1][0]-1
                    if self.curx>=self.limit_xcoord[1][0] and self.cury<self.limit_xcoord[1][1]:
                        dx=self.curx=self.limit_xcoord[1][0]
        #Limit graph with in the axis
        if self.curx not in range(50,305):
                if self.curx<50:
                    self.curx=dx=50
                else:
                    self.curx=dx=305
                    
        if self.cury not in range(20,275):
                if self.cury<20:
                    self.cury=dy=20
                else:
                    self.cury=dy=275
        #adding start,end points  
        new1points.append(self.startpoint)
        new1points.append(self.endpoint)
        #adding current point to list
        if (dx,dy) not in new1points and (dx,dy) not in [self.startpoint,self.endpoint]:
            new1points.append((dx,dy))
        #adding oldpoints to list
        if hasattr(self,"ox"):
          for op in self.oldpoints:
                if op not in new1points:
                    new1points.append(op)
        new1points.sort()
        #removing oval point that is on drag
        if hasattr(self,"curoval"):
            if self.curoval not in self.curovals:
               self.canvas.delete(self.curoval)
        self.canvas.delete(self.curline)
        #if points that start with 50 or 51  or 305,304 other than start ,end
        #points exists remove start or end points
        #remove ovals
        #finding oval for start point and endpoint
        for i in new1points:
            if i[0]==51 or i[0]==50:
                if i!=self.startpoint:
                    if self.startpoint in new1points:
                        new1points.remove(self.startpoint)
                        ###removing start point oval 
                        x = 50
                        y = 275
                        st_oval_1= self.canvas.find_enclosed(x-3,y-3,x+3,y+3)
                         
                        if st_oval_1:
                            for so in st_oval_1:
                                if so!=[]:
                                    st_oval=so
                                    st_oval_coords=self.canvas.coords(st_oval)
                                    if (int(st_oval_coords[0]+2),int(st_oval_coords[1]+2))==self.startpoint: 
                                        self.canvas.delete(st_oval)
                                        if st_oval in self.curovals:
                                            self.curovals.remove(st_oval)
                        
                        
        for i in new1points:   
            if i[0]==304 or i[0]==305: 
                if i!=self.endpoint :
                    if self.endpoint in new1points:
                        new1points.remove(self.endpoint)
                        ###removing end point oval
                        x = 305
                        y = 20
                        end_oval_1= self.canvas.find_enclosed(x-3,y-3,x+3,y+3)
                        if end_oval_1:
                            for eo in end_oval_1:
                                if eo!=[]:
                                    end_oval=eo
                                    end_oval_coords=self.canvas.coords(end_oval)
                                    if (int(end_oval_coords[0]+2),int(end_oval_coords[1]+2))==self.endpoint: 
                                        self.canvas.delete(end_oval)
                                        if end_oval in self.curovals:
                                            self.curovals.remove(end_oval)
        new1points.sort()
          
        for (x,y) in new1points:
            curve_points.append(x)
            curve_points.append(y)
        self.smoothened_points= self.smooth(curve_points)
        #drawing line
        if len(self.smoothened_points)>2: 
            if self.Smooth:
                self.curline=self.canvas.create_line(self.smoothened_points)
            else:
                self.curline=self.canvas.create_line(curve_points)
        else:
            
            if curve_points[0]==50 or 51:
                 if self.Smooth:
                    self.curline=self.canvas.create_line(curve_points,smooth=1)
                 else:
                    self.curline=self.canvas.create_line(curve_points)   
            else:    
                self.curline=self.canvas.create_line(self.startpoint,self.endpoint)
        
        ##Adding oval when start or end point in new1points
        coval_coords=[]
        for i in self.curovals:
            coval_coords.append(self.canvas.coords(i))
        if self.endpoint in new1points:
            co=self.canvas.create_oval(self.endpoint[0]-2,self.endpoint[-1]-2,self.endpoint[0]+2,self.endpoint[-1]+2,width=1,outline='black',fill='black')
            endco_coords =self.canvas.coords(co)
            if endco_coords not in  coval_coords:
                self.curovals.append(co)
        if self.startpoint in new1points:
            co=self.canvas.create_oval(self.startpoint[0]-2,self.startpoint[-1]-2,self.startpoint[0]+2,self.startpoint[-1]+2,width=1,outline='black',fill='black')          
            startco_coords=self.canvas.coords(co)
            if startco_coords not in  coval_coords:
                self.curovals.append(co)
        #drawing ovals
        if (self.curx,self.cury)!=self.endpoint:
            self.curoval=self.canvas.create_oval(self.curx-2,self.cury-2,self.curx+2,self.cury+2,width=1,outline='black',fill='black') 
            
        if (self.curx,self.cury)==self.endpoint and self.endpoint in new1points:
            self.curoval=self.canvas.create_oval(self.curx-2,self.cury-2,self.curx+2,self.cury+2,width=1,outline='black',fill='black')
        self.newd1ramp= self.caluculate_ramp()
        if self.continuous:
                self.callbacks.CallCallbacks(self.newd1ramp)
            
  
        
    ########  convert coordinates to ramp##################  
    def caluculate_ramp(self):
        """ 
        """
        dramp=[] 
        mypoints=[]
        mynewpoints=[]
        self.oldpoints.sort()
        calcpoints=[]
        #if self.continuous :
        if hasattr(self,"curx"):
                if (self.curx,self.cury) not in self.oldpoints and (self.curx,self.cury) not in [self.startpoint,self.endpoint]:
                    calcpoints.append((self.curx,self.cury))
        if len(self.oldpoints)!=0:
            for o in self.oldpoints:
                if o not in calcpoints:
                    calcpoints.append(o)
        if self.startpoint not in calcpoints:
                calcpoints.append(self.startpoint)
        if self.endpoint not in calcpoints:
                calcpoints.append(self.endpoint)
        
        calcpoints.sort()    
        length=len(calcpoints)
        for l in range(length):
            if l+1<=length-1:
                mypoints=[calcpoints[l],calcpoints[l+1]]
                if calcpoints[l] not in mynewpoints:
                    mynewpoints.append( calcpoints[l])   
                
                (x1,y1)=calcpoints[l]
                (x2,y2)=calcpoints[l+1]
                if x1>x2:
                    dcx=x1-x2
                    px=x1-1
                else:
                    dcx=x2-x1
                    px=x1+1
                if y1>y2:
                    dcy=y1-y2
                    if dcx>=1:
                        py=y1-float(dcy)/float(dcx)
                    else:
                        py=y1
                else:   
                    dcy=y2-y1
                    if dcx>=1:
                        py=y1+float(dcy)/float(dcx)
                    else:
                        py=y2
                mynewpoints.append( (px,int(round(py))))
                for dc in range(dcx-1):
                    
                    if x1>x2:
                        px=px-1
                    else:
                        px=px+1
                    if y1>y2:
                        if dcx>=1:
                            py=py-float(dcy)/float(dcx)
                        else:
                            py=y1
                    else:
                        if dcx>=1:
                            py=py+float(dcy)/float(dcx)
                        else:
                            py=y2
                    mynewpoints.append( (px,int(round(py))))
        ramp=[]
        for r in mynewpoints:
            #scale
            ra=float(275-r[1])
            if ra>=256:
                ra=255.0
            ramp.append(ra)
            dramp=Numeric.array(ramp,'f')
        if len(dramp)!=0:
            return   dramp     
        else:
            dramp=Numeric.arange(0,256,1,'f')
            return   dramp
        
    def get(self):
        if hasattr(self,"newd1ramp"):
           return self.newd1ramp
        else:    
            return self.caluculate_ramp()
        
    def configure(self, **kw):
        if 'type' in kw.keys(): # make sure type is set first
            self.setType(kw['type'])
            del kw['type']
            
        for key,value in kw.items():
            if key=='callback': 
                self.setCallbacks(value)
            
   
    def setCallbacks(self, cb):
        """Set widget callback. Must be callable function. Callback is called
every time the widget value is set/modified"""

        assert cb is None or callable(cb) or type(cb) is types.ListType,\
               "Illegal callback: must be either None or callable. Got %s"%cb
        if cb is None: return
        elif type(cb) is types.ListType:
            for func in cb:
                assert callable(func), "Illegal callback must be callable. Got %s"%func
                self.callbacks.AddCallback(func)
        else:
            self.callbacks.AddCallback(cb)
        self.callback = cb
        

    def invertGraph(self):
        """This function is for inverting graph by reverse computing controlpoints"""
        if self.history!=[]:
            if self.history[-1][1]!=self.d1scalewheel.get():
                self.history[-1]=(self.history[-1][0],self.d1scalewheel.get())
        invert_points=[]
        #self.oldpoints=[]
        points=self.getControlPoints()
        if len(points)<2:
            points=[self.startpoint,self.endpoint]
        for p in points:
            if p[1] in range(20,276):
                y=275 -(p[1]-20)
                invert_points.append((p[0],y))
        self.reset()
        ###################################################
        #Some times start and end points are not deleted 
        #So for deleting them canvas.find_enclosed points at
        #startpoint and endpoint are caluculated(returns alist of
        #canvas objects present there) and if the coords of
        #any canvas objects matches with start or end point that gets deleted
        #####################################################             
        x = 50
        y = 275
        st_oval_1= self.canvas.find_enclosed(x-3,y-3,x+3,y+3)
        if st_oval_1:
            for so in st_oval_1:
                if so!=[]:
                    st_oval=so
                    st_oval_coords=self.canvas.coords(st_oval)
                    if (int(st_oval_coords[0]+2),int(st_oval_coords[1]+2))==self.startpoint: 
                        self.canvas.delete(st_oval)
                        if st_oval in self.curovals:
                            self.curovals.remove(st_oval)
        x = 305
        y = 20
        end_oval_1= self.canvas.find_enclosed(x-3,y-3,x+3,y+3)
        if end_oval_1:
            for eo in end_oval_1:
                if eo!=[]:
                    end_oval=eo
                    end_oval_coords=self.canvas.coords(end_oval)
                    if (int(end_oval_coords[0]+2),int(end_oval_coords[1]+2))==self.endpoint: 
                        self.canvas.delete(end_oval)
                        if end_oval in self.curovals:
                            self.curovals.remove(end_oval)
        self.canvas.delete(self.curline)
        if self.Smooth:
              self.curline=self.canvas.create_line(invert_points,smooth=1)  
        else:
                self.curline=self.canvas.create_line(invert_points)
        self.oldpoints=invert_points
        for p in invert_points:
            self.curoval=self.canvas.create_oval(p[0]-2,p[1]-2,p[0]+2,p[1]+2,width=1,outline='black',fill='black')
            self.curovals.append(self.curoval)
        (self.curx,self.cury) =invert_points[-2]    
        if self.continuous or self.mousebuttonup:
            self.newd1ramp=self.caluculate_ramp()
            self.callbacks.CallCallbacks([self.newd1ramp])
        self.history.append((deepCopySeq(self.oldpoints),self.d1scalewheel.get()))
         

    def defaultcurve_cb(self):
        """draws curve with default points"""
        if self.history!=[]:
            if self.history[-1][1]!=self.d1scalewheel.get():
                self.history[-1]=(self.history[-1][0],self.d1scalewheel.get())
        points=[]
        self.default_points=[]
        self.oldpoints=[]
        self.d1scalewheel.set(0.013)
        self.default_points=[(50,275),(88, 238), (101, 150), (154, 78), (75, 271),(305,20)]
        self.reset()
        self.canvas.delete(self.curline)
        self.default_points.sort()
        if self.Smooth:
            self.curline=self.canvas.create_line(self.default_points,smooth=1) 
        else:
                self.curline=self.canvas.create_line(self.default_points)
        self.oldpoints=self.default_points
        for p in self.default_points:
            self.curoval=self.canvas.create_oval(p[0]-2,p[1]-2,p[0]+2,p[1]+2,width=1,outline='black',fill='black')
            self.curovals.append(self.curoval)
        (self.curx,self.cury) =self.default_points[-2]    
        if self.continuous or self.mousebuttonup:
            self.newd1ramp=self.caluculate_ramp()
            self.callbacks.CallCallbacks(self.newd1ramp)
        self.history.append((deepCopySeq(self.oldpoints),self.d1scalewheel.get()))
        self.default_ramp= self.newd1ramp
        
    def read_cb(self):
        fileTypes = [("Graph",'*_Graph.py'), ("any file",'*.*')]
        fileBrowserTitle = "Read Graph"
        fileName  = self.fileOpenAsk(types=fileTypes,
                                     title=fileBrowserTitle)
        if not fileName:
            return
        self.read(fileName)

    def read( self,fileName):
        if self.history!=[]:
            if self.history[-1][1]!=self.d1scalewheel.get():
                self.history[-1]=(self.history[-1][0],self.d1scalewheel.get())
        fptr=open(fileName,"r")
        data=fptr.readlines()
        cpoints=data[0][:-1]
        sensitivity=data[1]
        self.d1scalewheel.set(eval(sensitivity))
        if len(cpoints)==0:
            return
        else:
            points=cpoints
        self.oldpoints=[]        
        self.reset()
        if hasattr(self,"curline"):
            self.canvas.delete(self.curline)
            for c in self.curovals:
                self.canvas.delete(c)
                self.curovals.remove(c)
            if hasattr(self,"curoval"):
                self.canvas.delete(self.curoval)
            self.curovals=[]
            if self.Smooth:
                self.curline=self.canvas.create_line(eval(points),smooth=1)
            else:
                self.curline=self.canvas.create_line(eval(points))
            self.readpoints=self.oldpoints=eval(points)[1:-1]
            for p in eval(points)[1:-1]:
                self.curoval=self.canvas.create_oval(p[0]-2,p[1]-2,p[0]+2,p[1]+2,width=1,outline='black',fill='black')
                self.curovals.append(self.curoval)
            (self.curx,self.cury) =eval(points)[-2]
            self.history.append((deepCopySeq(self.oldpoints),self.d1scalewheel.get()))

    def fileOpenAsk(self, idir=None, ifile=None, types=None,
                    title='Open'):
        if types==None: types = [ ('All files', '*') ]
        file = tkFileDialog.askopenfilename( filetypes=types,
                                             initialdir=idir,
                                             initialfile=ifile,
                                             title=title)
        if file=='': file = None
        return file 
    
    def write_cb(self):
        fileTypes = [("Graph",'*_Graph.py'), ("any file",'*.*')]
        fileBrowserTitle = "Write Graph"
        fileName  = self.fileSaveAsk(types=fileTypes,
                                     title=fileBrowserTitle)
        if not fileName:
            return
        self.write(fileName)
    
    def write(self,fileName):
        fptr=open(fileName,"w")
        points= self.getControlPoints()        
        points.sort()
        fptr.write(str(points))
        fptr.write("\n")
        fptr.write(str(self.d1scalewheel.get()))
        fptr.close()
        
    def fileSaveAsk(self, idir=None, ifile=None, types = None,
                title='Save'):
        if types==None: types = [ ('All files', '*') ]
        file = tkFileDialog.asksaveasfilename( filetypes=types,
                                               initialdir=idir,
                                               initialfile=ifile,
                                               title=title)
        if file=='': file = None
        return file
        
    def reset(self):
        """This function deletes current line removes current ovals"""  
        self.canvas.delete(self.curline)
        self.oldpoints=[]
        for c in self.curovals:
            self.canvas.delete(c)
        if hasattr(self,"curoval"):
            self.canvas.delete(self.curoval)
        self.curovals=[]
        if hasattr(self,"curoval"):
            delattr(self,"curoval")
        if hasattr(self,"curx"):
            delattr(self,"curx")
        

    def resetAll_cb(self):
        """Resetting curve as slant line 0 to 255"""
        self.reset()    
        self.curline=self.canvas.create_line([self.startpoint,self.endpoint],width=1,fill='black')
        for p in [self.startpoint,self.endpoint]:
            self.curoval=self.canvas.create_oval(p[0]-2,p[1]-2,p[0]+2,p[1]+2,width=1,outline='black',fill='black') 
            self.curovals.append(self.curoval)
        self.oldpoints=[self.startpoint,self.endpoint]
        (self.curx,self.cury)=self.endpoint
        self.d1scalewheel.set(0.013)
        
        if self.continuous or self.mousebuttonup:
            self.newd1ramp=Numeric.arange(0,256,1,'f')
            self.callbacks.CallCallbacks(self.newd1ramp)
        #self.histvar.set(0)
        self.history=[]
        
        
    def stepBack_cb(self):
        """when stepBack button clicked previous step is displayed.History of
        all the steps done is remembered and when stepback clicked from history
        list previous step is shown and that step is removed from history list """
        
        if self.history!=[]:
            if len(self.history)==1:
                self.resetAll_cb()
            
            else:
                del self.history[-1]
                pns = self.history[-1][0]
                #deleting 
                self.oldpoints=pns
                self.canvas.delete(self.curline)
                for c in self.curovals:
                    self.canvas.delete(c)
                if hasattr(self,"curoval"):
                    self.canvas.delete(self.curoval)
                self.curovals=[]
                ###################################################
                #Some times start and end points are not deleted 
                #So for deleting them canvas.find_enclosed points at
                #startpoint and endpoint are caluculated(returns alist of
                #canvas objects present there) and if the coords of
                #any canvas objects matches with start or end point that gets deleted
                #####################################################             
                x = 50
                y = 275
                st_oval_1= self.canvas.find_enclosed(x-3,y-3,x+3,y+3)
                if st_oval_1:
                            for so in st_oval_1:
                                if so!=[]:
                                    st_oval=so
                                    st_oval_coords=self.canvas.coords(st_oval)
                                    if (int(st_oval_coords[0]+2),int(st_oval_coords[1]+2))==self.startpoint: 
                                        self.canvas.delete(st_oval)
                                        if st_oval in self.curovals:
                                            self.curovals.remove(st_oval)
                x = 305
                y = 20
                end_oval_1= self.canvas.find_enclosed(x-3,y-3,x+3,y+3)
                if end_oval_1:
                            for eo in end_oval_1:
                                if eo!=[]:
                                    end_oval=eo
                                    end_oval_coords=self.canvas.coords(end_oval)
                                    if (int(end_oval_coords[0]+2),int(end_oval_coords[1]+2))==self.endpoint: 
                                        self.canvas.delete(end_oval)
                                        if end_oval in self.curovals:
                                            self.curovals.remove(end_oval)
                
                
                
                pns.sort()
                #if no start or end points 
                if pns[0][0]>51 :
                    pns.insert(0,self.startpoint)
                l=len(pns) 
                if pns[-1][0]<304:
                    pns.insert(l,self.endpoint)
                
                #if start or endpoints and points with (50or 51) or (305or305)
                if self.startpoint in pns:
                    for p in pns:
                        if p!=self.startpoint:
                            if p[0]== 50 or p[0]==51:
                               pns.remove(self.startpoint)
                if self.endpoint in pns:
                    for p in pns:
                        if p!=self.endpoint:
                            if p[0]==305 or p[0]==304:
                                pns.remove(self.endpoint)
                      
                print pns
                if self.Smooth:
                    self.curline=self.canvas.create_line(pns,width=1,fill='black',smooth=1)
                else:
                    self.curline=self.canvas.create_line(pns,width=1,fill='black')
                
                for p in pns:
                    self.curoval=self.canvas.create_oval(p[0]-2,p[1]-2,p[0]+2,p[1]+2,width=1,outline='black',fill='black')
                    self.curovals.append(self.curoval)
                self.d1scalewheel.set(self.history[-1][1])
                
                 
                if self.continuous or self.mousebuttonup:
                    self.newd1ramp=Numeric.arange(0,256,1,'f')
                    self.callbacks.CallCallbacks(self.newd1ramp)
                (self.curx,self.cury)=self.endpoint    
    def getControlPoints(self):
        """fuction to get current control points of the curve"""
        if not self.oldpoints==[self.startpoint,self.endpoint]:
          for i in range(len(self.oldpoints)):
            if self.startpoint  in self.oldpoints:
                self.oldpoints.remove(self.startpoint)
            if self.endpoint in self.oldpoints:
                self.oldpoints.remove(self.endpoint)
        self.controlpoints=[]
        
        if hasattr(self,"curoval"):
           c=self.canvas.coords(self.curoval)
           if len(c)!=0:
                if (int(c[0]+2),int(c[1]+2)) not in self.oldpoints and (int(c[0]+2),int(c[1]+2)) not in [self.startpoint,self.endpoint]:
                    self.controlpoints.append((int(c[0]+2),int(c[1]+2))) 
        for op in self.oldpoints:
                self.controlpoints.append(op)
        self.controlpoints.sort()
        if len(self.controlpoints)>0:
            if self.controlpoints[0][0]==50 or self.controlpoints[0][0]==51 :
                pass
            else:
                self.controlpoints.append(self.startpoint)
            self.controlpoints.sort()
            if self.controlpoints[-1][0]==305 or self.controlpoints[-1][0]==304:
                pass
            else:
                self.controlpoints.append(self.endpoint)
            self.controlpoints.sort() 
        return  self.controlpoints   
                
    def setControlPoints(self,points):
        """function to set curve control points"""
        assert isinstance(points, types.ListType),"Illegal type for points"
        
        for (x,y) in points:
            assert  x in range(50,306),"coordinates are out of range,x should be in [50,305]"
            assert  y in range(20,276),"coordinates are out of range,y should be in [20,275]"
                
        self.oldpoints=[]
        self.controlpoints=[]
        self.reset()
        self.oldpoints=self.controlpoints=points
        self.controlpoints.sort()
        if self.controlpoints[0]!=self.startpoint:
            self.controlpoints.append(self.startpoint)
        if  self.controlpoints[-1]!=self.endpoint:
            self.controlpoints.append(self.endpoint)
        self.canvas.delete(self.curline)   
        self.controlpoints.sort()
        if self.Smooth:
            self.curline=self.canvas.create_line( self.controlpoints,smooth=1)
        else:
            self.curline=self.canvas.create_line( self.controlpoints)
        for p in self.controlpoints[1:-1]:
            self.curoval=self.canvas.create_oval(p[0]-2,p[1]-2,p[0]+2,p[1]+2,width=1,outline='black',fill='black')
            self.curovals.append(self.curoval)
        (self.curx,self.cury)=   self.controlpoints[-2] 
        self.history.append((deepCopySeq(self.oldpoints),self.d1scalewheel.get()))
        if self.continuous or self.mousebuttonup:
            self.newd1ramp=self.caluculate_ramp()
            self.callbacks.CallCallbacks(self.newd1ramp)
             
    def setSensitivity(self,val):
        self.d1scalewheel.set(val)
        
    def Update(self):
        if self.update==1:
            dramp=self.caluculate_ramp()
            self.newd1ramp=dramp
            self.callbacks.CallCallbacks(dramp)
                
    def dismiss_cb(self):
       
       try:
            if self.master.winfo_ismapped():
                self.master.withdraw()
       except:
            if self.master.master.winfo_ismapped():
                self.master.master.withdraw()

        
    #draw Histogram

    def removeHistogram(self):
        """removes Histograms"""
        for b in self.bars:
            self.canvas.delete(b)
        self.bars=[]
    
    def drawHistogram(self):
       """This function draws histogram from list of pixel counts ,one for
            each value in the source image"""        
       self.removeHistogram()
       if self.histvar.get():
           h=self.histvalues 
           if h==[]:
            return
           list_pixels_count=h
           c=[]
           maxc=max(list_pixels_count)
           if maxc==0:
            return
           if list_pixels_count.index(maxc):
                
                list_pixels_count.remove(maxc)
                list_pixels_count.insert(255,0)
           
           for i in list_pixels_count[:256]:
                max_list=max(list_pixels_count)
                if max_list==0:
                    return
                val=i*200/max_list
                c.append(val)
           
           for i in range(0,len(c)):
                x1=50+i
                x2=50+i
                y1=275-c[i]
                y2=275
                r=self.canvas.create_line([(x1,y1),(x2,y2)],fill="gray70",width=1)
                
                self.bars.append(r)
           #displaying line and ovals ontop
           self.canvas.tkraise(self.curline)
           for i in self.curovals:
                self.canvas.tkraise(i)
           if hasattr(self,"curoval"):
                self.canvas.tkraise(self.curoval)
           self.canvas.update()
           ##Update Histograms on Graphtool
           #if self.update!=1:
           #      prev_option=self.optionType.get()
           #      self.optionType.set(2)
           #      self.update=1
           #      self.Update()
           #      self.optionType.set(prev_option)
           #      self.update=0
           #      return
           
               
    ################SMOOTHING  CODE############################
    def addcurve(self,out, xy, steps):
        
        add = out.append
        for i in range(1, steps+1):
            t = float(i) / steps; t2 = t*t; t3 = t2*t
            u = 1.0 - t; u2 = u*u; u3 = u2*u
            add(xy[0]*u3 + 3*(xy[2]*t*u2 + xy[4]*t2*u) + xy[6]*t3)
            add(xy[1]*u3 + 3*(xy[3]*t*u2 + xy[5]*t2*u) + xy[7]*t3)

    def smooth(self,xy, steps=12):
        
        if not xy:
            return xy
        closed = xy[0] == xy[-2] and xy[1] == xy[-1]
        out = []
        if closed:
            # connect end segment to first segment
            control = (
            0.500*xy[-4] + 0.500*xy[0],
            0.500*xy[-3] + 0.500*xy[1],
            0.167*xy[-4] + 0.833*xy[0],
            0.167*xy[-3] + 0.833*xy[1],
            0.833*xy[0]  + 0.167*xy[2],
            0.833*xy[1]  + 0.167*xy[3],
            0.500*xy[0]  + 0.500*xy[2],
            0.500*xy[1]  + 0.500*xy[3],
            )
            out = [control[0], control[1]]
            self.addcurve(out, control, steps)
        else:
            out = [xy[0], xy[1]]
        for i in range(0, len(xy)-4, 2):
            if i == 0 and not closed:
                control = (xy[i],xy[i+1],0.333*xy[i] + 0.667*xy[i+2],0.333*xy[i+1] + 0.667*xy[i+3],)
            else:
                control = (
                0.500*xy[i] + 0.500*xy[i+2],
                0.500*xy[i+1] + 0.500*xy[i+3],
                0.167*xy[i] + 0.833*xy[i+2],
                0.167*xy[i+1] + 0.833*xy[i+3],
                )
            if i == len(xy)-6 and not closed:
                control = control + (
                0.667*xy[i+2] + 0.333*xy[i+4],
                0.667*xy[i+3] + 0.333*xy[i+5],
                xy[i+4],
                xy[i+5],
                )
            else:
                control = control + (
                0.833*xy[i+2] + 0.167*xy[i+4],
                0.833*xy[i+3] + 0.167*xy[i+5],
                0.500*xy[i+2] + 0.500*xy[i+4],
                0.500*xy[i+3] + 0.500*xy[i+5],
                )
            if ((xy[i] == xy[i+2] and xy[i+1] == xy[i+3]) or
                (xy[i+2] == xy[i+4] and xy[i+3] == xy[i+5])):
                out.append(control[6])
                out.append(control[7])
            else:
                self.addcurve(out, control, steps)
        return out 



