## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#
#$Id: spline.py,v 1.5 2007/08/08 02:03:49 vareille Exp $

#This class is to create a default ("natural") spline, simply use sp = Spline(points,slopes).
#where points=[(x1,y1),(x2,y2)..............(xn,yn)]
#where slopes=[(None,h1),(l2,h2)...............(ln,None)]
#(l1,h1)=(lowslope,highslope)
#sp.spline_points gives list of spline points

from numpy.oldnumeric import *
import math
from types import *
from types import TupleType,ListType

ArrayType = type(asarray(1.0))

def array_map(f, ar):
        "Apply an ordinary function to all values in an array."
        flat_ar = ravel(ar)
        out = zeros(len(flat_ar), flat_ar.dtype.char)
        for i in xrange(len(flat_ar)):
	        out[i] = f(flat_ar[i])
                out.shape = ar.shape
        return out

class Spline:
    def __init__(self, points, slopes):
        if len(points)!=len(slopes):
            return "Number of points should be equal to number of slopes"
        self.list_y2_vals=[]
        self.spline_points=[]
        self.s_points=[]
        self.slopeVectors=slopes
        if len(points)>=2:
            self.points = points
        else:
            return "insufficient data for points,no. of points should be greater than or equal to two"
        #checks for type of slopes if [(None,(x,y),((x1,y1),(x2,y2)),.......] converts to  [(None,s1),(s2,s3),................]
        #where s=y/x
        if len(slopes)>=2:
            if len(slopes[0])==2:  
                if type(slopes[0][1])==TupleType or type(slopes[0][1])== ListType:
                    self.slopes=[]
                    if slopes[0][0]==None:
                        _slope_low=None
                        _slope_high=slopes[0][1][1]/slopes[0][1][0]
                        self.slopes.append((_slope_low,_slope_high))
                    for s in slopes[1:-1]:
                        _slope_low=s[0][1]/s[0][0]
                        _slope_high=s[1][1]/s[1][0]
                        self.slopes.append((_slope_low,_slope_high))
                    if slopes[-1][1]==None:
                        _slope_low=slopes[-1][0][1]/slopes[-1][0][0]
                        _slope_high=None
                        self.slopes.append((_slope_low,_slope_high))
                else:   
                    self.slopes=slopes
            else:
                return "invalid data for slopes"
        else:
            return "insufficient data for points,no. of slopes should be greater than or equal to two "
        
        
        
        #if len(points)%4!=0:
        #    print "enter correct points"
        #    return
        
        #finding vector points
        new_points = []
        len_sl = len(self.slopeVectors)
        for (s,p) in zip(self.slopeVectors,points):
            
            sll,slr = s
            
            if sll:
                vx,vy = sll
                n  =  30/sqrt(vx*vx+vy*vy)
                vx *=  n
                vy *=  n
                new_points.append((p[0]+vx,p[1]+vy))
                new_points.append(p)
            if slr:
                vx,vy = slr
                n  =  30/sqrt(vx*vx+vy*vy)
                vx *=  n
                vy *=  n
                new_points.append(p)
                new_points.append((p[0]+vx,p[1]+vy))       
                    
                    
        
        _points=[] 
        time_steps=[]
        for i in range(len(new_points)):
          if i%4==0:
            if i+3 <=len(new_points)-1  :
                _points.append([new_points[i],new_points[i+1],new_points[i+2],new_points[i+3]])
        for j in range(100):
            time_steps.append(float(j)/100)
        for p in _points:
          curve_points=[]
          for  t in   time_steps:
                res = self.beizer_calc(t,p)
                curve_points.append(res)
          self.spline_points.append(curve_points)
               
               
    def beizer_calc(self,t,points):
         _t  =1-t
         p0,p1,p2,p3 = points
         (x0,y0) =p0
         (x1,y1) =p1
         (x2,y2) =p2
         (x3,y3) =p3
         c0 = ((t*x0+_t*x1),(t*y0+_t*y1))
         c1 = ((t*x1+_t*x2),(t*y1+_t*y2))
         c2 =((t*x2+_t*x3),(t*y2+_t*y3))
         
         d0 = ((t*c0[0]+_t*c1[0]),(t*c0[1]+_t*c1[1]))
         d1 = ((t*c1[0]+_t*c2[0]),(t*c1[1]+_t*c2[1]))
         R = ((t*d0[0]+_t*d1[0]),(t*d0[1]+_t*d1[1]))
         return R
        
        
        
    

    def getControlPoints(self):
        """function to get control points"""
        return self.points

    def getSlopesVectors(self):
        """function to getslopes"""
        if  len(self.slopeVectors)>0:
            if len(self.slopeVectors[0])==2:
                if type(self.slopeVectors[0][1])==TupleType or type(self.slopeVectors[0][1])==ListType:
                    return self.slopeVectors
                else:
                    print "not given in vector form"
                    
    def setControlPoints(self,points,slopes):
        """function to setcontrolpoints"""
        self.__init__(points,slopes)
  
    def insertControlPoint(self,point,slope,position):
        """function to insert a control point and its slope as
        (lowslope,highslope) at a position """
        points = self.getControlPoints()
        slopes=self.getSlopesVectors()
        points.insert(position,point)
        if len(slope)==2:
            if type(slope[1])==TupleType or type(slopes[1])==ListType:
                slopes.insert(position,slope)
        self.setControlPoints(points,slopes)
        
    def deleteControlPoint(self,point):
        """function to delete controlpoint"""
        points=self.getControlPoints()
        points=map(lambda p: ((int(round(p[0])),int(round(p[1])))) ,points)
        slopes=self.getSlopesVectors()
        if point in points:
            point_index=points.index(point)
            points.remove(point)
            slope=slopes[point_index]
            slopes.remove(slope)
            self.setControlPoints(points,slopes)
        else:
            return "point not in controlPoints"
        


#example
from numpy.oldnumeric import arange, cos, Float
import Tkinter
class FunctionGraph1(Tkinter.Canvas):

    def __init__(self, x, y,pn):
        
        sizex = 2*len(x)
        sizey = 2*len(y)
        self.canvas = Tkinter.Canvas(width=sizex+500, height=sizey+300)
        coords = []
        bb = min(x), min(y), max(x), max(y)
        points=[]
        for p in pn:
            points.append((2*p[0],sizey-2*p[1]))
            
        for pq in points:
            self.canvas.create_oval( pq[0]-2,pq[1]-2,pq[0]+2,pq[1]+2,fill="black" )
            self.canvas.create_text(pq[0]+10,pq[1]+10,text = (pq[0],pq[1]))
        for a,b in zip(x,y):
            coords.append(2*a)
            coords.append(sizey-2*b)
        print len(coords)
        #self.canvas.create_line( *coords )
        self.canvas.create_line( coords[0:200] )
        self.canvas.create_line( coords [200:400])
        self.canvas.pack()

    def destroy(self):
        self.canvas.master.destroy()
##p=[(0,0),(100.,100.),(200.,200.),(250.,150.),(300.,200.),(350.,250.)]
###p=[(0,0),(100,100),(200,200),(250,150),(300,200),(350,250)]
###s=[[None,1],[-1,1],[-1,1],[-1,-1],[-1,1],[-1,None]]
###s=[(None,1),(-1,1),(-1,1),(-1,-1),(-1,1),(-1,None)]
##s=[(None,(1,1)),((1,-1),(1,1)),((1,-1),(1,1)),((1,-1),(1,-1)),((1,-1),(1,1)),((1,-1),None)]
#p=[(20,0),(100.,20.),(200.,200.)]
##s=[(None,(30,30)),((20,-20),(14,16)),((13,-17),None)]
#s=[(None,(1,1)),((1,1),(1,1)),((1,1),None)]
#sp=Spline(p,s)
###sp=Spline([(0,0),(.5,.7),(1,1)],[ (None, (1,1)), ((-1,-1),(1,0)), ((1,0), None) ] )
#xc=[]
#yc=[]
#res = sp.spline_points
#for r in res:
#     
#     for i in r:
#        xc.append(i[0])
#        yc.append(i[1])
#
#FunctionGraph1(xc,yc,p)
        
                
