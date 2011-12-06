## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#
#
#
#$Id: animatorSpline.py,v 1.5 2007/07/24 17:30:40 vareille Exp $


from splineGUI import splineGUI
from spline import Spline
from mglutil.util.callback import CallbackFunction
from types import IntType,FloatType

class AnimatorSpline(splineGUI):

    def __init__(self, spline, canvas, viewport, uniqueId, scale = None):
        """This class is used in AnimatorGUI  for spline interpolation"""
        splineGUI.__init__(self,spline, canvas, viewport, uniqueId, scale = None)
        
    def doit(self,startvalue,endvalue,range_begin,range_end,slopes):
        self.configure(startvalue = startvalue,endvalue = endvalue,range_begin = range_begin,range_end = range_end,slopes = slopes)
        if self.startvalue<self.range_begin or self.endvalue>self.range_end:
            print "invalid input,startvalue should be >= range_begin and endvalue<=range_end"
            return
        p1 = (self.viewport[0],self.getCoordsForVal(self.startvalue))
        p2 = (self.viewport[1],self.getCoordsForVal(self.endvalue))
        self.spline.setControlPoints([p1,p2],self.slopes )
        self._line_On =1
        self._slopes_0n =1
        self.drawSpline()
        self.tagBind()
        self.tagUnBind()
         
    def configure(self,**kw):
        startvalue = kw.get("startvalue")
        if startvalue is not None:
           if isinstance(startvalue, IntType or type(startvalue) == FloatType):
                self.startvalue = startvalue
           else:
                raise TypeError("bad type for startvalue")
        endvalue = kw.get("endvalue")
        if endvalue is not None:
            if isinstance(endvalue, IntType or type(endvalue) == FloatType):
                self.endvalue = endvalue
            else:
                raise TypeError("bad type for endvalue")
        range_begin = kw.get("range_begin")
        if  range_begin is not None:
            if isinstance(range_begin, IntType or type(range_begin) == FloatType):
                self.range_begin = range_begin
            else:
                raise TypeError("bad type for range_begin")
        range_end =   kw.get("range_end")  
        if range_end is not None:
            if isinstance(range_end, IntType or type(range_end) == FloatType):
                self.range_end = range_end
            else:
                raise TypeError("bad type for range_end")
        slopes = kw.get("slopes")
        if slopes is not None:
            self.slopes = slopes

    def tagBind(self):
        self.canvas.bind("<Motion>", self.OnCanvas)
        lid = self.canvas.find_withtag("line_"+self.uniqueId)
        for i in lid:
            cb = CallbackFunction(self.OnLine,i)
            self.canvas.tag_bind(i,"<Motion>",cb)
        
    def tagUnBind(self):
        lid = self.canvas.find_withtag("line_"+self.uniqueId)
        for i in lid:
            self.canvas.tag_unbind(i,'<Double-Button-1>')
    

    def moveControlPoint(self, cid, index, event):
        splineGUI.moveControlPoint(self, cid, index, event)
        self.tagBind()
        self.tagUnBind()
        
    def moveSlopePoint(self ,sid, lid,tid, event):
        splineGUI.moveSlopePoint(self ,sid, lid,tid, event)
        self.tagBind()
        self.tagUnBind()

    def getValue(self,fraction):
        spline_points = self.spline_points
        p = self.getCoords(fraction)
        y = p[1]
        value = (self.endvalue)*(y/self.range_end)
        return val
        
    def getCoords(self,fraction):
        spline_points = self.spline_points
        fr = fraction*10
        fn = round(fr)
        return spline_points[int(fn)]
        
    def getCoordsForVal(self,val):
        x1 = self.viewport[0]
        x2 = self.viewport[1]
        y1 = self.viewport[1]
        y2 = self.viewport[0]
        y = (y2-y1)*(val/float(self.range_end))+y1
        return y

    def getValueforCoords(self,coords):
        """returns value for coords"""
        y = coords[-1]
        rng = self.range_end-(self.range_begin)
        result = self.range_begin+rng*((self.viewport[1]-float(y))/(self.viewport[1]-self.viewport[0]))
        return result

    def getValuesForAllFourPoints(self):
        """returns values for all four coords"""
        cids = self.canvas.find_withtag("ctrlPoints_"+self.uniqueId)
        values = []
        newlist = []
        points = self.getControlPoints()
        l,r = self.find_sids(cids[0])
        pr = self.canvas.coords(r)
        newlist.append(points[0][1])
        newlist.append(pr[1]+2) 
        l,r = self.find_sids(cids[1])
        pl = self.canvas.coords(l)
        newlist.append(points[1][1])
        newlist.append(pl[1]+2)
        rng = self.range_end-(self.range_begin)
        for p in newlist:
            result = self.range_begin+rng*((self.viewport[1]-p)/(self.viewport[1]-self.viewport[0]))
            if result<self.range_begin:
                result  = self.range_begin
            if result > self.range_end:
                result  = self.range_end
            values.append(result)
        return values
    
    
    def OnCanvas(self,event=None):
        x = event.x
        y = event.y
        c = self.canvas
        items = c.find_overlapping(x-2,y-2,x+2,y+2)
        if items == ():
            pn = c.find_withtag("point_online")
            if pn!=[]:
                for i in pn:
                    c.delete(i)
            tx = c.find_withtag("text_online")
            if tx!=[]:
                for i in tx:
                    c.delete(i)        
                   
    def OnLine(self,id,event=None):
        """when mouse is online creates oval and displays its value"""
        c = self.canvas    
        #del old
        pn = c.find_withtag("point_online")
        if pn!=[]:
            for i in pn:
                c.delete(i)
        tx = c.find_withtag("text_online")
        if tx!=[]:
            for i in tx:
                c.delete(i)
        x = event.x
        y = event.y
        self.canvas.create_oval(x-2,y-2,x+2,y+2,fill="white",tags=("point_online",))
        val = self.getValueforCoords((x,y))
        self.canvas.create_text((x+10,y+10),text=val,tags= ("text_online",))

if __name__ == "__main__":
    ###########
    import Tkinter
    sp = Spline([(0.0,0.0),(1.0,1.0)],[ (None, (1,-1)), ((1,1), None) ] )
    canvas  =  Tkinter.Canvas(width = 700, height = 700)
    canvas.pack()
    #canvas.configure(cursor="cross")
    canvas.create_rectangle([100,600,600,100])
    anim  =  AnimatorSpline(sp, canvas, (100, 600, 500, 500), 'abc')
    slopes = [ (None, (1,1)), ((1,1), None) ]
    anim.doit(50,80,5,100,slopes)

