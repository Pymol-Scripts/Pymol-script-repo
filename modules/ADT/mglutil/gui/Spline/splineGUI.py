
#
#$Id: splineGUI.py,v 1.21 2007/07/13 19:07:41 sowjanya Exp $


import Tkinter
from math import sqrt
from mglutil.util.callback import CallbackFunction
from mglutil.util.misc import deepCopySeq
from spline import Spline
import copy

class splineGUI:
    """Draw a spline that can be edited
the self.canvas on which to draw is provided as an argument
the viewport provides the coordinates on the self.canvas of the lower left corner
and the width and height of the drawing area.
uniqueID is a string that is used as a tag for self.canvas items for this GUI. It
should not be used by items on the self.canvas not created by this object.
doubleClick on control points to add or delete a point.
rightClick to tie or untie normals.
"""
    def __init__(self, spline, canvas, viewport, uniqueId, scale = None):
        assert isinstance(spline, Spline)
        self.spline  =  spline
        self.uniqueId  =  uniqueId
        
        assert isinstance(canvas, Tkinter.Canvas)
        self.canvas  =  canvas

        self.autoScale  =  False

        assert len(viewport) == 4
        self.setViewport(viewport)

        self.setScale(scale)
        
        self.pointRadius  =  2 # radius of a control point
        self.drawControlPoints()
        self._slopes_0n = 0
        self.history_cpts = []
        self.history_spline_points = []
        self.history_spline_slopes = []
        self.history_spline_slopes_lengths =[]
        self.slopes_lengths=[]
        self._var_ = 1
        self.right_dic={}
        self.left_dic={}
        self.canvas.configure(cursor="cross")
        self.drawSlopes()
        self.drawFunctionLine()

    def setScale(self, scale):
        # scale can be an (sx,sy) tuple or None in which case the scale
        # is computed automatically to contain the function
        if scale is None:
            self.doAutoScale()
        else:
            assert len(scale) == 2
            assert isinstance(scale[0], float)
            assert isinstance(scale[1], float)
            self.scale  =  scale

            
    def doAutoScale(self):
        # compute scaling factors to fit function in viewport
        self.autoScale  =  True
        ptx = []
        pty = []
        pts  =  self.spline.getControlPoints()
        for p in pts:
            ptx.append(p[0])
            pty.append(p[1])        
        if pty[-1]<=1.5:
            dx=1.0
            dy =1.1
        else:
             minx  =  min(ptx)
             maxx  =  max(ptx)
             dx  =  maxx-minx
             miny  =  min(pty)*1.1 # FIXME we should really evaluate the fucntion
             maxy  =  max(pty)*1.1 # to find real min and max for y
             dy  =  maxy-miny
        self.scale  = (self.width/float(dx), self.height/float(dy))
        
            

    def setViewport(self, viewport):
        # portion of the self.canvas on which the function will be drawn
        x,y,w,h  =  viewport
        self.originX  =  x
        self.originY  =  y
        self.width  =  w
        self.height  =  h
        self.viewport  =  viewport
        if self.autoScale:
            self.doAutoScale()

    
    def drawControlPoints(self):
        self.right_dic={}
        self.left_dic={}
        canvas  =  self.canvas
        pts  =  self.getControlPoints()
        ptRad  =  self.pointRadius
        sx, sy  =  self.scale
        ox  =  self.originX
        oy  =  self.originY
        uid  =  self.uniqueId
        tag  =  'ctrlPoints_'+uid
        index  =  0
        
        
        for (x, y) in pts:
         
         rx = x
         ry = y
                
               
            
         
            
         
         id  =  self.canvas.create_oval(rx-ptRad, ry-ptRad, rx+ptRad, ry+ptRad,
                                    fill = 'black', tags = (tag, uid,
                                                        'cp'+str(index)+uid))
         cb  =  CallbackFunction(self.startMoveControlPoint, id, index)
         self.canvas.tag_bind(id, '<ButtonPress-1>', cb)
         cb  =  CallbackFunction( self.endModeControlPoint, id)
         self.canvas.tag_bind(id, '<ButtonRelease-1>', cb)
         cb  =  CallbackFunction(self.deletePoint, id)
         self.canvas.tag_bind(id, '<Double-Button-1>', cb)
         cb = CallbackFunction(self.tieNormals, id)
         self.canvas.tag_bind(id, '<ButtonPress-3>', cb)
         index += 1
                      

    def removeControlPoints(self):
        c  =  self.canvas
        c.delete('ctrlPoints_'+self.uniqueId)
        
    def startEditingCurve(self):
        
        c  =  self.canvas
        c.itemconfigure('ctrlPoints_'+self.uniqueId, fill = 'red')
        self.drawSlopes()


    
    def drawSlopes(self):
        self._slopes_0n = 1
        
        canvas  =  self.canvas
        pts  =  self.getControlPoints()
        slopes  =  self.getSlopesVectors()
        ptRad  =  self.pointRadius
        sx, sy  =  self.scale
        ox  =  self.originX
        oy  =  self.originY
        uid  =  self.uniqueId
        tag  =  'slopes_'+uid
        
        self._points=[]        
        #index = 0 
        for p, slope in zip(pts, slopes):
            left_slope = None
            right_slope = None
            ind = slopes.index(slope)
            rx = p[0]
            ry = p[1]
            
            
                
            sll, slr  =  slope # unpack left and right slope
            if sll: # draw vector on left side of the point
                vx, vy  =  sll
                if len(self.slopes_lengths)<ind+1:
                    n  =  30/sqrt(vx*vx+vy*vy)
                else:
                     n  =  30/sqrt(vx*vx+vy*vy)
                     
                     if round(n)!=round(self.slopes_lengths[ind][0]):
                        n= self.slopes_lengths[ind][0]/sqrt(vx*vx+vy*vy)
                vx *=  n
                vy *=  n
                Ry = ry-vy
                Rx = rx-vx
                if Rx>rx:
                    Rx = rx+vx
                if Rx<self.originX:
                    Rx = self.originX
                elif Rx>self.originY:
                    Rx = self.originY
                if Ry<self.originX:
                    Ry = self.originX
                elif Ry>self.originY:
                    Ry=self.originY             
                    
                lid  =  self.canvas.create_line(Rx, Ry, rx, ry,
                                        fill = 'black', tags = (tag,
                                        uid,"left","ln"))
                tid = self.canvas.create_text((Rx+10, Ry+10),text = round(vy/vx,1), tags = (tag, uid,"left","tx"))
                sid  =  self.canvas.create_oval(Rx-ptRad, Ry-ptRad,
                                        Rx+ptRad, Ry+ptRad,
                                        fill = 'green', tags = (tag,
                                        uid,"left","pt"))
                cb  =  CallbackFunction(self.startMoveSlopePoint, sid, lid,tid)
                self.canvas.tag_bind(sid, '<ButtonPress-1>', cb)
                cb  =  CallbackFunction( self.endMoveSlopePoint, sid,lid,tid)
                self.canvas.tag_bind(sid, '<ButtonRelease-1>', cb)
                
                left_slope = n
                self._points.append((Rx,Ry))    
                
                self._points.append((rx,ry))
                
                    
            if slr:
                vx, vy  =  slr
                if len(self.slopes_lengths)<ind+1:
                    n  =  30/sqrt(vx*vx+vy*vy)
                else:
                    n  =  30/sqrt(vx*vx+vy*vy)
                    if round(n)!=round(self.slopes_lengths[ind][1]):
                        n= self.slopes_lengths[ind][1]/sqrt(vx*vx+vy*vy)
                    
                
                
                vx *=  n
                vy *=  n
                Ry = ry+vy
                Rx = rx+vx
                if Rx<rx:
                    Rx = rx-vx
                if Rx>rx:
                    Rx = rx+vx
                if Rx<self.originX:
                    Rx =self.originX
                elif Rx>self.originY:
                    Rx = self.originY
                if Ry<self.originX:
                    Ry =self.originX
                elif Ry>self.originY:
                    Ry=self.originY
                
                lid  =  self.canvas.create_line(rx, ry, Rx, Ry,
                                        fill = 'black', tags = (tag,
                                        uid,"right","ln"))
                tid = self.canvas.create_text((Rx+10, Ry+10),text = round(vy/vx,1), tags = (tag, uid,"right","tx"))
                
                sid  =  self.canvas.create_oval(Rx-ptRad, Ry-ptRad,
                                        Rx+ptRad, Ry+ptRad,
                                        fill = 'green', tags = (tag,
                                        uid,"right","pt"))
                cb  =  CallbackFunction(self.startMoveSlopePoint, sid, lid,tid)
                self.canvas.tag_bind(sid, '<ButtonPress-1>', cb)
                cb  =  CallbackFunction( self.endMoveSlopePoint, sid,lid,tid)
                right_slope = n
                self._points.append((rx,ry))
                self._points.append((Rx,Ry))
            if len(self.slopes_lengths)!=len(slopes):
                self.slopes_lengths.append((left_slope,right_slope))
        self.canvas.tag_raise("pt")
        self.canvas.tag_raise('ctrlPoints_'+self.uniqueId)
         

    def stopEditingCurve(self):
        self._slopes_0n = 0
        c  =  self.canvas
        c.itemconfigure('ctrlPoints_'+self.uniqueId, fill = 'black')
        c.delete('slopes_'+self.uniqueId)
    
    def endModeControlPoint(self, cid, index):
        c  =  self.canvas
        c.tag_unbind(cid, '<B1-Motion>')
        self.appendToHistory() 
         

    def findSplinePoints(self):
        c=self.canvas
        spline_points=[]
        self.spline_points=[]
        self.n_points=[] 
        time_steps=[]
         
                
        
        for i in range(len(self._points)):
          if i%4==0:
            if i+3 <=len(self._points)-1  :
                
                p0 = self._points[i]
                p1 = self._points[i+1]
                p2 = self._points[i+2]
                p3 = self._points[i+3]
             
                    
                
                self.n_points.append([p0,p1,p2,p3])
                
        for j in range(10):
            time_steps.append(float(j)/10)
        for p in self.n_points:
          
          p0,p1,p2,p3=p
          ind_p=self.n_points.index(p)
          curve_points=[]
          result_points=[]
          for  t in   time_steps:
                ind_t = time_steps.index(t)
                res = self.spline.beizer_calc(t,p)
                result_points.append(res)
                curve_points.append(res)
          if p0 not in curve_points:
            curve_points.append(p0)
          if hasattr(self,"id"):
            #right
            if self.id == "right":
                #for  i in   range(0,len(result_points)):
                #    if i+1 <= len(result_points)-1:
                      if result_points[4][0]<result_points[5][0] :
                          f = ind_p*11#(len(time_steps))
                          l = (ind_p+1)*11#(len(time_steps))
                          required_points = self.history_spline_points[-1][f:l]
                          required_points.reverse()
                          curve_points=[]
                          for q in required_points:
                              curve_points.append(q)
                          ind = self._points.index(p1)
                          self._points.remove(self._points[ind])
                          self._points.insert(ind,p1)
                          self.point_p1 =  self._points[ind]
                          if hasattr(self,"sid"):
                              sid = self.sid
                              if len(c.coords(sid))>0:
                                  [x1,y1,x2,y2] = c.coords(sid)
                          
                                  if self.point_p1[0]>x1 and self.point_p1[0]<x2:
                                      n_x =    self.point_p1[0]
                                      if self.point_p1[1]>y1 and self.point_p1[1]<y2:
                                          n_y =    self.point_p1[1]    
                                          self.right_dic[sid]=(n_x,n_y)   
                                               
                          #break
            #result_points.reverse()
            if self.id == "left":
            #    for  i in   range(0,len(result_points)):
            #         if i+1 <= len(result_points)-1:
                      if result_points[7][0]>result_points[6][0]:
                            
                          f = ind_p*11#(len(time_steps))
                          l = (ind_p+1)*11#(len(time_steps))
                          required_points = self.history_spline_points[-1][f:l]
                          required_points.reverse()
                          curve_points=[]
                          for q in required_points:
                              curve_points.append(q)   
                          ind = self._points.index(p2)
                          point_p2 =  self._points[ind]
                          if hasattr(self,"sid"):
                              sid = self.sid
                              if len(c.coords(sid))>0:
                                  [x1,y1,x2,y2] = c.coords(sid)
                                  if point_p2[0]>x1 and point_p2[0]<x2:
                                      n_x =    point_p2[0]
                                      if point_p2[1]>y1 and point_p2[1]<y2:
                                          n_y =    point_p2[1]    
                                          self.left_dic[sid]=(n_x,n_y)
                          #break
                        
            #result_points.reverse() 
          
          spline_points.append(curve_points)
        for r in spline_points:
            r.reverse()
            for i in r:
                self.spline_points.append(i)
        
        self.history_spline_points.append(deepCopySeq(self.spline_points))
        
        return  self.spline_points

    def drawFunctionLine(self):
        self._line_On = 1
        c = self.canvas
        uid  =  self.uniqueId
        tag  =  'ctrlPoints_'+uid
        tag1 =  'line_'+uid
        spline_points = self.findSplinePoints()
        self.line = id = c.create_line(spline_points,tags = (tag1, uid))
        cb  =  CallbackFunction(self.createPoint, id)
        self.canvas.tag_bind(id, '<Double-Button-1>', cb)
        self.canvas.tag_raise('ctrlPoints_'+self.uniqueId) 
        if self.history_cpts==[]:    
                self.appendToHistory()
        #if self.history_cpts==[]:
        #    self.history_cpts.append(deepCopySeq(self.getControlPoints()))
        #    self.history_spline_slopes.append((deepCopySeq(self.getSlopesVectors()),deepCopySeq(self.slopes_lengths))) 

    def removeFunctionLine(self):
        self._line_On = 0
        c  =  self.canvas
        c.delete('line_'+self.uniqueId)    
    
    
    ##
    ## callbacks
    ##
    def startMoveControlPoint(self, cid, index, event):
        # cid of control point and control point index
        
        c  =  self.canvas
        cb =  CallbackFunction(self.moveControlPoint, cid, index)
        c.tag_bind(cid, '<B1-Motion>', cb)
        self.origx = event.x
        self.origy = event.y
        cpts=self.getControlPoints()
        rounded_cpts = []
        ind_cpt=None 
        
        for p in cpts:
            rounded_cpts.append((round(p[0]),round(p[1])))
        
        if (round(c.coords(cid)[0]+2),round(c.coords(cid)[1]+2)) in rounded_cpts:
            ind_cpt = rounded_cpts.index((round(c.coords(cid)[0]+2),round(c.coords(cid)[1]+2)))
        if (round(c.coords(cid)[2]),round(c.coords(cid)[3])) in rounded_cpts:
            ind_cpt = rounded_cpts.index((round(c.coords(cid)[2]),round(c.coords(cid)[3])))
        if ind_cpt:
            if ind_cpt-1 in range(len(cpts)):
                self.p_cpt = rounded_cpts[ind_cpt -1]
            if ind_cpt+1 in range(len(cpts)):
                self.n_cpt = rounded_cpts[ind_cpt +1]
        
           
    def moveControlPoint(self, cid, index, event):
        
        c  =  self.canvas
        uid  =  self.uniqueId
        tag  =  'ctrlPoints_'+uid
        items = c.find_withtag(tag)
        x  =  c.canvasx(event.x)
        y  =  c.canvasy(event.y)
        coords = c.coords(cid)
        cx1,cy1=(coords[0],coords[1])
        #fixing to limit cpts to move with in Viewport
        if y>self.originY:
            y=self.originY
        if y<self.originX:
            y=self.originX
        #making start ,end points to move along y-axis    
        if index == 0 or index == len(items)-1:
            dx = 0
            if index == 0:
                x = self.originX
            else:
                x = self.originY
        else:
            dx = x-self.origx
        if hasattr(self,"p_cpt"):
            if x<self.p_cpt[0]:
                if cx1>self.p_cpt[0]:
                    x= self.p_cpt[0]+3
                    dx = x-cx1+2
                else:
                    dx = 0
        if hasattr(self,"n_cpt"):
            if x>self.n_cpt[0]:
                if cx1<self.n_cpt[0]:
                    x=self.n_cpt[0]-3
                    dx = x-cx1-2
                else:
                    dx = 0
        if y ==  self.originY or y ==  self.originX:
            dy = 0
        else:
            dy = y - self.origy
        c.move(cid,dx,dy)
        self.origx = event.x
        self.origy = event.y
        #setting splines Controlpoints with new points.Such that
        #spline.getControlPoints returns new control points
        cpts = []
        sx, sy  =  self.scale
        ox  =  self.originX
        oy  =  self.originY
        for i in items:
            cx = c.coords(i)[0]
            cy = c.coords(i)[1]
            cpts.append((cx+2,cy+2))
         
        self.spline.setControlPoints(cpts,self.getSlopesVectors())
        self.id=None
        self._drawSpline()
        cids = c.find_withtag("ctrlPoints_"+self.uniqueId)    
        for cid in cids:
            if "tied" in c.gettags(cid):
                lsid,rsid = self.find_sids(cid)
                if lsid!=None and rsid!=None:
                    c.itemconfigure(lsid, fill = 'blue')
                    c.itemconfigure(rsid, fill = 'blue')
        

    def startMoveSlopePoint(self, sid, lid,tid, event):
        
        c  =  self.canvas
        self.origx = event.x
        self.origy = event.y
        
        self.line_coords = c.coords(lid)    
        cb =  CallbackFunction(self.moveSlopePoint, sid,lid, tid)
        c.tag_bind(sid, '<B1-Motion>', cb)
        
        
        self._var_=1 
             
        
        
    def moveSlopePoint(self ,sid, lid,tid, event):
        
        c  =  self.canvas
        x  =  c.canvasx(event.x)
        y  =  c.canvasy(event.y)
        other_sid,cid =  self.find_othersid(sid)
        
        #when tied ,controlling other_side slope point not to cross viewport limits
        if "left" in c.gettags(sid):
                c_rid, c_cid= self.find_othersid(sid)
                cid = c_cid
                self.line_coords = (c.coords(sid)[0]+2,c.coords(sid)[1]+2,c.coords(c_cid)[0]+2,c.coords(c_cid)[1]+2)
                if c_cid:
                 if "tied" in c.gettags(c_cid):
                    if c_rid:
                        if c.coords(c_rid)[0]>=self.originY-2:
                            if x<=c.coords(sid)[0]:
                                x = c.coords(sid)[0]+2
                        if c.coords(c_rid)[1]<=self.originX+2:
                            if y>=c.coords(sid)[1]:
                               y = c.coords(sid)[1]+2 
                        if c.coords(c_rid)[1]>=self.originY-2:
                            if y<=c.coords(sid)[1]+2:
                               y = c.coords(sid)[1]+2 
                                
        if "right" in c.gettags(sid):
                c_lid , c_cid= self.find_othersid(sid)
                cid = c_cid
                self.line_coords = (c.coords(c_cid)[0]+2,c.coords(c_cid)[1]+2,c.coords(sid)[0]+2,c.coords(sid)[1]+2)   
                if c_cid:
                    if "tied" in c.gettags(c_cid):
                        if c_lid:
                            if c.coords(c_lid)[0]<=self.originX+2:
                                if x>=c.coords(sid)[0]:
                                    x = c.coords(sid)[0]+2
                            if c.coords(c_lid)[1]<=self.originX+2:
                                if y>=c.coords(sid)[1]:
                                    y = c.coords(sid)[1]+2 
                            if c.coords(c_lid)[1]>=self.originY-2:
                                if y<=c.coords(sid)[1]+2:
                                    y = c.coords(sid)[1]+2 
        c.delete(lid)
        c.delete(tid)
        if hasattr(self,"sid"):
                self.old_sid = self.sid
        if hasattr(self,"lid"):
           self.old_lid=self.lid
           self.old_tid =self.tid
           if self._var_ == 0:
                c.delete(self.lid)
                c.delete(self.tid)
        uid  =  self.uniqueId
        tag  =  'slopes_'+uid     
        tag1 =  'ctrlPoints_'+uid 
        current_points = self.getControlPoints()
        cpts = []
        for i in current_points:
            cx = i[0]
            cy = i[1]
            cpts.append((round(cx),round(cy)))
            
        
        slope_vectors = self.getSlopesVectors()
        cid_coords = []
        for i in c.find_withtag("ctrlPoints_"+self.uniqueId):
            cid_coords.append(c.coords(i))
        if "left" in c.gettags(sid):
            
                    
            
            self.id = "left"
            if self.left_dic!={}:
                
                if sid in self.left_dic.keys():
                    
                    (px,py) = self.left_dic[sid]
                    self.left_dic.pop(sid)
                    if (round(px),round(py))==(round(self.line_coords[0]),round(self.line_coords[1])):
                         
                         if (x,y)<(px,py): 
                            (x,y)=(px,py)
                        
             
            #limiting left slope in one direction
            if x>self.line_coords[2]:
                x = self.line_coords[2]
            #when slope is infinity
            if self.line_coords[2] - x == 0:
                
                x = self.line_coords[2]-2
            
            dx = x-self.line_coords[0]
            dy = y-self.line_coords[1] 
            
            
                    
                       
            #moving point
            #limiting slope points with in view port
            cx,cy = (c.coords(sid)[0]+2,c.coords(sid)[1]+2)
            if x<=self.originX:
                if cx>=self.originX :
                    x = self.originX+2
                    dx = x - cx
                     
                else:
                    dx = 0
                    x = self.originX+2
                    
            if y>=self.originY:
                if cy<=self.originY:
                    y = self.originY-2
                    
                    dy = y-cy
                else:
                    dy =0
                    
                    y = self.originY-2
            if y<=self.originX:
                if cy>=self.originX:
                    y = self.originX+2
                    
                    dy = y - cy
                else:
                    dy =0
                    y = self.originX+2
                    
            
            #if sid coords ==  cid coords
            for i in cid_coords:
                if (round(x),round(y)) == (round(i[0]+2),round(i[1]+2)):
                    x =c.coords(sid)[0]+2+4
                    y = c.coords(sid)[1]+2+4
                    dx=4
                    dy =4
                if (round(x),round(y)) == (round(i[0]),round(i[1])):    
                    x =c.coords(sid)[0]+2+4
                    y = c.coords(sid)[1]+2+4
                    dx = 4
                    dy = 4
                    
                            
                                
                    
                if (round(x),round(y)) == (round(i[0]+4),round(i[1]+4)):    
                    x =c.coords(sid)[0]+2+4
                    y = c.coords(sid)[1]+2+4
                    dx = 4
                    dy = 4
                
            #checking if tied other sid coords crosses viewport limits before
            #moving slope point
            if cid:
                if "tied" in c.gettags(cid):
                   coords_other_sid = c.coords(other_sid)
                   if (coords_other_sid[0]+2)-dx>=self.originY:
                        dx =  -1*(self.originY - (coords_other_sid[0]+2))
                        x = c.coords(sid)[0]+dx+2
                   if (coords_other_sid[1]+2)-dy<=self.originX:
                        dy =  ((coords_other_sid[1]+2)-self.originX)
                        y = c.coords(sid)[1]+dy+2 
                   if (coords_other_sid[1]+2)-dy>=self.originY:
                        dy =  -1*(self.originY - (coords_other_sid[1]+2))
                        y = c.coords(sid)[1]+dy+2  
                          
            if x>self.line_coords[2]-6 :
               if y in range(int(round(self.line_coords[3]-6)) ,int(round(self.line_coords[3]+6))):
                    x=self.line_coords[2]-6
                    dx = (self.line_coords[2]-6) - (c.coords(sid)[0]+2)
                    if dx%2!=0:
                        x =x+1
                        dx = dx+1
                        
                    #d = x - (self.line_coords[2]-6 )
                    #d = 6-diff
                    #if d!=0:
                    #    if d%2!=0:
                    #        x=c.coords(sid)[0]+2-(d+1)
                    #        dx = -(d+1)
                    #        
                    #    else:
                    #        x = c.coords(sid)[0]+2-d
                    #        dx=-d
            c.move(sid,dx,dy)
            #x,y = (c.coords(sid)[0]+2,c.coords(sid)[1]+2)
            #drawing line and text
            self.lid = c.create_line((x,y),(self.line_coords[2],self.line_coords[3]),tags = (tag,uid,"left","ln"))
            
            _current_slope_length = sqrt((self.line_coords[2]-x)*(self.line_coords[2]-x)+(self.line_coords[3]-y)*(self.line_coords[3]-y))
            _current_slope  =  (self.line_coords[3] - y)/(self.line_coords[2] - x)
            self.tid = self.canvas.create_text((x+10,y+10),text = round(_current_slope,1), tags = (tag, uid,"left","tx"))
            #considering radius of point
            if  (round(self.line_coords[2]),round(self.line_coords[3])) in cpts :
                _point_index  =  cpts.index((round(self.line_coords[2]),round(self.line_coords[3])))
            else:
                _point_index  =  cpts.index((round(self.line_coords[2]-2),round(self.line_coords[3]-2)))
            
            right_slope = slope_vectors[_point_index][1]
            slope_vectors.remove(slope_vectors[_point_index])
            
            cpt=(round(self.line_coords[2]),round(self.line_coords[3]))
            
            #updating self._points
            rx,ry = (round(self.line_coords[0]),round(self.line_coords[1]))
            new_p=[]
            
            for i in  self._points:
                new_p.append((round(i[0]),round(i[1])))
            
            if (rx,ry) in new_p : 
                ind = new_p.index((rx,ry))
            if (rx-2,ry-2)  in new_p :    
                ind = new_p.index((rx-2,ry-2))
            if (rx+2,ry+2)  in new_p :    
                ind = new_p.index((rx+2,ry+2))
            
            
            Ry =y 
            Rx = x  
            self._points.remove(self._points[ind])
            self._points.insert(ind,(Rx,Ry))
            
            #finding current cid        
            cids=c.find_withtag(tag1)
            coords_cids_rounded  = []
            for i in cids:  
                coords_cids = c.coords(i)
                coords_cids_rounded.append((round(coords_cids[0]+2),round(coords_cids[1]+2),round(coords_cids[2]),round(coords_cids[3])))
            current_id = None    
            #finding cid
            for j in coords_cids_rounded:
                
                if (j[0],j[1])==cpt:
                    cpt = (j[0]+2,j[1]+2)
                    ind = coords_cids_rounded.index(j)
                    current_id= cids[ind] 
                elif (j[2],j[3])==cpt:
                    cpt =(j[2],j[3])
                    ind = coords_cids_rounded.index(j)
                    if ind >0 and ind <len(coords_cids_rounded)-1:
                        current_id= cids[ind]
                  
            ##finding lid,sid,tid with same cpt
            for i in c.find_withtag("pt"):
                if i>sid:
                    right_sid = i
                    break
            
            

            if current_id:
              if "tied" in c.gettags(current_id): 
                #finding and deleting if right lid and right tid exists
                
                
                if sid+1 in c.find_all():
                    r_lid = sid+1
                    r_tid = sid+2
                    c.delete(r_lid)
                    c.delete(r_tid)
                
                right_sl = c.find_withtag("right")
                for l,t in zip(c.find_withtag("ln"),c.find_withtag("tx")):
                    if "right" in c.gettags(l)  and "right" in c.gettags(t):
                        
                        if (round(c.coords(l)[0]),round(c.coords(l)[1]))==(round(c.coords(current_id)[0]+2),round(c.coords(current_id )[1]+2)):
                            c.delete(l)
                            c.delete(t)
                            
                            
                            
                if hasattr(self,"old_lid"):
                    coords_l = c.coords(self.old_lid)
                    if coords_l!=[]:
                        if (round(coords_l[0]),round(coords_l[1])) == (round(c.coords(current_id)[0]+2),round(c.coords(current_id )[1]+2)):
                           c.delete(self.old_lid)
                           c.delete(self.old_tid)
                
                if hasattr(self,"right_lid"):
                    coords_l = c.coords(self.right_lid)
                    if coords_l!=[]:
                        if (round(coords_l[0]),round(coords_l[1])) == (round(c.coords(current_id)[0]+2),round(c.coords(current_id )[1]+2)):
                            
                            c.delete(self.right_lid)
                            c.delete(self.right_tid)
                
                right_sid_coords=c.coords(right_sid)
                
                #moving point
                #limiting slope point to move with in viewport
                   
                c.move(right_sid,-dx,-dy)
                new_coords = c.coords(right_sid)
                
                self.right_lid = c.create_line((self.line_coords[2],self.line_coords[3]),(new_coords[2]-2,new_coords[3]-2),tags = (tag,uid,"right","ln"))
                
                self.right_tid = self.canvas.create_text((new_coords[2]-2+10,new_coords[3]-2+10),text = round(_current_slope,1), tags = (tag, uid,"right","tx"))
                slope_vectors.insert(_point_index,((1,_current_slope),(1,_current_slope)))
                #r_slope = self.slopes_lengths[_point_index][1]
                self.slopes_lengths.remove(self.slopes_lengths[_point_index])
                self.slopes_lengths.insert(_point_index,(_current_slope_length,_current_slope_length))
                cb  =  CallbackFunction(self.startMoveSlopePoint, right_sid, self.right_lid,self.right_tid)
                self.canvas.tag_bind(right_sid, '<ButtonPress-1>', cb)
                cb  =  CallbackFunction( self.endMoveSlopePoint, right_sid,self.right_lid,self.right_tid)
                self.canvas.tag_bind(right_sid, '<ButtonRelease-1>', cb)
                
                #updating self._points
                if right_sid_coords:
                    rx,ry = (round(right_sid_coords[2]),round(right_sid_coords[3]))
                    new_p=[]
                    
                    for i in  self._points:
                        new_p.append((round(i[0]),round(i[1])))
                    
                    if (rx,ry) in new_p : 
                        ind1 = new_p.index((rx,ry))
                    if (rx-2,ry-2)  in new_p :    
                        ind1 = new_p.index((rx-2,ry-2))
                    if (rx-4,ry-4)  in new_p :    
                        ind1 = new_p.index((rx-4,ry-4))
                    if (rx-4,ry-2)  in new_p :    
                        ind1 = new_p.index((rx-4,ry-2))
                    if (rx-2,ry-4)  in new_p :    
                        ind1 = new_p.index((rx-2,ry-4))
                    Ry =new_coords[3]
                    Rx = new_coords[2]
                    self._points.remove(self._points[ind1])
                    self._points.insert(ind1,(Rx,Ry))
              else:
                
                slope_vectors.insert(_point_index,((1,_current_slope),right_slope)) 
                r_slope = self.slopes_lengths[_point_index][1]
                self.slopes_lengths.remove(self.slopes_lengths[_point_index])
                self.slopes_lengths.insert(_point_index,(_current_slope_length,r_slope))
                #c.itemconfigure(sid, fill = 'green')
                
            else:
                r_slope = self.slopes_lengths[_point_index][1]
                self.slopes_lengths.remove(self.slopes_lengths[_point_index])
                self.slopes_lengths.insert(_point_index,(_current_slope_length,r_slope))
                slope_vectors.insert(_point_index,((1,_current_slope),right_slope))  
            
            
            
        
        if "right" in c.gettags(sid):
            self.id = "right"
            
            if self.right_dic!={}:
                
                if sid in self.right_dic.keys():
                    
                    (px,py) = self.right_dic[sid]
                    self.right_dic.pop(sid)
                    if (round(px),round(py))==(round(self.line_coords[2]),round(self.line_coords[3])):
                         #if x>px:
                         
                         if (x,y)>(px,py): 
                            (x,y)=(px,py)
            
            #limiting right slope in one direction
            if x<self.line_coords[0]:
                x = self.line_coords[0]
            #when slope is infinity
            if x-self.line_coords[0] == 0:
                x = self.line_coords[0]+2
                
            
            dx  =  x-self.line_coords[2]
            dy  =  y-self.line_coords[3]    
            
            #moving point
            #limiting slope point to move with in viewport
            cx,cy = (c.coords(sid)[0]+2,c.coords(sid)[1]+2)
            if x>=self.originY:
                if cx<=self.originY :
                    x = self.originY-2
                    dx = x-cx 
                else:
                    dx = 0
                    x = self.originY-2
            if y>=self.originY:
                if cy<=self.originY:
                    y = self.originY-2
                    dy = y-cy
                else:
                    dy =0
                    y = self.originY-2
            if y<=self.originX:
                if cy>=self.originX:
                    y = self.originX+2
                    dy = y - cy
                else:
                    dy =0
                    y = self.originX+2
                    
            
            #if sid coords ==  cid coords        
            for i in cid_coords:
                if (round(x),round(y)) == (round(i[0]+2),round(i[1]+2)):
                    x =c.coords(sid)[0]+2+4
                    y = c.coords(sid)[1]+2+4
                    dx=4
                    dy =4
                     
                
                if (round(x),round(y)) == (round(i[0]),round(i[1])): 
                    x =c.coords(sid)[0]+2+4
                    y = c.coords(sid)[1]+2+4
                    dx = 4
                    dy = 4
                if (round(x),round(y)) == (round(i[0]+4),round(i[1]+4)):  
                    x =c.coords(sid)[0]+2+4
                    y = c.coords(sid)[1]+2+4
                    dx = 4
                    dy = 4
            
            #checking if tied other sid coords crosses viewport limits
            if cid:
               
                
                if "tied" in c.gettags(cid):
                   coords_other_sid = c.coords(other_sid)
                   if (coords_other_sid[0]+2)-dx<=self.originX:
                        dx =  ((coords_other_sid[0]+2)-self.originX)
                        x = c.coords(sid)[0]+dx+2
                   if (coords_other_sid[1]+2)-dy<=self.originX:
                        dy =  ((coords_other_sid[1]+2)-self.originX)
                        y = c.coords(sid)[1]+dy+2 
                   if (coords_other_sid[1]+2)-dy>=self.originY:
                        dy =  -1*(self.originY - (coords_other_sid[1]+2))
                        y = c.coords(sid)[1]+dy+2 
            
            
            
            #if x<self.line_coords[0]+6 :
            #   if y in range(int(round(self.line_coords[1]-6)) ,int(round(self.line_coords[1]+6))):
            #        x=self.line_coords[0]+6
            #        dx = (c.coords(sid)[0]+2) - (self.line_coords[0]+6) 
            #        if dx%2!=0:
            #            x=x+1
            #            dx = dx+1
                        
            c.move(sid,dx,dy)
            
            #drawing line and text
            self.lid  =  c.create_line((self.line_coords[0],self.line_coords[1]),(x,y),tags = (tag,uid,"right","ln"))
            _current_slope_length = sqrt((x-self.line_coords[0])*(x-self.line_coords[0])+(y-self.line_coords[1])*(y-self.line_coords[1]))
            _current_slope  =  (y-self.line_coords[1])/(x-self.line_coords[0])
            self.tid = c.create_text((x+10,y+10),text = round(_current_slope,1), tags = (tag, uid,"right","tx"))
            #considering radius of point
            if (round(self.line_coords[0]),round(self.line_coords[1])) in cpts:
                _point_index  =  cpts.index((round(self.line_coords[0]),round(self.line_coords[1])))
            else:    
                _point_index  =  cpts.index((round(self.line_coords[0]-2),round(self.line_coords[1]-2)))
            
            left_slope = slope_vectors[_point_index][0]
            slope_vectors.remove(slope_vectors[_point_index])
            
            cpt=(round(self.line_coords[0]),round(self.line_coords[1]))
            
            #updating self._points which is used to compute beizercurve
            rx,ry = (round(self.line_coords[2]),round(self.line_coords[3]))
            new_p=[]
            for i in  self._points:
                new_p.append((round(i[0]),round(i[1])))
             
            if (rx,ry) in new_p : 
                ind = new_p.index((rx,ry))
            if (rx-2,ry-2)  in new_p:
                ind = new_p.index((rx-2,ry-2))
            if (rx+2,ry+2)  in new_p:
                ind = new_p.index((rx+2,ry+2))
            Ry =y 
            Rx = x
            self._points.remove(self._points[ind])
            self._points.insert(ind,(Rx,Ry))          
            
            
            #finding current cid        
            cids=c.find_withtag(tag1)
            coords_cids_rounded  = []
            for i in cids:  
                coords_cids = c.coords(i)
                coords_cids_rounded.append((round(coords_cids[0]+2),round(coords_cids[1]+2),round(coords_cids[2]),round(coords_cids[3])))
            current_id = None    
            #finding cid
            for j in coords_cids_rounded:
                if (j[0],j[1])==cpt:
                    cpt = (j[0]+2,j[1]+2)
                    ind = coords_cids_rounded.index(j)
                    current_id= cids[ind] 
                elif (j[2],j[3])==cpt:
                    cpt =(j[2],j[3])
                    ind = coords_cids_rounded.index(j)
                    if ind >0 and ind <len(coords_cids_rounded)-1:
                        current_id= cids[ind]
                       
            ##finding lid,sid,tid with same cpt
            for i in c.find_withtag("pt"):
                if i == sid:
                    ind = list(c.find_withtag("pt")).index(sid)
                    left_sid=c.find_withtag("pt")[ind-1]
       

            if current_id:
              if "tied" in c.gettags(current_id): 
                #finding and deleting if left lid and left tid exists
                if sid-4 in c.find_all():
                    l_lid = sid-4
                    l_tid = sid-5
                    c.delete(l_lid)
                    c.delete(l_tid)
                for l,t in zip(c.find_withtag("ln"),c.find_withtag("tx")):
                    if "left" in c.gettags(l)  and "left" in c.gettags(t):
                        if (round(c.coords(l)[2]),round(c.coords(l)[3]))==(round(c.coords(current_id)[0]+2),round(c.coords(current_id )[1]+2)):
                            c.delete(l)
                            c.delete(t)
                if hasattr(self,"left_lid"):
                    coords_l = c.coords(self.left_lid)
                    if coords_l!=[]:
                        if (round(coords_l[2]),round(coords_l[3])) == (round(c.coords(current_id)[0]+2),round(c.coords(current_id )[1]+2)):
                            c.delete(self.left_lid)
                            c.delete(self.left_tid)
                if hasattr(self,"old_lid"):
                    coords_l = c.coords(self.old_lid)
                    if coords_l!=[]:
                        if (round(coords_l[2]),round(coords_l[3])) == (round(c.coords(current_id)[0]+2),round(c.coords(current_id )[1]+2)):
                            c.delete(self.old_lid)
                            c.delete(self.old_tid)
                
                #moving point
                #limiting slope points with in view port
                left_sid_coords = c.coords(left_sid)
                c.move(left_sid,-dx,-dy)
                new_coords = c.coords(left_sid)
                
                self.left_lid = c.create_line((new_coords[2]-2,new_coords[3]-2),(self.line_coords[0],self.line_coords[1]),tags = (tag,uid,"left","ln"))
            
                self.left_tid = self.canvas.create_text((new_coords[2]-2+10,new_coords[3]-2+10),text = round(_current_slope,1), tags = (tag, uid,"left","tx"))
                slope_vectors.insert(_point_index,((1,_current_slope),(1,_current_slope)))
                self.slopes_lengths.remove(self.slopes_lengths[_point_index])
                self.slopes_lengths.insert(_point_index,(_current_slope_length,_current_slope_length))
                cb  =  CallbackFunction(self.startMoveSlopePoint, left_sid, self.left_lid,self.left_tid)
                self.canvas.tag_bind(left_sid, '<ButtonPress-1>', cb)
                cb  =  CallbackFunction( self.endMoveSlopePoint, left_sid,self.left_lid,self.left_tid)
                self.canvas.tag_bind(left_sid, '<ButtonRelease-1>', cb) 
                
                #updating self._points
                if left_sid_coords:
                    rx,ry = (round(left_sid_coords[2]),round(left_sid_coords[3]))
                    new_p=[]
                    for i in  self._points:
                        new_p.append((round(i[0]),round(i[1])))
                    
                    if (rx,ry) in new_p : 
                        ind1 = new_p.index((rx,ry))
                    if (rx-2,ry-2)  in new_p :    
                        ind1 = new_p.index((rx-2,ry-2))
                    if (rx-4,ry-4)  in new_p :    
                        ind1 = new_p.index((rx-4,ry-4))
                    if (rx-4,ry-2)  in new_p :    
                        ind1 = new_p.index((rx-4,ry-2))
                    if (rx-2,ry-4)  in new_p :    
                        ind1 = new_p.index((rx-2,ry-4)) 
                    Ry =new_coords[3]
                    Rx = new_coords[2]
                    self._points.remove(self._points[ind1])
                    self._points.insert(ind1,(Rx,Ry))
                

              else:
                    l_slope = self.slopes_lengths[_point_index][0]
                    self.slopes_lengths.remove(self.slopes_lengths[_point_index])
                    self.slopes_lengths.insert(_point_index,(l_slope,_current_slope_length))
                    slope_vectors.insert(_point_index,(left_slope,(1,_current_slope)))
                    
            else:
                l_slope = self.slopes_lengths[_point_index][0]
                self.slopes_lengths.remove(self.slopes_lengths[_point_index])
                self.slopes_lengths.insert(_point_index,(l_slope,_current_slope_length))
                slope_vectors.insert(_point_index,(left_slope,(1,_current_slope)))
      
                
        cb  =  CallbackFunction(self.startMoveSlopePoint, sid, self.lid,self.tid)
        self.canvas.tag_bind(sid, '<ButtonPress-1>', cb)
        cb  =  CallbackFunction( self.endMoveSlopePoint, sid,self.lid,self.tid)
        self.canvas.tag_bind(sid, '<ButtonRelease-1>', cb)
        #finding slope and calling spline with new slope
        self.spline.setControlPoints(cpts,slope_vectors) 
        
        
        if self._line_On == 1:
            self.removeFunctionLine()
            self.drawFunctionLine()
        #self.line_coords = c.coords(self.lid)
        self._var_ = 0
        self.canvas.tag_raise('pt')
        self.sid = sid
        

    def endMoveSlopePoint(self, sid,lid, tid,index):
        c  =  self.canvas
        c.tag_unbind(sid, '<Any-Motion>')    
        self.appendToHistory()
        
        
    def createPoint(self,id,event = None):
        """function to insert point on spline""" 
        c = self.canvas
        old_cids = c.find_withtag("ctrlPoints_"+self.uniqueId)
        coords = []
        for cid in old_cids:
            if "tied" in c.gettags(cid):
                coords.append([round(c.coords(cid)[0]),round(c.coords(cid)[1]),round(c.coords(cid)[2]),round(c.coords(cid)[3])])
        x = event.x
        y = event.y
        slope = ((1,0),(1,0))
        vx1,vy1 = slope[0]
        slope_length=(30,30)
        pts = self.getControlPoints()
        pts.append((x,y))
        pts.sort()
        ind = pts.index((x,y))
        pos = ind
        self.slopes_lengths.insert(pos,slope_length)
        s_vectors = copy.deepcopy(self.getSlopesVectors())
        s_vectors.insert(pos,slope)
        self.spline.setControlPoints(pts,s_vectors)
        self.drawSpline()
        self.appendToHistory()
        cids = c.find_withtag("ctrlPoints_"+self.uniqueId)    
        
        for cid in cids:
                if  [round(c.coords(cid)[0]),round(c.coords(cid)[1]),round(c.coords(cid)[2]),round(c.coords(cid)[3])] in coords:
                    tags = c.gettags(cid)
                    tags += ("tied",)
                    if "tied" not in c.gettags(cid):
                        c.itemconfig(cid,tags = tags)
                    lsid,rsid = self.find_sids(cid)
                    if lsid!=None and rsid!=None:
                        
                        c.itemconfigure(lsid, fill = 'blue')
                        c.itemconfigure(rsid, fill = 'blue')
                    
    def deletePoint(self,id,event = None):
        """function to delete point on spline"""
        c = self.canvas
        old_cids = c.find_withtag("ctrlPoints_"+self.uniqueId)
        coords = []
        for cid in old_cids:
            if "tied" in c.gettags(cid):
                coords.append([round(c.coords(cid)[0]),round(c.coords(cid)[1]),round(c.coords(cid)[2]),round(c.coords(cid)[3])])
        p=None
        pts = self.getControlPoints()
        slopes=copy.deepcopy(self.getSlopesVectors())
        cur_point = c.coords(id)
        pts = map(lambda p: ((int(round(p[0])),int(round(p[1])))) ,pts)
        pn = (round(cur_point[0]+2),round(cur_point[1]+2))
        if pn in pts:
                p = pn
        
        if p:
            ind = pts.index(p)
            pts.remove(p)
            slope=slopes[ind]
            slopes.remove(slope)
            self.spline.setControlPoints(pts,slopes)
            self.slopes_lengths.remove(self.slopes_lengths[ind])
            self.drawSpline()
            self.appendToHistory()
            cids = c.find_withtag("ctrlPoints_"+self.uniqueId)    
            for cid in cids:
                if  [round(c.coords(cid)[0]),round(c.coords(cid)[1]),round(c.coords(cid)[2]),round(c.coords(cid)[3])] in coords:
                    tags = c.gettags(cid)
                    tags += ("tied",)
                    if "tied" not in c.gettags(cid):
                        c.itemconfig(cid,tags = tags)
                    lsid,rsid = self.find_sids(cid)
                    if lsid!=None and rsid!=None:
                        c.itemconfigure(lsid, fill = 'blue')
                        c.itemconfigure(rsid, fill = 'blue')
        else:
            return
    
    
    def untieNormals(self,cid,event=None):
        c=self.canvas
        if "tied" in c.gettags(cid):
            c.dtag(cid,"tied")
            ls,rs = self.find_sids(cid)
            if ls!=None and rs!=None:
                c.itemconfigure(ls, fill = 'green')
                c.itemconfigure(rs, fill = 'green')
            cb= CallbackFunction(self.tieNormals, cid)
            c.tag_bind(cid, '<ButtonPress-3>', cb)
            
    def tieNormals(self,cid,event=None):
         
        c = self.canvas
        cid_coords = c.coords(cid)
        cpts= self.getControlPoints()
        rounded_cpts=[]
        uid = self.uniqueId
        tag1 =  'ctrlPoints_'+uid
        #if cid is first or last cpts cannot be tied
        cs = c.find_withtag(tag1)
        if cid ==cs[0] or cid==cs[-1]:
           print "cannot tie first and last coords of spline"
           return
        for p in cpts:
            rounded_cpts.append((round(p[0]),round(p[1])))
        if (round(cid_coords[0]+2),round(cid_coords[1]+2)) in  rounded_cpts :
            ind = rounded_cpts.index((round(cid_coords[0]+2),round(cid_coords[1]+2)))
            slopes_vectors=self.getSlopesVectors()
            slopes_vectors.remove(slopes_vectors[ind])
            slopes_vectors.insert(ind,((1,0),(1,0)))
            cur_len = 21.6
            #if  self.slopes_lengths[ind][0]>self.slopes_lengths[ind][1]:
            #    cur_len = self.slopes_lengths[ind][0]
            #else:
            #    cur_len = self.slopes_lengths[ind][1]
            self.slopes_lengths.remove(self.slopes_lengths[ind])
            self.slopes_lengths.insert(ind,(cur_len,cur_len))
            self.spline.setControlPoints(cpts,slopes_vectors)
            if self._slopes_0n == 1:
                self.stopEditingCurve()
                self.startEditingCurve()
                cids = c.find_withtag(tag1)
                for ci in cids:  
                    if "tied" in c.gettags(ci):
                        ls,rs = self.find_sids(ci)
                        if ls!=None and rs!=None:
                            c.itemconfigure(ls, fill = 'blue')
                            c.itemconfigure(rs, fill = 'blue')    
            if self._line_On == 1:
                self.removeFunctionLine()
                self.drawFunctionLine()
            
            ls,rs = self.find_sids(cid)
            if ls!=None and rs!=None:
                c.itemconfigure(ls, fill = 'blue')
                c.itemconfigure(rs, fill = 'blue')
            tags = c.gettags(cid)
            tags += ("tied",)
            if "tied" not in c.gettags(cid):
                c.itemconfig(cid,tags = tags)
                
            cb= CallbackFunction(self.untieNormals, cid)
            c.tag_bind(cid, '<ButtonPress-3>', cb)        
        
    def getControlPoints(self):
        sx, sy  =  self.scale
        ox  =  self.originX
        oy  =  self.originY
        pts =self.spline.getControlPoints()       
        cpts = []
        for (x, y) in pts:
         
         rx = x
         ry = y
         if y <= 1.5:   
            
            rx  =  ox + sx*x
            ry  =  oy - sy*y
         cpts.append((rx,ry))
        self.controlpoints = cpts
        return cpts

    
    def find_sids(self,cid):
            
            #finding sids at cid
            c=self.canvas
            lids = c.find_withtag("ln")
            pids = c.find_withtag("pt")
            rsids = c.find_withtag("right")
            lsids = c.find_withtag("left")
            r_lids = []
            l_lids = []
            r_sids = []
            l_sids = []
            cid_coords = c.coords(cid)
            current_rsid =None
            current_lsid = None
            
            for  s in pids:
                if s in rsids:
                    r_sids.append(s)
                if s in lsids:
                    l_sids.append(s)
            
            for  s in lids:
                if s in rsids:
                    r_lids.append(s)
                if s in lsids:
                    l_lids.append(s)
            for r in r_lids:
                if (round(c.coords(r)[0]),round(c.coords(r)[1]))==(round(cid_coords[0]+2),round(cid_coords[1]+2)):
                    current_rlid = r
                    for i in r_sids:
                        if (round(c.coords(r)[2]),round(c.coords(r)[3]))==(round(c.coords(i)[0]+2),round(c.coords(i)[1]+2)):
                            current_rsid = i
            for l in l_lids:
                if (round(c.coords(l)[2]),round(c.coords(l)[3]))==(round(cid_coords[0]+2),round(cid_coords[1]+2)):
                    current_llid = l
                    for i in l_sids:
                        if (round(c.coords(l)[0]),round(c.coords(l)[1]))==(round(c.coords(i)[0]+2),round(c.coords(i)[1]+2)):
                            current_lsid = i
                            
            return current_lsid,current_rsid

    
    def find_othersid(self,sid):
        c=self.canvas
        lsids = c.find_withtag("left")
        rsids = c.find_withtag("right")
        lids = c.find_withtag("ln")
        pids = c.find_withtag("pt")
        cids = c.find_withtag("ctrlPoints_"+self.uniqueId)
        cids_coords = []
        for cid in cids:
            cids_coords.append(c.coords(cid))
        r_lids = []
        l_lids = []
        r_sids = []
        l_sids = []    
        current_id = None
        current_cid = None
        for  s in pids:
                if s in rsids:
                    r_sids.append(s)
                if s in lsids:
                    l_sids.append(s)
            
        
        for  s in lids:
                if s in rsids:
                    r_lids.append(s)
                if s in lsids:
                    l_lids.append(s)
        if "left" in c.gettags(sid):
            for l in l_lids:
                
                #llid=sid
                if (round(c.coords(l)[0])-2,round(c.coords(l)[1]-2)) == (round(c.coords(sid)[0]),round(c.coords(sid)[1])) or (round(c.coords(l)[0]),round(c.coords(l)[1])) == (round(c.coords(sid)[0]),round(c.coords(sid)[1])) :
                    current_lid = l
                    for cid in  cids_coords:
                        #llid=cid
                        if (round(c.coords(l)[2]),round(c.coords(l)[3])) == (round(cid[2]-2),round(cid[3]-2)):
                            ind = cids_coords.index(cid)
                            current_cid = cids[ind]
                            for r in r_lids:
                                #cid =rlid
                                if (round(cid[0]),round(cid[1]))==(round(c.coords(r)[0]-2),round(c.coords(r)[1]-2)):
                                    for p in r_sids:
                                        #rlid=rsid 
                                        if (round(c.coords(r)[2]+2),round(c.coords(r)[3]+2)) == (round(c.coords(p)[2]),round(c.coords(p)[3])):
                                            current_id = p
        if "right" in  c.gettags(sid):
            for r in r_lids:
                #rlid=rsid
                if (round(c.coords(r)[2]+2),round(c.coords(r)[3]+2)) == (round(c.coords(sid)[2]),round(c.coords(sid)[3])) or (round(c.coords(r)[2]),round(c.coords(r)[3])) == (round(c.coords(sid)[2]),round(c.coords(sid)[3])): 
                    for cid in  cids_coords:
                        #cid =rlid
                        if (round(cid[0]),round(cid[1]))==(round(c.coords(r)[0]-2),round(c.coords(r)[1]-2)):
                            ind = cids_coords.index(cid)
                            current_cid = cids[ind]
                            for l in l_lids:
                                #llid=cid
                                if (round(c.coords(l)[2]),round(c.coords(l)[3])) == (round(cid[2]-2),round(cid[3]-2)):
                                    for p in l_sids:
                                        #llid=sid
                                        if (round(c.coords(l)[0])-2,round(c.coords(l)[1]-2)) == (round(c.coords(p)[0]),round(c.coords(p)[1])):
                                            current_id = p
                    
            
              
        return current_id,current_cid

    def drawSpline(self):
        self.removeControlPoints()
        self.drawControlPoints()
        self._drawSpline()

    def _drawSpline(self):
        if self._slopes_0n == 1:
                self.stopEditingCurve()
                self.startEditingCurve()
        if self._line_On == 1:
                self.removeFunctionLine()
                self.drawFunctionLine()
                
                
    def invertSpline(self):
        """This function is for inverting graph by reverse computing controlpoints"""
        
        invert_points = []
        cpts = self.getControlPoints()
        slopes = self.getSlopesVectors()
        for p in cpts:
            y=self.originY -(p[1]-50)
            invert_points.append((p[0],y))
        self.spline.setControlPoints(invert_points,slopes)
        self.drawSpline()    
        self.appendToHistory()
         
        
    def stepBack(self):
       """when stepBack button clicked previous step is displayed.History of
        all the steps done is remembered and when stepback clicked from history
        list previous step is shown and that step is removed from history list """ 
        
       if len(self.history_cpts) == 1: 
            cpts = self.history_cpts[-1]    
            slopes = self.history_spline_slopes[-1]
            slopes_lengths = self.history_spline_slopes_lengths[-1]
       if len(self.history_cpts)>1:     
            self.history_cpts = self.history_cpts[:-1]
            self.history_spline_slopes = self.history_spline_slopes[:-1]
            self.history_spline_slopes_lengths = self.history_spline_slopes_lengths[:-1]
            cpts = self.history_cpts[-1]    
            slopes = self.history_spline_slopes[-1]
            slopes_lengths = self.history_spline_slopes_lengths[-1]
       self.slopes_lengths = [] 
       for j in  slopes_lengths :
           self.slopes_lengths.append(j)
        
       self.spline.setControlPoints( cpts, slopes)
       self.drawSpline()
             
    def appendToHistory(self):
        """function to append  to history"""
        self.history_cpts.append(copy.deepcopy(self.controlpoints))        
        self.history_spline_slopes_lengths.append(copy.deepcopy(self.slopes_lengths))
        self.history_spline_slopes.append(copy.deepcopy(self.slopesvectors))

    def getSlopesVectors(self):
        sls = self.spline.getSlopesVectors()
        self.slopesvectors = sls
        return sls

    def setControlPoints(self,pn,sl):
        self.spline.setControlPoints(pn,sl)

if __name__ == "__main__":

    from spline import Spline
    #sp  =  Spline( [0, .5, 1], [0, .7, 1], [ (None, (1,0)), ((1,1),(1,0)), ((1,0), None) ] )
    #sp = Spline([(0,0),(0.3,0.5),(.5,.7),(1,1)],[ (None, (1,-1)), ((1,-3),(1,-4)),((1,-2),(1,-1)), ((1,-3), None) ] )
    sp = Spline([(0.0,.5),(1.0,0.3)],[ (None, (1,-1)), ((1,-1), None) ] )
    #sp = Spline([(100,600),(200,400),(400,200),(600,100)],[ (None, (1,-4)),((1,-3),(1,-4)),((1,-2),(1,-1)),((2,3), None) ])
    canvas  =  Tkinter.Canvas(width = 700, height = 700)
    canvas.pack()
    #canvas.configure(cursor="cross")
    spg  =  splineGUI(sp, canvas, (100, 600, 500, 500), 'abc')
    #canvas.create_rectangle([100,600,600,100])
    spg.drawSlopes()
    spg.drawFunctionLine()









