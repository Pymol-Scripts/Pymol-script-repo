# axes.py
from pymol.cgo import *
from pymol import cmd
from pymol.vfont import plain
 
# create the axes object, draw axes with cylinders coloured red, green,
#blue for X, Y and Z
 
obj = [
   CYLINDER, 0., 0., 0., 10., 0., 0., 0.2, 1.0, 1.0, 1.0, 1.0, 0.0, 0.,
   CYLINDER, 0., 0., 0., 0., 10., 0., 0.2, 1.0, 1.0, 1.0, 0., 1.0, 0.,
   CYLINDER, 0., 0., 0., 0., 0., 10., 0.2, 1.0, 1.0, 1.0, 0., 0.0, 1.0,
   ]
 
# add labels to axes object (requires pymol version 0.8 or greater, I
# believe
 
cyl_text(obj,plain,[-5.,-5.,-1],'Origin',0.20,axes=[[3,0,0],[0,3,0],[0,0,3]])
cyl_text(obj,plain,[10.,0.,0.],'X',0.20,axes=[[3,0,0],[0,3,0],[0,0,3]])
cyl_text(obj,plain,[0.,10.,0.],'Y',0.20,axes=[[3,0,0],[0,3,0],[0,0,3]])
cyl_text(obj,plain,[0.,0.,10.],'Z',0.20,axes=[[3,0,0],[0,3,0],[0,0,3]])
 
# then we load it into PyMOL
cmd.load_cgo(obj,'axes')
