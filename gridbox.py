# Copyright (c) 2014 Chen Zhaoqiang DDDC SIMM
#qq 744891290 qqgroup 144539789
# A PyMOL script for drawing a CGO box by the coordinate of the center,and the size of box
#center_x,center_y,center_z,size_x,size_y,size_z,



import pymol
from pymol import cmd
from pymol.wizard import Wizard
from chempy import cpv
from pymol.cgo import *

def gridbox(center_x,center_y,center_z,size_x,size_y,size_z,name="gridbox",r1=0,g1=0,b1=1,trasp=0.2):
	"""
	DESCRIPTION
	Create a box from the center coordinate of the box and the size of box
	
	USAGE
	run gridbox.py
	1the simplest
	gridbox center_x,center_y,center_z,size_x,size_y,size_z
	2rename the box object
	gridbox center_x,center_y,center_z,size_x,size_y,size_z,name,
	3set the color of the box object
	gridbox center_x,center_y,center_z,size_x,size_y,size_z,name,r1,g1,b1
	4set the trasp of the box
	gridbox center_x,center_y,center_z,size_x,size_y,size_z,name,r1,g1,b1,trasp
	
	ps:the value of r1,g1,b1 trasp   range  is 0-1
	   trasp=1,no trasprent
	
	
	"""
	center_x=float(center_x)
	center_y=float(center_y)
	center_z=float(center_z)
	size_x=float(size_x)
	size_y=float(size_y)
	size_z=float(size_z)
	r1=float(r1)
	g1=float(g1)
	b1=float(b1)
	trasp=float(trasp)
	
	p1=[center_x-size_x/2,center_y-size_y/2,center_z+size_z/2]
	p2=[center_x-size_x/2,center_y+size_y/2,center_z+size_z/2]
	p3=[center_x+size_x/2,center_y+size_y/2,center_z+size_z/2]
	p4=[center_x+size_x/2,center_y-size_y/2,center_z+size_z/2]
	p5=[center_x-size_x/2,center_y-size_y/2,center_z-size_z/2]
	p6=[center_x-size_x/2,center_y+size_y/2,center_z-size_z/2]
	p7=[center_x+size_x/2,center_y+size_y/2,center_z-size_z/2]
	p8=[center_x+size_x/2,center_y-size_y/2,center_z-size_z/2]
	obj=[
			 ALPHA,trasp,
			 COLOR,r1,g1,b1,
			 BEGIN, TRIANGLE_STRIP,
			 VERTEX,p1[0],p1[1],p1[2],
			 VERTEX,p2[0],p2[1],p2[2],
			 VERTEX,p4[0],p4[1],p4[2],
			 VERTEX,p3[0],p3[1],p3[2],
			 END,
			 
			 BEGIN, TRIANGLE_STRIP,
#			 COLOR,1,0,0,
			 VERTEX,p1[0], p1[1], p1[2],
			 VERTEX,p5[0], p5[1], p5[2],
			 VERTEX,p4[0], p4[1], p4[2],
			 VERTEX,p8[0], p8[1], p8[2],
			 END,
			 
			 BEGIN, TRIANGLE_STRIP,
			 VERTEX,p4[0],p4[1],p4[2],
			 VERTEX,p8[0],p8[1],p8[2],
			 VERTEX,p3[0],p3[1],p3[2],
			 VERTEX,p7[0],p7[1],p7[2],
			 END,
			 
			 
			 BEGIN, TRIANGLE_STRIP,
			 VERTEX,p7[0],p7[1],p7[2],
			 VERTEX,p3[0],p3[1],p3[2],
			 VERTEX,p6[0],p6[1],p6[2],
			 VERTEX,p2[0],p2[1],p2[2],
			 END,
			 
			 BEGIN, TRIANGLE_STRIP,
#			 COLOR,0,1,0,
			 VERTEX,p2[0],p2[1],p2[2],
			 VERTEX,p6[0],p6[1],p6[2],
			 VERTEX,p1[0],p1[1],p1[2],
			 VERTEX,p5[0],p5[1],p5[2],
			 END,
			 
			 BEGIN, TRIANGLE_STRIP,
			 VERTEX,p1[0],p1[1],p1[2],
			 VERTEX,p5[0],p5[1],p5[2],
			 VERTEX,p4[0],p4[1],p4[2],
			 VERTEX,p8[0],p8[1],p8[2],
			 END,
			 
			 BEGIN, TRIANGLE_STRIP,
			 VERTEX,p5[0],p5[1],p5[2],
			 VERTEX,p8[0],p8[1],p8[2],
			 VERTEX,p6[0],p6[1],p6[2],
			 VERTEX,p7[0],p7[1],p7[2],
			 END,
			 ]
	cmd.load_cgo(obj, name)
	


cmd.extend('gridbox', gridbox)