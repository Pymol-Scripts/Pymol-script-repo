###########################################################
#
#  Pymol script copyright Matthew O'Meara and Xavier Ambroggio 2007
#
#  Last updated Nov 29, 2007
#
#  Draw an axis given a point and a direction.  Optionally give color,
#  length and width.
#
#  Usage: from the pymol menu click file->run...  then find this file.
#  Then at the prompt, type
#
#           draw_axis x,y,z, i,k,j
#
#  where (x,y,z) is the point and (i,k,j) is the direction
#
#
#  Also one can run the script as follows
#
#
#           draw_axis x,y,z, i,k,j, length r,g,b, width,
#
#  where (r,g,b) is the color for the line in red, green, colors,
#  width is the thickness of the line and length is the length of the
#  line.
#
#
#  For a fun example of the script run from the command line after the
#  script is loaded
#
#           draw_axis_example
#
#
 
from pymol.cgo import *    # get constants
from pymol import cmd
 
import math
 
class Counter:
   def __init__(self):
       self.state = 1
counter = Counter()
 
def draw_axis(x=None, y=None, z=None, i=None, j=None, k=None, length=20.0, r=1.0, g=1.0, b=1.0, width=1.0 ):
   if x == None or y == None or z == None or i == None or j == None or k== None :
       print 'Usage: draw_axis x,y,z, i,k,j, length, r,g,b, width'
       print 'draw a line centered at (x,y,z) with the direction vector (i,j,k)'
       print 'length, color (r,g,b), and width arguments are optional'
#        print 'For a fun example of the command, run draw_axis_example'
   else :
       x,y,z = float(x), float(y), float(z)
       i,j,k = float(i), float(j), float(k)
       r,g,b = float(r), float(g), float(b)
       width = float(width)
       length = float(length) / 2.0
 
       x1,y1,z1 = (x+i*length,y+j*length,z+k*length)
       x2,y2,z2 = (x-i*length,y-j*length,z-k*length)
 
       obj = [
           LINEWIDTH, width,
           BEGIN, LINES,
 
           COLOR,   r,  g,  b,
           VERTEX, x1, y1, z1,
           VERTEX, x2, y2, z2,
 
           END
           ]
 
       cmd.load_cgo(obj,'axis'+str(counter.state))
       counter.state += 1
 
cmd.extend("draw_axis", draw_axis)
 
# a simple example
#draw_line(x=18.232,  y=17.150,  z=9.488,
#          i=-.226639,j=0.708772,k=-.668039,
#          r=1,       b=1,       g=1,
#          width=1,   length=1)
 
 
 
 
# a more complex example
 
#import random
#def example1(n, f):
#    """draw a gradient field with n segments with the function f(x,y,z)=(i,j,k)"""
#    for i in range(n):
#        scale = 4
#        x,y,z = [random.random()*scale for i in range(3)]
#        i,j,k = f(x,y,z)
 
#        draw_axis(x,y,z,i,j,k,abs(i),abs(j),abs(k))
 
 
#def f(x,y,z):
#    return (2*x,pow(z,2)+x,y-z)
 
#cmd.extend("draw_axis_example", lambda :example1(1000,f))
