##############################################################################
#
#
#   Author: Sowjanya Karnati
#
#
###############################################################################
#
#
#
#
#$Id: test_graphtool.py,v 1.8 2007/09/27 21:00:50 annao Exp $
import sys,unittest
from Tkinter import Tk
from time import sleep
from mglutil.gui.BasicWidgets.Tk.graphtool import GraphApp
import numpy
import os
def pause():
    import time
    time.sleep(0.1)

class GrapToolTest(unittest.TestCase):
    
    def test_default_graphtool(self):
        """create and withdarw"""
        root=Tk()
        app=GraphApp(root)
        app.master.update()
        pause()
        app.master.destroy()


    def test_graphtool_default_curve(self):
        """testing default curve"""
        root=Tk()
        app=GraphApp(root)
        app.master.update()
        app.defaultcurve_cb()
        app.master.update()
        pause()
        self.assertEqual(app.getControlPoints()[1:-1],app.default_points)
        app.master.destroy()
    
    def test_graphtool_reset(self):
        """testing resetting graph"""
        root=Tk()
        app=GraphApp(root)
        app.master.update()
        app.resetAll_cb()
        app.master.update()
        self.assertEqual(app.getControlPoints(),[(50, 275), (305, 20)])
        self.assertEqual(len(app.curovals),2)
        self.assertEqual(app.history,[])
        self.assertEqual(app.oldpoints,[(50, 275), (305, 20)])
        pause()
        app.master.destroy()
    
    def test_graphtool_setControlPoints(self):
        """tests setControlPoints"""
        root=Tk()
        app=GraphApp(root)
        app.master.update()
        points=[(100,150),(150,200)]
        app.setControlPoints(points)
        app.master.update()
        pause()
        self.assertEqual(app.getControlPoints()[1:-1],points)
        self.assertEqual(len(app.curovals),len(points))
        pause()
        app.master.destroy()

    def test_stepback_cb(self):
        """testing setback function"""
        root=Tk()
        app=GraphApp(root)
        app.master.update()
        points=[(100,150),(150,200)]
        app.setControlPoints(points)
        app.master.update()
        npoints=[(200,200)]
        app.setControlPoints(npoints)
        app.master.update()
        app.stepBack_cb()
        app.master.update()
        self.assertEqual(app.getControlPoints()[1:-1],[(100,150),(150,200)])
        app.master.destroy()

    def test_graphtool_setSensitivity(self):
        """testing setSensitivity function"""
        root=Tk()
        app=GraphApp(root)
        app.master.update()
        app.setSensitivity(0.06)
        app.master.update()
        self.assertEqual(app.d1scalewheel.get(),0.06)
        app.master.destroy()

    def test_graphtool_MouseButtonUp(self):
        """testing when mousebuttonup var set to 1"""
        root=Tk()
        app=GraphApp(root)
        app.master.update()
        pause()
        val1=app.caluculate_ramp()
        app.mousebuttonup=1
        points=[(100,150),(150,200)]
        app.setControlPoints(points)
        app.master.update()
        val2=app.caluculate_ramp()
        self.assertFalse(numpy.alltrue(val1 == val2) )
        #self.assertEqual(val1!=val2,True)
        app.master.destroy()
        
    def test_graphtool_Continuous(self):
        """testing when continuous var set to 1"""
        root=Tk()
        app=GraphApp(root)
        app.master.update()
        val1=app.caluculate_ramp()
        app.continuous=1
        points=[(100,150),(150,200)]
        app.setControlPoints(points)
        app.master.update()
        val2=app.caluculate_ramp()
        self.assertFalse(numpy.alltrue(val1 == val2) )
        app.master.destroy()
        

    def test_graphtool_Update(self):
        """testing when update var set to 1"""
        root=Tk()
        app=GraphApp(root)
        app.master.update()
        points=[(130,150)]
        app.setControlPoints(points)
        val1=app.caluculate_ramp()
        app.update=1
        points=[(100,150),(150,200)]
        app.setControlPoints(points)
        app.Update()
        app.master.update()        
        val2=app.caluculate_ramp()
        #self.assertEqual(val1!=val2,True)
        self.assertFalse(numpy.alltrue(val1 == val2) )
        app.master.destroy()

    def test_graphtool_read_write(self):
        """testing read ,write functions"""    
        root=Tk()
        app=GraphApp(root)
        app.master.update()
        app.setControlPoints(app.default_points)
        app.master.update()
        app.write("new1_Graph.py")
        app.master.update()
        self.assertEqual(os.path.exists("new1_Graph.py"),True)
        fptr=open("new1_Graph.py")
        data=fptr.readlines()
        points=[(38, 37), (51, 125), (104, 197), (25, 4)]
             
        self.assertEqual(eval(data[0][:-1])[1:-1],app.default_points)
        app.read("new1_Graph.py")
        app.master.update()
        self.assertEqual(app.getControlPoints()[1:-1],app.default_points)
        os.system("rm -rf new1_Graph.py")
        app.master.destroy()
    
    ####TESTING ILLEGAL INPUT

    def test_graphtool_setControlpoints_illegal_input(self):
        """testing illegal input  for setControlpoints"""
        root=Tk()
        app=GraphApp(root)
        app.master.update()
        self.assertRaises(AssertionError,app.setControlPoints,"abc")
        app.master.destroy()
        
    def test_graphtool_setControlpoints_illegal_input_1(self):
        """testing illegal input coordinates out of range raises assertion error"""
        root=Tk()
        app=GraphApp(root)
        app.master.update()
        self.assertRaises(AssertionError,app.setControlPoints,[(0,200)])
        app.master.destroy()
    
if __name__ == '__main__':
   unittest.main()

