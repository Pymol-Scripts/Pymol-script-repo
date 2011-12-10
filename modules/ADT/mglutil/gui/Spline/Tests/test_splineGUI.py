#
#Author : Sowjanya Karnati
#
#$Id: test_splineGUI.py,v 1.1 2007/07/13 18:22:24 sowjanya Exp $


from mglutil.gui.Spline.spline import Spline
from mglutil.gui.Spline.splineGUI import splineGUI
import unittest
import Tkinter

class splineGUIBaseTest(unittest.TestCase):

    def test_defualt_args(self):
        """calling splineGUI with default args"""
        sp = Spline([(0.0,.5),(1.0,0.3)],[ (None, (1,-1)), ((1,-1), None) ] )
        canvas  =  Tkinter.Canvas(width = 700, height = 700)
        canvas.pack()
        spg  =  splineGUI(sp, canvas, (100, 600, 500, 500), 'abc')
        self.assertEqual(spg.getControlPoints(),[(100.0, 372.72727272727275), (600.0, 463.63636363636363)])


    def test_setControlPoints(self):
        """tests settting control points"""
        sp = Spline([(0.0,.5),(1.0,0.3)],[ (None, (1,-1)), ((1,-1), None) ] )
        canvas  =  Tkinter.Canvas(width = 700, height = 700)
        canvas.pack()
        spg  =  splineGUI(sp, canvas, (100, 600, 500, 500), 'abc')
        points = [(100,300),(600,300)]
        slopes = [(None,(1,-1)),((1,3),None)]
        spg.setControlPoints(points,slopes)    
        spg.drawSpline()
        self.assertEqual(spg.getControlPoints(),points)


    def test_invertSpline(self):
        """tests calling inverspline"""
        sp = Spline([(0.0,.5),(1.0,0.3)],[ (None, (1,-1)), ((1,-1), None) ] )
        canvas  =  Tkinter.Canvas(width = 700, height = 700)
        canvas.pack()
        spg  =  splineGUI(sp, canvas, (100, 600, 500, 500), 'abc')  
        points = spg.getControlPoints()
        spg.invertSpline()
        self.assertEqual(spg.getControlPoints()!=points,True)


    def test_stepBack(self):
        """tests stepBack function"""
        sp = Spline([(0.0,.5),(1.0,0.3)],[ (None, (1,-1)), ((1,-1), None) ] )
        canvas  =  Tkinter.Canvas(width = 700, height = 700)
        canvas.pack()
        spg  =  splineGUI(sp, canvas, (100, 600, 500, 500), 'abc') 
        first_points = spg.getControlPoints()
        points = [(100,300),(600,300)]
        slopes = [(None,(1,-1)),((1,3),None)]
        spg.setControlPoints(points,slopes)    
        spg.drawSpline()
        spg.stepBack()
        self.assertEqual(spg.getControlPoints()!=points,True)
        self.assertEqual(spg.getControlPoints(),first_points)

if __name__  ==  '__main__':
   unittest.main()
