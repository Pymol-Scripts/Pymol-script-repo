#########################################################################
#
# Date: Jun 2002 Authors: Daniel Stoffler, Michel Sanner
#
#    stoffler@scripps.edu
#    sanner@scripps.edu
#
# Copyright:  Daniel Stoffler, Michel Sanner, and TSRI
#
#########################################################################

# $Id: test_thumbwheel.py,v 1.12 2007/12/04 21:28:04 vareille Exp $

import sys,unittest,Tkinter
from time import sleep
tw = None
wasCalled = False
mySum = 0.
from mglutil.gui.BasicWidgets.Tk.thumbwheel import ThumbWheel
from mglutil.util.misc import ensureFontCase

class Dummyevent:
    
    def __init__(self,y = 'None',x = 'None'):
        self.y = y
        self.x = x
        
def pause():
    import time
    time.sleep(0.1)
class ThumbWheelBaseTest(unittest.TestCase,Dummyevent):

    def test_constructorOptions(self):
        # test all possible constructor options
        root = Tkinter.Tk()
        tw = ThumbWheel(width=100, height=26,wheelPad=4,master = root,
                        labcfg={'fg':'black', 'side':'left', 'text':'Test:'},
                        wheelLabcfg1={'font':(ensureFontCase('times'),14,'bold')},
                        wheelLabcfg2={'font':(ensureFontCase('times'),14,'bold')},
                        canvascfg={'bg':'blue'}) 
        tw.master.update()
        pause()
        tw.master.destroy()


    def test_setValue(self):
        # test setting of a value
        root = Tkinter.Tk()
        tw = ThumbWheel(width=100, height=26, value=1.0,master = root)
        tw.set(10)
        tw.master.update()
        #tw.master.update()
        pause()
        self.assertEqual(tw.value,10.0)
        tw.configure(type=float)
        tw.set(20.0)
        tw.master.update()
        #tw.master.update()
        pause()
        self.assertEqual(tw.value,20.0)
        tw.master.destroy()


    def test_callback_1(self):
        global wasCalled
        def foo(val):
            global wasCalled
            wasCalled=True
        root = Tkinter.Tk()
        tw = ThumbWheel(width=100, height=26,master = root)
        tw.callbacks.AddCallback(foo)
        # setValue(val) should NOT call callback
        wasCalled = 0
        tw.setValue(5.0)
        #tw.master.update()
        if wasCalled:
            raise RuntimeError
        # set(val) should call callback
        tw.set(6.0)
        if not wasCalled:
            raise RuntimeError
        wasCalled = False
        tw.set(7.0, update=0)
        if wasCalled:
            raise RuntimeError
        tw.master.update()
        pause()
        tw.master.destroy()


    def test_callback_2(self):
        global wasCalled
        def foo(val):
            global wasCalled
            wasCalled=True
        root = Tkinter.Tk()
        tw = ThumbWheel(width=100, height=26, master=root, callback=foo)
        self.assertEqual(tw.callbacks.callbacks, [foo,], "Expecting to have foo added to the callbackmanager callbacks list")
        #tw.callbacks.AddCallback(foo)
        # setValue(val) should NOT call callback
        wasCalled = 0
        tw.setValue(5.0)
        #tw.master.update()
        if wasCalled:
            raise RuntimeError
        # set(val) should call callback
        tw.set(6.0)
        if not wasCalled:
            raise RuntimeError
        wasCalled = False
        tw.set(7.0, update=0)
        if wasCalled:
            raise RuntimeError
        tw.master.update()
        pause()
        tw.master.destroy()

    def test_callback_3(self):
        global wasCalled1
        global wasCalled2
        
        def foo1(val):
            global wasCalled1
            wasCalled1=True

        def foo2(val):
            global wasCalled2
            wasCalled2=True
            
        root = Tkinter.Tk()
        tw = ThumbWheel(width=100, height=26, master=root, callback=[foo1, foo2])
        self.assertEqual(tw.callbacks.callbacks, [foo1,foo2], "Expecting to have foo added to the callbackmanager callbacks list")
        # setValue(val) should NOT call callback
        wasCalled1 = False
        wasCalled2 = False
        tw.setValue(5.0)
        if wasCalled1 and wasCalled2:
            raise RuntimeError

        # set(val) should call callback
        tw.set(6.0)
        if not wasCalled1 and not wasCalled2:
            raise RuntimeError

        wasCalled1 = False
        wasCalled2 = False
        tw.set(7.0, update=0)
        if wasCalled1 and wasCalled2:
            raise RuntimeError
        tw.master.update()
        pause()
        tw.master.destroy()


    def test_setType(self):
        # test setting of type
        tw = ThumbWheel(width=100, height=26, type='float')
        tw.set(100.0)
        # we can mix configure(type=float) and configure(type='float')
        tw.configure(type=float)
        self.assertEqual(tw.type,float)
        self.assertEqual(type(tw.get()),type(1.0))
        tw.configure(type='int')
        tw.master.update()
        pause()
        self.assertEqual(tw.type,int)
        self.assertEqual(type(tw.get()),type(1))
        tw.master.update()
        pause()
        tw.master.destroy()
    

    def test_setMin(self):
        global wasCalled
        # test that we do not call the callback when we set the minimum value
        def foo(val):
            global wasCalled
#           print "I should not be called"
            wasCalled=1

        from mglutil.gui.BasicWidgets.Tk.thumbwheel import ThumbWheel
    
        tw = ThumbWheel(width=100, height=26)
        tw.callbacks.AddCallback(foo)
        wasCalled=0
        # set the value without calling callback
        tw.setValue(4)
        # change minimum without affecting value ==> nocallback
        tw.configure(min=2)
        tw.master.update()
        pause()
        if wasCalled:
            raise RuntimeError
        # change minimum with affecting value ==> callback
        tw.configure(min=6)
        if wasCalled==0:
            raise RuntimeError
        tw.master.update()
        pause()
        tw.master.destroy()

    
    def test_setMax(self):
        global wasCalled
        # test that we do not call the callback when we set the maximum value
        def foo(val):
            global wasCalled
#           print "I should not be called"
            wasCalled=1
        tw = ThumbWheel(width=100, height=26)
        tw.callbacks.AddCallback(foo)
        tw.master.update()
        pause()
        wasCalled=0
        # set the value without calling callback
        tw.setValue(4)
        # change maximum without affecting value ==> nocallback
        tw.configure(max=6)
        if wasCalled:
            raise RuntimeError
        # change maximum with affecting value ==> callback
        tw.configure(max=2)
        if wasCalled==0:
            raise RuntimeError
        tw.master.update()
        pause()
        tw.master.destroy()


    def test_setIncrement(self):
        #THIS IS A BUG
        #it is not working correctly
        tw = ThumbWheel(width=100, height=26)
        tw.set(5.0)
        tw.master.update()
        pause()
        tw.configure(increment=10)
        tw.set(24.0)
        self.assertEqual(tw.increment,10.0)
        #24 is shown instead of 25
        self.assertEqual(tw.get(),24.0)
        tw.master.update()
        pause()
        tw.master.destroy()

    
    def test_setPrecision(self):
        
        tw = ThumbWheel(width=100, height=26)
        tw.configure(type='float')
        tw.configure(precision=5)
        # increase precision to 5 - this is only visual. value is not changed
        self.assertEqual(tw.precision,5)
        tw.set(4.123456789)
        tw.master.update()
        pause()
        # make sure that value did not change
        self.assertEqual(tw.value,4.123456789)
        # test that out-of-range values are set to 1 - 10
        tw.configure(precision=0)
        self.assertEqual(tw.precision,1)
        tw.configure(precision=11)
        self.assertEqual(tw.precision,10)
        tw.master.update()
        pause()
        tw.master.destroy()


    def test_setContinuous(self):
        """tests setContinuous
        """
        tw = ThumbWheel(width=100, height=26)
        tw.configure(continuous=1)
        tw.master.update()
        pause()
        self.assertEqual(tw.continuous,(1 or True))
        tw.configure(continuous=0)
        self.assertEqual(tw.continuous,(0 or False))
        tw.master.update()
        pause()
        tw.master.destroy()


    def test_setOneTurn(self):
        """tests setOneturn
        """
        tw = ThumbWheel(width=100, height=26)
        tw.configure(oneTurn=23.0)
        tw.master.update()
        self.assertEqual(tw.oneTurn,23.0)
        tw.master.update()
        pause()
        tw.master.destroy()


    def test_setShowLabel(self):
        """tests setShowLAbel
        """
        tw = ThumbWheel(width=100, height=26)
        tw.configure(showLabel=0)
        self.assertEqual(tw.showLabel,0)
        tw.configure(showLabel=1)
        tw.master.update()
        pause()
        self.assertEqual(tw.showLabel,1)
        tw.configure(showLabel=2)
        self.assertEqual(tw.showLabel,2)
        tw.master.update()
        pause()
        tw.master.destroy()

    def test_setOrient_horizontal_to_vertical(self):
        """tests setOrient vertical
        """
        tw = ThumbWheel(width=100, height=50)
        tw.setOrient('vertical')
        tw.set(100)
        tw.master.update()
        pause()
        myevent = Dummyevent(y = 100)
        tw.mouseDown(myevent)
        tw.mouseMove(myevent)
        self.assertEqual(tw.orient,'vertical')
        tw.master.update()
        pause()
        tw.master.destroy()
        
    def test_setOrient_vertical_to_horizontal(self):
        """tests setOrient horizontal
        """
        tw = ThumbWheel(width=100, height=50,orient = 'vertical')
        tw.setOrient('horizontal')
        tw.set(100)
        tw.master.update()
        pause()
        tw.lasty =100
        myevent = Dummyevent(y = 100,x =100)
        tw.mouseDown(myevent)
        tw.mouseMove(myevent)
        self.assertEqual(tw.orient,'horizontal')
        tw.master.update()
        pause()
        tw.master.destroy()
        
    def test_setReportDelta(self):
        """tests setReportDelta
        """
        global mySum
        mySum = 0.
        ### define callback method ###
        def foo(val):
            global mySum
            mySum = mySum + val
            self.assertEqual(val, 1.0)
        tw = ThumbWheel(width=100, height=26, reportDelta=1)
        tw.callbacks.AddCallback(foo)
        tw.master.update()
        pause()
        for i in xrange(10):
            tw.computeNewValue(1.) # emulate mouse rotation
        self.assertEqual(tw.value,10.0)
        #check that the get method returns deltaValue when setReportDelta is True
        self.assertEqual(tw.get(), 1.0)
        #mySum will be 1.0 because the value is deltaVal and never changes
        #so the callback is only called once
        self.assertEqual(mySum, 1.0)
        tw.master.update()
        pause()
        tw.master.destroy()

#invalid input
    def test_thumbwheel_invalid_type(self):
        """tests thumbwheel invalid type
        """
        root = Tkinter.Tk()
        self.assertRaises(AssertionError,ThumbWheel,type ='hai',master = root)
        root.destroy()
        
    def test_thumbwheel_invalid_value(self):
        """tests thumbwheel invlaid value
        """
        root = Tkinter.Tk()
        self.assertRaises(AssertionError,ThumbWheel,value = 'hai',master = root)
        root.destroy()
        
    def test_thumbwheel_invalid_oneTurn(self):
        """tests thumbwheel invalid oneturn
        """
        root = Tkinter.Tk()
        self.assertRaises(AssertionError,ThumbWheel,oneTurn = 'hai',master = root)
        root.destroy()
        
    def test_thumbwheel_invalid_height(self):
        """tests thumbwheel invalid size
        """
        root = Tkinter.Tk()
        self.assertRaises(AssertionError,ThumbWheel,height = 'hai',master = root)
        root.destroy()
        
    def test_thumbwheel_invalid_width(self):
        """tests thumbwheel invalid size
        """
        root = Tkinter.Tk()
        self.assertRaises(AssertionError,ThumbWheel,
                          width = 'hai',master = root)
        root.destroy()
        
    def test_thumbwheel_invalid_callbacks(self):
        """tests thumbwheel invalid callbacks
        """
        root = Tkinter.Tk()
        self.assertRaises(AssertionError,ThumbWheel, callback='jlsd',
                          master = root)
        
        self.assertRaises(AssertionError,ThumbWheel, callback=['jlsd', 'afd'],
                          master = root)
        root.destroy()
        
    def test_thumbwheel_invalid_master(self):
        """tests thumbwheel invalid master
        """
        root = Tkinter.Tk()
        self.assertRaises(AttributeError,ThumbWheel,master = 'asgjf')
        root.destroy()
    
    def test_thumbwheel_invalid_labcfg(self):
        """tests invalid labcfg
        """
        root = Tkinter.Tk()
        self.assertRaises(AttributeError,ThumbWheel,labCfg = 'abcd',master = root)
        root.destroy()
        
    def test_thumbwheel_invalid_min(self):
        """tests thumbwheel invalid min
        """
        root = Tkinter.Tk()
        self.assertRaises(AssertionError,ThumbWheel,min = 'hkdf',master = root)
        root.destroy()
        
    def test_thumbwheel_invalid_max(self):
        """tests thumbwheel invalid max
        """
        root = Tkinter.Tk()
        self.assertRaises(AssertionError,ThumbWheel,max = 'hkdf',master = root)
        root.destroy()
       
        
    def test_thumbwheel_invalid_increment(self):
        """tests invalid increment
        """
        root = Tkinter.Tk()
        self.assertRaises(AssertionError,ThumbWheel,increment = 'hkdf',master = root)
        root.destroy()
    
        
    def test_thumbwheel_invalid_continuous(self):
        """tests invalid continuous
        """
        root = Tkinter.Tk()
        self.assertRaises(AssertionError,ThumbWheel,continuous = 'hkdf')
        root.destroy()
       
    def test_thumbwheel_invalid_precision(self):
        """tests invalid precision
        """
        root = Tkinter.Tk()
        self.assertRaises(AssertionError,ThumbWheel,precision = 'jhgdfj',master =
        root) 
        root.destroy()

    def test_thumbwheel_options_widget_is_mapped(self):
        """tests widget is mapped
        """
        thumbwheel = ThumbWheel(size = 100)
        thumbwheel.opPanel.displayPanel(1)
        thumbwheel.master.update()
        pause()
        self.assertEqual(thumbwheel.opPanel.optionsForm.root.winfo_ismapped(),1)
        thumbwheel.opPanel.optionsForm.root.withdraw()
        thumbwheel.master.destroy() 

#invalid input for configure methods

    def test_invalid_setMin(self):
        """tests invalid input for setMin
        """
        root =Tkinter.Tk()
        tw =ThumbWheel(master = root)
        self.assertRaises(AssertionError,tw.setMin,'hai')
        tw.master.destroy()

    def test_invalid_setMax(self):
        """tests invalid input for setMax
        """
        root =Tkinter.Tk()
        tw =ThumbWheel(master = root)
        self.assertRaises(AssertionError,tw.setMax,'hai')
        tw.master.destroy()
        
    def test_invalid_setIncrement(self):
        """tests invalid input for setIncrement
        """
        root =Tkinter.Tk()
        tw =ThumbWheel(master = root)
        self.assertRaises(AssertionError,tw.setIncrement,'hai')    
        tw.master.destroy()
        
    def test_invalid_setValue(self):
        """tests invalid input for setValue
        """
        root =Tkinter.Tk()
        tw =ThumbWheel(master = root)
        self.assertRaises(AssertionError,tw.setValue,'hai')
        tw.master.destroy()
        
    def test_invalid_setType(self):
        """tests invalid input for setType
        """
        root =Tkinter.Tk()
        tw =ThumbWheel(master = root)
        self.assertRaises(AssertionError,tw.setType,'hai')    
        tw.master.destroy()
        
    def test_invalid_setShowLabel(self):
        """tests invalid input for setShowLabel
        """
        root =Tkinter.Tk()
        tw =ThumbWheel(master = root)
        self.assertRaises(AssertionError,tw.setShowLabel,"iurey")
        tw.master.destroy()
        
    def test_invalid_setLabel(self):
        """tests invalid input for setLabel
        """
        root =Tkinter.Tk()
        tw =ThumbWheel(master = root)
        self.assertRaises(AttributeError,tw.setLabel,'hai')    
        tw.master.destroy()
        
    def test_invalid_setPrecision(self):
        """tests invalid input for setPrecision
        """
        root =Tkinter.Tk()
        tw =ThumbWheel(master = root)
        self.assertRaises(AssertionError,tw.setPrecision,'hai')    
        tw.master.destroy()
        
    def test_invalid_setContinuous(self):
        """tests invalid input for setContinuous
        """
        root =Tkinter.Tk()
        tw =ThumbWheel(master = root)
        self.assertRaises(AssertionError,tw.setContinuous,"iurey")
        tw.master.destroy()
        
    def test_invalid_setOneTurn(self):
        """tests invalid input for setOneTurn
        """
        root =Tkinter.Tk()
        tw =ThumbWheel(master = root)
        self.assertRaises(AssertionError,tw.setOneTurn,'hai')    
        tw.master.destroy()
        
    def test_invalid_setOrient(self):
        """tests invalid input for setOrient
        """
        root =Tkinter.Tk()
        tw =ThumbWheel(master = root)
        self.assertRaises(AssertionError,tw.setOrient,'hai')    
        tw.master.destroy()

#lockmethods
    def test_tw_lockmin(self):
        """tests lockmin
        """
        tw = ThumbWheel(size = 100)
        tw.opPanel.displayPanel(1)
        tw.master.update()
        pause()
        self.assertEqual(tw.opPanel.min_entry.cget('state'),'disabled')
        tw.lockMinCB(1)
        self.assertEqual(tw.opPanel.min_entry.cget('state'),'normal')
        tw.opPanel.Dismiss_cb()
        tw.master.destroy()
        
    def test_tw_lockmax(self):
        """tests lockmax
        """
        tw = ThumbWheel(size = 100)
        tw.opPanel.displayPanel(1)
        tw.master.update()
        pause()
        self.assertEqual(tw.opPanel.max_entry.cget('state'),'disabled')
        tw.lockMaxCB(1)
        self.assertEqual(tw.opPanel.max_entry.cget('state'),'normal')    
        tw.opPanel.Dismiss_cb()
        tw.master.destroy()
        
    def test_tw_lockincrement(self):
        """tests lockincrement
        """
        tw = ThumbWheel(size = 100)
        tw.opPanel.displayPanel(1)
        tw.master.update()
        pause()
        self.assertEqual(tw.opPanel.incr_entry.cget('state'),'disabled')
        tw.lockIncrementCB(1)
        self.assertEqual(tw.opPanel.incr_entry.cget('state'),'normal')    
        tw.opPanel.Dismiss_cb()
        tw.master.destroy()
    
     
    def test_tw_lockOneTurnCB(self):
        """tests lockoneturn
        """
        tw = ThumbWheel(size = 100)
        tw.opPanel.displayPanel(1)
        tw.lockOneTurnCB(1)
        tw.master.update()
        pause()
        self.assertEqual(tw.opPanel.sens_entry.cget('state'),'disabled')
        tw.lockOneTurnCB(0)
        self.assertEqual(tw.opPanel.sens_entry.cget('state'),'normal')
        tw.opPanel.Dismiss_cb()
        tw.master.destroy()
        
    def test_tw_lockValueCB(self):
        """tests lockvalue 
        """
        tw = ThumbWheel(size = 100)
        tw.opPanel.displayPanel(1)
        tw.lockValueCB(1)
        tw.master.update()
        pause()
        self.assertEqual(tw.opPanel.val_entry.cget('state'),'disabled')
        tw.lockValueCB(0)
        self.assertEqual(tw.opPanel.val_entry.cget('state'),'normal')   
        tw.opPanel.Dismiss_cb()
        tw.master.destroy()
        

if __name__ == '__main__':
   unittest.main() 

    

