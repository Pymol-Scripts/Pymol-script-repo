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

# $Id: test_dial.py,v 1.8 2009/05/28 22:07:20 vareille Exp $

import sys,unittest,Tkinter
from time import sleep
from mglutil.gui.BasicWidgets.Tk.Dial import Dial
wasCalled = False
def MyCallback(val):
    global wasCalled
    print val
    wasCalled = True
def pause():
    import time
    time.sleep(0.1)
    

class DialBaseTest(unittest.TestCase):
   
    def test_setValue(self):
        # test setting of a value
        dial = Dial(size=50, value=10.0)
        dial.master.update()
        pause()
        self.assertEqual(dial.value == 10.0,True)
        dial.configure(type=float)
        dial.set(20.0)
        self.assertEqual(dial.value == 20.0,True)
        dial.master.update()
        pause()
        dial.master.withdraw()
        dial.master.destroy()



    def test_setValueWithIncrement(self):
        
        dial = Dial(size=50, value=5.0)
        dial.master.update()
        pause()
        self.assertEqual(dial.increment == 0.0,True) # increment is 0.0 if not set
        self.assertEqual(dial.value == 5.0,True)
        # now set increment and then set a new value
        dial.configure(increment=2.0)
        self.assertEqual(dial.increment == 2.0,True)
        dial.set(21.0)
        #set value is 21 but 23 is shown
        #print dial.value
        self.assertEqual(dial.value == 21.0,True)
        dial.master.update()
        pause()
        dial.master.withdraw()
        dial.master.destroy()

       
    def test_callback(self):
        # test callback
        global wasCalled 
        dial = Dial(size=50, value=10.0)
        dial.callbacks.AddCallback(MyCallback)
        dial.master.update()
        pause()
        self.assertEqual(len(dial.callbacks.callbacks),1)
        dial.callbacks.RemoveCallback(MyCallback)
        self.assertEqual(len(dial.callbacks.callbacks),0)
        # setValue(val) should NOT call callback
        dial.callbacks.AddCallback(MyCallback)
        dial.setValue(5.0)
        dial.master.update()
        pause()
        if wasCalled == True:
            raise RuntimeError
        # set(val) should call callback
        dial.set(6.0)
        dial.master.update()
        pause()
        if wasCalled == False:
            raise RuntimeError
        wasCalled = False
        dial.set(7.0, update=0)
        if wasCalled:
            raise RuntimeError
        dial.master.destroy()    
       
    def test_callback_2(self):
        global wasCalled1
        global wasCalled2
        
        def foo1(val):
            global wasCalled1
            wasCalled1=True

        def foo2(val):
            global wasCalled2
            wasCalled2=True
            
        dial = Dial(size=50, value=10.0, callback=[foo1,foo2])
        self.assertEqual(dial.callbacks.callbacks, [foo1,foo2], "Expecting to have foo added to the callbackmanager callbacks list")
        # setValue(val) should NOT call callback
        wasCalled1 = False
        wasCalled2 = False
        dial.setValue(5.0)
        if wasCalled1 and wasCalled2:
            raise RuntimeError

        # set(val) should call callback
        dial.set(6.0)
        if not wasCalled1 and not wasCalled2:
            raise RuntimeError

        wasCalled1 = False
        wasCalled2 = False
        dial.set(7.0, update=0)
        if wasCalled1 and wasCalled2:
            raise RuntimeError
        dial.master.update()
        pause()
        dial.master.destroy()

    def test_setType(self):
        # test setting of type
        dial = Dial(size=50, type='float')
        dial.set(100.0)
        dial.master.update()
        pause()
        # we can mix configure(type=float) and configure(type='float')
        dial.configure(type=float)
        self.assertEqual(dial.type == float,True)
        self.assertEqual(type(dial.get()) == type(1.0),True)
        dial.configure(type='int')
        self.assertEqual(dial.type == int,True)
        self.assertEqual(type(dial.get()) == type(1),True)
        dial.master.update()
        pause()
        dial.master.withdraw()
        dial.master.destroy()
    

    def test_setMin(self):
        global wasCalled
        # test that we do not call the callback when we set the minimum value
        dial = Dial(size=50)
        dial.callbacks.AddCallback(MyCallback)
        # set the value without calling callback
        dial.setValue(4)
        # change minimum without affecting value ==> nocallback
        dial.configure(min=2)
        dial.master.update()
        pause()
        if wasCalled == True:
            raise RuntimeError
        # change minimum with affecting value ==> callback
        dial.configure(min=6)
        if wasCalled == False:
            raise RuntimeError
        wasCalled = False
        dial.master.update()
        pause()
        dial.master.destroy()

    
    def test_setMax(self):
        global wasCalled
        # test that we do not call the callback when we set the maximum value
        dial = Dial(size=50)
        dial.callbacks.AddCallback(MyCallback)
        dial.master.update()
        pause()
        # set the value without calling callback
        dial.setValue(4)
        # change maximum without affecting value ==> nocallback
        dial.configure(max=6)
        if wasCalled == True:
            raise RuntimeError
        # change maximum with affecting value ==> callback
        dial.configure(max=2)
        if wasCalled == False:
            raise RuntimeError
        wasCalled = False
        dial.master.update()
        pause()
        dial.master.destroy()


    def test_setIncrement(self):
        """tests set increment 
        """
        dial = Dial(size=50)
        dial.configure(increment=5)
        self.assertEqual(dial.increment == 5,True)
        dial.master.update()
        pause()
        dial.master.destroy()

    
    def test_setPrecision(self):
        """tests set precision
        """
        dial = Dial(size=50)
        dial.configure(type='float')
        dial.configure(precision=5)
        dial.master.update()
        pause()
        # increase precision to 5 - this is only visual. value is not changed
        self.assertEqual(dial.precision == 5,True)
        dial.set(4.123456789)
        dial.master.update()
        pause()
        # make sure that value did not change
        self.assertEqual(dial.value == 4.123456789,True)
        dial.master.destroy()


    def test_setContinuous(self):
        """tests setcontinuous
        """
        dial = Dial(size=50)
        dial.configure(continuous=1)
        dial.master.update()
        pause()
        self.assertEqual(dial.continuous == 1,True)
        dial.configure(continuous=0)
        dial.master.update()
        pause()
        self.assertEqual(dial.continuous == (0 or None),True)
        dial.master.destroy()


    def test_setOneTurn(self):
        """tests setoneturn
        """
        dial = Dial(size=50)
        dial.configure(oneTurn=23.0)
        dial.master.update()
        pause()
        self.assertEqual(dial.oneTurn == 23.0,True)
        dial.master.destroy()


    def test_setShowLabel(self):
        """tests set showlabel
        """
        dial = Dial(size=50)
        dial.configure(showLabel=0)
        dial.master.update()
        pause()
        self.assertEqual(dial.showLabel == 0,True)
        dial.configure(showLabel=1)
        dial.master.update()
        pause()
        self.assertEqual(dial.showLabel == 1,True)
        dial.configure(showLabel=2)
        dial.master.update()
        pause()
        self.assertEqual(dial.showLabel == 2,True)
        dial.master.destroy()

    def test_setArrow(self):
        """tests set arrow
        """
        dial = Dial(size=50)
        dial.setArrow(200)
        dial.master.update()
        pause()
        #arrowLength = max(3,3*size/40)
        self.assertEqual(dial.arrowLength,15)
        dial.setArrow(160)
        dial.master.update()
        pause()
        self.assertEqual(dial.arrowLength,12)
        dial.master.destroy()
        
    def test_dial_default_parameters(self):
        """tests set default parameters
        """
        dial = Dial(master=None, type='float',
                 labCfg={'fg':'black','side':'left', 'text':None},
                 min=None, max=None, increment=0.0, precision=2,
                 showLabel=1, value=0.0, continuous=1, oneTurn=360.,
                 size=50, callback=None,
                 lockMin=0, lockBMin=0, lockMax=0, lockBMax=0,
                 lockIncrement=0, lockBIncrement=0,
                 lockPrecision=0, lockShowLabel=0, lockValue=0,
                 lockType=0, lockContinuous=0, lockOneTurn=0)
        self.assertEqual(str(dial.opPanel.type),"<type 'float'>")
        self.assertEqual(dial.opPanel.minInput.get(),'')
        self.assertEqual(dial.opPanel.maxInput.get(),'')
        self.assertEqual(dial.opPanel.incrInput.get(),'')
        dial.opPanel.displayPanel(1)
        dial.master.update()
        pause()
        self.assertEqual(dial.precision,2)
        self.assertEqual(dial.showLabel,1)
        self.assertEqual(dial.value,0.0)
        self.assertEqual(dial.continuous,1)
        self.assertEqual(dial.oneTurn,360)
        self.assertEqual(dial.size,50)
        dial.opPanel.Apply_cb()
        dial.master.destroy()

    def test_set_arguements(self):
        """tests set arguements
        """
        root =Tkinter.Tk()
        dial = Dial(master=root,type='int',labCfg={'fg':'blue','side':'left', 'text':'hai'},min = 2,max = 100, increment = 2,precision =4,showlabel =0,value =10,continuous =0,oneTurn =180,size = 100,callback = None)
        
        dial.opPanel.displayPanel(1)
        dial.master.update()
        pause()
        self.assertEqual(str(dial.opPanel.type),"<type 'float'>")
        self.assertEqual(dial.opPanel.minInput.get(),'2')
        self.assertEqual(dial.opPanel.maxInput.get(),'100')
        self.assertEqual(dial.opPanel.incrInput.get(),'2')
        self.assertEqual(dial.precision,4)
        self.assertEqual(dial.labCfg,{'fg':'blue','side':'left','text':'hai'})
        self.assertEqual(dial.value,10)
        self.assertEqual(dial.continuous,None)
        dial.opPanel.Apply_cb()
        dial.master.destroy()
        
    def test_dial_set_lock_arguements(self):
        """tests set lock arguements
        """
        dial = Dial(size = 100,lockMin=1, lockBMin=1, lockMax=1, lockBMax=1,lockIncrement=1, lockBIncrement=1,lockPrecision=1, lockShowLabel=1, lockValue=1,lockType=1, lockContinuous=1, lockOneTurn=1)
        dial.master.update()
        pause()
        self.assertEqual(dial.lockMin,1)
        self.assertEqual(dial.lockMax,1)
        self.assertEqual(dial.lockBMin,1)
        self.assertEqual(dial.lockBMax,1)
        self.assertEqual(dial.lockIncrement,1)
        self.assertEqual(dial.lockIncrement,1)
        self.assertEqual(dial.lockValue,1)
        self.assertEqual(dial.lockType,1)
        self.assertEqual(dial.lockContinuous,1)
        self.assertEqual(dial.lockOneTurn,1)
        dial.master.destroy()
        
#Invalid Inputs 

    def test_dial_invalid_type(self):
        """tests dial invalid type
        """
        root =Tkinter.Tk()
        self.assertRaises(AssertionError,Dial,type ='hai',master =root)
        root.destroy()
        
    def test_dial_invalid_value(self):
        """tests dial invlaid value
        """
        root =Tkinter.Tk()
        self.assertRaises(ValueError,Dial,value = 'hai',master =root)
        root.destroy()
        
    def test_dial_invalid_oneTurn(self):
        """tests dial invalid oneturn
        """
        root =Tkinter.Tk()
        self.assertRaises(AssertionError,Dial,oneTurn = 'hai',master =root)
        root.destroy()
        
    def test_dial_invalid_size(self):
        """tests dial invalid size
        """
        root =Tkinter.Tk()
        self.assertRaises(AssertionError,Dial,size = 'hai',master =root)
        root.destroy()
    
    def test_dial_invalid_callback(self):
        """tests dial invalid callback
        """
        root =Tkinter.Tk()
        self.assertRaises(AssertionError,Dial,callback ='jlsd',master =root)
        root.destroy()
        
    def test_dial_invalid_master(self):
        """tests dial invalid master
        """
        root =Tkinter.Tk()
        self.assertRaises(AttributeError,Dial,master = 'asgjf')
        root.destroy()
        
    def test_dial_invalid_labcfg(self):
        """tests invalid labcfg
        """
        root =Tkinter.Tk()
        self.assertRaises(AttributeError,Dial,labCfg = 'abcd',master =root)
        root.destroy()
        
    def test_dial_invalid_min(self):
        """tests dial invalid min
        """
        root =Tkinter.Tk()
        self.assertRaises(AssertionError,Dial,min = 'hkdf',master =root)
        root.destroy()
        
    def test_dial_invalid_max(self):
        """tests dial invalid max
        """
        root =Tkinter.Tk()
        
        self.assertRaises(AssertionError,Dial,max = 'hkdf',master =root)
        
        root.destroy()
        
    def test_dial_invalid_increment(self):
        """tests invalid increment
        """
        root =Tkinter.Tk()
        self.assertRaises(AssertionError,Dial,increment = 'hkdf',master =root)
        root.destroy()
        
    def test_dial_invalid_continuous(self):
        """tests invalid continuous:failed
        """
        root =Tkinter.Tk()
        self.assertRaises(AssertionError,Dial,continuous = 'hkdf',master =root)
        root.destroy()
        
    def test_dial_invalid_precision(self):
        """tests invalid precision
        """
        root =Tkinter.Tk()
        self.assertRaises(AssertionError,Dial,precision = 'jhgdfj',master =root) 
        root.destroy()

    def test_dial_options_widget_is_mapped(self):
        """tests widget is mapped
        """
        dial = Dial(size = 100)
        dial.opPanel.displayPanel(1)
        dial.master.update()
        pause()
        self.assertEqual(dial.opPanel.optionsForm.root.winfo_ismapped(),1)
        dial.opPanel.optionsForm.root.withdraw()
        dial.master.destroy() 

#lockmethods
    def test_dial_lockmin(self):
        """tests lockmin
        """
        dial = Dial(size = 100)
        dial.opPanel.displayPanel(1)
        dial.master.update()
        pause()
        self.assertEqual(dial.opPanel.min_entry.cget('state'),'disabled')
        dial.lockMinCB(1)
        self.assertEqual(dial.opPanel.min_entry.cget('state'),'normal')
        dial.opPanel.Dismiss_cb()
        dial.master.destroy()
        
    def test_dial_lockmax(self):
        """tests lockmax
        """
        dial = Dial(size = 100)
        dial.opPanel.displayPanel(1)
        dial.master.update()
        pause()
        self.assertEqual(dial.opPanel.max_entry.cget('state'),'disabled')
        dial.lockMaxCB(1)
        self.assertEqual(dial.opPanel.max_entry.cget('state'),'normal')    
        dial.opPanel.Dismiss_cb()
        dial.master.destroy()
        
    def test_dial_lockincrement(self):
        """tests lockincrement
        """
        dial = Dial(size = 100)
        dial.opPanel.displayPanel(1)
        dial.master.update()
        pause()
        self.assertEqual(dial.opPanel.incr_entry.cget('state'),'disabled')
        dial.lockIncrementCB(1)
        self.assertEqual(dial.opPanel.incr_entry.cget('state'),'normal')    
        dial.opPanel.Dismiss_cb()
        dial.master.destroy()
    
     
    def test_dial_lockOneTurnCB(self):
        """tests lockoneturn
        """
        dial = Dial(size = 100)
        dial.opPanel.displayPanel(1)
        dial.lockOneTurnCB(1)
        dial.master.update()
        pause()
        self.assertEqual(dial.opPanel.sens_entry.cget('state'),'disabled')
        dial.lockOneTurnCB(0)
        self.assertEqual(dial.opPanel.sens_entry.cget('state'),'normal')
        dial.opPanel.Dismiss_cb()
        dial.master.destroy()
        
    def test_dial_lockValueCB(self):
        """tests lockvalue 
        """
        dial = Dial(size = 100)
        dial.opPanel.displayPanel(1)
        dial.lockValueCB(1)
        dial.master.update()
        pause()
        self.assertEqual(dial.opPanel.val_entry.cget('state'),'disabled')
        dial.lockValueCB(0)
        self.assertEqual(dial.opPanel.val_entry.cget('state'),'normal')   
        dial.opPanel.Dismiss_cb()
        dial.master.destroy()
        
#invalid input for configure methods


    def test_invalid_setMin(self):
        """tests invalid input for setMin
        """
        root =Tkinter.Tk()
        dial =Dial(master = root)
        self.assertRaises(AssertionError,dial.setMin,'hai')
        dial.master.destroy()

    def test_invalid_setMax(self):
        """tests invalid input for setMax
        """
        root =Tkinter.Tk()
        dial =Dial(master = root)
        self.assertRaises(AssertionError,dial.setMax,'hai')
        dial.master.destroy()
        
    def test_invalid_setIncrement(self):
        """tests invalid input for setIncrement
        """
        root =Tkinter.Tk()
        dial =Dial(master = root)
        self.assertRaises(AssertionError,dial.setIncrement,'hai')    
        dial.master.destroy()
        
    def test_invalid_setValue(self):
        """tests invalid input for setValue
        """
        root =Tkinter.Tk()
        dial =Dial(master = root)
        self.assertRaises(ValueError,dial.setValue,'hai')
        dial.master.destroy()
        
    def test_invalid_setType(self):
        """tests invalid input for setType
        """
        root =Tkinter.Tk()
        dial =Dial(master = root)
        self.assertRaises(AssertionError,dial.setType,'hai')    
        dial.master.destroy()
        
    def test_invalid_setShowLabel(self):
        """tests invalid input for setShowLabel
        """
        root =Tkinter.Tk()
        dial =Dial(master = root)
        self.assertRaises(AssertionError,dial.setShowLabel,"iurey")
        dial.master.destroy()
        
    def test_invalid_setLabel(self):
        """tests invalid input for setLabel
        """
        root =Tkinter.Tk()
        dial =Dial(master = root)
        self.assertRaises(AttributeError,dial.setLabel,'hai')    
        dial.master.destroy()
        
    def test_invalid_setPrecision(self):
        """tests invalid input for setPrecision
        """
        root =Tkinter.Tk()
        dial =Dial(master = root)
        self.assertRaises(AssertionError,dial.setPrecision,'hai')    
        dial.master.destroy()
        
    def test_invalid_setContinuous(self):
        """tests invalid input for setContinuous
        """
        root =Tkinter.Tk()
        dial =Dial(master = root)
        self.assertRaises(AssertionError,dial.setContinuous,"iurey")
        dial.master.destroy()
        
    def test_invalid_setOneTurn(self):
        """tests invalid input for setOneTurn
        """
        root =Tkinter.Tk()
        dial =Dial(master = root)
        self.assertRaises(AssertionError,dial.setOneTurn,'hai')    
        dial.master.destroy()
        
if __name__ == '__main__':
   unittest.main() 
