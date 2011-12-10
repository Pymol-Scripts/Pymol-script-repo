#
# $Header: /opt/cvs/python/packages/share1.5/mglutil/gui/BasicWidgets/Tk/Tests/test_customizedWidgets.py,v 1.3 2005/06/21 16:56:46 sowjanya Exp $
#
# $Id: test_customizedWidgets.py,v 1.3 2005/06/21 16:56:46 sowjanya Exp $
#

import sys,unittest
from mglutil.regression import testplus
from time import sleep
inCB = False
widget = None
wasCalled = 0
def pause(sleepTime=0.2):
    widget.master.update()
    sleep(sleepTime)

#########################################################################
###  TEST KBCOMBOBOX
#########################################################################

class KBComboBox(unittest.TestCase):
    def test_kbComboBox(self):
        global widget
        from mglutil.gui.BasicWidgets.Tk.customizedWidgets import kbComboBox

        # default root widget
        widget = kbComboBox(label_text='Numbers:', labelpos='w')
        widget.pack()
        pause()
        widget.master.destroy()

        # with root widget
        import Tkinter
        root = Tkinter.Toplevel()
        widget = kbComboBox(root, label_text='Numbers:', labelpos='w')
        widget.pack()
        pause()
        widget.master.destroy()


    def test_kbComboBox_setValue(self):
        # test if we can set the value of a kbComboBox
        global widget, wasCalled
        from mglutil.gui.BasicWidgets.Tk.customizedWidgets import kbComboBox

        def chooseNumber(number):
            global wasCalled
            print 'A you chose:', number
            wasCalled=1


        namelist = map(str, range(30))
        widget = kbComboBox(label_text='Numbers:', labelpos='w',
                        selectioncommand = chooseNumber,
                        scrolledlist_items = namelist )
        widget.pack()
        pause()

        # set value and call command
        wasCalled = 0
        widget.set( 8 )
        pause()
        self.assertEqual(wasCalled == 1,True)
        self.assertEqual(widget.get() == '8',True)
        # set value and dio NOT call command
        wasCalled = 0
        widget.set( 2, run=0 )
        pause()
        self.assertEqual(wasCalled == 0,True)
        self.assertEqual(widget.get() == '2',True)
    
        widget.master.destroy()


#########################################################################
###  TEST SLIDERWIDGET
#########################################################################
class SLIDERWIDGETBaseTest(unittest.TestCase):
    def test_sliderwidget(self):
        global widget
        from mglutil.gui.BasicWidgets.Tk.customizedWidgets import SliderWidget
        # default root widget
        widget = SliderWidget()
        widget.pack()
        pause()
        widget.master.destroy()
        # with root widget
        import Tkinter
        root = Tkinter.Toplevel()
        widget = SliderWidget(root)
        widget.pack()
        pause()
        widget.master.destroy()

    def test_sliderwidget_set(self):
        global inCB
        inCB = False
        def slider_cb(val):
            global inCB
            inCB = True
        global widget
        from mglutil.gui.BasicWidgets.Tk.customizedWidgets import SliderWidget
        import Tkinter
        root = Tkinter.Toplevel()
        widget = SliderWidget(root, command=slider_cb)
        widget.pack()
        pause()
        widget.set(9.567, update=0)
        self.assertEqual(not inCB,True)
        val = widget.get()
        self.assertEqual(val == 9.567,True)
        widget.set(10.567, update=1)
        self.assertEqual(inCB,True)
        val = widget.get()
        self.assertEqual(val==10.567,True)
        # need to check if the cursorlabel has the proper format
        pause()
    


        widget.master.destroy()

    def test_sliderwidget_discrete(self):
        global widget
        from mglutil.gui.BasicWidgets.Tk.customizedWidgets import SliderWidget
        import Tkinter
        root = Tkinter.Toplevel()
        widget = SliderWidget(root,label='discrete slider',
                          labelsCursorFormat='%4d', immediate=0,
                          lookup=[10, 15, 25, 46, 78, 99] )
        widget.pack()
        # the widget min, max and val are indices in the lookup table.
        self.assertEqual(widget.min==0,True)
        pause()
        # set a non existent val the val must be the min.
        try:
            widget.set(9.567)
        except ValueError:
            val = widget.get()
            self.assertEqual(val == 10,True)

        pause()
        # set to a value in the lookup table. the widget.val must be the index of
        # this val.
        widget.set(78)
        val = widget.get()
        self.assertEqual(val == 78,True)
        self.assertEqual(widget.val == widget.lookup.index(78),True)
        pause()

        # set to a value greater than the last value of the look up table.
        try:
            widget.set(110)

        except ValueError:
            val = widget.get()
            self.assertEqual(val == 78,True)

        pause()

        # set the min value of the discrete slider
        widget.setMin(15)
        val = widget.get()
        self.assertEqual(val == 15,True)
        self.assertEqual(widget.min == 1,True)

        widget.set(10)
        val = widget.get()
        self.assertEqual(val == widget.lookup[widget.min],True)

        pval = widget.get()
        widget.setMax(78)
        val = widget.get()
        self.assertEqual(val == pval,True)
        self.assertEqual(widget.max == widget.lookup.index(78),True)
        widget.set(99)
        self.assertEqual(widget.get() == widget.lookup[widget.max],True)

        pause()

        widget.master.destroy()
    
#########################################################################
###  TEST EXTENDEDSLIDERWIDGET
#########################################################################
class EXTENDEDSLIDERWIDGETBaseTest(unittest.TestCase):
    def test_extendedsliderwidget_discrete(self):
        global widget
        from mglutil.gui.BasicWidgets.Tk.customizedWidgets import ExtendedSliderWidget
        import Tkinter
        root = Tkinter.Toplevel()
        widget = ExtendedSliderWidget(root,label='discrete extendedslider',
                                  labelsCursorFormat='%4d', immediate=0,
                                  lookup=[10, 15, 25, 46, 78, 99] )
        widget.pack()
        # the widget min, max and val are indices in the lookup table.
        self.assertEqual(widget.min==0,True)
        pause()
        # set a non existent val the val must be the min.
        try:
            widget.set(9.567)
        except ValueError:
            val = widget.get()
            self.assertEqual(val == 10,True)
        
        # need to check the value of the entry though.
        entryval = widget.entryContent.get()
        self.assertEqual(entryval == widget.labelsCursorFormat%val,True)
        pause()

        # set to a value in the lookup table. the widget.val must be the index of
        # this val.
        widget.set(78)
        val = widget.get()
        self.assertEqual(val == 78,True)
        self.assertEqual(widget.val == widget.lookup.index(78),True)
        # need to check the value of the entry though.
        entryval = widget.entryContent.get()
        self.assertEqual(entryval == widget.labelsCursorFormat%val,True)
        pause()

        # set to a value greater than the last value of the look up table.
        try:
            widget.set(110)
        except ValueError:
            val = widget.get()
            self.assertEqual(val == 78,True)
        #self.assertEqual(val == widget.lookup[widget.max]
        #self.assertEqual(widget.val == widget.max
        pause()

        # set the min value of the discrete slider
        widget.setMin(15)
        val = widget.get()
        self.assertEqual(val == 15,True)
        self.assertEqual(widget.min == 1,True)

        widget.set(10)
        val = widget.get()
        self.assertEqual(val == widget.lookup[widget.min],True)

        widget.setMax(78)
        self.assertEqual(widget.max == widget.lookup.index(78),True)
        widget.set(99)
        self.assertEqual(widget.get() == widget.lookup[widget.max],True)

        # TEST THE CALLBACK OF THE entry
        widget.entryContent.set(46)
        widget.setval('<Return>')
        self.assertEqual(widget.get() == 46,True)
        pause()

        try:
            widget.entryContent.set(110)
            widget.setval('<Return>')
        except ValueError:
            self.assertEqual(widget.get() == 46,True)

        pause()

        widget.master.destroy()

    def test_extendedsliderwidget(self):
        global widget
        from mglutil.gui.BasicWidgets.Tk.customizedWidgets import ExtendedSliderWidget
        import Tkinter
        root = Tkinter.Toplevel()
        widget = ExtendedSliderWidget(root)
        widget.pack()
        pause()
        # TEST SET method and the labelsCursorFormat
        widget.set(9.567)
        val = widget.get()
        self.assertEqual(val == 9.567,True)
        # need to check the value of the entry though.
        entryval = widget.entryContent.get()
        self.assertEqual(entryval == widget.labelsCursorFormat%val,True)
        pause()

        # TEST THE CALLBACK OF THE entry
        widget.entryContent.set(12.956)
        widget.setval('<Return>')
        self.assertEqual(widget.get() == 12.956,True)

        # TEST SETMIN
        widget.setMin(3.564)
        widget.set(0.0)
        val = widget.get()
        self.assertEqual(val == 3.564,True)
        pause()

        widget.setMin(5.)
        val = widget.get()
        self.assertEqual(val == 5.,True)
        pause()

        # TEST SETMAX
        widget.setMax(30.566)
        val = widget.get()
        self.assertEqual(val == 5.0,True)
        pause()

        widget.set(25.9)
        self.assertEqual(widget.get()==25.9,True)
        widget.setMax(20.)
        self.assertEqual(widget.get()==20.,True)
        pause()
    
        widget.master.destroy()
    

if __name__ == '__main__':
    unittest.main()
    
