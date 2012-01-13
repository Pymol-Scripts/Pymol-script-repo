#
#
#
#$Id: test_colorWidgets.py,v 1.3 2010/09/30 19:29:12 annao Exp $
#
###########################################################################
#
#   Authors : Sowjanya Karnati,Michel F Sanner
#
###########################################################################
#
#
#

import sys,unittest,Tkinter
from time import sleep
from mglutil.gui.BasicWidgets.Tk.colorWidgets import *
from mglutil.util.colorUtil import *

     
def pause():
    sleep(0.1)

class colorEditorTest(unittest.TestCase):

#####set and get method tests

    def test_colorEditor_set_HSV(self):
        """tests setting color when mode is HSV"""
        self.master = Tkinter.Toplevel()       
        editFrame = Tkinter.Frame(self.master,borderwidth=2,relief='ridge')
        ce = ColorEditor(editFrame)
        ce.pack()
        editFrame.pack()
        color=(0.5,0.5,0.5)
        ce.set(color=color,mode='HSV')
        rgbcolor=ToRGB(color)
        for i in range(0,len(rgbcolor)-1):
            self.assertEqual(rgbcolor[i],ce.get()[i])
        ce.master.update()
        pause()
        ce.master.master.destroy()
        
    def test_colorEditor_set_HEX(self):
        """tests setting color when mode is HEX"""
        self.master = Tkinter.Toplevel()       
        editFrame = Tkinter.Frame(self.master,borderwidth=2,relief='ridge')
        ce = ColorEditor(editFrame)
        ce.pack()
        editFrame.pack()
        color='#FFFF00'
        ce.set(color=color,mode='HEX')
        rgbcolor=ToRGB(color,mode="HEX")
        for i in range(0,len(rgbcolor)-1):
            self.assertEqual(rgbcolor[i],ce.get()[i])   
        ce.master.update()
        pause()
        ce.master.master.destroy()
        
    def test_colorEditor_set_RGB(self):
        """tests setting color when mode is RGB"""
        self.master = Tkinter.Toplevel()       
        editFrame = Tkinter.Frame(self.master,borderwidth=2,relief='ridge')
        ce = ColorEditor(editFrame)
        ce.pack()
        editFrame.pack()
        color=(1.0, 0.0, 0.0)
        ce.set(color=color,mode='RGB')
        for i in range(0,len(color)-1):    
            self.assertEqual(ce.get()[i],color[i])
        ce.master.update()
        pause()
        ce.master.master.destroy()
   
        

#####updateWidgetColor method tests


    def test_colorEditor_updateWidgetColor_HSV(self):
        """tests updateWidgetColor when mode is 'hsv' """
        self.master = Tkinter.Toplevel()       
        editFrame = Tkinter.Frame(self.master,borderwidth=2,relief='ridge')
        ce = ColorEditor(editFrame)
        ce.pack()
        editFrame.pack()
        old_hVal = ce.hVal.get()
        old_sVal = ce.sVal.get()
        old_vVal = ce.vVal.get()
        ce.master.update()
        pause()
        ce.updateWidgetsColor((0.0,1.0,0.0),who='hsv')
        ce.master.update()
        pause()
        new_hVal = ce.hVal.get()
        new_sVal = ce.sVal.get()
        new_vVal = ce.vVal.get()
        self.assertEqual(old_hVal!=new_hVal ,True)
        self.assertEqual(old_sVal!=new_sVal ,True)
        #self.assertEqual(old_vVal!=new_vVal ,True)
        ce.master.update()
        pause()
        ce.master.master.destroy() 

    def test_colorEditor_updateWidgetColor_RGB(self):
        """tests updateWidgetColor when mode is 'rgb' """
        self.master = Tkinter.Toplevel()       
        editFrame = Tkinter.Frame(self.master,borderwidth=2,relief='ridge')
        ce = ColorEditor(editFrame)
        ce.pack()
        editFrame.pack()
        old_rVal = ce.rVal.get()
        old_gVal = ce.gVal.get()
        old_bVal = ce.bVal.get()
        ce.updateWidgetsColor((0.5,0.5,0.5),who='rgb')
        new_rVal = ce.rVal.get()
        new_gVal = ce.gVal.get()
        new_bVal = ce.bVal.get()
        self.assertEqual(old_rVal!=new_rVal ,True)
        self.assertEqual(old_gVal!=new_gVal ,True)
        self.assertEqual(old_bVal!=new_bVal ,True)
        ce.master.update()
        pause()
        ce.master.master.destroy()

    def test_colorEditor_updateWidgetColor_HEX(self):
        """tests updateWidgetColor when mode is 'hex' """
        self.master = Tkinter.Toplevel()       
        editFrame = Tkinter.Frame(self.master,borderwidth=2,relief='ridge')
        ce = ColorEditor(editFrame)
        ce.pack()
        editFrame.pack()
        old_hexVal = ce.hexVal.get()
        ce.updateWidgetsColor((0.5,1.0,0.5),who='hex')
        new_hexVal = ce.hexVal.get()
        self.assertEqual(old_hexVal!=new_hexVal,True)
        ce.master.update()
        pause()
        ce.master.master.destroy()

###########Entry Val Tests################


    def test_colorEditor_rVal(self):
        """tests colorEditor values by setting rVal and invoking """
        self.master = Tkinter.Toplevel()       
        editFrame = Tkinter.Frame(self.master,borderwidth=2,relief='ridge')
        ce = ColorEditor(editFrame)
        ce.pack()
        editFrame.pack()    
        oldcol=ce.get()
        ce.rVal.setvalue(0.5)
        ce.master.update()
        pause()
        ce.rVal.invoke()
        newcol = ce.get()
        self.assertEqual(oldcol!=newcol,True)
        ce.master.update()
        pause()
        ce.master.master.destroy()
        
    def test_colorEditor_gVal(self):
        """tests colorEditor values by setting gVal and invoking """
        self.master = Tkinter.Toplevel()       
        editFrame = Tkinter.Frame(self.master,borderwidth=2,relief='ridge')
        ce = ColorEditor(editFrame)
        ce.pack()
        editFrame.pack()    
        oldcol=ce.get()
        ce.gVal.setvalue(0.5)
        ce.master.update()
        pause()
        ce.gVal.invoke()
        newcol = ce.get()
        self.assertEqual(oldcol!=newcol,True)
        ce.master.update()
        pause()
        ce.master.master.destroy()

    def test_colorEditor_bVal(self):
        """tests colorEditor values by setting bVal and invoking """
        self.master = Tkinter.Toplevel()       
        editFrame = Tkinter.Frame(self.master,borderwidth=2,relief='ridge')
        ce = ColorEditor(editFrame)
        ce.pack()
        editFrame.pack()    
        oldcol=ce.get()
        ce.bVal.setvalue(0.5)
        ce.master.update()
        pause()
        ce.bVal.invoke()
        newcol = ce.get()
        self.assertEqual(oldcol!=newcol,True)
        ce.master.update()
        pause()
        ce.master.master.destroy()

    def test_colorEditor_hVal(self):
        """tests colorEditor values by setting hVal and invoking """
        self.master = Tkinter.Toplevel()       
        editFrame = Tkinter.Frame(self.master,borderwidth=2,relief='ridge')
        ce = ColorEditor(editFrame)
        ce.pack()
        editFrame.pack()    
        oldcol=ce.get()
        #set vVal,sVal also since they shouldn't be 0.0 and 1.0
        ce.sVal.setvalue(0.32)
        ce.sVal.setvalue(0.32)
        ce.hVal.setvalue(0.32)
        ce.vVal.invoke()
        ce.sVal.invoke()
        ce.hVal.invoke()
        newcol = ce.get()
        self.assertEqual(oldcol!=newcol,True)
        ce.master.update()
        pause()
        ce.master.master.destroy()

    def test_colorEditor_sVal(self):
        """tests colorEditor values by setting sVal and invoking """
        self.master = Tkinter.Toplevel()       
        editFrame = Tkinter.Frame(self.master,borderwidth=2,relief='ridge')
        ce = ColorEditor(editFrame)
        ce.pack()
        editFrame.pack()    
        oldcol=ce.get()
        ce.sVal.setvalue(0.5)
        ce.master.update()
        pause()
        ce.sVal.invoke()
        newcol = ce.get()
        self.assertEqual(oldcol!=newcol,True)
        ce.master.update()
        pause()
        ce.master.master.destroy()

    def test_colorEditor_vVal(self):
        """tests colorEditor values by setting vVal and invoking """
        self.master = Tkinter.Toplevel()       
        editFrame = Tkinter.Frame(self.master,borderwidth=2,relief='ridge')
        ce = ColorEditor(editFrame)
        ce.pack()
        editFrame.pack()    
        oldcol=ce.get()
        ce.vVal.setvalue(0.5)
        ce.vVal.invoke()
        newcol = ce.get()
        self.assertEqual(oldcol!=newcol,True)
        ce.master.update()
        pause()
        ce.master.master.destroy()

    def test_colorEditor_hexVal(self):
        """tests colorEditor values by setting hexVal and invoking """
        self.master = Tkinter.Toplevel()       
        editFrame = Tkinter.Frame(self.master,borderwidth=2,relief='ridge')
        ce = ColorEditor(editFrame)
        ce.pack()
        editFrame.pack()    
        oldcol=ce.get()
        ce.hexVal.setvalue("#FFFF00")
        ce.hexVal.invoke()
        newcol = ce.get()
        self.assertEqual(oldcol!=newcol,True)
        ce.master.update()
        pause()
        ce.master.master.destroy()

    
        
class ColorChooserTest(unittest.TestCase):
    
     def test_colorChooser_1(self):
        """tests colorchooser by invoking color radio select buttons """
        self.master = Tkinter.Toplevel()       
        self.masterFrame = Tkinter.Frame(self.master,borderwidth=2,relief='ridge')
        self.menuBar = Pmw.MenuBar(self.masterFrame,
                                   hull_relief = 'raised',
                                   hull_borderwidth = 1)
        self.mainFrame = Tkinter.Frame(self.masterFrame,
                                       borderwidth=2, relief='ridge',
                                        width=150, height=200)
        cc = ColorChooser(self.mainFrame)
        cc.pack()
        self.mainFrame.pack()
        self.masterFrame.pack()
        pause()
        cc.master.update()
        #cc.editColor()
        #cc.colorChips.invoke(0)    
        
        cc.master.update()
        pause()
        color = cc.ce.get()
        self.assertEqual(color,[1.0, 1.0, 1.0])
        cc.master.master.master.destroy() 
        
  
if __name__ == '__main__':
    unittest.main()


        
