#########################################################################
#
# Date: Jun 2002 Authors: Daniel Stoffler, Michel Sanner
#
#    stoffler@scripps.edu
#    sanner@scripps.edu
#
# Copyright:  Daniel Stoffler, Sowjanya Karnati, Michel Sanner, and TSRI
#
#########################################################################

# $Id: test_optionsPanel.py,v 1.4 2005/06/20 17:32:03 sowjanya Exp $

import sys,unittest,Tkinter
from time import sleep
from mglutil.gui.BasicWidgets.Tk.optionsPanel import OptionsPanel,VectorOptionsPanel

widget = None
wasCalled = 0
def pause():
    import time
    time.sleep(0.1)
class OptionsPanelBaseTest(unittest.TestCase):

    def test_inpval(self):
        """checks value entry in options panel widget
        """
        oppanel  = OptionsPanel(title="My Options Panel")
        oppanel.displayPanel(create=1)
        oppanel.idf.entryByName['inpVal']['wcfg']['textvariable'].set(1)
        oppanel.master.update()
        pause()
        self.assertEqual(oppanel.valInput.get(),'1')
        oppanel.master.destroy()

    
    def test_inpMin(self):
        """checks minimum value entry in options panel widget
        """
        oppanel  = OptionsPanel(title="My Options Panel")
        oppanel.displayPanel(create=1)
        oppanel.bmin_entry.invoke()
        oppanel.idf.entryByName['inpMin']['wcfg']['textvariable'].set(1)
        oppanel.master.update()
        pause()
        self.assertEqual(oppanel.minInput.get(),'1')
        oppanel.master.destroy()
    
    
    def test_inpMax(self):
        """checks maximum value entry in options panel widget
        """
        oppanel  = OptionsPanel(title="My Options Panel")
        oppanel.displayPanel(create=1)
        oppanel.bmax_entry.invoke()
        oppanel.idf.entryByName['inpMax']['wcfg']['textvariable'].set(1)
        oppanel.master.update()
        pause()
        self.assertEqual(oppanel.maxInput.get(),'1')
        oppanel.master.destroy()
        
    def test_inputIncr(self):
        """checks increment value entry in options panel widget
        """
        oppanel  = OptionsPanel(title="My Options Panel")
        oppanel.displayPanel(create=1)
        oppanel.idf.entryByName['inpIncr']['wcfg']['textvariable'].set(1)
        oppanel.master.update()
        pause()
        self.assertEqual(oppanel.incrInput.get(),'1')
        oppanel.master.destroy()
    
    def test_togglemin(self):
        """checks toggle min button
        """
        oppanel  = OptionsPanel(title="My Options Panel")
        oppanel.displayPanel(create=1)
        oppanel.idf.entryByName['togMin']['widget'].invoke()
        oppanel.master.update()
        pause()
        self.assertEqual(oppanel.min_entry.cget('state'),'normal')
        oppanel.master.destroy()
        
    def test_togglemax(self):
        """checks toggle max button
        """
        oppanel  = OptionsPanel(title="My Options Panel")
        oppanel.displayPanel(create=1)
        oppanel.idf.entryByName['togMax']['widget'].invoke()
        oppanel.master.update()
        pause()
        self.assertEqual(oppanel.max_entry.cget('state'),'normal')
        oppanel.master.destroy()
        
    def test_toggleIncr(self):
        """checks toggleIncr button
        """
        oppanel  = OptionsPanel(title="My Options Panel")
        oppanel.displayPanel(create=1)
        oppanel.idf.entryByName['togIncr']['widget'].invoke()
        oppanel.master.update()
        pause()
        self.assertEqual(oppanel.incr_entry.cget('state'),'normal')
        oppanel.master.destroy()
        
   
    def test_Continuous(self):
        """checks continuous setting on,off
        """
        oppanel = OptionsPanel(title="My Options Panel")
        oppanel.displayPanel(create=1)
        pause()
        oppanel.idf.entryByName['togCont']['widget'].setitems(items=('on','off'), index=1)
        oppanel.master.update()
        pause()
        self.assertEqual(oppanel.idf.entryByName['togCont']['widget'].getvalue(),'off')
        
        oppanel.idf.entryByName['togCont']['widget'].setitems(items=('on','off'), index=0)
        oppanel.master.destroy()    
    
    def test_showLabel(self):
        """checks show label setting all 'never', 'always', 'move'
        """
        oppanel = OptionsPanel(title="My Options Panel")
        oppanel.displayPanel(create=1)
        pause()
        oppanel.lab_entry.setitems(('never', 'always', 'move'), index=0)
        oppanel.master.update()
        pause()
        self.assertEqual(oppanel.idf.entryByName['togLabel']['widget'].getvalue(),'never')
        oppanel.lab_entry.setitems(('never', 'always', 'move'), index=1)
        oppanel.master.update()
        pause()
        self.assertEqual(oppanel.idf.entryByName['togLabel']['widget'].getvalue(),'always')
        oppanel.lab_entry.setitems(('never', 'always', 'move'), index=2)
        oppanel.master.update()
        pause()
        self.assertEqual(oppanel.idf.entryByName['togLabel']['widget'].getvalue(),'move')   
        oppanel.master.destroy()
        
    def test_type(self):
        """checks setting type int,float
        """
        oppanel = OptionsPanel(title="My Options Panel")
        oppanel.displayPanel(create=1)
        oppanel.master.update()
        pause()
        oppanel.idf.entryByName['togIntFloat']['widget'].setitems(items=('float','int'), index=1)
        oppanel.master.update()
        pause()
        self.assertEqual(oppanel.idf.entryByName['togIntFloat']['widget'].getvalue(),'int')
        oppanel.idf.entryByName['togIntFloat']['widget'].setitems(items=('float','int'), index=0)
        oppanel.master.update()
        pause()
        self.assertEqual(oppanel.idf.entryByName['togIntFloat']['widget'].getvalue(),'float')
        oppanel.master.destroy()


    def test_Precision(self):
        """checks setting precision 1 to 10
        """
        oppanel = OptionsPanel(title="My Options Panel")
        oppanel.displayPanel(create=1)
        items = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10']
        for i in range(10):
            oppanel.idf.entryByName['selPrec']['widget'].setitems(items=('1', '2', '3', '4', '5', '6', '7', '8', '9', '10'),index=i)
            oppanel.master.update()
            pause()
            self.assertEqual(oppanel.idf.entryByName['selPrec']['widget'].getvalue(),items[i])
        oppanel.master.destroy()
   
    def test_setTitle(self):
        """checks setting title
        """
        oppanel = OptionsPanel(title="My Options Panel")
        oppanel.displayPanel(create=1)
        oppanel.setTitle("mypanel")
        self.assertEqual(oppanel.title ,'mypanel')
        oppanel.master.destroy()

    def test_setSensitivity(self):
        """checks setting sensitivity
        """
        oppanel = OptionsPanel(title="My Options Panel")
        oppanel.displayPanel(create=1)
        oppanel.idf.entryByName['inpSens']['wcfg']['textvariable'].set(300)
        oppanel.master.update()
        pause()
        self.assertEqual(oppanel.idf.entryByName['inpSens']['widget'].get(),'300')
        oppanel.master.destroy()



if __name__ == '__main__':
   unittest.main()
