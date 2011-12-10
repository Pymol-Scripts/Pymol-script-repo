# This module handles Splash Screen
# $Header: /opt/cvs/python/packages/share1.5/mglutil/splashregister/splashscreen.py,v 1.22 2008/04/16 18:19:20 annao Exp $
# $Date: 2008/04/16 18:19:20 $
# $Id: splashscreen.py,v 1.22 2008/04/16 18:19:20 annao Exp $
from mglutil.gui.BasicWidgets.Tk.progressBar import ProgressBar
from mglutil.util.packageFilePath import getResourceFolderWithVersion, findFilePath
from mglutil.util.misc import ensureFontCase
import Tkinter, os, sys, time

class SplashScreen:
    """
    package    : mglutil
    module     : splashregister.splashscreen
    class      : SplashScreen
    description:
        Provides splash screen
    """
    def __init__(self, about, noSplash=False):
        """Constructor for SplashScreen"""
        self.version = about.version
        rcWithVersion = getResourceFolderWithVersion()
        if rcWithVersion is not None:
            registration = rcWithVersion + os.sep + '.registration'
            self.timing_rc = rcWithVersion + os.sep + '.timing_' + about.title
            self.usage_rc = rcWithVersion + os.sep + '.usage_' + about.title
            open(self.usage_rc,'a').write(str(time.time())+"\n")
            n_usage = os.path.getsize(self.usage_rc)
        else:
            registration = None
            self.timing_rc = None
            self.usage_rc = None

        self.noSplash = noSplash
        self.registered = False # change this to False to check registration
        self.register = None
        
        # check if the user has registered
        if registration is None or os.path.exists(registration):
            self.registered = True
            #self.check_for_updates(registration)        
        else:
            # count number of use and ask to register
            lines = open(self.usage_rc).readlines()
            if len(lines) < 10:
                self.registered = True
        
        self.waitTk = Tkinter.IntVar()
        self.percent = 0
        if self.timing_rc is None:
            self.timing_data = []
        else:
            
            try:
                timing_file = open(self.timing_rc,'r+')
            except:
                timing_file = open(self.timing_rc,'w+')
            self.timing_data = timing_file.readlines()
            timing_file.close()
            ###self.percent = 0
            if self.timing_data:
                for data in self.timing_data:
                    self.percent += int(data)
                self.percent = int( self.percent/len(self.timing_data) ) 
                self.percent = 100 / self.percent        
        if self.percent == 0: self.percent = 1    
        self.counter = 0

        if self.registered and self.noSplash:
            return
        
        self.splash_win = Tkinter.Toplevel()
        self.splash_win.overrideredirect(1)
        self.splash_win.withdraw()        
 
        frame = Tkinter.Frame(self.splash_win,relief='raised', bd=3)
        frame.pack()
        try:
            about.gui(master=frame)
        except Exception, inst:
            print inst
        self.progressBar = ProgressBar(master=frame, labelside=None,
                                width=420, height=20, mode='percent')
        self.progressBar.setLabelText('Loading Modules...')
        self.progressBar.set(0)

        if not self.registered:
            text = """ Please Register! It helps us secure funding for 
supporting development and you won't have to 
click these buttons again in the future. Thanks."""
            Tkinter.Label(frame, text = text, font =(ensureFontCase('helvetica'), 14, 'bold') ).\
                                                                          pack()
            Tkinter.Button(frame, text='Register Now', bg = 'Green',
                            command=self.Register_Now).pack(side='left')
            self.Later_Button = Tkinter.Button(frame, text='Remind Me Later', 
                       state='disabled', command=self.Later )
            self.Later_Button.pack(side='right')
      
        self._updateprogressBar()
        self.splash_win.update_idletasks()

        width = self.splash_win.winfo_reqwidth()
        height = self.splash_win.winfo_reqheight()
        screenwidth = self.splash_win.winfo_screenwidth()
        if screenwidth > 2000:
            screenwidth = 1000
        x = (screenwidth - width) / 2 - self.splash_win.winfo_vrootx()
        y = (self.splash_win.winfo_screenheight() - height) / 4 - \
                                                    self.splash_win.winfo_vrooty()
        if x < 0:
            x = 0
        if y < 0:
            y = 0
        geometry = '%dx%d+%d+%d' % (width, height, x, y)
        self.splash_win.geometry(geometry)
        self.splash_win.update_idletasks()
        self.splash_win.deiconify()
        self.splash_win.update()
        self.splash_win_children = self.splash_win.children
        self.splash_win.children = [] # this is needed to avoid problems self.changeFont

    def  _updateprogressBar(self):
        """Updates Progress Bar"""
        abs_percent =  (self.counter + 1)*self.percent
        self.counter += 1 
        if abs_percent < 100:
            self.progressBar.set(abs_percent)
            self.splash_win.after(1, lambda:self._updateprogressBar())
        else:
            self.progressBar.setLabelText('Please wait ')        
            self.progressBar.set(100)

    def finish(self):
        """Used to remove Splash Screen"""
        if self.registered and self.noSplash:
            return
        if not self.registered:
            if self.splash_win.wm_state() == 'normal':
                self.Later_Button.config(state='normal', bg='red')
                self.splash_win.update()
                self.splash_win.wait_variable(self.waitTk)

        self.splash_win.children = self.splash_win_children
        self.splash_win.destroy()

        if len(self.timing_data) > 10:
            self.timing_data = self.timing_data[1:10]
        self.timing_data.append(str(self.counter)+'\n')
        try:
            timing_file = open(self.timing_rc,'w') 
            timing_file.writelines(self.timing_data)
            timing_file.close()
        except:
            pass


    def Later(self):
        """Command for self.Later_Button"""
        self.waitTk.set(1)
        
    def Register_Now(self):   
        if not self.register or not self.register.master:            
            from mglutil.splashregister.register import Register_User
            self.splash_win.withdraw()
            self.register = Register_User(self.version)
            self.waitTk.set(1)
            
    def check_for_updates(self,registration):
        """To be implemented """
        self.registered = True
        #from mglutil.splashregister.register import Update_User
        #Update_User(registration)
