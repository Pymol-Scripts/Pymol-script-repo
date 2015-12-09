#This script intends to make a slider gui for Pymol such that you can quickly grab it and slide through frames
#By Matthew Baumgartner
#University of Pittsburgh
#Laboratory of Carlos Camacho
#Department of Computational and Systems Biology
#mpb21@pitt.edu
#3-17-15 #woo! My birthday!

#It's pretty self explanatory, but here is the wiki page.
#http://www.pymolwiki.org/index.php/Frame_slider


#The MIT License (MIT)
#
#Copyright (c) 2015 Matthew Baumgartner 
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in
#all copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#THE SOFTWARE.

from __future__ import print_function

from pymol import cmd

try:
    import Tkinter
    from Tkinter import *
except ImportError:
    import tkinter as Tkinter
    from tkinter import *

DEBUG = 2

VERSION = 'v0.1'


class FrameSlider:
    ''' Main Pymol Plugin Class '''
    def __init__(self, app):

        parent = app.root
        self.parent = parent
        
        self.app = app
        

        #make a little window for the slider
        self.master = Tkinter.Tk()
        self.master.minsize(width=400, height=1)
        self.master.title("Frame Slider " + VERSION)


        #variable for the current frame that the slider is set to
        self.frame = IntVar()
        self.frame.set(1)
        

        self.max_frame = IntVar()
        self.max_frame.set(self.get_max_enabled_frame())


        #main slider object
        self.slider = Tkinter.Scale(self.master, 
                                    from_ = 1, 
                                    to = self.max_frame.get(), 
                                    orient=HORIZONTAL, 
                                    variable = self.frame, 
                                    command = self.updateframe, 
                                    #label = 'Frame', 
                                    length = 100,
                                    showvalue=False)
        
        self.slider.pack(side=LEFT, fill=BOTH, expand=YES) #
        
        #Make it so when you click on the slider, it rechecks the max number of frames
        #There may be an event that you can listen for when an object is clicked or unclicked in the object pane, but I haven't looked
        self.slider.bind("<Button-1>", self.update_max_frames)
        
        #In order for us to misuse the validate command of a Tkinter.Entry, we need to 'register' the function first. 
        #My first guess was to register it with the master object and it works, but this might not be the right thing to do.
        update_frame = self.master.register(self.update_frame_field)
        
        #make a field that you can type the frame into
        self.frame_entry = Tkinter.Entry(self.master, 
                                           textvariable = self.frame, 
                                           bg = 'black', 
                                           fg = 'white', 
                                           width = 10,
                                           validate = 'key',
                                           validatecommand = (update_frame, '%P'), #'%P' means update on any key press
                                           insertofftime = 200,
                                           insertontime = 200,
                                           insertborderwidth = 5
                                           )

        self.frame_entry.pack(side=LEFT, fill=BOTH)

    
    
    def update_frame_field(self, string):
        ''' 
        Check to see if what was put in the frame field is a int
        Also bastardize this validation code to also change the frame. 
        Entry doesn't have a damn 'command' attribute, so I can't figure out a 
        better way to do this.  
        
        '''
        
        if DEBUG > 3:
            print('Validating Text entry:', string)
        
        #if the string is blank, say that it's valid, so that you can delete out the frame
        if string == '':
            return True
        
        #say no if there is a space in it
        if ' ' in string:
            return False
        
        #say no if you can't convert it to a string
        try:
            f = int(string)
        except TypeError:
            if DEBUG > 3:
                print('Not valid')
            return False
        
        except ValueError:
            #delete the last charchter
            if DEBUG > 2:
                print('Validation found an empty string')
            return False
        
        
        #tell pymol to change frames and update the position of the slider
        self.updateframe(f)
        
        return True
         

    def updateframe(self, f):
        ''' When the slider is moved, run this.  ''' 

        #Check the max number of frames
        self.update_max_frames()
        
        #change to that frame
        cmd.frame(f)
    
    def update_max_frames(self, arg = None):
        ''' Call the couple of functions required to update the max frames 
        Any arguments are ignored.
        '''
        #figure out what the max frame is and save it as a variable
        self.max_frame.set(self.get_max_enabled_frame())
        #update the max position of the slider
        self.slider.config(to = self.max_frame.get())
    
    def get_max_enabled_frame(self):
        ''' Go through the enabled objects in Pymol and get the max frame'''
        fr_max = 1
        
        for obj in cmd.get_names(enabled_only=1):
            frms = cmd.count_states(obj)
            
            if DEBUG > 4:
                print('obj:', obj, 'has', frms, 'frames')
                
            #if it is higher, set it as the new max
            if frms >  fr_max:
                fr_max = frms
        if DEBUG > 3:
            print("Max frames:", fr_max)
            
        return fr_max
            
    

def __init__(self):
    self.menuBar.addmenuitem('Plugin', 'command', 'Frame Slider', label = 'Frame Slider', command = lambda s=self : FrameSlider(s))  
        
if __name__ == '__main__':
    print('Error: This should be installed as a Pymol Plugin. Plugins > Install...')
    sys.exit()

