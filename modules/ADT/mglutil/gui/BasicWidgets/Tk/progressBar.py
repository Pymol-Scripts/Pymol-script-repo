#########################################################################
#
# Date: Apr 2003 Authors: Daniel Stoffler
#
#    stoffler@scripps.edu
#
# Copyright: Daniel Stoffler and TSRI
#
#########################################################################

import Tkinter

class ProgressBarConf:
    """callable object to call the progress bar configure method. Used to
    overwrite methods of class instances
    (example: the MolKit/molecule.Molecule)"""
    
    def __init__(self, pb):
        self.pb = pb

    def __call__(self, **kw):
        apply(self.pb.configure, (), kw)
        

class ProgressBarUpd:
    """callable object to call the progress bar set method. Used to
    overwrite methods of class instances
    (example: the MolKit/molecule.Molecule)"""
    
    def __init__(self, pb):
        self.pb = pb

    def __call__(self, value=None):
        self.pb.set(value)
        

class ProgressBar:
    """This class provides a TKinter progress bar widget.

SUMMARY:
   To use this widget, use the configure() method and the set() method.
   use show() to display the widget, use hide() to hide the widget
   .
   Example:
       # create an instance of the progress bar widh an additional label,
       # the label is packed above (top) of the bar
       mybar = ProgressBar(master=None, labelside='top')
       # configure the progress bar with my own width and height, set it to
       # percent mode, initialize it to 0 (init=1), and add a label
       # which says 'Read...'
       mybar.configure(height=150, width=18, init=1, mode='percent',
                       'labeltext'='Read...')
       mybar.set(20) # this sets the bar to 20%
       # now change the mode
       mybar.configure(mode='absolute', max=200)
       mybar.set(100) # this sets the bar to 50%, because the mas is now 200

INSTANCIATION:
    simplest case: instanciate it without any options.
    Example: mybar = ProgressBar()
    Additional init options:
        master: the Tkinter master of this widget. If master is None, the
                master will be set to Tkinter.Toplevel()
        height: the height of the progress bar (not of the entire widget)
                in pixel
        width:  the width of the progress bar (not the entire widget)
        mode:   the widget runs in 3 different modes: 'percent', 'absolute',
                and 'increment'. Explanation below.
        max:    the maximum value. Explanation below.
        labelside: can be 'left', 'right', 'top', 'bottom' or None.
                Used to speicfy where to place the label. None will place
                the label inside the progress bar.

THE CONFIGURE METHOD:
    The progress bar can be configured with many options.
    calling the configure without any arguments returns a dict with the current
    widget configuration.
    Usage: mybar.configure(key1=value1, key2=value2, ...)
    options are:
        init: can be 0 or 1. If set to 1, the progress bar is reset to 0
        mode: 'percent' runs the progress bar in percent mode, i.e. setting
              the value to 0 means 0% and setting it to 100 means 100%
              'absolute' means the widget's set() method computes the current
              percent value based on the max value.
              when in 'increment' mode, the progress bar automatically
              increments its value by 1 and computes the corresponding percent
              value.
        max: the max value is used to compute the current percent value
             (i.e., setting max to 200 and then set the progress bar to 100
              will set it tp 50%)
        progressformat: can be either 'percent' or 'ratio'. 'Percent' displays
                        the current value as 'xxx%', 'ratio' displays the
                        current value as 'x/x, like 2/7'
        granularity: an integer describing at how many % the bar shall update.
                     Default value is 1 (i.e. every 1% the bar redraws)
                     This is onlu usefull when in 'increment' mode.
        width: see above
        height: see above
        labelformat: see above

THE SET METHOD
    The progress bar has a set method to set the current value,
    depending on the mode (see above) the proper percent value is computed,
    and then the bar is redrawn if necessary
    usage: mybar.set(value)

"""
    def __init__(self, master=None, width=150, height=30, max=None,
                 mode='percent', labelside='top'):

        if master is None:
            master = Tkinter.Toplevel()
            master.protocol('WM_DELETE_WINDOW', self.hide )
        self.master = master           # master of this widget

        self.mode = mode               # can be 'percent' or 'absolute'
        self.ONOFF = 'on'              # can be on or off. When off, the
                                       # set method immediately returns

        self.width = width             # total width of widget
        self.height = height           # total height of widget
        self.borderwidth = 1           # the canvas borderwidth
        self.frame = None              # the frame holding canvas & label
        self.canvas = None             # the canvas holding the progress bar
        self.progressBar = None        # the progress bar geometry
        self.backgroundBar = None      # the background bar geometry
        self.progressLabel = None      # the '%' label
        self.progressformat = 'percent'# the format of the progress label,
                                       # can be 'percent' or 'ratio'
        self.label = None              # the text label
        self.labeltext = None          # the text string of the label
        self.labelside = labelside     # where the text label shall be packed
        
        self.progress = 0              # the progress value

        if max is None:
            max = 100
        self.max = max                 # the max value used when mode='absol.'
        self.increment = None          # increment used when mode = 'incr.'
        self.incrementCounter = 0      # counter used for increment

        self.granularity = 1           # update progress bar every 1%
        
        self.createWidget()            # create the widget
        self.frame.pack()              # pack the frame
        self.set(0.0)                  # initialize to 0%


    def createWidget(self):
        self.frame = Tkinter.Frame(self.master, borderwidth=self.borderwidth,
                                   relief='sunken')
        self.canvas = Tkinter.Canvas(self.frame)

        # create the background bar geometry
        bw = self.borderwidth
        self.backgroundBar = self.canvas.create_rectangle( #x0 y0 x1 y1
                0, 0, self.width, self.height, fill='darkblue')
        # create the progress bar geometry
        self.progressBar = self.canvas.create_rectangle( #x0 y0 x1 y1
                0, 0, self.width, self.height, fill='blue')
        # set overall height and with of the widget
        self.setWidth()
        self.setHeight()
        # create the progress label geometry
        self.progressLabel = self.canvas.create_text(
            int(self.width/2), int((self.height/2)+1), text="0%")
        self.canvas.itemconfig(self.progressLabel, fill='white')
        # create the text label geometry
        self.label = Tkinter.Label(self.frame, text='', width=20)
        if self.labelside is not None:
            self.label.pack(side=self.labelside)
        self.canvas.pack()
        

    def setWidth(self, width=None):
        # set overall widget width
        if width is not None:
            self.width = width
        self.canvas.configure(width=self.width) # set width
        self.canvas.coords(self.backgroundBar,  # set background bar
                           0, 0, self.width, self.height)
        self.setBar() # update progress bar


    def setHeight(self, height=None):
        # set overall widget height
        if height is not None:
            self.height = height
        self.canvas.configure(height=self.height)
        self.canvas.coords(self.backgroundBar,  # set background bar
                           0, 0, self.width, self.height)
        self.setBar() # update progress bar
            

    def set(self, value=None):
        # set the progress bar to value
        # Note: value can be None, which can be true when the mode is increment
        # if the mode is 'percent', the value represents the percentage,
        # the values must range from 0 100
        # if the mode is 'absolute', the percentage has to be computed
        # based on self.max
        # if the mode is 'increment', set() is only executed when the new
        # value has reached another full percent

        if self.ONOFF == 'off': # no need to set and redraw if hidden
            return
        
        # test for mode
        if self.mode == 'percent':
            if value is None:
                return
            else:
                self.progress = value
                self.setBar()
                return
        
        elif self.mode == 'absolute':
            if value is None:
                return
            else:
                if self.max == 0:
                    self.progress = 100
                    self.setBar()
                else:
                    self.progress = value * 100.0 / self.max
                    self.setBar()
                return
        
        elif self.mode == 'increment':
            # note that the value is not used here!
            self.incrementCounter = self.incrementCounter + 1
            if self.max == 0:
                progress = 100
            else:
                progress = (100.*self.incrementCounter)/self.max
            if (progress-self.progress)> self.granularity:
                self.progress = progress
                self.setBar()
            return
##             if self.incrementCounter >= self.increment:
##                 if self.max == 0:
##                     self.progress = 100
##                     self.incrementCounter = 0
##                     self.setBar()
##                     return
##                 else:
##                     self.progress = self.progress + self.granularity +\
##                      ( (self.incrementCounter-self.increment)*100.0/self.max )
##                     self.incrementCounter = 0
##                     self.setBar()
##                     return
##             else:
##                 #print '.............',self.progress
##                 return

            
    def setBar(self):
        self.canvas.coords(self.progressBar,
                           0, 0, self.width * self.progress/100.0, self.height)

        # set the label
        if self.progressLabel:
            if self.progressformat == 'percent':
                #print self.progress
                text = str(int(self.progress))+'%'
                #print text
                
            elif self.progressformat == 'ratio':
                if self.mode == 'absolute' or self.mode == 'increment':
                    text = str(int(self.progress*self.max/100.0))+\
                           '/'+str(self.max)
                    #print text, self.progress
                elif self.mode == 'percent':
                    text = str(self.progress)+'/'+str(self.max)

            text = self.label.cget('text')+' '+text
            self.canvas.itemconfig(self.progressLabel,
                                   text=text)
            self.canvas.coords(self.progressLabel, int(self.width/2),
                               int(self.height/2)+1)
        # update idletasks
        self.canvas.update_idletasks()


    def get(self):
        return self.progress
    

    def reset(self):
        # initialize the bar
        oldmode = self.mode
        self.mode = 'percent'
        self.set(0)   # reset progress bar to 0%
        self.mode = oldmode
        if self.mode == 'increment':
            self.setIncrement()    # reset increment counter to 0
        

    def setMax(self, max):
        if max is None:
            if mode == 'percent':
                max = 100
            else:
                max = self.width # reset max
        self.max = max
        self.setIncrement()
        

    def setIncrement(self):
        self.increment = ( float(self.max) * self.granularity ) / 100.0
        self.incrementCounter = 0 # reset increment counter to 0


    def setMode(self, mode):
        if mode not in ['percent', 'absolute', 'increment']:
            mode = 'percent'
            print 'PROGRESSBAR Warning: illegal mode, set to "percent"'
        self.mode = mode
        if mode == 'increment':
            self.setIncrement() # this needs to be reset


    def setLabelText(self, label):
        if label is None:
            return
        self.label.configure(text=label)
        self.labeltext = label
        

    def setProgressFormat(self, format):
        if format not in ['percent','ratio']:
            format = 'percent'
            print 'PROGRESSBAR Warning: illegal progressformant, '+ \
                  'set to "percent"'
        self.progressformat = format

        
    def setGranularity(self, value):
        self.granularity = value
        self.setIncrement()
        

    def show(self):
        if isinstance(self.master, Tkinter.Toplevel):
            self.master.deiconify()
        else:
            self.frame.pack(side='left', pady=3)
            
        self.ONOFF = 'on'


    def hide(self):
        if isinstance(self.master, Tkinter.Toplevel):
            self.master.withdraw()
        else:
            self.frame.forget()
        self.ONOFF = 'off'
        
        
    def configure(self, **kw):
        if len(kw)==0: # return a dict of the current configuration
            cfg = {}
            cfg['width'] = self.width
            cfg['height'] = self.height
            cfg['max'] = self.max
            cfg['mode'] = self.mode
            cfg['labeltext'] = self.labeltext
            cfg['labelside'] = self.labelside
            cfg['progressformat'] = self.progressformat
            cfg['granularity'] = self.granularity
            return cfg

        else: # do a configure
            mode = None
            init = 0

            for key,value in kw.items():
                if key=='width':
                    self.setWidth(value)
                elif key=='height':
                    self.setHeight(value)
                elif key=='max':
                    self.setMax(value)
                elif key=='mode':
                    mode = value
                elif key=='labeltext':
                    self.setLabelText(value)
                elif key=='init':
                    init = value
                elif key=='progressformat':
                    self.setProgressFormat(value)
                elif key=='granularity':
                    self.setGranularity(value)

            # some things have to be set in a particular order
            if mode: # this needs to be called after setting the max
                self.setMode(mode)
            if init: # this needs to be called last
                self.reset()
            

###############################
if __name__ == '__main__':
    import time

    bar = ProgressBar(width=150, height=18)

    print 'type "loop()" and "loop2(max=value)" to run the progress bar'
    
    def loop(sleep=None):
        if sleep is None:
            sleep = 0.001
        bar.configure(mode='percent', labeltext='Read...',progressformat='ratio')
        for i in range(101):
            bar.set(i)
            if i == 33:
                bar.configure(labeltext='Parse...')
            elif i == 66:
                bar.configure(labeltext='Compute...')
            elif i == 100:
                bar.configure(labeltext='Done')
                
            time.sleep(sleep)


    def loop2(sleep=None, max=None):
        if sleep is None:
            sleep = 0.000001

        if max is None:
            max = 50
            
        bar.configure(init=1, mode='increment', max=max,
                      labeltext='Write...', progressformat='ratio')

        for i in range(max):
            bar.set()
            time.sleep(sleep)

        #bar.configure(labeltext='Done...')
    loop(0.1)
