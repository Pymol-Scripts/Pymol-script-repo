#########################################################################
#
# Date: Jun 2002 Author: Daniel Stoffler
#
#    stoffler@scripps.edu
#
# Copyright: Daniel Stoffler and TSRI
#
##########################################################################

import Tkinter
import Pmw
import string, time, os

# class ToolBarButton modified from 'Python and Tkinter Programming'
# handbook by John E. Grayson
class ToolBarButton(Tkinter.Label):
    """
    This class implements a toolbar button.
    - balloonmaster is the master the balloon help is bound to
    - master is the frame the toolbar button is added to
    - name is the name of the toolbarb button
    the toolbar button can have None, one or two icons which will represent
    the button in its disabled (icon1) and enabled (icon2) form.
    - command is the command bound to the button
    - balloonhelp is a text string to be presented as balloon help if the
      application supports balloon help
    - statushelp is a text string to be written in the status bar of the
      application if available
    - height is the height of the button
    - width is the width of the button
    - bd is the border in pixels of the raised button
    - activebackground is the background of the raised button
    - padx and pady are for padding the buddon
    - state can be 'normal' or 'disabled'
    - bg is the background color of the button
    - iconpath is a path to a directory
The master (frame) will get a dictionary 'toolbarButtonDict' which will store
the buttons. key is the button name, value is the button.  Thus, button
names need to be different.

    """

    def __init__(self, balloonmaster=None, master=None, name=None, icon1=None,
                 icon2=None, command=None, balloonhelp='',
                 statushelp='', height=20, width=21,
                 bd=1, activebackground='lightgrey', padx=0, pady=0,
                 state='normal', bg=None, iconpath=None):

        assert master is not None

        # this list stores the buttons
        if not hasattr(master, 'toolbarButtonDict'):
            master.toolbarButtonDict = {}
        # assert the button name doesn't exist already
        assert not master.toolbarButtonDict.has_key(name), name

        self.activebackground = activebackground
        self.name = name
        self.icommand = command
        self.command  = self.activate
        self.buttonPressed = 0  # 0 if button not pressed, 1 if pressed
        self.buttonFocus = 0    # 0 if cursor not over button, 1 if over button
        
        Tkinter.Label.__init__(self, master, height=height, width=width,
                       relief='flat', bd=bd, bg=bg)

        if iconpath is None:
            from mglutil.util.packageFilePath import findFilePath
            iconpath = findFilePath('Tk','mglutil.gui.BasicWidgets')
            iconpath = os.path.join(iconpath,'icons' )
            #iconpath = 'mglutil.gui.BasicWidgets.Tk.icons'
        if icon1 is not None:
            ICONPATH1 = os.path.join(iconpath, icon1)
            if string.splitfields(icon1, '.')[1] == 'bmp':
                self.icon1 = Tkinter.BitmapImage(file=ICONPATH1, master=master)
            else:
                self.icon1 = Tkinter.PhotoImage(file=ICONPATH1, master=master)
        else:
            self.icon1 = None

        if icon2 is not None:
            ICONPATH2 = os.path.join(iconpath, icon2)
            if string.splitfields(icon2, '.')[1] == 'bmp':
                self.icon2 = Tkinter.BitmapImage(file=ICONPATH2, master=master)
            else:
                self.icon2 = Tkinter.PhotoImage(file=ICONPATH2, master=master)
        else:
            self.icon2 = self.icon1

        self.config(image=self.icon1)

        # to prevent binding of the balloon overriding the binding of
        # Leave and Enter events, first the balloon is bound

        if balloonmaster is None:
            master.balloons = Pmw.Balloon(master, yoffset=0)
            balloonmaster = master

        if balloonhelp or statushelp:
            if hasattr(balloonmaster, 'balloons'):
                balloonmaster.balloons.bind(self, balloonhelp, statushelp)

        # a little trick: the '+' adds this binding, otherwise it might
        # be overwritten
        self.bind("<Enter>",           self.buttonEnter, '+')
        self.bind("<Leave>",           self.buttonLeave, '+')
        self.bind("<ButtonPress-1>",   self.buttonDown)
        self.bind("<ButtonRelease-1>", self.buttonUp)
        self.pack(side='left', anchor=Tkinter.NW, padx=padx, pady=pady)
        self.state = state    

        # add the button to the list stored in the master
        master.toolbarButtonDict[name] = self

        if bg is None:
            self.bg = self.configure()['background'][-1]
        else:
            self.bg = bg

      
    def activate(self):
        self.icommand(self.name)


    def buttonEnter(self, event):
        if self.state != 'disabled':
            self.buttonFocus = 1
            if self.buttonPressed == 1:
                self.config(relief='sunken', bg=self.bg)
            else:
                self.config(relief='raised', bg=self.bg)


    def buttonLeave(self, event):
        if self.state != 'disabled':
            self.config(relief='flat', bg=self.bg)
            self.buttonFocus = 0
        

    def buttonDown(self, event):
        if self.state != 'disabled':
            self.config(relief='sunken', bg=self.activebackground)
            self.buttonPressed = 1


    def buttonUp(self, event):
        if self.state != 'disabled':
            if self.command != None and self.buttonFocus:
                self.command()
            time.sleep(0.05)
            self.config(relief='flat', bg=self.bg)  
            self.buttonPressed = 0

    def enable(self):
        self.state = 'normal'
        self.config(image=self.icon2)


    def disable(self):
        self.state = 'disabled'
        self.config(image=self.icon1)


if __name__ == '__main__':
    
    tbframe = Tkinter.Frame()
    tbframe.balloons = Pmw.Balloon(tbframe)

    buttonFuncs = {}

    def foo():
        print 'You pressed a button.'

    def selectFunc(name):
        curFunc = buttonFuncs[name]
        if curFunc:
            curFunc()

    # add separator icon
    ToolBarButton(tbframe, tbframe, name='sep1', icon1='sep.gif', width=10,
                  state='disabled')
    # add 5 smilies
    for name, icon, func, balloon in [
        ('smiley1', 'smiley.gif', foo, 'This is icon1'),
        ('smiley2', 'smiley.gif', foo, 'This is icon2'),
        ('smiley3', 'smiley.gif', foo, 'This is icon3'),
        ('smiley4', 'smiley.gif', foo, 'This is icon4'),
        ('smiley5', 'smiley.gif', foo, 'This is icon5'),
        ]:
    
        ToolBarButton(tbframe, tbframe, name=name, icon1=icon,
                      command=selectFunc, balloonhelp=balloon)
        buttonFuncs[name] = func
        

    # add separator icon
    ToolBarButton(tbframe, tbframe, name='sep2', icon1='sep.gif', width=10,
                  state='disabled')
  

    # pack at the end for performance issue
    tbframe.pack()
