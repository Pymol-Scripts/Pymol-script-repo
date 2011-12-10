#############################################################################
#
# Author: Ruth HUEY
#
# Copyright: M. Sanner TSRI 2003
#
#############################################################################

#
#
# $Header: /opt/cvs/python/packages/share1.5/mglutil/gui/BasicWidgets/Tk/player.py,v 1.18.4.1 2011/03/24 23:42:50 annao Exp $
#
# $Id: player.py,v 1.18.4.1 2011/03/24 23:42:50 annao Exp $
#
#
#
#

import Tkinter, Pmw
import types, time, os
from mglutil.util.callback import CallBackFunction
from mglutil.gui.InputForm.Tk.gui import InputFormDescr,InputForm,evalString
from mglutil.util.callback import CallbackManager
from mglutil.gui.BasicWidgets.Tk.thumbwheel import ThumbWheel
from mglutil.util.packageFilePath import findFilePath
from mglutil.util.misc import ensureFontCase
import tkMessageBox

try:
    import pymedia
    pymediaFound = True
except:
    pymediaFound = False


class Player(Tkinter.Frame):
    """Widget to play sequence of frames.
    The GUI allows to specify direction and speed of playback.

    the nextFrame(number) methoid should be overridden and define how
    this player plays the animation
    
    the root constructor argument can be a Tkinter container in which to
    embde the control panel.
    
    required attributes:
        currentFrameIndex = 0
        startFrame = 0
        endFrame = 0
        maxFrame = 0
        stepSize = 1
        target = endFrame           #determines direction of play
        increment = 1               #implements decrementing vs incrementing
                                    #  in getNextFrameIndex
        playMode = 0                #playMode options:
                                    #   0   play once and stop
                                    #   1   play continuously in 1 direction
                                    #   2   play once in 2 directions
                                    #   3   play continuously in 2 directions
        framerate =15.              # number of frame per second to be display
        hasSlider = False           #adds Tkinter Scale Widget if True

    customization attributes:
        buttonMask is a dictionary that provides a button name as a key and
        a boolean as a value. For each button in the GUI if the name appears
        in self.buttonMask and the value is False the button willnot be
        addded. This only works when form2 is true. slider and counter
        are handle by the hasSlider and coutner keyword arguments

    required methods:
    #play methods:
        play                           #play according to play mode
        getNextFrameIndex              #returns index of next frame to
                                        play according to current index
                                        and playMode
        nextFrame                      #actually display the nextFrame
                                        ??possibly update widgets???

    #methods to regulate play
        Play_cb                        #play according to playMode
        PlayRev_cb                     #play backwards  
        FastForward_cb()               #play foward at max speed
        FastReverse_cb()               #play reverse at max speed
        Stop_cb                        #stop current play
        ?Pause_cb                      #stop at current frame??
        ?Loop_cb                       #???????


    #methods to set Frame to a specific frame
        SetState_cb                    #counter callback for random access
        GoToStart_cb                   #set to current startFrame
        GoToEnd_cb                     #set to current endFrame
        ?setCurrentFrameIndex          #set currentframe index???


    #methods for player gui
        showGUI                        #show the gui
        buildForm2                     #opens image-based gui
        buildForm                      #opens text-based gui
        Close_cb                       #withdraws gui

    #methods for pmw counter
        custom_validate                #used by pmw to check entry
        custom_counter                 #used by pmw for counter


    #methods for changing playMode
        SetMode_cb                     #opens playModeForm
        setPlayMode_cb                 #sets playMode, sets delay,
                                          startFrame, endFrame 
                                          AND closes playModeForm
        cancelPlayMode_cb              #closes playModeForm w/out changes

        #methods for end points:
            setStartFrame()            #set startFrame
            setEndFrame()              #set endFrame
            setStepSize()              #sets increment, def=1
        additional methods:
    """

    def __init__(self, master=None, root=None,
                        height=80,width=200,
                        currentFrameIndex=0,
                        startFrame=0, 
                        endFrame=0,
                        maxFrame=0,
                        stepSize=1, 
                        playMode=0,
                        ##afterDelay=50,
                        titleStr='Player',
                        gotoStartfile = 'go_to_start.gif',
                        gotoEndfile = 'go_to_end.gif',
                        ff_revfile = 'ff_rev.gif',
                        ff_fwdfile = 'ff_fwd.gif',
                        stopfile = 'stop.gif',
                        playfile = 'play_fwd.gif',
                        playRevfile = 'play_rev.gif',
                        chmodfile = 'chmod.gif',
                        closefile = 'close.gif',
                        iconpath = None,
                        counter = 1,
                        form2=1,
                        gui=1,
                        framerate=15.,
                        hasSlider=False,
                        buttonMask=None):

        # frame rate in fps
        self.framerate = framerate 
        self.gui = gui
        self.currentFrameIndex = currentFrameIndex
        self.startFrame = startFrame
        self.endFrame = endFrame
        self.targetFrame = self.endFrame
        self.maxFrame = maxFrame
        if not maxFrame:
            self.maxFrame = endFrame
        self.stepSize = stepSize
        #used for play once in 2 directions
        self.oneDirection = 0
        self.increment = 1
        self.target = self.endFrame
        self.hasSlider = hasSlider
        #used to coordinate switching icons
        #playMode options:
        #   0   play once and stop
        #   1   play continuously in 1 direction
        #   2   play once in 2 directions and stop
        #   3   play continuously in 2 directions
        self.playMode = playMode
        ### replace by framerate
        #amt of time to sleep
        #self.delay = delay


        #amt of time for after
        #self.afterDelay = afterDelay

        self.afterID = None
        self.stop = 1
        self.hasCounter = 0
        # gui variable
        #self.master = master
        #self.root = root
        self.form2 = form2
        self.hasCounter = counter

        if buttonMask is None:
            self.buttonMask = {}
            # self.buttonMask provides a button name as a key and a boolean
            # as a value. For each button in the GUI if the name appears in
            # self.buttonMask and the value is False the button willnot be
            # addded
        else:
            self.buttonMask = buttonMask
            
        if gui:
            self.showGUI(master=master,root=root,
                         titleStr=titleStr,height=height,
                         width=width,iconpath=iconpath)

        # make sure we have a master as we will call master.update in play()
        if hasattr(master, 'update'):
            self.masterForUpdate = self.master
        else:
            self.masterForUpdate = self.root


    def showGUI(self,master=None,root=None,
                width=200,height=80,titleStr='player',
                gotoStartfile = 'go_to_start.gif',
                gotoEndfile = 'go_to_end.gif',
                ff_revfile = 'ff_rev.gif',
                ff_fwdfile = 'ff_fwd.gif',
                stopfile = 'stop.gif',
                playfile = 'play_fwd.gif',
                playRevfile = 'play_rev.gif',
                chmodfile = 'chmod.gif',
                closefile = 'close.gif',iconpath=None):
        """ function to display the player gui."""

        if hasattr(self, 'form'):
	    if hasattr(self.form, "deiconify"):
                self.form.deiconify()
                return   
        self.master=master
        self.root=root

        if not self.form2:
            self.form = self.buildForm(titleStr)  # pass some arguments here
            self.form2=0
        else:
            if iconpath is None:
                iconpath = ('mglutil.gui.BasicWidgets.Tk','icons')
            ICONDIR = findFilePath(iconpath[1], iconpath[0])
            #if findFilePath failed, already returned
            gotoStartfile = os.path.join(ICONDIR,gotoStartfile)
            gotoEndfile = os.path.join(ICONDIR,gotoEndfile)
            ff_revfile = os.path.join(ICONDIR,ff_revfile)
            ff_fwdfile = os.path.join(ICONDIR,ff_fwdfile)
            stopfile = os.path.join(ICONDIR,stopfile)
            playfile = os.path.join(ICONDIR,playfile)
            playRevfile = os.path.join(ICONDIR,playRevfile)
            chmodfile = os.path.join(ICONDIR,chmodfile)
            closefile = os.path.join(ICONDIR,closefile)
            recordfile = os.path.join(ICONDIR,"record.gif")
            record1file = os.path.join(ICONDIR, "record1.gif")
            
            self.gotoStartIcon = Tkinter.PhotoImage(file=gotoStartfile, master=master)
            self.gotoEndIcon = Tkinter.PhotoImage(file=gotoEndfile, master=master)
            self.ff_revIcon = Tkinter.PhotoImage(file=ff_revfile, master=master)
            self.ff_fwdIcon = Tkinter.PhotoImage(file=ff_fwdfile, master=master)
            self.stopIcon = Tkinter.PhotoImage(file=stopfile, master=master)
            self.playIcon = Tkinter.PhotoImage(file=playfile, master=master)
            self.playRevIcon = Tkinter.PhotoImage(file=playRevfile, master=master)
            self.recIcon = Tkinter.PhotoImage(file=recordfile, master=master)
            self.rec1Icon = Tkinter.PhotoImage(file=record1file, master=master)
            self.chmodIcon = Tkinter.PhotoImage(file=chmodfile, master=master)
            self.closeIcon = Tkinter.PhotoImage(file=closefile, master=master)
            self.form = self.buildForm2(titleStr)  # pass some argument here

    #play methods:
    # play, waitTime, getNextFrameIndex, nextFrame
    def play(self, framerate, event=None):
        t1 = time.time() # previous frame time
        timestamp = 1./framerate # rate to update frame
        self.stop = 0
        ind = self.currentFrameIndex

        # if player at the end we resume from begining
        if ind>= self.endFrame:
            self.GoToStart_cb()
            
        while not self.stop: #this has to be more complex
            if self.stop: 
                print 'play stopped!'
                break
            
            #do something different here if ff_fwd or ff_rev
            if self.gui:
                #self.afterID = self.master.after(self.afterDelay, self.waitTime)
                self.masterForUpdate.update()

            t2 = time.time()  # current time
            t = t2 -t1        # time difference between current frame and previous

            if framerate > -1 and t < timestamp:
                pass
            else:
                id = self.getNextFrameIndex(self.currentFrameIndex)
                if id==None:
                    self.stop = 1
                    if self.gui:
                        self.Stop_cb()
                        break
                self.nextFrame(id) #maybe set entry then call nextFrame??
                self.currentFrameIndex = id
                t1 = t2
            
    def getNextFrameIndex(self, index):
        newFrame = index + self.increment*self.stepSize
        # check if  newFrame in current range:
        #       if incrementing, has to be <=self.endFrame
        #       if decrementing, has to be >=self.startFrame
        #   NB: incPos indicates whether currently incrementing or decrementing
        # FORCE newFrame into range 
        incPos = self.increment>0
        if incPos and newFrame>self.endFrame:
            newFrame = self.endFrame
        elif not incPos and newFrame<self.startFrame:
            newFrame = self.startFrame
        # check whether reached current targetFrame
        #       if so, action depends on current playMode 
        mode = self.playMode
        if index==self.targetFrame:
            # playModes 0 and 2 are play once and stop (in 1 or 2 directions)
            if not mode%2:
                if mode==0:
                    # None means stop
                    return None
                else:
                    # play once in 2 directions
                    #   check whether already reached halfway pt
                    if self.oneDirection:
                        return None
                    self.oneDirection = 1
                    # to reverse direction:
                    #   toggle increment
                    self.increment = -1 * self.increment
                    #   toggle targetFrame
                    if self.targetFrame==self.endFrame:
                        self.targetFrame = self.startFrame
                    else:
                        self.targetFrame = self.endFrame
                    return newFrame
            elif self.playMode==1:
                # play continuously in 1 direction
                # toggle targetFrame to the opposite end
                if self.targetFrame==self.endFrame:
                    return self.startFrame
                else:
                    return self.endFrame
            elif self.playMode==3:
                # loop continuously in 2 directions
                #   toggle increment
                self.increment = -1 * self.increment
                newFrame = self.targetFrame + self.stepSize * self.increment
                #   toggle targetFrame
                if self.targetFrame==self.endFrame:
                    self.targetFrame = self.startFrame
                else:
                    self.targetFrame = self.endFrame
                return newFrame
        return newFrame


##     def waitTime(self):
##         self.afterID = None
##         t1 = time.time()
##         delta = time.time() - t1
##         #curBut is either play or playRev button
##         hasCurBut = hasattr(self, 'curBut')
##         while delta < self.delay:
##             if hasCurBut:
##                 self.curBut.config(bg='red')
##             self.master.update()
##             time.sleep(self.delay/100)
##             delta = time.time() - t1
##         if hasCurBut:
##             self.curBut.config(bg='white')
##         if self.afterID is not None and hasCurBut:
##             self.curBut.after_cancel(self.afterID)
##             self.afterID = self.curBut.after(self.afterDelay, self.waitTime)
        

    def nextFrame(self, id):
        ##pass #must be overriden
        id = int(id)
        if id == self.currentFrameIndex: return        
        if self.hasCounter and self.gui:
            self.form.ent2.delete(0,'end')
            self.form.ent2.insert(0, str(id))
            if self.hasSlider:
                self.form.ifd.entryByName['slider']['widget'].set(id)
        print 'playing ', id
        self.currentFrameIndex = int(id)


    #methods to call play
    #Play_cb, PlayRev_cb, FastForward_cb, FastReverse_cb
    #Stop_cb and possibly Pause_cb, Loop_cb
    def startRecording_cb(self, event=None):
        pass
 

    def stopRecording_cb(self, event=None):
        pass

    def Play_cb(self, framerate=None,event=None):
        #print 'Play_cb'
        self.form.ifd.entryByName['playB']['widget'].grab_release()
        #this is a new call to Play_cb
        if framerate == None:
            framerate = self.framerate
        self.stop = 0
        self.oneDirection = 0
        self.targetFrame = self.endFrame
        self.increment = 1
        #print 'currentFrameIndex=', self.currentFrameIndex
        #possibly currently playing reverse
        if hasattr(self,'curBut') and self.curBut==self.form.playRevB:
            if self.form2:
                self.curBut.config(command=self.oldCmd,
                        image=self.oldImage,bg='white')
            else:
                self.curBut.config(command=self.oldCmd, text=self.oldText)
        self.oldTT = self.form.playTT
        self.oldTTtext = 'play forward according to current play mode'
        self.curBut = self.form.playB
        self.oldCmd = self.Play_cb
        self.oldImage = self.playIcon
        self.oldText = 'Play'
        if self.form2:
            self.curBut.config(command=self.Stop_cb, image=self.stopIcon)
        else:
            self.curBut.config(command=self.Stop_cb, text='Stop')
        self.oldTT.bind(self.curBut,'stop play')
        #self.master.update()
        self.play(framerate)
        self.stopRecording_cb()

    def PlayRev_cb(self, framerate=None,event=None):
        #figure out current number and then play backwards
        #print 'currentFrameIndex=', self.currentFrameIndex
        self.form.ifd.entryByName['playRevB']['widget'].grab_release()
        if framerate == None:
            framerate=self.framerate
        self.increment = -1
        self.targetFrame = self.startFrame
        self.oneDirection = 0
        self.stop = 0
        #possibly currently playing 
        if hasattr(self,'curBut') and self.curBut==self.form.playB:
            if self.form2:
                self.curBut.config(command=self.oldCmd, 
                        image=self.oldImage,bg='white')
            else:
                self.curBut.config(command=self.oldCmd, text=self.oldText)
        self.oldTT = self.form.playRevTT
        self.oldTTtext = 'play reverse according to current play mode'
        self.curBut = self.form.playRevB
        self.oldCmd = self.PlayRev_cb
        self.oldImage = self.playRevIcon
        self.oldText = 'Play Reverse'
        if self.form2:
            self.curBut.config(command=self.Stop_cb, image=self.stopIcon)
        else:
            self.curBut.config(command=self.Stop_cb, text='Stop')
        self.oldTT.bind(self.curBut,'stop play')
        self.play(framerate)

    def FastReverse_cb(self, event=None):
        #print 'FastReverse'
        #framerate = self.framerate * 2
        self.oneDirection = 0
        self.PlayRev_cb(framerate=-1)#framerate)
        

    def FastForward_cb(self, event=None):
        #print 'FastForward'
        self.oneDirection = 0
        #framerate = self.framerate * 2
        self.Play_cb(framerate=-1)
        


    def Stop_cb(self, event=None):
        self.stop = 1
        if hasattr(self, 'curBut'):
            if self.form2:
                self.curBut.config(command=self.oldCmd, image=self.oldImage)
            else:
                self.curBut.config(command=self.oldCmd, text=self.oldText)
                #FIX THIS!!!!
                #DECIDE if stop means reset to start or not
                #self.nextFrame(self.startFrame)
            self.oldTT.bind(self.curBut,self.oldTTtext)

        self.stopRecording_cb()

        #FIX THIS: does this go here?  or in nextFrame
        #clean up form
        #if self.hasCounter:
        #    self.form.ent2.delete(0,'end')
        #    self.form.ent2.insert(0, str(self.currentFrameIndex))


    #these may be superfluous: called by original form
    def PlayReturn_cb(self, event=None):
        #print 'PlayReturn_cb'
        #this should be superfluous
        self.playMode = 3
        self.Play_cb()


    def Loop_cb(self, event=None):
        #print 'Loop_cb'
        #this should be superfluous
        self.playMode = 3
        self.Play_cb()


    def Pause_cb(self, event=None):
        self.stop = 1


    #methods to set Frame to a specific frame
    #SetState_cb, GoToStart_cb, GoToEnd_cb, setCurrentFrameIndex
    def SetState_cb(self,  event=None):
        #do nothing if no counter
        if self.hasCounter:
            index = self.form.counter.get()
            self.nextFrame(index)


    def GoToStart_cb(self, event=None):
        #print "GoToStart_cb", self.startFrame, self.currentFrameIndex
        #self.currentFrameIndex = self.startFrame
        self.oneDirection = 0
        self.nextFrame(self.startFrame)


    def GoToEnd_cb(self, event=None):
        #print 'GoToEnd'
        self.oneDirection = 0
        #self.currentFrameIndex = self.endFrame
        self.nextFrame(self.endFrame)


    def setCurrentFrameIndex(self, index):
        #print 'setting currentFrameIndex to', index
        self.currentFrameIndex = index 


    #methods for player gui
    def Close_cb(self, event=None):
        self.stop = 1
        if hasattr(self,'form'):
            self.form.withdraw()
    def SetAnim_cb(self):
        """ function to be overwritten.
        Use to call a functoin to set the animatiom frame."""

        return

    def buildForm2(self, titleStr):
        self.stop = 1
        if hasattr(self, 'form'):
            if hasattr(self.form, 'deiconify'):
                self.form.deiconify()
                return
        maxval = self.endFrame
        ifd = InputFormDescr(title=titleStr)

        if self.buttonMask.get('gotoStartB', None) is not False:
            ifd.append({'name': 'gotoStartB',
                'widgetType': Tkinter.Button,
                #'text':'gotoStart',
                'tooltip':'sets frame to current startFrame',
                'wcfg':{'bd':4,
                        'image':self.gotoStartIcon,
                        'width':self.gotoStartIcon.width(),
                        'height':self.gotoStartIcon.height()
                },
                'gridcfg':{'sticky':'nesw'},
                'command':self.GoToStart_cb})

        if self.buttonMask.get('fastReverseB', None) is not False:
            ifd.append({'name': 'fastReverseB',
                'widgetType': Tkinter.Button,
                #'text':'fastReverse',
                'tooltip':'play reverse as fast as possible',
                'wcfg':{'bd':4,
                        'image':self.ff_revIcon,
                        'width':self.ff_revIcon.width(),
                        'height':self.ff_revIcon.height()
                },
                'gridcfg':{'sticky':'nesw', 'row':-1, 'column':1},
                'command':self.FastReverse_cb})

        if self.buttonMask.get('playRevB', None) is not False:
            ifd.append({'name': 'playRevB',
                'widgetType': Tkinter.Button,
                #'text':'Play Reverse',
                'tooltip':'play reverse according to current play mode',
                'wcfg':{'bd':4,
                        'image':self.playRevIcon,
                        'width':self.playRevIcon.width(),
                        'height':self.playRevIcon.height()},
                'gridcfg':{'sticky':'nesw','row':-1, 'column':2},
                'command':self.PlayRev_cb})

            if self.hasCounter:
                ifd.append({'widgetType':Pmw.Counter,
                        'name':'statesCounter',
                        'required':1,
                        'tooltip':'used to show frames via random access',
                        'wcfg':{#'labelpos': 'n,
                             #'label_text':'conformation:  ',
                            'autorepeat':0,
                            'entryfield_value':self.startFrame,
                            #'entryfield_value':self.idList[0],
                            'datatype': self.custom_counter,
                            'entry_width':3,
                            'entryfield_validate': self.custom_validate },
                         'gridcfg':{'sticky':'nesw', 'row':-1, 'column':3,'columnspan':2}})

        if self.buttonMask.get('playB', None) is not False:
            ifd.append({'name': 'playB',
                'widgetType': Tkinter.Button,
                #'text':'Play',
                'tooltip':'play forward according to current play mode',
                'wcfg':{'bd':4,
                        'image':self.playIcon,
                        'width':self.playIcon.width(),
                        'height':self.playIcon.height()
                },
                'gridcfg':{'sticky':'nesw', 'row':-1, 'column':5},
                'command':self.Play_cb})

        if self.buttonMask.get('fastForwardB', None) is not False:
            ifd.append({'name': 'fastForwardB',
                'widgetType': Tkinter.Button,
                #'text':'fastForward',
                'tooltip':'play as fast as possible',
                'wcfg':{'bd':4,
                        'image': self.ff_fwdIcon,
                        'width':self.ff_fwdIcon.width(),
                        'height':self.ff_fwdIcon.height()
                        },
                'gridcfg':{'sticky':'nesw','row':-1, 'column':6},
                'command':self.FastForward_cb})

        if self.buttonMask.get('gotoEndB', None) is not False:
            ifd.append({'name': 'gotoEndB',
                'widgetType': Tkinter.Button,
                #'text':'gotoEnd',
                'tooltip':'sets frame to current endFrame',
                'wcfg':{'bd':4,
                        'image':self.gotoEndIcon,
                        'width':self.gotoEndIcon.width(),
                        'height':self.gotoEndIcon.height()
                },
                'gridcfg':{'sticky':'nesw','row':-1, 'column':7},
                'command':self.GoToEnd_cb})

        if self.buttonMask.get('modeB', None) is not False:
            ifd.append({'name': 'modeB',
                'widgetType': Tkinter.Button,
                'text':'Change Mode',
                'tooltip':'opens panel to change play options',
                'wcfg':{'bd':4,
                        'image':self.chmodIcon,
                        'width':self.chmodIcon.width(),
                        'height':self.chmodIcon.height()
                        },
                'gridcfg':{'sticky':'nesw','row':-1, 'column':8},
                'command':self.SetMode_cb})

        if pymediaFound and self.buttonMask.get('recordB', None) is not False:
            ifd.append({'name': 'recordB',
                'widgetType': Tkinter.Checkbutton,
                'text':'Record',
                'tooltip':'record an mpeg movie into movie_####.mpg',
                'defaultValue':0,
                'wcfg':{'bd':4,
                        'variable':Tkinter.IntVar(),
                        'image':self.recIcon,
                        'width':self.recIcon.width(),
                        'height':self.recIcon.height(),
                        'indicatoron':0,
                        },
                'gridcfg':{'sticky':'nesw','row':-1, 'column':9},
                'command':self.startRecording_cb})

        if self.buttonMask.get('setanimB', None) is not False:
            ifd.append({'name': 'setanimB',
                'widgetType': Tkinter.Button,
                'text':'SetAnim',
                'tooltip':'Set Animation',
                'wcfg':{'bd':4},
                'gridcfg':{'sticky':'nesw','row':-1, 'column':10},
                'command':self.SetAnim_cb})
        
        if self.buttonMask.get('closeB', None) is not False:
            ifd.append({'name': 'closeB',
                'widgetType': Tkinter.Button,
                'text':'Close',
                'tooltip':'closes player',
                'wcfg':{'bd':4,
                        'image':self.closeIcon,
                        'width':self.closeIcon.width(),
                        'height':self.closeIcon.height(),
                        },
                'gridcfg':{'sticky':'nesw','row':-1, 'column':11},
                #'gridcfg':{'sticky':'nesw', 'columnspan':2},
                'command':self.Close_cb})

        if self.hasSlider:
            ifd.append({'name': 'slider',
                'widgetType': Tkinter.Scale,
                'wcfg':{'orient':'horizontal',
                        'from_':self.startFrame,
                        'to':self.maxFrame,
                        'showvalue':False},
                'gridcfg':{'sticky':'nesw','row':1, 'column':0, 'columnspan':12},
                'command':self.nextFrame
                })
        #form = self.vf.getUserInput(ifd, modal=0,blocking=0)
        form = InputForm(self.master, self.root,
                         descr = ifd,
                         modal = 0, blocking = 0,
                         closeWithWindow=1)
        form.ifd = ifd
        form.playB = form.ifd.entryByName['playB']['widget']
        form.playRevB = form.ifd.entryByName['playRevB']['widget']
        #set up link to balloon help which needs to change, also
        form.playTT = form.ifd.entryByName['playB']['balloon']
        form.playRevTT = form.ifd.entryByName['playRevB']['balloon']
        if self.hasCounter:
            ctr = ifd.entryByName['statesCounter']['widget']
            entF = ctr.component('entryfield')
            form.ent2 = entF._entryFieldEntry
            da = ctr.component('downarrow')
            ua = ctr.component('uparrow')
            for item in [da,ua]:
                item.bind('<ButtonPress-1>', self.SetState_cb, '+')
            form.ent2.bind('<Return>', self.SetState_cb, '+')
            form.counter = form.ifd.entryByName['statesCounter']['widget']
        if self.hasSlider:
            slider = form.ifd.entryByName['slider']['widget']
            slider.set(self.currentFrameIndex)
        #print 'returning form'
        return form


    def buildForm(self, titleStr):
        #??FIX THIS:
        self.stop = 1
        if hasattr(self, 'form'):
            self.form.deiconify()
            return
        maxval = self.endFrame
        ifd = InputFormDescr(title=titleStr)
        if self.hasCounter:
            ifd.append({'widgetType':Pmw.Counter,
                    'name':'statesCounter',
                    'required':1,
                    'tooltip':'used to show frames via random access',
                    'wcfg':{#'labelpos': 'n,
                         #'label_text':'conformation:  ',
                        'autorepeat':0,
                        'entryfield_value':self.startFrame,
                        #'entryfield_value':self.idList[0],
                        'datatype': self.custom_counter,
                        'entry_width':9,
                        'entryfield_validate': self.custom_validate },
                     'gridcfg':{'sticky':'nesw', 'columnspan':2}})
        ifd.append({'name': 'playB',
            'widgetType': Tkinter.Button,
            'text':'Play',
            'tooltip':'play sequence according to current play mode',
            'wcfg':{'bd':4},
            'gridcfg':{'sticky':'nesw','columnspan':1},
            'command':self.Play_cb})
        ifd.append({'name': 'playRevB',
            'widgetType': Tkinter.Button,
            'text':'Play Reverse',
            'wcfg':{'bd':4},
            'gridcfg':{'sticky':'nesw','row':-1, 'column':1},
            'command':self.PlayRev_cb})
        ifd.append({'name': 'playTB',
            'widgetType': Tkinter.Button,
            'text':'Play+Return',
            'wcfg':{'bd':4},
            'gridcfg':{'sticky':'nesw','columnspan':1},
            'command':self.PlayReturn_cb})
        ifd.append({'name': 'loopB',
            'widgetType': Tkinter.Button,
            'text':'Loop',
            'wcfg':{'bd':4},
            'gridcfg':{'sticky':'nesw','row':-1, 'column':1},
            'command':self.Loop_cb})
        ifd.append({'name': 'stopB',
            'widgetType': Tkinter.Button,
            'text':'Stop',
            'tooltip':'stop play',
            'wcfg':{'bd':4,
                    },
            'gridcfg':{'sticky':'nesw'},
            'command':self.Stop_cb})
        #add fastforward, fastrewind, thumbwheel for speed
        ifd.append({'name': 'pauseB',
            'widgetType': Tkinter.Button,
            'text':'Pause',
            'wcfg':{'bd':4},
            'gridcfg':{'sticky':'nesw','row':-1, 'column':1},
            'command':self.Pause_cb})

        ifd.append({'name': 'closeB',
            'widgetType': Tkinter.Button,
            'text':'Close',
            'wcfg':{'bd':4},
            #'gridcfg':{'sticky':'nesw','row':-1, 'column':1},
            'gridcfg':{'sticky':'nesw', 'columnspan':2},
            'command':self.Close_cb})
        form = InputForm(self.master, self.root,
                         descr = ifd, modal = 0,
                         blocking = 0,closeWithWindow=1)
        form.ifd = ifd
        if self.hasCounter:
            ctr = ifd.entryByName['statesCounter']['widget']
            entF = ctr.component('entryfield')
            form.ent2 = entF._entryFieldEntry
            da = ctr.component('downarrow')
            ua = ctr.component('uparrow')
            for item in [da,ua]:
                item.bind('<ButtonPress-1>', self.SetState_cb, '+')
            form.ent2.bind('<Return>', self.SetState_cb, '+')
            form.counter = form.ifd.entryByName['statesCounter']['widget']
        form.stopB = form.ifd.entryByName['stopB']['widget']
        form.playB = form.ifd.entryByName['playB']['widget']
        form.playRevB = form.ifd.entryByName['playRevB']['widget']
        #print 'returning form1'
        self.form = form
        return form


    #methods for changing playMode
    def SetMode_cb(self, event=None):
        #print 'SetMode'
        #playMode options:
        #   0   play once and stop
        #   1   play continuously in 1 direction
        #   2   play once in 2 directions
        #   3   play continuously in 2 directions
        #play framerate is frame/per second
        if not hasattr(self, 'playModeForm'):
            self.playModeList=[ 'once and stop', 
                                'continuously in 1 direction',
                                'once in 2 directions', 
                                'continuously in 2 directions']
            ifd2 = InputFormDescr(title='Set Play Mode')    

            ifd2.append( {'name': 'playModeLabel',
                    'widgetType':Tkinter.Label,
                    'wcfg':{'text':'play mode options:', 
                        'font':(ensureFontCase('helvetica'),12,'bold')},
                    'gridcfg':{'sticky':'w'}})

            ifd2.append({'name':'playMode',
                        'widgetType': Tkinter.Radiobutton,
                        'defaultValue': self.playModeList[self.playMode],
                        'listtext':self.playModeList,
                        'gridcfg':{'sticky':'w'}})

            ifd2.append({'name': 'framerateTW',
                    'widgetType':ThumbWheel,
                    'tooltip':
                       """Framerate to enforce during playback""",
                    'gridcfg':{'sticky':'we'},
                    'wcfg':{
                             'value':self.framerate, 'oneTurn':100., 
                             'type':'float',
                             'continuous': True,
                             'wheelPad':2,'width':145,'height':18,
                             'labCfg':{'text':'framerate:  '},
                        },
                         })

            ifd2.append({'name': 'startFrameTW',
                    'widgetType':ThumbWheel,
                    'tooltip':
                       """First frame used in playback""",
                    'gridcfg':{'sticky':'we'},
                    'wcfg':{
                             'value':self.startFrame, 'oneTurn':100, 
                             'type':'int',
                             'continuous': True,
                             'wheelPad':2,'width':145,'height':18,
                             'labCfg':{'text':'start frame: '},
                        },
                         })

            ifd2.append({'name': 'endFrameTW',
                    'widgetType':ThumbWheel,
                    'tooltip':
                       """Last frame used in playback""",
                    'gridcfg':{'sticky':'we'},
                    'wcfg':{
                             'value':self.endFrame, 'oneTurn':100, 
                             'type':'int',
                             'continuous': True,
                             'wheelPad':2,'width':145,'height':18,
                             'labCfg':{'text':' end frame: '},
                        },
                         })

            ifd2.append({'name': 'stepSizeTW',
                    'widgetType':ThumbWheel,
                    'tooltip':
                       """???""",
                    'gridcfg':{'sticky':'we'},
                    'wcfg':{
                             'value':self.stepSize, 'oneTurn':100, 
                             'type':'int',
                             'continuous': True,
                             'wheelPad':2,'width':145,'height':18,
                             'labCfg':{'text':' step size:  '},
                        },
                         })

            ifd2.append({'name':'acceptB',
                        'widgetType': Tkinter.Button,
                        'wcfg':{
                            'text': 'ok',
                            'command': self.setPlayMode_cb,
                        },
                        'gridcfg':{'sticky':'nesw'}})

            ifd2.append({'name':'cancelB',
                        'widgetType': Tkinter.Button,
                        'wcfg':{
                            'text': 'cancel',
                            'command': self.cancelPlayMode_cb,
                        },
                        'gridcfg':{'sticky':'nesw', 'row':-1, 'column':1}})
            if self.master is None:
                master = self.root
            else:
                master = Tkinter.Toplevel()
            self.playModeForm = InputForm(master, None,
                        descr = ifd2,
                        modal = 0, blocking = 0)
            self.playModeVar = self.playModeForm.descr.entryByName['playMode']['variable']
            self.framerateWidget = self.playModeForm.descr.entryByName['framerateTW']['widget']
            self.startFrameWidget = self.playModeForm.descr.entryByName['startFrameTW']['widget']
            self.endFrameWidget = self.playModeForm.descr.entryByName['endFrameTW']['widget']
            self.stepSizeWidget = self.playModeForm.descr.entryByName['stepSizeTW']['widget']
        else:
            self.playModeForm.deiconify()
        

    def setPlayMode_cb(self, event=None):
        curVal = self.playModeVar.get()
        self.playMode = self.playModeList.index(curVal)
        #print 'setting playMode to ', curVal
        self.framerate = round(self.framerateWidget.get(),4)
        #print 'setting self.framerate ->', self.framerate
        self.timestamp= 1./self.framerate

        self.startFrame = self.startFrameWidget.get()
        self.endFrame = self.endFrameWidget.get()
        self.stepSize = self.stepSizeWidget.get()
        self.cancelPlayMode_cb()
        self.oneDirection = 0
        

    def cancelPlayMode_cb(self, event=None):
        self.playModeForm.withdraw()


    #methods for counter
    def custom_validate(self, text):
        #print 'in custom_validate, text=', text 
        if not len(text):
            return -1

        if text in ['start', 'end']:
            return 1

        tt = int(text)
        okList = range(self.startFrame, self.endFrame+1)
        if tt in okList:
            return 1
        else:
            return -1


    def custom_counter(self, text, factor, increment, **kw):
        # text is current content of entry
        # factor is 1 for increment and -1 for decrement
        # increment is value of increment megawidget option
        ###if not text in self.idList:
            ###raise ValueError, text + ' not in current idList'
        #check whether ind+factor is in range
        newval = self.currentFrameIndex + factor
        #print 'newval=', newval
        if newval<0 or newval>self.endFrame:
            #print 'custom_counter returning ', text
            return text
        else:
            #print 'custom_counter returning ', newval
            return newval


if __name__ == '__main__':
    def foo(val):
        print val
    root = Tkinter.Tk()
    root.withdraw()
    pl = Player( master=root, endFrame=10, maxFrame=10, form2=1)


