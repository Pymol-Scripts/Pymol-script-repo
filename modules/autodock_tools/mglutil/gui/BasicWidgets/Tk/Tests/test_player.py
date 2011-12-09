#########################################################################
#
# Date: Jan 2003 Authors: Ruth Huey,  ????
#
#    rhuey@scripps.edu
#    ?????@scripps.edu
#
# Copyright:  Ruth Huey, and TSRI
#
#########################################################################

#
#
# $Header: /opt/cvs/python/packages/share1.5/mglutil/gui/BasicWidgets/Tk/Tests/test_player.py,v 1.3 2005/06/21 16:16:03 sowjanya Exp $
#
# $Id: test_player.py,v 1.3 2005/06/21 16:16:03 sowjanya Exp $
#
#
#
#



import sys, Tkinter,unittest
from mglutil.regression import testplus
from mglutil.gui.BasicWidgets.Tk.player import Player
from time import sleep


widget = None
wasCalled = 0


def pause(sleepTime=1):
    widget.master.update()
    sleep(sleepTime)

class PlayerBaseTest(unittest.TestCase):
    def test_constructor(self):
        # test if we can display a very basic Player
        global widget
        root = Tkinter.Toplevel()
        widget = Player( master=root,gui=1)
        widget.master.update()
        pause()
        widget.master.master.destroy()


    def test_constructorOptionsForm1(self):
        # test all possible constructor options
        global widget
        root = Tkinter.Toplevel()
        widget = Player(    master=root,
                        width=100, height=26,
                        currentFrameIndex=1,
                        startFrame=1, 
                        endFrame=10,
                        maxFrame=30,
                        stepSize=3, 
                        playMode=3,
                        titleStr='Player Widget',
                        gotoStartfile = 'stop.gif',
                        gotoEndfile = 'stop.gif',
                        ff_revfile = 'stop.gif',
                        ff_fwdfile = 'stop.gif',
                        stopfile = 'stop.gif',
                        playfile = 'stop.gif',
                        playRevfile = 'stop.gif',
                        chmodfile = 'stop.gif',
                        closefile = 'stop.gif',
                        iconpath = None,
                        counter = 0,
                        form2=0,gui=1)
        widget.master.update()
        pause()
        widget.master.master.destroy()


    def test_constructorOptionsForm2(self):
        # test all possible constructor options
        global widget
        root = Tkinter.Toplevel()
        widget = Player(    master=root,
                        width=100, height=26,
                        currentFrameIndex=1,
                        startFrame=1, 
                        endFrame=10,
                        maxFrame=30,
                        stepSize=3, 
                        playMode=3,
                        titleStr='Player Widget',
                        gotoStartfile = 'stop.gif',
                        gotoEndfile = 'stop.gif',
                        ff_revfile = 'stop.gif',
                        ff_fwdfile = 'stop.gif',
                        stopfile = 'stop.gif',
                        playfile = 'stop.gif',
                        playRevfile = 'stop.gif',
                        chmodfile = 'stop.gif',
                        closefile = 'stop.gif',
                        iconpath = None,
                        counter = 0,
                        form2=1,gui=1)

        pause()
        widget.master.master.destroy()


    def test_nextFrame(self):
        #THIS SHOULD BE OVERWRITTEN for other players
        global widget
        root = Tkinter.Toplevel()
        widget = Player(    master=root,
                        width=100, height=26,
                        currentFrameIndex=1,
                        startFrame=1, 
                        endFrame=10,
                        form2=1,gui=1)
        widget.nextFrame(5)
        self.assertEqual(widget.currentFrameIndex, 5)
        pause()
        widget.master.master.destroy()


    def test_getnextFrame(self):
        global widget
        root = Tkinter.Toplevel()
        widget = Player(    master=root,
                        width=100, height=26,
                        currentFrameIndex=1,
                        startFrame=1, 
                        endFrame=10,
                        form2=1,gui=1)
        frameIndex0 = widget.currentFrameIndex
        frameIndex1 = widget.getNextFrameIndex(frameIndex0)
        self.assertEqual(frameIndex1 == frameIndex0 + widget.stepSize*widget.increment,True)
        widget.stepSize = 2
        frameIndex2 = widget.getNextFrameIndex(frameIndex1)
        self.assertEqual(frameIndex2 == frameIndex1 + 2*widget.increment,True)
        widget.increment = -1
        frameIndex3 = widget.getNextFrameIndex(frameIndex2)
        self.assertEqual(frameIndex3 == frameIndex2 - 2,True)
        pause()
        widget.master.master.destroy()
    #add code to change increment, step size
    #add code to test what happens at the end of range
    #add code to test the different play modes


    ## def xtest_waitTime():
    ##     global widget
    ##     from mglutil.gui.BasicWidgets.Tk.player import Player
    ##     import time
    ##     root = Tkinter.Toplevel()
    ##      
    ##     widget = Player(    master=root,
    ##                         width=100, height=26,
    ##                         currentFrameIndex=1,
    ##                         startFrame=1, 
    ##                         endFrame=10,
    ##                         form2=1,gui=1)
    ##     t0 = time.time()
    ##     print t0
    ##     widget.waitTime()
    ##     t1 = time.time()
    ##     print t1
    ##     assert t1!= t0
    ##     pause()
    ##     widget.master.master.destroy()


    def test_Play_cb(self):
        global widget
        root = Tkinter.Toplevel()
        widget = Player(    master=root,
                        width=100, height=26,
                        currentFrameIndex=1,
                        startFrame=1, 
                        endFrame=10,
                        form2=1,gui=1)
        widget.Play_cb()
        self.assertEqual(widget.currentFrameIndex == widget.endFrame,True)
        pause()
        widget.master.master.destroy()


    def test_PlayRev_cb(self):
        global widget
        root = Tkinter.Toplevel()
        widget = Player(    master=root,
                        width=100, height=26,
                        currentFrameIndex=1,
                        startFrame=1, 
                        endFrame=10,
                        form2=1,gui=1)
        widget.nextFrame(widget.endFrame)
        widget.PlayRev_cb()
        self.assertEqual(widget.startFrame==widget.currentFrameIndex,True)
        pause()
        widget.master.master.destroy()


    def test_FastForward_cb(self):
        global widget
        root = Tkinter.Toplevel()
        widget = Player(    master=root,
                        width=100, height=26,
                        currentFrameIndex=1,
                        startFrame=1, 
                        endFrame=10,
                        form2=1,gui=1)
        widget.FastForward_cb()
        self.assertEqual(widget.currentFrameIndex == widget.endFrame,True)
        pause()
        widget.master.master.destroy()


    def test_FastReverse_cb(self):
        global widget
        root = Tkinter.Toplevel()
        widget = Player(    master=root,
                        width=100, height=26,
                        currentFrameIndex=1,
                        startFrame=1, 
                        endFrame=10,
                        form2=1,gui=1)
        widget.nextFrame(widget.endFrame)
        widget.FastReverse_cb()
        self.assertEqual(widget.startFrame==widget.currentFrameIndex,True)
        pause()
        widget.master.master.destroy()


    def test_Stop_cb(self):
        global widget
        root = Tkinter.Toplevel()
        widget = Player(    master=root,
                        width=100, height=26,
                        currentFrameIndex=1,
                        startFrame=1, 
                        endFrame=10,
                        form2=1,gui=1)
        widget.Stop_cb()
        self.assertEqual(widget.stop==1,True)
        pause()
        widget.master.master.destroy()


    def test_GoToStart_cb(self):
        global widget
        root = Tkinter.Toplevel()
        widget = Player(    master=root,
                        width=100, height=26,
                        currentFrameIndex=1,
                        startFrame=1, 
                        endFrame=10,
                        form2=1)
        widget.GoToStart_cb()
        self.assertEqual(widget.currentFrameIndex==widget.startFrame,True)
        pause()
        widget.master.master.destroy()


    def test_GoToEnd_cb(self):
        global widget
        root = Tkinter.Toplevel()
        widget = Player(    master=root,
                        width=100, height=26,
                        currentFrameIndex=1,
                        startFrame=1, 
                        endFrame=10,
                        form2=1,gui=1)
        widget.GoToEnd_cb()
        self.assertEqual(widget.currentFrameIndex==widget.endFrame,True)
        pause()
        widget.master.master.destroy()


    def test_SetState_cb(self):
        global widget
        root = Tkinter.Toplevel()
        widget = Player(    master=root,
                        width=100, height=26,
                        currentFrameIndex=1,
                        startFrame=1, 
                        endFrame=10,
                        form2=1,gui=1)
        widget.form.ent2.delete(0, 'end')
        widget.form.ent2.insert(0, '5')
        widget.SetState_cb()
        self.assertEqual(widget.currentFrameIndex==5,True)
        pause()
        pause()
        widget.master.master.destroy()


    def test_SetMode_cb(self):
        global widget
        root = Tkinter.Toplevel()
        widget = Player(    master=root,
                        width=100, height=26,
                        currentFrameIndex=1,
                        startFrame=1, 
                        endFrame=10,
                        form2=1,gui=1)
        widget.SetMode_cb()
        self.assertEqual(widget.playModeForm.root.winfo_ismapped(),1)
        pause()
        widget.master.master.destroy()
    #add tests of changing mode here


    def test_SetPlayMode_cb(self):
        global widget
        root = Tkinter.Toplevel()
        widget = Player(    master=root,
                        width=100, height=26,
                        currentFrameIndex=1,
                        startFrame=1, 
                        endFrame=10,
                        form2=1,gui=1)
        widget.SetMode_cb()
        widget.playModeVar.set('once in 2 directions')
        widget.framerateWidget.set(10.)
        #widget.playDelayWidget.set(.5)
        #widget.playAfterDelayWidget.set(57)
        widget.startFrameWidget.set(2)
        widget.endFrameWidget.set(8)
        widget.setPlayMode_cb()
        self.assertEqual(widget.playMode==2,True)
        self.assertEqual(widget.framerate==10.,True)
        #assert widget.delay==.5
        #assert widget.afterDelay==57
        self.assertEqual(widget.startFrame==2,True)
        self.assertEqual(widget.endFrame==8,True)
        self.assertEqual(widget.oneDirection==0,True)
        pause()
        widget.master.master.destroy()
    #add tests of changing mode here


    def test_cancelPlayMode_cb(self):
        global widget
        root = Tkinter.Toplevel()
        widget = Player(    master=root,
                        width=100, height=26,
                        currentFrameIndex=1,
                        startFrame=1, 
                        endFrame=10,
                        form2=1,gui=1)
        widget.SetMode_cb()
        
        widget.master.master.update()
        pause()
        widget.cancelPlayMode_cb()
        self.assertEqual(widget.playModeForm.root.winfo_ismapped()==0,True)
        pause()
        widget.master.master.destroy()
        #add tests of changing mode here


    def test_Close_cb(self):
        global widget
        root = Tkinter.Toplevel()
        widget = Player(    master=root,
                        width=100, height=26,
                        currentFrameIndex=1,
                        startFrame=1, 
                        endFrame=10,
                        form2=1,gui=1)
        widget.Close_cb()
        self.assertEqual(widget.form.root.winfo_ismapped()==0,True)
        pause()
        widget.master.master.destroy()


if __name__ == '__main__':
    unittest.main()
