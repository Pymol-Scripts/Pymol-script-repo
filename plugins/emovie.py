'''
See more at: http://www.pymolwiki.org/index.php/emovie

#eMovie plugin tool for PyMOL
#28-Dec-2006
#created by Eran Hodis during summer 2006
#at the Israel Structural Proteomics Center
#Weizmann Institute of Science, Rehovot, Israel

#editing functions for storyboard added 28-Jul-07 - EH

CHANGELOG

2015-01-23 Thomas Holder
    * Fixed "Make Morph" feature for PyMOL 1.6+
2017-06-02 Thomas Holder
    * Python 3 compatibility

'''

from pymol import cmd
import sys

if sys.version_info[0] < 3:
    from Tkinter import *
    import tkSimpleDialog
    import tkMessageBox
    import tkFileDialog
else:
    from tkinter import *
    import tkinter as Tkinter
    from tkinter import simpledialog as tkSimpleDialog
    from tkinter import messagebox as tkMessageBox
    from tkinter import filedialog as tkFileDialog

import os
import pickle
import time

#
# KR: Check whether Morphing is available
#
if 'rigimol' in cmd.keyword or 'morph' in cmd.keyword:
    MORPHING_OPTIONS = 1
else:
    MORPHING_OPTIONS = 0


class eMovie:

    # create a storyBoard list to hold a storyboard of tuples (frames, action, info,script)
    storyBoard = []
    morphList = []
    sceneList = []

    # this is important for editing actions using the storyboard view
    actionToEdit = -1

    def __init__(self, master):

        frame = Frame(master)
        frame.pack()

        self.scenesButton = Button(frame, text="Scenes (views, appearances)", command=self.scenes, width=25, relief=RIDGE, fg="black")
        self.scenesButton.pack()

        self.rotationButton = Button(frame, text="     Rotation     ", command=self.rotation, width=25, relief=RIDGE, fg="black")
        self.rotationButton.pack()

        self.zoomButton = Button(frame, text="     Zoom     ", command=self.zoom, width=25, relief=RIDGE, fg="black")
        self.zoomButton.pack()

        self.commandButton = Button(frame, text="     Add Command     ", command=self.command, width=25, relief=RIDGE, fg="black")
        self.commandButton.pack()

        self.fadingButton = Button(frame, text="     Fading     ", command=self.fading, width=25, relief=RIDGE, fg="black")
        self.fadingButton.pack()

        self.wormButton = Button(frame, text="     Backbone Trace     ", command=self.worm, width=25, relief=RIDGE, fg="black")
        self.wormButton.pack()

        self.pauseButton = Button(frame, text="     Add Pause     ", command=self.pause, width=25, relief=RIDGE, fg="black")
        self.pauseButton.pack()

        self.stopButton = Button(frame, text="     Add Stop     ", command=self.stop, width=25, relief=RIDGE, fg="black")
        self.stopButton.pack()

        # KR--altered by EH to include only makeMorph button:
        # morphing buttons will disappear if rigimol is not installed
        if MORPHING_OPTIONS:
            self.makeMorphButton = Button(frame, text="     Make Morph     ", command=self.makeMorph, width=25, relief=RIDGE, fg="darkblue")
            self.makeMorphButton.pack()

        self.addMorphButton = Button(frame, text="Add Morph to eMovie", command=self.addMorph, width=25, relief=RIDGE, fg="darkblue")
        self.addMorphButton.pack()

        self.loadMorphButton = Button(frame, text="Load Previously Made Morph", command=self.loadMorph, width=25, relief=RIDGE, fg="darkblue")
        self.loadMorphButton.pack()

        self.storyButton = Button(frame, text="     View Storyboard     ", command=self.story, width=25, fg="blue", relief=RIDGE)
        self.storyButton.pack()

        self.saveButton = Button(frame, text="     Save eMovie     ", command=self.save, width=25, fg="purple", relief=RIDGE)
        self.saveButton.pack()

        self.loadButton = Button(frame, text="     Load eMovie     ", command=self.load, width=25, fg="magenta", relief=RIDGE)
        self.loadButton.pack()

        self.exportButton = Button(frame, text="     Export eMovie     ", command=self.export, width=25, fg="red", relief=RIDGE)
        self.exportButton.pack()

    def insertScene(self):
        insertSceneWindow = InsertScene(root, title="Insert Scene")

    def scenes(self):
        scenesWindow = Scenes(root, title="Scenes")

    def createScene(self):
        createSceneWindow = CreateScene(root, title="Create Scene")

    def zoom(self):
        zoomWindow = Zoom(root, title="Zoom")

    def command(self):
        commandWindow = Command(root, title="Add Command")

    def rotation(self):
        rotationWindow = Rotation(root, title="Rotation")

    def fading(self):
        fadingWindow = Fading(root, title="Fading")

    def worm(self):
        wormWindow = Worm(root, title="Backbone Trace")

    def pause(self):
        pauseWindow = Pause(root, title="Add Pause")

    def stop(self):
        stopWindow = Stop(root, title="Add Stop")

    def makeMorph(self):
        makeMorphWindow = MakeMorph(root, title="Make Morph")

    def addMorph(self):
        addMorphWindow = AddMorph(root, title="Add Morph")

    def loadMorph(self):
        loadMorphWindow = LoadMorph(root, title="Load Morph")

    def story(self):
        storyWindow = Story(root, title="Storyboard")

    def save(self):
        saveWindow = Save(root, title="Save")

    def load(self):
        loadWindow = Load(root, title="Load")

    def export(self):
        exportWindow = Export(root, title="Export")

    def moveActions(self):
        moveActionsWindow = MoveActions(root, title="Move Actions")

    def editAction(self):
        editActionWindow = EditAction(root, title="Edit Action")

# need to initialize plugin:


def __init__(self):

    self.menuBar.addcascademenu('Plugin', 'eMovie', 'eMovie plugin', label='eMovie')
    self.menuBar.addmenuitem('eMovie', 'command', 'Open eMovie', label='Open eMovie', command=lambda: open_eMovie(self.root))

    # also initialize Kristian Rother's movie.py commands
    #
    # all movie commands
    #
    # define a set of new PyMOL commands for creation of movies
    # this part will be executed first upon the
    # 'run movie.py' command from within PyMOL.
    cmd.extend('movie', movie)
    cmd.extend('mv_clear', mv_clear)  # clear cached instructions
    cmd.extend('mv_rot', mv_rotate)     # rotate over frame range
    cmd.extend('mv_move', mv_move)   # move over frame range
    cmd.extend('mv_trans', mv_trans)   # translate over frame range
    cmd.extend('mv_turn', mv_turn)   # turn over frame range
    cmd.extend('mv_cmd', mv_cmd)     # do anything during a frame range
    cmd.extend('mv_set', mv_set)     # adjust parameter over frame range

    # downward compatibility section
    cmd.extend('mvClear', mv_clear)  # clear cached instructions
    cmd.extend('mvRot', mv_rotate)     # rotate over frame range
    cmd.extend('mvMove', mv_move)   # move over frame range
    cmd.extend('mvTurn', mv_turn)   # turn over frame range
    cmd.extend('mvCmd', mv_cmd)     # do anything during a frame range
    cmd.extend('mvSet', mv_set)     # adjust parameter over frame range

    # KR: finally, call open_eMovie:
    # UN-COMMENTED TO REMOVE AUTO OPEN
    # open_eMovie()


def open_eMovie(parent=None):
    global root
    root = Toplevel(parent)
    root.title(' eMovie 1.04 ')
    global emovie
    emovie = eMovie(root)

#
# The following are classes for our button actions
# they use the Dialog class defined within tkSimpleDialog.py
#


class InsertScene(tkSimpleDialog.Dialog):

    def body(self, master):

        Label(master, text="Scene name to insert :").grid(row=0, column=0)
        Label(master, text="At frame :").grid(row=1, column=0)

        self.e1 = Entry(master)
        self.e2 = Entry(master)

        self.e1.grid(row=0, column=1)
        self.e2.grid(row=1, column=1)

        return self.e1

    def apply(self):

        sceneName = self.e1.get()
        startFrame = int(self.e2.get())

        (startFrameState, endFrameState, blankCmd) = getStartEndFrameStates(startFrame)

        # add command to movie
        mv_cmd("%i:%i" % (startFrame, startFrameState), 'cmd.scene("%s","recall")' % (sceneName))

        # append to storyboard
        emovie.storyBoard.append(("%i" % (startFrame), "Scene Set", "Scene name: %s" % (sceneName), ['eMovie.mv_cmd("%i:%i","cmd.scene(\'%s\',\'recall\')")' % (startFrame, startFrameState, sceneName), blankCmd]))

        movie()
        cmd.do("mstop")


class Scenes(tkSimpleDialog.Dialog):

    def body(self, master):

        Label(master, text="Scene List").grid(row=0, column=0)

        if len(emovie.sceneList) > 10:
            HeightLb = 10
            verticalScrollFlag = True
        else:
            HeightLb = len(emovie.sceneList)
            verticalScrollFlag = False

        scrollbar = Scrollbar(master, orient=VERTICAL)
        self.lb1 = Listbox(master, yscrollcommand=scrollbar.set, relief=RIDGE, height=HeightLb)
        scrollbar.config(command=self.yview)
        if verticalScrollFlag:
            scrollbar.grid(row=1, column=1, sticky=N + S)

        emovie.sceneList.sort()

        for a in emovie.sceneList:
            self.lb1.insert(END, a)

        self.lb1.grid(row=1, column=0)

        b1 = Button(master, text="Create new scene/Save scene...", command=self.createScene, width=27)
        b1.grid(row=2, column=0)
        b2 = Button(master, text="Recall selected scene", command=self.recallScene, width=27)
        b2.grid(row=3, column=0)
        b3 = Button(master, text="Delete selected scene", command=self.deleteScene, width=27)
        b3.grid(row=4, column=0)
        b4 = Button(master, text="Insert a scene into movie", command=self.addScene, width=27)
        b4.grid(row=5, column=0)

    def yview(self, *args):
        self.lb1.yview(*args)

    def addScene(self):
        emovie.insertScene()

    def recallScene(self):
        selection = self.lb1.curselection()
        try:
            selection = list(map(int, selection))
        except ValueError:
            pass
        sceneToRecall = emovie.sceneList[selection[0]]
        cmd.scene("%s" % (sceneToRecall), "recall")

    def deleteScene(self):
        selection = self.lb1.curselection()
        try:
            selection = list(map(int, selection))
        except ValueError:
            pass
        sceneToDelete = emovie.sceneList[selection[0]]
        cmd.scene("%s" % (sceneToDelete), "clear")
        emovie.sceneList.remove(sceneToDelete)
        self.destroy()
        emovie.scenes()

    def createScene(self):
        self.destroy()
        emovie.createScene()

    def buttonbox(self):  # overide tkSimpleDialog button box in order to insert a help button
        '''add standard button box.

        override if you do not want the standard buttons
        '''

        box = Frame(self)

        w = Button(box, text="OK", width=10, command=self.ok, default=ACTIVE)
        w.pack(side=LEFT, padx=5, pady=5)

        w = Button(box, text="Cancel", width=10, command=self.cancel)
        w.pack(side=LEFT, padx=5, pady=5)

        # help button
        helpButton = Button(box, command=self.help, bitmap="question", width=20)
        helpButton.pack(side=LEFT, padx=5, pady=5)

        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)

        box.pack()

    def help(self):

        tkMessageBox.showinfo("Scenes Help", "Create new Scene:  Here you create a new PyMOL scene.  The current PyMOL view, colors, and/or representations, depending on your choices, are/is saved as a scene with the name of your choice.  This new scene will appear on the scene list, and can then be inserted to the movie at a specific frame.\n\nRecall Selected Scene: The selected scene is immediately recalled.  You can then edit it and save your changes to the scene by using the 'Create new scene/Save scene...' button.\n\nDelete Scene:  The selected scene is cleared.  Any instance of the scene inserted into the movie is not removed.  To do so, you must remove instances of the scene from the storyboard.\n\nInsert Scene into Movie: The scene of the given name is recalled at the given frame.")


class CreateScene(tkSimpleDialog.Dialog):

    def body(self, master):

        Label(master, text="New scene name (no spaces) :").grid(row=0, column=0)
        Label(master, text="Store view? (y/n) :").grid(row=1, column=0)
        Label(master, text="Store colors? (y/n) :").grid(row=2, column=0)
        Label(master, text="Store representations? (y/n) :").grid(row=3, column=0)

        self.e1 = Entry(master)
        self.e2 = Entry(master)
        self.e3 = Entry(master)
        self.e4 = Entry(master)

        self.e1.grid(row=0, column=1)
        self.e2.grid(row=1, column=1)
        self.e3.grid(row=2, column=1)
        self.e4.grid(row=3, column=1)

        self.e2.insert(END, "y")
        self.e3.insert(END, "y")
        self.e4.insert(END, "y")

        return self.e1

    def apply(self):

        sceneName = self.e1.get()
        viewYN = self.e2.get()
        colorsYN = self.e3.get()
        repYN = self.e4.get()

        if viewYN.lower() == "n" or viewYN.lower() == "no":
            viewFlag = 0
        else:
            viewFlag = 1
        if colorsYN.lower() == "n" or colorsYN.lower() == "no":
            colorsFlag = 0
        else:
            colorsFlag = 1
        if repYN.lower() == "n" or repYN.lower() == "no":
            repFlag = 0
        else:
            repFlag = 1

        cmd.scene(sceneName, "store", "", viewFlag, colorsFlag, 1, repFlag, 0)

        while emovie.sceneList.count(sceneName) > 0:  # remove all previous instances of sceneName
            emovie.sceneList.remove(sceneName)

        emovie.sceneList.append(sceneName)

        emovie.scenes()


class Zoom(tkSimpleDialog.Dialog):

    def body(self, master):

        Label(master, text="Angstroms to zoom (can be negative) :").grid(row=0)
        Label(master, text="Start frame :").grid(row=1)
        Label(master, text="Action length (in frames) :").grid(row=2)

        self.e1 = Entry(master)
        self.e2 = Entry(master)
        self.e3 = Entry(master)

        self.e1.grid(row=0, column=1)
        self.e2.grid(row=1, column=1)
        self.e3.grid(row=2, column=1)

        return self.e1  # initial focus

    def apply(self):

        # take the 3 parameters from the dialog box and execute the movie command for zoom

        amtZoom = int(self.e1.get())
        startFrame = int(self.e2.get())
        actionLength = int(self.e3.get())

        endFrame = startFrame + actionLength - 1

        (startFrameState, endFrameState, blankCmd) = getStartEndFrameStates(startFrame, endFrame)

        # add command to movie
        mv_move("%i-%i:%i-%i" % (startFrame, endFrame, startFrameState, endFrameState), "z", "%i" % (amtZoom), "linear")

        # append to storyBoard
        emovie.storyBoard.append(("%i-%i" % (startFrame, endFrame), "Zoom", "Angstroms: %i" % (amtZoom), ['mv_move %i-%i:%i-%i,z,%i,linear' % (startFrame, endFrame, startFrameState, endFrameState, amtZoom), blankCmd]))

        movie()
        cmd.do("mstop")

    def buttonbox(self):  # overide tkSimpleDialog button box in order to insert a help button
        '''add standard button box.

        override if you do not want the standard buttons
        '''

        box = Frame(self)

        w = Button(box, text="OK", width=10, command=self.ok, default=ACTIVE)
        w.pack(side=LEFT, padx=5, pady=5)

        w = Button(box, text="Cancel", width=10, command=self.cancel)
        w.pack(side=LEFT, padx=5, pady=5)

        # help button
        helpButton = Button(box, command=self.help, bitmap="question", width=20)
        helpButton.pack(side=LEFT, padx=5, pady=5)

        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)

        box.pack()

    def help(self):

        tkMessageBox.showinfo("Zoom Help", "Moves the camera a given number of angstroms along the z-axis.  The action length is the amount of frames it takes for the completion of the zoom.  Positive angstroms will zoom-in, while negative angstroms constitute zooming-out.")
        self.lift()


class Command(tkSimpleDialog.Dialog):

    def body(self, master):

        Label(master, text="PyMOL command :").grid(row=0)
        Label(master, text="Execute at frame :").grid(row=1)

        self.e1 = Entry(master)
        self.e2 = Entry(master)

        self.e1.grid(row=0, column=1)
        self.e2.grid(row=1, column=1)

        return self.e1  # initial focus

    def apply(self):

        # take the 2 parameters from the dialog box and execute the movie command for adding a command at a given frame

        inputCmd = self.e1.get()
        startFrame = int(self.e2.get())

        (startFrameState, endFrameState, blankCmd) = getStartEndFrameStates(startFrame)

        # add command to movie
        mv_cmd("%i:%i" % (startFrame, startFrameState), inputCmd)

        # append to storyBoard
        emovie.storyBoard.append(("%i" % (startFrame), "Command", inputCmd, ['eMovie.mv_cmd("%i:%i","%s")' % (startFrame, startFrameState, inputCmd), blankCmd]))

        movie()
        cmd.do("mstop")

    def buttonbox(self):  # overide tkSimpleDialog button box in order to insert a help button
        '''add standard button box.

        override if you do not want the standard buttons
        '''

        box = Frame(self)

        w = Button(box, text="OK", width=10, command=self.ok, default=ACTIVE)
        w.pack(side=LEFT, padx=5, pady=5)

        w = Button(box, text="Cancel", width=10, command=self.cancel)
        w.pack(side=LEFT, padx=5, pady=5)

        # help button
        helpButton = Button(box, command=self.help, bitmap="question", width=20)
        helpButton.pack(side=LEFT, padx=5, pady=5)

        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)

        box.pack()

    def help(self):

        tkMessageBox.showinfo("Command Help", "Enter any PyMOL command to be executed at a particular frame.  The PyMOL command can include commas.  For the execution of many commands at once, try inserting a PyMOL command to execute an external script ('@script.py' or 'run script.py').")
        self.lift()


class Rotation(tkSimpleDialog.Dialog):

    def body(self, master):

        Label(master, text="Axis :").grid(row=0, column=0)
        Label(master, text="Degrees :").grid(row=1, column=0)
        Label(master, text="Start frame :").grid(row=2, column=0)
        Label(master, text="Action length :").grid(row=3, column=0)

        self.e1 = Entry(master)
        self.e2 = Entry(master)
        self.e3 = Entry(master)
        self.e4 = Entry(master)

        self.e1.grid(row=0, column=1)
        self.e2.grid(row=1, column=1)
        self.e3.grid(row=2, column=1)
        self.e4.grid(row=3, column=1)

        return self.e1  # initial focus

    def apply(self):

        # take the 4 parameters from the dialog box and execute the movie command for adding a rotation at a given frame

        axis = self.e1.get()
        degrees = self.e2.get()
        startFrame = int(self.e3.get())
        actionLength = int(self.e4.get())

        endFrame = startFrame + actionLength - 1

        (startFrameState, endFrameState, blankCmd) = getStartEndFrameStates(startFrame, endFrame)

        # add command to movie
        mv_turn("%i-%i:%i-%i" % (startFrame, endFrame, startFrameState, endFrameState), axis, degrees, "linear")

        # append to storyBoard
        emovie.storyBoard.append(("%i-%i" % (startFrame, endFrame), "Rotation", "Axis: %s; Degrees: %s" % (axis, degrees), ['mv_turn %i-%i:%i-%i,%s,%s,linear' % (startFrame, endFrame, startFrameState, endFrameState, axis, degrees), blankCmd]))

        movie()
        cmd.do("mstop")

    def buttonbox(self):  # overide tkSimpleDialog button box in order to insert a help button
        '''add standard button box.

        override if you do not want the standard buttons
        '''

        box = Frame(self)

        w = Button(box, text="OK", width=10, command=self.ok, default=ACTIVE)
        w.pack(side=LEFT, padx=5, pady=5)

        w = Button(box, text="Cancel", width=10, command=self.cancel)
        w.pack(side=LEFT, padx=5, pady=5)

        # help button
        helpButton = Button(box, command=self.help, bitmap="question", width=20)
        helpButton.pack(side=LEFT, padx=5, pady=5)

        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)

        box.pack()

    def help(self):

        tkMessageBox.showinfo("Rotation Help", "Inserts a rotation of a given degree about a given axis at a given frame.  The action length is the amount of frames it takes for the completion of the rotation.  Rotations along an axis other than x, y, or z can be created by stacking combinations of rotations.  To set the origin of rotation, use the PyMOL command 'origin selection', where selection is the object or selection at which you wish to place the origin of rotation.")
        self.lift()


class Fading(tkSimpleDialog.Dialog):

    def body(self, master):

        Label(master, text="Molecule :").grid(row=0, column=0)
        Label(master, text="Shown as :").grid(row=1, column=0)
        Label(master, text="Fade from :").grid(row=2, column=0)
        Label(master, text="Fade to   :").grid(row=3, column=0)
        Label(master, text="Start frame :").grid(row=4, column=0)
        Label(master, text="Action length :").grid(row=5, column=0)

        Label(master, text="(stick, surface, cartoon, sphere)").grid(row=1, column=2, sticky=W)
        Label(master, text="percent visible").grid(row=2, column=2, sticky=W)
        Label(master, text="percent visible").grid(row=3, column=2, sticky=W)

        self.e1 = Entry(master)
        self.e2 = Entry(master)
        self.e3 = Entry(master)
        self.e4 = Entry(master)
        self.e5 = Entry(master)
        self.e6 = Entry(master)

        self.e1.grid(row=0, column=1)
        self.e2.grid(row=1, column=1)
        self.e3.grid(row=2, column=1)
        self.e4.grid(row=3, column=1)
        self.e5.grid(row=4, column=1)
        self.e6.grid(row=5, column=1)

        return self.e1  # initial focus

    def apply(self):

        # take the 4 parameters from the dialog box and execute the movie command for adding a rotation at a given frame

        molecule = self.e1.get()
        representation = self.e2.get()
        fadeFrom = float(self.e3.get())
        fadeTo = float(self.e4.get())
        startFrame = int(self.e5.get())
        actionLength = int(self.e6.get())

        endFrame = startFrame + actionLength - 1

        (startFrameState, endFrameState, blankCmd) = getStartEndFrameStates(startFrame, endFrame)

        # convert fadeFrom and fadeTo from percentage to fractions of a whole
        fadeFromMod = fadeFrom / 100
        fadeToMod = fadeTo / 100
        # convert to amount of transparency rather than amt visibility
        fadeFromMod = 1 - fadeFromMod
        fadeToMod = 1 - fadeToMod

        if startFrameState < endFrameState:
            anteEndFrameState = endFrameState - 1
        elif startFrameState > endFrameState:
            anteEndFrameState = endFrameState + 1
        else:
            anteEndFrameState = endFrameState

        # since surface transparency setting is transparency and all others are of the form TYPE_transparency, we do:
        if representation == "surface":
            appearance = ""
        else:
            appearance = "%s_" % (representation)

        # add command to movie
        mv_set("%i-%i:%i-%i" % (startFrame, endFrame - 1, startFrameState, anteEndFrameState), "%stransparency" % (appearance), "%f" % (fadeFromMod), "%f" % (fadeToMod), molecule, 'linear')
        # because there is a issue with mv_set's linear setting of values, we force the last value on the last frame for now
        mv_cmd("%i:%i" % (endFrame, endFrameState), "set %stransparency,%f,%s" % (appearance, fadeToMod, molecule))

        # append to storyBoard
        emovie.storyBoard.append(("%i-%i" % (startFrame, endFrame), "Fading", "Molecule: %s; Shown as: %s; Fade from: %.2f%%; Fade to: %.2f%%" % (molecule, representation, fadeFrom, fadeTo), ['mv_set %i-%i:%i-%i,%stransparency,%f,%f,%s,linear' % (startFrame, endFrame - 1, startFrameState, anteEndFrameState, appearance, fadeFromMod, fadeToMod, molecule), 'mv_cmd %i:%i,cmd.set("%stransparency","%f","%s")' % (endFrame, endFrameState, appearance, fadeToMod, molecule), blankCmd]))

        movie()
        cmd.do("mstop")

    def buttonbox(self):  # overide tkSimpleDialog button box in order to insert a help button
        '''add standard button box.

        override if you do not want the standard buttons
        '''

        box = Frame(self)

        w = Button(box, text="OK", width=10, command=self.ok, default=ACTIVE)
        w.pack(side=LEFT, padx=5, pady=5)

        w = Button(box, text="Cancel", width=10, command=self.cancel)
        w.pack(side=LEFT, padx=5, pady=5)

        # help button
        helpButton = Button(box, command=self.help, bitmap="question", width=20)
        helpButton.pack(side=LEFT, padx=5, pady=5)

        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)

        box.pack()

    def help(self):

        tkMessageBox.showinfo("Fading Help", "The representation of the molecule/selection will be faded from one visibility level to the other visibility level over the action length.  Make sure this representation is being shown; otherwise the fading will not be visible.  You may fade one or more representations of a molecule/selection while simultaneously showing other representations.")


class Worm(tkSimpleDialog.Dialog):

    def body(self, master):

        Label(master, text="Molecule :").grid(row=0, column=0)
        Label(master, text="Starting Residue :").grid(row=1, column=0)
        Label(master, text="Ending Residue :").grid(row=2, column=0)
        Label(master, text="Start Frame :").grid(row=3, column=0)

        self.e1 = Entry(master)
        self.e2 = Entry(master)
        self.e3 = Entry(master)
        self.e4 = Entry(master)

        self.e1.grid(row=0, column=1)
        self.e2.grid(row=1, column=1)
        self.e3.grid(row=2, column=1)
        self.e4.grid(row=3, column=1)

        return self.e1  # initial focus

    def apply(self):

        # take the parameters from the dialog box and execute the movie command for adding a worm feature at a given frame

        molecule = self.e1.get()
        startAA = int(self.e2.get())
        endAA = int(self.e3.get())
        startFrame = int(self.e4.get())

        (startFrameState, endFrameState, blankCmd) = getStartEndFrameStates(startFrame)

        lastFrame = wormFunction(molecule, startAA, endAA, startFrame, startFrameState)

        movie()
        cmd.do("mstop")

        emovie.storyBoard.append(("%i-%i" % (startFrame, lastFrame), "Worm", "Molecule: %s; Residues: %i-%i" % (molecule, startAA, endAA), ['eMovie.wormFunction("%s",%i,%i,%i,%i)' % (molecule, startAA, endAA, startFrame, startFrameState), blankCmd]))

        tkMessageBox.showinfo("Backbone Trace Feature", "Backbone trace inserted at frames %i to %i.\n" % (startFrame, lastFrame))

    def buttonbox(self):  # overide tkSimpleDialog button box in order to insert a help button
        '''add standard button box.

        override if you do not want the standard buttons
        '''

        box = Frame(self)

        w = Button(box, text="OK", width=10, command=self.ok, default=ACTIVE)
        w.pack(side=LEFT, padx=5, pady=5)

        w = Button(box, text="Cancel", width=10, command=self.cancel)
        w.pack(side=LEFT, padx=5, pady=5)

        # help button
        helpButton = Button(box, command=self.help, bitmap="question", width=20)
        helpButton.pack(side=LEFT, padx=5, pady=5)

        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)

        box.pack()

    def help(self):

        tkMessageBox.showinfo("Backbone Trace Help", "This function first draws the first 30 residues between the given starting amino acid and the given ending amino acid in stick representation.  Then a cartoon of the carbon backbone (default setting = cartoon tube) begins to appear starting at the given start amino acid and ending at the given end amino acid.  The cartoon adds one amino acid to its structure with each passing frame.  Meanwhile, the 30 residues shown as sticks fade out entirely.  The cartoon is colored as a rainbow spectrum.  While this trace occurs, the molecule/selection rotates about the y axis.  The user cannot specify an action length with this feature because the action length depends on the amount of residues in the backbone trace.")
        self.lift()


class Pause(tkSimpleDialog.Dialog):

    def body(self, master):

        Label(master, text="How many frames should the pause last? :").grid(row=0, column=0)
        Label(master, text="Start pause at frame :").grid(row=1, column=0)

        self.e1 = Entry(master)
        self.e2 = Entry(master)

        self.e1.grid(row=0, column=1)
        self.e2.grid(row=1, column=1)

        return self.e1  # initial focus

    def apply(self):

        actionLength = int(self.e1.get())
        startFrame = int(self.e2.get())

        endFrame = startFrame + actionLength - 1

        (startFrameState, endFrameState, blankCmd) = getStartEndFrameStates(startFrame, endFrame)

        # add pause to movie
        mv_cmd("%i-%i:%i-%i" % (startFrame, endFrame, startFrameState, endFrameState), "")

        # append to storyBoard
        emovie.storyBoard.append(("%i-%i" % (startFrame, endFrame), "Pause", "-", ['mv_cmd %i-%i:%i-%i' % (startFrame, endFrame, startFrameState, endFrameState), blankCmd]))

        movie()
        cmd.do("mstop")

    def buttonbox(self):  # overide tkSimpleDialog button box in order to insert a help button
        '''add standard button box.
        override if you do not want the standard buttons
        '''

        box = Frame(self)

        w = Button(box, text="OK", width=10, command=self.ok, default=ACTIVE)
        w.pack(side=LEFT, padx=5, pady=5)

        w = Button(box, text="Cancel", width=10, command=self.cancel)
        w.pack(side=LEFT, padx=5, pady=5)

        # help button
        helpButton = Button(box, command=self.help, bitmap="question", width=20)
        helpButton.pack(side=LEFT, padx=5, pady=5)

        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)

        box.pack()

    def help(self):

        tkMessageBox.showinfo("Pause Help", "A pause is added to the movie for the specified action length, at the specified frame.  This is not a hard-fast pause, as actions can be inserted into the pause.  Also, it is pointless to insert a pause in the middle of a movie, because simply by lack of inserted actions over a given frame range, a pause will occur.  The purpose of this button is to insert pauses to the end of movies.")
        self.lift()


class Stop(tkSimpleDialog.Dialog):

    def body(self, master):

        Label(master, text="Insert stop at frame :").grid(row=0, column=0)

        self.e1 = Entry(master)

        self.e1.grid(row=0, column=1)

        return self.e1

    def apply(self):

        startFrame = int(self.e1.get())

        (startFrameState, endFrameState, blankCmd) = getStartEndFrameStates(startFrame)

        # add pause to movie
        mv_cmd("%i:%i" % (startFrame, startFrameState), "mstop")

        # append to storyBoard
        emovie.storyBoard.append(("%i" % (startFrame), "Stop", "-", ['mv_cmd %i:%i,mstop' % (startFrame, startFrameState), blankCmd]))

        movie()
        cmd.do("mstop")

    def buttonbox(self):  # overide tkSimpleDialog button box in order to insert a help button
        '''add standard button box.
        override if you do not want the standard buttons
        '''

        box = Frame(self)

        w = Button(box, text="OK", width=10, command=self.ok, default=ACTIVE)
        w.pack(side=LEFT, padx=5, pady=5)

        w = Button(box, text="Cancel", width=10, command=self.cancel)
        w.pack(side=LEFT, padx=5, pady=5)

        # help button
        helpButton = Button(box, command=self.help, bitmap="question", width=20)
        helpButton.pack(side=LEFT, padx=5, pady=5)

        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)

        box.pack()

    def help(self):

        tkMessageBox.showinfo("Stop Help", "A stop is inserted into the movie at the specified frame.  When the movie plays, and reaches the frame with the stop, the movie will stop until the user presses the play button again.  Stops are useful for automatically stopping the movie at specific points during a presentation so that the presenter can explain, and then continue the movie at his or her leisure.  Warning: Remove all stop actions from the storyboard before exporting a movie.")
        self.lift()


class Story(tkSimpleDialog.Dialog):

    # override simpleDialog initialization to disable self.grab_set()

    def __init__(self, parent, title=None):
        '''Initialize a dialog.

        Arguments:

        parent -- a parent window (the application window)

        title -- the dialog title
        '''
        Toplevel.__init__(self, parent)
        self.transient(parent)

        if title:
            self.title(title)

        self.parent = parent

        self.result = None

        body = Frame(self)
        self.initial_focus = self.body(body)
        body.pack(padx=5, pady=5)

        self.buttonbox()

        # self.grab_set()  #comment this out so that user can keep storyboard open while using other eMovie functions

        if not self.initial_focus:
            self.initial_focus = self

        self.protocol("WM_DELETE_WINDOW", self.cancel)

        if self.parent is not None:
            self.geometry("+%d+%d" % (parent.winfo_rootx() + 50, parent.winfo_rooty() + 50))
            self.initial_focus.focus_set()

        self.wait_window(self)

    def body(self, master):

        # set the widths for the listboxes based on the entries in emovie.storyBoard
        # find the max length string in each entry in emovie.storyboard
        WidthLb1 = 1
        WidthLb2 = 1
        WidthLb3 = 1
        WidthLb4 = 1
        verticalScrollFlag = False
        horizontalScrollFlag = False

        i = 0
        for a in emovie.storyBoard:
            i = i + 1
            WidthLb1 = int(i / 10) + 1  # set this width to the number of digits in the number of items in the storyboard
            if len(a[0]) > WidthLb2:
                WidthLb2 = len(a[0])
            if len(a[1]) > WidthLb3:
                WidthLb3 = len(a[1]) + 1
            if len(a[2]) > WidthLb4:
                WidthLb4 = len(a[2])

        # set the heights of the listboxes so that they are 10 max, but adaptive from 0 to 10
        if i <= 10:
            HeightLb = i
        else:
            HeightLb = 10
            verticalScrollFlag = True

        # now we make sure the WidthLb3 doesn't exceed our max desired width
        if WidthLb4 > 45:
            WidthLb4 = 45
            horizontalScrollFlag = True

        scrollbar = Scrollbar(master, orient=VERTICAL)
        scrollbarInfo = Scrollbar(master, orient=HORIZONTAL)
        self.lb1 = Listbox(master, yscrollcommand=scrollbar.set, relief=FLAT, width=WidthLb1, height=HeightLb)
        self.lb2 = Listbox(master, yscrollcommand=scrollbar.set, relief=FLAT, width=WidthLb2, height=HeightLb)
        self.lb3 = Listbox(master, yscrollcommand=scrollbar.set, relief=FLAT, width=WidthLb3, height=HeightLb)
        self.lb4 = Listbox(master, yscrollcommand=scrollbar.set, xscrollcommand=scrollbarInfo.set, relief=FLAT, width=WidthLb4, height=HeightLb)
        scrollbar.config(command=self.yview)
        if verticalScrollFlag:
            scrollbar.grid(row=1, column=4, sticky=N + S)
        scrollbarInfo.config(command=self.xview)
        if horizontalScrollFlag:
            scrollbarInfo.grid(row=2, column=0, columnspan=5, sticky=E + W)
        self.lb1.grid(row=1, column=0)
        self.lb2.grid(row=1, column=1)
        self.lb3.grid(row=1, column=2)
        self.lb4.grid(row=1, column=3)

        Label(master, text="#", fg="red").grid(row=0, column=0)
        Label(master, text="FRAME(S)", fg="red").grid(row=0, column=1)
        Label(master, text="ACTION", fg="red").grid(row=0, column=2)
        Label(master, text="INFORMATION", fg="red").grid(row=0, column=3)

        # sort Storyboard so it is displayed in order of startFrame, using the sorting function cmpStoryBoard

        emovie.storyBoard.sort(cmpStoryBoard)

        i = 1
        for a in emovie.storyBoard:
            self.lb1.insert(END, i)
            self.lb2.insert(END, a[0])
            self.lb3.insert(END, a[1])
            self.lb4.insert(END, a[2])
            i = i + 1

        # deleting actions button:
        deleteActionFlag = False
        deleteButton = Button(master, text="Delete selected action", command=self.deleteAction).grid(row=3, column=2, columnspan=3, sticky=E + W)

        # refresh storyboard button
        refreshButton = Button(master, text="REFRESH", command=self.refresh).grid(row=3, column=0, columnspan=2, sticky=E + W)

        # edit an action
        editActionButton = Button(master, text="Edit action", command=self.editActionButton).grid(row=4, column=0, columnspan=2, sticky=E + W)

        # move a group of actions button
        moveActionsButton = Button(master, text="Move a group of actions", command=self.moveActionsButton).grid(row=4, column=2, columnspan=3, sticky=E + W)

    def editActionButton(self):
        # need to get the index number of the selection, the selection could be in any of the 4 listboxes
        actionNumber = -1
        if len(self.lb1.curselection()) == 1:
            # get the selection but in some versions of Tkinter, the selection is given as ints, and some as strings, so we check for both
            selection = self.lb1.curselection()
            try:
                selection = list(map(int, selection))
            except ValueError:
                pass
            actionNumber = selection[0]
        elif len(self.lb2.curselection()) == 1:
            selection = self.lb2.curselection()
            try:
                selection = list(map(int, selection))
            except ValueError:
                pass
            actionNumber = selection[0]
        elif len(self.lb3.curselection()) == 1:
            selection = self.lb3.curselection()
            try:
                selection = list(map(int, selection))
            except ValueError:
                pass
            actionNumber = selection[0]
        elif len(self.lb4.curselection()) == 1:
            selection = self.lb4.curselection()
            try:
                selection = list(map(int, selection))
            except ValueError:
                pass
            actionNumber = selection[0]
        emovie.actionToEdit = actionNumber
        emovie.editAction()
        self.destroy()
        emovie.story()

    def moveActionsButton(self):
        emovie.moveActions()
        self.destroy()
        emovie.story()

    def refresh(self):
        self.destroy()
        emovie.story()

    def deleteAction(self):
        deleteActionFlag = True
        self.apply(deleteActionFlag)
        self.destroy()
        emovie.story()

    def yview(self, *args):
        self.lb1.yview(*args)
        self.lb2.yview(*args)
        self.lb3.yview(*args)
        self.lb4.yview(*args)

    def xview(self, *args):
        self.lb4.xview(*args)

    def apply(self, deleteActionFlag=False):

        # need to get the index number of the selection, the selection could be in any of the 4 listboxes
        actionNumber = -1
        if len(self.lb1.curselection()) == 1:
            # get the selection but in some versions of Tkinter, the selection is given as ints, and some as strings, so we check for both
            selection = self.lb1.curselection()
            try:
                selection = list(map(int, selection))
            except ValueError:
                pass
            actionNumber = selection[0]
        elif len(self.lb2.curselection()) == 1:
            selection = self.lb2.curselection()
            try:
                selection = list(map(int, selection))
            except ValueError:
                pass
            actionNumber = selection[0]
        elif len(self.lb3.curselection()) == 1:
            selection = self.lb3.curselection()
            try:
                selection = list(map(int, selection))
            except ValueError:
                pass
            actionNumber = selection[0]
        elif len(self.lb4.curselection()) == 1:
            selection = self.lb4.curselection()
            try:
                selection = list(map(int, selection))
            except ValueError:
                pass
            actionNumber = selection[0]

        if (deleteActionFlag == True) and (actionNumber != -1):
            # delete the actionNumber-th tuple in emovie.storyBoard
            deletedAction = emovie.storyBoard.pop(actionNumber)
            tkMessageBox.showinfo("Deleted Action", "Action number %s deleted.  \nAction frames: %s . \nAction: %s. \nInformation: %s. \n" % (actionNumber + 1, deletedAction[0], deletedAction[1], deletedAction[2]))

            # now need to reload movie:
            # clear movie and add storyBoard commands to movie:
            mv_clear()

            cmd.do("from pmg_tk.startup import eMovie")
            for a in emovie.storyBoard:
                for b in a[3]:  # action commands are stored in 4th position of tuples in storyBoard list, as a list of cmds
                    cmd.do(b)

            time.sleep(2)

            movie()
            cmd.do("mstop")

    def buttonbox(self):  # overide tkSimpleDialog button box in order to insert a help button
        '''add standard button box.

        override if you do not want the standard buttons
        '''

        box = Frame(self)

        w = Button(box, text="OK", width=10, command=self.ok, default=ACTIVE)
        w.pack(side=LEFT, padx=5, pady=5)

        w = Button(box, text="Cancel", width=10, command=self.cancel)
        w.pack(side=LEFT, padx=5, pady=5)

        # help button
        helpButton = Button(box, command=self.help, bitmap="question", width=20)
        helpButton.pack(side=LEFT, padx=5, pady=5)

        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)

        box.pack()

    def help(self):

        tkMessageBox.showinfo("Storyboard Help", "A list of the actions that comprise the movie, organized according to the frame numbers in which they take place.  By selecting an action and clicking 'Delete selected action', the selected action is removed from the movie, and the movie is reloaded and the storyboard refreshed. 'Edit action' allows changing of a selected action's parameters. 'Move a group of actions' allows movement of a group of consecutive actions up or down by a specified number of frames.")
        self.lift()


class MoveActions(tkSimpleDialog.Dialog):
        # take every action between and including firstAction and lastAction and remove them from our emovie.storyboard
        # store each of the actions we removed in an array
        # clear and recompile our movie without the removed actions
        # one by one, get the values for each action in our new array,
        # then change the start and end Frame by the moveAmt
        # finally reinsert each moved action into the movie with insertByActionValues
        # recompiling the movie after insertion of each action

    def body(self, master):

        Label(master, text="Group starts at action number :").grid(row=0, column=0)
        Label(master, text="Group ends with action number :").grid(row=1, column=0)
        Label(master, text="Move how many frames \n (negative number for backward) :").grid(row=2, column=0)

        self.e1 = Entry(master)
        self.e2 = Entry(master)
        self.e3 = Entry(master)

        self.e1.grid(row=0, column=1)
        self.e2.grid(row=1, column=1)
        self.e3.grid(row=2, column=1)

        return self.e1

    def apply(self):

        firstAction = int(self.e1.get())
        lastAction = int(self.e2.get())
        moveAmt = int(self.e3.get())

        queue = []

        for i in range(firstAction - 1, lastAction):
            # remove each action from emovie.storyBoard and put into new array
            deletedAction = emovie.storyBoard.pop(firstAction - 1)
            queue.append(deletedAction)

        # now need to reload movie:
        # clear movie and add storyBoard commands to movie:
        mv_clear()

        cmd.do("from pmg_tk.startup import eMovie")
        for a in emovie.storyBoard:
            for b in a[3]:  # action commands are stored in 4th position of tuples in storyBoard list, as a list of cmds
                cmd.do(b)

        time.sleep(2)

        movie()
        cmd.do("mstop")

        for i in range(firstAction - 1, lastAction):
            actionToAdd = queue.pop(0)
            # get the values of the action
            currentValues = getActionValues(actionToAdd)
            # change the start and end frames of the action to add
            newStartFrame = int(currentValues[0]) + moveAmt
            newEndFrame = int(currentValues[1]) + moveAmt
            newFrameValues = newStartFrame, newEndFrame
            newValues = newFrameValues + currentValues[2:]
            # now reinsert the action
            insertActionByValues(newValues)


class EditAction(tkSimpleDialog.Dialog):
        # the user selected some action in the storyboard
        # from that selected action, we get the values using getActionValues function
        # depending on the type of action, we display the relevant values in input boxes for editing
        # the user clicks OK and we get the values from the input boxes
        # we then delete the action from the movie, recompile the movie, and then add the changed action

    def body(self, master):

        if emovie.actionToEdit == -1:
            Label(master, text="No action selected to edit.").grid(row=0, column=0)
        else:
            actionValues = getActionValues(emovie.storyBoard[emovie.actionToEdit])

            actionType = actionValues[2]

            # start with the values common to all actions: Action Type (uneditable) and Start frame
            Label(master, text="Action type :").grid(row=0, column=0)
            Label(master, text="Start frame :").grid(row=1, column=0)
            Label(master, text=actionType).grid(row=0, column=1, sticky=W)

            self.e1 = Entry(master)
            self.e1.grid(row=1, column=1)
            self.e1.insert(END, actionValues[0])

            # now show editable values based on the type of action chosen for editing
            if actionType == "Pause":

                Label(master, text="End frame :").grid(row=2, column=0)

                self.e2 = Entry(master)
                self.e2.grid(row=2, column=1)
                self.e2.insert(END, actionValues[1])

            elif actionType == "Scene Set":

                Label(master, text="Scene name :").grid(row=2, column=0)

                self.e2 = Entry(master)
                self.e2.grid(row=2, column=1)
                self.e2.insert(END, actionValues[3])

            elif actionType == "Zoom":

                Label(master, text="End frame :").grid(row=2, column=0)
                Label(master, text="Angstroms to zoom : \n (negative to zoom out)").grid(row=3, column=0)

                self.e2 = Entry(master)
                self.e3 = Entry(master)

                self.e2.grid(row=2, column=1)
                self.e3.grid(row=3, column=1)

                self.e2.insert(END, actionValues[1])
                self.e3.insert(END, actionValues[3])

            elif actionType == "Command":

                Label(master, text="Command to execute :").grid(row=2, column=0)

                self.e2 = Entry(master)
                self.e2.grid(row=2, column=1)
                self.e2.insert(END, actionValues[3])

            elif actionType == "Rotation":

                Label(master, text="End frame :").grid(row=2, column=0)
                Label(master, text="Axis :").grid(row=3, column=0)
                Label(master, text="Degrees :").grid(row=4, column=0)

                self.e2 = Entry(master)
                self.e3 = Entry(master)
                self.e4 = Entry(master)

                self.e2.grid(row=2, column=1)
                self.e3.grid(row=3, column=1)
                self.e4.grid(row=4, column=1)

                self.e2.insert(END, actionValues[1])
                self.e3.insert(END, actionValues[3])
                self.e4.insert(END, actionValues[4])

            elif actionType == "Fading":

                Label(master, text="End frame :").grid(row=2, column=0)
                Label(master, text="Molecule :").grid(row=3, column=0)
                Label(master, text="Representation :\n (stick, surface, cartoon, sphere)").grid(row=4, column=0)
                Label(master, text="Fade from (% vis) :").grid(row=5, column=0)
                Label(master, text="Fade to (% vis) :").grid(row=6, column=0)

                self.e2 = Entry(master)
                self.e3 = Entry(master)
                self.e4 = Entry(master)
                self.e5 = Entry(master)
                self.e6 = Entry(master)

                self.e2.grid(row=2, column=1)
                self.e3.grid(row=3, column=1)
                self.e4.grid(row=4, column=1)
                self.e5.grid(row=5, column=1)
                self.e6.grid(row=6, column=1)

                self.e2.insert(END, actionValues[1])
                self.e3.insert(END, actionValues[3])
                self.e4.insert(END, actionValues[4])
                self.e5.insert(END, actionValues[5])
                self.e6.insert(END, actionValues[6])

            elif actionType == "Worm":

                Label(master, text="Molecule :").grid(row=2, column=0)
                Label(master, text="Start amino acid :").grid(row=3, column=0)
                Label(master, text="End amino acid :").grid(row=4, column=0)

                self.e2 = Entry(master)
                self.e3 = Entry(master)
                self.e4 = Entry(master)

                self.e2.grid(row=2, column=1)
                self.e3.grid(row=3, column=1)
                self.e4.grid(row=4, column=1)

                self.e2.insert(END, actionValues[3])
                self.e3.insert(END, actionValues[4])
                self.e4.insert(END, actionValues[5])
            elif actionType == "Stop":
                pass
            else:
                print("Error: Action type of selected action not found.")

            return self.e1

    def apply(self):

        if emovie.actionToEdit == -1:  # do nothing
            self.destroy()
        else:
            # first remove the old action from the movie
            deletedAction = emovie.storyBoard.pop(emovie.actionToEdit)

            # now need to reload movie:
            # clear movie and add storyBoard commands to movie:
            mv_clear()

            cmd.do("from pmg_tk.startup import eMovie")
            for a in emovie.storyBoard:
                for b in a[3]:  # action commands are stored in 4th position of tuples in storyBoard list, as a list of cmds
                    cmd.do(b)

            time.sleep(2)

            movie()
            cmd.do("mstop")

            # e1 always holds the start frame
            startFrame = int(self.e1.get())

            actionValues = getActionValues(deletedAction)

            actionType = actionValues[2]

            # now append other values depending on the type of action
            if actionType == "Stop":

                endFrame = 0

                # create new editedValues tuple
                editedValues = startFrame, endFrame, actionType

            elif actionType == "Pause":

                endFrame = int(self.e2.get())

                editedValues = startFrame, endFrame, actionType

            elif actionType == "Scene Set":

                endFrame = 0
                sceneName = self.e2.get()

                editedValues = startFrame, endFrame, actionType, sceneName

            elif actionType == "Zoom":

                endFrame = int(self.e2.get())
                zoomAmt = self.e3.get()

                editedValues = startFrame, endFrame, actionType, zoomAmt

            elif actionType == "Command":

                endFrame = 0
                inputCmd = self.e2.get()

                editedValues = startFrame, endFrame, actionType, inputCmd

            elif actionType == "Rotation":

                endFrame = int(self.e2.get())
                axis = self.e3.get()
                degrees = self.e4.get()

                editedValues = startFrame, endFrame, actionType, axis, degrees

            elif actionType == "Fading":

                endFrame = int(self.e2.get())
                molecule = self.e3.get()
                representation = self.e4.get()
                fadeFrom = float(self.e5.get())
                fadeTo = float(self.e6.get())

                editedValues = startFrame, endFrame, actionType, molecule, representation, fadeFrom, fadeTo

            elif actionType == "Worm":

                endFrame = 0
                molecule = self.e2.get()
                startAA = self.e3.get()
                endAA = self.e4.get()

                editedValues = startFrame, endFrame, actionType, molecule, startAA, endAA

            else:
                print("Error: Action type of selected action not found.")

            # now reinsert the action
            insertActionByValues(editedValues)

            # and reset emovie.actionToEdit to -1
            emovie.actionToEdit = -1


# class Save(tkSimpleDialog.Dialog):
# KR: same as with Load
class Save:

    def __init__(self, root, title="Save"):
        #
        # KR: added this method to produce a dialog right away
        #
        filename = tkFileDialog.asksaveasfile(title=title, defaultextension=".emov", filetypes=[("eMovie files", ".emov"), ("All files", "*")])
        # added file type menu.
        if not filename:
            return
        fileName = filename.name

        # want to make sure the extension is .emov
        # split the string using "." as the delimiter string and check if last entry is "pse" or "emov" if so, remove them
        # if not, add emov to string

        splitFileName = fileName.split(".")
        if splitFileName[len(splitFileName) - 1] == "pse":
            fileName = fileName[:-4]
        elif splitFileName[len(splitFileName) - 1] == "emov":
            fileName = fileName[:-5]

        # save session

        mv_clear()  # delete the movie so that the sesson can be saved

        cmd.do("save %s.pse" % (fileName))

        f = open('%s.emov' % (fileName), 'w')

        # make x = storyBoard with morphList added as last item and sceneList as 2nd to last item (should do this just as a tuple(A,B,C) in the future)
        x = emovie.storyBoard
        x.append('0')
        x[len(x) - 1] = emovie.morphList
        x.append('0')
        x[len(x) - 1] = emovie.sceneList

        # dump x to file f
        pickle.dump(x, f)

        f.close()

        # must now load the movie again
        time.sleep(2)

        emovie.sceneList = x.pop()
        emovie.morphList = x.pop()
        emovie.storyBoard = x

        cmd.do("from pmg_tk.startup import eMovie")
        for a in emovie.storyBoard:
            for b in a[3]:  # action commands are stored in 4th position of tuples in storyBoard list, as a list of cmds
                cmd.do(b)

        time.sleep(2)

        movie()
        cmd.do("mstop")

    #
    # KR: again, some now dead code starts
    #
    def body(self, master):

        Label(master, text="File name (.emov) :").grid(row=0)

        self.e1 = Entry(master)

        self.e1.grid(row=0, column=1)

        browseButton = Button(master, text="Browse...", command=self.openFileDialog, relief=RIDGE).grid(row=1, column=1)

        return self.e1  # initial focus

    def openFileDialog(self):

        filename = tkFileDialog.asksaveasfile(title="Browse...", defaultextension=".emov")

        if filename:
            self.e1.delete(0, END)  # clear the entrybox
            self.e1.insert(END, filename.name)
            self.e1.xview(END)  # scroll so you can see the end of the filename

            filename.close()  # close the file that the tkFileDialog opened for writing
            os.remove(filename.name)  # delete the file that tkFileDialog created, we will make it ourselves later
            # now automatically press the OK button
            self.ok()
        else:
            self.lift()

    def apply(self):

        fileName = self.e1.get()

        # want to make sure the extension is .emov
        # split the string using "." as the delimiter string and check if last entry is "pse" or "emov" if so, remove them
        # if not, add emov to string

        splitFileName = fileName.split(".")
        if splitFileName[len(splitFileName) - 1] == "pse":
            fileName = fileName[:-4]
        elif splitFileName[len(splitFileName) - 1] == "emov":
            fileName = fileName[:-5]

        # save session

        mv_clear()  # delete the movie so that the sesson can be saved

        cmd.do("save %s.pse" % (fileName))

        f = open('%s.emov' % (fileName), 'w')

        # make x = storyBoard with morphList added as last item and sceneList as 2nd to last item (should do this just as a tuple(A,B,C) in the future)
        x = emovie.storyBoard
        x.append('0')
        x[len(x) - 1] = emovie.morphList
        x.append('0')
        x[len(x) - 1] = emovie.sceneList

        # dump x to file f
        pickle.dump(x, f)

        f.close()

        # must now load the movie again
        time.sleep(2)

        emovie.sceneList = x.pop()
        emovie.morphList = x.pop()
        emovie.storyBoard = x

        cmd.do("from pmg_tk.startup import eMovie")
        for a in emovie.storyBoard:
            for b in a[3]:  # action commands are stored in 4th position of tuples in storyBoard list, as a list of cmds
                cmd.do(b)

        time.sleep(2)

        movie()
        cmd.do("mstop")

    def buttonbox(self):  # overide tkSimpleDialog button box in order to insert a help button
        '''add standard button box.

        override if you do not want the standard buttons
        '''

        box = Frame(self)

        w = Button(box, text="OK", width=10, command=self.ok, default=ACTIVE)
        w.pack(side=LEFT, padx=5, pady=5)

        w = Button(box, text="Cancel", width=10, command=self.cancel)
        w.pack(side=LEFT, padx=5, pady=5)

        # help button
        helpButton = Button(box, command=self.help, bitmap="question", width=20)
        helpButton.pack(side=LEFT, padx=5, pady=5)

        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)

        box.pack()
    #
    # KR: end dead code for Save
    #

    def help(self):

        tkMessageBox.showinfo("Save eMovie Help", "Saves the eMovie.  Two files are created: filename.emov and filename.pse.  Both files are crucial for later loading of the eMovie.  In giving a file name for saving, one can leave out the extension, or alternatively use either .emov or .pse.  First the session file is saved, and then the eMovie actions are saved.")
        self.lift()

# class Load(tkSimpleDialog.Dialog):
# KR: made Load dialog a simple class so no extra window will pop up


class Load:

    def __init__(self, root, title="Load"):
        #
        # KR: added this method to produce a dialog right away
        #
        filename = tkFileDialog.askopenfile(title=title, defaultextension=".emov", filetypes=[("eMovie files", ".emov"), ("All files", "*")])
        # added file type menu.
        if not filename:
            return
        fileName = filename.name

        # here we want to look at the fileName
        # if it does not have an extension of either .pse or .emov, we will use it as is
        # if it has an extenstion of either .pse or .emov, we will strip away that extension

        splitFileName = fileName.split(".")

        if splitFileName[len(splitFileName) - 1] == "emov":
            fileName = fileName[:-5]
        elif splitFileName[len(splitFileName) - 1] == "pse":
            fileName = fileName[:-4]

        cmd.delete("all")

        cmd.load("%s.pse" % (fileName))

        time.sleep(2)

        f = open('%s.emov' % (fileName), 'r')

        x = pickle.load(f)

        if isinstance(x[-1], bytes):
            temp = x.pop()  # to handle any old versions of emovies that saved initialview, just throws away the initial view

        if isinstance(x[-1], list) and isinstance(x[-2], list):  # check if sceneList and morphList are present because earlier versions didnt save sceneList
            emovie.sceneList = x.pop()
            emovie.morphList = x.pop()
            emovie.storyBoard = x
        else:
            emovie.sceneList = []
            emovie.morphList = x.pop()
            emovie.storyBoard = x

        try:
            emovie.morphList.remove(("MORPHS (30 frames each)", 0, 0))  # to handle any versions of emovies that saved this within morphList
        except:
            pass

        # for compatibility with old versions of eMovie:
        for a in emovie.storyBoard:
            for b in a[3]:
                u = b.split(".", 1)
                if (u[0] == "eMovie_box_beta") or (u[0] == "movie"):
                    newAction = "eMovie." + u[1]
                    x = a[3]
                    index = x.index(b)
                    x[index] = newAction

        # clear movie and add scriptList commands to movie:
        mv_clear()

        cmd.do("from pmg_tk.startup import eMovie")

        for a in emovie.storyBoard:
            for b in a[3]:  # action commands are stored in 4th position of tuples in storyBoard list, as a list of cmds
                cmd.do(b)

        f.close()

        time.sleep(2)

        movie()
        cmd.do("mstop")

    #
    # KR: the following stuff is IMHO unnecessary now
    #

    def body(self, master):

        Label(master, text="File name (.emov) :").grid(row=0)

        self.e1 = Entry(master)
        self.e1.grid(row=0, column=1)

        browseButton = Button(master, text="Browse...", command=self.openFileDialog, relief=RIDGE).grid(row=1, column=1)

        return self.e1  # initial focus

    def openFileDialog(self):
        filename = tkFileDialog.askopenfile(title="Browse...", defaultextension=".emov")

        if filename:
            self.e1.delete(0, END)  # clear the entrybox
            self.e1.insert(END, filename.name)
            self.e1.xview(END)  # scroll so you can see the end of the filename

            filename.close()  # close the file that the tkFileDialog opened for reading
            self.ok()
        else:
            self.lift()

    def apply(self):

        fileName = self.e1.get()

        # here we want to look at the fileName
        # if it does not have an extension of either .pse or .emov, we will use it as is
        # if it has an extenstion of either .pse or .emov, we will strip away that extension

        splitFileName = fileName.split(".")

        if splitFileName[len(splitFileName) - 1] == "emov":
            fileName = fileName[:-5]
        elif splitFileName[len(splitFileName) - 1] == "pse":
            fileName = fileName[:-4]

        cmd.delete("all")

        cmd.load("%s.pse" % (fileName))

        time.sleep(2)

        f = open('%s.emov' % (fileName), 'r')

        x = pickle.load(f)

        if isinstance(x[-1], bytes):
            temp = x.pop()  # to handle any old versions of emovies that saved initialview, just throws away the initial view

        if isinstance(x[-1], list) and isinstance(x[-2], list):  # check if sceneList and morphList are present because earlier versions didnt save sceneList
            emovie.sceneList = x.pop()
            emovie.morphList = x.pop()
            emovie.storyBoard = x
        else:
            emovie.sceneList = []
            emovie.morphList = x.pop()
            emovie.storyBoard = x

        try:
            emovie.morphList.remove(("MORPHS (30 frames each)", 0, 0))  # to handle any versions of emovies that saved this within morphList
        except:
            pass

        # for compatibility with old versions of eMovie:
        for a in emovie.storyBoard:
            for b in a[3]:
                u = b.split(".", 1)
                if (u[0] == "eMovie_box_beta") or (u[0] == "movie"):
                    newAction = "eMovie." + u[1]
                    x = a[3]
                    index = x.index(b)
                    x[index] = newAction

        # clear movie and add scriptList commands to movie:
        mv_clear()

        cmd.do("from pmg_tk.startup import eMovie")

        for a in emovie.storyBoard:
            for b in a[3]:  # action commands are stored in 4th position of tuples in storyBoard list, as a list of cmds
                cmd.do(b)

        f.close()

        time.sleep(2)

        movie()
        cmd.do("mstop")

    def buttonbox(self):  # overide tkSimpleDialog button box in order to insert a help button
        '''add standard button box.

        override if you do not want the standard buttons
        '''

        box = Frame(self)

        w = Button(box, text="OK", width=10, command=self.ok, default=ACTIVE)
        w.pack(side=LEFT, padx=5, pady=5)

        w = Button(box, text="Cancel", width=10, command=self.cancel)
        w.pack(side=LEFT, padx=5, pady=5)

        # help button
        helpButton = Button(box, command=self.help, bitmap="question", width=20)
        helpButton.pack(side=LEFT, padx=5, pady=5)

        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)

        box.pack()

    def help(self):

        tkMessageBox.showinfo("Load eMovie Help", "Loads an eMovie.  Two files must be present to load from: filename.emov and filename.pse.  You can specify either file, and the loading will occur successfully.  Alternatively, you can specify the filename without the extension.  First the session file is loaded, and then the eMovie actions are loaded.")
        self.lift()

    #
    # KR: end of dead code in Load class
    #


class MakeMorph(tkSimpleDialog.Dialog):

    def body(self, master):

        Label(master, text="Morph name (no extension, no spaces) :").grid(row=0)
        Label(master, text="Morph from :").grid(row=1)
        Label(master, text="Morph to :").grid(row=2)
        Label(master, text="Refinement amount :\n (1=little;20=medium;100=a lot)").grid(row=3)
        Label(master, text="Place resulting morph on top of which molecule :").grid(row=4)

        self.e1 = Entry(master)
        self.e2 = Entry(master)
        self.e3 = Entry(master)
        self.e4 = Entry(master)
        self.e5 = Entry(master)

        self.e1.grid(row=0, column=1)
        self.e2.grid(row=1, column=1)
        self.e3.grid(row=2, column=1)
        self.e4.grid(row=3, column=1)
        self.e5.grid(row=4, column=1)

        return self.e1

    def apply(self):

        morphName = self.e1.get()
        morphFrom = self.e2.get()
        morphTo = self.e3.get()
        refinementAmt = int(self.e4.get())
        placeOn = self.e5.get()

        # perform alignment
        if placeOn == morphFrom:
            cmd.align(morphTo, morphFrom)
        elif placeOn == morphTo:
            cmd.align(morphFrom, morphTo)

        numberOfSteps = 30  # this cannot be changed except through altering the rigimol input file (.inp)

        try:
            cmd.morph
        except AttributeError:
            # PyMOL version < 1.6
            cmd.do("from ipymol import rigimol")

            # save morphFrom and morphTo as morphFrom.pdb and morphTo.pdb
            cmd.save("morphFrom.pdb", morphFrom, state=1)
            cmd.save("morphTo.pdb", morphTo, state=1)

            time.sleep(2)

            # need to insert a pause here to allow for saving to finish

            cmd.do("delete all")
            # call rigimol.inp because you can't run it except through an external script
            cmd.do("rigimol eMovie_rigimol.inp")

            # time.sleep(120)  #change to pymol API command sync?
            tkMessageBox.showinfo("Making Morph", "Morph is in the process of being made.\n Wait for the PyMOL Tcl/Tk GUI window to read: \n 'RigiMOL: normal program termination.' \n Before pressing 'OK'")

            # delete morphFrom.pdb and morphTo.pdb
            os.remove("morphFrom.pdb")
            os.remove("morphTo.pdb")

            # run Delano's refinement commands on eMovie_rigimol_morph.pdb
            # refine.py:
            # Copyright (C) 2004 DeLano Scientific LLC.  All Rights Reserved.

            # to run: "pymol -qc refine.py" (command mode)
            #         "pymol -xFil refine.py" (to watch it run)

            # This script improves the atomic geometry of RigiMOL interpolations
            # without losing the desirable smoothness of these interpolations.

            # NOTE: this script can take a while to run (up to dozens of minutes)

            # define input and output files

            input_file = "eMovie_rigimol_morph.pdb"
            output_file = "%s_eMorph.pdb" % (morphName)

            # how much refinement should be done?

            # 1   = almost nothing
            # 20  = a reasonable amount
            # 100 = a lot

            how_much_refinement = refinementAmt

            # ======================================
            # no changes usually required below here

            object_name = "refining"

            cmd.do("load %s,%s" % (input_file, object_name))
            # cmd.show("ribbon")
            #cmd.spectrum(selection="elem c")

            cmd.do('rigimol.refine("%i","%s")' % (how_much_refinement, object_name))

            # cmd.alter('all','type="HETATM"')

            cmd.save(output_file, object_name, state=-1)  # multimodel PDB file

            cmd.delete(object_name)
            os.remove("eMovie_rigimol_morph.pdb")
        else:
            # PyMOL version >= 1.6
            cmd.morph(morphName, morphFrom, morphTo, refinement=refinementAmt, steps=numberOfSteps)

            # save to disk and delete everything to stay compatible with eMovie protocol
            cmd.save("%s_eMorph.pdb" % (morphName), morphName, 0)
            cmd.delete("all")

        # append morphName.pdb into morphList

        if len(emovie.morphList) == 0:
            lastState = 0
        else:
            for a in emovie.morphList:
                lastState = a[2]

        emovie.morphList.append((morphName, lastState + 1, lastState + numberOfSteps, "%s_eMorph.pdb" % (morphName)))

        # load all morphs including new one
        for a in emovie.morphList:
            cmd.do('load %s,morph' % (a[3]))

        # append morphName.pdb into scriptList so it's saved with the script and loads automatically
        #emovie.scriptList.append('load %s,morph'%(output_file))

    def buttonbox(self):  # overide tkSimpleDialog button box in order to insert a help button
        '''add standard button box.

        override if you do not want the standard buttons
        '''

        box = Frame(self)

        w = Button(box, text="OK", width=10, command=self.ok, default=ACTIVE)
        w.pack(side=LEFT, padx=5, pady=5)

        w = Button(box, text="Cancel", width=10, command=self.cancel)
        w.pack(side=LEFT, padx=5, pady=5)

        # help button
        helpButton = Button(box, command=self.help, bitmap="question", width=20)
        helpButton.pack(side=LEFT, padx=5, pady=5)

        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)

        box.pack()

    def help(self):

        tkMessageBox.showinfo("Make Morph Help", "This feature only works with iPyMOL (incentive PyMOL), because iPyMOL is the only version to include the morphing tool RigiMOL.  A morph of 30 steps is created, the first and last steps being the structures the user specifies to morph from and to, and the 28 intermediate steps being RigiMOL's interpolations between the 2 specified structures.  The file 'eMovie_rigimol.inp' must be located in the current working directory.  The morph is generated and added to the eMovie session, but is also saved under the name 'morphname_eMorph.pdb' in the current working directory.  Refinement is important in creating real-looking motions/morphs, however it takes a lot of time.  If you specify a refinement value of 20, RigiMOL will cycle through your prepared morph 20 times, each time implementing its refinement algorithm.  A refinement value of 1 is pretty much no refinement, 20 is a medium amount, and 100 is a lot.  The entry entitled 'Place resulting morph on top of:' functions to align the molecules invovled in the morph; in this entry, specify the molecule (of the 2 being morphed) whose space you wish the resulting morph to occupy. If you have already performed alignment of the start and end states of the morph, leave this entry blank.  If you performing serial morphing of more than 2 morphs, it is suggested you do your own alignment using PyMOL and leave the 'Place resulting morph on top of:' entry blank.")
        self.lift()


class AddMorph(tkSimpleDialog.Dialog):

    # override simpleDialog initialization to disable self.grab_set()

    def __init__(self, parent, title=None):
        '''Initialize a dialog.

        Arguments:

        parent -- a parent window (the application window)

        title -- the dialog title
        '''
        Toplevel.__init__(self, parent)
        self.transient(parent)

        if title:
            self.title(title)

        self.parent = parent

        self.result = None

        body = Frame(self)
        self.initial_focus = self.body(body)
        body.pack(padx=5, pady=5)

        self.buttonbox()

        # self.grab_set()  #comment this out so that user can keep storyboard open while using other eMovie functions

        if not self.initial_focus:
            self.initial_focus = self

        self.protocol("WM_DELETE_WINDOW", self.cancel)

        if self.parent is not None:
            self.geometry("+%d+%d" % (parent.winfo_rootx() + 50, parent.winfo_rooty() + 50))
        self.initial_focus.focus_set()

        self.wait_window(self)

    def body(self, master):

        Label(master, text="MORPHS (30 frames each)", fg="purple").grid(row=0, column=0)

        if len(emovie.morphList) > 10:
            HeightLb = 10
            verticalScrollFlag = True
        else:
            HeightLb = len(emovie.morphList)
            verticalScrollFlag = False

        scrollbar = Scrollbar(master, orient=VERTICAL)
        self.lb1 = Listbox(master, yscrollcommand=scrollbar.set, relief=RIDGE, height=HeightLb)
        scrollbar.config(command=self.yview)
        if verticalScrollFlag:
            scrollbar.grid(row=1, column=1, sticky=N + S)

        for a in emovie.morphList:
            self.lb1.insert(END, a[0])

        self.lb1.grid(row=1, column=0)

        b = Button(master, text="Click to add a morph to movie", command=self.openAddMorphParam)
        b.grid(row=2, column=0)

    def yview(self, *args):
        self.lb1.yview(*args)

    def openAddMorphParam(self):

        addMorphParamWindow = AddMorphParam(root)

    def buttonbox(self):  # overide tkSimpleDialog button box in order to insert a help button
        '''add standard button box.

        override if you do not want the standard buttons
        '''

        box = Frame(self)

        w = Button(box, text="OK", width=10, command=self.ok, default=ACTIVE)
        w.pack(side=LEFT, padx=5, pady=5)

        w = Button(box, text="Cancel", width=10, command=self.cancel)
        w.pack(side=LEFT, padx=5, pady=5)

        # help button
        helpButton = Button(box, command=self.help, bitmap="question", width=20)
        helpButton.pack(side=LEFT, padx=5, pady=5)

        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)

        box.pack()

    def help(self):

        tkMessageBox.showinfo("Add Morph to eMovie Help", "Click on 'Click to add morph to movie' and then a new window pops up.  Enter the name of a morph from the list and play it forwards, backwards, or in a loop (with an optional pause in the middle) at a specified frame.")
        self.lift()


class AddMorphParam(tkSimpleDialog.Dialog):

    def body(self, master):

        Label(master, text="Morph name (from list) :").grid(row=0, column=0)
        Label(master, text="Start frame :").grid(row=1, column=0)
        Label(master, text="Direction :\n(forward,backward,loop)").grid(row=2, column=0)
        Label(master, text="Pause at morphed stage :\n(only for loop display)").grid(row=3, column=0)
        Label(master, text="(forward and backward morph takes 30 frames \n loop takes 60 frames plus your pause)").grid(row=4, column=0, columnspan=2)

        self.e1 = Entry(master)
        self.e2 = Entry(master)
        self.e3 = Entry(master)
        self.e4 = Entry(master)

        self.e1.grid(row=0, column=1)
        self.e2.grid(row=1, column=1)
        self.e3.grid(row=2, column=1)
        self.e4.grid(row=3, column=1)

        self.e4.insert(END, "0")

        return self.e1

    def apply(self):

        morphName = self.e1.get()
        startFrame = int(self.e2.get())
        direction = self.e3.get()
        pause = int(self.e4.get())

        blankCmd = ""

        # need to get the start and end states of the morph
        for a in emovie.morphList:
            if a[0] == morphName:
                startState = a[1]
                endState = a[2]

        # need to adjust for states present in movie
        if startFrame > cmd.count_frames():
            # get state of last Frame, use blank cmd to set state of blank frames until start Frame
            lastState = getFrameState(cmd.count_frames())
            mv_cmd("%i-%i:%i" % (cmd.count_frames(), startFrame - 1, lastState), "")
            blankCmd = 'eMovie.mv_cmd("%i-%i:%i","")' % (cmd.count_frames(), startFrame - 1, lastState)

        if direction == "loop":

            # add command to movie
            mv_cmd("%i-%i:%i-%i" % (startFrame, startFrame + 29, startState, endState), "")
            mv_cmd("%i-%i:%i" % (startFrame + 30, startFrame + 30 + pause - 1, endState), "")
            mv_cmd("%i-%i:%i-%i" % (startFrame + 30 + pause, startFrame + 30 + pause + 29, endState, startState), "")

            # append to storyBoard
            emovie.storyBoard.append(("%i-%i" % (startFrame, startFrame + 29), "Morph", "Name: %s; Direction: forward" % (morphName), ['mv_cmd %i-%i:%i-%i' % (startFrame, startFrame + 29, startState, endState), blankCmd]))
            emovie.storyBoard.append(("%i-%i" % (startFrame + 30, startFrame + 30 + pause - 1), "Morph", "Name: %s; Action: pause" % (morphName), ['mv_cmd %i-%i:%i' % (startFrame + 30, startFrame + 30 + pause - 1, endState)]))
            emovie.storyBoard.append(("%i-%i" % (startFrame + 30 + pause, startFrame + 59 + pause), "Morph", "Name: %s; Direction: backwards" % (morphName), ['mv_cmd %i-%i:%i-%i' % (startFrame + 30 + pause, startFrame + 30 + pause + 29, endState, startState)]))

        elif direction == "backward":

            # add command to movie
            mv_cmd("%i-%i:%i-%i" % (startFrame, startFrame + 29, endState, startState), "")

            # append to storyBoard
            emovie.storyBoard.append(("%i-%i" % (startFrame, startFrame + 29), "Morph", "Name: %s; Direction: %s" % (morphName, direction), ['mv_cmd %i-%i:%i-%i' % (startFrame, startFrame + 29, endState, startState), blankCmd]))

        else:
            # add command to movie
            mv_cmd("%i-%i:%i-%i" % (startFrame, startFrame + 29, startState, endState), "")

            # append to storyBoard
            emovie.storyBoard.append(("%i-%i" % (startFrame, startFrame + 29), "Morph", "Name: %s; Direction: %s" % (morphName, direction), ['mv_cmd %i-%i:%i-%i' % (startFrame, startFrame + 29, startState, endState), blankCmd]))

        movie()
        cmd.do("mstop")


class LoadMorph(tkSimpleDialog.Dialog):

    def body(self, master):

        Label(master, text="Enter filename of morph (.pdb) :").grid(row=0, column=0)
        Label(master, text="Enter name to give morph :").grid(row=1, column=0)
        Label(master, text="(by using an existing morph name, existing morph will be overwritten)").grid(row=2, column=0, columnspan=3)

        self.e1 = Entry(master)
        self.e2 = Entry(master)

        self.e1.grid(row=0, column=1)
        self.e2.grid(row=1, column=1)

        browseButton = Button(master, text="Browse...", command=self.openFileDialog, relief=RIDGE).grid(row=0, column=2)

        return self.e1

    def openFileDialog(self):
        filename = tkFileDialog.askopenfile(title="Browse...", defaultextension=".pdb")

        if filename:
            self.e1.delete(0, END)  # clear the entrybox
            self.e1.insert(END, filename.name)
            self.e1.xview(END)  # scroll so you can see the end of the filename

            filename.close()  # close the file that the tkFileDialog opened for reading

        self.lift()

    def apply(self):

        morphFile = self.e1.get()
        morphName = self.e2.get()

        numberOfSteps = 30  # this pertains to the number of states in the loaded morph, should find a way to autodetect this but as of right now all morphs are made with 30 steps

        # check if morphName matches an existing morph name in morphList, if so, and numberOfSteps of each matches, replace it, if not, print error
        morphNameExistsFlag = False
        for a in emovie.morphList:
            if a[0] == morphName:
                existingMorph = a
                morphNameExistsFlag = True
        if morphNameExistsFlag == True:
            if (int(existingMorph[2]) - int(existingMorph[1]) + 1) == numberOfSteps:
                b = (morphName, existingMorph[1], existingMorph[2], morphFile)
                index = emovie.morphList.index(existingMorph)
                emovie.morphList.insert(index, b)
                temp = emovie.morphList.pop(index + 1)
                cmd.do('load %s,morph,%i' % (morphFile, existingMorph[1]))
        else:
            if len(emovie.morphList) == 0:
                lastState = 0
            else:
                for a in emovie.morphList:
                    lastState = a[2]

            emovie.morphList.append((morphName, lastState + 1, lastState + numberOfSteps, morphFile))

            # load new morph
            cmd.do('load %s,morph' % (morphFile))

    def buttonbox(self):  # overide tkSimpleDialog button box in order to insert a help button
        '''add standard button box.

        override if you do not want the standard buttons
        '''

        box = Frame(self)

        w = Button(box, text="OK", width=10, command=self.ok, default=ACTIVE)
        w.pack(side=LEFT, padx=5, pady=5)

        w = Button(box, text="Cancel", width=10, command=self.cancel)
        w.pack(side=LEFT, padx=5, pady=5)

        # help button
        helpButton = Button(box, command=self.help, bitmap="question", width=20)
        helpButton.pack(side=LEFT, padx=5, pady=5)

        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)

        box.pack()

    def help(self):

        tkMessageBox.showinfo("Load Previously Made Morph Help", "Load a morph from a file into the eMovie session and into the list of available morphs for insertion to the movie.  The morphs must have 30 steps/states.  If you give a morph name that is already taken, the old morph will be replaced with the new morph.")
        self.lift()


class Export(tkSimpleDialog.Dialog):

    def body(self, master):

        Label(master, text="Save image sequence as :").grid(row=0, column=0)
        Label(master, text="Ray trace frames? (y/n) :").grid(row=1, column=0)
        Label(master, text="(set GUI viewing window to the size you wish for the movie)").grid(row=2, column=0, columnspan=3)

        self.e1 = Entry(master)
        self.e2 = Entry(master)

        self.e1.grid(row=0, column=1)
        self.e2.grid(row=1, column=1)

        browseButton = Button(master, text="Browse...", command=self.openFileDialog, relief=RIDGE).grid(row=0, column=2)

        return self.e1  # initial focus

    def openFileDialog(self):

        filename = tkFileDialog.asksaveasfile(title="Browse...")

        if filename:
            self.e1.delete(0, END)  # clear the entrybox
            self.e1.insert(END, filename.name)
            self.e1.xview(END)  # scroll so you can see the end of the filename

            filename.close()  # close the file that the tkFileDialog opened for writing
            os.remove(filename.name)  # delete the file that tkFileDialog created, we will make it ourselves later

        self.lift()

    def apply(self):

        prefix = self.e1.get()
        rayTrace = self.e2.get()

        if rayTrace == "y":
            rayTrace = "on"
        else:
            rayTrace = "off"

        cmd.do("set ray_trace_frames = %s" % (rayTrace))
        cmd.do("set cache_frames = off")
        cmd.do("mpng %s" % (prefix))
        cmd.do("set ray_trace_frames = off")

    def buttonbox(self):  # overide tkSimpleDialog button box in order to insert a help button
        '''add standard button box.

        override if you do not want the standard buttons
        '''

        box = Frame(self)

        w = Button(box, text="OK", width=10, command=self.ok, default=ACTIVE)
        w.pack(side=LEFT, padx=5, pady=5)

        w = Button(box, text="Cancel", width=10, command=self.cancel)
        w.pack(side=LEFT, padx=5, pady=5)

        # help button
        helpButton = Button(box, command=self.help, bitmap="question", width=20)
        helpButton.pack(side=LEFT, padx=5, pady=5)

        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)

        box.pack()

    def help(self):

        tkMessageBox.showinfo("Export eMovie Help", "Exports the current eMovie as a .png image sequence of the given name.  Each frame is exported to one image in the sequence.  The exported image size is determined by the size of the PyMOL viewport or GUI window; make sure to set this before you export the eMovie.  Ray-tracing can be used to enhance the images, but this takes a lot of time.  Warning: Remove all stop actions from the storyboard before exporting a movie.")
        self.lift()

# a function to figure out what state the movie is on during frame


def getFrameState(frame):
    cmd.frame(frame)
    frameState = cmd.get_state()
    return frameState


def wormFunction(molecule, startAA, endAA, startFrame, state):
    # Core Worm action Script

    mv_cmd("%i:%i" % (startFrame, state), "set stick_transparency=0.0, %s,%i" % (molecule, state))
    mv_cmd("%i:%i" % (startFrame, state), "hide everything, %s" % (molecule))
    mv_cmd("%i:%i" % (startFrame, state), "cartoon tube")

    # start to show residues one by one as ball and sticks
    currentAA = startAA

    mv_cmd("%i:%i" % (startFrame, state), "set stick_ball=on; set stick_radius=0.15; set stick_ball_ratio=2.0; util.cbag %s" % (molecule))
    for a in range(startFrame + 1, startFrame + 31):
        mv_cmd("%i:%i" % (a, state), "show sticks, %s & resi %i" % (molecule, currentAA))
        currentAA = currentAA + 1

    # build protein as a backbone trace, colored as rainbow
    mv_cmd("%i:%i" % (startFrame + 32, state), "spectrum count, selection=%s&e. c" % (molecule))
    currentAA = startAA
    mv_cmd("%i:%i" % (startFrame + 33, state), "set stick_transparency = 0.5,%s,%i" % (molecule, state))
    transparency = 0.5

    # now start to fade the first 30 residues while building up the worm
    for a in range(startFrame + 34, startFrame + 64):
        mv_cmd("%i:%i" % (a, state), "show cartoon, %s & resi %i" % (molecule, currentAA))
        mv_cmd("%i:%i" % (a, state), "set stick_transparency = %f, %s,%i" % (transparency, molecule, state))
        transparency = transparency + 0.017
        currentAA = currentAA + 1

    # when the sticks are totally transparent they can be blanked
    mv_cmd("%i:%i" % (startFrame + 65, state), "hide sticks, %s" % (molecule))

    range1 = startFrame + 66 + endAA - currentAA

    for a in range(startFrame + 66, range1):
        mv_cmd("%i:%i" % (a, state), "show cartoon, %s & resi %i" % (molecule, currentAA))
        currentAA = currentAA + 1

    mv_cmd("%i:%i" % (range1 + 1, state), "set stick_transparency = 0.0, %s,%i" % (molecule, state))

    mv_turn("%i-%i:%i" % (startFrame, range1 + 1, state), "y", 360)

    return (range1 + 1)


def getStartEndFrameStates(startFrame, endFrame=0):

    blankCmd = ""

    # resolve state issue
    # 3 cases we consider here:
    # A.) startFrame is out of range of movie:
    #		1.) get last frame's state
    #		2.) set every frame's state from last frame to startFrame to last frame's state
    #		3.) return startFrameState and endFrameState both equal to last frame's state
    # B.) startFrame is in range of movie but endFrame doesn't exist
    #		1.) set startFrameState = startFrame's state and endFrameState = startFrameState
    # C.) startFrame is in range of movie and endFrame exists
    #		1.) set startFrameState = startFrame's state and endFrameState = endFrame's state
    #		2.) do command with "startFrame-endFrame:startFrameState-endFrameState" (cycling through all states btwn start and end)
    if startFrame > cmd.count_frames():
        # get state of last Frame
        startFrameState = getFrameState(cmd.count_frames())

        # if startFrame = lastFrame +1 we need no action, otherwise we need to use blank cmd to set state of blank frames until startFrame
        if startFrame > (cmd.count_frames() + 1):
            mv_cmd("%i-%i:%i" % (cmd.count_frames() + 1, startFrame - 1, startFrameState), "")
            # want to append this to storyBoard later
            blankCmd = 'mv_cmd %i-%i:%i' % (cmd.count_frames() + 1, startFrame - 1, startFrameState)

        endFrameState = startFrameState
    elif endFrame == 0:
        startFrameState = getFrameState(startFrame)
        endFrameState = startFrameState
    else:
        startFrameState = getFrameState(startFrame)
        endFrameState = getFrameState(endFrame)

    result = (startFrameState, endFrameState, blankCmd)
    return result


# function for comparing entries in storyBoard by startframe
# returns positive if a comes before b, 0 if equal, negative if b comes before a
def cmpStoryBoard(a, b):

    # check if there is no frame assigned to a or b
    if a[0] == "-":
        if b[0] == "-":
            return 0
        return -1
    elif b[0] == "-":
        return 1

    x = a[0].split("-")
    y = b[0].split("-")

    xStartFrame = int(x[0])
    yStartFrame = int(y[0])

    if xStartFrame > yStartFrame:
        return 1
    elif xStartFrame < yStartFrame:
        return -1
    else:
        # use endframes to compare

        # first get endFrames
        if len(x) == 2:
            xEndFrame = int(x[1])
        else:
            xEndFrame = xStartFrame
        if len(y) == 2:
            yEndFrame = int(y[1])
        else:
            yEndFrame = yStartFrame

        # now compare
        if xEndFrame > yEndFrame:
            return 1
        elif xEndFrame < yEndFrame:
            return -1
        else:
            return 0


def getActionValues(action):
    # takes as input an element from emovie.storyboard and returns the type of action and all the
    # important values associated with it
    #
    # Gets from (along with action type):
    #
    # SCENE SET:
    #   startFrame, sceneName
    #
    # ZOOM:
    #	startFrame, endFrame, amtZoom
    #
    # COMMAND:
    #	startFrame, inputCmd
    #
    # ROTATION:
    #	startFrame, endFrame, axis, degrees
    #
    # FADING:
    #	startFrame, endFrame, molecule, representation, fadeFrom, fadeTo
    #
    # WORM:
    #	startFrame, molecule, startAA, endAA
    #
    # PAUSE:
    #	startFrame, endFrame
    #
    # STOP:
    #	startFrame
    #
    # MORPH (FORWARD):
    #	startFrame, endFrame
    #
    # MORPH (PAUSE):
    #	startFrame, endFrame
    #
    # MORPH (BACKWARD):
    #	startFrame, endFrame

    # first we get the startFrame and the endFrame
    frames = action[0]
    t = frames.split("-")
    startFrame = t[0]  # the first token is the start frame
    if len(t) == 2:
        endFrame = t[1]
    else:
        endFrame = 0

    # next, we get the actionType
    actionType = action[1]

    # next, we get the rest of the information that we need (dependent on type of action)
    #    and return the values
    if actionType == "Pause" or actionType == "Stop":
        result = (startFrame, endFrame, actionType)
        return result
    elif actionType == "Morph":
        info = action[2]
        q = info.split(";")
        p = q[0].split(" ")
        r = q[1].split(" ")
        morphName = p[len(p) - 1]
        direction = r[len(r) - 1]

        result = (startFrame, endFrame, actionType, morphName, direction)
        return result
    elif actionType == "Scene Set":
        info = action[2]
        q = info.split(" ")
        sceneName = q[len(q) - 1]

        result = (startFrame, endFrame, actionType, sceneName)
        return result
    elif actionType == "Zoom":
        info = action[2]
        q = info.split(" ")
        if len(q) == 2:
            amtZoom = q[1]
        else:
            amtZoom = 0

        result = (startFrame, endFrame, actionType, amtZoom)
        return result
    elif actionType == "Command":
        inputCmd = action[2]

        result = (startFrame, endFrame, actionType, inputCmd)
        return result
    elif actionType == "Rotation":
        info = action[2]
        q = info.split(";")
        p = q[0].split(" ")
        r = q[1].split(" ")
        axis = p[1]
        degrees = r[2]

        result = (startFrame, endFrame, actionType, axis, degrees)
        return result
    elif actionType == "Fading":
        info = action[2]
        q = info.split(";")
        p = q[0].split(" ")
        r = q[1].split(" ")
        s = q[2].split(" ")
        u = q[3].split(" ")
        molecule = p[len(p) - 1]
        representation = r[len(r) - 1]
        fadeFrom = s[len(s) - 1]
        fadeTo = u[len(u) - 1]

        fadeFromMod = fadeFrom[:-1]  # everything but the last character which is a percent sign
        fadeToMod = fadeTo[:-1]  # everything but the last char which is a percent sign

        result = (startFrame, endFrame, actionType, molecule, representation, fadeFromMod, fadeToMod)
        return result
    elif actionType == "Worm":
        info = action[2]
        q = info.split(";")
        p = q[0].split(" ")
        r = q[1].split(" ")
        s = r[-1].split("-")
        molecule = p[len(p) - 1]
        startAA = s[0]
        endAA = s[1]

        result = (startFrame, endFrame, actionType, molecule, startAA, endAA)
        return result
    else:
        print("Error: ActionType not found in getActionValues function.")
        return (0)


def insertActionByValues(actionValues):

    # this function takes the actionValues in the format of those obtained by the function getActionValues
    # and inserts an action into the movie with the values given
    #
    # most of the code is copied from the "apply" functions of the individual action window functions
    # and altered accordingly

    if actionValues[2] == "Pause":

        startFrame = int(actionValues[0])
        endFrame = int(actionValues[1])

        (startFrameState, endFrameState, blankCmd) = getStartEndFrameStates(startFrame, endFrame)

        # add pause to movie
        mv_cmd("%i-%i:%i-%i" % (startFrame, endFrame, startFrameState, endFrameState), "")

        # append to storyBoard
        emovie.storyBoard.append(("%i-%i" % (startFrame, endFrame), "Pause", "-", ['mv_cmd %i-%i:%i-%i' % (startFrame, endFrame, startFrameState, endFrameState), blankCmd]))

        movie()
        cmd.do("mstop")
    elif actionValues[2] == "Stop":

        startFrame = int(actionValues[0])

        (startFrameState, endFrameState, blankCmd) = getStartEndFrameStates(startFrame)

        # add pause to movie
        mv_cmd("%i:%i" % (startFrame, startFrameState), "mstop")

        # append to storyBoard
        emovie.storyBoard.append(("%i" % (startFrame), "Stop", "-", ['mv_cmd %i:%i,mstop' % (startFrame, startFrameState), blankCmd]))

        movie()
        cmd.do("mstop")
    elif actionValues[2] == "Morph":

        morphName = actionValues[3]
        startFrame = int(actionValues[0])
        endFrame = int(actionValues[1])
        direction = actionValues[4]

        blankCmd = ""

        # need to get the start and end states of the morph
        for a in emovie.morphList:
            if a[0] == morphName:
                startState = a[1]
                endState = a[2]

        # need to adjust for states present in movie
        if startFrame > cmd.count_frames():
            # get state of last Frame, use blank cmd to set state of blank frames until start Frame
            lastState = getFrameState(cmd.count_frames())
            mv_cmd("%i-%i:%i" % (cmd.count_frames(), startFrame - 1, lastState), "")
            blankCmd = 'eMovie.mv_cmd("%i-%i:%i","")' % (cmd.count_frames(), startFrame - 1, lastState)

        if direction == "pause":

            # add command to movie
            mv_cmd("%i-%i:%i" % (startFrame, endFrame, endState), "")

            # append to storyBoard
            emovie.storyBoard.append(("%i-%i" % (startFrame, endFrame), "Morph", "Name: %s; Action: pause" % (morphName), ['mv_cmd %i-%i:%i' % (startFrame, endFrame, endState)]))

        elif direction == "backwards":

            # add command to movie
            mv_cmd("%i-%i:%i-%i" % (startFrame, startFrame + 29, endState, startState), "")

            # append to storyBoard
            emovie.storyBoard.append(("%i-%i" % (startFrame, startFrame + 29), "Morph", "Name: %s; Direction: %s" % (morphName, direction), ['mv_cmd %i-%i:%i-%i' % (startFrame, startFrame + 29, endState, startState), blankCmd]))

        else:  # direction = forward

            # add command to movie
            mv_cmd("%i-%i:%i-%i" % (startFrame, startFrame + 29, startState, endState), "")

            # append to storyBoard
            emovie.storyBoard.append(("%i-%i" % (startFrame, startFrame + 29), "Morph", "Name: %s; Direction: %s" % (morphName, direction), ['mv_cmd %i-%i:%i-%i' % (startFrame, startFrame + 29, startState, endState), blankCmd]))

        movie()
        cmd.do("mstop")
    elif actionValues[2] == "Scene Set":

        sceneName = actionValues[3]
        startFrame = int(actionValues[0])

        (startFrameState, endFrameState, blankCmd) = getStartEndFrameStates(startFrame)

        # add command to movie
        mv_cmd("%i:%i" % (startFrame, startFrameState), 'cmd.scene("%s","recall")' % (sceneName))

        # append to storyboard
        emovie.storyBoard.append(("%i" % (startFrame), "Scene Set", "Scene name: %s" % (sceneName), ['eMovie.mv_cmd("%i:%i","cmd.scene(\'%s\',\'recall\')")' % (startFrame, startFrameState, sceneName), blankCmd]))

        movie()
        cmd.do("mstop")
    elif actionValues[2] == "Zoom":

        amtZoom = int(actionValues[3])
        startFrame = int(actionValues[0])
        endFrame = int(actionValues[1])

        (startFrameState, endFrameState, blankCmd) = getStartEndFrameStates(startFrame, endFrame)

        # add command to movie
        mv_move("%i-%i:%i-%i" % (startFrame, endFrame, startFrameState, endFrameState), "z", "%i" % (amtZoom), "linear")

        # append to storyBoard
        emovie.storyBoard.append(("%i-%i" % (startFrame, endFrame), "Zoom", "Angstroms: %i" % (amtZoom), ['mv_move %i-%i:%i-%i,z,%i,linear' % (startFrame, endFrame, startFrameState, endFrameState, amtZoom), blankCmd]))

        movie()
        cmd.do("mstop")
    elif actionValues[2] == "Command":

        inputCmd = actionValues[3]
        startFrame = int(actionValues[0])

        (startFrameState, endFrameState, blankCmd) = getStartEndFrameStates(startFrame)

        # add command to movie
        mv_cmd("%i:%i" % (startFrame, startFrameState), inputCmd)

        # append to storyBoard
        emovie.storyBoard.append(("%i" % (startFrame), "Command", inputCmd, ['eMovie.mv_cmd("%i:%i","%s")' % (startFrame, startFrameState, inputCmd), blankCmd]))

        movie()
        cmd.do("mstop")
    elif actionValues[2] == "Rotation":

        startFrame = int(actionValues[0])
        endFrame = int(actionValues[1])
        axis = actionValues[3]
        degrees = actionValues[4]

        (startFrameState, endFrameState, blankCmd) = getStartEndFrameStates(startFrame, endFrame)

        # add command to movie
        mv_turn("%i-%i:%i-%i" % (startFrame, endFrame, startFrameState, endFrameState), axis, degrees, "linear")

        # append to storyBoard
        emovie.storyBoard.append(("%i-%i" % (startFrame, endFrame), "Rotation", "Axis: %s; Degrees: %s" % (axis, degrees), ['mv_turn %i-%i:%i-%i,%s,%s,linear' % (startFrame, endFrame, startFrameState, endFrameState, axis, degrees), blankCmd]))

        movie()
        cmd.do("mstop")
    elif actionValues[2] == "Fading":

        startFrame = int(actionValues[0])
        endFrame = int(actionValues[1])
        molecule = actionValues[3]
        representation = actionValues[4]
        fadeFrom = float(actionValues[5])
        fadeTo = float(actionValues[6])

        (startFrameState, endFrameState, blankCmd) = getStartEndFrameStates(startFrame, endFrame)

        # convert fadeFrom and fadeTo from percentage to fractions of a whole
        fadeFromMod = fadeFrom / 100
        fadeToMod = fadeTo / 100
        # convert to amount of transparency rather than amt visibility
        fadeFromMod = 1 - fadeFromMod
        fadeToMod = 1 - fadeToMod

        if startFrameState < endFrameState:
            anteEndFrameState = endFrameState - 1
        elif startFrameState > endFrameState:
            anteEndFrameState = endFrameState + 1
        else:
            anteEndFrameState = endFrameState

        # since surface transparency setting is transparency and all others are of the form TYPE_transparency, we do:
        if representation == "surface":
            appearance = ""
        else:
            appearance = "%s_" % (representation)

        # add command to movie
        mv_set("%i-%i:%i-%i" % (startFrame, endFrame - 1, startFrameState, anteEndFrameState), "%stransparency" % (appearance), "%f" % (fadeFromMod), "%f" % (fadeToMod), molecule, 'linear')
        # because there is a issue with mv_set's linear setting of values, we force the last value on the last frame for now
        mv_cmd("%i:%i" % (endFrame, endFrameState), "set %stransparency,%f,%s" % (appearance, fadeToMod, molecule))

        # append to storyBoard
        emovie.storyBoard.append(("%i-%i" % (startFrame, endFrame), "Fading", "Molecule: %s; Shown as: %s; Fade from: %.2f%%; Fade to: %.2f%%" % (molecule, representation, fadeFrom, fadeTo), ['mv_set %i-%i:%i-%i,%stransparency,%f,%f,%s,linear' % (startFrame, endFrame - 1, startFrameState, anteEndFrameState, appearance, fadeFromMod, fadeToMod, molecule), 'mv_cmd %i:%i,cmd.set("%stransparency","%f","%s")' % (endFrame, endFrameState, appearance, fadeToMod, molecule), blankCmd]))

        movie()
        cmd.do("mstop")
    elif actionValues[2] == "Worm":

        startFrame = int(actionValues[0])
        molecule = actionValues[3]
        startAA = int(actionValues[4])
        endAA = int(actionValues[5])

        (startFrameState, endFrameState, blankCmd) = getStartEndFrameStates(startFrame)

        lastFrame = wormFunction(molecule, startAA, endAA, startFrame, startFrameState)

        movie()
        cmd.do("mstop")

        emovie.storyBoard.append(("%i-%i" % (startFrame, lastFrame), "Worm", "Molecule: %s; Residues: %i-%i" % (molecule, startAA, endAA), ['eMovie.wormFunction("%s",%i,%i,%i,%i)' % (molecule, startAA, endAA, startFrame, startFrameState), blankCmd]))

    else:
        print("Error: ActionType not found in insertActionByValues function.")
        return (0)


# Here I have added code from Kristian Rother's "movie.py" (with some small edits)

#"""
#--- movie: easy movie scripting in PyMOL ---
# Script  : movie
# Author  : The PyMOL community
# Date    : Oct 2004
# Version : 0.7
# Contact : kristian.r@gmx.de
#
# Copyright (C) 2004 Kristian Rother, Laurence Pearl, Seth Harris,
# Morri Feldmann, Tserk Wassenaar and Lieven Buts
#
# Movie is an extension for the PyMOL Molecular Graphics System.
#
# Movie for PyMOL facilitates creation of high-quality animations. Translation, rotation
# and other commands can be assigned to frame ranges, resulting in a much shorter scripting
# of animations. Additionaly, smooth movements are enabled.
#
# It is assumed that you have already gathered experience with the standard PyMOL animation
# commands, or are at least familiar with the 'frame' and 'state' concepts. Most movie
# commands will require a frame/state range assigned to a specific action. Thus, you
# dont have to worry about how to use the mset command correctly. Just say what you
# want to see.
#
# Frames can be specified to cover individual frames, states or ranges of both. The frame
# ranges of several commands may very well overlap (but not the states).
# The 'mode' parameter of several commands takes the values 'linear' and 'trigon'
#
# Examples for frames/states:
# 10         # frame 10 only
# 20-199     # frame 20 to 199
# 1-100:5    # state 5 only
# 1-100:1-10 # multiple states are looped through
#
#"""
class Movie:

    def __init__(self):
        """Stores data of what should appear in the movie."""
        self.movie = []
        self.maxframe = 1
        self.framestates = {}
        self.png_prefix = ""
        self.views = {}
        self.ray = 0

    def add(self, framestate, command):
        self.movie.append((framestate[0], framestate[1], command))

    def add_frame_states(self, framestates):
        """Stores list of (frame,state) tuples in a dictionary."""
        for fs in framestates:
            self.framestates[str(fs[0])] = fs[1]

    def get_state(self, frame):
        key = str(frame)
        if key in self.framestates:
            return self.framestates[key]
        else:
            return 1  # state one, if none specified


moviedata = Movie()


# ---------------------------------------------------------
# internal stuff

def get_linear_values(number, start=0.0, end=1.0):
    """
    Returns a list of 'number' float values, where the first
    value is 'start', the last one is 'end' and the others
    are linearly in between them.
    """
    if number == 1:
        return [end]
    values = []
    i = 0.0
    while i <= number:
        values.append(start + (end - start) * i / (number))
        i += 1
    return values


def get_trigon_values(number, start=0.0, end=1.0):
    """
    Returns a list of 'number' float values, where the first
    value is 'start', the last one is 'end' and the others
    are calculated using a trigonometric function.
    """
    if number == 1:
        return [end]
    values = []
    arcIncrement = pi / (number - 1)

    prev = 1.0
    sum = 0.0
    i = 0.0
    while i < number:
        arc = cos(i * arcIncrement)
        sum += (end - start) * abs(arc - prev) * 0.5
        values.append(start + sum)
        prev = arc
        i += 1
    return values


def get_values(number, start, end, mode):
    if mode == 'linear':
        values = get_linear_values(number, start, end)
    elif mode == 'trigon':
        values = get_trigon_values(number, start, end)
    else:
        print('movie.py: INVALID MODE %s' % mode)
        return []

    return values


def get_increment_values(number, start=0.0, end=1.0, mode='linear'):
    """
    Returns a list of float numbers that are either linear or
    trigonometric increments that sum up to the difference between
    start and end.
    """
    if mode == 'linear':
        incr = []
        increment = (end - start) / number
        incr = [increment] * number
        return incr
    else:
        values = get_values(number + 1, start, end, mode)
        incr = []
        act = values[0]
        for v in values[1:]:
            incr.append(v - act)
            act = v
        return incr


def get_frame_states(fstring):
    """Returns a list of (frame,state) tuples parsed from a string like
    1-10       # frame 1 to 10, state 1
    10         # frame 10 only
    1-100:1-10 # multiple states are looped through
    1-100:5    # one state only
    1-50:10-5  # reverse order is ok, too
    1:8        # state no.8 in frame no.1
    1:1-10     # this is crap, don't know what happens
    """

    t = fstring.split(":")
    frames = t[0]  # the first token is frame range
    if len(t) == 2:
        states = t[1]  # the second token is state range
    else:
        states = '1'  # only state number 1 is used

    # parse frame substring
    t = frames.split("-")
    firstFrame = int(t[0])  # first token is starting frame
    if len(t) == 2:
        lastFrame = int(t[1])  # second token is end frame
    else:
        lastFrame = firstFrame  # only one frame is used

    # parse state substring
    t = states.split("-")
    firstState = int(t[0])  # first token is starting state
    if len(t) == 2:
        lastState = int(t[1])  # second token is end state
    else:
        lastState = firstState  # only one state is used

    # compile list of frames/states
    framestates = []
    if lastFrame == firstFrame or lastFrame < firstFrame:
        framestates.append((firstFrame, firstState))
        lastFrame = firstFrame
    else:
        nframes = lastFrame - firstFrame + 1
        if lastState >= firstState:
            stateinc = (lastState - firstState + 1) * 1.0 / nframes
        else:
            stateinc = (lastState - firstState - 1) * 1.0 / nframes
        for i in range(nframes):
            frame = firstFrame + i
            state = firstState + int(stateinc * i)
            framestates.append((frame, state))

    # put values into mv for compiling the movie later
    if moviedata.maxframe < lastFrame:
        moviedata.maxframe = lastFrame
    moviedata.add_frame_states(framestates)

    return framestates

#
#
# Commands revised
#
#


def mv_clear():
    """Deletes the movie."""
    moviedata.movie = []
    cmd.mset("1")
    cmd.mclear()
    cmd.frame(1)
    moviedata.maxframe = 1
    moviedata.framestates = {}


def mv_cmd(frames="1", command=""):
    """
    mv_cmd(frames,command) - executes a command in all frames specified.
    """
    framestates = get_frame_states(frames)
    if command != "":
        for fs in framestates:
            moviedata.add(fs, command)


def mv_turn(frames="1", axis="z", angle='360', mode='linear'):
    """
    mv_turn(frames,axis,angle,selection,mode) - turns the camera over the given
            frame range, the rotation angle summing up to the angle given.
            """
    framestates = get_frame_states(frames)
    nFrames = len(framestates)
    angleIncrement = get_increment_values(nFrames, 0.0, float(angle), mode)

    for i in range(nFrames):
        moviedata.add(framestates[i], "turn %s,%f" % (axis, angleIncrement[i]))


def mv_rotate(frames="1", axis="z", angle='360', selection='all', mode='linear'):
    """
    mv_rotate(frames,axis,angle,selection,mode) - rotates the object over the given
            frame range, the rotation angle summing up to the angle given.
            """
    framestates = get_frame_states(frames)
    nFrames = len(framestates)
    angleIncrement = get_increment_values(nFrames, 0.0, float(angle), mode)

    for i in range(nFrames):
        moviedata.add(framestates[i], "rotate %s,%f,%s" % (axis, angleIncrement[i], selection))


def mv_move(frames="1", axis="x", distance="0", mode='linear'):
    """
    mv_move(frames,axis,distance,selection,mode) - moves the environment over the given
            frame range, the moved distance summing up to the distance given.
            """
    framestates = get_frame_states(frames)
    nFrames = len(framestates)
    distanceIncrement = get_increment_values(nFrames, 0.0, float(distance), mode)

    for i in range(nFrames):
        moviedata.add(framestates[i], "move %s,%f" % (axis, distanceIncrement[i]))


def mv_trans(frames="1", axis="x", distance="0", selection='all', mode='linear'):
    """
    mv_move(frames,axis,distance,selection,mode) - moves the selection object over the given
            frame range, the moved distance summing up to the distance given.
            """
    framestates = get_frame_states(frames)
    nFrames = len(framestates)
    distanceIncrement = get_increment_values(nFrames, 0.0, float(distance), mode)

    for i in range(nFrames):
        if axis == 'x':
            vector = "[%f,0,0]" % (distanceIncrement[i])
        elif axis == 'y':
            vector = "[0,%f,0]" % (distanceIncrement[i])
        elif axis == 'z':
            vector = "[0,0,%f]" % (distanceIncrement[i])
        moviedata.add(framestates[i], "translate %s,%s" % (vector, selection))


def mv_set(frames="1", variable="", start="0.0", end="1.0", selection='all', mode='linear'):
    """
    mv_set(frames,variable,start,end,selection,mode) - lets a PyMOL variable go through a gradient
    in the specified frame range. Great for fading effects!
    """
    framestates = get_frame_states(frames)
    nFrames = len(framestates)
    values = get_values(nFrames, float(start), float(end), mode)

    for i in range(nFrames):
        moviedata.add(framestates[i], "set %s,%f,%s" % (variable, values[i], selection))


def movie(png_prefix="", ray="0"):
    """
    Creates the movie and plays it.
    """
    if png_prefix:
        moviedata.png_prefix = png_prefix
    if ray != "0":
        moviedata.ray = 1
    if moviedata.png_prefix != "":
        # stop movie after first loop
        mv_cmd(str(moviedata.maxframe + 1), "mstop")
        mv_cmd(str(moviedata.maxframe + 2), "dummy")
        print("done")

    cmd.frame(1)  # reset frame counter
    nFrames = moviedata.maxframe

    # compile list of molecule states
    states = []
    for n in range(nFrames):
        states.append(moviedata.get_state(n + 1))

    # compile mset command string "1 2 3 4 x3 5" from state list [1,2,3,4,4,4,5]
    statelist = ""
    count = 0
    for e in range(nFrames):
        actual = states[e]
        next = 0
        if e < nFrames - 1:
            next = states[e + 1]
        if next == actual:
            count += 1
        else:
            if count > 0:
                statelist += "%s x%i " % (actual, count + 1)
            else:
                statelist += "%s " % (actual)
            count = 0

    cmd.mset(statelist)

    do = ["zero unused"]  # create empty frame-2do-list
    for i in range(nFrames):
        do.append("")

    for m in moviedata.movie:
        do[m[0]] += m[2] + ";"  # push all movie commands to the 2do-list

    # check for png output and raytracing:
    # if moviedata.png_prefix!="":
    #    cmd.set('ray_trace_frames',moviedata.ray)
    #    cmd.set('cache_frames',0)
    #    cmd.mpng(moviedata.png_prefix)
    if moviedata.png_prefix != "":
        for i in range(nFrames):
            if moviedata.ray:
                do[i] += "ray;"
            num = "0000" + str(i)
            num = num[-4:]
            do[i] += "png %s.%s;" % (moviedata.png_prefix, num)

    # now let action happen in the frames
    for i in range(nFrames):
        cmd.mdo(i + 1, do[i + 1])

    cmd.mplay()  # start the movie
