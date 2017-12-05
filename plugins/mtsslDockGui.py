"""
#-----------------------------------------------------------------------------------------
mtsslDock: distance constrained rigid body docking in PyMOL
Author	: Gregor Hagelueken
Date	: July 2013
Version : 1.0
Mail	: hagelueken'at'pc.uni-bonn.de
#-----------------------------------------------------------------------------------------
"""
import pymol
import numpy
import copy
import threading
import scipy.spatial.distance
import random
import time
import math
import os
from pymol import cmd
from pymol import util
from pymol import stored
from operator import attrgetter
import sys

if sys.version_info[0] < 3:
    import Tkinter
    import tkFileDialog
    import Queue
    import ttk
else:
    import tkinter as Tkinter
    import tkinter.filedialog as tkFileDialog
    import queue as Queue
    import tkinter.ttk as ttk

import Pmw
from tkintertable.Tables import TableCanvas
from tkintertable.TableModels import TableModel
import pickle


def __init__(self):
    self.menuBar.addmenuitem('Plugin', 'command',
                             'Launch Plugin',
                             label='mtsslDock',
                             command=lambda s=self: mtsslDockPlugin(s))

##########################################################################################
# custom GUI classes
##########################################################################################

#-----------------------------------------------------------------------------------------


class ImportDialog(Pmw.Dialog):

    def __init__(self, parent=None, **kw):
        # Initialise base class
        Pmw.Dialog.__init__(self, parent, **kw)
        interior = self.interior()
        self.importGroup = self.createcomponent('importGroup', (), None, Pmw.Group, (interior,), tag_text='Import')
        #importGroup = Pmw.Group(interior, tag_text='Import')
        self.importGroup.pack(fill='both', expand=2, padx=10, pady=5)

        # protein selection group
        #proteinSelectionGroup = Pmw.Group(self.importGroup.interior(), tag_text='Select proteins to be docked')
        self.proteinSelectionGroup = self.createcomponent('proteinSelectionGroup', (), None, Pmw.Group, (self.importGroup.interior(),), tag_text='Select proteins to be docked')
        self.proteinSelectionGroup.pack(fill='x', anchor='n', padx=10, pady=5)
        self.proteinASelector = self.createcomponent('proteinASelector', (), None, Pmw.OptionMenu, (self.proteinSelectionGroup.interior(),), labelpos='w', label_text='Protein A', items=([]),)
        self.proteinBSelector = self.createcomponent('proteinBSelector', (), None, Pmw.OptionMenu, (self.proteinSelectionGroup.interior(),), labelpos='w', label_text='Protein B', items=([]),)
        items = cmd.get_object_list('(all)')
        self.proteinASelector.setitems(items)
        self.proteinBSelector.setitems(items)
        entries = (self.proteinASelector, self.proteinBSelector)
        for entry in entries:
            entry.pack(fill='x', anchor='s', padx=4, pady=1)  # vertical

        # assignLabels Group
        self.assignLabelsGroup = self.createcomponent('assignLabelsGroup', (), None, Pmw.Group, (self.importGroup.interior(),), tag_text='Assign labels')
        self.assignLabelsGroup.pack(fill='both', expand=1, anchor='n', padx=10, pady=5)
        self.assignLabelsGroupScrolledFrame = self.createcomponent('assignLabelsGroupScrolledFrame', (), None, Pmw.ScrolledFrame, (self.assignLabelsGroup.interior(),),)
        self.assignLabelsGroupScrolledFrame.pack(fill='both', expand=1, anchor='s', padx=10, pady=5)
        items = []
        items = cmd.get_object_list('(*_label)')
        self.listOfLabelSelectors = []
        for pymolObject in items:
            labelSelector = ABselect(self.assignLabelsGroupScrolledFrame.interior())
            labelSelector.pack(fill='x', anchor='n', padx=4, pady=1)
            labelSelector.labelText.set(pymolObject)
            self.listOfLabelSelectors.append(labelSelector)
        self.initialiseoptions()

#-----------------------------------------------------------------------------------------


class SettingsDialog(Pmw.Dialog):

    def __init__(self, parent=None, **kw):
        # Initialise base class (after defining options).
        Pmw.Dialog.__init__(self, parent, **kw)
        interior = self.interior()
        self.scoreClashes = True
        #self.c2 = False
        self.symmetry = "None"
        self.settingsGroup = self.createcomponent('settingsGroup', (), None, Pmw.Group, (interior,), tag_text='Settings')
        self.settingsGroup.pack(fill='x', padx=10, pady=5)
        self.numberOfPopulations = self.createcomponent('numberOfPopulations', (), None, Pmw.EntryField, (self.settingsGroup.interior(),), labelpos='w', label_text='Populations', value=str(10), validate={'validator': 'integer', 'min': 1},)
        self.numberOfGenerations = self.createcomponent('numberOfGenerations', (), None, Pmw.EntryField, (self.settingsGroup.interior(),), labelpos='w', label_text='Generations', value=str(1000), validate={'validator': 'integer', 'min': 1},)
        self.numberOfChromosomes = self.createcomponent('numberOfChromosomes', (), None, Pmw.EntryField, (self.settingsGroup.interior(),), labelpos='w', label_text='Chromosomes', value=str(400), validate={'validator': 'integer', 'min': 1},)
        self.numberOfRigidBodyCycles = self.createcomponent('numberOfRigidBodyCycles', (), None, Pmw.EntryField, (self.settingsGroup.interior(),), labelpos='w', label_text='Rigid body cycles', value=str(50), validate={'validator': 'integer', 'min': 1},)
        # Create and pack a vertical RadioSelect widget, with checkbuttons.
        self.checkbuttons = self.createcomponent('checkbuttons', (), None, Pmw.RadioSelect, (self.settingsGroup.interior(),), buttontype='checkbutton', orient='vertical', labelpos='w', command=self.checkbuttoncallback, label_text='', hull_borderwidth=2, hull_relief='ridge',)

        # Add some buttons to the checkbutton RadioSelect
        self.checkbuttons.add("Score clashes")
        self.checkbuttons.add("Force C2 symmetry")
        self.checkbuttons.invoke("Score clashes")
        entries = (self.numberOfPopulations, self.numberOfGenerations, self.numberOfChromosomes, self.numberOfRigidBodyCycles, self.checkbuttons)
        for entry in entries:
            entry.pack(fill='x', expand=1, anchor='n', padx=4, pady=1)  # vertical
        self.initialiseoptions()

    def checkbuttoncallback(self, tag, state):
        if tag == "Score clashes":
            if state:
                self.scoreClashes = True
            else:
                self.scoreClashes = False
        if tag == "Force C2 symmetry":
            if state:
                #self.c2 = True
                self.symmetry = "C2"
            else:
                #self.c2 = False
                self.symmetry = "None"

#-----------------------------------------------------------------------------------------


class ABselect(Tkinter.Frame):

    def __init__(self, parent, **options):
        Tkinter.Frame.__init__(self, parent, **options)
        self.checkbuttons = Pmw.RadioSelect(self,
                                            buttontype='radiobutton',
                                            orient='horizontal',
                                            labelpos='w',
                                            label_text='')

        for text in ('A', 'B'):
            self.checkbuttons.add(text)
        self.labelText = Tkinter.StringVar()
        self.labelText.set("ObjectName")
        self.objectLabel = Tkinter.Label(self, textvariable=self.labelText)
        self.checkbuttons.pack(side=Tkinter.LEFT)
        self.objectLabel.pack(side=Tkinter.LEFT)

#-----------------------------------------------------------------------------------------


class DockingFrame(Tkinter.Frame):

    def __init__(self, parent, result, **options):
        Tkinter.Frame.__init__(self, parent, **options)

        # protein selection group
        proteinSelectionGroup = Pmw.Group(self, tag_text='Select proteins to be docked')
        self.proteinSelectionGroup.pack(fill='x', anchor='n', padx=10, pady=5)
        self.proteinASelector = Pmw.OptionMenu(self, labelpos='w', label_text='Protein A', items=([]),)
        self.proteinBSelector = Pmw.OptionMenu(self, labelpos='w', label_text='Protein B', items=([]),)
        items = cmd.get_object_list('(all)')
        self.proteinASelector.setitems(items)
        self.proteinBSelector.setitems(items)
        entries = (self.proteinASelector, self.proteinBSelector)
        for entry in entries:
            entry.pack(fill='x', anchor='s', padx=4, pady=1)  # vertical

        self.assignLabelsGroupScrolledFrame = self.createcomponent('assignLabelsGroupScrolledFrame', (), None, Pmw.ScrolledFrame, (self.assignLabelsGroup.interior(),),)
        self.assignLabelsGroupScrolledFrame.pack(fill='both', expand=1, anchor='s', padx=10, pady=5)
        items = []
        items = cmd.get_object_list('(*_label)')
        self.listOfLabelSelectors = []
        for pymolObject in items:
            labelSelector = ABselect(self.assignLabelsGroupScrolledFrame.interior())
            labelSelector.pack(fill='x', anchor='n', padx=4, pady=1)
            labelSelector.labelText.set(pymolObject)
            self.listOfLabelSelectors.append(labelSelector)
#-----------------------------------------------------------------------------------------


class ResultsFrame(Tkinter.Frame):

    def __init__(self, parent, result, **options):
        Tkinter.Frame.__init__(self, parent, **options)
        self.result = result

        self.Frame1 = Tkinter.Frame(self)
        self.Frame1.pack(fill='both', side='right', expand=1, anchor='nw', padx=10, pady=5)

        self.Frame2 = Tkinter.Frame(self)
        self.Frame2.pack(fill='both', side='left', anchor='ne', padx=10, pady=5)

        #########
        # Setup listbox
        #########
        self.box = Pmw.ScrolledListBox(self.Frame2,
                                       labelpos='nw',
                                       label_text='Solutions:',
                                       selectioncommand=self.selectionCommand,
                                       usehullsize=1,
                                       hull_width=200,
                                       hull_height=200,
                                       )
        self.box.pack(fill='both', expand=1, anchor='nw', padx=10, pady=5)
        for population in result['environment'].populations:
            print(population.chromosomes[0].name)
            self.box.insert('end', population.chromosomes[0].name)

        #########
        # Setup results group
        #########
        self.resultsGroup = Pmw.Group(self.Frame1, tag_text='Results:')
        self.resultsGroup.pack(fill='both', side='right', expand=1, anchor='ne', padx=10, pady=5)

        self.createLogTable('rotX rotY rotZ tX tY tZ rmsd chi2 clashes')

        #########
        # Setup radioselect
        #########
        self.radioSelect = Pmw.RadioSelect(self.resultsGroup.interior(),
                                           command=self.callback)
        self.radioSelect.pack(fill='x', side='bottom', padx=10, pady=10)

        # Add some buttons to the horizontal RadioSelect.
        for text in ('Constraints', 'Population log', 'Solution log'):
            self.radioSelect.add(text)
        self.radioSelect.invoke('Population log')

        #########
        # Setup settings group
        #########
        self.settingsGroup = Pmw.Group(self.Frame2, tag_text='Settings for this run:')
        self.settingsGroup.pack(fill='both', anchor='sw', padx=10, pady=5)
        populations = "Populations: %i" % self.result['settings']['numberOfPopulations']
        generations = "Generations: %i" % self.result['settings']['numberOfGenerations']
        chromosomes = "Chromosomes: %i" % self.result['settings']['numberOfChromosomes']
        rigidBodyCycles = "Rigid Body Cycles: %i" % self.result['settings']['numberOfRigidBodyCycles']
        scoreClashes = "Score Clashes: %r" % self.result['settings']['scoreClashes']
        symmetry = "Symmetry: %s" % self.result['settings']['symmetry']

        populationsLabel = Tkinter.Label(self.settingsGroup.interior(), text=populations)
        populationsLabel.pack(anchor='w')

        generationsLabel = Tkinter.Label(self.settingsGroup.interior(), text=generations)
        generationsLabel.pack(anchor='w')

        chromosomesLabel = Tkinter.Label(self.settingsGroup.interior(), text=chromosomes)
        chromosomesLabel.pack(anchor='w')

        rigidBodyCyclesLabel = Tkinter.Label(self.settingsGroup.interior(), text=rigidBodyCycles)
        rigidBodyCyclesLabel.pack(anchor='w')

        scoreClashesLabel = Tkinter.Label(self.settingsGroup.interior(), text=scoreClashes)
        scoreClashesLabel.pack(anchor='w')

        symmetryLabel = Tkinter.Label(self.settingsGroup.interior(), text=symmetry)
        symmetryLabel.pack(anchor='w')

#-----------------------------------------------------------------------------------------

    def selectionCommand(self):
        if self.radioSelect.getcurselection() == 'Solution log':
            self.displaySolutionLog()
        if self.radioSelect.getcurselection() == 'Constraints':
            self.displayConstraints()
        if self.radioSelect.getcurselection() == 'Population log':
            self.displayPopulationLog()

#-----------------------------------------------------------------------------------------

    def callback(self, tag):
        if tag == 'Solution log':
            self.displaySolutionLog()
        if tag == 'Population log':
            self.displayPopulationLog()
        if tag == 'Constraints':
            self.displayConstraints()


#-----------------------------------------------------------------------------------------

    def displaySolutionLog(self):
        if len(self.box.getcurselection()) > 0:
            self.removeLogTable()
            self.createLogTable('rotX rotY rotZ tX tY tZ rmsd chi2 clashes')
            selectedSolutionName = self.box.getcurselection()[0]
            selectedSolution = Chromosome("None")
            # find selected solution
            for population in self.result['environment'].populations:
                if population.chromosomes[0].name == selectedSolutionName:
                    selectedSolution = population.chromosomes[0]
            cmd.zoom(selectedSolution.name)

            # print loglines
            loglines = selectedSolution.log.splitlines()

            # datalines
            for logline in loglines:
                dataline = ''
                logColumns = logline.split()
                for column in logColumns:
                    data = '%-7s' % (column)
                    dataline += data + '   '
                dataline = dataline[:-3]
                dataline = dataline + '\n'
                self.logTable.insert('end', dataline)

            self.logTable.configure(text_state='disabled', Header_state='disabled')

#-----------------------------------------------------------------------------------------

    def displayPopulationLog(self):
        if len(self.box.getcurselection()) > 0:
            self.removeLogTable()
            self.createLogTable('rotX rotY rotZ tX tY tZ rmsd chi2 clashes')
            selectedSolutionName = self.box.getcurselection()[0]
            selectedSolution = Chromosome("None")
            # find selected solution
            for population in self.result['environment'].populations:
                if population.chromosomes[0].name == selectedSolutionName:
                    selectedSolution = population.chromosomes[0]
                    selectedPopulation = population
            cmd.zoom(selectedSolution.name)

            # print loglines
            loglines = selectedPopulation.log.splitlines()

            # datalines
            for logline in loglines:
                dataline = ''
                logColumns = logline.split()
                for column in logColumns:
                    data = '%-7s' % (column)
                    dataline += data + '   '
                dataline = dataline[:-3]
                dataline = dataline + '\n'
                self.logTable.insert('end', dataline)

            self.logTable.configure(text_state='disabled', Header_state='disabled')


#-----------------------------------------------------------------------------------------

    def removeLogTable(self):
        self.logTable.pack_forget()

#-----------------------------------------------------------------------------------------

    def createLogTable(self, headers):
        #########
        # Setup table
        #########
        fixedFont = Pmw.logicalfont('Fixed')
        self.logTable = Pmw.ScrolledText(self.resultsGroup.interior(),
                                         # borderframe = 1,
                                         labelpos='nw',
                                         label_text='Log:',
                                         columnheader=1,
                                         usehullsize=1,
                                         hull_width=400,
                                         hull_height=300,
                                         text_wrap='none',
                                         text_font=fixedFont,
                                         Header_font=fixedFont,
                                         Header_foreground='blue',
                                         )
        self.logTable.pack(fill='both', side='top', expand=1, anchor='ne', padx=10, pady=5)
        # columnheaders
        headers = headers.split()
        headerLine = ''
        for column in range(len(headers)):
            headerLine = headerLine + ('%-7s   ' % (headers[column],))
        headerLine = headerLine[:-3]
        self.logTable.component('columnheader').insert('0.0', headerLine)


#-----------------------------------------------------------------------------------------

    def displayConstraints(self):
        if len(self.box.getcurselection()) > 0:
            self.removeLogTable()
            self.createLogTable('Exp StdDev Docked Diff')

            selectedSolutionName = self.box.getcurselection()[0]
            for population in self.result['environment'].populations:
                if population.chromosomes[0].name == selectedSolutionName:
                    selectedSolution = population.chromosomes[0]
                    selectedPopulation = population

            try:
                cmd.zoom(selectedSolution.name)
            except:
                print("Can't find that object in PyMOL ...")

            rowheaders = self.result['environment'].constraintNames
            expDistances = self.result['environment'].expDistances
            expErrors = self.result['environment'].expErrors
            trialDistances = selectedSolution.trialDistances

            for distance, expError, trialDistance in zip(expDistances, expErrors, trialDistances):
                if distance != 0:
                    dataline = ''
                    data = '%-.*f' % (5, distance)
                    dataline += data + '  '
                    data = '%-.*f' % (5, expError)
                    dataline += data + '  '
                    data = '%-.*f' % (5, trialDistance)
                    dataline += data + '  '
                    data = '%-.*f' % (5, distance - trialDistance)
                    dataline += data + '  '
                    dataline = dataline + '\n'
                    self.logTable.insert('end', dataline)
            self.logTable.configure(text_state='disabled', Header_state='disabled')

#-----------------------------------------------------------------------------------------


class StatusBar(Tkinter.Frame):

    def __init__(self, parent, **options):
        Tkinter.Frame.__init__(self, parent, **options)
        #########
        # Setup messagebar
        #########
        self.messagebar = Pmw.MessageBar(self, entry_width=40, entry_relief='sunken', labelpos='w', label_text='Status:')
        self.messagebar.pack(side='left', fill='x')

        #########
        # Setup progressbar
        #########
        self.progressbarLabel = Tkinter.Label(self, text="Progress: ")
        self.progressbar = ttk.Progressbar(self, orient="horizontal", length=200, mode="determinate")
        self.progressbar.pack(side='right', fill='x')
        self.progressbarLabel.pack(side='right', fill='x')


##########################################################################################
# mtsslDock class
##########################################################################################

#-----------------------------------------------------------------------------------------

class mtsslDockPlugin:

    #-------------------------------------------------------------------------------------

    def __init__(self, app):
        self.reset()
        self.parent = app.root
        self.messages = {
            'help': 'Save current file',
            'userevent': 'Saving file "foo"',
            'busy': 'Busy deleting all files from file system ...',
            'systemevent': 'File "foo" saved',
            'usererror': 'Invalid file name "foo/bar"',
            'systemerror': 'Failed to save file: file system full',
        }

        self.stateMessages = {
            0: '',
            1: 'Database is down',
            2: 'Waiting for reply from database',
        }
        #########
        # Setup main window
        #########

        self.showImportDialog()
        self.dialog = Pmw.Dialog(self.parent,
                                 buttons=('Dock',),
                                 title = 'mtsslDock Plugin',
                                 command = self.execute)

        self.center_window(self.dialog, 845, 725)
        # self.dialog.geometry('1024x768')
        # self.dialog.withdraw()
        Pmw.setbusycursorattributes(self.dialog.component('hull'))

        # self.master.geometry('1024x768')
        #self.mainWindowFrame = Tkinter.Frame(self.master, padx=0, pady=0)
        #self.mainWindowFrame.pack(fill=Tkinter.BOTH, side = Tkinter.TOP, expand = 5)
        #self.dockButtonFrame = Tkinter.Frame(master, padx = 0, pady = 0)
        #self.dockButtonFrame.pack(fill=Tkinter.BOTH, side = Tkinter.BOTTOM)

        #########
        # Setup menubar
        #########
        self.interior = self.dialog.interior()
        self.balloon = Pmw.Balloon(self.interior)
        menuBar = Pmw.MenuBar(self.interior, hull_relief=Tkinter.RAISED, hull_borderwidth=1, balloon=self.balloon)
        menuBar.pack(fill='x')
        menuBar.addmenu('File', 'Load and save data.')
        menuBar.addmenu('PyMOL', 'Interact with PyMOL')
        menuBar.addmenu('Settings', 'Adjust Settings')
        menuBar.addmenu('Help', 'Help for this program', side='right')
        menuBar.addmenuitem('File', 'command', 'Load constraints from file.', label='Load constraints...', command=lambda s=self: s.table.load())
        menuBar.addmenuitem('File', 'command', 'Save constraints to file.', label='Save constraints...', command=lambda s=self: s.table.save())
        menuBar.addmenuitem('File', 'separator')
        menuBar.addmenuitem('File', 'command', 'Load results from file.', label='Load result...', command=lambda s=self: s.loadResult())
        menuBar.addmenuitem('File', 'command', 'Save results to file.', label='Save result...', command=lambda s=self: s.saveResult())
        menuBar.addmenuitem('File', 'separator')
        menuBar.addmenuitem('File', 'command', 'Save table data to file.', label='Exit', command=lambda s=self: s.quit())
        #menuBar.addmenuitem('PyMOL', 'command', 'Synchronize with PyMOL session.', font=('StingerLight', 14), label='Synchronize', command = lambda s=self: s.execute('Synchronize with PyMOL'))
        menuBar.addmenuitem('PyMOL', 'command', 'Import Labels', label='Import Labels', command=lambda s=self: s.showImportDialog())
        menuBar.addmenuitem('Settings', 'command', 'Settings', label='Docking settings', command=lambda s=self: s.showSettingsDialog())
        menuBar.addmenuitem('Help', 'command', 'About', label='About', command=lambda s=self: s.showAboutDialog())

        #########
        # Setup messagebar
        #########
        self.statusbar = StatusBar(self.dialog.interior())
        self.statusbar.pack(side='bottom', fill='x', padx=10, pady=10)

        #########
        # Setup notebook
        #########
        self.notebook = Pmw.NoteBook(self.dialog.interior())
        self.notebook.pack(fill='both', expand=1, padx=10, pady=10)
        # Add the "Docking" page to the notebook.
        notebookpage = self.notebook.add('Docking')
        self.notebookpages['Docking'] = notebookpage
        self.notebook.tab('Docking').focus_set()

        #########
        # Table
        #########
        tableGroup = Pmw.Group(self.notebookpages['Docking'], tag_text='Constraints')
        tableGroup.pack(fill='both', side=Tkinter.RIGHT, anchor='nw', expand=2, padx=10, pady=5)

        # Setup table frame and table
        self.tframe = Tkinter.Frame(tableGroup.interior())
        self.tframe.pack(fill='both', anchor='n', expand=1, padx=10, pady=5)
        self.statusbar.messagebar.message('state', 'Waiting for input...')

    #-------------------------------------------------------------------------------------

    def loadResult(self):
        # get filename
        file_opt = options = {}
        filename = tkFileDialog.askopenfilename(**file_opt)
        if filename:
            file = open(filename, 'rb')
            results = pickle.load(file)
            # print results
            newTitle = "LoadedRun-%i" % self.dockingRunNumber
            notebookpage = self.notebook.add(newTitle)
            self.notebookpages[newTitle] = notebookpage
            self.notebook.tab(newTitle).focus_set()
            resultsFrame = ResultsFrame(self.notebookpages[newTitle], results)
            resultsFrame.pack(fill='both', expand=1, anchor='w', padx=10, pady=5)

    #-------------------------------------------------------------------------------------

    def saveResult(self):
        file_opt = options = {}
        filename = tkFileDialog.asksaveasfilename(**file_opt)
        if filename:
            file = open(filename, 'wb')
            selectedPage = self.notebook.getcurselection()
            if selectedPage != "Docking":
                children = {}
                # print self.notebookpages[selectedPage].children
                children = self.notebookpages[selectedPage].children
                # print children
                for v in children.values():
                    # print v.result
                    pickle.dump(v.result, file)
                    file.close()
    #-------------------------------------------------------------------------------------

    def addResultsPage(self):
        key = "dockingRun_%i" % self.dockingRunNumber
        newTitle = "Run-%i" % self.dockingRunNumber
        notebookpage = self.notebook.add(newTitle)
        self.notebookpages[newTitle] = notebookpage
        self.notebook.tab(newTitle).focus_set()
        resultsFrame = ResultsFrame(self.notebookpages[newTitle], self.dockingRuns[key])
        resultsFrame.pack(fill='both', expand=1, anchor='w', padx=10, pady=5)

    #-------------------------------------------------------------------------------------

    def center_window(self, window, w, h):
        # get screen width and height
        # gets screwed up with secondary screen - don't know how to fix that...
        ws = self.parent.winfo_screenwidth()
        hs = self.parent.winfo_screenheight()
        # calculate position x, y
        x = (ws / 2) - (w / 2)
        y = (hs / 2) - (h / 2)
        window.geometry('%dx%d+%d+%d' % (w, h, x, y))

    #-------------------------------------------------------------------------------------

    def reset(self):
        self.idle = True
        self.notebookpages = {}
        self.dockingRuns = {}
        self.objectPrefix = "mD"
        self.dockingRunNumber = 0
        self.myQueue = Queue.Queue()
        self.dockingProgress = 0
        self.symmetry = "None"
        #self.c2 = True
        self.scoreClashes = True
        self.scoreVdw = False
        self.farAwayPenalty = False
        self.numberOfLabelsForA = 0
        self.numberOfLabelsForB = 0
        self.labelPositionsProteinA = []
        self.labelPositionsProteinB = []
        self.labelNamesProteinA = []
        self.labelNamesProteinB = []
        self.numberOfPopulations = 10
        self.numberOfGenerations = 1000
        self.numberOfChromosomes = 400
        self.numberOfRigidBodyCycles = 50

    #-------------------------------------------------------------------------------------

    def showImportDialog(self):
        self.importDialog = ImportDialog(command=self.importLabels, buttons=('Synchronize with PyMOL', 'OK', ), title = 'Import Labels',)
        self.center_window(self.importDialog, 640, 480)

    #-------------------------------------------------------------------------------------

    def showSettingsDialog(self):
        self.settingsDialog = SettingsDialog(command=self.changeSettings, title='Settings')
        self.center_window(self.settingsDialog, 320, 240)

    #-------------------------------------------------------------------------------------

    def execute(self, result):
        if self.idle:
            self.statusbar.progressbar["maximum"] = self.numberOfGenerations + self.numberOfRigidBodyCycles
            self.statusbar.progressbar["value"] = 0
            self.dockingThread = threading.Thread(target=self.dock)
            self.dockingThread.setDaemon(1)
            self.dockingThread.start()
            self.parent.after(50, self.check_myQueue)
            self.parent.after(50, self.check_thread)

        # while self.idle == False:
        #	print "Working"
        # self.dock()

    #-------------------------------------------------------------------------------------

    def check_myQueue(self):
        while True:
            try:
                x = self.myQueue.get_nowait()
            except Queue.Empty:
                self.parent.after(25, self.check_myQueue)
                break
            else:  # continue from the try suite
                self.statusbar.progressbar["value"] += x

    #-------------------------------------------------------------------------------------

    def check_thread(self):
        if self.dockingThread.is_alive():
            self.parent.after(50, self.check_thread)
            self.statusbar.messagebar.message('state', 'Docking...')
            self.idle = False
        else:
            self.statusbar.messagebar.message('state', 'Idle.')
            self.addResultsPage()
            self.idle = True

    #-------------------------------------------------------------------------------------

    def changeSettings(self, result):
        self.numberOfPopulations = int(self.settingsDialog.numberOfPopulations.getvalue())
        self.numberOfGenerations = int(self.settingsDialog.numberOfGenerations.getvalue())
        self.numberOfChromosomes = int(self.settingsDialog.numberOfChromosomes.getvalue())
        self.numberOfRigidBodyCycles = int(self.settingsDialog.numberOfRigidBodyCycles.getvalue())
        self.scoreClashes = self.settingsDialog.scoreClashes
        #self.c2 = self.settingsDialog.c2
        self.symmetry = self.settingsDialog.symmetry
        print(self.symmetry)
        print(self.scoreClashes)
        self.settingsDialog.withdraw()

    #-------------------------------------------------------------------------------------

    def resetTable(self):
        emptyModel = TableModel(rows=1, columns=1)
        self.table = TableCanvas(self.tframe, model=emptyModel, showkeynamesinheader=True, cellwidth=200, rowheaderwidth=200)
        self.table.createTableFrame()
        self.table.tablecolheader.thefont = 'Arial 10'

    #-------------------------------------------------------------------------------------

    def showAboutDialog(self):
        Pmw.aboutversion('1.0')
        Pmw.aboutcontact(
            '  Gregor Hagelueken\n' +
            '  email: hagelueken@pc.uni-bonn.de\n' +
            '  www.pymolwiki.org/mtsslDock')
        self.about = Pmw.AboutDialog(self.parent, applicationname='mtsslDock')
        self.about.show()
        self.center_window(self.about, 350, 240)

    #-------------------------------------------------------------------------------------

    def importLabels(self, result):
        self.reset()
        self.resetTable()
        if result == 'OK':
            for labelSelector in self.importDialog.listOfLabelSelectors:
                if labelSelector.checkbuttons.getcurselection() == "A":
                    if self.numberOfLabelsForA == 0:
                        #self.table.addRow("A-%i" %(self.numberOfLabelsForA))
                        self.table.addRow("%i-%s" % (self.numberOfLabelsForA, labelSelector.labelText.get()))
                        self.table.model.deleteRows([0])
                        self.table.redrawTable()
                    else:
                        #self.table.addRow("A-%i" %(self.numberOfLabelsForA))
                        self.table.addRow("%i-%s" % (self.numberOfLabelsForA, labelSelector.labelText.get()))
                    coord = cmd.get_model(labelSelector.labelText.get(), 1).get_coord_list()
                    self.labelPositionsProteinA.append(coord[0])
                    self.labelNamesProteinA.append("%i-%s" % (self.numberOfLabelsForA, labelSelector.labelText.get()))
                    self.numberOfLabelsForA += 1
                elif labelSelector.checkbuttons.getcurselection() == "B":
                    if self.numberOfLabelsForB == 0:
                        #self.table.addColumn("B-%i" %(self.numberOfLabelsForB))
                        self.table.addColumn("%i-%s" % (self.numberOfLabelsForB, labelSelector.labelText.get()))
                        self.table.model.deleteColumn(0)
                        self.table.redrawTable()
                    else:
                        #self.table.addColumn("B-%i" %(self.numberOfLabelsForB))
                        self.table.addColumn("%i-%s" % (self.numberOfLabelsForB, labelSelector.labelText.get()))
                    coord = cmd.get_model(labelSelector.labelText.get(), 1).get_coord_list()
                    self.labelPositionsProteinB.append(coord[0])
                    self.labelNamesProteinB.append("%i-%s" % (self.numberOfLabelsForB, labelSelector.labelText.get()))
                    self.numberOfLabelsForB += 1
            self.proteinApymolString = "%s" % self.importDialog.proteinASelector.getcurselection()
            self.proteinBpymolString = "%s" % self.importDialog.proteinBSelector.getcurselection()
            self.importDialog.withdraw()
        elif result == 'Synchronize with PyMOL':
            for view in self.importDialog.listOfLabelSelectors:
                view.pack_forget()
            self.importDialog.items = []
            items = self.getPymolObjects("*_label")
            self.importDialog.listOfLabelSelectors = []
            for pymolObject in items:
                labelSelector = ABselect(self.importDialog.assignLabelsGroupScrolledFrame.interior())
                labelSelector.pack(fill='x', anchor='n', padx=4, pady=1)
                labelSelector.labelText.set(pymolObject)
                self.importDialog.listOfLabelSelectors.append(labelSelector)
            items = self.getPymolObjects("all")
            self.importDialog.proteinASelector.setitems(items)
            self.importDialog.proteinBSelector.setitems(items)

    #-------------------------------------------------------------------------------------

    def showAppModal(self):
        self.dialog.show()

    #-------------------------------------------------------------------------------------

    def dock(self):
        self.dockingRunNumber += 1
        key = "dockingRun_%i" % self.dockingRunNumber
        self.dockingRuns[key] = {}
        results = {}
        settings = {}
        settings['numberOfPopulations'] = self.numberOfPopulations
        settings['numberOfGenerations'] = self.numberOfGenerations
        settings['numberOfChromosomes'] = self.numberOfChromosomes
        settings['numberOfRigidBodyCycles'] = self.numberOfRigidBodyCycles
        settings['scoreClashes'] = self.scoreClashes
        settings['symmetry'] = self.symmetry

        # extract data from table and put into numpy array
        expDistances = numpy.zeros((self.numberOfLabelsForA, self.numberOfLabelsForB))
        expErrors = numpy.zeros((self.numberOfLabelsForA, self.numberOfLabelsForB))
        constraintNames = []
        rowCounter = 0
        for row in self.labelNamesProteinA:
            colCounter = 0
            for col in self.labelNamesProteinB:
                meanAndError = self.table.model.data[row][col].split(";")
                expDistances[rowCounter][colCounter] = float(meanAndError[0])
                constraintNames.append("%s-%s" % (row, col))
                if len(meanAndError) == 2:
                    expErrors[rowCounter][colCounter] = float(meanAndError[1])
                elif len(meanAndError) == 1:
                    expErrors[rowCounter][colCounter] = 1
                else:
                    print("Something wrong with input data!")
                    return
                colCounter += 1
            rowCounter += 1
        # print constraintNames
        trialDistances = quickDist2(numpy.array(self.labelPositionsProteinA), self.labelPositionsProteinB)
        print(trialDistances)
        print(numpy.reshape(trialDistances, (self.numberOfLabelsForA, self.numberOfLabelsForB)))
        expDistances = numpy.reshape(expDistances, (-1, 1))
        expErrors = numpy.reshape(expErrors, (-1, 1))

        testing = False
        if testing:
            # For testing: choose random distances
            expDistances = numpy.zeros((self.numberOfLabelsForA, self.numberOfLabelsForB))
            expErrors = numpy.zeros((self.numberOfLabelsForA, self.numberOfLabelsForB))
            expErrors.fill(1)
            #expDistances = numpy.reshape(expDistances, (-1, 1))
            trialDistancesTest = numpy.reshape(trialDistances, (self.numberOfLabelsForA, self.numberOfLabelsForB))
            print(trialDistancesTest)

            pickedDistances = 1
            alreadyPickedA = numpy.zeros(self.numberOfLabelsForA)
            alreadyPickedB = numpy.zeros(self.numberOfLabelsForB)
            row = 0
            unique = True
            errors = True

            while pickedDistances <= 10:
                posA = random.randrange(0, self.numberOfLabelsForA)
                posB = random.randrange(0, self.numberOfLabelsForB)
                if unique and self.symmetry == "None":
                    if alreadyPickedA[posA] == 0 and alreadyPickedB[posB] == 0:
                        expDistances[posA][posB] = trialDistancesTest[posA][posB]
                        if errors:
                            expErrors[posA][posB] = random.randrange(1, 6)
                        alreadyPickedA[posA] = 1
                        alreadyPickedB[posB] = 1
                        pickedDistances += 1
                elif unique and self.symmetry == "C2":
                    if alreadyPickedA[posA] == 0:
                        expDistances[posA][posA] = trialDistancesTest[posA][posA]
                        if errors:
                            expErrors[posA][posB] = random.randrange(1, 6)
                        alreadyPickedA[posA] = 1
                        alreadyPickedB[posB] = 1
                        pickedDistances += 1
                elif not unique and self.symmetry == "None":
                    if alreadyPickedA[posA] == 0 or alreadyPickedB[posB] == 0:
                        expDistances[posA][posB] = trialDistancesTest[posA][posB]
                        if errors:
                            expErrors[posA][posB] = random.randrange(1, 6)
                        alreadyPickedA[posA] = 1
                        alreadyPickedB[posB] = 1
                        pickedDistances += 1

            expDistances = numpy.reshape(expDistances, (-1, 1))
            expErrors = numpy.reshape(expErrors, (-1, 1))
            print(expDistances)
            print(expErrors)
            ####

        results['constraints'] = expDistances
        # print expDistances
        noDataIndices = numpy.where(expDistances == 0)[0]
        noDataIndices = numpy.where(expErrors == 0)[0]
        # print numpy.array(self.labelPositionsProteinA), self.labelPositionsProteinB

        print("trial distances: ")
        print(trialDistances)

        print("exp distances: ")
        print(expDistances)

        myView = cmd.get_view()
        # cmd.delete("*solution*")
        cmd.delete("tmp*")
        cmd.delete("test*")
        zeroChromosome = Chromosome("None")
        zeroChromosome.genes = numpy.array([0, 0, 0, 0, 0, 0])

        # setup Proteins
        # RdxR label cogs
        createPseudoatom(self.labelPositionsProteinA, "labelsProteinA", 1)
        #self.proteinApymolString = "%s" %self.proteinASelector.getcurselection()
        stored.atoms = []
        cmd.iterate_state(1, "%s & name ca" % self.proteinApymolString, 'stored.atoms.append((x,y,z))')
        allPositionsProteinA = numpy.array(stored.atoms)
        proteinA = Protein(self.labelPositionsProteinA, allPositionsProteinA, self.proteinApymolString)

        # Rdx label cogs
        #labelPositionsProteinB = numpy.array([[-12.460000038146973, -20.09000015258789, 24.219999313354492], [-26.829999923706055, -42.65999984741211, 19.31999969482422], [0.10999999940395355, -35.5099983215332, 20.15999984741211]])
        createPseudoatom(self.labelPositionsProteinB, "labelsProteinB", 1)

        #proteinBpymolString = "%s" %self.proteinBSelector.getcurselection()
        stored.atoms = []
        cmd.iterate_state(1, "%s & name ca" % self.proteinBpymolString, 'stored.atoms.append((x,y,z))')
        allPositionsProteinB = numpy.array(stored.atoms)
        proteinB = Protein(self.labelPositionsProteinB, allPositionsProteinB, self.proteinBpymolString)

        # move both to origin
        proteinA.moveToOrigin(proteinA.labelAtomsCog)
        proteinB.moveToOrigin(proteinB.labelAtomsCog)

        #proteinA.moveInPymol("tmpA_origin", zeroChromosome, 1)
        #proteinB.moveInPymol("tmpB_origin", zeroChromosome, 1)

        #######
        # evolve
        #######

        # setup populations
        print("setting up populations...")
        populations = []
        for i in range(0, self.numberOfPopulations):
            if self.symmetry == "C2":
                population = Population(self.numberOfChromosomes, "C2")
            elif self.symmetry == "None":
                population = Population(self.numberOfChromosomes, "None")
            population.name = "%i" % (i + 1)
            populations.append(population)

        # put them into an environment and evolve
        environment1 = Environment(populations, proteinA, proteinB, expDistances, expErrors, False, False, False)
        environment1.constraintNames = constraintNames
        environment1.applySelectionPressure()
        for population in environment1.populations:
            population.log += population.chromosomes[0].printChromosomeWithoutClashes()

        # first evolution
        print("evolving...")
        for i in range(0, self.numberOfGenerations):
            for population in environment1.populations:
                population.sacrifice(0.05)
                population.produceOffspring(rigidBody=False)
            environment1.applySelectionPressure()
            for population in environment1.populations:
                population.log += population.chromosomes[0].printChromosomeWithoutClashes()
            # for population in environment1.populations:
            #	print ("Pop. %s:" %(population.name)),
            #	population.chromosomes[0].printChromosome()
            self.myQueue.put(1)
            #self.statusbar.progressbar["value"] += i

        # prepare for rigid body
        tmpBestSolutions = []
        for population in environment1.populations:
            population.sacrifice(0.95)
            population.resetFitness()
            tmpBestSolutions.append(population.chromosomes[0])

        # print "Best solutions before rigid body: "
        # for solution in tmpBestSolutions:
        #	print solution.printChromosomeWithoutClashes()

        # add best chromosome from each population to each other population
        # for i in range (0, len(tmpBestSolutions)):
        #	bestOfOthers = tmpBestSolutions[:i] + tmpBestSolutions[(i + 1):]
            # print bestOfOthers
        #	environment1.populations[i].chromosomes += bestOfOthers
            # print environment1.populations[i].chromosomes

        # for population in environment1.populations:
        #	print "Best solutions: "
        #	for solution in population.chromosomes:
        #		print solution.printChromosomeWithoutClashes()

        environment1.scoreClashes = self.scoreClashes
        environment1.farAwayPenalty = self.farAwayPenalty
        environment1.applySelectionPressure()
        for population in environment1.populations:
            population.log += population.chromosomes[0].printChromosomeWithoutClashes()
        environment1.scoreVdw = self.scoreVdw

        print("rigid body starts")
        for i in range(0, self.numberOfRigidBodyCycles):
            for population in environment1.populations:
                population.sacrifice(0.95)
                population.produceOffspring(rigidBody=True)
            environment1.applySelectionPressure()
            for population in environment1.populations:
                population.log += population.chromosomes[0].printChromosomeWithClashes()
            # print "Rigid Body cycle: %i" %i
            # for population in environment1.populations:
            #	print ("Pop. %s:" %(population.name)),
            #	population.chromosomes[0].printChromosome()
            self.myQueue.put(1)

        # create solutions
        print("")
        print("Solutions:")
        nonClashingSolution = 1
        clashingSolution = 1
        for population in environment1.populations:
            createPseudoatom(self.labelPositionsProteinB, "tmpSolution-labels", 1)
            tmpProtein = Protein(self.labelPositionsProteinB, self.labelPositionsProteinB, "tmpSolution-labels")
            solution = population.chromosomes[0]
            # print solution.printChromosomeWithClashes()
            if solution.clashes <= 5:
                nameOfSolution = "%s-%i_sol-%i" % (self.objectPrefix, self.dockingRunNumber, nonClashingSolution)
                solution.name = nameOfSolution

                proteinB.moveInPymol(nameOfSolution, solution, 1)
                tmpProtein.moveInPymol("%s-labels" % nameOfSolution, solution, 1)
                cmd.translate(list(proteinA.labelAtomsCog.reshape(-1,)), nameOfSolution, 1, 0, None)
                cmd.translate(list(proteinA.labelAtomsCog.reshape(-1,)), "%s-labels" % nameOfSolution, 1, 0, None)
                nonClashingSolution += 1
                # for testing
                rms_cur = cmd.rms_cur(nameOfSolution, "target")
                # print rms_cur
                print("%s, %1.2f" % (nameOfSolution, rms_cur))
                ###
            elif solution.clashes > 5:
                nameOfSolution = "%s-%i_clash-%i" % (self.objectPrefix, self.dockingRunNumber, clashingSolution)
                solution.name = nameOfSolution
                proteinB.moveInPymol(nameOfSolution, solution, 1)
                tmpProtein.moveInPymol("%s-labels" % nameOfSolution, solution, 1)
                cmd.translate(list(proteinA.labelAtomsCog.reshape(-1,)), nameOfSolution, 1, 0, None)
                cmd.translate(list(proteinA.labelAtomsCog.reshape(-1,)), "%s-labels" % nameOfSolution, 1, 0, None)
                clashingSolution += 1
        cmd.group("%s-%i" % (self.objectPrefix, self.dockingRunNumber), "%s-%i*" % (self.objectPrefix, self.dockingRunNumber))
        cmd.set_view(myView)

        print("Done!")

        # add results to dictionary
        self.dockingRuns[key]['settings'] = settings
        self.dockingRuns[key]['environment'] = environment1

    #-------------------------------------------------------------------------------------

    def getPymolObjects(self, selection):
        return cmd.get_object_list('(%s)' % selection)

    #-------------------------------------------------------------------------------------

    def quit(self):
        self.dialog.withdraw()


##########################################################################################
# Some helper functions
##########################################################################################

#-----------------------------------------------------------------------------------------

def createPseudoatom(coordinates, objectName, state):
    for coordinate in coordinates:
        x = float(coordinate[0])
        y = float(coordinate[1])
        z = float(coordinate[2])
        posString = "[%3.2f,%3.2f,%3.2f]" % (x, y, z)
        cmd.pseudoatom(pos=posString, object=objectName, state=state)

#-----------------------------------------------------------------------------------------


def quickClash(atomsProteinA, atomsProteinB, cutoff):
    dist = scipy.spatial.distance.cdist(atomsProteinA, atomsProteinB)
    clashes = len(numpy.nonzero(dist < cutoff)[0])
    return clashes

#-----------------------------------------------------------------------------------------


def quickDist2(atoms1, atoms2):
    # if there is only one atom it has to be duplicated for quick_dist2 to work
    duplicatedA = False
    duplicatedB = False
    if len(atoms1) == 1:
        duplicatedA = True
        atoms1 = numpy.tile(atoms1, (2, 1))
        # print atoms1
    if len(atoms2) == 1:
        duplicatedB = True
        atoms2 = numpy.tile(atoms2, (2, 1))
        # print atoms2
    dist = scipy.spatial.distance.cdist(atoms1, atoms2)
    # print dist
    # remove the duplication depending on which selection contained the single atom
    if duplicatedA and not duplicatedB:
        dist = numpy.reshape(dist[0, :], (-1, 1))

    elif duplicatedB and not duplicatedA:
        dist = numpy.reshape(dist[:, 0], (-1, 1))

    elif duplicatedA and duplicatedB:
        dist = numpy.reshape(dist[:1, 0], (-1, 1))
    else:
        dist = numpy.reshape(dist, (-1, 1))
    return dist

##########################################################################################
# classes for genetic algorithm
##########################################################################################

#-----------------------------------------------------------------------------------------


class Protein:

    #-------------------------------------------------------------------------------------

    def __init__(self, labelAtoms, allAtoms, pymolString):
        self.originalAllAtoms = numpy.copy(allAtoms)
        self.originalLabelAtoms = numpy.copy(labelAtoms)
        self.labelAtoms = numpy.copy(labelAtoms)
        self.allAtoms = numpy.copy(allAtoms)
        self.labelAtomsCog = self.getCog(self.labelAtoms)
        self.pymolString = pymolString

    #-------------------------------------------------------------------------------------

    def getCog(self, atoms):
        return numpy.average(atoms, axis=0)

    #-------------------------------------------------------------------------------------

    def moveInPymol(self, solutionPymolString, chromosome, state):
        # recreate protein B at solution coordinates
        cmd.create(solutionPymolString, self.pymolString, 1, state)
        # translate to origin with proteinBcog. IMPORTANT: set camera to "0" so that the translation is not done along the camera coordinate system!
        cmd.translate(list(-1 * self.labelAtomsCog.reshape(-1,)), solutionPymolString, state, 0, None)

        # rotate and translate according to solution
        translation = chromosome.genes[3:6]
        # print list(translation.reshape(-1,))
        rotX = chromosome.genes[0]
        # IMPORTANT: set camera to "0" so that the translation is not done along the camera coordinate system! Also set rotation origin to 0,0,0!
        cmd.rotate([1, 0, 0], rotX, solutionPymolString, state, 0, None, [0, 0, 0])
        rotY = chromosome.genes[1]
        cmd.rotate([0, 1, 0], rotY, solutionPymolString, state, 0, None, [0, 0, 0])
        rotZ = chromosome.genes[2]
        cmd.rotate([0, 0, 1], rotZ, solutionPymolString, state, 0, None, [0, 0, 0])
        cmd.translate(list(translation.reshape(-1,)), solutionPymolString, state, 0, None)
        # print chromosome.c2RotationAngles

        if chromosome.symmetry == "C2":
            # translate again for C2
            cmd.translate(list(translation.reshape(-1,)), solutionPymolString, state, 0, None)
            # align along X axis
            angleX = chromosome.c2RotationAngles[0]
            # IMPORTANT: set camera to "0" so that the translation is not done along the camera coordinate system! Also set rotation origin to 0,0,0!
            cmd.rotate(chromosome.c2CrossProduct.tolist(), angleX, solutionPymolString, state, 0, None, [0, 0, 0])

            # match local and global z Axis
            cmd.rotate([1, 0, 0], chromosome.zAngle, solutionPymolString, state, 0, None, [0, 0, 0])
            #rotY = chromosome.c2RotationAngles[1]
            #cmd.rotate([0, 1, 0], -rotY, solutionPymolString,state,0,None,[0,0,0])
            #rotZ = chromosome.c2RotationAngles[2]
            #cmd.rotate([0, 0, 1], -rotZ, solutionPymolString,state,0,None,[0,0,0])

    #-------------------------------------------------------------------------------------

    def moveToOrigin(self, cog):
        self.allAtoms = self.allAtoms - cog
        self.labelAtoms = self.labelAtoms - cog

    #-------------------------------------------------------------------------------------

    def calculatePositionFromChromosome(self, chromosome, onlyLabelAtoms):

        if onlyLabelAtoms:
            newPosition = numpy.copy(self.labelAtoms)
        else:
            newPosition = numpy.copy(self.allAtoms)

        orientationVector = numpy.array([[1, 0, 0]])
        zVector = numpy.array([[0, 0, 1]])
        #createPseudoatom (numpy.array([[0,0,0]]), "myOrigin", 1)
        #createPseudoatom (numpy.array([[1,0,0]]), "x", 1)
        #createPseudoatom (numpy.array([[0,1,0]]), "y", 1)
        #createPseudoatom (numpy.array([[0,0,1]]), "z", 1)

        #createPseudoatom (orientationVector, "orientationVector", 1)
        # print "This should be 0, 90, 90: ", self.getAngle(orientationVector, numpy.array([1,0,0])), self.getAngle(orientationVector, numpy.array([0,1,0])), self.getAngle(orientationVector, numpy.array([0,0,1]))
        # print chromosome.printChromosomeWithClashes()
        # rotate around x
        rotationMatrix = self.setupRotationMatrix(chromosome.genes[0], numpy.array([1, 0, 0]))
        newPosition = self.rotatePoints(newPosition, rotationMatrix)
        orientationVector = self.rotatePoints(orientationVector, rotationMatrix)
        zVector = self.rotatePoints(zVector, rotationMatrix)

        # rotate around y
        rotationMatrix = self.setupRotationMatrix(chromosome.genes[1], numpy.array([0, 1, 0]))
        newPosition = self.rotatePoints(newPosition, rotationMatrix)
        orientationVector = self.rotatePoints(orientationVector, rotationMatrix)
        zVector = self.rotatePoints(zVector, rotationMatrix)

        # rotate around z
        rotationMatrix = self.setupRotationMatrix(chromosome.genes[2], numpy.array([0, 0, 1]))
        newPosition = self.rotatePoints(newPosition, rotationMatrix)
        orientationVector = self.rotatePoints(orientationVector, rotationMatrix)
        zVector = self.rotatePoints(zVector, rotationMatrix)

        #createPseudoatom (orientationVector, "orientationVector", 2)
        # print "Now rotated: ", self.getAngle(orientationVector, numpy.array([1,0,0])), self.getAngle(orientationVector, numpy.array([0,1,0])), self.getAngle(orientationVector, numpy.array([0,0,1]))

        # translate along x, y, z
        newPosition = newPosition + numpy.array([chromosome.genes[3], 0, 0])
        #orientationVector = orientationVector + numpy.array([chromosome.genes[3],0,0])
        newPosition = newPosition + numpy.array([0, chromosome.genes[4], 0])
        #orientationVector = orientationVector + numpy.array([0,chromosome.genes[4],0])
        newPosition = newPosition + numpy.array([0, 0, chromosome.genes[5]])
        #orientationVector = orientationVector + numpy.array([0,0,chromosome.genes[5]])

        # return position now, if no symmetry
        if chromosome.symmetry == "None":
            if onlyLabelAtoms:
                self.labelAtoms = newPosition
            else:
                self.allAtoms = newPosition
            return

        # create C2 symmetric dimer
        elif chromosome.symmetry == "C2":
            # rotate 180 around c2 axis
            rotationMatrix = self.setupRotationMatrix(180, numpy.array([1, 0, 0]))
            newPosition = self.rotatePoints(newPosition, rotationMatrix)
            orientationVector = self.rotatePoints(orientationVector, rotationMatrix)
            zVector = self.rotatePoints(zVector, rotationMatrix)

            #createPseudoatom (orientationVector, "orientationVector", 3)
            # print "Rotated by 180: ", self.getAngle(orientationVector, numpy.array([1,0,0])), self.getAngle(orientationVector, numpy.array([0,1,0])), self.getAngle(orientationVector, numpy.array([0,0,1]))

            # measure angles between orientationVector and axes of global coordinate system
            angleX = self.getAngle(orientationVector, numpy.array([1, 0, 0]))
            angleY = self.getAngle(orientationVector, numpy.array([0, 1, 0]))
            angleZ = self.getAngle(orientationVector, numpy.array([0, 0, 1]))
            chromosome.c2RotationAngles = numpy.array([angleX, angleY, angleZ])
            # print "Angles", angleX, angleY, angleZ

            # rotate points back by 180 around c2 axis, but leave orientationVector untouched
            rotationMatrix = self.setupRotationMatrix(-180, numpy.array([1, 0, 0]))
            newPosition = self.rotatePoints(newPosition, rotationMatrix)
            #orientationVector = self.rotatePoints(orientationVector, rotationMatrix)
            # print "Rotated back by -180: ", self.getAngle(orientationVector, numpy.array([1,0,0])), self.getAngle(orientationVector, numpy.array([0,1,0])), self.getAngle(orientationVector, numpy.array([0,0,1]))

            # translate to dimer position by repeating translation
            newPosition = newPosition + numpy.array([chromosome.genes[3], 0, 0])
            #orientationVector = orientationVector + numpy.array([chromosome.genes[3],0,0])
            newPosition = newPosition + numpy.array([0, chromosome.genes[4], 0])
            #orientationVector = orientationVector + numpy.array([0,chromosome.genes[4],0])
            newPosition = newPosition + numpy.array([0, 0, chromosome.genes[5]])
            #orientationVector = orientationVector + numpy.array([0,0,chromosome.genes[5]])

            # rotate
            # rotate so that orientationVector is along x
            rotationAxis = numpy.cross(orientationVector[0], numpy.array([1, 0, 0]))
            chromosome.c2CrossProduct = rotationAxis
            rotationMatrix = self.setupRotationMatrix(angleX, rotationAxis)
            newPosition = self.rotatePoints(newPosition, rotationMatrix)
            zVector = self.rotatePoints(zVector, rotationMatrix)
            orientationVector = self.rotatePoints(orientationVector, rotationMatrix)
            #createPseudoatom (orientationVector, "orientationVector", 4)
            # print "orientationVector angle to X: ", self.getAngle(orientationVector, numpy.array([1,0,0]))

            # print "zVector: ", zVector
            #createPseudoatom (zVector, "zVector", 1)
            # print "zVectorAngle to X: ", self.getAngle(zVector, numpy.array([1,0,0]))
            # print "zVectorAngle to Y: ", self.getAngle(zVector, numpy.array([0,1,0]))
            # print "zVectorAngle to Z: ", self.getAngle(zVector, numpy.array([0,0,1]))
            # print "crossproduct zVector and z:", numpy.cross(zVector[0], numpy.array([0,0,1]))
            zVectorAngleToZ = self.getAngle(zVector, numpy.array([0, 0, 1]))
            # print "zVector angle to z: ", zVectorAngleToZ

            # rotate around x to match original position
            if numpy.cross(zVector[0], numpy.array([0, 0, 1]))[0] <= 0:
                angle = -zVectorAngleToZ
            else:
                angle = zVectorAngleToZ
            chromosome.zAngle = angle
            rotationMatrix = self.setupRotationMatrix(angle, numpy.array([1, 0, 0]))
            newPosition = self.rotatePoints(newPosition, rotationMatrix)
            zVector = self.rotatePoints(zVector, rotationMatrix)

            # print "zVector angle should be 0 again:", self.getAngle(zVector, numpy.array([0,0,1]))

            # print "Rotated back along x: ", self.getAngle(orientationVector, numpy.array([1,0,0])), self.getAngle(orientationVector, numpy.array([0,1,0])), self.getAngle(orientationVector, numpy.array([0,0,1]))
            # print "Rotated back to 1,0,0 xaxis: ", self.getAngle(orientationVector, numpy.array([1,0,0])), self.getAngle(orientationVector, numpy.array([0,1,0])), self.getAngle(orientationVector, numpy.array([0,0,1]))

            # rotate around y
            #rotationMatrix = self.setupRotationMatrix(angleY, numpy.cross(orientationVector[0], numpy.array([0,1,0])))
            #newPosition = self.rotatePoints(newPosition, rotationMatrix)
            #tmporientationVector = self.rotatePoints(orientationVector, rotationMatrix)

            # rotate around z
            #rotationMatrix = self.setupRotationMatrix(angleX, numpy.cross(orientationVector[0], numpy.array([0,0,1])))
            #newPosition = self.rotatePoints(newPosition, rotationMatrix)
            #tmporientationVector = self.rotatePoints(orientationVector, rotationMatrix)

            if onlyLabelAtoms:
                self.labelAtoms = newPosition
            else:
                self.allAtoms = newPosition
            return

    #-------------------------------------------------------------------------------------

    def getAngle(self, a, b):
        return numpy.rad2deg(numpy.arccos(numpy.dot(a / numpy.linalg.norm(a), b / numpy.linalg.norm(b))))
    #-------------------------------------------------------------------------------------

    def rotatePoints(self, points, rotationMatrix):
        rotatedPoints = []
        for point in points:
            # add 1 for multiplication with 4x4 matrix
            point = numpy.append(point, 1)
            rotatedPoint = numpy.dot(rotationMatrix, point)
            # remove 1 again
            rotatedPoint = numpy.delete(rotatedPoint, 3)
            rotatedPoints.append(rotatedPoint)
        return numpy.array(rotatedPoints)

    #-------------------------------------------------------------------------------------

    def setupRotationMatrix(self, angle, axisPoint):
        # print angle
        u = axisPoint[0]
        v = axisPoint[1]
        w = axisPoint[2]
        L = (u * u + v * v + w * w)
        angle = angle * numpy.pi / 180.0
        u2 = u * u
        v2 = v * v
        w2 = w * w
        rotationMatrix = numpy.zeros((4, 4))
        rotationMatrix[0][0] = (u2 + (v2 + w2) * numpy.cos(angle)) / L
        rotationMatrix[0][1] = (u * v * (1 - numpy.cos(angle)) - w * numpy.sqrt(L) * numpy.sin(angle)) / L
        rotationMatrix[0][2] = (u * w * (1 - numpy.cos(angle)) + v * numpy.sqrt(L) * numpy.sin(angle)) / L
        rotationMatrix[0][3] = 0.0

        rotationMatrix[1][0] = (u * v * (1 - numpy.cos(angle)) + w * numpy.sqrt(L) * numpy.sin(angle)) / L
        rotationMatrix[1][1] = (v2 + (u2 + w2) * numpy.cos(angle)) / L
        rotationMatrix[1][2] = (v * w * (1 - numpy.cos(angle)) - u * numpy.sqrt(L) * numpy.sin(angle)) / L
        rotationMatrix[1][3] = 0.0

        rotationMatrix[2][0] = (u * w * (1 - numpy.cos(angle)) - v * numpy.sqrt(L) * numpy.sin(angle)) / L
        rotationMatrix[2][1] = (v * w * (1 - numpy.cos(angle)) + u * numpy.sqrt(L) * numpy.sin(angle)) / L
        rotationMatrix[2][2] = (w2 + (u2 + v2) * numpy.cos(angle)) / L
        rotationMatrix[2][3] = 0.0

        rotationMatrix[3][0] = 0.0
        rotationMatrix[3][1] = 0.0
        rotationMatrix[3][2] = 0.0
        rotationMatrix[3][3] = 1.0
        return rotationMatrix

#-----------------------------------------------------------------------------------------


class Environment:

    #-------------------------------------------------------------------------------------

    def __init__(self, populations, proteinA, proteinB, expDistances, expErrors, scoreClashes, scoreVdw, farAwayPenalty):
        self.populations = populations
        self.proteinA = proteinA
        self.proteinB = proteinB
        self.expDistances = expDistances
        self.expErrors = expErrors
        self.scoreClashes = scoreClashes
        self.scoreVdw = scoreVdw
        self.farAwayPenalty = farAwayPenalty
        self.constraintNames = []

    #-------------------------------------------------------------------------------------

    def getFitness(self, trialDistances, trialAtomsProteinB):
        rmsd = 0
        chi2 = 0
        clashes = 0
        vdwContacts = 0
        rmsd = numpy.sqrt(numpy.mean(numpy.square(self.expDistances - trialDistances)))
        chi2 = numpy.sum(numpy.square((self.expDistances - numpy.sqrt(numpy.square(trialDistances))) / self.expErrors))
        #chi2red = chi2/(len(numpy.where(self.expDistances == 0)[0])-6)

        if self.scoreClashes:
            clashes = quickClash(self.proteinA.allAtoms, trialAtomsProteinB, 3.5)

        fitnessResults = [rmsd, chi2, clashes]
        return fitnessResults

    #-------------------------------------------------------------------------------------

    def applySelectionPressure(self):
        # Determine fitness of childs
        for population in self.populations:
            counter = 0
            for chromosome in population.chromosomes:
                if chromosome.fitness == 0:
                    # print "Calculating fitness"
                    tmpProteinB = copy.deepcopy(self.proteinB)
                    tmpProteinB.calculatePositionFromChromosome(chromosome, onlyLabelAtoms=True)

                    # remove distances that where not measured form trialDistances
                    noDataIndices = numpy.where(self.expDistances == 0)[0]
                    trialDistances = quickDist2(self.proteinA.labelAtoms, tmpProteinB.labelAtoms)
                    # print trialDistances
                    for index in noDataIndices:
                        trialDistances[index] = 0
                    # print trialDistances

                    if self.scoreClashes or self.scoreVdw or self.farAwayPenalty:
                        tmpProteinB.calculatePositionFromChromosome(chromosome, onlyLabelAtoms=False)
                        fitness = self.getFitness(trialDistances, tmpProteinB.allAtoms)
                    else:
                        fitness = self.getFitness(trialDistances, tmpProteinB.labelAtoms)
                    #chromosome.fitness = fitness[1]
                    chromosome.rmsd = fitness[0]
                    chromosome.chi2 = fitness[1]
                    chromosome.clashes = fitness[2]
                    # add clash penalty
                    chromosome.fitness = chromosome.chi2 * 1 + chromosome.clashes
                    if self.scoreClashes:
                        chromosome.log += chromosome.printChromosomeWithClashes()
                    else:
                        chromosome.log += chromosome.printChromosomeWithoutClashes()
                    chromosome.trialDistances = trialDistances
                counter += 1
            population.rank()

#-----------------------------------------------------------------------------------------


class Population:

    #-------------------------------------------------------------------------------------

    def __init__(self, size, symmetry):
        self.size = size
        self.symmetry = symmetry
        self.chromosomes = self.randomSpawn(symmetry)
        self.sacrificed = 0
        self.rank()
        self.name = ""
        self.log = ""

    #-------------------------------------------------------------------------------------

    def resetFitness(self):
        for chromosome in self.chromosomes:
            chromosome.fitness = 0

    #-------------------------------------------------------------------------------------

    def randomSpawn(self, symmetry):
        chromosomes = []
        for i in range(0, self.size):
            chromosome = Chromosome(symmetry)
            chromosomes.append(chromosome)
            # chromosome.printChromosome()
        return chromosomes

    #-------------------------------------------------------------------------------------

    def rank(self):
        self.chromosomes = sorted(self.chromosomes, key=attrgetter('fitness'))

    #-------------------------------------------------------------------------------------

    def sacrifice(self, mutationFrequency):
        numberToSacrifice = int(len(self.chromosomes) * mutationFrequency)
        self.chromosomes = self.chromosomes[0:len(self.chromosomes) - numberToSacrifice]
        self.sacrificed = numberToSacrifice

    #-------------------------------------------------------------------------------------

    def printChromosomes(self):
        for chromosome in self.chromosomes:
            chromosome.printChromosomeWithClashes()

    #-------------------------------------------------------------------------------------

    def getGeneticOperator(self):
        if random.choice([True, False]):
            if random.choice([True, False]):
                return "smallCreepMutation"
            else:
                return "randomMutation"
        else:
            if random.choice([True, False]):
                return "singlePointCrossover"
            else:
                return "exchangeCrossover"

    #-------------------------------------------------------------------------------------

    def smallCreepMutation(self, parent):
        childs = []
        #child = Chromosome()
        child = copy.deepcopy(parent)

        position = random.randrange(6)
        # print position
        oldValue = parent.genes[position]
        newValue = 0
        # rotation
        if position < 3:
            creep = 0
            while creep == 0:
                creep = random.uniform(-5.0, 5.0)
                # print creep
            newValue = oldValue + creep
            if newValue > 360:
                newValue -= 360
            elif newValue < 0:
                newValue += 360
        # translation
        else:
            creep = 0
            while creep == 0:
                creep = random.uniform(-2.0, 2.0)
                # print creep
            newValue = oldValue + creep
        child.genes[position] = newValue
        #child.log = ''
        child.fitness = 0
        childs.append(child)
        return childs

    #-------------------------------------------------------------------------------------

    def randomMutation(self, parent):
        childs = []
        #child = Chromosome()
        child = copy.deepcopy(parent)
        position = random.randrange(6)
        # print child.printChromosomeWithClashes()
        # print position, child.c2Axis, parent.c2Axis
        oldValue = parent.genes[position]
        newValue = 0
        # rotation
        if position < 3:
            newValue = random.randrange(360)
        # translation
        else:
            newValue = random.randrange(-50, 50)
        child.genes[position] = newValue
        #child.log = ''
        child.fitness = 0
        childs.append(child)
        return childs

    #-------------------------------------------------------------------------------------

    def singlePointCrossover(self, parent1, parent2):
        childs = []
        position = random.randrange(6)
        child1 = copy.deepcopy(parent1)
        child2 = copy.deepcopy(parent2)
        for i in range(0, 5):
            if i <= position:
                child1.genes[i] = parent1.genes[i]
                child2.genes[i] = parent2.genes[i]
            else:
                child1.genes[i] = parent2.genes[i]
                child2.genes[i] = parent1.genes[i]
        child1.fitness = 0
        #child1.log = ''
        child2.fitness = 0
        #child2.log = ''
        childs.append(child1)
        childs.append(child2)
        return childs

    #-------------------------------------------------------------------------------------

    def exchangeCrossover(self, parent1, parent2):
        childs = []
        position = random.randrange(6)
        child1 = copy.deepcopy(parent1)
        child2 = copy.deepcopy(parent2)
        tmpValue = child1.genes[position]
        child1.genes[position] = child2.genes[position]
        child2.genes[position] = tmpValue
        child1.fitness = 0
        #child1.log = ''
        child2.fitness = 0
        #child2.log = ''
        childs.append(child1)
        childs.append(child2)
        return childs

    #-------------------------------------------------------------------------------------

    def produceOffspring(self, rigidBody):
        j = 0
        while j < self.sacrificed:
            geneticOperator = ""

            # make sure smallCreepMutation is selected for rigid body refinements
            if not rigidBody:
                geneticOperator = self.getGeneticOperator()
            else:
                geneticOperator = "smallCreepMutation"
                # while not geneticOperator == "smallCreepMutation":

            childs = []
            if geneticOperator == "smallCreepMutation":
                if rigidBody:
                    # take fittest chromosome as parent for rigid body
                    # print "Choosing best one"
                    parent = self.chromosomes[0]  # random.randrange(3)]
                else:
                    parent = self.chromosomes[random.randrange(len(self.chromosomes))]
                childs = self.smallCreepMutation(copy.deepcopy(parent))
                j += 1

            elif geneticOperator == "randomMutation":
                parent = self.chromosomes[random.randrange(len(self.chromosomes))]
                childs = self.randomMutation(copy.deepcopy(parent))
                j += 1

            elif geneticOperator == "singlePointCrossover" and (self.sacrificed - j <= 2):
                parent1 = self.chromosomes[random.randrange(len(self.chromosomes))]
                parent2 = self.chromosomes[random.randrange(len(self.chromosomes))]
                childs = self.singlePointCrossover(copy.deepcopy(parent1), copy.deepcopy(parent2))
                j += 2

            elif geneticOperator == "exchangeCrossover" and (self.sacrificed - j <= 2):
                parent1 = self.chromosomes[random.randrange(len(self.chromosomes))]
                parent2 = self.chromosomes[random.randrange(len(self.chromosomes))]
                childs = self.exchangeCrossover(copy.deepcopy(parent1), copy.deepcopy(parent2))
                j += 2

            # force translation along twofold (global x-Axis) to 0 for C2 symmetry to avoid screw axis
            if self.symmetry == "C2":
                for child in childs:
                    child.genes[3] = 0

            self.chromosomes += childs

#-----------------------------------------------------------------------------------------


class Chromosome:

    #-------------------------------------------------------------------------------------

    def __init__(self, symmetry):
        self.maxTranslation = 100
        self.symmetry = symmetry
        #self.c2Axis = 0
        if self.symmetry == "None":
            self.genes = self.generateRandomGenes()
        elif self.symmetry == "C2":
            self.genes = self.generateC2Genes()

        self.fitness = 0
        self.chi2 = 0
        self.rmsd = 0
        self.clashes = 0
        self.chi2red = 0
        self.log = ""
        self.trialDistances = numpy.zeros((3, 3))
        self.name = ""
        self.c2RotationAngles = numpy.zeros(3)
        self.c2CrossProduct = numpy.zeros(3)
        self.zAngle = 0

    #-------------------------------------------------------------------------------------

    def generateRandomGenes(self):
        # translation between -maxTranslation and maxTranslation
        translation = self.maxTranslation * numpy.random.uniform(-1, 1, size=3)
        rotation = 360 * numpy.random.uniform(0, 1, size=3)
        return numpy.concatenate((rotation, translation), axis=1)

    #-------------------------------------------------------------------------------------

    def generateC2Genes(self):
        # translation between -maxTranslation and maxTranslation
        translation = self.maxTranslation * numpy.random.uniform(-1, 1, size=3)
        rotation = 360 * numpy.random.uniform(0, 1, size=3)
        #self.c2Axis = 0
        translation[0] = 0
        return numpy.concatenate((rotation, translation), axis=1)

    #-------------------------------------------------------------------------------------
    def printChromosomeWithClashes(self):
        #string = "rotX: %1.2f\t rotY: %1.2f\t rotZ: %1.2f\t tX: %1.2f\t tY: %1.2f\t tZ: %1.2f\t rmsd:%1.2f\t chi2:%1.2f\t clashes:%i\n" %(self.genes[0],self.genes[1],self.genes[2],self.genes[3],self.genes[4],self.genes[5], self.rmsd, self.chi2, self.clashes)
        string = "%1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %i\n" % (self.genes[0], self.genes[1], self.genes[2], self.genes[3], self.genes[4], self.genes[5], self.rmsd, self.chi2, self.clashes)
        # print string
        return string

    def printChromosomeWithoutClashes(self):
        #string = "rotX: %1.2f\t rotY: %1.2f\t rotZ: %1.2f\t tX: %1.2f\t tY: %1.2f\t tZ: %1.2f\t rmsd:%1.2f\t chi2:%1.2f\t clashes:%i\n" %(self.genes[0],self.genes[1],self.genes[2],self.genes[3],self.genes[4],self.genes[5], self.rmsd, self.chi2, self.clashes)
        string = "%1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f -\n" % (self.genes[0], self.genes[1], self.genes[2], self.genes[3], self.genes[4], self.genes[5], self.rmsd, self.chi2)
        # print string
        return string


class ConsoleText(Tkinter.Text):

    '''A Tkinter Text widget that provides a scrolling display of console
    stderr and stdout.'''

    class IORedirector(object):

        '''A general class for redirecting I/O to this Text widget.'''

        def __init__(self, text_area):
            self.text_area = text_area

    class StdoutRedirector(IORedirector):

        '''A class for redirecting stdout to this Text widget.'''

        def write(self, str):
            self.text_area.write(str, False)

    class StderrRedirector(IORedirector):

        '''A class for redirecting stderr to this Text widget.'''

        def write(self, str):
            self.text_area.write(str, True)

    def __init__(self, master=None, cnf={}, **kw):
        '''See the __init__ for Tkinter.Text for most of this stuff.'''

        tk.Text.__init__(self, master, cnf, **kw)

        self.started = False
        self.write_lock = threading.Lock()

        self.tag_configure('STDOUT', background='white', foreground='black')
        self.tag_configure('STDERR', background='white', foreground='red')

        self.config(state=tk.DISABLED)

    def start(self):

        if self.started:
            return

        self.started = True

        self.original_stdout = sys.stdout
        self.original_stderr = sys.stderr

        stdout_redirector = ConsoleText.StdoutRedirector(self)
        stderr_redirector = ConsoleText.StderrRedirector(self)

        sys.stdout = stdout_redirector
        sys.stderr = stderr_redirector

    def stop(self):

        if not self.started:
            return

        self.started = False

        sys.stdout = self.original_stdout
        sys.stderr = self.original_stderr

    def write(self, val, is_stderr=False):

        # Fun Fact:	The way Tkinter Text objects work is that if they're disabled,
        # you can't write into them AT ALL (via the GUI or programatically).	 Since we want them
        # disabled for the user, we have to set them to NORMAL (a.k.a. ENABLED), write to them,
        # then set their state back to DISABLED.

        self.write_lock.acquire()
        self.config(state=tk.NORMAL)

        self.insert('end', val, 'STDERR' if is_stderr else 'STDOUT')
        self.see('end')

        self.config(state=tk.DISABLED)
        self.write_lock.release()
