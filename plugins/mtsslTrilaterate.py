'''
Title   : mtsslTrilaterate
Version : 1.0
Author  : Dinar Abdullin
Date    : Feb 2013
Email   : abdullin@pc.uni-bonn.de
Website : http://www.schiemann.uni-bonn.de
Short description:
    Trilateration is a plugin for the PyMOL Molecular Graphics System. It provides a set of
    mathematical tools which can be used to locate spin centers in biomolecules by means of EPR 
    spectroscopy.
'''

####################################################################################################
'''Import libraries'''
import os
import sys
import threading
import string
import math
import numpy
from numpy import *
from numpy.linalg import *
from pymol import *
from pymol import stored
from pymol.cgo import *
import wx
import wx.grid
from wx.lib.pubsub import Publisher
import  wx.lib.mixins.listctrl  as  listmix

####################################################################################################
class trilat(wx.Frame):
    ''''''
    
    #-----------------------------------------------------------------------------------------------
    def __init__(self, parent, ID, title):
        ''''''
        wx.Frame.__init__(self, parent, ID, title, wx.DefaultPosition, wx.Size(700,525))
        self.VariablesInitialization()
        Publisher().subscribe(self.OnImportParam, ("parameters"))
        self.UserInterface()
        self.Centre()
        self.Show()

    #-----------------------------------------------------------------------------------------------
    def VariablesInitialization(self):
        '''Initialize variables'''
        # System parameters
        self.maxNumOfLabels = 10
        self.curLabel = 0
        self.inputAccepted = False
        # Data arrays
        self.labelNames = [None for col in range(self.maxNumOfLabels)]
        for index in range(self.maxNumOfLabels):
            self.labelNames[index] = "Label_" + str(index+1)
        self.labelCoordMean = [[None for col in range(3)] for row in range(self.maxNumOfLabels)]
        self.labelCoordStd = [[None for col in range(3)] for row in range(self.maxNumOfLabels)]
        self.distMean = [None for col in range(self.maxNumOfLabels)]
        self.distStd = [None for col in range(self.maxNumOfLabels)]
        self.numOfLabels = 0
        self.targetCoordMean = [[None for col in range(3)]]
        self.targetCoordStd = [[None for col in range(3)]]
        self.chiSquare = None
        self.numOfIter = None
        # Calculation parameters (Preferences menu)
        self.coordFactor = 1.0
        self.distFactor = 10.0
        self.calcMode = 0
        self.maxNumOfIter = 10000
        self.minChiSquare = 0.000001
        self.lambdaFirst = 0.001
        self.lambdaStep = 10
        self.ellipsoidColor = [1.00, 0.43, 0.00] # orange
        self.spheresColor = [0.00, 1.00, 1.00] # blue
        self.spheresTransparency = 0.3
        self.ellipsoidTransparency = 0.0
        self.confidenceLevel = 1
    
    #-----------------------------------------------------------------------------------------------
    def OnImportParam(self, msg):
        '''Import parameters fron Preferences menu'''
        self.calcMode = msg.data[0]
        self.maxNumOfIter = msg.data[1]
        self.minChiSquare = msg.data[2]
        self.lambdaFirst = msg.data[3]
        self.lambdaStep = msg.data[4]
        self.coordFactor = msg.data[5]
        self.distFactor = msg.data[6]
        self.spheresColor = msg.data[7]
        self.spheresTransparency = msg.data[8]
        self.ellipsoidColor = msg.data[9]
        self.ellipsoidTransparency = msg.data[10]
        self.confidenceLevel = msg.data[11]
        
    #-----------------------------------------------------------------------------------------------     
    def UserInterface(self):
        '''Graphical User Interface'''
        panel = wx.Panel(self, -1)
        panel.SetFont( wx.Font(9, wx.DEFAULT, wx.NORMAL, wx.NORMAL) )
        
        # MAIN SIZER
        sizer = wx.GridBagSizer(hgap=1, vgap=1)
        sizer.AddGrowableCol(1)
        panel.SetSizer(sizer)
        
        # INPUT BOX
        inputBox = wx.StaticBox(panel, label='')
        sizerI1 = wx.StaticBoxSizer(inputBox, wx.VERTICAL)
        sizer.Add(sizerI1, pos=(0, 0), flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.BOTTOM, border=5)
        # Header
        headerInput = wx.StaticText(panel, label="INPUT")
        headerInput.SetForegroundColour( wx.Color(112,5,0) )
        headerInput.SetFont(wx.Font(11, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        sizerI1.Add(headerInput, flag=wx.ALIGN_CENTER_HORIZONTAL|wx.TOP, border=10)
        # Input data table 
        self.table = wx.grid.Grid(panel, -1, size=(-1,-1))
        rowsNum = self.maxNumOfLabels + 2
        self.table.CreateGrid(rowsNum, 6)       
        # Hide initial column labels and row labels 
        colLabels = ['' for col in range(6)]
        rowLabels = ['' for row in range(rowsNum)]
        for col in range( len(colLabels) ):
            self.table.SetColLabelValue(col, colLabels[col])
        for row in range( len(rowLabels) ):
            self.table.SetRowLabelValue(row, rowLabels[row])
        self.table.SetRowLabelSize(1)
        self.table.SetColLabelSize(1)
        # Create new column labels and row labels
        self.table.SetCellSize(0, 0, 2, 1)
        self.table.SetCellValue(0, 0, 'Label name')  
        self.table.SetCellSize(0, 1, 1, 3)
        self.table.SetCellValue(0, 1, 'Label coordinates ' + u'(\u00c5)')
        self.table.SetCellValue(1, 1, 'x')
        self.table.SetCellValue(1, 2, 'y')
        self.table.SetCellValue(1, 3, 'z')  
        self.table.SetCellSize(0, 4, 1, 2)
        self.table.SetCellValue(0, 4, 'Distance ' + u'(\u00c5)')
        self.table.SetCellValue(1, 4, 'mean')
        self.table.SetCellValue(1, 5, 'std')      
        # Do not edit column labels and row labels
        for col in range(6):
            for row in range(2):
                self.table.SetReadOnly(row, col, isReadOnly=True)
        # Hide white borders
        panelcolour = wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNFACE) 
        self.table.SetDefaultCellBackgroundColour(panelcolour) 
        for row in range(2,rowsNum): 
            for col in range(6): 
                self.table.SetCellBackgroundColour(row, col, wx.WHITE) 
        self.table.Refresh()
        # Set default text format
        self.table.SetDefaultCellAlignment(wx.ALIGN_CENTRE,wx.ALIGN_CENTRE)
        self.table.SetLabelFont( wx.Font(9, wx.DEFAULT, wx.NORMAL, wx.NORMAL) )
        # Set column and row size
        self.table.SetColSize(0, 100)
        for col in range(1, 6):
            self.table.SetColSize(col, 55)
        for row in range(0, rowsNum):
            self.table.SetRowSize(row, 22)    
        # Preset names of labels
        for row in range(2, rowsNum):
            self.table.SetCellValue(row, 0, self.labelNames[row-2])  
        # Table commands
        self.table.Bind(wx.grid.EVT_GRID_CELL_LEFT_CLICK, self.OnInputTableCellSelect)
        self.table.Bind(wx.grid.EVT_GRID_CELL_LEFT_DCLICK, self.OnInputTableCellSelect) 
        self.table.Bind(wx.grid.EVT_GRID_CELL_CHANGE, self.OnInputTableCellChange)
        # Attach to sizer
        sizerI1.Add(self.table, proportion=1, flag=wx.TOP|wx.LEFT|wx.RIGHT, border=10)
        # Load list with labels
        sizerI11 = wx.BoxSizer(wx.HORIZONTAL)
        sizerI1.Add(sizerI11, proportion=0, flag=wx.LEFT|wx.RIGHT, border=10)
        importCoordFromPymol = wx.Button(panel, label='Import from PyMOL')
        importCoordFromPymol.Bind(wx.EVT_BUTTON, self.OnImportCoordFromPymol)
        sizerI11.Add(importCoordFromPymol, flag=wx.RIGHT, border=5)
        # Load coordinates from data file
        loadCoordFromFile = wx.Button(panel, label='Load coordinates')
        loadCoordFromFile.Bind(wx.EVT_BUTTON, self.OnLoadCoordFromFile)
        sizerI11.Add(loadCoordFromFile, flag=wx.RIGHT, border=30)
        # Load distances from data file
        loadDistFromFile = wx.Button(panel, label='Load distances')
        loadDistFromFile.Bind(wx.EVT_BUTTON, self.OnLoadDistFromFile)
        sizerI11.Add(loadDistFromFile)
        # Import from PyMOL
        sizerI12 = wx.BoxSizer(wx.HORIZONTAL)
        sizerI1.Add(sizerI12, flag=wx.ALL, border=10)
        sizerI12.Add(wx.StaticText(panel, label="List of labels:"), flag=wx.ALIGN_CENTER_VERTICAL|wx.RIGHT, border=5)
        self.labelList = wx.ComboBox(panel, choices=[], size = [200, 20], style=wx.CB_READONLY)
        sizerI12.Add(self.labelList, flag=wx.RIGHT, border=5)
        self.loadFromPymol = wx.Button(panel, label='Load')
        self.loadFromPymol.Bind(wx.EVT_BUTTON, self.OnLoadFromPymol)
        sizerI12.Add(self.loadFromPymol)
        self.labelList.Disable()
        self.loadFromPymol.Disable()
        # General commands 
        sizerI13 = wx.BoxSizer(wx.HORIZONTAL)
        acceptInput = wx.Button(panel, label='Accept')
        acceptInput.Bind(wx.EVT_BUTTON, self.OnAcceptInput)
        sizerI13.Add(acceptInput)
        clearInput = wx.Button(panel, label='Clear All')
        clearInput.Bind(wx.EVT_BUTTON, self.OnClearInput)
        sizerI13.Add(clearInput, flag=wx.LEFT, border=10)
        sizerI1.Add(sizerI13, proportion=0, flag=wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, border=10)
        
        # OUTPUT BOX
        outputBox = wx.StaticBox(panel, label='')
        sizerO1 = wx.StaticBoxSizer(outputBox, wx.VERTICAL)
        sizer.Add(sizerO1, pos=(0, 1), flag=wx.EXPAND|wx.ALL, border=5)
        # Header
        headerOutput = wx.StaticText(panel, label='OUTPUT')
        headerOutput.SetForegroundColour( wx.Color(112,5,0) )
        headerOutput.SetFont(wx.Font(11, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        sizerO1.Add(headerOutput, flag=wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, border=10)
        # Computation
        copmute = wx.Button(panel, label='Compute')
        copmute.Bind(wx.EVT_BUTTON, self.OnCompute)
        sizerO1.Add(copmute, flag=wx.ALIGN_CENTER_HORIZONTAL|wx.BOTTOM, border=14)
        # Otput data table
        self.table1 = wx.grid.Grid(panel, -1, size=(-1,-1))
        self.table1.CreateGrid(4, 4)       
        # Hide initial column labels and row labels 
        colLabels = ['' for col in range(4)]
        rowLabels = ['' for row in range(4)]
        for col in range( len(colLabels) ):
            self.table1.SetColLabelValue(col, colLabels[col])
        for row in range( len(rowLabels) ):
            self.table1.SetRowLabelValue(row, rowLabels[row])
        self.table1.SetRowLabelSize(1)
        self.table1.SetColLabelSize(1)
        # Create new column labels and row labels
        self.table1.SetCellSize(0, 0, 2, 1) 
        self.table1.SetCellSize(0, 1, 1, 3)
        self.table1.SetCellValue(0, 1, 'Target coordinates ' + u'(\u00c5)')
        self.table1.SetCellValue(1, 1, 'x')
        self.table1.SetCellValue(1, 2, 'y')
        self.table1.SetCellValue(1, 3, 'z')  
        self.table1.SetCellValue(2, 0, 'mean')
        self.table1.SetCellValue(3, 0, 'std')      
        # Do not edit column labels and row labels
        for col in range(4):
            for row in range(4):
                self.table1.SetReadOnly(row, col, isReadOnly=True)
        # Hide white borders
        panelcolour = wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNFACE) 
        self.table1.SetDefaultCellBackgroundColour(panelcolour) 
        for row in range(2,4): 
            for col in range(1,4): 
                self.table1.SetCellBackgroundColour(row, col, wx.WHITE) 
        self.table1.Refresh()
        # Set default text format
        self.table1.SetDefaultCellAlignment(wx.ALIGN_CENTRE,wx.ALIGN_CENTRE)
        self.table1.SetLabelFont( wx.Font(9, wx.DEFAULT, wx.NORMAL, wx.NORMAL) )
        # Set column and row size
        for col in range(4):
            self.table1.SetColSize(col, 55)
        for row in range(4):
            self.table1.SetRowSize(row, 22)    
        # Attach to sizer
        sizerO1.Add(self.table1, flag=wx.ALL, border=10)
        # Statistics
        sizerO11 = wx.GridBagSizer(hgap=5, vgap=5)
        sizerO1.Add(sizerO11, flag=wx.LEFT|wx.RIGHT|wx.BOTTOM, border=10)
        sizerO11.Add(wx.StaticText(panel, label='Chi-square value:'), pos=(0, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        sizerO11.Add(wx.StaticText(panel, label='Number of iterations:'), pos=(1, 0), flag=wx.ALIGN_CENTER_VERTICAL)
        self.chiSquareBox = wx.TextCtrl(panel, size=(80, 20), style=wx.TE_CENTRE)
        sizerO11.Add(self.chiSquareBox, pos=(0, 1))
        self.numOfIterBox = wx.TextCtrl(panel, size=(80, 20), style=wx.TE_CENTRE)
        sizerO11.Add(self.numOfIterBox, pos=(1, 1))
        # General commands  
        sizerO12 = wx.BoxSizer(wx.HORIZONTAL)
        sizerO1.Add(sizerO12, flag=wx.ALIGN_CENTER_HORIZONTAL|wx.TOP, border=10)
        exportToPymol = wx.Button(panel, label='Export to PyMOL')
        exportToPymol.Bind(wx.EVT_BUTTON, self.OnExportToPymol)
        sizerO12.Add(exportToPymol)
        clearOutput =  wx.Button(panel, label='Clear')
        clearOutput.Bind(wx.EVT_BUTTON, self.OnClearOutput)
        sizerO12.Add(clearOutput, flag=wx.LEFT, border=10)

        # CREATE MENU BAR
        menubar = wx.MenuBar()
        fileMenu = wx.Menu()
        openItem = wx.MenuItem(fileMenu, 11, '&Open\tCtrl+O')
        self.Bind(wx.EVT_MENU, self.OnOpen, id=11)
        fileMenu.AppendItem(openItem)
        saveItem = wx.MenuItem(fileMenu, 12, '&Save As ...\tCtrl+S')
        self.Bind(wx.EVT_MENU, self.OnSave, id=12)
        fileMenu.AppendItem(saveItem)
        closeItem = wx.MenuItem(fileMenu, 13, '&Close\tCtrl+C')
        self.Bind(wx.EVT_MENU, self.OnClose, id=13)
        fileMenu.AppendItem(closeItem)
        fileMenu.AppendSeparator()
        prefItem = wx.MenuItem(fileMenu, 14, '&Preferences\tCtrl+P')
        self.Bind(wx.EVT_MENU, self.OnOpenPref, id=14)
        fileMenu.AppendItem(prefItem)
        fileMenu.AppendSeparator()
        exitItem = wx.MenuItem(fileMenu, 15, '&Exit\tAlt+F4')
        self.Bind(wx.EVT_MENU, self.OnExit, id=15)
        fileMenu.AppendItem(exitItem)
        pymolMenu = wx.Menu()
        loadItem = wx.MenuItem(pymolMenu, 21 , '&Load structure to PyMOL')
        self.Bind(wx.EVT_MENU, self.OnLoad, id=21)
        pymolMenu.AppendItem(loadItem)
        removeItem = wx.MenuItem(pymolMenu, 22, '&Remove structure from PyMOL')
        self.Bind(wx.EVT_MENU, self.OnRemove, id=22)
        pymolMenu.AppendItem(removeItem)
        mtsslWizardItem = wx.MenuItem(pymolMenu, 23, '&Run MtsslWizard')
        self.Bind(wx.EVT_MENU, self.OnRunMtsslWizard, id=23)
        pymolMenu.AppendItem(mtsslWizardItem)
        helpMenu = wx.Menu()
        aboutItem = wx.MenuItem(helpMenu, 31, '&About\tCtrl+A')
        self.Bind(wx.EVT_MENU, self.OnAbout, id=31)
        helpMenu.AppendItem(aboutItem)
        helpMenu.AppendSeparator()
        helpItem = wx.MenuItem(helpMenu, 32, '&Help\tF1')
        #self.Bind(wx.EVT_MENU, self.OnHelp, id=32)
        helpMenu.AppendItem(helpItem)
        menubar.Append(fileMenu, '&File')
        menubar.Append(pymolMenu, '&PyMOL')
        menubar.Append(helpMenu, '&Help')
        self.SetMenuBar(menubar)
        
        # CREATE STATUS BAR
        self.statusbar = self.CreateStatusBar()
        
    #-----------------------------------------------------------------------------------------------
    def OnInputTableCellSelect(self, event):
        '''Reads the number of selected cell'''
        curRow = event.GetRow()
        if (curRow >= 2):
            # Highlight selected label
            self.table.SetCellBackgroundColour((self.curLabel + 2), 0, wx.WHITE)
            self.table.SetCellBackgroundColour(curRow, 0, wx.Colour(255,168,0))
            self.table.Refresh()
            # Set number of current label
            self.curLabel = curRow - 2
        
        event.Skip()
    
    #-----------------------------------------------------------------------------------------------
    def OnInputTableCellChange(self, event):
        '''Saves the data entered to the table in input arrays'''
        curRow = event.GetRow()
        curCol = event.GetCol()
        if (curCol == 0):
            self.labelNames[self.curLabel] = self.table.GetCellValue(curRow,curCol)
        elif (curCol == 1):
            try:
                self.labelCoordMean[self.curLabel][0] = float( self.table.GetCellValue(curRow,curCol) )
                self.statusbar.SetStatusText('Entering new input data... Click "Accept" button after all input data has been entered.')
            except:
                self.labelCoordMean[self.curLabel][0] = None
                self.statusbar.SetStatusText('Incorrect input data!')
        elif (curCol == 2):
            try:
                self.labelCoordMean[self.curLabel][1] = float( self.table.GetCellValue(curRow,curCol) )
                self.statusbar.SetStatusText('Entering new input data... Click "Accept" button after all input data has been entered.')
            except:
                self.labelCoordMean[self.curLabel][1] = None
                self.statusbar.SetStatusText('Incorrect input data!')
        elif (curCol == 3):
            try:
                self.labelCoordMean[self.curLabel][2] = float( self.table.GetCellValue(curRow,curCol) )
                self.statusbar.SetStatusText('Entering new input data... Click "Accept" button after all input data has been entered.')
            except:
                self.labelCoordMean[self.curLabel][2] = None
                self.statusbar.SetStatusText('Incorrect input data!')
        elif (curCol == 4): 
            try:
                self.distMean[self.curLabel] = float( self.table.GetCellValue(curRow,curCol) )
                self.statusbar.SetStatusText('Entering new input data... Click "Accept" button after all input data has been entered.')
            except:
                self.distMean[self.curLabel] = None
                self.statusbar.SetStatusText('Incorrect input data!')
        elif (curCol == 5):
            try:
                self.distStd[self.curLabel] = float( self.table.GetCellValue(curRow,curCol) )
                self.statusbar.SetStatusText('Entering new input data... Click "Accept" button after all input data has been entered.')
            except:
                self.distStd[self.curLabel] = None
                self.statusbar.SetStatusText('Incorrect input data!')
        # Update table 
        self.UpdateInputTable() 
        # Input data doesn't accepted
        self.inputAccepted = False   
        
        event.Skip()
    
    #-----------------------------------------------------------------------------------------------
    def OnLoadCoordFromFile(self, event):
        '''Loads label coordinates from data file'''
        dlg = wx.FileDialog(self, message="Choose a file", defaultDir=os.getcwd(), defaultFile="", wildcard="Text Files (.txt, .dat)|*.txt;*.dat", style=wx.OPEN | wx.CHANGE_DIR) 
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            dataFile = open(path, 'r')
            # Load the array of coordinates
            coord = []
            for (i, line) in enumerate(dataFile):
                coordLine = line.split()
                coord += [[float(coordLine[0])*self.coordFactor, float(coordLine[1])*self.coordFactor, float(coordLine[2])*self.coordFactor]]
            # Calculate average coordinates and their standard deviations
            self.labelCoordMean[self.curLabel], self.labelCoordStd[self.curLabel] = AverageCoordCalc(coord)
            # Update table
            self.UpdateInputTable() 
            # Input data doesn't accepted
            self.inputAccepted = False
            # Statusbar message
            self.statusbar.SetStatusText('Coordinates of label no. '+str(self.curLabel+1)+' were imported from data file! Click "Accept" button after all input data has been entered.')
            
        event.Skip()
    
    #-----------------------------------------------------------------------------------------------
    def OnLoadDistFromFile(self, event):
        '''Load the distance from data file'''
        dlg = wx.FileDialog(self, message="Choose a file", defaultDir=os.getcwd(), defaultFile="", wildcard="Text Files (.txt, .dat)|*.txt;*.dat", style=wx.OPEN | wx.CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            dataFile = open(path, 'r')
            # Load a distance distribution
            dist = []
            density = []
            for (i, line) in enumerate(dataFile):
                distrib = line.split()
                dist += [ float(distrib[0]) * self.distFactor ]
                density += [ float(distrib[1]) ]
            # Calculate distance and its standard deviation
            self.distMean[self.curLabel], self.distStd[self.curLabel] = AverageDistCalc(dist, density)
            # Update table
            self.UpdateInputTable() 
            # Input data doesn't accepted
            self.inputAccepted = False
            # Statusbar message
            self.statusbar.SetStatusText('Distance between target and label no. '+str(self.curLabel+1)+' was imported from data file! Click "Accept" button after all input data has been entered.')
                      
        event.Skip()
 
    #-----------------------------------------------------------------------------------------------
    def UpdateInputTable(self):
        '''Fills in the table by data'''
        for i in range(self.maxNumOfLabels):
            self.table.SetCellValue(i+2, 0, self.labelNames[i])
            if not ( self.labelCoordMean[i][0] == None):
                self.table.SetCellValue(i+2, 1, '{0:.2f}'.format(self.labelCoordMean[i][0]) )
            else:
                self.table.SetCellValue(i+2, 1, '')
            if not ( self.labelCoordMean[i][1] == None):
                self.table.SetCellValue(i+2, 2, '{0:.2f}'.format(self.labelCoordMean[i][1]) )
            else:
                self.table.SetCellValue(i+2, 2, '')
            if not ( self.labelCoordMean[i][2] == None):
                self.table.SetCellValue(i+2, 3, '{0:.2f}'.format(self.labelCoordMean[i][2]) )
            else:
                self.table.SetCellValue(i+2, 3, '')
            if not ( self.distMean[i] == None):
                self.table.SetCellValue(i+2, 4, '{0:.2f}'.format(self.distMean[i]) )
            else:
                self.table.SetCellValue(i+2, 4, '')
            if not ( self.distStd[i] == None):
                self.table.SetCellValue(i+2, 5, '{0:.2f}'.format(self.distStd[i])) 
            else:
                self.table.SetCellValue(i+2, 5, '') 
            # Highlight complete sets of data
            if not (( self.labelCoordMean[i][0] == None) | ( self.labelCoordMean[i][1] == None) | ( self.labelCoordMean[i][2] == None) | ( self.distMean[i] == None) | ( self.distStd[i] == None)):   
                for j in range(1,6):
                    self.table.SetCellBackgroundColour(i+2, j, wx.Colour(127, 255, 0))
            else:
                for j in range(1,6):
                    self.table.SetCellBackgroundColour(i+2, j, wx.WHITE)
    
    #-----------------------------------------------------------------------------------------------
    def OnAcceptInput(self, event):
        '''Accepts input data if input data is correct'''
        # Number of defined labels
        N = 0
        # Rearrange input data arrays & calculate number of defined labels
        for i in range(self.maxNumOfLabels):
            if not (( self.labelCoordMean[i][0] == None) | ( self.labelCoordMean[i][1] == None) | ( self.labelCoordMean[i][2] == None) | ( self.distMean[i] == None) | ( self.distStd[i] == None)):
                # Rearrange input data arrays
                tempLabelName = self.labelNames[N]
                tempLabelCoordMean = self.labelCoordMean[N]
                tempLabelCoordStd = self.labelCoordMean[N]
                tempDistMean = self.distMean[N]
                tempDistStd = self.distStd[N]
                
                self.labelNames[N] = self.labelNames[i]
                self.labelCoordMean[N] = self.labelCoordMean[i]
                self.labelCoordStd[N] = self.labelCoordStd[i]
                self.distMean[N] = self.distMean[i]
                self.distStd[N] = self.distStd[i]
                
                self.labelNames[i] = tempLabelName
                self.labelCoordMean[i] = tempLabelCoordMean
                self.labelCoordStd[i] = tempLabelCoordStd
                self.distMean[i] = tempDistMean
                self.distStd[i] = tempDistStd
                
                N += 1
        # Update table
        self.UpdateInputTable()
        # Check in the number of defined labels
        if (N >= 4):
            # Input data does accepted
            self.inputAccepted = True
            self.numOfLabels = N
            # Statusbar message
            self.statusbar.SetStatusText('Input data is prepared!')
        else:
            # Input data doesn't accepted
            self.inputAccepted = False
            self.numOfLabels = N
            # Statusbar message
            self.statusbar.SetStatusText('The Number of the determined labels should be more than 3 to start the calculations!')   
        event.Skip()
    
    #-----------------------------------------------------------------------------------------------
    def OnClearInput(self, event):
        '''Clears all input data '''
        # Clear input arrays
        self.labelNames = [None for col in range(self.maxNumOfLabels)]
        for index in range(self.maxNumOfLabels):
            self.labelNames[index] = "Label_" + str(index+1)
        self.labelCoordMean = [[None for col in range(3)] for row in range(self.maxNumOfLabels)]
        self.labelCoordStd = [[None for col in range(3)] for row in range(self.maxNumOfLabels)]
        self.distMean = [None for col in range(self.maxNumOfLabels)]
        self.distStd = [None for col in range(self.maxNumOfLabels)]
        # Update table
        self.UpdateInputTable()
        # Input data doesn't accepted
        self.inputAccepted = False 
             
        event.Skip() 
        
    #-----------------------------------------------------------------------------------------------  
    def OnCompute(self, event):
        '''Computes target coordinates'''
        if (self.inputAccepted == False):
            # Statusbar message
            self.statusbar.SetStatusText('Correct input data is required!')
        else:
            # Linear least squares, Singular value decomposition
            self.targetCoordMean,self.targetCoordStd,self.chiSquare = SingularValueDecomposition(self.numOfLabels, 
                                                                                                 self.labelCoordMean, 
                                                                                                 self.distMean, 
                                                                                                 self.distStd)
            self.numOfIter = None
            # Statusbar message
            self.statusbar.SetStatusText('Target coordinates were calculated by means of SVD!')
            # Fill in the otput fields
            self.table1.SetCellValue(2, 1, '{0:.1f}'.format(self.targetCoordMean[0][0]) )
            self.table1.SetCellValue(2, 2, '{0:.1f}'.format(self.targetCoordMean[0][1]) )
            self.table1.SetCellValue(2, 3, '{0:.1f}'.format(self.targetCoordMean[0][2]) )
            self.table1.SetCellValue(3, 1, '{0:.1f}'.format(self.targetCoordStd[0][0]) )
            self.table1.SetCellValue(3, 2, '{0:.1f}'.format(self.targetCoordStd[0][1]) )
            self.table1.SetCellValue(3, 3, '{0:.1f}'.format(self.targetCoordStd[0][2]) )
            self.chiSquareBox.SetValue( '{0:.3f}'.format(self.chiSquare) )
            self.numOfIterBox.SetValue( '1' )
            if self.calcMode == 0:
                self.targetCoordMean,self.targetCoordStd,self.chiSquare,self.numOfIter = InverseHessian(self.numOfLabels, 
                                                                                                        self.labelCoordMean, 
                                                                                                        self.distMean, 
                                                                                                        self.distStd,
                                                                                                        self.targetCoordMean,
                                                                                                        self.minChiSquare,
                                                                                                        self.maxNumOfIter)
            elif self.calcMode == 1:
                self.targetCoordMean,self.targetCoordStd,self.chiSquare,self.numOfIter = LevenbergMarquardt(self.numOfLabels, 
                                                                                                            self.labelCoordMean, 
                                                                                                            self.distMean, 
                                                                                                            self.distStd,
                                                                                                            self.targetCoordMean,
                                                                                                            self.minChiSquare,
                                                                                                            self.maxNumOfIter,
                                                                                                            self.lambdaFirst,
                                                                                                            self.lambdaStep)
            if self.numOfIter == self.maxNumOfIter:
                self.statusbar.SetStatusText('Can not find a solution for a given set of input data! (Try to increase number of iterations.)')
            else:
                self.statusbar.SetStatusText('Target coordinates were calculated by use of nonlinear least squares algorithm!')
            # Fill in output fields
            self.table1.SetCellValue(2, 1, '{0:.1f}'.format(self.targetCoordMean[0][0]) )
            self.table1.SetCellValue(2, 2, '{0:.1f}'.format(self.targetCoordMean[0][1]) )
            self.table1.SetCellValue(2, 3, '{0:.1f}'.format(self.targetCoordMean[0][2]) )
            self.table1.SetCellValue(3, 1, '{0:.1f}'.format(self.targetCoordStd[0][0]) )
            self.table1.SetCellValue(3, 2, '{0:.1f}'.format(self.targetCoordStd[0][1]) )
            self.table1.SetCellValue(3, 3, '{0:.1f}'.format(self.targetCoordStd[0][2]) )
            self.chiSquareBox.SetValue( '{0:.3f}'.format(self.chiSquare) )
            self.numOfIterBox.SetValue( '{0:d}'.format(self.numOfIter) )

        event.Skip()    
        
    #----------------------------------------------------------------------------------------------- 
    def OnClearOutput(self, event):
        '''Clear all output data '''
        self.targetCoordMean = [[None for col in range(3)]]
        self.targetCoordStd = [[None for col in range(3)]]
        self.chiSquare = None
        self.numOfIter = None
        self.table1.SetCellValue(2, 1, '' )
        self.table1.SetCellValue(2, 2, '' )
        self.table1.SetCellValue(2, 3, '' )
        self.table1.SetCellValue(3, 1, '' )
        self.table1.SetCellValue(3, 2, '' )
        self.table1.SetCellValue(3, 3, '' )
        self.chiSquareBox.SetValue( '' )
        self.numOfIterBox.SetValue( '' )
        # Remove figures from PyMOL
        try:
            cmd.delete('trilateration')
            cmd.delete('target')
        except:
            pass
        
        event.Skip()
    
    #-----------------------------------------------------------------------------------------------
    def OnImportCoordFromPymol(self, event):
        '''Import list of labels from current PyMOL session'''
        self.labelList.Enable()
        self.loadFromPymol.Enable()
        # Load a dictionary with labels and their coordinates from PyMOL
        items = getPymolObjects("*_M-T-S-S-L")
        # Write names of labels to the combobox list
        self.labelList.Clear()
        for key in sorted(items):
            self.labelList.Append(key)
        # Statusbar message
        Publisher().sendMessage(('change_statusbar'), 'Choose label and push "Load".')
        
        event.Skip()
    
    #-----------------------------------------------------------------------------------------------
    def OnLoadFromPymol(self, event):
        '''Import names and coordinates of labels from PyMOL'''
        #try:
        # Read the name of the label
        s = self.labelList.GetStringSelection()
        # Load the coordinates of epr-active center of the spin label
        coordO1 = []
        coordN1 = []
        cmd.iterate_state(0, s+' and name O1', 'coordO1.append([x,y,z])', space=locals(), atomic=0)
        cmd.iterate_state(0, s+' and name N1', 'coordN1.append([x,y,z])', space=locals(), atomic=0)
        coord = 0.5 * numpy.add(coordO1,coordN1)      
        # Save the name of spin label
        self.labelNames[self.curLabel] = s
        # Calculate average label coordinates and their standard deviations
        self.labelCoordMean[self.curLabel], self.labelCoordStd[self.curLabel] =  AverageCoordCalc(coord)
        # Update table
        self.UpdateInputTable() 
        # Statusbar message
        self.statusbar.SetStatusText('Coordinates of '+s+' label were imported from PyMOL. Click "Accept" button after all input data has been entered.')
        #except:
            # Statusbar message
            #self.statusbar.SetStatusText('No label was chosen!')
        # Input data doesn't accepted
        self.inputAccepted = False
        
        event.Skip()
    
    #-----------------------------------------------------------------------------------------------
    def OnExportToPymol(self, event):
        '''Creates a plot in PyMOL'''
        if (self.inputAccepted == False):
            # Plot spheres
            PlotSpheres(self.numOfLabels,self.labelCoordMean,self.distMean,
                        self.spheresColor, self.spheresTransparency)
        else:
            # Plot spheres
            PlotSpheres(self.numOfLabels,self.labelCoordMean,self.distMean,
                        self.spheresColor, self.spheresTransparency)
            # Plot ellipsoid
            try:
                PlotEllipsoid(self.targetCoordMean[0][0],self.targetCoordMean[0][1],self.targetCoordMean[0][2],
                              self.targetCoordStd[0][0]*self.confidenceLevel,self.targetCoordStd[0][1]*self.confidenceLevel,self.targetCoordStd[0][2]*self.confidenceLevel,
                              self.ellipsoidColor,self.ellipsoidTransparency)
            except:
                pass
        
        event.Skip()
    
    #----------------------------------------------------------------------------------------------- 
    def OnSave(self, event):
        '''Saves current session'''
        if (self.inputAccepted == False):
            self.statusbar.SetStatusText('Data set is not complete!')
        else:
            dialog = wx.FileDialog(None, message="Save file as ...", defaultDir=os.getcwd(), defaultFile="", wildcard="*.dat", style = wx.SAVE | wx.OVERWRITE_PROMPT)
            if dialog.ShowModal() == wx.ID_OK:
                path = dialog.GetPath()
                # Content of a data file
                title = "*** Trilateration program data file ***"    
                object = "* Protein name:"
                numOfLabels = "* Number of labels:"
                labels = '''* Label parameters:
col 1: Name
col 2,3,4: (x,y,z) coordinates
col 5: Label-target distance
col 6: Standard deviation of label-target distance'''
                settings = '''* Calculation settings:
line 1: Inverse-Hessian method vs Levenberg-Marquardt method, 1 = ON and 0 = OFF
line 2: Maximal number of iterations
line 3: Minimal chi-square value
line 4: First lambda value
line 5: Lambda increment'''
                target = '''* Target parameters:
line 1: Name
line 2: (x,y,z) coordinates
line 3: Standard deviation of (x,y,z) coordinates
line 4: Obtained chi-square value
line 5: Obtained number of iterations'''
                # Write data to file
                dataFile = open(path, 'w')
                dataFile.write(title + '\n\n')
                dataFile.write(object + '\n')
                dataFile.write('name' + '\n')
                dataFile.write(numOfLabels + '\n')
                dataFile.write('{0:d}'.format(self.numOfLabels) + '\n')
                dataFile.write(labels+'\n')  
                for i in range(self.numOfLabels):
                    dataFile.write('{0:16s}'.format(self.labelNames[i]) + ' ')
                    dataFile.write('{0:8.2f}'.format(self.labelCoordMean[i][0]) + ' ')
                    dataFile.write('{0:8.2f}'.format(self.labelCoordMean[i][1]) + ' ')
                    dataFile.write('{0:8.2f}'.format(self.labelCoordMean[i][2]) + ' ')
                    dataFile.write('{0:8.2f}'.format(self.distMean[i]) + ' ')
                    dataFile.write('{0:8.2f}'.format(self.distStd[i]) + '\n')
                dataFile.write(settings + '\n') 
                if (self.calcMode == 0):
                    dataFile.write('1 0' + '\n')
                else:
                    dataFile.write('0 1' + '\n')
                dataFile.write(str(self.maxNumOfIter) + '\n')
                dataFile.write(str(self.minChiSquare) + '\n')
                dataFile.write(str(self.lambdaFirst) + '\n')
                dataFile.write(str(self.lambdaStep) + '\n')
                dataFile.write(target + '\n')
                dataFile.write('name' + '\n')
                dataFile.write('{0:8.2f}'.format(self.targetCoordMean[0][0]) + ' ')
                dataFile.write('{0:8.2f}'.format(self.targetCoordMean[0][1]) + ' ')
                dataFile.write('{0:8.2f}'.format(self.targetCoordMean[0][2]) + '\n')
                dataFile.write('{0:8.2f}'.format(self.targetCoordStd[0][0]) + ' ')
                dataFile.write('{0:8.2f}'.format(self.targetCoordStd[0][1]) + ' ')
                dataFile.write('{0:8.2f}'.format(self.targetCoordStd[0][2]) + '\n')
                dataFile.write('{0:.3f}'.format(self.chiSquare) + '\n')
                dataFile.write('{0:d}'.format(self.numOfIter))
                dataFile.close()
                # Statusbar message
                self.statusbar.SetStatusText('Current session is saved!')
            dialog.Destroy()
    
    #-----------------------------------------------------------------------------------------------
    def OnOpen(self, event):
        '''Opens saved session'''
        dlg = wx.FileDialog(self, message="Choose a file", defaultDir=os.getcwd(), defaultFile="", wildcard="*.dat", style=wx.OPEN | wx.CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            dataFile = open(path, 'r')
            len = 8 # lenght of numerical symbols in the file
            for (i, line) in enumerate(dataFile):
                if (i == 5):
                    self.numOfLabels = int(line)
                for j in range(self.numOfLabels):
                    if (i == (11+j)):
                        self.labelNames[j] = line[0:2*len].replace(' ','')
                        self.labelCoordMean[j][0] = float( line[2*len+1:3*len+1].replace(' ','') )
                        self.labelCoordMean[j][1] = float( line[3*len+2:4*len+2].replace(' ','') )
                        self.labelCoordMean[j][2] = float( line[4*len+3:5*len+3].replace(' ','') )
                        self.distMean[j] = float( line[5*len+4:6*len+4].replace(' ','') )
                        self.distStd[j] = float( line[6*len+5:7*len+5].replace(' ','') )   
                if (i == (17+self.numOfLabels)):
                    if (line[0] == '1'):
                        self.calcMode = 0
                    else:
                        self.calcMode = 1
                if (i == (18+self.numOfLabels)):
                    self.maxNumOfIter = int(line.replace(' ',''))   
                if (i == (19+self.numOfLabels)):
                    self.minChiSquare = float(line.replace(' ',''))     
                if (i == (20+self.numOfLabels)):
                    self.lambdaFirst = float(line.replace(' ',''))  
                if (i == (21+self.numOfLabels)):
                    self.lambdaStep = float(line.replace(' ',''))  
                if (i == (29+self.numOfLabels)):
                    self.targetCoordMean[0][0] = float( line[0:len].replace(' ','') )
                    self.targetCoordMean[0][1] = float( line[len+1:2*len+1].replace(' ','') )
                    self.targetCoordMean[0][2] = float( line[2*len+2:3*len+2].replace(' ','') )
                    self.table1.SetCellValue(2, 1, '{0:.1f}'.format(self.targetCoordMean[0][0]) )
                    self.table1.SetCellValue(2, 2, '{0:.1f}'.format(self.targetCoordMean[0][1]) )
                    self.table1.SetCellValue(2, 3, '{0:.1f}'.format(self.targetCoordMean[0][2]) )
                if (i == (30+self.numOfLabels)):
                    self.targetCoordStd[0][0] = float( line[0:len].replace(' ','') )
                    self.targetCoordStd[0][1] = float( line[len+1:2*len+1].replace(' ','') )
                    self.targetCoordStd[0][2] = float( line[2*len+2:3*len+2].replace(' ','') )
                    self.table1.SetCellValue(3, 1, '{0:.1f}'.format(self.targetCoordStd[0][0]) )
                    self.table1.SetCellValue(3, 2, '{0:.1f}'.format(self.targetCoordStd[0][1]) )
                    self.table1.SetCellValue(3, 3, '{0:.1f}'.format(self.targetCoordStd[0][2]) )
                if (i == (31+self.numOfLabels)):
                    self.chiSquare = float( line.replace(' ','') )
                    self.chiSquareBox.SetValue( line.replace(' ','') )
                if (i == (32+self.numOfLabels)):    
                    self.numOfIter = float( line.replace(' ','') ) 
                    self.numOfIterBox.SetValue( line.replace(' ','') )
            # Update table
            self.UpdateInputTable()
            self.inputAccepted = True
            # Statusbar message
            self.statusbar.SetStatusText('New session was loaded from data file!')
            # Remove figures from PyMOL
            try:
                cmd.delete('trilateration')
                cmd.delete('target')
            except:
                pass
        else:
            # Statusbar message
            self.statusbar.SetStatusText('Nothing was selected!')
        dlg.Destroy()
    
    #----------------------------------------------------------------------------------------------- 
    def OnClose(self, event):
        '''Close current session'''
        # Clear input arrays
        self.labelNames = [None for col in range(self.maxNumOfLabels)]
        for index in range(self.maxNumOfLabels):
            self.labelNames[index] = "Label_" + str(index+1)
        self.labelCoordMean = [[None for col in range(3)] for row in range(self.maxNumOfLabels)]
        self.labelCoordStd = [[None for col in range(3)] for row in range(self.maxNumOfLabels)]
        self.distMean = [None for col in range(self.maxNumOfLabels)]
        self.distStd = [None for col in range(self.maxNumOfLabels)]
        # Update table
        self.UpdateInputTable()
        # Input data doesn't accepted
        self.inputAccepted = False
        # Clear output arrays
        self.targetCoordMean = [[None for col in range(3)]]
        self.targetCoordStd = [[None for col in range(3)]]
        self.chiSquare = None
        self.numOfIter = None
        # Clear output text fields
        self.table1.SetCellValue(2, 1, '' )
        self.table1.SetCellValue(2, 2, '' )
        self.table1.SetCellValue(2, 3, '' )
        self.table1.SetCellValue(3, 1, '' )
        self.table1.SetCellValue(3, 2, '' )
        self.table1.SetCellValue(3, 3, '' )
        self.chiSquareBox.SetValue( '' )
        self.numOfIterBox.SetValue( '' )
        # Remove figures from PyMOL
        try:
            cmd.delete('trilateration')
            cmd.delete('target')
        except:
            pass
        
        event.Skip()        
    
    #-----------------------------------------------------------------------------------------------
    def OnOpenPref(self, event):
        '''Preferences menu'''
        prefFrame = Preferences(self,
                                self.coordFactor, self.distFactor,
                                self.calcMode,self.maxNumOfIter,self.minChiSquare,self.lambdaFirst,self.lambdaStep,
                                self.spheresColor,self.spheresTransparency,self.ellipsoidColor,self.ellipsoidTransparency,
                                self.confidenceLevel)
        prefFrame.Show()
    
    #-----------------------------------------------------------------------------------------------
    def OnExit(self, event):
        '''Close program'''
        self.Close()
    
    #-----------------------------------------------------------------------------------------------
    def OnLoad(self, event):
        '''Load new structure to PyMOL'''
        dialog = wx.FileDialog (None, message = 'Set PDB file', style = wx.OPEN)
        if dialog.ShowModal() == wx.ID_OK:
            pdbfile = dialog.GetPath()
            cmd.load(pdbfile)
            # make PyMOL view nicer
            cmd.remove("solvent")
            cmd.show("cartoon")
            #util.cbag("all")
            # Statusbar message
            self.statusbar.SetStatusText('New structure was loaded to PyMOL!')
        else:
            # Statusbar message
            self.statusbar.SetStatusText('Nothing was selected.')
        
        dialog.Destroy()
            
    #-----------------------------------------------------------------------------------------------
    def OnRemove(self, event):
        '''Remove current structure from PyMOL'''
        cmd.delete("all")
        # Statusbar message
        self.statusbar.SetStatusText('Current structure was removeded from PyMOL!')
        
    #-----------------------------------------------------------------------------------------------
    def OnRunMtsslWizard(self, event):
        '''Run MtsslWizard'''
        try:
            sys.path.append(os.path.dirname(os.path.abspath(__file__)) +'\mtsslWizard.py')
            from mtsslWizard import MtsslWizard
            wiz = MtsslWizard()
            cmd.set_wizard(wiz)
        except:
            pass
    
    #-----------------------------------------------------------------------------------------------   
    def OnAbout(self, event):
        '''Programs info'''
        info = wx.AboutDialogInfo()
        info.SetIcon(wx.Icon(os.path.dirname(os.path.abspath(__file__)) + '/icons/about.png', wx.BITMAP_TYPE_PNG))
        info.SetName('mtsslTrilaterate')
        info.SetVersion('1.0')
        description = """
This program calculates coordinates of the point object which is defined by
distances from this object to few (more than 3) fixed points in 3D space."""
        info.SetDescription(description)
        info.SetCopyright('(C) 2012 Dinar Abdullin')
        info.SetWebSite('http://www.schiemann.uni-bonn.de')
        licence = """This program is free software. You can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation. See the GNU General Public License for more details."""
        info.SetLicence(licence)
        wx.AboutBox(info)

   
####################################################################################################
class ListCtrlLeft(wx.ListCtrl):
    '''Content of Preferences menu'''

    #-----------------------------------------------------------------------------------------------
    def __init__(self, parent, id):
        ''''''
        wx.ListCtrl.__init__(self, parent, id, style=wx.LC_REPORT|wx.LC_HRULES|wx.LC_NO_HEADER|wx.LC_SINGLE_SEL)
        self.parent = parent
        self.Bind(wx.EVT_SIZE, self.OnSize)
        self.Bind(wx.EVT_LIST_ITEM_SELECTED, self.OnSelect)
        self.InsertColumn(0, '')
        titles = ['Otput data', 'Calculations', 'Input data']
        for i in titles:
            self.InsertStringItem(0, i)
            
    #-----------------------------------------------------------------------------------------------            
    def OnSize(self, event):
        ''''''
        size = self.parent.GetSize()
        self.SetColumnWidth(0, size.x-5)
        event.Skip()
    
    #-----------------------------------------------------------------------------------------------
    def OnSelect(self, event):
        ''''''
        window = self.parent.GetGrandParent().FindWindowByName('ListControlOnRight')
        index = event.GetIndex()
        window.LoadData(index)

####################################################################################################
class ListCtrlRight(wx.grid.Grid):
    '''Names and values of parameters'''
    
    #-----------------------------------------------------------------------------------------------
    def __init__(self, parent, id,
                 coordFactor, distFactor,
                 calcMode,maxNumOfIter,minChiSquare,lambdaFirst,lambdaStep,
                 spheresColor,spheresTransparency,ellipsoidColor,ellipsoidTransparency,confidenceLevel):
        ''''''
        wx.grid.Grid.__init__(self, parent, id)
        # Create hollow grid
        self.CreateGrid(0, 2)
        self.SetRowLabelSize(0)
        self.SetColLabelSize(0)
        self.parent = parent
        size = self.parent.GetSize()
        self.SetColSize(0, size.x-100)
        self.SetColSize(1, 95)
        #self.EnableGridLines(False)
        self.Bind(wx.grid.EVT_GRID_CELL_CHANGE, self.OnTableCellChange)
        # Initialize parameters
        self.OnParamInit(coordFactor, distFactor,
                         calcMode,maxNumOfIter,minChiSquare,lambdaFirst,lambdaStep,
                         spheresColor,spheresTransparency,ellipsoidColor,ellipsoidTransparency,confidenceLevel)
    
    #-----------------------------------------------------------------------------------------------
    def OnParamInit(self,
                    coordFactor, distFactor,
                    calcMode,maxNumOfIter,minChiSquare,lambdaFirst,lambdaStep,
                    spheresColor,spheresTransparency,ellipsoidColor,ellipsoidTransparency,confidenceLevel):
        ''''''
        colorLine1 = ''
        colorLine2 = ''
        for i in range(3):
            colorLine1 += str( int(spheresColor[i]*255.0) ) + ' '
            colorLine2 += str( int(ellipsoidColor[i]*255.0) ) + ' '
        self.param = {
        'input': [['Multiply the coordinates loaded from data file by factor of', coordFactor],
                  ['Multiply the distances loaded from data file by factor of', distFactor]],
        'calc': [['Inverse-Hessian method (0) vs Lavenberg-Marquardt method (1)', calcMode],
                 ['Maximal number of iterations', maxNumOfIter],
                 ['Minimal chi-square value', minChiSquare],
                 ['Initial lambda (for Lavenberg-Marquardt method only)', lambdaFirst],
                 ['Lambda increment (for Lavenberg-Marquardt method only)', lambdaStep]],
        'output': [['Color of <trilateration> graphical object (RGB)', colorLine1],
                   ['Transparency of <trilateration> graphical object', spheresTransparency],
                   ['Color of <target> graphical object (RGB)', colorLine2],
                   ['Transparency of <target> graphical object', ellipsoidTransparency],
                   ['Confidence level', confidenceLevel]]
        }
    #-----------------------------------------------------------------------------------------------
    def LoadData(self, index):
        ''''''
        # Clear last grid
        self.ClearGrid()
        oldRows = self.GetNumberRows()
        self.DeleteRows(0,oldRows,True)
        # Select type of parameters
        self.index = index
        if self.index == 0:
            data = self.param['input']
        elif self.index == 1:
            data = self.param['calc']
        elif self.index == 2:
            data = self.param['output']
        # Create a new grid
        newRows = len(data)
        self.AppendRows(newRows, True)
        self.SetColLabelValue(0, "Name")
        self.SetColLabelValue(1, "Value")
        self.SetLabelFont( wx.Font(8, wx.DEFAULT, wx.NORMAL, wx.NORMAL) )
        self.SetRowLabelSize(0)
        self.SetColLabelSize(20)
        size = self.parent.GetSize()
        self.SetColSize(0, size.x-105)
        self.SetColSize(1, 100)
        # Fill in a new grid
        for i in range(len(data)):
            self.SetCellValue(i, 0, data[i][0])
            self.SetCellValue(i, 1, str(data[i][1]))
            self.SetReadOnly(i, 0, isReadOnly=True)
            
    #-----------------------------------------------------------------------------------------------
    def OnTableCellChange(self, event):
        ''''''
        row = event.GetRow()       
        if self.index == 0:
            if row == 0:
                self.param['input'][0][1] = float(self.GetCellValue(row,1))
            elif row == 1:
                self.param['input'][1][1] = float(self.GetCellValue(row,1))
        if self.index == 1:
            if row == 0:
                self.param['calc'][0][1] = int(self.GetCellValue(row,1))
            elif row == 1:
                self.param['calc'][1][1] = int(self.GetCellValue(row,1))
            elif row == 2:
                self.param['calc'][2][1] = float(self.GetCellValue(row,1))
            elif row == 3:
                self.param['calc'][3][1] = float(self.GetCellValue(row,1))
            elif row == 4:
                self.param['calc'][4][1] = float(self.GetCellValue(row,1))
        if self.index == 2:
            if row == 0:
                self.param['output'][0][1] = self.GetCellValue(row,1)
            elif row == 1:
                self.param['output'][1][1] = float(self.GetCellValue(row,1))
            elif row == 2:
                self.param['output'][2][1] = self.GetCellValue(row,1)
            elif row == 3:
                self.param['output'][3][1] = float(self.GetCellValue(row,1))
            elif row == 4:
                self.param['output'][4][1] = int(self.GetCellValue(row,1))
        
        event.Skip()

####################################################################################################        
class Preferences(wx.Frame):
    '''Preferences menu'''
    
    #-----------------------------------------------------------------------------------------------
    def __init__(self, parent, 
                 coordFactor, distFactor,
                 calcMode,maxNumOfIter,minChiSquare,lambdaFirst,lambdaStep,
                 spheresColor,spheresTransparency,ellipsoidColor,ellipsoidTransparency,confidenceLevel):
        wx.Frame.__init__(self, parent, -1, 'Preferences', size=(600,400))

        panelcolour = wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNFACE) 
        self.SetBackgroundColour(panelcolour) 

        mbox = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(mbox)
        
        splitter = wx.SplitterWindow(self, -1, style=wx.SP_LIVE_UPDATE|wx.SP_NOBORDER)
        mbox.Add(splitter, 1, wx.EXPAND)

        accept = wx.Button(self, label = 'OK')
        mbox.Add(accept, 0, flag=wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, border=5)
        accept.Bind(wx.EVT_BUTTON, self.OnReloadParam)

        panel1 = wx.Panel(splitter, -1)
        vbox1 = wx.BoxSizer(wx.VERTICAL)
        panel1.SetSizer(vbox1)

        panel11 = wx.Panel(panel1, -1, size=(-1, 40), style=wx.NO_BORDER)
        panel11.SetBackgroundColour('#53728c')
        st1 = wx.StaticText(panel11, -1, 'Bloks', (5, 5))
        st1.SetForegroundColour('WHITE')

        panel12 = wx.Panel(panel1, -1, style=wx.BORDER_SUNKEN)
        panel12.SetBackgroundColour('WHITE')
        vbox11 = wx.BoxSizer(wx.VERTICAL)
        panel12.SetSizer(vbox11)
        list1 = ListCtrlLeft(panel12, -1)
        vbox11.Add(list1, 1, wx.EXPAND)

        vbox1.Add(panel11, 0, wx.EXPAND)
        vbox1.Add(panel12, 1, wx.EXPAND)
        
        panel2 = wx.Panel(splitter, -1)
        vbox2 = wx.BoxSizer(wx.VERTICAL)
        panel2.SetSizer(vbox2)
        
        panel21 = wx.Panel(panel2, -1, size=(-1, 40), style=wx.NO_BORDER)
        panel21.SetBackgroundColour('#53728c')
        st2 = wx.StaticText(panel21, -1, 'Parameters', (5, 5))
        st2.SetForegroundColour('WHITE')

        panel22 = wx.Panel(panel2, -1, style=wx.BORDER_SUNKEN)
        vbox22 = wx.BoxSizer(wx.VERTICAL)
        panel22.SetSizer(vbox22)
        self.list2 = ListCtrlRight(panel22, -1,
                                   coordFactor, distFactor,
                                   calcMode,maxNumOfIter,minChiSquare,lambdaFirst,lambdaStep,
                                   spheresColor,spheresTransparency,ellipsoidColor,ellipsoidTransparency,confidenceLevel)
        self.list2.SetName('ListControlOnRight')
        vbox22.Add(self.list2, 1, wx.EXPAND)

        vbox2.Add(panel21, 0, wx.EXPAND)
        vbox2.Add(panel22, 1, wx.EXPAND)

        splitter.SplitVertically(panel1, panel2,150)
        self.Centre()
        self.Show(True)
    
    #-----------------------------------------------------------------------------------------------
    def OnReloadParam(self, event):
    
        calcMode = self.list2.param['calc'][0][1]
        maxNumOfIter = self.list2.param['calc'][1][1]
        minChiSquare = self.list2.param['calc'][2][1]
        lambdaFirst = self.list2.param['calc'][3][1]
        lambdaStep = self.list2.param['calc'][4][1]
        coordFactor = self.list2.param['input'][0][1]
        distFactor = self.list2.param['input'][1][1]
        spheresColor = []
        line = self.list2.param['output'][0][1].split()
        for i in range(3):
            spheresColor += [float(line[i]) / 255.0]
        spheresTransparency = self.list2.param['output'][1][1]
        ellipsoidColor = []
        line = self.list2.param['output'][2][1].split()
        for i in range(3):
            ellipsoidColor += [float(line[i]) / 255.0]
        ellipsoidTransparency = self.list2.param['output'][3][1]
        confidenceLevel = self.list2.param['output'][4][1]
        msg = [calcMode,
               maxNumOfIter,
               minChiSquare,
               lambdaFirst,
               lambdaStep,
               coordFactor,
               distFactor,
               spheresColor,
               spheresTransparency,
               ellipsoidColor,
               ellipsoidTransparency,
               confidenceLevel]
        Publisher().sendMessage(("parameters"), msg)
        self.Close()

####################################################################################################
class PlotSpheres:
    '''Plots spheres'''
    
    def __init__(self,num,centers,radii,color,transparency):
        ''''''
        obj = []
        cmd.set("cgo_sphere_quality", 4)
        for i in range(num):
            r,g,b = color
            obj.extend( [ ALPHA, 1-transparency] )
            obj.extend( [ COLOR, r, g, b  ] )
            obj.extend( [ SPHERE, float(centers[i][0]),float(centers[i][1]),float(centers[i][2]), float(radii[i]) ] )
        cmd.load_cgo(obj,'trilateration')

####################################################################################################         
class PlotEllipsoid:
    '''Plots ellipsoid'''
    
    #-----------------------------------------------------------------------------------------------
    def __init__(self, x, y, z, rx, ry, rz, color, transparency):
        ''''''
        cmd.load_cgo(self.makeEllipsoid(x, y, z, rx, ry, rz, color, transparency), 'target')
    
    #-----------------------------------------------------------------------------------------------
    def makeEllipsoid(self, x, y, z, a1, a2, a3, color, transparency):
        ''''''
        return self.makeSuperQuadricEllipsoid(x, y, z, a1, a2, a3, 1.0, 1.0, -pi / 2, pi / 2, -pi, pi, 50, 50, color, transparency)   
    
    #-----------------------------------------------------------------------------------------------
    def makeSuperQuadricEllipsoid(self, x, y, z, a1, a2, a3, n, event, u1, u2, v1, v2, u_segs, v_segs, color, transparency):
        ''''''
        r, g, b = color
        # Calculate delta variables
        dU = (u2 - u1) / u_segs
        dV = (v2 - v1) / v_segs
        obj = [ BEGIN, TRIANGLES ]
        U = u1
        for Y in range(0, u_segs):
                # Initialize variables for loop
                V = v1
                for X in range(0, v_segs):
                        # VERTEX #1 */
                        x1, y1, z1, n1x, n1y, n1z = self.sqEllipsoid(x, y, z, a1, a2, a3, U, V, n, event)
                        x2, y2, z2, n2x, n2y, n2z = self.sqEllipsoid(x, y, z, a1, a2, a3, U + dU, V, n, event)
                        x3, y3, z3, n3x, n3y, n3z = self.sqEllipsoid(x, y, z, a1, a2, a3, U + dU, V + dV, n, event)
                        x4, y4, z4, n4x, n4y, n4z = self.sqEllipsoid(x, y, z, a1, a2, a3, U, V + dV, n, event)
                        obj.extend([COLOR, r, g, b, ALPHA, 1-transparency, NORMAL, n1x, n1y, n1z, VERTEX, x1, y1, z1])
                        obj.extend([COLOR, r, g, b, ALPHA, 1-transparency, NORMAL, n2x, n2y, n2z, VERTEX, x2, y2, z2])
                        obj.extend([COLOR, r, g, b, ALPHA, 1-transparency, NORMAL, n4x, n4y, n4z, VERTEX, x4, y4, z4])
                        obj.extend([COLOR, r, g, b, ALPHA, 1-transparency, NORMAL, n2x, n2y, n2z, VERTEX, x2, y2, z2])
                        obj.extend([COLOR, r, g, b, ALPHA, 1-transparency, NORMAL, n3x, n3y, n3z, VERTEX, x3, y3, z3])
                        obj.extend([COLOR, r, g, b, ALPHA, 1-transparency, NORMAL, n4x, n4y, n4z, VERTEX, x4, y4, z4])
                        # Update variables for next loop
                        V += dV
                # Update variables for next loop
                U += dU
        obj.append(END)
        return obj
    
    #-----------------------------------------------------------------------------------------------
    def sqEllipsoid(self, x, y, z, a1, a2, a3, u, v, n, event):
        ''''''
        x = a1 * self.sqC(u, n) * self.sqC(v, event) + x
        y = a2 * self.sqC(u, n) * self.sqS(v, event) + y
        z = a3 * self.sqS(u, n) + z
        nx = self.sqC(u, 2 - n) * self.sqC(v, 2 - event) / a1
        ny = self.sqC(u, 2 - n) * self.sqS(v, 2 - event) / a2
        nz = self.sqS(u, 2 - n) / a3
        return x, y, z, nx, ny, nz
    
    #-----------------------------------------------------------------------------------------------
    def signOfFloat(self,f):
        ''''''
        if f < 0: return -1
        if f > 0: return 1
        return 0
 
    #-----------------------------------------------------------------------------------------------
    def sqC(self, v, n):
        ''''''
        return self.signOfFloat(cos(v)) *  pow(fabs(cos(v)), n)
 
    #-----------------------------------------------------------------------------------------------
    def sqS(self, v, n):
        ''''''
        return self.signOfFloat(sin(v)) * pow(fabs(sin(v)), n)

####################################################################################################

# SUPLEMENTARY FUNCTIONS
#---------------------------------------------------------------------------------------------------
def getPymolObjects(selection):
	return cmd.get_object_list('(%s)' %selection)
    
#---------------------------------------------------------------------------------------------------
def AverageCoordCalc(coord):
    '''Calculates average coordinates and their standard deviations for a set of coordinates'''
    # Average coordinates
    coordMean = numpy.mean(coord, axis=0)
    # Standard deviations
    coordStd = numpy.std(coord, axis=0)
    
    return coordMean, coordStd         

#---------------------------------------------------------------------------------------------------
def AverageDistCalc(dist,density):
    '''Calculates an average distance and its standard deviations for a distance distribution'''
    # Average distance
    
    distAve = dist[numpy.argmax(density)]
    # Standard deviation
    numerator = 0
    demominator = 0
    for i in range(len(dist)):
        numerator += density[i] * (dist[i]-distAve)**2
        demominator += density[i]
    distStd = sqrt( numerator / demominator )
    
    return distAve, distStd  

#---------------------------------------------------------------------------------------------------
def SquaredRadiusVector(vector1, vector2):
    ''''''
    srv = 0
    for i in range(3):
        srv += (vector1[i] - vector2[i])**2
    
    return srv 

#---------------------------------------------------------------------------------------------------
def ChiSquareCalc(numOfLabels, labelCoord, targetCoordMean, distMean, distStd):
    ''''''
    chiSquare = 0
    for i in range(numOfLabels):
        f = sqrt( SquaredRadiusVector(targetCoordMean[0], labelCoord[i]) ) - distMean[i]
        chiSquare += (1/distStd[i]**2) * f**2
    
    return chiSquare    

#---------------------------------------------------------------------------------------------------
def SingularValueDecomposition(numOfLabels, labelCoord, distMean, distStd):
    '''Linear Least Squares: Singular Value Decomposition algorithm'''
    # Compute A and b coefficients
    A = zeros( (numOfLabels-1,3) )
    b = zeros( (numOfLabels-1,1) )
    for i in range(1, numOfLabels):
        for j in range(3): 
            A[i-1][j] = (labelCoord[i][j]-labelCoord[0][j]) / distStd[i]
        b[i-1] = 0.5 * ( distMean[0]**2 - distMean[i]**2 + SquaredRadiusVector(labelCoord[i], labelCoord[0]) ) / distStd[i]
    # SVD decomposition
    U, S, Vh = linalg.svd(A, full_matrices=0)
    D = linalg.inv(diag(S))
    # Compute the target coordinates and their standard deviations
    targetCoordMean = zeros( (1,3) )
    targetCoordMean = dot(dot(dot(Vh.T, D), U.T), b).T
    for j in range(3): 
        targetCoordMean[0][j] += labelCoord[0][j]
    targetCoordStd = zeros( (1,3) )
    for j in range(3): 
        for i in range(size(A, axis=1)):
            targetCoordStd[0][j] += (Vh.T[j][i] / S[i])**2
        targetCoordStd[0][j] = sqrt(targetCoordStd[0][j])  
    # Compute the Chi-Square value
    chiSquare = ChiSquareCalc(numOfLabels, labelCoord, targetCoordMean, distMean, distStd)
    
    return targetCoordMean, targetCoordStd, chiSquare   

#---------------------------------------------------------------------------------------------------
def InverseHessian(numOfLabels, labelCoord, distMean, distStd, targetCoordInitial, chiSquareMin, iterMax):  
    '''Nonlinear Least Squares: Inverse-Hessian algorithm'''
    # Use previously determined values of target coordinates as initial data 
    targetCoordMean = targetCoordInitial
    # Calculate current Chi-Square value
    chiSquare = ChiSquareCalc(numOfLabels, labelCoord, targetCoordMean, distMean, distStd)
    # Start iterations
    iter = 0
    while (chiSquare > chiSquareMin) and (iter < iterMax):
        # Number of iteration
        iter += 1
        # Save previous Chi-Square value
        chiSquare_temp = chiSquare
        # Calculate the first order and the second order (main part of it) derivatives
        d1f = zeros( (3,1) )
        d2f = zeros( (3,3) )
        for i in range(numOfLabels):
            f = sqrt( SquaredRadiusVector(targetCoordMean[0], labelCoord[i]) ) - distMean[i]
            for j in range(3):
                d1f[j] += (2/distStd[i]**2) * (targetCoordMean[0][j]-labelCoord[i][j]) * f / (f+distMean[i])
                for k in range(3):
                    d2f[j][k] += (2/distStd[i]**2) * (targetCoordMean[0][j]-labelCoord[i][j]) * (targetCoordMean[0][k]-labelCoord[i][k]) / (f+distMean[i])**2
        # Calculate the target coordinates
        targetCoordMean -= dot(inv(d2f), d1f).T
        # Calculate new Chi-Square value
        chiSquare = ChiSquareCalc(numOfLabels, labelCoord, targetCoordMean, distMean, distStd)
        # Check in changes in Chi-Square value
        if (chiSquare < chiSquare_temp) and ((chiSquare_temp - chiSquare)/chiSquare_temp < 0.001): # HIDDEN PARAMETER!!!
            break
    # Calculate the standard deviation of target coordinates
    targetCoordStd = zeros( (1,3) )
    d2f = zeros( (3,3) )
    for i in range(numOfLabels):
            f = sqrt( SquaredRadiusVector(targetCoordMean[0], labelCoord[i]) ) - distMean[i]
            for j in range(3):
                for k in range(3):
                    d2f[j][k] += (2/distStd[i]**2) * (targetCoordMean[0][j]-labelCoord[i][j]) * (targetCoordMean[0][k]-labelCoord[i][k]) / (f+distMean[i])**2          
    correlationMatrix = inv(d2f)
    for j in range(3):
        targetCoordStd[0][j] = sqrt( abs(correlationMatrix[j][j]) ) 
        
    return targetCoordMean, targetCoordStd, chiSquare, iter  

#---------------------------------------------------------------------------------------------------
def LevenbergMarquardt(numOfLabels, labelCoord, distMean, distStd, targetCoordInitial, chiSquareMin, iterMax, lambdaFirst, lambdaStep):
    '''Nonlinear Least Squares: Levenberg-Marquardt algorithm'''
    # Use previously determined values of target coordinates as initial data 
    targetCoordMean = targetCoordInitial
    # Calculate current Chi-Square value
    chiSquare = ChiSquareCalc(numOfLabels, labelCoord, targetCoordMean, distMean, distStd)
    # Start iterations
    lambdaCurrent = lambdaFirst
    iter = 0
    while (chiSquare > chiSquareMin) and (iter < iterMax):
        # Number of iteration
        iter += 1
        # Save previous Chi-Square value and target coordinates
        chiSquare_temp = chiSquare
        targetCoord_temp = targetCoordMean
        # Calculate the first order and the second order (main part of it) derivatives
        d1f = zeros( (3,1) )
        d2f = zeros( (3,3) )
        for i in range(numOfLabels):
            f = sqrt( SquaredRadiusVector(targetCoordMean[0], labelCoord[i]) ) - distMean[i]
            for j in range(3):
                d1f[j] += (2/distStd[i]**2) * (targetCoordMean[0][j]-labelCoord[i][j]) * f / (f+distMean[i])
                for k in range(3):
                    d2f[j][k] += (2/distStd[i]**2) * (targetCoordMean[0][j]-labelCoord[i][j]) * (targetCoordMean[0][k]-labelCoord[i][k]) / (f+distMean[i])**2
        # Modify the second order derivatives matrix
        for j in range(3):
            d2f[j][j] *= (1 + lambdaCurrent)
        # Calculate the target coordinates
        targetCoordMean -= dot(inv(d2f), d1f).T
        # Calculate new Chi-Square value
        chiSquare = ChiSquareCalc(numOfLabels, labelCoord, targetCoordMean, distMean, distStd)
        # Check in changes in Chi-Square value
        if (chiSquare < chiSquare_temp):
            if ((chiSquare_temp - chiSquare)/chiSquare_temp < 0.001): # HIDDEN PARAMETER!!!
                break
            lambdaCurrent /= lambdaStep
        else:
            lambdaCurrent *= lambdaStep
            targetCoordMean = targetCoord_temp
    # Calculate the standard deviation of target coordinates
    targetCoordStd = zeros( (1,3) )
    d2f = zeros( (3,3) )
    for i in range(numOfLabels):
            f = sqrt( SquaredRadiusVector(targetCoordMean[0], labelCoord[i]) ) - distMean[i]
            for j in range(3):
                for k in range(3):
                    d2f[j][k] += (2/distStd[i]**2) * (targetCoordMean[0][j]-labelCoord[i][j]) * (targetCoordMean[0][k]-labelCoord[i][k]) / (f+distMean[i])**2
    correlationMatrix = inv(d2f)
    for j in range(3):
        targetCoordStd[0][j] = sqrt( abs(correlationMatrix[j][j]) ) 
        
    return targetCoordMean, targetCoordStd, chiSquare, iter

####################################################################################################
class MainApp(wx.App):

    def OnInit(self):
        frame = trilat(None, -1, "Trilateration")
        frame.Show(True)
        self.SetTopWindow(frame)
        return True
 

def run():
    app = MainApp(0)
    app.MainLoop()
    

       
def start():
    t = threading.Thread(target=run,args=())
    t.setDaemon(1),
    t.start() 

    
def __init__(self):
    self.menuBar.addmenuitem('Plugin','command','trilater',
                             label = 'mtsslTrilaterate',
                             command = lambda s=self : start() )