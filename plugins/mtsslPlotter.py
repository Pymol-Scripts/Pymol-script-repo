# Import libraries
#----------------------------------------------
import os
import sys
import threading
import warnings
warnings.filterwarnings('ignore')
import string
import math
import numpy
import pymol
from pymol import stored
import matplotlib
matplotlib.use('WXAgg', warn=False, force=False)
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
from matplotlib.ticker import MultipleLocator
from matplotlib import colors
import wx

# Main plotter class
#----------------------------------------------
class Plotter(wx.Frame):
	
	def __init__(self, parent, ID, title):
		wx.Frame.__init__(self, parent, ID, title, wx.DefaultPosition, wx.Size(600,600))
		self.initGui()
		self.initMenu()
		self.Centre()
		self.Show()
		self.graphs = []
		self.currentGraph = {}		  
	def initGui(self):
		# Create main panel
		panel = wx.Panel(self)
		panel.SetFont( wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL) )
		# Set sizers
		sizer0 = wx.BoxSizer(wx.VERTICAL)
		panel.SetSizer(sizer0)
		# Set gui elements
		sizer1 = wx.BoxSizer(wx.HORIZONTAL)
		loadData = wx.Button(panel, label='Load Graphs')
		loadData.Bind(wx.EVT_BUTTON, self.OnLoad)
		sizer1.Add(loadData, flag=wx.ALIGN_CENTER_VERTICAL)
		sizer1.Add((20, -1))
		sizer1.Add(wx.StaticText(panel, label='Current Graph'), flag=wx.ALIGN_CENTER_VERTICAL)
		sizer1.Add((4, -1))
		self.graphsList = wx.ComboBox(panel, choices=[], size = [300, 24], style=wx.CB_READONLY)
		self.graphsList.Bind(wx.EVT_COMBOBOX, self.OnSelectGraph)
		sizer1.Add(self.graphsList, proportion=1)
		sizer0.Add(sizer1, flag=wx.EXPAND|wx.ALIGN_CENTER_HORIZONTAL|wx.TOP|wx.LEFT|wx.RIGHT|wx.BOTTOM, border=20)
		# Matplotlib item
		panelcolour = wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNFACE).Get()
		facecolor = [channel/float(255) for channel in panelcolour]
		self.figure = Figure(figsize=(-1, -1), dpi=80, facecolor=facecolor, linewidth=1.0)
		self.figure.subplots_adjust(left=0.08, bottom=0.08, right=0.98, top=0.98)
		self.canvas = FigCanvas(panel, -1, self.figure)
		self.axes = self.figure.add_subplot(111) 
		self.axes.set_xlabel('X-Axis')
		self.axes.set_ylabel('Y-Axis')
		self.axes.set_xlim(0.0, 1.0)
		self.axes.set_ylim(0.0, 1.0)
		matplotlib.rcParams.update({'font.size': 12})
		self.axes.grid(False)
		sizer0.Add(self.canvas, 1, flag=wx.GROW|wx.LEFT|wx.RIGHT, border=20)
		self.toolbar = NavigationToolbar(self.canvas)
		sizer0.Add(self.toolbar, 0, flag=wx.ALIGN_CENTER_HORIZONTAL|wx.LEFT|wx.RIGHT, border=20)
		sizer0.Add((-1, 20))
		self.statusBar = wx.StatusBar(self, -1)
		self.statusBar.SetFieldsCount(1)
		self.SetStatusBar(self.statusBar)
		self.canvas.mpl_connect('motion_notify_event', self.UpdateStatusBar)
		self.canvas.Bind(wx.EVT_ENTER_WINDOW, self.ChangeCursor)
		# Redefine close event
		self.Bind(wx.EVT_CLOSE, self.OnClose)
	def initMenu(self):
		menubar = wx.MenuBar()
		# File menu
		fileMenu = wx.Menu()
		openItem = wx.MenuItem(fileMenu, 11, '&Open\tCtrl+O')
		self.Bind(wx.EVT_MENU, self.OnOpen, id=11)
		saveItem = wx.MenuItem(fileMenu, 12, '&Save\tCtrl+S')
		self.Bind(wx.EVT_MENU, self.OnSave, id=12)
		exitItem = wx.MenuItem(fileMenu, 13, '&Exit\tAlt+F4')
		self.Bind(wx.EVT_MENU, self.OnClose, id=13)
		fileMenu.AppendItem(openItem)
		fileMenu.AppendItem(saveItem)
		fileMenu.AppendSeparator()
		fileMenu.AppendItem(exitItem)
		# Edit menu
		editMenu = wx.Menu()
		loadItem = wx.MenuItem(editMenu, 21, '&Load data')
		self.Bind(wx.EVT_MENU, self.OnLoad, id=21)
		clearItem = wx.MenuItem(editMenu, 22, '&Clear memory')
		self.Bind(wx.EVT_MENU, self.OnClear, id=22)
		editMenu.AppendItem(loadItem)
		editMenu.AppendItem(clearItem)
		# Help menu
		helpMenu = wx.Menu()
		aboutItem = wx.MenuItem(helpMenu, 31, '&About\tCtrl+A')
		#self.Bind(wx.EVT_MENU, self.OnAbout, id=31)
		helpItem = wx.MenuItem(helpMenu, 32, '&Help\tF1')
		#self.Bind(wx.EVT_MENU, self.OnHelp, id=32)
		helpMenu.AppendItem(aboutItem)
		helpMenu.AppendSeparator()
		helpMenu.AppendItem(helpItem)
		menubar.Append(fileMenu, '&File')
		menubar.Append(editMenu, '&Edit')
		menubar.Append(helpMenu, '&Help')
		self.SetMenuBar(menubar)
	def ChangeCursor(self,event):
		self.canvas.SetCursor(wx.StockCursor(wx.CURSOR_BULLSEYE))
	def UpdateStatusBar(self,event):
		if event.inaxes:
			x, y = event.xdata, event.ydata
			self.statusBar.SetStatusText(( "x = " + str(numpy.rint(x)) + "; y = " +str(numpy.rint(y)) ), 0) 
	def OnLoad(self,event):
		# Load a dictionary with plots
		self.graphs = []
		self.graphs.extend(stored.plots)
		# Update a list of graphs
		#self.graphsList.Enable()
		self.graphsList.Clear()
		for graph in self.graphs:
			self.graphsList.Append(graph['Title'])
		event.Skip()
	def OnClear(self,event):
		# Clear a dictionary with plots
		self.graphs = []
		self.graphsList.Clear()
		event.Skip()
	def OnOpen(self,event):
		dlg = wx.FileDialog(self, message="Choose a file", defaultDir=os.getcwd(), defaultFile="", wildcard="*.dat", style=wx.OPEN | wx.CHANGE_DIR)
		if dlg.ShowModal() == wx.ID_OK:
			path = dlg.GetPath()
			dataFile = open(path, 'r')
			# Read data file
			# ...
		dlg.Destroy()
	def OnSave(self,event):
		if (bool(self.graphs) == True) and (bool(self.currentGraph) == True):
			dialog = wx.FileDialog(None, message="Save data as ...", defaultDir=os.getcwd(), defaultFile="", wildcard="*.dat", style = wx.SAVE | wx.OVERWRITE_PROMPT)
			if dialog.ShowModal() == wx.ID_OK:
				path = dialog.GetPath()
				# Create a file
				dataFile = open(path, 'w')
				# Fill in the content of a data file
				title  = 'Title: '	+ self.currentGraph['Title']  + '\n'
				type   = 'Type: '	+ self.currentGraph['Type']	  + '\n'
				xTitle = 'xTitle: ' + self.currentGraph['xTitle'] + '\n'
				yTitle = 'yTitle: ' + self.currentGraph['yTitle'] + '\n'
				axesTitles = '{0:>10s}'.format('x') + '{0:>10s}'.format('y') + '{0:>10s}'.format('z') + '\n'
				dataFile.write(title)
				dataFile.write(type)
				dataFile.write(xTitle)
				dataFile.write(yTitle)
				dataFile.write(axesTitles)
				# Write x,y,z to a data file
				graphType = self.currentGraph['Type']
				print self.currentGraph['xData']
				#print self.currentGraph['xData'].shape[0], self.currentGraph['yData'].shape[0], self.currentGraph['zData'].shape
				dim = self.currentGraph['xData'].shape[0]
				if (graphType == 'DistanceMap'):
					for i in range(dim):
						x = self.currentGraph['xData'][i]
						dataFile.write('{0:10.2f}'.format(x))
						y = self.currentGraph['yData'][i]
						dataFile.write('{0:10.2f}'.format(y))
						for j in range(dim):
							z = self.currentGraph['zData'][i][j]
							dataFile.write('{0:10.2f}'.format(z))
						dataFile.write('\n')  
				if (graphType == 'DifferenceMap'):
					pass
				if (graphType == 'DistancePlot'):
					pass
				if (graphType == 'AccessibilityPlot'):
					pass
				dataFile.close()
			dialog.Destroy()
	def OnClose(self,event):
		stored.mtsslplot = 0
		self.Destroy()
	def OnSelectGraph(self,event):
		# Read on the selected graph
		title = self.graphsList.GetStringSelection()
		for graph in self.graphs:
			if (graph['Title'] == title):
				self.currentGraph = graph
		# Plot the selected graph
		self.plotGraph()
		event.Skip()
	def plotGraph(self):
		#print "Selecting graph"
		self.refreshFigure()
		# Read on the type of a graph
		graphType = self.currentGraph['Type']
		if (graphType == 'DistanceMap'):
			self.plotDistanceMap()
		if (graphType == 'DifferenceMap'):
			self.plotDifferenceMap()
		if (graphType == 'DistancePlot'):
			self.plotDistancePlot()
		if (graphType == 'AccessibilityPlot'):
			self.plotAccessibilityPlot()
		if (graphType == 'DistanceDistribution'):
			#print "plotting distribution"
			self.plotDistanceDistribution()
	def refreshFigure(self):
		self.figure.clear()
		self.axes = self.figure.add_subplot(111)
		self.axes.set_xlabel('X-Axis')
		self.axes.set_ylabel('Y-Axis')
		self.axes.set_xlim(0.0, 1.0)
		self.axes.set_ylim(0.0, 1.0)
		matplotlib.rcParams.update({'font.size': 12})
	def plotDistanceMap(self):
		# Read on x,y,z
		x = self.currentGraph['xData'] - 0.5 * numpy.ones(len(self.currentGraph['xData']))
		y = self.currentGraph['yData'] - 0.5 * numpy.ones(len(self.currentGraph['yData']))
		X, Y = numpy.meshgrid(x, y)
		Z = self.currentGraph['zData']
		# Define colormap
		cmap = colors.ListedColormap(['blue', 'green', 'orange', 'red'])
		cmap.set_under('white')
		cmap.set_over('white')
		bounds = [1,15,50,80,100]
		norm = colors.BoundaryNorm(bounds, cmap.N)
		# Draw surface plot
		img = self.axes.pcolor(X, Y, Z, cmap=cmap, norm=norm)
		self.axes.set_xlim(x.min(), x.max())
		self.axes.set_ylim(y.min(), y.max())
		self.axes.set_xlabel(self.currentGraph['xTitle'])
		self.axes.set_ylabel(self.currentGraph['yTitle'])
		# Cosmetics
		#matplotlib.rcParams.update({'font.size': 12})
		xminorLocator = MultipleLocator(10)
		yminorLocator = MultipleLocator(10)
		self.axes.xaxis.set_minor_locator(xminorLocator)
		self.axes.yaxis.set_minor_locator(yminorLocator)
		self.axes.tick_params(direction='out', length=6, width=1)
		self.axes.tick_params(which='minor', direction='out', length=3, width=1)
		self.axes.xaxis.labelpad = 15
		self.axes.yaxis.labelpad = 15
		# Draw colorbar
		colorbar = self.figure.colorbar(img, boundaries = [0,1,15,50,80,100], 
										spacing = 'proportional',
										ticks = [15,50,80,100], 
										extend = 'both')
		colorbar.ax.set_xlabel('Angstrom')
		colorbar.ax.xaxis.set_label_position('top')
		colorbar.ax.xaxis.labelpad = 20
		self.figure.tight_layout()		
		self.canvas.draw()
	def plotDifferenceMap(self):
		# Read on x,y,z
		x = self.currentGraph['xData'] - 0.5 * numpy.ones(len(self.currentGraph['xData']))
		y = self.currentGraph['yData'] - 0.5 * numpy.ones(len(self.currentGraph['yData']))
		X, Y = numpy.meshgrid(x, y)
		Z = self.currentGraph['zData']	   
		# Draw surface plot
		img = self.axes.pcolor(X, Y, Z, cmap='jet', vmin=0.0, vmax=numpy.amax(Z))
		self.axes.set_xlim(x.min(), x.max())
		self.axes.set_ylim(y.min(), y.max())
		self.axes.set_xlabel(self.currentGraph['xTitle'])
		self.axes.set_ylabel(self.currentGraph['yTitle'])
		# Cosmetics
		#matplotlib.rcParams.update({'font.size': 18})
		xminorLocator = MultipleLocator(10)
		yminorLocator = MultipleLocator(10)
		self.axes.xaxis.set_minor_locator(xminorLocator)
		self.axes.yaxis.set_minor_locator(yminorLocator)
		self.axes.tick_params(direction='out', length=6, width=1)
		self.axes.tick_params(which='minor', direction='out', length=3, width=1)
		self.axes.xaxis.labelpad = 15
		self.axes.yaxis.labelpad = 15
		# Draw colorbar
		colorbar = self.figure.colorbar(img)
		colorbar.ax.set_xlabel('Angstrom')
		colorbar.ax.xaxis.set_label_position('top')
		colorbar.ax.xaxis.labelpad = 20
		self.figure.tight_layout()		
		self.canvas.draw()
	def plotDistancePlot(self):
		# Read on x,y
		x = self.currentGraph['xData']
		y = self.currentGraph['yData']
		# Remove data points with distances <= 0
		index=numpy.zeros(1)
		for i in range (0,y.shape[0]):
			if (y[i] <= 0):
				index = numpy.append(index,i)
				index = numpy.delete(index, 0)
		Y = numpy.delete(y, index)
		X = numpy.delete(x, index)		
		# Plot graph
		self.axes.plot(X, Y, 
					   linestyle='solid',
					   linewidth=0.5, 
					   color='black', 
					   marker='o',
					   markerfacecolor='white', 
					   markeredgecolor='black', 
					   markersize=6, 
					   markeredgewidth=1)
		Xmin = 0
		Xmax = X.max()+1
		Ymin = 0
		Ymax = Y.max() * 1.05
		self.axes.set_xlim(Xmin, Xmax)
		self.axes.set_ylim(Ymin, Ymax)
		self.axes.set_xlabel(self.currentGraph['xTitle'])
		self.axes.set_ylabel(self.currentGraph['yTitle'])
		# Color the distance ranges
		self.axes.fill_between([Xmin, Xmax], 15, Ymin, color='blue', alpha=0.4)
		self.axes.fill_between([Xmin, Xmax], 50, 15, color='green', alpha=0.4)
		self.axes.fill_between([Xmin, Xmax], 80, 50, color='orange', alpha=0.4)
		self.axes.fill_between([Xmin, Xmax], Ymax, 80, color='red', alpha=0.4)
		# Cosmetics
		#matplotlib.rcParams.update({'font.size': 18})
		xminorLocator = MultipleLocator(10)
		yminorLocator = MultipleLocator(10)
		self.axes.xaxis.set_minor_locator(xminorLocator)
		self.axes.yaxis.set_minor_locator(yminorLocator)
		self.axes.tick_params(direction='out', length=6, width=1)
		self.axes.tick_params(which='minor', direction='out', length=3, width=1)
		self.axes.xaxis.labelpad = 15
		self.axes.yaxis.labelpad = 15
		self.figure.tight_layout()		
		self.canvas.draw()
	def plotDistanceDistribution(self):
		# Read on x,y
		x = self.currentGraph['xData']
		y = self.currentGraph['yData']
		# Remove data points with distances <= 0
		index=numpy.zeros(1)
		for i in range (0,y.shape[0]):
			if (y[i] <= 0):
				index = numpy.append(index,i)
				index = numpy.delete(index, 0)
		Y = numpy.delete(y, index)
		X = numpy.delete(x, index)		
		# Plot graph
		self.axes.plot(X, Y, 
					   linestyle='solid',
					   linewidth=0.5, 
					   color='red', 
					   #marker='.',
					   markerfacecolor='white', 
					   markeredgecolor='black', 
					   markersize=6, 
					   markeredgewidth=1)
		Xmin = 0
		Xmax = X.max()+1
		Ymin = 0
		Ymax = Y.max() * 1.05
		self.axes.set_xlim(Xmin, Xmax)
		self.axes.set_ylim(Ymin, Ymax)
		self.axes.set_xlabel(self.currentGraph['xTitle'])
		self.axes.set_ylabel(self.currentGraph['yTitle'])
		# Color the distance ranges
		#self.axes.fill_between([Xmin, Xmax], 15, Ymin, color='blue', alpha=0.4)
		#self.axes.fill_between([Xmin, Xmax], 50, 15, color='green', alpha=0.4)
		#self.axes.fill_between([Xmin, Xmax], 80, 50, color='orange', alpha=0.4)
		#self.axes.fill_between([Xmin, Xmax], Ymax, 80, color='red', alpha=0.4)
		# Cosmetics
		#matplotlib.rcParams.update({'font.size': 18})
		xminorLocator = MultipleLocator(10)
		yminorLocator = MultipleLocator(10)
		self.axes.xaxis.set_minor_locator(xminorLocator)
		self.axes.yaxis.set_minor_locator(yminorLocator)
		self.axes.tick_params(direction='out', length=6, width=1)
		self.axes.tick_params(which='minor', direction='out', length=3, width=1)
		self.axes.xaxis.labelpad = 15
		self.axes.yaxis.labelpad = 15
		self.figure.tight_layout()		
		self.canvas.draw()
	def plotAccessibilityPlot(self):
		# Read on x,y
		x = self.currentGraph['xData']
		y = self.currentGraph['yData']
		X = x
		Y = y
		# Plot graph
		self.axes.plot(X, Y, 
					   linestyle='solid',
					   linewidth=0.5, 
					   color='red', 
					   marker='.',
					   markerfacecolor='white', 
					   markeredgecolor='black', 
					   markersize=6, 
					   markeredgewidth=1)
		Xmin = 0
		Xmax = X.max()+1
		Ymin = 0
		Ymax = Y.max() * 1.05
		self.axes.set_xlim(Xmin, Xmax)
		self.axes.set_ylim(Ymin, Ymax)
		self.axes.set_xlabel(self.currentGraph['xTitle'])
		self.axes.set_ylabel(self.currentGraph['yTitle'])
		# Cosmetics
		#self.figure.set_size_inches(30,3)
		#matplotlib.rcParams.update({'font.size': 18})
		xminorLocator = MultipleLocator(10)
		yminorLocator = MultipleLocator(10)
		self.axes.xaxis.set_minor_locator(xminorLocator)
		self.axes.yaxis.set_minor_locator(yminorLocator)
		self.axes.tick_params(direction='out', length=6, width=1)
		self.axes.tick_params(which='minor', direction='out', length=3, width=1)
		self.axes.xaxis.labelpad = 15
		self.axes.yaxis.labelpad = 15
		self.figure.tight_layout()		
		self.canvas.draw()
	
# Initialize App
#----------------------------------------------
class MainApp(wx.App):
	def OnInit(self):
		frame = Plotter(None, -1, "mtsslPlotter")
		frame.Show(True)
		self.SetTopWindow(frame)
		return True

def run():
	app = MainApp(0)
	app.MainLoop()

def start():
	if hasattr(stored, 'mtsslplot'):
		if (stored.mtsslplot == 1):
			print "mtsslPlotter is already opened!"
		else:
			t = threading.Thread(target=run,args=())
			t.setDaemon(0)
			t.start()
			stored.mtsslplot = 1
	else:
		t = threading.Thread(target=run,args=())
		t.setDaemon(0)
		t.start()
		stored.mtsslplot = 1

def __init__(self):
	self.menuBar.addmenuitem('Plugin','command','mtsslPlotterV3',
							 label = 'mtsslPlotterV3',
							 command = lambda s=self : start() ) 