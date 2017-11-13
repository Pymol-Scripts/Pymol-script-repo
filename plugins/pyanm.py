'''
   PyANM--Pymol Plugin for easier Anisotropic Network Model (ANM) building and visualizing
   
   Written by:
               Yuan Wang (yuanwang at iastate dot edu)
               Dr. Robert Jernigan Lab
               Iowa State University
               Sept, 2014
'''

# Copyright Notice
# ================
# 
# The PyMOL Plugin source code in this file is copyrighted, but you can
# freely use and copy it as long as you don't change or remove any of
# the copyright notices.
# 
# ----------------------------------------------------------------------
# This PyMOL Plugin is Copyright (C) 2014 by Yuan Wang <yuanwang at iastate dot edu>
# 
#                        All Rights Reserved
# 
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name(s) of the author(s) not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
# 
# THE AUTHOR(S) DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
# INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.  IN
# NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY SPECIAL, INDIRECT OR
# CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
# USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
# OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
# PERFORMANCE OF THIS SOFTWARE.
# ------------------------------------
import warnings
import os 
import sys
if sys.version_info[0] < 3:
    import Tkinter
    import tkSimpleDialog
    import tkMessageBox
    import tkFileDialog
    import tkColorChooser
else:
    import tkinter as Tkinter
    import tkinter.simpledialog as tkSimpleDialog
    import tkinter.messagebox as tkMessageBox
    import tkinter.filedialog as tkFileDialog
    import tkinter.colorchooser as tkColorChooser
import Pmw
from pymol.cgo import *

try:
    import numpy as np
except ImportError:
    raise ImportError('Failed to import numpy, which is required for building ANM')
    
try:
    # In MacPyMOL, even importing linalg causes a crash
    import platform
    if platform.system() == 'Darwin':
        haveSP = False
    else:
        import scipy.sparse.linalg
        haveSP = True
except ImportError:
    haveSP = False

try:
    from pymol import *
    havePymol = True
except ImportError:
    havePymol = False
    warnings.warn("Failed to load module pymol, functions might be incomplete")

#####################################
### Add PyANM in Pymol Plugin menu ##
#####################################
def __init__(self):
    
    self.menuBar.addmenuitem('Plugin','command','PyANM',label = 'PyANM', command = lambda s = self : PyANMPlugin(s))

#####################################
########### GUI Class ###############
#####################################
class PyANMPlugin:
    
    def __init__(self, app):
        '''
        this is the GUI class for the plugin
        calculations are done in class ANM and structure
        '''
        
        self.parent = app.root
        self.dialog = Pmw.Dialog(self.parent,
                                 buttons = ('Build ANM',
                                            'Show Springs',
                                            'Color by MSFs',
                                            'Make Movies',
                                            'Draw Arrows',
                                            'Export Data',
                                            'Exit'),
                                 title = 'PyANM -- ANM Plugin for Pymol',
                                 defaultbutton = 'Build ANM',
                                 command = self.execute)
        Pmw.setbusycursorattributes(self.dialog.component('hull'))
        
        self.structure = Tkinter.StringVar()
        self.cutoffValue = Tkinter.DoubleVar()
        self.useCAonly = Tkinter.BooleanVar()
        self.includeHETATM = Tkinter.BooleanVar()
        self.saveDir = Tkinter.StringVar()
        self.modesForMovies = Tkinter.StringVar()
        self.modesForArrows = Tkinter.StringVar()
        self.powerValue = Tkinter.DoubleVar()
        self.anmMethod = None
        self.anmBuilt = False
        self.arrowRGB = [1, 1, 1]
        self.arrowLoaded = False
        self.movieLoaded = False
        self.arrowSize = 5
        self.movieScales = 5
        self.exportE = Tkinter.BooleanVar()
        self.exportV = Tkinter.BooleanVar()
        self.exportMSF = Tkinter.BooleanVar()
        self.exportCX = Tkinter.BooleanVar()
        self.exportHess = Tkinter.BooleanVar()
        self.exportAll = Tkinter.BooleanVar()
        
        w = Tkinter.Label(self.dialog.interior(),
                          text = '\nPyANM \nYuan Wang \nIowa State University, 2014',
                          background = 'black', foreground = 'green')
        w.pack(expand = 1, fill = 'both', padx = 10, pady = 2)
        
        
        self.notebook = Pmw.NoteBook(self.dialog.interior())
        self.notebook.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        
######################################
########## ANM TAB ###################
######################################
        
        page = self.notebook.add('ANM')
        self.notebook.tab('ANM').focus_set()
        group_struc = Tkinter.LabelFrame(page, text = 'Structure Handling')
        group_struc.pack(fill = 'both', expand = True, padx = 10, pady = 5)
        
       
        
        structure_entry = Pmw.EntryField(group_struc,
                                         label_text = 'Select Structure:',
                                         labelpos = 'wn',
                                         entry_textvariable = self.structure)
                                         
        structure_button = Tkinter.Button(group_struc, text = 'Browse...',
                                          command = self.getStructurePath)
       
        
                                      
        useCAonly_checkbox = Tkinter.Checkbutton(group_struc,
                                    text = "Use only Alpha Carbons?",
                                    variable = self.useCAonly)
                                    
        includeHET_checkbox = Tkinter.Checkbutton(group_struc,
                                                  text = "Include HETATMs?",
                                                  variable = self.includeHETATM)
        useCAonly_checkbox.select()                                        
                        
        structure_entry.grid(sticky = 'w', row = 0, column = 0,
                             columnspan = 2, padx = 5, pady = 5)
        structure_button.grid(sticky = 'w', row = 0, column = 1,
                              padx = 5, pady = 1)
        
        useCAonly_checkbox.grid(sticky = 'w', row = 2, column = 0,
                       padx = 3, pady = 5)
        includeHET_checkbox.grid(sticky = 'w', row = 3, column = 0,
                                 padx= 3, pady = 5)
       
        
        self.group_struc2 = Tkinter.LabelFrame(page, text = 'ANM Parameters')
        self.group_struc2.pack(fill = 'both', expand = True, padx = 10, pady = 5)
        
        self.method_select = Pmw.ScrolledListBox(self.group_struc2,
                                            items = ('Cutoff Model', 'Parameter Free'),
                                            labelpos = 'ne',
                                            listbox_height = 2,
                                            label_text = 'Please select your ANM Method:',
                                            selectioncommand = self.askFurther)
                                            
        
        self.method_select.grid(sticky = 'w', row = 0, column = 0, padx = 3, pady = 1)
        
######################################
######## Control Panel ###############
######################################
        
        page = self.notebook.add('Control Panel')
        
        moviePanel = Tkinter.LabelFrame(page, text = 'Movie Controls')
        moviePanel.grid(sticky = 'we', row = 0, column = 0,
                        padx = 10, pady = 5)
        movieScaleLabel = Tkinter.Label(moviePanel, text = 'Scale of movies:')
        movieScaleLabel.grid(sticky = 'e', row = 0, column = 0, padx = 5, pady = 3)
        self.movieScale = Tkinter.Scale(moviePanel,
                                   from_ = 0.0,
                                   to = 10.0,
                                   resolution = 0.5,
                                   orient = Tkinter.HORIZONTAL,
                                   command = self.changeMoiveScale)
                                   
        self.movieScale.grid(stick = 'we', row = 0, column = 1, padx = 5, pady = 3)
        self.movieScale.set(5.0)
        arrowPanel = Tkinter.LabelFrame(page, text = 'Arrow Controls')
        arrowPanel.grid(sticky = 'we', row = 1, column = 0,
                        padx= 10, pady = 5)
        arrowColorLabel = Tkinter.Label(arrowPanel, text = 'Color of the arrows:')
        self.arrowColorButton = Tkinter.Button(arrowPanel, 
                                          bg = 'white',
                                          command = self.getColorRGB)
        arrowSizeLabel = Tkinter.Label(arrowPanel, text = 'Size of the arrows:')
        self.arrowSizeScale = Tkinter.Scale(arrowPanel,
                                            from_ = 0.0,
                                            to = 10.0,
                                            resolution = 0.5,
                                            orient = Tkinter.HORIZONTAL,
                                            command = self.changeArrowSize)
        self.arrowSizeScale.set(5.0)
        arrowColorLabel.grid(sticky = 'e', row = 0, column = 0, padx = 5, pady = 3)
        self.arrowColorButton.grid(sticky = 'we', row = 0, column = 1, padx = 5, pady = 3 )
        arrowSizeLabel.grid(sticky = 'e', row = 1, column = 0, padx = 5, pady = 3)
        self.arrowSizeScale.grid(sticky = 'we', row = 1, column = 1, padx = 5, pady = 3)      
        
        exportPanel = Tkinter.LabelFrame(page, text = 'Export Selections')
        exportPanel.grid(sticky = 'e', row = 0, column = 1, rowspan = 2,
                         padx = 10, pady = 5)
                         
        exportE_checkbox = Tkinter.Checkbutton(exportPanel,
                                               text = 'Export Eigenvalues',
                                               variable = self.exportE)
        exportE_checkbox.grid(sticky = 'w', row = 0, column = 0)
        exportE_checkbox.select()
        
        exportV_checkbox = Tkinter.Checkbutton(exportPanel,
                                               text = 'Export Eigenvectors',
                                               variable = self.exportV)
        exportV_checkbox.grid(sticky = 'w', row = 1, column = 0)
        exportV_checkbox.select()
        
        exportMSF_checkbox = Tkinter.Checkbutton(exportPanel,
                                                 text = 'Export MSF',
                                                 variable = self.exportMSF)
        exportMSF_checkbox.grid(sticky = 'w', row = 2, column = 0)
        
        exportCX_checkbox = Tkinter.Checkbutton(exportPanel,
                                                 text = 'Export Contact Matrix',
                                                 variable = self.exportCX)
        exportCX_checkbox.grid(sticky = 'w', row = 3, column = 0)
        
        exportHess_checkbox = Tkinter.Checkbutton(exportPanel,
                                                 text = 'Export Hessian Matrix',
                                                 variable = self.exportHess)
        exportHess_checkbox.grid(sticky = 'w', row = 4, column = 0)
        
        exportAll_checkbox = Tkinter.Checkbutton(exportPanel,
                                                 text = 'All of Above',
                                                 variable = self.exportAll,
                                                 command = self.selectAll)
        exportAll_checkbox.grid(sticky = 'w', row = 5, column = 0)
######################################
########## ABOUT TAB #################
######################################
        
        page = self.notebook.add('About')
        group_about = Tkinter.LabelFrame(page, text = 'About PyANM')
        group_about.grid(sticky = 'we', row = 0, column = 0, 
                         padx = 5, pady = 3)
        about_plugin = """PyANM is a Pymol Plugin which allows you to build 
                          Anisotropic Network Models (ANM), 
                          visulize springs and mode motions within Pymol.
                          Author: Yuan Wang
                          Iowa State University
                          Sept, 2014"""
        label_about = Tkinter.Label(group_about, text = about_plugin)
        label_about.grid(sticky = 'w', row = 0, column = 0,
                          padx = 5, pady = 10)
                          
        
        return
        
    def askFurther(self):
        '''
        choose the specific ANM (cutoff or pf)
        '''
        
        if self.method_select.getcurselection()[0] == 'Cutoff Model':
            self.cutoff_entry = Pmw.EntryField(self.group_struc2,
                                      label_text = 'Set cutoff value (A):',
                                      labelpos = 'wn',
                                      value = 12,
                                      entry_textvariable = self.cutoffValue)
            self.cutoff_entry.grid(sticky = 'w', row = 1, column = 0, 
                          padx = 5, pady = 1) 
       
        if self.method_select.getcurselection()[0] == 'Parameter Free':
            self.cutoff_entry = Pmw.EntryField(self.group_struc2,
                                      label_text = 'Set Power:',
                                      labelpos = 'wn',
                                      value = 6,
                                      entry_textvariable = self.powerValue)
            self.cutoff_entry.grid(sticky = 'w', row = 1, column = 0, 
                          padx = 5, pady = 1) 
                          
        self.anmMethod = self.method_select.getcurselection()[0]        
        
    def getStructurePath(self):
        '''
        get the path to pdb file if user chooses to use a local file
        '''
        
        filePath = tkFileDialog.askopenfilename(title = 'File Location',
                                                initialdir = '',
                                                filetypes = [('all', '*')],
                                                parent = self.parent)
        if filePath:
            self.structure.set(filePath)
            cmd.load(filePath)
        return
        
    def getColorRGB(self):
        '''
        get the rgb color for the arrows
        '''
        
        try:
            colorTuple, color = tkColorChooser.askcolor()
            if colorTuple is not None and color is not None:            
                self.arrowColorButton.config(bg = color)               
                self.arrowRGB = list()
                for i in range(3):
                    self.arrowRGB.append((float(colorTuple[i])/255))
        except Tkinter._tkinter.TclError:
            self.arrowRGB = [1,1,1]
            
        if self.arrowLoaded:
            for arrows in self.arrowObj:
                cmd.delete(arrows)
            
            self.arrowDrawer(modes = self.modesForArrows)
        
    def changeArrowSize(self, value):
        '''
        change the size of the arrows
        
        value: return value from the scroll bar for arrow sizes
        '''
        
        self.arrowSize = float(value)
        
        if self.arrowLoaded:
            for arrows in self.arrowObj:
                cmd.delete(arrows)
            
            self.arrowDrawer(modes = self.modesForArrows)
            
    def changeMoiveScale(self, value):
        '''
        change the scale for movies
        
        value: return value from the scroll bar for movie scales
        '''
        
        self.movieScales = float(value) 
        
        if self.movieLoaded:
            for movie in anm.movieFiles:
                cmd.delete(anm.movieFiles[movie].split('.')[0])
            
                
            self.movieMaker(modes = self.modesForMovies)
            
    def selectAll(self):
        '''
        select all button behavior
        '''        
        if self.exportAll.get():
            self.exportE.set(True)
            self.exportV.set(True)
            self.exportMSF.set(True)
            self.exportCX.set(True)
            self.exportHess.set(True)
                                                 
    def ANMbuilder(self):
        '''
        builds ANM using class ANM
        '''        
        
        global anm
        flag = False
        if self.anmMethod is None:
            tkMessageBox.showerror('Unable to Build ANM','Please select your ANM Model First')
        
        if self.anmMethod == 'Cutoff Model':
            flag = True
            anm = ANM(self.structure.get(), cutoff = self.cutoffValue.get(),
                      useCA = self.useCAonly.get())
        elif self.anmMethod == 'Parameter Free':
            flag = True            
            anm = ANM(self.structure.get(), method = 'pf', power = self.powerValue.get(),
                      useCA = self.useCAonly.get())
        if flag: 
            anm.buildCX()
            anm.buildHess()
            anm.calcModes()
            if anm.freqChecker(): 
                self.anmBuilt = True
                print("ANM Built successfully")
                tkMessageBox.showinfo("Built Successfully","ANM was built successfully",
                                      parent = self.parent)
            else:
                self.anmBuilt = True
                print("Not all of first 6 eigenvalues are zero")
                tkMessageBox.showwarning("Warning","Not all of first 6 eigenvalues are zero, calculations can continue but might not be accurate")
        
        
            
    def SpringsBuilder(self):
        '''
        makes the spring file to show springs in ANM
        '''        
        proteinName = os.path.basename(self.structure.get())
        filename = proteinName.split('.')[0] + '_springs.pdb'
        obj_name = proteinName.split('.')[0] + '_springs'
        obj_list = cmd.get_names()
        if obj_name in obj_list:
            cmd.delete(obj_name)
        struc.writePDB(filename)
        springFile = open(filename,'a')
        for i in range (0,len(anm.contactPair_x)):
            springFile.write('CONECT %4d %4d\n' %(struc.atom_number[anm.contactPair_x[i]], struc.atom_number[anm.contactPair_y[i]]))
        springFile.close()
        cmd.hide('everything')
        cmd.load(filename)
        os.remove(filename)
        
    def MSFcolorer(self):
        '''
        color structure based on MSF calculations
        '''
        struc.b_factor = anm.getMSF().tolist()[0]      
        proteinName = os.path.basename(self.structure.get())
        filename = proteinName.split('.')[0] + '_MSFs.pdb'
        obj_name = proteinName.split('.')[0] + '_MSFs'        
        obj_list = cmd.get_names()
        if obj_name in obj_list:
            cmd.delete(obj_name)
        struc.writePDB(filename)
        cmd.hide('everything')        
        cmd.load(filename)
        if struc.isCA():
            cmd.set('ribbon_trace',1)
        cmd.show('ribbon',obj_name)
        cmd.spectrum('b','blue_red', obj_name)
        
    def movieMaker(self, modes = None):
        '''
        builds the movie files for different modes
        '''
        
        if modes is None:    
            self.modesForMovies = tkSimpleDialog.askstring('Input Mode Number',
                                                           'Please enter the mode numbers:',
                                                           initialvalue = '1 2 3',
                                                           parent = self.parent)
            if self.modesForMovies is not None:
                self.modesForMovies = list(map(int, self.modesForMovies.split()))
                modes = self.modesForMovies
            else:
                return
                
        anm.modeAnimator(modes = modes, scaler = self.movieScales)
        cmd.hide('everything')
        if struc.isCA():
            cmd.set('ribbon_trace',1)
            
        for i in range(len(anm.movieFiles)):
            cmd.load(anm.movieFiles[i])
            cmd.show('ribbon',anm.movieFiles[i].split('.')[0])
        
        self.movieLoaded = True
        
    def arrowDrawer(self, modes = None):
        '''
        draws arrows to represent different modes
        '''
        
        if modes is None:
            self.modesForArrows = tkSimpleDialog.askstring('Input Mode Number',
                                                           'Please enter the mode numbers:',
                                                           initialvalue = '1 2 3',
                                                           parent = self.parent)
            if self.modesForArrows is not None:
                self.modesForArrows = list(map(int, self.modesForArrows.split()))
                modes = self.modesForArrows
            else:
                return
                
        '''    
        cmd.hide('everything')
        if struc.isCA():
            cmd.set('ribbon_trace',1)'''
        self.arrowObj = list()
        for i in modes:
            anm.arrowGenerator(mode = i, color = self.arrowRGB, size = self.arrowSize)
            arrowName = os.path.basename(self.structure.get()) + '_mode_' + str(i) + '_arrows'            
            self.arrowObj.append(arrowName)            
            cmd.load_cgo(anm.arrows,arrowName)
            self.arrowLoaded = True
    
    def dataExporter(self):
        '''
        export data to local directories
        '''        
        
        saveDirDialog = tkFileDialog.askdirectory(title = 'Save data in:',
                                                  initialdir = '',
                                                  parent = self.parent)
        if saveDirDialog:
            self.saveDir.set(saveDirDialog)
            os.chdir(self.saveDir.get())
            if self.exportE.get():
                filename = anm.proteinName + '_frequncies.txt'                
                np.savetxt(filename, anm.e, fmt = '%.3f')
            if self.exportV.get():
                filename = anm.proteinName + '_modes.txt'                
                np.savetxt(filename, anm.v, fmt = '%.3f')
            if self.exportCX.get() and anm.method == 'cutoff':
                filename = anm.proteinName + '_contactMatrix.txt'
                np.savetxt(filename, anm.getCX(), fmt = '%d')
            if self.exportCX.get() and anm.method == 'pf':
                filename = anm.proteinName + '_contactMatrix.txt'
                np.savetxt(filename, anm.getCX(), fmt = '%.3e')
            if self.exportHess.get():
                filename = anm.proteinName + '_Hessian.txt'
                np.savetxt(filename, anm.getHess(), fmt = '%.3f')
            if self.exportMSF.get():
                filename = anm.proteinName + '_MSF.txt'
                np.savetxt(filename, anm.getMSF(), fmt = '%.3f')
        
            message = 'Data has been saved in the following directory: ' + self.saveDir.get()
        
            tkMessageBox.showinfo("Data Saved", message, parent = self.parent)
    
    def execute(self, result):
        if result == 'Build ANM':
            self.ANMbuilder()
        elif result == 'Show Springs' and not self.anmBuilt:
            tkMessageBox.showerror('ANM not found','Please Build your ANM First')
        elif result == 'Show Springs' and self.anmBuilt:
            self.SpringsBuilder()
        elif result == 'Color by MSFs' and not self.anmBuilt:
            tkMessageBox.showerror('ANM not found','Please Build your ANM First')       
        elif result == 'Color by MSFs' and self.anmBuilt:
            self.MSFcolorer()
        elif result == 'Make Movies' and not self.anmBuilt:
            tkMessageBox.showerror('ANM not found','Please Build your ANM First')
        elif result == 'Make Movies' and self.anmBuilt:
            self.movieMaker()
        elif result == 'Draw Arrows' and not self.anmBuilt:
            tkMessageBox.showerror('ANM not found','Please Build your ANM First')        
        elif result == 'Draw Arrows' and self.anmBuilt:
            self.arrowDrawer()
        elif result == 'Export Data' and not self.anmBuilt:
            tkMessageBox.showerror('ANM not found','Please Build your ANM First')        
        elif result == 'Export Data' and self.anmBuilt:
            self.dataExporter()
        elif result == 'Exit':
            self.Quit()
            
    def Quit(self):
        self.dialog.destroy()

###################################
###### Simple Math Functions ######
###################################

def dis(cor_a, cor_b):
    '''
    returns distance between two points in 3 dimensions given the coordinates of the two points
    '''    
    summ = 0    
    for i in range (0,3):
        summ = summ + (cor_a[i]-cor_b[i]) ** 2
        
    distance = summ ** 0.5
    return distance

def mat2vec(np_matrix):
    '''
    transfers a matrix to a vector
    '''    
    vector = np.asarray(np_matrix).reshape(-1)
    return vector
    
def vec2mat(np_vector):
    '''
    transfers a vector to a 3D matrix
    '''
    
    length = len(np_vector)
    if length%3 != 0:
        raise ValueError('vector length cannot be divided by 3')
    else:
        matrix = np_vector.reshape((length/3,3))
    
    return matrix
    
##################################
########## ANM Class #############
##################################    

class ANM:
    def __init__(self, protein, cutoff = 12, useCA = True, method = 'cutoff', power = 4):
        '''
        initiation of class ANM
        '''        
        
        
        self.cx = None
        self.hess = None
        self.freq = None
        self.modes = None
        self.cutoff = cutoff
        self.e = None
        self.v = None
        self.proteinName = protein
        self.method = method
        self.power = power
        self.MSF = None
        self.arrows = list()
        global struc        
        try:
            struc = structure(protein)
        except IOError:
            if havePymol:
                obj_list = cmd.get_names()
                filename = protein + '_PyANM.pdb'
                if protein in obj_list:
                    cmd.save(filename, protein, state = -1 )
                else:
                    cmd.fetch(protein)
                    cmd.save(filename, protein, state = -1)
                struc = structure(filename)
                os.remove(filename)
            else:
                raise IOError('Failed to find object %s both in PDB and local directory' %protein)
              
        
        if useCA:
            struc.CA_only()
            if struc.isCA():
                pass
            else:
                raise IOError('Failed to extract all alpha carbons from %s' %protein)
        
        
    
    def getCX(self):
        if self.cx is None:
            warnings.warn("Contact Matrix not found, Rebuilding now...") 
            self.buildCX()
        
        return self.cx
            
    def resetCX(self):
        self.cx = None
        
    def getHess(self):
        if self.hess is None:
            warnings.warn("Hessian Matrix not found, Rebuilding now...")
            self.buildHess()
        
        return self.hess
        
    def resetHess(self):
        self.hess = None
        
    def getMSF(self):
        if self.MSF is None:
            self.calcMSF()
            
        return self.MSF
        
    def resetMSF(self):
        self.MSF = None
        
    def getE(self):
        if self.e is None:
            self.calcModes()
            
        return self.e
    
    def resetE(self):
        self.e = None
        
    def getV(self):
        if self.v is None:
            self.calcModes()
            
        return self.v
        
    def resetV(self):
        self.v = None
        
    def buildCX(self):
        '''
        build contact matrix
        '''
        if not isinstance(struc.cor, np.ndarray):        
            raise ValueError('Structure coordinates are not numpy arrays')
        if struc.cor.shape[1] != 3:
            raise ValueError('Dimension of coordinates does not equal 3')
            
        self.contactPair_x = list()
        self.contactPair_y = list()    
        length = struc.cor.shape[0]    
        self.cx = np.zeros((length,length),dtype = float)
        
        if self.method == 'cutoff':
            for i in range (1,length):
                for j in range (0,i):
                    temp_distance = dis(struc.cor[i,:], struc.cor[j,:])
                    if temp_distance < self.cutoff:
                        self.cx[i,j] = 1
                        self.contactPair_x.append(i)
                        self.contactPair_y.append(j)
        elif self.method == 'pf':
            for i in range (1,length):
                for j in range (0,i):
                    temp_distance = dis(struc.cor[i,:], struc.cor[j,:])
                    self.cx[i,j] = (1/temp_distance)**self.power
                    if temp_distance <= self.cutoff:
                        self.contactPair_x.append(i)
                        self.contactPair_y.append(j)
        else:
            raise IOError('Wrong ANM Method, should be either cutoff or pf')
                    
    def buildHess(self):
        '''
        build hessian matrix
        '''
        self.cx = self.getCX()
        length = self.cx.shape[0]
        contactNumber = len(self.contactPair_x)
        self.hess = np.zeros((3*length,3*length),dtype = float)
        for i in range (0,contactNumber):
            x = self.contactPair_x[i]
            y = self.contactPair_y[i]
            dR = struc.cor[y,:] - struc.cor[x,:]
            dR = dR/np.linalg.norm(dR)
            K33 = -dR[np.newaxis,:].T*dR
            self.hess[3*x:3*x+3,3*y:3*y+3] = K33
            self.hess[3*y:3*y+3,3*x:3*x+3] = K33
            self.hess[3*x:3*x+3,3*x:3*x+3] = self.hess[3*x:3*x+3,3*x:3*x+3] - K33
            self.hess[3*y:3*y+3,3*y:3*y+3] = self.hess[3*y:3*y+3,3*y:3*y+3] - K33
            
    def calcModes(self,numeig = None):
        '''
        calculat modes for ANM
        
        numeig: number of eigen values calculated
                default value is none
                all modes are calculated
        '''
        self.hess = self.getHess()
        if not isinstance(self.hess, np.ndarray):        
            raise ValueError('Hessian Matrix is not in the type of numpy arrays')
            
        if numeig is None:
            numeig = struc.cor.shape[0]*3
        
        if haveSP:
            try:
                [self.e,self.v] = scipy.sparse.linalg.eigsh(self.hess, numeig, which = 'SM', maxiter = 5000)
            except ValueError:
                [self.e,self.v] = scipy.sparse.linalg.eigsh(self.hess, numeig-1, which = 'SM', maxiter = 5000)
        else:
            [self.e,self.v] = np.linalg.eig(self.hess)
            self.v = self.v[:,self.e.argsort()]
            self.e.sort()
            self.v = self.v[:,0:numeig]
            self.e = self.e[0:numeig]
            
        if not self.freqChecker():
            warnings.warn("Not all eigenvalues converged, results might be inaccurate")
            print(self.e[0:7])
            
    
    def freqChecker(self):
        '''
        checks if the first 6 modes have a frequency of 0
        '''
        flag = False
        if sum(self.e[0:6]) < 1e-10 and sum(self.e[0:7]) > 1e-10:
            flag = True
        
        return flag
    
    def modeAnimator(self, modes = [1,2,3], useAllAtom = False, scaler = 5):
        '''
        make movies for different modes
        
        modes:        mode numbers to make movie for
                      default is first 3 modes
        scaler:       motion scales for movie
                      default is 5
        '''
        self.v = self.getV()
        self.movieFiles = {}
        numberTracker = 0

        if not useAllAtom:
            for i in modes:
                movieFile = self.proteinName.split('.')[0] + '_mode' + str(i) + '.pdb'
                self.movieFiles[numberTracker] = movieFile
                numberTracker += 1
                outfile = open(movieFile,'w')
                
                for j in range(11):
                    tempCor = struc.cor + scaler*j*vec2mat(self.v[:,i+5])
                    outfile.write('MODEL  %4d\n' %(j+1))
                    for k in range(struc.length):
                        outfile.write('%-6s%5d  %-4s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f           %-4s\n' %(struc.line_type[k], struc.atom_number[k], struc.atom_name[k], struc.residue_name[k], struc.chain_id[k], struc.residue_number[k],tempCor[k,0], tempCor[k,1], tempCor[k,2], struc.occupancy[k], struc.b_factor[k], struc.element[k]))
                    outfile.write('ENDMDL\n')
                
                for j in range(10):
                    tempCor = struc.cor + scaler*(10-j)*vec2mat(self.v[:,i+5])
                    outfile.write('MODEL  %4d\n' %(j+12))
                    for k in range(struc.length):
                        outfile.write('%-6s%5d  %-4s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f           %-4s\n' %(struc.line_type[k], struc.atom_number[k], struc.atom_name[k], struc.residue_name[k], struc.chain_id[k], struc.residue_number[k],tempCor[k,0], tempCor[k,1], tempCor[k,2], struc.occupancy[k], struc.b_factor[k], struc.element[k]))
                    outfile.write('ENDMDL\n')
                
                for j in range(11):
                    tempCor = struc.cor - scaler*j*vec2mat(self.v[:,i+5])
                    outfile.write('MODEL  %4d\n' %(j+22))
                    for k in range(struc.length):
                        outfile.write('%-6s%5d  %-4s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f           %-4s\n' %(struc.line_type[k], struc.atom_number[k], struc.atom_name[k], struc.residue_name[k], struc.chain_id[k], struc.residue_number[k],tempCor[k,0], tempCor[k,1], tempCor[k,2], struc.occupancy[k], struc.b_factor[k], struc.element[k]))
                    outfile.write('ENDMDL\n')
                    
                for j in range(10):
                    tempCor = struc.cor - scaler*(10-j)*vec2mat(self.v[:,i+5])
                    outfile.write('MODEL  %4d\n' %(j+33))
                    for k in range(struc.length):
                        outfile.write('%-6s%5d  %-4s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f           %-4s\n' %(struc.line_type[k], struc.atom_number[k], struc.atom_name[k], struc.residue_name[k], struc.chain_id[k], struc.residue_number[k],tempCor[k,0], tempCor[k,1], tempCor[k,2], struc.occupancy[k], struc.b_factor[k], struc.element[k]))
                    outfile.write('ENDMDL\n')
                
                outfile.close()
                
    def calcMSF(self):
        '''
        calculate the Mean Square Fluctuations from ANM
        '''
        self.cx = self.getCX()
        self.hess = self.getHess()
        
        if self.e is None:
            self.calcModes()
            
        self.MSF = np.zeros((1,struc.length), dtype = float)
        
        cv1 = range(0,struc.length*3,3)
        cv2 = range(1,struc.length*3,3)
        cv3 = range(2,struc.length*3,3)
        
        for i in range(6,struc.length*3):
            v_temp = self.v[:,i]**2
            self.MSF = self.MSF + (1/(self.e[i])) * (v_temp[cv1]+v_temp[cv2]+v_temp[cv3])
        
        T  = 300
        Kb = 1.38065e-23
        Gamma = 1e-20
        pi = 3.1415926
         
        self.MSF =  (8/3)*(pi**2)*Kb*(T/Gamma)*self.MSF
        
    def arrowGenerator(self, mode = 1, color = [1,1,1], size = 5):
        '''
        generate arrow objects for pymol
        
        mode:   the number of mode to draw arrows for
                default is the first mode
        color:  rgb color for the arrows
                default is white
        size:   size for the arrows
                default is 5
        '''
        self.v = self.getV()
        self.arrows = list()
        
        for i in range(struc.length):
            tail = [CYLINDER, struc.cor[i,0], struc.cor[i,1], struc.cor[i,2],\
                    struc.cor[i,0]+size*16*self.v[3*i,mode+5], struc.cor[i,1]+size*16*self.v[3*i+1,mode+5], struc.cor[i,2]+size*16*self.v[3*i+2,mode+5],\
                    0.3, color[0], color[1], color[2], color[0], color[1], color[2]]
            self.arrows.extend(tail)
            
            head = [CONE, struc.cor[i,0]+size*16*self.v[3*i,mode+5], struc.cor[i,1]+size*16*self.v[3*i+1,mode+5], struc.cor[i,2]+size*16*self.v[3*i+2,mode+5],\
                   struc.cor[i,0]+size*20*self.v[3*i,mode+5], struc.cor[i,1]+size*20*self.v[3*i+1,mode+5], struc.cor[i,2]+size*20*self.v[3*i+2,mode+5],\
                   1, 0.0, color[0], color[1], color[2], color[0], color[1], color[2], 1.0, 1.0]
            self.arrows.extend(head)
        #cmd.load_cgo(self.arrows,'test')
        
##################################
###### PDB Parsing Class #########
##################################

class structure:
    def __init__(self, filename, includeHET = False):
        '''
        A pdb parser for formal pdb format
        
        filename:   the pdb file to be parsed
        includeHET: if HETATM in pdb should be included, default is false (do not include)
        '''
        pdbfile = open(filename,'r')
        self.line_type = list()
        self.atom_number = list()
        self.atom_name = list()
        self.residue_name = list()
        self.chain_id = list()
        self.residue_number = list()
        self.cor_x = list()
        self.cor_y = list()
        self.cor_z = list()
        self.occupancy = list()       
        self.b_factor = list()
        self.element = list()
        
        for line in pdbfile:
            if includeHET:
                if (line[0:4] == 'ATOM') | (line[0:6] == 'HETATM'):
                    self.line_type.append(line[0:6].strip())                
                    self.atom_number.append(int(line[7:11].strip()))
                    self.atom_name.append(line[12:16].strip())
                    self.residue_name.append(line[17:21].strip())
                    self.chain_id.append(line[21])
                    self.residue_number.append(int(line[22:26].strip()))
                    self.cor_x.append(float(line[30:38].strip()))
                    self.cor_y.append(float(line[38:46].strip()))
                    self.cor_z.append(float(line[46:54].strip()))
                    self.occupancy.append(float(line[54:60].strip()))
                    self.b_factor.append(float(line[60:66].strip()))
                    self.element.append(line[77:80].strip())
            else:
                if (line[0:4] == 'ATOM'):
                    self.line_type.append(line[0:6].strip())                
                    self.atom_number.append(int(line[7:11].strip()))
                    self.atom_name.append(line[12:16].strip())
                    self.residue_name.append(line[17:21].strip())
                    self.chain_id.append(line[21])
                    self.residue_number.append(int(line[22:26].strip()))
                    self.cor_x.append(float(line[30:38].strip()))
                    self.cor_y.append(float(line[38:46].strip()))
                    self.cor_z.append(float(line[46:54].strip()))
                    self.occupancy.append(float(line[54:60].strip()))
                    self.b_factor.append(float(line[60:66].strip()))
                    self.element.append(line[77:80].strip())
        
        self.length = len(self.atom_name)
        self.cor = np.zeros((self.length,3),dtype=float)
        self.cor[:,0] = self.cor_x;
        self.cor[:,1] = self.cor_y;
        self.cor[:,2] = self.cor_z;
        pdbfile.close()
    
    def CA_only(self):
        '''
        parse only the alpha carbons
        '''
        temp_line_type = list()        
        temp_atom_number = list()
        temp_atom_name = list()
        temp_residue_name = list()
        temp_chain_id = list()
        temp_residue_number = list()
        temp_occupancy = list()
        temp_b_factor = list()
        temp_element = list()
        
        ca_index = list()
        for i in range(self.length):
            if self.atom_name[i] == 'CA':
                ca_index.append(i)
        
        temp_cor = self.cor[ca_index,:]
        ca_index.reverse()
        while ca_index:
            temp = ca_index.pop()
            temp_line_type.append(self.line_type[temp])            
            temp_atom_number.append(self.atom_number[temp])
            temp_atom_name.append(self.atom_name[temp])
            temp_residue_name.append(self.residue_name[temp])
            temp_chain_id.append(self.chain_id[temp])
            temp_residue_number.append(self.residue_number[temp])
            temp_occupancy.append(self.occupancy[temp])
            temp_b_factor.append(self.b_factor[temp])
            temp_element.append(self.element[temp])
        
        self.line_type = temp_line_type        
        self.atom_number = temp_atom_number
        self.atom_name = temp_atom_name
        self.residue_name = temp_residue_name
        self.chain_id = temp_chain_id
        self.residue_number = temp_residue_number
        self.occupancy = temp_occupancy
        self.b_factor = temp_b_factor
        self.element = temp_element
        self.cor = temp_cor
        self.length = len(self.atom_name)
        
    def isCA(self):
        '''
        returns true if the structure contains only alpha carbons
        '''
        ca_flag = True        
        for i in range(self.length):
            if self.atom_name[i] == 'CA':
                continue
            else:
                ca_flag = False
                break
        return ca_flag
        
    def writePDB(self,outname):
        '''
        write structure back to pdb format
        
        outname: output filename
        '''
        outfile = open(outname,'w')
        for i in range(self.length):
           outfile.write('%-6s%5d  %-4s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f           %-4s\n' %(self.line_type[i], self.atom_number[i], self.atom_name[i], self.residue_name[i], self.chain_id[i], self.residue_number[i],self.cor[i,0], self.cor[i,1], self.cor[i,2], self.occupancy[i], self.b_factor[i], self.element[i]))
       
        outfile.close()
        
        
