import sys
import Tkinter as tk
from Tkinter import *
#Note: By default, the widgets like button, frame, radiobutton are all ttk widgets.
#To implement the widgets from Tkinter, specify tk.Button, tk.Frame etc...!
import ttk
from ttk import *
import tkMessageBox
import tkFileDialog 
import distutils
import os
import subprocess

import logging
from logging import handlers
#log file is used to debug the code when necessary

from pymol import cmd, plugins

import time
import multiprocessing
import datetime

import shutil
 # Use for checking License Status- About Tab
import site
import distutils.core

from distutils.errors import DistutilsFileError
from calendar import month

import pymol.plugins

startupDirectories=pymol.plugins.get_startup_path()
for pyDir in startupDirectories:
    sys.path.append(pyDir)
    print "pydir1 = ", pyDir.encode('string-escape')
    print "pydir2 = ", pyDir

print "sys.path = ", sys.path

import lisica

configure=lisica.Configuration()
exe_filename=configure.exe_File()

lisicaFolder=lisica.LISICA_DIRECTORY
resultsFolder=os.path.abspath(os.curdir)
logFolder=os.path.join(lisicaFolder,"Log")
iconFolder=os.path.join(lisicaFolder,"Icons")

exe_path=os.path.normpath(os.path.join(lisicaFolder,"bin",exe_filename))
#process_id=os.getpid()
#process id may not be unique for the plugin instances. PyMOL uses multithreading

timestamp=datetime.datetime.fromtimestamp(time.time())
log_filename="log_"+str(timestamp.strftime('%m%d%H%M%S'))+".log"
## log file for each lisica process 

####################################################################
####                                                            ####
####      Setting up the log file for LiSiCA                    ####
####                                                            ####
####################################################################

### Configuring a logger named lisica
log=logging.getLogger("lisica")
### logger data is written to the file licisa.log in the LiSiCA folder in C:/Program files
handle=logging.FileHandler(os.path.join(logFolder,log_filename))
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
handle.setFormatter(formatter)
log.addHandler(handle)
log.setLevel(logging.DEBUG)



sys.path.append(lisicaFolder)
    

try:
    log.info("This application uses Tkinter Treectrl module")
    log.info("http://tkintertreectrl.sourceforge.net/")
    import TkTreectrl
    
except ImportError:
    #catches the exception raised if the TkTreectrl package (Tkinter wrapper of Tktreectrl tcl module) is incorrectly installed 
    #TkTreectrl package is not a standard python package and has to be installed separately
    #Installation normally done by running setup.py file
    log.critical("Error Message: '%s'" %ImportError.message)
    log.debug("TkTreectrl Import error")
    print "Import Error: Tkinter Treectrl not installed correctly"
    print("Copy and paste the TkTreectrl folder from the installation directory into the site-packages folder in Python's lib folder")
except:
    #catches all other exceptions
    #Possible exception: The tcl module tktreectrl is missing.
    #Installation of tktreectrl  normally done by running setup.py file
    log.critical("Error Message: '%s'" %Exception.message)
    log.debug("Possible: Tcl- treectrl error")
    print "Exception caught: Not Import Error"
    print "Possible: Tcl- treectrl error "
    print("Copy and paste the treectrl folder from the installation directory into the site-packages folder in Python's tcl folder")
    

####################################################################
####                                                            ####
####      Function to delete old files                          ####
####                                                            ####
####################################################################

def deleteOldFiles(folder):
    try:
        deleted=0
        if folder==logFolder:
            days=7
        
        now=time.time()
        for f in os.listdir(folder):
            f=os.path.join(folder, f)
        
            if os.stat(f).st_mtime < now - days * 86400:
                #if the contents of the file/folder is last modified before 'days' days
                if os.stat(f).st_atime < now - 5 * 86400:
                    #if the contents of the file/folder is last accessed before 5 days
                    if os.path.isfile(f):
                        os.remove(f) 
                        deleted+=1   
                    elif os.path.isdir(f):
                        os.rmdir(f)
                        deleted+=1
    except IOError:
        log.error("IO Error")
    except OSError:
        log.error("OSError")
    except:
        log.error("Some error while deleting old log files")
    else:
        log.info("Old log files removed")
    



class GUI(Frame):
    def __init__(self,parent):
        Frame.__init__(self, parent)
        
        self.plugin = parent
        ### Configuring the parent frame (Plugin Frame) of LiSiCA
        #~ self.plugin.geometry('{}x{}'.format(750,680))
        self.plugin.title("LiSiCA - Ligand Similarity using Clique Algorithm")
        
        try:
            self.plugin.iconbitmap(os.path.join(iconFolder,"lisica_icon.ico"))
        except TclError:
            log.error(TclError.message)
            log.debug("Icon file not found in the path specified")
        except:
            log.debug("Not TclError")
            #Unix/Linux-@name.xbm
            
        self.initUI()    

    def close(self):
        #Command to close LiSiCA plugin
        
        if  tkMessageBox.askyesno("Exit", "Exit LiSiCA"):
            try:   
                #If the lisica.exe process is running in the background, terminates it.
                self.proc.kill()
                log.info("lisica.exe terminated")
            except:
                pass
            self.plugin.destroy()
            
    def getRefFileName(self):
        self.ref_Entry
        self.Rfilename=tkFileDialog.askopenfilename(filetypes = [('Mol2 files','*.mol2')])
        if os.path.isfile(self.Rfilename):
            self.ref_Entry.delete(0,END)
            self.ref_Entry.insert(0,self.Rfilename)
        self.updateCommand()
    def getTarFileName(self):
    #get Target Ligand database file name 
        self.tar_Entry
        self.Tfilename=tkFileDialog.askopenfilename(filetypes = [('Mol2 files','*.mol2')])
        if os.path.isfile(self.Tfilename):
            self.tar_Entry.delete(0,END)
            self.tar_Entry.insert(0,self.Tfilename)
        self.updateCommand()
    
    def dim(self):
    #get the product graph dimension from the radio button
        if self.dimension.get()==3:
            self.update()
            self.parameter.set("Maximum allowed Spatial Distance Difference: ")
            self.units.set("?") #Not used unicode error
            # for 3D screening, display an extra parameter "No of conformations"
            self.conf_Label.grid(row=8,sticky=W,pady=10,padx=20) #label
            self.conformations.grid(row=8,column=2,sticky=W) #Entry widget
        else:
            self.update()
            self.parameter.set("Maximum allowed Shortest Path difference: ")#initial setting req-Default 2D?
            self.units.set("bonds")#initial setting req?
            # for 2D screening, remove the parameter "No of conformations"
            self.conf_Label.grid_remove()
            self.conformations.grid_remove()
        self.updateCommand()

        
            
    def inputValidation(self):
        msg=" "
        if os.path.isfile(self.ref_Entry.get()):
            if self.ref_Entry.get().endswith(".mol2"):
                log.info("Ref molecule file: '%s'",self.ref_Entry.get())
            else:
                log.error("Wrong file format")
                self.ref_Entry.config(highlightbackground='#B0171F')
                self.update()
        else:
            log.error("Ref molecule file path is invalid")
            msg="Invalid path for reference molecule\n"
            
            self.update()
        if os.path.isfile(self.tar_Entry.get()):
            if self.tar_Entry.get().endswith(".mol2"):
                log.info("Target molecule file: '%s'",self.tar_Entry.get())
            else:
                log.error("Wrong file format")
        else:
            log.error("Target molecule file path is invalid")
            msg=msg + "\nInvalid path for Target molecule"
            
        if os.path.isdir(self.saveResultsEntry.get()):
            log.info("The results will be written in the folder: '%s' " %self.saveResultsEntry.get())
        else:    
            log.debug("Folder name does not exist..!")
            self.saveResultsEntry.delete(0,END)
            self.saveResultsEntry.insert(0, resultsFolder) 
            msg=msg + "\nInvalid path for saving the results"
        return msg
    
    def createCommand(self):
        self.lisica_path=exe_path
        
        self.command=self.lisica_path\
                    +" -R "+os.path.normpath(self.ref_Entry.get())\
                    +" -T "+os.path.normpath(self.tar_Entry.get())\
                    +" -n "+str(self.CPU_Cores.get())\
                    +" -d "+str(self.dimension.get())\
                    +" -w "+self.w_Entry.get()\
                    +" -f "+self.saveResultsEntry.get()
                    #+" -h "+str(self.hydrogen.get())
        
        if self.hydrogen.get()==1:
            self.command=self.command+" -h "
        if self.dimension.get()==2:
            self.command=self.command+" -s "+str(self.max_Entry.get())
        elif self.dimension.get()==3:
            self.command=self.command\
                        +" -m "+str(self.max_Entry.get())\
                        +" -c "+str(self.conformations.get())
        else:
            log.error("Dimension is neither 2 or 3..!")
            
        self.command=self.command+" --plugin"        
        log.info("Command for lisica: '%s'" ,self.command)
        return self.command
                    
            
    def submitted(self):
        self.updateCommand()
        try:
            alert=self.inputValidation()
            if alert==" ":
                self.createCommand()
                self.go_Button.lower(self.input_Tab)
            
                self.progress=Progressbar(self.note,mode='indeterminate',length=500)
                self.progress.grid(row=17,columnspan=6,sticky=W+E,padx=(40,5),in_=self.input_Tab)
            
                print self.command
    
                self.startupinfo=None
                if os.name=='nt':
                    try:
                        self.startupinfo=subprocess.STARTUPINFO()
                        self.startupinfo.dwFlags|=subprocess.STARTF_USESHOWWINDOW
                    except:
                        self.startupinfo.dwFlags|=subprocess._subprocess.STARTF_USESHOWWINDOW
                
                os.chdir( lisicaFolder ) 
                self.proc = subprocess.Popen(self.command.split(),shell=False,startupinfo=self.startupinfo)
                #self.proc = subprocess.Popen(self.command.split(),shell=False)
               
                self.progress.start()
            
            
                while self.proc.poll() is None:
                    try:
                        self.update()
                    except:
                        log.info("Application has been terminated")
                self.progress.stop()
                self.progress.lower(self.input_Tab)
                self.go_Button.lift(self.input_Tab)
                if os.path.isfile(os.path.join(lisicaFolder,'done.txt')):
                    
                    #Read the result folder name from the done.txt file
                    done_File=open(os.path.join(lisicaFolder,'done.txt'))
                    self.timestamp_FolderName=done_File.readline()
                    print self.timestamp_FolderName
                    #timestamp_FolderName=self.timestamp_FolderName[:-1]
                    self.timestamp_Folder=os.path.join(self.saveResultsEntry.get(),self.timestamp_FolderName)
                    done_File.close()
                    
    
                    self.Ref=self.ref_Entry.get()
    
                    print "before hokuspokus =", self.timestamp_Folder
                    #~ self.timestamp_Folder=self.timestamp_Folder.encode('string-escape')
                    self.timestamp_Folder=self.timestamp_Folder[:-1]
                    print "after =", self.timestamp_Folder
                    self.addRefFileToResults()
    
                    
                    log.debug("Execution of liSiCA is complete")
                    log.info("Files with common structures are available in results directory")
                    self.note.select(self.output_Tab)
                    self.uploadResults()
                
                    try:
                        os.remove(os.path.join(lisicaFolder,'done.txt'))
                    except IOError:
                            log.exception(IOError.message)
                else:
                    tkMessageBox.showerror("Alert", "Abnormal termination of LiSiCA")
            
                
            else:
                tkMessageBox.showwarning("Alert", alert)
        except:
            print "note: something went wrong in the middle of calculation"

    def getResultsDir(self,arg):
        global resultsFolder
        rfolder=tkFileDialog.askdirectory(initialdir=resultsFolder)
        if os.path.isdir(rfolder):
            resultsFolder=rfolder #set new path globally
            if arg==0:
                self.saveResultsEntry.delete(0,END)
                self.saveResultsEntry.insert(0,resultsFolder)
            if arg==1:
                self.loadResultsEntry.delete(0,END)
                self.loadResultsEntry.insert(0,resultsFolder)
        self.updateCommand()
    def loadResults(self):
            
        self.timestamp_Folder=self.loadResultsEntry.get()
        self.Ref=os.path.join(self.timestamp_Folder,"Ref.mol2")
        self.note.select(self.output_Tab)
        self.uploadResults()

    def uploadResults(self):
        #Upload the first 100 results in lisica_result.txt to the output tab
        self.multiListBox1.delete(0,END)
        self.multiListBox2.delete(0,END)
        row=()
        compd_num=0
        try: 
            print self.w_Entry
            print self.w_Entry.get()
            '''
            print "before hokuspokus =", self.timestamp_Folder
            self.timestamp_Folder=self.timestamp_Folder.encode('string-escape')
            self.timestamp_Folder=self.timestamp_Folder[:-1]
            print "after =", self.timestamp_Folder
            '''
            result_File=os.path.normpath(os.path.join(self.timestamp_Folder,r"lisica_results.txt"))
            #result_File=result_File.encode('string-escape')
            print result_File
            Textfile=open(result_File,"r") 
        #the file location -dynamically allocate..!!!!
            log.info("No IO Error. File '%s' open for reading" %Textfile.name)
        
            for block in iter(lambda: Textfile.readline(), ""):
                compd_num=compd_num+1
                row=block.split()
                rankedList=(compd_num,row[1],row[0])
                self.multiListBox1.insert(END,*rankedList)
                #print compd_num
                if compd_num==int(self.w_Entry.get()):
                    print "Enters"
                    #if not self.w_Entry.get()==0:#displays only the first 100 results from the text file 
                    break
            log.info("Has finished uploading data")
            Textfile.close()
            self.multiListBox1.focus()
            self.multiListBox1.select_set(0, None)
            
        except IOError:
            log.error(IOError.errno)
            log.error(IOError.message)
            log.error(IOError.filename)
    
    def showCorr(self,ar):
        
        try:
            
            cmd.reinitialize()
            self.multiListBox1.focus()
            
            cmd.load(self.Ref,object="obj1")
            self.index=ar[0]
            self.selection_row=self.multiListBox1.get(self.index)
            self.selection_data=self.selection_row[0]
            self.rank=self.selection_data[0]
            self.zincId=self.selection_data[1]
            
            
            #Empty box2
            self.multiListBox2.delete(0,END)
            
            file_type=".mol2"
            result_Filename=self.rank+ "_" +self.zincId+file_type 
            log.info("Filename:'%s'" %result_Filename)
            i=0
            j=0
            k=0
            line=()
            #timestamp_Folder=os.path.join(resultsFolder,timestamp_Done)
            file_Mol=os.path.join(self.timestamp_Folder,result_Filename)
            print file_Mol
           
            self.ref=()
            self.tar=()
            
            try:  
                fh=open(file_Mol,"r") #handle possible exceptions
                cmd.load(file_Mol,object="obj2")
            except IOError:
                tkMessageBox.showinfo("No data found", "The result molecules were not written to the output")
            for block in iter(lambda: fh.readline(), ""):
                if i==1:
                    if j==1:
                        if not block=='\n':
                        
                            k=k+1
                            line=block.split()
                            self.ref=self.ref+(line[1],)
                            self.tar=self.tar+ (line[3],)
                            
                            self.multiListBox2.insert(END,*line)
                        else:
                            i=0
                            j=0
                if "LiSiCA  RESULTS" in block:
                    i=1
                if i==1:
                    if "--------------" in block:
                        j=1
            log.info("Displayed the atom correspondence of the selected atom")
                     
            
            if self.dimension.get()==2:
                cmd.set("grid_mode",value=1)
                cmd.set("grid_slot",value=1,selection="obj1")
                cmd.set("grid_slot",value=2,selection="obj2")
                cmd.align("obj1","obj2")
            elif self.dimension.get()==3:
                selection=[]
                for idx, item in enumerate(self.ref):
                    selection.append("'obj1////" + self.ref[idx] + "'")
                    selection.append("'obj2////" + self.tar[idx] + "'")
                fitCode="cmd.pair_fit(" + ",".join(selection) + ")"
                exec fitCode in globals(), locals()
            else:
                log.error("Dimension value is wrong")    
            
            cmd.orient()
            cmd.turn("z", 90)
        
        except:
            pass
        
    def addRefFileToResults(self):
        print self.timestamp_Folder
        #~ self.timestamp_Folder=self.timestamp_Folder.encode('string-escape')
        #~ self.timestamp_Folder=self.timestamp_Folder[:-1]
        print self.timestamp_Folder
        
        shutil.copyfile(self.ref_Entry.get(),os.path.normpath(os.path.join(self.timestamp_Folder,"Ref.mol2")))
        
        
    def showCorrAtoms(self,ar):
        try:
            self.multiListBox2.select_set(ar[0], None)
            self.sel_index=self.multiListBox2.get(ar[0])
            self.selected=self.sel_index[0]
            self.ref_Atom=self.selected[1]
            self.tar_Atom=self.selected[3]
            self.selected_Pair="obj1////"+self.ref_Atom+ " + obj2////"+self.tar_Atom
            cmd.select("ref",selection=self.selected_Pair)
            
        except:
            print("None selected from box2")

    def onEvent(self,event):

        self.updateCommand()
        
    def updateCommand(self):
        self.input_Tab.update()
		
        x=self.createCommand()
        self.display_Command.delete(1.0, END)
        self.display_Command.insert(1.0, x) 
    
    def retag(self,tag, *args):
        for widget in args:
            x0=widget.bindtags()[0]
            x1=widget.bindtags()[1]
            x2=widget.bindtags()[2]
            x3=widget.bindtags()[3] 
            widget.bindtags((x1,x0,tag, x2,x3))
            
            
    
####################################################################
####                                                            ####
####                GUI design                                  ####
####                                                            ####
####################################################################

    
    def initUI(self):
        #Frame to display the name and icon of the software
        self.frame=Frame(self.plugin)
        self.frame.pack(fill=BOTH)
        
        
        try:
            photo=PhotoImage(master=self,file=os.path.join(iconFolder,"lisica_icon.gif"))
            self.display=Label(self.frame,image=photo,background="#BFBFBF")
            self.display.image=photo
            self.display.pack(side=RIGHT)
        except:
            pass
        self.heading=Label(self.frame,text="LiSiCA",font=("Times",30, "bold"),
                           foreground="brown", background="#BFBFBF",anchor=CENTER)
        
        self.heading.pack(ipady=10, fill=BOTH, expand=1)
        
       
        #Notebook with 3 tabs
        self.note=Notebook(self.plugin)
        
        self.input_Tab=Frame(self.note)
        self.load_Project_Tab=Frame(self.note)
        
        self.output_Tab=Frame(self.note)
        self.about_Tab=Frame(self.note)
        
        
        
        self.note.add(self.input_Tab,text="   Inputs     ")
        self.note.add(self.load_Project_Tab,text=" Load Project ")
        self.note.add(self.output_Tab, text="   Outputs   ")
        self.note.add(self.about_Tab,text="   About     ")
        self.note.pack(padx=10,pady=10, ipadx=10,ipady=10,anchor=CENTER,fill=BOTH,expand=1,after=self.frame)
        
        #Button to close LiSiCA plugin
        #self.closeButton=Button(self.plugin, text='Close', command=self.close)
        #self.closeButton.pack(pady=10,after=self.note,ipady=2)

#Load Project Tab Design

        Label(self.load_Project_Tab, text="Choose the results directory :").grid(row=1,sticky=W,padx=20,pady=(40))
        self.loadResultsEntry=Entry(self.load_Project_Tab,width=50)
        self.loadResultsEntry.grid(row=1,column=2)
        self.loadResultsEntry.insert(0, resultsFolder) 
       
        Button(self.load_Project_Tab,text="Browse",command=lambda: self.getResultsDir(1)).grid(row=1,column=4)
        Button(self.load_Project_Tab,text="Load",command=self.loadResults).grid(row=2,column=4)
        


       
        
#About Tab Design
        #label_About_style=ttk.Style()
        #label_About_style.configure('AboutTabLabel.TLabel', font=('Helvetica', 12, 'Bold'),foreground="black", background="white")
        
        self.version_Frame=tk.LabelFrame(self.about_Tab,text="Product Version Information",labelanchor="nw",font=("Times", 12),relief="ridge",borderwidth=4)
        
        self.version_Frame.pack(fill=X,padx=(10,10),pady=(20,20))

        from lisica import UpgraderGitlab
        upgraderObj=UpgraderGitlab()
        upgraderObj.findCurrentVersion()
        upgraderObj.findLatestVersionGUI()

        self.currentVersionGUI=StringVar(master=self)
        self.currentVersionGUI.set(upgraderObj.currentVersionGUI)

        #~ upgraderObj.latestVersion = upgraderObj.latestVersion[1:]
        #~ print "lversion = ", upgraderObj.latestVersion
        
        Label(self.version_Frame,text="LiSiCA GUI Version :",font=("Times",11)).grid(row=7,columnspan=2,padx=(10,10),pady=(10,10),sticky=W)
        Label(self.version_Frame,textvariable=self.currentVersionGUI,font=("Times",11)).grid(row=7,column=2,columnspan=2,padx=(10,10),pady=(10,10),sticky=W)
        
        #Checking for latest upgrades
        from distutils.version import StrictVersion
        if StrictVersion(upgraderObj.currentVersionGUI) < StrictVersion(upgraderObj.latestVersionGUI):
            Label(self.version_Frame,text="Upgrade of LiSiCA GUI is available!",font=("Times",11,'bold')).grid(row=8,column=0,columnspan=2,padx=(10,10),pady=(10,10),sticky=W)

            def doUpgrade():
                upgraderObj.upgrade()
                self.close()
                
            self.upgrade_Button=Button(self.version_Frame,text="Upgrade to version " + upgraderObj.latestVersionGUI,command=lambda: doUpgrade())
            self.upgrade_Button.grid(row=9,column=3,sticky=E,pady=(5,5),padx=(5,5))
            
        
        
        self.contact_Frame=tk.LabelFrame(self.about_Tab,text="Contact Us",labelanchor="nw",font=("Times", 12),relief="ridge",borderwidth=4)
        self.contact_Frame.pack(fill=X,padx=(10,10),pady=(20,20))
        Label(self.contact_Frame,text="Contact Us at info@insilab.com",font=("Times",11)).grid(columnspan=2,padx=(10,10),pady=(10,10),sticky=W)
        Label(self.contact_Frame,text="Please feel free to visit our website : ",font=("Times",11)).grid(row=2,columnspan=2,padx=(10,5),pady=(10,10),sticky=W)
        
        
        self.website=Label(self.contact_Frame,text=r"http://www.insilab.com",font=("Times",11),foreground="blue",underline=True)
        import tkFont
        hyperlink_font=tkFont.Font(self.website,self.website.cget("font"))
        hyperlink_font.configure(underline=True)
        self.website.configure(font=hyperlink_font)
        
        self.website.grid(row=2,column=2,columnspan=2,padx=(0,10),pady=(10,10),sticky=W)
        self.website.bind("<Button-1>", self.openWebsite)
#Input Tab Design
        
        #Choose Mol2 file for Reference Ligand 
        Label(self.input_Tab,text="Reference Ligand: ").grid(row=1,sticky=W,pady=(20,10),padx=20)
        self.ref_Entry=Entry(self.input_Tab,width=50)
        
        self.ref_Entry.grid(row=1,column=2,pady=(20,10))
        
        self.button_BrowseRef=Button(self.input_Tab,text="Browse",command=self.getRefFileName)
        self.button_BrowseRef.grid(row=1,column=4,pady=(20,10))
        
        
        #Choose Mol2 file for Target Ligand(s)
        Label(self.input_Tab, text="Target Ligand(s):").grid(row=2,sticky=W,padx=20)
        self.tar_Entry=Entry(self.input_Tab,width=50)
        self.tar_Entry.grid(row=2,column=2)
        Button(self.input_Tab,text="Browse",command=self.getTarFileName).grid(row=2,column=4)
        
        #set variables for dimension parameters
        self.parameter=StringVar(master=self)
        self.units=StringVar(master=self)
        #set default for parameters
        self.parameter.set("Maximum allowed Shortest Path difference: ")
        self.units.set("bonds")
        
        global Conf_Label,conformations 
        #declaring these widget variables outside the function dim()
        #prevents repeated creation of widgets.

        self.conf_Label=Label(self.input_Tab, text="No of conformations: ")
        self.conformations=Entry(self.input_Tab, width=10)
        self.conformations.insert(0, "1")

        #2D or 3D Screening? Default is 2D- RadioButtons
        self.dimension=IntVar(master=self)
        self.dimension.set(2)
       
        self.d2=Radiobutton(self.input_Tab,text="2 Dimensional Screening", variable=self.dimension, value=2, command=self.dim)
        self.d3=Radiobutton(self.input_Tab,text="3 Dimensional Screening", variable=self.dimension, value=3, command=self.dim)
        Label(self.input_Tab,text="Product Graph Dimension:").grid(row=5,columnspan=2,sticky=W,pady=10,padx=20)
        self.d2.grid(row=5,column=2,columnspan=2,sticky=W+E)
        self.d3.grid(row=6,column=2,columnspan=2,sticky=W+E)
        
        #Depending on the dimension chosen 
        self.maxLabel=Label(self.input_Tab, textvariable=self.parameter)
        self.maxLabel.grid(row=7,sticky=W,pady=10,padx=20)
        self.maxLabel.update()
        self.max_Entry=Entry(self.input_Tab, width=10)
        self.max_Entry.grid(row=7,column=2,sticky=W)
        self.max_Entry.insert(0,"1")
       
        #Value to be acquired and passed
            #To display the unit, uncomment the following line
        #Label(self.input_Tab, textvariable=self.units).grid(row=7,column=2)
           
        
        self.w_Label=Label(self.input_Tab, text="No of highest ranked molecules to write to the output: ")
        self.w_Entry=Entry(self.input_Tab, width=10)
        
        
        self.w_Label.grid(row=9,sticky=W,pady=10,padx=20)
        self.w_Entry.grid(row=9,column=2,sticky=W)
        self.w_Entry.insert(0, "100")
        
        
        self.hydrogen=IntVar(master=self)
        
        
        Label(self.input_Tab, text="Consider Hydrogens :").grid(row=11,sticky=W,pady=10,padx=20)
        self.checkBox=Checkbutton(self.input_Tab,variable=self.hydrogen,command=self.updateCommand)
       
        self.checkBox.grid(row=11,column=2,sticky=W)
        
        #Get the no of CPU cores available in the system
        self.CPU_Cores=IntVar(master=self)
        self.CPU_Cores.set(multiprocessing.cpu_count())#Deafult=all CPU cores
        #Widget to set the CPU count manually- Slider
        Label(self.input_Tab, text="Number of CPU cores to be used: ").grid(row=12,sticky=W,pady=10,padx=20)
            #Defualt = All available CPU cores 
        self.slider=tk.Scale(self.input_Tab,
                          variable=self.CPU_Cores,
                          from_=1,
                          to=multiprocessing.cpu_count(),
                          length=100,
                          resolution=1,
                          orient=HORIZONTAL,
                          command=self.onEvent)
        self.slider.grid(row=12,column=2,pady=10,sticky=W)

        #Option to choose the directory for Results of LiSiCA
        Label(self.input_Tab, text="Save results in:").grid(row=13,sticky=W,padx=20)
        self.saveResultsEntry=Entry(self.input_Tab,width=50)
        self.saveResultsEntry.grid(row=13,column=2,pady=10)
        self.saveResultsEntry.insert(0, resultsFolder) 
        Button(self.input_Tab,text="Browse",command=lambda: self.getResultsDir(0)).grid(row=13,column=4)
        
        self.display_Command=Text(self.input_Tab,height=3)
	self.display_Command.grid(row=14,columnspan=6,sticky=W+E)
	
	self.retag("Parameters", self.ref_Entry, self.tar_Entry, self.conformations,self.w_Entry,self.max_Entry,self.saveResultsEntry,self.input_Tab)	
	self.bind_class("Parameters", "<Button-1>", self.onEvent)
	self.bind_class("Parameters", "<KeyPress>", self.onEvent)
	
	self.go_Button=Button(self.note,text="GO",command=self.submitted)
	self.go_Button.grid(row=17,column=4,pady=(5,5),in_=self.input_Tab)
    
    # Output Tab Design
        #Frame on the left half of Output Tab for displaying the Results
        self.frame_Result=Frame(self.output_Tab)
        self.frame_Result.pack(side=LEFT,fill=BOTH,expand=1)
        #Frame on the right half Output Tab for displaying the Atom correspondence
        self.frame_Corr=Frame(self.output_Tab)
        self.frame_Corr.pack(side=RIGHT,fill=BOTH,expand=1)

        #Design of the Result frame
        self.ref=()
        self.tar=()
        
            #Add a scroll bar named scrollbar1
        self.scrollbar1=Scrollbar(self.frame_Result)
        self.scrollbar1.pack(side=RIGHT,fill=Y)
            #Add a multi-column list box to the frame
            #The MultiListBox class belongs to TkTreectrl module
            #TkTreectrl package is an external package and is to be downloaded and installed  
        self.multiListBox1=TkTreectrl.MultiListbox(self.frame_Result,yscrollcommand=self.scrollbar1.set,
                                                   expandcolumns=(1,2),selectcmd=self.showCorr)
        self.multiListBox1.pack(side=LEFT,fill=BOTH,expand=1)
        self.multiListBox1.config(columns=('Rank','ZINC ID', "Tanimoto score"))
            #The scrollbar2 is linked to the multiListBox2
        self.scrollbar1.config(command=self.multiListBox1.yview)

        #Design of the Atoms' Correspondence displaying frame
            #Add a scroll bar named scrollbar2
        self.scrollbar2=Scrollbar(self.frame_Corr)
        self.scrollbar2.pack(side=RIGHT,fill=Y)
            #Add a multi-column list box to the frame
            #The MultiListBox class belongs to TkTreectrl module
            #TkTreectrl package is an external package and is to be downloaded and installed  
        self.multiListBox2=TkTreectrl.MultiListbox(self.frame_Corr,yscrollcommand=self.scrollbar2.set,expandcolumns=(0,1,2,3,4),selectcmd=self.showCorrAtoms)
        self.multiListBox2.pack(side=LEFT,fill=BOTH,expand=1)
        self.multiListBox2.config(columns=('Ref. Num','Ref. Atom','Tar. Num','Tar. Atom','Atom Type'))
            #The scrollbar2 is linked to the multiListBox2
        self.scrollbar2.config(command=self.multiListBox2.yview)
        
        
        self.pack()
        
    def openWebsite(self,event):
        import webbrowser
        try:
            webbrowser.open_new(r"http://www.insilab.com")
        except:
            print "Error. Could not open the Website"
#######################################################################################################################


def main():
    import platform

    #For each 
    if os.path.isfile(os.path.join(lisicaFolder,'lisica.log')):
        open(os.path.join(lisicaFolder,'lisica.log'),"w").close()
    
    
            
    log.info("*******************")
    log.info("LiSiCA: version 1.0") 
    log.info("Plugin_timestamp: '%s'" %timestamp.strftime("%Y-%m-%d %H: %M :%S"))
    root = Tk()

    treectrl=os.path.join(lisicaFolder,'treectrl2.4.1')
    print "system = ", platform.system()
    if platform.system()=='Linux':
        treectrl=os.path.join(lisicaFolder,'treectrl2.4.1-linux')

    root.tk.eval('lappend auto_path {'+treectrl+'}')
    lisica=GUI(root)
    root.protocol('WM_DELETE_WINDOW', lisica.close)
    root.mainloop()
#if __name__ == "__main__":
    
    #main()
    #deleteOldFiles(logFolder)

#def __init__(self):
    #main()
    #self.menuBar.addmenuitem('Plugin', 'command', 'LiSiCA', label = 'LiSiCA', command=lambda s=self: main())
