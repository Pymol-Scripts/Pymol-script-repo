# This modules handels updates for MGLTools
# Author: Sargis Dallakyan (sargis at scripps.edu)
# $Header: /opt/cvs/Support/update.py,v 1.37.4.1 2011/05/24 22:57:20 annao Exp $
# $Id: update.py,v 1.37.4.1 2011/05/24 22:57:20 annao Exp $
import os, sys, platform, urllib, tarfile, unzip, Tkinter, Pmw
from tkFileDialog import *
from tkMessageBox import *
import pickle, httplib, urllib
import webbrowser
base_url = 'http://mgltools.scripps.edu/downloads/tars/releases/nightly'
import mglutil
from mglutil.util.misc import ensureFontCase

#Paltform Independent Packages
share_packages = ['AutoDockTools', 'mglutil', 'NetworkEditor', 'PyAutoDock', 
                  'symserv', 'Vision', 'DejaVu', 'MolKit', 'Pmv', 'PyBabel', 
                  'ViewerFramework', 'Volume']
#Paltform Dependent Packages
packages = ['bhtree', 'cMolKit', 'gle', 'opengltk', 'QSlimLib', 'stride',
            'binaries', 'geomutils', 'mslib', 'pyglf', 'sff', 'UTpackages']

def rm_dirs(path):
    """source http://docs.python.org/lib/os-file-dir.html
    Delete everything reachable from the directory named in 'top',
    assuming there are no symbolic links.
    CAUTION:  This is dangerous!  For example, if top == '/', it
    could delete all your disk files."""
    for root, dirs, files in os.walk(path, topdown=False):
        for name in files:
            try:
                os.remove(os.path.join(root, name))
            except OSError:
                pass
        for name in dirs:
            try:
                os.rmdir(os.path.join(root, name))
            except OSError:
                pass   
    try:
        os.removedirs(path) # added to remove path itself
    except OSError:
        pass
        
    
class Update:
    """
    package    : Support
    module     : update
    class      : Update
    description:
        Provides GUI and command line interface for updating MGLTools
    """
    def __init__(self):
        if sys.platform == 'linux2':
            if sys.maxint> 2147483647: #os.popen('uname -m').read() == 'x86_64\n':
                self.sys_prefix = 'x86_64Linux2'
            else:
                self.sys_prefix = 'i86Linux2'
            self.update_file = self.sys_prefix+'.tar.gz'
        elif os.name == 'nt': #sys.platform == 'win32':
            self.sys_prefix = 'win32'
            self.update_file = self.sys_prefix+'.zip'
        elif sys.platform == 'darwin':
            uname  = os.uname()
            if uname[2].startswith('8'):
                name = '8'
            else:
                name = '9'
            if uname[-1] == "i386":
                self.sys_prefix = "i86Darwin"+name
            else:
                self.sys_prefix = "ppcDarwin"+name            
            self.update_file = self.sys_prefix+'.tar.gz'
        else:
            print "Sorry no nightly builds is available for your platform."
            print sys.platform
        from user import home
        from Support.version import __version__
        self.rc = home + os.sep + ".mgltools" + os.sep + __version__
        self.updates_rc_dir = self.rc + os.sep + 'update'
        if not os.path.isdir(self.updates_rc_dir):
            os.mkdir(self.updates_rc_dir)
        self.latest = 'tested' # or 'nightly'
        self.master = None
        self.updates_dir = None
        from Support.version import __version__
        self.Version = __version__
        self.baseURL = base_url +'/' + self.Version.split("rc")[0]
        self.date_latest = None
        self.date_tested = None
        self.latestURL = None
        self.testedURL = None
        
    def setUpdatesDir(self):
        """Gets update directory path from '.mgltools/update/'+self.latest.
        If this file doesn't exists, updates_dir is set to one directory up from
        PIL packages directory"""
        if not self.updates_dir:
            latest_rc = self.updates_rc_dir + os.sep + 'tested'    
            if os.path.exists(latest_rc):
                self.updates_dir = open(latest_rc).read()
                if self.updates_dir:
                    self.updates_dir = self.updates_dir.split('\t')[0]
                    self.updates_dir= os.path.split(self.updates_dir)[0]

            if not self.updates_dir:
                latest_rc = self.updates_rc_dir + os.sep + 'nightly'
                if os.path.exists(latest_rc):
                    self.updates_dir = open(latest_rc).read()
                    if self.updates_dir:
                        self.updates_dir = self.updates_dir.split('\t')[0]
                        self.updates_dir= os.path.split(self.updates_dir)[0]
                    
            if not self.updates_dir or not os.path.exists(self.updates_dir):
                path = os.path.split(os.path.split(os.path.abspath(__file__))[0])[0]
                # to find ../lib/python2.5/site-packages/
                self.updates_dir = path + os.sep + 'Updates'
    
    def testWritable(self):
        "Tests if we can open(mgltools.tar.gz,'w') in the Updates directory"
        if not os.path.isdir(self.updates_dir):
            try:
                os.mkdir(self.updates_dir)
            except Exception, inst:
                print inst    
                if self.master:
                    showinfo("Could not create "+self.updates_dir,
                             "Please select directory for downloading updates.")                
                    self.browseUpdatesDir()
                else:
                    self.updates_dir = raw_input("Could not create " +
                                             self.updates_dir +
                         "\nPlease enter directory path for saving updates\n")
                if not self.updates_dir:
                    return False
                else:
                    self.testWritable()
        try:
            open(self.updates_dir + os.sep + 'mgltools.tar.gz','w')
            os.remove(self.updates_dir + os.sep + 'mgltools.tar.gz')
        except Exception, inst:
            if self.master:
                showinfo("Could not create " + self.updates_dir + os.sep + 
                         'mgltools.tar.gz', 
                         "Please select directory for downloading updates.")                
                self.browseUpdatesDir()
            else:
                self.updates_dir = raw_input("Could not create "+
                                             self.updates_dir + os.sep + 
                                             "mgltools.tar.gz\n"
                             "Please enter directory path for saving updates\n")
            if not self.updates_dir:
                return False
            else:
                self.testWritable()    
        return True
    
    def getUpdates(self, packages=[]):
        "Downloads and unpacks updates"
        self.setUpdatesDir()
        if not self.testWritable():
            print "Updates are not installed!"
            return

        if not self.checkRegistration():
            print "Updates are not installed!"
            self.cancel()
            return            
        
        regDict = pickle.load(open(self.rc + os.sep + ".registration"))
        UserID = regDict['UserID']

        update_url = self.baseURL
        if not self.date_latest or not self.date_tested:
            self.getInfo()            
        date_tested = self.date_tested
        date_latest = self.date_latest
        
        if self.latest == 'tested':
            update_url = self.testedURL
            self.clearUpdates()
        else:
            if date_tested == date_latest:
                try:
                     self.clearUpdates()
                except:
                    pass
                update_url = self.testedURL
                self.latest = 'tested'
            else:
                latest_rc = self.updates_rc_dir + os.sep + self.latest
                if os.path.exists(latest_rc):                    
                    lines = open(latest_rc).readlines()
                    if lines:
                        date = lines[0].split('\t')[1]
                        if date < date_latest:
                            self.updates_dir += os.sep + 'nightly_'
                            self.updates_dir += date_latest.replace('/','_')
                            lines.insert(0,self.updates_dir + "\t" + 
                                         date_latest +'\n')
                            open(latest_rc,'w').writelines(lines)
                    else:
                        self.updates_dir += os.sep + 'nightly'
                        self.updates_dir += '_' + date_latest.replace('/','_')
                        open(latest_rc,'w').write(self.updates_dir + "\t" + 
                                                  date_latest + '\n')                          
                else:
                    self.updates_dir += os.sep + 'nightly'
                    self.updates_dir += '_' + date_latest.replace('/','_')
                    open(latest_rc,'w').write(self.updates_dir + "\t" + 
                                              date_latest + '\n')       
                update_url = self.latestURL

        if self.latest == 'tested':
            if not os.path.isdir(self.updates_dir):
                    os.mkdir(self.updates_dir)        
            self.updates_dir += os.sep + 'tested'
            
        if not os.path.isdir(self.updates_dir):
                os.mkdir(self.updates_dir)        

        tar_file = self.update_file
        download_file = open(tar_file,'wb')

        url_file =  urllib.FancyURLopener().open(update_url)
        if self.master:
            upload_size = url_file.headers['content-length']
            per_size = int(upload_size)/99
            self.progressBar.configure(labeltext="Progress ...")
            self.master.update()
            for i in range(101):
                if self.updates_dir == 'cancel':
                    return
                self.progressBar.set(i)
                self.master.update()
                tmp = url_file.read(per_size)
                if tmp:
                    download_file.write(tmp)
                else:
                    download_file.close()
                    url_file.close()
                    break
            self.progressBar.canvas.itemconfig(self.progressBar.progressLabel,
                                               text=" Please wait...")
            self.master.update()
        else:
            print "Downloading updates from\n" + update_url +"\nPlease wait..."
            download_file.write(url_file.read())
            download_file.close()

        if os.name != 'nt': #sys.platform != 'win32':
            tar = tarfile.open(tar_file)
            for tarinfo in tar:
                tar.extract(tarinfo, path=self.updates_dir)
            tar.close()
        else:
            unzipper = unzip.unzip()
            unzipper.extract(tar_file, self.updates_dir)
        os.remove(tar_file)
        if self.latest == 'tested':
            latest_rc = self.updates_rc_dir + os.sep + self.latest
            open(latest_rc,'w').write(self.updates_dir + "\t" + date_tested)       
        headers = {"Content-type": "application/x-www-form-urlencoded",
                                                         "Accept": "text/plain"}
        update_dict = {}
        update_dict['UserID'] = UserID
        update_dict['type'] = self.latest
        if self.latest == 'tested':
            update_dict['date'] = date_tested
        else:
            update_dict['date'] = date_latest
        params = urllib.urlencode(update_dict)
        conn = httplib.HTTPConnection("www.scripps.edu:80")
        conn.request("POST", "/cgi-bin/sanner/update_mgltools_user.py", 
                     params, headers)
        response = conn.getresponse()
        conn.close()
        if self.master:
            self.finishGUI()
            
    def gui(self):
        "GUI for MGLTools updates"
        import Tkinter
        self.master = Tkinter.Tk()
        self.master.lift()
        self.master.title("Update MGLTools")
        self.master.option_add('*font',"Times 12 bold")
        text = "Update Manager"
        Tkinter.Label(self.master, text=text, bg='white', font=(ensureFontCase('helvetica'), 16)
                      ).grid(column=0, row=0, columnspan=3, sticky='snew')
                       
        tested_rc = self.updates_rc_dir + os.sep + "tested"
        nightly_rc = self.updates_rc_dir + os.sep + "nightly"
        tested_dir = None
        nightly_dir = None

        if os.path.exists(tested_rc):
            tested_dir = open(tested_rc).read()
        
        if os.path.exists(nightly_rc):
            nightly_dir = open(nightly_rc).readlines()

        text = "You are running MGLTools " + self.Version
        if nightly_dir:
            tmp = nightly_dir[0]
            text += " with nightly update from " +tmp.split()[1]
            if len(nightly_dir)>1:
                text += "\nOlder installed updates:\n"
                for item in nightly_dir[1:]:
                    text += "    nightly " + item.split('\t')[1] +"\n"
                if tested_dir:
                    text += "    tested  " + tested_dir.split('\t')[1]
            elif tested_dir:
                text += "\nOlder installed updates:\n"
                text += "    tested " + tested_dir.split('\t')[1]
        elif tested_dir:
            text += " with tested update from " + tested_dir.split()[1]

        Tkinter.Label(self.master, text=text, bg='white', justify='left'
                      ).grid(row=1, column=0, columnspan=3, sticky='snew')

        self.notebook = Pmw.NoteBook(self.master)
        self.notebook.grid(row=2, column=0, columnspan=3, sticky='snew')
        webPage = self.notebook.add('Online')
        self.webPage = webPage
        self.updateCheckbutton = Tkinter.Checkbutton(webPage,
                                                     text="Check for Updates", 
                                                     command=self.checkUpdates)
        self.updateCheckbutton.grid(row=0, column=0, columnspan=3,sticky='nw')
        
        self.tk_latest = Tkinter.IntVar()

        from mglutil.gui.BasicWidgets.Tk.progressBar import ProgressBar
        self.frame = Tkinter.Frame(self.master, relief='groove')
        self.frame.grid(row=4, column=0, columnspan=3, sticky='ew')
        self.progressBar = ProgressBar(master=self.frame, labelside=None, 
                                       width=200, height=20, mode='percent')
        self.progressBar.setLabelText('Progress...')
        self.frame.grid_forget()
        self.waitTk = Tkinter.IntVar()
        if tested_dir or nightly_dir:
            Tkinter.Button(self.master, text="Revert to "+self.Version,
                       command=self.clearUpdatesGUI).grid(row=5, column=0)
        if nightly_dir:
            if len(nightly_dir) > 1:
                Tkinter.Button(self.master, text="Rollback to nightly " + \
                               nightly_dir[1].split('\t')[1],
                               command=self.rollback).grid(row=5, column=1)
            elif tested_dir:
                Tkinter.Button(self.master, text="Rollback to tested " + \
                               tested_dir.split('\t')[1],
                               command=self.rollback).grid(row=5, column=1)
                
                       
        Tkinter.Button(self.master, text="Cancel", command=self.cancel).\
                                                           grid(row=5, column=2)

        self.tested_dir = tested_dir
        self.nightly_dir = nightly_dir

        
        filePage = self.notebook.add('File')
        self.localLabel = Tkinter.Label(filePage, text="Update from Local File", 
                                        justify='left', fg='gray45')
        self.localLabel.grid()        
        self.localButton = Tkinter.Button(filePage, text="Browse...", 
                                     command=self.updateLocal)
        self.localButton.grid()
        self.notebook.setnaturalsize()
        return self.waitTk 
    
    def checkUpdates(self):
        """This command is called when Check for Updates button is pressed"""
        self.updateCheckbutton.configure(state='disabled')
        
        #check to see if self.Version is the latest
        version = urllib.urlopen(base_url+'/version').read().strip()
        if self.Version < version:
            txt = "New version (" + version + ") is available for download."
            Tkinter.Label(self.webPage, text=txt).\
                              grid(row=1, column=0, sticky='ew')
            
            l = Tkinter.Label(self.webPage, fg='Blue', cursor='hand1',
                          text='http://mgltools.scripps.edu/downloads')
            l.grid(row=2, column=0, sticky='ew')
            def openurl(evt=None):
                webbrowser.open('http://mgltools.scripps.edu/downloads')
            l.bind(sequence="<Button-1>", func=openurl)
            self.notebook.setnaturalsize()
        
        if not self.date_latest or not self.date_tested:
            self.getInfo()
        date_tested = self.date_tested
        date_latest = self.date_latest

        updatesLabel = Tkinter.Label(self.webPage,  
                                     text='Available updates:',
                                     justify='left')
        updatesLabel.grid(row=3, column=0, sticky='w')
        
        self.testedButton = None
        self.nightlyButton = None
        
        def getTested():
            webbrowser.open_new(self.testedURL)

        def getNightly():
            webbrowser.open_new(self.latestURL)
            
        if self.tested_dir and date_tested:
            if date_tested > self.tested_dir.split()[1]:
                testedGroup = Tkinter.LabelFrame(self.webPage, padx=5, pady=5,
                                                 text = "Updates - Tested Builds")
                testedGroup.grid(row=4, column=0, columnspan=3, sticky='ew')
                Tkinter.Label(testedGroup, text=date_tested +'   ',
                              justify='left').\
                              grid(row=0, column=0, sticky='w')
                self.testedButton = Tkinter.Button(testedGroup, 
                                text='  Install  ', command=self.updateTested)
                self.testedButton.grid(row=0, column=1, sticky='ew')
                self.downloadTested = Tkinter.Button(testedGroup, 
                                                     text='  Save to File  ',
                                                     command=getTested)
                self.downloadTested.grid(row=0, column=2, sticky='ew')
                               
        elif date_tested:
            testedGroup = Tkinter.LabelFrame(self.webPage, padx=5, pady=5,
                                             text = "Updates - Tested Builds")
            testedGroup.grid(row=4, column=0, columnspan=3, sticky='ew')
            Tkinter.Label(testedGroup, text=date_tested+'   ',justify='left').\
                          grid(row=0, column=0, sticky='w')
            self.testedButton = Tkinter.Button(testedGroup, 
                            text='  Install  ', command=self.updateTested)
            self.testedButton.grid(row=0, column=1, sticky='ew')
            self.downloadTested = Tkinter.Button(testedGroup, 
                                                 text='  Save to File  ',
                                                 command=getTested)
            self.downloadTested.grid(row=0, column=2, sticky='ew')
                          
                          
        if self.nightly_dir and date_latest:
            if date_latest > self.nightly_dir[0].split()[1]:
                nightlyGroup = Tkinter.LabelFrame(self.webPage, padx=5, pady=5,
                                             text = "Updates - Nightly Builds")
                nightlyGroup.grid(row=5, column=0, columnspan=3, sticky='ew')

                Tkinter.Label(nightlyGroup, text=date_latest,justify='left').\
                              grid(row=0, column=0, sticky='w')
                self.nightlyButton = Tkinter.Button(nightlyGroup, 
                             text='Install', command=self.updateNightly)
                self.nightlyButton.grid(row=0, column=1, sticky='ew')
                self.downloadNightly = Tkinter.Button(nightlyGroup, 
                                                      text='Save to File', 
                                                      command=getNightly)
                self.downloadNightly.grid(row=0, column=2, sticky='ew')

        elif date_latest:
            if date_latest > date_tested:
                nightlyGroup = Tkinter.LabelFrame(self.webPage, padx=5, pady=5,
                                             text = "Updates - Nightly Builds")
                nightlyGroup.grid(row=5, column=0, columnspan=3, sticky='ew')
    
                Tkinter.Label(nightlyGroup, text=date_latest,justify='left').\
                              grid(row=0, column=0, sticky='w')
                self.nightlyButton = Tkinter.Button(nightlyGroup, 
                             text='Install', command=self.updateNightly)
                self.nightlyButton.grid(row=0, column=1, sticky='ew')
                self.downloadNightly = Tkinter.Button(nightlyGroup, 
                                                      text='Save to File', 
                                                      command=getNightly)
                self.downloadNightly.grid(row=0, column=2, sticky='ew')

                # If test results summary is available for current platform (or at least for Linux),
                # add a button to show it.
                
                self.summary_file = "testsummary"+ sys.platform + ".html"
                if sys.platform == "darwin":
                    if os.uname()[-1] == 'Power Macintosh':
                        prt = "powerpc"
                    else:
                        prt = "i386"
                    self.summary_file = "testsummary"+sys.platform+ prt+".html"
                
                date_summary = self.getDate(self.baseURL+'/latest/', 
                              self.baseURL+'/latest/'+self.summary_file)
                if not date_summary: # no test summary for current platform
                    self.summary_file = "testsummarylinux2.html"
                    date_summary = self.getDate(self.baseURL+'/latest/', 
                              self.baseURL+'/latest/'+self.summary_file )
                #print "date_summary:", date_summary
                #print "self.summary_file:", self.summary_file
                if date_summary == date_latest:
                    # add a button to show the report summary    
                    Tkinter.Label(nightlyGroup,  text='Click to see test summary:',
                        justify='left', fg='gray45').grid(row=1, column=0, 
                                                          columnspan=2, 
                                                          sticky='w' )
                    self.summaryButton = Tkinter.Button(nightlyGroup,
                                                        text="Test Summary",
                                                        command=self.showTestSummary )
                    self.summaryButton.grid(row=1, column = 2, sticky='ew')
                    
        if not self.testedButton and not self.nightlyButton:
            updatesLabel.configure(text="No updates are available at this time.")
        self.notebook.setnaturalsize()
        
    def checkRegistration(self):
        if not self.rc:
            return False        
        regfile = self.rc + os.sep + ".registration"
        if not os.path.exists(regfile):
            from mglutil.splashregister.register import Register_User
            register = Register_User(self.Version)
            while register.master: 
                self.master.update()
            if not os.path.exists(regfile):
                return False
        return True
    
    def updateLocal(self):
        "Updates from local file"
        fileTypes = [('', '*.zip *.gz')]        
        file = askopenfilename(parent=self.master,title='Choose a file',
                              filetypes=fileTypes)
        if not file:
            return
        
        self.setUpdatesDir()
        if not self.testWritable():
            print "Updates are not installed!"
            return
        
        self.localButton.configure(state='disabled')
        self.localLabel.configure(text='Please Wait...')
        self.master.update()
        mtime = os.path.getmtime(file)
        import time
        mdate = time.localtime(mtime)
        date_latest = str(mdate[0]) + '/' + str(mdate[1]) + '/' + str(mdate[2])
        latest_rc = self.updates_rc_dir + os.sep + 'nightly'
        if os.path.exists(latest_rc):                    
            lines = open(latest_rc).readlines()
            if lines:
                self.updates_dir += os.sep + 'nightly_'
                self.updates_dir += date_latest.replace('/','_')
                lines.insert(0,self.updates_dir + "    " + 
                             date_latest +'\n')
                open(latest_rc,'w').writelines(lines)
            else:
                self.updates_dir += os.sep + 'nightly'
                self.updates_dir += '_' + date_latest.replace('/','_')
                open(latest_rc,'w').write(self.updates_dir + "    " + 
                                          date_latest + '\n')                          
        else:
            self.updates_dir += os.sep + 'nightly'
            self.updates_dir += '_' + date_latest.replace('/','_')
            open(latest_rc,'w').write(self.updates_dir + "    " + 
                                      date_latest + '\n')       
        if not os.path.isdir(self.updates_dir):
                os.mkdir(self.updates_dir)        
        if file.endswith('.gz'):
            tar = tarfile.open(file)
            for tarinfo in tar:
                tar.extract(tarinfo, path=self.updates_dir)
            tar.close()
        else:
            unzipper = unzip.unzip()
            unzipper.extract(file, self.updates_dir)
        if self.master:
            self.finishGUI()
    
    def browseUpdatesDir(self):
        self.updates_dir = askdirectory()

    def cancel(self):
        self.updates_dir = 'cancel'
        self.finishGUI()
        
    def finishGUI(self):
        self.waitTk.set(1)
        self.master.destroy()

    def rollback(self):
        nighlty_dir = None
        tested_dir = None
        
        nighlty_rc = self.updates_rc_dir + os.sep + 'nightly'
        if os.path.exists(nighlty_rc):
            nighlty_dir = open(nighlty_rc).readlines()

        tested_rc = self.updates_rc_dir + os.sep + 'tested'
        if os.path.exists(tested_rc):
            tested_dir = open(tested_rc).read()

        if nighlty_dir:
            if len(nighlty_dir) > 1:
                open(nighlty_rc,'w').writelines(nighlty_dir[1:])
            else:
                os.remove(nighlty_rc)
            rm_dirs(nighlty_dir[0].split()[0])
        elif tested_dir:
            self.clearUpdates()
        self.finishGUI()
            
    def updateTested(self):
        self.testedButton.configure(state='disabled')
        self.frame.grid(row=4, column=0, columnspan=3, sticky='ew')
        self.progressBar.setLabelText("Downloading Updates. Please wait...")
        self.master.update()
        self.latest = 'tested'
        self.getUpdates()

    def updateNightly(self):
        self.nightlyButton.configure(state='disabled')
        self.frame.grid(row=4, column=0, columnspan=3, sticky='ew')
        self.progressBar.setLabelText("Downloading Updates. Please wait...")
        self.master.update()
        self.latest = 'nightly'
        self.getUpdates()
            
    def clearUpdates(self):
        "Removes all updates files"
        updates_rc = []
        updates_rc.append(self.updates_rc_dir + os.sep + 'nightly')
        updates_rc.append(self.updates_rc_dir + os.sep + 'tested')
        for update in updates_rc:
            if os.path.exists(update):
                update_dir = open(update).read()
                os.remove(update)
                if update_dir:
                    update_dir = update_dir.split('\t')[0]
                    if update_dir.find('tested') == -1 and update_dir.find('nightly') == -1:
                        print "Refusing to delete " + update_dir
                        print "Please delete " +update_dir +" to clear updates"
                        return
                    if os.path.isdir(update_dir):
                        mgl_dirs = os.listdir(update_dir)
                        for dir in mgl_dirs:
                            dir_path = update_dir + os.sep  + dir
                            if os.path.isdir(dir_path):
                                rm_dirs(dir_path)
                            else:
                                os.remove(dir_path)
                        os.rmdir(update_dir)
                                
    def showTestSummary(self):
        webbrowser.open_new(self.baseURL+'/latest/'+self.summary_file)
           
    def getInfo(self):
        "Gets info about nighly and tested builds"
        nightly_url = self.baseURL + '/latest/'
        f = urllib.urlopen(nightly_url)
        html = f.read()
        start = html.find(self.update_file)
        if start > 0:
            html = html[start:].split()
            self.date_latest = html[10].split('</tt>')[0]
            self.latestURL = html[6].split("\"")[1]
        f.close()

        nightly_url = self.baseURL + '/latest-tested/'
        f = urllib.urlopen(nightly_url)
        html = f.read()
        start = html.find(self.update_file)
        if start > 0:
            html = html[start:].split()
            self.date_tested = html[10].split('</tt>')[0]
            self.testedURL = html[6].split("\"")[1]
        f.close()

    def getDate(self, nightly_url, update_url):
        "This function is used to extract date from remote directory listing"
        f = urllib.urlopen(nightly_url)
        html = f.read()
        start = html.find(update_url)
        if start < 0:
            return None
        html = html[start:].split()
        date = html[10].split('</tt>')[0]
        f.close()
        return date


    def clearUpdatesGUI(self):    
        try:
            self.clearUpdates()
        except:
            pass
        self.finishGUI()
        
if __name__ == "__main__":
    u = Update()
    u.gui()
    u.master.mainloop()
    print u.updates_dir
