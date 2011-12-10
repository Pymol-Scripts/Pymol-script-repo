#$Header: /opt/cvs/python/packages/share1.5/AutoDockTools/WebServices.py,v 1.21 2009/08/12 21:40:05 lclement Exp $
#$Id: WebServices.py,v 1.21 2009/08/12 21:40:05 lclement Exp $
# Author: Sargis Dallakyan (sargis@scripps.edu)

import Tkinter, Pmw, os, httplib, webbrowser, urllib, re
from ViewerFramework.VFCommand import CommandGUI, Command
from Pmv.mvCommand import MVCommand
from autostartCommands import menuText
from mglutil.gui.InputForm.Tk.gui import InputFormDescr
from mglutil.util.packageFilePath import getResourceFolderWithVersion
from tkFileDialog import *
from tkMessageBox import *
from mglutil.gui.BasicWidgets.Tk.progressBar import ProgressBar
from mglutil.web.services.AppService_client import AppServiceLocator, launchJobRequest, \
getOutputsRequest, queryStatusRequest
from mglutil.web.services.AppService_types import ns0
import os

class WebServices(MVCommand):
    def __init__(self):
        MVCommand.__init__(self)
        rc =  getResourceFolderWithVersion() + os.sep + 'ws' + os.sep
        if not os.path.exists(rc):
            os.mkdir(rc)
        self.proxy_gama = rc + 'proxy_gama'
        self.rc_ad = rc + "rc_ad"
        self.login = False
        if hasattr(self, 'vf.GUI.ROOT'):
            self.dpf = Tkinter.StringVar(self.vf.GUI.ROOT)
            self.gpf = Tkinter.StringVar(self.vf.GUI.ROOT)
            self.prev_dir = Tkinter.StringVar(self.vf.GUI.ROOT)
            self.ad_radio = Tkinter.IntVar(self.vf.GUI.ROOT)
        else:
            self.dpf = Tkinter.StringVar()
            self.gpf = Tkinter.StringVar()
            self.prev_dir = Tkinter.StringVar()
            self.ad_radio = Tkinter.IntVar()
        self.current_job = None
            
    def guiCallback(self, event=None):
        mainform = self.showForm('default', modal=0, blocking=1.,
                                  initFunc=self.initForm)
         
    def buildFormDescr(self, formName):
        ifd = InputFormDescr(title = "AutoGrid/AutoDock Web Services")
        #Web Services Login
        ifd.append({'name':"LoginGroup", 'widgetType':Pmw.Group,
                     'container':{'LoginGroup':'w.interior()'},
                     'wcfg':{'tag_text':'Web Services Location'},
                     'gridcfg':{'sticky':'nswe'}
                    })
        ifd.append({'widgetType':Pmw.ComboBox, 'name':'WS_address',
                    'parent':'LoginGroup',
                    'wcfg':{'scrolledlist_items':
                            ('http://ws.nbcr.net/opal2/services',),
                            'listheight':50, 'dropdown':1, 'history':1, 
                            'autoclear':1},
                     'gridcfg':{'sticky':'ew', 'row':0, 'column':0,
                                'columnspan':3}
                     })

#        ifd.append({'widgetType':Tkinter.Label, 'name':'New_User',
#                    'parent':'LoginGroup', 'wcfg':{'text':'   New Users?',
#                                           'fg':'Blue','cursor':'hand1'},
#                    'gridcfg':{'sticky':'w', 'row':1, 'column':0}
#                    })
#        ifd.append({'widgetType':Tkinter.Label, 'name':'UserName_Label',
#                    'parent':'LoginGroup', 'wcfg':{'text':'User Name'},
#                    'gridcfg':{'sticky':'e', 'row':1, 'column':1}
#                    })
#        ifd.append({'widgetType':Tkinter.Entry, 'name':'UserName_Entry',
#                    'parent':'LoginGroup','wcfg':{},
#                    'gridcfg':{'sticky':'ew', 'row':1, 'column':2}
#                    })
#        ifd.append({'widgetType':Tkinter.Label, 'name':'Password_Label',
#                    'parent':'LoginGroup', 'wcfg':{'text':'Password'},
#                    'gridcfg':{'sticky':'e', 'row':2, 'column':1}
#                    })
#        ifd.append({'widgetType':Tkinter.Entry, 'name':'Password_Entry',
#                    'parent':'LoginGroup', 'wcfg':{'show':'*'},
#                    'gridcfg':{'sticky':'ew', 'row':2, 'column':2}
#                    })
#        ifd.append({'widgetType':Tkinter.Label, 'name':'Remember_Label',
#                    'parent':'LoginGroup',
#                    'wcfg':{'text':'Remember User Name and Password'},
#                    'gridcfg':{'sticky':'e', 'row':3, 'column':0,'columnspan':2}
#                    })
#        self.RememberLogin_var = Tkinter.BooleanVar()
#        ifd.append({'widgetType':Tkinter.Checkbutton, 'name':'Remember_Checkbutton',
#                    'parent':'LoginGroup', 'variable':self.RememberLogin_var,
#                    'gridcfg':{'sticky':'w', 'row':3, 'column':2}
#                    })       
        #AutoGrid group
        ifd.append({'name':"AutoGrid", 'widgetType':Pmw.Group,
                    'container':{'AutoGrid':'w.interior()'},
                    'wcfg':{'tag_text':'AutoGrid'},
                    'gridcfg':{'sticky':'nswe'}
                    })
        
        ifd.append({'widgetType':Tkinter.Button, 'name':'Run_autogrid',
                    'parent':'AutoGrid', 
                    'wcfg':{'text':'Run AutoGrid ',
                            'command':self.startAutogrid},
                    'gridcfg':{'sticky':'w', 'row':0, 'column':0}
                    })

        ifd.append( {'name': 'gpf_entry', 'parent':'AutoGrid',
                     'widgetType':Tkinter.Entry, 
                     'wcfg':{'width':30,'textvariable':self.gpf},
                     'gridcfg':{'sticky':'w','row':0,'column':1}
                     })

        ifd.append({'name': 'browse_gpf', 'widgetType': Tkinter.Button,
                    'parent':'AutoGrid', 'text':'Browse',
                    'command':self.browse_gpf,
                    'gridcfg':{'sticky':'w','row':0, 'column':2}
                    })
        #AutoDock group
        ifd.append({'name':"AutoDock", 'widgetType':Pmw.Group,
                    'container':{'AutoDock':'w.interior()'},
                    'wcfg':{'tag_text':'AutoDock'},
                    'gridcfg':{'sticky':'nswe'}
                    })

        ifd.append({'widgetType':Tkinter.Button, 'name':'Run_autodock',
                    'parent':'AutoDock', 
                    'wcfg':{'text':'Run AutoDock',
                            'command':self.startAutodock},
                    'gridcfg':{'sticky':'w', 'row':0, 'column':0}
                    })

        ifd.append( {'name': 'dpf_entry', 'parent':'AutoDock',
                     'widgetType':Tkinter.Entry, 
                     'wcfg':{'width':30,'textvariable':self.dpf},
                     'gridcfg':{'sticky':'w','row':0,'column':1}
                     })

        ifd.append({'name': 'browse_dpf', 'widgetType': Tkinter.Button,
                    'parent':'AutoDock', 'text':'Browse',
                    'command':self.browse_dpf,
                    'gridcfg':{'sticky':'w','row':0, 'column':2}
                    })

        ifd.append({'name': 'ag_local', 'widgetType': Tkinter.Radiobutton,
                    'parent':'AutoDock', 'text':'Use local grids',
                    'tooltip':"This option sends locally stored grid files with Web Services request",
                    'wcfg':{'variable':self.ad_radio,'value':0},
                    'gridcfg':{'sticky':'w','row':1, 'column':0,'columnspan':2}
                    })

#        ifd.append({'name': 'ag_before', 'widgetType': Tkinter.Radiobutton,
#                    'parent':'AutoDock', 'text':'Run AutoGrid first',
#                    'tooltip':"This option runs AutoGrid Web Services and uses resulting map files for AutoDock",
#                    'wcfg':{'variable':self.ad_radio,'value':1,'state':'disabled'},
#                    'gridcfg':{'sticky':'w','row':2, 'column':0,'columnspan':2}
#                    })

        ifd.append({'name': 'use_remote', 'widgetType': Tkinter.Radiobutton,
                    'parent':'AutoDock', 'text':'Use grids from server directory',
                    'tooltip':"This option copies map files from previous AutoGrid run",
                    'wcfg':{'variable':self.ad_radio,'value':2,},
                    'gridcfg':{'sticky':'w','row':3, 'column':0,'columnspan':2}
                    })

        ifd.append( {'name': 'remote_dir', 'parent':'AutoDock',
                     'widgetType':Tkinter.Entry, 
                     'wcfg':{'width':23,'textvariable':self.prev_dir},
                     'gridcfg':{'sticky':'e','row':3,'column':1,'columnspan':2}
                     })

        #Status
        ifd.append({'name':"StatusGroup", 'widgetType':Pmw.Group,
                    'container':{'StatusGroup':'w.interior()'},
                    'wcfg':{'tag_text':'Web Services Status'},
                    'gridcfg':{'sticky':'nswe'}
                    })
        
        ifd.append({'widgetType':Tkinter.Label, 'name':'status0',
                    'parent':'StatusGroup', 
                    'wcfg':{'text':'   ',},
                    'gridcfg':{'sticky':'w', 'row':0, 'column':0}
                    })

        ifd.append({'widgetType':Tkinter.Label, 'name':'status1',
                    'parent':'StatusGroup', 
                    'wcfg':{'text':'   ',},
                    'gridcfg':{'sticky':'w', 'row':1, 'column':0}
                    })

        ifd.append({'name':'WS_ProgressBar', 'widgetType':Tkinter.Frame,
                    'parent':'StatusGroup', 'wcfg':{'height':30},
                     'gridcfg':{'sticky':'ew', 'row':2,'column':0}
                    })
        
        ifd.append({'widgetType':Tkinter.Label, 'name':'down_label',
                    'parent':'StatusGroup', 
                    'wcfg':{'text':'   ',},
                    'gridcfg':{'sticky':'w', 'row':3, 'column':0}
                    })
        
        return ifd
    
    def browse_gpf(self):
        filename = askopenfilename(filetypes=[('Grid Parameter File','*.gpf')],\
                                   title="Please Select Grid Parameter File",
                                   parent=self.cmdForms['default'].root)
        if filename:
            self.gpf.set(filename)
            self.cmdForms['default'].descr.entryByName['Run_autogrid']['widget'].configure(state='normal')
            #self.cmdForms['default'].descr.entryByName['ag_before']['widget'].configure(state='normal')

    def browse_dpf(self):
        filename = askopenfilename(filetypes=[('Dock Parameter File','*.dpf')],\
                                   title="Please Select Dock Parameter File",
                                   parent=self.cmdForms['default'].root)
        if filename:
            self.dpf.set(filename)
            self.cmdForms['default'].descr.entryByName['Run_autodock']['widget'].configure(state='normal')
           
    def initForm(self, cmdForm=None):
        cmdForm.descr.entryByName['WS_address']['widget'].selectitem(0)
#        if not os.path.exists(self.rc_ad):
#            open(self.rc_ad,'w')
#        else:
#            file = open(self.rc_ad)
#            text = file.read()
#            text = text.split()
#            for line in text:
#                tmp_line = line.split('User:')
#                if len(tmp_line) > 1:
#                    cmdForm.descr.entryByName['UserName_Entry']['wcfg']\
#                    ['textvariable'].set(tmp_line[1])
#                tmp_line = line.split('Password:')
#                if len(tmp_line) > 1:
#                    cmdForm.descr.entryByName['Password_Entry']['wcfg']\
#                    ['textvariable'].set(tmp_line[1])
#            file.close()
            
#        def openurl(event):
#            webbrowser.open('https://nbcr.net:8443/worksphere/start?cid=apply')
#        cmdForm.descr.entryByName['New_User']['widget'].bind(sequence="<Button-1>",
#                                                                   func=openurl)
        if hasattr(self.vf,'dpo') and self.vf.dpo.dpf_filename:
            self.dpf.set(self.vf.dpo.dpf_filename)
            cmdForm.descr.entryByName['Run_autodock']['widget'].configure(state='normal')
        else:
            if not self.dpf.get():
                cmdForm.descr.entryByName['Run_autodock']['widget'].configure(state='disabled')

        if hasattr(self.vf,'gpo') and self.vf.gpo.gpf_filename:
            self.gpf.set(self.vf.gpo.gpf_filename)
            cmdForm.descr.entryByName['Run_autogrid']['widget'].configure(state='normal')
            #cmdForm.descr.entryByName['ag_before']['widget'].configure(state='normal')
        else:
            if not self.gpf.get():
                cmdForm.descr.entryByName['Run_autogrid']['widget'].configure(state='disabled')

        self.progressBar = ProgressBar(
                          cmdForm.descr.entryByName['WS_ProgressBar']['widget'], 
                          labelside=None, width=200, height=20, mode='percent')

        self.progressBar.setLabelText('Progress...')
        self.progressBar.set(0)
        cmdForm.descr.entryByName['WS_ProgressBar']['widget'].grid_forget()        

    def startAutogrid(self):
        self.cmdForms['default'].descr.entryByName['Run_autogrid']['widget']\
                                                    .configure(state='disabled')
        gpf_file = self.gpf.get()
        if not os.path.exists(gpf_file):
            self.cmdForms['default'].descr.entryByName['status0']['widget'].\
            configure(text = 'ERROR: gpf file ' + gpf_file + ' does not exist!')
            return
        
        self.host = self.cmdForms['default'].descr.entryByName['WS_address']['widget'].get()
        
#        if not self.login :
#            self.cmdForms['default'].descr.entryByName['status0']['widget'].\
#                configure(text='Connecting to '+ self.host + ". Please wait...")
#            self.vf.GUI.ROOT.update()
#            f = self.validate_login()
#            if f == "Failed":
#                return
        
        self.appLocator = AppServiceLocator()
        self.req = launchJobRequest()
        input_file = os.path.basename(gpf_file)
        options = '-p ' +  input_file + ' -l ' + os.path.splitext(input_file)[0] + '.glg'
        self.req._argList = options
        
        #input_gpf = ns0.InputFileType_Def('inputFile')
        #input_gpf._name = input_file
        gpfFile = open(gpf_file, 'r')
        gpfFileString = gpfFile.read()
        gpfFile.close()
        #input_gpf._contents = gpfFileString
        
        gpfFileString = gpfFileString.split('\n')
        for line in gpfFileString:
            if line[0:9] == 'receptor ':
                pdbqs = line.split()[1]
        
        #input_pdbqs = ns0.InputFileType_Def('inputFile')
        #input_pdbqs._name = pdbqs
        pdbqs = os.path.join(os.path.split(gpf_file)[0],pdbqs)
        #pdbqsFile = open(pdbqs, 'r')
        #pdbqsFileString = pdbqsFile.read()
        #pdbqsFile.close()
        #input_pdbqs._contents = pdbqsFileString
        
        inputFiles = []
        
        #inputFiles.append(input_gpf)
        #inputFiles.append(input_pdbqs)
        inputFiles.append(self.uploadFile(gpf_file))
        inputFiles.append(self.uploadFile(pdbqs))
        
        self.req._inputFile = inputFiles
        
        self.appServicePort = self.appLocator.getAppServicePort(
                                 self.host+'/AutogridOpalService')
        
        resp = self.appServicePort.launchJob(self.req)

        self.JobID = resp._jobID
        self.cmdForms['default'].descr.entryByName['status0']['widget'].\
                                configure(text = 'Running Autogrid Job ID: ' + self.JobID)
        self.vf.GUI.ROOT.update()
        self.vf.GUI.ROOT.after(5, self.checkStatus)
        self.cmdForms['default'].descr.entryByName['Run_autogrid']['widget'].configure(state='normal')
        self.prev_dir.set(self.JobID)
        self.cmdForms['default'].descr.entryByName['use_remote']['widget'].configure(state='normal')
        
    def startAutodock(self):

        self.cmdForms['default'].descr.entryByName['Run_autodock']['widget']\
                                                    .configure(state='disabled')
        dpf_file = self.dpf.get()
        if not os.path.exists(dpf_file):
            self.cmdForms['default'].descr.entryByName['status0']['widget'].\
            configure(text = 'ERROR: dpf file ' + fpf_file + ' does not exist!')
            return
        
        self.host = self.cmdForms['default'].descr.entryByName['WS_address']['widget'].get()
        
#        if not self.login :
#            self.cmdForms['default'].descr.entryByName['status0']['widget'].\
#                configure(text='Connecting to '+ self.host + ". Please wait...")
#            self.vf.GUI.ROOT.update()
#            f = self.validate_login()
#            if f == "Failed":
#                return

        self.appLocator = AppServiceLocator()
        self.req = launchJobRequest()
        input_file = os.path.basename(dpf_file)
        options = '-p ' +  input_file + ' -l ' + os.path.splitext(input_file)[0] + '.dlg'
        self.req._argList = options
        
        #input_dpf = ns0.InputFileType_Def('inputFile')
        #input_dpf._name = input_file
        dpfFile = open(dpf_file, 'r')
        dpfFileString = dpfFile.read()
        dpfFile.close()
        #input_dpf._contents = dpfFileString
        #DPF file
        inputFiles = []       
        inputFiles.append(self.uploadFile(dpf_file))
        
        
        run_option = self.ad_radio.get()
        if run_option == 0: # sends locally stored grid files
            inputs = re.findall("\w*.\w*\.map ",dpfFileString)
            inputs.extend(re.findall("\w*\.maps.fld",dpfFileString))
            inputs.extend(re.findall("\w*.pdbq[t]*",dpfFileString))
            for input in inputs:
                input = input.strip()
                #ws_input = ns0.InputFileType_Def('inputFile')
                #ws_input._name = input
                input_full_name = os.path.join(os.path.split(dpf_file)[0],input)
                #inputFile = open(input_full_name, 'r')
                #inputFileString = inputFile.read()
                #inputFile.close()
                #ws_input._contents = inputFileString
                inputFiles.append(self.uploadFile(input_full_name))
                
        elif run_option == 2: # runs AutoGrid first
            prev_dir = self.prev_dir.get()
            inputs = re.findall("\w*.\w*\.map ",dpfFileString)
            inputs.extend(re.findall("\w*\.maps.fld",dpfFileString))
            
            host = 'http://'+self.host.split('/')[2]
            for input in inputs:
                self.req._argList += " " +host+"/"+prev_dir+"/"+input
            pdbq_input = re.findall("\w*.pdbq[t]*",dpfFileString)
            pdbq_input = pdbq_input[0].strip()

            #ws_input = ns0.InputFileType_Def('inputFile')
            #ws_input._name = pdbq_input
            input_full_name = os.path.join(os.path.split(dpf_file)[0],pdbq_input)
            #inputFile = open(input_full_name, 'r')
            #inputFileString = inputFile.read()
            #inputFile.close()
            #ws_input._contents = inputFileString
            inputFiles.append(self.uploadFile(input_full_name))
        self.req._inputFile = inputFiles
        self.vf.GUI.ROOT.update()
        self.appServicePort = self.appLocator.getAppServicePort(
                                 self.host+'/AutodockOpalService',)
        resp = self.appServicePort.launchJob(self.req)

        self.JobID = resp._jobID
        self.cmdForms['default'].descr.entryByName['status0']['widget'].\
                                configure(text = 'Running Autodock Job ID: ' + self.JobID)
        self.vf.GUI.ROOT.update()
        self.vf.GUI.ROOT.after(5, self.checkStatus)
        self.cmdForms['default'].descr.entryByName['Run_autodock']['widget'].configure(state='normal')

    def uploadFile(self, path):
        """
        this function given a string containing a path creates a 
        InputFileType to be used with jobLaunch
        """
        inputFile = ns0.InputFileType_Def('inputFile')
        inputFile._name = os.path.basename(path)
        if self.isOpal2():
            #use attachment this is opal2 server
            inputFile._attachment = open(path, "r")
        else:
            #it's not a opal2 server don't user attachment
            infile = open(path, "r")
            inputFile._contents = infile.read()
            infile.close()
        return inputFile

        
    def isOpal2(self):
        """return True if we are using Opal2"""
        print "self.host is: " + self.host
        if self.host.find("/opal2/") != -1:
            return True
        else:
            return False

    def checkStatus(self):
        resp = self.appServicePort.queryStatus(queryStatusRequest(self.JobID))
        if resp._code == 8: # 8 = GramJob.STATUS_DONE
            descr = self.cmdForms['default'].descr
            descr.entryByName['status0']['widget'].configure(text=resp._message)       
            webbrowser.open(resp._baseURL)
            descr.entryByName['status1']['widget'].configure(text=resp._baseURL,
                                                       fg='Blue',cursor='hand1')
            def openurl(event):
                webbrowser.open(resp._baseURL)
            descr.entryByName['status1']['widget'].bind(sequence="<Button-1>",
                                                                   func=openurl)
            self.resp = self.appServicePort.getOutputs(getOutputsRequest(self.JobID))
            descr.entryByName['WS_ProgressBar']['widget'].grid(sticky='ew', 
                                                                row=2, column=0)

            self.opener = urllib.FancyURLopener(cert_file=self.proxy_gama, 
                                                       key_file=self.proxy_gama)
            self.download_finished = False
            self.new_download = True
            self.file_counter = -1
            inputs = [x for x in self.resp._outputFile if x._name[-3:] !='dlg']
            if len(inputs) != len(self.resp._outputFile):
                for input in inputs:
                    self.resp._outputFile.remove(input)
            self.download()                
            return
        else:
            self.cmdForms['default'].descr.entryByName['status0']['widget'].\
                                    configure(text = "Status: " + resp._message)       
            self.cmdForms['default'].descr.entryByName['status1']['widget'].\
                                    configure(text = "")
            
        self.vf.GUI.ROOT.after(5000, self.checkStatus)

    def download(self):
        if self.new_download:
            self.file_counter += 1
            if self.file_counter > self.resp._outputFile.__len__() - 1 :
                self.cmdForms['default'].descr.entryByName['WS_ProgressBar']\
                                                        ['widget'].grid_forget()
                self.cmdForms['default'].descr.entryByName['down_label']\
                                                 ['widget'].configure(text = "")
                self.cmdForms['default'].descr.entryByName['Run_autogrid']\
                                            ['widget'].configure(state='normal')
                return
            
            self.progressBar.configure(progressformat='percent',
                                        labeltext='Progress ... ', max =100)                   
            self.progressBar.set(0)
            remote_file = self.resp._outputFile[self.file_counter]
            
            self.cmdForms['default'].descr.entryByName['down_label']['widget'].\
                configure(text = "Downloading " + remote_file._name + "      " +
                                              str(self.file_counter+1) +" of " + 
                                           str(self.resp._outputFile.__len__()))

            self._url = self.opener.open(remote_file._url)
            try:
                self._out = open(remote_file._name,"w")
            except IOError:
                showerror("Download Failed!", 
                          "Permission denied: " +os.path.join(os.getcwd(),remote_file._name),
                          parent = self.cmdForms['default'].root)
                return             
            bytes = int(self._url.headers.dict['content-length'])
            self._progress_counter = 0
            self._download_bytes = bytes/100
            if self._download_bytes == 0: self._download_bytes = 1
            self.new_download = False
            self.vf.GUI.ROOT.after(1, self.download)
            return
        else:
            self._progress_counter += 1
            if self._progress_counter >  100:
                self._progress_counter =  100
            self.progressBar.set(self._progress_counter)
            tmp = self._url.read(self._download_bytes)
            if tmp:
                self._out.write(tmp)
            else:
                self._url.close()
                self._out.close()
                self.new_download = True
                
            self.vf.GUI.ROOT.after(50, self.download)
        
    def validate_login(self):
        self.login = False
        from mglutil.web.services.SecuritymyproxyloginImplService_services import \
                                          loginUserMyProxyRequestWrapper, \
                                      SecuritymyproxyloginImplServiceLocator
        gamaLoginLocator = SecuritymyproxyloginImplServiceLocator()
        gamaLoginService = gamaLoginLocator.getSecuritymyproxyloginImpl(
                                    ssl=1,transport=httplib.HTTPSConnection)
        req = loginUserMyProxyRequestWrapper()
        username =  self.cmdForms['default'].descr.\
                               entryByName['UserName_Entry']['widget'].get()
        passwd =  self.cmdForms['default'].descr.\
                               entryByName['Password_Entry']['widget'].get()
        if not username or not passwd:
            showerror("Username or Password is missing", 
                      "Login failed. Please type your User Name and Password,\
 or click on New User?", parent = self.cmdForms['default'].root)
            return "Failed"
        req._username = username
        req._passwd = passwd
        resp = gamaLoginService.loginUserMyProxy(req)
        f = open(self.proxy_gama, "w")
        f.write(resp._loginUserMyProxyReturn)
        f.close()
        if self.RememberLogin_var.get():
            file = open(self.rc_ad,'w')
            user = self.cmdForms['default'].descr.entryByName\
                                          ['UserName_Entry']['widget'].get()
            passwd = self.cmdForms['default'].descr.entryByName\
                                          ['Password_Entry']['widget'].get()
            file.write("User:%s\nPassword:%s\n"%(user,passwd))
        self.login = True
        
WebServicesGUI=CommandGUI()
WebServicesGUI.addMenuCommand('AutoToolsBar', menuText['StartMB'], "Web Services...")

commandList = [{'name':'ADweb_services','cmd':WebServices(),'gui':WebServicesGUI}]  

WebServices4GUI=CommandGUI()
WebServices4GUI.addMenuCommand('AutoTools4Bar', menuText['StartMB'], "Web Services...")


def initModule(viewer):
    if not hasattr(viewer, 'ADweb_services') and hasattr(viewer, 'GUI')\
        and hasattr(viewer.GUI, 'currentADTBar'):
        viewer.addCommand(WebServices(),'ADweb_services',WebServices4GUI)
    #else:
    #    for _dict in commandList:
    #        viewer.addCommand(_dict['cmd'],_dict['name'],_dict['gui'])
