# This module handles User Registration
# $Author: sargis $
# $Header: /opt/cvs/python/packages/share1.5/mglutil/splashregister/register.py,v 1.21.2.2 2011/06/06 17:19:38 sargis Exp $
# $Date: 2011/06/06 17:19:38 $
# $Id: register.py,v 1.21.2.2 2011/06/06 17:19:38 sargis Exp $
import Tkinter, tkMessageBox, Pmw
import os, sys, socket, httplib, urllib, pickle, time, shutil
from mglutil.util.packageFilePath import getResourceFolderWithVersion, getResourceFolder

class Register_User:
    """Opens TopLevel Dialog for User Registration"""
    def __init__(self, version = None):
        self.version = version
        master = Tkinter.Toplevel()
        self.master = master
        font = self.master.option_get('font', '*')
        self.master.option_add('*font',"Times 12 bold")
        x = self.master.winfo_screenwidth()/6
        y = self.master.winfo_screenheight()/6        
        geometry = '+%d+%d' % ( x, y)
        self.master.geometry(geometry)
        #self.master.config(font="Helvetica 10 bold italic")
        master.title("MGLTools Registration Form")
        Tkinter.Label(master, text="MGLTools Registration Form").\
                                                        grid(row=0,columnspan=5)
        text = """This data will be kept confidential and will not be made available to any third party."""
        Tkinter.Label(master, text=text).grid(row=1,columnspan=5)
        
        Tkinter.Label(master, text="Fields in Red are Required",fg="Red" ).\
                                                        grid(row=2,columnspan=5)            
        
        Tkinter.Label(master, text="First Name", fg='Red').grid(row=3,sticky='W')
        Tkinter.Label(master, text="Last Name", fg='Red').grid(row=4,sticky='W')
        Tkinter.Label(master, text="Email", fg='Red').grid(row=5,sticky='W')
        Tkinter.Label(master, text="Institution", fg='Red').grid(row=6,sticky='W')
        Tkinter.Label(master, text="Institution Type", fg='Red').grid(row=7,sticky='W')
        Tkinter.Label(master, text="Position").grid(row=8,sticky='W')    
        Tkinter.Label(master, text="Department", ).grid(row=9,sticky='W')
    
        Tkinter.Label(master, text="      ").grid(row=6,column=2) #placeholder
        Tkinter.Label(master, text="Address").grid(row=3,column=3,sticky='W')
        Tkinter.Label(master, text="City" ).grid(row=4,column=3,sticky='W')
        Tkinter.Label(master, text="State" ).grid(row=5,column=3,sticky='W')
        Tkinter.Label(master, text="PostalCode").grid(row=6,column=3,sticky='W')
        Tkinter.Label(master, text="Country").grid(row=7,column=3,sticky='W')
        Tkinter.Label(master, text="Phone").grid(row=8,column=3,sticky='W')
        Tkinter.Label(master, text="Fax").grid(row=9,column=3,sticky='W')
        
        self.e_First_Name = Tkinter.Entry(master,bg="White")
        self.e_Last_Name = Tkinter.Entry(master,bg="White")
        self.e_Email = Tkinter.Entry(master,bg="White")
        self.e_Institution = Tkinter.Entry(master,bg="White")
        self.e_Institution_Type = Pmw.OptionMenu(master, items = 
                  ['Academic','Government','Commercial'], menubutton_width = 16)
        self.e_Position = Tkinter.Entry(master,bg="White")
        self.e_Department = Tkinter.Entry(master,bg="White")
        self.e_Address = Tkinter.Entry(master,bg="White")
        self.e_City = Tkinter.Entry(master,bg="White")
        self.e_State = Tkinter.Entry(master,bg="White")
        self.e_PostalCode = Tkinter.Entry(master,bg="White")
        self.e_Country = Tkinter.Entry(master,bg="White")
        self.e_Phone = Tkinter.Entry(master,bg="White")
        self.e_Fax = Tkinter.Entry(master,bg="White")
        self.e_First_Name.grid(row=3, column=1)
        self.e_Last_Name.grid(row=4, column=1)
        self.e_Email.grid(row=5, column=1)
        self.e_Institution.grid(row=6, column=1)
        self.e_Institution_Type.grid(row=7, column=1)
        self.e_Position.grid(row=8, column=1)
        self.e_Department.grid(row=9, column=1)
        self.e_Address.grid(row=3, column=4)
        self.e_City.grid(row=4, column=4)
        self.e_State.grid(row=5, column=4)
        self.e_PostalCode.grid(row=6, column=4)            
        self.e_Country.grid(row=7, column=4)            
        self.e_Phone.grid(row=8, column=4)
        self.e_Fax.grid(row=9, column=4)

        mgl_group = Tkinter.LabelFrame(master, text="Which of the "+ 
"following program(s) are you planning to use?",labelanchor='n')
        mgl_group.grid(row=11,columnspan=8,sticky='WESN',pady=5) 
        label_message = Tkinter.Label(mgl_group, text="Check all that apply")
        label_message.grid(row=12,columnspan=5) 
        self.adt_var = Tkinter.IntVar()
        checkbutton = Tkinter.Checkbutton(mgl_group, 
               text="AutoDockTools      ", variable=self.adt_var)
        checkbutton.grid(row=13, column=0)
        
        self.pmv_var = Tkinter.IntVar()
        checkbutton = Tkinter.Checkbutton(mgl_group, 
                               text="PMV       ", variable=self.pmv_var)
        checkbutton.grid(row=13,column=2)

        self.vision_var = Tkinter.IntVar()
        checkbutton = Tkinter.Checkbutton(mgl_group, text="Vision      ", 
                                           variable=self.vision_var)
        checkbutton.grid(row=13,column=4)

        self.DejaVu_var = Tkinter.IntVar()
        checkbutton = Tkinter.Checkbutton(mgl_group, text="DejaVu", 
                                           variable=self.DejaVu_var)
        checkbutton.grid(row=13,column=5)

        bin_source = Tkinter.LabelFrame(master, text="Did you install MGLTools "+ 
"from Binary and/or Source?",labelanchor='n')
        bin_source.grid(row=14,columnspan=8,sticky='WESN',pady=5) 
        
        self.bin_var = Tkinter.IntVar()
        checkbutton = Tkinter.Checkbutton(bin_source, 
               text="Binary     ", variable=self.bin_var)
        checkbutton.grid(row=15, column=0)
        
        self.source_var = Tkinter.IntVar()
        checkbutton = Tkinter.Checkbutton(bin_source, 
                               text="Source       ", variable=self.source_var)
        checkbutton.grid(row=15,column=2)
        self.label_message = Tkinter.Label(master, text="")
        self.label_message.grid(row=16,columnspan=5) #placeholder
        
        self.regButton =  Tkinter.Button(master, text='Register',
                                                          command=self.Register)
        self.regButton.grid(row=17,column=1,columnspan=2)
    
        Tkinter.Button(master, text='Cancel',
                       command=self.Cancel).grid(row=17, column=3, columnspan=2,
                                                                         pady=2)
        self.master.option_add('*font', font) 
        self.preFill()
        
    def Register(self):
        self.regButton.configure(state='disabled')
        form_dict = {}
        First_Name  = self.e_First_Name.get()
        First_Name = First_Name.strip()
        if not First_Name:
            self.label_message.configure(text='First Name is missing', fg='Red')
            self.e_First_Name.configure(bg='Red')
            self.regButton.configure(state='normal')
            return
        else:
            self.e_First_Name.configure(bg='White')
            self.label_message.configure(text = '')
        form_dict['First_Name'] = First_Name
        Last_Name  = self.e_Last_Name.get()
        Last_Name = Last_Name.strip()
        if not Last_Name:
            self.label_message.configure(text='Last Name is missing', fg='Red')
            self.e_Last_Name.configure(bg='Red')
            self.regButton.configure(state='normal')
            return
        else:
            self.e_Last_Name.configure(bg='White')
            self.label_message.configure(text='')
        form_dict['Last_Name'] = Last_Name
        Email  = self.e_Email.get()
        Email = Email.strip()
        if not Email or Email.find('@') == -1:
            self.label_message.configure(text='Please provide a valid Email address', fg = 'Red')
            self.e_Email.configure(bg='Red')
            self.regButton.configure(state='normal')
            return
        else:
            self.e_Email.configure(bg='White')
            self.label_message.configure(text='')
        form_dict['Email'] = Email
        Institution  = self.e_Institution.get()
        Institution = Institution.strip()
        if not Institution:
            self.label_message.configure(text = 'Institution is missing', fg = 'Red')
            self.e_Institution.configure(bg='Red')
            self.regButton.configure(state='normal')            
            return
        else:
            self.e_Institution.configure(bg='White')
            self.label_message.configure(text = '')
        form_dict['Institution'] = Institution
        Institution_Type  = self.e_Institution_Type.getvalue()
        if Institution_Type not in ['Academic','Government','Commercial']:
            self.label_message.configure(text = 'Institution Type should be\
            Academic,Government or Commercial', fg = 'Red')
            self.regButton.configure(state='normal')            
            return
        else:
            self.label_message.configure(text = '')
        form_dict['Institution_Type'] = Institution_Type
        Position  = self.e_Position.get()
        Position = Position.strip()
        if not Position:
            Position = ' '
        form_dict['Position'] = Position
        Department = self.e_Department.get()
        Department = Department.strip()
        form_dict['Department'] = Department
        Address = self.e_Address.get()
        Address = Address.strip()
        form_dict['Address'] = Address
        City = self.e_City.get()
        City = City.strip()
        form_dict['City'] = City
        State = self.e_State.get()
        State = State.strip()
        form_dict['State'] = State
        PostalCode = self.e_PostalCode.get()
        PostalCode = PostalCode.strip()
        form_dict['PostalCode'] = PostalCode
        Country = self.e_Country.get()
        Country = Country.strip()
        form_dict['Country'] = Country
        Phone = self.e_Phone.get()
        Phone = Phone.strip()
        form_dict['Phone'] = Phone
        Fax = self.e_Fax.get()
        Fax = Fax.strip()
        form_dict['Fax'] = Fax
        planning_to_use = ''
        if self.adt_var.get():
            planning_to_use = 'ADT'
        if self.pmv_var.get():
            planning_to_use += ',PMV'
        if self.vision_var.get():
            planning_to_use += ',Vision'
        if self.DejaVu_var.get():
            planning_to_use += ',DejaVu'
            
        form_dict['PlanningToUse'] = planning_to_use
        build_from = ''
        if self.bin_var.get():
            build_from = 'Binary'
        if self.source_var.get():
            build_from += ',Source'
        form_dict['BuildFrom'] = build_from
        self.label_message.configure(text = 'Submitting the Registration Form')
        form_dict['version'] = self.version.split('(')[0]
        os_name = os.name
        form_dict['os_name'] = os_name
        sys_platfrom = sys.platform
        form_dict['sys_platfrom'] = sys_platfrom
        sys_version = sys.version
        form_dict['sys_version'] = sys_version.replace('\n','')
        try:
            hostname = gethostname()
            form_dict['hostname'] = hostname
        except:
            form_dict['hostname'] = "problem"
        params = urllib.urlencode(form_dict)
        self.label_message.configure(text = 'Please Wait', fg = 'Red')
        headers = {"Content-type": "application/x-www-form-urlencoded",
                                                         "Accept": "text/plain"}
        conn = httplib.HTTPConnection("www.scripps.edu:80")
        try:
            conn.request("POST", "/cgi-bin/sanner/register_mgltools.py", params, headers)
            response = conn.getresponse()
        except Exception, inst:
            from traceback import print_exc
            print_exc()
            return            

        self.master.update()
        if response.status == 200:
            getResourceFolder
            reg_file =  os.path.join(getResourceFolder(),  '.registration')
            UserID = response.read()
            if UserID:
                form_dict['UserID'] = UserID
                file = open(reg_file,'w')
                pickle.dump(form_dict, file)
                file.close()
                c_reg_file = os.path.join(getResourceFolderWithVersion(),  '.registration')
                shutil.copy(reg_file, c_reg_file)
            else:
                tkMessageBox.showerror("ERROR", "Registration failed to create User." +
                                       "\nPlease contact mgltools@scripps.edu")
                self.regButton.configure(state='normal')
                return
        else:
            tkMessageBox.showerror("ERROR", "Unable to connect to Registration Database" +
            "\nPlease try again")
            self.regButton.configure(state='normal')
            return
        conn.close()        
        self.Cancel()
        tkMessageBox.showinfo("Thank You", "Thank you for registering MGLTools!")

    def Cancel(self):
        self.master.destroy()
        self.master = None

    def preFill(self):
        old_rc = getResourceFolder()
        regfile = os.path.join(old_rc, ".registration")
        if os.path.exists(regfile):
            form_dict =  pickle.load(open(regfile, 'rb'))
            if form_dict.has_key("First_Name"):
                self.e_First_Name.insert(0, form_dict['First_Name'])
            if form_dict.has_key("Last_Name"):
                self.e_Last_Name.insert(0, form_dict['Last_Name'])
            if form_dict.has_key("Email"):
                self.e_Email.insert(0, form_dict['Email'])
            if form_dict.has_key("Institution"):
                self.e_Institution.insert(0, form_dict['Institution'])
            if form_dict.has_key("Institution_Type"):
                self.e_Institution_Type.setvalue(form_dict['Institution_Type'])
            if form_dict.has_key("Position"):
                self.e_Position.insert(0, form_dict['Position'])
            if form_dict.has_key("Department"):
                self.e_Department.insert(0, form_dict['Department'])
            if form_dict.has_key("Address"):
                self.e_Address.insert(0, form_dict['Address'])
            if form_dict.has_key("City"):
                self.e_City.insert(0, form_dict['City'])
            if form_dict.has_key("State"):
                self.e_State.insert(0, form_dict['State'])
            if form_dict.has_key("PostalCode"):
                self.e_PostalCode.insert(0, form_dict['PostalCode'])
            if form_dict.has_key("Country"):
                self.e_Country.insert(0, form_dict['Country'])
            if form_dict.has_key("Phone"):
                self.e_Phone.insert(0, form_dict['Phone'])
            if form_dict.has_key("Fax"):
                self.e_Fax.insert(0, form_dict['Fax'])
            
        
    def Update(self):
        regfile = None
        old_rc = getResourceFolder()
        rcWithVersion = getResourceFolderWithVersion()
        regfile = os.path.join(old_rc, ".registration")
        if os.path.exists(old_rc + os.sep + ".registration"):
            regfile = old_rc + os.sep + ".registration"
        else:
            dirlist = os.listdir(old_rc)
            for item in dirlist:
                tmpRegFile = old_rc+os.sep+item+os.sep + ".registration"
                if os.path.exists(tmpRegFile):
                    regfile = tmpRegFile
                    break
        
        regDict = pickle.load(open(regfile))
        regDict['Version'] = Version
        form_dict = {}
        form_dict['UserID'] = regDict['UserID']
        form_dict['Version'] = regDict['Version']
        import httplib, urllib
        params = urllib.urlencode(form_dict)
        headers = {"Content-type": "application/x-www-form-urlencoded",
                                                 "Accept": "text/plain"}
        conn = httplib.HTTPConnection("www.scripps.edu:80")
        conn.request("POST", "/cgi-bin/sanner/update_mgltools_version.py", 
                     params, headers)
        response = conn.getresponse()
        if response.status == 200:
            reg = open(lRessourceFolder + os.sep + ".registration", 'w')
            pickle.dump(regDict,reg)
            reg.close()
        conn.close()

        
def getdnsnames(name):
    d = socket.gethostbyaddr(name)
    names = [ d[0] ] + d[1] + d[2]
    return names

def resolve(name):
    names = getdnsnames(name)
    for dnsname in names:
        if '.' in dnsname:
            fullname = dnsname
            break
    else:
        fullname = name
    return fullname

def gethostname():
    fullname = socket.gethostname()
    if '.' not in fullname:
        fullname = resolve(fullname)
    return fullname 

def Update_User(registration):
    """To be implemented """
    pass

if __name__ == '__main__':
    m = Tkinter.Tk()
    Register_User('sd')
    m.mainloop()
