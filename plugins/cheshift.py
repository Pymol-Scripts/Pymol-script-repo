'''
cheshift plugin-in: Validate your protein model with PyMOL
Described at PyMOL wiki: http://www.pymolwiki.org/index.php/cheshift

Author : Osvaldo Martin
email: aloctavodia@gmail.com
Date: August 2013
License: GNU General Public License
Version 3.0
'''

import Tkinter
from Tkinter import *
import tkFileDialog
import Pmw
from pymol import cmd
import tempfile
import time
import re
import os
try:
    import mechanize
except:
    print '<'*80 + '\n\nCheShift plug-in needs the mechanize library, please follows the instructions at\nhttp://www.pymolwiki.org/index.php/Pycheshift to install this library.\n\n' + '>'*80


try:
    import gconf
except:
    pass


def __init__(self):
    """Add this Plugin to the PyMOL menu"""
    self.menuBar.addmenuitem('Plugin', 'command',
                            'Cheshift',
                            label = 'Cheshift',
                            command = lambda : mainDialog())

def cheshift(pdb_path,cs_path):
    """Conects to CheShift to send and receive data"""
    pdb_filename = pdb_path.split('/')[-1]
    pdb_filenamenoext = re.sub(r'[^A-Za-z0-9-.]', '', pdb_filename.split('.')[0])
    cs_filename = cs_path.split('/')[-1]
# Send data to CheShift.com and captures the name of the new web page.
    br = mechanize.Browser()
    br.set_handle_refresh(False)
    br.set_proxies(proxy_values)
    try:
        br.open('http://www.cheshift.com/visual')
        br.select_form(nr=0)
        br.form.add_file(open(pdb_path), "chemical/x-pdb", pdb_filename,
        name = 'uploaded')
        br.form.add_file(open(cs_path), "text/plain", cs_filename, name='CS_file')
        br.form.set_value([ref_value.get()],name='radiogroup')
        br.form.set_all_readonly(False)
        response = br.submit().geturl()
        subfix = response.split('/')[-1].split('.')[0]
    except:
        Pmw.MessageDialog(title='Error',message_text=('Something went wrong while processing your request\n Please check your files and try again\n\nIf you could not find the problem please send and e-mail to aloctavodia@gmail.com'))

# Checks if the colored PDB file was created or if something went wrong
    test = 0
    while True:
        try:# If the colored PDB was created, creates a local temporal PDB file
# and loads that file into PyMOL.
            result0 = br.open('http://www.cheshift.com/jobs/%s/%s_Ca.pdb?r=%d' % (subfix, pdb_filenamenoext, test)).read()
            result1 = br.open('http://www.cheshift.com/jobs/%s/%s_Cb.pdb?r=%d' % (subfix, pdb_filenamenoext, test)).read()
            result2 = br.open('http://www.cheshift.com/jobs/%s/%s_CaCb.pdb?r=%d' % (subfix, pdb_filenamenoext, test)).read()
            cmd.reinitialize()
            fd = tempfile.NamedTemporaryFile(bufsize=0,delete=False)
            fd.write(result0)
            fd.close()
            pdb_tmp = fd.name
            rename = ('%s_Ca.pdb' % (pdb_filenamenoext))
            cmd.load(pdb_tmp, rename)
            fd = tempfile.NamedTemporaryFile(bufsize=0,delete=False)
            fd.write(result1)
            fd.close()
            pdb_tmp = fd.name
            rename = ('%s_Cb.pdb' % (pdb_filenamenoext))
            cmd.load(pdb_tmp, rename)
            fd = tempfile.NamedTemporaryFile(bufsize=0,delete=False)
            fd.write(result2)
            fd.close()
            pdb_tmp = fd.name
            rename = ('%s_CaCb.pdb' % (pdb_filenamenoext))
            cmd.load(pdb_tmp, rename)
            colorize()
            result = br.open('http://www.cheshift.com/jobs/%s/%s.zip?r=%d' % (subfix, pdb_filenamenoext, test)).read()
            details_path = pdb_path.split('.')[0]
            fd = open('%s.zip' % details_path, 'w')
            fd.write(result)
            fd.close()
            print '<'*80 + '\nCheShift-2 Validation Report saved at\n%s.zip\n' % details_path + '>'*80
            break
        except:
            test += 1
            try:# check for posible errors
                result = br.open('http://www.cheshift.com/jobs/%s/error.txt?r=%d' % (subfix, test)).read().split() # behind a proxy you need to add a "dummy" modifier to the requested file (?r=%d'). Otherwise you will just be checking the cached file and not the actual file.
                if int(result[0]) == 0:
                    Pmw.MessageDialog(title = 'Error',message_text=('The residue at position %s is missing from your PDB file\n Please fix the problem and try again' % (result[1])))
                elif int(result[0]) == 1:
                    Pmw.MessageDialog(title = 'Error',message_text = ('The residue %s%s in your PDB file does not match\n with residue %s%s in your chemical shift file.\n\n Please check your files and try again' % (result[1], result[3], result[2], result[3])))
                else:
                    Pmw.MessageDialog(title = 'Error',message_text=('Something went wrong while processing your request\n Please check your files and try again\n\nIf you could not find the problem please write to us'))
                break
            except:
                pass
        time.sleep(3)

def mainDialog():
    """ Creates the GUI """
    global ref_value, conection, httpProxy, port, username, password, entry_proxy, entry_port, check_authentication, entry_username, entry_password, authentication
    master = Tk()
    master.title(' CheShift ')
    w = Tkinter.Label(master, text='\nCheShift: Validate your protein model with PyMOl\n - www.cheshift.com -\n',
                                background = 'black',
                                foreground = 'white')
    w.pack(expand=1, fill = 'both', padx=4, pady=4)
############################ NoteBook #########################################
    Pmw.initialise()
    nb = Pmw.NoteBook(master, hull_width=430, hull_height=250)
    p1 = nb.add('Run CheShift')
    p2 = nb.add(' Color code ')
    p3 = nb.add('    Proxy   ')
    p4 = nb.add('    About   ')
    nb.pack(padx=5, pady=5, fill=BOTH, expand=1)
############################ RUN CheShift TAB #################################
# select files
    group = Pmw.Group(p1,tag_text='Select your files')
    group.pack(fill='both', expand=1, padx=5, pady=5)
    Button(group.interior(), text='        PDB file        ', command=retrieve_pdb).pack()
    Button(group.interior(), text='Chemical Shift file', command=retrieve_cs).pack()
# reference compound
    group = Pmw.Group(p1,tag_text='Select a reference compound')
    group.pack(fill='both', expand=1, padx=5, pady=5)
    ref_value = StringVar(master=group.interior())
    ref_value.set('1.70')
    Radiobutton(group.interior(), text="DSS", variable=ref_value, value='1.70').pack(side=LEFT)
    Radiobutton(group.interior(), text="TSP", variable=ref_value, value='1.82').pack(side=LEFT)
    Radiobutton(group.interior(), text="TMS", variable=ref_value, value='0.00').pack(side=LEFT)
# submit
    Button(p1, text="Submit", command=submit).pack(side=BOTTOM)
############################ COLOR TAB ########################################
    Label(p2, text =u"""
Colors indicate the difference between predicted and
observed 13C\u03B1 and 13C\u03B2 chemical shifts values
averaged over all uploaded conformers.
Green, yellow and red colors represent small, medium
and large differences, respectively. 
White is used if either the prediction fail or the
observed value is missing
CheShift-2 provied alternative rotamers for blue residues.
""",justify=LEFT).pack()
    Button(p2, text="Reset View", command=colorize).pack(side=BOTTOM)
############################ PROXY TAB #####################################
    group1 = Pmw.Group(p3,tag_text='Proxy setting')
    group1.pack(fill='both', expand=1)
# auto-detect vs manual
    conection = StringVar(master=group1.interior())
    conection.set('auto')
    Radiobutton(group1.interior(), text="Auto-detect proxy configuration", variable=conection, value='auto', command=disable_entry).grid(row=0, columnspan=3)
    Radiobutton(group1.interior(), text="Manual proxy configuration       ", variable=conection, value='manual', command=enable_entry).grid(row=1, columnspan=3, ipady=3)
# http proxy and Port
    Label(group1.interior(), text='http proxy').grid(row=2, column=0)
    httpProxy = StringVar(master=group1.interior())
    entry_proxy = Entry(group1.interior(),textvariable=httpProxy, width=20)
    entry_proxy.grid(row=2, column=1)
    entry_proxy.configure(state='disabled')
    entry_proxy.update()
    Label(group1.interior(), text='Port').grid(row=3, column=0)
    port = StringVar(master=group1.interior())
    entry_port = Entry(group1.interior(),textvariable=port, width=20)
    entry_port.grid(row=3, column=1)
    entry_port.configure(state='disabled')
    entry_port.update()
# enable authentication
    authentication = IntVar(master=group1.interior())
    check_authentication=Checkbutton(group1.interior(),text='authentication', variable=authentication, command=enable_entry_authentication, width=20)
    check_authentication.grid(row=5, rowspan=2, column=2)
    check_authentication.configure(state='disabled')
    check_authentication.update()
# username and password
    Label(group1.interior(), text='username').grid(row=5, column=0)
    username = StringVar(master=group1.interior())
    entry_username=Entry(group1.interior(),textvariable=username, width=20)
    entry_username.grid(row=5, column=1)
    entry_username.configure(state='disabled')
    entry_username.update()
    Label(group1.interior(), text='password').grid(row=6, column=0)
    password = StringVar(master=group1.interior())
    entry_password = Entry(group1.interior(),textvariable=password, show='*', width=20)
    entry_password.grid(row=6, column=1)
    entry_password.configure(state='disabled')
    entry_password.update()
############################ About TAB ########################################
    Label(p4, text = """
If you find CheShift useful please cite:

Martin O.A. Arnautova Y.A. Icazatti A.A. 
Scheraga H.A. and Vila J.A. 
A Physics-Based Method to Validate and 
Repair Flaws in Protein Structures. 
Proc Natl Acad Sci USA 2013. (in press). 
""",justify=CENTER).pack()
    master.mainloop()

def retrieve_pdb():
    """Loads a PDB file provided by the user"""
    global pdb_path
    pdb_path = tkFileDialog.askopenfilename(title="Open PDB file", filetypes=[("PDB files", ".pdb"),("All files",".*")])
    if len(pdb_path) == 0:
        del pdb_path
    else:
        cmd.reinitialize()
        cmd.load(pdb_path)


def retrieve_cs():
    """Loads a Chemical Shift file provided by the user"""
    global cs_path
    cs_path = tkFileDialog.askopenfilename(title="Open chemical shift file", filetypes=[("All files","*")])
    if len(cs_path) == 0:
        del cs_path

def enable_entry():
    """enables the fields for proxy and port"""
    entry_port.configure(state='normal')
    entry_port.update()
    entry_proxy.configure(state='normal')
    entry_proxy.update()
    check_authentication.configure(state='normal')
    check_authentication.update()

def disable_entry():
    """disables all the fields related to the proxy tab"""
    entry_port.configure(state='disabled')
    entry_port.update()
    entry_proxy.configure(state='disabled')
    entry_proxy.update()
    check_authentication.configure(state='disabled')
    check_authentication.update()
    entry_username.configure(state='disabled')
    entry_username.update()
    entry_password.configure(state='disabled')
    entry_password.update()

def enable_entry_authentication():
    """enables or disable the fields username and password"""
    if authentication.get() == 1:
        entry_username.configure(state='normal')
        entry_username.update()
        entry_password.configure(state='normal')
        entry_password.update()
    else:
        entry_username.configure(state='disabled')
        entry_username.update()
        entry_password.configure(state='disabled')
        entry_password.update()

def colorize():
    """color according to the b-factor using the cheshift-code"""
    try:
        cmd.spectrum('b', 'red_yellow_green', minimum='-1.0', maximum='0.0')
        cmd.select('missing', 'b = -2.0')
        cmd.color('white','missing')
        cmd.delete('missing')
        cmd.select('fixable', 'b = 2.0')
        cmd.color('blue','fixable')
        cmd.delete('fixable')
        cmd.hide()
        cmd.show('cartoon')
    except:
        pass

def guess_proxy_settings():
    """try to guess the correct proxy settings on Linux (KDE and Gnome)"""
# Mechanize can detect the proxy setting on Windows, MacOsX and sometimes on Linux. On Linux there is not standard way to set the proxy settings, Mechanize can only checks for the existence of http_proxy environment variable (but proxy setting can be stored on other "locations"). This function uses two methods one restricted to gnome the other to KDE. I don`t test the KDE method, because I don`t have a KDE machine.
    proxy_values={}
    osystem = sys.platform
    if osystem.startswith('linux'):
        try:
            mode = gconf.client_get_default().get_string('/system/proxy/mode')
            if mode != None:
                g_httpProxy = gconf.client_get_default().get_string('/system/http_proxy/host')
                g_port = gconf.client_get_default().get_int('/system/http_proxy/port')
                g_username = gconf.client_get_default().get_string('/system/http_proxy/authentication_user')
                g_password = gconf.client_get_default().get_string('/system/http_proxy/authentication_password')
                g_authentication2 = gconf.client_get_default().get_bool('/system/http_proxy/use_authentication')
                if g_authentication2 is True:
                    proxy_values['http'] = g_username+':'+g_password+'@'+g_httpProxy+':'+str(g_port)
                else:
                    proxy_values['http'] = g_httpProxy+':'+str(g_port)
        except:
            try:
                result = os.popen(r'kreadconfig --file kioslaverc --group "Proxy Settings" --key httpProxy').read()
                proxy_values['http']=result.split('//')[2]
            except:
                pass
    return proxy_values


def user_proxy_settings():
    """save the proxy settings provided by the user"""
    if conection.get() == 'manual':
        if authentication.get() == 1:
            if len(username.get()) != 0 and len(password.get()) != 0 and len(httpProxy.get()) != 0 and len(port.get()) != 0:
                proxy_values['http'] = username.get()+':'+password.get()+'@'+httpProxy.get()+':'+port.get()
        else:
            if len(httpProxy.get()) != 0 and len(port.get()) != 0:
                proxy_values['http'] = httpProxy.get()+':'+port.get()
#################### Saves the proxy setting for future use ###################
    if len(proxy_values) != 0:
        osystem = sys.platform
        if osystem.startswith('linux'):
            fd = open('%s/.proxy_cache' % (os.getenv('HOME')),'w')
            fd.write('%s' % proxy_values['http'])
            fd.close()
        elif osystem == 'win32':
            fd = open('%s\proxy_cache' % (os.getcwd()),'w')
            fd.write('%s' % proxy_values['http'])
            fd.close()
        elif osystem == 'darwin':
            fd = open('%s/Lybrary/.proxy_cache' % (os.path.expanduser("~")),'w')
            fd.write('%s' % proxy_values['http'])
            fd.close()
        else:
            pass
    return proxy_values

def user_proxy_file():
    proxy_values={}
    try:
        osystem = sys.platform
        if osystem.startswith('linux'):
            for line in open('%s/.proxy_cache' % os.getenv('HOME')).readlines():
                proxy_values['http'] = line
        elif osystem == 'win32':
            for line in open('%s\proxy_cache' % os.getcwd()).readlines():
                proxy_values['http'] = line
        elif osystem == 'darwin':
            for line in open('%s/Lybrary/.proxy_cache' % os.path.expanduser("~")).readlines():
                proxy_values['http'] = line
        return proxy_values
    except:
        return proxy_values

def test_internet_conection():
    br = mechanize.Browser()
    br.set_proxies(proxy_values)
    try:
        response = br.open('http://www.cheshift.com/visual.html', timeout=2)
        return True
    except:
        pass
        return False

def submit():
    """Checks if files were provided, call functions related to the proxy
configuration and if everything seems Ok calls the cheshift function"""
    global proxy_values
    proxy_values = {}
    pdb = 0
    cs = 0
    try:
        pdb_path
        pdb = 1
    except:
        Pmw.MessageDialog(title = 'Error',message_text = 'Please choose a\n PDB file')
    if pdb == 1:
        try:
            cs_path
            cs = 1
        except:
                Pmw.MessageDialog(title = 'Error',message_text = 'Please choose a\n Chemical Shift file')
    if pdb == 1 and cs == 1:
####################### checks the Internet connection ########################
        if test_internet_conection() == True:
            cheshift(pdb_path, cs_path)
        else:
            proxy_values = user_proxy_file()
            if test_internet_conection() == True:
                    cheshift(pdb_path, cs_path)
            elif conection.get() == 'manual': # uses the values provided by the user
                    proxy_values = user_proxy_settings()
                    if test_internet_conection() == True:
                        cheshift(pdb_path, cs_path)
                    else:
                        Pmw.MessageDialog(title='Error', message_text=('I can`t connect to the Internet,\n  please check your Internet connection\n and your proxy configuration' ))
            else:
                proxy_values = guess_proxy_settings()
                if test_internet_conection() == True:
                   cheshift(pdb_path, cs_path)
                else:
                    Pmw.MessageDialog(title='Error', message_text=('I can`t connect to the Internet,\n  please check your Internet connection\n  and your proxy configuration' ))



#cmd.extend('cheshift', cheshift)
