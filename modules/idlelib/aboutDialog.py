"""About Dialog for IDLE"""
import os
import sys
from Tkinter import Toplevel, Frame, Button, Label, TkVersion
from Tkconstants import LEFT, NSEW, SUNKEN, EW, W, BOTH, TOP, BOTTOM

import idlever
import textView
from stylist import PoorManStyle
from configHandler import idleConf

TTK = idleConf.GetOption('main', 'General', 'use-ttk', type='int')
if TTK:
    from ttk import Frame, Button, Label, Style

class AboutDialog(Toplevel):
    """Modal about dialog for idle

    """
    def __init__(self,parent,title):
        Toplevel.__init__(self, parent)
        self.configure(borderwidth=5)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+30,
                                  parent.winfo_rooty()+30))
        self.bg = "#707070"
        self.fg = "#ffffff"
 
        self.SetupStyles()
        self.CreateWidgets()
        self.resizable(height=False, width=False)
        self.title(title)
        self.transient(parent)
        self.grab_set()
        self.protocol("WM_DELETE_WINDOW", self.Ok)
        self.parent = parent
        self.buttonOk.focus_set()
        self.bind('<Return>',self.Ok) #dismiss dialog
        self.bind('<Escape>',self.Ok) #dismiss dialog
        self.wait_window()

    def SetupStyles(self):
        if TTK:
            style = Style(self.master)
            style.configure('Color.TLabel', foreground=self.fg,
                            background=self.bg)
            style.configure('Color.TFrame', background=self.bg)
            self.ttkstyle = style
            self.style = lambda w, style: w.configure(style=style)
        else:
            self.style = PoorManStyle(self,
                styles={'Color.TLabel': {'fg': self.fg, 'bg': self.bg},
                        'Color.TFrame': {'bg': self.bg}}).style_it

    def CreateWidgets(self):
        frameMain = Frame(self, borderwidth=2, relief=SUNKEN)
        frameButtons = Frame(self)
        frameButtons.pack(side=BOTTOM, pady=3)
        frameMain.pack(side=TOP, expand=True, fill=BOTH)
        self.buttonOk = Button(frameButtons, text='Close', command=self.Ok)
        self.buttonOk.pack()
        frameBg = Frame(frameMain)
        frameBg.pack(expand=True, fill=BOTH)
        labelTitle = Label(frameBg, text='IDLE', font=('courier', 24, 'bold'))
        labelTitle.grid(row=0, column=0, sticky=W, padx=10, pady=10)
        byline = "Python's Integrated DeveLopment Environment" + 5*'\n'
        labelDesc = Label(frameBg, text=byline, justify=LEFT)
        labelDesc.grid(row=2, column=0, sticky=W, columnspan=3, padx=10, pady=5)
        labelEmail = Label(frameBg, text='email:  idle-dev@python.org',
                           justify=LEFT)
        labelEmail.grid(row=6, column=0, columnspan=2,
                        sticky=W, padx=10, pady=0)
        labelWWW = Label(frameBg, text='www:  http://www.python.org/idle/',
                         justify=LEFT)
        labelWWW.grid(row=7, column=0, columnspan=2, sticky=W, padx=10, pady=0)
        fbg = Frame(frameBg, borderwidth=1, relief=SUNKEN,  height=2)
        fbg.grid(row=8, column=0, sticky=EW, columnspan=3, padx=5, pady=5)
        labelPythonVer = Label(frameBg, text='Python version:  ' + \
                               sys.version.split()[0])
        labelPythonVer.grid(row=9, column=0, sticky=W, padx=10, pady=0)
        # handle weird tk version num in windoze python >= 1.6 (?!?)
        tkVer = repr(TkVersion).split('.')
        tkVer[len(tkVer)-1] = str('%.3g' % (float('.'+tkVer[len(tkVer)-1])))[2:]
        if tkVer[len(tkVer)-1] == '':
            tkVer[len(tkVer)-1] = '0'
        tkVer = '.'.join(tkVer)
        labelTkVer = Label(frameBg, text='Tk version:  '+ tkVer)
        labelTkVer.grid(row=9, column=1, sticky=W, padx=2, pady=0)
        py_button_f = Frame(frameBg)
        py_button_f.grid(row=10, column=0, columnspan=2, sticky=NSEW)
        buttonLicense = Button(py_button_f, text='License', width=8,
                               command=self.ShowLicense)
        buttonLicense.pack(side=LEFT, padx=10, pady=10)
        buttonCopyright = Button(py_button_f, text='Copyright', width=8,
                                 command=self.ShowCopyright)
        buttonCopyright.pack(side=LEFT, padx=10, pady=10)
        buttonCredits = Button(py_button_f, text='Credits', width=8,
                               command=self.ShowPythonCredits)
        buttonCredits.pack(side=LEFT, padx=10, pady=10)
        fbg2 = Frame(frameBg, borderwidth=1, relief=SUNKEN, height=2)
        fbg2.grid(row=11, column=0, sticky=EW, columnspan=3, padx=5, pady=5)
        idle_v = Label(frameBg, text='IDLE version:   ' + idlever.IDLE_VERSION)
        idle_v.grid(row=12, column=0, sticky=W, padx=10, pady=0)
        idle_button_f = Frame(frameBg)
        idle_button_f.grid(row=13, column=0, columnspan=3, sticky=NSEW)
        idle_about_b = Button(idle_button_f, text='README', width=8,
                                command=self.ShowIDLEAbout)
        idle_about_b.pack(side=LEFT, padx=10, pady=10)
        idle_news_b = Button(idle_button_f, text='NEWS', width=8,
                                command=self.ShowIDLENEWS)
        idle_news_b.pack(side=LEFT, padx=10, pady=10)
        idle_credits_b = Button(idle_button_f, text='Credits', width=8,
                                command=self.ShowIDLECredits)
        idle_credits_b.pack(side=LEFT, padx=10, pady=10)

        s = self.style
        s(frameButtons, 'RootColor.TFrame')
        s(frameBg, 'Color.TFrame')
        s(labelTitle, 'Color.TLabel')
        s(labelDesc, 'Color.TLabel')
        s(labelEmail, 'Color.TLabel')
        s(labelWWW, 'Color.TLabel')
        s(fbg, 'Color.TFrame')
        s(labelPythonVer, 'Color.TLabel')
        s(labelTkVer, 'Color.TLabel')
        s(py_button_f, 'Color.TFrame')
        s(fbg2, 'Color.TFrame')
        s(idle_v, 'Color.TLabel')
        s(idle_button_f, 'Color.TFrame')

    def ShowLicense(self):
        self.display_printer_text('About - License', license)

    def ShowCopyright(self):
        self.display_printer_text('About - Copyright', copyright)

    def ShowPythonCredits(self):
        self.display_printer_text('About - Python Credits', credits)

    def ShowIDLECredits(self):
        self.display_file_text('About - Credits', 'CREDITS.txt', 'iso-8859-1')

    def ShowIDLEAbout(self):
        self.display_file_text('About - Readme', 'README.txt')

    def ShowIDLENEWS(self):
        self.display_file_text('About - NEWS', 'NEWS.txt')

    def display_printer_text(self, title, printer):
        printer._Printer__setup()
        text = '\n'.join(printer._Printer__lines)
        textView.view_text(self, title, text)

    def display_file_text(self, title, filename, encoding=None):
        fn = os.path.join(os.path.abspath(os.path.dirname(__file__)), filename)
        textView.view_file(self, title, fn, encoding)

    def Ok(self, event=None):
        self.destroy()

if __name__ == '__main__':
    # test the dialog
    from Tkinter import Tk
    root = Tk()
    def run():
        AboutDialog(root, 'About')
    Button(root, text='Dialog', command=run).pack()
    root.mainloop()
