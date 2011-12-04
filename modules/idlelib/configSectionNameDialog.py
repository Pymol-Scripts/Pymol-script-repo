"""
Dialog that allows user to specify a new config file section name.
Used to get new highlight theme and keybinding set names.
"""
from Tkinter import Toplevel, Entry, Frame, Label, Button, StringVar
from Tkconstants import TOP, BOTTOM, RIGHT, BOTH, SUNKEN, X
import tkMessageBox

from configHandler import idleConf

TTK = idleConf.GetOption('main', 'General', 'use-ttk', type='int')
if TTK:
    from ttk import Entry, Frame, Label, Button

class GetCfgSectionNameDialog(Toplevel):
    def __init__(self, parent, title, message, usedNames):
        """
        message - string, informational message to display
        usedNames - list, list of names already in use for validity check
        """
        Toplevel.__init__(self, parent)
        self.configure(borderwidth=5)
        self.resizable(height=False, width=False)
        self.title(title)
        self.transient(parent)
        self.grab_set()
        self.protocol("WM_DELETE_WINDOW", self.Cancel)
        self.parent = parent
        self.message = message
        self.usedNames = usedNames
        self.result = ''
        self.CreateWidgets()
        self.withdraw() #hide while setting geometry
        self.update_idletasks()
        self.geometry("+%d+%d" %
            ((parent.winfo_rootx() + ((parent.winfo_width() / 2)
                - (self.winfo_reqwidth() / 2)),
              parent.winfo_rooty() + ((parent.winfo_height() / 2)
                - (self.winfo_reqheight() / 2)) )) ) #centre dialog over parent
        self.deiconify() #geometry set, unhide
        self.wait_window()

    def CreateWidgets(self):
        self.name = StringVar(self)
        self.fontSize = StringVar(self)
        self.frameMain = Frame(self, borderwidth=2, relief=SUNKEN)
        self.messageInfo = Label(self.frameMain, text=self.message)
        entryName = Entry(self.frameMain, textvariable=self.name, width=30)
        frameButtons = Frame(self)
        self.buttonOk = Button(frameButtons, text='Ok', command=self.Ok)
        self.buttonCancel = Button(frameButtons, text='Cancel',
                command=self.Cancel)

        entryName.focus_set()

        self.frameMain.pack(side=TOP, expand=True, fill=BOTH)
        self.messageInfo.pack(padx=5, pady=5)
        entryName.pack(padx=5, pady=5)
        frameButtons.pack(side=BOTTOM, fill=X)
        self.buttonOk.pack(padx=1, pady=5, side=RIGHT)
        self.buttonCancel.pack(pady=5, padx=5, side=RIGHT)

        if TTK:
            self.messageInfo['padding'] = 5
            frameButtons['style'] = 'RootColor.TFrame'
        else:
            self.messageInfo.configure(padx=5, pady=5)

    def NameOk(self):
        #simple validity check for a sensible
        #ConfigParser file section name
        nameOk = 1
        name = self.name.get()
        name.strip()

        if not name: #no name specified
            tkMessageBox.showerror(title='Name Error',
                message='No name specified.', parent=self)
            nameOk = 0
        elif len(name) > 30: #name too long
            tkMessageBox.showerror(title='Name Error',
                message=('Name too long. It should be no more than '
                         '30 characters.'), parent=self)
            nameOk=0
        elif name in self.usedNames:
            tkMessageBox.showerror(title='Name Error',
                    message='This name is already in use.', parent=self)
            nameOk=0

        return nameOk

    def Ok(self, event=None):
        if self.NameOk():
            self.result = self.name.get().strip()
            self.destroy()

    def Cancel(self, event=None):
        self.result = ''
        self.destroy()

if __name__ == '__main__':
    from Tkinter import Tk
    #test the dialog
    def run():
        keySeq = ''
        dlg = GetCfgSectionNameDialog(root, 'Get Name',
            'The information here should need to be word wrapped. Test.', [])
        print dlg.result

    root=Tk()
    Button(root, text='Dialog', command=run).pack()
    root.mainloop()
