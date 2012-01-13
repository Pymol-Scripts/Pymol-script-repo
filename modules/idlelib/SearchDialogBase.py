from Tkinter import Toplevel, Frame, Label, Entry, Button, Checkbutton, \
                    Radiobutton

from configHandler import idleConf

if idleConf.GetOption('main', 'General', 'use-ttk', type='int'):
    from ttk import Frame, Label, Entry, Button, Checkbutton, Radiobutton

class SearchDialogBase:

    title = "Search Dialog"
    icon = "Search"
    needwrapbutton = 1
    bottom_btns = None

    def __init__(self, root, engine):
        self.root = root
        self.engine = engine
        self.ttop = None

    def open(self, text, searchphrase=None):
        self.text = text
        if not self.ttop:
            self.create_widgets()
        else:
            self.ttop.deiconify()
            self.ttop.tkraise()
        if searchphrase:
            self.ent.delete(0, "end")
            self.ent.insert("end", searchphrase)
        self.ent.focus_set()
        self.ent.selection_range(0, "end")
        self.ent.icursor(0)
        self.ttop.grab_set()

    def close(self, event=None):
        if self.ttop:
            self.ttop.grab_release()
            self.ttop.withdraw()

    def create_widgets(self):
        top = Toplevel(self.root)
        top.bind("<Return>", self.default_command)
        top.bind("<Escape>", self.close)
        top.protocol("WM_DELETE_WINDOW", self.close)
        top.wm_title(self.title)
        top.wm_iconname(self.icon)
        top.resizable(height=False, width=False)
        self.ttop = top
        self.top = Frame(top)

        self.row = 0
        self.top.grid(sticky='news')

        self.create_entries()
        self.create_option_buttons()
        self.create_other_buttons()
        self.create_command_buttons()


    def make_entry(self, label, var):
        l = Label(self.top, text=label)
        l.grid(row=self.row, column=0, sticky="ne", padx=6, pady=6)
        e = Entry(self.top, textvariable=var, exportselection=0)
        e.grid(row=self.row, column=1, sticky="nwe", padx=6, pady=6)
        self.row = self.row + 1
        return e

    def make_frame(self,labeltext=None):
        if labeltext:
            l = Label(self.top, text=labeltext)
            l.grid(row=self.row, column=0, sticky="ne", padx=6, pady=6)
        f = Frame(self.top)
        f.grid(row=self.row, column=1, columnspan=1, sticky="nwe",
               padx=6, pady=6 if labeltext else 0)
        self.row = self.row + 1
        return f

    def create_entries(self):
        self.ent = self.make_entry("Find", self.engine.patvar)

    def create_option_buttons(self):
        f = self.make_frame("Options")

        btn = Checkbutton(f, variable=self.engine.revar,
                          text="Regular expression")
        btn.pack(side="left", fill="both")
        if self.engine.isre():
            btn.invoke()

        btn = Checkbutton(f, variable=self.engine.casevar, text="Match case")
        btn.pack(side="left", fill="both")
        if self.engine.iscase():
            btn.invoke()

        btn = Checkbutton(f, variable=self.engine.wordvar, text="Whole word")
        btn.pack(side="left", fill="both")
        if self.engine.isword():
            btn.invoke()

        if self.needwrapbutton:
            btn = Checkbutton(f, variable=self.engine.wrapvar,
                              text="Wrap around")
            btn.pack(side="left", fill="both")
            if self.engine.iswrap():
                btn.invoke()

    def create_other_buttons(self):
        f = self.make_frame("Direction")

        btn = Radiobutton(f, variable=self.engine.backvar, value=1, text="Up")
        btn.pack(side="left")
        if self.engine.isback():
            btn.invoke()

        btn = Radiobutton(f, variable=self.engine.backvar, value=0, text="Down")
        btn.pack(side="left")
        if not self.engine.isback():
            btn.invoke()

    def create_command_buttons(self):
        self.bottom_btns = self.bottom_btns or []
        f = Frame(self.top)
        f.grid(row=self.row, column=0, columnspan=len(self.bottom_btns) + 1,
               pady=6)

        column = 0
        b = Button(f, text="Close", command=self.close)
        b.grid(row=self.row, column=column, padx=6, pady=6)
        column += 1

        btns = {}
        for tbtn in self.bottom_btns:
            opts = {'text': tbtn[0], 'command': getattr(self, tbtn[1])}
            if len(tbtn) == 3:
                opts['default'] = tbtn[2] and 'active' or 'normal'

            btns[opts['text']] = Button(f, **opts).grid(row=self.row, padx=6,
                                                        pady=6, column=column)
            column += 1
