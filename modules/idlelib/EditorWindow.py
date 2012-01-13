import os
import re
import sys
import traceback
import webbrowser
import tkMessageBox
import tkSimpleDialog
from Tkinter import Menu, Scrollbar, TclError, BooleanVar
from Tkconstants import INSERT, END, RIGHT, BOTTOM, TOP, X, Y, BOTH

import macosxSupport
import Bindings
import WindowList
from editorpage import EditorPage, classifyws, filename_to_unicode
from tabbedpages import get_tabbedpage
from configHandler import idleConf
from MultiStatusBar import MultiStatusBar

TabbedPageSet = get_tabbedpage()

TTK = idleConf.GetOption('main', 'General', 'use-ttk', type='int')
if TTK:
    from ttk import Style, Scrollbar

# The default tab setting for a Text widget, in average-width characters.
TK_TABWIDTH_DEFAULT = 8

class EditorWindow(object):
    from ColorDelegator import ColorDelegator # overridden by PyShell
    from UndoDelegator import UndoDelegator   # overridden by PyShell

    help_url = None
    menu_specs = [
        ("file", "_File"),
        ("edit", "_Edit"),
        ("format", "F_ormat"),
        ("run", "_Run"),
        ("options", "_Options"),
        ("windows", "_Windows"),
        ("help", "_Help"),
    ]

    if macosxSupport.runningAsOSXApp():
        del menu_specs[-3]
        menu_specs[-2] = ("windows", "_Window")

    def __init__(self, flist=None, filename=None, key=None, root=None,
       start_page=EditorPage):
        if EditorWindow.help_url is None:
            dochome =  os.path.join(sys.prefix, 'Doc', 'index.html')
            if sys.platform.count('linux'):
                # look for html docs in a couple of standard places
                pyver = 'python-docs-' + '%s.%s.%s' % sys.version_info[:3]
                if os.path.isdir('/var/www/html/python/'):  # "python2" rpm
                    dochome = '/var/www/html/python/index.html'
                else:
                    basepath = '/usr/share/doc/'  # standard location
                    dochome = os.path.join(basepath, pyver,
                                           'Doc', 'index.html')
            elif sys.platform[:3] == 'win':
                chmfile = os.path.join(sys.prefix, 'Doc',
                                       'Python%d%d.chm' % sys.version_info[:2])
                if os.path.isfile(chmfile):
                    dochome = chmfile

            elif macosxSupport.runningAsOSXApp():
                # documentation is stored inside the python framework
                dochome = os.path.join(sys.prefix,
                        'Resources/English.lproj/Documentation/index.html')

            dochome = os.path.normpath(dochome)
            if os.path.isfile(dochome):
                EditorWindow.help_url = dochome
                if sys.platform == 'darwin':
                    # Safari requires real file:-URLs
                    EditorWindow.help_url = 'file://' + EditorWindow.help_url
            else:
                EditorWindow.help_url = "http://www.python.org/doc/current"

        self.flist = flist
        root = root or flist.root
        self.root = root
        try:
            sys.ps1
        except AttributeError:
            sys.ps1 = '>>> '
        self.menubar = Menu(root)
        self.top = WindowList.ListedToplevel(root, menu=self.menubar)
        if flist:
            self.tkinter_vars = flist.vars
            # self.top.instance_dict makes flist.inversedict avalable to
            # configDialog.py so it can access all EditorWindow instaces
            self.top.instance_dict = flist.inversedict
        else:
            self.tkinter_vars = {}  # keys: Tkinter event names
                                    # values: Tkinter variable instances
            self.top.instance_dict = {}
        self.recent_files_path = os.path.join(idleConf.GetUserCfgDir(),
            'recent-files.lst')

        if flist:
            flist.inversedict[self] = key
            if key:
                flist.dict[key] = self

        self.menudict = None

        # create a Notebook where the text pages for this EditorWindow will
        # reside
        self.text_notebook = TabbedPageSet(self.top)
        self.text_notebook.pack(fill=BOTH, expand=True)
        self.text_notebook.bind('<<NotebookTabChanged>>', self._update_controls)
        self.new_tab(filename=filename, load_ext=False, ptype=start_page)
        self.text = self.current_page.text # XXX
        self.top.focused_widget = self.text
        self.top.bind('<<tab-closed>>', self._post_tab_close)

        # The following "width" attribute is used by PyShell, so keep it here
        self.width = idleConf.GetOption('main', 'EditorPage', 'width')

        self.top.protocol("WM_DELETE_WINDOW", self.close)
        self.top.bind("<<close-window>>", self.close_event)

        self._create_statusbar()
        self.top.after_idle(self.set_line_and_column)

        # usetabs true  -> literal tab characters are used by indent and
        #                  dedent cmds, possibly mixed with spaces if
        #                  indentwidth is not a multiple of tabwidth,
        #                  which will cause Tabnanny to nag!
        #         false -> tab characters are converted to spaces by indent
        #                  and dedent cmds, and ditto TAB keystrokes
        # Although use-spaces=0 can be configured manually in config-main.def,
        # configuration of tabs v. spaces is not supported in the configuration
        # dialog.  IDLE promotes the preferred Python indentation: use spaces!
        usespaces = idleConf.GetOption('main', 'Indent', 'use-spaces',
            type='bool')
        self.usetabs = not usespaces

        # tabwidth is the display width of a literal tab character.
        # CAUTION:  telling Tk to use anything other than its default
        # tab setting causes it to use an entirely different tabbing algorithm,
        # treating tab stops as fixed distances from the left margin.
        # Nobody expects this, so for now tabwidth should never be changed.
        self.tabwidth = 8    # must remain 8 until Tk is fixed.

        # indentwidth is the number of screen characters per indent level.
        # The recommended Python indentation is four spaces.
        self.indentwidth = self.tabwidth
        self.set_notabs_indentwidth()

        # If context_use_ps1 is true, parsing searches back for a ps1 line;
        # else searches for a popular (if, def, ...) Python stmt.
        self.context_use_ps1 = False

        # When searching backwards for a reliable place to begin parsing,
        # first start num_context_lines[0] lines back, then
        # num_context_lines[1] lines back if that didn't work, and so on.
        # The last value should be huge (larger than the # of lines in a
        # conceivable file).
        # Making the initial values larger slows things down more often.
        self.num_context_lines = 50, 500, 5000000

        if hasattr(self, 'ispythonsource'): # PyShell
            self.set_indentation_params(self.ispythonsource(filename))
        else:
            self.set_indentation_params(
                self.current_page.ispythonsource(filename))

        self.extensions = {}
        self._load_extensions()

        menu = self.menudict.get('windows')
        if menu:
            end = menu.index("end")
            if end is None:
                end = -1
            if end >= 0:
                menu.add_separator()
                end = end + 1
            self.wmenu_end = end
            WindowList.register_callback(self.postwindowsmenu)

        # Some abstractions so IDLE extensions are cross-IDE
        self.askyesno = tkMessageBox.askyesno
        self.askinteger = tkSimpleDialog.askinteger
        self.showerror = tkMessageBox.showerror

    @property
    def current_page(self):
        """Return the active EditorPage in EditorWindow."""
        if not self.text_notebook.pages: # no pages available
            return None
        curr_tab = self.text_notebook.select()
        if not curr_tab:
            return None

        if TTK:
            page = self.text_notebook.pages[self.text_notebook.tab(
                curr_tab)['text']].editpage
        else:
            page = self.text_notebook.pages[curr_tab].editpage
        return page

    def remove_tab_controls(self):
        """Remove tab area and most tab bindings from this window."""
        if TTK:
            self.text_notebook['style'] = 'PyShell.TNotebook'
            style = Style(self.top)
            style.layout('PyShell.TNotebook.Tab', [('null', '')])
        else:
            self.text_notebook._tab_set.grid_forget()

        # remove commands related to tab
        if 'file' in self.menudict:
            menu = self.menudict['file']
            curr_entry = None
            i = 0
            while True:
                last_entry, curr_entry = curr_entry, menu.entryconfigure(i)
                if last_entry == curr_entry:
                    # no more menu entries
                    break

                if 'label' in curr_entry and 'Tab' in curr_entry['label'][-1]:
                    if 'Close' not in ' '.join(curr_entry['label'][-1]):
                        menu.delete(i)
                i += 1

        self.current_page.text.unbind('<<new-tab>>')
        # close-tab is still available!

    def short_title(self):
        # overriden by PyShell
        return self.current_page.short_title()

    def next_tab(self, event):
        """Show next tab if not in the last tab already."""
        index = self.text_notebook.index(self.text_notebook.select())
        if index == len(self.text_notebook.tabs()) - 1:
            return
        self.text_notebook.select(index + 1)

    def prev_tab(self, event):
        """Show the previous tab if not in the first tab already."""
        index = self.text_notebook.index(self.text_notebook.select())
        if index == 0:
            return
        self.text_notebook.select(index - 1)

    def new_tab(self, event=None, filename=None, load_ext=True, ptype=None):
        """Create a new EditorPage and insert it into the notebook."""
        page_title = "#%d" % (len(self.text_notebook.pages) + 1)
        page = self.text_notebook.add_page(page_title)

        vbar = Scrollbar(page.frame, name='vbar')
        hbar = Scrollbar(page.frame, name='hbar', orient='horizontal')
        hbar.set(0, 0)
        vbar.set(0, 0)
        ptype = ptype or EditorPage
        page.editpage = ptype(page.frame, self, title=page_title,
            name='text', padx=5, wrap='none')

        firstpage = False # don't update window's title
        if self.menudict is None:
            # This EditorWindow is being created now, perform the following
            # tasks before.
            firstpage = True # will cause window's title to be updated
            self.menudict = {}
            self._createmenubar(page.editpage.text)
            # Create the recent files submenu
            self.recent_files_menu = Menu(self.menubar)
            self.menudict['file'].insert_cascade(3, label='Recent Files',
                underline=0, menu=self.recent_files_menu)
            self.update_recent_files_list()

        # pack widgets
        text = page.editpage.text
        vbar['command'] = text.yview
        hbar['command'] = text.xview
        text['yscrollcommand'] = vbar.set
        text['xscrollcommand'] = hbar.set
        vbar.pack(side=RIGHT, fill=Y)
        hbar.pack(side=BOTTOM, fill=X)
        fontWeight = 'normal'
        if idleConf.GetOption('main', 'EditorPage', 'font-bold', type='bool'):
            fontWeight = 'bold'
        text.config(font=(idleConf.GetOption('main', 'EditorPage', 'font'),
            idleConf.GetOption('main', 'EditorPage', 'font-size'),
            fontWeight))
        text.pack(side=TOP, fill=BOTH, expand=1)
        text.focus_set()

        self.apply_bindings(tab=page)
        if load_ext:
            self._load_extensions()

        # select the just created page
        self.text_notebook.select(len(self.text_notebook.pages) - 1)

        page.editpage.post_init(filename=filename,
            update_window_title=firstpage)

        self.top.event_generate('<<tab-created>>')
        return "break"

    def new_callback(self, event, page):
        dirname, basename = page.io.defaultfilename()
        self.flist.new(dirname)
        return "break"

    def set_line_and_column(self, event=None):
        # Used by PyShell too
        curr_page = self.current_page
        if not curr_page:
            return

        line, column = curr_page.text.index(INSERT).split('.')
        self.status_bar.set_label('column', 'Col: %s' % column)
        self.status_bar.set_label('line', 'Ln: %s' % line)

    def postwindowsmenu(self):
        # Only called when Windows menu exists
        menu = self.menudict['windows']
        end = menu.index("end")
        if end is None:
            end = -1
        if end > self.wmenu_end:
            menu.delete(self.wmenu_end+1, end)
        WindowList.add_windows_to_menu(menu)

    def newline_and_indent_event(self, event):
        """Call newline_and_indent_event on current EditorPage."""
        self.current_page.newline_and_indent_event(event)

    def get_selection_indices(self):
        """Call get_selection_indices on current EditorPage."""
        return self.current_page.get_selection_indices()

    def build_char_in_string_func(self, startindex):
        """Call build_char_in_string_func on current EditorPage."""
        return self.current_page.build_char_in_string_func(startindex)

    def gotoline(self, lineno):
        page = self.current_page
        text = page.text

        if lineno is not None and lineno > 0:
            text.mark_set("insert", "%d.0" % lineno)
            text.tag_remove("sel", "1.0", "end")
            text.tag_add("sel", "insert", "insert +1l")
            page.center()

    def close_hook(self):
        if self.flist:
            self.flist.unregister_maybe_terminate(self)
            self.flist = None

    def set_close_hook(self, close_hook):
        self.close_hook = close_hook

    def set_theme(self, ttkstyle):
        # called from configDialog.py
        ttkstyle.theme_use(idleConf.GetOption('main', 'Theme', 'displaytheme'))

    def ResetColorizer(self):
        "Update the colour theme"
        # Called from self.filename_change_hook and from configDialog.py
        for page in self.text_notebook.pages.itervalues():
            page.editpage.reset_colorizer()

    def ResetFont(self):
        "Update the text widgets' font if it is changed"
        # Called from configDialog.py
        fontWeight = 'normal'
        if idleConf.GetOption('main', 'EditorPage', 'font-bold', type='bool'):
            fontWeight = 'bold'

        for page in self.text_notebook.pages.itervalues():
            text = page.editpage.text
            text.config(font=(idleConf.GetOption('main', 'EditorPage', 'font'),
                idleConf.GetOption('main', 'EditorPage', 'font-size'),
                fontWeight))

    def RemoveKeybindings(self):
        "Remove the keybindings before they are changed."
        # Called from configDialog.py
        Bindings.default_keydefs = keydefs = idleConf.GetCurrentKeySet()

        for page in self.text_notebook.pages.itervalues():
            text = page.editpage.text
            for event, keylist in keydefs.items():
                text.event_delete(event, *keylist)

        for extensionName in self._get_standard_extension_names():
            xkeydefs = idleConf.GetExtensionBindings(extensionName)
            if xkeydefs:
                for page in self.text_notebook.pages.itervalues():
                    text = page.editpage.text
                    for event, keylist in xkeydefs.items():
                        text.event_delete(event, *keylist)

    def ApplyKeybindings(self):
        "Update the keybindings after they are changed"
        # Called from configDialog.py
        Bindings.default_keydefs = keydefs = idleConf.GetCurrentKeySet()
        self.apply_bindings()
        for extensionName in self._get_standard_extension_names():
            xkeydefs = idleConf.GetExtensionBindings(extensionName)
            if xkeydefs:
                self.apply_bindings(xkeydefs)
        #update menu accelerators
        menuEventDict = {}
        for menu in Bindings.menudefs:
            menuEventDict[menu[0]] = {}
            for item in menu[1]:
                if item:
                    menuEventDict[menu[0]][prepstr(item[0])[1]] = item[1]
        for menubarItem in self.menudict.keys():
            menu = self.menudict[menubarItem]
            end = menu.index(END) + 1
            for index in range(0, end):
                if menu.type(index) == 'command':
                    accel = menu.entrycget(index, 'accelerator')
                    if accel:
                        itemName = menu.entrycget(index, 'label')
                        event = ''
                        if menuEventDict.has_key(menubarItem):
                            if menuEventDict[menubarItem].has_key(itemName):
                                event = menuEventDict[menubarItem][itemName]
                        if event:
                            accel = get_accelerator(keydefs, event)
                            menu.entryconfig(index, accelerator=accel)

    def reset_help_menu_entries(self):
        "Update the additional help entries on the Help menu"
        help_list = idleConf.GetAllExtraHelpSourcesList()
        helpmenu = self.menudict['help']
        # first delete the extra help entries, if any
        helpmenu_length = helpmenu.index(END)
        if helpmenu_length > self.base_helpmenu_length:
            helpmenu.delete((self.base_helpmenu_length + 1), helpmenu_length)
        # then rebuild them
        if help_list:
            helpmenu.add_separator()
            for entry in help_list:
                cmd = self.__extra_help_callback(entry[1])
                helpmenu.add_command(label=entry[0], command=cmd)
        # and update the menu dictionary
        self.menudict['help'] = helpmenu

    def set_notabs_indentwidth(self):
        "Update the indentwidth if changed and not using tabs in this window"
        # Called from configDialog.py
        if not self.usetabs:
            self.indentwidth = idleConf.GetOption('main', 'Indent','num-spaces',
                                                  type='int')

    def update_recent_files_list(self, new_file=None):
        "Load and update the recent files list and menus"
        # IOBinding calls this
        rf_list = []
        if os.path.exists(self.recent_files_path):
            rf_list_file = open(self.recent_files_path,'r')
            try:
                rf_list = rf_list_file.readlines()
            finally:
                rf_list_file.close()
        if new_file:
            new_file = os.path.abspath(new_file) + '\n'
            if new_file in rf_list:
                rf_list.remove(new_file)  # move to top
            rf_list.insert(0, new_file)
        # clean and save the recent files list
        bad_paths = []
        for path in rf_list:
            if '\0' in path or not os.path.exists(path[0:-1]):
                bad_paths.append(path)
        rf_list = [path for path in rf_list if path not in bad_paths]
        ulchars = "1234567890ABCDEFGHIJK"
        rf_list = rf_list[0:len(ulchars)]
        rf_file = open(self.recent_files_path, 'w')
        try:
            rf_file.writelines(rf_list)
        finally:
            rf_file.close()
        # for each edit window instance, construct the recent files menu
        for instance in self.top.instance_dict.keys():
            menu = instance.recent_files_menu
            menu.delete(1, END)  # clear, and rebuild:
            for i, file in enumerate(rf_list):
                file_name = file[0:-1]  # zap \n
                # make unicode string to display non-ASCII chars correctly
                ufile_name = filename_to_unicode(file_name)
                callback = instance.__recent_file_callback(file_name)
                menu.add_command(label=ulchars[i] + " " + ufile_name,
                                 command=callback,
                                 underline=0)

    def get_geometry(self):
        "Return (width, height, x, y)"
        geom = self.top.wm_geometry()
        m = re.match(r"(\d+)x(\d+)\+(-?\d+)\+(-?\d+)", geom)
        tuple = (map(int, m.groups()))
        return tuple

    def close_event(self, event):
        self.close()

    def close(self):
        to_check = self.text_notebook.pages.copy()

        while to_check:
            curr_tab = self.text_notebook.select()
            if TTK:
                page_name = self.text_notebook.tab(curr_tab)['text']
            else:
                page_name = curr_tab
            page = to_check.pop(page_name)
            editpage = page.editpage
            reply = editpage.close_tab()
            if reply == "cancel":
                break

    def _close(self):
        WindowList.unregister_callback(self.postwindowsmenu)
        self._unload_extensions()
        self.tkinter_vars = None

        for page in self.text_notebook.pages.itervalues():
            page.editpage.close()

        self.top.destroy()
        if self.close_hook:
            # unless override: unregister from flist, terminate if last window
            self.close_hook()

    def apply_bindings(self, keydefs=None, tab=None):
        if keydefs is None:
            keydefs = Bindings.default_keydefs

        if tab:
            iter_over = [tab]
        else:
            iter_over = self.text_notebook.pages.itervalues()

        for page in iter_over:
            text = page.editpage.text
            text.keydefs = keydefs
            for event, keylist in keydefs.items():
                if keylist:
                    text.event_add(event, *keylist)

    def getvar(self, name):
        var = self.get_var_obj(name)
        if var:
            value = var.get()
            return value
        else:
            raise NameError, name

    def setvar(self, name, value, vartype=None):
        var = self.get_var_obj(name, vartype)
        if var:
            var.set(value)
        else:
            raise NameError, name

    def get_var_obj(self, name, vartype=None, text=None):
        var = self.tkinter_vars.get(name)
        if not var and vartype:
            # create a Tkinter variable object with self.text as master:
            self.tkinter_vars[name] = var = vartype(text or self.text)
        return var

    # Tk implementations of "virtual text methods" -- each platform
    # reusing IDLE's support code needs to define these for its GUI's
    # flavor of widget.

    # Return the text widget's current view of what a tab stop means
    # (equivalent width in spaces).
    def get_tabwidth(self): # XXX depends on self.text
        current = self.text['tabs'] or TK_TABWIDTH_DEFAULT
        return int(current)

    # Set the text widget's current view of what a tab stop means.
    def set_tabwidth(self, newtabwidth): # XXX depends on self.text
        text = self.text
        if self.get_tabwidth() != newtabwidth:
            pixels = text.tk.call("font", "measure", text["font"],
                                  "-displayof", text.master,
                                  "n" * newtabwidth)
            text.configure(tabs=pixels)

    # If ispythonsource and guess are true, guess a good value for
    # indentwidth based on file content (if possible), and if
    # indentwidth != tabwidth set usetabs false.
    # In any case, adjust the Text widget's view of what a tab
    # character means.
    def set_indentation_params(self, ispythonsource, guess=True):
        if guess and ispythonsource:
            i = self.guess_indent()
            if 2 <= i <= 8:
                self.indentwidth = i
            if self.indentwidth != self.tabwidth:
                self.usetabs = False
        self.set_tabwidth(self.tabwidth)

    # Guess indentwidth from text content.
    # Return guessed indentwidth.  This should not be believed unless
    # it's in a reasonable range (e.g., it will be 0 if no indented
    # blocks are found).
    def guess_indent(self): # XXX depends on self.text
        opener, indented = IndentSearcher(self.text, self.tabwidth).run()
        if opener and indented:
            raw, indentsmall = classifyws(opener, self.tabwidth)
            raw, indentlarge = classifyws(indented, self.tabwidth)
        else:
            indentsmall = indentlarge = 0
        return indentlarge - indentsmall

    # Private methods/attributes

    # extensions won't have more than one instance per window
    _unique_extensions = ['CodeContext', 'ScriptBinding', 'FormatParagraph']

    def _unload_extensions(self):
        for ins in self.extensions.values():
            if hasattr(ins, "close"):
                ins.close()
        self.extensions = {}

    def _load_extension(self, name, tab):
        ext_loaded = self.extensions.get(name)

        try:
            mod = __import__(name, globals(), locals(), [])
        except ImportError:
            print "\nFailed to import extension: ", name
            return

        keydefs = idleConf.GetExtensionBindings(name)

        if name not in self._unique_extensions or not ext_loaded:
            # create a new instance
            cls = getattr(mod, name)
            ins = cls(tab.editpage)
            self.extensions.setdefault(name, []).append(ins)
            if not ext_loaded:
                # create new items in menu only if this is the first time this
                # extension is being loaded in this window
                if hasattr(cls, "menudefs"):
                    self._fill_menus(cls.menudefs, keydefs)
        elif name in self._unique_extensions and ext_loaded:
            # get an existing instance
            ins = self.extensions[name][0]

        if keydefs:
            self.apply_bindings(keydefs, tab)
            for vevent in keydefs.keys():
                methodname = vevent.replace("-", "_")
                while methodname[:1] == '<':
                    methodname = methodname[1:]
                while methodname[-1:] == '>':
                    methodname = methodname[:-1]
                methodname = methodname + "_event"
                if hasattr(ins, methodname):
                    tab.editpage.text.bind(vevent, getattr(ins, methodname))

    def _load_extensions(self):
        self._load_standard_extensions(self.text_notebook.last_page())

    def _load_standard_extensions(self, tab):
        for name in self._get_standard_extension_names():
            try:
                self._load_extension(name, tab)
            except:
                print "Failed to load extension", repr(name)
                traceback.print_exc()

    def _get_standard_extension_names(self):
        return idleConf.GetExtensions(editor_only=True)

    def _post_tab_close(self, event):
        if not self.current_page:
            # no tabs now, close window
            self._close()
            return

    def _update_controls(self, event):
        curr_page = self.current_page
        if not curr_page:
            return

        self.text = curr_page.text
        curr_page.saved_change_hook(True, False) # update window title
        curr_page.text.focus_set()
        self.set_line_and_column()

        # update references in extensions that are loaded only once
        for ext in self._unique_extensions:
            if ext not in self.extensions:
                continue
            ext = self.extensions[ext][0]
            ext.editpage = curr_page

    def _create_statusbar(self):
        self.status_bar = MultiStatusBar(self.top)
        if macosxSupport.runningAsOSXApp():
            # Insert some padding to avoid obscuring some of the statusbar
            # by the resize widget.
            self.status_bar.set_label('_padding1', '    ', side=RIGHT)
        self.status_bar.set_label('column', 'Col: ?', side=RIGHT)
        self.status_bar.set_label('line', 'Ln: ?', side=RIGHT)
        self.status_bar.pack(side=BOTTOM, fill=X)

    def _createmenubar(self, text):
        mbar = self.menubar
        menudict = self.menudict
        for name, label in self.menu_specs:
            underline, label = prepstr(label)
            menudict[name] = menu = Menu(mbar, name=name, tearoff=0)
            mbar.add_cascade(label=label, menu=menu, underline=underline)

        if sys.platform == 'darwin' and '.framework' in sys.executable:
            # Insert the application menu
            menudict['application'] = menu = Menu(mbar, name='apple')
            mbar.add_cascade(label='IDLE', menu=menu)

        self._fill_menus(text=text)
        self.base_helpmenu_length = self.menudict['help'].index(END)
        self.reset_help_menu_entries()

    def _fill_menus(self, menudefs=None, keydefs=None, text=None):
        """Add appropriate entries to the menus and submenus

        Menus that are absent or None in self.menudict are ignored.
        """
        if menudefs is None:
            menudefs = Bindings.menudefs
        if keydefs is None:
            keydefs = Bindings.default_keydefs
        menudict = self.menudict
        for mname, entrylist in menudefs:
            menu = menudict.get(mname)
            if not menu:
                continue
            for entry in entrylist:
                if not entry:
                    menu.add_separator()
                else:
                    label, eventname = entry
                    checkbutton = (label[:1] == '!')
                    if checkbutton:
                        label = label[1:]
                    underline, label = prepstr(label)
                    accelerator = get_accelerator(keydefs, eventname)
                    def command(eventname=eventname):
                        self.text.event_generate(eventname)
                    if checkbutton:
                        var = self.get_var_obj(eventname, BooleanVar, text)
                        menu.add_checkbutton(label=label, underline=underline,
                            command=command, accelerator=accelerator,
                            variable=var)
                    else:
                        menu.add_command(label=label, underline=underline,
                            command=command, accelerator=accelerator)

    def __recent_file_callback(self, file_name):
        def open_recent_file(fn_closure=file_name):
            self.current_page.io.open(editFile=fn_closure)
        return open_recent_file

    def __extra_help_callback(self, helpfile):
        "Create a callback with the helpfile value frozen at definition time"
        def display_extra_help(helpfile=helpfile):
            if not helpfile.startswith(('www', 'http')):
                url = os.path.normpath(helpfile)
            if sys.platform[:3] == 'win':
                os.startfile(helpfile)
            else:
                webbrowser.open(helpfile)
        return display_extra_help

# Look at the leading whitespace in s.
# Return pair (# of leading ws characters,
#              effective # of leading blanks after expanding
#              tabs to width tabwidth)

import tokenize
_tokenize = tokenize
del tokenize

class IndentSearcher(object):

    # .run() chews over the Text widget, looking for a block opener
    # and the stmt following it.  Returns a pair,
    #     (line containing block opener, line containing stmt)
    # Either or both may be None.

    def __init__(self, text, tabwidth):
        self.text = text
        self.tabwidth = tabwidth
        self.i = self.finished = 0
        self.blkopenline = self.indentedline = None

    def readline(self):
        if self.finished:
            return ""
        i = self.i = self.i + 1
        mark = repr(i) + ".0"
        if self.text.compare(mark, ">=", "end"):
            return ""
        return self.text.get(mark, mark + " lineend+1c")

    def tokeneater(self, type, token, start, end, line,
                   INDENT=_tokenize.INDENT,
                   NAME=_tokenize.NAME,
                   OPENERS=('class', 'def', 'for', 'if', 'try', 'while')):
        if self.finished:
            pass
        elif type == NAME and token in OPENERS:
            self.blkopenline = line
        elif type == INDENT and self.blkopenline:
            self.indentedline = line
            self.finished = 1

    def run(self):
        save_tabsize = _tokenize.tabsize
        _tokenize.tabsize = self.tabwidth
        try:
            try:
                _tokenize.tokenize(self.readline, self.tokeneater)
            except _tokenize.TokenError:
                # since we cut off the tokenizer early, we can trigger
                # spurious errors
                pass
        finally:
            _tokenize.tabsize = save_tabsize
        return self.blkopenline, self.indentedline

### end autoindent code ###

def prepstr(s):
    # Helper to extract the underscore from a string, e.g.
    # prepstr("Co_py") returns (2, "Copy").
    i = s.find('_')
    if i >= 0:
        s = s[:i] + s[i+1:]
    return i, s


keynames = {
 'bracketleft': '[',
 'bracketright': ']',
 'slash': '/',
}

def get_accelerator(keydefs, eventname):
    keylist = keydefs.get(eventname)
    if not keylist:
        return ""
    s = keylist[0]
    s = re.sub(r"-[a-z]\b", lambda m: m.group().upper(), s)
    s = re.sub(r"\b\w+\b", lambda m: keynames.get(m.group(), m.group()), s)
    s = re.sub("Key-", "", s)
    s = re.sub("Cancel","Ctrl-Break",s)   # dscherer@cmu.edu
    s = re.sub("Control-", "Ctrl-", s)
    s = re.sub("-", "+", s)
    s = re.sub("><", " ", s)
    s = re.sub("<", "", s)
    s = re.sub(">", "", s)
    return s


def fixwordbreaks(root):
    # Make sure that Tk's double-click and next/previous word
    # operations use our definition of a word (i.e. an identifier)
    tk = root.tk
    tk.call('tcl_wordBreakAfter', 'a b', 0) # make sure word.tcl is loaded
    tk.call('set', 'tcl_wordchars', '[a-zA-Z0-9_]')
    tk.call('set', 'tcl_nonwordchars', '[^a-zA-Z0-9_]')


def test():
    from Tkinter import Tk
    root = Tk()
    fixwordbreaks(root)
    root.withdraw()
    if sys.argv[1:]:
        filename = sys.argv[1]
    else:
        filename = None
    edit = EditorWindow(root=root, filename=filename)
    edit.set_close_hook(root.quit)
    edit.text.bind("<<close-all-windows>>", edit.close_event)
    root.mainloop()
    root.destroy()

if __name__ == '__main__':
    test()
