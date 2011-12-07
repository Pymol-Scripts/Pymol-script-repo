import os
import sys
import imp
import webbrowser
import tkMessageBox
import tkSimpleDialog
from Tkinter import Text, Menu, TclError

import utils
import textView
import aboutDialog
import configDialog
import macosxSupport
import PyParse
import IOBinding
import GrepDialog
import PathBrowser
import ClassBrowser
import SearchDialog
import ReplaceDialog
from configHandler import idleConf
from MultiCall import MultiCallCreator
from Percolator import Percolator

def classifyws(s, tabwidth):
    raw = effective = 0
    for ch in s:
        if ch == ' ':
            raw = raw + 1
            effective = effective + 1
        elif ch == '\t':
            raw = raw + 1
            effective = (effective // tabwidth + 1) * tabwidth
        else:
            break
    return raw, effective

def index2line(index):
    """"line.col" -> line, as an int"""
    return int(float(index))

def filename_to_unicode(filename):
    """Convert filename to unicode in order to display it in Tk"""
    if isinstance(filename, unicode) or not filename:
        return filename
    else:
        try:
            return filename.decode(IOBinding.filesystemencoding)
        except UnicodeDecodeError:
            # XXX
            try:
                return filename.decode(IOBinding.encoding)
            except UnicodeDecodeError:
                # byte-to-byte conversion
                return filename.decode('iso8859-1')

def _find_module(fullname, path=None):
    """Version of imp.find_module() that handles hierarchical module names"""
    file = None

    for tgt in fullname.split('.'):
        if file is not None:
            file.close()            # close intermediate files
        (file, filename, descr) = imp.find_module(tgt, path)
        if descr[2] == imp.PY_SOURCE:
            break                   # find but not load the source file
        module = imp.load_module(tgt, file, filename, descr)
        try:
            path = module.__path__
        except AttributeError:
            raise ImportError('No source for module %s' % module.__name__)

    return file, filename, descr

class EditorPage(object):
    rmenu = None
    rmenu_specs = [
        ("Set Breakpoint", "<<set-breakpoint-here>>"),
        ("Clear Breakpoint", "<<clear-breakpoint-here>>")
    ]

    def __init__(self, parent_frame, editwin, title=None, **kwargs):
        self.editwin = editwin
        self.title = title
        self.tab_initialized = False
        kwargs.setdefault('width', idleConf.GetOption('main', 'EditorPage',
            'width'))
        kwargs.setdefault('height', idleConf.GetOption('main', 'EditorPage',
            'height'))

        self.text = MultiCallCreator(Text)(parent_frame, **kwargs)
        self.color = None # initialized in reset_colorizer
        self.per = Percolator(self.text)
        self.undo = self.editwin.UndoDelegator()
        self.per.insertfilter(self.undo)
        self.text.undo_block_start = self.undo.undo_block_start
        self.text.undo_block_stop = self.undo.undo_block_stop
        self.io = IOBinding.IOBinding(self)

        self.undo.set_saved_change_hook(self.saved_change_hook)
        self.io.set_filename_change_hook(self.filename_change_hook)
        self.reset_colorizer()
        self._setup_bindings()

    def post_init(self, filename=None, update_window_title=False):
        if filename:
            if os.path.exists(filename) and not os.path.isdir(filename):
                self.io.loadfile(filename)
            else:
                self.io.set_filename(filename)
        self.saved_change_hook(update_window_title=update_window_title)
        self.tab_initialized = True

    def close_tab(self, event=None):
        """Close current tab, if no more tabs present, close the window."""
        if hasattr(self.editwin, 'interp'):
            # this is a PyShell, don't ask to save
            reply = 'yes'
        else:
            reply = str(self.maybesave())
        if reply != "cancel":
            if self.io.filename:
                self.editwin.update_recent_files_list(new_file=self.io.filename)
            self.close()
            self.editwin.text_notebook.remove_page(self.title)
            self.editwin.top.event_generate('<<tab-closed>>')

        return reply

    def close(self):
        """Perform necessary cleanup for this page before closing it."""
        self.io.close()
        self.io = None

        self.undo = None

        if self.color:
            self.color.close(False)
            self.color = None

        self.per.close()
        self.per = None
        self.text = None

    # XXX (1) mark where these functions are used
    def saved_change_hook(self,update_window_title=False,update_tab_title=True):
        short = self.editwin.short_title()
        long = self.long_title()

        if short and long:
            title = short + " - " + long
            tabtitle = os.path.split(long)[-1]
        elif short:
            title = short
            tabtitle = short
        elif long:
            title = long
            tabtitle = os.path.split(long)[-1]
        else:
            title = tabtitle = "Untitled"
        icon = short or long or title
        if not self.get_saved():
            title = "*%s*" % title
            tabtitle = "*%s*" % tabtitle
            icon = "*%s" % icon

        if update_tab_title:
            self.editwin.text_notebook.update_tabtitle(self, tabtitle)
        if update_window_title or (
           update_window_title is None and self.tab_initialized):
            self.editwin.top.wm_title(title)
            self.editwin.top.wm_iconname(icon)

    def get_saved(self):
        return self.undo.get_saved()

    def set_saved(self, flag):
        self.undo.set_saved(flag)

    def filename_change_hook(self):
        if self.editwin.flist:
            self.editwin.flist.filename_changed_edit(self, self.editwin)
        self.saved_change_hook(self.tab_initialized)
        self.editwin.top.update_windowlist_registry(self.editwin)
        self.reset_colorizer()

    def reset_undo(self):
        self.undo.reset_undo()

    def reset_colorizer(self):
        "Update the colour theme"
        # Called from self.filename_change_hook and from configDialog.py
        self.__rmcolorizer()
        self.__addcolorizer()
        theme = idleConf.GetOption('main','Theme','name')
        normal_colors = idleConf.GetHighlight(theme, 'normal')
        cursor_color = idleConf.GetHighlight(theme, 'cursor', fgBg='fg')
        select_colors = idleConf.GetHighlight(theme, 'hilite')

        self.text.config(
            foreground=normal_colors['foreground'],
            background=normal_colors['background'],
            insertbackground=cursor_color,
            selectforeground=select_colors['foreground'],
            selectbackground=select_colors['background'])

    def short_title(self):
        filename = self.io.filename
        if filename:
            filename = os.path.basename(filename)
        # return unicode string to display non-ASCII chars correctly
        return filename_to_unicode(filename)

    def long_title(self):
        # return unicode string to display non-ASCII chars correctly
        return filename_to_unicode(self.io.filename or "")

    def maybesave(self):
        if self.io:
            if not self.get_saved():
                if self.editwin.top.state()!= 'normal':
                    self.editiwn.top.deiconify()
                self.editwin.top.lower()
                self.editwin.top.lift()
            return self.io.maybesave()
    # XXX (1) end

    def center(self, mark="insert"):
        # Used by EditorWindow.gotoline
        text = self.text
        top, bot = self._getwindowlines()
        lineno = self._getlineno(mark)
        height = bot - top
        newtop = max(1, lineno - height//2)
        text.yview(float(newtop))

    def ispythonsource(self, filename):
        if not filename or os.path.isdir(filename):
            return True
        base, ext = os.path.splitext(os.path.basename(filename))
        if os.path.normcase(ext) in (".py", ".pyw"):
            return True
        try:
            f = open(filename)
            line = f.readline()
            f.close()
        except IOError:
            return False
        return line.startswith('#!') and line.find('python') >= 0

    def newline_and_indent_event(self, event):
        # Used by EditorWindow.newline_and_indent_event which is used by PyShell
        text = self.text
        first, last = self.get_selection_indices()
        text.undo_block_start()
        try:
            if first and last:
                text.delete(first, last)
                text.mark_set("insert", first)
            line = text.get("insert linestart", "insert")
            i, n = 0, len(line)
            while i < n and line[i] in " \t":
                i = i+1
            if i == n:
                # the cursor is in or at leading indentation in a continuation
                # line; just inject an empty line at the start
                text.insert("insert linestart", '\n')
                return "break"
            indent = line[:i]
            # strip whitespace before insert point unless it's in the prompt
            i = 0
            last_line_of_prompt = sys.ps1.split('\n')[-1]
            while line and line[-1] in " \t" and line != last_line_of_prompt:
                line = line[:-1]
                i = i+1
            if i:
                text.delete("insert - %d chars" % i, "insert")
            # strip whitespace after insert point
            while text.get("insert") in " \t":
                text.delete("insert")
            # start new line
            text.insert("insert", '\n')

            # adjust indentation for continuations and block
            # open/close first need to find the last stmt
            lno = index2line(text.index('insert'))
            #print self.editwin.indentwidth, self.editwin.tabwidth
            y = PyParse.Parser(self.editwin.indentwidth, self.editwin.tabwidth)
            if not self.editwin.context_use_ps1:
                for context in self.editwin.num_context_lines:
                    startat = max(lno - context, 1)
                    startatindex = `startat` + ".0"
                    rawtext = text.get(startatindex, "insert")
                    y.set_str(rawtext)
                    bod = y.find_good_parse_start(
                              self.editwin.context_use_ps1,
                              self.build_char_in_string_func(startatindex))
                    if bod is not None or startat == 1:
                        break
                y.set_lo(bod or 0)
            else:
                r = text.tag_prevrange("console", "insert")
                if r:
                    startatindex = r[1]
                else:
                    startatindex = "1.0"
                rawtext = text.get(startatindex, "insert")
                y.set_str(rawtext)
                y.set_lo(0)

            c = y.get_continuation_type()
            if c != PyParse.C_NONE:
                # The current stmt hasn't ended yet.
                if c == PyParse.C_STRING_FIRST_LINE:
                    # after the first line of a string; do not indent at all
                    pass
                elif c == PyParse.C_STRING_NEXT_LINES:
                    # inside a string which started before this line;
                    # just mimic the current indent
                    text.insert("insert", indent)
                elif c == PyParse.C_BRACKET:
                    # line up with the first (if any) element of the
                    # last open bracket structure; else indent one
                    # level beyond the indent of the line with the
                    # last open bracket
                    self.__reindent_to(y.compute_bracket_indent())
                elif c == PyParse.C_BACKSLASH:
                    # if more than one line in this stmt already, just
                    # mimic the current indent; else if initial line
                    # has a start on an assignment stmt, indent to
                    # beyond leftmost =; else to beyond first chunk of
                    # non-whitespace on initial line
                    if y.get_num_lines_in_stmt() > 1:
                        text.insert("insert", indent)
                    else:
                        self.__reindent_to(y.compute_backslash_indent())
                else:
                    assert 0, "bogus continuation type %r" % (c,)
                return "break"

            # This line starts a brand new stmt; indent relative to
            # indentation of initial line of closest preceding
            # interesting stmt.
            indent = y.get_base_indent_string()
            text.insert("insert", indent)
            if y.is_block_opener():
                self._smart_indent_event(event)
            elif indent and y.is_block_closer():
                self._smart_backspace_event(event)
            return "break"
        finally:
            text.see("insert")
            text.undo_block_stop()

    # If a selection is defined in the text widget, return (start,
    # end) as Tkinter text indices, otherwise return (None, None)
    def get_selection_indices(self):
        # Used by EditorWindow.get_selection_indices which is used by
        # FormatParagraph
        try:
            first = self.text.index("sel.first")
            last = self.text.index("sel.last")
            return first, last
        except TclError:
            return None, None

    # Our editpage provides a __is_char_in_string function that works
    # with a Tk text index, but PyParse only knows about offsets into
    # a string. This builds a function for PyParse that accepts an
    # offset.
    def build_char_in_string_func(self, startindex):
        # Used by EditorWindow.build_char_in_string_func which is used by
        # HyperParser
        def inner(offset, _startindex=startindex,
                  _icis=self.__is_char_in_string):
            return _icis(_startindex + "+%dc" % offset)
        return inner

    def _setup_bindings(self):
        text = self.text
        def bind_them(to_bind, prefix='_%s'):
            for tb in to_bind:
                prefix_size = tb.count('<')
                method_name = tb[prefix_size:-prefix_size].replace('-', '_')
                text.bind(tb, getattr(self, prefix % method_name.lower()))
            
        actions = ('<<help>>', '<<python-docs>>',
            '<<about-idle>>', '<<open-config-dialog>>', '<<open-module>>',
            '<<cut>>', '<<copy>>', '<<paste>>', '<<select-all>>',
            '<<remove-selection>>', '<<del-word-left>>', '<<del-word-right>>',
            '<<beginning-of-line>>')
        events = ('<<find>>', '<<center-insert>>', '<<find-again>>',
            '<<find-in-files>>', '<<find-selection>>', '<<replace>>',
            '<<goto-line>>', '<<smart-backspace>>', '<<smart-indent>>',
            '<<indent-region>>', '<<dedent-region>>', '<<comment-region>>',
            '<<tabify-region>>', '<<untabify-region>>', '<<toggle-tabs>>',
            '<<change-indentwidth>>')
        parent_actions = ('<<new-tab>>', '<<next-tab>>', '<<prev-tab>>')

        bind_them(actions)
        bind_them(events, prefix="_%s_event")
        for action in parent_actions:
            prefix_size = action.count('<')
            method_name = action[prefix_size:-prefix_size].replace('-', '_')
            text.bind(action, getattr(self.editwin, method_name))

        text.bind('<<close-tab>>', self.close_tab)
        text.bind('<<newline-and-indent>>', self.newline_and_indent_event)
        text.bind("<<do-nothing>>", lambda event: "break")
        text.bind("<Left>", self._move_at_edge_if_selection(0))
        text.bind("<Right>", self._move_at_edge_if_selection(1))
        text.bind("<3>", self._right_menu)
        text.bind('<<set-line-and-column>>', self.editwin.set_line_and_column)
        text.event_add("<<set-line-and-column>>",
                                    "<KeyRelease>", "<ButtonRelease>")

        if self.editwin.flist:
            text.bind("<<open-new-window>>",
                utils.callback(self.editwin.new_callback, self))
            text.bind("<<close-all-windows>>",
                self.editwin.flist.close_all_callback)
            text.bind("<<open-class-browser>>", self._open_class_browser)
            text.bind("<<open-path-browser>>", self._open_path_browser)

        if macosxSupport.runningAsOSXApp():
            # Command-W on editorwindows doesn't work without this.
            text.bind('<<close-window>>', self.editwin.close_event)

    def _help(self, event=None):
        fn = os.path.join(os.path.abspath(os.path.dirname(__file__)),
            'help.txt')
        textView.view_file(self.text, 'Help', fn)

    def _python_docs(self, event=None):
        if sys.platform[:3] == 'win':
            os.startfile(self.editwin.help_url)
        else:
            webbrowser.open(self.editwin.help_url)
        return "break"

    def _about_idle(self, event=None):
        aboutDialog.AboutDialog(self.text, 'About IDLE')

    def _open_class_browser(self, event=None):
        filename = self.io.filename
        if not filename:
            tkMessageBox.showerror("No filename",
                "This buffer has no associated filename",
                master=self.text)
            self.text.focus_set()
            return None
        head, tail = os.path.split(filename)
        base, ext = os.path.splitext(tail)
        ClassBrowser.ClassBrowser(self.editwin.flist, base, [head])

    def _open_path_browser(self, event=None):
        PathBrowser.PathBrowser(self.editwin.flist)

    def _open_config_dialog(self, event=None):
        # When changing colors and saving it, it requires the attribute
        # instance_dict making necessary to pass self.editwin.top as the
        # parent
        configDialog.ConfigDialog(self.editwin.top, 'Settings')

    def _open_module(self, event=None):
        try:
            name = self.text.get("sel.first", "sel.last")
        except TclError:
            name = ""
        else:
            name = name.strip()

        name = tkSimpleDialog.askstring("Module",
            "Enter the name of a Python module\n"
            "to search on sys.path and open:",
            parent=self.text, initialvalue=name)

        if name:
            name = name.strip()
        if not name:
            return
        # XXX Ought to insert current file's directory in front of path
        try:
            (f, file, (suffix, mode, type)) = _find_module(name)
        except (NameError, ImportError), msg:
            tkMessageBox.showerror("Import error", str(msg), parent=self.text)
            return

        if type != imp.PY_SOURCE:
            tkMessageBox.showerror("Unsupported type",
                "%s is not a source module" % name, parent=self.text)
            return
        if f:
            f.close()
        if self.editwin.flist:
            if idleConf.GetOption('main', 'EditorWindow', 'file-in-tab',
               default=1, type='bool'):
                self.editwin.flist.open(file)
        else:
            self.io.loadfile(file)

    def _find_event(self, event):
        SearchDialog.find(self.text)
        return "break"

    def _find_again_event(self, event):
        SearchDialog.find_again(self.text)
        return "break"

    def _find_selection_event(self, event):
        SearchDialog.find_selection(self.text)
        return "break"

    def _find_in_files_event(self, event):
        GrepDialog.grep(self.text, self.io, self.editwin.flist)
        return "break"

    def _replace_event(self, event):
        ReplaceDialog.replace(self.text)
        return "break"

    def _center_insert_event(self, event):
        self.center()

    def _getwindowlines(self):
        text = self.text
        top = self._getlineno("@0,0")
        bot = self._getlineno("@0,65535")
        if top == bot and text.winfo_height() == 1:
            # Geometry manager hasn't run yet
            height = int(text['height'])
            bot = top + height - 1
        return top, bot

    def _getlineno(self, mark="insert"):
        return int(float(self.text.index(mark)))

    def _goto_line_event(self, event):
        text = self.text
        lineno = tkSimpleDialog.askinteger("Goto", "Go to line number:",
            parent=text)

        if lineno is None:
            return "break"

        if lineno <= 0:
            text.bell()
            return "break"

        text.mark_set("insert", "%d.0" % lineno)
        text.see("insert")

    def _cut(self, event):
        self.text.event_generate("<<Cut>>")
        return "break"

    def _copy(self, event):
        if not self.text.tag_ranges("sel"):
            # There is no selection, so do nothing and maybe interrupt.
            return

        self.text.event_generate("<<Copy>>")
        return "break"

    def _paste(self, event):
        self.text.event_generate("<<Paste>>")
        self.text.see("insert")
        return "break"

    def _select_all(self, event=None):
        self.text.tag_add("sel", "1.0", "end-1c")
        self.text.mark_set("insert", "1.0")
        self.text.see("insert")
        return "break"

    def _remove_selection(self, event=None):
        self.text.tag_remove("sel", "1.0", "end")
        self.text.see("insert")

    def _del_word_left(self, event):
        self.text.event_generate('<Meta-Delete>')
        return "break"

    def _del_word_right(self, event):
        self.text.event_generate('<Meta-d>')
        return "break"

    def _move_at_edge_if_selection(self, edge_index):
        """Cursor move begins at start or end of selection
        When a left/right cursor key is pressed create and return to Tkinter a
        function which causes a cursor move from the associated edge of the
        selection.
        """
        text_index = self.text.index
        text_mark_set = self.text.mark_set
        edges_table = ("sel.first+1c", "sel.last-1c")

        def move_at_edge(event):
            if (event.state & 5) == 0: # no shift(==1) or control(==4) pressed
                try:
                    text_index("sel.first")
                    text_mark_set("insert", edges_table[edge_index])
                except TclError:
                    pass

        return move_at_edge

    def _beginning_of_line(self, event):
        if (event.state & 12) != 0 and event.keysym == "Home":
            # state&1==shift, state&4==control, state&8==alt
            return # <Modifier-Home>; fall back to class binding

        text = self.text

        if text.index("iomark") and \
           text.compare("iomark", "<=", "insert lineend") and \
           text.compare("insert linestart", "<=", "iomark"):
            insertpt = int(text.index("iomark").split(".")[1])
        else:
            line = text.get("insert linestart", "insert lineend")
            for insertpt in xrange(len(line)):
                if line[insertpt] not in (' ','\t'):
                    break
            else:
                insertpt=len(line)

        lineat = int(text.index("insert").split('.')[1])

        if insertpt == lineat:
            insertpt = 0

        dest = "insert linestart+%sc" % str(insertpt)

        if (event.state & 1) == 0:
            # shift not pressed
            text.tag_remove("sel", "1.0", "end")
        else:
            if not text.index("sel.first"):
                text.mark_set("anchor", "insert")

            first = text.index(dest)
            last = text.index("anchor")

            if text.compare(first, ">", last):
                first, last = last, first

            text.tag_remove("sel", "1.0", "end")
            text.tag_add("sel", first, last)

        text.mark_set("insert", dest)
        text.see("insert")
        return "break"

    def _smart_backspace_event(self, event):
        text = self.text
        first, last = self.get_selection_indices()
        if first and last:
            text.delete(first, last)
            text.mark_set("insert", first)
            return "break"
        # Delete whitespace left, until hitting a real char or closest
        # preceding virtual tab stop.
        chars = text.get("insert linestart", "insert")
        if chars == '':
            if text.compare("insert", ">", "1.0"):
                # easy: delete preceding newline
                text.delete("insert-1c")
            else:
                text.bell()     # at start of buffer
            return "break"
        if  chars[-1] not in " \t":
            # easy: delete preceding real char
            text.delete("insert-1c")
            return "break"
        # Ick.  It may require *inserting* spaces if we back up over a
        # tab character!  This is written to be clear, not fast.
        tabwidth = self.editwin.tabwidth
        have = len(chars.expandtabs(tabwidth))
        assert have > 0
        want = ((have - 1) // 
            self.editwin.indentwidth) * self.editwin.indentwidth
        # Debug prompt is multilined....
        last_line_of_prompt = sys.ps1.split('\n')[-1]
        ncharsdeleted = 0
        while 1:
            if chars == last_line_of_prompt:
                break
            chars = chars[:-1]
            ncharsdeleted = ncharsdeleted + 1
            have = len(chars.expandtabs(tabwidth))
            if have <= want or chars[-1] not in " \t":
                break
        text.undo_block_start()
        text.delete("insert-%dc" % ncharsdeleted, "insert")
        if have < want:
            text.insert("insert", ' ' * (want - have))
        text.undo_block_stop()
        return "break"

    def _smart_indent_event(self, event):
        # if intraline selection:
        #     delete it
        # elif multiline selection:
        #     do indent-region
        # else:
        #     indent one level
        text = self.text
        first, last = self.get_selection_indices()
        text.undo_block_start()
        try:
            if first and last:
                if index2line(first) != index2line(last):
                    return self._indent_region_event(event)
                text.delete(first, last)
                text.mark_set("insert", first)
            prefix = text.get("insert linestart", "insert")
            raw, effective = classifyws(prefix, self.editwin.tabwidth)
            if raw == len(prefix):
                # only whitespace to the left
                self.__reindent_to(effective + self.editwin.indentwidth)
            else:
                # tab to the next 'stop' within or to right of line's text:
                if self.editwin.usetabs:
                    pad = '\t'
                else:
                    effective = len(prefix.expandtabs(self.editwin.tabwidth))
                    n = self.editwin.indentwidth
                    pad = ' ' * (n - effective % n)
                text.insert("insert", pad)
            text.see("insert")
            return "break"
        finally:
            text.undo_block_stop()

    def _indent_region_event(self, event):
        head, tail, chars, lines = self.__get_region()
        for pos in range(len(lines)):
            line = lines[pos]
            if line:
                raw, effective = classifyws(line, self.editwin.tabwidth)
                effective = effective + self.editwin.indentwidth
                lines[pos] = self.__make_blanks(effective) + line[raw:]
        self.__set_region(head, tail, chars, lines)
        return "break"

    def _dedent_region_event(self, event):
        head, tail, chars, lines = self.__get_region()
        for pos in range(len(lines)):
            line = lines[pos]
            if line:
                raw, effective = classifyws(line, self.editwin.tabwidth)
                effective = max(effective - self.editwin.indentwidth, 0)
                lines[pos] = self.__make_blanks(effective) + line[raw:]
        self.__set_region(head, tail, chars, lines)
        return "break"

    def _comment_region_event(self, event):
        head, tail, chars, lines = self.__get_region()
        for pos in range(len(lines) - 1):
            line = lines[pos]
            lines[pos] = '##' + line
        self.__set_region(head, tail, chars, lines)

    def _uncomment_region_event(self, event):
        head, tail, chars, lines = self.__get_region()
        for pos in range(len(lines)):
            line = lines[pos]
            if not line:
                continue
            if line[:2] == '##':
                line = line[2:]
            elif line[:1] == '#':
                line = line[1:]
            lines[pos] = line
        self.__set_region(head, tail, chars, lines)

    def _tabify_region_event(self, event):
        head, tail, chars, lines = self.__get_region()
        tabwidth = self.__asktabwidth()
        for pos in range(len(lines)):
            line = lines[pos]
            if line:
                raw, effective = classifyws(line, tabwidth)
                ntabs, nspaces = divmod(effective, tabwidth)
                lines[pos] = '\t' * ntabs + ' ' * nspaces + line[raw:]
        self.__set_region(head, tail, chars, lines)

    def _untabify_region_event(self, event):
        head, tail, chars, lines = self.__get_region()
        tabwidth = self.__asktabwidth()
        for pos in range(len(lines)):
            lines[pos] = lines[pos].expandtabs(tabwidth)
        self.__set_region(head, tail, chars, lines)

    def _toggle_tabs_event(self, event):
        if self.editwin.askyesno(
              "Toggle tabs",
              "Turn tabs " + ("on", "off")[self.editwin.usetabs] +
              "?\nIndent width " +
              ("will be", "remains at")[self.editwin.usetabs] + " 8." +
              "\n Note: a tab is always 8 columns",
              parent=self.text):
            self.editwin.usetabs = not self.editwin.usetabs
            # Try to prevent inconsistent indentation.
            # User must change indent width manually after using tabs.
            self.editwin.indentwidth = 8
        return "break"

    def _change_indentwidth_event(self, event):
        new = self.editwin.askinteger(
                  "Indent width",
                  "New indent width (2-16)\n(Always use 8 when using tabs)",
                  parent=self.text,
                  initialvalue=self.editwin.indentwidth,
                  minvalue=2,
                  maxvalue=16)
        if new and new != self.editwin.indentwidth and not self.editwin.usetabs:
            self.editwin.indentwidth = new
        return "break"

    def _right_menu(self, event):
        self.text.tag_remove("sel", "1.0", "end")
        self.text.mark_set("insert", "@%d,%d" % (event.x, event.y))

        if not self.rmenu:
            self.__make_rmenu()

        rmenu = self.rmenu
        self.event = event
        iswin = sys.platform[:3] == 'win'
        if iswin:
            self.text.config(cursor="arrow")
        rmenu.tk_popup(event.x_root, event.y_root)
        if iswin:
            self.text.config(cursor="ibeam")

    def __make_rmenu(self):
        rmenu = Menu(self.text, tearoff=False)

        for label, eventname in self.rmenu_specs:
            def command(text=self.text, eventname=eventname):
                text.event_generate(eventname)
            rmenu.add_command(label=label, command=command)

        self.rmenu = rmenu

    def __rmcolorizer(self):
        if not self.color:
            return
        self.color.removecolors()
        self.per.removefilter(self.color)
        self.color = None

    def __addcolorizer(self):
        if self.color:
            return

        if hasattr(self.editwin, 'ispythonsource'): # PyShell window
            pychecker = self.editwin.ispythonsource
        else:
            pychecker = self.ispythonsource
        if pychecker(self.io.filename):
            self.color = self.editwin.ColorDelegator()

        # can add more colorizers here...
        if self.color:
            self.per.removefilter(self.undo)
            self.per.insertfilter(self.color)
            self.per.insertfilter(self.undo)

    def __asktabwidth(self):
        return self.editwin.askinteger(
            "Tab width",
            "Columns per tab? (2-16)",
            parent=self.text,
            initialvalue=self.editwin.indentwidth,
            minvalue=2,
            maxvalue=16) or self.editwin.tabwidth

    # Make string that displays as n leading blanks.
    def __make_blanks(self, n):
        if self.editwin.usetabs:
            ntabs, nspaces = divmod(n, self.editwin.tabwidth)
            return '\t' * ntabs + ' ' * nspaces
        else:
            return ' ' * n

    def __get_region(self):
        text = self.text
        first, last = self.get_selection_indices()
        if first and last:
            head = text.index(first + " linestart")
            tail = text.index(last + "-1c lineend +1c")
        else:
            head = text.index("insert linestart")
            tail = text.index("insert lineend +1c")
        chars = text.get(head, tail)
        lines = chars.split("\n")
        return head, tail, chars, lines

    def __set_region(self, head, tail, chars, lines):
        text = self.text
        newchars = "\n".join(lines)
        if newchars == chars:
            text.bell()
            return
        text.tag_remove("sel", "1.0", "end")
        text.mark_set("insert", head)
        text.undo_block_start()
        text.delete(head, tail)
        text.insert(head, newchars)
        text.undo_block_stop()
        text.tag_add("sel", head, "insert")

    # Delete from beginning of line to insert point, then reinsert
    # column logical (meaning use tabs if appropriate) spaces.
    def __reindent_to(self, column):
        text = self.text
        text.undo_block_start()
        if text.compare("insert linestart", "!=", "insert"):
            text.delete("insert linestart", "insert")
        if column:
            text.insert("insert", self.__make_blanks(column))
        text.undo_block_stop()

    # Tk implementations of "virtual text methods" -- each platform
    # reusing IDLE's support code needs to define these for its GUI's
    # flavor of widget.

    # Is character at text_index in a Python string?  Return 0 for
    # "guaranteed no", true for anything else.  This info is expensive
    # to compute ab initio, but is probably already known by the
    # platform's colorizer.
    def __is_char_in_string(self, text_index):
        if self.color:
            # Return true iff colorizer hasn't (re)gotten this far
            # yet, or the character is tagged as being in a string
            return self.text.tag_prevrange("TODO", text_index) or \
                   "STRING" in self.text.tag_names(text_index)
        else:
            # The colorizer is missing: assume the worst
            return 1


import re

class OutputPage(EditorPage):
    rmenu_specs = [
        ("Go to file/line", "<<goto-file-line>>"),
    ]

    file_line_pats = [
        r'file "([^"]*)", line (\d+)',
        r'([^\s]+)\((\d+)\)',
        r'([^\s]+):\s*(\d+):',
    ]

    file_line_progs = None

    def __init__(self, parent_frame, editwin, title=None, **kwargs):
        EditorPage.__init__(self, parent_frame, editwin, title, **kwargs)
        self.text.bind("<<goto-file-line>>", self.goto_file_line)

    def ispythonsource(self, filename):
        # No colorization needed
        return 0

    def short_title(self):
        return "Output"

    def maybesave(self):
        # Override base class method -- don't ask any questions
        if self.get_saved():
            return "yes"
        else:
            return "no"

    def goto_file_line(self, event):
        text = self.text
        if self.file_line_progs is None:
            l = []
            for pat in self.file_line_pats:
                l.append(re.compile(pat, re.IGNORECASE))
            self.file_line_progs = l

        line = text.get("insert linestart", "insert lineend")
        result = self._file_line_helper(line)
        if not result:
            # Try the previous line.  This is handy e.g. in tracebacks,
            # where you tend to right-click on the displayed source line
            line = text.get("insert -1line linestart", "insert -1line ineend")
            result = self._file_line_helper(line)
            if not result:
                tkMessageBox.showerror("No special line",
                    "The line you point at doesn't look like "
                    "a valid file name followed by a line number.",
                    master=text)
                return
        filename, lineno = result
        self._open_file_at_line(filename, lineno)

    def _file_line_helper(self, line):
        for prog in self.file_line_progs:
            m = prog.search(line)
            if m:
                break
        else:
            return None

        filename, lineno = m.group(1, 2)
        try:
            f = open(filename, "r")
            f.close()
        except IOError:
            return None
        try:
            return filename, int(lineno)
        except TypeError:
            return None

    def _open_file_at_line(self, filename, lineno):
        edit = self.editwin.flist.open(filename)
        edit.gotoline(lineno)
