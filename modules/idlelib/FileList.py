import os
import tkMessageBox

from configHandler import idleConf

def _canonize(filename):
    if not os.path.isabs(filename):
        try:
            pwd = os.getcwd()
        except os.error:
            pass
        else:
            filename = os.path.join(pwd, filename)
    return os.path.normpath(filename)

class FileList:

    from EditorWindow import EditorWindow  # class variable, may be overridden
                                           # e.g. by PyShellFileList

    def __init__(self, root):
        self.root = root
        self.dict = {}
        self.inversedict = {}
        self.vars = {} # For EditorWindow.getrawvar (shared Tcl variables)

    def open(self, filename, action=None):
        assert filename
        filename = _canonize(filename)
        if os.path.isdir(filename):
            # This can happen when bad filename is passed on command line:
            tkMessageBox.showerror(
                "File Error",
                "%r is a directory." % (filename,),
                master=self.root)
            return None
        key = os.path.normcase(filename)
        if self.dict.has_key(key):
            edit = self.dict[key]
            edit.top.wakeup()
            return edit

        if action:
            # Don't create window, perform 'action', e.g. open in same window
            return action(filename=filename)

        elif idleConf.GetOption('main', 'EditorWindow', 'file-in-tab',
           default=1, type='bool'):
            # check if there is some PyShellEditorWindow running to open this
            # new tab
            for entry in self.inversedict:
                if hasattr(entry, 'set_breakpoint'):
                    # PyShellEditorWindow, good
                    entry.new_tab(filename=filename)
                    # select the last page created
                    entry.text_notebook.select(len(entry.text_notebook.pages)-1)
                    return entry
            else:
                # no PyShellEditorWindows, create one
                return self.EditorWindow(self, filename, key)

        else:
            return self.EditorWindow(self, filename, key)

    def gotofileline(self, filename, lineno=None):
        edit = self.open(filename)
        if edit is not None and lineno is not None:
            edit.gotoline(lineno)

    def new(self, filename=None):
        return self.EditorWindow(self, filename)

    def close_all_callback(self, event):
        for edit in self.inversedict.keys():
            reply = edit.close()
            if reply == "cancel":
                break
        return "break"

    def unregister_maybe_terminate(self, edit):
        try:
            key = self.inversedict[edit]
        except KeyError:
            print "Don't know this EditorWindow object.  (close)"
            return
        if key:
            del self.dict[key]
        del self.inversedict[edit]
        if not self.inversedict:
            self.root.quit()

    def filename_changed_edit(self, page, editwin):
        page.saved_change_hook(page.tab_initialized)
        try:
            key = self.inversedict[editwin]
        except KeyError:
            print "Don't know this EditorWindow object.  (rename)"
            return
        filename = page.io.filename
        if not filename:
            if key:
                del self.dict[key]
            self.inversedict[editwin] = None
            return
        filename = _canonize(filename)
        newkey = os.path.normcase(filename)
        if newkey == key:
            return
        if self.dict.has_key(newkey):
            conflict = self.dict[newkey]
            self.inversedict[conflict] = None
            tkMessageBox.showerror(
                "Name Conflict",
                "You now have multiple edit windows open for %r" % (filename,),
                master=self.root)
        self.dict[newkey] = editwin
        self.inversedict[editwin] = newkey
        if key:
            try:
                del self.dict[key]
            except KeyError:
                pass


def _test():
    import sys
    from Tkinter import Tk
    from EditorWindow import fixwordbreaks
    root = Tk()
    fixwordbreaks(root)
    root.withdraw()
    flist = FileList(root)
    if sys.argv[1:]:
        for filename in sys.argv[1:]:
            flist.open(filename)
    else:
        flist.new()
    if flist.inversedict:
        root.mainloop()

if __name__ == '__main__':
    _test()
