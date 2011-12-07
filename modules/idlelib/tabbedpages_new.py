"""Classes exported:

TabbedPageSet -- A custom ttk.Notebook used by IDLE.
"""
from ttk import Frame, Notebook

from tabbedpages import InvalidNameError, AlreadyExistsError

class FramePage(object):
    def __init__(self, notebook):
        self.frame = Frame(notebook)

class TabbedPageSet(Notebook):
    """
    Pages may be accessed through the 'pages' attribute, which is a dictionary
    of pages, using the name given as the key. A page is an instance of a
    subclass of ttk's Frame widget.

    Pages may be added or removed at any time using the add_page() and
    remove_page() methods.
    """

    def __init__(self, master, page_names=None, **kw):
        """Constructor arguments:

        page_names -- A list of strings, each will be the dictionary key to a
        page's widget, and the name displayed on the page's tab. Should be
        specified in the desired page order. The first page will be the default
        and first active page. If page_names is None or empty, the
        TabbedPageSet will be initialized empty.
        """
        Notebook.__init__(self, master, **kw)

        self.pages = {}
        page_names = page_names or ()
        for name in page_names:
            self.add_page(name)

    def update_tabtitle(self, tab, newtitle):
        """Update tab title to newtitle."""
        currpage = self.pages[tab.title].frame
        old = tab.title

        # resolve title duplicate
        if newtitle in self.pages and currpage != self.pages[newtitle].frame:
            # newtitle is already present, and the current tab is not the
            # one who owns it
            count = 1
            temptitle = newtitle
            while temptitle in self.pages:
                if currpage == self.pages[temptitle].frame:
                    break
                temptitle = "%s #%d" % (newtitle, count)
                count += 1
            newtitle = temptitle

        tab.title = newtitle
        self.pages[newtitle] = self.pages.pop(old)
        self.tab(currpage, text=newtitle)

    def add_page(self, page_name):
        """Add a new page with the name given in page_name."""
        if not page_name:
            raise InvalidNameError("Invalid TabPage name: '%s'" % page_name)
        if page_name in self.pages:
            raise AlreadyExistsError(
                "TabPage named '%s' already exists" % page_name)

        fpage = FramePage(self)
        self.pages[page_name] = fpage
        self.add(fpage.frame, text=page_name, padding=6)

        # workaround for bug #1878298 at tktoolkit sf bug tracker
        self.event_generate('<Expose>')

        return fpage

    def remove_page(self, page_name):
        """Remove page_name from the notebook."""
        if not page_name in self.pages:
            raise KeyError("No such TabPage: '%s" % page_name)

        self.forget(self.index(self.pages[page_name].frame))
        del self.pages[page_name]

        # workaround for bug #1878298 at tktoolkit sf bug tracker
        self.event_generate('<Expose>')

    def last_page(self):
        """Return the last page in the notebook."""
        return self.pages[self.tab(self.index('end') - 1)['text']]

if __name__ == '__main__':
    from Tkinter import Tk
    from Tkconstants import TOP, BOTH
    from ttk import Label, Entry, Button, Style
    # test dialog
    root=Tk()
    style = Style()
    style.configure('C.TLabel', padding=20)
    tabPage=TabbedPageSet(root, page_names=['Foobar','Baz'])
    tabPage.pack(side=TOP, expand=True, fill=BOTH)
    Label(tabPage.pages['Foobar'].frame, text='Foo', style='C.TLabel').pack()
    Label(tabPage.pages['Foobar'].frame, text='Bar', style='C.TLabel').pack()
    Label(tabPage.pages['Baz'].frame, text='Baz').pack()
    entryPgName=Entry(root)
    buttonAdd=Button(root, text='Add Page',
            command=lambda:tabPage.add_page(entryPgName.get()))
    buttonRemove=Button(root, text='Remove Page',
            command=lambda:tabPage.remove_page(entryPgName.get()))
    labelPgName=Label(root, text='name of page to add/remove:')
    buttonAdd.pack(padx=5, pady=5)
    buttonRemove.pack(padx=5, pady=5)
    labelPgName.pack(padx=5)
    entryPgName.pack(padx=5)
    root.mainloop()
