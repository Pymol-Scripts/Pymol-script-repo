'''
Described at: http://www.pymolwiki.org/index.php/rendering_plugin

Plugin contributed by Michael Lerner (mglerner@gmail.com)
'''

from __future__ import print_function

import sys
if sys.version_info[0] < 3:
    import Tkinter
else:
    import tkinter as Tkinter

import Pmw

try:
    import pymol
    REAL_PYMOL = True
except ImportError:
    REAL_PYMOL = False

    class pymol:

        class cmd:

            def ray(*args):
                print("ray", args)

            def png(*args, **kwargs):
                print("png", args, kwargs)

            def draw(*args, **kwargs):
                print("draw", args, kwargs)
            ray = staticmethod(ray)
            png = staticmethod(png)
            draw = staticmethod(draw)


def __init__(self):
    self.menuBar.addmenuitem('Plugin', 'command',
                             'Launch Rendering Plugin',
                             label='Rendering...',
                             command=lambda s=self: RenderPlugin(s))


class RenderPlugin:

    def __init__(self, app):
        self.parent = app.root
        self.dialog = Pmw.Dialog(self.parent,
                                 buttons=('Ray', 'Draw', 'Exit'),
                                 title = 'PyMOL Rendering Plugin',
                                 command = self.execute)
        self.dialog.withdraw()
        Pmw.setbusycursorattributes(self.dialog.component('hull'))

        w = Tkinter.Label(self.dialog.interior(),
                          text='PyMOL Rendering Plugin\nMichael Lerner, 2006 - www.umich.edu/~mlerner/PyMOL\nNOTE: make sure this window is not on top of the PyMOL window.',
                          background='black',
                          foreground='white',
                          )
        w.pack(expand=1, fill='both', padx=4, pady=4)

        #self.notebook = Pmw.NoteBook(self.dialog.interior())
        # self.notebook.pack(fill='both',expand=1,padx=10,pady=10)

        # Set up the Main page
        #page = self.notebook.add('Main')
        #group = Pmw.Group(page,tag_text='Main options')

        group = Pmw.Group(self.dialog.interior())

        group.pack(fill='both', expand=1, padx=10, pady=5)
        self.filename = Pmw.EntryField(group.interior(),
                                       labelpos='w',
                                       label_text='Filename',
                                       value='picture.png',
                                       )
        self.height = Pmw.EntryField(group.interior(), labelpos='w',
                                     label_text='Height:',
                                     value=str(4.0),
                                     validate={'validator': 'real',
                                                 'min': 0, }
                                     )

        self.width = Pmw.EntryField(group.interior(), labelpos='w',
                                    label_text='Width:',
                                    value=str(4.0),
                                    validate={'validator': 'real',
                                                'min': 0, }
                                    )
        self.units = Pmw.OptionMenu(group.interior(), labelpos='w',
                                    label_text='Units',
                                    items=('inch', 'cm',),
                                    )
        self.dpi = Pmw.EntryField(group.interior(), labelpos='w',
                                  label_text='DPI:',
                                  value=str(300),
                                  validate={'validator': 'real',
                                              'min': 0, }
                                  )
        entries = (self.height, self.width, self.units, self.filename, self.dpi)
        for entry in entries:
            # entry.pack(side='left',fill='both',expand=1,padx=4) # side-by-side
            entry.pack(fill='x', expand=1, padx=4, pady=1)  # vertical
        # self.notebook.setnaturalsize()
        self.showAppModal()

    def showAppModal(self):
        #self.dialog.activate(geometry = 'centerscreenalways', globalMode = 'nograb')
        self.dialog.show()
        #self.dialog.activate(geometry = 'centerscreenalways')

    def execute(self, result):
        dpi = float(self.dpi.getvalue())
        if self.units.getvalue() == 'cm':
            dp_unit = dpi / 2.54
        else:
            dp_unit = dpi
        h = int(float(self.height.getvalue()) * dp_unit)
        w = int(float(self.width.getvalue()) * dp_unit)
        if result == 'Ray':
            pymol.cmd.ray(w, h)
            pymol.cmd.png(self.filename.getvalue(), dpi=dpi)
        elif result == 'Draw':
            pymol.cmd.draw(w, h)
            pymol.cmd.png(self.filename.getvalue(), dpi=dpi)
        else:
            #
            # Doing it this way takes care of clicking on the x in the top of the
            # window, which as result set to None.
            #
            if __name__ == '__main__':
                #
                # dies with traceback, but who cares
                #
                self.parent.destroy()
            else:
                self.dialog.withdraw()


# Create demo in root window for testing.
if __name__ == '__main__':
    class App:

        def my_show(self, *args, **kwargs):
            pass
    app = App()
    app.root = Tkinter.Tk()
    Pmw.initialise(app.root)
    app.root.title('Some Title')

    widget = RenderPlugin(app)
    exitButton = Tkinter.Button(app.root, text='Exit', command=app.root.destroy)
    exitButton.pack()
    app.root.mainloop()
