import Tkinter, Pmw
from mglutil.util.callback import CallbackManager, CallBackFunction
from mglutil.gui.InputForm.Tk.moveTk import MovableWidget
from mglutil.gui.BasicWidgets.Tk import customizedWidgets
cw = customizedWidgets

tkinterWidgets = [('Button',), ('Checkbuttons',), ('Entry',), ('Label',),
                  ('Text',), ('ListBox',), ('Frame',), ('Scale',)]
pmwWidgets = [('ComboBox',), ( 'EntryField',), ( 'RadioSelect',),
              ( 'ButtonBox',), ('Counter',), ( 'Group',),
              ( 'LabeledWidget',), ( 'ScrolledCanvas',), ('ScrolledField',),
              ('ScrolledListBox',), ( 'ScrolledText',), ('TimeCounter',)]
customizedWidgets = [('ListChooser',), ('ExtendedSliderWidget',),
                     ( 'SliderWidget',)]
              
class InteractiveInputFormBuilder(Tkinter.Frame):
    """ This class implements an interactive inputform builder that provides:
    - a frame inside of which you can widgets and edit them.
    - a set of widgets that can be added to the frame
    - mechanism to modify these widgets
    - mechanism to save the created inputform.
    """
    def __init__(self, name='NoName', master = None, **kw):
        """
        instance <- InteractiveInputFormBuilder(name='NoName, master=None, **kw'
        kw can contain any """
        # Editor Name
        self.name = name
        # List of the ItemDescr object for each widgets added to the
        # inputform
        self.widgets = []
        self.containers = {}


##          self.actionCallback = {
##              'onAddNode' : CallbackManager(),
##              'onDeleteNodes' : CallbackManager(),
##              'onAddConnection' : CallbackManager(),
##              'onDeleteConnections' : CallbackManager(),
##              'onSelectNodes' : CallbackManager(),
##              'onDeselectNodes' : CallbackManager(),
##              'onSelectConnections' : CallbackManager(),
##              'onDeselectConnections' : CallbackManager(),
##              }

        Tkinter.Frame.__init__(self, master)
        Tkinter.Pack.config(self, expand=1, fill=Tkinter.BOTH)

        mainFrame = Tkinter.Frame(bg='white',width=400, height=400)
        self.containers['mainFrame'] = mainFrame
        mainFrame.pack(side=Tkinter.LEFT)

        self.createMenus()
        self.widgetArea = Tkinter.Frame(borderwidth=2, relief = 'sunken')
##          self.widgetArea.pack(side=Tkinter.RIGHT, expand=1, fill ='both')

    def createMenus(self):
        self.mBar = Tkinter.Frame(self, relief=Tkinter.RAISED, borderwidth=2)
        self.mBar.pack(fill=Tkinter.X)
        self.menuButtons = {}
        self.makeFileMenu()
        self.makeEditMenu()
        apply( self.mBar.tk_menuBar, self.menuButtons.values() )
	self.title = Tkinter.Label(self.mBar, text=self.name)
	self.title.pack(side=Tkinter.RIGHT)

    def makeFileMenu(self):
        File_button = Tkinter.Menubutton(self.mBar, text='File', underline=0)
        self.menuButtons['File'] = File_button
        File_button.pack(side=Tkinter.LEFT, padx="1m")
        File_button.menu = Tkinter.Menu(File_button)
##          File_button.menu.add_command(label='New...', underline=0, 
##                                       command=self.loadFile)
        
        File_button.menu.add_command(label='Load...', underline=0, 
                                     command=self.loadFile)
        File_button.menu.entryconfig(1, state=Tkinter.DISABLED)
        File_button.menu.add_command(label='Save...', underline=0, 
                                     command=self.saveFile)
        File_button.menu.entryconfig(2, state=Tkinter.DISABLED)
        File_button.menu.add_command(label='Quit', underline=0, 
                                     command=self.exit)
        File_button['menu'] = File_button.menu


    def exit(self, event=None):
        self.master.destroy()

        
    def makeEditMenu(self):
        Edit_button = Tkinter.Menubutton(self.mBar, text='Edit', underline=0)
        self.menuButtons['Edit'] = Edit_button
        Edit_button.pack(side=Tkinter.LEFT, padx="1m")
        Edit_button.menu = Tkinter.Menu(Edit_button)

        Edit_button.menu.add('command', label="Undo", command=self.undo)
        Edit_button.menu.entryconfig(1, state=Tkinter.DISABLED)

        Edit_button.menu.add('command', label="Add widget",
                             command=self.buildWidgetsFrame)
        Edit_button.menu.entryconfig(2)
        
        
        # and these are just for show. No "command" callbacks attached.
        Edit_button.menu.add_command(label="Cut")
        Edit_button.menu.entryconfig(3, state=Tkinter.DISABLED)
        Edit_button.menu.add_command(label="Copy")
        Edit_button.menu.entryconfig(4, state=Tkinter.DISABLED)
        Edit_button.menu.add_command(label="Paste")
        Edit_button.menu.entryconfig(5, state=Tkinter.DISABLED)
        Edit_button.menu.add_command(label="Delete", command=self.delete)
        Edit_button.menu.entryconfig(5, state=Tkinter.DISABLED)

        # set up a pointer from the file menubutton back to the file menu
        Edit_button['menu'] = Edit_button.menu

    def addWidget(self,event=None):
        wtype = self.guiTK.get()+'.'+self.chosenWidgets.get()[0]
        container = self.containers[self.chosenContainers.get()[0]]
        widget = apply( eval(wtype), (container,) )
        b = MovableWidget(widget, 0, 0)
        self.widgets.append(b)
    
    def buildWidgetsFrame(self, event=None):
        self.widgetsFrame = Tkinter.Frame(master = self.widgetArea)
        self.widgetsFrame.pack()
        self.choices = {'Tkinter':tkinterWidgets,'Pmw':pmwWidgets,
                   'customizedWidgets':customizedWidgets}
        # Combox with the different widgets libraries.
        self.guiTK = Pmw.ComboBox(self.widgetsFrame,
                                  label_text = 'Choose a Gui ToolKit',
                                  labelpos = 'nw',
                                  selectioncommand = self.guiType,
                                  scrolledlist_items = self.choices.keys() )
        
        self.guiTK.grid(padx = 8, pady = 8)
        self.guiTK.selectitem('Tkinter', setentry=1)
        self.guiTK.grid(row=1, column=0, padx = 8,
                        pady = 8)
        # ListChooser with the widgets contained in the chosen library.
        self.chosenWidgets = cw.ListChooser(self.widgetsFrame, mode='single',
                                      title='Choose a widget',
                                      entries = self.choices['Tkinter'],
                                      lbwcfg ={'exportselection':0})
        self.chosenWidgets.grid(row=2, column=0, padx = 8,
                        pady = 8)
        # ListChooser with the container already created
        entries = map(lambda x: (x,), self.containers.keys())

        self.chosenContainers = cw.ListChooser(self.widgetsFrame,
                                               mode='single',
                                         title='Choose a container',
                                         entries = entries,
                                         lbwcfg ={'exportselection':0})
        self.chosenContainers.grid(row=2, column=1, padx = 8,
                                   pady = 8)

##          # Posx EntryField
##          self.posxW = Pmw.EntryField(self.widgetsFrame,
##                                     label_text = 'Enter a x position',
##                                     labelpos = 'w',
##                                     validate = {'validator':'integer'})
        
##          self.posxW.grid(row=3, column=0, padx = 8,
##                          pady = 8)
##          # Posy EntryField
##          self.posyW = Pmw.EntryField(self.widgetsFrame,
##                                     label_text = 'Enter a y position',
##                                     labelpos = 'w',
##                                     validate = {'validator':'integer'})
##          self.posyW.grid(row=3, column=1, padx = 8,
##                          pady = 8)
                             
        # Add button to add the selected widget in the frame window.
        add = Tkinter.Button(self.widgetsFrame,
                             text = 'Add', command = self.addWidget)
        add.grid(row=4, column=0)

        # Edit button:
        self.editVar = Tkinter.IntVar()
        edit = Tkinter.Checkbutton(self.widgetsFrame,
                                   text = 'Edit Widgets',
                                   variable=self.editVar)
        edit.bind('<ButtonPress>', self.edit_cb)
        edit.grid(row=4, column=1)
        
        self.widgetArea.pack(side=Tkinter.RIGHT, expand=1, fill ='both')
        
        

    def edit_cb(self, event=None):
        print self.editVar.get()
        if self.editVar.get() == 0:
            for b in self.widgets:
                b.widget.bind('<Right>', b.moveRight)
                b.widget.bind('<Left>', b.moveLeft)
                b.widget.bind('<Up>', b.moveUp)
                b.widget.bind('<Down>', b.moveDown)
                b.widget.bind('<ButtonPress-1>', b.recordPosition)
                b.widget.bind('<B1-Motion>', b.moveButton)
                b.widget.bind('<ButtonRelease-1>', b.endMove)
        else:
            for b in self.widgets:
                b.widget.unbind('<Right>')
                b.widget.unbind('<Left>')
                b.widget.unbind('<Up>')
                b.widget.unbind('<Down>')
                b.widget.unbind('<ButtonPress-1>')
                b.widget.unbind('<B1-Motion>')
                b.widget.unbind('<ButtonRelease-1>')
            
    
    def guiType(self, event=None):
        self.chosenWidgets.clear()
        map(self.chosenWidgets.add, self.choices[self.guiTK.get()])
        

    def loadFile(self, event=None):
        pass


    def saveFile(self, event=None):
        pass

    def newFile(sele, event=None):
        pass

    def delete(self, event=None):
        pass


    def undo(self, event =None):
        pass
    
