import Pmw, Tkinter
# Pmw_0_8_5

class VFRadioSelect(Pmw.RadioSelect):
    """ This class is a wrapper for the Pmw.RadioSlect class.
    It implements :
    - get method: to get the current value.
    - set method: to set the radioslect to a given value.

    Pmw.RadioSelect constructor arguments:
        buttontype='button',
        command=None,
        labelmargin=0,
        labelpos=None,
        orient='vertical',
        padx=5, pady=5,
        selectmode='single'

    VFPmw.VFRadioSelect new arguments:
        buttonText = [] list of the label for each button contained in the
                        RadioSelect widget.
    """
    def __init__(self, parent, buttonText=[], **kw):
        """ New constructor argument:
        - buttonText: label for the button in the RadioSelect"""

        # 1- Define optionDefs
        optiondefs = (())
        # 2 Call the self.defineOptions method
        self.defineoptions(kw, optiondefs)
        # 3- Call the constructor of the base class Pmw.EntryField
        Pmw.RadioSelect.__init__(self, parent)
        # 4- Call the initialiseoption method
        self.initialiseoptions(VFRadioSelect)

        for text in buttonText:
            self.add(text)

    def set(self, buttonName):
        """ Check the button corresponding the buttonName"""
        self.invoke(buttonName)

    def get(self):
        """Get the name of button curent selected and return it"""
        curbuttonName = self.getcurselection()
        return curbuttonName


class VFButtonBox(Pmw.ButtonBox):
    """
    This class is a wrapper for the Pmw.ButtonBox class.
    the __init__method takes a new argument buttonText which is the list of
    the button text contained in the ButtonBox.
    It implements :
    - set method: to set the default button in the button box to a given
    value.

    Pmw.ButtonBox constructor arguments:
    buttontype (default Tkinter.Button)
    label_text: default ''
    labelmargin=0
    labelpos=None
    orient='vertical'
    padx=3, pady=3
    """
    def __init__(self, parent, buttonText=[], **kw):

        # 1- Define optionDefs
        optiondefs = (())
        # 2 Call the self.defineOptions method
        self.defineoptions(kw, optiondefs)
        # NO NEED of the new label_text argument already existing !
        # 3- Call the constructor of the base class Pmw.EntryField
        Pmw.ButtonBox.__init__(self, parent)
        # 4- Call the initialiseoption method
        self.initialiseoptions(VFButtonBox)

        # Here buttonText should be a list of tuple:
        # (buttonName, callBack)
        for text in buttonText:
            self.add(text)

    def set(self, buttonName):
        """
        call the command_callBack bound to the button buttonName.
        Problem if no command
        """
        # What happens here when setdefault is called should also call
        # the callback
        self.setdefault(buttonName)



class VFEntryField(Pmw.EntryField):
    """
    This class is a wrapper for the Pmw.EntryField class.
    Pmw.EntryField arguments/options and their default value:
        command          : None,
        errorbackground  : 'pink',
        extravalidators  : {},
        invalidcommand   : None,
        label_text        : '' ????
        labelmargin      : 0,
        labelpos         : None,
        modifiedcommand  : None,
        validate         : {},
        value            : ''

    New methods:
    - set method: to set the default entry in the counterentry field
      to a given value.
    - get method: to get the current value in the entry field
    """
    def __init__(self, parent=None, **kw):
        # 1- Define optionDefs
        optiondefs = (())
        # 2 Call the self.defineOptions method
        self.defineoptions(kw, optiondefs)
        # NO NEED of the new label_text argument already existing !
        # 3- Call the constructor of the base class Pmw.EntryField
        Pmw.EntryField.__init__(self, parent)
        # 4- Call the initialiseoption method
        self.initialiseoptions(VFEntryField)

        
    def set(self, defaultvalue):
        self.setentry(defaultvalue)



    def get(self):
        curValue = self.invoke()
        return curValue

class VFNoteBook(Pmw.NoteBook):
    """ The VFNoteBookR class is a wrapper for the Pmw.NoteBookR class."""
    def __init__(self, parent, listLevel, pagesDescr= {}, **kw):
        #self.listLevel = listLevel
        apply(Pmw.NoteBook.__init__, (self,parent), kw)
        print pagesDescr
        # Creates the pages
        self.pages = {}
        for level in listLevel:
            self.pages[level] = self.add(level)
            if not pagesDescr.has_key(level):
                continue
            descr = pagesDescr[level]
            widget = apply(descr['widgetType'], (self.pages[level],),
                          descr['wcfg'])
            apply(widget.grid, (), descr['gridcfg'])
        #print self.pages

        
