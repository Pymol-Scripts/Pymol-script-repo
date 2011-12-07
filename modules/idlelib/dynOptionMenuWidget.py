"""
OptionMenu widget modified to allow dynamic menu reconfiguration
and setting of highlightthickness
"""
from Tkinter import OptionMenu, Menu
from Tkinter import _setit
import copy

from configHandler import idleConf
TTK = idleConf.GetOption('main', 'General', 'use-ttk', type='int')
if TTK:
    from ttk import OptionMenu

class DynOptionMenu(OptionMenu):
    """Unlike OptionMenu, our kwargs can include highlightthickness"""
    def __init__(self, master, variable, value, *values, **kwargs):
        #get a copy of kwargs before OptionMenu.__init__ munges them
        kwargsCopy=copy.copy(kwargs)
        if 'highlightthickness' in kwargs.keys():
            del(kwargs['highlightthickness'])
        self.command=kwargs.get('command')
        self.variable=variable

        OptionMenu.__init__(self, master, variable, value, *values, **kwargs)
        self.config(highlightthickness=kwargsCopy.get('highlightthickness'))

    def SetMenu(self, valueList, value=None):
        """
        clear and reload the menu with a new set of options.
        valueList - list of new options
        value - initial value to set the optionmenu's menubutton to
        """
        if TTK:
            self.set_menu(value, *valueList)
        else:
            menu = self['menu']
            menu.delete(0,'end')
            for item in valueList:
                menu.add_command(label=item,
                    command=_setit(self.variable,item,self.command))
            if value:
                self.variable.set(value)
