#
# Last modified on Mon Nov 26 19:00:44 PST 2001 by lindy
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/autopilotCommands.py,v 1.3 2005/01/12 21:53:52 rhuey Exp $
#

"""
This Module facilitates the launching of large numbers of AutoDock jobs
"""


import string, Tkinter
from ViewerFramework.VFCommand import CommandGUI
from Pmv.mvCommand import MVCommand
from mglutil.gui.InputForm.Tk.gui import InputFormDescr

##  class ADParamSearch(MVCommand):
##      """facilitate parameter search"""

##      def onAddCmdToViewer(self):
##          self.dpo = DockingParameters.DockingParameters()
##          self.num_params = Tkinter.StringVar()
##      self.num_params.set('3')
##          self.base_dpf_name = Tkinter.StringVar()
##      self.parameters = []
##      self.first_values = []
##      self.last_value = []
##      self.step_value = []
##      for i in range(10):
##          self.first_values.append(Tkinter.StringVar())
##          self.last_values.append(Tkinter.StringVar())
##          self.step_values.append(Tkinter.StringVar())
##          self.parameter.append(Tkinter.StringVar())


##      def guiCallback(self):
##          #called each time the param_search button is pressed
##          print "called guiCallback in ADParamSearch"
##          if not hasattr(self, 'form'):
##              # create the ifd
##              self.buildGUI()
##          else:
##              # use the old ifd
##              self.form.deiconify()

##      def dismiss_cb(self):
##          self.form.withdraw()

##      def browse_cb(self):
##          base_dpf_name = self.vf.askFileOpen(types=[('base docking parameter file:', '*.dpf')], title = 'Select DPF')
##      if base_dpf_name:
##          # os.path.basename???
##          self.base_dpf_name.set(base_dpf_name)
##          #instantiate DPO
##          #read DPF
##          #set parameter lists accordingly


##      def doit_cb(self):
##          print "num_params:",  self.num_params.get()


##      def buildParamWids(self, num):
##      oldlen = len(self.paramWids)
##      diff = num - oldlen
##      if diff < 0:
##          for i in range(diff, 0, -1):
##          self.paramWids[i].grid_forget()
##      elif diff > 0:
##          self.paramWids[i].grid(self.grids[i])


##      def buildMenu(self, mB, nameList, varDict, oldvarDict, cmd):
##          for i in nameList:
##              if i not in varDict.keys():
##                  varDict[i]=Tkinter.IntVar()
##                  oldvarDict[i]=0
##          if hasattr(mB, 'menu'):
##              mB.menu.delete(1,'end')
##          else:
##              mB.menu = Tkinter.Menu(mB)
##              mB['menu']=mB.menu
##          #PACK all the entries
##          for i in varDict.keys():
##              mB.menu.add_radiobutton(label=i, var=varDict[i], command=cmd)

##      def buildParamMenu(self, paramButton, num, event=None):
##          #molMenubutton = self.ifd.entryByName['Mol List']['widget']
##          self.menuVars[num]={}
##          self.oldmenuVars[num]={}
##          paramNames = self.dpo.names
##      cmd = CallbackFunction(self.getParamVal, num)
##          self.buildMenu(paramButton,paramNames,self.menuVars[num],self.oldmolVar,cmd)

##      def buildEntry(self, i):
##          print 'in buildEntry'
    

##      def getEntries(self, num):
##      #first forget all the entries
##      for i in self.paramEntries:
##          i.grid_forget()
##      #then build new ones
##      for i in range(num):
##          print 'build ', i, 'th entry here'
##          self.buildEntry(i)

    
    

##      def buildGUI(self):
##          """Create the Parameter Search input form
##          """
##      self.ifd = ifd = InputFormDescr(title='Parameter Search')
##      # append new dicts to the ifd list
##      # one-widget-one-dict (generally)
##      ifd.append({ 'widgetType'   : Tkinter.Entry,
##               'label'        : "Number of Paramters",
##               'textvariable' : self.num_params,
##               'gridcfg' : { 'sticky' : Tkinter.W,
##                     'columnspan' : 2 }
##               })
##      ifd.append({ 'widgetType'   : Tkinter.Entry,
##               'label'        : "Base DPF",
##               'textvariable' : self.base_dpf_name,
##               'gridcfg' : { 'sticky' : Tkinter.E,
##                     'columnspan' : 2,
##                     'row' : -1,
##                     'column' : 2 }
##               })
##      ifd.append({ 'widgetType'   : Tkinter.Button,
##               'text'         : "Browse",
##               'command'      : self.browse_cb,
##               'gridcfg' : { 'sticky' : Tkinter.E,
##                     'columnspan' : 1,
##                     'row' : -1,
##                     'column' : 4 }
##               })
##      #
##      # Throw up some labels for the widgets to follow
##      #
##      ifd.append({ 'widgetType'   : Tkinter.Label,
##               'text'         : "Parameter",
##               'gridcfg' : { 'sticky' : Tkinter.W + Tkinter.E,
##                     'columnspan' : 1 }
##               })
##      ifd.append({ 'widgetType'   : Tkinter.Label,
##               'text'         : "First Value",
##               'gridcfg' : { 'sticky' : Tkinter.W + Tkinter.E,
##                     'columnspan' : 1,
##                     'row' : -1,
##                     'column' : 1 }
##               })
##      ifd.append({ 'widgetType'   : Tkinter.Label,
##               'text'         : "Last Value",
##               'gridcfg' : { 'sticky' : Tkinter.W + Tkinter.E,
##                     'columnspan' : 1,
##                     'row' : -1,
##                     'column' : 2 }
##               })
##      ifd.append({ 'widgetType'   : Tkinter.Label,
##               'text'         : "Step Value",
##               'gridcfg' : { 'sticky' : Tkinter.W + Tkinter.E,
##                     'columnspan' : 1,
##                     'row' : -1,
##                     'column' : 3 }
##               })
##      #
##      # these widgets define the parameter and range
##      #
##      ifd.append({ 'widgetType'   : Pmw.ScrolledListBox,
##               'textvariable' : self.parameters,
##               'wcfg' : {'items' : self.dpo.docking_parameter_list,
##                     'listbox_height' : 1},
##               'gridcfg' : { 'sticky' : Tkinter.W,
##                     'columnspan' : 1 }
##               })
##      ifd.append({ 'widgetType'   : Tkinter.Entry,
##               'textvariable' : self.first_value,
##               'gridcfg' : { 'sticky' : Tkinter.W,
##                     'columnspan' : 1,
##                     'row' : -1,
##                     'column' : 1 }
##               })
##      ifd.append({ 'widgetType'   : Tkinter.Entry,
##               'textvariable' : self.last_value,
##               'gridcfg' : { 'sticky' : Tkinter.W,
##                     'columnspan' : 1,
##                     'row' : -1,
##                     'column' : 2 }
##               })
##      ifd.append({ 'widgetType'   : Tkinter.Entry,
##               'textvariable' : self.step_value,
##               'gridcfg' : { 'sticky' : Tkinter.W,
##                     'columnspan' : 1,
##                     'row' : -1,
##                     'column' : 3 }
##               })

##      #
##      # 
##      #
##      ifd.append({ 'widgetType'   : Tkinter.Button,
##               'text'         : "Doit",
##               'command'      : self.doit_cb,
##               'gridcfg' : { 'sticky' : Tkinter.W + Tkinter.E,
##                     'columnspan' : 1 }
##               })
##      ifd.append({ 'widgetType'   : Tkinter.Button,
##               'text'         : "Dismiss",
##               'command'      : self.dismiss_cb,
##               'gridcfg' : { 'sticky' : Tkinter.W + Tkinter.E,
##                     'columnspan' : 1,
##                     'row' : -1,
##                     'column' : 1 }
##               })

##      # after the dict is built, do this
##          self.form = self.vf.getUserInput(self.ifd, modal= 0, blocking= 0)


##      def doit():
##          print 'in doit'

##  ADParamSearchGUI = CommandGUI()
##  ADParamSearchGUI.addMenuCommand('AutoToolsBar', 'AutoPilot', 'Parameter Search')

#
# Menu text
#
menuText = {}
menuText['AutoPilotMB'] = 'Pilot'
menuText['TemplateMB'] = 'Template'
menuText['DatabaseMB'] = 'Database'
menuText['ClusterMB'] = 'Cluster'

#
# Cluster command and GUICommand
#
class ADpilot_Cluster(MVCommand):
    """
    """

    #
    # the signiture for __call__ must be compatible with that of doit!!!
    #
    def __call__(self, **kw):
        """This method is invoked by saying self.vf.ADpilot_Cluster
        in the Python Shell.
        """
        apply(self.doitWrapper, (), kw)


    def onAddCmdToViewer(self):
        """Initialization.

        This method is called by vf.addCommand (below)
        """
        pass


    def doit(self):
        """The real work gets done here
        """
        print 'in doit'


class ADpilot_ClusterGUICommand(MVCommand):
    """A template GUI command"""

    #
    # signiture for __call__ must be compatible with that of doit!!!
    #
    def __call__(self, **kw):
        """This method is invoked by saying self.vf.ADpilot_ClusterGUICommand
        in the Python Shell.
        """
        apply(self.doitWrapper, (), kw)


    def onAddCmdToViewer(self):
        """Initialization.

        This method is called by vf.addCommand (below)
        """
        pass


    def buildFormDescr(self, formName):
        if formName == 'template':
            idf = self.idf = InputFormDescr(title = self.name)
            # append to idf here
            idf.append( {'name': 'template',
                         'widgetType': Tkinter.Button,
                         'wcfg': { 'text': 'template button label'},
                         'gridcfg': { 'sticky' : 'we'}})
            # only return idf for template
            return idf


    def guiCallback(self):
        """
        """
        idf_dict = self.showForm('template')
    
        # get the kw key:value pairs from the gui
        apply(self.doitWrapper, (), idf_dict)


    def doit(self, **kw):
        """The GUI command calls the real command here
        """
        # process keyword args here
        self.vf.ADpilot_Cluster()

ADpilot_ClusterGUI = CommandGUI()
ADpilot_ClusterGUI.addMenuCommand('AutoToolsBar',
                                   menuText['AutoPilotMB'],
                                   menuText['ClusterMB'],
                                   cascadeName = menuText['DatabaseMB'] )


#
# Template command and GUICommand
#
class ADpilot_Template(MVCommand):
    """A template command"""
    #
    # the signiture for __call__ must be compatible with that of doit!!!
    #
    def __call__(self, **kw):
        """This method is invoked by saying self.vf.ADpilot_Template
        in the Python Shell.
        """
        apply(self.doitWrapper, (), kw)


    def onAddCmdToViewer(self):
        """Initialization.

        This method is called by vf.addCommand (below)
        """
        pass



    def doit(self):
        """The real work gets done here
        """
        print 'in doit'



class ADpilot_TemplateGUICommand(MVCommand):
    """A template GUI command"""

    #
    # signiture for __call__ must be compatible with that of doit!!!
    #
    def __call__(self, **kw):
        """This method is invoked by saying self.vf.ADpilot_TemplateGUICommand
        in the Python Shell.
        """
        apply(self.doitWrapper, (), kw)


    def onAddCmdToViewer(self):
        """Initialization.

        This method is called by vf.addCommand (below)
        """
        pass


    def buildFormDescr(self, formName):
        if formName == 'template':
            idf = self.idf = InputFormDescr(title = self.name)
            # append to idf here
            idf.append( {'name': 'template',
                         'widgetType': Tkinter.Button,
                         'wcfg': { 'text': 'template button label'},
                         'gridcfg': { 'sticky' : 'we'}})
            # only return idf for template
            return idf


    def guiCallback(self):
        """
        """
        idf_dict = self.showForm('template')
    
        # get the kw key:value pairs from the gui
        apply(self.doitWrapper, (), idf_dict)


    def doit(self, **kw):
        """The GUI command calls the real command here
        """
        # process keyword args here
        self.vf.ADpilot_Template()

ADpilot_TemplateGUI = CommandGUI()
ADpilot_TemplateGUI.addMenuCommand('AutoToolsBar',
                                   menuText['AutoPilotMB'],
                                   menuText['TemplateMB'])


# A convenient list of dictionaries with arguments for vf.addCommand
commandList = [
    # Cluster Command
    {'name': 'ADpilot_Cluster',
     'cmd': ADpilot_Cluster(),
     'gui': None},
    {'name': 'ADpilot_ClusterGC',
     'cmd': ADpilot_ClusterGUICommand(),
     'gui': ADpilot_ClusterGUI},

    # Template Command
    {'name': 'ADpilot_Template',
     'cmd': ADpilot_Template(),
     'gui': None},
    {'name': 'ADpilot_TemplateGC',
     'cmd': ADpilot_TemplateGUICommand(),
     'gui': ADpilot_TemplateGUI},
    ]


def initModule(vf):
    """Add each command in the modules's commandList into ViewerFramework
    """
    for dict in commandList:
        vf.addCommand(dict['cmd'], dict['name'], dict['gui'])

# To load this module, type into the adt python interpreter:
# self.loadModule('autopilotCommands', 'AutoDockTools', log = 0)

