#!/usr/bin/env python

'''
See more here: http://www.pymolwiki.org/index.php/caver2

# CAVER Copyright Notice
# ============================
#
'''
from __future__ import division
from __future__ import generators

import os,math,platform
import Tkinter
from Tkinter import *
import Pmw
import distutils.spawn # used for find_executable
from pymol import cmd,selector
import sys,subprocess
from pymol.cmd import _feedback,fb_module,fb_mask,is_list,_cmd
from pymol.cgo import *
from chempy.models import Indexed
from chempy import Bond, Atom
import threading

#
# Global config variables
#

#win/linux ========================
# 1 for windows, 0 for linux

if sys.platform.startswith('win'):
    WINDOWZ = 1
elif sys.platform.startswith('linux'):
    WINDOWZ = 0
else:
    WINDOWZ = 0

VERS = "2"
JOPTS = "-Xmx400m"

VERSION = ".%s" % (VERS,)

#win/linux
if WINDOWZ:
#  OUTPUT_LOCATION = "C:\\Caver2_1\\Output"
  OUTPUT_LOCATION = os.getcwd()
else: #linux:
  OUTPUT_LOCATION = os.getcwd()

if WINDOWZ:
    if 'PYMOL_GIT_MOD' in os.environ:
        PYMOL_LOCATION = os.path.join(os.environ['PYMOL_GIT_MOD'],"Caver2_1_2","windows")
    else:
        PYMOL_LOCATION = "C:\\Program Files\\DeLano Scientific\\PyMol12"
else: #linux:
    if 'PYMOL_GIT_MOD' in os.environ:
        PYMOL_LOCATION = os.path.join(os.environ['PYMOL_GIT_MOD'],"Caver2_1_2","linux_mac")
    else:
        PYMOL_LOCATION = "directory/where/jar/with/plugin/is/located"

if WINDOWZ:
  LABEL_TEXT = "PyMOL location:"
else:
  LABEL_TEXT = "Directory with plugin .jar:"

# Python backward-compatibility...
try:
    True
except:
    True = 1
try:
    False
except:
    False = 0
#
# Cheap hack for testing purposes
#
try:
    import pymol
    REAL_PYMOL = True
except ImportError:
    REAL_PYMOL = False
    class pymol:
        class cmd:
            def load(self,name,sel=''):
                pass
            def get_names(self):
                return ['mol1','mol2','map1','map2']
            def get_type(self,thing):
                if thing.startswith('mol'):
                    return 'object:molecule'
                else:
                    return 'object:map'
                f.close()
        cmd = cmd()
    pymol = pymol()


# pridani do menu
def __init__(self):
        lbb = "Caver 2.1%s" % (VERSION,)
        self.menuBar.addmenuitem('Plugin', 'command',
                                 'Launch Caver 2.1'  + VERSION,
                                 label=lbb,
                                 command = lambda s=self: AnBeKoM(s))

defaults = {
    "startingpoint": ('-5','-5','-5'),
    "compute_command": 'Compute tunnels',
    "exit_command": 'Exit',
    "default_tunnels": '3',
    "surroundings" : '',
    "startingacids":('117','283','54'),
    "default_block": '10.0'
    }

class FileDialogButtonClassFactory:
    def get(fn,filter='*'):
        """This returns a FileDialogButton class that will
        call the specified function with the resulting file.
        """
        class FileDialogButton(Tkinter.Button):
            # This is just an ordinary button with special colors.

            def __init__(self, master=None, cnf={}, **kw):
                '''when we get a file, we call fn(filename)'''
                self.fn = fn
                self.__toggle = 0
                apply(Tkinter.Button.__init__, (self, master, cnf), kw)
                self.configure(command=self.set)
            def set(self):
                fd = PmwFileDialog(self.master,filter=filter)
                fd.title('Please choose a file')
                n=fd.askfilename()
                if n is not None:
                    self.fn(n)
        return FileDialogButton
    get = staticmethod(get)

class MyThread ( threading.Thread ):
  def run ( self ):
    os.system("start http://loschmidt.chemi.muni.cz/caver")

class AnBeKoM:
    def setBinaryLocation(self,value):
        self.binlocation.setvalue(value)
    def setPymolLocation(self,value):
        self.pymollocation.setvalue(value)
    def __init__(self,app):
        parent = app.root
        self.parent = parent

        #by default select all
        self.whichModelSelect = 'all';
        # backup the site
        self.originalX = -5
        self.originalY = -5
        self.originalZ = -5

        self.optimizeNearValue = StringVar()
        self.optimizeNearValue.set("4.0")

        # Create the dialog.
        self.dialog = Pmw.Dialog(parent,
                                 buttons = (defaults["compute_command"], defaults["exit_command"]),
                                 #defaultbutton = 'Run CAVER',
                                 title = 'Caver 2.1' + VERSION,
                                 command = self.execute)
        self.dialog.withdraw()
        #Pmw.setbusycursorattributes(self.dialog.component('hull'))
        lbb = "Caver 2.1%s" % (VERSION,)

        #labelfont = ('-weight bold')
        w = Tkinter.Label(self.dialog.interior(),
                                text = lbb ,
                                background = 'orange',
                                foreground = 'white',
                                #padx = 100,
                                )
        #w.config(font=labelfont)
        w.pack(expand = 1, fill = 'both', padx = 4, pady = 4)
        #ww = Tkinter.Label(self.dialog.interior(), background = 'orange', foreground="white", text = 'WWW and Help: http://loschmidt.chemi.muni.cz/caver')
        ww = Tkinter.Button(self.dialog.interior(), text = 'WWW and Help', command = self.launchHelp)
        ww.pack()

        self.stdam_list = [ 'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'ASX', 'CYX', 'GLX', 'HI0', 'HID', 'HIE', 'HIM', 'HIP', 'MSE', 'ACE', 'ASH', 'CYM', 'GLH', 'LYN', 'NME']

        #group = Pmw.Group(self.dialog.interior(),tag_text='Main options')
	#group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)

        def quickFileValidation(s):
            if s == '': return Pmw.PARTIAL
            elif os.path.isfile(s): return Pmw.OK
            elif os.path.exists(s): return Pmw.PARTIAL
            else: return Pmw.PARTIAL

        if (1):
          self.pymollocation = Pmw.EntryField(self.dialog.interior(),
                                       labelpos='w',
                                       value = PYMOL_LOCATION,
                                       label_text = LABEL_TEXT)

          self.pymollocation.pack(fill='x',padx=4,pady=1) # vertical
#win/linux

        self.binlocation = Pmw.EntryField(self.dialog.interior(),
                                     labelpos='w',
                                     value = OUTPUT_LOCATION,
                                     label_text = 'Output directory:')

        self.binlocation.pack(fill='x',padx=4,pady=1) # vertical

        self.tunnels = Pmw.EntryField(self.dialog.interior(),
                                     labelpos='w',
                                     value = defaults["default_tunnels"],
                                     label_text = 'Number of tunnels:')

        self.tunnels.pack(fill='x',padx=4,pady=1) # vertical

        #self.block = Pmw.EntryField(group.interior(),
        #                             labelpos='w',
        #                             value = defaults["default_block"],
        #                             label_text = 'Blocking radius:')


        #self.block.pack(fill='x',padx=4,pady=1) # vertical

	#labframe0 = Tkinter.Frame(group.interior())
	#labframe0.pack(fill='x',padx=4,pady=2)

        labframe0 = Tkinter.Frame(self.dialog.interior())
        labframe0.pack(fill='x',padx=4,pady=2)

        self.varremovewater = IntVar()

        self.removewaterbutton = Checkbutton(labframe0, text="Ignore waters", variable=self.varremovewater)
        self.varremovewater.set(1)
        #self.removewaterbutton.pack(side='left',expand='yes',)

        self.inModelGroup = Pmw.Group(self.dialog.interior(), tag_text='Input model:');
        self.listbox1 = Tkinter.Listbox(self.inModelGroup.interior(), width=25, height=6,exportselection=0)
        self.listbox1.bind('<<ListboxSelect>>',self.inputAnalyseWrap)
        yscroll1 = Tkinter.Scrollbar(self.inModelGroup.interior(),command=self.listbox1.yview, orient=Tkinter.VERTICAL)
        self.listbox1.pack(side=LEFT)
        yscroll1.pack(side=LEFT, fill='y')
        self.listbox1.configure(yscrollcommand=yscroll1.set)
        self.reloadListButton = Tkinter.Button(self.inModelGroup.interior(), text = 'Reload', command = self.updateList)
        self.reloadListButton.pack(side=LEFT)
        self.inModelGroup.pack()

        self.filterGroup = Pmw.Group(self.dialog.interior(), tag_text='Input atoms:');
        self.filterGroup.pack()
        self.checklist = []
        self.buttonlist = []

        self.updateList()
        #fill with data
        #self.listbox1.insert(0,"all")
        #self.listbox1.selection_set(0, 0) # Default sel
        #tindex = 1;
        #for item in cmd.get_object_list():
        #  self.listbox1.insert(tindex,str(item))
        #  tindex = tindex + 1;

        self.s = dict()
        self.s["AA"] = IntVar();
        self.reinitialise()
        self.inputAnalyse()

	groupstart = Pmw.Group(self.dialog.interior(),tag_text='Starting point')

#	labframe = Tkinter.Frame(groupstart.interior())
#	labframe.pack(fill='x',padx=4,pady=2)
#	label1 = Tkinter.Label(labframe,
#	                      justify=LEFT,
#                              text = "Specify x,y,z of starting point directly or use the conversion\nto get x,y,z from average centers of selection",
#			      )
#	label1.pack(side='left')

#-------------
#        radiogroups = []
        self.surroundingsvar = Tkinter.IntVar()
#        self.surroundingsvar.set(1)

        radioframe = Tkinter.Frame(groupstart.interior())
        group1 = Pmw.Group(radioframe,
                tag_text='Convert surroundings to x,y,x coordinates of starting point')

        group1.pack(side='top',expand = 'yes',fill='x')
#        cw = Tkinter.Frame(group1.interior())
#        cw.pack(padx = 2, pady = 2, expand='yes', fill='both')
#        radiogroups.append(group1)
	try: defaults['surroundings'] = cmd.get_names('objects')[0]
	except: defaults['surroundings'] = ''
	self.selectionlist = Pmw.EntryField(group1.interior(),
                                  labelpos='w',
                                  label_text='Specify selection: ',
                                  value=defaults['surroundings'],
                                  entry_width=50
                                  )
	self.selectionlist.pack(fill='x',expand='yes',padx=4,pady=1) # vertical

	self.convertButton = Tkinter.Button(group1.interior(), text = 'Convert to x,y,z', command = self.convert)
	self.convertButton.pack(fill='x',expand='yes',padx=4,pady=1)

        group2 = Pmw.Group(radioframe,
                tag_text='x, y, z coordinates of starting point')
        group2.pack(fill = 'x', expand = 1, side='top')
        radioframe.pack(side='left',expand='yes',fill='x')
#-------------
	groupstart.pack(side='left',padx=4,pady=1,expand='yes',fill='x')

	self.xlocvar=DoubleVar()
        self.xlocvar.set(float(defaults["startingpoint"][0]))
	self.ylocvar=DoubleVar()
        self.ylocvar.set(float(defaults["startingpoint"][1]))
	self.zlocvar=DoubleVar()
        self.zlocvar.set(float(defaults["startingpoint"][2]))

	self.xlocfr = Tkinter.Frame(group2.interior())
	labX = Label(self.xlocfr,text="x");
	self.xlocation = Entry(self.xlocfr,textvariable=self.xlocvar);
	self.scrX=Scrollbar(self.xlocfr,orient="horizontal",command=self.changeValueX)

        self.ylocfr = Tkinter.Frame(group2.interior())
	labY = Label(self.ylocfr,text="y");
	self.ylocation = Entry(self.ylocfr,textvariable=self.ylocvar);
	self.scrY=Scrollbar(self.ylocfr,orient="horizontal",command=self.changeValueY)

        self.zlocfr = Tkinter.Frame(group2.interior())
	labZ = Label(self.zlocfr,text="z");
	self.zlocation = Entry(self.zlocfr,textvariable=self.zlocvar);
	self.scrZ=Scrollbar(self.zlocfr,orient="horizontal",command=self.changeValueZ)

	labX.pack(side=LEFT)
	self.xlocation.pack(side=LEFT)
	self.scrX.pack(side=LEFT)
	self.xlocfr.pack(fill='x',padx=4,pady=1) # vertical
	labY.pack(side=LEFT)
	self.ylocation.pack(side=LEFT)
	self.scrY.pack(side=LEFT)
	self.ylocfr.pack(fill='x',padx=4,pady=1) # vertical
	labZ.pack(side=LEFT)
	self.zlocation.pack(side=LEFT)
	self.scrZ.pack(side=LEFT)
	self.zlocfr.pack(fill='x',padx=4,pady=1) # vertical

	self.OpGroup = Pmw.Group(radioframe,tag_text = "Optimize starting point")
	self.OpGroup.pack(fill='x');
	self.optimizeLabel = Tkinter.Label(self.OpGroup.interior(),text = 'Neighbourhood:')
	self.optimizeLabel.pack(side=LEFT)

	self.optimizeNear = Tkinter.Entry(self.OpGroup.interior(),textvariable=self.optimizeNearValue,justify='right')
	self.optimizeNear.pack(side=LEFT,padx=4,pady=1)
	self.angstrom = Tkinter.Label(self.OpGroup.interior(),text="(A)")
	self.angstrom.pack(side=LEFT, padx=0, pady=1)
	self.optimizeButton = Tkinter.Button(self.OpGroup.interior(), text = 'Optimize', command = self.optimize)
	self.optimizeButton.pack(side=LEFT,padx=5,pady=1)
	self.UoptimizeButton = Tkinter.Button(self.OpGroup.interior(), text = 'Undo', command = self.uoptimize)
	self.UoptimizeButton.pack(side=LEFT,padx=1,pady=1)

# blocking method =========================
#	groupMethod = Pmw.Group(self.dialog.interior(),tag_text='Choose blocking method')

# by default, set 4: branch blocking
        self.methodvar = Tkinter.IntVar()
        self.methodvar.set(4)

#        methodframe = Tkinter.Frame(groupMethod.interior())
#        methodframe.pack(expand='yes')
#        metody=[('Simple normal', 0), ('Simple restrictive', 1), ('Decreasing normal', 2),('Decreasing restrictive', 3), ('Branching normal (default)', 4), ('Branching restrictive', 5), ('Narrowest', 6)]
#
#        for text, value in metody:
#          Radiobutton(methodframe, text=text, value=value, variable=self.methodvar).pack(anchor=W)


#        groupMethod.pack(fill = 'both', expand = 1, padx = 10, pady = 5)

# eoblock method ==========================
        self.showAppModal()

    def showCrisscross(self):
	    startpoint=(float(self.xlocvar.get()),float(self.ylocvar.get()),float(self.zlocvar.get()))
            cmd.delete("crisscross")
	    self.crisscross(startpoint[0],startpoint[1],startpoint[2],0.5,"crisscross")

#win/linux
    if WINDOWZ:
      def changeValueX(self,obj,a,b=0,c=0):
        val=float(self.xlocvar.get())+float(a)*0.2
        self.xlocvar.set(val)
        self.showCrisscross()
      def changeValueY(self,obj,a,b=0,c=0):
        val=float(self.ylocvar.get())+float(a)*0.2
        self.ylocvar.set(val)
        self.showCrisscross()
      def changeValueZ(self,obj,a,b=0,c=0):
        val=float(self.zlocvar.get())+float(a)*0.2
        self.zlocvar.set(val)
        self.showCrisscross()
    else:
      def changeValueX(self,a):
        val=float(self.xlocvar.get())+float(a)*0.2
        self.xlocvar.set(val)
        self.showCrisscross()
      def changeValueY(self,a):
        val=float(self.ylocvar.get())+float(a)*0.2
        self.ylocvar.set(val)
        self.showCrisscross()
      def changeValueZ(self,a):
        val=float(self.zlocvar.get())+float(a)*0.2
        self.zlocvar.set(val)
        self.showCrisscross()

    def showAppModal(self):
        self.dialog.show()

    def updateList(self):
      self.listbox1.delete(0, Tkinter.END)
      #fill with data
      self.listbox1.insert(0,"all")
      self.listbox1.selection_set(0, 0) # Default sel
      tindex = 1;
      for item in cmd.get_object_list():
        self.listbox1.insert(tindex,str(item))
        tindex = tindex + 1;
      self.inputAnalyse()

    def launchHelp(self):
      if WINDOWZ:
        t = MyThread();
        t.start()
      else:
        import webbrowser
        webbrowser.open("http://loschmidt.chemi.muni.cz/caver")
        Pmw.MessageDialog(self.parent,title = 'Information',message_text = "see http://loschmidt.chemi.muni.cz/caver")

    def testBinary(self):
            #test
            if WINDOWZ:
              #anbekomloc = "%s\\modules\\Caver2_1_%s\\Caver2_1.jar" % (self.pymollocation.getvalue(), VERS,)
              anbekomloc = "%s\\Caver2_1.jar" % (self.pymollocation.getvalue(), )
              if not os.path.isfile(anbekomloc):
                          error_dialog = Pmw.MessageDialog(self.parent,title = 'Error',message_text = "ERROR: PyMOL directory incorrectly specified or Caver2_1 not installed, check %s" % anbekomloc ,)
                          junk = error_dialog.activate()
                          return 0
            #test
            if not WINDOWZ:
              anbekomloc = "%s/Caver2_1.jar" % (self.pymollocation.getvalue(), )
              if not os.path.isfile(anbekomloc):
                          error_dialog = Pmw.MessageDialog(self.parent,title = 'Error',message_text = "ERROR: Directory incorrectly specified -plugin .jar not found, check %s" % anbekomloc ,)
                          junk = error_dialog.activate()
                          return 0
            return 1
    def execute(self, result):

	if result == defaults["compute_command"]:
            if self.testBinary() == 0:
              return
            self.showCrisscross()
            #input
            sel1index = self.listbox1.curselection()[0]
            sel1text = self.listbox1.get(sel1index)
            self.whichModelSelect = sel1text;
            print 'selected ' + self.whichModelSelect
            sel=cmd.get_model(self.whichModelSelect)
            cnt=0
            for a in sel.atom:
               cnt+=1
            print cnt

            if cnt == 0:
                error_dialog = Pmw.MessageDialog(self.parent,title = 'Error',
                                   message_text = 'ERROR: No molecule loaded.',)
                junk = error_dialog.activate()
                return
            outdir = self.binlocation.getvalue()
            if os.path.isfile(outdir):
               error_dialog = Pmw.MessageDialog(self.parent,title = 'Error',
                                   message_text = 'ERROR: Output directory is file.',)
               junk = error_dialog.activate()
               return
            elif not os.path.exists(outdir):
               self.CreateDirectory(outdir)
            self.stdamString = string.join(self.stdam_list, "+")
            # jen to zaskrtnute
            generatedString = ""
            for key in self.s:
              if self.s[key].get() == 1:
                # pak pouzit do vyberu:
                if key == "AA":
                  generatedString = generatedString + "+" + self.stdamString
                else:
                  generatedString = generatedString + "+" + key

            generatedString = generatedString[1:]
            print "Checked: " + generatedString

            mmodel = cmd.get_model(self.whichModelSelect)
            print self.whichModelSelect + " asize: " + str(len(mmodel.atom))
            newmodel = Indexed()
            for matom in mmodel.atom:
              if generatedString.find(matom.resn) > -1:
                #print matom.resn
                newmodel.atom.append(matom)
            cmd.load_model(newmodel,"tmpCaverModel")
            #cmd.label("example","name")
            #fix outdir slashes
            outdir = outdir.replace("\\","/")
            if (outdir.endswith("/")):
				outdir = outdir[:-1]
            input = "%s/out.pdb" % (outdir)
            #cmd.save(input, self.whichModelSelect) # to by ulozilo cely model whichModelSelect.
            cmd.save(input, "tmpCaverModel")
            cmd.delete("tmpCaverModel")
            cesta = os.getcwd()
            # set ignore waters to false -- the model is already filtered by input model and aminos
            self.varremovewater.set(0)
            if WINDOWZ:
              #commandXYZ = "java %s -jar \"%s/modules/Caver2_1_%s/Caver2_1.jar\" \"%s\" %f %f %f %s \"%s\" %d %d %d %s" % (JOPTS, self.pymollocation.getvalue(), VERS, input, float(self.xlocvar.get()), float(self.ylocvar.get()), float(self.zlocvar.get()), self.tunnels.getvalue(), outdir, self.varremovewater.get(), 0,self.methodvar.get(), "tun_" + self.whichModelSelect)
              commandXYZ = "java %s -jar \"%s\Caver2_1.jar\" \"%s\" %f %f %f %s \"%s\" %d %d %d %s" % (JOPTS, self.pymollocation.getvalue(),input, float(self.xlocvar.get()), float(self.ylocvar.get()), float(self.zlocvar.get()), self.tunnels.getvalue(), outdir, self.varremovewater.get(),0,self.methodvar.get(), "tun_" + self.whichModelSelect)
            else:
              commandXYZ = "java %s -jar \"%s/Caver2_1.jar\" \"%s\" %f %f %f %s \"%s\" %d %d %d %s" % (JOPTS, self.pymollocation.getvalue(),input, float(self.xlocvar.get()), float(self.ylocvar.get()), float(self.zlocvar.get()), self.tunnels.getvalue(), outdir, self.varremovewater.get(),0,self.methodvar.get(), "tun_" + self.whichModelSelect)
            tunnels = int(self.tunnels.getvalue())
            # ted vymazat v output dir vsechny soubory s path_*.py os.path.isfile(string)
            for i in range(tunnels):
              tunpy = "%s/path_%i.py" % (outdir,i)
              if (os.path.isfile(tunpy)):
                os.remove(tunpy);
            print commandXYZ
            #os.system(commandXYZ)
            status = subprocess.call(commandXYZ, shell=True)
            for i in range(tunnels):
               pathpy = "%s/path_%i.py" % (outdir, i)
               if os.access(pathpy,os.F_OK):
                  view = cmd.get_view()
                  execfile(pathpy)
                  cmd.set_view(view)
                  if i != 0:
                     cmd.disable("tunnel%i" % i);
               else:
                  if i == 0:
                    error_dialog = Pmw.MessageDialog(self.parent,title = 'Error',
                                       message_text = 'Error: No tunnel was found, the starting point is probably outside of the molecule.',)

                    junk = error_dialog.activate()
                    break
                  else:
                    messa = "No more tunnels than %i were found" % (i)
                    error_dialog = Pmw.MessageDialog(self.parent,title = 'Info',
                                       message_text = messa,)
                    junk = error_dialog.activate()
                    break
	    #pass
	    #self.deleteTemporaryFiles()
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
                #self.dialog.deactivate(result)
                global CAVER_BINARY_LOCATION
                CAVER_BINARY_LOCATION = self.binlocation.getvalue()
                self.dialog.withdraw()

    def CreateDirectory(self,dir):
      if os.path.isdir(dir):
       return;
      parent, base = os.path.split(dir)
      self.CreateDirectory(parent)
      os.mkdir(dir)

    def convert(self):
      sel=cmd.get_model('(all)')
      cnt=0
      for a in sel.atom:
          cnt+=1
      if cnt == 0:
          error_dialog = Pmw.MessageDialog(self.parent,title = 'Error',message_text = 'ERROR: No molecule loaded.',)
      try:
		startpoint=[]
		startpoint=self.computecenter(self.selectionlist.getvalue())
		self.xlocvar.set(startpoint[0])
		self.ylocvar.set(startpoint[1])
		self.zlocvar.set(startpoint[2])
		self.crisscross(startpoint[0],startpoint[1],startpoint[2],0.5,"crisscross")
      except:
		error_dialog = Pmw.MessageDialog(self.parent,title = 'Error',message_text = 'ERROR: Invalid selection name',)

    def uoptimize(self):
      self.xlocvar.set(self.originalX)
      self.ylocvar.set(self.originalY)
      self.zlocvar.set(self.originalZ)
      self.showCrisscross()

    def optimize(self):
      if self.testBinary() == 0:
        return
      self.originalX = self.xlocvar.get()
      self.originalY = self.ylocvar.get()
      self.originalZ = self.zlocvar.get()
      outdir = self.binlocation.getvalue()
      if os.path.isfile(outdir):
         error_dialog = Pmw.MessageDialog(self.parent,title = 'Error',
                             message_text = 'ERROR: Output directory is file.',)
         junk = error_dialog.activate()
         return
      elif not os.path.exists(outdir):
         self.CreateDirectory(outdir)

      #fix outdir slashes
      outdir = outdir.replace("\\","/")
      if (outdir.endswith("/")):
			outdir = outdir[:-1]
      input = "%s/outOpti.pdb" % (outdir)
      #input
      sel1index = self.listbox1.curselection()[0]
      sel1text = self.listbox1.get(sel1index)
      self.whichModelSelect = sel1text;

      #save selected
      cmd.save(input, self.whichModelSelect)
      commandOPT =""
      if WINDOWZ:
        commandOPT = "java -classpath \"%s/modules/Caver2_1_%s/Caver2_1.jar\" AnBeKoM.ActiveSiteOptimize \"%s\" %f %f %f \"%s\" %f" % (self.pymollocation.getvalue(), VERS, input, float(self.xlocvar.get()), float(self.ylocvar.get()), float(self.zlocvar.get()),  outdir,float(self.optimizeNearValue.get()))
      else:
        commandOPT = "java -classpath \"%s/Caver2_1.jar\" AnBeKoM.ActiveSiteOptimize \"%s\" %f %f %f \"%s\" %f" % (self.pymollocation.getvalue(),input, float(self.xlocvar.get()), float(self.ylocvar.get()), float(self.zlocvar.get()), outdir, float(self.optimizeNearValue.get()))

      outAsite = "%s/startingPointOptimized.txt" % outdir
      if os.path.exists(outAsite):
        os.remove(outAsite);

      print commandOPT
      #os.system(commandOPT)
      status = subprocess.call(commandOPT, shell=True)
      if not os.path.exists(outAsite):
         error_dialog = Pmw.MessageDialog(self.parent,title = 'Error',
                             message_text = 'ERROR: Optimize failed. Check starting point above',)
         junk = error_dialog.activate()
         return
      try:
        optAsite = "";
        file = open(outAsite)
        while 1:
            line = file.readline()
            if not line:
                break
            if line.startswith('%'):
                break
            optAsite = line.rstrip('\n').split(' ');
        self.xlocvar.set(optAsite[0]);
        self.ylocvar.set(optAsite[1]);
        self.zlocvar.set(optAsite[2]);
        self.crisscross(self.xlocvar.get(),self.ylocvar.get(),self.zlocvar.get(),0.5,"crisscross")
        self.showCrisscross()
      except Exception:
        Pmw.MessageDialog(self.parent,title = 'Information',message_text = 'Optimize failed, probably outside of: ' + self.whichModelSelect)

    def containsValue(self, array, value):
	    for v in array:
	      if (v == value):
	        return 1
	    return 0

    def stdamMessage(self):
	    Pmw.MessageDialog(self.parent,title = 'Information',message_text = 'AA: Standard amino acids: \n ' + string.join(self.stdam_list, ", "))

    def inputAnalyseWrap(self, args):
	      #print self.listbox1.curselection()[0] # aby to fungovalo, musi byt bindnute na <<ListboxSelect>>
	      self.inputAnalyse()

    def inputAnalyse(self):
	      sel1index = self.listbox1.curselection()[0]
	      sel1text = self.listbox1.get(sel1index)
	      self.whichModelSelect = sel1text;
	      sel=cmd.get_model(self.whichModelSelect)
	      #pripravit kontrolni strukturu pro nalezene
	      self.s = dict()
	      self.s.clear()
	      self.s["AA"] = IntVar()
	      #cntr = 0
	      for a in sel.atom:
	        if not (self.s.has_key(a.resn)):
	          if (self.containsValue(self.stdam_list,a.resn)):
	             self.s["AA"].set(1)
	          else:
	            self.s[a.resn] = IntVar()
	            self.s[a.resn].set(0)  # default only standard checked

	      self.reinitialise()

    def reinitialise(self):
	        #TODO: pouzivat zde uz setrideny
	        ksorted = sorted(self.s.keys())
	        ksorted.remove("AA")
	        ksorted = sorted(ksorted)
	        stdstrings = ["AA"]
	        #na zacatek dam standardni a pak setrideny zbytek
	        ksorted = stdstrings + ksorted
	        #print "calling initialise"
	        for xs in self.checklist:
	          xs.grid_remove()
	        self.checklist = []

	        for xs in self.buttonlist:
	          xs.grid_remove()
	        self.buttonlist = []

	        cntr = 0
	        # tady uz setrideny, se STDAM a STDRNA na zacatku
	        for key in ksorted:
	          #print "{" + key + "}{" + str(self.s[key].get()) + "}"
	          #
	          if cntr == 1:
	            cntr = cntr + 4
	          tmpButton = Tkinter.Checkbutton(self.filterGroup.interior(), text=key, variable=self.s[key])
	          tmpButton.var = self.s[key]
	          tmpButton.grid(sticky=W, row = int(cntr/5), column = (cntr % 5))
	          self.checklist.append(tmpButton)

	          # zaridit balloon help -- neni potreba kdyz je tlacitko
	          #if key == "AA":
	          #  balloon = Pmw.Balloon(self.parent)
	          #  balloon.bind(tmpButton, 'STanDard AMino acids: \n ' + string.join(self.stdam_list, ", "), 'STanDard AMino acids')

	          # kdyz poprve, je tam pridano STDAM, vlozit tam  tedy i napovedu
	          if cntr == 0:
	            xButton = Tkinter.Button(self.filterGroup.interior(), text='?', command=self.stdamMessage, width = 5)
	            xButton.grid(sticky=W, row = 0, column=1) # 0,1 = stdam, 0,2 = help

	          cntr = cntr + 1
	        ##print "size: " + str(len(self.checklist))
	        #aButton = Tkinter.Button(self.filterGroup.interior(), text='Reload', command=self.inputAnalyse)
	        ## zarovnat
	        #if not cntr % 3 == 0:
	        #  cntr = cntr + 3 - (cntr % 3)
	        ##analyse,save,load
	        #aButton.grid(row = int(cntr/3), column = (cntr % 3))
	        #self.buttonlist.append(aButton)
    def computecenter(self,selection="(all)"):
	gcentx=0
	gcenty=0
	gcentz=0
	gcnt=0
	for selstr in selection.split():
		sel=cmd.get_model(selstr)
		centx=0
		centy=0
		centz=0
		cnt=len(sel.atom)
		for a in sel.atom:
			centx+=a.coord[0]
			centy+=a.coord[1]
			centz+=a.coord[2]
		centx/=cnt
		centy/=cnt
		centz/=cnt
	#	fmttext="%lf\t%lf\t%lf\n" % (centx,centy,centz)
#		print centx,centy,centz
		gcentx+=centx
		gcenty+=centy
		gcentz+=centz
		gcnt+=1
	gcentx/=gcnt
	gcenty/=gcnt
	gcentz/=gcnt
	return (gcentx,gcenty,gcentz)

    def TestProgram(self,cmd):

#	f1 = os.popen(cmd, 'r');
#	content=f1.read()
#	status = f1.close()
#	if status is None:
#	   return (content,0)
#	else:
#	   return ("",1)
    #status = os.system(cmd)
        status = subprocess.call(cmd, shell=True)
	return status

    def box(self,x1,y1,z1,x2,y2,z2,name="box"):

		obj = [
		BEGIN, LINE_STRIP,
		VERTEX, float(x1), float(y1), float(z1),
		VERTEX, float(x2), float(y1), float(z1),
		VERTEX, float(x2), float(y2), float(z1),
		VERTEX, float(x1), float(y2), float(z1),
		VERTEX, float(x1), float(y1), float(z1),
		END,

		BEGIN, LINE_STRIP,
		VERTEX, float(x1), float(y1), float(z2),
		VERTEX, float(x2), float(y1), float(z2),
		VERTEX, float(x2), float(y2), float(z2),
		VERTEX, float(x1), float(y2), float(z2),
		VERTEX, float(x1), float(y1), float(z2),
		END,

		BEGIN, LINES,
		VERTEX, float(x1), float(y1), float(z1),
		VERTEX, float(x1), float(y1), float(z2),
		END,

		BEGIN, LINES,
		VERTEX, float(x2), float(y1), float(z1),
		VERTEX, float(x2), float(y1), float(z2),
		END,

		BEGIN, LINES,
		VERTEX, float(x1), float(y2), float(z1),
		VERTEX, float(x1), float(y2), float(z2),
		END,

		BEGIN, LINES,
		VERTEX, float(x2), float(y2), float(z1),
		VERTEX, float(x2), float(y2), float(z2),
		END

		]

		cmd.load_cgo(obj,name)


    def crisscross(self,x,y,z,d,name="crisscross"):
		obj = [
		LINEWIDTH, 3,

		BEGIN, LINE_STRIP,
		VERTEX, float(x-d), float(y), float(z),
		VERTEX, float(x+d), float(y), float(z),
		END,

		BEGIN, LINE_STRIP,
		VERTEX, float(x), float(y-d), float(z),
		VERTEX, float(x), float(y+d), float(z),
		END,

		BEGIN, LINE_STRIP,
		VERTEX, float(x), float(y), float(z-d),
		VERTEX, float(x), float(y), float(z+d),
		END

		]
		view = cmd.get_view()
		cmd.load_cgo(obj,name)
		cmd.set_view(view)

#
# The classes PmwFileDialog and PmwExistingFileDialog and the _errorpop function
# are taken from the Pmw contrib directory.  The attribution given in that file
# is:
################################################################################
# Filename dialogs using Pmw
#
# (C) Rob W.W. Hooft, Nonius BV, 1998
#
# Modifications:
#
# J. Willem M. Nissink, Cambridge Crystallographic Data Centre, 8/2002
#    Added optional information pane at top of dialog; if option
#    'info' is specified, the text given will be shown (in blue).
#    Modified example to show both file and directory-type dialog
#
# No Guarantees. Distribute Freely.
# Please send bug-fixes/patches/features to <r.hooft@euromail.com>
#
################################################################################
import os,fnmatch,time
import Tkinter,Pmw
#Pmw.setversion("0.8.5")

def _errorpop(master,text):
    d=Pmw.MessageDialog(master,
                        title="Error",
                        message_text=text,
                        buttons=("OK",))
    d.component('message').pack(ipadx=15,ipady=15)
    d.activate()
    d.destroy()

class PmwFileDialog(Pmw.Dialog):
    """File Dialog using Pmw"""
    def __init__(self, parent = None, **kw):
	# Define the megawidget options.
	optiondefs = (
	    ('filter',    '*',              self.newfilter),
	    ('directory', os.getcwd(),      self.newdir),
	    ('filename',  '',               self.newfilename),
	    ('historylen',10,               None),
	    ('command',   None,             None),
            ('info',      None,             None),
	    )
	self.defineoptions(kw, optiondefs)
        # Initialise base class (after defining options).
	Pmw.Dialog.__init__(self, parent)

	self.withdraw()

        # Create the components.
	interior = self.interior()

        if self['info'] is not None:
            rowoffset=1
            dn = self.infotxt()
            dn.grid(row=0,column=0,columnspan=2,padx=3,pady=3)
        else:
            rowoffset=0

	dn = self.mkdn()
	dn.grid(row=0+rowoffset,column=0,columnspan=2,padx=3,pady=3)
	del dn

	# Create the directory list component.
	dnb = self.mkdnb()
	dnb.grid(row=1+rowoffset,column=0,sticky='news',padx=3,pady=3)
	del dnb

	# Create the filename list component.
	fnb = self.mkfnb()
	fnb.grid(row=1+rowoffset,column=1,sticky='news',padx=3,pady=3)
	del fnb

	# Create the filter entry
	ft = self.mkft()
	ft.grid(row=2+rowoffset,column=0,columnspan=2,padx=3,pady=3)
	del ft

	# Create the filename entry
	fn = self.mkfn()
	fn.grid(row=3+rowoffset,column=0,columnspan=2,padx=3,pady=3)
	fn.bind('<Return>',self.okbutton)
	del fn

	# Buttonbox already exists
	bb=self.component('buttonbox')
	bb.add('OK',command=self.okbutton)
	bb.add('Cancel',command=self.cancelbutton)
	del bb

	Pmw.alignlabels([self.component('filename'),
			 self.component('filter'),
			 self.component('dirname')])

    def infotxt(self):
        """ Make information block component at the top """
        return self.createcomponent(
                'infobox',
                (), None,
                Tkinter.Label, (self.interior(),),
                width=51,
                relief='groove',
                foreground='darkblue',
                justify='left',
                text=self['info']
            )

    def mkdn(self):
        """Make directory name component"""
        return self.createcomponent(
	    'dirname',
	    (), None,
	    Pmw.ComboBox, (self.interior(),),
	    entryfield_value=self['directory'],
	    entryfield_entry_width=40,
            entryfield_validate=self.dirvalidate,
	    selectioncommand=self.setdir,
	    labelpos='w',
	    label_text='Directory:')

    def mkdnb(self):
        """Make directory name box"""
        return self.createcomponent(
	    'dirnamebox',
	    (), None,
	    Pmw.ScrolledListBox, (self.interior(),),
	    label_text='directories',
	    labelpos='n',
	    hscrollmode='none',
	    dblclickcommand=self.selectdir)

    def mkft(self):
        """Make filter"""
        return self.createcomponent(
	    'filter',
	    (), None,
	    Pmw.ComboBox, (self.interior(),),
	    entryfield_value=self['filter'],
	    entryfield_entry_width=40,
	    selectioncommand=self.setfilter,
	    labelpos='w',
	    label_text='Filter:')

    def mkfnb(self):
        """Make filename list box"""
        return self.createcomponent(
	    'filenamebox',
	    (), None,
	    Pmw.ScrolledListBox, (self.interior(),),
	    label_text='files',
	    labelpos='n',
	    hscrollmode='none',
	    selectioncommand=self.singleselectfile,
	    dblclickcommand=self.selectfile)

    def mkfn(self):
        """Make file name entry"""
        return self.createcomponent(
	    'filename',
	    (), None,
	    Pmw.ComboBox, (self.interior(),),
	    entryfield_value=self['filename'],
	    entryfield_entry_width=40,
            entryfield_validate=self.filevalidate,
	    selectioncommand=self.setfilename,
	    labelpos='w',
	    label_text='Filename:')

    def dirvalidate(self,string):
        if os.path.isdir(string):
            return Pmw.OK
        else:
            return Pmw.PARTIAL

    def filevalidate(self,string):
        if string=='':
            return Pmw.PARTIAL
        elif os.path.isfile(string):
            return Pmw.OK
        elif os.path.exists(string):
            return Pmw.PARTIAL
        else:
            return Pmw.OK

    def okbutton(self):
	"""OK action: user thinks he has input valid data and wants to
           proceed. This is also called by <Return> in the filename entry"""
	fn=self.component('filename').get()
	self.setfilename(fn)
	if self.validate(fn):
	    self.canceled=0
	    self.deactivate()

    def cancelbutton(self):
	"""Cancel the operation"""
	self.canceled=1
	self.deactivate()

    def tidy(self,w,v):
	"""Insert text v into the entry and at the top of the list of
           the combobox w, remove duplicates"""
	if not v:
	    return
	entry=w.component('entry')
	entry.delete(0,'end')
	entry.insert(0,v)
	list=w.component('scrolledlist')
	list.insert(0,v)
	index=1
	while index<list.index('end'):
	    k=list.get(index)
	    if k==v or index>self['historylen']:
		list.delete(index)
	    else:
		index=index+1
        w.checkentry()

    def setfilename(self,value):
	if not value:
	    return
	value=os.path.join(self['directory'],value)
	dir,fil=os.path.split(value)
	self.configure(directory=dir,filename=value)

	c=self['command']
	if callable(c):
	    c()

    def newfilename(self):
	"""Make sure a newly set filename makes it into the combobox list"""
	self.tidy(self.component('filename'),self['filename'])

    def setfilter(self,value):
	self.configure(filter=value)

    def newfilter(self):
	"""Make sure a newly set filter makes it into the combobox list"""
	self.tidy(self.component('filter'),self['filter'])
	self.fillit()

    def setdir(self,value):
	self.configure(directory=value)

    def newdir(self):
	"""Make sure a newly set dirname makes it into the combobox list"""
	self.tidy(self.component('dirname'),self['directory'])
	self.fillit()

    def singleselectfile(self):
	"""Single click in file listbox. Move file to "filename" combobox"""
	cs=self.component('filenamebox').curselection()
	if cs!=():
	    value=self.component('filenamebox').get(cs)
            self.setfilename(value)

    def selectfile(self):
	"""Take the selected file from the filename, normalize it, and OK"""
        self.singleselectfile()
	value=self.component('filename').get()
        self.setfilename(value)
        if value:
	    self.okbutton()

    def selectdir(self):
	"""Take selected directory from the dirnamebox into the dirname"""
	cs=self.component('dirnamebox').curselection()
	if cs!=():
	    value=self.component('dirnamebox').get(cs)
	    dir=self['directory']
	    if not dir:
		dir=os.getcwd()
	    if value:
		if value=='..':
		    dir=os.path.split(dir)[0]
		else:
		    dir=os.path.join(dir,value)
	    self.configure(directory=dir)
	    self.fillit()

    def askfilename(self,directory=None,filter=None):
	"""The actual client function. Activates the dialog, and
	   returns only after a valid filename has been entered
           (return value is that filename) or when canceled (return
           value is None)"""
	if directory!=None:
	    self.configure(directory=directory)
	if filter!=None:
	    self.configure(filter=filter)
	self.fillit()
        self.canceled=1 # Needed for when user kills dialog window
	self.activate()
	if self.canceled:
	    return None
	else:
	    return self.component('filename').get()

    lastdir=""
    lastfilter=None
    lasttime=0
    def fillit(self):
	"""Get the directory list and show it in the two listboxes"""
        # Do not run unnecesarily
        if self.lastdir==self['directory'] and self.lastfilter==self['filter'] and self.lasttime>os.stat(self.lastdir)[8]:
            return
        self.lastdir=self['directory']
        self.lastfilter=self['filter']
        self.lasttime=time.time()
	dir=self['directory']
	if not dir:
	    dir=os.getcwd()
	dirs=['..']
	files=[]
        try:
            fl=os.listdir(dir)
            fl.sort()
        except os.error,arg:
            if arg[0] in (2,20):
                return
            raise
	for f in fl:
	    if os.path.isdir(os.path.join(dir,f)):
		dirs.append(f)
	    else:
		filter=self['filter']
		if not filter:
		    filter='*'
		if fnmatch.fnmatch(f,filter):
		    files.append(f)
	self.component('filenamebox').setlist(files)
	self.component('dirnamebox').setlist(dirs)

    def validate(self,filename):
	"""Validation function. Should return 1 if the filename is valid,
           0 if invalid. May pop up dialogs to tell user why. Especially
           suited to subclasses: i.e. only return 1 if the file does/doesn't
           exist"""
	return 1

class PmwExistingFileDialog(PmwFileDialog):
    def filevalidate(self,string):
        if os.path.isfile(string):
            return Pmw.OK
        else:
            return Pmw.PARTIAL

    def validate(self,filename):
        if os.path.isfile(filename):
            return 1
        elif os.path.exists(filename):
            _errorpop(self.interior(),"This is not a plain file")
            return 0
        else:
            _errorpop(self.interior(),"Please select an existing file")
            return 0

# Create demo in root window for testing.
if __name__ == '__main__':
    class App:
        def my_show(self,*args,**kwargs):
            pass
    app = App()
    app.root = Tkinter.Tk()
    Pmw.initialise(app.root)
    app.root.title('Some Title')

    widget = AnBeKoM(app)
    exitButton = Tkinter.Button(app.root, text = 'Exit', command = app.root.destroy)
    exitButton.pack()
    app.root.mainloop()
