#!/usr/bin/env python

# MOLE Ver. 1.2; 14. Nov. 2007
#
# MOLE Copyright Notice
# ============================
#

from __future__ import division
from __future__ import generators
from __future__ import print_function
from __future__ import absolute_import

import os
import math
import platform
import Pmw
import distutils.spawn  # used for find_executable
import pymol
from pymol import cmd, selector
import sys
import subprocess
from pymol.cmd import _feedback, fb_module, fb_mask, is_list, _cmd
from pymol.cgo import *
from chempy.models import Indexed
from chempy import Bond, Atom

if sys.version_info[0] < 3:
    import Tkinter
    from Tkinter import *
else:
    import tkinter as Tkinter
    from tkinter import *

#
# Global config variables
#
MOLE_OUTPUT = os.getcwd()

if 'PYMOL_GIT_MOD' in os.environ:
    pymol_env = os.environ.copy()
    #pymol_env['PYTHONPATH'] = os.path.join(os.environ['PYMOL_GIT_MOD'],"Mole")
    sys.path.append(os.path.join(os.environ['PYMOL_GIT_MOD'], "Mole"))
    if sys.platform.startswith('linux'):
        MOLE_BINARY_LOCATION = os.path.join(os.environ['PYMOL_GIT_MOD'], "Mole", "linux", "bin", "Mole.exe")
        qhull_dir = os.path.join(os.environ['PYMOL_GIT_MOD'], "Mole", "linux", "bin")
        mole_dir = os.path.join(os.environ['PYMOL_GIT_MOD'], "Mole", "linux")
        pymol_env['PATH'] = pymol_env['PATH'] + ":" + qhull_dir
        pymol_env['MOLEDIR'] = mole_dir
    elif sys.platform.startswith('darwin'):
        MOLE_BINARY_LOCATION = os.path.join(os.environ['PYMOL_GIT_MOD'], "Mole", "mac", "bin", "Mole.exe")
        qhull_dir = os.path.join(os.environ['PYMOL_GIT_MOD'], "Mole", "mac", "bin")
        mole_dir = os.path.join(os.environ['PYMOL_GIT_MOD'], "Mole", "mac")
        pymol_env['PATH'] = pymol_env['PATH'] + ":" + qhull_dir
        pymol_env['MOLEDIR'] = mole_dir
    elif sys.platform.startswith('win'):
        MOLE_BINARY_LOCATION = os.path.join(os.environ['PYMOL_GIT_MOD'], "Mole", "win32", "bin", "Mole.exe")
        qhull_dir = os.path.join(os.environ['PYMOL_GIT_MOD'], "Mole", "win32", "bin")
        mole_dir = os.path.join(os.environ['PYMOL_GIT_MOD'], "Mole", "win32")
        pymol_env['PATH'] = pymol_env['PATH'] + ":" + qhull_dir
        pymol_env['MOLEDIR'] = mole_dir
    else:
        MOLE_BINARY_LOCATION = None
        qhull_dir = ""
        mole_dir = ""
else:
    MOLE_BINARY_LOCATION = None
    qhull_dir = ""
    mole_dir = ""

if MOLE_OUTPUT is None:
    if 'MOLEDIR' in os.environ:
        MOLE_OUTPUT = os.path.join(os.environ['MOLEDIR'], "tmp")
    else:
        if 'TEMP' in os.environ:
            MOLE_OUTPUT = os.environ['TEMP']
        else:
            MOLE_OUTPUT = ''

if MOLE_BINARY_LOCATION is None:
    if 'MOLEDIR' in os.environ:
        MOLE_BINARY_LOCATION = os.path.join(os.environ['MOLEDIR'], "bin", "Mole.exe")
    else:
        MOLE_BINARY_LOCATION = distutils.spawn.find_executable('mole')
        if MOLE_BINARY_LOCATION is None:
            MOLE_BINARY_LOCATION = ''

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

            def load(self, name, sel=''):
                pass

            def get_names(self):
                return ['mol1', 'mol2', 'map1', 'map2']

            def get_type(self, thing):
                if thing.startswith('mol'):
                    return 'object:molecule'
                else:
                    return 'object:map'
                f.close()
        cmd = cmd()
    pymol = pymol()


def __init__(self):
    self.menuBar.addmenuitem('Plugin', 'command',
                             'Launch MOLE Tools',
                             label='MOLE Tools...',
                             command=lambda s=self: MOLETools(s))

defaults = {
    "residues": '(sele)',
    "startingpoint": ('0.0', '0.0', '0.0'),
    "outdir": MOLE_OUTPUT,
    "numtunnels": 3,
    "numinterp": 200,
    "activesiteradius": 3.0,
}

#   "atomselection": "all",


class FileDialogButtonClassFactory:

    def get(fn, filter='*'):
        """This returns a FileDialogButton class that will
        call the specified function with the resulting file.
        """
        class FileDialogButton(Tkinter.Button):
            # This is just an ordinary button with special colors.

            def __init__(self, master=None, cnf={}, **kw):
                '''when we get a file, we call fn(filename)'''
                self.fn = fn
                self.__toggle = 0
                Tkinter.Button.__init__(self, master, cnf, **kw)
                self.configure(command=self.set)

            def set(self):
                fd = PmwFileDialog(self.master, filter=filter)
                fd.title('Please choose a file')
                n = fd.askfilename()
                if n is not None:
                    self.fn(n)
        return FileDialogButton
    get = staticmethod(get)


class MOLETools:

    def setBinaryLocation(self, value):
        self.binlocation.setvalue(value)

    def __init__(self, app):
        parent = app.root
        self.parent = parent

        # Create the dialog.
        self.dialog = Pmw.Dialog(parent,
                                 buttons=('Show starting point', 'Run MOLE', 'Exit MOLE tools'),
                                 title = 'PyMOL MOLE Tools',
                                 command = self.execute)
        self.dialog.withdraw()
        Pmw.setbusycursorattributes(self.dialog.component('hull'))

        w = Tkinter.Label(self.dialog.interior(),
                          text='PyMOL MOLE Tools\nMartin Petrek, 2007, http://mole.chemi.muni.cz',
                          background='navy',
                          foreground='white',
                          #pady = 20,
                          )
        w.pack(expand=1, fill='both', padx=4, pady=4)

        self.notebook = Pmw.NoteBook(self.dialog.interior())
        self.notebook.pack(fill='both', expand=1, padx=10, pady=10)

        def quickFileValidation(s):
            if s == '':
                return Pmw.PARTIAL
            elif os.path.isfile(s):
                return Pmw.OK
            elif os.path.exists(s):
                return Pmw.PARTIAL
            else:
                return Pmw.PARTIAL

        # Set up the Main page
        page = self.notebook.add('Main')
        group = Pmw.Group(page, tag_text='Main options')
        group.pack(fill='both', expand=1, padx=10, pady=5)

        self.selection = Pmw.EntryField(group.interior(),
                                        labelpos='w',
                                        label_text='Selection to use: ',
                                        value='(all)',
                                        )

        self.binlocation = Pmw.EntryField(group.interior(),
                                          labelpos='w',
                                          label_pyclass=FileDialogButtonClassFactory.get(self.setBinaryLocation),
                                          validate={'validator': quickFileValidation, },
                                          value=MOLE_BINARY_LOCATION,
                                          label_text='Mole binary location:')

        self.numbertunnels = Pmw.EntryField(group.interior(),
                                            labelpos='w',
                                            label_text='Number of tunnels: ',
                                            value=str(defaults['numtunnels']),
                                            )

        for entry in (self.selection, self.binlocation, self.numbertunnels):
            entry.pack(fill='x', padx=4, pady=1)  # vertical

        labframe = Tkinter.Frame(group.interior())
        labframe.pack(fill='x', padx=4, pady=2)

        label1 = Tkinter.Label(labframe, justify=LEFT, text="Starting point specification:",)
        label1.pack(side='left')

        radiogroups = []
        self.var = Tkinter.IntVar()
        self.var.set(1)
        radioframe = Tkinter.Frame(group.interior())

        w = Pmw.Group(radioframe,
                      tag_pyclass=Tkinter.Radiobutton,
                      tag_text='Use average point from centers of given selections \n(e.g. /Molecule///GLY`37/N  /Molecule///PHE`151)',
                      tag_value=0,
                      tag_variable=self.var)
        w.pack(fill='x', expand=1, side='top')
        cw = Tkinter.Frame(w.interior())
        cw.pack(padx=2, pady=2, expand='yes', fill='both')
        radiogroups.append(w)
        self.selectionlist = Pmw.EntryField(w.interior(),
                                            labelpos='w',
                                            label_text='Selections List: ',
                                            value=defaults['residues'],
                                            #                                  value='(sele)',
                                            )

        self.selectionlist.pack(fill='x', padx=4, pady=1)  # vertical

        w = Pmw.Group(radioframe,
                      tag_pyclass=Tkinter.Radiobutton,
                      tag_text='Use X, Y, Z coordinates to specify starting point',
                      tag_value=1,
                      tag_variable=self.var)
        w.pack(fill='x', expand=1, side='top')
        cw = Tkinter.Frame(w.interior())
        cw.pack(padx=2, pady=2, expand='yes', fill='both')
        radiogroups.append(w)

# set starting point to center of molecule
        startpoint = self.computecenter1()
#        print startpoint
        self.xlocvar = DoubleVar()
        self.ylocvar = DoubleVar()
        self.zlocvar = DoubleVar()

        self.xlocvar.set(float("%.2f" % (startpoint[0])))
        self.ylocvar.set(float("%.2f" % (startpoint[1])))
        self.zlocvar.set(float("%.2f" % (startpoint[2])))
        self.var.set(1)

        self.xlocfr = Tkinter.Frame(w.interior())
        labX = Label(self.xlocfr, text="X:")
        self.xlocation = Entry(self.xlocfr, textvariable=self.xlocvar)
        self.scrX = Scrollbar(self.xlocfr, orient="horizontal", command=self.changeValueX)

        self.ylocfr = Tkinter.Frame(w.interior())
        labY = Label(self.ylocfr, text="Y:")
        self.ylocation = Entry(self.ylocfr, textvariable=self.ylocvar)
        self.scrY = Scrollbar(self.ylocfr, orient="horizontal", command=self.changeValueY)

        self.zlocfr = Tkinter.Frame(w.interior())
        labZ = Label(self.zlocfr, text="Z:")
        self.zlocation = Entry(self.zlocfr, textvariable=self.zlocvar)
        self.scrZ = Scrollbar(self.zlocfr, orient="horizontal", command=self.changeValueZ)

        labX.pack(side=LEFT)
        self.xlocation.pack(side=LEFT)
        self.scrX.pack(side=LEFT)
        self.xlocfr.pack(fill='x', padx=4, pady=1)  # vertical
        labY.pack(side=LEFT)
        self.ylocation.pack(side=LEFT)
        self.scrY.pack(side=LEFT)
        self.ylocfr.pack(fill='x', padx=4, pady=1)  # vertical
        labZ.pack(side=LEFT)
        self.zlocation.pack(side=LEFT)
        self.scrZ.pack(side=LEFT)
        self.zlocfr.pack(fill='x', padx=4, pady=1)  # vertical

        self.varremovewater = IntVar()
        self.removewaterbutton = Checkbutton(group.interior(), text="Remove water (HOH,WAT,H2O)", variable=self.varremovewater)

        radioframe.pack(padx=6, pady=6, expand='yes', fill='both')
        Pmw.aligngrouptags(radiogroups)

#        self.removewaterbutton.pack(padx = 6, pady = 6, expand='yes', fill='both')
        self.removewaterbutton.pack(side='left', expand='yes', padx=6, pady=6, fill='both')

# Other Configuration card
        page = self.notebook.add('Other Configuration')
        group = Pmw.Group(page, tag_text='Locations')
        group.pack(fill='both', expand=1, padx=10, pady=5)

        self.pymol_outdir = Pmw.EntryField(group.interior(),
                                           labelpos='w',
                                           label_text='Output directory: ',
                                           value=str(defaults['outdir']),
                                           )
        self.pymol_outdir.pack(fill='x', padx=20, pady=10)
        self.pymol_numinterp = Pmw.EntryField(group.interior(),
                                              labelpos='w',
                                              label_text='Number of point for tunnel axis interpolation: ',
                                              value=str(defaults['numinterp']),
                                              )
        self.pymol_numinterp.pack(fill='x', padx=20, pady=10)

        self.pymol_activesiteradius = Pmw.EntryField(group.interior(),
                                                     labelpos='w',
                                                     label_text='Radius for starting point optimization: ',
                                                     value=str(defaults['activesiteradius']),
                                                     )
        self.pymol_activesiteradius.pack(fill='x', padx=20, pady=10)

#        self.pymol_atomselection = Pmw.EntryField(group.interior(),
#                                                          labelpos = 'w',
#                                                          label_text = 'Selection of atoms to be used for calculation: ',
#                                                          value = str(defaults['atomselection']),
#                                                          )
#        self.pymol_atomselection.pack(fill = 'x', padx = 20, pady = 10)

# About card
        page = self.notebook.add('About')
        group = Pmw.Group(page, tag_text='About PyMOL MOLE Tools')
        group.pack(fill='both', expand=1, padx=10, pady=5)
        text = """
This plugin integrates PyMOL (http://pymol.org/) with MOLE.
Please cite the following reference when reporting the results using MOLE:

 Petrek M., Kosinova P., Koca J., Otyepka M.:
 MOLE: A Voronoi Diagram-Based Explorer of Molecular Channels,
 Pores, and Tunnels. Structure 2007, 15, 1357-1363

It should be fairly self-explanatory.  In the simplest case,

1) Load a protein structure into PyMOL.
2) Select Residues or Atoms making the starting point.
3) Start this plugin.
4) Make sure that the path to the MOLE binary is correct on the "MOLE Location" entry field.
5) Specify starting point by geometric center of selected atoms/residues
     (e.g. "/yourprotein///GLU`78/OE1  /yourprotein///GLN`35")
   You can also select residues by hand first and use pymol''s generated name (e.g. "(sele)"  )
6) Click "Show starting point" button. You should see small crisscross object in the PyMOL Viewer window.
7) Now you can switch to the choice "Use X, Y, Z coordinates to specify starting point" and you can state starting point
   more precisely by clicking arrows. You should move crisscross object in an empty space near active place in the molecule.
8) Click the "Run MOLE" button.
9) Wait some time.
10) See and enjoy results.

11) If no tunnel was found, you can try to remove waters from the structure by checking the checkbox.

For any bug, please send an e-mail with description to the author.

Details about used algorithm, tutorial for PyMOL plugin, examples you can find on the web page <http://mole.chemi.muni.cz/>.

Citation: see the web <http://mole.chemi.muni.cz/>.

Many thanks to
 - Warren DeLano for everything involving PyMOL
 - Michael Lerner <http://www.umich.edu/~mlerner/Pymol> for inspiration for Pymol Plugin design

Created by Martin Petrek petrek@chemi.muni.cz
LCC Group, NCBR Brno <http://ncbr.chemi.muni.cz/>"""

        lfre = Frame(group.interior())
        bar = Scrollbar(lfre,)
        ll = Text(lfre, yscrollcommand=bar.set, background="#ddddff", font="Times 12", width=60, height=15,)
        bar.config(command=ll.yview)

        ll.insert(END, text)
        ll.pack(side=LEFT, expand="yes", fill="both")
        bar.pack(side=LEFT, expand="no", fill="y")
        lfre.pack(expand="no", fill="both")

        self.notebook.setnaturalsize()
        self.showAppModal()

    def show_error1(self, message):
        error_dialog1 = Pmw.MessageDialog(self.parent, title='Error', message_text=message)
        junk1 = error_dialog1.activate()

    def showCrisscross(self):
        startpoint = (float(self.xlocvar.get()), float(self.ylocvar.get()), float(self.zlocvar.get()))
        cmd.delete("crisscross")
        self.crisscross(startpoint[0], startpoint[1], startpoint[2], 0.5, "crisscross")

    def changeValueX(self, obj, a, b=0, c=0):
        val = float(self.xlocvar.get()) + float(a) * 0.2
        self.xlocvar.set(val)
        self.showCrisscross()

    def changeValueY(self, obj, a, b=0, c=0):
        val = float(self.ylocvar.get()) + float(a) * 0.2
        self.ylocvar.set(val)
        self.showCrisscross()

    def changeValueZ(self, obj, a, b=0, c=0):
        val = float(self.zlocvar.get()) + float(a) * 0.2
        self.zlocvar.set(val)
        self.showCrisscross()

    def showAppModal(self):
        #self.dialog.activate(geometry = 'centerscreenalways', globalMode = 'nograb')
        self.dialog.show()
        #self.dialog.activate(geometry = 'centerscreenalways')

    def computecenter1(self, selection="(all)"):
        sel = cmd.get_model(selection)
        centx = 0
        centy = 0
        centz = 0
        cnt = len(sel.atom)
        if (cnt == 0):
            return (0, 0, 0)
        for a in sel.atom:
            centx += a.coord[0]
            centy += a.coord[1]
            centz += a.coord[2]
        centx /= cnt
        centy /= cnt
        centz /= cnt
#	fmttext="%lf\t%lf\t%lf\n" % (centx,centy,centz)
        return (centx, centy, centz)

    def computecenter(self, selection="(all)"):
        gcentx = 0
        gcenty = 0
        gcentz = 0
        gcnt = 0
        for selstr in selection.split():
            sel = cmd.get_model(selstr)

            centx = 0
            centy = 0
            centz = 0
            cnt = len(sel.atom)
            if (cnt == 0):
                self.show_error1("Bad starting point selection")
                return (0, 0, 0)
            for a in sel.atom:
                centx += a.coord[0]
                centy += a.coord[1]
                centz += a.coord[2]
            centx /= cnt
            centy /= cnt
            centz /= cnt
#	fmttext="%lf\t%lf\t%lf\n" % (centx,centy,centz)
#		print centx,centy,centz
            gcentx += centx
            gcenty += centy
            gcentz += centz
            gcnt += 1

        if (cnt == 0):
            self.show_error1("Bad starting point selection")
            return (0, 0, 0)

        gcentx /= gcnt
        gcenty /= gcnt
        gcentz /= gcnt
        return (gcentx, gcenty, gcentz)

    def prepareinput(self, selection='(all)', startpoint='(0,0,0)'):
        if self.varremovewater.get():
            selection = '(%s) and not resn HOH and not resn H2O and not resn WAT' % (selection)

        print(selection)
        outdir = self.pymol_outdir.getvalue()
        tmppdb = os.path.join(outdir, "tmp.pdb")
        cmd.save(tmppdb, selection)

    def removewater(self, selection='(all)'):
        seloutwater = Indexed()
        for a in selection.atom:
            if (a.resn != "HOH") and (a.resn != "H2O") and (a.resn != "WAT"):
                seloutwater.atom.append(a)
        return seloutwater

    def exportvdw(self, selection='(all)'):
        sel = cmd.get_model(selection)
        ext = cmd.get_extent(selection)
        mincorner = ext[0]
        maxcorner = ext[1]

        cnt = 0
        f = open("vystup.dat", "w")
#	f.write("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n" % (mincorner[0],mincorner[1],mincorner[2],maxcorner[0],maxcorner[1],maxcorner[2]))
        f.write("%d\n" % len(sel.atom))
        for a in sel.atom:
            x = a.coord[0]
            y = a.coord[1]
            z = a.coord[2]
            rad = a.vdw
            name = a.name
            text = "%s\t%lf\t%lf\t%lf\t%lf\n" % (name, x, y, z, rad)
            f.write(text)
            cnt += 1
        f.close()

    def deleteFileIfExist(self, path):
        if os.access(path, os.F_OK):
            os.remove(path)

    def deleteTemporaryFiles(self):
        #	self.deleteFileIfExist(self.pymol_generated_config_filename.getvalue())
        pass

    def execute(self, result):
        if result == 'Show starting point':
            choice = int(self.var.get())
            if (choice == 0):
                startpoint = self.computecenter(self.selectionlist.getvalue())
                self.xlocvar.set(float("%.2f" % (startpoint[0])))
                self.ylocvar.set(float("%.2f" % (startpoint[1])))
                self.zlocvar.set(float("%.2f" % (startpoint[2])))
                self.var.set(1)
            else:
                if (choice == 1):
                    startpoint = (float(self.xlocvar.get()), float(self.ylocvar.get()), float(self.zlocvar.get()))
            cmd.delete("crisscross")

            self.crisscross(startpoint[0], startpoint[1], startpoint[2], 0.5, "crisscross")
            return

        if result == 'Run MOLE':
            if self.checkInput(False) == False:
                return

#           self.exportvdw(self.selection.getvalue())
            choice = int(self.var.get())
            if (choice == 0):
                startpoint = self.computecenter(self.selectionlist.getvalue())
            else:
                if (choice == 1):
                    startpoint = (float(self.xlocvar.get()), float(self.ylocvar.get()), float(self.zlocvar.get()))

            print("Starting point of MOLE: %f %f %f\n" % startpoint)
            cmd.delete("crisscross")

            self.crisscross(startpoint[0], startpoint[1], startpoint[2], 0.5, "crisscross")

            self.deleteTemporaryFiles()
            outdir = self.pymol_outdir.getvalue()

            for i in range(1, int(self.numbertunnels.getvalue()) + 1):
                pathpy = os.path.join(outdir, "found_001_%03d.py" % i)
                self.deleteFileIfExist(pathpy)
                cmd.delete("tunnel%03d" % i)
# Now, prepare files for MOLE
            self.prepareinput(self.selection.getvalue(), startpoint)
            tmppdb = os.path.join(outdir, "tmp.pdb")
            mole_stdout = os.path.join(outdir, "mole.stdout")
            mole_stderr = os.path.join(outdir, "mole.stderr")
            if sys.platform.startswith('win'):
                command = '%s --fi pdb --input "%s" --numinterp %d --activesiteradius %f --noclustering --selection "all" --outdir "%s" --vmd --pymol --numtun %d --x %f --y %f --z %f >"%s" 2>"%s"' % (self.binlocation.getvalue(), tmppdb, int(self.pymol_numinterp.getvalue()), float(self.pymol_activesiteradius.getvalue()), self.pymol_outdir.getvalue(), int(self.numbertunnels.getvalue()), startpoint[0], startpoint[1], startpoint[2], mole_stdout, mole_stderr)
                print('\nRunning command:', command)
                status = subprocess.call(command)
            else:
                #command = [self.binlocation.getvalue(),'--fi','pdb','--input',tmppdb,'--numinterp',str(int(self.pymol_numinterp.getvalue())),'--activesiteradius',str(float(self.pymol_activesiteradius.getvalue())),'--noclustering','--selection','all','--outdir',self.pymol_outdir.getvalue(),'--vmd','--pymol','--numtun',str(int(self.numbertunnels.getvalue())),'--x',str(startpoint[0]),'--y',str(startpoint[1]),'--z',str(startpoint[2])]
                # print '\nRunning command:',command
                # print '\n pyenv is:', pymol_env
                #callfunc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=pymol_env)
                #child_stdout, child_stderr = callfunc.communicate()
                # callfunc.wait()
                #status= callfunc.returncode
                # print(child_stdout)
                # print(child_stderr)

                #command = '%s --fi pdb --input "%s" --numinterp %d --activesiteradius %f --noclustering --selection "all" --outdir "%s" --vmd --pymol --numtun %d --x %f --y %f --z %f >"%s" 2>"%s"'  % (self.binlocation.getvalue(),tmppdb,int(self.pymol_numinterp.getvalue()),float(self.pymol_activesiteradius.getvalue()),self.pymol_outdir.getvalue(),int(self.numbertunnels.getvalue()),startpoint[0],startpoint[1],startpoint[2],mole_stdout, mole_stderr);
                command = '%s --fi pdb --input %s --numinterp %d --activesiteradius %f --noclustering --selection "all" --outdir %s --vmd --pymol --numtun %d --x %f --y %f --z %f' % (self.binlocation.getvalue(), tmppdb, int(self.pymol_numinterp.getvalue()), float(self.pymol_activesiteradius.getvalue()), self.pymol_outdir.getvalue(), int(self.numbertunnels.getvalue()), startpoint[0], startpoint[1], startpoint[2])
                print('\nRunning command:', command)
                status = subprocess.call(command, shell=True)
            print("Result from execution was: %s" % status)
            view = cmd.get_view()
            for i in range(1, int(self.numbertunnels.getvalue()) + 1):
                pathpy = os.path.join(outdir, "found_001_%03d.py" % i)
                if os.access(pathpy, os.F_OK):
                    cmd.do("run %s" % pathpy)
                    # execfile(pathpy)
                    #callfunc = subprocess.Popen(pathpy, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=pymol_env)
                    # pymol.finish_launching()
                    #child_stdout, child_stderr = callfunc.communicate()
                    #result = callfunc.returncode
                    #print("Result from execution: %s"%result)
                    # print(child_stdout)
                    # print(child_stderr)
                else:
                    print("Error: Tunnel number %d wasn't found by MOLE." % i)
                    if sys.platform.startswith('win'):
                        print("Did you set your variables??")
                        print("PATH = PATH + ;%s" % qhull_dir)
                        print("MOLEDIR = %s" % mole_dir)
                    else:
                        print("Did you set your variables??")
                        print("PATH=$PATH:%s" % qhull_dir)
                        print("export PATH")
                        print("MOLEDIR=%s" % mole_dir)
                        print("export MOLEDIR")
            cmd.set_view(view)
            pass

#	    self.deleteTemporaryFiles()
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
                # self.dialog.deactivate(result)
                global MOLE_BINARY_LOCATION
                MOLE_BINARY_LOCATION = self.binlocation.getvalue()
                self.dialog.withdraw()

    def TestProgram(self, cmd):
        pathtoprog = os.path.join(os.path.split(MOLE_BINARY_LOCATION)[0], cmd)
        print('\nTesting command:', pathtoprog)
        if sys.platform.startswith('win'):
            status = subprocess.call(pathtoprog)
        else:
            callfunc = subprocess.Popen(pathtoprog, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=pymol_env)
            child_stdout, child_stderr = callfunc.communicate()
            status = callfunc.returncode
            # print(child_stdout)
            # print(child_stderr)
        print("Result from test was: %s" % status)
        return status

    def box(self, x1, y1, z1, x2, y2, z2, name="box"):

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
        cmd.load_cgo(obj, name)

    def crisscross(self, x, y, z, d, name="crisscross"):

        obj = [
            LINEWIDTH, 3,

            BEGIN, LINE_STRIP,
            VERTEX, float(x - d), float(y), float(z),
            VERTEX, float(x + d), float(y), float(z),
            END,

            BEGIN, LINE_STRIP,
            VERTEX, float(x), float(y - d), float(z),
            VERTEX, float(x), float(y + d), float(z),
            END,

            BEGIN, LINE_STRIP,
            VERTEX, float(x), float(y), float(z - d),
            VERTEX, float(x), float(y), float(z + d),
            END

        ]
        view = cmd.get_view()
        cmd.load_cgo(obj, name)
        cmd.set_view(view)

    def checkInput(self, silent=False):
        """If silent is True, we'll just return a True/False value
        """
        if not silent:
            def show_error(message):
                error_dialog = Pmw.MessageDialog(self.parent,
                                                 title='Error',
                                                 message_text=message,
                                                 )
                junk = error_dialog.activate()
        else:
            def show_error(message):
                pass

        #
        # First, check to make sure we have valid locations for MOLE
        #
        if not self.binlocation.valid():
            show_error('Please set the MOLE binary location')
            return False
        #
        # Now check if 'qconvex' can by run from MOLE as system('qconvex < ... > ...')
        #
        if sys.platform.startswith('win'):
            err = self.TestProgram("qhull.exe")
        else:
            err = self.TestProgram("qhull -")
        if err:
            print("Result from execution was: %s" % err)
            if sys.platform.startswith('win'):
                messgl = "\nDid you set your variables??\nPATH=PATH;%s\nMOLEDIR = %s" % (qhull_dir, mole_dir)
            else:
                messgl = "\nDid you set your variables??\nPATH=$PATH:%s\nexport PATH\nMOLEDIR=%s\nexport MOLEDIR" % (qhull_dir, mole_dir)
            errormessg = ('Can\'t run qhull (needed by MOLE). Please, add qhull in your PATH variable.%s' % messgl)
            print(errormessg)
            show_error(errormessg)
            return False
#
# Now check the temporary filenames
#

#        if not self.pymol_generated_config_filename.getvalue():
#            show_error('Please choose a name for the PyMOL\ngenerated config file')
#            return False
        return True

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
import os
import fnmatch
import time
import Pmw
# Pmw.setversion("0.8.5")


def _errorpop(master, text):
    d = Pmw.MessageDialog(master,
                        title="Error",
                        message_text=text,
                        buttons=("OK",))
    d.component('message').pack(ipadx=15, ipady=15)
    d.activate()
    d.destroy()


class PmwFileDialog(Pmw.Dialog):

    """File Dialog using Pmw"""

    def __init__(self, parent=None, **kw):
        # Define the megawidget options.
        optiondefs = (
            ('filter', '*', self.newfilter),
            ('directory', os.getcwd(), self.newdir),
            ('filename', '', self.newfilename),
            ('historylen', 10, None),
            ('command', None, None),
            ('info', None, None),
        )
        self.defineoptions(kw, optiondefs)
        # Initialise base class (after defining options).
        Pmw.Dialog.__init__(self, parent)

        self.withdraw()

        # Create the components.
        interior = self.interior()
        if self['info'] is not None:
            rowoffset = 1
            dn = self.infotxt()
            dn.grid(row=0, column=0, columnspan=2, padx=3, pady=3)
        else:
            rowoffset = 0

        dn = self.mkdn()
        dn.grid(row=0 + rowoffset, column=0, columnspan=2, padx=3, pady=3)
        del dn

        # Create the directory list component.
        dnb = self.mkdnb()
        dnb.grid(row=1 + rowoffset, column=0, sticky='news', padx=3, pady=3)
        del dnb

        # Create the filename list component.
        fnb = self.mkfnb()
        fnb.grid(row=1 + rowoffset, column=1, sticky='news', padx=3, pady=3)
        del fnb

        # Create the filter entry
        ft = self.mkft()
        ft.grid(row=2 + rowoffset, column=0, columnspan=2, padx=3, pady=3)
        del ft

        # Create the filename entry
        fn = self.mkfn()
        fn.grid(row=3 + rowoffset, column=0, columnspan=2, padx=3, pady=3)
        fn.bind('<Return>', self.okbutton)
        del fn

        # Buttonbox already exists
        bb = self.component('buttonbox')
        bb.add('OK', command=self.okbutton)
        bb.add('Cancel', command=self.cancelbutton)
        del bb

        Pmw.alignlabels([self.component('filename'),
                         self.component('filter'),
                         self.component('dirname')])

        self.lastdir = ""
        self.lastfilter = None
        self.lasttime = 0

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

    def dirvalidate(self, string):
        if os.path.isdir(string):
            return Pmw.OK
        else:
            return Pmw.PARTIAL

    def filevalidate(self, string):
        if string == '':
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
        fn = self.component('filename').get()
        self.setfilename(fn)
        if self.validate(fn):
            self.canceled = 0
            self.deactivate()

    def cancelbutton(self):
        """Cancel the operation"""
        self.canceled = 1
        self.deactivate()

    def tidy(self, w, v):
        """Insert text v into the entry and at the top of the list of
           the combobox w, remove duplicates"""
        if not v:
            return
        entry = w.component('entry')
        entry.delete(0, 'end')
        entry.insert(0, v)
        list = w.component('scrolledlist')
        list.insert(0, v)
        index = 1
        while index < list.index('end'):
            k = list.get(index)
            if k == v or index > self['historylen']:
                list.delete(index)
            else:
                index = index + 1
        w.checkentry()

    def setfilename(self, value):
        if not value:
            return
        value = os.path.join(self['directory'], value)
        dir, fil = os.path.split(value)
        self.configure(directory=dir, filename=value)

        c = self['command']
        if callable(c):
            c()

    def newfilename(self):
        """Make sure a newly set filename makes it into the combobox list"""
        self.tidy(self.component('filename'), self['filename'])

    def setfilter(self, value):
        self.configure(filter=value)

    def newfilter(self):
        """Make sure a newly set filter makes it into the combobox list"""
        self.tidy(self.component('filter'), self['filter'])
        self.fillit()

    def setdir(self, value):
        self.configure(directory=value)

    def newdir(self):
        """Make sure a newly set dirname makes it into the combobox list"""
        self.tidy(self.component('dirname'), self['directory'])
        self.fillit()

    def singleselectfile(self):
        """Single click in file listbox. Move file to "filename" combobox"""
        cs = self.component('filenamebox').curselection()
        if cs != ():
            value = self.component('filenamebox').get(cs)
            self.setfilename(value)

    def selectfile(self):
        """Take the selected file from the filename, normalize it, and OK"""
        self.singleselectfile()
        value = self.component('filename').get()
        self.setfilename(value)
        if value:
            self.okbutton()

    def selectdir(self):
        """Take selected directory from the dirnamebox into the dirname"""
        cs = self.component('dirnamebox').curselection()
        if cs != ():
            value = self.component('dirnamebox').get(cs)
            dir = self['directory']
            if not dir:
                dir = os.getcwd()
            if value:
                if value == '..':
                    dir = os.path.split(dir)[0]
                else:
                    dir = os.path.join(dir, value)
            self.configure(directory=dir)
            self.fillit()

    def askfilename(self, directory=None, filter=None):
        """The actual client function. Activates the dialog, and
        returns only after a valid filename has been entered
        (return value is that filename) or when canceled (return
        value is None)"""
        if directory != None:
            self.configure(directory=directory)
        if filter != None:
            self.configure(filter=filter)
        self.fillit()
        self.canceled = 1  # Needed for when user kills dialog window
        self.activate()
        if self.canceled:
            return None
        else:
            return self.component('filename').get()

    def fillit(self):
        """Get the directory list and show it in the two listboxes"""
        # Do not run unnecesarily
        if self.lastdir == self['directory'] and self.lastfilter == self['filter'] and self.lasttime > os.stat(self.lastdir)[8]:
            return
        self.lastdir = self['directory']
        self.lastfilter = self['filter']
        self.lasttime = time.time()
        dir = self['directory']
        if not dir:
            dir = os.getcwd()
        dirs = ['..']
        files = []
        try:
            fl = os.listdir(dir)
            fl.sort()
        except os.error as arg:
            if arg[0] in (2, 20):
                return
            raise
        for f in fl:
            if os.path.isdir(os.path.join(dir, f)):
                dirs.append(f)
            else:
                filter = self['filter']
                if not filter:
                    filter = '*'
                if fnmatch.fnmatch(f, filter):
                    files.append(f)
        self.component('filenamebox').setlist(files)
        self.component('dirnamebox').setlist(dirs)

    def validate(self, filename):
        """Validation function. Should return 1 if the filename is valid,
           0 if invalid. May pop up dialogs to tell user why. Especially
           suited to subclasses: i.e. only return 1 if the file does/doesn't
           exist"""
        return 1


class PmwExistingFileDialog(PmwFileDialog):

    def filevalidate(self, string):
        if os.path.isfile(string):
            return Pmw.OK
        else:
            return Pmw.PARTIAL

    def validate(self, filename):
        if os.path.isfile(filename):
            return 1
        elif os.path.exists(filename):
            _errorpop(self.interior(), "This is not a plain file")
            return 0
        else:
            _errorpop(self.interior(), "Please select an existing file")
            return 0

# Create demo in root window for testing.
if __name__ == '__main__':
    class App:

        def my_show(self, *args, **kwargs):
            pass

    app = App()
    app.root = Tkinter.Tk()
    Pmw.initialise(app.root)
    app.root.title('Some Title')

    widget = MOLETools(app)
    exitButton = Tkinter.Button(app.root, text='Exit', command=app.root.destroy)
    exitButton.pack()
    app.root.mainloop()
