""" 2010_04_15: Hongbo Zhu
      PyMOL plugin for displaying MSMS surface.
"""

# Copyright Notice
# ================
#
# The PyMOL Plugin source code in this file is copyrighted, but you can
# freely use and copy it as long as you don't change or remove any of
# the copyright notices.
#
# ----------------------------------------------------------------------
#               This PyMOL Plugin is Copyright (C) 2010 by
#                 Hongbo Zhu <macrozhu at gmail dot com>
#
#                        All Rights Reserved
#
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name(s) of the author(s) not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
#
# THE AUTHOR(S) DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
# INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.  IN
# NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY SPECIAL, INDIRECT OR
# CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
# USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
# OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
# PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------

from __future__ import absolute_import
from __future__ import print_function

# python lib
import os
import sys
import platform
import subprocess
import time
if sys.version_info[0] < 3:
    import Tkinter
    import tkMessageBox
    import tkFileDialog
    import tkColorChooser
else:
    import tkinter as Tkinter
    import tkinter.messagebox as tkMessageBox
    import tkinter.filedialog as tkFileDialog
    import tkinter.colorchooser as tkColorChooser

# pymol lib
try:
    from pymol import cmd
    from pymol.cgo import *
except ImportError:
    print('Warning: pymol library cmd not found.')
    sys.exit(1)

# external lib
try:
    import Pmw
except ImportError:
    print('Warning: failed to import Pmw. Exit ...')
    sys.exit(1)

VERBOSE = True

#################
# here we go
#################


def __init__(self):
    """ MSMS plugin for PyMol
    """
    self.menuBar.addmenuitem('Plugin', 'command',
                             'MSMS', label='MSMS',
                             command=lambda s=self: MSMSPlugin(s))


#################
# GUI related
#################
class MSMSPlugin:

    def __init__(self, app):
        self.parent = app.root
        self.dialog = Pmw.Dialog(self.parent,
                                 buttons=('Run MSMS', 'Display Mesh',
                                            'Display Vertices', 'Exit'),
                                 title = 'MSMS Plugin for PyMOL',
                                 command = self.execute)
        Pmw.setbusycursorattributes(self.dialog.component('hull'))

        # parameters used by MSMS
        self.pdb_fn = Tkinter.StringVar()
        self.pymol_sel = Tkinter.StringVar()
        self.msms_bin = Tkinter.StringVar()
##         self.pdb2xyzr_bin  = Tkinter.StringVar()
        self.pdb2xyzrn_bin = Tkinter.StringVar()
        self.tmp_dir = Tkinter.StringVar()

        self.cleanup_saved_pymol_sel = Tkinter.BooleanVar()
        self.cleanup_saved_pymol_sel.set(True)  # by default, clean up

        self.pdb_fn.set('')
        if 'MSMS_BIN' not in os.environ and 'PYMOL_GIT_MOD' in os.environ:
            if sys.platform.startswith('linux') and platform.machine() == 'x86_32':
                initialdir_msms = os.path.join(os.environ['PYMOL_GIT_MOD'], "MSMS", "i86Linux2", "msms.i86Linux2.2.6.1")
                os.environ['MSMS_BIN'] = initialdir_msms
            elif sys.platform.startswith('linux') and platform.machine() == 'x86_64':
                initialdir_msms = os.path.join(os.environ['PYMOL_GIT_MOD'], "MSMS", "i64Linux2", "msms.x86_64Linux2.2.6.1")
                os.environ['MSMS_BIN'] = initialdir_msms
            elif sys.platform.startswith('darwin'):
                initialdir_msms = os.path.join(os.environ['PYMOL_GIT_MOD'], "MSMS", "universalDarwin", "msms.MacOSX.2.6.1")
                os.environ['MSMS_BIN'] = initialdir_msms
            elif sys.platform.startswith('win'):
                initialdir_msms = os.path.join(os.environ['PYMOL_GIT_MOD'], "MSMS", "win32", "msms_win32_6.2.1", "msms.exe")
                os.environ['MSMS_BIN'] = initialdir_msms
            else:
                pass
        if 'MSMS_BIN' in os.environ:
            if VERBOSE:
                print('Found MSMS_BIN in environmental variables', os.environ['MSMS_BIN'])
            self.msms_bin.set(os.environ['MSMS_BIN'])
        else:
            if VERBOSE:
                print('MSMS_BIN not found in environmental variables.')
            self.msms_bin.set('')
# self.pdb2xyzr_bin.set('')
        if 'PDB2XYZRN' not in os.environ and 'PYMOL_GIT_MOD' in os.environ:
            if sys.platform.startswith('linux') and platform.machine() == 'x86_32':
                initialdir_msms = os.path.join(os.environ['PYMOL_GIT_MOD'], "MSMS", "i86Linux2", "pdb_to_xyzrn")
                os.environ['PDB2XYZRN'] = initialdir_msms
            elif sys.platform.startswith('linux') and platform.machine() == 'x86_64':
                initialdir_msms = os.path.join(os.environ['PYMOL_GIT_MOD'], "MSMS", "i64Linux2", "pdb_to_xyzrn")
                os.environ['PDB2XYZRN'] = initialdir_msms
            elif sys.platform.startswith('darwin'):
                initialdir_msms = os.path.join(os.environ['PYMOL_GIT_MOD'], "MSMS", "universalDarwin", "pdb_to_xyzrn")
                os.environ['PDB2XYZRN'] = initialdir_msms
            elif sys.platform.startswith('win'):
                initialdir_msms = os.path.join(os.environ['PYMOL_GIT_MOD'], "MSMS", "win32", "msms_win32_6.2.1", "pdb2xyzrn.py")
                os.environ['PDB2XYZRN'] = initialdir_msms
            else:
                pass
        if 'PDB2XYZRN' in os.environ:
            self.pdb2xyzrn_bin.set(os.environ['PDB2XYZRN'])
        else:
            self.pdb2xyzrn_bin.set('')
        cmd.select("protpolymer", "polymer")
        cmd.disable("protpolymer")
        self.pymol_sel.set(cmd.get_names('selections')[-1])
        # self.tmp_dir.set('/tmp')
        self.tmp_dir.set(os.getcwd())
        self.cleanup_msms_output = Tkinter.BooleanVar()
        self.cleanup_msms_output.set(True)  # by default, clean up msms output

        # MSMS parameters
        self.probe_radius = Tkinter.DoubleVar()
        self.density = Tkinter.DoubleVar()
        self.hdensity = Tkinter.DoubleVar()
        self.noh = Tkinter.BooleanVar()  # ignore hydrogen atoms
        self.allcpn = Tkinter.BooleanVar()  # consider all surface components
        self.probe_radius.set(1.5)
        self.density.set(1.0)
        self.hdensity.set(3.0)
        self.noh.set(True)      # by default ignore hydrogens
        self.allcpn.set(False)  # by default consider all components

        # MSMS output
        self.msms_vert_fn = None  # external surface
        self.msms_face_fn = None
        self.msms_cpn_vert_fn_list = []  # internal components
        self.msms_cpn_face_fn_list = []

        # MSMSSurfPymol object
        self.msp = MSMSSurfPymol()
        self.cpn_msp_list = []  # MSMSSurfPymol objects for internal components

        # MSMS visualization color
        self.mesh_col = '#ffffff'
        self.vert_col = '#ffffff'
        self.norm_col = '#ffb432'
        self.mesh_col_R = 255  # mesh color
        self.mesh_col_G = 255
        self.mesh_col_B = 255
        self.vert_col_R = 255  # vertex color
        self.vert_col_G = 255
        self.vert_col_B = 255
        self.norm_col_R = 255  # normal vector color (orange)
        self.norm_col_G = 180
        self.norm_col_B = 50

        self.vert_rad = Tkinter.DoubleVar()  # radius for spheres representing vertices
        self.norm_len = Tkinter.DoubleVar()  # length of normal vectors
        self.vert_rad.set(0.2)
        self.norm_len.set(1.0)

        w = Tkinter.Label(self.dialog.interior(),
                          text='\nMSMS Plugin for PyMOL\nHongbo Zhu, 2010.\n\nDisplaying protein surface calculated by MSMS.',
                          background='black', foreground='green'
                          )
        w.pack(expand=1, fill='both', padx=10, pady=5)

        # make a few tabs within the dialog
        self.notebook = Pmw.NoteBook(self.dialog.interior())
        self.notebook.pack(fill='both', expand=1, padx=10, pady=10)

        ######################
        # Tab : Structure Tab
        ######################
        page = self.notebook.add('Structure')
        self.notebook.tab('Structure').focus_set()
        group_struc = Tkinter.LabelFrame(page, text='Structure')
        group_struc.pack(fill='both', expand=True, padx=10, pady=5)

        pymol_sel_ent = Pmw.EntryField(group_struc,
                                       label_text='PyMOL selection:',
                                       labelpos='wn',
                                       entry_textvariable=self.pymol_sel
                                       )
        clean_cb = Tkinter.Checkbutton(group_struc,
                                       text='Clean up tmp pdb (saved PyMOL selection) in the temp dir.',
                                       variable=self.cleanup_saved_pymol_sel,
                                       onvalue=True, offvalue=False)
        label = Tkinter.Label(group_struc, text='or')

        pdb_fn_ent = Pmw.EntryField(group_struc,
                                    label_text='PDB file:', labelpos='wn',
                                    entry_textvariable=self.pdb_fn)
        pdb_fn_but = Tkinter.Button(group_struc, text='Browse...',
                                    command=self.getPDBFile)

        # arrange widgets using grid
        pymol_sel_ent.grid(sticky='we', row=0, column=0,
                           columnspan=2, padx=5, pady=5)
        clean_cb.grid(sticky='w', row=1, column=0,
                      columnspan=2, padx=1, pady=1)
        label.grid(sticky='we', row=2, column=0, columnspan=2, padx=5, pady=10)
        pdb_fn_ent.grid(sticky='we', row=3, column=0, padx=5, pady=5)
        pdb_fn_but.grid(sticky='we', row=3, column=1, padx=5, pady=5)
        group_struc.columnconfigure(0, weight=9)
        group_struc.columnconfigure(1, weight=1)

        ######################
        # Tab : MSMS Configuration
        ######################
        page = self.notebook.add('MSMS Configuration')

        group_loc = Tkinter.LabelFrame(page, text='Locations')
        group_msms_param = Tkinter.LabelFrame(page, text='Parameters')
        group_loc.grid(sticky='eswn', row=0, column=0, columnspan=2, padx=10, pady=5)
        group_msms_param.grid(sticky='eswn', row=0, column=2, padx=10, pady=5)
        page.columnconfigure(0, weight=2)
        page.columnconfigure(1, weight=1)

        msms_bin_ent = Pmw.EntryField(group_loc,
                                      label_text='MSMS binary:', labelpos='wn',
                                      entry_textvariable=self.msms_bin,
                                      entry_width=20)
        msms_bin_but = Tkinter.Button(group_loc, text='Browse...',
                                      command=self.getMsmsBin)

# pdb2xyzr_bin_ent = Pmw.EntryField(group_loc,
# label_text='pdb2xyzr binary:', labelpos='wn',
# entry_textvariable=self.pdb2xyzr_bin,
# entry_width=20)
# pdb2xyzr_bin_but = Tkinter.Button(group_loc, text = 'Browse...',
# command = self.getPdb2xyzrBin)
        pdb2xyzrn_bin_ent = Pmw.EntryField(group_loc,
                                           label_text='pdb2xyzrn binary:', labelpos='wn',
                                           entry_textvariable=self.pdb2xyzrn_bin,
                                           entry_width=20)
        pdb2xyzrn_bin_but = Tkinter.Button(group_loc, text='Browse...',
                                           command=self.getPdb2xyzrnBin)

        tmp_dir_ent = Pmw.EntryField(group_loc,
                                     label_text='Temporary dir:', labelpos='wn',
                                     entry_textvariable=self.tmp_dir,
                                     entry_width=20)
        tmp_dir_but = Tkinter.Button(group_loc, text='Browse...',
                                     command=self.getTmpDir)
        ko_cb = Tkinter.Checkbutton(group_loc,
                                    text='Clean up MSMS output (.vert and .face files) in the temp dir.',
                                    variable=self.cleanup_msms_output,
                                    onvalue=True, offvalue=False)

        msms_bin_ent.grid(sticky='we', row=0, column=0, padx=5, pady=1)
        msms_bin_but.grid(sticky='we', row=0, column=1, padx=5, pady=1)
        pdb2xyzrn_bin_ent.grid(sticky='we', row=1, column=0, padx=5, pady=1)
        pdb2xyzrn_bin_but.grid(sticky='we', row=1, column=1, padx=5, pady=1)
        tmp_dir_ent.grid(sticky='we', row=2, column=0, padx=5, pady=1)
        tmp_dir_but.grid(sticky='we', row=2, column=1, padx=5, pady=1)
        ko_cb.grid(sticky='w', row=3, column=0, padx=1, pady=1)
        group_loc.columnconfigure(0, weight=9)
        group_loc.columnconfigure(1, weight=1)

        pr_ent = Pmw.EntryField(group_msms_param, labelpos='wn',
                                label_text='Probe radius:',
                                value=self.probe_radius.get(),
                                validate={'validator': 'real', 'min': 0},
                                entry_textvariable=self.probe_radius,
                                entry_width=10
                                )
        # density can be set to be as low as 0.25.
        # MSMS fails to generate output for density values lower than 0.25
        # (at least for some PDB files like 1BID).
        den_ent = Pmw.EntryField(group_msms_param, labelpos='wn',
                                 label_text='Density:', value=self.density.get(),
                                 validate={'validator': 'real', 'min': 0.25},
                                 entry_textvariable=self.density,
                                 entry_width=10
                                 )
        hden_ent = Pmw.EntryField(group_msms_param, labelpos='wn',
                                  label_text='High density:', value=self.hdensity.get(),
                                  validate={'validator': 'real', 'min': 0.25},
                                  entry_textvariable=self.hdensity,
                                  entry_width=10
                                  )
        noh_cb = Tkinter.Checkbutton(group_msms_param,
                                     text='Ignore hydrogens.',
                                     variable=self.noh,
                                     onvalue=True, offvalue=False)
        allcpn_cb = Tkinter.Checkbutton(group_msms_param,
                                        text='Consider all surface components.',
                                        variable=self.allcpn,
                                        onvalue=True, offvalue=False)
        pr_ent.grid(sticky='we', row=0, column=0, padx=5, pady=1)
        den_ent.grid(sticky='we', row=1, column=0, padx=5, pady=1)
        hden_ent.grid(sticky='we', row=2, column=0, padx=5, pady=1)
        noh_cb.grid(sticky='w', row=3, column=0, padx=0, pady=1)
        allcpn_cb.grid(sticky='w', row=4, column=0, padx=0, pady=1)

        group_msms_param.columnconfigure(0, weight=1)
        group_msms_param.columnconfigure(1, weight=1)

        ######################
        # Tab : Visualization Tab
        ######################
        page = self.notebook.add('Visualization')

        group_vis_msms = Tkinter.LabelFrame(page, text='MSMS Surface')
        group_vis_msms.grid(sticky='eswn', row=0, column=0, padx=10, pady=5)

        # colors for MSMS surface
        mesh_col_lab = Tkinter.Label(group_vis_msms, text='Surface mesh color:')
        self.mesh_col_but = Tkinter.Button(group_vis_msms,
                                           bg=self.mesh_col,
                                           activebackground=self.mesh_col,
                                           command=self.custermizeMeshColor)
        vert_col_lab = Tkinter.Label(group_vis_msms, text='Surface vertex color:')
        self.vert_col_but = Tkinter.Button(group_vis_msms,
                                           bg=self.vert_col,
                                           activebackground=self.vert_col,
                                           command=self.custermizeVertColor)

        norm_col_lab = Tkinter.Label(group_vis_msms, text='Normal vector color:')
        self.norm_col_but = Tkinter.Button(group_vis_msms,
                                           bg=self.norm_col,
                                           activebackground=self.norm_col,
                                           command=self.custermizeNormColor)
        vert_rad_ent = Pmw.EntryField(group_vis_msms, labelpos='wn',
                                      label_text='Surface vertex radius:',
                                      value=self.vert_rad.get(),
                                      validate={'validator': 'real', 'min': 0.00},
                                      entry_textvariable=self.vert_rad
                                      )

        norm_len_ent = Pmw.EntryField(group_vis_msms, labelpos='wn',
                                      label_text='Normal vector length:',
                                      value=self.norm_len.get(),
                                      validate={'validator': 'real', 'min': 0.00},
                                      entry_textvariable=self.norm_len
                                      )

        mesh_col_lab.grid(sticky='e', row=0, column=0, padx=5, pady=1)
        self.mesh_col_but.grid(sticky='we', row=0, column=1, padx=5, pady=1)
        vert_col_lab.grid(sticky='e', row=1, column=0, padx=5, pady=1)
        self.vert_col_but.grid(sticky='we', row=1, column=1, padx=5, pady=1)
        norm_col_lab.grid(sticky='e', row=2, column=0, padx=5, pady=1)
        self.norm_col_but.grid(sticky='we', row=2, column=1, padx=5, pady=1)
        vert_rad_ent.grid(sticky='e', row=3, column=0, columnspan=2, padx=5, pady=1)
        norm_len_ent.grid(sticky='e', row=4, column=0, columnspan=2, padx=5, pady=1)

        ######################
        # Tab : About Tab
        ######################
        page = self.notebook.add('About')
        group_about = Tkinter.LabelFrame(page, text='About MSMS Plugin for PyMOL')
        group_about.grid(sticky='we', row=0, column=0, padx=10, pady=5)
        about_plugin = """ This plugin provides a GUI for running MSMS and displaying its results in PyMOL.
Created by Hongbo Zhu <hongbo.zhu.cn@googlemail.com>, Biotechnology Center (BIOTEC), TU Dresden.

Please cite this plugin if you use it in a publication.
Citation for this plugin:
Hongbo Zhu. MSMS plugin for PyMOL, 2010, Biotechnology Center (BIOTEC), TU Dresden.
Citation for PyMOL can be found at:
http://pymol.sourceforge.net/faq.html#CITE
"""

        label_about = Tkinter.Label(group_about, text=about_plugin)
        label_about.grid(sticky='we', row=0, column=0, padx=5, pady=10)

        self.notebook.setnaturalsize()

        return

    def getPDBFile(self):
        file_name = tkFileDialog.askopenfilename(
            title='PDB File', initialdir='',
            filetypes=[('pdb files', '*.pdb *.ent'), ('all files', '*')],
            parent=self.parent)
        if file_name:
            self.pdb_fn.set(file_name)

    def getMsmsBin(self):
        msms_bin_fname = tkFileDialog.askopenfilename(
            title='MSMS Binary', initialdir='',
            filetypes=[('all', '*')], parent=self.parent)
        if msms_bin_fname:
            self.msms_bin.set(msms_bin_fname)

# def getPdb2xyzrBin(self):
# pdb2xyzr_bin_fname = tkFileDialog.askopenfilename(
# title='pdb2xyzr Binary', initialdir='',
# filetypes=[('all','*')], parent=self.parent)
# self.pdb2xyzr_bin.set(pdb2xyzr_bin_fname)

    def getPdb2xyzrnBin(self):
        pdb2xyzrn_bin_fname = tkFileDialog.askopenfilename(
            title='pdb2xyzrn Binary', initialdir='',
            filetypes=[('all', '*')], parent=self.parent)
        if pdb2xyzrn_bin_fname:
            self.pdb2xyzrn_bin.set(pdb2xyzrn_bin_fname)

    def getTmpDir(self):
        tmp_dir = tkFileDialog.askdirectory()
        if tmp_dir:
            self.tmp_dir.set(tmp_dir)

    def getStrucPDBFname(self):
        """ get the PDB file name for the structure to work on
            if the structure is specified by a pymol selection, save it in the temp dir;
            if the structure is specified by a separate PDB file, use it.
        """
        pdb_fn = None
        sel = self.pymol_sel.get()

        if len(sel) > 0:  # if any pymol selection is specified
            # save the pymol selection in the tmp dir
            all_sel_names = cmd.get_names('selections')  # get names of all selections

            tmp_dir = self.tmp_dir.get()
            if tmp_dir[-1] == '/' or tmp_dir[-1] == '\\':
                tmp_dir = tmp_dir[:-1]
                self.tmp_dir.set(tmp_dir)

            if sel in all_sel_names:
                # pdb_fn = '%s/pymol_sele_%s_%s.pdb' % (self.tmp_dir.get(), sel,
                #                                      str(time.time()).replace('.',''))
                pdb_fn = os.path.join(self.tmp_dir.get(), "pymol_sele_%s_%s.pdb" % (sel, str(time.time()).replace('.', '')))
                cmd.save(filename=pdb_fn, selection=sel)
                if VERBOSE:
                    print('Selection %s saved to %s.' % (sel, pdb_fn))
            else:  # sel is unknown
                err_msg = 'Unknown selection name: %s' % (sel,)
                print('ERROR: %s' % (err_msg,))
                tkMessageBox.showinfo(title='ERROR', message=err_msg)

        elif len(self.pdb_fn.get()) > 0:  # if no selection specified, use specified PDB file
            pdb_fn = self.pdb_fn.get()
            if not os.path.isfile(pdb_fn):
                err_msg = 'PDB file does not exist: %s' % (pdb_fn,)
                print('ERROR: %s' % (err_msg,))
                tkMessageBox.showinfo(title='ERROR', message=err_msg)

        else:   # what structure do you want MSMS to work on?
            err_msg = 'Neither PyMOL selection nor PDB file specified!'
            print('ERROR: %s' % (err_msg,))
            tkMessageBox.showinfo(title='ERROR', message=err_msg)

        return pdb_fn

    def cleanMSMSOutput(self):

        if os.path.isfile(self.msms_vert_fn):
            if VERBOSE:
                print('Cleaning msms vert file', self.msms_vert_fn)
            os.remove(self.msms_vert_fn)
            self.msms_vert_fn = None
        if os.path.isfile(self.msms_face_fn):
            if VERBOSE:
                print('Cleaning msms face file', self.msms_face_fn)
            os.remove(self.msms_face_fn)
            self.msms_face_fn = None

        for i in range(len(self.msms_cpn_vert_fn_list)):
            vfn = self.msms_cpn_vert_fn_list[i]
            ffn = self.msms_cpn_face_fn_list[i]
            if os.path.isfile(vfn):
                if VERBOSE:
                    print('Cleaning msms face file', vfn)
                os.remove(vfn)
            if os.path.isfile(ffn):
                if VERBOSE:
                    print('Cleaning msms face file', ffn)
                os.remove(ffn)

        self.msms_cpn_vert_fn_list = []
        self.msms_cpn_face_fn_list = []

        return

    def runMSMS(self):
        """
            @return: whether MSMS has been executed successfully
            @rtype: boolean
        """
        # clean up old results, which might be from previous execution
        self.msms_vert_fn = None  # external surface
        self.msms_face_fn = None
        self.msms_cpn_vert_fn_list = []  # internal components
        self.msms_cpn_face_fn_list = []
        self.msp = MSMSSurfPymol()
        self.cpn_msp_list = []  # MSMSSurfPymol objects for internal components

        tmp_dir = self.tmp_dir.get()
        if tmp_dir[-1] == '/' or tmp_dir[-1] == '\\':
            tmp_dir = tmp_dir[:-1]
            self.tmp_dir.set(tmp_dir)

        if VERBOSE:
            print('MSMS bin  =', self.msms_bin.get())
# print self.pdb2xyzr_bin.get()
            print('pdb2xyzrn =', self.pdb2xyzrn_bin.get())
            print('tmp dir   =', self.tmp_dir.get())

        pdb_fn = self.getStrucPDBFname()
        if pdb_fn is None:
            return None

        print('Running MSMS ...')
        if VERBOSE:
            print('Probe raidus            =', self.probe_radius.get())
            print('Surface vertex density  =', self.density.get())
            print('Surface vertex hdensity =', self.hdensity.get())
            print('Ignore hydrogen atoms   =', str(self.noh.get()))
            print('Consider all surface components =', str(self.allcpn.get()))

        msms = Msms(msms_bin=self.msms_bin.get(),
                    # pdb2xyzr_bin=self.pdb2xyzr_bin.get(),
                    pdb2xyzrn_bin=self.pdb2xyzrn_bin.get(),
                    pr=self.probe_radius.get(),
                    den=self.density.get(),
                    hden=self.hdensity.get(),
                    noh=self.noh.get(), all_components=self.allcpn.get(),
                    output_dir=self.tmp_dir.get()
                    )
        msms.run(pdb_fn)
        print('done!')

        # remove temp file (saved pymol selection)
        if self.cleanup_saved_pymol_sel and \
                len(self.pymol_sel.get()) > 0 and os.path.isfile(pdb_fn):
            if VERBOSE:
                print('Cleaning temp file(s)', pdb_fn)
            # !os.remove(pdb_fn)  # clean up (remove pdb file of the pymol selection)

        fn_list = msms.getOutputFiles()
        self.msms_vert_fn = fn_list[0]
        self.msms_face_fn = fn_list[1]

        # make copies. do not use reference
        [self.msms_cpn_vert_fn_list.append(fni) for fni in fn_list[2]]
        [self.msms_cpn_face_fn_list.append(fni) for fni in fn_list[3]]

        if VERBOSE:
            print('MSMS .vert file:', self.msms_vert_fn)
            print('MSMS .face file:', self.msms_face_fn)

        return msms

    def custermizeMeshColor(self):
        try:
            color_tuple, color = tkColorChooser.askcolor(color=self.mesh_col)
            if color_tuple is not None and color is not None:
                self.mesh_col_R, self.mesh_col_G, self.mesh_col_B = color_tuple
                self.mesh_col = color
                self.mesh_col_but['bg'] = self.mesh_col
                self.mesh_col_but['activebackground'] = self.mesh_col
                self.mesh_col_but.update()
        except Tkinter._tkinter.TclError:
            print('Old color (%s) will be used.' % (self.mesh_col))

        return

    def custermizeVertColor(self):
        try:
            color_tuple, color = tkColorChooser.askcolor(color=self.vert_col)
            if color_tuple is not None and color is not None:
                self.vert_col_R, self.vert_col_G, self.vert_col_B = color_tuple
                self.vert_col = color
                self.vert_col_but['bg'] = self.vert_col
                self.vert_col_but['activebackground'] = self.vert_col
                self.vert_col_but.update()
        except Tkinter._tkinter.TclError:
            print('Old color (%s) will be used.' % (self.vert_col))

        return

    def custermizeNormColor(self):
        try:
            color_tuple, color = tkColorChooser.askcolor(color=self.norm_col)
            if color_tuple is not None and color is not None:
                self.norm_col_R, self.norm_col_G, self.norm_col_B = color_tuple
                self.norm_col = color
                self.norm_col_but['bg'] = self.norm_col
                self.norm_col_but['activebackground'] = self.norm_col
                self.norm_col_but.update()
        except Tkinter._tkinter.TclError:
            print('Old color (%s) will be used.' % (self.norm_col))

        return

    def execute(self, cmd):
        """ Run the cmd represented by the botton clicked by user.
        """
        if cmd == 'OK':
            print('is everything OK?')

        elif cmd == 'Run MSMS':

            if self.runMSMS() is not None:  # msms has been executed successfully

                if VERBOSE:
                    print('Parsing MSMS output ...')
                self.msp.parseVertFile(self.msms_vert_fn)
                self.msp.parseFaceFile(self.msms_face_fn)

                if self.allcpn:  # all componenents of surface
                    for i in range(len(self.msms_cpn_vert_fn_list)):
                        vfn = self.msms_cpn_vert_fn_list[i]
                        ffn = self.msms_cpn_face_fn_list[i]
                        cpn_msp = MSMSSurfPymol()
                        cpn_msp.parseVertFile(vfn)
                        cpn_msp.parseFaceFile(ffn)
                        self.cpn_msp_list.append(cpn_msp)

                if VERBOSE:
                    print('done!')

                if self.cleanup_msms_output.get():  # clean up
                    self.cleanMSMSOutput()

        elif cmd == 'Display Mesh':
            if len(self.msp.vert_coords) == 0:
                err_msg = 'Please execute MSMS first.'
                print('ERROR: %s' % (err_msg,))
                tkMessageBox.showinfo(title='ERROR', message=err_msg)
            else:
                self.msp.displayMsmsSurfMesh(mesh_cgo_color=(self.mesh_col_R / 255.0,
                                              self.mesh_col_G / 255.0,
                                              self.mesh_col_B / 255.0))
                if self.allcpn:
                    for i in range(len(self.cpn_msp_list)):
                        self.cpn_msp_list[i].displayMsmsSurfMesh(
                            mesh_cgo_name='msms_surf_mesh_%d' % (i + 1),
                            mesh_cgo_color=(self.mesh_col_R / 255.0,
                                             self.mesh_col_G / 255.0,
                                             self.mesh_col_B / 255.0))

        elif cmd == 'Display Vertices':
            if len(self.msp.vert_coords) == 0:
                err_msg = 'Please execute MSMS first.'
                print('ERROR: %s' % (err_msg,))
                tkMessageBox.showinfo(title='ERROR', message=err_msg)

            else:
                print('Vertex sphere radius =', self.vert_rad.get())
                print('Normal vector length =', self.norm_len.get())
                print('Vertex color = (%.2f, %.2f, %.2f)' % \
                      (self.vert_col_R / 255.0, self.vert_col_G / 255.0, self.vert_col_B / 255.0))
                print('Normal vector color = (%.2f, %.2f, %.2f)' % \
                      (self.norm_col_R / 255.0, self.norm_col_G / 255.0, self.norm_col_B / 255.0))
                self.msp.displayMsmsSurfVert(r=self.vert_rad.get(),
                                             norm_len=self.norm_len.get(),
                                             vert_cgo_color=(self.vert_col_R / 255.0,
                                                             self.vert_col_G / 255.0,
                                                             self.vert_col_B / 255.0),
                                             norm_cgo_color=(self.norm_col_R / 255.0,
                                                             self.norm_col_G / 255.0,
                                                             self.norm_col_B / 255.0)
                                             )
                if self.allcpn:
                    for i in range(len(self.cpn_msp_list)):
                        self.cpn_msp_list[i].displayMsmsSurfVert(
                            r=self.vert_rad.get(),
                            norm_len=self.norm_len.get(),
                            vert_cgo_color=(self.vert_col_R / 255.0,
                                             self.vert_col_G / 255.0,
                                             self.vert_col_B / 255.0),
                            norm_cgo_color= (self.norm_col_R / 255.0,
                                             self.norm_col_G / 255.0,
                                             self.norm_col_B / 255.0),
                            vert_cgo_name='msms_surf_vert_%d' % (i + 1),
                            norm_cgo_name='msms_surf_nrom_%d' % (i + 1)
                        )

        elif cmd == 'Exit':
            print('Exiting MSMS Plugin ...')
            if __name__ == '__main__':
                self.parent.destroy()
            else:
                self.dialog.withdraw()
            print('Done.')
        else:
            print('Exiting MSMS Plugin ...')
            self.dialog.withdraw()
            print('Done.')

    def quit(self):
        self.dialog.destroy()


#########################
#
# MSMS related
#
#########################
class Msms:

    def __init__(self, msms_bin,  # pdb2xyzr_bin,
                 pdb2xyzrn_bin,
                 pr=1.4, den=3.0, hden=3.0, noh=False, all_components=False,
                 output_dir='/tmp'):
        """ set MSMS binary and parameters

            @param msms_bin: MSMS binary
            @param type: string

            @param pdb2xyzr_bin: script pdb_to_xyzr
            @param type: string

            @param pdb2xyzrn_bin: script pdb_to_xyzrn
            @param type: string

            @param pr: probe parameter (default in MSMS is 1.5)
            @param type: float

            @param den: triangulation density (default in MSMS is 1.0)
              recommended by M. Sanner to set to 1.0 for large
              molecules (>1000 atoms) and 3.0 for small molecules
            @param type: float

            @param hden: triangulation high density (default in MSMS is 3.0)
            @param type: float

            @param noh: Whether hydrogens are considered
            @param type: boolean

            @param all_components: Whether MSMS should find all components.
              If False, only external component is found.
            @param type: boolean

            @param output_dir: where .xyzr and .xyzrn files are stored
            @param type: string
        """
        self.msms_bin = os.path.abspath(os.path.expanduser(msms_bin))
        # msms_wd: MSMS work dir. The file atmtypenumbers is stored there.
        self.msms_wd = os.path.dirname(self.msms_bin)

##         self.pdb2xyzr_bin  = os.path.abspath(os.path.expanduser(pdb2xyzr_bin))
        self.pdb2xyzrn_bin = os.path.abspath(os.path.expanduser(pdb2xyzrn_bin))

        self.param_pr = pr   # probe radius
        self.param_den = den  # surface vertex density
        self.param_hden = hden  # surface vertex density
        self.noh = noh  # whether hydrogens are considered
        self.all_components = all_components

        self.output_dir = os.path.abspath(os.path.expanduser(output_dir))

##         self.output_xyzr_fn  = None
##         self.output_xyzrn_fn = None
        self.output_vert_fn = None
        self.output_face_fn = None
        self.output_cpn_vert_fn = []  # component .vert file name (if all_compoents == True)
        self.output_cpn_face_fn = []  # component .face file name (if all_compoents == True)

        return

    def setMsmsBin(self, msms_bin):
        """ set the executable of MSMS program.
        """
        self.msms_bin = os.path.abspath(os.path.expanduser(msms_bin))
        self.msms_wd = os.path.dirname(self.msms_bin)
        return

# def setPdb2xyzrBin(self, pdb2xyzr_bin):
##         self.pdb2xyzr_bin = os.path.abspath(os.path.expanduser(pdb2xyzr_bin))
# return

    def setPdb2xyzrnBin(self, pdb2xyzrn_bin):
        self.pdb2xyzrn_bin = os.path.abspath(os.path.expanduser(pdb2xyzrn_bin))
        return

    def setMsmsParam(self, pr=1.4, den=3.0, hden=3.0, noh=False, all_components=False):
        """ set parameters *probe radius* used by the MSMS program

            @param pr: probe parameter
              (default in MSMS is 1.5, here we use 1.4 by default)
            @param type: float

            @param den: triangulation density
              (default in MSMS is 1.0, here we use 3.0 by default)
              recommended by M. Sanner to set to 1.0 for large
              molecules (>1000 atoms) and 3.0 for small molecules
            @param type: float

        """
        self.param_pr = pr
        self.param_den = den
        self.param_hden = hden
        self.noh = noh
        self.all_components = all_components
        return

    def setOutputDir(self, output_dir):
        self.output_dir = os.path.abspath(os.path.expanduser(output_dir))
        return

    def printMsmsParam(self):
        print('Probe radius               = %.2f' % (self.param_pr,))
        print('Triangulation density      = %.2f' % (self.param_den,))
        print('Triangulation high density = %.2f' % (self.param_hden,))
        print('No hydrogen    = %s' % (self.noh,))
        print('All components = %s' % (self.all_components,))
        return

    def run(self, pdb_fn, ofn_root=None):
        """ Run MSMS on given pdb file. Output file names are stored.

            @param ofn_root: output file name root.
              The output of MSMS will be ofn_root.vert and ofn_root.face.
            @param type: string

        """
        # convert .pdb file to .xyzr and .xyzrn files
        #fname_root = '.'.join(self.pdb_fn.split('/')[-1].split('.')[:-1])
        fname_root = os.path.splitext(os.path.split(pdb_fn)[-1])[0]
##         xyzr_fname  = '%s/%s.xyzr' % (self.output_dir, fname_root)
        #xyzrn_fname = '%s/%s.xyzrn' % (self.output_dir, fname_root)
        xyzrn_fname = os.path.join(self.output_dir, "%s.xyzrn" % fname_root)

        if ofn_root is None:
            ofn_root = '%s_surface' % (fname_root,)

        old_cwd = os.getcwd()
        os.chdir(self.msms_wd)
##         cmd = '%s %s > %s' % (self.pdb2xyzr_bin, pdb_fn, xyzr_fname)
# os.system(cmd)
        cmd = '%s %s > %s' % (self.pdb2xyzrn_bin, pdb_fn, xyzrn_fname)
        # os.system(cmd)
        print(cmd)
        if sys.platform.startswith('win') and 'PYMOL_GIT_MOD' in os.environ:
            pymol_env = os.environ.copy()
            callfunc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=pymol_env)
            child_stdout, child_stderr = callfunc.communicate()
            print(child_stdout)
            print(child_stderr)
            retval = callfunc.returncode
            print("pdb2xyzrn's mainCommand returned", retval)
        else:
            status = subprocess.call(cmd, shell=True)
        # read in .xyzr and .xyzrn data
# try:
##             xyzr_fh = open(xyzr_fname)
##             self.xyzr_fd = xyzr_fh.readlines()
# xyzr_fh.close()
# except IOError:
# print 'ERROR: pdb2xyzr failed to convert pdb file to xyzr file!'
# print '       pdb2xyzr = %s' % (self.pdb2xyzr_bin,)
# print '       pdb file = %s' % (pdb_fn,)
# sys.exit()

        try:
            xyzrn_fh = open(xyzrn_fname)
            self.xyzrn_fd = xyzrn_fh.readlines()
            xyzrn_fh.close()
        except IOError:
            print('ERROR: pdb2xyzrn failed to convert pdb file to xyzrn file!')
            print('       pdb2xyzrn = %s' % (self.pdb2xyzrn_bin,))
            print('       pdb file  = %s' % (pdb_fn,))
            sys.exit()

        #output_root = '%s/%s' % (self.output_dir, ofn_root)
        output_root = os.path.join(self.output_dir, ofn_root)

        # run MSMS on .xyzrn file
        msms_bin_str = '\"%s\"' % (self.msms_bin,)  # there may be whitespace in path
        cmd = '%s -if %s -probe_radius %f -density %f -hdensity %f -no_area -of %s' % \
              (msms_bin_str, xyzrn_fname,
               self.param_pr, self.param_den, self.param_hden,
               output_root)

        if self.noh:               # ignore hydrogen atoms
            cmd += ' -noh'

        if self.all_components:    # force MSMS to search all surface components
            cmd += ' -all_components'

        if VERBOSE:
            print('command line for running msms:')
            print(cmd)

        # os.system(cmd)
        status = subprocess.call(cmd, shell=True)
        os.chdir(old_cwd)

##         self.output_xyzr_fn  = xyzr_fname
##         self.output_xyzrn_fn = xyzrn_fname

        # clean up intermediate files
# if os.path.isfile(xyzr_fname):
# os.remove(xyzr_fname)
        #!if os.path.isfile(xyzrn_fname):
        #!    os.remove(xyzrn_fname)

        self.output_vert_fn = '%s.vert' % (output_root,)
        self.output_face_fn = '%s.face' % (output_root,)
        if self.all_components:
            fn_idx = 1
            component_vert_fn = '%s_%d.vert' % (output_root, fn_idx)
            component_face_fn = '%s_%d.face' % (output_root, fn_idx)

            while os.path.isfile(component_vert_fn) and \
                    os.path.isfile(component_face_fn):
                self.output_cpn_vert_fn.append(component_vert_fn)
                self.output_cpn_face_fn.append(component_face_fn)

                fn_idx += 1
                component_vert_fn = '%s_%d.vert' % (output_root, fn_idx)
                component_face_fn = '%s_%d.face' % (output_root, fn_idx)

        return

    def getOutputFiles(self):
        return [self.output_vert_fn, self.output_face_fn,
                self.output_cpn_vert_fn, self.output_cpn_face_fn]


class MsmsSurfaceVertex:

    def __init__(self, vid, coord, norm, fnum, snum, vtype, lip=0.0):
        """ @param lip: lipophilicity value of the vertex
            @param type: float
        """
        self.vid = vid
        self.coord = coord
        self.norm = norm
        self.fnum = fnum
        self.snum = snum  # sphere number, or the number of the atom
        self.vtype = vtype
        self.lip = lip
        return


class MsmsSurfaceFace:

    def __init__(self, fid, tri_vid, ftype, fnum):
        self.fid = fid
        self.tri_vid = tri_vid  # a tuple of 3 vertices of the face triangle
        self.ftype = ftype
        self.fnum = fnum
        return


class MsmsSurfacePatch:

    def __init__(self, pid, center, r, atom, vert, face, sid_map, vid_map):
        self.pid = pid
        self.center = center
        self.r = r
        self.atom = atom
        self.vert = vert
        self.face = face
        self.sid_map = sid_map  # each face is related to a sphere, or an atom
        # the atom is specified by an id in the full list of atoms
        # after split, atoms ids are different
        # this is the map from original id to new id in the patch
        self.vid_map = vid_map  # the face is defined by vertex ids
        # the vertex id are 1-based in the total list of id
        # since they are split into patches, the ids are different
        # this is the map from original id to new id in the patch
        return


class MsmsOutputParser:

    """ This parser reads the output .vert and .face files generated by Msms,
        and return MsmsSurface instance.
    """

    def __init__(self):
        return

    def parseMsmsSurface(self, pdb_fn, vert_fn, face_fn, rm_dup=False):
        """ Parse the main component, i.e., the external surface.
            the original pdb file is required for mapping the surface
            vertices to atoms in the PDB.

            @param rm_dup: whether to remove vertices of same
              coordinates (but different normal vectors)
            @param type: boolean
        """
        surf = MsmsSurface()

        ##############
        # read ATOM and HETATM entries in the PDB file
        ##############
        pdb_fh = open(pdb_fn)  # let the exception raise
        buf = pdb_fh.readlines()
        pdb_fh.close()
        surf.pdb_fn = pdb_fn

        # the ATOM and HETATM lines are kept
        surf.pdb_atom = [line for line in buf \
                         if (line.startswith('ATOM  ') or line.startswith('HETATM')) \
                         and len(line) > 52]

        ##############
        # read .vert and .face files
        ##############
        fh = open(vert_fn)
        surf.surf_vert_fd = fh.readlines()
        fh.close()
        surf.surf_vert_fn = vert_fn

        fh = open(face_fn)
        surf.surf_face_fd = fh.readlines()
        fh.close()
        surf.surf_face_fn = face_fn

        # NOTE: MSMS sometimes generates duplicate vertices, i.e.
        #       vertices on surface with exactly the same coordinates
        #       This might result from the limited precision in .vert files.
        #       Nevertheless, i remove duplicate coordinates!
        #       Meanwhile, I also update edges in the patch to reflect
        #       the change of vertices.

        ##############
        # parse vertices from .vert file data
        ##############
        buf_v = surf.surf_vert_fd[2].split()
        surf.surf_vert_num = int(buf_v[0])
        surf.total_sphere_num = int(buf_v[1])
        surf.surf_tri_den = float(buf_v[2])
        surf.probe_radius = float(buf_v[3])

        vid = 0         # NOTE: vid is 1-based in MSMS
        dup_vdict = {}  # correspondance between old vid and new vid
        coord_dict = {}  # dictionary for coordinates
        new_vid_dict = {}  # vid in faces should be updated as some vertices are removed

        if not rm_dup:

            for line in surf.surf_vert_fd[3:]:
                vid += 1
                ck = line[0:29]
                coord = (float(line[0:9]), float(line[10:19]), float(line[20:29]))
                norm = (float(line[30:39]), float(line[40:49]), float(line[50:59]))
                fnum = int(line[60:67])
                snum = int(line[68:75])
                vtype = int(line[76:78])
                vert = MsmsSurfaceVertex(vid, coord, norm, fnum, snum, vtype)
                surf.surf_vert.append(vert)

            surf.surf_vert_num -= len(dup_vdict)  # update number of vertex

        else:  # dup should be removed

            for line in surf.surf_vert_fd[3:]:

                vid += 1
                ck = line[0:29]

                if ck in coord_dict:  # coord has been reported before

                    print("this vertex found before:", vid, 'dup to', coord_dict[ck])
                    print(ck)
                    dup_vdict[vid] = coord_dict[ck]  # point this v to its dup
                    print(dup_vdict)
                    new_vid_dict[vid] = coord_dict[ck]
                    continue

                else:  # only if the coord of the vertex is new

                    coord_dict[ck] = vid - len(dup_vdict)
                    new_vid_dict[vid] = vid - len(dup_vdict)
                    coord = (float(line[0:9]), float(line[10:19]), float(line[20:29]))
                    norm = (float(line[30:39]), float(line[40:49]), float(line[50:59]))
                    fnum = int(line[60:67])
                    snum = int(line[68:75])
                    vtype = int(line[76:78])
                    # note: vid should be updated
                    vert = MsmsSurfaceVertex(vid - len(dup_vdict),
                                              coord, norm, fnum, snum, vtype)
                    surf.surf_vert.append(vert)

            surf.surf_vert_num -= len(dup_vdict)  # update number of vertex
            print("INFO: Remove %d surface vertices due to duplicate coordinates." % (len(dup_vdict)))

        # In the MSMS output .vert file, each vertice is related to a sphere,
        #    which is an atom. The surface atoms can be extracted.
        # get close sphere id (1-based)
        sid_list = [int(line[68:75]) for line in surf.surf_vert_fd[3:]]
        sid_list = set(sid_list)
        for sid in sid_list:
            surf.surf_atom.append(surf.pdb_atom[sid - 1])

        ##############
        # parse faces from .face file data
        ##############
        buf_f = surf.surf_face_fd[2].split()
        surf.surf_face_num = int(buf_f[0])
        assert(surf.total_sphere_num == int(buf_f[1]))
        assert(surf.surf_tri_den == float(buf_f[2]))
        assert(surf.probe_radius == float(buf_f[3]))
        fid = 0  # NOTE: fid is 1-based
        ignored_face_num = 0  # how many faces have been ignored because of removd v
        face_dict = {}  # since faces are re-organized for the removed vertices
        # I shall check no duplicate faces are parsed

        if not rm_dup:
            for line in surf.surf_face_fd[3:]:
                fid += 1
                vid1, vid2, vid3 = int(line[0:6]), int(line[7:13]), int(line[14:20])
                vid_tri = [vid1, vid2, vid3]
                vid_tri.sort()
                vid_tri_key = '%d_%d_%d' % (vid_tri[0], vid_tri[1], vid_tri[2])
                face_dict[vid_tri_key] = 1
                ftype = int(line[21:23])
                fnum = int(line[24:30])
                face = MsmsSurfaceFace(fid, (vid1, vid2, vid3), ftype, fnum)
                surf.surf_face.append(face)

        else:  # dup vertices are removed, faces need to be updated

            for line in surf.surf_face_fd[3:]:
                fid += 1
                vid1, vid2, vid3 = int(line[0:6]), int(line[7:13]), int(line[14:20])
                vid1 = new_vid_dict[vid1]  # update vid
                vid2 = new_vid_dict[vid2]
                vid3 = new_vid_dict[vid3]
                if vid1 == vid2 or vid1 == vid3 or vid2 == vid3:
                    ignored_face_num += 1
                    continue

                vid_tri = [vid1, vid2, vid3]
                vid_tri.sort()
                vid_tri_key = '%d_%d_%d' % (vid_tri[0], vid_tri[1], vid_tri[2])

                if vid_tri_key in face_dict:
                    ignored_face_num += 1
                    continue
                else:
                    face_dict[vid_tri_key] = 1
                    ftype = int(line[21:23])
                    fnum = int(line[24:30])
                    face = MsmsSurfaceFace(fid - ignored_face_num,
                                            (vid1, vid2, vid3), ftype, fnum)
                    surf.surf_face.append(face)

            print("INFO: Remove %d surface faces due to duplicate coordinates." % (ignored_face_num))

        return surf

    def parseMsmsSurfaceAllComponents(self, pdb_fn,
                                      vert_fn, face_fn,
                                      cpn_vert_fns, cpn_face_fns):
        """ Parse all components, including both external and internal surfaces.
        """
        # parse the external surface
        extl_surf = self.parseMsmsSurface(pdb_fn, vert_fn, face_fn)
        # parse internal surfaces
        itnl_surf = [self.parseMsmsSurface(pdb_fn, vf, ff) for vf, ff in zip(cpn_vert_fns, cpn_face_fns)]
        return extl_surf, itnl_surf


class MsmsSurface:

    def __init__(self):  # , msms, pdb_fn, all_components=False):
        """
            @param msms: msms instance
            @param type: an instance of Msms class
        """
        self.pdb_fn = None
        self.pdb_atom = []
        self.pdb_atom_lip = []  # lipophilicity of atoms in the PDB

##         self.xyzr_fd  = []
        self.xyzrn_fd = []

        self.surf_vert_fn = None  # .vert file name
        self.surf_face_fn = None  # .face file name
        self.surf_vert_fd = []   # .vert file data
        self.surf_face_fd = []   # .face file data

        self.surf_vert_lip = []  # lipophilicity values for each surface vertex

        self.surf_atom = []  # surface atoms
        self.surf_vert = []  # surface vertices
        self.surf_face = []  # surface faces

        self.total_sphere_num = 0  # total number of shperes in the pdb
        self.surf_vert_num = 0  # number of vertices on the surface
        self.surf_face_num = 0  # number of faces
        self.surf_sphere_num = 0  # number of surface spheres
        self.surf_tri_den = 0  # triangulation density
        self.surf_tri_hden = 0  # triangulation high density
        self.probe_radius = 0  # probe radius

        self.surf_patch = []
        self.surf_patch_num = 0
        return

    def saveSurfaceAsObj(self, obj_fn):
        """ Read output files from MSMS and organize them into a .obj file.
            This is need by parameterization.

        """
        if not os.path.isfile(self.surf_vert_fn) or \
                not os.path.isfile(self.surf_face_fn):
            # TODO: exception should be thrown out here
            return False

        try:
            fh = open(self.surf_vert_fn)
            vert_fdata = fh.readlines()
            fh.close()
        except IOError:
            print('File not found:', self.surf_vert_fn)

        try:
            fh = open(self.surf_face_fn)
            face_fdata = fh.readlines()
            fh.close()
        except IOError:
            print('File not found:', self.surf_face_fn)

        # the first three lines are not needed
        vert_fdata = vert_fdata[3:]
        # each line contains one vertex, the first three fields
        # are the coord of vertices: vert1,vert2,vert3
        vertice_lines = ['v %.6f %.6f %.6f' % \
                         (float(line.split()[0]),
                          float(line.split()[1]),
                          float(line.split()[2])) \
                         for line in vert_fdata]

        # the first three lines are not needed
        face_fdata = face_fdata[3:]
        # Face file contains triangles composed of vertices (1-based index):
        # vert_index1,vert_index2,vert_index3,toric_flag,face_num(of analytical surf)
        face_lines = ['f %d/%d %d/%d %d/%d' % \
                      (int(line.split()[0]), int(line.split()[0]),
                       int(line.split()[1]), int(line.split()[1]),
                       int(line.split()[2]), int(line.split()[2])) \
                      for line in face_fdata]

        obj_buf = vertice_lines
        obj_buf.extend(face_lines)

        fh = open(obj_fn, 'w')
        fh.writelines(obj_buf)
        fh.close()

        return True


class MSMSSurfPymol:

    """ This class parse MSMS output files (.vert, .face)
        and convert them into PyMOL objects
    """

    def __init__(self, vert_fn=None, face_fn=None):

        self.vert_fn = vert_fn
        self.face_fn = face_fn

        self.vert_coords = []
        self.vert_norms = []
        self.vert_indicator = []  # 0 for v on surface, nega for v on edge
        self.vert_sidx = []  # closest sphere index
        self.vert_feature = []  # v feature

        self.face_tri = []
        self.face_indicator = []
        self.face_ana_num = []

    def parseVertFile(self, vert_fn=None):
        """ Parse MSMS .vert file. Read
                http://www.scripps.edu/~sanner/html/msms_man.html
            for the specification of the file format.
        """
        if vert_fn is None:
            vert_fn = self.vert_fn

        try:
            fh = open(vert_fn, 'r')
            fd = fh.readlines()
            fh.close()
        except IOError:
            print('Error: MSMS .vert file not found:', self.vert_fn)
            return

        self.vert_coords = []
        self.vert_norms = []
        self.vert_indicator = []  # 0 for v on surface, nega for v on edge
        self.vert_sidx = []  # closest sphere index
        self.vert_feature = []  # v feature

        if fd[0].startswith('#'):  # remove info lines
            fd = fd[3:]

        for line in fd:
            self.vert_coords.append([float(line[:9]),
                                     float(line[10:19]),
                                     float(line[20:29])])
            self.vert_norms.append([float(line[30:39]),
                                    float(line[40:49]),
                                    float(line[50:59])])
            self.vert_indicator.append(int(line[60:67]))
            self.vert_sidx.append(int(line[68:75]))
            self.vert_feature.append(int(line[76:78]))

        print('Number of vertices =', len(self.vert_coords))
        print('Number of normal vectors =', len(self.vert_norms))

        return

    def parseFaceFile(self, face_fn=None):

        if face_fn is None:
            face_fn = self.face_fn

        try:
            fh = open(face_fn, 'r')
            fd = fh.readlines()
            fh.close()
        except IOError:
            print('Error: MSMS .face file not found:', self.face_fn)
            return

        self.face_tri = []
        self.face_indicator = []
        self.face_ana_num = []

        if fd[0].startswith('#'):  # remove info lines
            fd = fd[3:]

        for line in fd:
            self.face_tri.append([int(line[0:6]), int(line[7:13]), int(line[14:20])])
            self.face_indicator.append(int(line[21:23]))
            self.face_ana_num.append(int(line[24:30]))

        print('Number of surface triangles =', len(self.face_tri))

        return

    def displayMsmsSurfVert(self,
                            r=0.5, display_norm=True,
                            norm_len=1.0,
                            vert_cgo_color=(1.0, 1.0, 1.0),
                            norm_cgo_color = (1.0, 0.7, 0.2),
                            vert_cgo_name='msms_surf_vert',
                            norm_cgo_name='msms_surf_norm'):
        """
            @param r: the radius of spheres representing vertices
            @param type: float

            @param norm_len: normal vector length
            @param type: float
        """
        vert_cgo = [COLOR,
                    vert_cgo_color[0], vert_cgo_color[1], vert_cgo_color[2]]
        [vert_cgo.extend([SPHERE, c[0], c[1], c[2], r]) for c in self.vert_coords]
        # vert_cgo.append(END)
        cmd.load_cgo(vert_cgo, vert_cgo_name)
        print('Number of spheres for vertices =', len(self.vert_coords))

        if display_norm and len(self.vert_norms) > 0:
            norm_cgo = [BEGIN, LINES,
                        COLOR,
                        norm_cgo_color[0], norm_cgo_color[1], norm_cgo_color[2]]
            [norm_cgo.extend([VERTEX, c[0], c[1], c[2], VERTEX,
                              c[0] + n[0] * norm_len,
                              c[1] + n[1] * norm_len,
                              c[2] + n[2] * norm_len]) \
             for c, n, in zip(self.vert_coords, self.vert_norms)]
            norm_cgo.append(END)
            cmd.load_cgo(norm_cgo, norm_cgo_name)
            print('Number of lines for normal vectors =', len(self.vert_norms))

        return

    def displayMsmsSurfMesh(self,
                            mesh_cgo_name='msms_surf_mesh',
                            mesh_cgo_color=(1.0, 1.0, 1.0)):
        """ Render MSMS surface as mesh (triangles).
        """
        line_dict = {}
        mesh_cgo = [BEGIN, LINES,
                    COLOR,
                    mesh_cgo_color[0], mesh_cgo_color[1], mesh_cgo_color[2]]

        line_num = 0
        for t in self.face_tri:

            if t[0] < t[1]:
                k = '%d_%d' % (t[0], t[1])
            else:
                k = '%d_%d' % (t[1], t[0])

            # NOTE: vertices indices are 1-based
            if k not in line_dict:
                line_dict[k] = 1
                mesh_cgo.extend([VERTEX,
                                 self.vert_coords[t[0] - 1][0],
                                 self.vert_coords[t[0] - 1][1],
                                 self.vert_coords[t[0] - 1][2],
                                 VERTEX,
                                 self.vert_coords[t[1] - 1][0],
                                 self.vert_coords[t[1] - 1][1],
                                 self.vert_coords[t[1] - 1][2]])
                line_num += 1

            if t[1] < t[2]:
                k = '%d_%d' % (t[1], t[2])
            else:
                k = '%d_%d' % (t[2], t[1])
            if k not in line_dict:
                line_dict[k] = 1
                mesh_cgo.extend([VERTEX,
                                 self.vert_coords[t[1] - 1][0],
                                 self.vert_coords[t[1] - 1][1],
                                 self.vert_coords[t[1] - 1][2],
                                 VERTEX,
                                 self.vert_coords[t[2] - 1][0],
                                 self.vert_coords[t[2] - 1][1],
                                 self.vert_coords[t[2] - 1][2]])
                line_num += 1

            if t[0] < t[2]:
                k = '%d_%d' % (t[0], t[2])
            else:
                k = '%d_%d' % (t[2], t[0])
            if k not in line_dict:
                line_dict[k] = 1
                mesh_cgo.extend([VERTEX,
                                 self.vert_coords[t[0] - 1][0],
                                 self.vert_coords[t[0] - 1][1],
                                 self.vert_coords[t[0] - 1][2],
                                 VERTEX,
                                 self.vert_coords[t[2] - 1][0],
                                 self.vert_coords[t[2] - 1][1],
                                 self.vert_coords[t[2] - 1][2]])
                line_num += 1

        mesh_cgo.append(END)
        cmd.load_cgo(mesh_cgo, mesh_cgo_name)
        print('Number of lines in the mesh =', line_num)

        return


class ObjPymol:

    """ Parse .obj files and convert to PyMOL cgo objects
    """

    def __init__(self, obj_fn=None):

        self.obj_fn = obj_fn

        self.vert_coords = []
        self.vert_norms = []
        self.vert_texture = []

        self.face_tri = []
        self.face_vtexture = []

    def parseObjFile(self, obj_fn=None):

        if obj_fn is None:
            obj_fn = self.obj_fn

        try:
            fh = open(obj_fn, 'r')
            fd = fh.readlines()
            fh.close()
        except IOError:
            print('Error: MSMS .vert file not found:', self.obj_fn)
            return

        for line in fd:

            if line.startswith('v '):      # vertices
                buf = line.strip().split()
                self.vert_coords.append((float(buf[1]), float(buf[2]), float(buf[3])))

            elif line.startswith('vn '):   # normal vectors
                buf = line.strip().split()
                self.vert_norms.append((float(buf[1]), float(buf[2]), float(buf[3])))

            elif line.startswith('f '):    # faces (in general, faces are not
                                           # necessarily triangles in .obj file)
                buf = line.strip().split()
                vs = [int(v.split('/')[0]) for v in buf[1:]]
                self.face_tri.append(vs)
                if len(buf[1].split('/')) > 1:
                    vt = [int(v.split('/')[1]) for v in buf[1:]]  # vertex texture
                    self.face_vtexture.append(vt)

        print('Number of vertices =', len(self.vert_coords))
        print('Number of normal vectors =', len(self.vert_norms))

        return

    def displayVert(self,
                    r=0.5, display_norm=True,
                    norm_len=1.0,
                    vert_cgo_color=(1.0, 1.0, 1.0),
                    norm_cgo_color = (1.0, 0.7, 0.2),
                    vert_cgo_name='patch_vert',
                    norm_cgo_name='patch_vert_norm'):
        """
        """
        vert_cgo = [COLOR,
                    vert_cgo_color[0], vert_cgo_color[1], vert_cgo_color[2]]
        [vert_cgo.extend([SPHERE, c[0], c[1], c[2], r]) for c in self.vert_coords]
        cmd.load_cgo(vert_cgo, vert_cgo_name)
        print('Number of spheres for vertices =', len(self.vert_coords))

        if display_norm and len(self.vert_norms) > 0:
            norm_cgo = [BEGIN, LINES,
                        COLOR,
                        norm_cgo_color[0], norm_cgo_color[1], norm_cgo_color[2]]
            [norm_cgo.extend([VERTEX, c[0], c[1], c[2], VERTEX,
                              c[0] + n[0] * norm_len,
                              c[1] + n[1] * norm_len,
                              c[2] + n[2] * norm_len]) \
             for c, n, in zip(self.vert_coords, self.vert_norms)]
            norm_cgo.append(END)
            cmd.load_cgo(norm_cgo, norm_cgo_name)
            print('Number of lines for normal vectors =', len(self.vert_norms))

        return


#############################################
#
#
# Create demo in root window for testing.
#
#
##############################################
if __name__ == '__main__':

    class App:

        def my_show(self, *args, **kwargs):
            pass

    app = App()
    app.root = Tkinter.Tk()
    Pmw.initialise(app.root)
    app.root.title('It seems to work!')

    widget = MSMSPlugin(app)
    app.root.mainloop()
