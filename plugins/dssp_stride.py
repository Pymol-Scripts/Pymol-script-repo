""" 2011_04_11: Hongbo Zhu
      PyMOL plugin for running DSSP and display SS using different colors.

      The DSSP codes for secondary structure are: 
   
      H  -  Alpha helix (4-12) 
      G  -  3-10 helix 
      I  -  pi helix 

      E  -  Strand 
      B  -  Isolated beta-bridge residue 

      T  -  Turn 
      S  -  Bend 
      -  -  None

      2011_04_24: add stride
      STRIDE code (taken from stride.doc) is nearly the same as DSSP
      H         Alpha helix
      G         3-10 helix
      I         PI-helix
      E         Extended conformation
      B or b    Isolated bridge
      T         Turn
      C         Coil (none of the above)

"""

# Copyright Notice
# ================
#
# The PyMOL Plugin source code in this file is copyrighted, but you can
# freely use and copy it as long as you don't change or remove any of
# the copyright notices.
#
# ----------------------------------------------------------------------
#               This PyMOL Plugin is Copyright (C) 2011 by
#            Hongbo Zhu <hongbo.zhu.cn at googlemail dot com>
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

from __future__ import print_function

# python lib
import os
import sys
import platform
if sys.version_info >= (2, 4):
    import subprocess  # subprocess is introduced in python 2.4
import math
import random
import tempfile


if sys.version_info[0] > 2:
    import tkinter.simpledialog as tkSimpleDialog
    import tkinter.messagebox as tkMessageBox
    import tkinter.filedialog as tkFileDialog
    import tkinter.colorchooser as tkColorChooser
    import tkinter as Tkinter
    from shutil import which
else:
    import tkSimpleDialog
    import tkMessageBox
    import tkFileDialog
    import tkColorChooser
    import Tkinter
    which = lambda cmd: None

# pymol lib
try:
    from pymol import cmd
    from pymol.cgo import *
except ImportError:
    print('Warning: pymol library cmd not found.')
    sys.exit(1)

# external lib
import Pmw

VERBOSE = True

#################
# here we go
#################


def __init__(self):
    """ DSSP and Stride plugin for PyMol
    """
    self.menuBar.addmenuitem('Plugin', 'command',
                             'DSSP & Stride', label='DSSP & Stride',
                             command=lambda s=self: DSSPPlugin(s))


#################
# GUI related
#################
class DSSPPlugin:

    def __init__(self, app):

        self.parent = app.root
        self.dialog = Pmw.Dialog(self.parent,
                                 buttons=('Run DSSP',
                                            'Run Stride',
                                            'Update Color',
                                            'Update ss',
                                            'Exit'),
                                 title = 'DSSP and Stride Plugin for PyMOL',
                                 command = self.execute)
        Pmw.setbusycursorattributes(self.dialog.component('hull'))

        # parameters used by DSSP
        self.pymol_sel = Tkinter.StringVar(value="all")
        self.dssp_bin = Tkinter.StringVar()
        self.stride_bin = Tkinter.StringVar()
        self.dssp_rlt_dict = {}
        self.stride_rlt_dict = {}
        self.ss_asgn_prog = None  # which program is used to assign ss

        self.sel_obj_list = []  # there may be more than one seletion or object
        # defined by self.pymol_sel
        # treat each selection and object separately
        if 'DSSP_BIN' not in os.environ and 'PYMOL_GIT_MOD' in os.environ:
            if sys.platform.startswith('linux') and platform.machine() == 'x86_32':
                initialdir_dssp = os.path.join(os.environ['PYMOL_GIT_MOD'], "DSSP", "i86Linux2", "dssp-2")
                os.environ['DSSP_BIN'] = initialdir_dssp
            elif sys.platform.startswith('linux') and platform.machine() == 'x86_64':
                initialdir_dssp = os.path.join(os.environ['PYMOL_GIT_MOD'], "DSSP", "ia64Linux2", "dssp-2")
                os.environ['DSSP_BIN'] = initialdir_dssp
            elif sys.platform.startswith('win'):
                initialdir_dssp = os.path.join(os.environ['PYMOL_GIT_MOD'], "DSSP", "win32", "dssp.exe")
                os.environ['DSSP_BIN'] = initialdir_dssp
            else:
                pass
        if 'DSSP_BIN' in os.environ:
            if VERBOSE:
                print('Found DSSP_BIN in environmental variables', os.environ['DSSP_BIN'])
            self.dssp_bin.set(os.environ['DSSP_BIN'])
        else:
            if VERBOSE:
                print('DSSP_BIN not found in environmental variables.')
            self.dssp_bin.set(which('mkdssp') or '')
        if 'STRIDE_BIN' not in os.environ and 'PYMOL_GIT_MOD' in os.environ:
            if sys.platform.startswith('linux') and platform.machine() == 'x86_32':
                initialdir_stride = os.path.join(os.environ['PYMOL_GIT_MOD'], "Stride", "i86Linux2", "stride")
                os.environ['STRIDE_BIN'] = initialdir_stride
            elif sys.platform.startswith('linux') and platform.machine() == 'x86_64':
                initialdir_stride = os.path.join(os.environ['PYMOL_GIT_MOD'], "Stride", "ia64Linux2", "stride")
                os.environ['STRIDE_BIN'] = initialdir_stride
            elif sys.platform.startswith('win'):
                initialdir_stride = os.path.join(os.environ['PYMOL_GIT_MOD'], "Stride", "win32", "stride.exe")
                os.environ['STRIDE_BIN'] = initialdir_stride
            else:
                pass
        if 'STRIDE_BIN' in os.environ:
            if VERBOSE:
                print('Found STRIDE_BIN in environmental variables', os.environ['STRIDE_BIN'])
            self.stride_bin.set(os.environ['STRIDE_BIN'])
        else:
            if VERBOSE:
                print('STRIDE_BIN not found in environmental variables.')
            self.stride_bin.set(which('stride') or '')

        # DSSP visualization color
        # - H        Alpha helix (4-12)
        # - G        3-10 helix
        # - I        pi helix
        #
        # - E        Extended strand
        # - B        Isolated beta-bridge residue
        #
        # - T        Turn
        # - S        Bend
        # - -        None
        # STRIDE does not have S and None, but it has two more
        # codes: 'b' same as 'B' and 'C' for coil
        self.DSSP_SSE_list = ['H', 'G', 'I', 'E', 'B', 'T', 'S', '-']  # for the record of the order
        self.STRIDE_SSE_list = ['H', 'G', 'I', 'E', 'B', 'b', 'T', 'C']  # for the record of the order
        self.SSE_map = {'H': 'H',
                        'G': 'H',
                        'I': 'H',
                        'E': 'S',
                        'B': 'S',
                        'T': 'L',
                        'S': 'L',
                        '-': 'L',

                        'b': 'S',
                        'C': 'L'
                        }
        self.SSE_name = {
            'H': 'Alpha helix',
            'G': '3-10 helix',
            'I': 'Pi helix',

            'E': 'Extended strand',
            'B': 'Isolated beta-bridge',

            'T': 'Turn',
            'S': 'Bend',
            '-': 'None',

            'b': 'Isolated beta-bridge',
            'C': 'Coil'
        }

        self.SSE_col_RGB = {
            'H': (255, 0, 0),
            'G': (255, 0, 128),
            'I': (255, 170, 170),

            'E': (255, 255, 0),
            'B': (153, 119, 85),

            'T': (0, 255, 255),
            'S': (153, 255, 119),
            '-': (179, 179, 179),

            'b': (153, 119, 85),  # same as 'B'
            'C': (0, 0, 255)
        }
        self.SSE_col = {}
        for sse in self.SSE_col_RGB.keys():
            self.SSE_col[sse] = '#%s%s%s' % (hex(self.SSE_col_RGB[sse][0])[2:].zfill(2),
                                             hex(self.SSE_col_RGB[sse][1])[2:].zfill(2),
                                             hex(self.SSE_col_RGB[sse][2])[2:].zfill(2))

        self.SSE_res_dict = {}
        self.SSE_sel_dict = {}

        w = Tkinter.Label(self.dialog.interior(),
                          #                          text = '\nDSSP Plugin for PyMOL\nHongbo Zhu, 2011.\n\nColor proteins according to the secondary structure determined by DSSP.',
                          text='\nDSSP and Stride Plugin for PyMOL\nby Hongbo Zhu, 2011\n',
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
                                       label_text='PyMOL selection/object:',
                                       labelpos='wn',
                                       entry_textvariable=self.pymol_sel
                                       )

        dssp_bin_ent = Pmw.EntryField(group_struc,
                                      label_text='DSSP binary:', labelpos='wn',
                                      entry_textvariable=self.dssp_bin,
                                      entry_width=20)
        dssp_bin_but = Tkinter.Button(group_struc, text='Browse...',
                                      command=self.getDSSPBin)

        stride_bin_ent = Pmw.EntryField(group_struc,
                                        label_text='Stride binary:', labelpos='wn',
                                        entry_textvariable=self.stride_bin,
                                        entry_width=20)
        stride_bin_but = Tkinter.Button(group_struc, text='Browse...',
                                        command=self.getStrideBin)

        # arrange widgets using grid
        pymol_sel_ent.grid(sticky='we', row=0, column=0,
                           columnspan=2, padx=5, pady=5)
        dssp_bin_ent.grid(sticky='we', row=1, column=0, padx=5, pady=1)
        dssp_bin_but.grid(sticky='we', row=1, column=1, padx=5, pady=1)
        stride_bin_ent.grid(sticky='we', row=2, column=0, padx=5, pady=1)
        stride_bin_but.grid(sticky='we', row=2, column=1, padx=5, pady=1)
        group_struc.columnconfigure(0, weight=9)
        group_struc.columnconfigure(1, weight=1)

        ######################
        # Tab : Color Tab
        ######################
        # H = alpha helix
        # G = 3-helix (3/10 helix)
        # I = 5 helix (pi helix)

        # E = extended strand, participates in beta ladder
        # B = residue in isolated beta-bridge

        # T = hydrogen bonded turn
        # S = bend
        #
        #  H        Alpha helix (4-12)
        #  G        3-10 helix
        #  I        pi helix
        #  E        Strand
        #  B        Isolated beta-bridge residue
        #  T        Turn
        #  S        Bend
        #  -        None

        page = self.notebook.add('Color')

        group_sse_color = Tkinter.LabelFrame(page, text='Secondary Structure Element Color')
        group_sse_color.grid(sticky='eswn', row=0, column=0, padx=10, pady=5)

        # colors for DSSP surface
        H_col_lab = Tkinter.Label(group_sse_color, text='Alpha helix (H):')
        self.H_col_but = Tkinter.Button(group_sse_color,
                                        bg=self.SSE_col['H'],
                                        activebackground=self.SSE_col['H'],
                                        command=self.custermizeHColor)
        G_col_lab = Tkinter.Label(group_sse_color, text='3-10 helix (G):')
        self.G_col_but = Tkinter.Button(group_sse_color,
                                        bg=self.SSE_col['G'],
                                        activebackground=self.SSE_col['G'],
                                        command=self.custermizeGColor)
        I_col_lab = Tkinter.Label(group_sse_color, text='PI helix (I):')
        self.I_col_but = Tkinter.Button(group_sse_color,
                                        bg=self.SSE_col['I'],
                                        activebackground=self.SSE_col['G'],
                                        command=self.custermizeIColor)

        E_col_lab = Tkinter.Label(group_sse_color, text='Extended strand (E):')
        self.E_col_but = Tkinter.Button(group_sse_color,
                                        bg=self.SSE_col['E'],
                                        activebackground=self.SSE_col['E'],
                                        command=self.custermizeEColor)
        B_col_lab = Tkinter.Label(group_sse_color, text='Isolated beta-bridge (B):')
        self.B_col_but = Tkinter.Button(group_sse_color,
                                        bg=self.SSE_col['B'],
                                        activebackground=self.SSE_col['B'],
                                        command=self.custermizeBColor)

        T_col_lab = Tkinter.Label(group_sse_color, text='Turn (T):')
        self.T_col_but = Tkinter.Button(group_sse_color,
                                        bg=self.SSE_col['T'],
                                        activebackground=self.SSE_col['T'],
                                        command=self.custermizeTColor)
        S_col_lab = Tkinter.Label(group_sse_color, text='Bend (DSSP S):')
        self.S_col_but = Tkinter.Button(group_sse_color,
                                        bg=self.SSE_col['S'],
                                        activebackground=self.SSE_col['S'],
                                        command=self.custermizeSColor)
        N_col_lab = Tkinter.Label(group_sse_color, text='None (DSSP NA):')
        self.N_col_but = Tkinter.Button(group_sse_color,
                                        bg=self.SSE_col['-'],
                                        activebackground=self.SSE_col['-'],
                                        command=self.custermizeNColor)

        b_col_lab = Tkinter.Label(group_sse_color, text='Isolated beta-bridge (Stride b):')
        self.b_col_but = Tkinter.Button(group_sse_color,
                                        bg=self.SSE_col['b'],
                                        activebackground=self.SSE_col['b'],
                                        command=self.custermizebColor)

        C_col_lab = Tkinter.Label(group_sse_color, text='Coil (Stride C):')
        self.C_col_but = Tkinter.Button(group_sse_color,
                                        bg=self.SSE_col['C'],
                                        activebackground=self.SSE_col['C'],
                                        command=self.custermizeCColor)

        H_col_lab.grid(sticky='e', row=0, column=0, padx=5, pady=3)
        self.H_col_but.grid(sticky='we', row=0, column=1, padx=5, pady=3)
        G_col_lab.grid(sticky='e', row=1, column=0, padx=5, pady=3)
        self.G_col_but.grid(sticky='we', row=1, column=1, padx=5, pady=3)
        I_col_lab.grid(sticky='e', row=2, column=0, padx=5, pady=3)
        self.I_col_but.grid(sticky='we', row=2, column=1, padx=5, pady=3)

        E_col_lab.grid(sticky='e', row=0, column=2, padx=5, pady=3)
        self.E_col_but.grid(sticky='we', row=0, column=3, padx=5, pady=3)
        B_col_lab.grid(sticky='e', row=1, column=2, padx=5, pady=3)
        self.B_col_but.grid(sticky='we', row=1, column=3, padx=5, pady=3)

        T_col_lab.grid(sticky='e', row=0, column=4, padx=5, pady=3)
        self.T_col_but.grid(sticky='we', row=0, column=5, padx=5, pady=3)
        S_col_lab.grid(sticky='e', row=1, column=4, padx=5, pady=3)
        self.S_col_but.grid(sticky='we', row=1, column=5, padx=5, pady=3)
        N_col_lab.grid(sticky='e', row=2, column=4, padx=5, pady=3)
        self.N_col_but.grid(sticky='we', row=2, column=5, padx=5, pady=3)

        b_col_lab.grid(sticky='e', row=2, column=2, padx=5, pady=3)
        self.b_col_but.grid(sticky='we', row=2, column=3, padx=5, pady=3)

        C_col_lab.grid(sticky='e', row=3, column=4, padx=5, pady=3)
        self.C_col_but.grid(sticky='we', row=3, column=5, padx=5, pady=3)
        ######################
        # Tab : About Tab
        ######################
        page = self.notebook.add('About')
        group_about = Tkinter.LabelFrame(page, text='About DSSP and Stride Plugin for PyMOL')
        group_about.grid(sticky='we', row=0, column=0, padx=5, pady=3)
        about_plugin = """ Assign and color secondary structures using DSSP or Stride.
by Hongbo Zhu <hongbo.zhu.cn .at. googlemail.com>
Please cite this plugin if you use it in a publication.
Hongbo Zhu. DSSP and Stride plugin for PyMOL, 2011, BIOTEC, TU Dresden.
"""

        label_about = Tkinter.Label(group_about, text=about_plugin)
        label_about.grid(sticky='we', row=0, column=0, padx=5, pady=10)

        self.notebook.setnaturalsize()

        return

    def getDSSPBin(self):
        dssp_bin_fname = tkFileDialog.askopenfilename(
            title='DSSP Binary', initialdir='',
            filetypes=[('all', '*')], parent=self.parent)
        if dssp_bin_fname:  # if nonempty
            self.dssp_bin.set(dssp_bin_fname)
        return

    def getStrideBin(self):
        stride_bin_fname = tkFileDialog.askopenfilename(
            title='Stride Binary', initialdir='',
            filetypes=[('all', '*')], parent=self.parent)
        if stride_bin_fname:  # if nonempty
            self.stride_bin.set(stride_bin_fname)
        return

    def runDSSPOneObj(self, one_obj_sel):
        """ Run DSSP on only one object.

            @param one_obj_sel: the selection/object involving only one object
            @param type: string
        """
        SSE_res = {
            'H': {}, 'G': {}, 'I': {},
            'E': {}, 'B': {},
            'T': {}, 'S': {}, '-': {}
        }
        SSE_sel = {
            'H': None, 'G': None, 'I': None,
            'E': None, 'B': None,
            'T': None, 'S': None, '-': None
        }

        pdb_fn = None
        pdb_os_fh, pdb_fn = tempfile.mkstemp(suffix='.pdb')  # file os handle, file name
        os.close(pdb_os_fh)
        # DSSP 2.0.4 ignores all residues after 1st TER in the same chain
        v = cmd.get(name='pdb_use_ter_records')
        if v:
            cmd.set(name='pdb_use_ter_records', value=0)  # do not insert TER into the pdb
        cmd.save(filename=pdb_fn, selection=one_obj_sel)
        if v:
            cmd.set(name='pdb_use_ter_records', value=v)  # restore old value

        if VERBOSE:
            print('Selection %s saved to %s.' % (one_obj_sel, pdb_fn))

        if pdb_fn is None:
            print('WARNING: DSSP has no pdb file to work on!')
            return None

        print('Running DSSP for %s ...' % (one_obj_sel,))
        dssp_sse_dict = {}
        if sys.version_info >= (2, 4):
            dssp_proc = subprocess.Popen([self.dssp_bin.get(), pdb_fn],
                                         stdout=subprocess.PIPE,
                                         stderr=subprocess.STDOUT)
            dssp_stdout, dssp_stderr = dssp_proc.communicate()
        else:  # use os.system + tempfile
            dssp_tmpout_os_fh, dssp_tmpout_fn = tempfile.mkstemp(suffix='.dssp')
            os.close(dssp_tmpout_os_fh)
            dssp_cmd = '%s %s > %s' % (self.dssp_bin.get(), pdb_fn, dssp_tmpout_fn)
            os.system(dssp_cmd)
            fh = open(dssp_tmpout_fn)
            dssp_stdout = ''.join(fh.readlines())
            fh.close()

        sse_started = False
        for line in dssp_stdout.splitlines():
            if not isinstance(line, str):  # Python 3
                line = line.decode('ascii', 'ignore')
            if line.startswith('  #  RESIDUE'):
                sse_started = True
                continue
            elif line.startswith(' !!!'):
                sse_started = False
                continue
            elif sse_started:
                if len(line) < 10 or line[9] == ' ':
                    continue
                ch, resname = line[11], line[13]
                residen, sscode = line[5:11].strip(), line[16]  # residen = resnum+icode, col 10 is for icode
                if sscode == ' ':
                    sscode = '-'
                k = (ch, resname, residen)
                dssp_sse_dict[k] = sscode

        self.dssp_rlt_dict[one_obj_sel] = dssp_sse_dict
        print('Got SSE for %d residues.' % (len(self.dssp_rlt_dict[one_obj_sel]),))

        # group residues according to their SSE, and chain name
        for k in self.dssp_rlt_dict[one_obj_sel].keys():
            # res = '/%s//%s/%d%s/' % (sel,k[0],k[1][1], k[1][2].strip()) # sel name, chain ID, res serial num, icode
            sse = self.dssp_rlt_dict[one_obj_sel][k]
            chn = k[0]
            res = k[2]
            if VERBOSE:
                print('(%s) and \"%s\"/%s/  sse=%s' % (str(one_obj_sel),
                                                               chn.strip(), res, sse))
            SSE_res[sse].setdefault(chn, []).append(res)

        self.SSE_res_dict[one_obj_sel] = SSE_res

        for sse in self.DSSP_SSE_list:
            sse_sel_name = self.selectSSE(one_obj_sel, sse)
            SSE_sel[sse] = sse_sel_name

        self.SSE_sel_dict[one_obj_sel] = SSE_sel

        print('\nNumber of residues with SSE element:')
        for sse in self.DSSP_SSE_list:
            num = sum([len(SSE_res[sse][chn]) for chn in SSE_res[sse]])
            print('%20s (%s) : %5d' % (self.SSE_name[sse], sse, num))
        print('')

        # clean up pdb_fn and dssp_tmpout_fn created by tempfile.mkstemp()
        if os.path.isfile(pdb_fn):
            os.remove(pdb_fn)
        if sys.version_info < (2, 4) and os.path.isfile(dssp_tmpout_fn):
            os.remove(dssp_tmpout_fn)

        return

    def runDSSP(self):
        """
            @return: whether DSSP has been executed successfully
            @rtype: boolean
        """
        # delete old results
        self.sel_obj_list = []
        self.dssp_rlt_dict = {}
        self.SSE_res_dict = {}
        self.SSE_sel_dict = {}

        pdb_fn = None
        sel_name = None
        sel = self.pymol_sel.get()

        if len(sel) > 0:  # if any pymol selection/object is specified
            # save the pymol selection/object in the tmp dir
            all_sel_names = cmd.get_names('all')  # get names of all selections
            if sel in all_sel_names:
                if cmd.count_atoms(sel) == 0:
                    err_msg = 'ERROR: The selection %s is empty.' % (sel,)
                    print('ERROR: %s' % (err_msg,))
                    tkMessageBox.showinfo(title='ERROR', message=err_msg)
                    return False
                else:
                    sel_name = sel
            # no selection/object with the input name is found
            # we assume either a single-word selector or
            # some other selection-expression is used
            # NOTE: if more than one selection is selected, they should be
            #       saved separately and SSE should be calculated for each
            #       of the selections.
            else:
                print('The selection/object you specified is not found.')
                print('Your input will be interpreted as a selection-expression.')
                # check whether the selection is empty
                tmpsel = self.randomSeleName(prefix='your_sele_')
                cmd.select(tmpsel, sel)
                if cmd.count_atoms(tmpsel) == 0:
                    cmd.delete(tmpsel)
                    err_msg = 'ERROR: The selection %s is empty.' % (sel,)
                    print('ERROR: %s' % (err_msg,))
                    tkMessageBox.showinfo(title='ERROR', message=err_msg)
                    return False
                else:
                    sel_name = tmpsel
        else:   # what structure do you want DSSP to work on?
            err_msg = 'No PyMOL selection/object specified!'
            print('ERROR: %s' % (err_msg,))
            tkMessageBox.showinfo(title='ERROR', message=err_msg)
            return False

        # each object in the selection is treated as an independent struc
        objlist = cmd.get_object_list(sel_name)
        self.ss_asgn_prog = 'DSSP'
        print('Starting %s ...' % (self.ss_asgn_prog, ))

        for objname in objlist:
            self.sel_obj_list.append('%s and %s' % (sel_name, objname))
            self.runDSSPOneObj(self.sel_obj_list[-1])

# cmd.delete(tmpsel)

        return True

    def runStrideOneObj(self, one_obj_sel):
        """ Run Stride on only one object.

            @param one_obj_sel: the selection/object involving only one object
            @param type: string
        """
        SSE_res = {
            'H': {}, 'G': {}, 'I': {},
            'E': {}, 'B': {}, 'b': {},
            'T': {}, 'C': {}
        }
        SSE_sel = {
            'H': None, 'G': None, 'I': None,
            'E': None, 'B': None, 'b': None,
            'T': None, 'C': None
        }

        pdb_fn = None
        pdb_os_fh, pdb_fn = tempfile.mkstemp(suffix='.pdb')  # file os handle, file name
        os.close(pdb_os_fh)
        cmd.save(filename=pdb_fn, selection=one_obj_sel)
        if VERBOSE:
            print('Selection %s saved to %s.' % (one_obj_sel, pdb_fn))

        if pdb_fn is None:
            print('WARNING: Stride has no pdb file to work on!')
            return None

        print('Running Stride for %s ...' % (one_obj_sel,))
        stride_sse_dict = {}
        if sys.version_info >= (2, 4):
            stride_proc = subprocess.Popen([self.stride_bin.get(), pdb_fn],
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.STDOUT
                                           )
            stride_stdout, stride_stderr = stride_proc.communicate()
        else:  # use os.system + tempfile
            stride_tmpout_os_fh, stride_tmpout_fn = tempfile.mkstemp(suffix='.stride')
            os.close(stride_tmpout_os_fh)
            stride_cmd = '%s %s > %s' % (self.stride_bin.get(), pdb_fn, stride_tmpout_fn)
            os.system(stride_cmd)
            fh = open(stride_tmpout_fn)
            stride_stdout = ''.join(fh.readlines())
            fh.close()

        for line in stride_stdout.splitlines():
            if not isinstance(line, str):  # Python 3
                line = line.decode('ascii', 'ignore')
            if line.startswith('ASG'):
                resname, ch = line[5:8], line[9]
                # according to stride doc, col 12-15 are used for residue number
                # actually, resnum+icode is used in 12-15.
                # In addition, if residue number is 4-digit or longer,
                # the field 12-16 or more are used.

                if line[15] == ' ':  # col 16 is not occupied
                    residen, sscode = line[11:15].strip(), line[24]  # residen = resnum+icode
                else:
                    shift = line[15:].find(' ')  # find where residen stops
                    residen, sscode = line[11:15 + shift].strip(), line[24 + shift]

                k = (ch, resname, residen)
                stride_sse_dict[k] = sscode

        self.stride_rlt_dict[one_obj_sel] = stride_sse_dict
        print('Got SSE for %d residues.' % (len(self.stride_rlt_dict[one_obj_sel]),))

        # group residues according to their SSE, and chain name
        for k in self.stride_rlt_dict[one_obj_sel].keys():
            sse = self.stride_rlt_dict[one_obj_sel][k]
            chn = k[0]  # this can be a space!
            res = k[2]
            if VERBOSE:
                print('(%s) and \"%s\"/%s/  sse=%s' % (str(one_obj_sel),
                                                               chn.strip(), res, sse))
            SSE_res[sse].setdefault(chn, []).append(res)

        self.SSE_res_dict[one_obj_sel] = SSE_res

        for sse in self.STRIDE_SSE_list:
            sse_sel_name = self.selectSSE(one_obj_sel, sse)
            SSE_sel[sse] = sse_sel_name
        self.SSE_sel_dict[one_obj_sel] = SSE_sel

        print('\nNumber of residues with SSE element:')
        for sse in self.STRIDE_SSE_list:
            num = sum([len(SSE_res[sse][chn]) for chn in SSE_res[sse]])
            print('%20s (%s) : %5d' % (self.SSE_name[sse], sse, num))
        print('')

        # clean up pdb_fn and dssp_tmpout_fn created by tempfile.mkstemp()
        if os.path.isfile(pdb_fn):
            os.remove(pdb_fn)
        if sys.version_info < (2, 4) and os.path.isfile(stride_tmpout_fn):
            os.remove(stride_tmpout_fn)

        return True

    def runStride(self):
        """
        """
        # delete old results
        self.sel_obj_list = []
        self.stride_rlt_dict = {}
        self.SSE_res_dict = {}
        self.SSE_sel_dict = {}

        pdb_fn = None
        sel_name = None
        sel = self.pymol_sel.get()

        if len(sel) > 0:  # if any pymol selection/object is specified
            all_sel_names = cmd.get_names('all')  # get names of all selections
            if sel in all_sel_names:
                if cmd.count_atoms(sel) == 0:
                    err_msg = 'ERROR: The selection %s is empty.' % (sel,)
                    print('ERROR: %s' % (err_msg,))
                    tkMessageBox.showinfo(title='ERROR', message=err_msg)
                    return False
                else:
                    sel_name = sel
            # no selection/object with the input name is found
            # we assume either a single-word selector or
            # some other selection-expression is uesd
            else:
                print('The selection/object you specified is not found.')
                print('Your input will be interpreted as a selection-expression.')
                tmpsel = self.randomSeleName(prefix='your_sele_')
                cmd.select(tmpsel, sel)
                if cmd.count_atoms(tmpsel) == 0:
                    cmd.delete(tmpsel)
                    err_msg = 'ERROR: The selection %s is empty.' % (sel,)
                    print('ERROR: %s' % (err_msg,))
                    tkMessageBox.showinfo(title='ERROR', message=err_msg)
                    return False
                else:
                    sel_name = tmpsel

        else:   # what structure do you want Stride to work on?
            err_msg = 'No PyMOL selection/object specified!'
            print('ERROR: %s' % (err_msg,))
            tkMessageBox.showinfo(title='ERROR', message=err_msg)
            return False

        # each object in the selection is treated as an independent struc
        objlist = cmd.get_object_list(sel_name)
        self.ss_asgn_prog = 'Stride'
        print('Starting %s ...' % (self.ss_asgn_prog, ))

        for objname in objlist:
            self.sel_obj_list.append('%s and %s' % (sel_name, objname))
            self.runStrideOneObj(self.sel_obj_list[-1])

        return True

    def randomSeleName(self, prefix='sele', suffix=''):
        """ generate a random selection name.
        """
        sel_list = cmd.get_names('all')
        sel_dict = dict(zip(sel_list, range(len(sel_list))))
        sel_name = '%s%d%s' % (prefix, random.randint(1000, 9999), suffix)
        while(sel_name in sel_dict):
            sel_name = '%s%d%s' % (prefix, random.randint(1000, 9999), suffix)

        return sel_name

    def selectSSE(self, sel, sse):
        """ generate selector for selecting all residues having the given sse.
            return the selection name.
        """
        #sel = self.pymol_sel.get()
        sel_list_chn = []

        if VERBOSE:
            print('\nSelecting SSE %s ... \n' % (sse))

        for chn in self.SSE_res_dict[sel][sse]:  # color one chain at a time
            if chn == ' ':
                chn_str = '-'
            else:
                chn_str = chn
            if VERBOSE:
                print('Selecting SSE %s on chain %s ... \n' % (sse, chn))
            limit = 150  # color every 150 residues
            sel_name_chn = self.randomSeleName(prefix='%s_%s_%s_' % ('_'.join(sel.split()), chn_str, sse))

            if len(self.SSE_res_dict[sel][sse][chn]) < limit:
                #sel_expr = '/%s//%s/%s/' % (sel,chn.strip(), '+'.join(self.SSE_res[sse][chn]))
                # always quote chain name in case it is empty (otherwise sel misinterpreted)
                sel_expr = '(%s) and \"%s\"/%s/' % (sel, chn.strip(), '+'.join(self.SSE_res_dict[sel][sse][chn]))
                cmd.select(sel_name_chn, sel_expr)
                if VERBOSE:
                    print('select %s, %s' % (sel_name_chn, sel_expr))
                sel_list_chn.append(sel_name_chn)
            else:
                rn = len(self.SSE_res_dict[sel][sse][chn])
                print('total number of res with %s = %d' % (sse, rn))
                sz = int(math.ceil(rn / float(limit)))
                sel_list_seg = []
                for i in range(sz):
                    s, e = i * limit, min((i + 1) * limit, rn)
                    print(s, e)
                    #sel_expr = '/%s//%s/%s/' % (sel,chn.strip(), '+'.join(self.SSE_res[sse][chn][s:e]))
                    # always quote chain name in case it is empty (otherwise sel misinterpreted)
                    sel_expr = '(%s) and \"%s\"/%s/' % (sel, chn.strip(),
                                                        '+'.join(self.SSE_res_dict[sel][sse][chn][s:e]))
                    sel_name_seg = self.randomSeleName(prefix='%s_%s_%s_tmp_' % ('_'.join(sel.split()), chn_str, sse))
                    cmd.select(sel_name_seg, sel_expr)
                    if VERBOSE:
                        print('select %s, %s' % (sel_name_seg, sel_expr))
                    sel_list_seg.append(sel_name_seg)

                sel_expr = ' or '.join(sel_list_seg)
                cmd.select(sel_name_chn, sel_expr)
                if VERBOSE:
                    print('select %s, %s' % (sel_name_chn, sel_expr))
                [cmd.delete(asel) for asel in sel_list_seg]
                sel_list_chn.append(sel_name_chn)

        if len(sel_list_chn) > 0:
            sel_name = self.randomSeleName(prefix='%s_%s_%s_' % ('_'.join(sel.split()), sse, self.ss_asgn_prog))
            sel_expr = ' or '.join(sel_list_chn)
            cmd.select(sel_name, sel_expr)
            [cmd.delete(asel) for asel in sel_list_chn]
        else:
            print('INFO: No residues are assigned to SSE \'%s\'.' % (sse,))
            sel_name = None

        return sel_name

    def updateColor(self):
        if self.ss_asgn_prog is None:
            err_msg = 'Run DSSP or Stride to assign secondary structures first!'
            print('ERROR: %s' % (err_msg,))
            tkMessageBox.showinfo(title='ERROR', message=err_msg)
        else:
            print('Update color for %s' % (self.pymol_sel.get()), end=' ')
            print('using secondary structure assignment by %s' % (self.ss_asgn_prog,))

            if self.ss_asgn_prog == 'DSSP':
                SSE_list = self.DSSP_SSE_list
            elif self.ss_asgn_prog == 'Stride':
                SSE_list = self.STRIDE_SSE_list

            for sse in SSE_list:  # give color names
                cmd.set_color('%s_color' % (sse,), self.SSE_col_RGB[sse])
            for sse in SSE_list:  # color each SSE
                for sel_obj in self.sel_obj_list:
                    if self.SSE_sel_dict[sel_obj][sse] is not None:
                        cmd.color('%s_color' % (sse,), self.SSE_sel_dict[sel_obj][sse])
                        print('color', self.SSE_sel_dict[sel_obj][sse], ',', self.SSE_col_RGB[sse])
                    else:
                        print('No residues with SSE \'%s\' to color.' % (sse,))

        return

    def updateSS(self):
        if self.ss_asgn_prog is None:
            err_msg = 'Run DSSP or Stride to assign secondary structures first!'
            print('ERROR: %s' % (err_msg,))
            tkMessageBox.showinfo(title='ERROR', message=err_msg)
        else:
            print('Update secondary structures for %s' % (self.pymol_sel.get()), end=' ')
            print('using secondary structure assignment by %s' % (self.ss_asgn_prog,))
            print('SSE mapping: (H,G,I) ==> H; (E,B,b) ==> S; (T,S,-,C) ==> L')

            if self.ss_asgn_prog == 'DSSP':
                SSE_list = self.DSSP_SSE_list
            elif self.ss_asgn_prog == 'Stride':
                SSE_list = self.STRIDE_SSE_list

            for sse in SSE_list:
                for sel_obj in self.sel_obj_list:
                    if self.SSE_sel_dict[sel_obj][sse] is not None:
                        cmd.alter(self.SSE_sel_dict[sel_obj][sse], 'ss=\'%s\'' % (self.SSE_map[sse],))
                        print('alter %s, ss=%s' % (self.SSE_sel_dict[sel_obj][sse], self.SSE_map[sse]))
                    else:
                        print('No residue with SSE %s to update ss.' % (sse,))

            # show cartoon for the input selection, and rebuild
            cmd.show('cartoon', self.pymol_sel.get())
            cmd.rebuild(self.pymol_sel.get())
            return

    def custermizeHColor(self):
        self.custermizeSSEColor('H')

    def custermizeGColor(self):
        self.custermizeSSEColor('G')

    def custermizeIColor(self):
        self.custermizeSSEColor('I')

    def custermizeEColor(self):
        self.custermizeSSEColor('E')

    def custermizeBColor(self):
        self.custermizeSSEColor('B')

    def custermizeTColor(self):
        self.custermizeSSEColor('T')

    def custermizeSColor(self):
        self.custermizeSSEColor('S')

    def custermizeNColor(self):
        self.custermizeSSEColor('-')

    def custermizebColor(self):
        self.custermizeSSEColor('b')

    def custermizeCColor(self):
        self.custermizeSSEColor('C')

    def custermizeSSEColor(self, sse):
        SSE_col_but = {
            'H': self.H_col_but,
            'G': self.G_col_but,
            'I': self.I_col_but,

            'E': self.E_col_but,
            'B': self.B_col_but,

            'T': self.T_col_but,
            'S': self.S_col_but,
            '-': self.N_col_but,

            'b': self.b_col_but,
            'C': self.C_col_but,
        }
        try:
            color_tuple, color = tkColorChooser.askcolor(color=self.SSE_col[sse])
            if color_tuple is not None and color is not None:
                self.SSE_col_RGB[sse] = color_tuple
                self.SSE_col[sse] = color

                SSE_col_but[sse]['bg'] = self.SSE_col[sse]
                SSE_col_but[sse]['activebackground'] = self.SSE_col[sse]
                SSE_col_but[sse].update()
        except Tkinter._tkinter.TclError:
            print('Old color (%s) will be used.' % (self.mesh_col))

    def execute(self, butcmd):
        """ Run the cmd represented by the botton clicked by user.
        """
        if butcmd == 'OK':
            print('is everything OK?')

        elif butcmd == 'Run DSSP':
            rtn = self.runDSSP()
            if rtn and VERBOSE:
                print('Done with DSSP!')

        elif butcmd == 'Run Stride':
            rtn = self.runStride()
            if rtn and VERBOSE:
                print('Done with Stride!')

        elif butcmd == 'Update Color':
            self.updateColor()

        elif butcmd == 'Update ss':
            self.updateSS()

        elif butcmd == 'Exit':
            print('Exiting DSSP and Stride Plugin ...')
            if __name__ == '__main__':
                self.parent.destroy()
            else:
                self.dialog.withdraw()
            print('Done.')
        else:
            print('Exiting DSSP and Stride Plugin because of unknown button click ...')
            self.dialog.withdraw()
            print('Done.')

    def quit(self):
        self.dialog.destroy()


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

    widget = DSSPPlugin(app)
    app.root.mainloop()
