# Autodock/Vina plugin  Copyright Notice
# ============================
#
# The Autodock/Vina plugin source code is copyrighted, but you can freely use and
# copy it as long as you don't change or remove any of the copyright
# notices.
#
# ----------------------------------------------------------------------
# Autodock/Vina plugin is Copyright (C) 2009 by Daniel Seeliger
#
#                        All Rights Reserved
#
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name of Daniel Seeliger not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
#
# DANIEL SEELIGER DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
# SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS.  IN NO EVENT SHALL DANIEL SEELIGER BE LIABLE FOR ANY
# SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
# RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
# CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------

#================================================================
import sys,os,math,re, fnmatch, shutil
from os import stat
from os.path import abspath
from stat import ST_SIZE
from time import sleep, time
import Tkinter
import Tkinter,Pmw
from Tkinter import *
import tkMessageBox, tkFileDialog
import Pmw
from threading import Thread
from commands import getstatusoutput
from pymol import cmd,selector
from pymol.cmd import _feedback,fb_module,fb_mask,is_list,_cmd
from pymol.cgo import *
from pymol import stored
from numpy import *
import tkColorChooser
from pymol.vfont import plain
from glob import glob

__version__ = "2.1.1"
#=============================================================================
#
#     INITIALISE PLUGIN
#


def __init__(self):
    self.menuBar.addmenuitem('Plugin', 'command',
    'Launch Autodock',
    label='Autodock/Vina',
    command = lambda s=self: Autodock(s))
        
    cmd.set("retain_order") # keep atom ordering


#=============================================================================
#
#     SOME BASIC THINGIES
#

intro_text = """
Welcome to the updated version of the PyMOL Autodock Plugin.The current version has been
entirely rewritten and extended. Now you can set up a complete docking run from the very
beginning, perform the docking runs from within PyMOL and directly load the results into
the viewer. The plugin now supports the novel Autodock spawn \"VINA\" which is orders of
magnitudes faster than Autodock4. To get this plugin to work properly you should have the
mgltools (http://mgltools.scripps.edu/) installed and tell the plugin on this page where
to find the scripts (usually somewhere in /python/site-packages/AutoDockTools/Utilities24/ )
and make sure that they work. Additionally you should tell the plugin where to find the autogrid4,
autodock4, and vina executables. If you do this once and press the \"Save Configuration File\"
button this information is present the next time you use the plugin.

A complete setup and execution of a docking run or a virtual screening is basically a walk
from the left to the right along the notebook structure if this plugin. Each page is more
or less intuitive or contains a description of what to do. The basic workflow is:
Load Structure -> Define binding site -> Receptor preparation -> Ligand preparation ->
Docking -> Analysis.

Have fun and contact me (dseelig@gwdg.de) if you find bugs or have any suggestions for future
versions.
Daniel
Reference: J. Comput.-Aided Mol. Des. 24:417-422 (2010)
"""

receptor_prep_text = """
Here you can define receptors from PyMOL selections and setup 
docking runs with flexible sidechains. Pick a selection from the 
selection list and press the Generate Receptor ->" button to prepare 
the protein for a docking run. If you want to use flexible sidechains 
within the binding site just create a PyMOL selection containing the 
residues you want to be flexible. Then use the "Import Selections" 
button and andselect the group in the left list. Then press the 
Select as Flexible->" button. The receptor definition has now changed. 
You can see which residues are flexible in the right list. (Proline, 
Glycine and Alanine residues are never defined as flexible)\n"""

ligand_prep_text = """
On this page you can prepare ligands for docking.The first way 
to do so is to load the ligand into PyMOL and select it in the left 
list (use the "Import Selection" button to synchronize the list with 
PyMOL). The plugin saves your ligand in as .pdb file and puts it 
into the ligand list in the middle of the page. If you press the 
"Generate Ligand" button your ligand will be prepared for docking 
and appears in the right list. Alternatively to loading a ligand via 
PyMOL you can choose an entire directory containing compounds 
in .pdb or .mol2 format. Use the "Import" button to read all 
compounds in the directory and use the "Generate All" button to 
prepare each compound for the docking run.\n"""

docking_text = """
This is the docking page. You have to select a recetor and 
whether you want to use flexible sidechains. You may either 
dock a single ligand from the list above or all ligands you 
have prepared. If you use Autodock4 you have to calculate 
a grid before you can start the actual docking run. Just 
press "Run AutoGrid" and wait until its finished and then 
press the "Run AutoDock" button. The genrated grid maps 
can be loaded into PyMOL on the "Grid Maps" page. If you 
use vina you don't have to calculate a grid. Just press 
"Run VINA" and wait until it is finished.You can load the 
generated docking poses into PyMOL on the "View Poses" page\n"""


pdb_format="%6s%5d %-4s%1s%3s%2s%4d %11.3f %7.3f %7.3f %5.2f %5.2f"

# get a temporary file directory
if not sys.platform.startswith('win'):
    home = os.environ.get('HOME')
else:
    home = os.environ.get('PYMOL_PATH')

tmp_dir = os.path.join(home,'.ADplugin')
if not os.path.isdir(tmp_dir):
    os.mkdir(tmp_dir)
    print "Created temporary files directory:  %s" % tmp_dir
    
default_settings = {
    "grid_spacing" : '0.375',
    "n_points_X":'60',
    "n_points_Y":'60',
    "n_points_Z":'60',
    "grid_center_selection":'(all)',
    "grid_center_X":'0',
    "grid_center_Y":'0',
    "grid_center_Z":'0',
    "dx":'1.0',
    "dy":'1.0',
    "dz":'1.0',
    "box_cylinder_size":'0.2',
    "box_mesh_line_width":'1',
    "box_mesh_grid_size":'1',
    "box_file_name":'box.dat',
    "gpf_file_name":'grid.gpf',
    "config_file_name":'config.txt',
    "rank_dat_file_name":'scores.dat',
    "rank_csv_file_name":'scores.csv',
    "rank_pose_file_name":'poses.pdb',
    "dlg_input_file":'docked.pdbqt',
    "map_input_file":'receptor.C.map',
    "map_threshold":5.
    }

BOX_AS_BOX = 0
BOX_AS_WIREBOX = 1
GRID_CENTER_FROM_SELECTION = 0
GRID_CENTER_FROM_COORDINATES = 1



#==========================================================================
#
#    THREAD CLASSES FOR SPAWNING DOCKING JOBS
        

class Thread_run(Thread):
    def __init__ (self,command, previous = None, status_line = None, log_text = None):
        Thread.__init__(self)
        self.command = command
        self.status = -1
        self.previous = previous
        self.status_line = status_line
        self.log_text = log_text
    def run(self):
        if self.previous:
            self.previous.join()
        if self.status_line:
            self.status_line.configure(text = self.log_text)
        self.status = os.system(self.command)

class Thread_log(Thread):
    def __init__(self, logfile, page):
        Thread.__init__(self)
        self.page = page
        self.logfile = logfile
    def run(self):
        if not os.environ.has_key('ADPLUGIN_NO_OUTPUT_REDIRECT'):
            t = Tail(self.logfile)
            line = t.nextline()
            self.page.insert('end',"%s" % line)
            while line:
                line = t.nextline()
                self.page.insert('end',"%s" % line)
                self.page.yview('moveto', 1.0)#, 'page')
        else:
            line = 'LOG FILE OUTPUT NOT REDIRECTED'
            self.page.insert('end',"%s" % line)
            self.page.yview('moveto', 1.0)#, 'page')
            
#==========================================================================
#
#    CLASSES FOR HANDLING AUTODOCK FILES


class ADModel:
    """ STORAGE CLASS FOR DOCKED LIGANDS """
    
    def __init__(self, lst = None):
        self.atomlines = []
        self.energy = 0.
        self.name = ''
        self.poseN = 0
        self.info = ''
        self.lst = []
        self.num = 0
        self.as_string = ''
        if lst is not None:
            self.from_list(lst)
            
    def from_list(self, lst):
        self.lst = lst
        for line in lst:
            self.info+=line
            if line.startswith('ATOM') or \
                   line.startswith('HETATM'):
                self.atomlines.append(line)
                self.as_string+='ATOM  '+line[6:67]+'\n'
            elif line.startswith('USER'):
                if 'Free Energy of Binding' in line:
                    entr = line.split('=')[1]
                    self.energy = float(entr.split()[0])
            elif line.startswith('REMARK'):
                if 'VINA RESULT' in line:
                    entr = line.split(':')[1]
                    self.energy = float(entr.split()[0])
                    
    def as_pdb_string(self):
        return self.as_string
    def info_string(self):
        return self.info


#---------------------------------------------------------------------------

class ADGridMap:
    """ CLASS FOR HANDLING AUTODOCK GRID MAP FILES"""
    
    def __init__(self, fp = None, name = 'map'):

        self.name = ''
        self.npts = [0,0,0]
        self.n = [0,0,0]
        self.center = [0,0,0]
        self.origin = [0,0,0]
        self.nelem = 0
        self.spacing = 0.
        self.values = []
        self.datafile = ''
        self.molecule = ''
        self.paramfile = ''
        self.precision = 0.0001
        if fp is not None:
            self.read(fp,name)
            
    def read(self,fp,name='map'):
        self.name = name
        for i in range(6):
            line = fp.readline()
            if i == 0:
                self.paramfile = line.split()[1]
            elif i == 1:
                self.datafile = line.split()[1]
            elif i == 2:
                self.molecule = line.split()[1]
            elif i == 3:
                self.spacing = float(line.split()[1])
            elif i == 4:
                self.npts = [int(x) for x in line.split()[1:]]
            elif i == 5:
                self.center = [float(x) for x in line.split()[1:]]
        for i in range(3):
            self.n[i] = self.npts[i]+1
        self.nelem=self.n[0]*self.n[1]*self.n[2]
        i = 0
        while i < self.nelem:
            val = float(fp.readline())
            self.values.append(val)
            i+=1
        for i in range(3):
            self.origin[i] = self.center[i]-self.npts[i]/2*self.spacing
            
    def meta(self):
        s= 'GRID_PARAMETER_FILE %s\n' % self.paramfile + \
           'GRID_DATA_FILE %s\n' % self.datafile +\
           'MACROMOLECULE %s\n' % self.molecule +\
           'SPACING %4.3f\n' % self.spacing +\
           'NELEMENTS %d %d %d\n' % (self.npts[0],self.npts[1],self.npts[2]) +\
           'CENTER %5.3f %5.3f %5.3f\n' % (self.center[0],self.center[1],self.center[2]) 
        return s

    def write(self,fp):
        print >>fp, 'GRID_PARAMETER_FILE %s' % self.paramfile
        print >>fp, 'GRID_DATA_FILE %s' % self.datafile
        print >>fp, 'MACROMOLECULE %s' % self.molecule
        print >>fp, 'SPACING %4.3f' % self.spacing
        print >>fp, 'NELEMENTS %d %d %d' % (self.npts[0],self.npts[1],self.npts[2])
        print >>fp, 'CENTER %5.3f %5.3f %5.3f' % (self.center[0],self.center[1],self.center[2])
        for x in self.values:
            if abs(x) < self.precision:
                print >>fp, '0.'
            else:
                print >>fp, '%.3f' % x

    def writeDX(self,fname):
        fp = open(fname,'w')
        nx = self.n[0]
        ny = self.n[1]
        nz = self.n[2]
        ori = self.origin
        spacing = self.spacing
        vals = self.values

        print >>fp,'#=================================='
        print >>fp,'# AutoGrid Map File: %s' % self.name
        print >>fp,'# Receptor File Name: %s' % self.molecule
        print >>fp,'#=================================='
        print >>fp,'object 1 class gridpositions counts %d %d %d' % (nx,ny,nz)
        print >>fp,'origin %12.5E %12.5E %12.5E' % (ori[0],ori[1],ori[2])
        print >>fp,'delta %12.5E %12.5E %12.5E' % (spacing,0,0)
        print >>fp,'delta %12.5E %12.5E %12.5E' % (0,spacing,0)
        print >>fp,'delta %12.5E %12.5E %12.5E' % (0,0,spacing)
        print >>fp,'object 2 class gridconnections counts %d %d %d' % (nx,ny,nz)
        print >>fp,'object 3 class array type double rank 0 items %d data follows' % len(vals)
        for k in range(nz):
            col=0;
            for j in range(ny):
                for i in range(nx):
                    fp.write(" %12.5E" % vals[i*ny*nz + j*nz + k])
                    col+=1;
                    if col==3:
                        print >>fp
                        col=0
        print >>fp,'attribute \"dep\" string \"positions\"'
        print >>fp,'object \"regular positions regular connections\" class field'
        print >>fp,'component \"positions\" value 1'
        print >>fp,'component \"connections\" value 2'
        print >>fp,'component \"data\" value 3'
        fp.close()

#==========================================================================
#
#    CLASSES FOR INTERNAL HANDLING OF RECEPTORS AND LIGANDS

class Receptor:
    """CONTAINS ALL INFORMATION ABOUT A DEFINED RECEPTOR"""    
    def __init__(self):
        self.selection = ''
        self.pdb_file = ''
        self.receptor_pdbqt = ''
        self.receptor_rigid = ''
        self.receptor_flexible = ''
        self.flexible_residues = {}
        self.resi_dic = {}
        
    def flex_res_string(self):
        flex_res_str = ''
        for key, val in self.flexible_residues.items():
            s = os.path.basename(self.receptor_pdbqt)[:-6]
            lst = []
            for idx, resname in val:
                lst.append( resname + str(idx) )
            s+=':'+key+':'+'_'.join(lst)
            flex_res_str+=s
        return flex_res_str
##         lst = []
##         for idx, resname in self.flexible_residues:
##             lst.append( resname+str(idx) )
##         return '_'.join(lst)

    def info(self):
        s= '#===============================================\n'
        s+=' > Receptor                 : "%s"\n' % self.selection
        s+=' > Generated from pdb file  : "%s"\n' % self.pdb_file
        s+=' > Receptor file            : "%s"\n' % self.receptor_pdbqt
        s+=' > Receptor rigid           : "%s"\n' % self.receptor_rigid
        s+=' > Receptor flexible        : "%s"\n' % self.receptor_flexible
#        s+=' > Number of flex. residues : %d\n' % len(self.flexible_residues)
        s+='#===============================================\n'
        return s

class Ligand:
    """CONTAINS ALL INFORMATION ABOUT A LIGAND"""
    def __init__(self):
        self.name = ''
        self.selection = ''
        self.input_file = ''
        self.ligand_pdbqt = ''
        self.outfile_poses = ''

    def info(self):
        s= '#=============================================\n'
        s+=' > Ligand                  : %s\n' % self.name
        s+=' > Generated from file     : %s\n' % self.input_file
        s+=' > Ligand pdbqt file       : %s\n' % self.ligand_pdbqt
        s+=' > Poses output file       : %s\n' % self.outfile_poses
        s+='#=============================================\n'
        return s
    
#==========================================================================
#
#    THE MAJOR, PRETTY BIG PLUGIN CLASS

class Autodock:
    """ THE MAJOR PLUGIN CLASS """
    
    def __init__(self,app):
        parent = app.root
        self.parent = parent
        # receptors and ligands
        self.receptor_dic = {}
        self.ligand_dic = {}
        # directory with ligands
        self.ligand_dir = StringVar()
        
        # box display settings
        self.box_display_mode = IntVar()
        self.box_display_mode.set(BOX_AS_BOX)
        self.box_color = [1.,1.,1.]
        self.box_is_on_display = False
        self.box_display_cylinder_size = DoubleVar()
        self.box_display_cylinder_size.set(default_settings['box_cylinder_size'])
        self.box_display_line_width = DoubleVar()
        self.box_display_line_width.set(default_settings['box_mesh_line_width'])
        self.box_display_mesh_grid = DoubleVar()
        self.box_display_mesh_grid.set(default_settings['box_mesh_grid_size'])
        self.box_size = []
        
        # grid definition

        self.grid_spacing = DoubleVar()
        self.grid_spacing.set(default_settings['grid_spacing'])

        self.n_points_X = DoubleVar()
        self.n_points_X.set(default_settings['n_points_X'])
        self.n_points_Y = DoubleVar()
        self.n_points_Y.set(default_settings['n_points_Y'])
        self.n_points_Z = DoubleVar()
        self.n_points_Z.set(default_settings['n_points_Z'])

        self.grid_center = [DoubleVar(), DoubleVar(), DoubleVar()]
        self.grid_center[0].set(default_settings['grid_center_X'])
        self.grid_center[1].set(default_settings['grid_center_Y'])
        self.grid_center[2].set(default_settings['grid_center_Z'])

        # paths to executables
        self.config_settings = {}
        self.autodock_tools_path = StringVar()
        self.autogrid_exe = StringVar()
        self.autodock_exe = StringVar()
        self.vina_exe = StringVar()

        # ligand display settings

        self.ligand_display_mode = {
                'lines':True,
                'sticks':False,
                'spheres':False,
                'surface':False,
                'mesh':False
                }
        # keep in mind what threads are running
        self.current_thread = None

        # build main window
        
        self.dialog = Pmw.Dialog(parent,
                                 buttons = ('Exit',),
                                 title = 'PyMOL Autodock/Vina Plugin',
                                 command = self.button_pressed)
        self.dialog.withdraw()
        Pmw.setbusycursorattributes(self.dialog.component('hull'))
        self.status_line = Label(self.dialog.interior(), 
                                 relief='sunken',
                                 font='helvetica 12', anchor='w',fg='yellow',bg='black')
        self.status_line.pack(side=BOTTOM,fill='x', expand=1, padx=0, pady=0)


        self.dialog.geometry('650x780')
        self.dialog.bind('<Return>',self.button_pressed)

        # the title

        self.title_label = Tkinter.Label(self.dialog.interior(),
                                         text = 'PyMOL Autodock/Vina Plugin\nDaniel Seeliger\n<http://wwwuser.gwdg.de/~dseelig>',
                                         background = 'navy',
                                         foreground = 'white',
                                         )
        self.title_label.pack(expand = 0, fill = 'both', padx = 4, pady = 4)



        # the basic notebook

        self.notebook = Pmw.NoteBook(self.dialog.interior())
        self.notebook.pack(fill='both',expand=1,padx=3,pady=3)



        # build pages
        self.configuration_page = self.notebook.add('Configuration')

        self.grid_definition_page = self.notebook.add('Grid Settings')
        self.receptor_preparation_page = self.notebook.add('Receptor')
        self.ligand_preparation_page = self.notebook.add('Ligands')
        self.docking_page = self.notebook.add('Docking')
        self.pose_viewer_page = self.notebook.add('View Poses')
        self.rank_page = self.notebook.add('Score/Rank')
        self.map_viewer_page = self.notebook.add('Grid Maps')

        #---------------------------------------------------------------
        # GRID DEFINITION PAGE
        self.grid_page_main_group = Pmw.Group(self.grid_definition_page, tag_text='Grid Definition')
        self.grid_page_main_group.pack(fill = 'both', expand = 0, padx=10, pady=5)

        # grid parameters on the left
        self.grid_page_left_side = Pmw.Group(self.grid_page_main_group.interior(),tag_text = 'Parameters')
        self.grid_page_left_side.pack(side = LEFT, fill = 'both',expand = 0, padx = 10, pady = 3)

        # grid spacing entry
        self.grid_spacing_frame = Tkinter.Frame(self.grid_page_left_side.interior())
        self.grid_spacing_label = Label(self.grid_spacing_frame, text='Spacing:', width = 10)
        self.grid_spacing_location = Entry(self.grid_spacing_frame, textvariable= self.grid_spacing, bg='black',fg='green', width = 10)
        self.grid_spacing_scrollbar = Scrollbar(self.grid_spacing_frame, orient = 'horizontal', command = self.grid_spacing_changed)
        self.grid_spacing_label.pack(side = LEFT, anchor='w')
        self.grid_spacing_location.pack(side = LEFT, anchor='w')
        self.grid_spacing_scrollbar.pack(side = LEFT, anchor='w')
        self.grid_spacing_frame.pack(fill = 'x', padx = 4, pady=1)

        # n grid points entries
        self.n_points_X_frame = Tkinter.Frame(self.grid_page_left_side.interior())
        self.n_points_X_label = Label(self.n_points_X_frame, text='X-points:', width = 10)
        self.n_points_X_location = Entry(self.n_points_X_frame, textvariable= self.n_points_X, bg='black',fg='green', width = 10)
        self.n_points_X_scrollbar = Scrollbar(self.n_points_X_frame, orient = 'horizontal', command = self.n_points_X_changed)
        self.n_points_X_label.pack(side = LEFT, anchor='w')
        self.n_points_X_location.pack(side = LEFT, anchor='w')
        self.n_points_X_scrollbar.pack(side = LEFT, anchor='w')
        self.n_points_X_frame.pack(fill = 'x', padx = 4, pady=1)

        self.n_points_Y_frame = Tkinter.Frame(self.grid_page_left_side.interior())
        self.n_points_Y_label = Label(self.n_points_Y_frame, text='Y-points:', width = 10)
        self.n_points_Y_location = Entry(self.n_points_Y_frame, textvariable= self.n_points_Y, bg='black',fg='green', width = 10)
        self.n_points_Y_scrollbar = Scrollbar(self.n_points_Y_frame, orient = 'horizontal', command = self.n_points_Y_changed)
        self.n_points_Y_label.pack(side = LEFT, anchor='w')
        self.n_points_Y_location.pack(side = LEFT, anchor='w')
        self.n_points_Y_scrollbar.pack(side = LEFT, anchor='w')
        self.n_points_Y_frame.pack(fill = 'x', padx = 4, pady=1)

        self.n_points_Z_frame = Tkinter.Frame(self.grid_page_left_side.interior())
        self.n_points_Z_label = Label(self.n_points_Z_frame, text='Z-points:', width = 10)
        self.n_points_Z_location = Entry(self.n_points_Z_frame, textvariable= self.n_points_Z, bg='black',fg='green', width = 10)
        self.n_points_Z_scrollbar = Scrollbar(self.n_points_Z_frame, orient = 'horizontal', command = self.n_points_Z_changed)
        self.n_points_Z_label.pack(side = LEFT, anchor='w')
        self.n_points_Z_location.pack(side = LEFT, anchor='w')
        self.n_points_Z_scrollbar.pack(side = LEFT, anchor='w')
        self.n_points_Z_frame.pack(fill = 'x', padx = 4, pady=1)

        Pmw.alignlabels( [self.grid_spacing_label,
                          self.n_points_X_label,
                          self.n_points_Y_label,
                          self.n_points_Z_label
                          ])

        Pmw.alignlabels( [self.grid_spacing_location,
                          self.n_points_X_location,
                          self.n_points_Y_location,
                          self.n_points_Z_location
                          ])

        # display option buttons
        self.display_button_box = Pmw.ButtonBox(self.grid_page_main_group.interior(), padx=0, pady=1,orient='vertical')
        self.display_button_box.pack(side=LEFT)
        self.display_button_box.add('Show Box',command = self.show_box)
        self.display_button_box.add('Hide Box',command = self.hide_box)
        self.display_button_box.add('Change Box Color',command = self.change_box_color)
        

        # display options on the right
        self.grid_page_right_side = Pmw.Group(self.grid_page_main_group.interior(), tag_text = 'Display Options')
        self.grid_page_right_side.pack(side = LEFT, fill = 'both', expand = 0, padx = 10, pady = 3)

        self.box_display_radiogroups = []
        self.box_display_radioframe = Tkinter.Frame(self.grid_page_right_side.interior())

        self.box_display_cylinder_frame = Pmw.Group(self.box_display_radioframe,
                                                    tag_pyclass = Tkinter.Radiobutton,
                                                    tag_text = 'Cylindric Box',
                                                    tag_value = 0,
                                                    tag_variable = self.box_display_mode
                                                    )
        self.box_display_cylinder_frame.pack(fill='x',expand = 1,side=TOP)
        self.box_display_radiogroups.append(self.box_display_cylinder_frame)

        self.box_display_cylinder_size_frame = Tkinter.Frame(self.box_display_cylinder_frame.interior())
        self.box_display_cylinder_size_label = Label(self.box_display_cylinder_size_frame, text = 'Size:', width=10)
        self.box_display_cylinder_size_location = Entry(self.box_display_cylinder_size_frame,
                                                        textvariable = self.box_display_cylinder_size,
                                                        bg = 'black',fg = 'green', width = 10)
        self.box_display_cylinder_size_scrollbar = Scrollbar(self.box_display_cylinder_size_frame,
                                                             orient = 'horizontal',command = self.box_display_cylinder_size_changed)
        

        self.box_display_cylinder_size_label.pack(side = LEFT)
        self.box_display_cylinder_size_location.pack(side = LEFT)
        self.box_display_cylinder_size_scrollbar.pack(side = LEFT)
        self.box_display_cylinder_size_frame.pack(fill='x',padx = 4, pady=1)
        

        self.box_display_wire_frame = Pmw.Group(self.box_display_radioframe,
                                                tag_pyclass = Tkinter.Radiobutton,
                                                tag_text = 'Wired Box',
                                                tag_value = 1,
                                                tag_variable = self.box_display_mode
                                                )
        self.box_display_wire_frame.pack(fill='x',expand = 1)
        self.box_display_radiogroups.append(self.box_display_wire_frame)

        self.box_display_mesh_line_width_frame = Tkinter.Frame(self.box_display_wire_frame.interior())
        self.box_display_mesh_line_width_label = Label(self.box_display_mesh_line_width_frame, text = 'Line Width:', width=10)
        self.box_display_mesh_line_width_location = Entry(self.box_display_mesh_line_width_frame,
                                                          textvariable = self.box_display_line_width,
                                                          bg='black', fg='green',width = 10)
        self.box_display_mesh_line_width_scrollbar = Scrollbar(self.box_display_mesh_line_width_frame,
                                                               orient = 'horizontal',command = self.box_display_line_width_changed)
        self.box_display_mesh_line_width_label.pack(side = LEFT)
        self.box_display_mesh_line_width_location.pack(side = LEFT)
        self.box_display_mesh_line_width_scrollbar.pack(side = LEFT)
        self.box_display_mesh_line_width_frame.pack(fill='x',padx=4,pady=1)

        
        self.box_display_mesh_grid_frame = Tkinter.Frame(self.box_display_wire_frame.interior())
        self.box_display_mesh_grid_label = Label(self.box_display_mesh_grid_frame, text = 'Grid Size:', width=10)
        self.box_display_mesh_grid_location = Entry(self.box_display_mesh_grid_frame,
                                                    textvariable = self.box_display_mesh_grid,
                                                    bg='black', fg='green',width = 10)
        self.box_display_mesh_grid_scrollbar = Scrollbar(self.box_display_mesh_grid_frame,
                                                         orient = 'horizontal',command = self.box_display_mesh_grid_changed)
        self.box_display_mesh_grid_label.pack(side = LEFT)
        self.box_display_mesh_grid_location.pack(side = LEFT)
        self.box_display_mesh_grid_scrollbar.pack(side = LEFT)
        self.box_display_mesh_grid_frame.pack(fill='x',padx=4,pady=1)

        self.box_display_radioframe.pack(padx = 6, pady = 6, expand='yes', fill='both')
        Pmw.aligngrouptags( self.box_display_radiogroups )
        
        # grid center definition
        
        self.grid_center_radiogroups = []
        self.grid_center_selection_mode = IntVar()
        self.grid_center_selection_mode.set(GRID_CENTER_FROM_SELECTION)
        self.grid_center_radioframe = Tkinter.Frame(self.grid_definition_page)

        self.grid_center_radio_button_pymol_selection = Pmw.Group(self.grid_center_radioframe,
                                                                  tag_pyclass = Tkinter.Radiobutton,
                                                                  tag_text = 'Calculate Grid Center by Selection',
                                                                  tag_value = GRID_CENTER_FROM_SELECTION,
                                                                  tag_variable = self.grid_center_selection_mode
                                                                  )
        self.grid_center_radio_button_pymol_selection.pack(fill = 'x', expand = 1, side = TOP)

        self.grid_center_radiogroups.append(self.grid_center_radio_button_pymol_selection)
        
        self.grid_center_selection_entry = Pmw.EntryField(self.grid_center_radio_button_pymol_selection.interior(),
                                                          labelpos = 'w',
                                                          label_text = 'Selection',
                                                          value = default_settings['grid_center_selection'],
                                                          command = self.grid_center_from_selection_changed
                                                          )
        self.grid_center_selection_entry.pack(fill='x',padx=4,pady=1,expand=0)

        
        self.grid_center_radio_button_coordinates = Pmw.Group(self.grid_center_radioframe,
                                                              tag_pyclass = Tkinter.Radiobutton,
                                                              tag_text = 'Grid Center Coordinates',
                                                              tag_value = GRID_CENTER_FROM_COORDINATES,
                                                              tag_variable = self.grid_center_selection_mode
                                                              )                                                              
        self.grid_center_radio_button_coordinates.pack(fill = 'x', expand = 1, side = TOP)
        
        self.grid_center_radiogroups.append(self.grid_center_radio_button_coordinates)

        self.grid_center_radioframe.pack(padx = 6, pady = 6, expand='yes', fill='both')
        Pmw.aligngrouptags(self.grid_center_radiogroups)
        
        
        self.grid_center_X_frame = Tkinter.Frame(self.grid_center_radio_button_coordinates.interior())
        self.grid_center_X_label = Label(self.grid_center_X_frame, text = 'X:')
        self.grid_center_X_location = Entry(self.grid_center_X_frame, textvariable = self.grid_center[0], bg='black', fg='green', width=10)
        self.grid_center_X_scrollbar = Scrollbar(self.grid_center_X_frame,orient='horizontal',command = self.grid_center_X_changed)

        self.grid_center_Y_frame = Tkinter.Frame(self.grid_center_radio_button_coordinates.interior())
        self.grid_center_Y_label = Label(self.grid_center_Y_frame, text = 'Y:')
        self.grid_center_Y_location = Entry(self.grid_center_Y_frame, textvariable = self.grid_center[1], bg='black', fg='green', width=10)
        self.grid_center_Y_scrollbar = Scrollbar(self.grid_center_Y_frame,orient='horizontal',command = self.grid_center_Y_changed)

        self.grid_center_Z_frame = Tkinter.Frame(self.grid_center_radio_button_coordinates.interior())
        self.grid_center_Z_label = Label(self.grid_center_Z_frame, text = 'Z:')
        self.grid_center_Z_location = Entry(self.grid_center_Z_frame, textvariable = self.grid_center[2], bg='black', fg='green', width=10)
        self.grid_center_Z_scrollbar = Scrollbar(self.grid_center_Z_frame,orient='horizontal',command = self.grid_center_Z_changed)

        self.grid_center_X_label.pack(side = LEFT)
        self.grid_center_X_location.pack(side=LEFT)
        self.grid_center_X_scrollbar.pack(side=LEFT)
        self.grid_center_X_frame.pack(side=LEFT,padx=4,pady=1)

        self.grid_center_Y_label.pack(side = LEFT)
        self.grid_center_Y_location.pack(side=LEFT)
        self.grid_center_Y_scrollbar.pack(side=LEFT)
        self.grid_center_Y_frame.pack(side=LEFT,padx=4,pady=1)

        self.grid_center_Z_label.pack(side = LEFT)
        self.grid_center_Z_location.pack(side=LEFT)
        self.grid_center_Z_scrollbar.pack(side=LEFT)
        self.grid_center_Z_frame.pack(side=LEFT,padx=4,pady=1)

        self.select_binding_site_button_box = Pmw.ButtonBox(self.grid_center_radio_button_coordinates.interior(),orient='horizontal', padx=0,pady=0)
        self.select_binding_site_button_box.add('Select binding site',command = self.select_atoms_within_binding_site)
        self.select_binding_site_button_box.pack(side=TOP,expand = 1, padx = 3, pady = 3)
        
        # load/write gpf

        self.gpf_file_io = Pmw.Group(self.grid_definition_page, tag_text='GPF File')
        self.gpf_file_io.pack(side = TOP,expand=1, fill='x')

        self.gpf_file_location = Pmw.EntryField(self.gpf_file_io.interior(),
                                                labelpos = 'w',
                                                label_pyclass = FileDialogButtonClassFactory.get(self.set_gpf_filename,mode='w',filter=[("Grid Parameter File","*.gpf")]),                                                
                                                validate = {'validator':quickFileValidation,},
                                                value = default_settings['gpf_file_name'],
                                                label_text = 'Autodock GPF File:')
        self.gpf_file_location.pack(side=LEFT,fill = 'x', expand = 1, padx = 1, pady = 5)

        self.gpf_button_box = Pmw.ButtonBox(self.gpf_file_io.interior(),orient='horizontal', padx=0,pady=0)
        self.gpf_button_box.add('Load',command = self.load_gpf_file)
        self.gpf_button_box.add('Save',command = self.save_gpf_file)
        self.gpf_button_box.pack(side=BOTTOM,expand = 1, padx = 10, pady = 5)

        # load/write vina config file
        self.config_file_io = Pmw.Group(self.grid_definition_page, tag_text='Config File')
        self.config_file_io.pack(side = TOP,expand=1, fill='x')

        self.config_file_location = Pmw.EntryField(self.config_file_io.interior(),
                                                labelpos = 'w',
                                                label_pyclass = FileDialogButtonClassFactory.get(self.set_config_filename,mode='w',filter=[("Vina Config File","*.txt")]),                                                
                                                validate = {'validator':quickFileValidation,},
                                                value = default_settings['config_file_name'],
                                                label_text = 'VINA config File:')
        self.config_file_location.pack(side=LEFT,fill = 'x', expand = 1, padx = 1, pady = 5)

        self.config_button_box = Pmw.ButtonBox(self.config_file_io.interior(), padx=0, pady=0,orient='horizontal')
        self.config_button_box.pack(side=BOTTOM,expand = 1, padx = 10, pady = 5)
        self.config_button_box.add('Load',command = self.load_config_file)
        self.config_button_box.add('Save',command = self.save_config_file)


        #------------------------------------------------------------------
        #
        #     PLUGIN CONFIGURATION PAGE

        self.configuration_top_group = Pmw.Group(self.configuration_page,tag_text='General Notes')
        self.configuration_top_group.pack(fill = 'both', expand = 0, padx = 10, pady = 5)

        self.text_field = Tkinter.Label(self.configuration_top_group.interior(),
                                         text = intro_text,
                                         background = 'black',
                                         foreground = 'yellow',
                                        justify = LEFT,
                                         )
        self.text_field.pack(expand = 0, fill = 'both', padx = 4, pady = 4)


        self.configuration_group = Pmw.Group(self.configuration_page,tag_text='Scripts and Program Paths')
        self.configuration_group.pack(fill = 'both', expand = 0, padx = 10, pady = 5)

        self.config_settings = self.read_plugin_config_file()
        
        self.autodock_tools_location = Pmw.EntryField(self.configuration_group.interior(),
                                         labelpos='w',
                                         label_pyclass = DirDialogButtonClassFactory.get(self.set_autodock_tools_path),
                                        value = self.config_settings['autodock_tools_path'],
                                         label_text = 'AutoDockTools:')


        self.autogrid_location = Pmw.EntryField(self.configuration_group.interior(),
                                         labelpos='w',
                                         label_pyclass = FileDialogButtonClassFactory.get(self.set_autogrid_location),
                                         value = self.config_settings['autogrid_exe'],
                                         label_text = 'autogrid4 executable:')

        self.autodock_location = Pmw.EntryField(self.configuration_group.interior(),
                                         labelpos='w',
                                         label_pyclass = FileDialogButtonClassFactory.get(self.set_autodock_location),
                                         value = self.config_settings['autodock_exe'],
                                         label_text = 'autodock4 executable:')

        self.vina_location = Pmw.EntryField(self.configuration_group.interior(),
                                         labelpos='w',
                                         label_pyclass = FileDialogButtonClassFactory.get(self.set_vina_location),
                                         value = self.config_settings['vina_exe'],
                                         label_text = 'vina executable:')

        self.work_path_location = Pmw.EntryField(self.configuration_group.interior(),
                                                 labelpos='w',
                                                 label_pyclass = DirDialogButtonClassFactory.get(self.set_work_path_location),
                                                 value = os.path.abspath(os.curdir),
                                                 label_text = 'Working Directory:')

        for x in  [self.autodock_tools_location,
                   self.autogrid_location,
                   self.autodock_location,
                   self.vina_location,
                   self.work_path_location
                   ]:
            x.pack(fill = 'both', expand = 1, padx = 10, pady = 5)

        Pmw.alignlabels(  [self.autodock_tools_location,
                           self.autogrid_location,
                           self.autodock_location,
                           self.vina_location,
                           self.work_path_location
                           ] )


        self.config_button_box = Pmw.ButtonBox(self.configuration_page, padx=0, pady=0,orient='horizontal')
        self.config_button_box.pack(side=BOTTOM,expand = 1, padx = 10, pady = 0)
        self.config_button_box.add('Save Plugin Configuration File',command = self.save_plugin_config_file)



        #------------------------------------------------------------------
        #
        #     RECEPTOR PREPARATION PAGE

        self.receptor_preparation_top_group = Pmw.Group(self.receptor_preparation_page, tag_text='Selections')
        self.receptor_preparation_top_group.pack(fill = 'both', expand = 0, padx=10, pady=5)
        
        self.receptor_import_selection_button_box = Pmw.ButtonBox(self.receptor_preparation_top_group.interior(),orient='vertical', padx=0,pady=0)
        self.receptor_import_selection_button_box.add('Import Selections',command = self.import_selections)
        self.receptor_import_selection_button_box.pack(side=LEFT, fill='x', padx = 0, pady = 3)


        self.receptor_pdbqt_location = Pmw.EntryField(self.receptor_preparation_top_group.interior(),
                                                labelpos = 'w',
                                                label_pyclass = FileDialogButtonClassFactory.get(self.set_receptor_pdbqt_location,filter=[("PDBQT File","*.pdbqt")]),

                                                validate = {'validator':quickFileValidation,},
                                                value = '',
                                                label_text = 'Receptor:')
        self.receptor_pdbqt_location.pack(side=LEFT,fill = 'x', expand = 1, padx = 1, pady = 5)

        self.receptor_button_box = Pmw.ButtonBox(self.receptor_preparation_top_group.interior(), padx=0, pady=0,orient='horizontal')
        self.receptor_button_box.pack(side=BOTTOM,expand = 1, padx = 10, pady = 5)
        self.receptor_button_box.add('Load',command = self.load_receptor_pdbqt)

        

        

        
        self.receptor_preparation_center_group = Pmw.Group(self.receptor_preparation_page, tag_text='Receptor Preparation')
        self.receptor_preparation_center_group.pack(fill = 'both', expand = 0, padx=10, pady=5)

        
        self.selection_list = Pmw.ComboBox(self.receptor_preparation_center_group.interior(),
                                                  scrolledlist_items=cmd.get_names("selections")+cmd.get_names(),
                                                  labelpos='nw',
                                                  label_text='PyMOL Selections',
                                                  listbox_height = 10,
                                                  selectioncommand=self.selectionCommand,
                                           dropdown=False
                                           )

        self.receptor_list = Pmw.ComboBox(self.receptor_preparation_center_group.interior(),
                                          scrolledlist_items=[],
                                          labelpos='nw',
                                          label_text='Receptors',
                                          listbox_height = 10,
                                          selectioncommand=self.selected_receptor,
                                          dropdown=False
                                          )
 
        self.flexible_residues_list = Pmw.ScrolledListBox(self.receptor_preparation_center_group.interior(),
                                        items=[],
                                        labelpos='nw',
                                        label_text='Flexible Residues',
                                        listbox_height = 11,
                                        selectioncommand=self.delete_residue,
                                        )
        
        self.selection_list.pack(side=LEFT, padx=0, anchor='n')
        self.receptor_conversion_button_box = Pmw.ButtonBox(self.receptor_preparation_center_group.interior(),orient='vertical', padx=0,pady=0)
        self.receptor_conversion_button_box.add('Generate Receptor ->',command = self.generate_receptor)
        self.receptor_conversion_button_box.add('Select as Flexible ->',command = self.select_flexible_residues)
        self.receptor_conversion_button_box.add('Remove Receptor',command = self.remove_receptor)
        self.receptor_conversion_button_box.add('Remove Flexible',command = self.remove_flexible_residues)
        self.receptor_conversion_button_box.add('Remove All',command = self.remove_all_receptors)
        self.receptor_conversion_button_box.pack(side = LEFT, expand = 0, padx = 0, pady = 12, anchor='n')
        self.receptor_list.pack(side=LEFT, padx=3, anchor='n')
        self.flexible_residues_list.pack(side=LEFT, padx=3, anchor='n')


        self.receptor_preparation_bottom_group = Pmw.Group(self.receptor_preparation_page, tag_text='Log')
        self.receptor_preparation_bottom_group.pack(fill = 'both', expand = 0, padx=10, pady=5)

        self.receptor_page_log_text = Pmw.ScrolledText(self.receptor_preparation_bottom_group.interior(),
                                                        borderframe=5, 
                                                        vscrollmode='dynamic',
                                                        hscrollmode='dynamic',
                                                        labelpos='n',
                                                        text_width=150, text_height=15,
                                                        text_wrap='none',
                                                        text_background='#000000',
                                                        text_foreground='green'
                                                        )
        self.receptor_page_log_text.pack(side=LEFT, anchor='n',pady=0)

        self.receptor_page_log_text.insert('end',receptor_prep_text)


        #===============================================
        #
        #      LIGAND PREPARATION PAGE



        self.ligand_preparation_top_group = Pmw.Group(self.ligand_preparation_page, tag_text='Selections')
        self.ligand_preparation_top_group.pack(fill = 'both', expand = 0, padx=10, pady=5)
        
        self.ligand_import_selection_button_box = Pmw.ButtonBox(self.ligand_preparation_top_group.interior(),orient='vertical', padx=0,pady=0)
        self.ligand_import_selection_button_box.add('Import Selections',command = self.import_selections)
        self.ligand_import_selection_button_box.pack(side=LEFT, fill='x', padx = 0, pady = 3)

        self.ligand_dir_location = Pmw.EntryField(self.ligand_preparation_top_group.interior(),
                                                labelpos = 'w',
                                                label_pyclass = DirDialogButtonClassFactory.get(self.set_ligand_dir_location),                                                
                                                validate = {'validator':quickFileValidation,},
                                                value = '',
                                                label_text = 'Ligands:')
        self.ligand_dir_location.pack(side=LEFT,fill = 'x', expand = 1, padx = 1, pady = 5)

        self.ligand_button_box = Pmw.ButtonBox(self.ligand_preparation_top_group.interior(), padx=0, pady=0,orient='horizontal')
        self.ligand_button_box.pack(side=BOTTOM,expand = 1, padx = 10, pady = 5)
        self.ligand_button_box.add('Import',command = self.import_ligands)

        

        
        self.ligand_preparation_center_group = Pmw.Group(self.ligand_preparation_page, tag_text='Ligand Preparation')
        self.ligand_preparation_center_group.pack(fill = 'both', expand = 0, padx=10, pady=5)

        
        self.ligand_selection_list = Pmw.ComboBox(self.ligand_preparation_center_group.interior(),
                                                  scrolledlist_items=cmd.get_names("selections")+cmd.get_names(),
                                                  labelpos='nw',
                                                  label_text='PyMOL Selections',
                                                  listbox_height = 10,
                                                  selectioncommand=self.save_as_ligand_pdb,
                                                  dropdown=False
#                                                  vscrollmode='dynamic',
#                                                  hscrollmode='dynamic',

        )

        self.ligand_list = Pmw.ComboBox(self.ligand_preparation_center_group.interior(),
                                          scrolledlist_items=[],
                                          labelpos='nw',
                                          label_text='Ligand List',
                                          listbox_height = 10,
                                          selectioncommand=self.ligand_info,
                                          dropdown=False
#                                                 vscrollmode='dynamic',
#                                                 hscrollmode='dynamic',
        )

        self.ligand_pdbqt_list = Pmw.ComboBox(self.ligand_preparation_center_group.interior(),
                                              scrolledlist_items=[],
                                              labelpos='nw',
                                              label_text='Prepared Ligands',
                                              listbox_height = 10,
                                              selectioncommand=self.ligand_info,
                                              dropdown=False,
                                                          #                                                   dropdown=False
#                                                          vscrollmode='dynamic',
#                                                          hscrollmode='dynamic',
        )
        self.ligand_selection_list.pack(side=LEFT, padx=0, anchor='n')
        self.ligand_list.pack(side=LEFT, padx=3, anchor='n')
        self.ligand_conversion_button_box = Pmw.ButtonBox(self.ligand_preparation_center_group.interior(),orient='vertical', padx=0,pady=0)
        self.ligand_conversion_button_box.add('Generate Ligand ->',command = self.generate_ligand)
        self.ligand_conversion_button_box.add('Remove Ligand',command = self.remove_ligand)
        self.ligand_conversion_button_box.add('Display Ligand',command = self.display_ligand)
        self.ligand_conversion_button_box.add('Generate All',command = self.generate_all_ligands)
        self.ligand_conversion_button_box.add('Remove All',command = self.remove_all_ligands)
        self.ligand_conversion_button_box.pack(side = LEFT, expand = 0, padx = 0, pady = 12, anchor='n')
        self.ligand_pdbqt_list.pack(side=LEFT, padx=3, anchor='n')




        self.ligand_preparation_bottom_group = Pmw.Group(self.ligand_preparation_page, tag_text='Log')
        self.ligand_preparation_bottom_group.pack(fill = 'both', expand = 0, padx=10, pady=5)

        self.ligand_page_log_text = Pmw.ScrolledText(self.ligand_preparation_bottom_group.interior(),
                                                        borderframe=5, 
                                                        vscrollmode='dynamic',
                                                        hscrollmode='dynamic',
                                                        labelpos='n',
#                                                        label_text='Log',
                                                        text_width=150, text_height=15,
                                                        text_wrap='none',
                                                        text_background='#000000',
                                                        text_foreground='green'
                                                        )
        self.ligand_page_log_text.pack(side=LEFT, anchor='n',pady=0)
        self.ligand_page_log_text.insert('end',ligand_prep_text)
        


        #------------------------------------------------------------------
        #
        #     DOCKING PAGE


        
        self.docking_top_group = Pmw.Group(self.docking_page, tag_text='Docking')
        self.docking_top_group.pack(fill = 'both', expand = 0, padx=10, pady=5)

        
        self.docking_receptor_rigid_list = Pmw.ComboBox(self.docking_top_group.interior(),
                                                  scrolledlist_items=[],
                                                  labelpos='nw',
                                                  label_text='Receptor',
                                                  listbox_height = 10,
                                                  selectioncommand=self.selectionCommand,
                                                  dropdown=True

        )
        self.docking_receptor_flexible_list = Pmw.ComboBox(self.docking_top_group.interior(),
                                                           scrolledlist_items=['No','Yes'],
                                                           labelpos='nw',
                                                           label_text='Use flexible sidechains',
                                                           listbox_height = 2,
                                                           selectioncommand=self.selectionCommand,
                                                           dropdown=True

        )

        self.docking_ligand_list = Pmw.ComboBox(self.docking_top_group.interior(),
                                                scrolledlist_items=['All'],
                                                labelpos='nw',
                                                label_text='Ligands',
                                                listbox_height = 10,
                                                selectioncommand=self.selectionCommand,
                                                dropdown=True
        )

        self.docking_nposes_list = Pmw.ComboBox(self.docking_top_group.interior(),
                                                scrolledlist_items=range(1,101),
                                                labelpos='nw',
                                                label_text='# Poses',
                                                listbox_height = 10,
                                                selectioncommand=self.selectionCommand,
                                                dropdown=True
        )
        self.docking_nposes_list.selectitem(9)
        self.docking_receptor_flexible_list.selectitem('No')
        self.docking_ligand_list.selectitem('All')
        self.docking_receptor_rigid_list.pack(side=LEFT, padx=0, anchor='n')
        self.docking_receptor_flexible_list.pack(side=LEFT, padx=0, anchor='n')
        self.docking_ligand_list.pack(side=LEFT, padx=0, anchor='n')
        self.docking_nposes_list.pack(side=LEFT, padx=0, anchor='n')





        self.docking_center_group = Pmw.Group(self.docking_page, tag_text='AutoDock')
        self.docking_center_group.pack(fill = 'both', expand = 0, padx=10, pady=5)
        self.docking_center_group2 = Pmw.Group(self.docking_page, tag_text='VINA')
        self.docking_center_group2.pack(fill = 'both', expand = 0, padx=10, pady=5)
        
        self.docking_button_box = Pmw.ButtonBox(self.docking_center_group.interior(),orient='horizontal', padx=0,pady=0)
        self.docking_button_box.add('Run AutoGrid',command = self.run_autogrid)
        self.docking_button_box.add('Run AutoDock',command = self.run_autodock)
        self.docking_button_box.add('Write AutoDock Input File(s)',command = self.write_autodock_input_files)

        self.docking_button_box2 = Pmw.ButtonBox(self.docking_center_group2.interior(),orient='horizontal', padx=0,pady=0)
        self.docking_button_box2.add('Run Vina',command = self.run_vina)
        self.docking_button_box2.add('Write Vina Input File(s)',command = self.write_vina_input_files)
        self.docking_button_box.pack(side=LEFT, fill='x', padx = 0, pady = 3)
        self.docking_button_box2.pack(side=LEFT, fill='x', padx = 0, pady = 3)
        self.docking_button_box.alignbuttons()
        self.docking_button_box2.alignbuttons()
 
        self.docking_bottom_group = Pmw.Group(self.docking_page, tag_text='Log')
        self.docking_bottom_group.pack(fill = 'both', expand = 0, padx=10, pady=5)

        self.docking_page_log_text = Pmw.ScrolledText(self.docking_bottom_group.interior(),
                                                        borderframe=5, 
                                                        vscrollmode='dynamic',
                                                        hscrollmode='dynamic',
                                                        labelpos='n',
                                                        text_width=150, text_height=25,
                                                        text_wrap='none',
                                                        text_background='#000000',
                                                        text_foreground='green'
                                                        )
        self.docking_page_log_text.pack(side=LEFT, anchor='n',pady=0)
        self.docking_page_log_text.insert('end',docking_text)
        

        #------------------------------------------------------------------
        #
        #     POSE VIEWER PAGE



        self.pose_viewer_page_top_group = Pmw.Group(self.pose_viewer_page,tag_text='File')
        self.pose_viewer_page_top_group.pack(fill = 'both', expand = 0, padx = 10, pady = 5)


        
        self.pose_viewer_page_stucts = Pmw.Group(self.pose_viewer_page,tag_text='Poses')
        self.pose_viewer_page_stucts.pack(fill = 'both', expand = 1, padx = 10, pady = 0)

        self.pose_viewer_page_display = Pmw.Group(self.pose_viewer_page,tag_text='Display Options')
        self.pose_viewer_page_display.pack(fill = 'x', expand = 1, padx = 10, pady = 0)

        self.pose_viewer_notebook = Pmw.NoteBook(self.pose_viewer_page_stucts.interior())
        self.pose_viewer_notebook.pack(fill='both',expand=1,padx=3,pady=3)
        self.pose_viewer_pages = {}
        self.pose_viewer_ligand_dic = {}        
        self.pose_file = StringVar()
        self.pose_file.set(default_settings['dlg_input_file'])
        self.pose_file_location = Pmw.EntryField(self.pose_viewer_page_top_group.interior(),
                                                 labelpos='w',
                                                 label_pyclass = FileDialogButtonClassFactory.get(self.set_pose_filename,filter=[("PDBQT File","*.pdbqt"),("DLG File","*.dlg")]),
                                                 validate = {'validator':quickFileValidation,},
                                                 value = default_settings['dlg_input_file'],
                                                 label_text = 'Browse:')

#        self.pose_file_location.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        self.pose_file_location.pack(side=LEFT,fill = 'x', expand = 1, padx = 1, pady = 5)
        
        self.load_pose_file_buttonbox = Pmw.ButtonBox(self.pose_viewer_page_top_group.interior(), padx=0)
        self.load_pose_file_buttonbox.pack(side=BOTTOM,expand = 1, padx = 10, pady = 5)
        self.load_pose_file_buttonbox.add('Load',command=self.load_ligand_file)
        self.load_pose_file_buttonbox.add('Load All',command=self.load_all_ligand_files)


        self.pose_viewer_ligand_display_radio = Pmw.RadioSelect(self.pose_viewer_page_display.interior(),
                                     selectmode='multiple',
                                     buttontype='checkbutton',
                                     labelpos='w',
                                     label_text='Display mode',
                                     orient='horizontal',
                                     frame_relief='ridge',
                                     command=self.ligand_display_mode_changed)
        self.pose_viewer_ligand_display_radio.pack(side=TOP, padx=10, anchor='w')
        for entry in ('lines', 'sticks','spheres','surface','mesh'):
            self.pose_viewer_ligand_display_radio.add(entry)
        for entry in ('lines', 'sticks','spheres','surface','mesh'):
            if self.ligand_display_mode[entry]:
                self.pose_viewer_ligand_display_radio.invoke(entry)



        self.pose_viewer_radiobuttons = Pmw.RadioSelect(self.pose_viewer_page_display.interior(),
                                                        buttontype = 'radiobutton',
                                                        orient = 'horizontal',
                                                        labelpos = 'w',
                                                        )
        for text in ('Show Selected',
                     'Hide Selected'):
                self.pose_viewer_radiobuttons.add(text)
                self.pose_viewer_radiobuttons.setvalue('Show Selected')
        self.pose_viewer_radiobuttons.pack(padx=4,pady=1,side=TOP)

        self.pose_viewer_ligand_pages = {}
        
        #---------------------------------------------------------------
        # SCORE/RANK PAGE
        
        self.score_table = ScoreTable(self.rank_page)
        self.score_table.pack(pady=20)

        self.score_table_radiobuttons = Pmw.RadioSelect(self.rank_page,
                                                        buttontype = 'radiobutton',
                                                        orient = 'horizontal',
                                                        labelpos = 'w',
                                                        command = self.update_score_table
                                                        )
        for text in ('Show All Poses',
                     'Show Only Best Pose'):
                self.score_table_radiobuttons.add(text)
                self.score_table_radiobuttons.setvalue('Show All Poses')
        self.score_table_radiobuttons.pack(padx=4,pady=1,side=TOP)

        self.rank_dat_file_io = Pmw.Group(self.rank_page, tag_text='Export scores as data file')
        self.rank_dat_file_io.pack(side = TOP,expand=1, fill='x')

        self.rank_dat_file_location = Pmw.EntryField(self.rank_dat_file_io.interior(),
                                                labelpos = 'w',
                                                label_pyclass = FileDialogButtonClassFactory.get(self.set_rank_dat_filename,mode='w',filter=[("Data File","*.dat")]),                                                
                                                validate = {'validator':quickFileValidation,},
                                                value = default_settings['rank_dat_file_name'],
                                                label_text = 'Filename:')
        self.rank_dat_file_location.pack(side=LEFT,fill = 'x', expand = 1, padx = 1, pady = 0)

        self.rank_dat_button_box = Pmw.ButtonBox(self.rank_dat_file_io.interior(), padx=0, pady=0,orient='horizontal')
        self.rank_dat_button_box.pack(side=BOTTOM,expand = 1, padx = 10, pady = 0)
        self.rank_dat_button_box.add('Export',command = self.export_score_dat_file)

        self.rank_csv_file_io = Pmw.Group(self.rank_page, tag_text='Export scores as CSV file')
        self.rank_csv_file_io.pack(side = TOP,expand=1, fill='x')

        self.rank_pose_file_io = Pmw.Group(self.rank_page, tag_text='Export poses as PDB file')
        self.rank_pose_file_io.pack(side = TOP,expand=1, fill='x')


        self.rank_csv_file_location = Pmw.EntryField(self.rank_csv_file_io.interior(),
                                                labelpos = 'w',
                                                label_pyclass = FileDialogButtonClassFactory.get(self.set_rank_csv_filename,mode='w',filter=[("CSV File","*.csv")]),                                                
                                                validate = {'validator':quickFileValidation,},
                                                value = default_settings['rank_csv_file_name'],
                                                label_text = 'Filename:')
        self.rank_csv_file_location.pack(side=LEFT,fill = 'x', expand = 1, padx = 1, pady = 0)

        self.rank_csv_button_box = Pmw.ButtonBox(self.rank_csv_file_io.interior(), padx=0, pady=0,orient='horizontal')
        self.rank_csv_button_box.pack(side=BOTTOM,expand = 1, padx = 10, pady = 0)
        self.rank_csv_button_box.add('Export',command = self.export_score_csv_file)


        self.rank_pose_file_location = Pmw.EntryField(self.rank_pose_file_io.interior(),
                                                labelpos = 'w',
                                                label_pyclass = FileDialogButtonClassFactory.get(self.set_rank_pose_filename,mode='w',filter=[("PDB File","*.pdb")]),                                                
                                                validate = {'validator':quickFileValidation,},
                                                value = default_settings['rank_pose_file_name'],
                                                label_text = 'Filename:')
        self.rank_pose_file_location.pack(side=LEFT,fill = 'x', expand = 1, padx = 1, pady = 0)

        self.rank_pose_button_box = Pmw.ButtonBox(self.rank_pose_file_io.interior(), padx=0, pady=0,orient='horizontal')
        self.rank_pose_button_box.pack(side=BOTTOM,expand = 1, padx = 10, pady = 0)
        self.rank_pose_button_box.add('Export',command = self.export_score_pose_file)

        #------------------------------------------------------------------
        #
        #     MAP VIEWER PAGE


        self.map_threshold = {}
        self.map_meta = {}

        
        self.map_viewer_page_top_group = Pmw.Group(self.map_viewer_page,tag_text='Grid Map')
        self.map_viewer_page_top_group.pack(fill = 'both', expand = 0, padx = 10, pady = 5)

        # the maps card
        
        self.map_viewer_page_center_group = Pmw.Group(self.map_viewer_page,tag_text='Maps')
        self.map_viewer_page_center_group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        self.map_viewer_notebook = Pmw.NoteBook(self.map_viewer_page_center_group.interior())
        self.map_viewer_notebook.pack(fill='both',expand=1,padx=3,pady=3)
        self.map_pages = {}
        self.map_dic = {}        
        self.mapfile = StringVar()
        self.mapfile.set(default_settings['map_input_file'])
        self.map_file_location = Pmw.EntryField(self.map_viewer_page_top_group.interior(),
                                          labelpos='w',
                                          label_pyclass = FileDialogButtonClassFactory.get(self.set_mapfilename,filter=[("Autodock Map File","*.map")]),
                                          validate = {'validator':quickFileValidation,},
                                          value = default_settings['map_input_file'],
                                          label_text = 'Browse:')

        self.map_file_location.pack(side = LEFT,fill = 'x', expand = 1, padx = 10, pady = 5)

        self.load_map_buttonbox = Pmw.ButtonBox(self.map_viewer_page_top_group.interior(), padx=0)
        self.load_map_buttonbox.pack(side=BOTTOM,expand = 1, padx = 10, pady = 5)
        self.load_map_buttonbox.add('Load',command=self.load_grid_map)




        #------------------------------------------------------------------
        ##################################################################
        # DONE PAGES
        self.notebook.setnaturalsize()
        self.dialog.show()
        self.status_line.configure(text ="Ready..... (version: %s)" % __version__)
        #------------------------------------------------------------------
        ##################################################################




    def button_pressed(self, result):
        if hasattr(result,'keycode'):
            if result.keycode == 36:
                if self.notebook.getcurselection()=='Grid Settings':
                    self.show_box()
                elif self.notebook.getcurselection()=='View Poses':
                    self.load_ligand_file()
        elif result == 'Exit' or result == None:
            self.dialog.withdraw()

    #------------------------------------------------------------------
    # grid settings functions
    
    def grid_spacing_changed(self, x):
        val =  float(self.grid_spacing.get())+float(x)*0.005
        self.grid_spacing.set(val)
        if self.box_is_on_display:
            self.grid_center_selection_mode.set(GRID_CENTER_FROM_COORDINATES)
            
    def calculate_grid_center(self):
        if self.grid_center_selection_mode.get() == GRID_CENTER_FROM_SELECTION:
            sel = self.grid_center_selection_entry.get()
            if sel:
                stored.xyz = []
                cmd.iterate_state(1,sel,"stored.xyz.append([x,y,z])")
                xx = average(map(lambda a: a[0], stored.xyz))
                yy = average(map(lambda a: a[1], stored.xyz))
                zz = average(map(lambda a: a[2], stored.xyz))
                self.grid_center[0].set(round(xx,2))
                self.grid_center[1].set(round(yy,2))
                self.grid_center[2].set(round(zz,2))
            else:
                self.grid_center_selection_mode.set(GRID_CENTER_FROM_COORDINATES)
                    
    
    def n_points_X_changed(self, x):
        val = int(self.n_points_X.get())+int(x)
        self.n_points_X.set(val)
        self.grid_center_selection_mode.set(GRID_CENTER_FROM_COORDINATES)
        self.show_box()
        

    def n_points_Y_changed(self, x):
        val = int(self.n_points_Y.get())+int(x)
        self.n_points_Y.set(val)
        self.grid_center_selection_mode.set(GRID_CENTER_FROM_COORDINATES)
        self.show_box()

    def n_points_Z_changed(self, x):
        val = int(self.n_points_Z.get())+int(x)
        self.n_points_Z.set(val)
        self.grid_center_selection_mode.set(GRID_CENTER_FROM_COORDINATES)
        self.show_box()

    def grid_center_from_selection_changed(self):
        self.grid_center_selection_mode.set(GRID_CENTER_FROM_SELECTION)
        self.show_box()
        self.grid_center_selection_entry.clear()

    def grid_center_X_changed(self, x):
        self.grid_center_selection_mode.set(GRID_CENTER_FROM_COORDINATES)
        val=float(self.grid_center[0].get())+float(x)*1.0
        self.grid_center[0].set(val)
        self.show_box()

        
    def grid_center_Y_changed(self, x):
        self.grid_center_selection_mode.set(GRID_CENTER_FROM_COORDINATES)
        val=float(self.grid_center[1].get())+float(x)*1.0
        self.grid_center[1].set(val)
        self.show_box()

    def grid_center_Z_changed(self, x):
        self.grid_center_selection_mode.set(GRID_CENTER_FROM_COORDINATES)
        val=float(self.grid_center[2].get())+float(x)*1.0
        self.grid_center[2].set(val)
        self.show_box()

    
    def select_atoms_within_binding_site(self):
        m = cmd.get_model("polymer")
        xmin, xmax = self.box_coords[0]
        ymin, ymax = self.box_coords[1]
        zmin, zmax = self.box_coords[2]
        lst = filter(lambda a: a.coord[0] >= xmin and \
                     a.coord[0] <= xmax and \
                     a.coord[1] >= ymin and \
                     a.coord[1] <= ymax and \
                     a.coord[2] >= zmin and \
                     a.coord[2] <= zmax, m.atom)
        by_id = map(lambda a: a.id, lst)
        if len(by_id) > 1:
            cmd.select("binding_site", "ID %d" % by_id[0])
            for idx in by_id[1:]:
                cmd.select("binding_site", "binding_site or ID %d" % idx)
            self.status_line.configure(text = "Selector 'binding_site' created with %d atoms" % len(by_id))
            self.import_selections()
        
    def box_display_cylinder_size_changed(self, x):
        val=float(self.box_display_cylinder_size.get())+float(x)*0.1
        self.box_display_cylinder_size.set(val)
        self.grid_center_selection_mode.set(GRID_CENTER_FROM_COORDINATES)
        self.show_box()

    def box_display_line_width_changed(self, x):
        val=float(self.box_display_line_width.get())+float(x)*0.1
        self.box_display_line_width.set(val)
        self.grid_center_selection_mode.set(GRID_CENTER_FROM_COORDINATES)
        self.show_box()


    def box_display_mesh_grid_changed(self, x):
        val=float(self.box_display_mesh_grid.get())+float(x)*0.1
        self.box_display_mesh_grid.set(val)
        self.grid_center_selection_mode.set(GRID_CENTER_FROM_COORDINATES)
        self.show_box()


    def show_box(self):
        self.calculate_grid_center()
        self.show_crisscross()        
        self.calculate_box()

    def hide_box(self):
        cmd.delete("box")
        cmd.delete("wirebox")
        cmd.delete("grid_center")
        self.box_is_on_display = False
        
    def change_box_color(self):
        color = tkColorChooser.Chooser(
            initialcolor='white',title='Choose box color').show()
        if color[0] is not None:
            self.box_color = [color[0][0]/100.,
                              color[0][1]/100.,
                              color[0][2]/100.]
        self.show_box()
        
    def set_box_filename(self, filename):
        self.box_file_location.setvalue(filename)

    def set_gpf_filename(self, filename):
        self.gpf_file_location.setvalue(filename)

    def set_config_filename(self, filename):
        self.config_file_location.setvalue(filename)

    def set_rank_dat_filename(self, filename):
        self.rank_dat_file_location.setvalue(filename)

    def set_rank_csv_filename(self, filename):
        self.rank_csv_file_location.setvalue(filename)

    def set_rank_pose_filename(self, filename):
        self.rank_pose_file_location.setvalue(filename)

    def load_gpf_file(self):
        filename = self.gpf_file_location.get()
        fp = self.fileopen(filename,'r')
        if not fp:
            return
        lst = fp.readlines()
        new = []
        for line in lst:
            if line.strip():
                new.append(line.strip())
        lst = new
        for line in lst:
            entr = line.split()
            if entr[0] == 'npts':
                n_points_X = int(entr[1])
                n_points_Y = int(entr[2])
                n_points_Z = int(entr[3])
                self.n_points_X.set(n_points_X)
                self.n_points_Y.set(n_points_Y)
                self.n_points_Z.set(n_points_Z)
            elif entr[0] == 'spacing':
                spacing = float(entr[1])
                self.grid_spacing.set(spacing)
            elif entr[0] == 'gridcenter':
                if entr[1]!='auto':
                    grid_X = float(entr[1])
                    grid_Y = float(entr[2])
                    grid_Z = float(entr[3])
                    self.grid_center[0].set(grid_X)
                    self.grid_center[1].set(grid_Y)
                    self.grid_center[2].set(grid_Z)
        self.status_line.configure(text= 'Reading box info from %s' % filename)
        self.grid_center_selection_mode.set(GRID_CENTER_FROM_COORDINATES)
        self.calculate_box()
            
    def save_gpf_file(self):
        filename = self.gpf_file_location.get()
        fp = self.fileopen(filename,'w')
        if not fp:
                return
        n_points_X = self.n_points_X.get()
        n_points_Y = self.n_points_Y.get()
        n_points_Z = self.n_points_Z.get()
        spacing = self.grid_spacing.get()
        center_X = self.grid_center[0].get()
        center_Y = self.grid_center[1].get()
        center_Z = self.grid_center[2].get()
        print >>fp, 'npts %d %d %d' % (n_points_X, n_points_Y, n_points_Z)
        print >>fp, 'spacing %5.3f' % spacing
        print >>fp, 'gridcenter  %8.3f %8.3f %8.3f' % (center_X, center_Y, center_Z)
        fp.close()
        self.status_line.configure(text= 'Wrote box info to %s' % filename)


            
    def load_config_file(self):
        filename = self.config_file_location.get()
        fp = self.fileopen(filename,'r')
        spacing = self.grid_spacing.get()
        if not fp:
            return
        lst = fp.readlines()
        new = []
        for line in lst:
            if line.strip():
                new.append(line.strip())
        lst = new
        for line in lst:
            entr = line.split()
            if entr[0] == 'center_x':
                center_x = float(entr[2])
                self.grid_center[0].set(center_x)
            elif entr[0] == 'center_y':
                center_y = float(entr[2])
                self.grid_center[1].set(center_y)
            elif entr[0] == 'center_z':
                center_z = float(entr[2])
                self.grid_center[2].set(center_z)
            elif entr[0] == 'size_x':
                size_x = float(entr[2])
                n_points_X = size_x/spacing
                self.n_points_X.set(n_points_X)
            elif entr[0] == 'size_y':
                size_y = float(entr[2])
                n_points_Y = size_y/spacing
                self.n_points_Y.set(n_points_Y)
            elif entr[0] == 'size_z':
                size_z = float(entr[2])
                n_points_Z = size_z/spacing
                self.n_points_Z.set(n_points_Z)
        self.status_line.configure(text= 'Reading box info from %s' % filename)
        self.grid_center_selection_mode.set(GRID_CENTER_FROM_COORDINATES)
        self.calculate_box()
                            



                            
    def save_config_file(self):
        filename = self.config_file_location.get()
        fp = self.fileopen(filename,'w')
        if not fp:
            return
        n_points_X = self.n_points_X.get()
        n_points_Y = self.n_points_Y.get()
        n_points_Z = self.n_points_Z.get()
        spacing = self.grid_spacing.get()
        size_x = n_points_X*spacing
        size_y = n_points_Y*spacing
        size_z = n_points_Z*spacing
        center_x = self.grid_center[0].get()
        center_y = self.grid_center[1].get()
        center_z = self.grid_center[2].get()
        print >>fp, "size_x = %6.2f" % size_x
        print >>fp, "size_y = %6.2f" % size_y
        print >>fp, "size_z = %6.2f" % size_z
        print >>fp, "center_x = %6.2f" % center_x
        print >>fp, "center_y = %6.2f" % center_y
        print >>fp, "center_z = %6.2f" % center_z
        fp.close()
        self.status_line.configure(text= 'Wrote box info to %s' % filename)
        
    def show_crisscross(self):
        center = [float(self.grid_center[0].get()),
                  float(self.grid_center[1].get()),
                  float(self.grid_center[2].get())
                  ]
        cmd.delete("grid_center")
        self.crisscross(center[0], center[1], center[2], 0.5, "grid_center")
        
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

    def calculate_box(self):
        x = float(self.grid_center[0].get())
        y = float(self.grid_center[1].get())
        z = float(self.grid_center[2].get())
        xpts = int(self.n_points_X.get())
        ypts = int(self.n_points_Y.get())
        zpts = int(self.n_points_Z.get())
        spacing = float(self.grid_spacing.get())
        cylinder_size = float(self.box_display_cylinder_size.get())

        size = [xpts*spacing, ypts*spacing, zpts*spacing]
        xmax = x + size[0]/2.
        xmin = x - size[0]/2.
        ymax = y + size[1]/2.
        ymin = y - size[1]/2.
        zmax = z + size[2]/2.
        zmin = z - size[2]/2.
        box_edge_x = [xmin,xmax]
        box_edge_y = [ymin,ymax]
        box_edge_z = [zmin,zmax]
        self.box_coords  = [box_edge_x,box_edge_y,box_edge_z]
        cmd.delete('box')
        if self.box_display_mode.get()==BOX_AS_BOX:
            self.display_box(self.box_coords,cylinder_size)
        elif self.box_display_mode.get()==BOX_AS_WIREBOX:
            self.display_wire_box(self.box_coords)
        self.box_is_on_display = True

    def display_box(self, box, cylinder_size):
        view = cmd.get_view()
        name = "box"
        obj = []
        # build cgo object
        color = self.box_color
        for i in range(2):
            for k in range (2):
                for j in range(2):
                    if i != 1:
                        obj.append(CYLINDER)
                        obj.extend([box[0][i],box[1][j],box[2][k]])
                        obj.extend([box[0][i+1],box[1][j],box[2][k]])
                        obj.append(cylinder_size)
                        obj.extend(color)
                        obj.extend(color)
                        obj.append(COLOR)
                        obj.extend(color)
                        obj.append(SPHERE)
                        obj.extend([box[0][i],box[1][j],box[2][k],cylinder_size])
                        
                    if j != 1:
                        obj.append(CYLINDER)
                        obj.extend([box[0][i],box[1][j],box[2][k]])
                        obj.extend([box[0][i],box[1][j+1],box[2][k]])
                        obj.append(cylinder_size)
                        obj.extend(color)
                        obj.extend(color)
                        obj.append(COLOR)
                        obj.extend(color)
                        obj.append(SPHERE)
                        obj.extend([box[0][i],box[1][j+1],box[2][k],cylinder_size])
                    if k != 1:
                        obj.append(CYLINDER)
                        obj.extend([box[0][i],box[1][j],box[2][k]])
                        obj.extend([box[0][i],box[1][j],box[2][k+1]])
                        obj.append(cylinder_size)
                        obj.extend(color)
                        obj.extend(color)
                        obj.append(COLOR)
                        obj.extend(color)
                        obj.append(SPHERE)
                        obj.extend([box[0][i],box[1][j],box[2][k+1],cylinder_size])
        axes = [[2.0,0.0,0.0],[0.0,2.0,0.0],[0.0,0.0,2.0]]
        xpos = [box[0][1]+(box[0][1]-box[0][0])/5.,box[1][0],box[2][0]]
        cyl_text(obj,plain,xpos,'X',0.10,axes=axes)
        ypos = [box[0][0],box[1][1]+(box[1][1]-box[1][0])/5,box[2][0]]
        cyl_text(obj,plain,ypos,'Y',0.10,axes=axes)
        zpos = [box[0][0],box[1][0],box[2][1]+(box[2][1]-box[2][0])/5]
        cyl_text(obj,plain,zpos,'Z',0.10,axes=axes)
        cmd.load_cgo(obj,name)
        cmd.set_view(view)
        

    def display_wire_box(self, box):

        cmd.delete("wirebox")
        color = self.box_color
        view = cmd.get_view()
        spacing = float(self.box_display_mesh_grid.get())
        lwidth = float(self.box_display_line_width.get())
        xpts = int(round((box[0][1]-box[0][0])/spacing))+1
        ypts = int(round((box[1][1]-box[1][0])/spacing))+1
        zpts = int(round((box[2][1]-box[2][0])/spacing))+1
        obj = []
        for i in range(xpts):
            for k in range (ypts):
                obj.append(BEGIN)
                obj.append(LINE_STRIP)
                obj.append(COLOR)
                obj.extend(color)
                
                for j in range(zpts):
                    
                    obj.append(VERTEX)
                    obj.extend([box[0][0]+spacing*i,box[1][0]+spacing*k,\
                                box[2][0]+spacing*j])

                obj.append(END)
        for i in range(xpts):
            for j in range (zpts):
                obj.append(BEGIN)
                obj.append(LINE_STRIP)
                obj.append(COLOR)
                obj.extend(color)
                for k in range(ypts):
                    obj.append(VERTEX)
                    obj.extend([box[0][0]+spacing*i,box[1][0]+spacing*k,\
                                box[2][0]+spacing*j])
                obj.append(END)
        for j in range(zpts):
            for i in range (xpts):
                obj.append(BEGIN)
                obj.append(LINE_STRIP)
                obj.append(COLOR)
                obj.extend(color)
                for k in range(ypts):
                    obj.append(VERTEX)
                    obj.extend([box[0][0]+spacing*i,box[1][0]+spacing*k,\
                                box[2][0]+spacing*j])
                obj.append(END)
        for j in range(zpts):
            for k in range (ypts):
                obj.append(BEGIN)
                obj.append(LINE_STRIP)
                obj.append(COLOR)
                obj.extend(color)
                for i in range(xpts):
                    obj.append(VERTEX)
                    obj.extend([box[0][0]+spacing*i,box[1][0]+spacing*k,\
                                box[2][0]+spacing*j])
                obj.append(END)
        cmd.load_cgo(obj,"wirebox")
        cmd.set("cgo_line_width",lwidth)
        cmd.set_view(view)
        
    #---------------------------------------------------------------------
    # config functions
    def set_autodock_tools_path(self, dirname):
        self.autodock_tools_location.setvalue(dirname)
        self.autodock_tools_path.set(dirname)
        self.config_settings['autodock_tools_path'] = dirname
        
    def set_autogrid_location(self, filename):
        self.autogrid_location.setvalue(filename)
        self.autogrid_exe.set(filename)
        self.config_settings['autogrid_exe'] = filename
        
    def set_autodock_location(self, filename):
        self.autodock_location.setvalue(filename)
        self.autodock_exe.set(filename)
        self.config_settings['autodock_exe'] = filename
        
    def set_vina_location(self, filename):
        self.vina_location.setvalue(filename)
        self.vina_exe.set(filename)
        self.config_settings['vina_exe'] = filename

    def set_work_path_location(self, dirname):
        self.work_path_location.setvalue(dirname)

    def work_dir(self):
        return self.work_path_location.getvalue()
    
    def read_plugin_config_file(self):
        config_file_name = os.path.join(tmp_dir,"pymol_autodock_plugin.conf")
        self.config_settings = {}
        self.config_settings['autodock_tools_path'] = ''
        self.config_settings['autogrid_exe'] = ''
        self.config_settings['autodock_exe'] = ''
        self.config_settings['vina_exe'] = ''
        if os.path.isfile(config_file_name):
            self.status_line.configure(text = 'Reading configuration file: %s' % config_file_name)
            lst = self.fileopen(config_file_name,'r').readlines()
            for line in lst:
                if line[0]!='#':
                    entr = line.split('=')
                    self.config_settings[entr[0].strip()] = entr[1].strip()
            self.autogrid_exe.set(self.config_settings['autogrid_exe'])
            self.autodock_exe.set(self.config_settings['autodock_exe'])
            self.vina_exe.set(self.config_settings['vina_exe'])
            self.autodock_tools_path.set(self.config_settings['autodock_tools_path'])
        else:
            self.status_line.configure(text = 'Configuration file not found')
        return self.config_settings
    
    def save_plugin_config_file(self):
        config_file_name = os.path.join(tmp_dir,"pymol_autodock_plugin.conf")
        fp = self.fileopen(config_file_name,'w')
        print >>fp, '#========================================'
        print >>fp, '# Autodock/Vina Plugin configuration file'
        self.config_settings['autogrid_exe'] = self.autogrid_location.getvalue()
        self.config_settings['autodock_exe'] = self.autodock_location.getvalue()
        self.config_settings['vina_exe'] = self.vina_location.getvalue()
        self.config_settings['autodock_tools_path'] = self.autodock_tools_location.getvalue()
        print 'ADDD', self.autodock_location.getvalue()
        for key, val in self.config_settings.items():
            print >>fp, key, '=', val
        fp.close()
        self.status_line.configure(text = 'Wrote configuration file %s' % config_file_name)

    def ligand_display_mode_changed(self, button_name, pressed):
        if pressed:
            self.ligand_display_mode[button_name] = True
            action = 'Enabled'
        else:
            self.ligand_display_mode[button_name] = False
            action = 'Disabled'
        txt = action+' ligand display mode <'+button_name+'>'
        self.status_line.configure(text=txt)

    #----------------------------------------------------------------------------
    # receptors

    def delete_residue(self):
        sel = self.flexible_residues_list.getcurselection()

    
    def selectionCommand(self, value):
        sels = self.selection_list.getcurselection()
        ## if len(sels) == 0:
##             print 'No selection'
##         else:
##             print 'Selection:', sels[0]

    def selected_receptor(self, value):
        sel = self.receptor_list.get()
        receptor_object = self.receptor_dic[sel]
        lst = []
        for key, val in receptor_object.flexible_residues.items():
            for resn, resi in val:
                lst.append("%s:%4d %s" % (key, resn, resi))
        self.flexible_residues_list.setlist(lst)
        self.receptor_page_log_text.insert('end',receptor_object.info())
        self.receptor_page_log_text.yview('moveto',1.0)
        self.docking_receptor_rigid_list.selectitem(value)



    def import_selections(self):
        lst = cmd.get_names("selections")+cmd.get_names()
        if 'grid_center' in lst:
            lst.remove('grid_center')
        if 'box' in lst:
            lst.remove('box')
        self.selection_list.setlist(lst)
        self.ligand_selection_list.setlist(lst)
        

    def generate_receptor(self):
        print self.work_dir()
        sel = self.selection_list.getcurselection()
        tmp_rec_pdb = os.path.join(self.work_dir(),"receptor.%s.pdb" % sel[0])
        print tmp_rec_pdb
        cmd.save(tmp_rec_pdb,sel[0])
        util_program = os.path.join(self.autodock_tools_path.get(),"prepare_receptor4.py")
        outfilename = os.path.join(self.work_dir(), "receptor.%s.pdbqt" % sel[0])
        command = "%s -r %s -o %s -A checkhydrogens" % (util_program,tmp_rec_pdb,outfilename)
        self.receptor_page_log_text.insert('end',"Batch: %s\n" % command)
        self.receptor_page_log_text.yview('moveto',1.0)

        result, output = getstatusoutput(command)
        if result == 0:
            self.receptor_list.insert('end',sel[0])
            self.status_line.configure(text="Successfully generated receptor file receptor.%s.pdbqt" % sel[0])
            #self.receptor_page_log_text.insert('end',"Successfully generated receptor file receptor.%s.pdbqt\n" % sel[0])
            self.receptor_page_log_text.insert('end',output)
            self.receptor_page_log_text.yview('moveto',1.0)
            r = Receptor()
            r.selection = sel[0]
            r.pdb_file = tmp_rec_pdb
            r.receptor_pdbqt = outfilename
            stored.list = []
            cmd.iterate(sel[0]+' and name ca',"stored.list.append([resi,resn])")
            for resi, resn in stored.list:
                r.resi_dic[int(resi)] = resn
            self.receptor_dic[sel[0]] = r
            self.receptor_list.selectitem(-1)
            self.docking_receptor_rigid_list.insert('end',sel[0])
            self.docking_receptor_rigid_list.selectitem(-1)
            
        else:
            self.status_line.configure(text="An error occured while preparing receptor from selection %s" % sel[0])
            self.receptor_page_log_text.insert('end',output)
            
    def select_flexible_residues(self):
        sel = self.selection_list.get()
        rec = self.receptor_list.get()
        stored.list = []
        cmd.iterate(sel+' and name ca',"stored.list.append([chain, resn,resi])")
        receptor_object = self.receptor_dic[rec]
        receptor_object.flexible_residues = []
        chains = {}
        for chain, resn, resi in stored.list:
            if resn not in ['ALA','GLY','PRO']:
                if chains.has_key(chain):
                    chains[chain].append([int(resi), resn])
                else:
                    chains[chain] = [ [int(resi), resn] ]
        receptor_object.flexible_residues = chains

        #                receptor_object.flexible_residues.append( [int(resi), resn] )
#        self.status_line.configure(text = "Selected %d flexible residues for receptor %s" \
#                                   % (len(receptor_object.flexible_residues), rec))
        util_program = os.path.join(self.autodock_tools_path.get(),"prepare_flexreceptor4.py")
        flex_res_string = receptor_object.flex_res_string()
        receptor_filename = receptor_object.receptor_pdbqt
        rec_rigid = receptor_filename[:-6]+'.rigid.pdbqt'
        rec_flexible = receptor_filename[:-6]+'.flexible.pdbqt'
        command = "%s -r %s -s %s -g %s -x %s" % \
                  (util_program, receptor_filename, flex_res_string, rec_rigid, rec_flexible)
        self.receptor_page_log_text.insert('end',"Batch: %s\n" % command)
        self.receptor_page_log_text.yview('moveto',1.0)
#        result = os.system(command)
        result, output = getstatusoutput(command)
        if result == 0:
            receptor_object.receptor_rigid = rec_rigid
            receptor_object.receptor_flexible = rec_flexible
            self.status_line.configure(text="Successfully generated flexible receptor files for %s(%s)" % (rec,sel))
#            self.receptor_page_log_text.insert('end',"Successfully generated flexible receptor files for %s(%s)\n" % (rec,sel))
            self.receptor_page_log_text.insert('end',output)
            self.selected_receptor(rec)
            self.docking_receptor_flexible_list.selectitem('Yes')
        else:
            self.status_line.configure(text="An error occured while preparing a flexible receptor from selection %s" % sel)
            self.receptor_page_log_text.insert('end',output)
        self.receptor_page_log_text.yview('moveto',1.0)
           

    def set_receptor_pdbqt_location(self, filename):
        self.receptor_pdbqt_location.setvalue(filename)

    def load_receptor_pdbqt(self):
        fn = self.receptor_pdbqt_location.getvalue()
        outfile = os.path.join(self.work_dir(), os.path.basename(fn).split('.')[0]+'_pdb.pdb')
        util_program = os.path.join(self.autodock_tools_path.get(),"pdbqt_to_pdb.py")
        command = "%s -f %s -o %s " % (util_program, fn, outfile)
        self.receptor_page_log_text.insert('end',"Batch: %s\n" % command)
        self.receptor_page_log_text.yview('moveto',1.0)
        result, output = getstatusoutput(command)
        if result == 0:
            self.status_line.configure(text = "Loading receptor %s" % fn)
            self.receptor_page_log_text.insert('end',output)
            self.receptor_page_log_text.yview('moveto',1.0)

            name = os.path.basename(fn).split('.')[0]
            cmd.load(outfile,name)
            r = Receptor()
            r.selection = name
            r.pdb_file = outfile
            r.receptor_pdbqt = os.path.abspath(fn)
            stored.list = []
            cmd.iterate(name+' and name ca',"stored.list.append([resi,resn])")
            for resi, resn in stored.list:
                r.resi_dic[int(resi)] = resn
            self.receptor_dic[name] = r
            self.receptor_list.insert('end',name)
            self.receptor_list.selectitem(-1)
            self.docking_receptor_rigid_list.insert('end',name)
            self.docking_receptor_rigid_list.selectitem(-1)
            
        else:
            self.status_line.configure(text="An error occured while loading receptor  %s" % fn)
            self.receptor_page_log_text.insert('end',output)
            

        
    def remove_receptor(self):
        rec = self.receptor_list.get()
        del self.receptor_dic[rec]
        self.status_line.configure(text="Removed receptor %s" % rec)
        index = list(self.receptor_list.get(0, 'end')).index(rec)
        self.receptor_list.delete(index)
        self.docking_receptor_list.delete(index)
        try:
            self.receptor_list.selectitem(0)
            self.docking_receptor_list.selectitem(0)
        except:
            self.receptor_list.clear()
            self.docking_receptor_list.clear()
            
    def remove_flexible_residues(self):
        rec = self.receptor_list.get()
        rec_object = self.receptor_dic[rec]
        rec_object.flexible_residues = []
        rec_object.receptor_rigid = ""
        rec_object.receptor_flexible = ""
        self.status_line.configure(text="Removed flexible residues from receptor %s" % rec)
        self.selected_receptor(rec)
        
    def remove_all_receptors(self):
        self.receptor_dic = {}
        self.status_line.configure(text="Deleted all receptor objects")
        self.receptor_list.clear()
        self.flexible_residues_list.clear()
    
    #------------------------------------------------------------------------------
    # ligands
    
    def save_as_ligand_pdb(self, name):
        pdb_name = os.path.join(self.work_dir(), name+'.ligand.pdb')
        cmd.save(pdb_name, name)
        self.status_line.configure(text="Saving ligand pdb file %s from selection %s" % (pdb_name, name))
        self.ligand_list.insert('end',name)
        l = Ligand()
        l.input_file = pdb_name
        l.selection = name
        l.name = name
        self.ligand_dic[name] = l
        self.ligand_list.selectitem(-1)

    def generate_ligand(self):
        sel = self.ligand_list.get()
        filename = self.ligand_dic[sel].input_file
        outfile = os.path.join(self.work_dir(), os.path.basename(filename).split('.')[0]+'.pdbqt')
        util_program = os.path.join(self.autodock_tools_path.get(),"prepare_ligand4.py")
        command = "%s -l %s -o %s -A checkhydrogens" % (util_program,filename,outfile)
        self.ligand_page_log_text.insert('end',"Batch: %s\n" % command)
#        result = os.system(command)
        result, output = getstatusoutput(command)
        if result == 0:
            self.ligand_pdbqt_list.insert('end',sel)
            self.docking_ligand_list.insert('end',sel)
            self.docking_ligand_list.selectitem(-1)
            self.status_line.configure(text="Successfully generated ligand file %s" % outfile)
            #self.ligand_page_log_text.insert('end',"Successfully generated ligand file %s\n" % outfile)
            self.ligand_page_log_text.insert('end',output)
            self.ligand_pdbqt_list.selectitem(-1)
            self.ligand_dic[sel].ligand_pdbqt = outfile
        else:
            self.status_line.configure(text="An error occured while preparing ligand file from %s" % sel)
            self.ligand_page_log_text.insert('end',output)
        self.ligand_page_log_text.yview('moveto',1.0)

 
    def ligand_info(self, sel):
        try:
            sel = self.ligand_pdbqt_list.getcurselection()[0]
        except:
            sel = self.ligand_list.getcurselection()[0]
        self.ligand_page_log_text.insert('end',self.ligand_dic[sel].info())
        self.ligand_page_log_text.yview('moveto',1.0)
        

    def remove_ligand(self):
        lig = self.ligand_pdbqt_list.get()
        try:
            del self.ligand_dic[lig]
        except:
            lig = self.ligand_list.get()
            del self.ligand_dic[lig]
        try:
            index = list(self.ligand_pdbqt_list.get(0, 'end')).index(lig)
            self.ligand_pdbqt_list.delete(index)
            self.docking_ligand_list.delete(index+1)
        except:
            pass
        index = list(self.ligand_list.get(0, 'end')).index(lig)
        self.ligand_list.delete(index)
        self.status_line.configure(text="Removed ligand %s" % lig)
        try:
            self.ligand_list.selectitem(0)
        except:
            self.ligand_list.clear()
        try:
            self.ligand_pdbqt_list.selectitem(0)
        except:
            self.ligand_pdbqt_list.clear()
        try:
            self.docking_ligand_list.selectitem(0)
        except:
            self.docking_ligand_list.clear()

    def display_ligand(self):
        lig = self.ligand_list.get()
        ligand_input = self.ligand_dic[lig].input_file
        string = open(ligand_input).read()
        cmd.load(ligand_input)
        self.ligand_page_log_text.insert('end',string)
        self.ligand_page_log_text.yview('moveto',1.0)
            
        
    def remove_all_ligands(self):
        self.ligand_dic = {}
        self.status_line.configure(text="Deleted all ligand objects")
        self.ligand_list.clear()
        self.ligand_pdbqt_list.clear()


    def set_ligand_dir_location(self, dirname):
        pth = self.ligand_dir_location.setvalue(dirname)
        self.ligand_dir.set(dirname)
        
    def import_ligands(self):
        pth = self.ligand_dir.get()
        lst = glob(os.path.join(pth,"*.pdb"))\
              +glob(os.path.join(pth,"*.mol2"))\
              +glob(os.path.join(pth,"*.pdbqt"))
        self.ligand_dic['VS_DIR'] = pth
        self.multiple_ligands = True

        for f in lst:
            l = Ligand()
            l.input_file = f
            name = os.path.basename(f).split('.')[0]
            filetype = os.path.basename(f).split('.')[1]
            if filetype == 'pdbqt':
                if os.path.abspath(os.path.dirname(f)) != os.path.dirname('.'):
                    os.symlink(f,os.path.basename(f))
                    l.ligand_pdbqt = os.path.basename(f)
                else:
                    l.ligand_pdbqt = f
                self.ligand_pdbqt_list.insert('end',name)
                self.docking_ligand_list.insert('end',name)
                self.docking_ligand_list.selectitem(-1)
                self.ligand_pdbqt_list.selectitem(-1)

            self.ligand_dic[name] = l
            self.ligand_dic[name].name = name
            self.ligand_page_log_text.insert('end',"Importing ligand  %s\n" % f)
            self.ligand_page_log_text.yview('moveto',1.0)
            self.ligand_list.insert('end',name)
        self.status_line.configure(text="Loaded %d ligands" % len(lst))

    def generate_all_ligands(self):
        self.multiple_ligands = True
        for key, ligand in self.ligand_dic.items():
            if key !='VS_DIR':
                self.ligand_list.selectitem(key)
                self.generate_ligand()
        self.docking_ligand_list.selectitem(0)

    #-------------------------------------------------------------------
    # docking

    def run_autogrid(self):
        rec = self.docking_receptor_rigid_list.get()
        use_flex =  self.docking_receptor_flexible_list.get()
        receptor_object = self.receptor_dic[rec]
        if use_flex == 'No':
            rigid_receptor_file = receptor_object.receptor_pdbqt
            flex_receptor_file = None
        elif use_flex == 'Yes':
            if receptor_object.receptor_flexible:
                rigid_receptor_file = receptor_object.receptor_rigid
                flex_receptor_file = receptor_object.receptor_flexible
            else:
                self.status_line.configure(text="No flexible receptor defined! Cannot do flexible docking!")
                return
        self.docking_page_log_text.insert('end','#---------------------------------------\n')
        self.docking_page_log_text.insert('end','# SETTING UP GRID CALCULATION\n')
        self.docking_page_log_text.insert('end',' > RECEPTOR RIGID    : %s\n' % rigid_receptor_file)
        self.docking_page_log_text.insert('end',' > RECEPTOR FLEXIBLE : %s\n' % flex_receptor_file)

        ligands = self.docking_ligand_list.get()
        if ligands == 'All':
            pth = self.ligand_dic['VS_DIR']
        else:
            lig_pdbqt = self.ligand_dic[ligands].ligand_pdbqt
        self.docking_page_log_text.insert('end',' > LIGAND(s)         : %s\n' % ligands)
        if ligands == 'All':
            self.docking_page_log_text.insert('end',' > LIGAND VS-DIR      : %s\n' % pth)

        n_points_X = self.n_points_X.get()
        n_points_Y = self.n_points_Y.get()
        n_points_Z = self.n_points_Z.get()
        spacing = self.grid_spacing.get()
        center_X = self.grid_center[0].get()
        center_Y = self.grid_center[1].get()
        center_Z = self.grid_center[2].get()

        self.docking_page_log_text.insert('end',' > GRID POINTS      : %d %d %d\n' % (n_points_X, n_points_Y, n_points_Z))
        self.docking_page_log_text.insert('end',' > GRID CENTER      : %8.3f %8.3f %8.3f\n' % (center_X, center_Y, center_Z))
        self.docking_page_log_text.insert('end',' > GRID SPACING     : %8.3f \n' % spacing)

        
        template_gpf = os.path.join(self.work_dir(),'template.gpf')
        outfile_gpf = rec+'.gpf'
        fp = self.fileopen(template_gpf,'w')
        print >>fp, 'npts %d %d %d' % (n_points_X, n_points_Y, n_points_Z)
        print >>fp, 'spacing %5.3f' % spacing
        print >>fp, 'gridcenter  %8.3f %8.3f %8.3f' % (center_X, center_Y, center_Z)
        fp.close()
        util_program = os.path.join(self.autodock_tools_path.get(),"prepare_gpf4.py")
        command = "%s -r %s -i %s -o %s" % (util_program, rigid_receptor_file, template_gpf, outfile_gpf)
        if flex_receptor_file is not None:
            command+=' -x %s' % flex_receptor_file
        if hasattr(self, "multiple_ligands"):
            command+=' -d %s' % self.ligand_dic['VS_DIR']
        else:
            command+=' -l %s' % lig_pdbqt
            
        self.docking_page_log_text.insert('end',"Batch: %s\n" % command)
        self.docking_page_log_text.yview('moveto',1.0)
        result, output = getstatusoutput(command)
        if result == 0:
            self.status_line.configure(text="Successfully generated AutoGrid input file")
            self.docking_page_log_text.insert('end',output)
            self.docking_page_log_text.insert('end',"Running AutoGrid.....\n")
            autogrid = self.autogrid_exe.get()
            outfile_log = outfile_gpf.split('.')[0]+'.glg'
            receptor_object.autogrid_gpf = outfile_gpf
            receptor_object.autogrid_log = outfile_log
            command = '%s -p %s -l %s' % (autogrid, outfile_gpf, outfile_log)
            self.status_line.configure(text="Running AutoGrid....")
            # dirty
            if os.path.isfile(outfile_log):
                shutil.move(outfile_log,outfile_log+'~')
            os.system('touch %s' % outfile_log)
            r = Thread_run(command)
            r.start()
            sleep(1)
            ll = Thread_log(outfile_log, self.docking_page_log_text)
            ll.start()
        else:
            self.status_line.configure(text="An error occured while trying to run Autogrid....")
            self.docking_page_log_text.insert('end',output)
        self.docking_page_log_text.yview('moveto',1.0)            



    def run_autodock(self, write_only = False):
        rec = self.docking_receptor_rigid_list.get()
        use_flex =  self.docking_receptor_flexible_list.get()
        receptor_object = self.receptor_dic[rec]
        if use_flex == 'No':
            rigid_receptor_file = receptor_object.receptor_pdbqt
            flex_receptor_file = None
        elif use_flex == 'Yes':
            if receptor_object.receptor_flexible:
                rigid_receptor_file = receptor_object.receptor_rigid
                flex_receptor_file = receptor_object.receptor_flexible
            else:
                self.status_line.configure(text="No flexible receptor defined! Cannot do flexible docking!")
                return
        self.docking_page_log_text.insert('end','#---------------------------------------\n')
        self.docking_page_log_text.insert('end','# SETTING UP DOCKING RUN\n')
        self.docking_page_log_text.insert('end',' > RECEPTOR RIGID    : %s\n' % rigid_receptor_file)
        self.docking_page_log_text.insert('end',' > RECEPTOR FLEXIBLE : %s\n' % flex_receptor_file)

        ligands = self.docking_ligand_list.get()
        if ligands == 'All':
            pth = self.ligand_dic['VS_DIR']
        else:
            lig_pdbqt = self.ligand_dic[ligands].ligand_pdbqt
        self.docking_page_log_text.insert('end',' > LIGAND(s)         : %s\n' % ligands)
        if ligands == 'All':
            self.docking_page_log_text.insert('end',' > LIGAND VS-DIR      : %s\n' % pth)

        nposes = int(self.docking_nposes_list.get())

        self.docking_page_log_text.insert('end',' > # POSES      : %d \n' % nposes)
        
        template_dpf = os.path.join(self.work_dir(),'template.dpf')

        fp = self.fileopen(template_dpf,'w')
        print >>fp, 'ga_run %d ' % (nposes)
        fp.close()

        util_program = os.path.join(self.autodock_tools_path.get(),"prepare_dpf4.py")
        if ligands!='All':
            outfile_dpf = ligands+'.dpf'
            command = "%s -r %s -i %s -o %s " % (util_program, rigid_receptor_file, template_dpf, outfile_dpf)
            if flex_receptor_file is not None:
                command+=' -x %s' % flex_receptor_file
            else:
                command+=' -l %s' % lig_pdbqt
            self.docking_page_log_text.insert('end',"Batch: %s\n" % command)
            self.docking_page_log_text.yview('moveto', 1.0)
            result, output = getstatusoutput(command)
            if result == 0:
                self.status_line.configure(text="Wrote AutoDock input file for ligand: %s" % ligands)
#                self.docking_page_log_text.insert('end',"Running AutoDock.....\n")
                self.docking_page_log_text.insert('end',output)
                self.docking_page_log_text.yview('moveto', 1.0)
                if write_only:
                    return 
                autodock = self.autodock_exe.get()
                outfile_poses = outfile_dpf.split('.')[0]+'.dlg'
                self.ligand_dic[ligands].outfile_poses = outfile_poses 
                command = '%s -p %s -l %s' % (autodock, outfile_dpf, outfile_poses)
#                self.status_line.configure(text="Now docking ligand: %s...." % ligands)
                if os.path.isfile(outfile_poses):
                    shutil.move(outfile_poses,outfile_poses+'~')
                os.system('touch %s' % outfile_poses)
                r = Thread_run(command, self.current_thread, self.status_line, "Now docking ligand: %s...." % ligands)
                r.start()
                self.current_thread = r
                ll = Thread_log(outfile_poses, self.docking_page_log_text)
                ll.start()
            else:
                self.status_line.configure(text="Error while generating run input file for ligand: %s" % ligands)
                self.docking_page_log_text.insert('end',output)
                self.docking_page_log_text.yview('moveto', 1.0)
                
        else:
            ligand_list = []
            for lig in self.ligand_dic.values():
                if hasattr(lig,"ligand_pdbqt") and  lig.ligand_pdbqt!='':
                                    ligand_list.append(lig)
            for idx in range(len(ligand_list)):
                self.docking_ligand_list.selectitem(idx+1)
                self.run_autodock(write_only = write_only)
            

    def run_vina(self, write_only = False):

        rec = self.docking_receptor_rigid_list.get()
        use_flex =  self.docking_receptor_flexible_list.get()
        receptor_object = self.receptor_dic[rec]
        if use_flex == 'No':
            rigid_receptor_file = receptor_object.receptor_pdbqt
            flex_receptor_file = None
        elif use_flex == 'Yes':
            if receptor_object.receptor_flexible:
                rigid_receptor_file = receptor_object.receptor_rigid
                flex_receptor_file = receptor_object.receptor_flexible
            else:
                self.status_line.configure(text="No flexible receptor defined! Cannot do flexible docking!")
                return
        self.docking_page_log_text.insert('end','#---------------------------------------\n')
        self.docking_page_log_text.insert('end','# SETTING UP VINA RUN\n')
        self.docking_page_log_text.insert('end',' > RECEPTOR RIGID    : %s\n' % rigid_receptor_file)
        self.docking_page_log_text.insert('end',' > RECEPTOR FLEXIBLE : %s\n' % flex_receptor_file)


        n_points_X = self.n_points_X.get()
        n_points_Y = self.n_points_Y.get()
        n_points_Z = self.n_points_Z.get()
        spacing = self.grid_spacing.get()
        center_X = self.grid_center[0].get()
        center_Y = self.grid_center[1].get()
        center_Z = self.grid_center[2].get()

        size_X = n_points_X*spacing
        size_Y = n_points_Y*spacing
        size_Z = n_points_Z*spacing

        self.docking_page_log_text.insert('end',' > CENTER X : %s\n' % center_X)
        self.docking_page_log_text.insert('end',' > CENTER Y : %s\n' % center_Y)
        self.docking_page_log_text.insert('end',' > CENTER Z : %s\n' % center_Z)
        self.docking_page_log_text.insert('end',' > SIZE X   : %s\n' % str(size_X))
        self.docking_page_log_text.insert('end',' > SIZE Y   : %s\n' % str(size_Y))
        self.docking_page_log_text.insert('end',' > SIZE Z   : %s\n' % str(size_Z))


        ligands = self.docking_ligand_list.get()
        if ligands == 'All':
            pth = self.ligand_dic['VS_DIR']
        else:
            lig_pdbqt = self.ligand_dic[ligands].ligand_pdbqt
        self.docking_page_log_text.insert('end',' > LIGAND(s)         : %s\n' % ligands)
        if ligands == 'All':
            self.docking_page_log_text.insert('end',' > LIGAND VS-DIR      : %s\n' % pth)

        nposes = int(self.docking_nposes_list.get())

        self.docking_page_log_text.insert('end',' > # POSES      : %d \n' % nposes)
        self.docking_page_log_text.yview('moveto', 1.0)



        if ligands!='All':
            outfile_conf = os.path.join(self.work_dir(), ligands+'.vina_config.txt')
            outfile_log = os.path.join(self.work_dir(),ligands+'.vina.log')
            outfile_poses = os.path.join(self.work_dir(),ligands+'.docked.pdbqt')
            self.ligand_dic[ligands].outfile_poses = outfile_poses
            self.ligand_dic[ligands].outfile_log = outfile_log
            self.ligand_dic[ligands].vina_conf = outfile_conf
            fp = self.fileopen(outfile_conf,'w')
            print >>fp, "receptor = %s" % rigid_receptor_file
            if flex_receptor_file is not None:
                    print >>fp, "flex = %s" % flex_receptor_file
            print >>fp, "ligand = %s" % lig_pdbqt
            print >>fp, "center_x = %s" % center_X
            print >>fp, "center_y = %s" % center_Y
            print >>fp, "center_z = %s" % center_Z
            print >>fp, "size_x = %6.2f" % size_X
            print >>fp, "size_y = %6.2f" % size_Y
            print >>fp, "size_z = %6.2f" % size_Z
            print >>fp, "out = %s" % outfile_poses
            print >>fp, "log = %s" % outfile_log
            print >>fp, "num_modes = %d" % nposes
            fp.close()
            self.status_line.configure(text="Wrote VINA input file for ligand: %s" % ligands)
#            self.docking_page_log_text.insert('end',"Running VINA.....\n")
            self.docking_page_log_text.yview('moveto', 1.0)
            if write_only:
                return 
            vina = self.vina_exe.get()
            command = '%s --config %s' % (vina, outfile_conf)
#            self.status_line.configure(text="Now docking ligand: %s...." % ligands)
            if os.path.isfile(outfile_log):
                shutil.move(outfile_log,outfile_log+'~')
            os.system('touch %s' % outfile_log)
            r = Thread_run(command, self.current_thread, self.status_line, "Now docking ligand: %s...." % ligands)
            r.start()
            self.current_thread = r
            ll = Thread_log(outfile_log, self.docking_page_log_text)
            ll.start()

        else:
            ligand_list = []
            for lig in self.ligand_dic.values():
                if hasattr(lig,"ligand_pdbqt") and  lig.ligand_pdbqt!='':
                    ligand_list.append(lig)
            for idx in range(len(ligand_list)):
                self.docking_ligand_list.selectitem(idx+1)
                self.run_vina(write_only = write_only)
                        
                                
    def write_autodock_input_files(self):
        self.run_autodock(write_only = True)
    def write_vina_input_files(self):
        self.run_vina(write_only = True)
            


    #-------------------------------------------------------------------
    # view poses
    

    def set_pose_filename(self, filename):
        self.pose_file_location.setvalue(filename)

    def load_ligand_file(self, update = True):
        filename = self.pose_file_location.get()
        ext = os.path.basename(filename).split('.')[-1]
        if ext == 'pdbqt':
            #print 'loading pdbqt'
            self.load_pdbqt(filename)
        elif ext == 'dlg':
            self.load_dlg(filename)
        if update:
            self.make_complete_ligand_list()
            self.score_table.clearTable()
            self.score_table.updateView(self.all_ligands)

    def load_all_ligand_files(self):
        lst = []
        for name, ligand in self.ligand_dic.items():
            if name != 'VS_DIR' and ligand.outfile_poses!='':
                lst.append(ligand.outfile_poses)
        for f in lst:
            self.pose_file_location.setvalue(f)
            self.load_ligand_file(update=False)
        self.make_complete_ligand_list()
        self.score_table.clearTable()
        self.score_table.updateView(self.all_ligands)

    def load_dlg(self, filename):
#        print 'loading dlg file'
        name = '.'.join(os.path.basename(filename).split('.')[:-1])
        lst = self.fileopen(filename,'r').readlines()
        filt = []
        for line in lst:
            if line.startswith('DOCKED:'):
                filt.append(line[8:])
        modlist = []
        for i, line in enumerate(filt):
            if line.startswith('MODEL'):
                newmod = []
                for k, line2 in enumerate(filt[i:]):
                    if not line2.startswith('ENDMDL'):
                        newmod.append(line2)
                    else:
                        modlist.append(newmod)
                        break
        pose_list = []
            
        for m in modlist:
            model = ADModel(lst = m)
            model.name = name
            pose_list.append(model)
            
        pose_list.sort(lambda a, b: cmp(a.energy, b.energy))
        self.pose_viewer_ligand_dic[name] = {}
        for i in range(len(pose_list)):
            pose_list[i].poseN = i+1
            self.pose_viewer_ligand_dic[name].update({name+'::%d'%(i+1):pose_list[i]})
        self.update_combo(name)

    def load_pdbqt(self, filename):

        lst = self.fileopen(filename,'r').readlines()
        name = '.'.join(os.path.basename(filename).split('.')[:-1])
        modlist = []
        for i, line in enumerate(lst):
            if line.startswith('MODEL'):
                newmod = []
                for k, line2 in enumerate(lst[i:]):
                    if not line2.startswith('ENDMDL'):
                        newmod.append(line2)
                    else:
                        modlist.append(newmod)
                        break
        pose_list = []

        for m in modlist:
            model = ADModel(lst = m)
            model.name = name
#            print model.energy
            pose_list.append(model)
           # print 'done pose list'
        pose_list.sort(lambda a, b: cmp(a.energy, b.energy))
            
        self.pose_viewer_ligand_dic[name] = {}
        for i in range(len(pose_list)):
            pose_list[i].poseN = i+1
            self.pose_viewer_ligand_dic[name].update({name+'::%d'%(i+1):pose_list[i]})
            
        self.update_combo(name)


    def update_combo(self,name):
        try:
            self.pose_viewer_notebook.delete(name)
        except:
            pass
        self.pose_viewer_ligand_pages[name] = {'name':self.pose_viewer_notebook.add(name)}
        pose_list = self.pose_viewer_ligand_dic[name].keys()

        pose_list.sort(lambda a,b: cmp(int(a.split('::')[1]), int(b.split('::')[1])))
        self.pose_viewer_ligand_pages[name].update( {'poses':pose_list} )

        self.pose_viewer_buttonbox = Pmw.ButtonBox(self.pose_viewer_ligand_pages[name]['name'], padx=3)
        self.pose_viewer_buttonbox.add('Show best 10' ,command=self.show_best_poses)
        self.pose_viewer_buttonbox.add('Show all' ,command=self.show_all_poses)
        self.pose_viewer_buttonbox.add('Hide all' ,command=self.hide_all_poses)
        self.pose_viewer_buttonbox.add('Delete',command=self.delete_ligand)        
        self.pose_viewer_buttonbox.pack(fill='x',side=TOP)
        self.pose_viewer_buttonbox.alignbuttons()
        self.pose_viewer_ligand_pages[name]['combo'] = Pmw.ComboBox(self.pose_viewer_ligand_pages[name]['name'],
                                                             label_text='Docked',
                                                             labelpos='nw',
                                                             scrolledlist_items= pose_list,
                                                             selectioncommand=self.ligand_combo_box_selected,
                                                             listbox_height=10,
                                                             listbox_width=1,
                                                             dropdown=False)
        self.pose_viewer_ligand_pages[name]['combo'].pack(side='left', padx=3, anchor='n')

        self.pose_viewer_ligand_pages[name]['text'] = Pmw.ScrolledText(self.pose_viewer_ligand_pages[name]['name'],
                                                                borderframe=5, 
                                                                vscrollmode='dynamic',
                                                                hscrollmode='dynamic',
                                                                labelpos='n',
                                                                label_text=name,
                                                                text_width=150, text_height=15,
                                                                text_wrap='none',
                                                                text_background='#000000',
                                                                text_foreground='green'
                                                    )
        self.pose_viewer_ligand_pages[name]['text'].pack()
        self.pose_viewer_notebook.selectpage(name)
        self.status_line.configure(text ='Loading %s' % name)


    def ligand_combo_box_selected(self, value):
        name = value.split('::')[0]
        if self.pose_viewer_radiobuttons.getvalue()=='Show Selected':
            view = cmd.get_view()
            ligand = self.pose_viewer_ligand_dic[name][str(value)]
            cmd.read_pdbstr(ligand.as_pdb_string(),str(value))
            pmlname = value.replace(' ','_').replace('::','_')
            for key, val in self.ligand_display_mode.items():
        #        print key, val
                if val == True:
                    cmd.show(key,pmlname)
                else:
                    cmd.hide(key,pmlname)
            cmd.set_view(view)
                    #            cmd.zoom(str(value))
            text = 'Docked Energy: %8.2f kcal/mol' % ligand.energy
            self.status_line.configure(text=text)
                
            self.pose_viewer_ligand_pages[name]['text'].clear()
            self.pose_viewer_ligand_pages[name]['text'].insert('end',ligand.info_string())
        else:
            pmlname = value.replace(' ','_').replace('::','_')
            cmd.delete(pmlname)
            self.pose_viewer_ligand_pages[name]['text'].clear()
            self.status_line.configure(text='')

    def show_all_poses(self):
        name = self.pose_viewer_notebook.getcurselection()
        for pose in self.pose_viewer_ligand_pages[name]['poses']:
            self.ligand_combo_box_selected(pose)
        self.status_line.configure(text = 'Showing all poses %s' % name)

    def show_best_poses(self):
        name = self.pose_viewer_notebook.getcurselection()
        for pose in self.pose_viewer_ligand_pages[name]['poses'][:10]:
            self.ligand_combo_box_selected(pose)
        self.status_line.configure(text = 'Showing best 10 poses %s' % name)

    def hide_all_poses(self):
        name = self.pose_viewer_notebook.getcurselection()
#           pmlname = name.replace(' ','_')
#           cmd.delete(pmlname+'.*')
        cmd.delete(name+'_*')
        self.status_line.configure(text = 'Deleted all poses %s' % name)
           
    def delete_ligand(self):
        name = self.pose_viewer_notebook.getcurselection()
        cmd.delete(name+'_*')
        self.pose_viewer_notebook.delete(name)
        del self.pose_viewer_ligand_pages[name]
        del self.pose_viewer_ligand_dic[name]
        self.status_line.configure(text = 'Deleted %s' % name)
        how = self.score_table_radiobuttons.getvalue()
        self.update_score_table(how)

    def make_complete_ligand_list(self):
        self.all_ligands = []
        for key, val in self.pose_viewer_ligand_dic.items():
            self.all_ligands.extend(val.values())
        self.sort_complete_ligand_list()
        
    def sort_complete_ligand_list(self):
        self.all_ligands.sort(lambda a, b: cmp(a.energy, b.energy))
        
    def select_best_poses(self):
        self.sort_complete_ligand_list()
        self.best_ligand_list = []
        done = []
        for ligand in self.all_ligands:
            if ligand.name not in done:
                self.best_ligand_list.append(ligand)
                done.append(ligand.name)

    def update_score_table(self, value):
        if value=='Show All Poses':
            self.make_complete_ligand_list()
            self.score_table.resetTable()
            self.score_table.updateView(self.all_ligands)
        elif value=='Show Only Best Pose':
            self.select_best_poses()
            self.score_table.resetTable()
            self.score_table.updateView(self.best_ligand_list)
            
    def export_score_dat_file(self):
        what = self.score_table_radiobuttons.getvalue()
        filename = self.rank_dat_file_location.getvalue()
        fp = self.fileopen(filename,'w')
        if what == 'Show All Poses':
            ligand_list = self.all_ligands
        elif what == 'Show Only Best Pose':
            ligand_list = self.best_ligand_list
        print >>fp,"# RANK         NAME              POSE #   SCORE"
        for i, ligand in enumerate(ligand_list):
            print >>fp, "%8d %20s %5d %8.3f" % \
                  (i+1, ligand.name, ligand.poseN, ligand.energy)
        self.status_line.configure(text="Exported docking results to %s" % filename)
            
        
    def export_score_csv_file(self):
        what = self.score_table_radiobuttons.getvalue()
        filename = self.rank_csv_file_location.getvalue()
        fp = self.fileopen(filename,'w')
        if what == 'Show All Poses':
            ligand_list = self.all_ligands
        elif what == 'Show Only Best Pose':
            ligand_list = self.best_ligand_list
        print >>fp, "RANK,      NAME,       POSE,      SCORE"
        for i, ligand in enumerate(ligand_list):
            print >>fp, "%8d, %20s, %5d, %8.3f" % \
                  (i+1, ligand.name, ligand.poseN, ligand.energy)
        self.status_line.configure(text="Exported docking results to %s" % filename)

    def export_score_pose_file(self):
        what = self.score_table_radiobuttons.getvalue()
        filename = self.rank_pose_file_location.getvalue()
        fp = self.fileopen(filename,'w')
        if what == 'Show All Poses':
            ligand_list = self.all_ligands
        elif what == 'Show Only Best Pose':
            ligand_list = self.best_ligand_list
        print >>fp, "TITLE      DOCKING POSES    "
        for i, ligand in enumerate(ligand_list):
            print >>fp, "MODEL%5d" % (i+1)
            print >>fp, "REMARK    ligand: %s" % ligand.name
            print >>fp, "REMARK    pose #: %d" % ligand.poseN
            print >>fp, "REMARK    energy: %8.3f" % ligand.energy
            print >>fp, ligand.as_string.rstrip()
            print >>fp, "ENDMDL"
        self.status_line.configure(text="Exported docking poses to %s" % filename)
            

    #-------------------------------------------
    # maps

    def set_mapfilename(self,filename):
        self.map_file_location.setvalue(filename)



    def load_grid_map(self):
        view = cmd.get_view()
        fn = self.map_file_location.get()
        fp = self.fileopen(fn,'r')
        if not fp:
            return None
        try:
            name = fn#.split('.')[-2]
        except:
            name = fn
        mp = ADGridMap(fp,name)
        fname = os.path.basename(name)
        
        self.status_line.configure(text = 'Loading map %s' % fn)
        tmpn = os.path.join(self.work_dir(),fname+'.dx')
        mp.writeDX(tmpn)
        cmd.load(tmpn, fname)
        self.status_line.configure(text = 'Loading map %s' % tmpn)
        map_key = fname+'.5.0'
        self.map_color = [1,0,0]
        self.map_meta[fname] = mp.meta()
        self.map_dic[fname] = {map_key:{'color':self.map_color,
                                        'displ':'isosurface',
                                        'thresh':5,
                                        'meta':mp.meta()
                                       }
                              }
        info = '#====================================\n'
        info+= '# DISPLAY SETTINGS\n'
        info+= 'THRESHOLD: 5\n'
        info+= 'COLOR    : %4.3f %4.3f %4.3f\n' % (self.map_color[0], self.map_color[1], self.map_color[2])
        info+= 'MODE     : isosurface\n'
        info+= '#====================================\n'
        info+= '# GRID MAP INFO\n'
        info+= self.map_meta[fname]
        self.map_meta[map_key] = info
        self.map_threshold[fname] = DoubleVar()
        self.map_threshold[fname].set(float(default_settings["map_threshold"]))
        self.display_map(fname, self.map_dic[fname][map_key])
        self.update_mapcombo(fname, self.map_dic[fname])
        self.map_pages[fname]['combo'].selectitem(-1)
        self.status_line.configure(text ='Loading Map %s' % fn)
        cmd.set_view(view)
        self.status_mapbox(map_key)
        
    def display_map(self,mapname,map_dic):
        view = cmd.get_view()
        surfname = mapname+'.'+"%2.1f" % map_dic['thresh']
        cmd.delete(surfname)
        if map_dic['displ'] == 'isosurface':
            cmd.isosurface(surfname,mapname,map_dic['thresh'])
            cmd.set_color("mycol_"+surfname,map_dic['color'])
            cmd.color('mycol_'+surfname,surfname)
        elif map_dic['displ'] == 'isomesh':
            cmd.isomesh(surfname,mapname,map_dic['thresh'])
            cmd.set_color("mycol_"+surfname,map_dic['color'])
            cmd.color('mycol_'+surfname,surfname)

        info = '#====================================\n'
        info+= '# DISPLAY SETTINGS\n'
        info+= 'THRESHOLD: %4.2f\n' % map_dic['thresh']
        info+= 'COLOR    : %4.3f %4.3f %4.3f\n' % (map_dic['color'][0], map_dic['color'][1], map_dic['color'][2])
        info+= 'MODE     : %s\n' % map_dic['displ']
        info+= '#====================================\n'
        info+= '# GRID MAP INFO\n'
        info+= self.map_meta[mapname]
        self.map_meta[surfname] = info
        cmd.set_view(view)

        
    def create_surface(self):
        name = self.map_viewer_notebook.getcurselection()
        thresh = float(self.map_threshold[name].get())
        surfname = name+'.'+str(thresh)
        cmd.delete(surfname)
        if self.map_radiobuttons.getvalue()=='isosurface':
            mode = 'isosurface'
        elif self.map_radiobuttons.getvalue()=='isomesh':
            mode = 'isomesh'
        if self.map_dic[name].has_key(surfname):
            mmp = self.map_dic[name][surfname]
        else:
            self.map_dic[name].update({surfname:{
                                'color':self.map_color,
                                'displ':mode,
                                'thresh':thresh
                                }})
            mmp = self.map_dic[name][surfname]
        mmp['thresh'] = float(self.map_threshold[name].get())
        mmp['color'] = self.map_color
        self.display_map(name, mmp)
        self.update_mapcombo(name,self.map_dic[name])
        self.map_radiobuttons.setvalue(mode)
        self.map_pages[name]['combo'].selectitem(-1)
        self.status_mapbox(surfname)
            
    def delete_surface(self):
        name = self.map_viewer_notebook.getcurselection()
        s = self.map_pages[name]['combo'].getcurselection()
        for mp in s:
            self.status_line.configure(text = "Deleting map %s" % mp)
            cmd.delete(mp)
            del self.map_dic[name][mp]
        self.update_mapcombo(name,self.map_dic[name])
        try:
            self.map_pages[name]['combo'].selectitem(0)
        except:
            pass
        
    def select_map_color(self):
        color = tkColorChooser.Chooser(
            initialcolor='red',title='Choose map color').show()
        if color[0] is not None:
            self.map_color = [color[0][0]/255.,
                                color[0][1]/255.,
                                color[0][2]/255.]


    def delete_map(self):
        name = self.map_viewer_notebook.getcurselection()
        #print name
        for mp in self.map_dic[name].keys():
            cmd.delete(mp)
        cmd.delete(name+'*')
        del  self.map_pages[name]
        del self.map_meta[name]
        del self.map_dic[name]
        self.map_viewer_notebook.delete(name)

    def status_mapbox(self,key):
        name = self.map_viewer_notebook.getcurselection()
        self.map_pages[name]['text'].clear()
        self.map_pages[name]['text'].insert('end',self.map_meta[key])


    def update_mapcombo(self,name, map_dic):
        try:
            self.map_viewer_notebook.delete(name)
        except:
            pass
        self.map_pages[name] = {'name':self.map_viewer_notebook.add(name)}
#        self.map_pages[name].update({'maps':self.map_dic.keys()})

        self.upper_part = Tkinter.Frame(self.map_pages[name]['name'])
        self.lower_part = Tkinter.Frame(self.map_pages[name]['name'])
        self.upper_part.pack(side=TOP)
        self.lower_part.pack(side=BOTTOM)

        self.upper_part1 = Tkinter.Frame(self.upper_part)
        self.upper_part2 = Tkinter.Frame(self.upper_part)
        self.upper_part1.pack(side=TOP, fill='x')
        self.upper_part2.pack(side=BOTTOM, fill='x')

        self.map_thresholdfr = Tkinter.Frame(self.upper_part1)
        labmap_threshold = Label(self.map_thresholdfr,text="Threshold:");
        self.map_thresholdloc = Entry(self.map_thresholdfr,textvariable=self.map_threshold[name],bg='black',fg='green',width=15);
        self.scrmap_threshold=Scrollbar(self.map_thresholdfr,orient="horizontal",command=self.change_map_threshold)

        labmap_threshold.pack(side=LEFT)
        self.map_thresholdloc.pack(side=LEFT)
        self.scrmap_threshold.pack(side=LEFT)
        self.map_thresholdfr.pack(fill='x',padx=4,pady=3, side=LEFT) # vertical

        
        self.map_buttonbox = Pmw.ButtonBox(self.upper_part2, padx=3)

        self.map_buttonbox.add('Create Surface',command=self.create_surface)
        self.map_buttonbox.add('Delete Surface',command=self.delete_surface)        
        self.map_buttonbox.add('Color',command=self.select_map_color)        
        self.map_buttonbox.add('Delete Map',command=self.delete_map)        
        self.map_buttonbox.pack(fill='x',side=LEFT,pady=3)
        self.map_buttonbox.alignbuttons()
        self.map_pages[name]['combo'] = Pmw.ComboBox(self.lower_part,
                                                     label_text='Representations',
                                                     labelpos='nw',
                                                     scrolledlist_items= self.map_dic[name].keys(),
                                                     selectioncommand=self.status_mapbox,
                                                     listbox_height=15,
                                                     listbox_width=1,
                                                     dropdown=False)
        self.map_pages[name]['combo'].pack(side=LEFT, padx=3, anchor='n', fill='y')


        self.map_pages[name]['text'] = Pmw.ScrolledText(self.lower_part,
                                                        borderframe=5, 
                                                        vscrollmode='dynamic',
                                                        hscrollmode='dynamic',
                                                        labelpos='n',
                                                        label_text=name,
                                                        text_width=60, text_height=15,
                                                        text_wrap='none',
                                                        text_background='#000000',
                                                        text_foreground='green'
                                                    )
        self.map_pages[name]['text'].pack(fill='y')

        self.map_viewer_notebook.selectpage(name)

        self.map_radiobuttons = Pmw.RadioSelect(self.lower_part,
                                                buttontype = 'radiobutton',
                                                orient = 'horizontal',
                                                labelpos = 'w',
                                                )
        for text in ('isosurface',
                     'isomesh'):
            self.map_radiobuttons.add(text)
            self.map_radiobuttons.setvalue('isosurface')
        self.map_radiobuttons.pack(padx=4,pady=1,side='top')


    def change_map_threshold(self,a):
        current = self.map_viewer_notebook.getcurselection()
        val=float(self.map_threshold[current].get())+float(a)*0.2
        self.map_threshold[current].set(val)



    def fileopen(self, filename, mode):
        if mode=='w' and os.path.isfile(filename):
            p = os.path.abspath(filename)
            b = os.path.basename(p)
            pa,n= p.split(b)
            tmp = '#'+b+'#'
            fn = os.path.join(pa,tmp)
            os.rename(filename,fn)
            self.status_line.configure(text='Backing up %s to %s' % (filename,fn))
        try:
            fp = open(filename,mode)
            return fp
        except:
            tkMessageBox.showerror('Error','Could not open file %s' % filename)
            return None               

    
#==========================================================
#
#     FOREIGN DIALOG CLASSES

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

def quickFileValidation(s):
    if s == '': return Pmw.PARTIAL
    elif os.path.isfile(s): return Pmw.OK
    elif os.path.exists(s): return Pmw.PARTIAL
    else: return Pmw.PARTIAL


class FileDialogButtonClassFactory:
    def get(fn,mode = 'r',filter=[("Executable",'*')]):
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
                if mode == 'r':
                    n = MyFileDialog(types = filter).getopenfile()
                elif mode == 'w':
                    n = MyFileDialog(types = filter).getsavefile()
                    
#                n = MyFileDialog().get()
#                fd = PmwFileDialog(self.master,filter=filter)
#                fd.title('Please choose a file')
#                n=fd.askfilename()
                if n is not None:
                    self.fn(n)
        return FileDialogButton
    get = staticmethod(get)

class DirDialogButtonClassFactory:
    def get(fn):
        """This returns a FileDialogButton class that will
        call the specified function with the resulting file.
        """
        class DirDialogButton(Tkinter.Button):
            # This is just an ordinary button with special colors.

            def __init__(self, master=None, cnf={}, **kw):
                '''when we get a file, we call fn(filename)'''
                self.fn = fn
                self.__toggle = 0
                apply(Tkinter.Button.__init__, (self, master, cnf), kw)
                self.configure(command=self.set)
            def set(self):
                fd = PmwDirDialog(self.master)
                fd.title('Please choose a directory')
                n=fd.askfilename()
                if n is not None:
                    self.fn(n)
        return DirDialogButton
    get = staticmethod(get)

class MyFileDialog:

    def __init__(self,types = [("Executable","*")]):
        self.types = types

    def getopenfile(self):
        result = tkFileDialog.askopenfilename(filetypes=self.types)
        if result == "":
            return None
        else:
            return result
    def getsavefile(self):
        result = tkFileDialog.asksaveasfilename(filetypes=self.types)
        if result == "":
            return None
        else:
            return result
        
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
        self.lasttime=time()
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





class PmwDirDialog(PmwFileDialog):
    """Directory Dialog using Pmw"""
    def __init__(self, parent = None, **kw):
	# Define the megawidget options.
	optiondefs = (
	    ('directory', os.getcwd(),      self.newdir),
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
	dn.grid(row=1+rowoffset,column=0,columnspan=2,padx=3,pady=3)
	dn.bind('<Return>',self.okbutton)
	del dn

	# Create the directory list component.
	dnb = self.mkdnb()
	dnb.grid(row=0+rowoffset,column=0,columnspan=2,sticky='news',padx=3,pady=3)
	del dnb

	# Buttonbox already exists
	bb=self.component('buttonbox')
	bb.add('OK',command=self.okbutton)
	bb.add('Cancel',command=self.cancelbutton)
	del bb



    lastdir=""
    def fillit(self):
	"""Get the directory list and show it in the two listboxes"""
        # Do not run unnecesarily
        if self.lastdir==self['directory']:
            return
        self.lastdir=self['directory']
	dir=self['directory']
	if not dir:
	    dir=os.getcwd()
	dirs=['..']
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
	self.component('dirnamebox').setlist(dirs)

    def okbutton(self):
	"""OK action: user thinks he has input valid data and wants to
           proceed. This is also called by <Return> in the dirname entry"""
	fn=self.component('dirname').get()
	self.configure(directory=fn)
	if self.validate(fn):
	    self.canceled=0
	    self.deactivate()
    
    def askfilename(self,directory=None):
	"""The actual client function. Activates the dialog, and
	   returns only after a valid filename has been entered 
           (return value is that filename) or when canceled (return 
           value is None)"""
	if directory!=None:
	    self.configure(directory=directory)
	self.fillit()
	self.activate()
	if self.canceled:
	    return None
	else:
	    return self.component('dirname').get()

    def dirvalidate(self,string):
        if os.path.isdir(string):
            return Pmw.OK
        elif os.path.exists(string):
            return Pmw.PARTIAL
        else:
            return Pmw.OK

    def validate(self,filename):
	"""Validation function. Should return 1 if the filename is valid, 
           0 if invalid. May pop up dialogs to tell user why. Especially 
           suited to subclasses: i.e. only return 1 if the file does/doesn't 
           exist"""
        if filename=='':
            _errorpop(self.interior(),"Empty filename")
            return 0
        if os.path.isdir(filename) or not os.path.exists(filename):
            return 1
        else:
            _errorpop(self.interior(),"This is not a directory")
            return 0        





class Tail(object):
    """The Tail monitor object."""
    
    def __init__(self, path, only_new = False,
                 min_sleep = 1,
                 sleep_interval = 1,
                 max_sleep = 60):
        """Initialize a tail monitor.
             path: filename to open
             only_new: By default, the tail monitor will start reading from
               the beginning of the file when first opened. Set only_new to
               True to have it skip to the end when it first opens, so that
               you only get the new additions that arrive after you start
               monitoring. 
             min_sleep: Shortest interval in seconds to sleep when waiting
               for more input to arrive. Defaults to 1.0 second.
             sleep_interval: The tail monitor will dynamically recompute an
               appropriate sleep interval based on a sliding window of data
               arrival rate. You can set sleep_interval here to seed it
               initially if the default of 1.0 second doesn't work for you
               and you don't want to wait for it to converge.
             max_sleep: Maximum interval in seconds to sleep when waiting
               for more input to arrive. Also, if this many seconds have
               elapsed without getting any new data, the tail monitor will
               check to see if the log got truncated (rotated) and will
               quietly reopen itself if this was the case. Defaults to 60.0
               seconds.
        """

        # remember path to file in case I need to reopen
        self.path = abspath(path)
        self.f = open(self.path,"r")
        self.min_sleep = min_sleep * 1.0
        self.sleep_interval = sleep_interval * 1.0
        self.max_sleep = max_sleep * 1.0
        if only_new:
            # seek to current end of file
            file_len = stat(path)[ST_SIZE]
            self.f.seek(file_len)
        self.pos = self.f.tell()        # where am I in the file?
        self.last_read = time()         # when did I last get some data?
        self.queue = []                 # queue of lines that are ready
        self.window = []                # sliding window for dynamically
                                        # adjusting the sleep_interval

    def _recompute_rate(self, n, start, stop):
        """Internal function for recomputing the sleep interval. I get
        called with a number of lines that appeared between the start and
        stop times; this will get added to a sliding window, and I will
        recompute the average interarrival rate over the last window.
        """
        self.window.append((n, start, stop))
        purge_idx = -1                  # index of the highest old record
        tot_n = 0                       # total arrivals in the window
        tot_start = stop                # earliest time in the window
        tot_stop = start                # latest time in the window
        for i, record in enumerate(self.window):
            (i_n, i_start, i_stop) = record
            if i_stop < start - self.max_sleep:
                # window size is based on self.max_sleep; this record has
                # fallen out of the window
                purge_idx = i
            else:
                tot_n += i_n
                if i_start < tot_start: tot_start = i_start
                if i_stop > tot_stop: tot_stop = i_stop
        if purge_idx >= 0:
            # clean the old records out of the window (slide the window)
            self.window = self.window[purge_idx+1:]
        if tot_n > 0:
            # recompute; stay within bounds
            self.sleep_interval = (tot_stop - tot_start) / tot_n
            if self.sleep_interval > self.max_sleep:
                self.sleep_interval = self.max_sleep
            if self.sleep_interval < self.min_sleep:
                self.sleep_interval = self.min_sleep

    def _fill_cache(self):
        """Internal method for grabbing as much data out of the file as is
        available and caching it for future calls to nextline(). Returns
        the number of lines just read.
        """
        old_len = len(self.queue)
        line = self.f.readline()
        while line != "":
            self.queue.append(line)
            line = self.f.readline()
        # how many did we just get?
        num_read = len(self.queue) - old_len
        if num_read > 0:
            self.pos = self.f.tell()
            now = time()
            self._recompute_rate(num_read, self.last_read, now)
            self.last_read = now
        return num_read

    def _dequeue(self):
        """Internal method; returns the first available line out of the
        cache, if any."""
        if len(self.queue) > 0:
            line = self.queue[0]
            self.queue = self.queue[1:]
            return line
        else:
            return None

    def _reset(self):
        """Internal method; reopen the internal file handle (probably
        because the log file got rotated/truncated)."""
        self.f.close()
        self.f = open(self.path, "r")
        self.pos = self.f.tell()
        self.last_read = time()

    def nextline(self):
        """Return the next line from the file. Blocks if there are no lines
        immediately available."""

        # see if we have any lines cached from the last file read
        line = self._dequeue()
        if line:
            return line
        
        # ok, we are out of cache; let's get some lines from the file
        if self._fill_cache() > 0:
            # got some
            return self._dequeue()

        # hmm, still no input available
        while True:
            sleep(self.sleep_interval)
            if self._fill_cache() > 0:
                return self._dequeue()
            now = time()
            if (now - self.last_read > self.max_sleep):
                # maybe the log got rotated out from under us?
                if stat(self.path)[ST_SIZE] < self.pos:
                    # file got truncated and/or re-created
                    self._reset()
                    if self._fill_cache() > 0:
                        return self._dequeue()

    def close(self):
        """Close the tail monitor, discarding any remaining input."""
        self.f.close()
        self.f = None
        self.queue = []
        self.window = []

    def __iter__(self):
        """Iterator interface, so you can do:

        for line in filetail.Tail('log.txt'):
            # do stuff
            pass
        """
        return self

    def next(self):
        """Kick the iterator interface. Used under the covers to support:

        for line in filetail.Tail('log.txt'):
            # do stuff
            pass
        """
        return self.nextline()



class TableCell( Label ):


    params = {"relief":RIDGE, "bd":2, "bg":"#ffffff",
              "anchor":W,"font":("Helvetica",12),
              "width":1}

    def __init__(self, *args, **kwargs):

        apply(Label.__init__,((self,)+args), self.params)
        self.bind("<Enter>", self.on_enter)
        self.bind("<Leave>", self.on_leave)
        self.bind("<ButtonRelease-1>", self.on_clicked)
        self.table = kwargs["command"]
        self.entry = kwargs["entry"]
        self.row = kwargs["row"]
        self.col = kwargs["col"]

        self.HIGHLIGHT = "#ccccff"
        self.SELECTED = "#bbccdd"
        self.UNSELECTED = "#ffffff"

    def on_enter(self, event):

        event.widget["bg"] = self.HIGHLIGHT

    def on_leave(self, event):

        if self.table.selected != self.row:
            event.widget["bg"] = self.UNSELECTED
        else:
            event.widget["bg"] = self.SELECTED


    def on_clicked(self, event):

        for col in xrange(0,len(ScoreTable.fields)):
            self.entry[self.row][col]["bg"] = self.SELECTED
            if self.table.selected != -1:
                self.entry[self.table.selected][col]["bg"] = self.UNSELECTED

        self.table.on_select(self.row)
        


class ScoreTable(Frame):

    fields = {0:"Rank",1:"Ligand",2:"Pose #",3:"Score"}

    width = [10,30,10,20] 

    def __init__(self, parent, rows=16, cols = len(fields), command = None):

        Frame.__init__(self, parent, relief=SUNKEN, bd=2)
        self.windowSize = rows
        self.cols = cols
        self.rows = rows
        self.command = command
        self.table = Frame(self, relief=FLAT, bd=0)
        self.table.pack(side=LEFT, expand = YES, fill=BOTH)
        self.yscrollbar = Scrollbar(self, jump=0, command = self.on_scroll)
        self.yscrollbar.pack(side=RIGHT, fill=Y)
        self.clearTable()
        self.entry={}
        self.createTable(rows, cols)

    def clearTable(self):

        self.data = {}
        self.selected = -1
        self.rows = 0
        self.firstVisible = 0
        self.lastVisible = self.windowSize - 1
        self.yscrollbar.set(0.0,1.0)
#        self.createTable(self.rows, self.cols)
        
    def updateScrollbar(self):

        top = float(self.firstVisible)/float(self.rows)
        bottom = float(self.lastVisible+1)/float(self.rows)
        self.yscrollbar.set(top, bottom)

    def createTable(self, rows, cols):

        for col in xrange(0, len(self.fields)):
            Label(self.table, relief=SUNKEN, bd=2,
                  width=self.width[col],
                  font=("Helvetica",12), text=self.fields[col]).grid(row=0, column=col, sticky=W+E)

        for row in xrange(0, rows):
            self.entry[row] = {}
            for colx in xrange(0, cols):
                self.entry[row][colx] = TableCell(self.table, row=row,
                                                 col=colx, entry=self.entry,
                                                 command=self)
                self.entry[row][colx].grid(row=row+1,column=colx,sticky=N+W+S+E)

    def updateScreen(self, firstVisible):
        if firstVisible < 0:
            firstVisible = 0
        self.firstVisible = firstVisible
        self.lastVisible = firstVisible+self.windowSize-1
        if self.lastVisible < len(self.data.values()):
            for row in xrange(0, self.windowSize):
                for col in xrange(0, self.cols):
                    try:
                        self.entry[row][col]["text"] = self.data[row+firstVisible][col]
                    except:
                        pass


            self.updateScrollbar()

    def on_scroll(self, *args):
        if args[0] == "moveto":
            self.updateScreen(int(float(args[1])*self.rows))
        elif args[0] == "scroll":
            offset = int(args[1])
            if args[2] == "pages":
                offset*=self.windowSize
            newFirstVisible = self.firstVisible + offset
            if newFirstVisible < 0: newFirstVisible = 0
            if newFirstVisible > self.rows - self.windowSize:
                newFirstVisible = self.rows - self.windowSize
            self.updateScreen(newFirstVisible)

    def on_select(self, row):

        self.selected = row

    def append_row(self, header = None):

        self.data[self.rows] = {}
        self.rows+=1
        self.updateScrollbar()

    def set_cell_value(self, row, col, text):
        self.data[row][col] = str(text)
        if row<=self.lastVisible:
            self.entry[row-self.firstVisible][col]["text"] = text

    def add_ligand(self, ligand, rank):
        self.append_row()
        self.set_cell_value(self.rows-1, 0, rank)
        self.set_cell_value(self.rows-1, 1, ligand.name)
        self.set_cell_value(self.rows-1, 2, ligand.poseN)
        self.set_cell_value(self.rows-1, 3, ligand.energy)
        

    def updateView(self, ligand_list):
        for i, ligand in enumerate(ligand_list):
            self.add_ligand(ligand, i+1)

            
    def resetTable(self):
        for i in range(self.rows):
            for k in range(4):
                self.set_cell_value(i,k,"")
        self.clearTable()
