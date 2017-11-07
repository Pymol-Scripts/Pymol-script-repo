'''
This PyMOL plugin is described at: http://www.pymolwiki.org/index.php/Pytms
################################################################################

Author : Andreas Warnecke
email: 4ndreas.warneck3@gmail.com
Date: July 2014
License: GNU General Public License, version 2 (GPL-2.0)
Citation: Warnecke et al.: PyTMs: a useful PyMOL plugin for modeling common post-translational modifications. BMC Bioinformatics 2014 15:370.
Version: 1.2

Plugin contributed by Andreas Warnecke
(andreas.warnecke@ki.se, 4ndreas.warneck3@gmail.com)

Feel free to contact me in case of feedback
(bug reports, suggestions, comments, or questions).


################################################################################
                                     pytms.py
                               (c) Andreas Warnecke
################################################################################
ABOUT
    This is plugin designed to facilitate the modeling of
    post-translational modifications (PTMs) in silco using PyMOL.

##################################################
CITATION
    Please cite the use of PyTMs as follows:
    Warnecke et al.: PyTMs: a useful PyMOL plugin for modeling common post-translational modifications. BMC Bioinformatics 2014 15:370.

##################################################
HELP
    The respective functions can be called
        * in the PyMOL API directly
        * using Python code
        * using the interactive Plugin GUI

    The GUI offers a help button for each respectuve function.
    In PyMOL type e.g. "help mda_modify" to call detailed
    information on functions (cf. list below).
    More information is available in the citation
    and on the PyMOL wiki:
    http://www.pymolwiki.org/Pytms

    ##############################################
    CURRENTLY COVERED PTMS:
        Acetylation............acetylate
        Carbamylation:.........carbamylate
        Citrullination:........citrullinate
        Cysteine oxidation:....oxidize_cys
        Malondialdehyde:.......mda_modify
        Methionine oxidation:..oxidize_met
        Methylation............methylate
        Nitration:.............nitrate
        S-Nitrosylation:.......nitrosylate
        Proline hydroxylation:.hydroxy_pro
        Phosphorylation:.......phosphorylate
    ##############################################

##################################################
VERSION NOTES
    0.90 pre release testing
    1.00 first release
        * minor changes and typo correction
        * fixed non-incentive users receiving
          error messages related to the 'alter' command
        * new feature: integrated surface selection
        * new feature: integrated display of vdW clashes for all PTMs
        * new function: display of vdW clashes indpendent of modification
    1.1 Minor fixes
        * changed boolean processing of 'delocalized' keyword
    1.2 Additional function and minor improvements
        * fixed crashed related to random selections
        * new function: nitrosylate (Cysteine S-Nitrosylation)
        * updated usage descriptions
        * the user interface window is now resizable
        * changes in nitrate-function:
            * the torsion angle can now be defined for the nitro group and has
              a default of ~22.352 (slight angle);
              this feature is currently only supported for Tyrosines!
            * selecting a surface cutoff in conjuction with random placement
              of CE1 or CE2 for Nitro-tyrosines will now select the most
              accessible CE atom rather than a random one
        * font size can now be adjusted from the Main menu

##################################################
DISCLAIMER
    The scipts contained in this file are provided 'as is'
    without any warranty or guarantee. User discretion is advised when
    applying the scripts.
################################################################################
################################################################################
'''

from __future__ import print_function
from __future__ import absolute_import

Author = 'Andreas Warnecke'
email = '4ndreas.warneck3@gmail.com'
Date = 'October 2015'
License = 'GNU General Public License, version 2 (GPL-2.0)'
Citation = ('Citation: Warnecke et al.'
            ': PyTMs: a useful PyMOL plugin for modeling common '
            'post-translational modifications. BMC Bioinformatics 2014 15:370.')
Version = '1.2'


################################################################################
import pymol
from pymol import cmd
from chempy import cpv
from pymol import stored
import glob
import datetime
import math
import copy
import random
import sys
if sys.version_info[0] < 3:
    from Tkinter import *
else:
    from tkinter import *
import Pmw
################################################################################

################################################################################
global custom_font_regular
custom_font_regular = ('Arial', 12, 'normal')
global custom_font_bold
custom_font_bold = ('Arial', 12, 'bold')
global custom_font_unispace
custom_font_unispace = ('Consolas', 12, 'normal')
########## ##########
global keywords_default
global keywords_set
# default arguments
keywords_default={
    # basic
    'selection'         :       'all',
    'surface_cutoff'    :       0,
    'mode'              :       {
                                'Cysteine oxidation'    :    1,
                                'Methionine oxidation'  :    1,
                                'Methylation'           :    1,
                                'Nitration'             :    0,
                                'Phosphorylation'       :    0,
                                'Proline hydroxylation' :    1
                                },
    'position'          :       0,
    'position_tyr'      :       1,
    'position_trp'      :       6,
    'gct'               :       '253',
    'group'             :       2,
    'confomer'          :       5,
    'type'              :       3,
    'color_mod'         :       '',
    'color_base'        :       '',
    'hydrogens'         :       0,
    'optimize'          :       {
                                'Malondialdehyde adducts'   :    1,
                                'Phosphorylation'           :    0
                                },
    'show_clashes'      :       0,
    'include_SEC'       :       0,
    'disulfides'        :       0,
    'convert_MSE'       :       1,
    'delocalized'       :       1,

    # optimize
    'states'            :       0,
    'interpolate'       :       0,
    'local_radius'      :       10,
    'base_strain_limit' :       300,
    'intervals'         :       [45,45,45,45,45],
    'interval'          :       30,
    'remove_radius'     :       5,
    'optimize_ignore'   :       'hetatm',
    'optimize_add'      :       'none',

    # advanced
    'protonate'         :       -1,
    'torsions'          :       [106.5,180,0,-33,33,180,-180,0,120,60],
    'angles'            :       [118.35,120,120,110.5,120],
    'angle'             :       22.3519420624,

    # final
    'quiet'             :       0 # for menu, is 1 as default in script function
    } # keywords_default
# Copy for read/writing
keywords_set=copy.deepcopy(keywords_default)


################################################################################

class pytms:
    '''Main dialog for PyTMs'''
    def __init__ (self, master):
        ########## function dictionary (used in menu creation) ##########
        global funcname
        funcname=''
        global function_names
        function_names={
                'Acetylation'               : 'acetylate',
                'Carbamylation'             : 'carbamylate',
                'Citrullination'            : 'citrullinate',
                'Cysteine oxidation'        : 'oxidize_cys',
                'Malondialdehyde adducts'   : 'mda_modify',
                'Methionine oxidation'      : 'oxidize_met',
                'Methylation'               : 'methylate',
                'Nitration'                 : 'nitrate',
                'S-Nitrosylation'           : 'nitrosylate',
                'Proline hydroxylation'     : 'hydroxy_pro',
                'Phosphorylation'           : 'phosphorylate',
                ' Display vdW strain'        : 'pytms_show_clashes'
        } # function_names


        global functions_exec
        functions_exec={
                'Acetylation'               : acetylate,
                'Carbamylation'             : carbamylate,
                'Citrullination'            : citrullinate,
                'Cysteine oxidation'        : oxidize_cys,
                'Malondialdehyde adducts'   : mda_modify,
                'Methionine oxidation'      : oxidize_met,
                'Methylation'               : methylate,
                'Nitration'                 : nitrate,
                'S-Nitrosylation'           : nitrosylate,
                'Proline hydroxylation'     : hydroxy_pro,
                'Phosphorylation'           : phosphorylate,
                ' Display vdW strain'        : pytms_show_clashes
        } # function_exec


    ############################################################################
        # USER INTERFACE (PLUGIN)
    ############################################################################
        '''
        The basic window layout is:
            * a banner/ menubar at the top
            * a PTM selector to the left (thin)
            * an options menu to the right (with subframes)
            * a bottom message line
        '''
        ########## create frames ##########
        self.pytmbanner = Frame(pytmsroot, padx=5, pady=5, bg='steelblue')
        self.pytmselector = Frame(pytmsroot, padx=5, pady=5, bg='steelblue')
        self.pytmrightframe = Frame(pytmsroot, padx=5, pady=5, bg='steelblue')
        self.pytmbottom = Frame(pytmsroot, padx=5, pady=5, bg='steelblue')
        self.pytmoptions = Frame(pytmsroot, padx=5, pady=5)
        ########## Create menu bar and categories ##########
        Pmw.initialise()
        self.pytmballoon = Pmw.Balloon(pytmsroot)
        self.pytmbar = Pmw.MenuBar(pytmsroot,
                hull_relief = 'raised',
                hull_borderwidth = 1,
                balloon = self.pytmballoon)
        self.pytmbar.pack(fill = X)

        ##### main menubar entries
        self.pytmbar.addmenu('Main',
            'PyTMs', font=custom_font_regular)
        self.pytmbar.addmenu('About',
            'Info about PyTMs', font=custom_font_regular)

        ##### menubar entries in menubar/'Main'
        #FONTS
        self.pytmbar.addmenuitem('Main', 'command',
            'Increase font size',
            label='Increase font size',
            font=custom_font_regular,
            command= self.font_increase)
        self.pytmbar.addmenuitem('Main', 'command',
            'Decrease font size',
            label='Decrease font size',
            font=custom_font_regular,
            command= self.font_decrease)
        # Quit
        self.pytmbar.addmenuitem('Main', 'command',
            'Quit',
            label='Quit',
            font=custom_font_regular,
            command= lambda: pytmsroot.destroy())

        ##### menubar entries in menubar/'About'
        self.pytmbar.addmenuitem('About', 'command',
            'About info',
            label='Info about PyTMs',
            font=custom_font_regular,
            command = lambda: infotext(' About PyTMs', __doc__))
        ########## ##########

        #start layout procedure
        self.layout_general()
        self.layout_internal()
    ############################################################################
    def font_increase(self):
        global custom_font_regular
        global custom_font_bold
        global custom_font_unispace
        custom_font_regular = ('Arial', custom_font_regular[1]+1, 'normal')
        custom_font_bold = ('Arial', custom_font_bold[1]+1, 'bold')
        custom_font_unispace = ('Consolas', custom_font_unispace[1]+1, 'normal')
        pytmsroot.destroy()
        open_pytms()

    def font_decrease(self):
        global custom_font_regular
        global custom_font_bold
        global custom_font_unispace
        if (custom_font_regular[1]>1):
            custom_font_regular = ('Arial', custom_font_regular[1]-1, 'normal')
        if (custom_font_bold[1]>1):
            custom_font_bold = ('Arial', custom_font_bold[1]-1, 'bold')
        if (custom_font_unispace[1]>1):
            custom_font_unispace = ('Consolas', custom_font_unispace[1]-1, 'normal')
        pytmsroot.destroy()
        open_pytms()

    ############################################################################
    def layout_general(self):
        ########## ##########
        # main menu frames
        self.pytmbanner.pack(side=TOP, anchor=NW, fill=X, expand=1)
        self.pytmbottom.pack(side=BOTTOM, fill=X, expand=1)
        self.pytmselector.pack(side=LEFT, fill=BOTH, expand=0) # LEFT side
        self.pytmrightframe.pack(side=RIGHT, fill=Y, expand=0)
        self.pytmoptions.pack(side=LEFT, fill=BOTH, expand=1) # RIGHT side

        try:
            self.pytmoptionsmenu.destroy()
            self.pytmoptionsoptimize.destroy()
            self.pytmoptionsadvanced.destroy()
            self.pytmoptionsbar.destroy()
        except: pass

        # (re-)create option subframes
        self.pytmoptionsmenu = Frame(self.pytmoptions, padx=5, pady=5)
        self.pytmoptionsoptimize = Frame(self.pytmoptions, padx=5, pady=5)
        self.pytmoptionsadvanced = Frame(self.pytmoptions, padx=5, pady=5)
        self.pytmoptionsbar = Frame(self.pytmoptions, padx=0, pady=0)

        # Pack Options menu
        self.pytmoptionsmenu.pack(side=TOP, fill=BOTH, expand=1, anchor=W) # RIGHT top
        self.pytmoptionsbar.pack(side=BOTTOM, fill=BOTH, expand=1, anchor=S) # RIGHT bot.
        self.pytmoptionsoptimize.pack(side=LEFT, fill=BOTH, expand=0, anchor=W) # R. mid.
        self.pytmoptionsadvanced.pack(side=LEFT, fill=BOTH, expand=1, anchor=W) # R. mid.


    ############################################################################
    def layout_internal(self):
        # Internal Layouts
        self.layout_top()
        self.layout_selector()
        self.layout_options()
        self.layout_bottom()


    ############################################################################
    def layout_top(self):
        '''Layout for TOP frame (self.pytmbanner)'''
        ##### add banner
        self.headbanner = Label(self.pytmbanner, text=(
            'PyTMs: modeling post-translational modifications using PyMOL'),
            font = custom_font_bold,
            background = 'SteelBlue',
            foreground = 'white')
        self.headbanner.pack(expand=1, fill = X, padx=0, pady=0)


    ############################################################################
    def layout_selector(self):
        '''Layout for LEFT frame (self.pytmselector)'''
        ########## Left: PTM selector radiobuttons ##########
        Label(self.pytmselector,
            text='SELECT PTM:',
            font=custom_font_bold,
            background = 'SteelBlue',
            foreground = 'white')

        ########## Right: corresponding frame
        Label(self.pytmrightframe,
            text='',
            font=custom_font_bold,
            background = 'SteelBlue',
            foreground = 'white').pack(fill=Y)

        # variable for selected function/options
        self.opt_pytmselected = StringVar(master=self.pytmselector)
        sortednames=list(function_names)
        sortednames.sort()
        self.rb_selector = []
        for button in sortednames:
            self.rb_selector.append(
                Radiobutton(
                self.pytmselector,
                text = '%s'%button,
                variable = self.opt_pytmselected,
                value = '%s'%button,
                command = self.layout_options, #uses 'keywords_set' by default
                indicatoron=0,
                padx=5,
                pady=5,
                font=custom_font_regular,
                fg='white',
                bg='steelblue',
                selectcolor='mediumblue'
                ) #radiobutton
            ) # append

        # pack PTM selector radiobuttons
        for button in self.pytmselector.winfo_children():
            button.pack(fill=X)


    ############################################################################
    def layout_options(self, settings=keywords_set):
        '''
        Layout for RIGHT frame (self.pytmoptions)
        - Rebuilds the options menu based on the selector!
        '''
        global funcname

        # before changing menu, save settings
        self.store_settings()

        # which PTM is selected?
        funcname=self.opt_pytmselected.get()

        # kill all sub-widgets
        for widget in self.pytmoptionsmenu.winfo_children():
            widget.destroy()
        for widget in self.pytmoptionsoptimize.winfo_children():
            widget.destroy()
        for widget in self.pytmoptionsadvanced.winfo_children():
            widget.destroy()
        for widget in self.pytmoptionsbar.winfo_children():
            widget.destroy()
        # re-adjust layout
        self.layout_general()

        # standard startup message (no PTM selected)
        if not funcname:
            x=Label(
                self.pytmoptionsmenu,
                text = ('Welcome to PyTMs!'
                '\n\n<-- Start by selecting a PTM'),
                font = custom_font_bold)
            x.grid(row=2, column=2)
            return

        # Options menu #########################################################
        # keep track of rows for grid management
        widgetrow=0
        # cycle through function keywords
        # selection ############################################################
        for keyword in functions_exec[funcname].__code__.co_varnames:

            if keyword=='confomer': continue # see 'group'
            if keyword=='type': continue # see 'group'
            if keyword=='color_mod': continue # see color_base
            if keyword=='selection':
                Label(
                    self.pytmoptionsmenu,
                    text = 'Selection: ',
                    font = custom_font_bold
                    ).grid(
                        row=widgetrow,
                        column=0,
                        sticky=W)
                Label(
                    self.pytmoptionsmenu,
                    text = '(',
                    font = custom_font_regular
                    ).grid(
                        row=widgetrow,
                        column=1,
                        sticky=E)

                self.opt_selection = StringVar(master=self.pytmoptionsmenu)
                self.opt_selection.set(settings['selection'])
                Entry(
                    self.pytmoptionsmenu,
                    textvariable = self.opt_selection,
                    font = custom_font_regular,
                    justify=CENTER
                    ).grid(
                        row=widgetrow,
                        column=2,
                        columnspan=3,
                        sticky=W+E)
                Label(
                    self.pytmoptionsmenu,
                    text = ')',
                    font = custom_font_regular
                    ).grid(
                        row=widgetrow,
                        column=5,
                        sticky=W)

                options=cmd.get_names()
                options.extend(cmd.get_names('public_selections'))
                for selection in ['all', self.opt_selection.get()]:
                    if selection not in options: options.append(selection)
                options.sort()

                # alternative: selectable from names
                widgetrow=widgetrow+1
                Pmw.OptionMenu(
                    self.pytmoptionsmenu,
                            labelpos = 'w',
                            label_text = 'define above or choose:',
                            label_font = custom_font_regular,
                            menubutton_textvariable = self.opt_selection,
                            items = options,
                            menubutton_font = custom_font_regular
                    ).grid(
                        row=widgetrow,
                        column=2,
                        columnspan=2,
                        sticky=N+S+W+E)

                #spacer
                widgetrow=self.spacer(self.pytmoptionsmenu, widgetrow)+1
                continue # next keyword

            # surface_cutoff ###################################################
            if keyword=='surface_cutoff':
                self.opt_surface_cutoff = DoubleVar(master=self.pytmoptionsmenu)
                Label(
                    self.pytmoptionsmenu,
                    text = 'surface selection cutoff (A^2): ',
                    font = custom_font_bold
                    ).grid(
                        row=widgetrow,
                        column=0,
                        sticky=W)
                Entry(
                    self.pytmoptionsmenu,
                    textvariable = self.opt_surface_cutoff,
                    font=custom_font_regular,
                    justify=CENTER
                    ).grid(
                        row=widgetrow,
                        column=2,
                        columnspan=1,
                        sticky=W+E)
                widgetrow=widgetrow+1
                self.opt_surface_cutoff.set(settings[keyword])
                continue # next keyword

            # mode #############################################################
            if keyword=='mode':
                self.opt_mode = IntVar(master=self.pytmoptionsmenu)
                Label(
                    self.pytmoptionsmenu,
                    text = 'Mode: ',
                    font = custom_font_bold
                    ).grid(
                        row=widgetrow,
                        column=0,
                        sticky=W)
                break_at=3
                if funcname=='Cysteine oxidation':
                    options=[1,2,3,4]
                    break_at=2
                    names={
                        1: 's-hydroxy-cysteine (S-OH, CSO)',
                        2: 's-oxy-cysteine (S=O, CSX)',
                        3: 'cysteine-s-dioxide,\nR-configuration (SOOH, CSW)',
                        4: 'cysteine-s-dioxide,\nS-configuration (SOOH, CSW)'
                        }
                if funcname=='Methionine oxidation':
                    options=[1,2,3]
                    names={
                        1: 'methionine-R-sulfoxide',
                        2: 'methionine-S-sulfoxide',
                        3: 'methionine-sulfone'}
                if funcname=='Methylation':
                    options=[0,1,2,3]
                    break_at=1
                    names={
                        0: 'random (once)',
                        1: 'N-(mono)methylation',
                        2: 'N-dimethylation',
                        3: 'N-trimethylation'}
                if funcname=='Nitration':
                    options=[0,1,2]
                    names={
                        0: 'Tyrosines only',
                        1: 'Tryptophans only',
                        2: 'Both'}
                if funcname=='Phosphorylation':
                    options=[0,4,1,2,3]
                    break_at=2
                    names={
                        0: 'all residues: S/T/Y',
                        1: 'Serines only',
                        2: 'Threonines only',
                        3: 'Tyrosines only',
                        4: 'Serines and Threonines\n(no Tyrosines)'}
                if funcname=='Proline hydroxylation':
                    options=[1,2]
                    names={
                        1: '4R\n (biological - enzymatic)',
                        2: '4S\n (NB! not biological)'}

                widgetcolumn=2
                for button in options:
                    name='%s'%(names[button])
                    Radiobutton(
                        self.pytmoptionsmenu,
                        text = name,
                        variable = self.opt_mode,
                        value = button,
                        indicatoron=1,
                        padx=5,
                        pady=5,
                        font=custom_font_regular,
                        command=self.toggle_mode
                    ).grid(
                        row=widgetrow,
                        column=widgetcolumn,
                        sticky=W)
                    if widgetcolumn-1==break_at:
                        break_at=3
                        widgetcolumn=1
                        widgetrow=widgetrow+1
                    widgetcolumn=widgetcolumn+1
                widgetrow=widgetrow+1
                self.opt_mode.set(settings[keyword][funcname])
                continue # next keyword

            # position #########################################################
            if keyword=='position':
                self.opt_position = IntVar(master=self.pytmoptionsmenu)
                Label(
                    self.pytmoptionsmenu,
                    text = 'Position: ',
                    font = custom_font_bold
                    ).grid(
                        row=widgetrow,
                        column=0,
                        sticky=W)
                widgetcolumn=2
                for button in [0,1,2]:
                    name='%s'%(['Lysines only','N-termini only','Both'][button])
                    Radiobutton(
                        self.pytmoptionsmenu,
                        text = name,
                        variable = self.opt_position,
                        value = button,
                        indicatoron=1,
                        padx=5,
                        pady=5,
                        font=custom_font_regular,
                        command=self.toggle_position
                    ).grid(
                        row=widgetrow,
                        column=widgetcolumn,
                        sticky=W)
                    widgetcolumn=widgetcolumn+1
                widgetrow=widgetrow+1
                self.opt_position.set(settings[keyword])
                continue # next keyword

            # position_tyr #####################################################
            if keyword=='position_tyr':
                self.opt_position_tyr = IntVar(master=self.pytmoptionsmenu)
                Label(
                    self.pytmoptionsmenu,
                    text = 'Position, Tyrosines: ',
                    font = custom_font_bold
                    ).grid(
                        row=widgetrow,
                        column=0,
                        sticky=W)
                widgetcolumn=2

                self.rb_position_tyr = []
                for button in [0,1,2]:
                    name='%s'%({
                        0    :    '3: random or auto (if surface_cutoff is used)',
                        1    :    '3: CE1',
                        2    :    '3: CE2',
                    }[button])
                    self.rb_position_tyr.append(
                        Radiobutton(
                            self.pytmoptionsmenu,
                            text = name,
                            variable = self.opt_position_tyr,
                            value = button,
                            indicatoron=1,
                            padx=5,
                            pady=5,
                            font=custom_font_regular,
                    ))
                    self.rb_position_tyr[-1].grid(
                        row=widgetrow,
                        column=widgetcolumn,
                        sticky=W)
                    widgetcolumn=widgetcolumn+1
                    if name=='5: CZ3':
                        widgetrow=widgetrow+1
                        widgetcolumn=2
                widgetrow=widgetrow+1
                self.opt_position_tyr.set(settings[keyword])
                continue # next keyword
            # position_trp #####################################################
            if keyword=='position_trp':
                #spacer
                widgetrow=self.spacer(self.pytmoptionsmenu, widgetrow)
                self.opt_position_trp = IntVar(master=self.pytmoptionsmenu)
                Label(
                    self.pytmoptionsmenu,
                    text = 'Position, Trptophans: ',
                    font = custom_font_bold
                    ).grid(
                        row=widgetrow,
                        column=0,
                        sticky=W)
                widgetcolumn=2
                self.rb_position_trp = []
                for button in [0,4,5,6,7,1]:
                    name='%s'%({
                        0    :    '0: random (per residue)',
                        1    :    '1: CD1',
                        4    :    '4: CE3',
                        5    :    '5: CZ3',
                        6    :    '6: CH2',
                        7    :    '7: CZ2',
                    }[button])
                    self.rb_position_trp.append(
                        Radiobutton(
                            self.pytmoptionsmenu,
                            text = name,
                            variable = self.opt_position_trp,
                            value = button,
                            indicatoron=1,
                            padx=5,
                            pady=5,
                            font=custom_font_regular
                    ))
                    self.rb_position_trp[-1].grid(
                        row=widgetrow,
                        column=widgetcolumn,
                        sticky=W)
                    widgetcolumn=widgetcolumn+1
                    if name=='5: CZ3':
                        widgetrow=widgetrow+1
                        widgetcolumn=2
                widgetrow=widgetrow+1
                self.opt_position_trp.set(settings[keyword])
                continue # next keyword

            # group, confomer, type ############################################
            if keyword=='group':
                #spacer
                widgetrow=self.spacer(self.pytmoptionsmenu, widgetrow)
                Label(
                    self.pytmoptionsmenu,
                    text = 'Group/ confomer/ type: \n (ID: GCT)',
                    font = custom_font_bold
                    ).grid(
                        row=widgetrow,
                        column=0,
                        sticky=W)
                self.opt_gct = StringVar(master=self.pytmoptionsmenu)
                options=[
                'Random group, confomer and type (once) / ID=000',
                'Random group, confomer and type (by residue) / ID=111',
                'DHP (MDA/FA) - 2,6-dihyrdo-3,5-dicarbaldehyde-pyridine - cis-boat - (Hydro) MFA / ID=222',
                'DHP (MDA/FA) - 2,6-dihyrdo-3,5-dicarbaldehyde-pyridine - trans-boat - (Hydro) MFA / ID=242',
                'DHP (MDA/AA) - 4-methyl-2,6-dihyrdo-3,5-dicarbaldehyde-pyridine - cis-boat-alpha - Methyl (MAA) / ID=223',
                'DHP (MDA/MDA) - 4-ethenol-2,6-dihyrdo-3,5-dicarbaldehyde-pyridine - cis-boat-alpha - Ethenol, O in cis / ID=224',
                'DHP (MDA/MDA) - 4-ethenol-2,6-dihyrdo-3,5-dicarbaldehyde-pyridine - cis-boat-alpha - Ethenol, O in trans / ID=225',
                'DHP (MDA/MDA) - 4-ethanal-2,6-dihyrdo-3,5-dicarbaldehyde-pyridine - cis-boat-alpha - Ethanal, O in cis / ID=226',
                'DHP (MDA/MDA) - 4-ethanal-2,6-dihyrdo-3,5-dicarbaldehyde-pyridine - cis-boat-alpha - Ethanal, O in trans / ID=227',
                'DHP (MDA/AA) - 4-methyl-2,6-dihyrdo-3,5-dicarbaldehyde-pyridine - cis-boat-beta - Methyl (MAA) / ID=233',
                'DHP (MDA/MDA) - 4-ethenol-2,6-dihyrdo-3,5-dicarbaldehyde-pyridine - cis-boat-beta - Ethenol, O in cis / ID=234',
                'DHP (MDA/MDA) - 4-ethenol-2,6-dihyrdo-3,5-dicarbaldehyde-pyridine - cis-boat-beta - Ethenol, O in trans / ID=235',
                'DHP (MDA/MDA) - 4-ethanal-2,6-dihyrdo-3,5-dicarbaldehyde-pyridine - cis-boat-beta - Ethanal, O in cis / ID=236',
                'DHP (MDA/MDA) - 4-ethanal-2,6-dihyrdo-3,5-dicarbaldehyde-pyridine - cis-boat-beta - Ethanal, O in trans / ID=237',
                'DHP (MDA/AA) - 4-methyl-2,6-dihyrdo-3,5-dicarbaldehyde-pyridine - trans-boat-alpha - Methyl (MAA) / ID=243',
                'DHP (MDA/MDA) - 4-ethenol-2,6-dihyrdo-3,5-dicarbaldehyde-pyridine - trans-boat-alpha - Ethenol, O in cis / ID=244',
                'DHP (MDA/MDA) - 4-ethenol-2,6-dihyrdo-3,5-dicarbaldehyde-pyridine - trans-boat-alpha - Ethenol, O in trans / ID=245',
                'DHP (MDA/MDA) - 4-ethanal-2,6-dihyrdo-3,5-dicarbaldehyde-pyridine - trans-boat-alpha - Ethanal, O in cis / ID=246',
                'DHP (MDA/MDA) - 4-ethanal-2,6-dihyrdo-3,5-dicarbaldehyde-pyridine - trans-boat-alpha - Ethanal, O in trans / ID=247',
                'DHP (MDA/AA) - 4-methyl-2,6-dihyrdo-3,5-dicarbaldehyde-pyridine - trans-boat-beta - Methyl (MAA) / ID=253',
                'DHP (MDA/MDA) - 4-ethenol-2,6-dihyrdo-3,5-dicarbaldehyde-pyridine - trans-boat-beta - Ethenol, O in cis / ID=254',
                'DHP (MDA/MDA) - 4-ethenol-2,6-dihyrdo-3,5-dicarbaldehyde-pyridine - trans-boat-beta - Ethenol, O in trans / ID=255',
                'DHP (MDA/MDA) - 4-ethanal-2,6-dihyrdo-3,5-dicarbaldehyde-pyridine - trans-boat-beta - Ethanal, O in cis / ID=256',
                'DHP (MDA/MDA) - 4-ethanal-2,6-dihyrdo-3,5-dicarbaldehyde-pyridine - trans-boat-beta - Ethanal, O in trans / ID=257',
                'MDA - N-propenal - R1-N cis/ C1-C2 cis - enamine, O in cis / ID=322',
                'MDA - N-propenal - R1-N cis/ C1-C2 cis - enamine, O in trans / ID=323',
                'MDA - N-iminopropenol - R1-N cis/ C1-C2 cis - imine-enol, O in cis / ID=324',
                'MDA - N-iminopropenol - R1-N cis/ C1-C2 cis - imine-enol, O in trans / ID=325',
                'MDA - N-iminopropanal - R1-N cis/ C1-C2 cis - imine-al, O in cis / ID=326',
                'MDA - N-iminopropanal - R1-N cis/ C1-C2 cis - imine-al, O in trans / ID=327',
                'MDA - N-propenal - R1-N cis/ C1-C2 trans - enamine, O in cis / ID=332',
                'MDA - N-propenal - R1-N cis/ C1-C2 trans - enamine, O in trans / ID=333',
                'MDA - N-iminopropenol - R1-N cis/ C1-C2 trans - imine-enol, O in cis / ID=334',
                'MDA - N-iminopropenol - R1-N cis/ C1-C2 trans - imine-enol, O in trans / ID=335',
                'MDA - N-iminopropanal - R1-N cis/ C1-C2 trans - imine-al, O in cis / ID=336',
                'MDA - N-iminopropanal - R1-N cis/ C1-C2 trans - imine-al, O in trans / ID=337',
                'MDA - N-propenal - R1-N trans/ C1-C2 cis - enamine, O in cis / ID=342',
                'MDA - N-propenal - R1-N trans/ C1-C2 cis - enamine, O in trans / ID=343',
                'MDA - N-iminopropenol - R1-N trans/ C1-C2 cis - imine-enol, O in cis / ID=344',
                'MDA - N-iminopropenol - R1-N trans/ C1-C2 cis - imine-enol, O in trans / ID=345',
                'MDA - N-iminopropanal - R1-N trans/ C1-C2 cis - imine-al, O in cis / ID=346',
                'MDA - N-iminopropanal - R1-N trans/ C1-C2 cis - imine-al, O in trans / ID=347',
                'MDA - N-propenal - R1-N trans/ C1-C2 trans - enamine, O in cis / ID=352',
                'MDA - N-propenal - R1-N trans/ C1-C2 trans - enamine, O in trans / ID=353',
                'MDA - N-iminopropenol - R1-N trans/ C1-C2 trans - imine-enol, O in cis / ID=354',
                'MDA - N-iminopropenol - R1-N trans/ C1-C2 trans - imine-enol, O in trans / ID=355',
                'MDA - N-iminopropanal - R1-N trans/ C1-C2 trans - imine-al, O in cis / ID=356',
                'MDA - N-iminopropanal - R1-N trans/ C1-C2 trans - imine-al, O in trans / ID=357',
                'FAAB - 2-formyl-3-(alkylamino)butanal - R configuration - delocalized bonds, closed / ID=422',
                'FAAB - 2-formyl-3-(alkylamino)butanal - R configuration - delocalized bonds, open / ID=423',
                'FAAB - 2-formyl-3-(alkylamino)butanal - R configuration - regular bonds (open) / ID=424',
                'FAAB - 2-formyl-3-(alkylamino)butanal - S configuration - delocalized bonds, closed / ID=432',
                'FAAB - 2-formyl-3-(alkylamino)butanal - S configuration - delocalized bonds, open / ID=433',
                'FAAB - 2-formyl-3-(alkylamino)butanal - S configuration - regular bonds (open) / ID=434'
                ]
                Pmw.OptionMenu(
                    self.pytmoptionsmenu,
                            labelpos = 'w',
                            label_text = 'Choose adduct or enter an adduct ID:',
                            label_font = custom_font_regular,
                            menubutton_textvariable = self.opt_gct,
                            items = options,
                            menubutton_font = custom_font_regular
                    ).grid(
                        row=widgetrow,
                        column=2,
                        columnspan=2,
                        sticky=N+S+W+E)
                Entry(
                    self.pytmoptionsmenu,
                    textvariable = self.opt_gct,
                    font = custom_font_regular,
                    justify=CENTER
                    ).grid(
                        row=widgetrow,
                        column=4,
                        columnspan=1,
                        sticky=W+E)
                widgetrow=widgetrow+1
                self.opt_gct.set(settings['gct'])
                #spacer
                widgetrow=self.spacer(self.pytmoptionsmenu, widgetrow)
                continue # next keyword

            # color_base, color_mod ############################################
            if keyword=='color_base':
                widgetrow=self.spacer(self.pytmoptionsmenu, widgetrow)
                Label(
                    self.pytmoptionsmenu,
                    text = 'Coloring (optional)\n Base / PTM',
                    font = custom_font_bold
                    ).grid(
                        row=widgetrow,
                        column=0,
                        sticky=W)

                self.opt_color_base = StringVar(master=self.pytmoptionsmenu)
                self.opt_color_mod = StringVar(master=self.pytmoptionsmenu)
                options = get_colors()
                options.append('')
                options.sort()
                Entry(
                    self.pytmoptionsmenu,
                    textvariable = self.opt_color_base,
                    font=custom_font_regular,
                    justify=CENTER
                    ).grid(
                        row=widgetrow,
                        column=2,
                        columnspan=1,
                        sticky=W+E)
                Entry(
                    self.pytmoptionsmenu,
                    textvariable = self.opt_color_mod,
                    font=custom_font_regular,
                    justify=CENTER
                    ).grid(
                        row=widgetrow,
                        column=3,
                        columnspan=1,
                        sticky=W+E)
                widgetrow=widgetrow+1
                Pmw.OptionMenu(
                    self.pytmoptionsmenu,
                            labelpos = 'w',
                            label_text = '',
                            label_font=custom_font_regular,
                            menubutton_textvariable = self.opt_color_base,
                            items = options,
                            menubutton_font = custom_font_regular
                    ).grid(
                        row=widgetrow,
                        column=2,
                        columnspan=1,
                        sticky=N+S+W+E)
                Pmw.OptionMenu(
                    self.pytmoptionsmenu,
                            labelpos = 'w',
                            label_text = '',
                            label_font=custom_font_regular,
                            menubutton_textvariable = self.opt_color_mod,
                            items = options,
                            menubutton_font = custom_font_regular
                    ).grid(
                        row=widgetrow,
                        column=3,
                        columnspan=1,
                        sticky=N+S+W+E)
                widgetrow=widgetrow+1
                self.opt_color_base.set(settings['color_base'])
                self.opt_color_mod.set(settings['color_mod'])
                #spacer
                widgetrow=self.spacer(self.pytmoptionsmenu, widgetrow)
                continue # next keyword

            # hydrogens ########################################################
            if keyword=='hydrogens':
                self.opt_hydrogens = IntVar(master=self.pytmoptionsmenu)
                Label(
                    self.pytmoptionsmenu,
                    text = 'Hydrogens: ',
                    font = custom_font_bold
                    ).grid(
                        row=widgetrow,
                        column=0,
                        sticky=W)
                widgetcolumn=2
                for button in [-1,0,1]:
                    name='%s'%({
                        -1   :    'remove hydrogens',
                        0    :    'as is (detect)',
                        1    :    'add hydrogens'}[button])
                    Radiobutton(
                        self.pytmoptionsmenu,
                        text = name,
                        variable = self.opt_hydrogens,
                        value = button,
                        indicatoron=1,
                        padx=5,
                        pady=5,
                        font=custom_font_regular
                    ).grid(
                        row=widgetrow,
                        column=widgetcolumn,
                        sticky=W)
                    widgetcolumn=widgetcolumn+1
                widgetrow=widgetrow+1
                self.opt_hydrogens.set(settings[keyword])
                continue # next keyword
            # optimize #########################################################
            if keyword=='optimize':
                self.opt_optimize = IntVar(master=self.pytmoptionsmenu)
                Label(
                    self.pytmoptionsmenu,
                    text = 'Optimization level: ',
                    font = custom_font_bold
                    ).grid(
                        row=widgetrow,
                        column=0,
                        sticky=W)
                widgetcolumn=2
                self.rb_optimize = []
                if funcname=='Malondialdehyde adducts': options=[0,1,2,3,4,5]
                if funcname=='Phosphorylation': options=[0,1]
                for button in options:
                    self.rb_optimize.append(
                        Radiobutton(
                            self.pytmoptionsmenu,
                            text = '%d' %button,
                            variable = self.opt_optimize,
                            value = button,
                            indicatoron=1,
                            padx=5,
                            pady=5,
                            font=custom_font_regular,
                            command=self.toggle_optimize
                    ))
                    self.rb_optimize[-1].grid(
                        row=widgetrow,
                        column=widgetcolumn,
                        sticky=W)
                    if widgetcolumn==4:
                        widgetcolumn=1
                        widgetrow=widgetrow+1
                    widgetcolumn=widgetcolumn+1
                widgetrow=widgetrow+1
                self.opt_optimize.set(settings['optimize'][funcname])
                continue # next keyword

            # include_sec  #####################################################
            if keyword=='include_sec':
                self.opt_include_sec = IntVar(master=self.pytmoptionsmenu)
                Label(
                    self.pytmoptionsmenu,
                    text = 'Include selenocysteins?: ',
                    font = custom_font_bold
                    ).grid(
                        row=widgetrow,
                        column=0,
                        sticky=W)
                widgetcolumn=2
                for button in [1,0]:
                    Radiobutton(
                        self.pytmoptionsmenu,
                        text = ['No','Yes'][button],
                        variable = self.opt_include_sec,
                        value = button,
                        indicatoron=1,
                        padx=5,
                        pady=5,
                        font=custom_font_regular
                    ).grid(
                        row=widgetrow,
                        column=widgetcolumn,
                        sticky=W)
                    if widgetcolumn==4:
                        widgetcolumn=1
                        widgetrow=widgetrow+1
                    widgetcolumn=widgetcolumn+1
                widgetrow=widgetrow+1
                self.opt_include_sec.set(settings[keyword])
                continue # next keyword

            # disulfides   #####################################################
            if keyword=='disulfides':
                self.opt_disulfides = IntVar(master=self.pytmoptionsmenu)
                Label(
                    self.pytmoptionsmenu,
                    text = 'Modify disulfides?: ',
                    font = custom_font_bold
                    ).grid(
                        row=widgetrow,
                        column=0,
                        sticky=W)
                widgetcolumn=2
                for button in [1,0]:
                    Radiobutton(
                        self.pytmoptionsmenu,
                        text = ['No','Yes'][button],
                        variable = self.opt_disulfides,
                        value = button,
                        indicatoron=1,
                        padx=5,
                        pady=5,
                        font=custom_font_regular
                    ).grid(
                        row=widgetrow,
                        column=widgetcolumn,
                        sticky=W)
                    if widgetcolumn==4:
                        widgetcolumn=1
                        widgetrow=widgetrow+1
                    widgetcolumn=widgetcolumn+1
                widgetrow=widgetrow+1
                self.opt_disulfides.set(settings[keyword])
                continue # next keyword

            # convert_MSE  #####################################################
            if keyword=='convert_MSE':
                self.opt_convert_MSE = IntVar(master=self.pytmoptionsmenu)
                Label(
                    self.pytmoptionsmenu,
                    text = 'Convert Selenomethionines?: ',
                    font = custom_font_bold
                    ).grid(
                        row=widgetrow,
                        column=0,
                        sticky=W)
                widgetcolumn=2
                for button in [1,0]:
                    Radiobutton(
                        self.pytmoptionsmenu,
                        text = ['No','Yes'][button],
                        variable = self.opt_convert_MSE,
                        value = button,
                        indicatoron=1,
                        padx=5,
                        pady=5,
                        font=custom_font_regular
                    ).grid(
                        row=widgetrow,
                        column=widgetcolumn,
                        sticky=W)
                    if widgetcolumn==4:
                        widgetcolumn=1
                        widgetrow=widgetrow+1
                    widgetcolumn=widgetcolumn+1
                widgetrow=widgetrow+1
                self.opt_convert_MSE.set(settings[keyword])
                continue # next keyword

            # convert_MSE  #####################################################
            if keyword=='delocalized':
                self.opt_delocalized = IntVar(master=self.pytmoptionsmenu)
                Label(
                    self.pytmoptionsmenu,
                    text = 'Format as delocalized?: ',
                    font = custom_font_bold
                    ).grid(
                        row=widgetrow,
                        column=0,
                        sticky=W)
                widgetcolumn=2
                for button in [1,0]:
                    Radiobutton(
                        self.pytmoptionsmenu,
                        text = ['No','Yes'][button],
                        variable = self.opt_delocalized,
                        value = button,
                        indicatoron=1,
                        padx=5,
                        pady=5,
                        font=custom_font_regular
                    ).grid(
                        row=widgetrow,
                        column=widgetcolumn,
                        sticky=W)
                    if widgetcolumn==4:
                        widgetcolumn=1
                        widgetrow=widgetrow+1
                    widgetcolumn=widgetcolumn+1
                widgetrow=widgetrow+1
                self.opt_delocalized.set(settings[keyword])
                continue # next keyword

            # show_clashes #####################################################
            if keyword=='show_clashes':
                self.opt_show_clashes = IntVar(master=self.pytmoptionsmenu)
                Label(
                    self.pytmoptionsmenu,
                    text = 'Visualize clashes?: ',
                    font = custom_font_bold
                    ).grid(
                        row=widgetrow,
                        column=0,
                        sticky=W)
                widgetcolumn=2
                for button in [1,0]:
                    Radiobutton(
                        self.pytmoptionsmenu,
                        text = ['No','Yes'][button],
                        variable = self.opt_show_clashes,
                        value = button,
                        indicatoron=1,
                        padx=5,
                        pady=5,
                        font=custom_font_regular
                    ).grid(
                        row=widgetrow,
                        column=widgetcolumn,
                        sticky=W)
                    if widgetcolumn==4:
                        widgetcolumn=1
                        widgetrow=widgetrow+1
                    widgetcolumn=widgetcolumn+1
                widgetrow=widgetrow+1
                self.opt_show_clashes.set(settings[keyword])
                continue # next keyword

            # NB! widgetrow is used throughout and is not altered between widgets,
            # the corresponding empty rows are not gridded
            # states ###########################################################
            if keyword=='states':
                self.opt_states = IntVar(master=self.pytmoptionsoptimize)
                Label(
                    self.pytmoptionsoptimize,
                    text = 'States: ',
                    font = custom_font_bold
                    ).grid(
                        row=widgetrow,
                        column=0,
                        sticky=W)
                Entry(
                    self.pytmoptionsoptimize,
                    textvariable = self.opt_states,
                    font=custom_font_regular,
                    justify=CENTER
                    ).grid(
                        row=widgetrow,
                        column=2,
                        columnspan=1,
                        sticky=W+E)
                widgetrow=widgetrow+1
                self.opt_states.set(settings[keyword])
                continue # next keyword

            # interpolate  #####################################################
            if keyword=='interpolate':
                self.opt_interpolate = IntVar(master=self.pytmoptionsoptimize)
                Label(
                    self.pytmoptionsoptimize,
                    text = 'Interpolate optima?: ',
                    font = custom_font_bold
                    ).grid(
                        row=widgetrow,
                        column=0,
                        sticky=W)
                widgetcolumn=2
                for button in [1,0]:
                    Radiobutton(
                        self.pytmoptionsoptimize,
                        text = ['No','Yes'][button],
                        variable = self.opt_interpolate,
                        value = button,
                        indicatoron=1,
                        padx=5,
                        pady=5,
                        font=custom_font_regular
                    ).grid(
                        row=widgetrow,
                        column=widgetcolumn,
                        sticky=W)
                    if widgetcolumn==4:
                        widgetcolumn=1
                        widgetrow=widgetrow+1
                    widgetcolumn=widgetcolumn+1
                widgetrow=widgetrow+1
                self.opt_interpolate.set(settings[keyword])
                continue # next keyword

            # local_radius #####################################################
            if keyword=='local_radius':
                self.opt_local_radius = DoubleVar(master=self.pytmoptionsoptimize)
                Label(
                    self.pytmoptionsoptimize,
                    text = 'Local radius: ',
                    font = custom_font_bold
                    ).grid(
                        row=widgetrow,
                        column=0,
                        sticky=W)
                Entry(
                    self.pytmoptionsoptimize,
                    textvariable = self.opt_local_radius,
                    font=custom_font_regular,
                    justify=CENTER
                    ).grid(
                        row=widgetrow,
                        column=2,
                        columnspan=1,
                        sticky=W+E)
                widgetrow=widgetrow+1
                self.opt_local_radius.set(settings[keyword])
                continue # next keyword

            # base_strain_limit ################################################
            if keyword=='base_strain_limit':
                self.opt_base_strain_limit = DoubleVar(master=self.pytmoptionsoptimize)
                Label(
                    self.pytmoptionsoptimize,
                    text = 'Base strain limit: ',
                    font = custom_font_bold
                    ).grid(
                        row=widgetrow,
                        column=0,
                        sticky=W)
                Entry(
                    self.pytmoptionsoptimize,
                    textvariable = self.opt_base_strain_limit,
                    font=custom_font_regular,
                    justify=CENTER
                    ).grid(
                        row=widgetrow,
                        column=2,
                        columnspan=1,
                        sticky=W+E)
                widgetrow=widgetrow+1
                self.opt_base_strain_limit.set(settings[keyword])
                continue # next keyword

            # interval #########################################################
            if keyword=='interval':
                self.opt_interval = DoubleVar(master=self.pytmoptionsoptimize)
                Label(
                    self.pytmoptionsoptimize,
                    text = 'Interval: ',
                    font = custom_font_bold
                    ).grid(
                        row=widgetrow,
                        column=0,
                        sticky=W)
                Entry(
                    self.pytmoptionsoptimize,
                    textvariable = self.opt_interval,
                    font=custom_font_regular,
                    justify=CENTER
                    ).grid(
                        row=widgetrow,
                        column=2,
                        columnspan=1,
                        sticky=W+E)
                widgetrow=widgetrow+1
                self.opt_interval.set(settings[keyword])
                continue # next keyword

            # intervals ########################################################
            if keyword=='intervals':
                self.opt_intervals = StringVar(master=self.pytmoptionsoptimize)
                Label(
                    self.pytmoptionsoptimize,
                    text = 'Intervals: ',
                    font = custom_font_bold
                    ).grid(
                        row=widgetrow,
                        column=0,
                        sticky=W)
                Entry(
                    self.pytmoptionsoptimize,
                    textvariable = self.opt_intervals,
                    font=custom_font_regular,
                    justify=CENTER
                    ).grid(
                        row=widgetrow,
                        column=2,
                        columnspan=1,
                        sticky=W+E)
                widgetrow=widgetrow+1
                self.opt_intervals.set(settings[keyword])
                continue # next keyword

            # remove_radius ####################################################
            if keyword=='remove_radius':
                self.opt_remove_radius = DoubleVar(master=self.pytmoptionsoptimize)
                Label(
                    self.pytmoptionsoptimize,
                    text = 'Removal radius: ',
                    font = custom_font_bold
                    ).grid(
                        row=widgetrow,
                        column=0,
                        sticky=W)
                Entry(
                    self.pytmoptionsoptimize,
                    textvariable = self.opt_remove_radius,
                    font=custom_font_regular,
                    justify=CENTER
                    ).grid(
                        row=widgetrow,
                        column=2,
                        columnspan=1,
                        sticky=W+E)
                widgetrow=widgetrow+1
                self.opt_remove_radius.set(settings[keyword])
                continue # next keyword

            # optimize_ignore ##################################################
            if keyword=='optimize_ignore':
                self.opt_optimize_ignore = StringVar(master=self.pytmoptionsoptimize)
                Label(
                    self.pytmoptionsoptimize,
                    text = 'Ignore during optimization: ',
                    font = custom_font_bold
                    ).grid(
                        row=widgetrow,
                        column=0,
                        sticky=W)
                Entry(
                    self.pytmoptionsoptimize,
                    textvariable = self.opt_optimize_ignore,
                    font=custom_font_regular,
                    justify=CENTER
                    ).grid(
                        row=widgetrow,
                        column=2,
                        columnspan=1,
                        sticky=W+E)
                widgetrow=widgetrow+1
                self.opt_optimize_ignore.set(settings[keyword])
                continue # next keyword

            # optimize_add #####################################################
            if keyword=='optimize_add':
                self.opt_optimize_add = StringVar(master=self.pytmoptionsoptimize)
                Label(
                    self.pytmoptionsoptimize,
                    text = 'Include for optimization: ',
                    font = custom_font_bold
                    ).grid(
                        row=widgetrow,
                        column=0,
                        sticky=W)
                Entry(
                    self.pytmoptionsoptimize,
                    textvariable = self.opt_optimize_add,
                    font=custom_font_regular,
                    justify=CENTER
                    ).grid(
                        row=widgetrow,
                        column=2,
                        columnspan=1,
                        sticky=W+E)
                widgetrow=widgetrow+1
                self.opt_optimize_add.set(settings[keyword])
                continue # next keyword

            # protonate ########################################################
            if keyword=='protonate':
                self.opt_protonate = IntVar(master=self.pytmoptionsadvanced)
                Label(
                    self.pytmoptionsadvanced,
                    text = 'Protonation: ',
                    font = custom_font_bold
                    ).grid(
                        row=widgetrow,
                        column=0,
                        sticky=W)
                widgetcolumn=2
                for button in [-1,0,1]:
                    names={
                    -1 : 'deprotonated',
                    0 :  'unchanged',
                    1 :  'forced charge'
                    }[button]
                    Radiobutton(
                        self.pytmoptionsadvanced,
                        text = names,
                        variable = self.opt_protonate,
                        value = button,
                        indicatoron=1,
                        padx=5,
                        pady=5,
                        font=custom_font_regular
                    ).grid(
                        row=widgetrow,
                        column=widgetcolumn,
                        sticky=W)
                    widgetrow=widgetrow+1
                widgetrow=widgetrow+1
                self.opt_protonate.set(settings[keyword])
                continue # next keyword

            # torsions  ########################################################
            if keyword=='torsions':
                # master button to toggle
                self.opt_advanced = IntVar(master=self.pytmoptionsadvanced)
                self.rb_advanced=Checkbutton(
                    self.pytmoptionsadvanced,
                    text="set advanced options",
                    variable=self.opt_advanced,
                    onvalue=1,
                    offvalue=0,
                    indicatoron=0,
                    font=custom_font_regular,
                    command=self.toggle_advanced
                    )
                self.rb_advanced.grid(
                        row=widgetrow,
                        column=0,
                        sticky=W)
                self.opt_advanced.set(0)
                widgetrow=widgetrow+1

                self.opt_torsions = StringVar(master=self.pytmoptionsadvanced)
                Label(
                    self.pytmoptionsadvanced,
                    text = 'Torsions: ',
                    font = custom_font_bold
                    ).grid(
                        row=widgetrow,
                        column=0,
                        sticky=W)
                widgetcolumn=2
                Entry(
                    self.pytmoptionsadvanced,
                    textvariable = self.opt_torsions,
                    font=custom_font_regular,
                    justify=CENTER
                    ).grid(
                        row=widgetrow,
                        column=2,
                        columnspan=3,
                        sticky=W+E)
                widgetrow=widgetrow+1
                self.opt_torsions.set(settings[keyword])
                continue # next keyword

            # angles ###########################################################
            if keyword=='angles':
                self.opt_angles = StringVar(master=self.pytmoptionsadvanced)
                Label(
                    self.pytmoptionsadvanced,
                    text = 'Angles: ',
                    font = custom_font_bold
                    ).grid(
                        row=widgetrow,
                        column=0,
                        sticky=W)
                widgetcolumn=2
                Entry(
                    self.pytmoptionsadvanced,
                    textvariable = self.opt_angles,
                    font=custom_font_regular,
                    justify=CENTER
                    ).grid(
                        row=widgetrow,
                        column=2,
                        columnspan=3,
                        sticky=W+E)
                widgetrow=widgetrow+1
                self.opt_angles.set(settings[keyword])
                continue # next keyword

            # angle ###################################################
            if keyword=='angle':
                self.opt_angle = DoubleVar(master=self.pytmoptionsmenu)
                Label(
                    self.pytmoptionsmenu,
                    text = 'Torsion angle CZ-CE1/2-NN-O1: ',
                    font = custom_font_bold
                    ).grid(
                        row=widgetrow,
                        column=0,
                        sticky=W)
                Entry(
                    self.pytmoptionsmenu,
                    textvariable = self.opt_angle,
                    font=custom_font_regular,
                    justify=CENTER
                    ).grid(
                        row=widgetrow,
                        column=2,
                        columnspan=1,
                        sticky=W+E)
                widgetrow=widgetrow+1
                self.opt_angle.set(settings[keyword])
                continue # next keyword

        # quiet ############################################################
            if keyword=='quiet':
                self.opt_quiet = IntVar()
                Label(
                    self.pytmoptionsmenu,
                    text = 'Verbosity: ',
                    font = custom_font_bold
                    ).grid(
                        row=widgetrow,
                        column=0,
                        sticky=W)
                Checkbutton(
                    self.pytmoptionsmenu,
                    text="quiet (no output on progress etc.)",
                    variable=self.opt_quiet,
                    onvalue=1,
                    offvalue=0,
                    font=custom_font_regular
                    ).grid(
                        row=widgetrow,
                        column=2,
                        sticky=W)
                widgetrow=widgetrow+1
                self.opt_quiet.set(settings[keyword])
                # final notice #################################################
                if funcname== ' Display vdW strain':
                    displaymessage=('This function is not for any PTM, but can be used to\n'
                               'visualize van-der-Waals clashes independently, i.e.\n'
                               'retrospectively, or in unmodified objects.\n'
                               'Note that hydrogens significantly affect the \n'
                               'vdW strain display and calculation.')
                    Label(
                        self.pytmoptionsmenu,
                        text = 'Note: ',
                        font = custom_font_bold
                        ).grid(
                            row=widgetrow,
                            column=0,
                            sticky=N)
                    Label(
                        self.pytmoptionsmenu,
                        text = '%s'%displaymessage,
                        justify = 'left',
                        font = custom_font_regular
                        ).grid(
                            row=widgetrow,
                            column=2,
                            sticky=W)
                    self.opt_quiet.set(settings[keyword])
                break # quiet is the last keyword

        # equalize colums
        for c in [2,3,4]:
            self.pytmoptionsmenu.columnconfigure(c, weight=1, pad=1, minsize=1, uniform='a')
        for c in [2,3,4]:
            self.pytmoptionsadvanced.columnconfigure(c, weight=1, pad=1, minsize=1, uniform='b')

        # Function Buttons (Bottom Frame) ######################################
        # adjust name of button
        if (funcname != ' Display vdW strain'):
            buttonexecstring='Modify: %s!'%funcname
        else:
            buttonexecstring='%s!'%funcname

        execbutton = Button(
            self.pytmoptionsbar,
            text = '%s'%buttonexecstring,
            command = self.run_function,
            font=custom_font_regular,
            background='steelblue',
            foreground='white')
        execbutton.pack(side=RIGHT, expand=0)
        helpbutton = Button(
            self.pytmoptionsbar,
            text = 'Help',
            command = lambda: infotext(funcname, functions_exec[funcname].__doc__),
            font=custom_font_regular,
            background='steelblue',
            foreground='white')
        helpbutton.pack(side=LEFT, expand=0)
        defaultbutton = Button(
            self.pytmoptionsbar,
            text = 'Reset defaults',
            command = lambda: self.layout_options(settings=keywords_default),
            font=custom_font_regular,
            background='steelblue',
            foreground='white')
        defaultbutton.pack(side=LEFT, padx=10, expand=0)
        ########################################################################

        # Toggles

        self.toggle_position()
        self.toggle_optimize()
        self.toggle_mode()
        self.toggle_advanced()

        # update label
        self.layout_bottom()
    ############################################################################

    ############################################################################
    def layout_bottom(self):
        '''Layout for BOTTOM frame (self.pytmbottom)'''
        global funcname
        try: self.pytmmessage.destroy()
        except: pass
        self.pytmmessage = Label(
            self.pytmbottom,
            font=custom_font_regular,
            background = 'SteelBlue',
            foreground = 'white')

        # which PTM is selected?
        funcname=self.opt_pytmselected.get()

        if not funcname:
            self.pytmmessage.config(text='To start, please select a PTM!')
        else:
            self.pytmmessage.config(text='Selected: %s'%funcname)

        self.pytmmessage.pack(fill=BOTH, expand=1)


    ############################################################################
    def spacer(self, widget, widgetrow):
        '''
        #inserts an empty row
        '''
        Label(
            widget,
            text = ' ',
            font = custom_font_bold
            ).grid(
                row=widgetrow,
                column=0,
                sticky=W)
        return widgetrow+1

    ############################################################################
    def toggle_position(self):
        '''enables/ disables options depending on position'''
        try:
            position = self.opt_position.get()
            for button in self.rb_optimize:
                button.config(state=NORMAL)
            if position==1: # N-termini only
                for button in self.rb_optimize[2:]:
                    button.config(state=DISABLED)
        except: return

    ########################################################################
    def toggle_optimize(self):
        '''enables/ disables options depending on optimize'''
        try:
            optimize = self.opt_optimize.get()
            if optimize==0:
                for button in self.pytmoptionsoptimize.winfo_children():
                    button.config(state=DISABLED)
            else:
                for button in self.pytmoptionsoptimize.winfo_children():
                    button.config(state=NORMAL)
        except: return

    ########################################################################
    def toggle_mode(self):
        '''enables/ disables options depending on mode'''
        try:
            mode = self.opt_mode.get()
            # Nitrations
            if mode==0: # 'Tyrosines only'
                for button in self.rb_position_tyr:
                    button.config(state=NORMAL)
                for button in self.rb_position_trp:
                    button.config(state=DISABLED)
            elif mode==1: # 'Tryptophans only'
                for button in self.rb_position_tyr:
                    button.config(state=DISABLED)
                for button in self.rb_position_trp:
                    button.config(state=NORMAL)
            else: # Both
                for button in self.rb_position_tyr:
                    button.config(state=NORMAL)
                for button in self.rb_position_trp:
                    button.config(state=NORMAL)
        except: return

    ########################################################################
    def toggle_advanced(self):
        '''enables/ disables advanced options'''
        try:
            advanced = self.opt_advanced.get()
            for button in self.pytmoptionsadvanced.winfo_children():
                button.config(state=DISABLED)
            if advanced==1:
                for button in self.pytmoptionsadvanced.winfo_children():
                    button.config(state=NORMAL)
            # always enable toggle button
            self.rb_advanced.config(state=NORMAL)
        except: return

    ############################################################################
    def store_settings(self):
        '''stores all altered settings'''
        global funcname
        if not funcname: return # skip unless function was defined
        try: keywords_set['selection'] = self.opt_selection.get()
        except: pass
        try: keywords_set['surface_cutoff'] = self.opt_surface_cutoff.get()
        except: pass
        try: keywords_set['mode'][funcname] = self.opt_mode.get()
        except: pass
        try: keywords_set['position'] = self.opt_position.get()
        except: pass
        try: keywords_set['position_tyr'] = self.opt_position_tyr.get()
        except: pass
        try: keywords_set['position_trp'] = self.opt_position_trp.get()
        except: pass
        try:
            gct=self.opt_gct.get() # ends with 3-digit ID
            keywords_set['gct'] = gct
            keywords_set['group'] = gct[-3]
            keywords_set['confomer'] = gct[-2]
            keywords_set['type'] = gct[-1]
        except: pass
        try:
            keywords_set['color_base'] = self.opt_color_base.get()
            keywords_set['color_mod'] = self.opt_color_mod.get()
        except: pass
        try: keywords_set['hydrogens'] = self.opt_hydrogens.get()
        except: pass
        try: keywords_set['optimize'][funcname] = self.opt_optimize.get()
        except: pass
        try: keywords_set['show_clashes'] = self.opt_show_clashes.get()
        except: pass
        try: keywords_set['include_SEC'] = self.include_SEC.get()
        except: pass
        try: keywords_set['disulfides'] = self.opt_disulfides.get()
        except: pass
        try: keywords_set['convert_MSE'] = self.opt_convert_MSE.get()
        except: pass
        try: keywords_set['delocalized'] = self.opt_delocalized.get()
        except: pass
        try: keywords_set['states'] = self.opt_states.get()
        except: pass
        try: keywords_set['interpolate'] = self.opt_interpolate.get()
        except: pass
        try: keywords_set['local_radius'] = self.opt_local_radius.get()
        except: pass
        try: keywords_set['base_strain_limit'] = self.opt_base_strain_limit.get()
        except: pass
        try: keywords_set['intervals'] = self.opt_intervals.get()
        except: pass
        try: keywords_set['interval'] = self.opt_interval.get()
        except: pass
        try: keywords_set['remove_radius'] = self.opt_remove_radius.get()
        except: pass
        try: keywords_set['optimize_ignore'] = self.opt_optimize_ignore.get()
        except: pass
        try: keywords_set['optimize_add'] = self.opt_optimize_add.get()
        except: pass
        try: keywords_set['protonate'] = self.opt_protonate.get()
        except: pass
        try: keywords_set['torsions'] = self.opt_torsions.get()
        except: pass
        try: keywords_set['angles'] = self.opt_angles.get()
        except: pass
        try: keywords_set['angle'] = self.opt_angle.get()
        except: pass
        try: keywords_set['quiet'] = self.opt_quiet.get()
        except: pass


    ############################################################################
    def run_function(self):
        '''
        passes keyword arguments and executes selectd function
        '''
        global funcname

        # save settings
        self.store_settings()
        # which PTM is selected?
        funcname=self.opt_pytmselected.get()

        # create and append kwarg dictionary
        keywordarguments={}
        for keyword in functions_exec[funcname].__code__.co_varnames:
            if keyword in ['mode', 'optimize']:
                keywordarguments[keyword]=(keywords_set[keyword][funcname])
            else:
                keywordarguments[keyword]=(keywords_set[keyword])
            if keyword=='quiet': break # final keyword

        # Finally, run it:
        functions_exec[funcname](**keywordarguments)

################################################################################
# Startup and menu functions
################################################################################
def get_colors(selection='', quiet=1):
    '''an integrated version of get_colors'''
    import pymol
    pymol_color_list=[]
    for tuplepair in pymol.querying.get_color_indices(selection):
        pymol_color_list.append(tuplepair[0])
    pymol_color_list.sort()
    return pymol_color_list

################################################################################
def infotext(title=' INFO ',messagetext=''):
    '''Displays an independent infoscreen'''
    inforoot = Tk()
    inforoot.title(title)

    # Text
    info = Text(
        inforoot,
        bg='white',
        fg='black',
        font=custom_font_unispace,
        wrap=NONE
        )
    info.config(state=NORMAL)
    info.insert(INSERT, messagetext)
    info.config(state=DISABLED) #supresses editing

    #Scrollbars
    info_xscrollbar = Scrollbar(inforoot, orient=HORIZONTAL)
    info_yscrollbar = Scrollbar(inforoot)
    info_xscrollbar.pack(side=BOTTOM, fill=X)
    info_yscrollbar.pack(side=RIGHT, fill=Y)

    # attach
    info.config(
        xscrollcommand=info_xscrollbar.set,
        yscrollcommand=info_yscrollbar.set)

    # pack window
    info_xscrollbar.config(command=info.xview)
    info_yscrollbar.config(command=info.yview)
    info.pack(expand=1, fill=BOTH)

    inforoot.mainloop()

################################################################################
def __init__(self):
    '''Add PyTMs to PyMOL Plugin menu'''
    # trigger PyMOL 2.0 legacy initialization
    self.root

    self.menuBar.addmenuitem('Plugin', 'separator')
    self.menuBar.addmenuitem('Plugin', 'command',
                             'PyTMs',
                             label = 'PyTMs',
                             command = open_pytms)
    self.menuBar.addmenuitem('Plugin', 'separator')

################################################################################
def open_pytms():
    '''opens main dialog'''
    global pytmsroot
    global pytmsvar
    pytmsroot = Tk()
    pytmsroot.title(' PyTMs ')
    pytmsroot.resizable(1,1) # prevent/allow rezise
    Pmw.initialise()
    pytmsvar = pytms(pytmsroot)

################################################################################
# DONE WITH MENU-ASSOCIATED FUNCTIONS: BELOW ARE FUNCTIONS RELATED TO PTMs
################################################################################





################################################################################
#Repeatedly-used functions
################################################################################


################################################################################
# selection and startup variables
def pytms_get_sele(selection='none', surface_cutoff=0):
    '''
    method startup: tests PyMOl version,
    verifies selection and defines temporary variables
    If surface_cutoff is set>0.0, then the selection will be processed to include
    surface exposed atoms.
    '''
    # prevent resi-splitting bug
    cmd.set('retain_order', 0)

    # Global variables
    global objects
    global names
    global temp_sel
    global temp_obj
    global temp_clash_A
    global version_ok

    ##### Get version info
    try:
        version_ok=cmd.get_version()
        if ((version_ok[1]>=1.7) and (pymol.invocation.options.incentive_product==1)):
            version_ok=True
        else: version_ok=False
    except:
        version_ok=False
        try:
            cmd.get_unused_name('temp_obj')
            # Doesn't work in ancient PyMOL versions
        except:
            raise Exception(("Outdated PyMOL!\n#####\n"
            "Sorry, PyTMs does not support this PyMOL version!\n"
            "Please install at least PyMOL v1.3 or higher!\n#####"))


    # temporary selections/ object names
    temp_obj=cmd.get_unused_name('temp_obj')
    temp_obj=cmd.get_unused_name('temp_obj')
    temp_clash_A=cmd.get_unused_name('temp_clash_A')
    temp_sel=cmd.get_unused_name('temp_sel')

    ##### Selection
    # fixed occasional problem with single objects requiring "model"
    # affects later selecetion using selection
    try:
        # group selection with bracketing
        selection='('+str(selection)+')'
        # checks functional selection
        cmd.count_atoms(selection)
    except:
        # attempt to fix (works for single model)
        selection=selection.replace('(', '(model ',1)
        try:
            if (len(cmd.get_object_list(selection))!=1):
                # not single object
                raise Exception(("Malformed selection\n#####\n"
                "Sorry, there is a problem with the input selection.\n"
                "Some objects may require adding 'model'\n#####"))
        except:
            raise Exception(("Malformed selection\n#####\n"
            "Sorry, there is a problem with the input selection.\n#####"))
    # otherwise proceed

    #Check input for states and alternative coordinates
    if (cmd.count_states(selection)>1):
        raise Exception(("Unsupported input\n"
        "#####\nThe selection contains more than one state.\n"
        "Provide object copies reduced to one state (see wiki for help)!\n"
        "try,e.g.:\ncmd.create('newobject', 'oldobject', '1', '1') "
        "# copy with only state 1\n#####"))

    #Check input for states and alternative coordinates
    if (cmd.count_atoms('%s and ((alt *) and not (alt ""))' %selection)>0):
        raise Exception(("Unsupported input\n"
        "#####\nThe selection contains alternative atom coordinates.\n"
        "Provide object copies reduced to one alternative (see wiki for help)!\n"
        "try,e.g.:\ncmd.create('newobject','oldobject and (alt \"\" or alt A)') "
        "# common and alt A\n"
        "cmd.alter('newobject', 'alt=\"\"') # reset alt\n#####"))

    # SELECTION LIST ###########################################################
    # variable for the selections
    # get the names of the proteins in the selection
    objects = cmd.get_object_list(selection)

    # INTEGRATED SURFACE SELECTION ####
    # cf. findSurfaceResidues by Jason Vertrees (PyMOL wiki)
    if surface_cutoff!=0:
        if (cmd.get('dot_solvent') in [0,'off']):
            print("PyTMs: NB! setting 'dot_solvent' to 1 for surface selection.")
            cmd.set('dot_solvent', 1)

        tmpObj = cmd.get_unused_name("_tmp")
        cmd.create(tmpObj, "(byobj (%s)) and not hetatm"%selection, zoom=0)
        cmd.get_area(selection=tmpObj, load_b=1)

        # threshold on what one considers an "exposed" atom (in A**2):
        cmd.remove('%s and b<%f' %(tmpObj,surface_cutoff))

        cmd.select('pytms_surface_input', '(byobj (%s)) in %s' %(selection, tmpObj))
        cmd.delete(tmpObj)

        # modify input selection to include area cutoff selection
        selection='(%s) and (pytms_surface_input)'%selection
    #####

    #include chains
    # subselect chains
    names=[]
    for p in objects:
        for chain in cmd.get_chains('model %s'%p) or ['']:
            names.append("(model %s and chain '%s')"%(p, chain))
    ##### SELECTION LIST DONE #####

    #proceed and return selection
    return selection
################################################################################


################################################################################
def verify_sele(selection='none'):
    '''returns valid selection'''
    try:
        # group selection with bracketing
        selection='('+str(selection)+')'
        # checks functional selection
        cmd.count_atoms(selection)
    except:
        raise Exception(("Malformed selection:\n#####\n"
        "%s\n"
        "is invalid. Please check this selection.\n#####"%selection))
    #proceed and return selection
    return selection
################################################################################


################################################################################
def pytms_rebuild():
    '''unpicks, rebuilds and deletes temporary objects'''
    #delete temporary selections and objects
    cmd.delete(temp_sel)
    cmd.delete(temp_obj)
    cmd.delete(temp_clash_A)
    # update PyMOL
    cmd.unpick()
    cmd.rebuild()
################################################################################


################################################################################
# immature function to select defined parts
def select_part(name,obj,start,limit):
    '''selects "name" from "object" and "start" but not beyond "limit"'''
    cmd.select(name,'%s and (%s)' %(obj,start))
    last_count=0
    increment=1
    while increment>0:
        cmd.select(name,
        ('((%s) extend 1) and (%s) '
        'and not (%s)' %(name,obj,limit)))
        increment=(cmd.count_atoms(name)-last_count)
        last_count=last_count+increment
################################################################################


##### START OF SET_ANGLE #######################################################
# immature function to set angles during editing
def set_angle(
angle=90,
atom1='pk1',
atom2='pk2',
atom3='pk3',
selection='',
axis=''):
    '''
    immature function for use during editing:
    rotates atom 3 and eveything distal or specified selection
    around perpendicular axis fromed by
    atom2-atom1 and atom2-atom3 to a set angle
    axis calculation is overridden by setting 'axis'
    '''
    from pymol import cmd
    from chempy import cpv
    import math

    def get_coord(v):
        '''gets coordinates'''
        if not isinstance(v, str):
            return v
        if v.startswith('['):
            try:
                return cmd.safe_list_eval(v)
            except:
                v=[round(y,6) for y in v]
                return cmd.safe_list_eval(v)
        return cmd.get_atom_coords(v)

    # get vectors
    ori = get_coord(atom2)
    ori = [round(y,6) for y in ori]
    ax1 = cpv.sub(get_coord(atom1),ori)
    ax2 = cpv.sub(get_coord(atom3),ori)

    # angle in rad
    try:
        rotang = (
        (((-1*cpv.get_angle(ax1, ax2))*180/math.pi)+float(angle)))
    except:
        rotang = (float(angle))
    rotang=round(rotang,6)

    # rotation axis
    if not axis:
        ax3 = cpv.get_system2(ax1, ax2)[2]
    else:
        ax3 = get_coord(axis)
    if ax3==[0.0, 0.0, 0.0]:
        print("Cannot calculate rotation axis!")
        ax3=False
        return ax3
    ax3=[round(y,6) for y in ax3]

    # refine selection
    if not selection:
        selection=cmd.get_unused_name('temp')
        cmd.select(selection,atom3)
        increment=1
        last_count=cmd.count_atoms(selection)
        if last_count==0:
            print("atom 3 does not define a selection! - aborting")
            return False
        # select everything distal including atom3
        while increment>0:
            cmd.select(selection,
            '((%s) extend 1) and not (%s)' %(selection,atom2))
            increment=(cmd.count_atoms(selection)-last_count)
            last_count=last_count+increment
    else:
        selection='('+str(selection)+')'

    cmd.rotate(axis=ax3, angle=rotang, camera='0',
    selection=selection, origin=ori)
    cmd.delete(selection)
    return ax3
################################################################################


# immature functions for angle conversion ######################################
def to360(v):
    '''convert angle to 360 degrees'''
    while v>360:
        v=v-360
    while v<0:
        v=v+360
    return v

def to180(v):
    '''convert angle to +/- 180 degrees'''
    while v>180:
        v=v-360
    while v<-180:
        v=v+360
    return v
################################################################################


##### START OF GET_STRAIN ######################################################
# Function for strain calculation and clashing
def get_strain(obj, objname, optimize=0):
    '''
    # credits to show_bumps.py
    '''
    optimize=bool(int(optimize))
    cmd.delete(objname)
    cmd.create(objname, '(%s)' %obj, zoom=0)
    cmd.sculpt_activate(objname)
    cmd.show_as('cgo', objname)
    cmd.set('sculpt_vdw_vis_mode', 1, objname)
    cmd.set('sculpt_field_mask', 0x020) # cSculptVDW
    strain=cmd.sculpt_iterate(objname, cycles=1)
    if strain==0:
        if ((optimize) and (cmd.count_atoms('(%s) and elem H'%obj)!=0)):
            raise Exception(("VDW strain acquisition failed:\n#####\n"
            "Sorry, an error (zero strain) occured during optimization or strain calulation!\n"
            "Try restarting PyMOL (or disable optimization / display of clashes).\n#####"))
        else:
            print('Warning! Unexpected zero strain occured during VdW strain calculation!')
    return strain
################################################################################


##### Logging Progress #########################################################
def log_pytms_prog(output=False, message='', count='', total='', last=''):
    '''log time, message and progress'''
    now=datetime.datetime.now()
    message=str(message)

    if output:
        if ((count!='') and (total!='')):
            # log message and progress
            count=float(count)
            total=float(total)
            percent=count/total*100.00
            print ("[%s]: %s | %.2f%% [%.0f of %.0f]"
            %(now, message, percent, count, total))
            if last!='':
                # also log ETA
                # based on last step
                try:
                    finish=(
                    datetime.timedelta(seconds=(100-percent)
                    *(now-last).total_seconds()/percent)+now)
                except AttributeError:
                    #PyMOL 1.3 uses python 2.5 and has no .total_seconds()
                    td=(now-last) #timedelta
                    # conversion (avoiding ctypes dependency)
                    td=(float((td.microseconds +
                      (td.seconds + td.days * 24 * 3600) * 10**6)) / 10**6)
                    finish=(
                    datetime.timedelta(seconds=(100-percent)
                    *td/percent)+now)
                # print finish
                print("ETA of completion: [%s]"%finish)
        else:
            # log time and message only
            print("[%s]: %s" %(now,message))
    # can be used to update time, even without output
    return now
################################################################################


################################################################################
def get_resi_macro_name(sele):
    '''
    will return the resi name of the last iteration
    accepts non-macro selection style only
    '''
    from pymol import stored
    # name macro
    try:
        sele='first (%s)'%sele # reduce to one atom
        stored.macroname=''
        cmd.iterate(sele, 'stored.macroname="/%s/%s/%s/%s`%s/" %(model, segi, chain, resn, resi)')
        return stored.macroname
    except:
        raise Exception(("Invalid selection:\n#####\n"
            "Error in get_resi_macro_name\n#####"))
################################################################################


########## Evaluate hydrogens ##################################################
def eval_hydrogens(selection, hydrogens):
    try: hydrogens=int(hydrogens)
    except: hydrogens=0
    if (not int(hydrogens)):
        hydrogens=True #default
        if (cmd.count_atoms('((%s) and not hetatm) and elem H'%selection)==0):
            hydrogens=False # if no H present
    elif hydrogens<0: hydrogens=False
    elif hydrogens>0: hydrogens=True
    return hydrogens
################################################################################

################################################################################
# END OF REPEATEDLY-USED FUNCTIONS: BELOW ARE PTM FUNCTIONS
################################################################################





################################################################################
# Integrated display of vdW clashes (possible for unmodified)
################################################################################
def pytms_show_clashes(
selection='all',
hydrogens=0,
quiet=1
):


    '''
DESCRIPTION

    An integrated version of the 'show_bumps.py' script, originally by Thomas Holder
    Displays the vdW clashed between atoms using colored discs.
    This function is intended to allow display these clashes independently of
    modification as an integrated part of PyTMs.

EXAMPLE

    frag ARG
    pytms_show_clashes

USAGE

    pytms_show_clashes [ selection [hydrogens [, quiet ]]]

ARGUMENTS

    selection: selection to be used {default: 'all'}
               Note that the vdW clashes will be shown and calculated objectwise
               and ignore inter-object clashes (e.g. in case of overlapping models)

    hydrogens: <int>     toggles if the object will have hydrogens or not
                         0  : as is (detection) {default}
                         1  : adds hydrogens;
                         -1 : no hydrogens;
    BE AWARE THAT THE HYDROGENS AFFECT THE VDW CALCULATION SIGNIFICANTLY!!!

    quiet: <boolean> toggles output {default: quiet=1}
                     updates on operation and remaining time
                     * NB! print appears first in the console

    '''
    ##### BEGINNING OF SCRIPT #####

    last_time=start_time=datetime.datetime.now()
    selection=pytms_get_sele(selection, 0)

    # argument settings
    try:
        # transform boolean to reduce checking loops
        hydrogens=eval_hydrogens(selection, hydrogens)
        show_clashes=True
        output=bool(not int(quiet))
    except:
        raise Exception("Input error!\n Please check the input parameters!")

    # SELECTION LIST ###########################################################
    # VDW PREP
    # OBJCHI (dictionary) structure
    # obj: [strain]
    OBJCHIS = {}
    if objects==[]:
        print('PyTMS: no objects in selection!')
        return OBJCHIS

    if show_clashes:
        # get unmodified VDW strain
        log_pytms_prog(output, 'Calculating VdW strain!')
        for p in objects:
            if (hydrogens):
                cmd.h_add(p)
            else:
                cmd.remove('((%s) and elem H)' %p)
            OBJCHIS[p] = [get_strain('(%s)'%(p), temp_clash_A, 0)]
            cmd.set_name(temp_clash_A,'%s_clashes'%p)
        if output:
            print('PyTMs: STRAIN REPORT:')
            print('OBJECT','vdW_STRAIN')
            for p in objects:
                print('%s'%p,OBJCHIS[p][0])
    # rebuild
    pytms_rebuild()

    #exit
    last_time=log_pytms_prog(output, 'vdW calculation and display done!')
    return OBJCHIS

cmd.extend( "pytms_show_clashes", pytms_show_clashes );

################################################################################
################################################################################





################################################################################
# Automated in silico citrullination (Model)
################################################################################
def citrullinate(
selection='all',
surface_cutoff=0,
show_clashes=0,
color_mod='',
color_base='',
hydrogens=0,
quiet=1
):


    '''
DESCRIPTION

    Modifies the Arginines of a selection to Citrullines

EXAMPLE

    frag ARG
    citrullinate

USAGE

    citrullinate [ selection [, surface_cutoff [, show_clashes [, color_mod [, color_base [, hydrogens [, quiet ]]]]]]]

ARGUMENTS

    selection: selection to be modified {default: 'all'}
               Arginines are automatically sub-selected!

    surface_cutoff: <float> variable for integrated selection of surface atoms.
                    If set >0, PyTMs will automatically calculate the solvent-
                    accessible surface area, and sub-select all atoms above this
                    cutoff.
                    Note that this operation will create a reserved selection:
                    'pytms_surface_input', but only modify residues if the target atoms
                    are also part of the original selection.
                    see also: PyTMs wiki page or findSurfaceResidues (PyMOL wiki)

    show_clashes: <boolean> toggles if clashes will be visualized {0}
    # consider removing potentially clashing heteroatoms

    color_base, color_mod: color names for
                           object selection and citrullination, respectively
                           {default: ''} = off

    hydrogens: <int>     toggles if the object will have hydrogens or not
                         after modification
                         0  : as is (detection) {default}
                         1  : adds hydrogens;
                         -1 : no hydrogens;

    quiet: <boolean> toggles output {default: quiet=1}
                     updates on operation and remaining time
                     * NB! print appears first in the console

    '''
    ##### BEGINNING OF SCRIPT #####

    last_time=start_time=datetime.datetime.now()
    try:
        surface_cutoff=abs(float(surface_cutoff))
    except:
        raise Exception("Input error!\n Illegal value for surface cutoff!")
    selection=pytms_get_sele(selection, surface_cutoff)

    # argument settings
    try:
        color_base, color_mod = str(color_base), str(color_mod)
        if color_base: cmd.color(color_base,selection)
        if color_mod: cmd.color(color_mod,''+selection+' and name O7')
        # transform boolean to reduce checking loops
        hydrogens=eval_hydrogens(selection, hydrogens)
        show_clashes=bool(int(show_clashes))
        output=bool(not int(quiet))
    except:
        raise Exception("Input error!\n Please check the input parameters!")

    # SELECTION LIST ###########################################################

    # create a empty list for appending
    stored.resi_list = []
    stored.temp = []

    # for each object and chain --> create selection lists
    for p in names:
        # for each object and chain get resi of arginine residues
        cmd.iterate((
        '(%s) and (resn ARG) '
        'and (name CA)' %(p)),
        'stored.resi_list.append("(%s and resi "+str(resi)+")")'  %p
        )

    # kick out non-selected or missing
    for p in list(stored.resi_list):
        if (cmd.count_atoms('(%s) and (%s) and name CZ' %(selection,p))==0):
            stored.resi_list.remove(p)
    ##### finished creating selection lists #####

    # premature exit if nothing to process
    if stored.resi_list==[]:
        print("PyTMs: No modifyable residues found in selection!")
        return False

    # VDW PREP
    # OBJCHI (dictionary) structure
    # obj: [umod. strain, mod. strain.]
    OBJCHIS = {}

    if show_clashes:
        # get unmodified VDW strain
        log_pytms_prog(output, 'Calculating base VdW strain!')
        for p in objects:
            if (hydrogens):
                cmd.h_add(p)
            else:
                cmd.remove('((%s) and elem H)' %p)

            OBJCHIS[p] = [get_strain('(%s)'%(p), temp_clash_A, 0), 0]

    last_time=log_pytms_prog(output, 'Initialized citrullination!')

    ##### CIT BUILDER #####

    # cycles through residues
    listcount=-1
    for p in stored.resi_list:
        listcount=listcount+1
        cmd.unpick()

        # get info on residue being processed
        stored.residue=''
        cmd.iterate('%s and name CA' %p, 'stored.residue=str(resn)')

        # build citrulline
        cmd.remove('%s and (((name NH1) extend 1) and not (name CZ))' %p)
        cmd.edit('%s and name CZ' %p)
        cmd.attach('O','2','2')
        cmd.alter('%s and (neighbor (name CZ)) and (elem O)' %p,'name="O7"')
        cmd.unbond('pk1', '%s and name O7' %p)
        cmd.bond('pk1', '%s and name O7' %p,'2')

        ## change name
        cmd.alter(p,'resn="CIR"')

        #Log
        last_time=log_pytms_prog(output,
        'Modified: %s'%get_resi_macro_name(p),listcount+1,len(stored.resi_list))
        cmd.unpick()
    # End of resi cycle
    if (hydrogens):
        for p in objects: cmd.h_add(p)
    else:
        for p in objects:
            cmd.remove('((%s) and elem H)' %p)

    # clashes
    if show_clashes:
        # get modified strain
        log_pytms_prog(output, 'Calculating modified VdW strain!')
        if output:
            print('STRAIN REPORT:')
            print('OBJECT','NATIVE_STRAIN','MODIFIED_STRAIN', 'DIFFERENCE')
        for p in objects:
            OBJCHIS[p][1] = get_strain('(%s)'%(p), temp_clash_A, 0)
            cmd.set_name(temp_clash_A,'%s_clashes'%p)
            if output:
                print('%s'%p,OBJCHIS[p][0],OBJCHIS[p][1],OBJCHIS[p][1]-OBJCHIS[p][0])
    # rebuild
    pytms_rebuild()

    # coloring
    if color_base:
        cmd.color(color_base,
        '%s or ((%s extend 1) and hydrogens)' %(selection, selection))
    for p in stored.resi_list:
        # enable selection "p.PTM in citrullination"
        if version_ok:
            cmd.alter('(%s and name O7)' %p,
            'p.PTM="citrullination"')
        if color_base:
            cmd.color(color_base,
            '(%s)' %(p))
        if color_mod:
            cmd.color(color_mod,
            '(%s and (name O7))' %(p))

    #exit
    last_time=log_pytms_prog(output, 'Citrullination done!')
    return [stored.resi_list, OBJCHIS]

cmd.extend( "citrullinate", citrullinate );

################################################################################
################################################################################






################################################################################
# Automated in silico carbamylation (Model)
################################################################################
def carbamylate(
selection='all',
surface_cutoff=0,
position=0,
show_clashes=0,
color_mod='',
color_base='',
hydrogens=0,
quiet=1
):


    '''
DESCRIPTION

    Modifies the Lysines of a selection to Carbamyl-Lysines (Homocitrullines)

EXAMPLE

    frag LYS
    carbamylate

USAGE

    carbamylate [ selection [, surface_cutoff [, position [, show_clashes
    [, color_mod [, color_base [, hydrogens [, quiet ]]]]]]]]

ARGUMENTS

    selection: selection to be modified {default: 'all'}
               Arginines are automatically sub-selected!

    surface_cutoff: <float> variable for integrated selection of surface atoms.
                    If set >0, PyTMs will automatically calculate the solvent-
                    accessible surface area, and sub-select all atoms above this
                    cutoff.
                    Note that this operation will create a reserved selection:
                    'pytms_surface_input', but only modify residues if the target atoms
                    are also part of the original selection.
                    see also: PyTMs wiki page or findSurfaceResidues (PyMOL wiki)

    position: <int> toggles the sites to be modified
                0: Lysines only {default}
                1: N-termini only
                2: Lysines and N-termini

    show_clashes: <boolean> toggles if clashes will be visualized {0}
    # consider removing potentially clashing heteroatoms

    color_base, color_mod: color names for
                           object selection and modification, respectively
                           {default: ''} = off

    hydrogens: <int>     toggles if the object will have hydrogens or not
                         after modification
                         0  : as is (detection) {default}
                         1  : adds hydrogens;
                         -1 : no hydrogens;

    quiet: <boolean> toggles output {default: quiet=1}
                     updates on operation and remaining time
                     * NB! print appears first in the console

    '''
    ##### BEGINNING OF SCRIPT ##################################################

    last_time=start_time=datetime.datetime.now()
    try:
        surface_cutoff=abs(float(surface_cutoff))
    except:
        raise Exception("Input error!\n Illegal value for surface cutoff!")
    selection=pytms_get_sele(selection, surface_cutoff)

    # argument settings
    try:
        color_base, color_mod = str(color_base), str(color_mod)
        if color_base: cmd.color(color_base,selection)
        if color_mod: cmd.color(color_mod,
        (''+selection+' and ((name OC1) or '+
        '((name NC1) extend 1) or (name CC1))'))
        position=abs(int(position))
        if (position not in range(0,3)):
            raise Exception("position out of range!")
        # transform boolean to reduce checking loops
        show_clashes=bool(int(show_clashes))
        hydrogens=eval_hydrogens(selection, hydrogens)
        output=bool(not int(quiet))
    except:
        raise Exception("Input error!\n Please check the input parameters!")

    ##### SELECTION LIST #######################################################

    # create a empty list for appending
    stored.resi_list = []
    stored.temp_k = []
    stored.temp_n = []
    instance_list = []

    # for each object and chain --> create selection lists
    for p in names:
        if (position!=1):
            # for each object and chain get resi of lysine residues
            cmd.iterate((
            '(%s) and ((resn LYS) or (resn LYN)) '
            'and (name NZ)' %(p)),
            'stored.temp_k.append("(%s and resi "+str(resi)+")")'  %p
            )

        # only append N-term if selected
        if (position!=0):
            # for each object and chain get resi of N-term
            # NB! not hetatm was used instead of polymer to allow 'frag'
            # amino acids to be selected
            cmd.iterate(
            '(byres (first %s)) and (not (hetatm)) '
            'and (name N) and not ((resn PRO) or (resn HYP))' %(p),
            'stored.temp_n.append("(%s and resi "+str(resi)+")")'  %p
            )
            # verify real N-terminus
            for c in list(stored.temp_n):
                if (cmd.count_atoms(('(neighbor (%s and (name N))) '
                'and (elem C)'%c))!=1):
                    stored.temp_n.remove(c)

    # now there are lists for K and N-terms which will be cleaned and merged
    # kick out non-selected or missing (not in selection or already modified)
    for p in list(stored.temp_k):
        if (cmd.count_atoms('(%s) and (%s) and (name NZ)' %(selection,p))==0):
            stored.temp_k.remove(p)
        elif (cmd.count_atoms('(%s) and elem C and (neighbor (name NZ))' %(p))>1):
            # modified otherwise
            stored.temp_k.remove(p)
    for p in list(stored.temp_n):
        if (cmd.count_atoms('(%s) and (%s) and (name N)' %(selection,p))==0):
            stored.temp_n.remove(p)
        elif (cmd.count_atoms(('(neighbor (%s and (name N))) '
            'and (elem C)'%p))!=1):
            # modified otherwise
            stored.temp_n.remove(p)

    # instance will be
    # 4: non n-terminal Lysine (target: NZ)
    # 3: n-terminal Lysine (target: NZ)
    # 2: n-terminal Lysine (target: N)
    # 1: n-term (target: N)
    # 0: skipped (should not occur)

    instance_list = list(stored.temp_k)
    for p in range(0,len(stored.temp_k)):
        # overwrite entry
        instance_list[p]=0
        if (stored.temp_k[p] in stored.temp_n):
            # n-terminal K
            instance_list[p]=3
        else:
            # non n-terminal K
            instance_list[p]=4
        stored.resi_list.append(stored.temp_k[p])

    for p in range(0,len(stored.temp_n)):
        if (stored.temp_n[p] in stored.temp_k):
            # n-terminal K
            instance_list.append(2)
        else:
            # non n-terminal K
            instance_list.append(1)
        stored.resi_list.append(stored.temp_n[p])

    ##### finished creating selection lists #####

    # premature exit if nothing to process
    if stored.resi_list==[]:
        print("PyTMs: No modifyable residues found in selection!")
        return False

    # VDW PREP
    # OBJCHI (dictionary) structure
    # obj: [umod. strain, mod. strain.]
    OBJCHIS = {}

    if show_clashes:
        # get unmodified VDW strain
        log_pytms_prog(output, 'Calculating base VdW strain!')
        for p in objects:
            if (hydrogens):
                cmd.h_add(p)
            else:
                cmd.remove('((%s) and elem H)' %p)

            OBJCHIS[p] = [get_strain('(%s)'%(p), temp_clash_A, 0),0]

    last_time=log_pytms_prog(output, 'Initialized carbamylation!')

    # cycles through residues (first lysines then N-terms) #####################
    ##### PRE_SETTINGS #####
    listcount=-1
    for p in list(stored.resi_list):
        cmd.unpick()
        listcount=listcount+1
        #jump skipped
        if instance_list[listcount]==0:continue

        #pick starting Nitrogen
        if (instance_list[listcount] in [3,4]):
            # handled as lysine
            cmd.remove('(neighbor (%s and name NZ)) and hydrogens' %p)
            cmd.edit('%s and name NZ' %p)
        else:
            # handled as n-term
            cmd.remove('(neighbor (%s and name N)) and hydrogens' %p)
            cmd.edit('%s and (name N)' %p)

        # adjust (assumed) geometry and valence of nitrogen
        cmd.set_geometry('pk1', 3, 3)
        cmd.alter('pk1', 'formal_charge=0')

        # CARB BUILDER #########################################################

        # build carbamyl (pick starting atom)
        cmd.remove('%s and (neighbor (pk1)) and hydrogens' %p)

        cmd.attach('C','3','3')
        # temporary unique name for selection (in case of double modification!)
        cmd.alter(('(%s and (neighbor (pk1)) and '
                   '(elem C)) and (not (name CE+CA))' %p),
        'name="CXXX"')

        cmd.unbond('pk1', '%s and name CXXX' %p)
        cmd.bond('pk1', '%s and name CXXX' %p,2)

        cmd.edit('%s and name CXXX' %p)
        cmd.attach('O','2','2')
        cmd.alter('(%s and (neighbor (name CXXX)) and elem O)' %p,
        'name="OC1"')

        select_part(temp_sel,'%s'%p,'pk1','(name CB)')
        # fix bonds used for temporary sp2
        cmd.unbond('%s and name NZ' %temp_sel, '%s and name CXXX' %temp_sel)
        cmd.bond('%s and name NZ' %temp_sel, '%s and name CXXX' %temp_sel,1)
        cmd.unbond('%s and name CXXX' %temp_sel, '%s and name OC1' %temp_sel)
        cmd.bond('%s and name CXXX' %temp_sel, '%s and name OC1' %temp_sel,2)

        cmd.attach('N','3','3')
        cmd.alter(('(%s and (neighbor (pk1)) and '
                   '(elem N)) and (not (name N+NZ))' %p),
        'name="NC1"')

        select_part(temp_sel,'%s'%p,'pk1','(name CB)')

        cmd.edit('%s and name NC1' %temp_sel)
        cmd.attach('H','1','1')
        cmd.attach('H','1','1')

        # rename CXXX
        cmd.alter('(%s and (name CXXX))'%temp_sel,
        'name="CC1"')

        ## change name if LYS
        if instance_list[listcount] in  [3,4]:
            cmd.alter(p,'resn="CAK"')
        else:
            select_part(temp_sel,'%s'%p,'(name CA)','(name CB)')
            # adjust resn, resi and resv
            stored.temp_n=''
            cmd.iterate('(%s) and name N' %(temp_sel),
            'stored.temp_n=int(resv)')
            stored.temp_n=stored.temp_n-1
            cmd.alter(('(%s and ((name NC1) extend 2)) '
            'and not (name N)' %temp_sel),
            'resn="NCA"')
            #N-Carbamyl
            cmd.alter(('(%s and ((name NC1) extend 2)) '
            'and not (name N)' %temp_sel),
            'resv=%s' %stored.temp_n)
            cmd.alter(('(%s and ((name NC1) extend 2)) '
            'and not (name N)' %temp_sel),
            'resi=%s' %stored.temp_n)
            # replace resi_list entry
            stored.resi_list[listcount]=str(
            "((model "+str(cmd.get_object_list(temp_sel)[0])+
            " and chain '"+str((cmd.get_chains(temp_sel) or [''])[0])+
            "' and resi "+str(stored.temp_n)+") extend 3) "+
            "and not ((name CB) or ((neighbor (name CB)) and (elem H)))"
            )

        #Log
        last_time=log_pytms_prog(output,
        'Modified: %s'%get_resi_macro_name(p),listcount+1,len(stored.resi_list))

        cmd.unpick()
    # End of resi cycle ########################################################

    # fixes ordering of resis
    for p in objects:
        cmd.remove('((%s) and elem H)' %p)
        cmd.h_add(p)

    # reinstate hydrogens
    if (hydrogens):
        for p in objects: cmd.h_add(p)
    else:
        for p in objects:
            cmd.remove('((%s) and elem H)' %p)

    # clashes
    if show_clashes:
        # get modified strain
        log_pytms_prog(output, 'Calculating modified VdW strain!')
        if output:
            print('STRAIN REPORT:')
            print('OBJECT','NATIVE_STRAIN','MODIFIED_STRAIN', 'DIFFERENCE')
        for p in objects:
            OBJCHIS[p][1] = get_strain('(%s)'%(p), temp_clash_A, 0)
            cmd.set_name(temp_clash_A,'%s_clashes'%p)
            if output:
                print('%s'%p,OBJCHIS[p][0],OBJCHIS[p][1],OBJCHIS[p][1]-OBJCHIS[p][0])

    # rebuild
    pytms_rebuild()

    # coloring
    if color_base:
        cmd.color(color_base,
        '%s or ((%s extend 1) and hydrogens)' %(selection, selection))
    for p in stored.resi_list:
        # allow selection "p.PTM in carbamylation"
        if version_ok:
            cmd.alter(('(%s and ((name OC1) or '+
            '((name NC1) extend 1) or (name CC1)))') %p,
            'p.PTM="carbamylation"')
        if color_base:
            cmd.color(color_base,
            '(%s)' %(p))
        if color_mod:
            cmd.color(color_mod,
            ('(%s and ((name OC1) or '+
            '((name NC1) extend 1) or (name CC1)))') %p)

    # exit
    last_time=log_pytms_prog(output, 'Carbamylation done!')
    return [stored.resi_list, OBJCHIS]

cmd.extend( "carbamylate", carbamylate );
################################################################################
################################################################################





################################################################################
# Automated in silico 3-Nitrotyrosinylation (Model)
################################################################################
def nitrate(
selection='all',
surface_cutoff=0,
mode=0,
position_tyr=1,
position_trp=6,
angle=22.3519420624,
show_clashes=0,
color_mod='',
color_base='',
delocalized=1,
hydrogens='',
quiet=1
):


    '''
DESCRIPTION

    Modifies the Tyrosines (and/or Tryptophans, depending on mode)
    of a selection by nitration

EXAMPLE

    frag TYR
    nitrate

USAGE

    nitrate [ selection [, surface_cutoff [, mode
    [, position_tyr [, position_trp [, show_clashes
    [, color_mod [, color_base [, delocalized [, hydrogens
    [, quiet ]]]]]]]]]]]


ARGUMENTS

    selection: selection to be modified {default: 'all'}
               Tyrosines (and/or Tryptophans) are automatically sub-selected!

    surface_cutoff: <float> variable for integrated selection of surface atoms.
                    If set >0, PyTMs will automatically calculate the solvent-
                    accessible surface area, and sub-select all atoms above this
                    cutoff.
                    Note that this operation will create a reserved selection:
                    'pytms_surface_input', but only modify residues if the target atoms
                    are also part of the original selection.
                    see also: PyTMs wiki page or findSurfaceResidues (PyMOL wiki)

    mode: <int> Toggles the residues to be nitrated
                0: Tyrosines only {default}
                1: Tryptophans only
                2: Both

    position_tyr: <int> toggles the ortho-atom to be modified for TYR
                    0: random/ per residue (if surface_cutoff is set >0 it will select the most accessible)
                    1: CE1 (position 3) {default}
                    2: CE2 (position 3)


    position_trp: <int> toggles the carbon to be modified for TRP
                    0: random/ per residue (excluding 1)
                    4: CE3
                    5: CZ3
                    6: CH2  {default}
                    7: CZ2
                    1: CD1 (irrelevant)

    angle: <float> torsion angle for CZ-CE(1/2)-NN-O1, i.e. the rotation of the nitro group
                    note that negating this vlaue may yield different results,
                    also depending on the choice of CE1 vs CE2
                    the default is: 22.3519420624 as taken from PDB: 4NDA

    show_clashes: <boolean> toggles if clashes will be visualized {0}
    # consider removing potentially clashing heteroatoms

    color_base, color_mod: color names for
                           object selection and modification, respectively
                           {default: ''} = off

    delocalized: <boolean> toggles whether the bonds in the nitro group
                           are delocalized or not {default: 1}

    hydrogens: <boolean> toggles if the object will have hydrogens or not
                         after modification
                         ''    : as is (detection);
                         False : no hydrogens;
                         any other entry=True : hydrogens
                         {default: ''}

    quiet: <boolean> toggles output {default: quiet=1}
                     updates on operation and remaining time
                     * NB! print appears first in the console

    '''
    ##### BEGINNING OF SCRIPT #####

    last_time=start_time=datetime.datetime.now()
    try:
        surface_cutoff=abs(float(surface_cutoff))
    except:
        raise Exception("Input error!\n Illegal value for surface cutoff!")
    selection=pytms_get_sele(selection, surface_cutoff)

    # argument settings
    try:
        mode=abs(int(mode))
        position_tyr=abs(int(position_tyr))
        position_trp=abs(int(position_trp))
        angle=float(angle)
        color_base, color_mod = str(color_base), str(color_mod)
        if color_base: cmd.color(color_base,selection)
        if color_mod: cmd.color(color_mod,
        ''+selection+' and ((name NN) or (name O1) or (name O2))')
        # transform boolean to reduce checking loops
        show_clashes=bool(int(show_clashes))
        delocalized=bool(int(delocalized))
        hydrogens=eval_hydrogens(selection, hydrogens)
        output=bool(not int(quiet))
    except:
        raise Exception("Input error!\n Please check the input parameters!")

    # range check
    if (mode not in range(0,3)):
        raise Exception("mode out of range!")
    if (position_tyr not in range(0,3)):
        raise Exception("position_tyr out of range!")
    if (position_trp not in [0,4,5,6,7,1]):
        raise Exception("position_trp out of range!")

    ##### SELECTION LIST #####

    # create a empty list for appending
    stored.resi_list = []
    stored.temp = []

    # for each object and chain --> create selection lists
    for p in names:
        # for each object and chain get resi of residues
        cmd.iterate((
        '(%s) and ((resn TYR) or (resn PTR) '
        'or (resn TRP)) and (name CA)' %(p)),
        'stored.resi_list.append("(%s and resi "+str(resi)+")")'  %p
        )
    # kick out non-selected or missing
    for p in list(stored.resi_list):
        if (cmd.count_atoms('(%s) and (%s) and ((name CE1) or (name CH2))' %(selection,p))==0):
            stored.resi_list.remove(p)
    for p in list(stored.resi_list):
        # get info on residue being processed
        stored.residue=''
        cmd.iterate('%s and name CA' %p, 'stored.residue=str(resn)')
        if ((mode==0) and (stored.residue in ['TRP'])):
            # W-based AA, but only Y processed
            stored.resi_list.remove(p)
        if ((mode==1) and (stored.residue in ['PTR', 'TYR'])):
            # is Y-based AA, but only W is processed
            stored.resi_list.remove(p)

    ##### finished creating selection lists #####

    # premature exit if nothing to process:
    if stored.resi_list==[]:
        print("PyTMs: No modifyable residues found in selection!")
        return False

    # VDW PREP
    # OBJCHI (dictionary) structure
    # obj: [umod. strain, mod. strain.]
    OBJCHIS = {}

    if show_clashes:
        # get unmodified VDW strain
        log_pytms_prog(output, 'Calculating base VdW strain!')
        for p in objects:
            if (hydrogens):
                cmd.h_add(p)
            else:
                cmd.remove('((%s) and elem H)' %p)

            OBJCHIS[p] = [get_strain('(%s)'%(p), temp_clash_A, 0),0]

    last_time=log_pytms_prog(output, 'Initialized nitration!')


    ##### NIT BUILDER #####

    # cycles through residues
    listcount=-1
    for p in stored.resi_list:
        listcount=listcount+1
        cmd.unpick()

        # get info on residue being processed
        stored.residue=''
        cmd.iterate('%s and name CA' %p, 'stored.residue=str(resn)')

        ## change name
        if (stored.residue=='PTR'):
            # is phosphorylated
            cmd.alter(p,'resn="PNIY"')
        elif (stored.residue=='TYR'):
            # regular Tyrosine
            cmd.alter(p,'resn="NIY"')
        elif (stored.residue=='TRP'):
            # W
            cmd.alter(p,'resn="NIW"')
        else:
            # unknown
            continue

        # W or Y?
        if (cmd.count_atoms('%s and name NE1'%p)==1):
            #W (==1)
            posname='CH2'
            if position_trp==1: posname='CD1'
            if position_trp==4: posname='CE3'
            if position_trp==5: posname='CZ3'
            if position_trp==6: posname='CH2'
            if position_trp==7: posname='CZ2'
            if position_trp==0:
                posname=random.choice(['CE3','CZ3','CH2','CZ2'])
        else:
            #Y (==0)
            posname='CE1'
            if position_tyr==1: posname='CE1'
            if position_tyr==2: posname='CE2'
            if position_tyr==0:
                if surface_cutoff>0:
                    stored.b1=0
                    stored.b2=0
                    cmd.iterate('%s and (name CE1)' %(p), '"stored.b1=b"')
                    cmd.iterate('%s and (name CE2)' %(p), '"stored.b2=b"')
                    posname='CE1'
                    if stored.b1<stored.b2: posname='CE2'
                else:
                    posname=random.choice(['CE1','CE2'])

        # build Nitro group
        cmd.remove('%s and ((neighbor (name %s)) and elem H)' %(p,posname))
        cmd.edit('%s and (name %s)' %(p,posname))
        cmd.attach('N','3','4')
        cmd.alter(('%s and ((neighbor (name %s)) and (elem N) '
        'and not (name NE1))' %(p,posname)),
        'name="NN"')
        cmd.edit('%s and (name NN)' %p)
        cmd.attach('O','2','2')
        cmd.alter('%s and (neighbor (name NN)) and elem O' %p,'name="O1"')
        cmd.attach('O','2','2')
        cmd.alter(('%s and (neighbor (name NN)) '
                   'and (elem O) and not (name O1)' %p),
        'name="O2"')
        cmd.unbond('%s and name NN' %p, '%s and name O1' %p)
        cmd.unbond('%s and name NN' %p, '%s and name O2' %p)

        if delocalized:
            cmd.bond('%s and name NN' %p, '%s and name O1' %p,'4')
            cmd.bond('%s and name NN' %p, '%s and name O2' %p,'4')
        else:
            cmd.bond('%s and name NN' %p, '%s and name O1' %p,'2')
            cmd.bond('%s and name NN' %p, '%s and name O2' %p,'1')

        #fix angle in tyrosines
        if (cmd.count_atoms('%s and name NE1'%p)==0):
            #Y (==1)
            cmd.set_dihedral(
            '%s and name CZ' %p,
            '%s and name %s' %(p,posname),
            '%s and name NN' %p,
            '%s and name O1' %p,
            angle) # default value from PDB: 4NDA

        #Log
        last_time=log_pytms_prog(output,
        'Modified: %s'%get_resi_macro_name(p),listcount+1,len(stored.resi_list))

        cmd.unpick()
    # End of resi cycle
    if (hydrogens):
        for p in objects: cmd.h_add(p)
    else:
        for p in objects:
            cmd.remove('((%s) and elem H)' %p)

    # clashes
    if show_clashes:
        # get modified strain
        log_pytms_prog(output, 'Calculating modified VdW strain!')
        if output:
            print('STRAIN REPORT:')
            print('OBJECT','NATIVE_STRAIN','MODIFIED_STRAIN', 'DIFFERENCE')
        for p in objects:
            OBJCHIS[p][1] = get_strain('(%s)'%(p), temp_clash_A, 0)
            cmd.set_name(temp_clash_A,'%s_clashes'%p)
            if output:
                print('%s'%p,OBJCHIS[p][0],OBJCHIS[p][1],OBJCHIS[p][1]-OBJCHIS[p][0])

    # rebuild
    pytms_rebuild()

    # coloring
    if color_base:
        cmd.color(color_base,
        '%s or ((%s extend 1) and hydrogens)' %(selection, selection))
    for p in stored.resi_list:
        # enable selection "p.PTM in nitration"
        if version_ok:
            cmd.alter('((name NN) or (name O1) or (name O2)) and %s' %p,
            'p.PTM="nitration"')
        if color_base:
            cmd.color(color_base,
            '(%s)' %(p))
        if color_mod:
            cmd.color(color_mod,
            '(%s and ((name NN) or (name O1) or (name O2)))' %(p))

    last_time=log_pytms_prog(output, 'Nitration done!')
    return [stored.resi_list, OBJCHIS]

cmd.extend( "nitrate", nitrate );

################################################################################
################################################################################





################################################################################
# Automated in silico Cysteine S-Nitrosylation (Model)
################################################################################
def nitrosylate(
selection='all',
surface_cutoff=0,
include_SEC=0,
disulfides=0,
show_clashes=0,
color_mod='',
color_base='',
hydrogens=0,
quiet=1
):


    '''
DESCRIPTION

    Modifies the Cysteins of a selection by S-Nitrosylation

EXAMPLE

    frag CYS
    nitrosylate

USAGE

    nitrosylate [ selection [, surface_cutoff [, include_SEC [, disulfides
    [, show_clashes [, color_mod [, color_base [, hydrogens [, quiet ]]]]]]]]]

ARGUMENTS

    selection: selection to be modified {default: 'all'}
               Cysteines etc. are automatically sub-selected
               (see also: include_SEC and disulfides)!

    surface_cutoff: <float> variable for integrated selection of surface atoms.
                    If set >0, PyTMs will automatically calculate the solvent-
                    accessible surface area, and sub-select all atoms above this
                    cutoff.
                    Note that this operation will create a reserved selection:
                    'pytms_surface_input', but only modify residues if the target atoms
                    are also part of the original selection.
                    see also: PyTMs wiki page or findSurfaceResidues (PyMOL wiki)

    include_SEC: <boolean> toggles if the selenocysteines (SEC/CSE) will be included or not
                 {default: 0}

    disulfides: <boolean> toggles if disulfide bridges will be affected or not
                {default: 0} NB! setting this True may give erroneous results

    show_clashes: <boolean> toggles if clashes will be visualized {0}
    # consider removing potentially clashing heteroatoms

    color_base, color_mod: color names for
                           object selection and modification, respectively
                           {default: ''} = off

    hydrogens: <int>     toggles if the object will have hydrogens or not
                         after modification
                         0  : as is (detection) {default}
                         1  : adds hydrogens;
                         -1 : no hydrogens;

    quiet: <boolean> toggles output {default: quiet=1}
                     updates on operation and remaining time
                     * NB! print appears first in the console

    '''
    ##### BEGINNING OF SCRIPT #####

    last_time=start_time=datetime.datetime.now()
    try:
        surface_cutoff=abs(float(surface_cutoff))
    except:
        raise Exception("Input error!\n Illegal value for surface cutoff!")
    selection=pytms_get_sele(selection, surface_cutoff)

    # argument settings
    try:
        color_base, color_mod = str(color_base), str(color_mod)
        if color_base: cmd.color(color_base,selection)
        if color_mod:
            cmd.color(color_mod,
            ('(%s and (resn CNO)) and '
            '((name NO) or (name ON))' %selection))
        # transform boolean to reduce checking loops
        show_clashes=bool(int(show_clashes))
        include_SEC=bool(int(include_SEC))
        disulfides=bool(int(disulfides))
        hydrogens=eval_hydrogens(selection, hydrogens)
        output=bool(not int(quiet))
    except:
        raise Exception("Input error!\n Please check the input parameters!")

    # SELECTION LIST ###########################################################

    # rename to allow later conversion
    if include_SEC:
        cmd.alter('byres (%s and %s and (resn SCE or resn SEC))' %(p,selection),
        'resn="CYS"')

    # create a empty list for appending
    stored.resi_list = []
    stored.temp = []

    # for each object and chain --> create selection lists
    for p in names:
        # for each object and chain get resi of residues
        cmd.iterate((
        '(%s) and (resn CYS) '
        'and (name CA)' %(p)),
        'stored.resi_list.append("(%s and resi "+str(resi)+")")'  %p
        )

    # kick out non-selected or missing
    for p in list(stored.resi_list):
        if (cmd.count_atoms('(%s) and (%s) and (elem S or elem SE)' %(selection,p))==0):
            stored.resi_list.remove(p)

    # kick out disulfides (if not set!)
    if not disulfides:
        for p in list(stored.resi_list):
            if (cmd.count_atoms((
            '(neighbor ((%s) and (elem S or elem SE))) '
            'and (elem S or elem SE)' %(p)))!=0):
                # not regular cysteine
                stored.resi_list.remove(p)


    ##### finished creating selection lists #####

    # premature exit if nothing to process:
    if stored.resi_list==[]:
        print("PyTMs: No modifyable residues found in selection!")
        return False

    # VDW PREP
    # OBJCHI (dictionary) structure
    # obj: [umod. strain, mod. strain.]
    OBJCHIS = {}

    if show_clashes:
        # get unmodified VDW strain
        log_pytms_prog(output, 'Calculating base VdW strain!')
        for p in objects:
            if (hydrogens):
                cmd.h_add(p)
            else:
                cmd.remove('((%s) and elem H)' %p)

            OBJCHIS[p] = [get_strain('(%s)'%(p), temp_clash_A, 0),0]

    last_time=log_pytms_prog(output,
    'Initialized Cysteine S-Nitrosylation!')

    ##### BUILDER #####

    # cycles through residues
    listcount=-1
    for p in stored.resi_list:
        cmd.unpick()
        listcount=listcount+1

        # get info on residue being processed
        stored.residue=''
        cmd.iterate('%s and name CA' %p, 'stored.residue=str(resn)')


        ## change name
        cmd.alter(p,'resn="CNO"')

        # get rid of bound non-chain neighbors
        # but conserve potential S
        cmd.remove(('%s and (neighbor (elem SE or elem S)) '
        'and not (elem C or elem S)' %p))
        # edit
        cmd.edit('%s and (elem SE or elem S)' %p)
        cmd.attach('N','2','3')
        cmd.alter('%s and (neighbor (elem SE or elem S)) and (elem N)' %p,'name="NO"')
        cmd.edit('%s and (name NO)' %p)
        cmd.attach('O','2','2')
        cmd.alter('%s and (neighbor (name NO)) and (elem O)' %p,'name="ON"')
        cmd.unbond('pk1', '%s and name ON' %p)
        cmd.bond('pk1', '%s and name ON' %p,'2')

        #Log
        last_time=log_pytms_prog(output,
        'Modified: %s'%get_resi_macro_name(p),listcount+1,len(stored.resi_list))

        cmd.unpick()
    # End of resi cycle
    if (hydrogens):
        for p in objects: cmd.h_add(p)
    else:
        for p in objects:
            cmd.remove('((%s) and elem H)' %p)

    # clashes
    if show_clashes:
        # get modified strain
        log_pytms_prog(output, 'Calculating modified VdW strain!')
        if output:
            print('STRAIN REPORT:')
            print('OBJECT','NATIVE_STRAIN','MODIFIED_STRAIN', 'DIFFERENCE')
        for p in objects:
            OBJCHIS[p][1] = get_strain('(%s)'%(p), temp_clash_A, 0)
            cmd.set_name(temp_clash_A,'%s_clashes'%p)
            if output:
                print('%s'%p,OBJCHIS[p][0],OBJCHIS[p][1],OBJCHIS[p][1]-OBJCHIS[p][0])

    # rebuild
    pytms_rebuild()

    # coloring
    if color_base:
        cmd.color(color_base,
        '%s or ((%s extend 1) and hydrogens)' %(selection, selection))
    for p in stored.resi_list:
        # enable selection "p.PTM in s_nitrosylation"
        if version_ok:
            cmd.alter('%s and ((name NO) or (name ON))' %p,
            'p.PTM="s_nitrosylation"')
        if color_base:
            cmd.color(color_base,
            '(%s)' %(p))
        if color_mod:
            cmd.color(color_mod,
            ('(%s and (resn CNO)) and '
            '((name NO) or (name ON))' %p))
    #exit
    last_time=log_pytms_prog(output,
    'Cysteine S-Nitrosylation complete!')
    return [stored.resi_list, OBJCHIS]

cmd.extend( "nitrosylate", nitrosylate );
################################################################################
################################################################################



################################################################################
# Automated in silico MDA modification (Model)
################################################################################
def mda_modify(
selection='all',
surface_cutoff=0,
position=0,
group=2,
confomer=5,
type=3,
torsions=[106.5,180,0,-33,33,180,-180,0,120,60],
angles=[118.35,120,120,110.5,120],
protonate=-1,
optimize=1,
interpolate=0,
local_radius=10,
base_strain_limit=300,
intervals=[45,45,45,45,45],
states=0,
remove_radius=5,
optimize_ignore='hetatm',
optimize_add='none',
show_clashes=0,
color_mod='',
color_base='',
hydrogens=0,
quiet=1
):


    '''
DESCRIPTION

    Modifies the free amines of a selection with Malondialdehyde (MDA)

EXAMPLE

    frag LYS
    mda_modify

    frag ARG
    mda_modify position=1

USAGE

    mda_modify [ selection [, surface_cutoff [, position
    [, group [, confomer [, type [, torsions [, angles [, protonate
    [, optimize [, interpolate [, local_radius [, base_strain_limit
    [, intervals [, states [, remove_radius [, optimize_ignore [, optimize_add
    [, show_clashes [, color_mod [, color_base [, hydrogens
    [, quiet ]]]]]]]]]]]]]]]]]]]]]]]

NOTES

    * MDA-Modified Lysines change name to MMK (Malodialdehyde-modified K)
    * modified N-terminals are treated as additional residue
    * The MDA adducts are selectable by "p.PTM in MDA" or
      "(((name CM*) or (name OM*)) extend 1) and not (elem N)"

ARGUMENTS

    selection: selection to be modified {default: 'all'}
               Lysines are automatically sub-selected!

    surface_cutoff: <float> variable for integrated selection of surface atoms.
                    If set >0, PyTMs will automatically calculate the solvent-
                    accessible surface area, and sub-select all atoms above this
                    cutoff.
                    Note that this operation will create a reserved selection:
                    'pytms_surface_input', but only modify residues if the target atoms
                    are also part of the original selection.
                    see also: PyTMs wiki page or findSurfaceResidues (PyMOL wiki)

    position: <int> toggles the sites to be modified {default=0}
                0: Lysines only
                1: N-termini only
                2: both Lysines and N-termini

    group, confomer, type:  <int> define MDA adduct type and positioning, see below:
    #defaults: group=2; confomer=5; type=3; (MAA adduct/ trans/ beta)
        group: group of adduct
                0: Random (once for all)
                1: Random (by residue)
                2: DHP (HDD or EDD)
                3: MDA
                4: FAAB

        confomer: DHP-TYPE                  MDA-TYPE                FAAB
                  0: Random (once for all)  0: Rand.                    0: Rand.
                  1: Random (by residue)    1: Rand.                    1: Rand.
                  2: cis-boat and alpha     2: R-N cis / C1-C2 cis      2: R
                  3: cis-boat and beta      3: R-N cis / C1-C2 trans    3: S
                  4: trans-boat and alpha   4: R-N trans / C1-C2 cis   (4: R)
                  5: trans-boat and beta    5: R-N trans / C1-C2 trans (5: S)

            type: MDD/EDD-TYPE                  MDA-TYPE
                  0: Random (once for all)  0: Rand.
                  1: Random (by residue)    1: Rand.
                  2: Hydro (HDD)            2: enamine, O in cis
                  3: Methyl (MDD)           3: enamine, O in trans
                  4: Ethenol, O in cis      4: imine-enol, O in cis
                  5: Ethenol, O in trans    5: imine-enol, O in trans
                  6: Ethanal, O in cis      6: imine-al, O in cis (pseudo)
                  7: Ethanal, O in trans    7: imine-al, O in trans (pseudo)

                  FAAB
                  0: Rand.
                  1: Rand.
                  2 or 5: delocalized bonds, closed
                  3 or 6: delocalized bonds, open
                  4 or 7: regular bonds (open)

    torsions/angles:
     <list> advanced setting for DHP-type (and FAAB) adducs,
            set the (starting) torsion/triangular angles for:
            # assume default DHP settings
            #torsions=[106.5,180,0,-33,33,180,-180,0,120,60]
            #angles=[118.35,120,120,110.5,120]
             TORSIONS                ANGLES
             [0]: CD-CE-NZ-CM2       [0]: CE-NZ-CM2
             [1]: CE-NZ-CM2-CM3      [1]: NZ-CM2-CM3
             [2]: NZ-CM2-CM3-CM4     [2]: CM2-CM3-CM4
             [3]: CM2-CM3-CM4-CM5    [3]: CM3-CM4-CM5
             [4]: CM3-CM4-CM5-CM6    [4]: CM4-CM5-CM6
             [5]: OMC3
             [6]: OMC5
             [7]: EDD: 4' adduct orientation
             [8]: FAAB moiety: orientation (CE-NZ-CM2-CM3)
             [9]: FAAB moiety: orientation (NZ-CM2-CM3-CM4A)
             ([1],[3],[4] will be auto-adjusted for cis boats)
             ([8] and [9] will be inverted for S-confomers)
             # note that the values are inter-dependent
             # changing them may distort the shape of the ring
             # for a planar MAA use, e.g.:
                 #torsions=[90,180,0,-0,0,180,-180,0,120,60]
                 #angles=[118.75,120,120,120,120]
                 #and confomer=4 or 5
                 #(become equivalent, whereas confomer 2 or 3 are invalid)

             # examine code for details

    protonate:  <int> toggles protonation of nitrogen
               {default : -1}
                -1: deprotonated (uncharged)
                0: no change (as is)
                1: forced protanation (charged)

    ##### Parameters below affect/take effect depending on 'optimize' ##########
    optimize:  <int> toggles VDW optimization by
               rotation of chi torsion angles
               (avoids clashing) NB! calculation time increases exponentially!
               ID : BOND
                0 : off
                1 : CE-NZ (CA-N for Nterm) {default}
                2 : CE-NZ, CD-CE
                3 : CE-NZ, CD-CE, CG-CD
                4 : CE-NZ, CD-CE, CG-CD, CB-CG
                5 : CE-NZ, CD-CE, CG-CD, CB-CG, CA-CB

    interpolate: <boolean>
                {default=0} toggles whether the script will
                attempt to average the angles obtained from equally minimal
                VDW strains. If off, the closest minimum will be used.

    local_radius: <float>
                        cutoff distance for local refinement
                       defines a local area around each modified residue
                       used for optimizing; {default=10}
                       smaller values speed up calculation,
                       but may result in misplacement
                       (increase if this is the case)

    base_strain_limit <float>:
        limits the max. baseline strain increase for the lysine residues during
        optimization, i.e. strain excluding MDA adduct. Rotamers with a base
        strain above this limit will not be considered for further optimization.
        Smaller values may accelerate calculation time, but may exclude
        potentially more optimal orientations for the adduct.
        (Note that if the adduct can only be optimized at the expense of high
        increase in base strain, it is either unlikely to occur or will impinge
        on the protein's conformation)
        {default: 300}
        # the average strain of all lysine rotamers is ~461 (range 0 to ~1444)
        # values>1500 will deactivate this option entirely

    intervals: <float> defines the step (in deg) for each angle in 'optimize':
                       lower values may yield better accuracy but will
                       increase the required time exponentially!
                       {default=[45,45,45,45,45]}

    states: <int> number of states used to animate the VDW optimization
                  only takes effect if "optimize" is performed
                  {default: 0} (set 0 or 1 for only final result)

    remove_radius: <float> toggles if heteroatoms are removed
               after optimization {default: 5.0}
               will remove all atoms (defined by optimize_ignore) surrounding
               modified residues (within <value> angstrom)
               Set to remove_radius=0 to supress any removal
               Optional: manually remove heteroatoms before/after running the
               script, if preferred (see also: optimize_ignore)

    optimize_ignore: <str> a string selection defining atoms to
                     be ignored during optimization {default: 'hetatm'}
                     can be used to ignore only specific heteroatoms
                     (e.g. 'resn HOH') or parts of objects
                     # Note that these will be removed, depending on
                     remove_radius

    optimize_add:   <str> a string selection defining additional atoms to
                    be accounted for during optimization other than the parent
                    object {default: 'none'}
                    (useful to account for adjacent objects)

    ##### Parameters above will affect/take effect depending on 'optimize' #####

    show_clashes: <boolean> toggles if clashes will be visualized {0}
    # consider removing potentially clashing heteroatoms

    color_base, color_mod: color names for
                           object selection and MDA adducts, respectively
                           {default: ''} = off

    hydrogens: <int>     toggles if the object will have hydrogens or not
                         after modification
                         0  : as is (detection) {default}
                         1  : adds hydrogens;
                         -1 : no hydrogens;
                         # will affect the VDW strain and optimization

    quiet: <boolean> toggles output {default: quiet=1}
                     updates on operation and remaining time
                     * NB! print appears first in the console
                     prints optimization values


    '''
    ##### BEGINNING OF SCRIPT #####

    last_time=start_time=datetime.datetime.now()
    try:
        surface_cutoff=abs(float(surface_cutoff))
    except:
        raise Exception("Input error!\n Illegal value for surface cutoff!")
    selection=pytms_get_sele(selection, surface_cutoff)

    # argument settings
    try:
        color_base, color_mod = str(color_base), str(color_mod)
        if color_base: cmd.color(color_base,selection)
        if color_mod: cmd.color(color_mod,
        ('('+selection+' and (((name CM*) or (name OM*)) '+
        'extend 1) and not (elem N))'))
        group=abs(int(group))
        confomer=abs(int(confomer))
        type=abs(int(type))
        torsions=eval(str(torsions))
        for c in range(0,10):
            torsions[c]=float(torsions[c])
        angles=eval(str(angles))
        for c in range(0,5):
            angles[c]=float(angles[c])
        local_radius=abs(float(local_radius))
        base_strain_limit=abs(float(base_strain_limit))
        states=abs(int(states))
        remove_radius=float(remove_radius)
        optimize_ignore=verify_sele(optimize_ignore)
        optimize_add=verify_sele(optimize_add)
        protonate=int(protonate)
        optimize=abs(int(optimize))
        intervals=eval(str(intervals))
        for c in range(0,5):
            intervals[c]=abs(float(intervals[c]))
        position=abs(int(position))
        # transform boolean to reduce checking loops
        interpolate=bool(int(interpolate))
        show_clashes=bool(int(show_clashes))
        hydrogens=eval_hydrogens(selection, hydrogens)
        output=bool(not int(quiet))
        if optimize<2: intervals[1]=361.0
        if optimize<3: intervals[2]=361.0
        if optimize<4: intervals[3]=361.0
        if optimize<5: intervals[4]=361.0
    except:
        raise Exception("Input error!\n Please check the input parameters!")

    # range check
    if (group not in range(0,5)):
        raise Exception("group out of range!")
    if (confomer not in range(0,6)):
        raise Exception("confomer out of range!")
    if (type not in range(0,8)):
        raise Exception("type out of range!")
    if (position not in range(0,3)):
        raise Exception("position out of range!")
    if (protonate not in [-1,0,1]):
        raise Exception("optimize out of range!")
    if (optimize not in range(0,6)):
        raise Exception("optimize out of range!")
    if (len(intervals)!=5):
        raise Exception("too many or too few intervals provided!")
    if (len(torsions)!=10):
        raise Exception("torsions out of range!")
    for c in range(0,10):
        if ((torsions[c]<-180.0) or (torsions[c]>180.0)):
            raise Exception("torsion value out of range!")
    if (len(angles)!=5):
        raise Exception("angles out of range!")

    # set G_ID, C_ID, T_ID and random modification (if == 0)
    # (if==1 --> defined later)
    G_ID, C_ID, T_ID = group, confomer, type
    if group==0:
        G_ID=random.choice([2,3,4])
    if confomer==0:
        C_ID=random.choice([2,3,4,5])
    if type==0:
        T_ID=random.choice([2,3,4,5,6,7])

    # SELECTION LIST ###########################################################

    # create a empty list for appending
    stored.resi_list = []
    stored.temp_k = []
    stored.temp_n = []
    instance_list = []

    # for each object and chain --> create selection lists
    for p in names:
        if (position!=1):
            # for each object and chain get resi of lysine residues
            cmd.iterate((
            '(%s) and ((resn LYS) or (resn LYN)) '
            'and (name NZ)' %(p)),
            'stored.temp_k.append("(%s and resi "+str(resi)+")")'  %p
            )

        # only append N-term if selected
        if (position!=0):
            # for each object and chain get resi of N-term
            # NB! not hetatm was used instead of polymer to allow 'frag'
            # amino acids to be selected
            cmd.iterate(
            '(byres (first %s)) and (not (hetatm)) '
            'and (name N) and not ((resn PRO) or (resn HYP))' %(p),
            'stored.temp_n.append("(%s and resi "+str(resi)+")")'  %p
            )
            # verify real N-terminus
            for c in list(stored.temp_n):
                if (cmd.count_atoms(('(neighbor (%s and (name N))) '
                'and (elem C)'%c))!=1):
                    stored.temp_n.remove(c)

    # now there are lists for K and N-terms which will be cleaned and merged
    # kick out non-selected or missing (not in selection or already modified)
    for p in list(stored.temp_k):
        if (cmd.count_atoms('(%s) and (%s) and (name NZ)' %(selection,p))==0):
            stored.temp_k.remove(p)
        elif (cmd.count_atoms('(%s) and elem C and (neighbor (name NZ))' %(p))>1):
            # modified otherwise
            stored.temp_k.remove(p)
    for p in list(stored.temp_n):
        if (cmd.count_atoms('(%s) and (%s) and (name N)' %(selection,p))==0):
            stored.temp_n.remove(p)
        elif (cmd.count_atoms(('(neighbor (%s and (name N))) '
            'and (elem C)'%p))!=1):
            # modified otherwise
            stored.temp_n.remove(p)

    # instance will be
    # 4: non n-terminal Lysine (target: NZ)
    # 3: n-terminal Lysine (target: NZ)
    # 2: n-terminal Lysine (target: N)
    # 1: n-term (target: N)
    # 0: skipped (should not occur)

    instance_list = list(stored.temp_k)
    for p in range(0,len(stored.temp_k)):
        # overwrite entry
        instance_list[p]=0
        if (stored.temp_k[p] in stored.temp_n):
            # n-terminal K
            instance_list[p]=3
        else:
            # non n-terminal K
            instance_list[p]=4
        stored.resi_list.append(stored.temp_k[p])

    for p in range(0,len(stored.temp_n)):
        if (stored.temp_n[p] in stored.temp_k):
            # n-terminal K
            instance_list.append(2)
        else:
            # non n-terminal K
            instance_list.append(1)
        stored.resi_list.append(stored.temp_n[p])


    ##### finished creating selection lists #####

    # premature exit if nothing to process:
    if stored.resi_list==[]:
        print("PyTMs: No modifyable residues found in selection!")
        return False

    if (output):
        print("----------------------------------------")
        print("--------- mda_modify output ------------")
        print("----------------------------------------")
    last_time=log_pytms_prog(output, 'Initialized MDA modification!')

    # VDW PREP
    # OBJCHI (dictionary) structure
    # obj: [umod. strain, mod. strain, opt. strain],
    # + [n*([Res,N-type,
    #       [[Start angles],[opt. angles],[diff]],
    #       [G_ID,C_ID,T_ID],[start local strain, optimal local strain]])]
    OBJCHIS = {}

    # get unmodified VDW strain
    for p in objects:
        if (hydrogens):
            cmd.h_add(p)
        else:
            cmd.remove('((%s) and elem H)' %p)

        OBJCHIS[p] = ([[get_strain('(%s or %s) and not (%s)'%(p,optimize_add,optimize_ignore),
        temp_clash_A,optimize),0,0],[]])

    ##### MDA BUILDER #####

    ##### PRE_SETTINGS #####
    # cycles through residues (first lysines then N-terms)
    MDA_ID = []
    listcount=-1
    for p in list(stored.resi_list):
        cmd.unpick()
        listcount=listcount+1
        #jump skipped
        if instance_list[listcount]==0:continue

        # set random G_ID, if selected
        if group==1:
            G_ID=random.choice([2,3,4])
        if confomer==1:
            C_ID=random.choice([2,3,4,5])
        if type==1:
            T_ID=random.choice([2,3,4,5,6,7])
        MDA_ID.append([G_ID,C_ID,T_ID])

        # get info on residue being processed
        stored.residue=''
        cmd.iterate('%s and name CA' %p, 'stored.residue=str(resn)')

        #pick starting Nitrogen
        if (instance_list[listcount] in [3,4]):
            # handled as lysine
            cmd.edit('%s and name NZ' %p)
        else:
            # handled as n-term
            cmd.edit('%s and (name N)' %p)

        # get valence (iterate is not supported for valence in older PyMOL versions)
        #make sure to fill hydrogens
        cmd.h_fill()
        val=cmd.count_atoms('neighbor pk1')

        cmd.remove('%s and (neighbor pk1) and hydrogens' %p)

        # adjust (assumed) geometry and valence of nitrogen
        geom=4

        if ((G_ID==3) and (T_ID>3)):
            # double bond at nitrogen, assumes sp2
            geom=3
        if (protonate==-1):
            val=3
            cmd.alter('pk1', 'formal_charge=0')
        if (protonate==1):
            val=4
            cmd.alter('pk1', 'formal_charge=1')
        cmd.set_geometry('pk1', geom, val)


        ##### MDD/EDD TYPE ADDUCTS START #####
        # selected = N
        if G_ID==2:
            # attach C-ring
            for c in ['CM2','CM3','CM4','CM5','CM6']:
                if c=='CM4':
                    cmd.attach('C','4','4')
                else:
                    cmd.attach('C','3','3')
                cmd.alter('(neighbor pk1) and (elem C) '+
                'and not ((name CE) or (name CA) or (name CM*))',
                'name="%s"' %c)

                #double bonds
                if c in ['CM3','CM6']:
                    cmd.unbond('pk1',
                    '%s and (pk1 around 2) and (name %s)' %(p, c))
                    cmd.bond('pk1',
                    '%s and (pk1 around 2) and (name %s)' %(p, c),'2')

                # fix angle
                x=angles[['CM2','CM3','CM4','CM5','CM6'].index(c)]
                set_angle(angle=x,
                atom1='(neighbor pk1) and not (name %s)' %c,
                atom2='pk1',
                atom3='(neighbor pk1) and (name %s)' %c)

                # pick added
                cmd.edit('(neighbor pk1) and (name %s)' %c)

            # Done adding ring C atoms

            # carbaldehydes
            select_part(temp_sel,'%s'%p,'pk1','(name CB)')
            for c in [3,5]:
                cmd.edit('%s and (name CM%s)' %(temp_sel,str(c)))
                cmd.attach('C','3','3')
                cmd.alter('(neighbor pk1) and not (name CM*)',
                'name="CMC%s"' %(str(c)))
                cmd.edit('%s and (pk1 extend 2) and (name CMC%s)' %(p,str(c)))
                cmd.attach('O','2','2')
                cmd.alter('(neighbor pk1) and (elem O) ',
                'name="OMC%s"' %(str(c)))
                # double bond
                cmd.unbond('pk1',
                '%s and (pk1 around 2) and (name OMC%s)' %(p,str(c)))
                cmd.bond('pk1',
                '%s and (pk1 around 2) and (name OMC%s)' %(p,str(c)),'2')
                set_angle(angle=120,
                atom1='(neighbor pk1) and not (name OMC%s)'%c,
                atom2='pk1',
                atom3='(neighbor pk1) and (name OMC%s)'%c)

            # adjust torsion of ring before closing
            #re-pick N
            cmd.edit('(%s) and (name N or name NZ)' %temp_sel)

            if (C_ID in [2,3]):
                #cis boat
                x=[['(name CD or name C)','(name CE or name CA)',
                    '(name N or name NZ)','name CM2',torsions[0]],
                   ['(name CE or name CA)','(name N or name NZ)',
                    'name CM2','name CM3',torsions[1]-90],
                   ['(name N or name NZ)','name CM2',
                    'name CM3','name CM4',torsions[2]],
                   ['name CM2','name CM3','name CM4',
                    'name CM5',torsions[4]],
                   ['name CM3','name CM4','name CM5',
                    'name CM6',torsions[3]]]
            else:
                #trans boat
                x=[['(name CD or name C)','(name CE or name CA)',
                    '(name N or name NZ)','name CM2',torsions[0]],
                   ['(name CE or name CA)','(name N or name NZ)',
                    'name CM2','name CM3',torsions[1]],
                   ['(name N or name NZ)','name CM2',
                    'name CM3','name CM4',torsions[2]],
                   ['name CM2','name CM3','name CM4','name CM5',torsions[3]],
                   ['name CM3','name CM4','name CM5','name CM6',torsions[4]]]


            for c in x:
                cmd.set_dihedral(
                '((%s) and (%s))'  %(temp_sel,c[0]),
                '((%s) and (%s))'  %(temp_sel,c[1]),
                '((%s) and (%s))'  %(temp_sel,c[2]),
                '((%s) and (%s))' %(temp_sel,c[3]),
                c[4]
                )

            #close ring
            cmd.bond('%s and name CM6' %temp_sel,
            '%s and elem N' %temp_sel, '1')


            ##### MDD or EDD position 4 #####
            cmd.edit('%s and (pk1 extend 3) and name CM4' %p)

            cmd.attach('H','1','1')
            ##### MDD #####
            if (T_ID==3):
                cmd.attach('C','4','4')
                cmd.alter('(neighbor pk1) and (elem C) '+
                'and not (name CM*)',
                'name="CMM1"')
                # adjust torsion to avoid sterical clashing of H atoms
                cmd.edit('(neighbor pk1) and (elem C) '+
                'and (name CMM1)')
                cmd.attach('H','1','1')
                select_part(temp_sel,'%s'%p,'pk1','(name CB)')
                cmd.attach('H','1','1')
                cmd.attach('H','1','1')
                cmd.set_dihedral(
                '%s and ((neighbor (name CM4)) and (elem H))' %temp_sel,
                '%s and (name CM4)' %temp_sel,
                '%s and (name CMM1)' %temp_sel,
                '%s and ((neighbor (name CMM1)) and (elem H))' %temp_sel,
                180)

            ##### EDD #####
            if (T_ID in [4,5,6,7]):
                cmd.attach('C','3','3')
                cmd.alter('(neighbor pk1) and (elem C) '+
                'and not (name CM3 or name CM5)',
                'name="CME1"')
                cmd.edit('%s and (pk1 extend 2) and name CME1' %p)
                cmd.attach('C','3','3')
                cmd.alter('(neighbor pk1) and (elem C) '+
                'and not (name CM4)',
                'name="CME2"')
                # double bonds
                if (T_ID in [4,5]):
                    #Ethenol
                    cmd.unbond('pk1',
                    '%s and (pk1 around 2) and (name CME2)' %p)
                    cmd.bond('pk1',
                    '%s and (pk1 around 2) and (name CME2)' %p,'2')

                cmd.edit('%s and (pk1 extend 2) and name CME2' %p)
                cmd.attach('O','2','2')
                cmd.alter('(neighbor pk1) and (elem O) ',
                'name="OME1"')
                set_angle(angle=120,
                atom1='(neighbor pk1) and not (name OME1)',
                atom2='pk1',
                atom3='(neighbor pk1) and (name OME1)')

                if (T_ID in [6,7]):
                    # Ethanal
                    cmd.unbond('pk1',
                    '%s and (pk1 around 2) and (name OME1)' %p)
                    cmd.bond('pk1',
                    '%s and (pk1 around 2) and (name OME1)' %p,'2')

            ## final hydration
            # fill hydrogens
            MDA= (
                  ['CM2','CM3','CM4','CM5','CM6','CMC3','CMC5']+
                  ['OMC3','OMC5']
                 )
            # CMM1 is omitted to avoid resetting of rotation!

            if (T_ID in [4,5,6,7]):
                MDA= (
                  ['CM2','CM3','CM4','CM5','CM6','CMC3','CMC5']+
                  ['OMC3','OMC5','CME1','CME2','OME1']
                     )
            for c in MDA:
                try:
                    cmd.edit('%s and (name %s) and (pk1 extend 7)' %(p, c))
                    cmd.h_fill()
                except:
                    continue

            # proximal/distal oxgen rotamers
            select_part(temp_sel,'%s'%p,'pk1','(name CB)')

            cmd.set_dihedral(
            '%s and (name CM4)' %temp_sel,
            '%s and (name CM3)' %temp_sel,
            '%s and (name CMC3)' %temp_sel,
            '%s and (name OMC3)' %temp_sel,
            torsions[5])

            cmd.set_dihedral(
            '%s and (name CM4)' %temp_sel,
            '%s and (name CM5)' %temp_sel,
            '%s and (name CMC5)' %temp_sel,
            '%s and (name OMC5)' %temp_sel,
            torsions[6])

            # alpha_beta
            if (T_ID>2):
                ang_0=cmd.get_angle(
                atom1=('(%s) and (neighbor (name CM4)) '
                       'and (elem H)' %(temp_sel)),
                atom2='(%s) and (name CM4)' %(temp_sel),
                atom3='(%s) and ((name CME1) or (name CMM1))' %(temp_sel))
                # move to overlapping pos and get axis
                ang_1=set_angle(
                angle=0,
                atom1=('(%s) and (neighbor (name CM4)) '
                       'and (elem H)' %(temp_sel)),
                atom2='(%s) and (name CM4)' %(temp_sel),
                atom3='(%s) and ((name CME1) or (name CMM1))' %(temp_sel))
                # note: in trans boats alpha and beta get reversed!
                # thus [2,5] and not [2,4]

                if (C_ID in [2,5]):
                    # alpha
                    set_angle(angle=ang_0,
                    atom1=('(%s) and (neighbor (name CM4)) '
                           'and (elem H)' %(temp_sel)),
                    atom2='(%s) and (name CM4)' %(temp_sel),
                    atom3=('(%s) and ((name CME1) or '
                           '(name CMM1))' %(temp_sel)),
                    axis=ang_1)
                else:
                    # beta
                    set_angle(angle=ang_0,
                    atom1=('(%s) and ((name CME1) or '
                           '(name CMM1))' %(temp_sel)),
                    atom2='(%s) and (name CM4)' %(temp_sel),
                    atom3=('(%s) and (neighbor (name CM4)) '
                           'and (elem H)' %(temp_sel)),
                    axis=ang_1)

            # EDD adduct
            if (T_ID in [4,5,6,7]):
                # oxygen cis/trans
                t_angle=0
                if (T_ID in [5,7]):
                    t_angle=180

                cmd.set_dihedral(
                '%s and (name CM4)' %temp_sel,
                '%s and (name CME1)' %temp_sel,
                '%s and (name CME2)' %temp_sel,
                '%s and (name OME1)' %temp_sel,
                t_angle)

                if (T_ID in [4,5]):
                    cmd.set_dihedral(
                    ('%s and ((neighbor (name CME2)) '
                     'and (elem H))' %temp_sel),
                    '%s and (name CME2)' %temp_sel,
                    '%s and (name OME1)' %temp_sel,
                    ('%s and ((neighbor (name OME1)) '
                    'and (elem H))' %temp_sel),
                    0)

                # adduct orientation
                t_angle=torsions[7]
                cmd.set_dihedral(
                ('%s and ((neighbor (name CM4)) '
                'and (elem H))' %temp_sel),
                '%s and (name CM4)' %temp_sel,
                '%s and (name CME1)' %temp_sel,
                '%s and (name CME2)' %temp_sel,
                t_angle)
            # Done rotating EDD

        ##### MDD/EDD TYPE ADDUCTS END #####

        ##### MDA TYPE ADDUCTS START #####
        if G_ID==3:
            ## add C atoms using loop
            for c in range(1,4):
                # attach planar C
                # here we assume a average sp2-like planar
                # state for all involved atoms
                cmd.attach('C','3','3')
                cmd.alter('(neighbor pk1) and (elem C) '+
                'and not (name CE or name CM1 or name CA or name CD)',
                'name="CM%s"' %(str(c)))
                if (not ((c==1) and (T_ID in [2,3]) and
                    (protonate!=-1) and (stored.residue=='LYS'))):
                    # it is not a protonated lysine with single bond at N-C
                    set_angle(angle=120,
                    atom1='(neighbor pk1) and not (name CM%s)' %(str(c)),
                    atom2='pk1',
                    atom3='(neighbor pk1) and (name CM%s)' %(str(c)))
                # double bonds
                if (
                    ((c==2) and (T_ID in [2,3])) or
                    ((c!=2) and (T_ID in [4,5])) or
                    ((c==1) and (T_ID in [6,7]))
                   ):
                    cmd.unbond('pk1',
                    '%s and (pk1 around 2) and (name CM%s)' %(p,str(c)))
                    cmd.bond('pk1',
                    '%s and (pk1 around 2) and (name CM%s)' %(p,str(c)),'2')
                # pick new to edit
                cmd.edit('%s and (name CM%s) and (pk1 around 2)' %(p,str(c)))

            # add oxygen
            cmd.attach('O','2','2')
            cmd.alter('(neighbor pk1) and not (name CM2)','name="OM1"')
            set_angle(angle=120,
            atom1='(neighbor pk1) and (name CM2)',
            atom2='pk1',
            atom3='(neighbor pk1) and (name OM1)')
            # double bonded for "-al types"
            if (not (T_ID in [4,5])):
                cmd.unbond('pk1', '%s and (pk1 around 2) and (name OM1)' %p)
                cmd.bond('pk1', '%s and (pk1 around 2) and (name OM1)' %p,'2')

            ## final hydration of c-atoms
            select_part(temp_sel,'%s'%p,'pk1','(name CB)')
            MDA=['CM1','CM2','CM3','OM1']
            for c in MDA:
                cmd.edit('%s and (name %s)' %(temp_sel, c))
                cmd.h_fill()

            # rotamers/isomers
            select_part(temp_sel,'%s'%p,'pk1','(name CB)')
            # adduct cis(E)/trans(Z) decision - adduct to R1-N
            t_angle=0
            if (C_ID>3):
                t_angle=180
            cmd.set_dihedral(
            '%s and ((name CA) or (name CE))' %temp_sel,
            '%s and ((name N) or (name NZ))' %temp_sel,
            '%s and (name CM1)' %temp_sel,
            '%s and (name CM2)' %temp_sel,
            t_angle)

            # adduct cis/trans decision - adduct to double bond
            t_angle=0
            if (C_ID in [3,5]):
                t_angle=180
            cmd.set_dihedral(
            '%s and ((name N) or (name NZ))' %temp_sel,
            '%s and (name CM1)' %temp_sel,
            '%s and (name CM2)' %temp_sel,
            '%s and (name CM3)' %temp_sel,
            t_angle)

            # adduct cis/trans decision - Oxygen
            t_angle=0
            if (T_ID in [3,5,7]):
                t_angle=180
            cmd.set_dihedral(
            '%s and (name CM1)' %temp_sel,
            '%s and (name CM2)' %temp_sel,
            '%s and (name CM3)' %temp_sel,
            '%s and (name OM1)' %temp_sel,
            t_angle)

            # force OH to trans (sterically favorable)
            if (T_ID in [4,5]):
                t_angle=180
                cmd.set_dihedral(
                '%s and (name CM2)' %temp_sel,
                '%s and (name CM3)' %temp_sel,
                '%s and (name OM1)' %temp_sel,
                '%s and ((neighbor (name OM1)) and (elem H))' %temp_sel,
                t_angle)

            #pymol may incorrectly position H-atom at CM1 using cmd.h_fill()
            cmd.edit('(%s) and ((neighbor (name CM1)) and (elem H))' %temp_sel)
            cmd.remove('pk1')
            cmd.edit('(%s) and (name CM1)' %temp_sel)
            cmd.attach('H','1','1')
            # N is hydrated later
        ##### MDA TYPE ADDUCTS END #####


        ##### FAAB TYPE ADDUCTS START #####
        if G_ID==4:
            cmd.attach('C','4','4')
            cmd.alter('(neighbor pk1) and (elem C) '+
            'and not (name CE or name CA)',
            'name="CM2"')
            cmd.edit('(neighbor pk1) and (name CM2)')

            cmd.attach('C','4','4')
            cmd.alter('(neighbor pk1) and (elem C) ',
            'name="CMA"')
            cmd.attach('C','3','3')
            cmd.alter('(neighbor pk1) and (elem C) and not (name CMA)',
            'name="CM3"')
            cmd.attach('C','3','3')
            cmd.alter('(neighbor pk1) and (elem C) and not (name CMA or name CM3)',
            'name="CMB"')

            # R or S?
            # CM2 now has 4 neighbors, NZ and CM3 being fixed, one of the
            # remaining two will be removed to yield R or S
            # picture an orientation looking at CM2 with the nitrogen facing up and
            # CM3 to the front
            # in R, the methyl group will be on the left (smaller dihedral)
            # in S, the methyl group will be on the right (larger dihedral)

            select_part(temp_sel,'%s'%p,'pk1','(name CB)')

            # zero plane
            dihedCM3=to360(cmd.get_dihedral(
            '%s and ((name CE) or (name CA))' %temp_sel,
            '%s and ((name NZ) or (name N))' %temp_sel,
            '%s and (name CM2)' %temp_sel,
            '%s and (name CM3)' %temp_sel))

            # first oxygen
            dihedCMA=to360(cmd.get_dihedral(
            '%s and ((name CE) or (name CA))' %temp_sel,
            '%s and ((name NZ) or (name N))' %temp_sel,
            '%s and (name CM2)' %temp_sel,
            '%s and (name CMA)' %temp_sel))
            dihedCMA=to360(dihedCMA-dihedCM3)

            # other oxygen
            dihedCMB=to360(cmd.get_dihedral(
            '%s and ((name CE) or (name CA))' %temp_sel,
            '%s and ((name NZ) or (name N))' %temp_sel,
            '%s and (name CM2)' %temp_sel,
            '%s and (name CMB)' %temp_sel))
            dihedCMB=to360(dihedCMB-dihedCM3)

            # default to delete: CMB assuming R
            if C_ID in [2,4]:
                # R
                dihedCM3='CMB'
            else:
                # S
                dihedCM3='CMA'
            # test angles
            if dihedCMB<dihedCMA:
                # invert default
                if C_ID in [2,4]:
                    # R
                    dihedCM3='CMA'
                else:
                    # S
                    dihedCM3='CMB'
            # delete extra
            cmd.remove('%s and (name %s)' %(temp_sel,dihedCM3))
            # rename to CM1
            cmd.alter('%s and (name CMA or name CMB)' %temp_sel,
            'name="CM1"')
            # R/S is now set!

            # fake/temp sp2 for CM3; pk1=CM2
            select_part(temp_sel,'%s'%p,'pk1','(name CB)')
            cmd.unbond('(%s) and (name CM2)' %temp_sel,
            '(%s) and (name CM3)' %temp_sel)
            cmd.bond('(%s) and (name CM2)' %temp_sel,
            '(%s) and (name CM3)' %temp_sel,2)

            for c in [['CM4A','OMA'],['CM4B','OMB']]:
                cmd.edit('(neighbor pk1) and name CM3')
                cmd.attach('C','3','3')
                cmd.alter('(neighbor pk1) and not (name CM*)',
                'name="%s"' %(c[0]))
                set_angle(angle=120,
                atom1='(neighbor pk1) and not (name CM4*)',
                atom2='pk1',
                atom3='(neighbor pk1) and (name %s)' %(c[0]))
                #attach Oxygen
                cmd.edit('(neighbor pk1) and (name %s)' %(c[0]))
                cmd.attach('O','2','2')
                cmd.alter('(neighbor pk1) and (elem O)',
                'name="%s"' %(c[1]))
                set_angle(angle=120,
                atom1='(neighbor pk1) and not (name %s)' %(c[1]),
                atom2='pk1',
                atom3='(neighbor pk1) and (name %s)' %(c[1]))

            #Hydrogen HF
            cmd.edit('(pk1 extend 3) and (name OMA)')
            cmd.attach('H','1','1')
            cmd.alter('(neighbor pk1) and (elem H)',
            'name="HF"')
            set_angle(angle=120,
            atom1='(neighbor pk1) and not (name HF)',
            atom2='pk1',
            atom3='(neighbor pk1) and (name HF)')

            #first delocalized double bonds
            bonds=2
            if T_ID in [2,3,5,6]:bonds=4

            select_part(temp_sel,'%s'%p,'pk1','(name CB)')
            x=[['CM3','CM4A'],['CM4B','OMB']]
            if (not (T_ID in [4,7])): x.append(['OMA','HF'])
            for c in x:
                cmd.unbond('(%s) and (name %s)' %(temp_sel,c[0]),
                '(%s) and (name %s)' %(temp_sel,c[1]))
                cmd.bond('(%s) and (name %s)' %(temp_sel,c[0]),
                '(%s) and (name %s)' %(temp_sel,c[1]),bonds)

            # undo fake DB at CM3
            cmd.unbond('(%s) and (name CM2)' %temp_sel,
            '(%s) and (name CM3)' %temp_sel)
            cmd.bond('(%s) and (name CM2)' %temp_sel,
            '(%s) and (name CM3)' %temp_sel,1)

            # final hydrations in chain
            MDA=['CM1','CM2','CM3','CM4A','CM4B']
            # omitted 'OMA','OMB'
            for c in MDA:
                cmd.edit('(%s) and (name %s)' %(temp_sel, c))
                cmd.h_fill()

            # fix torsions
            x=[['(name CE) or (name CA)','(name NZ) or (name N)',
                'name CM2','name CM3',torsions[8]],
               ['(name NZ) or (name N)','name CM2',
                'name CM3','name CM4A',torsions[9]],
               ['name CM2','name CM3','name CM4A','name OMA',180.0],
               ['name CM2','name CM3','name CM4B','name OMB',180.0],
               ['name CM3','name CM4A','name OMA','name HF',0.0]]

            if C_ID in [3,5]:
                x[0][4]=-1*torsions[8]
                x[1][4]=-1*torsions[9]

            for c in x:
                cmd.set_dihedral(
                '%s and (%s)' %(temp_sel,c[0]),
                '%s and (%s)' %(temp_sel,c[1]),
                '%s and (%s)' %(temp_sel,c[2]),
                '%s and (%s)' %(temp_sel,c[3]),
                c[4])

            # final delocalized bonds
            x=[['CM3','CM4B'],['CM4A','OMA']]
            if T_ID in [2,5]: x.append(['OMB','HF'])
            if (not (T_ID in [4,7])):
                bonds=4
                for c in x:
                    cmd.unbond('(%s) and (name %s)' %(temp_sel,c[0]),
                    '(%s) and (name %s)' %(temp_sel,c[1]))
                    cmd.bond('(%s) and (name %s)' %(temp_sel,c[0]),
                    '(%s) and (name %s)' %(temp_sel,c[1]),bonds)

        ##### FAAB TYPE ADDUCTS END #####


        ##### final hydration and protonations (all adducts) #####
        if ((stored.residue=='LYS') or (stored.residue=='LYN')):
            cmd.edit('%s and name NZ' %(p))
        else:
            cmd.edit('%s and name N' %(p))

        # complete hydration
        cmd.h_fill()
        ## change name if LYS
        if instance_list[listcount] in  [3,4]:
            cmd.alter(p,'resn="MMK"')
        else:
            # adjust resn, resi and resv
            stored.temp_n=''
            cmd.iterate('(%s) and name N' %(temp_sel),
            'stored.temp_n=int(resv)')
            stored.temp_n=stored.temp_n-1
            cmd.alter(('(%s and (((name CM*) or '
            '(name OM*)) extend 1)) and not (name N)' %temp_sel),
            'resn="MDA"')
            cmd.alter(('(%s and (((name CM*) or '
            '(name OM*)) extend 1)) and not (name N)' %temp_sel),
            'resv=%s' %stored.temp_n)
            cmd.alter(('(%s and (((name CM*) or '
            '(name OM*)) extend 1)) and not (name N)' %temp_sel),
            'resi=%s' %stored.temp_n)
            # replace resi_list entry
            stored.resi_list[listcount]=str(
            "((model "+str(cmd.get_object_list(temp_sel)[0])+
            " and chain '"+str((cmd.get_chains(temp_sel) or [''])[0])+
            "' and resi "+str(stored.temp_n)+") extend 3) "+
            "and not ((name CB) or ((neighbor (name CB)) and (elem H)))"
            )

        #Log
        last_time=log_pytms_prog(output,
        'Modified: %s'%get_resi_macro_name(p),listcount+1,len(stored.resi_list))

    cmd.unpick()
    # end of p in names
    #FINISHED BUILDING MDA ADDUCTS #####

    # fixes ordering of resis
    for p in objects:
        cmd.remove('((%s) and elem H)' %p)
        cmd.h_add(p)

    # re-adjust hydrogens
    if (hydrogens):
        for p in objects: cmd.h_add(p)
    else:
        for p in objects:
            cmd.remove('((%s) and elem H)' %p)

    # coloring
    if color_base:
        cmd.color(color_base,
        '%s or ((%s extend 1) and hydrogens)' %(selection, selection))
    for p in stored.resi_list:
        #enable selection "p.PTM in MDA"
        if version_ok:
            cmd.alter(('%s and (((name CM*) '+
            'or (name OM*)) extend 1) '+
            'and not (elem N)') %p,
            'p.PTM="MDA"')
        if color_base:
            cmd.color(color_base,
            '(%s)' %(p))
        if color_mod:
            cmd.color(color_mod,
            ('(%s and (((name CM*) or (name OM*)) '+
            'extend 1) and not (elem N))') %(p))

    #######
    #     #        Congratulations!
    # :-) # You've found the hidden smiley!
    #     #     Save it for a bad day!
    #######
    last_time=log_pytms_prog(output, 'MDA modifications complete!')

    ##### VDW OPTIMIZATION #####

    # get modified strain
    for p in objects:
        OBJCHIS[p][0][1]=OBJCHIS[p][0][2]= get_strain(
        '(%s or %s) and not (%s)'%(p,optimize_add,optimize_ignore),
        temp_clash_A,optimize)


    if optimize >0:
        last_time=log_pytms_prog(output, 'Starting optimization (be patient ;-))!')

        # Chi angles of LYSINE (dictionary)
        CHIS = [
                ['name N','name CA','name CB','name CG'],
                ['name CA','name CB','name CG','name CD'],
                ['name CB','name CG','name CD','name CE'],
                ['name CG','name CD','name CE','name NZ'],
                ['name CD','name CE','name NZ',
                 '(neighbor (name NZ)) and (name CM*) '+
                 'and not (name CM6)'],
                ['name C','name CA','name N',
                 '(neighbor (name N)) and (name CM*) '+
                 'and not (name CM6)']
               ]

        ##### START OF P IN OBJECTS

        # progress will increment by residue
        progress=0.0
        for p in objects:
            listcount=-1
            for q in list(stored.resi_list):
                listcount=listcount+1
                # pass if residue not in object or unknown
                if (instance_list[listcount]==0): continue
                if (cmd.count_atoms('%s and %s'%(p,q))==0): continue

                #reset VALS with starting values
                VALS=([q, instance_list[listcount],
                      [[0.0,0.0,0.0,0.0,0.0],
                      [0.0,0.0,0.0,0.0,0.0],
                      [0.0,0.0,0.0,0.0,0.0]],
                             MDA_ID[listcount],[0,0]])


                if (instance_list[listcount]<3):
                    # not treated as lysine
                    VALS[2][0]=VALS[2][1]=(
                    cmd.get_dihedral(
                    '%s and %s' %(q, CHIS[5][0]),
                    '%s and %s' %(q, CHIS[5][1]),
                    '%s and %s' %(q, CHIS[5][2]),
                    '%s and %s' %(q, CHIS[5][3])
                    ))
                else:
                    # lysine
                    for  c in range(0,5):
                        VALS[2][0][c]=VALS[2][1][c]=(
                        cmd.get_dihedral(
                        '%s and %s' %(q, CHIS[c][0]),
                        '%s and %s' %(q, CHIS[c][1]),
                        '%s and %s' %(q, CHIS[c][2]),
                        '%s and %s' %(q, CHIS[c][3])
                        ))
                OBJCHIS[p][1].append(list(VALS))
            # done writing starting values

            # optimizing
            for q in OBJCHIS[p][1]:
                # is not skipped resi
                if q[1]==0:
                    continue

                x=OBJCHIS[p][1].index(q)
                optstrain=0
                # get info
                stored.residue=''
                cmd.iterate('%s and name CA' %q[0], 'stored.residue=str(resi)')
                # define local and get starting strain
                # base selection
                optimize_local= ('((byres (%s expand %f)) and (%s or %s))'
                %(q[0], local_radius, p,optimize_add))
                # without ignored
                optimize_local= ('(%s and not (%s))'
                %(optimize_local, optimize_ignore))
                # without adducts optimized later
                optimize_local= ('(%s and not (((((name CM*) or (name OM*)) extend 1)'
                ' and not (elem N)) and %s and ((resi %s-) and not %s)))'
                %(optimize_local, p, stored.residue, q[0]))
                # make selection
                cmd.select(temp_sel,optimize_local)

                optstrain=get_strain(optimize_local,temp_clash_A,optimize)
                OBJCHIS[p][1][x][4][0]=OBJCHIS[p][1][x][4][1]=optstrain


                if q[1]<3:
                    # not treated as lysine

                    first=ang_4=q[2][0]
                    last=first

                    while (ang_4 < ((q[2][0])+360)):
                        # set new angles
                        cmd.set_dihedral(
                        '%s and %s' %(q[0], CHIS[5][0]),
                        '%s and %s' %(q[0], CHIS[5][1]),
                        '%s and %s' %(q[0], CHIS[5][2]),
                        '%s and %s' %(q[0], CHIS[5][3]),
                        ang_4
                        )

                        # get new vdw strain
                        strain=get_strain(temp_sel,
                        temp_clash_A,optimize)

                        # check if better
                        if strain<optstrain:
                            optstrain=strain
                            first = ang_4
                        if strain==optstrain:
                            last=ang_4
                        if optstrain==0:
                            #absolute min. (exit calc.)
                            ang_4=ang_4+361.0
                            print ('!!! PyTMs Warning: '
                                   'local zero strain '
                                   'occured for: %s'%q[0])
                        ang_4=ang_4+intervals[0]
                    # END OF N-Term optimization

                    # set to optimal and write optimal to dictionary
                    if interpolate:
                        ang_4=to180((to360(first)+to360(last))/2)
                    else:
                        if (
                        (abs(to180(to360(first)-to360(q[2][0])))<=
                         abs(to180(to360(last)-to360(q[2][0]))))
                        ):
                            ang_4=to180(first)
                        else:
                            ang_4=to180(last)

                    cmd.set_dihedral(
                    '%s and %s' %(q[0], CHIS[5][0]),
                    '%s and %s' %(q[0], CHIS[5][1]),
                    '%s and %s' %(q[0], CHIS[5][2]),
                    '%s and %s' %(q[0], CHIS[5][3]),
                    ang_4
                    )

                    if interpolate:
                        # check that average is good
                        strain=get_strain(temp_sel,temp_clash_A,optimize)
                        if strain>optstrain:
                            if (
                            (abs(to180(to360(first)-to360(q[2][0])))<=
                             abs(to180(to360(last)-to360(q[2][0]))))
                            ):
                                ang_4=to180(first)
                            else:
                                ang_4=to180(last)
                            cmd.set_dihedral(
                            '%s and %s' %(q[0], CHIS[5][0]),
                            '%s and %s' %(q[0], CHIS[5][1]),
                            '%s and %s' %(q[0], CHIS[5][2]),
                            '%s and %s' %(q[0], CHIS[5][3]),
                            ang_4
                            )

                    # write to dictionary
                    OBJCHIS[p][1][x][2][1]=ang_4
                    OBJCHIS[p][1][x][2][2]=(
                    to180(to360(ang_4)-to360(OBJCHIS[p][1][x][2][0]))
                    )
                    # optimal local strain
                    OBJCHIS[p][1][x][4][1]=optstrain

                else:
                    # Lysines
                    # optimal base strain = base strain - lysine only
                    # lysine only:
                    optimize_base=('(%s) and not ((%s) and '
                    '((((name CM*) or (name OM*)) extend 1) '
                    'and not (elem N)))'%(q[0],q[0]))
                    optbasestrain=get_strain(optimize_base,temp_clash_A)

                    # base selection = local-resi.MDA
                    optimize_base=('(%s) and not ((%s) and '
                    '((((name CM*) or (name OM*)) extend 1) '
                    'and not (elem N)))'%(temp_sel,q[0]))

                    # optimal base strain = base strain - lysine only
                    optbasestrain=get_strain(optimize_base,temp_clash_A)-optbasestrain
                    # allow for limit
                    optbasestrain=optbasestrain+base_strain_limit

                    first=[q[2][0][0],
                           q[2][0][1],
                           q[2][0][2],
                           q[2][0][3],
                           q[2][0][4]]
                    last=list(first)
                    # cycle through chi angles
                    ang_0=q[2][0][0]
                    while (ang_0 < ((q[2][0][0])+360)):
                        cmd.set_dihedral(
                        '%s and %s' %(q[0], CHIS[0][0]),
                        '%s and %s' %(q[0], CHIS[0][1]),
                        '%s and %s' %(q[0], CHIS[0][2]),
                        '%s and %s' %(q[0], CHIS[0][3]),
                        ang_0
                        )
                        ang_1=q[2][0][1]
                        while (ang_1 < ((q[2][0][1])+360)):
                            cmd.set_dihedral(
                            '%s and %s' %(q[0], CHIS[1][0]),
                            '%s and %s' %(q[0], CHIS[1][1]),
                            '%s and %s' %(q[0], CHIS[1][2]),
                            '%s and %s' %(q[0], CHIS[1][3]),
                            ang_1
                            )

                            ang_2=q[2][0][2]
                            while (ang_2 < ((q[2][0][2])+360)):
                                cmd.set_dihedral(
                                '%s and %s' %(q[0], CHIS[2][0]),
                                '%s and %s' %(q[0], CHIS[2][1]),
                                '%s and %s' %(q[0], CHIS[2][2]),
                                '%s and %s' %(q[0], CHIS[2][3]),
                                ang_2
                                )

                                ang_3=q[2][0][3]
                                while (ang_3 < ((q[2][0][3])+360)):
                                    cmd.set_dihedral(
                                    '%s and %s' %(q[0], CHIS[3][0]),
                                    '%s and %s' %(q[0], CHIS[3][1]),
                                    '%s and %s' %(q[0], CHIS[3][2]),
                                    '%s and %s' %(q[0], CHIS[3][3]),
                                    ang_3
                                    )

                                    ang_4=q[2][0][4]
                                    base_tested=False
                                    if base_strain_limit>1500: base_tested=True
                                    if ((ang_0==q[2][0][0]) and
                                        (ang_1==q[2][0][1]) and
                                        (ang_2==q[2][0][2]) and
                                        (ang_3==q[2][0][3])):
                                        # don't repeat start strain
                                        ang_4=ang_4+intervals[0]

                                    while (ang_4 < ((q[2][0][4])+360)):
                                        # set new angles
                                        cmd.set_dihedral(
                                        '%s and %s' %(q[0], CHIS[4][0]),
                                        '%s and %s' %(q[0], CHIS[4][1]),
                                        '%s and %s' %(q[0], CHIS[4][2]),
                                        '%s and %s' %(q[0], CHIS[4][3]),
                                        ang_4
                                        )

                                        # unfavorable base strain?
                                        # from lysine and surrounding
                                        if (not base_tested):
                                            base_tested=True
                                            basestrain=get_strain(optimize_base,
                                            temp_clash_A,optimize)
                                            if basestrain>(optbasestrain):
                                                #jump over loop
                                                break

                                        # get new vdw strain
                                        strain=get_strain(temp_sel,
                                        temp_clash_A,optimize)

                                        # check if better
                                        if strain<optstrain:
                                            optstrain=strain
                                            first=(
                                            [ang_0,ang_1,ang_2,ang_3,ang_4])
                                        if (strain==optstrain):
                                            last=(
                                            [ang_0,ang_1,ang_2,ang_3,ang_4])
                                        if optstrain==0:
                                            #absolute min. (exit calc.)
                                            ang_4=ang_4+361.0
                                            ang_3=ang_3+361.0
                                            ang_2=ang_2+361.0
                                            ang_1=ang_1+361.0
                                            ang_0=ang_0+361.0
                                            print ('!!! PyTMs Warning: '
                                                   'local zero strain '
                                                   'occured for: %s'%q[0])
                                        # run this max. 1 full turn
                                        # and only change angles
                                        #specified by 'optimize'
                                        ang_4=ang_4+intervals[0]
                                    ang_3=ang_3+intervals[1]
                                ang_2=ang_2+intervals[2]
                            ang_1=ang_1+intervals[3]
                        ang_0=ang_0+intervals[4]
                    # END OF lys optimization

                    # set to optimal and write optimal to dictionary
                    #lysine
                    ang_4=[0,0,0,0,0]
                    if interpolate:
                        for c in range(0,5):
                            ang_4[c]=to180((to360(first[c])+to360(last[c]))/2)
                    else:
                        sum_first=0
                        sum_last=0
                        for c in range(0,5):
                            sum_first=sum_first+(c+1)*abs(to180(to360(first[c])-to360(q[2][0][c])))
                            sum_last=sum_last+(c+1)*abs(to180(to360(last[c])-to360(q[2][0][c])))
                            #smallest total (factored with c+1 to give more impact to proximal chi angles)
                        if (sum_first<=sum_last):
                            for c in range(0,5):
                                ang_4[c]=to180(first[c])
                        else:
                            for c in range(0,5):
                                ang_4[c]=to180(last[c])

                    for c in range(0,5):
                        cmd.set_dihedral(
                        '%s and %s' %(q[0], CHIS[c][0]),
                        '%s and %s' %(q[0], CHIS[c][1]),
                        '%s and %s' %(q[0], CHIS[c][2]),
                        '%s and %s' %(q[0], CHIS[c][3]),
                        (ang_4[c])
                        )

                    if interpolate:
                        # check that minimum is OK
                        strain=get_strain(temp_sel,temp_clash_A,optimize)
                        if strain>optstrain:
                            sum_first=0
                            sum_last=0
                            for c in range(0,5):
                                sum_first=sum_first+(c+1)*abs(to180(to360(first[c])-to360(q[2][0][c])))
                                sum_last=sum_last+(c+1)*abs(to180(to360(last[c])-to360(q[2][0][c])))
                                #smallest total (factored with c+1 to give more impact to proximal chi angles)
                            if (sum_first<=sum_last):
                                for c in range(0,5):
                                    ang_4[c]=to180(first[c])
                            else:
                                for c in range(0,5):
                                    ang_4[c]=to180(last[c])

                            for c in range(0,5):
                                ang_4[c]=to180(first[c])
                                cmd.set_dihedral(
                                '%s and %s' %(q[0], CHIS[c][0]),
                                '%s and %s' %(q[0], CHIS[c][1]),
                                '%s and %s' %(q[0], CHIS[c][2]),
                                '%s and %s' %(q[0], CHIS[c][3]),
                                ang_4[c])

                    # write to dictionary
                    for c in range(0,5):
                        OBJCHIS[p][1][x][2][1][c]=ang_4[c]

                        OBJCHIS[p][1][x][2][2][c]=(
                        to180((
                        (to360(OBJCHIS[p][1][x][2][1][c]))-
                        (to360(OBJCHIS[p][1][x][2][0][c])))
                        ))
                    # optimal local strain
                    OBJCHIS[p][1][x][4][1]=optstrain

                #Residue optimization complete

                # progress
                progress=progress+1
                log_pytms_prog(output,
                'Optimized: %s'%get_resi_macro_name(q[0]), progress, len(stored.resi_list),last_time)
            # END OF Q IN OBJCHIS[p][1]
            # write optimized strain
            OBJCHIS[p][0][2]= get_strain(
            '(%s or %s) and not (%s)'%(p,optimize_add,optimize_ignore),
            temp_clash_A,optimize)

        ##### END OF P IN OBJECTS

        last_time=log_pytms_prog(output,
        'Optimization complete!')
    ##### VDW OPTIMIZATION END #####

    # remove atoms near modified residues
    for p in objects:
        for q in OBJCHIS[p][1]:
            optimize_local=('(((byres (%s expand %f)) and (%s or %s)) '
            'and (%s))' %(q[0], remove_radius, p,optimize_add, optimize_ignore))
            cmd.remove(optimize_local)


    ##### ANIMATE VDW OPTIMIZATION START #####
    progress=0.0
    if ((states>0) and (optimize>0)):
        last_time=log_pytms_prog(output,
        'Initiating animation morphing!')
        for p in objects:
            cmd.delete(temp_obj)
            cmd.delete(temp_clash_A)

            # check for N-terminal MMKs
            # animation needs to be set to distinct due to PyMOl error
            # mixed atom positions every second state!
            # obsolete with N terminus being treated as separate residue
            discrete=0
            for q in OBJCHIS[p][1]:
                if (q[1] in [2,3]):
                    discrete=1
                    break

            for n in range(1,states+1):
                cmd.set('state',1,'%s'%p)
                frac=(float(n-1)/float(states-1))
                #float is due to intervals during animation:
                #int would be rounded
                for q in OBJCHIS[p][1]:
                    # skipped resis
                    if q[1]==0: continue

                    if q[1]<3:
                        # N-term
                        # recalculate angle
                        ang_4=to180(to360(q[2][0])+(q[2][2]*frac))

                        cmd.set_dihedral(
                        '%s and %s' %(q[0], CHIS[5][0]),
                        '%s and %s' %(q[0], CHIS[5][1]),
                        '%s and %s' %(q[0], CHIS[5][2]),
                        '%s and %s' %(q[0], CHIS[5][3]),
                        ang_4
                        )

                    else:
                        # LYS
                        # recalculate angle
                        ang_0=to180((q[2][0][0])+(q[2][2][0]*frac))
                        ang_1=to180((q[2][0][1])+(q[2][2][1]*frac))
                        ang_2=to180((q[2][0][2])+(q[2][2][2]*frac))
                        ang_3=to180((q[2][0][3])+(q[2][2][3]*frac))
                        ang_4=to180((q[2][0][4])+(q[2][2][4]*frac))

                        for c in range(0,5):
                            cmd.set_dihedral(
                            '%s and %s' %(q[0], CHIS[c][0]),
                            '%s and %s' %(q[0], CHIS[c][1]),
                            '%s and %s' %(q[0], CHIS[c][2]),
                            '%s and %s' %(q[0], CHIS[c][3]),
                            ([ang_0,ang_1,ang_2,ang_3,ang_4][c])
                            )
                    #Done with residues
                # update temporary models
                cmd.create(temp_obj,p,1,n,discrete=discrete)
            #Done with states

            # when done with all states --> rename object
            cmd.delete(p)
            cmd.create(p,temp_obj)

            # clashes
            if (show_clashes):
                optimize_local=('((%s or %s)'
                'and not %s)' %(p,optimize_add, optimize_ignore))
                cmd.create('%s_clashes'%p, optimize_local)
                # activate clash visibility
                for n in range(1,states+1):
                    cmd.set('state',n)
                    cmd.sculpt_activate('%s_clashes'%p)
                    cmd.show_as('cgo', '%s_clashes'%p)
                    cmd.set('sculpt_vdw_vis_mode', 1, '%s_clashes'%p)
                    cmd.sculpt_iterate('%s_clashes'%p, cycles=1)

            progress=progress+1
            log_pytms_prog(output,
            'Animating objects!', progress, len(objects), last_time)

    elif (show_clashes):
        # no animation but clashes
        for p in objects:
            optimize_local=('((%s or %s) '
            'and not %s)' %(p,optimize_add, optimize_ignore))
            get_strain(optimize_local, temp_clash_A, optimize)
            cmd.set_name(temp_clash_A,'%s_clashes'%p)

    ##### ANIMATE VDW OPTIMIZATION END #####
    #OUTPUT
    if (output):
        CHIS=[ 'CA-CB', 'CB-CG', 'CG-CD', 'CD-CE', 'CE-NZ']
        print("----------------------------------------")
        print("VDW STRAIN REPORT:")
        print("----------------------------------------")
        for p in objects:
            print("Protein:",p)
            print("     Original:", OBJCHIS[p][0][0])
            print("     Modified:", OBJCHIS[p][0][1])
            if optimize==0:
                print("Optimization off!")
            else:
                print("    Optimized:", OBJCHIS[p][0][2])
            print("--------------------")

            #only print if optimize was set
            for q in OBJCHIS[p][1]:
                print("  Residue:", q[0])
                if q[1]==0:
                    print("not part of selection - skipped!")
                    print("--------------------")
                    continue
                x=('    Group: %d'
                   ' Confomer: %d'
                   ' Type: %d' %(q[3][0],q[3][1],q[3][2]))
                print(x)
                if q[1]==4: print("Lysine")
                if q[1]==3: print("N-term-LYS-NZ")
                if q[1]==2: print("N-term-LYS-N")
                if q[1]==1: print("N-term")
                print("--Torsion angles--")
                print("Bond".ljust(6)+\
                      "Original".ljust(14)+\
                      " Optimal".ljust(14))
                if q[1]>2:
                    x=0
                    for c in range(4,0,-1):
                        x=x+1
                        print(CHIS[c].ljust(6)+\
                              str(round(q[2][0][c],3)).ljust(14)+\
                              " "+str(round(q[2][1][c],3)).ljust(14))
                        if x==optimize: break
                else:
                    print("CA-N".ljust(6)+\
                          str(round(q[2][0],3)).ljust(14)+\
                          " "+str(round(q[2][1],3)).ljust(14))
                print("--Local strain--")
                print("Start".ljust(10)+\
                              "Optimal".ljust(10)+\
                              "Difference".ljust(10))
                print ("%.2f".ljust(10)%(q[4][0])+
                       "%.2f".ljust(10)%(q[4][1])+
                       "%.2f".ljust(10)%(q[4][1]-q[4][0])
                      )
                print("--------------------")
        print("----------------------------------------")
        print("required time:",datetime.datetime.now()-start_time)
        print("----------------------------------------")
        #print OBJCHIS

    # rebuild
    pytms_rebuild()

    last_time=log_pytms_prog(output,
    'MDA modification(s) done!')
    # export dictionary if run in python
    if optimize>0: return OBJCHIS
    else: return stored.resi_list

cmd.extend( "mda_modify", mda_modify );
################################################################################
################################################################################





################################################################################
# Automated in silico phosphorylation (Model)
################################################################################
def phosphorylate(
selection='all',
surface_cutoff=0,
mode=0,
color_mod='',
color_base='',
optimize=0,
interpolate=0,
local_radius=10,
interval=30,
states=0,
remove_radius=5,
optimize_ignore='hetatm',
optimize_add='none',
show_clashes=0,
hydrogens=0,
quiet=1
):


    '''
DESCRIPTION

    Introduces Phosphorylation of Serines, Threonines and Tyrosines, toggled by mode

EXAMPLE

    frag SER
    phosphorylate

USAGE

    phosphorylate [ selection [, surface_cutoff [, mode
    [, color_mod [, color_base [, optimize [, interpolate [, local_radius
    [, interval [, states [, remove_radius [, optimize_ignore
    [, optimize_add [, show_clashes [, hydrogens [, quiet ]]]]]]]]]]]]]]]]

ARGUMENTS

    selection: selection to be modified {default: 'all'}
               Residues are automatically sub-selected (see also: mode)!

    surface_cutoff: <float> variable for integrated selection of surface atoms.
                    If set >0, PyTMs will automatically calculate the solvent-
                    accessible surface area, and sub-select all atoms above this
                    cutoff.
                    Note that this operation will create a reserved selection:
                    'pytms_surface_input', but only modify residues if the target atoms
                    are also part of the original selection.
                    see also: PyTMs wiki page or findSurfaceResidues (PyMOL wiki)

    mode: <int> toggles the sub-selection {default=0}
                0: all residues: S/T/Y
                1: Serines only
                2: Threonines only
                3: Tyrosines only
                4: Serines and Threonines (no Tyrosines)

    color_base, color_mod: color names for
                           object selection and modification, respectively
                           {default: ''} = off

    optimize: <boolean> toggles if the dihedral angle of the phospho-group
                        will be adjusted to avoid vdw-clashing;
                        if off the phosphate will be positioned with a default
                        dihedral angle of 180 for SER/THR and 90 for TYR residues
                        {default: False}

    interpolate: <boolean> {default=0} toggles whether the script will
                           attempt to average the angles obtained from equally minimal
                           VDW strains. If off, the closest minimum will be used.

    local_radius: <float> cutoff distance for local refinement
                          defines a local area around each modified residue
                          used for optimizing; {default=10}
                          smaller values speed up calculation,
                          but may result in misplacement
                          (increase if this is the case)

    interval: <float> defines the step (in deg) for the angle in 'optimize':
                       lower values may yield better accuracy but will
                       increase calculation time!
                       {default=30}

    states: <int> number of states used to animate the VDW optimization
                  only takes effect if "optimize" is performed
                  {default: 0} (set 0 or 1 for only final result)

    remove_radius: <float> toggles if heteroatoms are removed
                   after optimization {default: 5.0}
                   will remove all atoms (defined by optimize_ignore) surrounding
                   modified residues (within <value> angstrom)
                   Set to remove_radius=0 to supress any removal
                   Optional: manually remove heteroatoms before/after running the
                   script, if preferred (see also: optimize_ignore)

    optimize_ignore: <str> a string selection defining atoms to
                     be ignored during optimization
                     {default: 'hetatm'}
                     can be used to ignore only specific heteroatoms
                     (e.g. 'resn HOH') or parts of objects
                     # Note that these will be removed, depending on
                     remove_radius

    optimize_add:    <str> a string selection defining additional atoms to
                     be accounted for during optimization other than the parent
                     object {default: 'none'} (useful to account for adjacent objects)

    show_clashes: <boolean> toggles if clashes will be visualized {0}
    # consider removing potentially clashing heteroatoms

    hydrogens: <int>     toggles if the object will have hydrogens or not
                         after modification
                         0  : as is (detection) {default}
                         1  : adds hydrogens;
                         -1 : no hydrogens;
                         #NB! affects vdW optimization!

    quiet: <boolean> toggles output {default: quiet=1}
                     updates on operation and remaining time
                     * NB! print appears first in the console

    '''
    ##### BEGINNING OF SCRIPT #####

    last_time=start_time=datetime.datetime.now()
    try:
        surface_cutoff=abs(float(surface_cutoff))
    except:
        raise Exception("Input error!\n Illegal value for surface cutoff!")
    selection=pytms_get_sele(selection, surface_cutoff)

    # argument settings
    try:
        color_base, color_mod = str(color_base), str(color_mod)
        if color_base: cmd.color(color_base,selection)
        if color_mod: cmd.color(color_mod,''+selection+' and (name OP1 extend 2)')
        mode=abs(int(mode))
        if (mode not in [0,1,2,3,4]):
            raise Exception("mode out of range!")
        local_radius=abs(float(local_radius))
        states=abs(int(states))
        remove_radius=abs(float(remove_radius))
        interval=abs(float(interval))
        # transform boolean to reduce checking loops
        interpolate=bool((str(interpolate)!='False'))
        hydrogens=eval_hydrogens(selection, hydrogens)
        optimize=bool(int(optimize))
        show_clashes=bool(int(show_clashes))
        output=bool(not int(quiet))
    except:
        raise Exception("Input error!\n Please check the input parameters!")

    if ((interval<0.1) or (interval>360)):
        raise Exception("interval value out of range!")

    # SELECTION LIST ###########################################################

    # create a empty list for appending
    stored.resi_list = []
    stored.temp = []

    # for each object and chain --> create selection lists
    for p in names:
        # for each object and chain get resi of residues
        cmd.iterate((
        '(%s) and ((resn SER) or (resn THR) or (resn TYR) or '
        '(resn NIY)) and (name CA)' %(p)),
        'stored.resi_list.append("(%s and resi "+str(resi)+")")'  %p
        )


    # kick out non-selected or missing
    for p in list(stored.resi_list):
        select_part(temp_sel,'%s'%p,'(name CB)','(name CA)')
        # select the one oxygen, and not those in a possible nitration
        if (cmd.count_atoms((
        '((%s) and (%s) and (%s) and (elem O) )and '
        'not (neighbor (elem N))' %(selection,p,temp_sel)))==0):
            stored.resi_list.remove(p)

    ##### finished creating selection lists #####

    # premature exit if nothing to process:
    if stored.resi_list==[]:
        print("PyTMs: No modifyable residues found in selection!")
        return False

    last_time=log_pytms_prog(output,
    'Initialized phosphorylation!')

    # VDW PREP
    # OBJCHI (dictionary) structure
    # obj: [umod. strain, mod. strain. opt. strain],
    # + [n*([Res,type,
    #       [[Start angles],[opt. angles],[diff]],
    #        [local stain (start), local stain (opt)]])]
    OBJCHIS = {}

    # get unmodified VDW strain
    for p in objects:
        if (hydrogens):
            cmd.h_add(p)
        else:
            cmd.remove('((%s) and elem H)' %p)
        OBJCHIS[p] = ([[get_strain('(%s or %s) and not (%s)'%(p,optimize_add,optimize_ignore),
        temp_clash_A,optimize),0,0],[]])


    ##### BUILDER #####
    isused_list=list(stored.resi_list)
    # isused_list is <int>
    # 0= not used; 1: serine; 2: threonine; 3: tyrosine

    # cycles through residues
    listcount=-1
    for p in stored.resi_list:
        cmd.unpick()
        listcount=listcount+1
        isused_list[listcount]=0

        stored.residue=''
        cmd.iterate('%s and name CA' %p, 'stored.residue=str(resn)')

        ## Skip or not? and change of name
        skip_resi=0
        if (stored.residue=='SER'):
            if mode not in [0,1,4]: continue
            cmd.alter(p,'resn="SEP"')
            isused_list[listcount]=1
        elif (stored.residue=='THR'):
            if mode not in [0,2,4]: continue
            cmd.alter(p,'resn="TPO"')
            isused_list[listcount]=2
        elif (stored.residue=='TYR'):
            if mode not in [0,3]: continue
            cmd.alter(p,'resn="PTR"')
            isused_list[listcount]=3
        elif (stored.residue=='NIY'):
            if mode not in [0,3]: continue
            cmd.alter(p,'resn="PNIY"')
            isused_list[listcount]=3
        else:
            print ('Skipped unidentified residue [%s]!' %stored.residue)
            print ('[%s]!' %p)
            continue

        select_part(temp_sel,'%s'%p,'(name CB)','(name CA)')

        # select the one oxygen, and not those in a possible nitration
        cmd.edit('(%s) and (elem O) and not (neighbor (elem N))' %temp_sel)

        stored.oname=''
        cmd.iterate('(%s) and (elem O) and not (neighbor (elem N))' %temp_sel, 'stored.oname=str(name)')

	# remove potential H atoms
        cmd.remove('(%s) and ((neighbor pk1) and (elem H))' %p)

        cmd.attach('P','3','4')
        cmd.alter('%s and (elem P)' %p, 'name="P"')
        cmd.edit('%s and (elem P)' %p)
        cmd.attach('O','2','2')
        cmd.alter(('%s and (neighbor pk1) and not (name %s or '
        '(name O1P) or (name O2P) or (name O3P))' %(p,stored.oname)),
        'name="O3P"')
        for c in ['O1P','O2P']:
            cmd.attach('O','2','1')
            cmd.alter(('%s and (neighbor pk1) and not (name %s or '
            '(name O1P) or (name O2P) or (name O3P))' %(p,stored.oname)),
            'name="%s"' %c)
            if version_ok:
                cmd.alter('%s and (neighbor pk1) and (name %s)' %(p,c),
                'formal_charge=-1')
        cmd.edit('%s and (name O3P)' %p)
        cmd.unbond('pk1', '%s and elem P' %p)
        cmd.bond('pk1', '%s and elem P' %p, '2')


        # default in trans unless optimize==True
        if optimize==False:
            x=180
            if (stored.oname=='OH'): x=90
            cmd.select(temp_sel, '%s and ((elem p) extend 4)'%p)
            cmd.set_dihedral(
            '%s and ((name CA) or (name CE1))' %temp_sel,
            '%s and ((name CB) or (name CZ))' %temp_sel,
            '%s and (name %s)' %(temp_sel, stored.oname),
            '%s and (elem P)' %temp_sel,
            x)

        #Log
        last_time=log_pytms_prog(output,
        'Modified: %s'%get_resi_macro_name(p),listcount+1,len(stored.resi_list))

        cmd.unpick()
    # End of resi cycle
    if (hydrogens):
        for p in objects: cmd.h_add(p)
    else:
        for p in objects:
            cmd.remove('((%s) and elem H)' %p)

    # rebuild
    pytms_rebuild()

    if (not version_ok):
        print("NB! Charge assignment is incorrect due do outdated PyMOL version!")
        if hydrogens:
            print("hydrogens were added to oxygens that should have negative charge!")
        print("Recommended: PyMOL Version >=1.7")

    # coloring
    if color_base:
        cmd.color(color_base,
        '%s or ((%s extend 1) and hydrogens)' %(selection, selection))
    for p in stored.resi_list:
        # enable selection "p.PTM in phosphorylation"
        if version_ok:
            cmd.alter('(%s) and ((name P) or (name O1P) or (name O2P) or (name O3P))' %p,
            'p.PTM="phosphorylation"')
        if color_base:
            cmd.color(color_base,
            '(%s)' %(p))
        if color_mod:
            cmd.color(color_mod,
           '(%s) and ((name P) or (name O1P) or (name O2P) or (name O3P))' %p)


    last_time=log_pytms_prog(output,
    'Modifications complete!')

    ##### VDW OPTIMIZATION #####

    # get modified strain
    for p in objects:
        OBJCHIS[p][0][1]=OBJCHIS[p][0][2]= get_strain(
        '(%s or %s) and not (%s)'%(p,optimize_add,optimize_ignore),
        temp_clash_A,optimize)

    if optimize:
        last_time=log_pytms_prog(output,
        'Starting optimization (be patient ;-))!')

        # Chi angles of residues
        CHIS = [
                ['name CA','name CB','name OG','elem P'],
                ['name CA','name CB','name OG1','elem P'],
                ['name CE1','name CZ','name OH','elem P']
               ]

        progress=0.0
        ##### START OF P IN OBJECTS

        for p in objects:
            listcount=-1
            for q in stored.resi_list:
                listcount=listcount+1
                # pass if residue not in object or unknown
                if isused_list[listcount]==0: continue
                if (cmd.count_atoms('%s and %s'%(p,q))==0): continue

                #reset VALS with starting values
                VALS=([q, isused_list[listcount], [0.0,0.0,0.0], [0,0]])

                # gets the starting dihedral
                VALS[2][0]=VALS[2][1]=(
                cmd.get_dihedral(
                '%s and %s' %(q, CHIS[VALS[1]-1][0]),
                '%s and %s' %(q, CHIS[VALS[1]-1][1]),
                '%s and %s' %(q, CHIS[VALS[1]-1][2]),
                '%s and %s' %(q, CHIS[VALS[1]-1][3])
                ))

                OBJCHIS[p][1].append(list(VALS))
            # done writing starting values
            # end q in stored.resi_list

            # optimizing
            for q in OBJCHIS[p][1]:
                # is not skipped resi
                if q[1]==0:
                    continue

                x=OBJCHIS[p][1].index(q)
                optstrain=0
                # get info
                stored.residue=''
                cmd.iterate('%s and name CA' %q[0], 'stored.residue=str(resi)')
                # define local and get starting strain
                # base selection
                optimize_local= ('((byres (%s expand %f)) and (%s or %s))'
                %(q[0], local_radius, p,optimize_add))
                # without ignored
                optimize_local= ('(%s and not (%s))'
                %(optimize_local, optimize_ignore))
                # without adducts optimized later
                optimize_local= ('(%s and not ((name P+O1P+O2P+O3P) '
                'and %s and ((resi %s-) and not %s)))'
                %(optimize_local, p, stored.residue, q[0]))
                # make selection
                cmd.select(temp_sel,optimize_local)

                optstrain=get_strain(temp_sel,temp_clash_A,optimize)
                # optimal local strain
                OBJCHIS[p][1][x][3][0]=OBJCHIS[p][1][x][3][1]=optstrain

                first=ang_4=q[2][0]
                last=first

                while (ang_4 < ((q[2][0])+360)):
                    # set new angles
                    cmd.set_dihedral(
                    '%s and %s' %(q[0], CHIS[q[1]-1][0]),
                    '%s and %s' %(q[0], CHIS[q[1]-1][1]),
                    '%s and %s' %(q[0], CHIS[q[1]-1][2]),
                    '%s and %s' %(q[0], CHIS[q[1]-1][3]),
                    ang_4
                    )

                    # get new vdw strain
                    strain=get_strain(temp_sel,temp_clash_A,optimize)

                    # check if better
                    if strain<optstrain:
                        optstrain=strain
                        first = ang_4
                    if strain==optstrain:
                        last=ang_4
                    if optstrain==0:
                        #absolute min. (exit calc.)
                        ang_4=ang_4+361.0
                        print ('!!! PyTMs Warning: '
                               'local zero strain '
                               'occured for: %s'%q[0])
                    ang_4=ang_4+interval

                # set to optimal and write optimal to dictionary
                if interpolate:
                    ang_4=to180((to360(first)+to360(last))/2)
                else:
                    if (
                        (abs(to180(to360(first)-to360(q[2][0])))<=
                         abs(to180(to360(last)-to360(q[2][0]))))
                    ):
                        ang_4=to180(first)
                    else:
                        ang_4=to180(last)

                cmd.set_dihedral(
                '%s and %s' %(q[0], CHIS[q[1]-1][0]),
                '%s and %s' %(q[0], CHIS[q[1]-1][1]),
                '%s and %s' %(q[0], CHIS[q[1]-1][2]),
                '%s and %s' %(q[0], CHIS[q[1]-1][3]),
                ang_4
                )

                if interpolate:
                # check that average is good
                    strain=get_strain(temp_sel,temp_clash_A,optimize)
                    if strain>optstrain:
                        if (
                            (abs(to180(to360(first)-to360(q[2][0])))<=
                             abs(to180(to360(last)-to360(q[2][0]))))
                        ):
                            ang_4=to180(first)
                        else:
                            ang_4=to180(last)
                        cmd.set_dihedral(
                        '%s and %s' %(q[0], CHIS[q[1]-1][0]),
                        '%s and %s' %(q[0], CHIS[q[1]-1][1]),
                        '%s and %s' %(q[0], CHIS[q[1]-1][2]),
                        '%s and %s' %(q[0], CHIS[q[1]-1][3]),
                        ang_4
                        )

                # write angle difference to dictionary
                OBJCHIS[p][1][x][2][1]=ang_4
                OBJCHIS[p][1][x][2][2]=(
                to180(to360(ang_4)-to360(OBJCHIS[p][1][x][2][0]))
                )
                # optimal local strain
                OBJCHIS[p][1][x][3][1]=optstrain

                # progress
                progress=progress+1
                log_pytms_prog(output,
                'Optimized: %s'%get_resi_macro_name(q[0]), progress, len(stored.resi_list),last_time)

            # END OF Q IN OBJCHIS[p][1]
            # write optimized strain
            OBJCHIS[p][0][2]= get_strain(
            '(%s or %s) and not (%s)'%(p,optimize_add,optimize_ignore),
            temp_clash_A,optimize)

        ##### END OF P IN OBJECTS
        last_time=log_pytms_prog(output,
        'Optimization complete!')
    ##### VDW OPTIMIZATION END #####

    # remove atoms near modified residues
    for p in objects:
        for q in OBJCHIS[p][1]:
            optimize_local=('(((byres (%s expand %f)) and (%s or %s)) '
            'and (%s))' %(q[0], remove_radius, p,optimize_add, optimize_ignore))
            cmd.remove(optimize_local)

    ##### ANIMATE VDW OPTIMIZATION START #####
    progress=0.0
    if ((states>0) and (optimize)):
        last_time=log_pytms_prog(output,
        'Initiating animation calulation!')
        for p in objects:
            cmd.delete(temp_obj)
            cmd.delete(temp_clash_A)

            for n in range(1,states+1):
                cmd.set('state',1,'%s'%p)
                frac=(float(n-1)/float(states-1))
                #float is due to intervals during animation:
                #int would be rounded
                for q in OBJCHIS[p][1]:
                    # skipped resis
                    if q[1]==0: continue

                    # recalculate angle
                    ang_4=to180(to360(q[2][0])+(q[2][2]*frac))

                    cmd.set_dihedral(
                    '%s and %s' %(q[0], CHIS[q[1]-1][0]),
                    '%s and %s' %(q[0], CHIS[q[1]-1][1]),
                    '%s and %s' %(q[0], CHIS[q[1]-1][2]),
                    '%s and %s' %(q[0], CHIS[q[1]-1][3]),
                    ang_4
                    )

                #Done with residues
                # update temporary models
                cmd.create(temp_obj,p,1,n)
            #Done with states

            # when done with all states --> rename object
            cmd.delete(p)
            cmd.create(p,temp_obj)

            # clashes
            if (show_clashes):
                optimize_local=('((%s or %s)'
                'and not %s)' %(p,optimize_add, optimize_ignore))
                cmd.create('%s_clashes'%p, optimize_local)
                # activate clash visibility
                for n in range(1,states+1):
                    cmd.set('state',n)
                    cmd.sculpt_activate('%s_clashes'%p)
                    cmd.show_as('cgo', '%s_clashes'%p)
                    cmd.set('sculpt_vdw_vis_mode', 1, '%s_clashes'%p)
                    cmd.sculpt_iterate('%s_clashes'%p, cycles=1)

            progress=progress+1
            log_pytms_prog(output,
            'Animating objects!', progress, len(objects), last_time)

    elif (show_clashes):
        # no animation but clashes
        for p in objects:
            optimize_local=('((%s or %s) '
            'and not %s)' %(p,optimize_add, optimize_ignore))
            get_strain(optimize_local, temp_clash_A, optimize)
            cmd.set_name(temp_clash_A,'%s_clashes'%p)

    ##### ANIMATE VDW OPTIMIZATION END #####
    #OUTPUT
    if (output):
        print("----------------------------------------")
        print("VDW STRAIN REPORT:")
        print("----------------------------------------")
        for p in objects:
            print("Protein:",p)
            print("     Original:", OBJCHIS[p][0][0])
            print("     Modified:", OBJCHIS[p][0][1])
            if not optimize:
                print("Optimization off!")
            else:
                print("    Optimized:", OBJCHIS[p][0][2])
            print("--------------------")

            #only prints if optimize was set
            for q in OBJCHIS[p][1]:
                print("  Residue:", q[0])
                if q[1]==0:
                    print("not part of selection - skipped!")
                    print("--------------------")
                    continue
                if q[1]==1: print("Type: SERINE")
                if q[1]==2: print("Type: THREONINE")
                if q[1]==3: print("Type: TYROSINE")
                print("--Torsion angles--")
                print("Bond".ljust(6)+\
                      "Original".ljust(14)+\
                      " Optimal".ljust(14))
                print("C-O".ljust(6)+\
                      str(round(q[2][0],3)).ljust(14)+\
                      " "+str(round(q[2][1],3)).ljust(14))
                print("--Local strain--")
                print("Start".ljust(10)+\
                              "Optimal".ljust(10)+\
                              "Difference".ljust(10))
                print ("%.2f".ljust(10)%(q[3][0])+
                       "%.2f".ljust(10)%(q[3][1])+
                       "%.2f".ljust(10)%(q[3][1]-q[3][0])
                      )
                print("--------------------")
            print("----------------------------------------")
            print("required time:",datetime.datetime.now()-start_time)
            print("----------------------------------------")
            #print OBJCHIS

    # rebuild and remove temporary
    pytms_rebuild()

    last_time=log_pytms_prog(output,
    'Phosphorylation done!')
    if optimize: return OBJCHIS
    else: return stored.resi_list

cmd.extend( "phosphorylate", phosphorylate );
################################################################################
################################################################################





################################################################################
# Automated in silico Cysteine oxidation (Model)
################################################################################
def oxidize_cys(
selection='all',
surface_cutoff=0,
mode=1,
include_SEC=0,
disulfides=0,
show_clashes=0,
color_mod='',
color_base='',
hydrogens=0,
quiet=1
):


    '''
DESCRIPTION

    Modifies the Cysteins of a selection by oxidation

EXAMPLE

    frag CYS
    oxidize_cys

USAGE

    oxidize_cys [ selection [, surface_cutoff [, mode [, include_SEC
    [, disulfides [, show_clashes [, color_mod [, color_base [, hydrogens
    [, quiet ]]]]]]]]]]

ARGUMENTS

    selection: selection to be modified {default: 'all'}
               Cysteines etc. are automatically sub-selected
               (see also: include_SEC and disulfides)!

    surface_cutoff: <float> variable for integrated selection of surface atoms.
                    If set >0, PyTMs will automatically calculate the solvent-
                    accessible surface area, and sub-select all atoms above this
                    cutoff.
                    Note that this operation will create a reserved selection:
                    'pytms_surface_input', but only modify residues if the target atoms
                    are also part of the original selection.
                    see also: PyTMs wiki page or findSurfaceResidues (PyMOL wiki)

    mode: <int> toggles the oxidation type and configuration {default=1}
                1: s-hydroxy-cysteine (S-OH, CSO)
                2: s-oxy-cysteine (S=O, CSX)
                3: cysteine-s-dioxide, R-configuration (SOOH, CSW)
                4: cysteine-s-dioxide, S-configuration (SOOH, CSW)

    include_SEC: <boolean> toggles if the selenocysteines (SEC/CSE) will be included or not
                 {default: 0}

    disulfides: <boolean> toggles if disulfide bridges will be affected or not
                {default: 0} NB! setting this True may give erroneous results

    show_clashes: <boolean> toggles if clashes will be visualized {0}
    # consider removing potentially clashing heteroatoms

    color_base, color_mod: color names for
                           object selection and modification, respectively
                           {default: ''} = off

    hydrogens: <int>     toggles if the object will have hydrogens or not
                         after modification
                         0  : as is (detection) {default}
                         1  : adds hydrogens;
                         -1 : no hydrogens;

    quiet: <boolean> toggles output {default: quiet=1}
                     updates on operation and remaining time
                     * NB! print appears first in the console

    '''
    ##### BEGINNING OF SCRIPT #####

    last_time=start_time=datetime.datetime.now()
    try:
        surface_cutoff=abs(float(surface_cutoff))
    except:
        raise Exception("Input error!\n Illegal value for surface cutoff!")
    selection=pytms_get_sele(selection, surface_cutoff)

    # argument settings
    try:
        color_base, color_mod = str(color_base), str(color_mod)
        if color_base: cmd.color(color_base,selection)
        if color_mod:
            cmd.color(color_mod,
            ('(%s and (resn CS*)) and '
            '((name OD*) or ((neighbor (name OD*)) and elem H))' %selection))
        mode=abs(int(mode))
        if (mode not in [1,2,3,4]):
            raise Exception("mode out of range!")
        # transform boolean to reduce checking loops
        show_clashes=bool(int(show_clashes))
        include_SEC=bool(int(include_SEC))
        disulfides=bool(int(disulfides))
        hydrogens=eval_hydrogens(selection, hydrogens)
        output=bool(not int(quiet))
    except:
        raise Exception("Input error!\n Please check the input parameters!")

    # SELECTION LIST ###########################################################

    # rename to allow later conversion
    if include_SEC:
        cmd.alter('byres (%s and %s and (resn SCE or resn SEC))' %(p,selection),
        'resn="CYS"')

    # create a empty list for appending
    stored.resi_list = []
    stored.temp = []

    # for each object and chain --> create selection lists
    for p in names:
        # for each object and chain get resi of residues
        cmd.iterate((
        '(%s) and (resn CYS) '
        'and (name CA)' %(p)),
        'stored.resi_list.append("(%s and resi "+str(resi)+")")'  %p
        )

    # kick out non-selected or missing
    for p in list(stored.resi_list):
        if (cmd.count_atoms('(%s) and (%s) and (elem S or elem SE)' %(selection,p))==0):
            stored.resi_list.remove(p)

    # kick out disulfides (if not set!)
    if not disulfides:
        for p in list(stored.resi_list):
            if (cmd.count_atoms((
            '(neighbor ((%s) and (elem S or elem SE))) '
            'and (elem S or elem SE)' %(p)))!=0):
                # not regular cysteine
                stored.resi_list.remove(p)


    ##### finished creating selection lists #####

    # premature exit if nothing to process:
    if stored.resi_list==[]:
        print("PyTMs: No modifyable residues found in selection!")
        return False

    # VDW PREP
    # OBJCHI (dictionary) structure
    # obj: [umod. strain, mod. strain.]
    OBJCHIS = {}

    if show_clashes:
        # get unmodified VDW strain
        log_pytms_prog(output, 'Calculating base VdW strain!')
        for p in objects:
            if (hydrogens):
                cmd.h_add(p)
            else:
                cmd.remove('((%s) and elem H)' %p)

            OBJCHIS[p] = [get_strain('(%s)'%(p), temp_clash_A, 0),0]

    last_time=log_pytms_prog(output,
    'Initialized cysteine oxidation!')

    ##### BUILDER #####

    resnames=['','CSO','CSX','CSW','CSW']

    # cycles through residues
    listcount=-1
    for p in stored.resi_list:
        cmd.unpick()
        listcount=listcount+1

        # get info on residue being processed
        stored.residue=''
        cmd.iterate('%s and name CA' %p, 'stored.residue=str(resn)')


        ## change name
        cmd.alter(p,'resn="%s"' %(resnames[mode]))

        # get rid of bound non-chain neighbors
        # but conserve potential S
        cmd.remove(('%s and (neighbor (elem SE or elem S)) '
        'and not (elem C or elem S)' %p))
        # edit
        cmd.edit('%s and (elem SE or elem S)' %p)
        if mode<3:
            cmd.attach('O','2','2')
            cmd.alter('%s and (neighbor (elem SE or elem S)) and (elem O)' %p,'name="OD"')
            if mode==1:
                cmd.edit('%s and (name OD)' %p)
                cmd.attach('H','1','1')
            if mode==2:
                cmd.unbond('pk1', '%s and name OD' %p)
                cmd.bond('pk1', '%s and name OD' %p,'2')
        else:
            # oxidize 3x
            for c in ['A','B','2']:
                cmd.attach('O','2','2')
                cmd.alter(('(%s and (neighbor (elem SE or elem S)) and (elem O)) '
                ' and not (name OD*)' %p),
                'name="OD%s"' %c)
            # double bond
            cmd.unbond('pk1', '%s and name OD2' %p)
            cmd.bond('pk1', '%s and name OD2' %p,'2')
            # R or S?
            # one oxygen will be removed
            # get dihedrals
            # the R oxygen is in front of the plane formed by CB-S-OD2 (smaller angle)
            # the S oxygen is behind the plane formed by CB-S-OD2 (larger angle)

            # zero plane
            dihedOD2=to360(cmd.get_dihedral(
            '%s and (name CA)' %p,
            '%s and (name CB)' %p,
            '%s and (elem S or elem SE)' %p,
            '%s and (name OD2)' %p))

            # first oxygen
            dihedODA=to360(cmd.get_dihedral(
            '%s and (name CA)' %p,
            '%s and (name CB)' %p,
            '%s and (elem S or elem SE)' %p,
            '%s and (name ODA)' %p))
            dihedODA=to360(dihedODA-dihedOD2)

            # other oxygen
            dihedODB=to360(cmd.get_dihedral(
            '%s and (name CA)' %p,
            '%s and (name CB)' %p,
            '%s and (elem S or elem SE)' %p,
            '%s and (name ODB)' %p))
            dihedODB=to360(dihedODB-dihedOD2)

            if mode==3: dihedOD2='ODA'
            else: dihedOD2='ODB'
            if dihedODA<dihedODB:
                # is in opposite configuration
                if mode==3: dihedOD2='ODB'
                else: dihedOD2='ODA'
            # S: deletes lower
            # R: deletes upper
            cmd.remove('%s and (name %s)' %(p,dihedOD2))
            # rename to OD1
            cmd.alter(('%s and (neighbor (name SG)) and '
            '(elem O) and not (name OD2)' %p),
            'name="OD1"')
            # hydration of OD1
            cmd.edit('%s and (name OD1)' %p)
            cmd.attach('H','1','1')
        # END IF MODE>=3

        #in PyMol valence corresponds more to
        #the number of neighbors than number of bonds
        if mode<3:
            cmd.set_geometry('%s and (name SG)' %p,'4','2')
        else:
            cmd.set_geometry('%s and (name SG)' %p,'4','3')
        cmd.h_add('%s' %p)

        #Log
        last_time=log_pytms_prog(output,
        'Modified: %s'%get_resi_macro_name(p),listcount+1,len(stored.resi_list))

        cmd.unpick()
    # End of resi cycle
    if (hydrogens):
        for p in objects: cmd.h_add(p)
    else:
        for p in objects:
            cmd.remove('((%s) and elem H)' %p)

    # clashes
    if show_clashes:
        # get modified strain
        log_pytms_prog(output, 'Calculating modified VdW strain!')
        if output:
            print('STRAIN REPORT:')
            print('OBJECT','NATIVE_STRAIN','MODIFIED_STRAIN', 'DIFFERENCE')
        for p in objects:
            OBJCHIS[p][1] = get_strain('(%s)'%(p), temp_clash_A, 0)
            cmd.set_name(temp_clash_A,'%s_clashes'%p)
            if output:
                print('%s'%p,OBJCHIS[p][0],OBJCHIS[p][1],OBJCHIS[p][1]-OBJCHIS[p][0])

    # rebuild
    pytms_rebuild()

    # coloring
    if color_base:
        cmd.color(color_base,
        '%s or ((%s extend 1) and hydrogens)' %(selection, selection))
    for p in stored.resi_list:
        # enable selection "p.PTM in oxidation_cys"
        if version_ok:
            cmd.alter('%s and ((name OD*) or ((neighbor (name OD*)) and elem H))' %p,
            'p.PTM="oxidation_cys"')
        if color_base:
            cmd.color(color_base,
            '(%s)' %(p))
        if color_mod:
            cmd.color(color_mod,
            ('(%s and (resn CS*)) and '
            '((name OD*) or ((neighbor (name OD*)) and elem H))' %p))
    #exit
    last_time=log_pytms_prog(output,
    'Cysteine oxidation complete!')
    return [stored.resi_list, OBJCHIS]

cmd.extend( "oxidize_cys", oxidize_cys );
################################################################################
################################################################################





################################################################################
# Automated in silico Methionine oxidation (Model)
################################################################################
def oxidize_met(
selection='all',
surface_cutoff=0,
mode=2,
convert_MSE=1,
show_clashes=0,
color_mod='',
color_base='',
hydrogens=0,
quiet=1
):


    '''
DESCRIPTION

    Modifies the Methionines of a selection by oxidation

EXAMPLE

    frag MET
    oxidize_met

USAGE

    oxidize_met [ selection [, surface_cutoff [, mode [, convert_MSE
    [, show_clashes [, color_mod [, color_base [, hydrogens [, quiet ]]]]]]]]]

ARGUMENTS

    selection: selection to be modified {default: 'all'}
               Methionines/Selenomethionines are automatically sub-selected (see also. convert_MSE)!

    surface_cutoff: <float> variable for integrated selection of surface atoms.
                    If set >0, PyTMs will automatically calculate the solvent-
                    accessible surface area, and sub-select all atoms above this
                    cutoff.
                    Note that this operation will create a reserved selection:
                    'pytms_surface_input', but only modify residues if the target atoms
                    are also part of the original selection.
                    see also: PyTMs wiki page or findSurfaceResidues (PyMOL wiki)

    mode: <int> toggles the oxidation type and configuration {default=2}
                1: methionine-R-sulfoxide
                2: methionine-S-sulfoxide
                3: methionine-sulfone

    convert_MSE: <boolean> toggles if the selenomethionines will be converted
                 NB!: this removes the Selen and introduces a Sulfur into the model
                 {default: 1}

    show_clashes: <boolean> toggles if clashes will be visualized {0}
    # consider removing potentially clashing heteroatoms

    color_base, color_mod: color names for
                           object selection and modification, respectively
                           {default: ''} = off

    hydrogens: <int>     toggles if the object will have hydrogens or not
                         after modification
                         0  : as is (detection) {default}
                         1  : adds hydrogens;
                         -1 : no hydrogens;

    quiet: <boolean> toggles output {default: quiet=1}
                     updates on operation and remaining time
                     * NB! print appears first in the console

    '''
    ##### BEGINNING OF SCRIPT #####

    last_time=start_time=datetime.datetime.now()
    try:
        surface_cutoff=abs(float(surface_cutoff))
    except:
        raise Exception("Input error!\n Illegal value for surface cutoff!")
    selection=pytms_get_sele(selection, surface_cutoff)

    # argument settings
    try:
        color_base, color_mod = str(color_base), str(color_mod)
        if color_base: cmd.color(color_base,selection)
        if color_mod: cmd.color(color_mod,''+selection+' and (resn SME) and name OS*')
        mode=abs(int(mode))
        if (mode not in [1,2,3]):
            raise Exception("mode out of range!")
        # transform boolean to reduce checking loops
        show_clashes=bool(int(show_clashes))
        convert_MSE=bool(int(convert_MSE))
        hydrogens=eval_hydrogens(selection, hydrogens)
        output=bool(not int(quiet))
    except:
        raise Exception("Input error!\n Please check the input parameters!")

    # SELECTION LIST ###########################################################

    # MSE to MET rename to allow later conversion
    if convert_MSE:
        cmd.alter('byres (%s and resn MSE)' %(selection), 'resn="MET"')


    # create a empty list for appending
    stored.resi_list = []
    stored.temp = []

    # for each object and chain --> create selection lists
    for p in names:
        # for each object and chain get resi of methionine residues
        cmd.iterate((
        '(%s) and (resn MET) '
        'and (name CA)' %(p)),
        'stored.resi_list.append("(%s and resi "+str(resi)+")")'  %p
        )

    # kick out non-selected or missing
    for p in list(stored.resi_list):
        if (cmd.count_atoms('(%s) and (%s) and (elem S or elem SE)' %(selection,p))==0):
            stored.resi_list.remove(p)

    ##### finished creating selection lists #####

    # premature exit if nothing to process:
    if stored.resi_list==[]:
        print("PyTMs: No modifyable residues found in selection!")
        return False

    # VDW PREP
    # OBJCHI (dictionary) structure
    # obj: [umod. strain, mod. strain.]
    OBJCHIS = {}

    if show_clashes:
        # get unmodified VDW strain
        log_pytms_prog(output, 'Calculating base VdW strain!')
        for p in objects:
            if (hydrogens):
                cmd.h_add(p)
            else:
                cmd.remove('((%s) and elem H)' %p)

            OBJCHIS[p] = [get_strain('(%s)'%(p), temp_clash_A, 0),0]

    last_time=log_pytms_prog(output,
    'Initialized methionine oxidation!')

    ##### OXM BUILDER #####
    # cycles through residues
    listcount=-1
    for p in stored.resi_list:
        cmd.unpick()
        listcount=listcount+1

        # get info on residue being processed
        ismse=False
        if cmd.count_atoms('%s and name SE' %p)==1: ismse=True
        # continue if MSE and convert_MSE is off
        if ((ismse) and not (convert_MSE)):
            continue

        #else: go on
        stored.residue=''
        cmd.iterate('%s and name CA' %p, 'stored.residue=str(resn)')


        ## change name
        # is MET or MSE
        cmd.alter(p,'resn="SME"')

        # build SME
        # first convert MSE (already skipped if convert is off)
        # replace SE or S to S with 4 or 8 valences
        x=cmd.get_dihedral(
        '%s and (name CA)' %p,
        '%s and (name CB)' %p,
        '%s and (name CG)' %p,
        '%s and (elem SE or elem S)' %p)
        cmd.remove('%s and (elem SE or elem S)' %p)
        cmd.edit('%s and (name CG)' %p)

        # temporary valence!
        cmd.attach('S','4','4')
        cmd.alter('%s and (neighbor (name CG)) and (elem S)' %p,'name="SD"')
        # inherit color of CG!
        stored.residue=''
        cmd.iterate('%s and name CG' %p, 'stored.residue=str(color)')
        cmd.color(stored.residue,'%s and name SD' %p)

        # reset dihedral to be sure
        cmd.set_dihedral(
        '%s and (name CA)' %p,
        '%s and (name CB)' %p,
        '%s and (name CG)' %p,
        '%s and (name SD)' %p,
        x)
        cmd.bond('%s and (name SD)' %p, '%s and (name CE)' %p,'1')
        # the S is now in position
        # R or S first (full oxidation, then removal)
        cmd.edit('%s and name SD' %p)
        cmd.attach('O','2','2')
        cmd.alter('%s and (neighbor (name SD)) and (elem O)' %p,'name="OS1"')
        cmd.attach('O','2','2')
        cmd.alter('%s and (neighbor (name SD)) and (elem O) and not (name OS1)' %p,
        'name="OS2"')
        # double bonds
        cmd.unbond('pk1', '%s and name OS1' %p)
        cmd.bond('pk1', '%s and name OS1' %p,'2')
        cmd.unbond('pk1', '%s and name OS2' %p)
        cmd.bond('pk1', '%s and name OS2' %p,'2')
        # it is now a sulfone (mode==3)
        if mode!=3:
            # one oxygen will be removed
            # get dihedrals
            dihedCE=to360(cmd.get_dihedral(
            '%s and (name CB)' %p,
            '%s and (name CG)' %p,
            '%s and (name SD)' %p,
            '%s and (name CE)' %p))
            dihedOS1=to360(cmd.get_dihedral(
            '%s and (name CB)' %p,
            '%s and (name CG)' %p,
            '%s and (name SD)' %p,
            '%s and (name OS1)' %p))
            dihedOS1=to360(dihedOS1-dihedCE)
            dihedOS2=to360(cmd.get_dihedral(
            '%s and (name CB)' %p,
            '%s and (name CG)' %p,
            '%s and (name SD)' %p,
            '%s and (name OS2)' %p))
            dihedOS2=to360(dihedOS2-dihedCE)
            # determine which oxygen to delete
            # picture an orientation looking at the methionine sulfur
            # the oxygens aligned vertically and facing towards the viewer
            # and main chain to the left
            # upper=S (smaller angle), lower=R (larger angle)
            # mode==1 (R) --> default OS2 (delete upper)
            if mode==1: dihedCE='OS2'
            else: dihedCE='OS1'
            if dihedOS2<dihedOS1:
                # OS1&2 in opposite configuration
                if mode==1: dihedCE='OS1'
                else: dihedCE='OS2'
            # S: deletes lower
            # R: deletes upper
            cmd.remove('%s and (name %s)' %(p,dihedCE))
            # rename to OS1
            cmd.alter('%s and (neighbor (name SD)) and (elem O)' %p,
            'name="OS1"')

        if mode==3:
            cmd.set_geometry('%s and (name SD)' %p,'4','4')
        else:
            cmd.set_geometry('%s and (name SD)' %p,'4','3')
        cmd.h_add('%s' %p)

        #Log
        last_time=log_pytms_prog(output,
        'Modified: %s'%get_resi_macro_name(p),listcount+1,len(stored.resi_list))

        cmd.unpick()
    # End of resi cycle
    if (hydrogens):
        for p in objects: cmd.h_add(p)
    else:
        for p in objects:
            cmd.remove('((%s) and elem H)' %p)

    # clashes
    if show_clashes:
        # get modified strain
        log_pytms_prog(output, 'Calculating modified VdW strain!')
        if output:
            print('STRAIN REPORT:')
            print('OBJECT','NATIVE_STRAIN','MODIFIED_STRAIN', 'DIFFERENCE')
        for p in objects:
            OBJCHIS[p][1] = get_strain('(%s)'%(p), temp_clash_A, 0)
            cmd.set_name(temp_clash_A,'%s_clashes'%p)
            if output:
                print('%s'%p,OBJCHIS[p][0],OBJCHIS[p][1],OBJCHIS[p][1]-OBJCHIS[p][0])

    # rebuild
    pytms_rebuild()

    # coloring
    if color_base:
        cmd.color(color_base,
        '%s or ((%s extend 1) and hydrogens)' %(selection, selection))
    for p in stored.resi_list:
        # enable selection "p.PTM oxidation_met"
        if version_ok:
            cmd.alter('(%s and (resn SME) and (name OS*))' %p,
            'p.PTM="oxidation_met"')
        if color_base:
            cmd.color(color_base,
            '(%s)' %(p))
        if color_mod:
            cmd.color(color_mod,
            '(%s and (resn SME) and (name OS*))' %(p))

    #exit
    last_time=log_pytms_prog(output,
    'Methionine oxidation complete!')
    return [stored.resi_list, OBJCHIS]

cmd.extend( "oxidize_met", oxidize_met );
################################################################################
################################################################################





################################################################################
# Automated in silico Proline hydroxylation (Model)
################################################################################
def hydroxy_pro(
selection='all',
surface_cutoff=0,
mode=1,
show_clashes=0,
color_mod='',
color_base='',
hydrogens=1,
quiet=1
):


    '''
DESCRIPTION

    Modifies the Prolines of a selection by hydroxylation

EXAMPLE

    frag PRO
    hydroxy_pro

USAGE

    hydroxy_pro [ selection [, surface_cutoff [, mode [, show_clashes
    [, color_mod [, color_base [, hydrogens [, quiet ]]]]]]]]

ARGUMENTS

    selection: selection to be modified {default: 'all'}
               Prolines are automatically sub-selected!

    surface_cutoff: <float> variable for integrated selection of surface atoms.
                    If set >0, PyTMs will automatically calculate the solvent-
                    accessible surface area, and sub-select all atoms above this
                    cutoff.
                    Note that this operation will create a reserved selection:
                    'pytms_surface_input', but only modify residues if the target atoms
                    are also part of the original selection.
                    see also: PyTMs wiki page or findSurfaceResidues (PyMOL wiki)

    mode: <int> toggles the configuration {default=1}
                1: 4R
                2: 4S

    show_clashes: <boolean> toggles if clashes will be visualized {0}
    # consider removing potentially clashing heteroatoms

    color_base, color_mod: color names for
                           object selection and modification, respectively
                           {default: ''} = off

    hydrogens: <int>     toggles if the object will have hydrogens or not
                         after modification
                         0  : as is (detection) {default}
                         1  : adds hydrogens;
                         -1 : no hydrogens;

    quiet: <boolean> toggles output {default: quiet=1}
                     * NB! print appears first in the console

    '''
    ##### BEGINNING OF SCRIPT #####

    last_time=start_time=datetime.datetime.now()
    try:
        surface_cutoff=abs(float(surface_cutoff))
    except:
        raise Exception("Input error!\n Illegal value for surface cutoff!")
    selection=pytms_get_sele(selection, surface_cutoff)

    # argument settings
    try:
        color_base, color_mod = str(color_base), str(color_mod)
        if color_base: cmd.color(color_base,selection)
        if color_mod: cmd.color(color_mod,
        ''+selection+' and (name OD1 or (neighbor (name OD1) and elem H))')
        mode=abs(int(mode))
        if (mode not in [1,2]):
            raise Exception("mode out of range!")
        # transform boolean to reduce checking loops
        show_clashes=bool(int(show_clashes))
        hydrogens=eval_hydrogens(selection, hydrogens)
        output=bool(not int(quiet))

    except:
        raise Exception("Input error!\n Please check the input parameters!")

    # SELECTION LIST ###########################################################
    # create a empty list for appending
    stored.resi_list = []
    stored.temp = []

    # for each object and chain --> create selection lists
    for p in names:
        # for each object and chain get resi of residues
        cmd.iterate((
        '(%s) and ((resn PRO)) '
        'and (name CA)' %(p)),
        'stored.resi_list.append("(%s and resi "+str(resi)+")")'  %p
        )

    # kick out non-selected or missing
    for p in list(stored.resi_list):
        if (cmd.count_atoms('(%s) and (%s) and (name CG)' %(selection,p))==0):
            stored.resi_list.remove(p)

    ##### finished creating selection lists #####

    # premature exit if nothing to process:
    if stored.resi_list==[]:
        print("PyTMs: No modifyable residues found in selection!")
        return False

    # VDW PREP
    # OBJCHI (dictionary) structure
    # obj: [umod. strain, mod. strain.]
    OBJCHIS = {}

    if show_clashes:
        # get unmodified VDW strain
        log_pytms_prog(output, 'Calculating base VdW strain!')
        for p in objects:
            if (hydrogens):
                cmd.h_add(p)
            else:
                cmd.remove('((%s) and elem H)' %p)

            OBJCHIS[p] = [get_strain('(%s)'%(p), temp_clash_A, 0),0]

    last_time=log_pytms_prog(output,
    'Initialized hydroxylation!')

    ##### HYP BUILDER #####

    # cycles through residues
    listcount=-1
    for p in stored.resi_list:
        cmd.unpick()
        listcount=listcount+1

        stored.residue=''
        cmd.iterate('%s and name CA' %p, 'stored.residue=str(resn)')

        ## change name
        cmd.alter(p,'resn="HYP"')

        # build HYP
	# remove potential H atoms
        cmd.remove('%s and ((neighbor (name CG)) and elem H)' %p)

        cmd.edit('%s and name CG' %p)
        cmd.attach('O','2','2')
        cmd.alter('%s and (neighbor (name CG)) and (elem O)' %p,'name="OD1"')
        cmd.attach('O','2','2')
        cmd.alter('%s and (neighbor (name CG)) and (elem O) and not (name OD1)' %p,
        'name="OD2"')
        # now two oxygens are attached, depending on mode, one will be removed
        # imagine looking at the proline from outside the ring and
        # oxygens aligned horizontally and the Nitogen facing up
        dihedOD1=cmd.get_dihedral(
        '%s and (name N)' %p,
        '%s and (name CD)' %p,
        '%s and (name CG)' %p,
        '%s and (name OD1)' %p)

        if mode==1:
            # 4R --> delete left (positive dihedral)
            if dihedOD1>0:
                cmd.remove('%s and (name OD2)' %p)
            else:
                cmd.remove('%s and (name OD1)' %p)
        else:
            # 4S --> delete right (negative dihedral)
            if dihedOD1>0:
                cmd.remove('%s and (name OD1)' %p)
            else:
                cmd.remove('%s and (name OD2)' %p)

        # rename to OS1
        cmd.alter('%s and (neighbor (name CG)) and (elem O)' %p,
        'name="OD1"')

        # fill hydrogens
        cmd.h_add('%s' %p)

        #Log
        last_time=log_pytms_prog(output,
        'Modified: %s'%get_resi_macro_name(p),listcount+1,len(stored.resi_list))
        cmd.unpick()
    # End of resi cycle
    if (hydrogens):
        for p in objects: cmd.h_add(p)
    else:
        for p in objects:
            cmd.remove('((%s) and elem H)' %p)

    # clashes
    if show_clashes:
        # get modified strain
        log_pytms_prog(output, 'Calculating modified VdW strain!')
        if output:
            print('STRAIN REPORT:')
            print('OBJECT','NATIVE_STRAIN','MODIFIED_STRAIN', 'DIFFERENCE')
        for p in objects:
            OBJCHIS[p][1] = get_strain('(%s)'%(p), temp_clash_A, 0)
            cmd.set_name(temp_clash_A,'%s_clashes'%p)
            if output:
                print('%s'%p,OBJCHIS[p][0],OBJCHIS[p][1],OBJCHIS[p][1]-OBJCHIS[p][0])

    # rebuild
    pytms_rebuild()

    # coloring
    if color_base:
        cmd.color(color_base,
        '%s or ((%s extend 1) and hydrogens)' %(selection, selection))
    for p in stored.resi_list:
        # enable selection "p.PTM in hydroxylation_pro"
        if version_ok:
            cmd.alter(('(%s and (resn HYP) and '
            '(name OD1 or ((neighbor (name OD1)) and elem H)))' %p),
            'p.PTM="hydroxylation_pro"')
        if color_base:
            cmd.color(color_base,
            '(%s)' %(p))
        if color_mod:
            cmd.color(color_mod,
            ('(%s and (resn HYP) and '
            '(name OD1 or ((neighbor (name OD1)) and elem H)))' %p))

    #exit
    last_time=log_pytms_prog(output,
    'Proline hydroxylation done!')
    return [stored.resi_list, OBJCHIS]

cmd.extend( "hydroxy_pro", hydroxy_pro );
################################################################################
################################################################################





################################################################################
# Automated in silico methylation (Model)
################################################################################
def methylate(
selection='all',
surface_cutoff=0,
mode=1,
position=0,
show_clashes=0,
color_mod='',
color_base='',
hydrogens=0,
quiet=1
):


    '''
DESCRIPTION

    Modifies the selection by methylation

EXAMPLE

    frag Lys
    methylate

USAGE

    methylate [ selection [, surface_cutoff [, mode [, position
    [, show_clashes [, color_mod [, color_base [, hydrogens [, quiet ]]]]]]]]]

ARGUMENTS

    selection: selection to be modified {default: 'all'}
               Residues are automatically sub-selected!

    surface_cutoff: <float> variable for integrated selection of surface atoms.
                    If set >0, PyTMs will automatically calculate the solvent-
                    accessible surface area, and sub-select all atoms above this
                    cutoff.
                    Note that this operation will create a reserved selection:
                    'pytms_surface_input', but only modify residues if the target atoms
                    are also part of the original selection.
                    see also: PyTMs wiki page or findSurfaceResidues (PyMOL wiki)

    mode: <int> toggles the type of methylation
                0: random (once)
                1: N-(mono)methylation {default}
                2: N-dimethylation
                3: N-trimethylation

    position: <int> toggles the sites to be modified
                0: Lysines only {default}
                1: N-termini only (likely irrelevant)
                2: Lysines and N-termini (likely irrelevant)

    show_clashes: <boolean> toggles if clashes will be visualized {0}
    # consider removing potentially clashing heteroatoms

    color_base, color_mod: color names for
                           object selection and modification, respectively
                           {default: ''} = off

    hydrogens: <int>     toggles if the object will have hydrogens or not
                         after modification
                         0  : as is (detection) {default}
                         1  : adds hydrogens;
                         -1 : no hydrogens;

    quiet: <boolean> toggles output {default: quiet=1}
                     * NB! print appears first in the console

    '''
    ##### BEGINNING OF SCRIPT #####

    last_time=start_time=datetime.datetime.now()
    try:
        surface_cutoff=abs(float(surface_cutoff))
    except:
        raise Exception("Input error!\n Illegal value for surface cutoff!")
    selection=pytms_get_sele(selection, surface_cutoff)

    # argument settings
    try:
        color_base, color_mod = str(color_base), str(color_mod)
        if color_base: cmd.color(color_base,selection)
        if color_mod: cmd.color(color_mod,
            ('(%s and (resn MLZ+MLY+M3L+MLN)) and '
            '(((name CM*) or (name CH*)) or '
            '((neighbor ((name CM*) or (name CH*))) and (elem H)))' %selection))
        mode=abs(int(mode))
        if mode==0: mode=random.choice([1,2,3])
        if (mode not in range(1,4)):
            raise Exception("mode out of range!")
        position=abs(int(position))
        if (position not in range(0,3)):
            raise Exception("position out of range!")
        # transform boolean to reduce checking loops
        show_clashes=bool(int(show_clashes))
        hydrogens=eval_hydrogens(selection, hydrogens)
        output=bool(not int(quiet))

    except:
        raise Exception("Input error!\n Please check the input parameters!")

    # create a empty list for appending
    stored.resi_list = []
    stored.temp_k = []
    stored.temp_n = []
    instance_list = []

    # for each object and chain --> create selection lists
    for p in names:
        if (position!=1):
            # for each object and chain get resi of lysine residues
            cmd.iterate((
            '(%s) and ((resn LYS) or (resn LYN)) '
            'and (name NZ)' %(p)),
            'stored.temp_k.append("(%s and resi "+str(resi)+")")'  %p
            )

        # only append N-term if selected
        if (position!=0):
            # for each object and chain get resi of N-term
            # NB! not hetatm was used instead of polymer to allow 'frag'
            # amino acids to be selected
            cmd.iterate(
            '(byres (first %s)) and (not (hetatm)) '
            'and (name N) and not ((resn PRO) or (resn HYP))' %(p),
            'stored.temp_n.append("(%s and resi "+str(resi)+")")'  %p
            )
            # verify real N-terminus
            for c in list(stored.temp_n):
                if (cmd.count_atoms(('(neighbor (%s and (name N))) '
                'and (elem C)'%c))!=1):
                    stored.temp_n.remove(c)

    # now there are lists for K and N-terms which will be cleaned and merged
    # kick out non-selected or missing (not in selection or already modified)
    for p in list(stored.temp_k):
        if (cmd.count_atoms('(%s) and (%s) and (name NZ)' %(selection,p))==0):
            stored.temp_k.remove(p)
        elif (cmd.count_atoms('(%s) and elem C and (neighbor (name NZ))' %(p))>1):
            # modified otherwise
            stored.temp_k.remove(p)
    for p in list(stored.temp_n):
        if (cmd.count_atoms('(%s) and (%s) and (name N)' %(selection,p))==0):
            stored.temp_n.remove(p)
        elif (cmd.count_atoms(('(neighbor (%s and (name N))) '
            'and (elem C)'%p))!=1):
            # modified otherwise
            stored.temp_n.remove(p)

    # instance will be
    # 4: non n-terminal Lysine (target: NZ)
    # 3: n-terminal Lysine (target: NZ)
    # 2: n-terminal Lysine (target: N)
    # 1: n-term (target: N)
    # 0: skipped (should not occur)

    instance_list = list(stored.temp_k)
    for p in range(0,len(stored.temp_k)):
        # overwrite entry
        instance_list[p]=0
        if (stored.temp_k[p] in stored.temp_n):
            # n-terminal K
            instance_list[p]=3
        else:
            # non n-terminal K
            instance_list[p]=4
        stored.resi_list.append(stored.temp_k[p])

    for p in range(0,len(stored.temp_n)):
        if (stored.temp_n[p] in stored.temp_k):
            # n-terminal K
            instance_list.append(2)
        else:
            # non n-terminal K
            instance_list.append(1)
        stored.resi_list.append(stored.temp_n[p])

    ##### finished creating selection lists #####

    # premature exit if nothing to process:
    if stored.resi_list==[]:
        print("PyTMs: No modifyable residues found in selection!")
        return False

    # VDW PREP
    # OBJCHI (dictionary) structure
    # obj: [umod. strain, mod. strain.]
    OBJCHIS = {}

    if show_clashes:
        # get unmodified VDW strain
        log_pytms_prog(output, 'Calculating base VdW strain!')
        for p in objects:
            if (hydrogens):
                cmd.h_add(p)
            else:
                cmd.remove('((%s) and elem H)' %p)

            OBJCHIS[p] = [get_strain('(%s)'%(p), temp_clash_A, 0),0]

    last_time=log_pytms_prog(output,
    'Initialized methylation!')

    ##### PRE_SETTINGS #####
    # cycles through residues (first lysines then N-terms)
    listcount=-1
    for p in list(stored.resi_list):
        cmd.unpick()
        listcount=listcount+1
        #jump skipped
        if instance_list[listcount]==0:continue

        #pick starting Nitrogen
        if (instance_list[listcount] in [3,4]):
            # handled as lysine
            cmd.remove('(neighbor (%s and name NZ)) and hydrogens' %p)
            cmd.edit('%s and name NZ' %p)
        else:
            # handled as n-term
            cmd.remove('(neighbor (%s and name N)) and hydrogens' %p)
            cmd.edit('%s and (name N)' %p)

        # adjust (assumed) geometry and valence of nitrogen
        geom=4
        val=4
        cmd.set_geometry('pk1', 4, 4)
        cmd.alter('pk1', 'formal_charge=1')

        ##### BUILDER #####

        c_names = [
            [''],
            ['','CM'],
            ['','CH1','CH2'],
            ['','CM1','CM2','CM3']
        ]

        for c in range(1,mode+1):
            cmd.attach('C','4','4')
            cmd.alter(('(neighbor pk1) and (elem C) '
            'and not (name CE or name CA or name CM* or name CH*)'),
            'name="%s"' %(c_names[mode][c]))
        # Hydrogens will be added later using h_add

        # adjust dihedral
        if instance_list[listcount] in  [3,4]:
            cmd.set_dihedral(
                '%s and name CD'%p,
                '%s and name CE'%p,
                '%s and name NZ'%p,
                '%s and (neighbor (name NZ)) and (name %s)'%(p, c_names[mode][1]),
                -149.880)
        else:
            # NB! N-terminal methylation is implemented,
            # but biologically probably not relevant;
            # especially not multiple methyl groups (steric clashes!)
            cmd.set_dihedral(
                '%s and name C'%p,
                '%s and name CA'%p,
                '%s and name N'%p,
                '%s and (neighbor (name N)) and (name %s)'%(p, c_names[mode][1]),
                -60)

        ## change name if LYS
        if instance_list[listcount] in  [3,4]:
            cmd.alter(p,'resn="%s"' %(['','MLZ','MLY','M3L'][mode]))
        else:
            select_part(temp_sel,'%s'%p,'(name CA)','(name CB)')
            # adjust resn, resi and resv
            stored.temp_n=''
            cmd.iterate('(%s) and name N' %(temp_sel),
            'stored.temp_n=int(resv)')
            stored.temp_n=stored.temp_n-1
            cmd.alter(('(%s and (((name CM*) or '
            '(name CH*)) extend 1)) and not (name N)' %temp_sel),
            'resn="MLN"')
            cmd.alter(('(%s and (((name CM*) or '
            '(name CH*)) extend 1)) and not (name N)' %temp_sel),
            'resv=%s' %stored.temp_n)
            cmd.alter(('(%s and (((name CM*) or '
            '(name CH*)) extend 1)) and not (name N)' %temp_sel),
            'resi=%s' %stored.temp_n)
            # replace resi_list entry
            stored.resi_list[listcount]=str(
            "((model "+str(cmd.get_object_list(temp_sel)[0])+
            " and chain '"+str((cmd.get_chains(temp_sel) or [''])[0])+
            "' and resi "+str(stored.temp_n)+") extend 3) "+
            "and not ((name CB) or ((neighbor (name CB)) and (elem H)))"
            )

        #Log
        last_time=log_pytms_prog(output,
        'Modified: %s'%get_resi_macro_name(p),listcount+1,len(stored.resi_list))

    cmd.unpick()
    # end of p in names
    #FINISHED BUILDING methyl group carbons #####

    # fixes ordering of resis and adds hydrogens
    for p in objects:
        cmd.remove('((%s) and elem H)' %p)
        cmd.h_add(p)

    # re-adjust hydrogens
    if (hydrogens):
        for p in objects: cmd.h_add(p)
    else:
        for p in objects:
            cmd.remove('((%s) and elem H)' %p)

    # clashes
    if show_clashes:
        # get modified strain
        log_pytms_prog(output, 'Calculating modified VdW strain!')
        if output:
            print('STRAIN REPORT:')
            print('OBJECT','NATIVE_STRAIN','MODIFIED_STRAIN', 'DIFFERENCE')
        for p in objects:
            OBJCHIS[p][1] = get_strain('(%s)'%(p), temp_clash_A, 0)
            cmd.set_name(temp_clash_A,'%s_clashes'%p)
            if output:
                print('%s'%p,OBJCHIS[p][0],OBJCHIS[p][1],OBJCHIS[p][1]-OBJCHIS[p][0])

    # rebuild
    pytms_rebuild()

    # coloring
    if color_base:
        cmd.color(color_base,
        '%s or ((%s extend 1) and hydrogens)' %(selection, selection))
    for p in stored.resi_list:
        # enable selection "p.PTM in methylation"
        if version_ok:
            cmd.alter(('(%s and (resn MLZ+MLY+M3L+MLN)) and '
            '(((name CM*) or (name CH*)) or '
            '((neighbor ((name CM*) or (name CH*))) and (elem H)))' %p),
            'p.PTM="methylation"')
        if color_base:
            cmd.color(color_base,
            '(%s)' %(p))
        if color_mod:
            cmd.color(color_mod,
            ('(%s and (resn MLZ+MLY+M3L+MLN)) and '
            '(((name CM*) or (name CH*)) or '
            '((neighbor ((name CM*) or (name CH*))) and (elem H)))' %p))

    #exit
    last_time=log_pytms_prog(output, 'Methylation complete!')
    return [stored.resi_list, OBJCHIS]

cmd.extend( "methylate", methylate );
################################################################################
################################################################################




################################################################################
# Automated in silico acetylation (Model)
################################################################################
def acetylate(
selection='all',
surface_cutoff=0,
position=0,
show_clashes=0,
color_mod='',
color_base='',
hydrogens=0,
quiet=1
):


    '''
DESCRIPTION

    Modifies the selection by acetylation

EXAMPLE

    frag LYS
    acetylate

USAGE

    acetylate [ selection [, surface_cutoff [, position [, show_clashes
    [, color_mod [, color_base [, hydrogens [, quiet ]]]]]]]]

ARGUMENTS

    selection: selection to be modified {default: 'all'}
               Residues are automatically sub-selected!

    surface_cutoff: <float> variable for integrated selection of surface atoms.
                    If set >0, PyTMs will automatically calculate the solvent-
                    accessible surface area, and sub-select all atoms above this
                    cutoff.
                    Note that this operation will create a reserved selection:
                    'pytms_surface_input', but only modify residues if the target atoms
                    are also part of the original selection.
                    see also: PyTMs wiki page or findSurfaceResidues (PyMOL wiki)

    position: <int> toggles the sites to be modified
                0: Lysines only {default}
                1: N-termini only
                2: Lysines and N-termini

    show_clashes: <boolean> toggles if clashes will be visualized {0}
    # consider removing potentially clashing heteroatoms

    color_base, color_mod: color names for
                           object selection and modification, respectively
                           {default: ''} = off

    hydrogens: <int>     toggles if the object will have hydrogens or not
                         after modification
                         0  : as is (detection) {default}
                         1  : adds hydrogens;
                         -1 : no hydrogens;

    quiet: <boolean> toggles output {default: quiet=1}
                     * NB! print appears first in the console

    '''
    ##### BEGINNING OF SCRIPT #####

    last_time=start_time=datetime.datetime.now()
    try:
        surface_cutoff=abs(float(surface_cutoff))
    except:
        raise Exception("Input error!\n Illegal value for surface cutoff!")
    selection=pytms_get_sele(selection, surface_cutoff)

    # argument settings
    try:
        color_base, color_mod = str(color_base), str(color_mod)
        if color_base: cmd.color(color_base,selection)
        if color_mod: cmd.color(color_mod,
            ('%s and (resn ACE+ALY) and ((name CH3) extend 2) '
            'and not (elem N)' %selection))
        position=abs(int(position))
        if (position not in range(0,3)):
            raise Exception("position out of range!")
        # transform boolean to reduce checking loops
        show_clashes=bool(int(show_clashes))
        hydrogens=eval_hydrogens(selection, hydrogens)
        output=bool(not int(quiet))
    except:
        raise Exception("Input error!\n Please check the input parameters!")

    # create a empty list for appending
    stored.resi_list = []
    stored.temp_k = []
    stored.temp_n = []
    instance_list = []

    # for each object and chain --> create selection lists
    for p in names:
        if (position!=1):
            # for each object and chain get resi of lysine residues
            cmd.iterate((
            '(%s) and ((resn LYS) or (resn LYN)) '
            'and (name NZ)' %(p)),
            'stored.temp_k.append("(%s and resi "+str(resi)+")")'  %p
            )

        # only append N-term if selected
        if (position!=0):
            # for each object and chain get resi of N-term
            # NB! not hetatm was used instead of polymer to allow 'frag'
            # amino acids to be selected
            cmd.iterate(
            '(byres (first %s)) and (not (hetatm)) '
            'and (name N) and not ((resn PRO) or (resn HYP))' %(p),
            'stored.temp_n.append("(%s and resi "+str(resi)+")")'  %p
            )
            # verify real N-terminus
            for c in list(stored.temp_n):
                if (cmd.count_atoms(('(neighbor (%s and (name N))) '
                'and (elem C)'%c))!=1):
                    stored.temp_n.remove(c)

    # now there are lists for K and N-terms which will be cleaned and merged
    # kick out non-selected or missing (not in selection or already modified)
    for p in list(stored.temp_k):
        if (cmd.count_atoms('(%s) and (%s) and (name NZ)' %(selection,p))==0):
            stored.temp_k.remove(p)
        elif (cmd.count_atoms('(%s) and elem C and (neighbor (name NZ))' %(p))>1):
            # modified otherwise
            stored.temp_k.remove(p)
    for p in list(stored.temp_n):
        if (cmd.count_atoms('(%s) and (%s) and (name N)' %(selection,p))==0):
            stored.temp_n.remove(p)
        elif (cmd.count_atoms(('(neighbor (%s and (name N))) '
            'and (elem C)'%p))!=1):
            # modified otherwise
            stored.temp_n.remove(p)

    # instance will be
    # 4: non n-terminal Lysine (target: NZ)
    # 3: n-terminal Lysine (target: NZ)
    # 2: n-terminal Lysine (target: N)
    # 1: n-term (target: N)
    # 0: skipped (should not occur)

    instance_list = list(stored.temp_k)
    for p in range(0,len(stored.temp_k)):
        # overwrite entry
        instance_list[p]=0
        if (stored.temp_k[p] in stored.temp_n):
            # n-terminal K
            instance_list[p]=3
        else:
            # non n-terminal K
            instance_list[p]=4
        stored.resi_list.append(stored.temp_k[p])

    for p in range(0,len(stored.temp_n)):
        if (stored.temp_n[p] in stored.temp_k):
            # n-terminal K
            instance_list.append(2)
        else:
            # non n-terminal K
            instance_list.append(1)
        stored.resi_list.append(stored.temp_n[p])

    ##### finished creating selection lists #####

    # premature exit if nothing to process:
    if stored.resi_list==[]:
        print("PyTMs: No modifyable residues found in selection!")
        return False

    # VDW PREP
    # OBJCHI (dictionary) structure
    # obj: [umod. strain, mod. strain.]
    OBJCHIS = {}

    if show_clashes:
        # get unmodified VDW strain
        log_pytms_prog(output, 'Calculating base VdW strain!')
        for p in objects:
            if (hydrogens):
                cmd.h_add(p)
            else:
                cmd.remove('((%s) and elem H)' %p)

            OBJCHIS[p] = [get_strain('(%s)'%(p), temp_clash_A, 0),0]

    last_time=log_pytms_prog(output, 'Initialized acetylation!')

    ##### PRE_SETTINGS #####
    # cycles through residues (first lysines then N-terms)
    listcount=-1
    for p in list(stored.resi_list):
        cmd.unpick()
        listcount=listcount+1
        #jump skipped
        if instance_list[listcount]==0:continue

        #pick starting Nitrogen
        if (instance_list[listcount] in [3,4]):
            # handled as lysine
            cmd.remove('(neighbor (%s and name NZ)) and hydrogens' %p)
            cmd.edit('%s and name NZ' %p)
        else:
            # handled as n-term
            cmd.remove('(neighbor (%s and name N)) and hydrogens' %p)
            cmd.edit('%s and (name N)' %p)

        # adjust (assumed) geometry and valence of nitrogen
        cmd.set_geometry('pk1', 3, 3)
        cmd.alter('pk1', 'formal_charge=0')

        ##### BUILDER #####

        cmd.attach('C','3','3')
        cmd.alter(('(neighbor pk1) and (elem C) '
        'and not (name CE or name CA or name C or name CA)'),
        'name="CXXX"')
        # fake DB for sp2
        cmd.unbond('pk1','(%s) and (name CXXX)' %p)
        cmd.bond('pk1','(%s) and (name CXXX)' %p,'2')
        # attach others
        cmd.edit('(neighbor pk1) and (name CXXX)')
        cmd.attach('C','4','4')
        cmd.alter(('(neighbor pk1) and (elem C) '
        'and not (name CE or name CA or name C or name CA)'),
        'name="CH3"')
        cmd.attach('O','2','2')
        cmd.alter('(neighbor pk1) and (elem O) ',
        'name="OXXX"')

        select_part(temp_sel,'%s'%p,'pk1','(name CB)')
        # fix bonds
        cmd.unbond('(%s) and (name CXXX)' %temp_sel,
        '(%s) and (name OXXX)' %temp_sel)
        cmd.bond('(%s) and (name CXXX)' %temp_sel,
        '(%s) and (name OXXX)' %temp_sel,'2')
        cmd.unbond('(%s) and (name CXXX)' %temp_sel,
        '(%s) and ((name N) or (name NZ))' %temp_sel)
        cmd.bond('(%s) and (name CXXX)' %temp_sel,
        '(%s) and ((name N) or (name NZ))' %temp_sel,'1')
        #fix names
        if instance_list[listcount] in  [3,4]:
            cmd.alter('%s and name CXXX' %p, 'name="CH"')
            cmd.alter('%s and name OXXX' %p, 'name="OH"')
        else:
            cmd.alter('%s and name CXXX' %p, 'name="C"')
            cmd.alter('%s and name OXXX' %p, 'name="O"')
        # Hydrogens will be added later using h_add

        ## change name if LYS
        if instance_list[listcount] in  [3,4]:
            cmd.alter(p,'resn="ALY"')
            #default dihedral
            # optimization may be implemented in future version
            cmd.set_dihedral(
            '%s and name CD'%p,
            '%s and name CE'%p,
            '%s and name NZ'%p,
            '%s and (neighbor name NZ) and name CH'%p,
            -173)
            cmd.set_dihedral(
            '%s and name CE'%p,
            '%s and name NZ'%p,
            '%s and (neighbor name NZ) and name CH'%p,
            '((%s and name NZ) extend 2) and name CH3'%p,
            180)
        else:
            # default dihedral (may depend on backbone)
            cmd.set_dihedral(
            '%s and (neighbor name CA) and name C'%p,
            '%s and name CA'%p,
            '%s and name N'%p,
            '%s and (neighbor (name CH3)) and name C'%p,
            -152)
            cmd.set_dihedral(
            '%s and name CA'%p,
            '%s and name N'%p,
            '%s and (neighbor (name CH3)) and name C'%p,
            '((%s and name N) extend 2) and name CH3'%p,
            180)
            # optimization may be implemented in future versions
            select_part(temp_sel,'%s'%p,'(name CA)','(name CB)')
            # adjust resn, resi and resv
            stored.temp_n=''
            cmd.iterate('(%s) and name N' %(temp_sel),
            'stored.temp_n=int(resv)')
            stored.temp_n=stored.temp_n-1
            cmd.alter(('(%s and ((name CH3) extend 2)) '
            'and not (elem N)' %temp_sel),
            'resn="ACE"')
            cmd.alter(('(%s and ((name CH3) extend 2)) '
            'and not (elem N)' %temp_sel),
            'resv=%s' %stored.temp_n)
            cmd.alter(('(%s and ((name CH3) extend 2)) '
            'and not (elem N)' %temp_sel),
            'resi=%s' %stored.temp_n)
            # replace resi_list entry
            stored.resi_list[listcount]=str(
            "((model "+str(cmd.get_object_list(temp_sel)[0])+
            " and chain '"+str((cmd.get_chains(temp_sel) or [''])[0])+
            "' and resi "+str(stored.temp_n)+") extend 3) "+
            "and not ((name CB) or ((neighbor (name CB)) and (elem H)))"
            )

        #Log
        last_time=log_pytms_prog(output,
        'Modified: %s'%get_resi_macro_name(p),listcount+1,len(stored.resi_list))

    cmd.unpick()
    # end of p in stored resi list
    #FINISHED BUILDING ACEs #####

    # fixes ordering of resis and adds hydrogens
    for p in objects:
        cmd.remove('((%s) and elem H)' %p)
        cmd.h_add(p)

    # re-adjust hydrogens
    if (hydrogens):
        for p in objects: cmd.h_add(p)
    else:
        for p in objects:
            cmd.remove('((%s) and elem H)' %p)

    # clashes
    if show_clashes:
        # get modified strain
        log_pytms_prog(output, 'Calculating modified VdW strain!')
        if output:
            print('STRAIN REPORT:')
            print('OBJECT','NATIVE_STRAIN','MODIFIED_STRAIN', 'DIFFERENCE')
        for p in objects:
            OBJCHIS[p][1] = get_strain('(%s)'%(p), temp_clash_A, 0)
            cmd.set_name(temp_clash_A,'%s_clashes'%p)
            if output:
                print('%s'%p,OBJCHIS[p][0],OBJCHIS[p][1],OBJCHIS[p][1]-OBJCHIS[p][0])

    # rebuild
    pytms_rebuild()

    # coloring
    if color_base:
        cmd.color(color_base,
        '%s or ((%s extend 1) and hydrogens)' %(selection, selection))
    for p in stored.resi_list:
        # enable selection "p.PTM in acetylation"
        if version_ok:
            cmd.alter(('(%s and ((name CH3) extend 2)) '
            'and not (elem N)' %p),
            'p.PTM="acetylation"')
        if color_base:
            cmd.color(color_base,
            '(%s)' %(p))
        if color_mod: cmd.color(color_mod,
            ('%s and ((name CH3) extend 2) '
            'and not (elem N)' %p))

    #exit
    last_time=log_pytms_prog(output, 'Acetylation complete!')
    return [stored.resi_list, OBJCHIS]

cmd.extend( "acetylate", acetylate );
################################################################################
################################################################################

