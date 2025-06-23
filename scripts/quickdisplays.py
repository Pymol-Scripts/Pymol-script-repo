'''
quickdisplays.py
Described at: http://www.pymolwiki.org/quickdisplays
##################################################

Bundled commands for quick default representations.
Type disp_list for listing of functions.

##################################################
Plugin contributed by Andreas Warnecke
(andreas.warnecke@ki.se, 4ndreas.warneck3@gmail.com)
##################################################
VERSION NOTES:
    1.0    2014     First release
    1.1    2016     disabled disp_surf from setting surface_carve_selection
                    corrected print output
'''
#-------------------------------------------------------------------------------
# IMPORT
from __future__ import print_function
from pymol import cmd
from pymol import stored
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# PRINT LIST OF FUNCTIONS
#-------------------------------------------------------------------------------
def disp_list():
    print(
        '''
quickdisplay - list of functions.
Enter, e.g. "help disp_ss" for more specific info.

    PURPOSE                          FUNCTION
    secondary structure cartoon......disp_ss
    ball and stick representation....disp_ball_stick
    mesh display.....................disp_mesh
    surface display..................disp_surf
    Putty b-factor sausage...........disp_putty
'''
    )
cmd.extend("disp_list", disp_list)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# PERCENTILE LIMITS (used below)
#-------------------------------------------------------------------------------
def get_b_limits(input='[0,100]', selection='all'):
    from pymol import cmd, stored
    import math
    try:
        limits = str(input)
        if limits.startswith('['):
            # list
            limits = eval(limits)[:2]
        else:
            # integer
            selection = str(selection)
            limits = abs(float(limits))
            if limits > 50:
                return False
            stored.temp = []
            cmd.iterate(selection, 'stored.temp.append(b)')
            stored.temp.sort()
            kl = (len(stored.temp) - 1) * limits / 100
            fl = math.floor(kl)
            cl = math.ceil(kl)
            kr = (len(stored.temp) - 1) * (100 - limits) / 100
            fr = math.floor(kr)
            cr = math.ceil(kr)
            limits = [0, 0]
            if fl == cl:
                limits[0] = stored.temp[int(kl)]
            else:
                limits[0] = stored.temp[int(kl)] * (cl - kl) + stored.temp[int(kl)] * (kl - fl)
            if fr == cr:
                limits[1] = stored.temp[int(kr)]
            else:
                limits[1] = stored.temp[int(kr)] * (cr - kr) + stored.temp[int(kr)] * (kr - fr)
    except:
        return False
    return limits
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# DISPLAY: SECONDARY STRUCTURE
#-------------------------------------------------------------------------------
def disp_ss(
        selection='all',
        colors='marine red white',
        only=False):
    '''
DESCRIPTION

    Formats the passed object into secondary structure cartoon

USAGE

    disp_ss [ selection [, colors [, only ]]]

PARAMETERS

    NAME=DEFAULT               TYPE    FUNCTION
    selection='all'            <str>   input selection
    colors='marine red white'  <str>   any three colors for: sheets, helices and loops
                                       e.g. 'marine red white'
                                       can also be set to either util.cbc, util.rainbow,
                                       or util.chainbow (alone)
                                       set='False' to supress coloring altogether, or enter False
                                       for the coloring to be omitted, e.g. 'marine False green'
    only                       <bool>  if True will use show_as; else show
    '''

    try:
        selection = '(' + selection + ')'
        colors = str(colors)
        only = bool(str(only) != 'False')
    except:
        print("Input error")
        return False

    if colors == 'False':
        color_s = color_h = color_l = False
    elif colors.startswith('util.'):
        util_choices = ['util.cbc', 'util.rainbow', 'util.chainbow', 'util.ss']
        if colors in util_choices:
            color_s = color_h = color_l = '%s ' % util_choices[util_choices.index(colors)]
        else:
            print("Input error! Please check the color setting using util. is one of:", util_choices)
            return False
    else:
        try:
            color_s, color_h, color_l = colors.split()[:3]
            if color_s != 'False':
                cmd.color(color_s, None)
                color_s = 'color %s, ' % color_s
            else:
                color_s = False
            if color_h != 'False':
                cmd.color(color_h, None)
                color_h = 'color %s, ' % color_h
            else:
                color_h = False
            if color_l != 'False':
                cmd.color(color_l, None)
                color_l = 'color %s, ' % color_l
            else:
                color_l = False
        except:
            print("Input error! Please check that three valid colors (or False) are provided")
            return False

    for p in cmd.get_object_list(selection):
        cmd.set('cartoon_discrete_colors', 'on', p)
        # settings
        cmd.cartoon('rectangle', 'ss s and %s and %s' % (selection, p))
        cmd.cartoon('dumbbell', 'ss h and %s and %s' % (selection, p))
        cmd.cartoon('loop', 'ss l+"" and %s and %s' % (selection, p))
        # sheets
        if color_s:
            print(cmd.do(color_s + '(ss s and %s and %s)' % (selection, p)))
        cmd.set('cartoon_rect_length', 1.5, p)
        cmd.set('cartoon_rect_width', 0.25, p)
        # a-helices
        if color_h:
            print(cmd.do(color_h + '(ss h and %s and %s)' % (selection, p)))
        cmd.set('cartoon_dumbbell_length', 1.5, p)
        cmd.set('cartoon_dumbbell_width', 0.25, p)
        cmd.set('cartoon_dumbbell_radius', 0.2, p)
        # loops
        if color_l:
            print(cmd.do(color_l + '(ss l+"" and %s and %s)' % (selection, p)))
        cmd.set('cartoon_loop_radius', 0.25, p)

        if only:
            cmd.show_as('cartoon', '%s and %s' % (selection, p))
        else:
            cmd.show('cartoon', '%s and %s' % (selection, p))

cmd.extend("disp_ss", disp_ss)
cmd.auto_arg[0]['disp_ss'] = [lambda: cmd.Shortcut(['all']), 'selection=', ',']
cmd.auto_arg[1]['disp_ss'] = [lambda: cmd.Shortcut(['default', 'blue', 'yellow']), 'color_s=', ',']
cmd.auto_arg[2]['disp_ss'] = [lambda: cmd.Shortcut(['default', 'red', 'blue']), 'color_h=', ',']
cmd.auto_arg[3]['disp_ss'] = [lambda: cmd.Shortcut(['color_l=default, only=False']), 'remaining (defaults)...', '']
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# DISPLAY: BALL AND STICK
#-------------------------------------------------------------------------------
def disp_ball_stick(selection='all', hydrogens=0, only=False):
    '''
DESCRIPTION

    Formats the passed object into ball and stick

USEAGE

    disp_ball_stick [ selection [, hydrogens [, only ]]]

EXAMPLE

    fetch 1hpv, async=0
    disp_ball_stick
    util.cbaw

PARAMETERS

    NAME=DEFAULT       TYPE    FUNCTION
    selection='all'    <str>   input selection
    hydrogens          <int>   -1: remove; 1: add; else: as is
    only=False         <bool>  if True will use show_as; else show

    '''
    try:
        selection = '(' + selection + ')'
        hydrogens = int(hydrogens)
        only = bool(str(only) != 'False')
    except:
        print("Input error")
        return False

    if hydrogens == 1:
        cmd.h_add('%s' % selection)
    if hydrogens == -1:
        cmd.remove('%s and elem H' % selection)

    for p in cmd.get_object_list(selection):
        cmd.set('valence', 'on', p)
        cmd.set('stick_ball', 'on', p)
        cmd.set('stick_ball_ratio', 3, p)
        cmd.set('stick_radius', 0.12, p)

    if only:
        cmd.show_as('sticks', '%s' % (selection))
    else:
        cmd.show('sticks', '%s' % (selection))


cmd.extend("disp_ball_stick", disp_ball_stick)
cmd.auto_arg[0]['disp_ball_stick'] = [lambda: cmd.Shortcut(['all']), 'selection=', ',']
cmd.auto_arg[1]['disp_ball_stick'] = [lambda: cmd.Shortcut(['-1', '0', '1']), 'hydrogens=', ',']
cmd.auto_arg[2]['disp_ball_stick'] = [lambda: cmd.Shortcut(['True', 'False']), 'only=', '']
#-------------------------------------------------------------------------------


def disp_stick_ball(selection='all', hydrogens=0, only=False):
    '''
    see help disp_stick_ball
    '''
    #disp_stick_ball - redirected
    disp_ball_stick(selection, hydrogens, only)
cmd.extend("disp_stick_ball", disp_stick_ball)
cmd.auto_arg[0]['disp_stick_ball'] = [lambda: cmd.Shortcut(['all']), 'selection=', ',']
cmd.auto_arg[1]['disp_stick_ball'] = [lambda: cmd.Shortcut(['-1', '0', '1']), 'hydrogens=', ',']
cmd.auto_arg[2]['disp_stick_ball'] = [lambda: cmd.Shortcut(['True', 'False']), 'only=', '']
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# DISPLAY: MESH
#-------------------------------------------------------------------------------
def disp_mesh(selection='all', color_m='default', hydrogens=0, only=False, limits=5):
    '''
DESCRIPTION

    Adds a mesh to the object
    Has advanced coloring options and automatically accounts for the hydrogens

USEAGE

    disp_mesh [ selection [, color_m [, hydrogens [, only [, limits]]]]]
    disp_mesh selection=all, color_m=default
    disp_mesh selection=all, color_m=white
    disp_mesh selection=all, color_m=putty

PARAMETERS

    NAME=DEFAULT       TYPE    FUNCTION
    selection='all'    <str>   input selection
    color_m='default'  <str>   'default': as current
                               'name': colors by color or ramp called name
                               'putty': b-factor on surface
    hydrogens=0        <int>   -1: remove; 1: add; else: as is
    only=False         <bool>  if True will use show_as; else show
    limits=5           <list or flaot>
                               applies only if color_m=='putty'
                               sets the b-factor range limits
                               <list> [min,max] # absolute values
                               <float> percentile cutoff (both sides) # relative for each protein
    '''

    try:
        selection = '(' + selection + ')'
        color_m = str(color_m)
        hydrogens = int(hydrogens)
        only = bool(str(only) != 'False')
    except:
        print("Input error")
        return False

    if hydrogens == 1:
        cmd.h_add('%s' % selection)
    if hydrogens == -1:
        cmd.remove('%s and elem H' % selection)

    for p in cmd.get_object_list(selection):
        cmd.set('mesh_width', 0.25, p)
        if (color_m == 'putty'):
            limits = get_b_limits(limits, p)
            if not limits:
                print("Input error (limits must be <list> or <float (<=50)>)!")
                return False

            cmd.set('mesh_color', 'default', p)
            cmd.spectrum('b', 'rainbow', '(not hetatm) and %s' % p, minimum='%f' % limits[0], maximum='%f' % limits[1], byres=0)
            print("disp_mesh:", p, "displayed in putty mode - mesh as putty - limits=[%.4f,%.4f]" % (limits[0], limits[1]))
        else:
            cmd.set('mesh_color', color_m, p)
            print("disp_mesh: regular mode - mesh - " + p)

        if only:
            cmd.show_as('mesh', '%s and %s' % (selection, p))
        else:
            cmd.show('mesh', '%s and %s' % (selection, p))
        cmd.rebuild()

cmd.extend("disp_mesh", disp_mesh)
cmd.auto_arg[0]['disp_mesh'] = [lambda: cmd.Shortcut(['all']), 'selection=', ',']
cmd.auto_arg[1]['disp_mesh'] = [lambda: cmd.Shortcut(['default', 'red', 'green', 'blue', 'yellow']), 'color_m=', ',']
cmd.auto_arg[2]['disp_mesh'] = [lambda: cmd.Shortcut(['-1', '0', '1']), 'hydrogens=', ',']
cmd.auto_arg[3]['disp_mesh'] = [lambda: cmd.Shortcut(['only=False, limits=5']), 'remaining (defaults)...', '']
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# DISPLAY: SURFACE AND SURFACE AS PUTTY
#-------------------------------------------------------------------------------
def disp_surf(
        selection='all',
        color_s='default',
        transparency=0,
        hydrogens=0,
        solvent=0,
        ramp_above=1,
        only=False,
        limits=5):
    '''
DESCRIPTION

    Advanced surface representation (cf. examples)

USAGE

    disp_surf [ selection [, color_s [, transparency [, hydrogens [, solvent [, ramp_above [, only [, limits]]]]]]]]

EXAMPLES

    disp_surf # opaque surface with default colors
    disp_surf all, white, 0.5 # half-transparent white surface
    disp_surf all, putty # b-factor on surface

PARAMETERS

    NAME=DEFAULT       TYPE    FUNCTION
    selection='all'    <str>   input selection
    color_s='default'  <str>   'default': as current
                               'name': colors by color or ramp called name
                               'putty': b-factor on surface (by resi)
    transparency=0     <float> set surface transparency
    hydrogens=0        <int>   -1: remove; 1: add; else: as is
    solvent=0          <int>   defines 'surface_solvent'
    ramp_above=1       <int>   defines 'surface_ramp_above_mode'
    only=False         <bool>  if True will use show_as; else show
    limits=5           <list or flaot>
                               applies only if color_s=='putty'
                               sets the b-factor range limits
                               <list> [min,max] # absolute values
                               <float> percentile cutoff (both sides) # relative for each protein
    '''
    try:
        selection = '(' + selection + ')'
        color_s = str(color_s)
        transparency = float(transparency)
        hydrogens = int(hydrogens)
        solvent = int(solvent)
        ramp_above = int(ramp_above)
        only = bool(str(only) != 'False')
    except:
        print("Input error")
        return False

    for p in cmd.get_object_list(selection):
        if hydrogens == -1:
            cmd.remove('%s and elem H' % p)
        if hydrogens == 1:
            cmd.h_add(p)
        # if hydrogens==0: as is

        # some defaults (settings can be changed later, too)
        cmd.set('surface_carve_cutoff', '4.5', p)
        cmd.set('surface_quality', '2', p)
        cmd.set('solvent_radius', '1.5', p)
        cmd.set('cavity_cull', '2', p)
        #cmd.set('surface_carve_selection', p)
        # this causes issues when using the command on different objects
        # it is better to define surface_carve_selection manually

        # defined
        cmd.set('surface_solvent', solvent, p)
        cmd.set('surface_ramp_above_mode', ramp_above, p)
        cmd.set('transparency', transparency, p)

        if (color_s == 'putty'):
            limits = get_b_limits(limits, p)
            if not limits:
                print("Input error (limits must be <list> or <float (<=50)>)!")
                return False

            cmd.set('surface_color', 'default', p)
            cmd.spectrum('b', 'rainbow',
                         '(not hetatm) and %s' % p, minimum='%f' % limits[0],
                         maximum='%f' % limits[1], byres=0)
            print("disp_surf:", p, "displayed in putty mode - surface as putty - limits=[%.4f,%.4f]" % (limits[0], limits[1]))
        else:
            cmd.set('surface_color', color_s, selection)
            print("disp_surf:", p, "displayed as regular surface")
        if only:
            cmd.show_as('surface', '(%s and %s)' % (selection, p))
        else:
            cmd.show('surface', '(%s and %s)' % (selection, p))
cmd.extend("disp_surf", disp_surf)
cmd.auto_arg[0]['disp_surf'] = [lambda: cmd.Shortcut(['all']), 'selection=', ',']
cmd.auto_arg[1]['disp_surf'] = [lambda: cmd.Shortcut(['default', 'red', 'green', 'blue', 'yellow']), 'color_s=', ',']
cmd.auto_arg[2]['disp_surf'] = [lambda: cmd.Shortcut(['0', '0.10', '0.20', '0.30', '0.40', '0.50', '0.60', '0.70', '0.80', '0.90', '1.00']), 'transparency=', ',']
cmd.auto_arg[3]['disp_surf'] = [lambda: cmd.Shortcut(['hydrogens=0, solvent=0, ramp_above=1, only=False, limits=5']), 'remaining (deafults)...', ',']
#cmd.auto_arg[4]['disp_surf']=[lambda: cmd.Shortcut(['0']), 'solvent=', ',']
#cmd.auto_arg[5]['disp_surf']=[lambda: cmd.Shortcut(['1']), 'ramp_above=', ',']
#cmd.auto_arg[6]['disp_surf']=[lambda: cmd.Shortcut(['True','False']), 'only=', '']
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# DISPLAY: PUTTY
#-------------------------------------------------------------------------------
def disp_putty(selection='all', limits=10, only=True):
    '''
DESCRIPTION

    Formats the passed object into a Putty b-factor sausage

USEAGE

    disp_putty [ selection ]
    selection    <str>    input selection
    limits=10    <list or flaot>
                          applies only if color_m=='putty'
                          sets the b-factor range limits (by protein)
                          <list> [min,max]
                          <float> percentile cutoff (both sides)
    only=True             <bool>  if True will use show_as; else show
    '''

    try:
        selection = '(' + selection + ')'
        only = bool(str(only) != 'False')
    except:
        print("Input error")
        return False

    for p in cmd.get_object_list(selection):
        limits = get_b_limits(limits, p)
        if not limits:
            print("Input error (limits must be <list> or <float (<=50)>)!")
            return False
        # settings
        cmd.set('cartoon_discrete_colors', 'off', p)
        cmd.set('cartoon_putty_scale_max', 5, p)
        cmd.set('cartoon_putty_scale_min', 1, p)
        # normalized nonlinear scaling
        cmd.set('cartoon_putty_transform', 0, p)
        cmd.spectrum('b', 'rainbow', '(not hetatm) and %s' % p, minimum='%f' % limits[0], maximum='%f' % limits[1], byres=1)
        cmd.cartoon('putty', '(%s and %s)' % (selection, p))
        if only:
            cmd.show_as('cartoon', '(%s and %s)' % (selection, p))
        else:
            cmd.show('cartoon', '(%s and %s)' % (selection, p))
cmd.extend("disp_putty", disp_putty)
cmd.auto_arg[0]['disp_putty'] = [lambda: cmd.Shortcut(['all', 'all and visible']), 'selection=', ',']
cmd.auto_arg[1]['disp_putty'] = [lambda: cmd.Shortcut('0', '5', '10', '[10,50]'), 'limits=', ',']
cmd.auto_arg[2]['disp_putty'] = [lambda: cmd.Shortcut(['True', 'False']), 'only=', '']
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
