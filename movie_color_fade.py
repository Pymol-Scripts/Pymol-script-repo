'''
movie_color_fade.py
Described at: http://www.pymolwiki.org/movie_color_fade
Version 1.0 (2014)
##################################################

    Fades the color of representations in movies
    #NB!: Defines and uses new color names using the selection name and frame numbers

##################################################
Plugin contributed by Andreas Warnecke
(andreas.warnecke@ki.se, 4ndreas.warneck3@gmail.com)
##################################################
VERSION NOTES:
    1.0    2014    First release
'''
#-------------------------------------------------------------------------------
# IMPORT
from __future__ import print_function
from pymol import cmd
from pymol import stored
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# MOVIE COLOR FADING
#-------------------------------------------------------------------------------
def movie_color_fade(
    startframe='',
    startcolor='red',
    endframe='',
    endcolor='green',
    selection='all'
):
    '''
DESCRIPTION

    Fades the color of representations in movies
    #NB!: Defines and uses new color names using the selection name and frame numbers

USE

    movie_color_fade startframe='', startcolor=red, endframe='', endcolor=green, selection=all
    #defaults indicated

PARAMETERS

   startframe, endframe = beginning and end movie frame for fading
   startcolor, endcolor = coloring at start and end
   selection: target selection

   NB! startframe and endframe can be omitted or set='' to assign current and last frame respectively

EXAMPLES

    ##### 1. #####
    fetch 1hpv, async=0
    as cartoon
    orient
    color yellow
    mset 1x120
    movie_color_fade 1, yellow, 60, blue
    movie_color_fade 60, blue, 120, yellow
    #####

    ##### 2. #####
    #repeat command and specify 'selection' to change multiple colors
    fetch 1hpv, async=0
    as cartoon
    orient
    color white
    mset 1x60
    movie_color_fade auto,white,auto,skyblue,ss s
    movie_color_fade auto,white,auto,red,ss h
    movie_color_fade auto,white,auto,grey,ss l+""
    #####

SEE ALSO

    mdo, mappend, set, movie_fade
    '''
    selection = '(' + str(selection) + ')'

    try:
        startframe = int(startframe)
    except:
        startframe = int(cmd.get('frame'))
    try:
        endframe = int(endframe)
    except:
        endframe = int(cmd.count_frames())

    if endframe == 0:
        endframe = 1

    if startframe == endframe:
        print("start == end")
        return False

    if startframe > endframe:
        startframe, endframe = endframe, startframe
        startcolor, endcolor = endcolor, startcolor

    # color RGBs
    try:
        startcolor = cmd.get_color_tuple(startcolor)
        endcolor = cmd.get_color_tuple(endcolor)
        diffcolor = [b - a for a, b in zip(startcolor, endcolor)]
    except:
        print("Input error - please provide regular colors")
        return False

    for frame in range(startframe, endframe + 1):
        # calculate intermediates
        frac = float(frame - startframe) / (endframe - startframe)
        endcolor = [a * frac for a in diffcolor]
        endcolor = list(map(sum, list(zip(startcolor, endcolor))))
        colorname = selection + "_" + str(frame)
        # define new color
        cmd.set_color(colorname, endcolor)

        cmd.mappend(frame, "/cmd.color(%s, %s)" % (repr(colorname), repr(selection)))

cmd.extend("movie_color_fade", movie_color_fade)
cmd.auto_arg[0]['movie_color_fade'] = [lambda: cmd.Shortcut(['auto', '1']), 'startframe=', ',']
cmd.auto_arg[1]['movie_color_fade'] = [lambda: cmd.Shortcut(['red', 'green', 'blue', 'yellow']), 'startcolor=', ',']
cmd.auto_arg[2]['movie_color_fade'] = [lambda: cmd.Shortcut(['auto']), 'endframe=', ',']
cmd.auto_arg[3]['movie_color_fade'] = [lambda: cmd.Shortcut(['red', 'green', 'blue', 'yellow']), 'endcolor=  ... more see "help movie_color_fade"', ',']
#cmd.auto_arg[4]['movie_color_fade']=[lambda: cmd.Shortcut(['all']), 'selection=', '']
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
