'''
get_colors.py
Described at: http://www.pymolwiki.org/get_colors
Version 1.0 (2014)
##################################################

Includes two functions:
    get_colors: returns a list of defined pymol colors
    get_random_color: returns a random color

##################################################
Plugin contributed by Andreas Warnecke
(andreas.warnecke@ki.se, 4ndreas.warneck3@gmail.com)
##################################################
VERSION NOTES:
    1.0    2014    First release
'''
#-------------------------------------------------------------------------------
from pymol import cmd

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


def get_colors(selection='', output=True):
    '''
DESCRIPTION:
    returns a list of available pymol colors

USAGE:
    get_colors [ selection [, output ]]

EXAMPLES:
    get_colors # basic colors
    get colors all # larger range with intermediates
    '''
    # there may be a more elegant way of doing this?
    import pymol
    pymol_color_list=[]
    for tuplepair in pymol.querying.get_color_indices(selection): 
        pymol_color_list.append(tuplepair[0])
    pymol_color_list.sort()
    # the str() is due to pymol passing aruments as strings, 
    # thus True becomes 'True'
    if str(output)=='True': print pymol_color_list
    return pymol_color_list
cmd.extend('get_colors',get_colors)
cmd.auto_arg[0]['get_colors']=[lambda: cmd.Shortcut(['""','all']), 'selection=', ',']
cmd.auto_arg[1]['get_colors']=[lambda: cmd.Shortcut(['True']), 'output=', '']
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def get_random_color(selection='', output=False):
    '''
DESCRIPTION:
    returns a random color name available in pymol
    ! Requires get_colors !Indended mostly for use in Python

USAGE:
    get_random_color [ selection [, output ]]

EXAMPLES:
    # print a random color name: 
    get_random_color output=True
    # color object randomly:
    fetch 1hpv, async=0
    cmd.color(get_random_color())
    '''
    import random
    randomcolor=random.choice(get_colors(selection, False))
    if str(output)!='False': print randomcolor
    return randomcolor
cmd.extend('get_random_color',get_random_color)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
