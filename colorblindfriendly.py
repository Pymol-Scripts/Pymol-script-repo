'''
More information and examples can be found at:
http://www.pymolwiki.org/index.php/colorblindfriendly

DESCRIPTION

    Certain colors are indistinguishable to people with the various forms of
    color blindness, and therefore are better not used in figures intended for
    public viewing.

    This script generates a palette of named colors for PyMOL that are
    unambiguous both to colorblind and non-colorblind people.

    The colors listed here are defined according to recommendations found at
    http://jfly.iam.u-tokyo.ac.jp/html/color_blind/.  This website is a good
    reference to consult when making all kinds of figures, not just those made
    using PyMOL.

    The colors are:

    * cb_black
    * cb_orange
    * cb_sky_blue (also: cb_skyblue, cb_light_blue, cb_lightblue)
    * cb_bluish_green (also: cb_bluishgreen, cb_green)
    * cb_yellow
    * cb_blue
    * cb_vermillion (also: cb_red, cb_redorange, cb_red_orange)
    * cb_reddish_purple (also: cb_rose, cb_violet, cb_magenta)

USAGE

    import colorblindfriendly as cbf

    # Add the new colors
    cbf.set_colors()
    color myObject, cb_red

    # Replace built-in colors with cbf ones
    cbf.set_colors(replace=True)
    color myOtherObject, yellow   # actually cb_yellow

    # Add a `cb_colors` menu item to the OpenGL GUI ([C] menu in the right panel)
    cbf.add_menu()

REQUIREMENTS

    The cb_colors menu (`add_menu()` function) requires PyMOL 1.6.0 or later.

AUTHOR

    Jared Sampson
    Github: @jaredsampson

LICENSE

Copyright (c) 2014-2021 Jared Sampson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

CHANGELOG

    0.3.0   Generalize the way colors and menus are efined and added, to
            enable the use of additional color palettes. Add Viridis and Magma
            palettes (contributed by Yehudi Bloch).  [2021-10-27]

    0.2.0   Complete overhaul for PyMOL 2.0 with conversion to module format.
            Now, setting the new `cb_*` color values requires a call to the
           `set_colors()` function after import.  You can also now add a
           `cb_colors` menu to the OpenGL GUI via the `add_menu()` function.
            [10/26/2017]

'''
from __future__ import print_function
import math

__author__ = 'Jared Sampson'
__version__ = '0.3.0'

import pymol
from pymol import cmd


# Color blind-friendly color list based on information found at:
# http://jfly.iam.u-tokyo.ac.jp/html/color_blind/#pallet
CB_COLORS = [
    {
        'name': 'red',
        'rgb': [213, 94, 0],
        'alt': ['vermillion', 'red_orange', 'redorange'],
    },
    {
        'name': 'orange',
        'rgb': [230, 159, 0],
        'alt': None,
    },
    {
        'name': 'yellow',
        'rgb': [240, 228, 66],
        'alt': None,
    },
    {
        'name': 'green',
        'rgb': [0, 158, 115],
        'alt': ['bluish_green', 'bluishgreen'],
    },
    {
        'name': 'light_blue',
        'rgb': [86, 180, 233],
        'alt': ['sky_blue', 'skyblue', 'lightblue'],
    },
    {
        'name': 'blue',
        'rgb': [0, 114, 178],
        'alt': None,
    },
    {
        'name': 'violet',
        'rgb': [204, 121, 167],
        'alt': ['reddish_purple', 'reddishpurple', 'rose', 'magenta'],
    },
    {
        'name': 'black',
        'rgb': [0, 0, 0],
        'code': '222',
        'alt': None,
    },
]

# Viridis and Magma palettes contributed by Yehudi Bloch, originally
# developed by Stéfan van der Walt and Nathaniel Smith for matplotlib.
# https://matplotlib.org/stable/users/prev_whats_new/whats_new_1.5.html
VIRIDIS_COLORS = [
    {'name':  'viridis1', 'rgb': [253, 231,  36], 'alt': None},    # noqa: E241
    {'name':  'viridis2', 'rgb': [186, 222,  39], 'alt': None},    # noqa: E241
    {'name':  'viridis3', 'rgb': [121, 209,  81], 'alt': None},    # noqa: E241
    {'name':  'viridis4', 'rgb': [ 66, 190, 113], 'alt': None},    # noqa: E241
    {'name':  'viridis5', 'rgb': [ 34, 167, 132], 'alt': None},    # noqa: E241
    {'name':  'viridis6', 'rgb': [ 32, 143, 140], 'alt': None},    # noqa: E241
    {'name':  'viridis7', 'rgb': [ 41, 120, 142], 'alt': None},    # noqa: E241
    {'name':  'viridis8', 'rgb': [ 52,  94, 141], 'alt': None},    # noqa: E241
    {'name':  'viridis9', 'rgb': [ 64,  67, 135], 'alt': None},    # noqa: E241
    {'name': 'viridis10', 'rgb': [ 72,  35, 116], 'alt': None},    # noqa: E241
    {'name': 'viridis11', 'rgb': [ 68,   1,  84], 'alt': None},    # noqa: E241
]

MAGMA_COLORS = [
    {'name':  'magma1', 'rgb': [251, 252, 191], 'alt': None},      # noqa: E241
    {'name':  'magma2', 'rgb': [253, 205, 114], 'alt': None},      # noqa: E241
    {'name':  'magma3', 'rgb': [253, 159, 108], 'alt': None},      # noqa: E241
    {'name':  'magma4', 'rgb': [246, 110,  91], 'alt': None},      # noqa: E241
    {'name':  'magma5', 'rgb': [221,  73, 104], 'alt': None},      # noqa: E241
    {'name':  'magma6', 'rgb': [181,  54, 121], 'alt': None},      # noqa: E241
    {'name':  'magma7', 'rgb': [140,  41, 128], 'alt': None},      # noqa: E241
    {'name':  'magma8', 'rgb': [ 99,  25, 127], 'alt': None},      # noqa: E241
    {'name':  'magma9', 'rgb': [ 59,  15, 111], 'alt': None},      # noqa: E241
    {'name': 'magma10', 'rgb': [ 20,  13,  53], 'alt': None},      # noqa: E241
    {'name': 'magma11', 'rgb': [  0,   0,   3], 'alt': None},      # noqa: E241
]

PALETTES = {
    'cb_colors': {
        'colors': CB_COLORS,
        'prefix': 'cb_',
    },
    'viridis': {
        'colors': VIRIDIS_COLORS,
        'prefix': '',
    },
    'magma': {
        'colors': MAGMA_COLORS,
        'prefix': '',
    },
}


def _get_palettes(palette=None):
    '''Return the desired palettes dict.'''
    if palette is None:
        palettes = PALETTES
    else:
        palettes = {palette: PALETTES[palette]}
    return palettes


def set_colors(palette=None, replace=False):
    '''Add the color blind-friendly colors to PyMOL.'''
    palettes = _get_palettes(palette)
    for pname, p in palettes.items():
        added_colors = []
        for c in p['colors']:
            # RGB tuple shortcut
            rgb = c['rgb']

            # Get the primary and alternate color names into a single list
            names = [c['name']]
            if c['alt']:
                names.extend(c['alt'])

            # Set the colors
            for name in names:
                try:
                    use_name = p['prefix'] + name
                except KeyError:
                    use_name = name
                cmd.set_color(use_name, rgb)

                # Optionally replace built-in colors
                if replace:
                    cmd.set_color(name, rgb)
                    spacer = (20 - len(name)) * ' '
                    added_colors.append('    {}{}{}'.format(name, spacer, use_name))
                else:
                    added_colors.append('    {}'.format(use_name))

        # Notify user of newly available colors
        print(f'These {pname} colors are now available:')
        print('\n'.join(added_colors))


def _add_palette_menu(name, palette, replace=False):
    '''Add a color palette to the PyMOL OpenGL menu.'''

    # Make sure cb_colors are installed.
    print(f'Checking for {name} colors...')
    try:
        for c in palette['colors']:
            if cmd.get_color_index(c['name']) == -1:
                # mimic pre-1.7.4 behavior
                raise pymol.CmdException
    except pymol.CmdException:
        print(f'Adding {name} palette colors...')
        set_colors(palette=name)

    # Abort if PyMOL is too old.
    try:
        from pymol.menu import all_colors_list
    except ImportError:
        print('PyMOL version too old for cb_colors menu. Requires 1.6.0 or later.')
        return

    # Add the menu
    print(f'Adding {name} menu...')
    # mimic pymol.menu.all_colors_list format
    # first color in list is used for menu item color

    def _get_color_code(rgb):
        '''Return a 3-digit string approximating the RGB color.'''
        return ''.join([str(math.floor(x / 256 * 10)) for x in rgb])

    # Menu item for each color should be a tuple in the form
    #    ('999', 'color_name')
    # where '999' is a string representing the 0-255 RGB color converted to
    # a 0-9 integer RGB format (i.e. 1000 colors).
    color_tuples = []
    for c in palette['colors']:
        try:
            # Allow code to be set explicitly in palette definition. This is
            # helpful for very dark colors, to allow contrast against the black
            # menu background.
            color_code = c['code']
        except KeyError:
            color_code = _get_color_code(c['rgb'])

        color_tuples.append((color_code, palette['prefix'] + c['name']))

    menu_colors = (name, color_tuples)

    # First `pymol` is the program instance, second is the Python module
    all_colors_list = pymol.pymol.menu.all_colors_list
    if menu_colors in all_colors_list:
        print(f'  - Menu for {name} was already added!')
    else:
        all_colors_list.append(menu_colors)
    print('    done.\n')


def add_menu(palette=None, replace=False):
    '''Add the specified color palettes to the PyMOL OpenGL menu.'''
    palettes = _get_palettes(palette)
    for name, pal in palettes.items():
        _add_palette_menu(name, pal, replace=replace)


def remove_menu(palette=None):
    '''Remove the color palette menu(s).'''
    palettes = _get_palettes(palette)
    all_colors_list = pymol.pymol.menu.all_colors_list
    for p in palettes.keys():
        for i, menu in enumerate(all_colors_list):
            if menu[0] == p:
                del(all_colors_list[i])
                print(f'Deleted menu for {p} palette.')


if __name__ == "pymol":
    add_menu()
