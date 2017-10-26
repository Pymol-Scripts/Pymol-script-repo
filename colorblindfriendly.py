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

Copyright (c) 2014-2017 Jared Sampson

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

    0.2.0   Complete overhaul for PyMOL 2.0 with conversion to module format.
            Now, setting the new `cb_*` color values requires a call to the
           `set_colors()` function after import.  You can also now add a
           `cb_colors` menu to the OpenGL GUI via the `add_menu()` function.
            [10/26/2017]

'''
from __future__ import print_function

__author__ = 'Jared Sampson'
__version__ = '0.2.0'

import pymol
from pymol import cmd


# Color blind-friendly color list based on information found at:
# http://jfly.iam.u-tokyo.ac.jp/html/color_blind/#pallet
CB_COLORS = {
    'black': {
        'rgb': [0, 0, 0],
        'alt': None,
    },
    'orange': {
        'rgb': [230, 159, 0],
        'alt': None,
    },
    'sky_blue': {
        'rgb': [86, 180, 233],
        'alt': ['skyblue', 'light_blue', 'lightblue'],
    },
    'bluish_green': {
        'rgb': [0, 158, 115],
        'alt': ['bluishgreen', 'green'],
    },
    'yellow': {
        'rgb': [240, 228, 66],
        'alt': None,
    },
    'blue': {
        'rgb': [0, 114, 178],
        'alt': None,
    },
    'vermillion': {
        'rgb': [213, 94, 0],
        'alt': ['red', 'red_orange', 'redorange'],
    },
    'reddish_purple': {
        'rgb': [204, 121, 167],
        'alt': ['reddishpurple', 'rose', 'violet', 'magenta'],
    },
}


def set_colors(replace=False):
    '''Add the color blind-friendly colors to PyMOL.'''
    # Track the added colors
    added_colors = []

    for color, properties in CB_COLORS.items():
        # RGB tuple shortcut
        rgb = properties['rgb']

        # Get the primary and alternate color names into a single list
        names = [color]
        if properties['alt']:
            names.extend(properties['alt'])

        # Set the colors
        for name in names:
            # Set the cb_color
            cb_name = 'cb_{}'.format(name)
            cmd.set_color(cb_name, rgb)

            # Optionally replace built-in colors
            if replace:
                cmd.set_color(name, rgb)
                spacer = (20 - len(name)) * ' '
                added_colors.append('    {}{}{}'.format(name, spacer, cb_name))
            else:
                added_colors.append('    {}'.format(cb_name))

    # Notify user of newly available colors
    print('\nColor blind-friendly colors are now available:')
    print('\n'.join(added_colors))
    print('')


def add_menu():
    '''Add a color blind-friendly list of colors to the PyMOL OpenGL menu.'''

    # Make sure cb_colors are installed.
    print('Checking for colorblindfriendly colors...')
    try:
        if cmd.get_color_index('cb_red') == -1:
            # mimic pre-1.7.4 behavior
            raise pymol.CmdException
    except pymol.CmdException:
        print('Adding colorblindfriendly colors...')
        set_colors()

    # Abort if PyMOL is too old.
    try:
        from pymol.menu import all_colors_list
    except ImportError:
        print('PyMOL version too old for cb_colors menu. Requires 1.6.0 or later.')
        return

    # Add the menu
    print('Adding cb_colors menu...')
    # mimic pymol.menu.all_colors_list format
    # first color in list is used for menu item color
    cb_colors = ('cb_colors', [
        ('830', 'cb_red'),
        ('064', 'cb_green'),
        ('046', 'cb_blue'),
        ('882', 'cb_yellow'),
        ('746', 'cb_magenta'),
        ('368', 'cb_skyblue'),
        ('860', 'cb_orange'),
    ])
    # First `pymol` is the program instance, second is the Python module
    all_colors_list = pymol.pymol.menu.all_colors_list
    if cb_colors in all_colors_list:
        print('Menu was already added!')
    else:
        all_colors_list.append(cb_colors)
    print('  done.')


def remove_menu():
    '''Remove the cb_colors menu.'''
    all_colors_list = pymol.pymol.menu.all_colors_list
    if all_colors_list[-1][0] == 'cb_colors':
        all_colors_list.pop()
        print('The `cb_colors` menu has been removed.')
    else:
        print('The `cb_colors` menu was not found! Aborting.')


if __name__ == "pymol":
    add_menu()
