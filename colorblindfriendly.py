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

    The "colorblind" color palette includes:

    * cb_black
    * cb_orange
    * cb_sky_blue (also: cb_skyblue, cb_light_blue, cb_lightblue)
    * cb_bluish_green (also: cb_bluishgreen, cb_green)
    * cb_yellow
    * cb_blue
    * cb_vermillion (also: cb_red, cb_redorange, cb_red_orange)
    * cb_reddish_purple (also: cb_rose, cb_violet, cb_magenta)

    Also added are two palettes from matplotlib, "viridis" and "magma", which
    are designed to be perceptually uniform in both color and black-and-white
    printouts.  These are available as "viridis[1-11]", "magma[1-11]".

USAGE

    With the PyMOL Script Repo installed and importable, import the module and
    set the colors:

    ```
    import colorblindfriendly as cbf

    # Add the new colors
    cbf.set_colors()
    color myObject, cb_red

    # Replace built-in colors of same names with cbf ones
    cbf.set_colors(replace=True)
    color myOtherObject, yellow   # actually cb_yellow

    # Add a `cb_colors` menu item to the OpenGL GUI ([C] menu in the right panel)
    cbf.add_menu()
    ```

    Or, to use without installing, run the script directly from GitHub.  This
    will add the colors and install GUI palette menus for  all three default
    color palettes:

    ```
    run https://github.com/Pymol-Scripts/Pymol-script-repo/blob/master/colorblindfriendly.py
    color myObject, cb_red
    ```

REQUIREMENTS

    The cb_colors menu (`add_menu()` function) requires PyMOL 1.6.0 or later.

AUTHOR

    Jared Sampson
    Github: @jaredsampson

LICENSE

Copyright (c) 2014-2025 Jared Sampson

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

    0.4.0   Add Palette and PaletteColor NamedTuples for cleaner declaration of
            color palettes.  [2025-02-10]

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
from typing import NamedTuple, Optional

__author__ = 'Jared Sampson'
__version__ = '0.4.0'

import pymol
from pymol import cmd


class PaletteColor(NamedTuple):
    '''Named tuple for storing color information.'''
    name: str
    rgb: tuple[int, int, int]
    alt_names: Optional[list[str]] = None
    # Allow code to be set explicitly in palette definition. This is helpful
    # for very dark colors, to allow contrast against the dark menu background.
    short_code: Optional[str] = None  # for GUI menu

    def all_names(self):
        '''Return a list of all names for this color.'''
        names = [self.name]
        if self.alt_names:
            names.extend(self.alt_names)
        return names

    def get_short_code(self):
        '''Return a 3-digit string approximating the RGB color.'''
        if self.short_code:
            return self.short_code
        return ''.join([str(math.floor(x / 256 * 10)) for x in self.rgb])


class Palette(NamedTuple):
    '''Named tuple for storing palette information.'''
    name: str
    colors: list[PaletteColor]
    prefix: str = ''

    def install(self):
        '''Install the palette, adding colors and the GUI menu.'''
        PALETTES_MAP[self.name] = self
        add_menu(self.name)


# Color blind-friendly color list based on information found at:
# http://jfly.iam.u-tokyo.ac.jp/html/color_blind/#pallet
CB_COLORS = [
    PaletteColor('red', (213, 94, 0),
                 ['vermillion', 'red_orange', 'redorange']),
    PaletteColor('orange', (230, 159, 0)),
    PaletteColor('yellow', (240, 228, 66)),
    PaletteColor('green', (0, 158, 115),
                 ['bluish_green', 'bluishgreen']),
    PaletteColor('light_blue', (86, 180, 233),
                 ['lightblue', 'sky_blue', 'skyblue']),
    PaletteColor('blue', (0, 114, 178)),
    PaletteColor('violet', (204, 121, 167),
                 ['reddish_purple', 'reddishpurple', 'rose', 'magenta']),
    PaletteColor('black', (0, 0, 0), short_code='222'),
]
CB_PALETTE = Palette('colorblind', CB_COLORS, prefix='cb_')

# Viridis and Magma palettes contributed by Yehudi Bloch, originally
# developed by Stéfan van der Walt and Nathaniel Smith for matplotlib.
# https://matplotlib.org/stable/users/prev_whats_new/whats_new_1.5.html
VIRIDIS_COLORS = [
    PaletteColor('viridis1',  (253, 231,  36)),
    PaletteColor('viridis2',  (186, 222,  39)),
    PaletteColor('viridis3',  (121, 209,  81)),
    PaletteColor('viridis4',  ( 66, 190, 113)),
    PaletteColor('viridis5',  ( 34, 167, 132)),
    PaletteColor('viridis6',  ( 32, 143, 140)),
    PaletteColor('viridis7',  ( 41, 120, 142)),
    PaletteColor('viridis8',  ( 52,  94, 141)),
    PaletteColor('viridis9',  ( 64,  67, 135)),
    PaletteColor('viridis10', ( 72,  35, 116)),
    PaletteColor('viridis11', ( 68,   1,  84)),
]
VIRIDIS_PALETTE = Palette('viridis', VIRIDIS_COLORS)

MAGMA_COLORS = [
    PaletteColor('magma1',  (251, 252, 191)),
    PaletteColor('magma2',  (253, 205, 114)),
    PaletteColor('magma3',  (253, 159, 108)),
    PaletteColor('magma4',  (246, 110,  91)),
    PaletteColor('magma5',  (221,  73, 104)),
    PaletteColor('magma6',  (181,  54, 121)),
    PaletteColor('magma7',  (140,  41, 128)),
    PaletteColor('magma8',  ( 99,  25, 127)),
    PaletteColor('magma9',  ( 59,  15, 111)),
    PaletteColor('magma10', ( 20,  13,  53)),
    PaletteColor('magma11', (  0,   0,   3)),
]
MAGMA_PALETTE = Palette('magma', MAGMA_COLORS)

PALETTES_MAP = {
    CB_PALETTE.name: CB_PALETTE,
    VIRIDIS_PALETTE.name: VIRIDIS_PALETTE,
    MAGMA_PALETTE.name: MAGMA_PALETTE,
}


def _get_palettes(palette_name: Optional[str] = None):
    '''Return the desired Palette(s).'''
    if palette_name is None:
        return PALETTES_MAP.values()
    if palette_name not in PALETTES_MAP:
        raise ValueError(f'Palette "{palette_name}" not found.')
    else:
        return [PALETTES_MAP[palette_name]]


def set_colors(palette=None, replace=False):
    '''Add the color blind-friendly colors to PyMOL.'''
    palettes = _get_palettes(palette)
    for palette in palettes:
        added_colors = []
        for color in palette.colors:
            # RGB tuple shortcut
            rgb = color.rgb

            # Set the colors
            for name in color.all_names():
                if palette.prefix:
                    use_name = f'{palette.prefix}{name}'
                else:
                    use_name = name
                cmd.set_color(use_name, rgb)

                # Optionally replace built-in colors
                if replace:
                    cmd.set_color(name, rgb)
                    # FIXME hard-coded column width
                    spacer = (20 - len(name)) * ' '
                    added_colors.append(f'    {name}{spacer}{use_name}')
                else:
                    added_colors.append('    {}'.format(use_name))

        # Notify user of newly available colors
        print(f'These {palette.name} colors are now available:')
        print('\n'.join(added_colors))


def _add_palette_menu(palette: Palette):
    '''Add a color palette to the PyMOL OpenGL menu.'''

    # Make sure cb_colors are installed.
    print(f'Checking for {palette.name} colors...')
    try:
        for color in palette.colors:
            if cmd.get_color_index(color.name) == -1:
                # mimic pre-1.7.4 behavior
                raise pymol.CmdException
    except pymol.CmdException:
        print(f'Adding {palette.name} palette colors...')
        set_colors(palette=palette.name)

    # Abort if PyMOL is too old.
    try:
        from pymol.menu import all_colors_list
    except ImportError:
        print('PyMOL version too old for cb_colors menu. Requires 1.6.0 or later.')
        return

    # Add the menu
    print(f'Adding {palette.name} menu...')
    # mimic pymol.menu.all_colors_list format
    # first color in list is used for menu item color

    # Menu item for each color in the menu should be a tuple in the form
    #    ('999', 'color_name')
    # where '999' is a string representing the 0-255 RGB color converted to
    # a 0-9 integer RGB format (i.e. 1000 colors).
    color_tuples = [
        (color.get_short_code(), palette.prefix + color.name)
        for color in palette.colors
    ]
    menu_colors = (palette.name, color_tuples)

    # First `pymol` is the program instance, second is the Python module
    all_colors_list = pymol.pymol.menu.all_colors_list
    if menu_colors in all_colors_list:
        print(f'  - Menu for {palette.name} was already added!')
    else:
        all_colors_list.append(menu_colors)
    print('    done.\n')


def add_menu(palette_name=None):
    '''Add the specified color palettes to the PyMOL OpenGL menu.'''
    palettes = _get_palettes(palette_name)
    for palette in palettes:
        _add_palette_menu(palette)


def remove_menu(palette_name=None):
    '''Remove the color palette menu(s).'''
    palettes = _get_palettes(palette_name)
    all_colors_list = pymol.pymol.menu.all_colors_list
    for palette in palettes:
        initial_length = len(all_colors_list)
        all_colors_list[:] = [
            color_menu for color_menu in all_colors_list
            if color_menu[0] != palette.name
        ]
        if len(all_colors_list) == initial_length:
            print(f'No menu for {palette.name} palette found. Nothing deleted.')
        else:
            print(f'Deleted menu for {palette.name} palette.')


if __name__ == "pymol":
    add_menu()
