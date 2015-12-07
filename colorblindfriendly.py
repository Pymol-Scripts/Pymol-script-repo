'''
More information and examples can be found at:
http://www.pymolwiki.org/index.php/color_blind_friendly

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

    import colorblindfriendly
    color myObject, cb_red
    color mySel, cb_yellow

REQUIREMENTS

    None.

AUTHOR

    Jared Sampson, NYU Langone Medical Center, 2014

LICENSE

Copyright (c) 2014 Jared Sampson

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

'''
from __future__ import print_function
__author__ = 'Jared Sampson'
__version__ = '0.1'

from pymol import cmd

# Color blind friendly color list based on information found at:
# http://jfly.iam.u-tokyo.ac.jp/html/color_blind/#pallet
# The RGB percentage values given on that page are less precise than the 0-255
# values, so the 0-255 values are converted here (e.g. 230/255 = 0.902).
cb_colors = (
    ("black", (0.000, 0.000, 0.000),                   # (  0,   0,   0)
            ()),
    ("orange", (0.902, 0.624, 0.000),                   # (230, 159,   0)
     ()),
    ("sky_blue", (0.337, 0.706, 0.914),                   # ( 86, 180, 233)
     ("skyblue", "light_blue", "lightblue")),
    ("bluish_green", (0.000, 0.620, 0.451),                   # (  0, 158, 115)
     ("bluishgreen", "green")),
    ("yellow", (0.941, 0.894, 0.259),                   # (240, 228,  66)
     ()),
    ("blue", (0.000, 0.447, 0.698),                   # (  0, 114, 178)
     ()),
    ("vermillion", (0.835, 0.369, 0.000),                   # (213,  94,   0)
     ("red", "red_orange", "redorange")),
    ("reddish_purple", (0.800, 0.475, 0.655),                   # (204, 121, 167)
     ("reddishpurple", "rose", "violet", "magenta")),
)

for c in cb_colors:
    # main name
    cmd.set_color("cb_%s" % c[0], c[1])
    print("Set color: cb_%s" % c[0])

    # alternate names
    for alt in c[2]:
        cmd.set_color("cb_%s" % alt, c[1])
        print("           cb_%s" % alt)
