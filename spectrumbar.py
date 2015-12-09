from __future__ import print_function

from pymol.cgo import *
from math import *
from pymol import cmd
import re


def spectrumbar(*args, **kwargs):
    """
    Author Sean M. Law
    University of Michigan
    seanlaw_(at)_umich_dot_edu

    USAGE

    While in PyMOL

    run spectrumbar.py

    spectrumbar (RGB_Colors,radius=1.0,name=spectrumbar,head=(0.0,0.0,0.0),tail=(10.0,0.0,0.0),length=10.0, ends=square)

    Parameter     Preset         Type     Description
    RGB_Colors    [1.0,1.0,1.0]  N/A      RGB colors can be specified as a
                                          triplet RGB value or as PyMOL
                                          internal color name (i.e. red)
    radius        1.0            float    Radius of cylindrical spectrum bar
    name          spectrumbar    string   CGO object name for spectrum bar
    head          (0.0,0.0,0.0)  float    Starting coordinate for spectrum bar
    tail          (10.0,0.0,0.0) float    Ending coordinate for spectrum bar
    length        10.0           float    Length of spectrum bar
    ends          square         string   For rounded ends use ends=rounded

    Examples:

    spectrumbar red, green, blue
    spectrumbar 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0

    The above two examples produce the same spectrumbar!

    spectrumbar radius=5.0
    spectrumbar length=20.0

    """

    rgb = [1.0, 1.0, 1.0]
    name = "spectrumbar"
    radius = 1.0
    ends = "square"
    x1 = 0
    y1 = 0
    z1 = 0
    x2 = 10
    y2 = 0
    z2 = 0
    num = re.compile('[0-9]')
    abc = re.compile('[a-z]')

    for key in kwargs:
        if (key == "radius"):
            radius = float(kwargs["radius"])
        elif (key == "name"):
            name = kwargs["name"]
        elif (key == "head"):
            head = kwargs["head"]
            head = head.strip('" []()')
            x1, y1, z1 = list(map(float, head.split(',')))
        elif (key == "tail"):
            tail = kwargs["tail"]
            tail = tail.strip('" []()')
            x2, y2, z2 = list(map(float, tail.split(',')))
        elif (key == "length"):
            if (abc.search(kwargs["length"])):
                print("Error: The length must be a value")
                return
            else:
                x2 = float(kwargs["length"])
        elif (key == "ends"):
            ends = kwargs["ends"]
        elif (key != "_self"):
            print("Ignoring unknown option \"" + key + "\"")
        else:
            continue

    args = list(args)
    if (len(args) >= 1):
        rgb = []
    while (len(args) >= 1):
        if (num.search(args[0]) and abc.search(args[0])):
            if (str(cmd.get_color_tuple(args[0])) != "None"):
                rgb.extend(cmd.get_color_tuple(args.pop(0)))
            else:
                return
        elif (num.search(args[0])):
            rgb.extend([float(args.pop(0))])
        elif (abc.search(args[0])):
            if (str(cmd.get_color_tuple(args[0])) != "None"):
                rgb.extend(cmd.get_color_tuple(args.pop(0)))
            else:
                return
        else:
            print("Error: Unrecognized color format \"" + args[0] + "\"")
            return

    if (len(rgb) % 3):
        print("Error: Missing RGB value")
        print("Please double check RGB values")
        return

    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1
    if (len(rgb) == 3):
        rgb.extend([rgb[0]])
        rgb.extend([rgb[1]])
        rgb.extend([rgb[2]])
    t = 1.0 / (len(rgb) / 3.0 - 1)
    c = len(rgb) / 3 - 1
    s = 0
    bar = []

    while (s < c):
        if (len(rgb) > 0):
            r = rgb.pop(0)
            g = rgb.pop(0)
            b = rgb.pop(0)
        if (s == 0 and ends == "rounded"):
            bar.extend([COLOR, float(r), float(g), float(b), SPHERE, x1 + (s * t) * dx, y1 + (s * t) * dy, z1 + (s * t) * dz, radius])
        bar.extend([CYLINDER])
        bar.extend([x1 + (s * t) * dx, y1 + (s * t) * dy, z1 + (s * t) * dz])
        bar.extend([x1 + (s + 1) * t * dx, y1 + (s + 1) * t * dy, z1 + (s + 1) * t * dz])
        bar.extend([radius, float(r), float(g), float(b)])
        if (len(rgb) >= 3):
            bar.extend([float(rgb[0]), float(rgb[1]), float(rgb[2])])
            r = rgb[0]
            g = rgb[1]
            b = rgb[2]
        else:
            bar.extend([float(r), float(g), float(b)])
        if (s == c - 1 and ends == "rounded"):
            bar.extend([COLOR, float(r), float(g), float(b), SPHERE, x1 + (s + 1) * t * dx, y1 + (s + 1) * t * dy, z1 + (s + 1) * t * dz, radius])
        s = s + 1

    cmd.delete(name)
    cmd.load_cgo(bar, name)

    return
cmd.extend("spectrumbar", spectrumbar)
