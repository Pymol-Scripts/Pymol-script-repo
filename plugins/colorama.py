'''
See more here: http://www.pymolwiki.org/index.php/Colorama

--- COLORAMA: Coloring Widget for PyMOL --- 
Author  : Gregor Hagelueken
Program : Color_select
Date    : Oct 2007
Version : 0.1.1
Mail    : gha@helmholtz-hzi.de

COLORAMA is a plugin for the PyMOL Molecular Graphics System. 
It allows to color molecules using RGB or HSV colors which can be manually adjusted. 
Alternatively, a user defined color gradient can be applied to the molecule.
The program works properly with PyMOL versions >=1.0.
 
The program uses a modified version of the color_b program by Robert L. Campbell & James Stroud
for the gradient calculation and the RGBToHTMLColor function by Paul Winkler.
 
Literature:
 DeLano, W.L. The PyMOL Molecular Graphics System (2002) DeLano Scientific, San Carlos, CA, USA. http://www.pymol.org
'''

from __future__ import print_function
from __future__ import absolute_import

import colorsys
import sys
from pymol import cmd, stored

if sys.version_info[0] < 3:
    from Tkinter import *
else:
    from tkinter import *


class Colorama:

    def __init__(self, master):
        # create frames
        self.F1 = Frame(roota, padx=5, pady=5, bg='red')
        self.F2 = Frame(roota, padx=5, pady=5, bg='green')
        self.F3 = Frame(roota, padx=5, pady=5, bg='blue')
        self.F4 = Frame(self.F1, padx=5, pady=5, bg='yellow')
        self.F5 = Frame(self.F1, padx=5, pady=5, bg='white')
        self.F6 = Frame(self.F1, padx=5, pady=5, bg='pink')

        # color system radiobuttons
        self.Radiocolorsystem = IntVar()
        self.RGB = Radiobutton(self.F3, text='RGB', indicatoron=0, variable=self.Radiocolorsystem, value=1, command=self.Scalergb)
        self.HSV = Radiobutton(self.F3, text='HSV', indicatoron=0, variable=self.Radiocolorsystem, value=2, command=self.Scalehsv)

        # mono/gradient and Farbe1/Farbe2 radiobuttons
        self.RadioMonoGradient = IntVar()
        self.RadioFarbe12 = IntVar()
        self.Monobutton = Radiobutton(self.F3, text='M', indicatoron=0, variable=self.RadioMonoGradient, value=1, command=self.Mono)
        self.Gradbutton = Radiobutton(self.F3, text='G', indicatoron=0, variable=self.RadioMonoGradient, value=2, command=self.Grad)
        self.Farbe1button = Radiobutton(self.F3, text='C1', indicatoron=0, variable=self.RadioFarbe12, value=1, command=self.Farbe1)
        self.Farbe2button = Radiobutton(self.F3, text='C2', indicatoron=0, variable=self.RadioFarbe12, value=2, command=self.Farbe2)

        # preselect RGB and mono
        self.RGB.select()
        self.Monobutton.select()
        self.Farbe1button.select()
        self.monograd = 'mono'
        self.colorsystem = 'rgb'
        self.farbe12 = 'farbe1'

        # initialize the scales
        self.Scales()

        # other GUI elements
        self.selectionentry = Entry(master=self.F5, font=('Arial', 10))
        self.selectionentry.insert(0, "")
        self.selectionbutton = Button(master=self.F5, text='Set', command=self.setselection)
        self.setgradientbutton = Button(master=self.F5, text='Set Gradient', command=self.setgradient)
        self.label = Label(master=self.F4, text="None", font=('Arial', 10))
        self.selectionlabel = Label(master=self.F4, text="Active:", font=('Arial', 10))
        self.inputlabel = Label(master=self.F5, text="Object:", font=('Arial', 10))
        self.colorfield1 = Label(master=self.F3, width=3, height=10)
        self.colorfield2 = Label(master=self.F3, width=3, height=10)

        self.selection = ""
        self.setselection()

        # start layout procedure
        self.layout()

    def layout(self):
        self.F1.pack(side=TOP, anchor=NW)
        self.F4.pack(side=BOTTOM, fill=X, anchor=W)
        self.F5.pack(side=TOP)
        self.F2.pack(side=RIGHT, fill=Y)
        self.F3.pack(side=LEFT, fill=X)

        #entry and buttons
        self.setgradientbutton.pack(side=RIGHT, fill=X, anchor=NE)
        self.selectionbutton.pack(side=RIGHT, anchor=N)
        self.selectionentry.pack(side=RIGHT, fill=X, anchor=NE)

        # labels
        self.inputlabel.pack(side=TOP, anchor=NW)
        self.selectionlabel.pack(side=LEFT, anchor=W)
        self.label.pack(side=LEFT)

        # colorfields
        self.colorfield2.pack(side=RIGHT)
        self.colorfield1.pack(side=RIGHT)

        # scales
        self.ScaleRed.pack(side=RIGHT, fill=Y)
        self.ScaleGreen.pack(side=RIGHT, fill=Y)
        self.ScaleBlue.pack(side=RIGHT, fill=Y)

        # radiobuttons
        self.RGB.pack(side=TOP, fill=X)
        self.HSV.pack(side=TOP, fill=X)
        self.Monobutton.pack(side=TOP, fill=X)
        self.Gradbutton.pack(side=TOP, fill=X)
        self.Farbe1button.pack(side=TOP, fill=X)
        self.Farbe2button.pack(side=TOP, fill=X)

    def Scales(self):
        self.ScaleRed = Scale(master=self.F2, label='R', length='3c',
                              from_=0, to=255,
                              # set(startred),
                              command=self.setzeFarbe)
        self.ScaleGreen = Scale(master=self.F2, label='G', length='3c',
                                from_=0, to=255,
                                # set(startgreen),
                                command=self.setzeFarbe)
        self.ScaleBlue = Scale(master=self.F2, label='B', length='3c',
                               from_=0, to=255,
                               # set(startblue),
                               command=self.setzeFarbe)

    def Scalergb(self):
        if (self.colorsystem == 'hsv'):
            h = float(self.ScaleRed.get())
            s = float(self.ScaleGreen.get())
            v = float(self.ScaleBlue.get())
            rgbcolor = colorsys.hsv_to_rgb(h, s, v)
            r = 255 * rgbcolor[0]
            g = 255 * rgbcolor[1]
            b = 255 * rgbcolor[2]
            self.ScaleRed.config(label='R', from_=0, to=255, resolution=1)
            self.ScaleGreen.config(label='G', from_=0, to=255, resolution=1)
            self.ScaleBlue.config(label='B', from_=0, to=255, resolution=1)
            self.ScaleRed.set(r)
            self.ScaleGreen.set(g)
            self.ScaleBlue.set(b)
            self.colorsystem = 'rgb'

    def Scalehsv(self):
        if (self.colorsystem == 'rgb'):
            r = float(self.ScaleRed.get()) / 255
            g = float(self.ScaleGreen.get()) / 255
            b = float(self.ScaleBlue.get()) / 255
            hsvcolor = colorsys.rgb_to_hsv(r, g, b)
            h = hsvcolor[0]
            s = hsvcolor[1]
            v = hsvcolor[2]
            self.ScaleRed.config(label='H', from_=0, to=1, resolution=0.01)
            self.ScaleGreen.config(label='S', from_=0, to=1, resolution=0.01)
            self.ScaleBlue.config(label='V', from_=0, to=1, resolution=0.01)
            self.ScaleRed.set(h)
            self.ScaleGreen.set(s)
            self.ScaleBlue.set(v)
            self.colorsystem = 'hsv'

    def Mono(self):
        self.monograd = 'mono'

    def Grad(self):
        self.monograd = 'grad'

    def Farbe1(self):
        # Let the scales know which color is to be changed
        self.farbe12 = 'farbe1'
        # set scales to farbe1
        if (self.monograd == 'grad'):
            if (self.colorsystem == 'rgb'):
                startred = self.farbe1[0]
                startgreen = self.farbe1[1]
                startblue = self.farbe1[2]
                self.ScaleRed.set(startred)
                self.ScaleGreen.set(startgreen)
                self.ScaleBlue.set(startblue)
            elif (self.colorsystem == 'hsv'):
                hsvcolor = colorsys.rgb_to_hsv(self.farbe1[0], self.farbe1[1], self.farbe1[2])
                h = hsvcolor[0]
                s = hsvcolor[1]
                v = hsvcolor[2]
                self.ScaleRed.set(h)
                self.ScaleGreen.set(s)
                self.ScaleBlue.set(v)

    def Farbe2(self):
        # Let the scales know which color is to be changed
        self.farbe12 = 'farbe2'
        # set scales to farbe1
        if (self.monograd == 'grad'):
            if (self.colorsystem == 'rgb'):
                startred = self.farbe2[0]
                startgreen = self.farbe2[1]
                startblue = self.farbe2[2]
                self.ScaleRed.set(startred)
                self.ScaleGreen.set(startgreen)
                self.ScaleBlue.set(startblue)
            elif (self.colorsystem == 'hsv'):
                hsvcolor = colorsys.rgb_to_hsv(self.farbe2[0], self.farbe2[1], self.farbe2[2])
                h = hsvcolor[0]
                s = hsvcolor[1]
                v = hsvcolor[2]
                self.ScaleRed.set(h)
                self.ScaleGreen.set(s)
                self.ScaleBlue.set(v)

    def setselection(self):
        if (self.selectionentry.get() != ""):
            self.selection = self.selectionentry.get()

            # Color of each residue is stored in  stored.colorlist to check if the molecule has a colorgradient
            stored.colorlist = []
            cmd.iterate(self.selection + " & name CA", "stored.colorlist.append(int(color))")

            if (len(stored.colorlist) == 0):
                # for other objects (e.g. density...)
                stored.colorlist.append(cmd.get_object_color_index(self.selection))
                stored.colorlist.append(cmd.get_object_color_index(self.selection))

            initialcolornterm = cmd.get_color_tuple(stored.colorlist[0])
            initialcolorcterm = cmd.get_color_tuple(stored.colorlist[len(stored.colorlist) - 1])
            self.farbe1 = initialcolornterm[0] * 255, initialcolornterm[1] * 255, initialcolornterm[2] * 255
            self.farbe2 = initialcolorcterm[0] * 255, initialcolorcterm[1] * 255, initialcolorcterm[2] * 255

            # Set active object to label
            self.label.config(text=self.selection)

            # check if there is a gradient and adjust Mono/Gradbutton
            if (initialcolornterm == initialcolorcterm):
                self.Monobutton.select()
                self.Mono()

            elif (initialcolornterm != initialcolorcterm):
                self.Gradbutton.select()
                self.Grad()

            # adjust colorfields
            self.colorfield1.config(bg=self.RGBToHTMLColor(self.farbe1))
            self.colorfield2.config(bg=self.RGBToHTMLColor(self.farbe2))
            self.Farbe1button.select()
            self.Farbe1()

            # Set scales to initialcolor of the new object
            if (self.colorsystem == 'rgb'):
                startred = 255 * initialcolornterm[0]
                startgreen = 255 * initialcolornterm[1]
                startblue = 255 * initialcolornterm[2]
                self.ScaleRed.set(startred)
                self.ScaleGreen.set(startgreen)
                self.ScaleBlue.set(startblue)
            elif (self.colorsystem == 'hsv'):
                hsvcolor = colorsys.rgb_to_hsv(initialcolornterm[0], initialcolornterm[1], initialcolornterm[2])
                h = hsvcolor[0]
                s = hsvcolor[1]
                v = hsvcolor[2]
                self.ScaleRed.set(h)
                self.ScaleGreen.set(s)
                self.ScaleBlue.set(v)

    def setzeFarbe(self, event):
        if ((self.selection != "") & (self.monograd == 'mono')):
            if (self.colorsystem == 'rgb'):
                col = []
                # read RGB values from scales
                r = int(self.ScaleRed.get())
                g = int(self.ScaleGreen.get())
                b = int(self.ScaleBlue.get())
                rgbcolor = r, g, b
                # Prepare a rgb tupel
                col.append(rgbcolor)
                # hexcolor for colorfields
                hexcolor = self.RGBToHTMLColor(rgbcolor)
                self.colorfield1.config(bg=hexcolor)
                self.colorfield2.config(bg=hexcolor)
                cmd.delete(self.selection + "_color")
                cmd.set_color(self.selection + "_color", col[0])
                cmd.color(self.selection + "_color", self.selection)
                del col[0]
            elif (self.colorsystem == 'hsv'):
                col = []
                # read HSV values from scales
                h = float(self.ScaleRed.get())
                s = float(self.ScaleGreen.get())
                v = float(self.ScaleBlue.get())

                # HSV to RGB and change from 1.0, 1.0, 1.0 format to 255,255,255 format
                rgbcolor = colorsys.hsv_to_rgb(h, s, v)
                r = 255 * rgbcolor[0]
                g = 255 * rgbcolor[1]
                b = 255 * rgbcolor[2]
                # as above
                rgbcolor = r, g, b
                col.append(rgbcolor)
                # hexcolor for colorfields
                hexcolor = self.RGBToHTMLColor(rgbcolor)
                self.colorfield1.config(bg=hexcolor)
                self.colorfield2.config(bg=hexcolor)
                cmd.delete(self.selection + "_color")
                cmd.set_color(self.selection + "_color", col[0])
                cmd.color(self.selection + "_color", self.selection)
                del col[0]
        elif ((self.selection != "") & (self.monograd == 'grad')):

            if (self.colorsystem == 'rgb'):
                col = []
                # read RGB values from scales
                r = int(self.ScaleRed.get())
                g = int(self.ScaleGreen.get())
                b = int(self.ScaleBlue.get())
                rgbcolor = r, g, b
                # Prepare a rgb tupel
                col.append(rgbcolor)
                # hexcolor for colorfields
                hexcolor = self.RGBToHTMLColor(rgbcolor)
                if (self.farbe12 == 'farbe1'):
                    self.colorfield1.config(bg=hexcolor)
                    self.farbe1 = rgbcolor
                elif (self.farbe12 == 'farbe2'):
                    self.colorfield2.config(bg=hexcolor)
                    self.farbe2 = rgbcolor

            elif (self.colorsystem == 'hsv'):
                col = []
                # read HSV values from scales
                h = float(self.ScaleRed.get())
                s = float(self.ScaleGreen.get())
                v = float(self.ScaleBlue.get())

                # HSV to RGB and change from 1.0, 1.0, 1.0 format to 255,255,255 format
                rgbcolor = colorsys.hsv_to_rgb(h, s, v)
                r = 255 * rgbcolor[0]
                g = 255 * rgbcolor[1]
                b = 255 * rgbcolor[2]
                # as above
                rgbcolor = r, g, b
                col.append(rgbcolor)
                # hexcolor for colorfields
                hexcolor = self.RGBToHTMLColor(rgbcolor)

                if (self.farbe12 == 'farbe1'):
                    self.colorfield1.config(bg=hexcolor)
                    self.farbe1 = rgbcolor
                elif (self.farbe12 == 'farbe2'):
                    self.colorfield2.config(bg=hexcolor)
                    self.farbe2 = rgbcolor

    def setgradient(self):

        stored.residuelist = []
        cmd.iterate(self.selection, "stored.residuelist.append(int(resi))")
        firstresidue = min(stored.residuelist)
        lastresidue = max(stored.residuelist)
        rs = float(self.farbe1[0]) / float(255)
        gs = float(self.farbe1[1]) / float(255)
        bs = float(self.farbe1[2]) / float(255)
        re = float(self.farbe2[0]) / float(255)
        ge = float(self.farbe2[1]) / float(255)
        be = float(self.farbe2[2]) / float(255)
        hsvcolorstart = colorsys.rgb_to_hsv(rs, gs, bs)
        hs = hsvcolorstart[0]
        ss = hsvcolorstart[1]
        vs = hsvcolorstart[2]
        hsvcolorend = colorsys.rgb_to_hsv(re, ge, be)
        he = hsvcolorend[0]
        se = hsvcolorend[1]
        ve = hsvcolorend[2]
        color_grad(selection=self.selection, minimum=firstresidue, maximum=lastresidue, hs=hs, he=he, ss=ss, se=se, vs=vs, ve=ve)

    def RGBToHTMLColor(self, rgb_tuple):
            # by Paul Winkler
        """ convert an (R, G, B) tuple to #RRGGBB """
        hexcolor = '#%02x%02x%02x' % tuple(map(int, rgb_tuple))
        # that's it! '%02x' means zero-padded, 2-digit hex values
        return hexcolor


def __init__(self):
    self.menuBar.addmenuitem('Plugin', 'command',
                             'Colorama',
                             label='Colorama',
                             command=lambda s=self: open_Colorama(s.root))


def open_Colorama(master):
    # initialize window (roota)
    global roota
    roota = Toplevel(master)
    roota.title(' COLORAMA by gha')
    global colorama
    colorama = Colorama(roota)


def color_grad(selection='', item='b', mode='hist', gradient='bgr', nbins=11, sat=1, value=1, minimum='1', maximum='1', dummy='dummy_all', hs=1, he=1, ss=1, se=1, vs=1, ve=1, colorname='init'):
    """
      --- color_grad: color gradient tool for PyMOL --- 
      Author  : Gregor Hagelueken
      Program : Color_grad
      Date    : Oct 2007
      Version : 0.1.0
      Mail    : gha@helmholtz-hzi.de




      This is a modified version of the color_b program by Robert L. Campbell & James Stroud

      Literature:
      DeLano, W.L. The PyMOL Molecular Graphics System (2002) DeLano Scientific, San Carlos, CA, USA. http://www.pymol.org

      ----------------------------------------------------------------------
      ----------------------------------------------------------------------
    """

    nbins = int(nbins)
    sat = float(sat)
    value = float(value)
    hs = float(hs)
    he = float(he)
    ss = float(ss)
    se = float(se)
    vs = float(vs)
    ve = float(ve)
    colorname = 'color_' + selection

    nbins = int(maximum) - int(minimum) + 2
    dummy = "dummy-" + selection
    colname = "col" + selection


# make sure sat and value are in the range 0-1.0
    sat = min(sat, 1.0)
    sat = max(sat, 0.0)
    value = min(value, 1.0)
    value = max(value, 0.0)

# make sure lowercase
    gradient.lower()
    mode.lower()

# Sanity checking
    if nbins == 1:
        print("\n     WARNING: You specified nbins=1, which doesn't make sense...resetting nbins=11\n")
        nbins = 11

    if mode not in ('hist', 'ramp'):
        print("\n     WARNING: Unknown mode ", mode, "    ----->   Nothing done.\n")
        return
    if selection == '':
        print("\n USAGE: color_grad dimB, minimum=380, maximum=531, hs=0.3, he=0.25,ss=0.7,se=0.2,vs=1,ve=0.5\n")
        return
    elif gradient not in ('bgr', 'rgb', 'rainbow', 'reverserainbow', 'bwr', 'rwb',
                          'bmr', 'rmb', 'rw', 'wr', 'gw', 'wg', 'bw', 'wb', 'gy', 'yg', 'gray', 'grey', 'reversegray', 'reversegrey'):
        print("\n     WARNING: Unknown gradient: ", gradient, "    ----->   Nothing done.\n")
        return

    print("MODE, GRADIENT, NBINS:", mode, gradient, nbins)

# get list of B-factors from selection
    m = cmd.get_model(selection)
    sel = []
    b_list = []

    if len(m.atom) == 0:
        print("Sorry, no atoms selected")

    else:
        if item == 'b':
            for i in range(len(m.atom)):
                m.atom[i].b = m.atom[i].resi
                b_list.append(m.atom[i].b)

        elif item == 'q':
            for i in range(len(m.atom)):
                b_list.append(m.atom[i].q)

        else:
            print("Not configured to work on item %s" % item)
            return

        cmd.load_model(m, dummy)

        print(selection)
        max_b = maximum
        min_b = minimum
        print("Minimum and Maximum B-values: ", min_b, max_b)
        #nbins = (max_b - min_b)

        if mode == 'hist':

            # check if minimum or maximum was specified and use the entered values
            if minimum != '':
                min_b = int(minimum) - 1
            if maximum != '':
                max_b = int(maximum) + 1
            # histogram:
            # color in bins of equal B-value ranges
            # subtract 0.1 from the lowest B in order to ensure that the single
            # atom with the lowest B value doesn't get omitted
            bin_width = (max_b - min_b) / nbins
            sel.append(selection + " and (%s = %4.4g" % (item, min_b + bin_width) + ")")
            for j in range(1, nbins):
                #sel.append(selection + " and %s > %4.4g" % (item,min_b + j*bin_width))
                sel.append(dummy + " and %s = %4.4g" % (item, min_b + j * bin_width))


# call the function to create the gradient which returns a list of colours
        colours = make_gradient(sel, gradient, nbins, sat, value, hs, he, ss, se, vs, ve, colorname)

# do the colouring now
        for j in range(nbins):
            print("Color select: ", sel[j])
            cmd.color(colours[j], sel[j])
    sel = []
    colours = []
# function for creating the gradient


def make_gradient(sel, gradient, nbins, sat, value, hs, he, ss, se, vs, ve, colorname):
    if gradient == 'bgr' or gradient == 'rainbow':
        col = []
        coldesc = []
        for j in range(nbins):
            # must append the str(sel[j]) to the color name so that it is unique
            # for the selection
            coldesc.append(colorname + str(j))
            # coldesc.append('col' + str(sel[j]) + str(j))

            # create colors using hsv scale (fractional) starting at blue(.6666667)
            # through red(0.00000) in intervals of .6666667/(nbins -1) (the "nbins-1"
            # ensures that the last color is, in fact, red (0)
            # rewrote this to use the colorsys module to convert hsv to rgb
            hsv = (hs - (hs - he) * float(j) / (nbins - 1), ss - (ss - se) * float(j) / (nbins - 1), vs - (vs - ve) * float(j) / (nbins - 1))
            # convert to rgb and append to color list
            rgb = colorsys.hsv_to_rgb(hsv[0], hsv[1], hsv[2])

            col.append(rgb)
            # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
            cmd.set_color(colorname + str(j), col[j])

            # cmd.color(,resi[j])

    # return the gradient as a list of colors named by their index (i.e. col0,col1,col2,col3,...)
    return coldesc
