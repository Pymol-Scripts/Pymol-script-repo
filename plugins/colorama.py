'''
See more here: http://www.pymolwiki.org/index.php/Colorama

--- COLORAMA: Coloring Widget for PyMOL ---
Author  : Gregor Hagelueken
Program : Color_select
Date    : Oct 2007
Version : 1.0.0 (Updated from the original Tkinter version 0.1.1)
Mail    : gha@helmholtz-hzi.de

Updated for PyMOL 3.x compatibility (Qt-based GUI)
Update Date: 2026 February 2
Update Author: Converted from Tkinter to PyQt5/PyQt6
Update Author: Blaine Mooers

COLORAMA is a plugin for the PyMOL Molecular Graphics System.
It allows coloring molecules using RGB or HSV colors which can be manually adjusted.
Alternatively, a user defined color gradient can be applied to the molecule.
This version works with PyMOL versions >=3.0.

The program uses a modified version of the color_b program by Robert L. Campbell & James Stroud
for the gradient calculation and the RGBToHTMLColor function by Paul Winkler.

Literature:
 DeLano, W.L. The PyMOL Molecular Graphics System (2002) DeLano Scientific, San Carlos, CA, USA. http://www.pymol.org
'''

from __future__ import print_function
from __future__ import absolute_import

import colorsys
from pymol import cmd, stored

# Import Qt - try PyQt6 first, then PyQt5
try:
    from PyQt6.QtWidgets import (
        QWidget, QDialog, QVBoxLayout, QHBoxLayout, QGridLayout,
        QLabel, QLineEdit, QPushButton, QRadioButton, QButtonGroup,
        QSlider, QFrame, QSizePolicy
    )
    from PyQt6.QtCore import Qt, pyqtSignal
    from PyQt6.QtGui import QColor, QPalette
    PYQT_VERSION = 6
except ImportError:
    from PyQt5.QtWidgets import (
        QWidget, QDialog, QVBoxLayout, QHBoxLayout, QGridLayout,
        QLabel, QLineEdit, QPushButton, QRadioButton, QButtonGroup,
        QSlider, QFrame, QSizePolicy
    )
    from PyQt5.QtCore import Qt, pyqtSignal
    from PyQt5.QtGui import QColor, QPalette
    PYQT_VERSION = 5


class ColorField(QLabel):
    """A colored label widget that displays a solid color."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumSize(30, 80)
        self.setAutoFillBackground(True)
        self.set_color("#808080")

    def set_color(self, hex_color):
        """Set the background color using a hex color string."""
        palette = self.palette()
        palette.setColor(QPalette.ColorRole.Window, QColor(hex_color))
        self.setPalette(palette)


class LabeledSlider(QWidget):
    """A vertical slider with a label above it."""

    valueChanged = pyqtSignal(int)

    def __init__(self, label_text, min_val=0, max_val=255, parent=None):
        super().__init__(parent)
        self.layout = QVBoxLayout(self)
        self.layout.setContentsMargins(5, 5, 5, 5)

        self.label = QLabel(label_text)
        self.label.setAlignment(Qt.AlignmentFlag.AlignCenter)

        self.slider = QSlider(Qt.Orientation.Vertical)
        self.slider.setMinimum(min_val)
        self.slider.setMaximum(max_val)
        self.slider.setValue(0)
        self.slider.setMinimumHeight(100)

        self.value_label = QLabel("0")
        self.value_label.setAlignment(Qt.AlignmentFlag.AlignCenter)

        self.layout.addWidget(self.label)
        self.layout.addWidget(self.slider)
        self.layout.addWidget(self.value_label)

        self.slider.valueChanged.connect(self._on_value_changed)

        # For HSV mode, we need to track if we are using float values
        self._is_float_mode = False
        self._resolution = 1

    def _on_value_changed(self, value):
        if self._is_float_mode:
            float_val = value * self._resolution
            self.value_label.setText(f"{float_val:.2f}")
        else:
            self.value_label.setText(str(value))
        self.valueChanged.emit(value)

    def set_label(self, text):
        self.label.setText(text)

    def set_range(self, min_val, max_val, resolution=1):
        """Set the range of the slider."""
        self._resolution = resolution
        if resolution < 1:
            # Float mode for HSV
            self._is_float_mode = True
            self.slider.setMinimum(int(min_val / resolution))
            self.slider.setMaximum(int(max_val / resolution))
        else:
            self._is_float_mode = False
            self.slider.setMinimum(int(min_val))
            self.slider.setMaximum(int(max_val))

    def get_value(self):
        """Return the current value, accounting for resolution."""
        if self._is_float_mode:
            return self.slider.value() * self._resolution
        return self.slider.value()

    def set_value(self, value):
        """Set the slider value, accounting for resolution."""
        if self._is_float_mode:
            self.slider.setValue(int(value / self._resolution))
        else:
            self.slider.setValue(int(value))


class Colorama(QDialog):
    """Main Colorama dialog window for PyMOL 3.x."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("COLORAMA by gha (Qt version)")
        self.setMinimumSize(400, 300)

        # Initialize state variables
        self.selection = ""
        self.monograd = 'mono'
        self.colorsystem = 'rgb'
        self.farbe12 = 'farbe1'
        self.farbe1 = (128, 128, 128)
        self.farbe2 = (128, 128, 128)

        self._setup_ui()

    def _setup_ui(self):
        """Set up the user interface."""
        main_layout = QVBoxLayout(self)

        # Top section: Object selection
        selection_layout = QHBoxLayout()
        self.input_label = QLabel("Object:")
        self.selection_entry = QLineEdit()
        self.selection_entry.setPlaceholderText("Enter object name")
        self.set_button = QPushButton("Set")
        self.set_button.clicked.connect(self.setselection)
        self.set_gradient_button = QPushButton("Set Gradient")
        self.set_gradient_button.clicked.connect(self.setgradient)

        selection_layout.addWidget(self.input_label)
        selection_layout.addWidget(self.selection_entry)
        selection_layout.addWidget(self.set_button)
        selection_layout.addWidget(self.set_gradient_button)
        main_layout.addLayout(selection_layout)

        # Active selection label
        active_layout = QHBoxLayout()
        self.selection_label = QLabel("Active:")
        self.active_label = QLabel("None")
        active_layout.addWidget(self.selection_label)
        active_layout.addWidget(self.active_label)
        active_layout.addStretch()
        main_layout.addLayout(active_layout)

        # Middle section: Controls and sliders
        middle_layout = QHBoxLayout()

        # Left side: Radio buttons
        radio_layout = QVBoxLayout()

        # Color system radio buttons (RGB/HSV)
        self.color_system_group = QButtonGroup(self)
        self.rgb_button = QRadioButton("RGB")
        self.hsv_button = QRadioButton("HSV")
        self.color_system_group.addButton(self.rgb_button, 1)
        self.color_system_group.addButton(self.hsv_button, 2)
        self.rgb_button.setChecked(True)
        self.rgb_button.clicked.connect(self.Scalergb)
        self.hsv_button.clicked.connect(self.Scalehsv)

        radio_layout.addWidget(self.rgb_button)
        radio_layout.addWidget(self.hsv_button)

        # Separator
        separator1 = QFrame()
        separator1.setFrameShape(QFrame.Shape.HLine)
        radio_layout.addWidget(separator1)

        # Mono/Gradient radio buttons
        self.mono_grad_group = QButtonGroup(self)
        self.mono_button = QRadioButton("Mono")
        self.grad_button = QRadioButton("Gradient")
        self.mono_grad_group.addButton(self.mono_button, 1)
        self.mono_grad_group.addButton(self.grad_button, 2)
        self.mono_button.setChecked(True)
        self.mono_button.clicked.connect(self.Mono)
        self.grad_button.clicked.connect(self.Grad)

        radio_layout.addWidget(self.mono_button)
        radio_layout.addWidget(self.grad_button)

        # Separator
        separator2 = QFrame()
        separator2.setFrameShape(QFrame.Shape.HLine)
        radio_layout.addWidget(separator2)

        # C1/C2 radio buttons (for gradient colors)
        self.farbe_group = QButtonGroup(self)
        self.farbe1_button = QRadioButton("Color 1")
        self.farbe2_button = QRadioButton("Color 2")
        self.farbe_group.addButton(self.farbe1_button, 1)
        self.farbe_group.addButton(self.farbe2_button, 2)
        self.farbe1_button.setChecked(True)
        self.farbe1_button.clicked.connect(self.Farbe1)
        self.farbe2_button.clicked.connect(self.Farbe2)

        radio_layout.addWidget(self.farbe1_button)
        radio_layout.addWidget(self.farbe2_button)
        radio_layout.addStretch()

        middle_layout.addLayout(radio_layout)

        # Color fields
        color_field_layout = QVBoxLayout()
        color_field_layout.addWidget(QLabel("C1"))
        self.colorfield1 = ColorField()
        color_field_layout.addWidget(self.colorfield1)
        color_field_layout.addWidget(QLabel("C2"))
        self.colorfield2 = ColorField()
        color_field_layout.addWidget(self.colorfield2)
        color_field_layout.addStretch()
        middle_layout.addLayout(color_field_layout)

        # Sliders
        slider_layout = QHBoxLayout()
        self.scale_red = LabeledSlider("R", 0, 255)
        self.scale_green = LabeledSlider("G", 0, 255)
        self.scale_blue = LabeledSlider("B", 0, 255)

        self.scale_red.valueChanged.connect(self.setzeFarbe)
        self.scale_green.valueChanged.connect(self.setzeFarbe)
        self.scale_blue.valueChanged.connect(self.setzeFarbe)

        slider_layout.addWidget(self.scale_red)
        slider_layout.addWidget(self.scale_green)
        slider_layout.addWidget(self.scale_blue)

        middle_layout.addLayout(slider_layout)
        main_layout.addLayout(middle_layout)

    def Scalergb(self):
        """Switch to RGB color mode."""
        if self.colorsystem == 'hsv':
            h = self.scale_red.get_value()
            s = self.scale_green.get_value()
            v = self.scale_blue.get_value()
            rgbcolor = colorsys.hsv_to_rgb(h, s, v)
            r = 255 * rgbcolor[0]
            g = 255 * rgbcolor[1]
            b = 255 * rgbcolor[2]

            self.scale_red.set_label('R')
            self.scale_green.set_label('G')
            self.scale_blue.set_label('B')
            self.scale_red.set_range(0, 255, 1)
            self.scale_green.set_range(0, 255, 1)
            self.scale_blue.set_range(0, 255, 1)
            self.scale_red.set_value(r)
            self.scale_green.set_value(g)
            self.scale_blue.set_value(b)
            self.colorsystem = 'rgb'

    def Scalehsv(self):
        """Switch to HSV color mode."""
        if self.colorsystem == 'rgb':
            r = float(self.scale_red.get_value()) / 255
            g = float(self.scale_green.get_value()) / 255
            b = float(self.scale_blue.get_value()) / 255
            hsvcolor = colorsys.rgb_to_hsv(r, g, b)
            h = hsvcolor[0]
            s = hsvcolor[1]
            v = hsvcolor[2]

            self.scale_red.set_label('H')
            self.scale_green.set_label('S')
            self.scale_blue.set_label('V')
            self.scale_red.set_range(0, 1, 0.01)
            self.scale_green.set_range(0, 1, 0.01)
            self.scale_blue.set_range(0, 1, 0.01)
            self.scale_red.set_value(h)
            self.scale_green.set_value(s)
            self.scale_blue.set_value(v)
            self.colorsystem = 'hsv'

    def Mono(self):
        """Set mode to mono (single color)."""
        self.monograd = 'mono'

    def Grad(self):
        """Set mode to gradient."""
        self.monograd = 'grad'

    def Farbe1(self):
        """Select color 1 for editing."""
        self.farbe12 = 'farbe1'
        if self.monograd == 'grad':
            if self.colorsystem == 'rgb':
                self.scale_red.set_value(self.farbe1[0])
                self.scale_green.set_value(self.farbe1[1])
                self.scale_blue.set_value(self.farbe1[2])
            elif self.colorsystem == 'hsv':
                hsvcolor = colorsys.rgb_to_hsv(
                    self.farbe1[0] / 255.0,
                    self.farbe1[1] / 255.0,
                    self.farbe1[2] / 255.0
                )
                self.scale_red.set_value(hsvcolor[0])
                self.scale_green.set_value(hsvcolor[1])
                self.scale_blue.set_value(hsvcolor[2])

    def Farbe2(self):
        """Select color 2 for editing."""
        self.farbe12 = 'farbe2'
        if self.monograd == 'grad':
            if self.colorsystem == 'rgb':
                self.scale_red.set_value(self.farbe2[0])
                self.scale_green.set_value(self.farbe2[1])
                self.scale_blue.set_value(self.farbe2[2])
            elif self.colorsystem == 'hsv':
                hsvcolor = colorsys.rgb_to_hsv(
                    self.farbe2[0] / 255.0,
                    self.farbe2[1] / 255.0,
                    self.farbe2[2] / 255.0
                )
                self.scale_red.set_value(hsvcolor[0])
                self.scale_green.set_value(hsvcolor[1])
                self.scale_blue.set_value(hsvcolor[2])

    def setselection(self):
        """Set the active selection from the entry field."""
        entry_text = self.selection_entry.text().strip()
        if entry_text != "":
            self.selection = entry_text

            # Color of each residue is stored to check if the molecule has a color gradient
            stored.colorlist = []
            try:
                cmd.iterate(self.selection + " & name CA", "stored.colorlist.append(int(color))")
            except Exception:
                pass

            if len(stored.colorlist) == 0:
                # For other objects (e.g. density...)
                try:
                    color_idx = cmd.get_object_color_index(self.selection)
                    stored.colorlist.append(color_idx)
                    stored.colorlist.append(color_idx)
                except Exception:
                    stored.colorlist = [0, 0]

            try:
                initialcolornterm = cmd.get_color_tuple(stored.colorlist[0])
                initialcolorcterm = cmd.get_color_tuple(stored.colorlist[len(stored.colorlist) - 1])
            except Exception:
                initialcolornterm = (0.5, 0.5, 0.5)
                initialcolorcterm = (0.5, 0.5, 0.5)

            self.farbe1 = (
                int(initialcolornterm[0] * 255),
                int(initialcolornterm[1] * 255),
                int(initialcolornterm[2] * 255)
            )
            self.farbe2 = (
                int(initialcolorcterm[0] * 255),
                int(initialcolorcterm[1] * 255),
                int(initialcolorcterm[2] * 255)
            )

            # Set active object to label
            self.active_label.setText(self.selection)

            # Check if there is a gradient and adjust Mono/Grad button
            if initialcolornterm == initialcolorcterm:
                self.mono_button.setChecked(True)
                self.Mono()
            else:
                self.grad_button.setChecked(True)
                self.Grad()

            # Adjust color fields
            self.colorfield1.set_color(self.RGBToHTMLColor(self.farbe1))
            self.colorfield2.set_color(self.RGBToHTMLColor(self.farbe2))
            self.farbe1_button.setChecked(True)
            self.Farbe1()

            # Set scales to initial color of the new object
            if self.colorsystem == 'rgb':
                self.scale_red.set_value(255 * initialcolornterm[0])
                self.scale_green.set_value(255 * initialcolornterm[1])
                self.scale_blue.set_value(255 * initialcolornterm[2])
            elif self.colorsystem == 'hsv':
                hsvcolor = colorsys.rgb_to_hsv(
                    initialcolornterm[0],
                    initialcolornterm[1],
                    initialcolornterm[2]
                )
                self.scale_red.set_value(hsvcolor[0])
                self.scale_green.set_value(hsvcolor[1])
                self.scale_blue.set_value(hsvcolor[2])

    def setzeFarbe(self, event=None):
        """Update color based on slider values."""
        if self.selection != "" and self.monograd == 'mono':
            if self.colorsystem == 'rgb':
                r = int(self.scale_red.get_value())
                g = int(self.scale_green.get_value())
                b = int(self.scale_blue.get_value())
                rgbcolor = (r, g, b)
                hexcolor = self.RGBToHTMLColor(rgbcolor)
                self.colorfield1.set_color(hexcolor)
                self.colorfield2.set_color(hexcolor)
                cmd.delete(self.selection + "_color")
                cmd.set_color(self.selection + "_color", (r / 255.0, g / 255.0, b / 255.0))
                cmd.color(self.selection + "_color", self.selection)
            elif self.colorsystem == 'hsv':
                h = float(self.scale_red.get_value())
                s = float(self.scale_green.get_value())
                v = float(self.scale_blue.get_value())
                rgbcolor_float = colorsys.hsv_to_rgb(h, s, v)
                r = int(255 * rgbcolor_float[0])
                g = int(255 * rgbcolor_float[1])
                b = int(255 * rgbcolor_float[2])
                rgbcolor = (r, g, b)
                hexcolor = self.RGBToHTMLColor(rgbcolor)
                self.colorfield1.set_color(hexcolor)
                self.colorfield2.set_color(hexcolor)
                cmd.delete(self.selection + "_color")
                cmd.set_color(self.selection + "_color", rgbcolor_float)
                cmd.color(self.selection + "_color", self.selection)

        elif self.selection != "" and self.monograd == 'grad':
            if self.colorsystem == 'rgb':
                r = int(self.scale_red.get_value())
                g = int(self.scale_green.get_value())
                b = int(self.scale_blue.get_value())
                rgbcolor = (r, g, b)
                hexcolor = self.RGBToHTMLColor(rgbcolor)
                if self.farbe12 == 'farbe1':
                    self.colorfield1.set_color(hexcolor)
                    self.farbe1 = rgbcolor
                elif self.farbe12 == 'farbe2':
                    self.colorfield2.set_color(hexcolor)
                    self.farbe2 = rgbcolor
            elif self.colorsystem == 'hsv':
                h = float(self.scale_red.get_value())
                s = float(self.scale_green.get_value())
                v = float(self.scale_blue.get_value())
                rgbcolor_float = colorsys.hsv_to_rgb(h, s, v)
                r = int(255 * rgbcolor_float[0])
                g = int(255 * rgbcolor_float[1])
                b = int(255 * rgbcolor_float[2])
                rgbcolor = (r, g, b)
                hexcolor = self.RGBToHTMLColor(rgbcolor)
                if self.farbe12 == 'farbe1':
                    self.colorfield1.set_color(hexcolor)
                    self.farbe1 = rgbcolor
                elif self.farbe12 == 'farbe2':
                    self.colorfield2.set_color(hexcolor)
                    self.farbe2 = rgbcolor

    def setgradient(self):
        """Apply a color gradient to the selection."""
        if self.selection == "":
            return

        stored.residuelist = []
        try:
            cmd.iterate(self.selection, "stored.residuelist.append(int(resi))")
        except Exception:
            print("Error: Could not iterate over selection")
            return

        if len(stored.residuelist) == 0:
            print("Error: No residues found in selection")
            return

        firstresidue = min(stored.residuelist)
        lastresidue = max(stored.residuelist)

        rs = float(self.farbe1[0]) / 255.0
        gs = float(self.farbe1[1]) / 255.0
        bs = float(self.farbe1[2]) / 255.0
        re = float(self.farbe2[0]) / 255.0
        ge = float(self.farbe2[1]) / 255.0
        be = float(self.farbe2[2]) / 255.0

        hsvcolorstart = colorsys.rgb_to_hsv(rs, gs, bs)
        hs = hsvcolorstart[0]
        ss = hsvcolorstart[1]
        vs = hsvcolorstart[2]

        hsvcolorend = colorsys.rgb_to_hsv(re, ge, be)
        he = hsvcolorend[0]
        se = hsvcolorend[1]
        ve = hsvcolorend[2]

        color_grad(
            selection=self.selection,
            minimum=firstresidue,
            maximum=lastresidue,
            hs=hs, he=he,
            ss=ss, se=se,
            vs=vs, ve=ve
        )

    def RGBToHTMLColor(self, rgb_tuple):
        """Convert an (R, G, B) tuple to #RRGGBB."""
        hexcolor = '#%02x%02x%02x' % tuple(map(int, rgb_tuple))
        return hexcolor


# Global reference to keep dialog alive
_colorama_dialog = None


def open_colorama():
    """Open the Colorama dialog."""
    global _colorama_dialog
    if _colorama_dialog is None or not _colorama_dialog.isVisible():
        _colorama_dialog = Colorama()
    _colorama_dialog.show()
    _colorama_dialog.raise_()
    _colorama_dialog.activateWindow()


def __init_plugin__(app=None):
    """
    PyMOL 3.x plugin initialization.
    This function is called by PyMOL when the plugin is loaded.
    """
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('Colorama', open_colorama)


def color_grad(selection='', item='b', mode='hist', gradient='bgr', nbins=11,
               sat=1, value=1, minimum='1', maximum='1', dummy='dummy_all',
               hs=1, he=1, ss=1, se=1, vs=1, ve=1, colorname='init'):
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

    # Make sure sat and value are in the range 0-1.0
    sat = min(sat, 1.0)
    sat = max(sat, 0.0)
    value = min(value, 1.0)
    value = max(value, 0.0)

    # Make sure lowercase
    gradient = gradient.lower()
    mode = mode.lower()

    # Sanity checking
    if nbins == 1:
        print("\n     WARNING: You specified nbins=1, which does not make sense...resetting nbins=11\n")
        nbins = 11

    if mode not in ('hist', 'ramp'):
        print("\n     WARNING: Unknown mode ", mode, "    ----->   Nothing done.\n")
        return

    if selection == '':
        print("\n USAGE: color_grad dimB, minimum=380, maximum=531, hs=0.3, he=0.25,ss=0.7,se=0.2,vs=1,ve=0.5\n")
        return
    elif gradient not in ('bgr', 'rgb', 'rainbow', 'reverserainbow', 'bwr', 'rwb',
                          'bmr', 'rmb', 'rw', 'wr', 'gw', 'wg', 'bw', 'wb', 'gy', 'yg',
                          'gray', 'grey', 'reversegray', 'reversegrey'):
        print("\n     WARNING: Unknown gradient: ", gradient, "    ----->   Nothing done.\n")
        return

    print("MODE, GRADIENT, NBINS:", mode, gradient, nbins)

    # Get list of B-factors from selection
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

        if mode == 'hist':
            # Check if minimum or maximum was specified and use the entered values
            if minimum != '':
                min_b = int(minimum) - 1
            if maximum != '':
                max_b = int(maximum) + 1

            # Histogram: color in bins of equal B-value ranges
            bin_width = (max_b - min_b) / nbins
            sel.append(selection + " and (%s = %4.4g" % (item, min_b + bin_width) + ")")
            for j in range(1, nbins):
                sel.append(dummy + " and %s = %4.4g" % (item, min_b + j * bin_width))

        # Call the function to create the gradient which returns a list of colors
        colours = make_gradient(sel, gradient, nbins, sat, value, hs, he, ss, se, vs, ve, colorname)

        # Do the coloring now
        for j in range(nbins):
            print("Color select: ", sel[j])
            cmd.color(colours[j], sel[j])

    sel = []
    colours = []


def make_gradient(sel, gradient, nbins, sat, value, hs, he, ss, se, vs, ve, colorname):
    """Create a color gradient and return a list of color names."""
    if gradient == 'bgr' or gradient == 'rainbow':
        col = []
        coldesc = []
        for j in range(nbins):
            coldesc.append(colorname + str(j))

            # Create colors using HSV scale
            hsv = (
                hs - (hs - he) * float(j) / (nbins - 1),
                ss - (ss - se) * float(j) / (nbins - 1),
                vs - (vs - ve) * float(j) / (nbins - 1)
            )
            # Convert to RGB and append to color list
            rgb = colorsys.hsv_to_rgb(hsv[0], hsv[1], hsv[2])
            col.append(rgb)
            cmd.set_color(colorname + str(j), col[j])

    # Return the gradient as a list of colors named by their index
    return coldesc


# Register the command with PyMOL
cmd.extend('colorama', open_colorama)
