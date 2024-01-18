'''
See more here: http://www.pymolwiki.org/index.php/dynoplot

###############################################
#  File:          dynoplot.py
#  Author:        Dan Kulp
#  Creation Date: 8/29/05
#
#  Modified 2011-11-17 by Thomas Holder
#  Ported to PyQt 2024 by Thomas Holder
#
#  Notes:
#  Draw plots that display interactive data.
#   Phi,Psi plot shown.
###############################################
'''

from typing import List, Optional, Tuple

from pymol import cmd
from pymol.Qt import QtCore, QtGui, QtWidgets

Qt = QtCore.Qt
QPoint = QtCore.QPoint
QColor = QtGui.QColor


class SimplePlot(QtWidgets.QWidget):

    # Class variables
    mark_size = 4
    xlabels = tuple(range(-180, 181, 60))
    ylabels = tuple(range(-180, 181, 30))

    def __init__(self, width: int = 320, height: int = 320):
        super().__init__()
        self.resize(width, height)

        self.name = ""
        self.spacingx = 0   # spacing in x direction
        self.spacingy = 0
        self.xmin = 0       # min value from each axis
        self.ymin = 0
        self.isdown = 0    # flag for mouse pressed
        self.item: Optional[dict] = None    # items array used for clickable events
        self.shapes: List[dict] = []  # store plot data, x,y etc..
        self.idx2resn = {}  # residue name mapping
        self.symbols = 0    # 0: amino acids, 1: secondary structure

    def closeEvent(self, event) -> None:
        cmd.delete(self.name)

    def mousePressEvent(self, event) -> None:
        button = event.button()
        if button == Qt.LeftButton:
            self.down(event)
        elif button == Qt.RightButton or button == Qt.MiddleButton:
            self.pickWhich(event)

    def mouseReleaseEvent(self, event) -> None:
        button = event.button()
        if button == Qt.LeftButton:
            self.up(event)

    def mouseMoveEvent(self, event) -> None:
        self.drag(event)

    def paintEvent(self, event) -> None:
        painter = QtGui.QPainter()
        painter.begin(self)
        # painter.setRenderHint(QtGui.QPainter.Antialiasing)
        self.paintAxes(painter, event.rect())
        self.paintColorDots(painter)
        painter.end()

    def paintAxes(self, painter, rect) -> None:
        xmin = rect.left() + 40
        xmax = rect.right() - 19
        ymin = rect.top() + 10
        ymax = rect.bottom() - 29
        xint = (ymax - ymin) // 2 + ymin
        yint = xmin

        xlabels = self.xlabels
        ylabels = self.ylabels

        # Store variables in self object
        self.spacingx = (xmax - xmin) / (len(xlabels) - 1)
        self.spacingy = (ymax - ymin) / (len(ylabels) - 1)
        self.xmin = xmin
        self.ymin = ymin

        pen = painter.pen()

        # Create axis lines
        pen.setWidth(3)
        painter.setPen(pen)
        painter.drawLine(xmin, xint, xmax, xint)
        painter.drawLine(yint, ymin, yint, ymax)

        # Text box size
        tsize = 30

        # Create tick marks and labels
        pen.setWidth(2)
        painter.setPen(pen)
        nextspot = xmin
        for label in xlabels:
            inextspot = round(nextspot)
            painter.drawLine(inextspot, xint + 5, inextspot, xint - 5)
            painter.drawText(inextspot - tsize // 2, xint - 10 - tsize,
                             tsize, tsize, Qt.AlignHCenter | Qt.AlignBottom,
                             str(label))

            if len(xlabels) == 1:
                nextspot = xmax
            else:
                nextspot += (xmax - xmin) / (len(xlabels) - 1)

        nextspot = ymax
        for label in ylabels:
            inextspot = round(nextspot)
            painter.drawLine(yint + 5, inextspot, yint - 5, inextspot)
            painter.drawText(yint - 10 - tsize, inextspot - tsize // 2,
                             tsize, tsize, Qt.AlignVCenter | Qt.AlignRight,
                             str(label))
            if len(ylabels) == 1:
                nextspot = ymin
            else:
                nextspot -= (ymax - ymin) / (len(ylabels) - 1)

    # Plot a point
    def plot(self, xp, yp, meta):

        resn, color, ss = self.idx2resn.get(meta)

        if self.symbols == 0:
            # symbols by amino acid (G/P/other)
            mark = {'GLY': 'Tri', 'PRO': 'Rect'}.get(resn, 'Oval')
        else:
            # symbols by secondary structure
            mark = {'H': 'Oval', 'S': 'Rect'}.get(ss, 'Tri')

        if color >= 0x40000000:
            qrgb = color | 0xFF000000
        else:
            qrgb = 0xFF
            for v in cmd.get_color_tuple(color):
                qrgb = (qrgb << 8) | int(0xFF * v)

        self.shapes.append({
            "phi": xp,
            "psi": yp,
            "mark": mark,
            "color": qrgb,
            "meta": meta,
            "dragged": False,
        })

    def paintColorDots(self, painter) -> None:
        pen = painter.pen()
        pen.setWidth(1)
        painter.setPen(pen)

        size = self.mark_size

        for point in self.shapes:
            # Convert from 'label' space to 'pixel' space
            x = round(self.convertToPixel("X", point["phi"]))
            y = round(self.convertToPixel("Y", point["psi"]))

            painter.setBrush(QColor.fromRgb(point["color"]))

            mark = point["mark"]
            if mark == 'Oval':
                painter.drawEllipse(x - size, y - size, 2 * size, 2 * size)
            elif mark == 'Tri':
                polygon = QtGui.QPolygon([
                    QPoint(x, y - size),
                    QPoint(x + size, y + size),
                    QPoint(x - size, y + size),
                ])
                painter.drawConvexPolygon(polygon)
            else:
                painter.drawRect(x - size, y - size, 2 * size, 2 * size)

    # Convert from pixel space to label space
    def convertToLabel(self, axis, value):

        # Defaultly use X-axis info
        label0 = self.xlabels[0]
        label1 = self.xlabels[1]
        spacing = self.spacingx
        min = self.xmin

        # Set info for Y-axis use
        if axis == "Y":
            label0 = self.ylabels[0]
            label1 = self.ylabels[1]
            spacing = self.spacingy
            min = self.ymin

        pixel = value - min
        label = pixel / spacing
        label = label0 + label * abs(label1 - label0)

        if axis == "Y":
            label = - label

        return label

    # Converts value from 'label' space to 'pixel' space
    def convertToPixel(self, axis, value):

        # Defaultly use X-axis info
        label0 = self.xlabels[0]
        label1 = self.xlabels[1]
        spacing = self.spacingx
        min = self.xmin

        # Set info for Y-axis use
        if axis == "Y":
            label0 = self.ylabels[0]
            label1 = self.ylabels[1]
            spacing = self.spacingy
            min = self.ymin

        # Get axis increment in 'label' space
        inc = abs(label1 - label0)

        # 'Label' difference from value and smallest label (label0)
        diff = float(value - label0)

        # Get whole number in 'label' space
        whole = int(diff / inc)

        # Get fraction number in 'label' space
        part = float(float(diff / inc) - whole)

        # Return 'pixel' position value
        pixel = whole * spacing + part * spacing

        # Reverse number by subtracting total number of pixels - value pixels
        if axis == "Y":
            tot_label_diff = float(self.ylabels[-1] - label0)
            tot_label_whole = int(tot_label_diff / inc)
            tot_label_part = float(float(tot_label_diff / inc) - tot_label_whole)
            tot_label_pix = tot_label_whole * spacing + tot_label_part * spacing

            pixel = tot_label_pix - pixel

        # Add min edge pixels
        pixel = pixel + min

        return pixel

    def findClosestShape(self, event) -> Optional[dict]:
        closest: Tuple[float, Optional[dict]] = ((self.mark_size * 4)**2, None)

        ex = event.x()
        ey = event.y()

        for point in self.shapes:
            # Convert from 'label' space to 'pixel' space
            x = self.convertToPixel("X", point["phi"])
            y = self.convertToPixel("Y", point["psi"])
            dist_sq = (x - ex)**2 + (y - ey)**2
            if dist_sq < closest[0]:
                closest = dist_sq, point

        return closest[1]

    # Print out which data point you just clicked on..
    def pickWhich(self, event):

        # Find closest data point
        point = self.findClosestShape(event)

        # Print the shape's meta information corresponding with the shape that was picked
        if point is not None:
            cmd.select('sele', '(%s`%d)' % point["meta"])
            cmd.iterate('sele', 'print(" You clicked /%s/%s/%s/%s`%s/%s (DynoPlot)" %'
                        ' (model, segi, chain, resn, resi, name))')
            cmd.center('byres sele', animate=1)

    # Mouse Down Event
    def down(self, event):
        # Find the currently selected item
        self.item = self.findClosestShape(event)

        # Identify that the mouse is down
        self.isdown = 1

    # Mouse Up Event
    def up(self, event):

        # Get label space version of x,y
        labelx = self.convertToLabel("X", event.x())
        labely = self.convertToLabel("Y", event.y())

        # Convert new position into label space..
        if self.item is not None:
            self.item["phi"] = labelx
            self.item["psi"] = labely
            self.item["dragged"] = True

        # Reset Flags
        self.item = None
        self.isdown = 0

    # Mouse Drag(Move) Event
    def drag(self, event):

        # Check that mouse is down and item clicked is a valid data point
        if self.isdown and self.item is not None:
            labelx = self.convertToLabel("X", event.x())
            labely = self.convertToLabel("Y", event.y())
            self.item["phi"] = labelx
            self.item["psi"] = labely
            self.update()  # repaint


def set_phipsi(model, index, phi, psi, state=-1):
    atsele = [
        'first ((%s`%d) extend 2 and name C)' % (model, index),  # prev C
        'first ((%s`%d) extend 1 and name N)' % (model, index),  # this N
        '(%s`%d)' % (model, index),                             # this CA
        'last ((%s`%d) extend 1 and name C)' % (model, index),  # this C
        'last ((%s`%d) extend 2 and name N)' % (model, index),  # next N
    ]
    try:
        cmd.set_dihedral(atsele[0], atsele[1], atsele[2], atsele[3], phi, state)
        cmd.set_dihedral(atsele[1], atsele[2], atsele[3], atsele[4], psi, state)
    except:
        print(' DynoPlot Error: cmd.set_dihedral failed')

# New Callback object, so that we can update the structure when phi,psi points are moved.


class DynoRamaObject:

    def __init__(self, selection=None, name=None, symbols='', state=-1):
        if name is None:
            try:
                name = cmd.get_unused_name('DynoRama')
            except AttributeError:
                name = 'DynoRamaObject'

        self.symbols = symbols
        self.selection = selection
        self.name = name
        self.lock = 0
        self.state = state

        self._init_canvas()

        if name != 'none':
            auto_zoom = cmd.get('auto_zoom')
            cmd.set('auto_zoom', 0)
            cmd.load_callback(self, name)
            cmd.set('auto_zoom', auto_zoom)

    def _init_canvas(self):
        canvas = SimplePlot()
        canvas.setWindowTitle(' Dynamic Angle Plotting ')
        canvas.destroyed.connect(self.close_callback)

        if self.symbols == 'ss':
            canvas.symbols = 1

        canvas.name = self.name
        self.canvas = canvas

        if self.selection is not None:
            self.start(self.selection)

        canvas.show()

    def __getstate__(self):
        state = dict(self.__dict__)
        state.pop("canvas", None)
        return state

    def __setstate__(self, state: dict):
        self.__dict__.update(state)
        self._init_canvas()

    def close_callback(self, obj=None) -> None:
        cmd.delete(self.name)

    def start(self, sel):
        self.lock = 1
        cmd.iterate('(%s) and name CA' % sel, 'idx2resn[model,index] = (resn, color, ss)',
                    space={'idx2resn': self.canvas.idx2resn})
        for model_index, (phi, psi) in cmd.get_phipsi(sel, self.state).items():
            print(" Plotting Phi,Psi: %8.2f,%8.2f" % (phi, psi))
            self.canvas.plot(phi, psi, model_index)
        self.lock = 0

    def __call__(self):
        if self.lock:
            return

        # Loop through each item on plot to see if updated
        for value in self.canvas.shapes:
            # Look for update flag...
            if value["dragged"]:
                # Set residue's phi,psi to new values
                model, index = value["meta"]
                print(" Re-setting Phi,Psi: %8.2f,%8.2f" % (value["phi"], value["psi"]))
                set_phipsi(model, index, value["phi"], value["psi"], self.state)
                value["dragged"] = False


def rama(sel='(all)', name=None, symbols='aa', filename=None, state=-1):
    '''
DESCRIPTION

    Ramachandran Plot
    http://pymolwiki.org/index.php/DynoPlot

ARGUMENTS

    sel = string: atom selection {default: all}

    name = string: name of callback object which is responsible for setting
    angles when canvas points are dragged, or 'none' to not create a callback
    object {default: DynoRamaObject}

    symbols = string: aa for amino acid or ss for secondary structure {default: aa}

    filename = string: filename for postscript dump of canvas {default: None}
    '''
    dyno = DynoRamaObject(sel, name, symbols, int(state))
    if filename is not None:
        dyno.canvas.postscript(file=filename)

# Extend these commands
cmd.extend('ramachandran', rama)
cmd.auto_arg[0]['ramachandran'] = cmd.auto_arg[0]['zoom']

# Add to plugin menu


def __init_plugin__(self):
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('Rama Plot', lambda: DynoRamaObject('(enabled)'))

# vi:expandtab:smarttab
