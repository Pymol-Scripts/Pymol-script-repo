"""
Representation Outliner
https://pymolwiki.org/index.php/Outline
Author: Jarrett Johnson (Schrodinger, Inc.)
"""

from pymol import cmd
from pymol.Qt import QtCore
from pymol.Qt import QtGui
from pymol.Qt import QtWidgets

import os
from io import BytesIO

from PIL import Image
from PIL import ImageChops
from PIL import ImageDraw
from PIL import ImageFilter
from PIL import ImageOps

VERSION = "0.1"


def __init_plugin__(app=None) -> None:
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('Outliner', run_plugin_gui)


def run_plugin_gui() -> None:
    '''
    Open our custom dialog
    '''
    global dialog

    if dialog is None:
        dialog = make_dialog()

    dialog.show()


# global reference to avoid garbage collection
dialog = None


def _create_clean_overlay(img: Image, target_color: tuple) -> Image:
    """
    Create a clean overlay of the image by removing inner edges
    :param img: Image
    :param target_color: Outline color
    """
    BLACK = (0, 0, 0, 255)
    WHITE = (255, 255, 255, 255)
    TRANSPARENT = (0, 0, 0, 0)

    # Convert to grayscale
    img = img.convert('L')

    # Floodfill approach adopted from
    # https://learnopencv.com/filling-holes-in-an-image-using-opencv-python-c/

    # Threshold
    threshold_val = 255 / 2
    threshold = img.point(lambda p: 255 if p > threshold_val else 0)
    threshold = threshold.convert("RGB")

    # Inverted threshold
    inverted = threshold.copy()
    seed = (0, 0)
    ImageDraw.floodfill(inverted, seed, WHITE[0:3], thresh=200)
    inverted = ImageOps.invert(inverted)

    # Combine
    bitwise_or_add = ImageChops.add(inverted, threshold)

    # Find edges
    edges = bitwise_or_add.filter(ImageFilter.FIND_EDGES).convert("RGBA")

    def _replace_color(img: Image, target_color: tuple,
                       replace_color: tuple) -> None:
        """
        Replace color in image
        :param img: Image
        :param target_color: Color to replace
        :param replace_color: Color to replace with
        """
        data = img.getdata()
        img.putdata(
            [item if item != target_color else replace_color for item in data])

    _replace_color(edges, BLACK, TRANSPARENT)
    _replace_color(edges, WHITE, target_color)

    return edges


def _outline(outline_sele: str, outline_color: tuple, reps: tuple) -> None:
    """
    Outline a selection's surface with a specific color.
    :param outline_sele: Selection to outline
    :param outline_color: Color to outline with
    :param reps: Representations to outline
    """
    tmp_scene = "tmp_scene"

    cmd.scene(tmp_scene, "store", quiet=1)

    # Render what we have
    base_bytes = cmd.png(filename=None, ray=1)

    # Render only whats outlined
    cmd.hide('everything')
    for rep in reps:
        cmd.show(rep, outline_sele)
    ray_trace_mode = cmd.get('ray_trace_mode')

    # Ray trace edges (we'll remove inner edges later)
    cmd.set('ray_trace_mode', 2)

    ray_opaque_background = cmd.get('ray_opaque_background')
    cmd.set('ray_opaque_background', 0)

    overlay_bytes = cmd.png(filename=None, ray=1)

    base = Image.open(BytesIO(base_bytes))
    overlay = Image.open(BytesIO(overlay_bytes))
    overlay = _create_clean_overlay(overlay, outline_color)

    composite = Image.composite(overlay, base, overlay)

    # TODO: load_png doesn't seem to take a filename so we have to save to disk
    tmp_composite_png = "_tmp_outline_comp.png"
    composite.save(tmp_composite_png)
    cmd.load_png(tmp_composite_png, quiet=1)

    # Revert scene and clean up
    cmd.scene(tmp_scene, "recall", quiet=1)
    cmd.scene(tmp_scene, "clear", quiet=1)
    cmd.set('ray_trace_mode', ray_trace_mode)
    cmd.set('ray_opaque_background', ray_opaque_background)
    os.remove(tmp_composite_png)


class StringListSelectorWidgetItem(QtWidgets.QWidget):
    """
    Widget for selecting a string from a list

    :cvar remove_item_signal: Signal for removing item
    :var_type remove_item_signal: QtCore.Signal
    """
    remove_item_signal = QtCore.Signal(QtWidgets.QWidget)

    def __init__(self, str_list, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.layout = QtWidgets.QHBoxLayout(self)
        self.remove_btn = QtWidgets.QPushButton('-')
        self.remove_btn.setFixedSize(20, 20)
        self.string_box = QtWidgets.QComboBox()
        self.string_box.addItems(str_list)
        self.layout.addWidget(self.remove_btn)
        self.layout.addWidget(self.string_box)

        self.remove_btn.clicked.connect(
            lambda: self.remove_item_signal.emit(self))


class StringListSelectorWidget(QtWidgets.QWidget):
    """
    Widget for selecting a string from a list
    """

    def __init__(self, add_label: str, str_list: list, *args,
                 **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.str_list = str_list
        self.layout = QtWidgets.QVBoxLayout(self)
        self.add_btn = QtWidgets.QPushButton(add_label)
        self.add_btn.clicked.connect(self._add_item)

        self.layout.addWidget(self.add_btn)
        self._add_item()

    def _add_item(self) -> None:
        """
        Add item to layout
        """
        item = StringListSelectorWidgetItem(self.str_list)
        self.layout.addWidget(item)
        item.remove_item_signal.connect(self._remove_item)

    def _remove_item(self, item: QtWidgets.QWidget) -> None:
        """
        Remove item from layout
        :param item: Item to remove
        """
        item.setParent(None)
        item.deleteLater()

    def get_rep_list(self) -> list:
        """
        :return: List of representations
        """
        rep_list = []
        for i in range(self.layout.count()):
            widget = self.layout.itemAt(i).widget()
            if isinstance(widget, StringListSelectorWidgetItem):
                rep_list.append(widget.string_box.currentText())
        return rep_list


class SurfaceOutlineDialog(QtWidgets.QDialog):
    """
    Surface Outline Dialog that allows the user to outline a selection's
        surface with a specific color.

    :cvar REP_LIST: List of outlinable representations
    :var_type REP_LIST: list
    """

    REP_LIST = [
        'surface', 'cartoon', 'mesh', 'dots', 'spheres', 'lines', 'nonbonded'
    ]

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.setWindowTitle("Outliner")
        self.default_sele = 'all'

        self.layout = QtWidgets.QVBoxLayout(self)

        # Combobox to hold Object & Selections
        self.combobox = QtWidgets.QComboBox()
        self._refreshCombobox()

        self.rep_list = StringListSelectorWidget(
            add_label='Add Representation', str_list=self.REP_LIST)

        # Outline Button to start Outlining
        self.outline_button = QtWidgets.QPushButton(self.default_sele)

        # Color Picker
        default_color = QtGui.QColor(255, 255, 0)
        self.color_dialogue = QtWidgets.QColorDialog(default_color)
        self.color_dialogue_btn = QtWidgets.QPushButton('Color')

        self._updateCol()

        # Brief note
        self.note = QtWidgets.QLabel(
            "Note: Selection must be completely enveloped by window.")

        self._connectSignals()

        self.layout.addWidget(QtWidgets.QLabel("Selection:"))
        self.layout.addWidget(self.combobox)
        self.layout.addWidget(self.rep_list)
        self.layout.addWidget(self.color_dialogue_btn)
        self.layout.addWidget(self.outline_button)
        self.layout.addWidget(self.note)

    def _updateCol(self) -> None:
        """
        Update color button
        """
        col = self.color_dialogue.currentColor().getRgb()
        self.color_dialogue_btn.setStyleSheet(
            f"background-color: rgb({col[0]},{col[1]},{col[2]});")

    def _connectSignals(self) -> None:
        """
        Connect signals
        """

        # Connect Combobox Signals
        def onComboChanged():
            self.outline_button.setText(
                f"Outline: {self.combobox.currentText()}")

        self.combobox.currentIndexChanged.connect(onComboChanged)

        # Connect Outline Button Signals
        def onOutlineClicked():
            col = self.color_dialogue.currentColor().getRgb()
            _outline(self.combobox.currentText(), col,
                     self.rep_list.get_rep_list())

        self.outline_button.clicked.connect(onOutlineClicked)

        # Connect Color Dialog Signals
        def onColorDialogClicked():
            self.color_dialogue.exec()
            self._updateCol()

        self.color_dialogue_btn.clicked.connect(onColorDialogClicked)

    def _refreshCombobox(self) -> None:
        self.combobox.clear()
        self.combobox.addItem(self.default_sele)
        for obj in cmd.get_names('all'):
            self.combobox.addItem(obj)

    def showEvent(self, event: QtGui.QShowEvent) -> None:
        self._refreshCombobox()
        super().showEvent(event)


def make_dialog() -> SurfaceOutlineDialog:
    """
    Create a Surface Outline Dialog
    """
    dialog = SurfaceOutlineDialog()
    return dialog