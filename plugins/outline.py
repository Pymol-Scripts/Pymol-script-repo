"""
Representation Outliner
https://pymolwiki.org/index.php/Outline
Author: Jarrett Johnson (Schrodinger, Inc.)
"""

from __future__ import annotations

from pymol import cmd
from pymol.Qt import QtCore
from pymol.Qt import QtGui
from pymol.Qt import QtWidgets

from dataclasses import dataclass

import os
from io import BytesIO
from pathlib import Path

from PIL import Image
from PIL import ImageChops
from PIL import ImageDraw
from PIL import ImageFilter
from PIL import ImageOps

__version__ = "0.2"


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


def _create_clean_overlay(img: Image, target_color: tuple,
                          outline_width: int) -> Image:
    """
    Create a clean overlay of the image by removing inner edges
    :param img: Image
    :param target_color: Outline color
    :param outline_width: Outline width
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
    edges = bitwise_or_add.\
        filter(ImageFilter.FIND_EDGES).\
        filter(ImageFilter.MaxFilter(outline_width)).\
        convert("RGBA")

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


@dataclass
class Extent2D:
    width: int = 0
    height: int = 0


def _outline(outline_sele: str, outline_color: tuple, outline_width: int,
             scale: int, reps: tuple, image_extent: Extent2D,
             save_file: bool) -> None:
    """
    Outline a selection's representations with a specific color.
    :param outline_sele: Selection to outline
    :param outline_color: Color to outline with
    :param outline_width: Width of outline
    :param scale: Scale factor for antialiasing
    :param reps: Representations to outline
    :param image_extent: Image extent
    :param save_file: Should save to file
    """
    try:
        tmp_scene = "tmp_scene"

        cmd.scene(tmp_scene, "store", quiet=1)

        # Render what we have
        width, height = image_extent.width, image_extent.height
        base_bytes = cmd.png(filename=None, ray=1, width=width, height=height)

        # Render only whats outlined
        cmd.hide('everything')
        for rep in reps:
            cmd.show(rep, outline_sele)

        # Ray trace edges (we'll remove inner edges later)
        ray_trace_mode = cmd.get('ray_trace_mode')
        cmd.set('ray_trace_mode', 2)

        ray_trace_color = cmd.get('ray_trace_color')
        cmd.set('ray_trace_color', 'white')

        bg_color = cmd.get('bg_rgb')
        cmd.bg_color('black')

        ray_opaque_background = cmd.get('ray_opaque_background')
        cmd.set('ray_opaque_background', 0)

        ray_antialias = cmd.get('antialias')
        cmd.set('antialias', 0)
        overlay_bytes = cmd.png(filename=None, ray=1, width=width*scale,
                                height=height*scale)

        base = Image.open(BytesIO(base_bytes))
        overlay = Image.open(BytesIO(overlay_bytes))
        overlay = _create_clean_overlay(overlay, outline_color, outline_width)
        overlay = overlay.resize((width, height), Image.Resampling.LANCZOS)

        composite = Image.composite(overlay, base, overlay)

        composition_file = Path.cwd() / "_tmp_outline_comp.png"
        composition_file_name = str(composition_file)
        composite.save(composition_file_name)
        # TODO: load_png doesn't take raw bytes so we have to save to disk
        cmd.load_png(composition_file_name, quiet=1)
        if save_file:
            new_name = QtWidgets.QFileDialog.getSaveFileName(
                None, "Save File", str(Path.cwd()), "PNG Files (*.png)")[0]
            # rename file
            new_path = Path(new_name)
            if new_path.exists():
                new_path.unlink()
            composition_file.rename(new_path)
            composition_file_name = str(new_path)
            msg = QtWidgets.QMessageBox()
            msg.setWindowTitle("File Saved")
            msg.setText(f"Saved to {composition_file_name}")
            msg.setStandardButtons(QtWidgets.QMessageBox.Ok)
            msg.exec_()

    finally:
        # Revert scene and clean up
        cmd.scene(tmp_scene, "recall", quiet=1)
        cmd.scene(tmp_scene, "clear", quiet=1)
        cmd.set('antialias', ray_antialias)
        cmd.set('ray_trace_mode', ray_trace_mode)
        cmd.set('ray_trace_color', ray_trace_color)
        cmd.set('bg_rgb', bg_color)
        cmd.set('ray_opaque_background', ray_opaque_background)
        if not save_file:
            composition_file.unlink()


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


class ButtonGroup(QtWidgets.QButtonGroup):
    """
    Helper class for creating a group of radio buttons
    """
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.setExclusive(True)

    def addButton(self, label: str) -> QtWidgets.QRadioButton:
        """
        Add button to group
        :param label: Button label
        :return: Button
        """
        button = QtWidgets.QRadioButton(label)
        super().addButton(button)
        return button


class ImageExtentSelectionGroup(QtWidgets.QWidget):
    """
    Two spinboxes for selecting image width and height
    """
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        width, height = cmd.get_viewport()
        self.width = width
        self.height = height
        self.width_box = QtWidgets.QSpinBox()
        self.height_box = QtWidgets.QSpinBox()
        self.width_box.setRange(32, 8192)
        self.height_box.setRange(32, 8192)
        self.width_box.setValue(self.width)
        self.height_box.setValue(self.height)
        self.layout = QtWidgets.QHBoxLayout(self)
        self.layout.addWidget(QtWidgets.QLabel("Img Width:"))
        self.layout.addWidget(self.width_box)
        self.layout.addWidget(QtWidgets.QLabel("Img Height:"))
        self.layout.addWidget(self.height_box)

    def get_extent(self) -> Extent2D:
        """
        :return: Image extent
        """
        return Extent2D(self.width_box.value(), self.height_box.value())


class RepresentationOutlineDialog(QtWidgets.QDialog):
    """
    Representation Outline Dialog that allows the user to outline a selection's
        representations with a specific color.

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
        self.default_kernel_size = 2

        self.layout = QtWidgets.QVBoxLayout(self)

        # Combobox to hold Object & Selections
        self.combobox_layout = QtWidgets.QHBoxLayout()
        self.combobox = QtWidgets.QComboBox()
        self.combobox_refresh_btn = QtWidgets.QPushButton()
        icon = self.style().standardIcon(QtWidgets.QStyle.SP_BrowserReload)
        self.combobox_refresh_btn.setIcon(icon)
        self.combobox_refresh_btn.setFixedSize(25, 25)
        self.combobox_refresh_btn.clicked.connect(self._refreshCombobox)
        self.combobox_layout.addWidget(self.combobox)
        self.combobox_layout.addWidget(self.combobox_refresh_btn)
        self._refreshCombobox()

        self.rep_list = StringListSelectorWidget(
            add_label='(+) Add Representation', str_list=self.REP_LIST)

        # Outline Button to start Outlining
        self.outline_button = QtWidgets.QPushButton(self.default_sele)

        # Color Picker
        default_color = QtGui.QColor(255, 255, 0)
        self.color_dialogue = QtWidgets.QColorDialog(default_color)
        self.color_dialogue_btn = QtWidgets.QPushButton('Color')

        # Width slider
        self.slider_layout = QtWidgets.QHBoxLayout()
        self.width_slider = self._createWidthSlider()
        width_min = self._kernelToWidth(self.width_slider.minimum())
        width_max = self._kernelToWidth(self.width_slider.maximum())
        self.width_label = QtWidgets.QLabel("Width: ")
        self.width_min = QtWidgets.QLabel(str(width_min))
        self.width_max = QtWidgets.QLabel(str(width_max))
        self.slider_layout.addWidget(self.width_label)
        self.slider_layout.addWidget(self.width_min)
        self.slider_layout.addWidget(self.width_slider)
        self.slider_layout.addWidget(self.width_max)

        # Antialias Radiobutton
        self.antialias_layout = QtWidgets.QHBoxLayout()
        self.antialias_group = ButtonGroup()
        self.no_aa = self.antialias_group.addButton(
            "1x Antialias (Fast; Jagged)")
        self.no_aa.setChecked(True)
        self.low_aa = self.antialias_group.addButton(
            "2x Antialias (Slow; Crisp)")
        self.hi_aa = self.antialias_group.addButton(
            "4x Antialias (Very Slow; Beautiful)")
        self.antialias_layout.addWidget(self.no_aa)
        self.antialias_layout.addWidget(self.low_aa)
        self.antialias_layout.addWidget(self.hi_aa)

        # Image Extent Spinbox group
        self.image_extent_selection_group = ImageExtentSelectionGroup()

        # File save field
        self.save_file_checkbox = QtWidgets.QCheckBox("Save to file")

        self._updateCol()

        # Brief note
        note = "Note: Selection must be completely enveloped by window.\n" +\
            "Avoid touching the window while outline is in process."
        self.note = QtWidgets.QLabel(note)

        self._connectSignals()

        self.layout.addWidget(QtWidgets.QLabel("Selection:"))
        self.layout.addLayout(self.combobox_layout)
        self.layout.addWidget(self.rep_list)
        self.layout.addWidget(self.color_dialogue_btn)
        self.layout.addLayout(self.slider_layout)
        self.layout.addLayout(self.antialias_layout)
        self.layout.addWidget(self.image_extent_selection_group)
        self.layout.addWidget(self.save_file_checkbox)
        self.layout.addWidget(self.outline_button)
        self.layout.addWidget(self.note)

    def _kernelToWidth(self, kernel_size: int) -> int:
        """
        Convert kernel size to width
        :param kernel_size: Kernel size
        :return: Width
        """
        return kernel_size * 2 - 1

    def _createWidthSlider(self) -> QtWidgets.QSlider:
        width_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        kernel_min = 1
        kernel_max = 5
        width_slider.setRange(kernel_min, kernel_max)
        width_slider.setValue(self.default_kernel_size)
        width_slider.setTickPosition(QtWidgets.QSlider.TicksBelow)
        return width_slider

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
            scale = int(self.antialias_group.checkedButton().text()[0])
            outline_width = self._kernelToWidth(self.width_slider.value() * scale)
            image_extent = self.image_extent_selection_group.get_extent()
            _outline(self.combobox.currentText(), col,
                     outline_width, scale, self.rep_list.get_rep_list(),
                     image_extent,
                     self.save_file_checkbox.isChecked())

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


def make_dialog() -> RepresentationOutlineDialog:
    """
    Create a Representation Outline Dialog
    """
    dialog = RepresentationOutlineDialog()
    return dialog
