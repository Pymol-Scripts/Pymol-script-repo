"""
    = vs.py =

    This plugin enables small scale virtual screening with the AutoDock Vina
    software stack. It uses Meeko and Scrubber to prepare molecular ligands,
    and PLIP to analyze the results.
    
    It was tested on PyMOL 3.0 with Python 3.10. Currently supports only
    Linux and probably Mac.

    @author Pedro Sousa Lacerda
    @email pslacerda@gmail.com
    @license BSD-2-Clause
"""

#
# SETUP FPOCKET
#

import platform
import subprocess
import sys
import os.path
from pymol import Qt


#
# SETUP SOME PIP PACKAGES
#

try:
    import numpy as np
    import pandas as pd
    from scipy.spatial import distance_matrix, distance
    from scipy.cluster.hierarchy import dendrogram, linkage
    from scipy.stats import pearsonr
    from matplotlib import pyplot as plt
    import seaborn as sb
    from strenum import StrEnum
    
except ImportError:
    subprocess.check_call(
        [
            "python",
            "-m",
            "pip",
            "--disable-pip-version-check",
            "install",
            "numpy",
            "scipy",
            "jinja2",
            "matplotlib",
            "seaborn",
            "pandas",
            "openpyxl",
            "StrEnum",
        ],
    )


#
# SETUP VINA & MEEKO & PLIP
#

try:
    from vina import Vina
    import meeko
    import openbabel
    import plip

except ImportError:
    subprocess.check_call(
        [
            "conda",
            "install",
            "-y",
            "-c",
            "conda-forge",
            "numpy",
            "swig",
            "boost-cpp",
            "libboost",
            "sphinx",
            "sphinx_rtd_theme",
            "meeko",
            "vina",
            "prody",
            "openbabel",
            "plip"
        ],
    )


#
# SCRUBBER
#

import shutil
import urllib.request

if not shutil.which('scrub.py'):
    zip = "%s/scrubber.zip" % data_dir
    url = "https://github.com/forlilab/scrubber/archive/refs/heads/develop.zip"
    urllib.request.urlretrieve(url, zip)

    subprocess.check_call(
        [
            "python",
            "-m",
            "pip",
            "--disable-pip-version-check",
            "install",
            zip     
        ]
    )


#
# CODE STARTS HERE
#

import os
from os.path import (
    expanduser,
    dirname,
    splitext,
    basename,
)
from glob import glob
import itertools
from operator import itemgetter
import shutil
import shlex
import textwrap
import subprocess
import json
from contextlib import contextmanager
import tempfile

import pymol
import pymol.gui
from pymol import cmd
from pymol.cgo import CYLINDER, SPHERE, COLOR
import numpy as np

QWidget = pymol.Qt.QtWidgets.QWidget
QFileDialog = pymol.Qt.QtWidgets.QFileDialog
QFormLayout = pymol.Qt.QtWidgets.QFormLayout
QPushButton = pymol.Qt.QtWidgets.QPushButton
QSpinBox = pymol.Qt.QtWidgets.QSpinBox
QDoubleSpinBox = pymol.Qt.QtWidgets.QDoubleSpinBox
QDockWidget = pymol.Qt.QtWidgets.QDockWidget
QLineEdit = pymol.Qt.QtWidgets.QLineEdit
QCheckBox = pymol.Qt.QtWidgets.QCheckBox
QApplication = pymol.Qt.QtWidgets.QApplication
QMessageBox = pymol.Qt.QtWidgets.QMessageBox
QVBoxLayout = pymol.Qt.QtWidgets.QVBoxLayout
QTextEdit = pymol.Qt.QtWidgets.QTextEdit
QDialog = pymol.Qt.QtWidgets.QDialog
QDialogButtonBox = pymol.Qt.QtWidgets.QDialogButtonBox
QDesktopWidget = pymol.Qt.QtWidgets.QDesktopWidget
QProgressBar = pymol.Qt.QtWidgets.QProgressBar
QAction = pymol.Qt.QtWidgets.QAction
QComboBox = pymol.Qt.QtWidgets.QComboBox
QTableWidget = pymol.Qt.QtWidgets.QTableWidget
QTableWidgetItem = pymol.Qt.QtWidgets.QTableWidgetItem
QHeaderView = Qt.QtWidgets.QHeaderView

LeftDockWidgetArea = pymol.Qt.QtCore.Qt.LeftDockWidgetArea
QRegExp = pymol.Qt.QtCore.QRegExp
QtCore = pymol.Qt.QtCore
QThread = pymol.Qt.QtCore.QThread
pyqtSignal = pymol.Qt.QtCore.Signal

QRegExpValidator = pymol.Qt.QtGui.QRegExpValidator
QPalette = pymol.Qt.QtGui.QPalette
QTextDocument = pymol.Qt.QtGui.QTextDocument
QIntValidator = pymol.Qt.QtGui.QIntValidator
QTextCursor = pymol.Qt.QtGui.QTextCursor


###############################################
#                Utils                        #
###############################################


def run(command):
    ret = subprocess.run(
        shlex.split(command), stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )
    output = ret.stdout.decode()
    success = ret.returncode == 0
    return output, success


@contextmanager
def chdir(dir):
    cwd = os.getcwd()
    os.chdir(dir)
    yield
    os.chdir(cwd)


class BaseThread(QThread):

    numSteps = pyqtSignal(int)
    currentStep = pyqtSignal(int)

    logEvent = pyqtSignal(str)
    logCodeEvent = pyqtSignal(str)
    logRawEvent = pyqtSignal(str)

    done = pyqtSignal(bool)

    def __init__(self, *args, parent=None):
        super().__init__(parent)
        self.args = args


def display_box_sel(name, sel, margin):
    coords = cmd.get_coords(sel)
    max = np.max(coords, axis=0) + margin
    min = np.min(coords, axis=0) - margin
    display_box(name, max, min)


def display_box(name, max_coords, min_coords):
    #
    # From the original AutoDock plugin
    #

    box = [
        [max_coords[0], min_coords[0]],
        [max_coords[1], min_coords[1]],
        [max_coords[2], min_coords[2]],
    ]
    cylinder_size = 0.2
    color = [1.0, 1.0, 1.0]

    view = cmd.get_view()
    obj = []

    cmd.delete("_box")

    # box_color
    for i in range(2):
        for k in range(2):
            for j in range(2):
                if i != 1:
                    obj.append(CYLINDER)
                    obj.extend([box[0][i], box[1][j], box[2][k]])
                    obj.extend([box[0][i + 1], box[1][j], box[2][k]])
                    obj.append(cylinder_size)
                    obj.extend(color)
                    obj.extend(color)
                    obj.append(COLOR)
                    obj.extend(color)
                    obj.append(SPHERE)
                    obj.extend([box[0][i], box[1][j], box[2][k], cylinder_size])

                if j != 1:
                    obj.append(CYLINDER)
                    obj.extend([box[0][i], box[1][j], box[2][k]])
                    obj.extend([box[0][i], box[1][j + 1], box[2][k]])
                    obj.append(cylinder_size)
                    obj.extend(color)
                    obj.extend(color)
                    obj.append(COLOR)
                    obj.extend(color)
                    obj.append(SPHERE)
                    obj.extend([box[0][i], box[1][j + 1], box[2][k], cylinder_size])
                if k != 1:
                    obj.append(CYLINDER)
                    obj.extend([box[0][i], box[1][j], box[2][k]])
                    obj.extend([box[0][i], box[1][j], box[2][k + 1]])
                    obj.append(cylinder_size)
                    obj.extend(color)
                    obj.extend(color)
                    obj.append(COLOR)
                    obj.extend(color)
                    obj.append(SPHERE)
                    obj.extend([box[0][i], box[1][j], box[2][k + 1], cylinder_size])
    axes = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]]
    cmd.load_cgo(obj, name)
    cmd.set_view(view)


###############################################
#          Load Result Pannel                 #
###############################################


def parse_vina_pdbqt(filename):
    name = basename(filename)
    name = name.rsplit('.', maxsplit=2)[0]
        
    with open(filename) as pdbqt_file:
        for line in pdbqt_file:
            if line.startswith("MODEL"):
                _, mode_txt = line.split()
                mode = int(mode_txt)
            elif line.startswith("REMARK VINA RESULT:"):
                parts = line.split()
                affinity = float(parts[3])
                yield {
                    "name": name,
                    "filename": filename,
                    "affinity": affinity,
                    "mode": mode
                }


class ResultsWidget(QWidget):

    class ResultsTableWidget(QTableWidget):
        def __init__(self, project_data, interactions_check):
            super().__init__()
            self.project_data = project_data
            self.props = ["Name", "Mode", "Affinity"]
            self.interactions_check = interactions_check

            self.setSelectionBehavior(QTableWidget.SelectRows)
            self.setSelectionMode(QTableWidget.SingleSelection)
            self.setColumnCount(3)
            self.setHorizontalHeaderLabels(self.props)
            header = self.horizontalHeader()
            for idx in range(len(self.props)):
                header.setSectionResizeMode(
                    idx, QHeaderView.ResizeMode.ResizeToContents
                )

            @self.itemClicked.connect
            def itemClicked(item):
                name = self.item(item.row(), 0).text()
                mode = self.item(item.row(), 1).text()
                filename = self.project_data['results_dir'] + f'/output/{name}.out.pdbqt'
                cmd.delete('Vina.lig')
                cmd.load(filename, 'Vina.lig', multiplex=True, zoom=False)
                cmd.set_name(f'Vina.lig_{mode.zfill(4)}', 'Vina.lig')
                cmd.delete('Vina.lig_*')
                cmd.group('Vina', 'Vina.lig')
                cmd.alter('Vina.lig', 'chain="Z"')
                cmd.alter('Vina.lig', 'resn="LIG"')
                cmd.alter('Vina.lig', 'resi=1')

                if 'Vina.prot' not in cmd.get_object_list():
                    if self.project_data['flexible']:
                        filename = self.project_data['rigid_pdbqt']
                    else:
                        filename = self.project_data['target_pdbqt']
                    cmd.load(filename, 'Vina.prot')
                    cmd.group('Vina', 'Vina.prot')
                        
                if self.interactions_check:
                    with tempfile.TemporaryDirectory() as tempdir:
                        # tempdir = '/tmp/testarr'
                        pdb_fname = f"{tempdir}/prot_lig.pdb"
                        pse_fname = f'{tempdir}/PLIP/PROT_LIG_PROTEIN_LIG_Z_1.pse'
                        cmd.save(pdb_fname, selection='Vina.*')
                        run(f"plip -f {pdb_fname} -q -s -y --nohydro -o {tempdir}/PLIP")
                        cmd.load(pse_fname)

                        

        def hideEvent(self, evt):
            self.clearSelection()
    
    class SortableItem(QTableWidgetItem):
        def __init__(self, obj):
            super().__init__(str(obj))
            self.setFlags(self.flags() & ~QtCore.Qt.ItemIsEditable)

        def __lt__(self, other):
            try:
                return float(self.text()) < float(other.text())
            except ValueError:
                return self.text() < other.text()

    def __init__(self, project_data, max_load, max_rank, interactions_check):
        super().__init__()
        self.project_data = project_data
        self.max_load = max_load
        self.max_rank = max_rank
        self.load_protein()

        layout = QVBoxLayout()
        self.setLayout(layout)

        self.table = self.ResultsTableWidget(project_data, interactions_check)
        layout.addWidget(self.table)

        
    def load_protein(self):
        cmd.delete('Vina.prot')
        if self.project_data['flexible']:
            filename = self.project_data['rigid_pdbqt']
        else:
            filename = self.project_data['target_pdbqt']
        cmd.load(filename, 'Vina.prot')
        cmd.group('Vina', 'Vina.prot')

    def showEvent(self, event):
        self.refresh()
        super().showEvent(event)

    def refresh(self):
        self.table.setSortingEnabled(False)

        # remove old rows
        while self.table.rowCount() > 0:
            self.table.removeRow(0)

        # append new rows
        results_dir = self.project_data["results_dir"]
        results = itertools.chain.from_iterable(
            map(parse_vina_pdbqt, glob(f"{results_dir}/output/*.out.pdbqt"))
        )
        results = list(sorted(results, key=itemgetter("affinity")))
        for idx, pose in enumerate(results):
            if idx >= self.max_load:
                break
            if pose['mode'] <= self.max_rank:
                self.appendRow(pose)

        self.table.setSortingEnabled(True)

    def appendRow(self, pose):
        self.table.insertRow(self.table.rowCount())
        line = self.table.rowCount() - 1

        self.table.setItem(line, 0, self.SortableItem(pose['name']))
        self.table.setItem(line, 1, self.SortableItem(pose['mode']))
        self.table.setItem(line, 2, self.SortableItem(pose['affinity']))


def new_load_results_widget():
    dockWidget = QDockWidget()
    dockWidget.setWindowTitle("Vina results")

    widget = QWidget()
    layout = QFormLayout(widget)
    widget.setLayout(layout)
    dockWidget.setWidget(widget)

    #
    # Max number of total loaded poses
    #
    max_load_spin = QSpinBox(widget)
    max_load_spin.setRange(1, 99999999)
    max_load_spin.setValue(15)
    max_load_spin.setGroupSeparatorShown(True)

    #
    # Only the best poses of each ligand
    #
    max_rank_spin = QSpinBox(widget)
    max_rank_spin.setRange(1, 20)
    max_rank_spin.setValue(9)
    max_rank_spin.setGroupSeparatorShown(True)

    #
    # Show interactions contacts
    #
    interactions_check = QCheckBox()
    interactions_check.setChecked(False)

    #
    # Choose output folder
    #
    project_button = QPushButton("Load docking...", widget)

    @project_button.clicked.connect
    def load_results():
        nonlocal results_widget
        docking_file = str(
            QFileDialog.getOpenFileName(
                project_button,
                "Docking file",
                expanduser("~"),
                "Docking file (docking.json)",
            )[0]
        )
        if not docking_file:
            return
        
        with open(docking_file, 'r') as file:
            project_data = json.load(file)

        results_widget = ResultsWidget(
            project_data,
            max_load_spin.value(),
            max_rank_spin.value(),
            interactions_check.isChecked(),
        )
        layout.addWidget(results_widget)
    
    #
    # Results Table
    #
    results_widget = None

    #
    # Setup form
    #
    layout.addRow("Max load:", max_load_spin)
    layout.addRow("Max rank:", max_rank_spin)
    layout.addRow("Run PLIP:", interactions_check)
    layout.addWidget(project_button)
    widget.setLayout(layout)

    return dockWidget


###############################################
#          Run Docking Pannel                 #
###############################################


class VinaThreadDialog(QDialog):
    def __init__(self, *vina_args, parent=None):
        super().__init__(parent)
        self.vina = VinaThread(*vina_args)
        self.vina.done.connect(self._done)

        # Setup window
        self.setModal(True)
        self.resize(QDesktopWidget().availableGeometry(self).size() * 0.7)
        self.setWindowFlags(self.windowFlags() | QtCore.Qt.CustomizeWindowHint)
        self.setWindowFlags(self.windowFlags() & ~QtCore.Qt.WindowCloseButtonHint)

        self.layout = QVBoxLayout(self)

        # Setup progress bar
        self.progress = QProgressBar()
        self.layout.addWidget(self.progress)
        self.progress.setValue(0)
        self.vina.numSteps.connect(self.progress.setMaximum)
        self.vina.currentStep.connect(self.progress.setValue)

        # Rich text output
        self.text = QTextEdit(self)
        self.layout.addWidget(self.text)
        self.text.setReadOnly(True)
        self.vina.logEvent.connect(self._appendHtml)
        self.vina.logCodeEvent.connect(self._appendCodeHtml)

        # Ok / Cancel buttons
        self.button_box = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Abort, QtCore.Qt.Horizontal, self
        )
        self.layout.addWidget(self.button_box)
        self.button_box.accepted.connect(self._start)
        self.button_box.rejected.connect(self._abort)

    def _appendHtml(self, html):
        self.text.moveCursor(QTextCursor.End)
        self.text.insertHtml(self._prepareHtml(html))

    def _appendCodeHtml(self, html):
        self.text.moveCursor(QTextCursor.End)
        self.text.insertHtml("<pre>" + self._prepareHtml(html) + "</pre>")

    def _start(self):
        ok_button = self.button_box.button(QDialogButtonBox.Ok)
        ok_button.setDisabled(True)
        self.vina.start()

    def _abort(self):
        self.vina.terminate()
        self.done(QDialog.Rejected)

    def _done(self, success):
        ok_button = self.button_box.button(QDialogButtonBox.Ok)
        abort_button = self.button_box.button(QDialogButtonBox.Abort)

        ok_button.setDisabled(False)
        abort_button.setDisabled(True)

        self.button_box.accepted.disconnect(self._start)

        @self.button_box.accepted.connect
        def _done():
            if success:
                self.accept()
            else:
                self.reject()

    @staticmethod
    def _prepareHtml(html):
        return textwrap.dedent(html)


#
# Run docking software
#
# Assumes that all arguments are ok. For instance the target_sel returns the
# atoms that will be at the final target.pdb file. Only input files and subprocess
# commands will be checked.
#


class VinaThread(BaseThread):
    def run(self):
        (
            results_dir,
            ligands_file,
            target_sel,
            delete_residue_sel,
            flex_sel,
            box_sel,
            box_margin,
            allow_errors,
            ph,
            exhaustiveness,
            num_modes,
            energy_range,
            cpu,
            seed,
        ) = self.args

        self.logEvent.emit("<h2>Preparation</h2>")

        #
        # Check previous output
        #
        if os.listdir(results_dir):
            self.logEvent.emit(f"""
                <br/>
                <font color="red">
                    <b>The output folder is not empty!</b>
                </font>
            """)

        #
        # Prepare target
        #
        target_pdb = f"{results_dir}/target.pdb"
        target_basename = f"{results_dir}/target"
        cmd.save(target_pdb, target_sel)
        command = (
            f"mk_prepare_receptor.py --read_pdb {target_pdb} -o {target_basename} -p"
        )
        if allow_errors:
            command = f"{command} -a"
        if delete_residue_sel != "":
            delete_residues = set()
            for atom in cmd.get_model(delete_residue_sel).atom:
                delete_residues.add(f"{atom.chain}:{atom.resi}")
            delete_residues = ",".join(delete_residues)
            command = f"{command} -d {delete_residues}"
        if flex_sel != "":
            flex_residues = set()
            for atom in cmd.get_model(flex_sel).atom:
                flex_residues.add(f"{atom.chain}:{atom.resi}")
            flex_residues = ",".join(flex_residues)
            command = f"{command} -f {flex_residues}"
        self.logEvent.emit(f"""
            <br/>
            <br/><b>Preparing target.</b>
            <br/><b>Command:</b> {command}
            <br/>
        """)
        output, success = run(command)
        self.logCodeEvent.emit(output)
        if not success:
            self.done.emit(False)
            return

        #
        # Scrubbe isomers
        #
        ligands_sdf = results_dir + "/ligands.sdf"
        command = (
            f"scrub.py -o {ligands_sdf} --ph {ph} --cpu {cpu} {ligands_file}"
        )
        self.logEvent.emit(
            f"""
                <br/>
                <br/><b>Scrubbing ligands.</b>
                <br/><b>Command:</b> {command}
                <br/>
            """
        )
        output, success = run(command)
        self.logCodeEvent.emit(output)
        if not success:
            self.done.emit(False)
            return

        #
        # Convert into PDBQT
        #

        ligands_pdbqt = results_dir + "/ligands_pdbqt"
        command = (
            f"mk_prepare_ligand.py -i {ligands_sdf} --multimol_outdir {ligands_pdbqt}"
        )
        self.logEvent.emit(f"""
            <br/>
            <br/><b>Converting ligands to PDBQT.</b>
            <br/><b>Command:</b> {command}
            <br/>
        """)
        output, success = run(command)
        self.logCodeEvent.emit(output)
        if not success:
            self.done.emit(False)
            return

        #
        # The number of dockings to do
        #
        count = int(output.split('\n')[-5].split(':')[1])
        n_ligands = count
        self.numSteps.emit(count)

        #
        # Compute box variables
        #
        box_coords = cmd.get_coords(box_sel)

        max = np.max(box_coords, axis=0)
        min = np.min(box_coords, axis=0)

        half_size = (max - min) / 2
        center = min + half_size

        size_x, size_y, size_z = (half_size + box_margin) * 2
        center_x, center_y, center_z = center

        size_x, size_y, size_z = (
            round(float(size_x), 2),
            round(float(size_y), 2),
            round(float(size_z), 2),
        )

        center_x, center_y, center_z = (
            round(float(center_x), 2),
            round(float(center_y), 2),
            round(float(center_z), 2),
        )
        
        #
        # Create Vina results directory
        #
        output_dir = f"{results_dir}/output"
        try:
            os.mkdir(output_dir)
        except FileExistsError:
            pass

        #
        # Project data
        #
        project_file = results_dir + "/docking.json"
        project_data = {}
        project_data.update(
            {
                "program": "vina",
                "results_dir": results_dir,
                "ligands_pdbqt": ligands_pdbqt,
                "output_dir": output_dir,
                "size_x": size_x,
                "size_y": size_y,
                "size_z": size_z,
                "center_x": center_x,
                "center_y": center_y,
                "center_z": center_z,
            }
        )
        if flex_sel == "":
            project_data.update(
                {"flexible": False, "target_pdbqt": f"{target_basename}.pdbqt"}
            )
        else:
            project_data.update(
                {
                    "flexible": True,
                    "rigid_pdbqt": f"{target_basename}_rigid.pdbqt",
                    "flex_pdbqt": f"{target_basename}_flex.pdbqt",
                }
            )
        #
        # Prompt for user confirmation
        #

        base_command = (
            f"vina"
            f" --center_x {center_x}"
            f" --center_y {center_y}"
            f" --center_z {center_z}"
            f" --size_x {size_x}"
            f" --size_y {size_y}"
            f" --size_z {size_z}"
            f" --cpu {cpu}"
            f" --seed {seed}"
            f" --exhaustiveness {exhaustiveness}"
            f" --num_modes {num_modes}"
            f" --energy_range {energy_range}"
        )
        self.logEvent.emit(
            f"""
            <br/>
            <h2>Docking</h2>
            <br/>
            <b>Vina base command:</b> {base_command}
        """
        )

        fail_count = 0
        for i, ligand_pdbqt in enumerate(glob(f"{ligands_pdbqt}/*.pdbqt")):
            name, _ = splitext(basename(ligand_pdbqt))
            output_pdbqt = f"{output_dir}/{name}.out.pdbqt"
            log_txt = f"{output_dir}/{name}.log"

            command = base_command + (
                f' --ligand "{ligand_pdbqt}"'
                f' --out "{output_pdbqt}"'
            )
            if project_data["flexible"]:
                rigid_pdbqt = project_data["rigid_pdbqt"]
                flex_pdbqt = project_data["flex_pdbqt"]
                command += f' --receptor "{rigid_pdbqt}"' f' --flex "{flex_pdbqt}"'
            else:
                target_pdbqt = project_data["target_pdbqt"]
                command += f' --receptor "{target_pdbqt}"'

            output, success = run(command)
            self.currentStep.emit(i + 1)
            if not success:
                fail_count += 1
                if fail_count <= 10:
                    self.logEvent.emit(
                        f"""
                        <br/>
                        <font color="red">
                            <b>Vina command failed:</b> {command}
                            <br/>
                            <pre>{output}</pre>
                        </font>
                    """
                    )
                elif fail_count == 11:
                    self.logEvent.emit(
                        f"""
                        <br/>
                        <h3>
                            <font color="red">
                                Too many errors. Omitting output.
                            </font>
                        <h3>f
                    """
                    )

        done_ligands = len(glob(f"{output_dir}/*.out.pdbqt"))

        self.logEvent.emit("<br/><h2>Summary</h2>")
        summary = f"""
            <br/><b>Total expected runs:</b> {n_ligands}
            <br/><b>Total failures:</b> {fail_count}
            <br/><b>Total found PDBQT files:</b> {done_ligands}
        """
        if done_ligands < n_ligands or fail_count > 0:
            self.logEvent.emit(f"<font color='red'>{summary}</font>")
        else:
            self.logEvent.emit(f"{summary}")

        with open(results_dir + f"/docking.json", "w") as docking_file:
            json.dump(project_data, docking_file, indent=4)
        self.done.emit(True)


def new_run_docking_widget():
    dockWidget = QDockWidget()
    dockWidget.setWindowTitle("AutoDock Vina: Run")

    widget = QWidget()

    layout = QFormLayout(widget)
    widget.setLayout(layout)
    dockWidget.setWidget(widget)

    #
    # Receptor selection
    #
    target_sel = QLineEdit("", widget)

    @target_sel.textEdited.connect
    def validate(text):
        validate_target_sel()

    def validate_target_sel():
        text = target_sel.text()
        palette = QApplication.palette(target_sel)
        palette.setColor(QPalette.Base, QtCore.Qt.white)
        valid = True
        try:
            if cmd.count_atoms(f"{text}") == 0:
                raise
        except:
            palette.setColor(QPalette.Base, QtCore.Qt.red)
            valid = False
        target_sel.setPalette(palette)
        return valid

    #
    # Delete residues selection
    #
    delete_residue_sel = QLineEdit("", widget)

    @delete_residue_sel.textEdited.connect
    def validate(text):
        validate_delete_residues_sel()

    def validate_delete_residues_sel():
        text = delete_residue_sel.text()
        palette = QApplication.palette(delete_residue_sel)
        palette.setColor(QPalette.Base, QtCore.Qt.white)
        valid = True

        if text.strip() == "":
            palette.setColor(QPalette.Base, QtCore.Qt.white)
            return True
        try:
            if cmd.count_atoms(f"({text}) and ({target_sel.text()})") == 0:
                raise
            palette.setColor(QPalette.Base, QtCore.Qt.white)
            valid = True
        except:
            palette.setColor(QPalette.Base, QtCore.Qt.red)
            valid = False
        delete_residue_sel.setPalette(palette)
        return valid


    #
    # Flexible residues selection
    #
    flex_sel = QLineEdit("", widget)

    @flex_sel.textEdited.connect
    def validate(text):
        validate_flex_sel()

    def validate_flex_sel():
        text = flex_sel.text()
        palette = QApplication.palette(flex_sel)
        palette.setColor(QPalette.Base, QtCore.Qt.white)
        valid = True

        if text.strip() == "":
            palette.setColor(QPalette.Base, QtCore.Qt.white)
            return True
        try:
            if cmd.count_atoms(f"({text}) and ({target_sel.text()})") == 0:
                raise
            palette.setColor(QPalette.Base, QtCore.Qt.white)
            valid = True
        except:
            palette.setColor(QPalette.Base, QtCore.Qt.red)
            valid = False
        flex_sel.setPalette(palette)
        return valid

    #
    # Box selection
    #
    box_sel = QLineEdit("", widget)
    @box_sel.textEdited.connect
    def validate(text):
        validate_box_sel()
    def validate_box_sel():
        text = box_sel.text()
        palette = QApplication.palette(box_sel)
        palette.setColor(QPalette.Base, QtCore.Qt.white)
        try:
            if cmd.count_atoms(text) == 0:
                raise
        except:
            palette.setColor(QPalette.Base, QtCore.Qt.red)
            box_sel.setPalette(palette)
            cmd.delete("box")
            return False
        display_box_sel("box", text, box_margin_spin.value())
        box_sel.setPalette(palette)
        return True

    #
    # Miscellaneous options
    #
    box_margin_spin = QDoubleSpinBox(widget)
    box_margin_spin.setRange(0.0, 10.0)
    box_margin_spin.setValue(3.0)
    box_margin_spin.setSingleStep(0.1)
    box_margin_spin.setDecimals(1)
    @box_margin_spin.valueChanged.connect
    def display_box(margin):
        cmd.delete("box")
        display_box_sel("box", box_sel.text(), margin)

    allow_errors_check = QCheckBox(widget)
    allow_errors_check.setChecked(False)

    ph_spin = QDoubleSpinBox(widget)
    ph_spin.setRange(0.0, 14.0)
    ph_spin.setValue(7.0)
    ph_spin.setSingleStep(0.1)
    ph_spin.setDecimals(1)

    exhaustiveness_spin = QSpinBox(widget)
    exhaustiveness_spin.setRange(1, 50)
    exhaustiveness_spin.setValue(8)

    num_modes_spin = QSpinBox(widget)
    num_modes_spin.setRange(1, 20)
    num_modes_spin.setValue(9)

    energy_range_spin = QDoubleSpinBox(widget)
    energy_range_spin.setRange(1.0, 10.0)
    energy_range_spin.setValue(3.0)

    cpu_count = QThread.idealThreadCount()
    cpu_spin = QSpinBox(widget)
    cpu_spin.setRange(1, cpu_count)
    cpu_spin.setValue(cpu_count)

    seed_spin = QSpinBox(widget)
    seed_spin.setRange(0, 10000)
    seed_spin.setValue(0)

    #
    # Choose ligand files and run docking
    #
    ligands_file = None
    ligands_button = QPushButton("Choose file...", widget)

    @ligands_button.clicked.connect
    def choose_ligands():
        nonlocal ligands_file
        ligands_file = str(
            QFileDialog.getOpenFileName(
                ligands_button, "Ligand files", expanduser("~"), "SMILES (*.smi *.txt)"
            )[0]
        )
        if not ligands_file:
            return

        ligands_button.setText(basename(ligands_file))

    #
    # Choose output folder
    #
    results_dir = None
    results_button = QPushButton("Choose folder...", widget)

    @results_button.clicked.connect
    def choose_results_dir():
        nonlocal results_dir
        results_dir = str(
            QFileDialog.getExistingDirectory(
                results_button,
                "Output folder",
                expanduser("~"),
                QFileDialog.ShowDirsOnly,
            )
        )
        if not results_dir:
            return

        results_button.setText(basename(results_dir))

    button = QPushButton("Run", widget)

    @button.clicked.connect
    def run():
        if not (validate_target_sel() & validate_flex_sel() & validate_box_sel()):
            return

        if not (results_dir and ligands_file):
            return

        dialog = VinaThreadDialog(
            results_dir,
            ligands_file,
            target_sel.text(),
            delete_residue_sel.text(),
            flex_sel.text(),
            box_sel.text(),
            box_margin_spin.value(),
            allow_errors_check.isChecked(),
            ph_spin.value(),
            exhaustiveness_spin.value(),
            num_modes_spin.value(),
            energy_range_spin.value(),
            cpu_spin.value(),
            seed_spin.value(),
        )
        dialog.exec_()

    layout.addRow("Target:", target_sel)
    layout.addRow("Delete residues:", delete_residue_sel)
    layout.addRow("Flexible residues:", flex_sel)
    layout.addRow("Box:", box_sel)
    layout.addRow("Box margin:", box_margin_spin)
    layout.addRow("Allow errors:", allow_errors_check)
    layout.addRow("Ligand pH:", ph_spin)
    layout.addRow("Exhaustiveness:", exhaustiveness_spin)
    layout.addRow("Number of modes:", num_modes_spin)
    layout.addRow("Energy range:", energy_range_spin)
    layout.addRow("Number of CPUs:", cpu_spin)
    layout.addRow("Seed number:", seed_spin)
    layout.addRow("Ligand file:", ligands_button)
    layout.addRow("Output folder:", results_button)
    layout.addWidget(button)
    widget.setLayout(layout)

    return dockWidget


def __init_plugin__(app=None):

    window = pymol.gui.get_qtwindow()
    menu_bar = window.menuBar()
    vina_menu = menu_bar.addMenu("Vina")

    run_docking_widget = new_run_docking_widget()
    run_docking_action = vina_menu.addAction("Run docking")
    window.addDockWidget(LeftDockWidgetArea, run_docking_widget)
    run_docking_widget.hide()

    @run_docking_action.triggered.connect
    def toggle():
        run_docking_widget.show()

    load_results_widget = new_load_results_widget()
    load_results_action = vina_menu.addAction("Load docking")
    window.addDockWidget(LeftDockWidgetArea, load_results_widget)
    load_results_widget.hide()

    @load_results_action.triggered.connect
    def toggle():
        load_results_widget.show()


if __name__ in ["pymol", "pmg_tk.startup.XDrugPy"]:
    __init_plugin__()
