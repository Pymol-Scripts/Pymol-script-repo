"""
    Vina Plugin for PyMOL
    
    Pedro Sousa Lacerda <pslacerda@gmail.com>
    LaBiMM: Laboratório de Bioinformática e Modelagem Molecular
    
    Be sure to adapt the DEFAULT_PREFS for your needs. It is defaulted to
    Debian/Ubuntu Linux.
"""

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
import subprocess
from tempfile import mkstemp
import shutil
import shlex
import textwrap
import subprocess
import json
from contextlib import contextmanager
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pymol
import pymol.gui
from pymol import cmd
from pymol.cgo import CYLINDER, SPHERE, COLOR, cyl_text
from pymol.vfont import plain
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

LeftDockWidgetArea = pymol.Qt.QtCore.Qt.LeftDockWidgetArea
QRegExp = pymol.Qt.QtCore.QRegExp
QtCore = pymol.Qt.QtCore
QThread = pymol.Qt.QtCore.QThread
pyqtSignal = pymol.Qt.QtCore.pyqtSignal

QRegExpValidator = pymol.Qt.QtGui.QRegExpValidator
QPalette = pymol.Qt.QtGui.QPalette
QTextDocument = pymol.Qt.QtGui.QTextDocument
QIntValidator = pymol.Qt.QtGui.QIntValidator
QTextCursor = pymol.Qt.QtGui.QTextCursor


#
# Default preferences
#


DEFAULT_PREFS = {
    'DOCKING_VINA': '/usr/bin/vina',
    'DOCKING_OBABEL': '/usr/bin/obabel',
    'DOCKING_ADT_PYTHON': '/usr/bin/python2.7',
    'DOCKING_PREPARE_RECEPTOR': '/usr/lib/python2.7/dist-packages/AutoDockTools/Utilities24/prepare_receptor4.py',
    'DOCKING_PREPARE_FLEXRECEPTOR': '/usr/lib/python2.7/dist-packages/AutoDockTools/Utilities24/prepare_flexreceptor4.py',
}



###############################################
#                Utils                        #
###############################################

def run(command):
    ret = subprocess.run(
        shlex.split(command),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT
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



def display_box_sel(name, sel):
    coords = cmd.get_coords(sel)
    max = np.max(coords, axis=0)
    min = np.min(coords, axis=0)
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
    color = [1., 1., 1.]

    view = cmd.get_view()
    obj = []
    
    cmd.delete('_box')
    
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


def parse_vina_log(filename):
    with open(filename) as pdbqt_file:
        for line in pdbqt_file:
            if line.startswith('MODEL'):
                _, mode_txt = line.split()
                mode = int(mode_txt)
                state = 'IN_MODE'
            elif line.startswith('REMARK VINA RESULT:'):
                parts = line.split()
                affinity = float(parts[3])
            elif line.startswith('REMARK  Name ='):
                parts = line.split()
                name = parts[3]
                yield {
                    'name': name,
                    'affinity': affinity,
                    'mode': mode,
                    'filename': filename
                }


def load_vina_results(project_file, group, max_load, max_rank, interactions_check):
    
    # Load project data
    with open(project_file) as _project_file:
        project_data = json.load(_project_file)
    
    # Load target
    target_name = f'{group}.target'
    if project_data['flexible']:
        cmd.load(project_data['rigid_pdbqt'], target_name)
    else:
        cmd.load(project_data['target_pdbqt'], target_name)
        
    cmd.group(group)
    cmd.group(group, target_name)
    
    # Show box
    box_name = f'{group}.box'
    display_box(
        box_name,
        (
            project_data['center_x'] + project_data['size_x'] / 2,
            project_data['center_y'] + project_data['size_y'] / 2,
            project_data['center_z'] + project_data['size_z'] / 2,
        ),
        (
            project_data['center_x'] - project_data['size_x'] / 2,
            project_data['center_y'] - project_data['size_y'] / 2,
            project_data['center_z'] - project_data['size_z'] / 2,
        )
    )
    cmd.group(group, box_name)
    
    # Parse results
    results_dir = project_data['results_dir']
    results = itertools.chain.from_iterable(map(
        parse_vina_log,
        glob(f"{results_dir}/poses/*.pdbqt")
    ))
    results = sorted(results, key=itemgetter('affinity'))
    
    cache = set()
    objects = set()
    count = 0
    for pose in results:
        # Ignore poses which mode is greater than max
        if pose["mode"] > max_rank:
            continue
        
        # Load molecule into cache
        cache_name = cmd.get_legal_name(pose['filename'].replace('.', '_'))
        if cache_name not in cache:
            cmd.load(
                pose['filename'],
                cache_name
            )
            cache.add(cache_name)
        
        # Compute object names
        score = int(-10 * pose["affinity"])
        state = pose['mode']
        base_name = f'{group}.{pose["name"]}_{pose["mode"]}_{score}'
        obj_name = f'{base_name}.mol'
        polar_name = f'{base_name}.polar'
        
        # Create group
        cmd.group(base_name)
        
        # Create molecule object
        cmd.create(
            obj_name,
            cache_name,
            state,
            1
        )
        cmd.group(base_name, obj_name)
        
        
        if interactions_check:
            cmd.distance(
                polar_name,
                target_name,
                obj_name,
                2
            )
            cmd.group(base_name, polar_name)
        
        objects.add(obj_name)
        count += 1
        if count >= max_load:
            break
    cmd.delete('delete ' + ' '.join(cache))


def new_load_results_widget():
    dockWidget = QDockWidget()
    dockWidget.setWindowTitle("AutoDock Vina: Load Results")
    
    widget = QWidget()
    
    layout = QFormLayout(widget)
    widget.setLayout(layout)
    dockWidget.setWidget(widget)
    
    #
    # Group name
    #
    group_title = QLineEdit("GroupTitle", widget)
    group_title.setValidator(QRegExpValidator(QRegExp('^[0-9A-Za-z_-]+$')))

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
    interactions_check.setChecked(True)
    
    #
    # Choose output folder
    #
    project_file = None
    project_button = QPushButton("Open docking...", widget)
    @project_button.clicked.connect
    def choose_results_dir():
        nonlocal project_file
        project_file = str(
            QFileDialog.getOpenFileName(
                project_button,
                "Docking file",
                expanduser("~"),
                "Docking project file (docking.json)"
            )[0]
        )
        if not project_file:
            return
        project_button.setText(basename(project_file))
    
    load_button = QPushButton("Load", widget)
    @load_button.clicked.connect
    def run():
        if not project_file:
            return
        
        load_vina_results(
            project_file,
            group_title.text(),
            max_load_spin.value(),
            max_rank_spin.value(),
            interactions_check.isChecked()
        )
    
    #
    # Setup form
    #
    layout.addRow('Title:', group_title)
    layout.addRow('Max load:', max_load_spin)
    layout.addRow('Max rank:', max_rank_spin)
    layout.addRow('Polar contacts:', interactions_check)
    layout.addWidget(project_button)
    layout.addWidget(load_button)
    widget.setLayout(layout)
    
    return dockWidget



###############################################
#          Run Docking Pannel                 #
###############################################


def get_flex_selection(flex_sel, target_sel):
    return f'{flex_sel} and {target_sel}'


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
            QDialogButtonBox.Ok | QDialogButtonBox.Abort,
            QtCore.Qt.Horizontal,
            self
        )
        self.layout.addWidget(self.button_box)
        self.button_box.accepted.connect(self._start)
        self.button_box.rejected.connect(self._abort)
    
    def _appendHtml(self, html):
        self.text.moveCursor(QTextCursor.End)
        self.text.insertHtml(
            self._prepareHtml(html)
        )
    
    def _appendCodeHtml(self, html):
        self.text.moveCursor(QTextCursor.End)
        self.text.insertHtml(
            "<pre>" + self._prepareHtml(html) + "</pre>"
        )
        
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
        return (
            textwrap.dedent(html)
        )
    
#
# Run docking software
#
# Assumes that all arguments are ok. For instance the target_sel returns the
# atoms that will be at the final target.pdb file. Only input files and subprocess
# commands will be checked.
#

class VinaThread(BaseThread):
    
    def run(self):
        (results_dir, ligands_file, target_sel, flex_sel, box_sel, ph,
            exhaustiveness, num_modes, energy_range, cpu, seed) = self.args
        
        self.logEvent.emit("<h2>Preparation</h2>")
        
        #
        # Check if the output
        #
        if os.listdir(results_dir):
            self.logEvent.emit(f"""
                <br/>
                <font color="red">
                    <b>The output folder is not empty!</b>
                </font>
            """)
        
        #
        # Create ligand directory
        #
        
        ligands_dir = results_dir + "/ligands"
        try:
            os.mkdir(ligands_dir)
        except FileExistsError:
            shutil.rmtree(ligands_dir)
            os.mkdir(ligands_dir)
        
        
        #
        # Convert SMILES file into PDBQT
        #
        obabel = pymol.plugins.pref_get('DOCKING_OBABEL')
        command = (
            f'"{obabel}" -i smi "{ligands_file}"'
            f' -ph {ph} --gen3d -m'
            f' -O "{ligands_dir}/.pdbqt"'
        )
        output, success = run(command)
        if success:
            self.logEvent.emit(f"""
                <br/>
                <br/><b>Ligands converted to PDBQT.</b>
                <br/><b>OpenBabel command:</b> {command}
            """)
            self.logCodeEvent.emit(output)
        else:
            self.logEvent.emit(f"""
                <br/>
                <br/><b>Ligands conversion to PDBQT failed.</b>
                <br/><b>OpenBabel command:</b> {command}
            """)
            self.logCodeEvent.emit(output)
            self.done.emit(False)
            return
        
        #
        # Rename PDBQT files accordingly to SMILES
        # Be aware that not every SMILES file has a name column
        #
        count = 0
        lineno = 0
        has_names = True
        with open(ligands_file) as smi:
            for line in smi:
                lineno += 1
                
                # skip empty lines
                if line.strip() == "":
                    continue
                
                count += 1
                if has_names:
                    try:
                        # it really has names
                        smiles, name = line.split()
                    except:
                        if count != 1 and has_names:
                            self.logEvent.emit(f"""
                                <br/>
                                <br/><b>Inconsistent SMILES naming at molecule line #{lineno}.</b>
                                <br/><b><i>Please check you SMILES file.</i></b>
                            """)
                            self.done.emit(False)
                            return
                        
                        # first line don't have name
                        # don't rename files
                        has_names = False
                        continue
                    shutil.move(
                        f"{ligands_dir}/{count}.pdbqt",
                        f"{ligands_dir}/{name}.pdbqt"
                    )
        
        if len(glob(f"{ligands_dir}/*.pdbqt")) != count:
            # The number of generated ligands and SMILES differ
            self.logEvent.emit(f"""
                <br/>
                <br/><b>Number of generated PDBQT files and SMILES molecules differ.</b>
                <br/><b>Please check you SMILES file.</b>
            """)
            self.done.emit(False)
            return
        
        
        #
        # The number of dockings to do
        #
        n_ligands = count
        self.numSteps.emit(count)
        
        
        #
        # Prepare rigid target
        #
        
        target_pdb = f'{results_dir}/target.pdb'
        cmd.save(target_pdb, target_sel)
        
        with chdir(dirname(target_pdb)):
            adt_python = pymol.plugins.pref_get('DOCKING_ADT_PYTHON')
            prepare_target = pymol.plugins.pref_get('DOCKING_PREPARE_RECEPTOR')
            command = (
                f'"{adt_python}"'
                f' "{prepare_target}" -r "{target_pdb}"'
            )
            output, success = run(command)
            if success:
                self.logEvent.emit(f"""
                    <br/>
                    <br/><b>Rigid target prepared.</b>
                    <br/><b>AutoDock command:</b> {command}
                """)
                self.logCodeEvent.emit(output)
            else:
                self.logEvent.emit(f"""
                    <br/>
                    <br/><b>Rigid target preparation failed.</b>
                    <br/><b>AutoDock command:</b> {command}
                """)
                self.logCodeEvent.emit(output)
                self.done.emit(False)
                return
        
        #
        # Prepare flexible target
        #
        if flex_sel != "":
            #
            # Construct residues string
            #
            flex_residues = set()
            for atom in cmd.get_model(flex_sel).atom:
                flex_residues.add(f'{atom.chain}:{atom.resn}{atom.resi}')
            
            flex_residues = ",".join(
                f"target:{res}"
                for res in flex_residues
            )
            
            #
            # Run AutoDock command
            #
            
            target_pdbqt = f'{results_dir}/target.pdbqt'
            with chdir(dirname(target_pdb)):
                adt_python = pymol.plugins.pref_get('DOCKING_ADT_PYTHON')
                prepare_flexreceptor = \
                    pymol.plugins.pref_get('DOCKING_PREPARE_FLEXRECEPTOR')
                command = (
                    f'"{adt_python}"'
                    f'"{prepare_flexreceptor}"'
                    f' -r "{target_pdbqt}"'
                    f' -s {flex_residues}'
                )
                output, success = run(command)
                if success:
                    self.logEvent.emit(f"""
                        <br/>
                        <br/><b>Flexible target prepared.</b>
                        <br/><b>AutoDock command:</b> {command}
                    """)
                    self.logCodeEvent.emit(output)
                else:
                    self.logEvent.emit(f"""
                        <br/>
                        <br/><b>Flexible target preparation failed.</b>
                        <br/><b>AutoDock command:</b> {command}
                    """)
                    self.logCodeEvent.emit(output)
                    self.done.emit(False)
                    return
        
        #
        # Create Vina results directory
        #
        
        output_dir = f'{results_dir}/poses'
        try:
            os.mkdir(output_dir)
        except FileExistsError:
            pass
        
        #
        # Compute box variables
        #
        box_coords = cmd.get_coords(box_sel)
        
        max = np.max(box_coords, axis=0)
        min = np.min(box_coords, axis=0)
        
        half_size = (max - min) / 2
        center = min + half_size
        
        size_x, size_y, size_z = half_size * 2
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
        # Project data
        #
        
        project_file = results_dir + "/docking.json"
        project_data = {}
        
        project_data.update({
            'program': 'vina',
            'results_dir': results_dir,
            'ligands_dir': ligands_dir,
            'output_dir': output_dir,
            'size_x': size_x,
            'size_y': size_y,
            'size_z': size_z,
            'center_x': center_x,
            'center_y': center_y,
            'center_z': center_z,
        })
        
        if flex_sel == "":
            project_data.update({
                'flexible': False,
                'target_pdbqt': f'{results_dir}/target.pdbqt'
            })
        else:
            project_data.update({
                'flexible': True,
                'rigid_pdbqt': f'{results_dir}/target_rigid.pdbqt',
                'flex_pdbqt': f'{results_dir}/target_flex.pdbqt'
            })
        print(project_data)
        #
        # Prompt for user confirmation
        #
        
        base_command = (
            f'vina'
            f' --center_x {center_x}'
            f' --center_y {center_y}'
            f' --center_z {center_z}'
            f' --size_x {size_x}'
            f' --size_y {size_y}'
            f' --size_z {size_z}'
            f' --cpu {cpu}'
            f' --seed {seed}'
            f' --exhaustiveness {exhaustiveness}'
            f' --num_modes {num_modes}'
            f' --energy_range {energy_range}'
        )
        self.logEvent.emit(f"""
            <br/>
            <h2>Docking</h2>
            <br/>
            <b>Vina base command:</b> {base_command}
        """)
        
        fail_count = 0
        for i, ligand_pdbqt in enumerate(glob(f"{ligands_dir}/*.pdbqt")):
            name, _ = splitext(basename(ligand_pdbqt))
            output_pdbqt = f'{output_dir}/{name}.out.pdbqt'
            log_txt = f'{output_dir}/{name}.log'
            
            command = base_command + (
                f' --ligand "{ligand_pdbqt}"'
                f' --out "{output_pdbqt}"'
                f' --log "{log_txt}"'
            )
            if project_data['flexible']:
                rigid_pdbqt = project_data['rigid_pdbqt']
                flex_pdbqt = project_data['flex_pdbqt']
                command += (
                    f' --receptor "{rigid_pdbqt}"'
                    f' --flex "{flex_pdbqt}"'
                )
            else:
                target_pdbqt = project_data['target_pdbqt']
                command += (
                    f' --receptor "{target_pdbqt}"'
                )
            
            output, success = run(command)
            self.currentStep.emit(i+1)
            if not success:
                fail_count += 1
                if fail_count <= 10:
                    self.logEvent.emit(f"""
                        <br/>
                        <font color="red">
                            <b>Vina command failed:</b> {command}
                            <br/>
                            <pre>{output}</pre>
                        </font>
                    """)
                elif fail_count == 11:
                    self.logEvent.emit(f"""
                        <br/>
                        <h3>
                            <font color="red">
                                Too many errors. Omitting output.
                            </font>
                        <h3>f
                    """)
        
        done_ligands = len(glob(f'{output_dir}/*.out.pdbqt'))
        
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
        
        with open(results_dir + f'/docking.json', 'w') as docking_file:
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
            if cmd.count_atoms(f'{text}') == 0:
                raise
        except:
            palette.setColor(QPalette.Base, QtCore.Qt.red)
            valid = False
        target_sel.setPalette(palette)
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
            return True
        try:
            if cmd.count_atoms(f'({text}) and ({target_sel.text()})') == 0:
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
        valid = True
        try:
            if cmd.count_atoms(text) == 0:
                raise
        except:
            palette.setColor(QPalette.Base, QtCore.Qt.red)
            box_sel.setPalette(palette)
            cmd.delete('box')
            return False
        display_box_sel('box', text)
        box_sel.setPalette(palette)
        return True
    
    #
    # Miscellaneous options
    #
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
                ligands_button,
                "Ligand files",
                expanduser("~"),
                "SMILES (*.smi *.txt)"
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
                QFileDialog.ShowDirsOnly
            )
        )
        if not results_dir:
            return
        
        results_button.setText(basename(results_dir))
    
    button = QPushButton("Run", widget)
    @button.clicked.connect
    def run():
        if not (validate_target_sel()
                    & validate_flex_sel()
                    & validate_box_sel()):
            return
        
        if not (results_dir and ligands_file):
            return
        
        dialog = VinaThreadDialog(
            results_dir,
            ligands_file,
            target_sel.text(),
            flex_sel.text(),
            box_sel.text(),
            ph_spin.value(),
            exhaustiveness_spin.value(),
            num_modes_spin.value(),
            energy_range_spin.value(),
            cpu_spin.value(),
            seed_spin.value()
        )
        dialog.exec_()
    
    layout.addRow('Target:', target_sel)
    layout.addRow('Flexible chains:', flex_sel)
    layout.addRow('Box:', box_sel)
    layout.addRow('Ligand pH:', ph_spin)
    layout.addRow('Exhaustiveness:', exhaustiveness_spin)
    layout.addRow('Number of modes:', num_modes_spin)
    layout.addRow('Energy range:', energy_range_spin)
    layout.addRow('Number of CPUs:', cpu_spin)
    layout.addRow('Seed number:', seed_spin)
    layout.addRow("Ligand file:", ligands_button)
    layout.addRow("Output folder:", results_button)
    layout.addWidget(button)
    widget.setLayout(layout)
    
    return dockWidget


def __init_plugin__(app=None):
    for pref in DEFAULT_PREFS:
#        if not pymol.plugins.pref_get(pref):
            pymol.plugins.pref_set(pref, DEFAULT_PREFS[pref])

    #cmd.set('group_auto_mode', 2)
    
    window = pymol.gui.get_qtwindow()
    menu_bar = window.menuBar()
    vina_menu = menu_bar.addMenu("Vina")
    
    run_docking_widget = new_run_docking_widget()
    run_docking_action = vina_menu.addAction("Run docking")
    window.addDockWidget(
        LeftDockWidgetArea,
        run_docking_widget
    )
    run_docking_widget.hide()
    @run_docking_action.triggered.connect
    def toggle():
        run_docking_widget.show()
    
    
    load_results_widget = new_load_results_widget()
    load_results_action = vina_menu.addAction("Load docking")
    window.addDockWidget(
        LeftDockWidgetArea,
        load_results_widget
    )
    load_results_widget.hide()
    @load_results_action.triggered.connect
    def toggle():
        load_results_widget.show()
