"""
    Set and get views stored on the .pse session.
    
    Pedro Sousa Lacerda <pslacerda@gmail.com>
    LaBiMM / UFBA: Laboratório de Bioinformática e Modelagem Molecular

    It is the last option ("Manage views") on the "Scene" menu. Double-click to
    rename a view. Erase it's name to delete it.

    In order to persist data in the PSE session it stores data into the
    mesh_clear_selection setting. It was the only way I found, if you use such
    option, try other string settings (some are persisted).
"""

import json

import pymol.gui
from pymol import cmd
from pymol.Qt import QtWidgets, QtCore, QtGui


DATA_SETTING_KEY = "mesh_clear_selection"


def new_view_widget():
    dockWidget = QtWidgets.QDockWidget()
    dockWidget.setWindowTitle("Views")

    widget = QtWidgets.QWidget()
    layout = QtWidgets.QVBoxLayout(widget)
    dockWidget.setWidget(widget)
    widget.setLayout(layout)

    viewsModel = QtGui.QStandardItemModel(0, 1, widget)
    def persist():
        records = []
        for i in range(viewsModel.rowCount()):
            item = viewsModel.itemFromIndex(viewsModel.index(i, 0))
            text = item.text()
            view = item.data()
            records.append((text, view))
        cmd.set(DATA_SETTING_KEY, json.dumps(records))

    viewsModel.rowsInserted.connect(persist)
    #viewsModel.rowsMoved.connect(persist)
    viewsModel.rowsRemoved.connect(persist)
    @viewsModel.dataChanged.connect
    def onDataChanged(index):
        item = viewsModel.itemFromIndex(index)
        if item.text() == "":
            viewsModel.removeRow(index.row())
        else:
            persist()

    @dockWidget.visibilityChanged.connect
    def onVisibilityChanged(visible):
        if not visible:
            return
        viewsModel.clear()
        try:
            records = json.loads(cmd.get(DATA_SETTING_KEY))
        except:
            records = []
        
        for (text, view) in records:
            item = QtGui.QStandardItem()
            item.setText(text)
            item.setData(view)
            viewsModel.appendRow(item)
    
    tree = QtWidgets.QTreeView(widget)
    tree.setModel(viewsModel)
    tree.setHeaderHidden(True)
    @tree.clicked.connect
    def onClicked(index):
        view = viewsModel.itemFromIndex(index).data()
        cmd.set_view(view)
    
    storeBtn = QtWidgets.QPushButton("Store", widget)
    @storeBtn.clicked.connect
    def onClicked():
        count = viewsModel.rowCount()
        item = QtGui.QStandardItem()
        item.setText(f"view{count}")
        item.setData(cmd.get_view())
        viewsModel.appendRow(item)
    
    layout.addWidget(tree)
    layout.addWidget(storeBtn)
    return dockWidget


def __init_plugin__(app=None):
    window = pymol.gui.get_qtwindow()
    menu = window.menudict['Scene']
    menu.addSeparator()
    action = menu.addAction("Manage views")
    
    widget = new_view_widget()
    window.addDockWidget(
        QtCore.Qt.LeftDockWidgetArea,
        widget
    )
    widget.hide()
    @action.triggered.connect
    def toggle():
        widget.show()
