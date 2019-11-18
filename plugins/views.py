"""
    Manage views and scenes.
    
    Pedro Sousa Lacerda <pslacerda@gmail.com>
    LaBiMM / UFBA: Laboratório de Bioinformática e Modelagem Molecular

    Click on the "Manage views" option on the "Scene" menu (the last one).
    Double-click an entry to rename the view. Erase it's name to delete it.
"""

import pymol
import pymol.gui
from pymol import cmd
from pymol import _cmd
from pymol.Qt import QtWidgets, QtCore, QtGui


class ViewManager:
    
    @staticmethod
    def recall(key):
        cmd.view(key, 'recall')

    @staticmethod
    def store(key):
        cmd.view(key, 'store')
    
    @staticmethod
    def rename(old_key, new_key):
        curr_view = cmd.get_view()
        cmd.view(old_key, 'recall', False)
        cmd.view(old_key, 'clear')
        cmd.view(new_key, 'store')
        cmd.set_view(curr_view)
    
    @staticmethod
    def clear(key):
        cmd.view(key, 'clear')
    
    @staticmethod
    def get_keys():
        return list(pymol._view_dict)
    
    @staticmethod
    def next_key():
        count = len(ViewManager.get_keys())
        return f'view{count}'


def new_manager_widget(manager, title):
    dockWidget = QtWidgets.QDockWidget()
    dockWidget.setWindowTitle(title)

    widget = QtWidgets.QWidget()
    layout = QtWidgets.QVBoxLayout(widget)
    dockWidget.setWidget(widget)
    widget.setLayout(layout)
    
    model = QtGui.QStandardItemModel(0, 1, widget)
    
    tree = QtWidgets.QTreeView(widget)
    tree.setModel(model)
    tree.setHeaderHidden(True)
        
    current_key = None
    
    def update_treeview():
        model.clear()
        for key in manager.get_keys():
            item = QtGui.QStandardItem()
            item.setText(key)
            model.blockSignals(True)
            model.appendRow(item)
            model.blockSignals(False)
    
    @model.rowsInserted.connect
    def onRowsInserted(index):
        key = model.itemFromIndex(index).text()
        store_view(key)
        update_treeview()
        
    @model.dataChanged.connect
    def onDataChanged(index):
        new_key = model.itemFromIndex(index).text()
        if new_key == "":
            manager.clear(current_key)
        else:
            manager.rename(current_key, new_key)
        update_treeview()
    
    @dockWidget.visibilityChanged.connect
    def onVisibilityChanged(visible):
        if not visible:
            return
        update_treeview()
    
    @tree.clicked.connect
    def onClicked(index):
        nonlocal current_key

        key = model.itemFromIndex(index).text()
        current_key = key
        
        manager.recall(key)
    
    storeBtn = QtWidgets.QPushButton("Store", widget)
    @storeBtn.clicked.connect
    def onClicked():
        count = len(manager.get_keys())
        manager.store(f"view{count}")
        update_treeview()
    
    layout.addWidget(tree)
    layout.addWidget(storeBtn)
    return dockWidget


def __init_plugin__(app=None):
    window = pymol.gui.get_qtwindow()
    menu = window.menudict['Scene']
    menu.addSeparator()
    
    action = menu.addAction("Manage views")
    widget = new_manager_widget(ViewManager, "Views")
    window.addDockWidget(
        QtCore.Qt.LeftDockWidgetArea,
        widget
    )
    widget.hide()
    @action.triggered.connect
    def toggle():
        widget.show()
