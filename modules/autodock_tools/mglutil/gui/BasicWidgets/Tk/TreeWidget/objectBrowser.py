#########################################################################
#
# Date: Jun 2004  Author: Daniel Stoffler
#
#    stoffler@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Daniel Stoffler and TSRI
#
#########################################################################

import Tkinter
import types
import inspect

from mglutil.gui.BasicWidgets.Tk.TreeWidget.tree import TreeView

class ObjectBrowser:
    """Introspect an object with a tree-browser GUI.
Click on the [+] or double-click on the folder icons to expand 1 level.
At the top of this GUI 3 Radiobuttons are visible. Default is 'Attributes'
which means we do not include callable methods. 'Methods' will only display
callable methods. 'All' will display everything.

Required constructor argument: the objetc to be browsed.
Optional constructor arguments:
    rootName: a string used as the name for the root node. Default: 'root'
    title: a string used as window title. Default: 'Python Introspector'
    refresh: this has to be either None or a callable method that returns
             an object to be introspected. If refresh is not None a new
             button is added to the left of the Dismiss button. Clicking this
             button calls the callable method (which must return a new object)
             and rebuilds the tree. If refresh is None, we do not add the
             'Refresh' button. Default: None

API: show() displays the GUI
     hide() hides the GUI (this does not destroy the window)
     rebuild() destroys the current tree and rebuilds it

To get an idea of how this thing works, start Python, import this class, then
type:
browser = ObjectBrowser( __builtins__ )
"""


    def __init__(self, object, rootName='root', title='Python Introspector',
                 refresh=None):

        self.object = object       # the object to be introspeced
        self.rootName = rootName   # name the root node
        self.title = title         # title of the GUI

        assert refresh is None or callable(refresh)
        self.refresh = refresh # used to rebuild tree with new data
        self.busyRebuilding = False

        self.buildGUI()
        

    def buildGUI(self):
        self.root = Tkinter.Toplevel()
        self.root.title(self.title)
        self.root.protocol('WM_DELETE_WINDOW', self.hide )
        self.topFrame = Tkinter.Frame(self.root, relief='raised',
                           bd=4)
        self.topFrame.pack(fill="x", side='top', )

        choices = ['All', 'Attributes', 'Methods']
        self.choicesVarTk = Tkinter.StringVar()
        self.choicesVarTk.set('Attributes')
        for c in choices:
            b = Tkinter.Radiobutton(
                self.topFrame,
                variable = self.choicesVarTk,
                text=c, value=c, command=self.rebuild)
            b.pack(side='left', expand=1, fill='x')

        self.frame = Tkinter.Frame(self.root)
        self.frame.pack(expand=1, fill="both")#, side='bottom')

        frame = Tkinter.Frame(self.root, bd=3)
        frame.pack(fill="x", side='bottom')

        # add Refresh button if specified
        if self.refresh is not None:
            button1 = Tkinter.Button(frame, text='Refresh',
                                     command=self.refresh_cb)
            button1.pack(expand=1, fill="x", side='left')

            button2 = Tkinter.Button(frame, text='Dismiss', command=self.hide)
            button2.pack(expand=1, fill="x", side='left')

        else:
            button2 = Tkinter.Button(frame, text='Dismiss', command=self.hide)
            button2.pack(expand=1, fill="x")
        
        # add root node        
        self.createTree()
        

    def show(self, event=None):
        """show GUI"""
        self.root.deiconify()
        if self.refresh is not None:
            data = self.refresh()
            if data != self.object:
                self.object = data
                self.rebuild()
                

    def hide(self, event=None):
        """hide GUI, note: this does not destroy the GUI"""
        self.root.withdraw()


    def createTree(self):
        """build a TreeView object, add it to the GUI"""
        self.tree = TreeView(master=self.frame, nohistory=True,
                             displayValue=True)
        # I want white background, so:
        self.tree.canvas.configure(bg='white')

        hasChildren = not self.isLeafNode(self.object)

        node = self.tree.addNode(parent=None, name=self.rootName,
                                 object=self.object, hasChildren=hasChildren,
                                 firstExpand_cb=self.expandTreeNode_cb)

    
    def expandTreeNode_cb(self, node=None, object=None):
        """expand the given object by 1 level (i.e. all its children)"""

        if type(object) in [types.ListType, types.TupleType]:
            children = []
            i = 0
            for o in object:
                children.append( (str(i), o) )
                i = i + 1
                
        elif type(object) == types.DictType:
            children = []
            for k, v in object.items():
                children.append( (k, v) )

        elif type(object) in [ types.BooleanType, types.FloatType,
                               types.IntType, types.LongType, types.NoneType,
                               types.StringType ]:
            children = []

        else:
            if self.choicesVarTk.get() == 'All':
                children = inspect.getmembers(object)
            elif self.choicesVarTk.get() == 'Attributes':
                children = inspect.getmembers(object,
                       lambda x: type(x) is not types.MethodType)
            elif self.choicesVarTk.get() == 'Methods':
                if node.parent is None: # only at the first level do we
                    # distinguish between methods, attributes, etc
                    children = inspect.getmembers(object, inspect.ismethod)
                else: #second and deeper levels: here we want to see everything
                    children = inspect.getmembers(object)

        from time import time
        t1 = time()
        nameList=[]
        objList=[]
        hasChildrenList=[]      
        for (name, data) in children:
            hasChildren = not self.isLeafNode(data)
            hasChildrenList.append(hasChildren)
            nameList.append(str(name))
            objList.append(data)
            
        n = self.tree.addNodeSet(parent=node, name=nameList, object=objList,
                                 hasChildren=hasChildrenList,
                                 firstExpand_cb=self.expandTreeNode_cb)
            
#            self.tree.canvas.update()
        #print "firstExpand_cb :", time()-t1


    def isLeafNode(self, x):
        """Helper method: returns True if this object does not have children
        (e.g., int, float, boolean, etc)"""
        
        if type(x) in [
            types.BooleanType, types.FloatType, types.IntType,
            types.LongType, types.NoneType, types.StringType]:
            return True
        elif type(x) in [types.DictionaryType, types.ListType] and \
             len(x)==0:
            return True
        else:
            return False
                       

    def rebuild(self, event=None):
        """destroy old tree, delete frame, create new tree and add it to the
        GUI"""
        self.tree.destroy()
        self.tree = None
        self.createTree()
        self.busyRebuilding = False


    def refresh_cb(self, event=None):
        """rebuild tree with new data (calls rebuild)"""
        
        # Note: on slow machines, if one repeatetly clicks the 'Refresh'
        # button very, very fast, we can reach a state where the canvas has
        # been destroyed while we try to add an icon on it
        if self.busyRebuilding is True:
            return

        if self.refresh is not None:
            self.busyRebuilding = True
            # call the callback to get a new object to introspect
            self.object = self.refresh()
            # and rebuild tree
            self.rebuild()



if __name__ == "__main__":
    # define a callback method
    def foo():
        return __builtins__

    # instanciate the browser, passing the Python moduel __builtins__ to be
    # browsed, and a callback function foo that will add the 'Refresh'
    # button
    browser = ObjectBrowser( __builtins__, refresh=foo )



