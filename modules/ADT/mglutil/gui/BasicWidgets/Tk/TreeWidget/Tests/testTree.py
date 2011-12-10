from time import sleep
from mglutil.gui.BasicWidgets.Tk.TreeWidget.tree import TreeView
import Tkinter

def pause(sleepTime=None):
    if sleepTime is None:
        from mglutil.gui.BasicWidgets.Tk.TreeWidget.Tests import pauseLength as sleepTime
    sleep(sleepTime)

def test_addNode():
    from mglutil.gui.BasicWidgets.Tk.TreeWidget.tree import TreeView
    tv = TreeView()
    # addNode(nodename, parentname = None)
    tv.addNode('protein_1')
    tv.addNode('residue_11',    parent='protein_1')
    tv.addNode('AminoAcid',     parent='protein_1|residue_11')
    tv.addNode('A',             parent='protein_1|residue_11|AminoAcid')
    tv.addNode('H',             parent='protein_1|residue_11|AminoAcid')
    tv.addNode('protein_2')
    tv.addNode('protein_3')
    tv.addNode('residue_21',    parent='protein_2')
    tv.addNode('residue_25',    parent='protein_2')
    tv.addNode('basdfe',        parent='protein_2|residue_21')
    tv.addNode('AminoAcid',     parent='protein_2|residue_21')
    tv.addNode('etc',       parent='protein_1|residue_11')
    tv.addNode('color',     parent='protein_1|residue_11|etc')
    tv.addNode('density',   parent='protein_1|residue_11|etc')
    tv.addNode('residue_12',parent='protein_1')
    tv.addNode('2', parent='protein_2|residue_21')
    tv.addNode('3', parent='protein_2|residue_21')
    tv.addNode('4', parent='protein_2|residue_21')
    tv.addNode('L', parent='protein_2|residue_21|AminoAcid')
    for a in range(10):
        name = 'A' + str(a)
        tv.addNode(name,  parent='protein_2|residue_21|AminoAcid')
    tv.addNode('protein_4')
    tv.addNode('residue_22', parent='protein_2')

    # the addNode function check the node before adding to the parentname
    # if already in the child list, the new node will NOT be added.
    tv.addNode('residue_22', parent='protein_2')

    # addNode returns the handle of the noded added
    node = tv.addNode('S', parent='protein_2|residue_21|AminoAcid')
    print node.name, "was added to" , node.parentFullname

    tv.destroy()

def test_createTree():
    from mglutil.gui.BasicWidgets.Tk.TreeWidget.tree import TreeView
    tv = TreeView()
    # paint canvas red so we see it disappear when distroyed
    tv.canvas.configure(bg='red')
    Tkinter._default_root.update()
    pause(0.2)
    tv.destroy()
    
def test_createWithMaster():
    from mglutil.gui.BasicWidgets.Tk.TreeWidget.tree import TreeView
    import Tkinter
    master = Tkinter.Toplevel()
    tv = TreeView(master)
    tv.canvas.configure(bg='red')
    tv.canvas.update()
    pause(1)
    tv.destroy()

def test_hideShow():
    from mglutil.gui.BasicWidgets.Tk.TreeWidget.tree import TreeView
    tv = TreeView()
    # paint canvas red so we see it disappear when distroyed
    tv.canvas.configure(bg='red')
    tv.hide()
    tv.canvas.update()
    pause(1)
    tv.show()
    tv.canvas.update()
    pause(1)
    tv.destroy()

def test_deleteNode():
    from mglutil.gui.BasicWidgets.Tk.TreeWidget.tree import TreeView
    tv = TreeView()
    #Add some nodes 
    tv.addNode('protein_1')
    tv.addNode('residue_11',    parent='protein_1')
    tv.addNode('AminoAcid',     parent='protein_1|residue_11')
    tv.addNode('A',             parent='protein_1|residue_11|AminoAcid')
    node = tv.addNode('H',      parent='protein_1|residue_11|AminoAcid')
    tv.addNode('protein_2')
    
    # Now deleteing
    tv.deleteNode_byName('A', parent='protein_1|residue_11|AminoAcid')
    tv.deleteNode_byName('protein_2')
    
    tv.deleteNode(node)
    
    
    # the following returns error message, 
    # since AminoAcid is NOT the child of 'protein_1'
    tv.deleteNode_byName('AminoAcid',     parent='protein_1')
    
    # This should work
    tv.deleteNode_byName('AminoAcid',     parent='protein_1|residue_11')
    
    tv.destroy()

def test_Expand_or_Collaps():
    from mglutil.gui.BasicWidgets.Tk.TreeWidget.tree import TreeView
    tv = TreeView()
    #Add some nodes 
    tv.addNode('protein_1')
    node = tv.addNode('residue_11',    parent='protein_1')
    tv.addNode('AminoAcid',     parent='protein_1|residue_11')
    tv.addNode('A',             parent='protein_1|residue_11|AminoAcid')
    tv.addNode('H',             parent='protein_1|residue_11|AminoAcid')
    tv.addNode('protein_2')
    
    # equivalent to 
    print "Is node", node.name, "expanded?: ",
    print node.expanded
    
    tv.ExpandNode('protein_1|residue_11')
    print "After expanding, is node", node.name, "expanded?: ",
    print node.expanded
    
    tv.CollapseNode('protein_1|residue_11')
    print "After collapsing, is node", node.name, "expanded?: ",
    print node.expanded
    
    tv.destroy()


def test_SelectNode_deSelectNode():
    from mglutil.gui.BasicWidgets.Tk.TreeWidget.tree import TreeView
    tv = TreeView()
    #Add some nodes 
    tv.addNode('protein_1')
    node = tv.addNode('residue_11',    parent='protein_1')
    tv.addNode('AminoAcid',     parent='protein_1|residue_11')
    tv.addNode('A',             parent='protein_1|residue_11|AminoAcid')
    tv.addNode('H',             parent='protein_1|residue_11|AminoAcid')
    tv.addNode('protein_2')

    selection = tv.GetSelected()
    if selection:
        print "Now", selection.name, "is selected"
    else:
        print "Nothing is selected"

    # Now select a node 
    tv.Select("protein_1|residue_11|AminoAcid|A")
    
    selection = tv.GetSelected()
    if selection:
        print "***", selection.name, "is selected ***"
    else:
        print "Nothing is selected"
    
    tv.destroy()


def test_Add_to_History():
    from mglutil.gui.BasicWidgets.Tk.TreeWidget.tree import TreeView
    tv = TreeView()
    tv.addNode('protein_1')
    tv.addNode('residue_11',    parent='protein_1')
    tv.addNode('AminoAcid',     parent='protein_1|residue_11')
    tv.addNode('A',             parent='protein_1|residue_11|AminoAcid')
    tv.addNode('H',             parent='protein_1|residue_11|AminoAcid')
    tv.addNode('protein_2')

    tv.Select("protein_1|residue_11|AminoAcid|A")
    
    selection = tv.GetSelected()
    if selection:
        print "Adding", selection.name, "to the history list"
    else:
        print "Nothing is selected"

    # add the current selected node to the list
    tv.AddToHistoryList()
    
    tv.destroy()

def test_Delete_from_History():
    from mglutil.gui.BasicWidgets.Tk.TreeWidget.tree import TreeView
    tv = TreeView()
    tv.addNode('protein_1')
    tv.addNode('residue_11',    parent='protein_1')
    tv.addNode('AminoAcid',     parent='protein_1|residue_11')
    tv.addNode('A',             parent='protein_1|residue_11|AminoAcid')
    tv.addNode('H',             parent='protein_1|residue_11|AminoAcid')
    tv.addNode('protein_2')

    tv.Select("protein_1|residue_11|AminoAcid|A")
    tv.AddToHistoryList()
    
    tv.Select("protein_2")
    tv.AddToHistoryList()
    
    tv.Select("protein_1|residue_11")
    tv.AddToHistoryList()
    
    tv.Select("protein_1|residue_11|AminoAcid|H")
    tv.AddToHistoryList()
    
    # delete 
    tv.DeleteFromHistoryList(2)
    list=[0,2]
    tv.historyList.Delete_multi(list)
    tv.destroy()

def test_Lock_Unlock_History():
    from mglutil.gui.BasicWidgets.Tk.TreeWidget.tree import TreeView
    tv = TreeView()
    tv.addNode('protein_1')
    tv.addNode('residue_11',    parent='protein_1')
    tv.addNode('AminoAcid',     parent='protein_1|residue_11')
    tv.addNode('A',             parent='protein_1|residue_11|AminoAcid')
    tv.addNode('H',             parent='protein_1|residue_11|AminoAcid')
    tv.addNode('protein_2')

    tv.Select("protein_1|residue_11|AminoAcid|A")
    tv.AddToHistoryList()
    tv.Select("protein_2")
    tv.AddToHistoryList()
    tv.Select("protein_1|residue_11")
    tv.AddToHistoryList()
    tv.Select("protein_1|residue_11|AminoAcid|H")
    tv.AddToHistoryList()
    
    # lock 
    tv.LockHistoryItem(2)
    tv.LockHistoryItem(3)
    # unlock        
    tv.UnLockHistoryItem(2)
    
    tv.destroy()

def foo(item): # called by test_SetAction()
    print item.name

def test_SetAction():
    from mglutil.gui.BasicWidgets.Tk.TreeWidget.tree import TreeView
    tv = TreeView()
    tv.addNode('protein_1')
    tv.addNode('residue_11',    parent='protein_1')
    tv.addNode('AminoAcid',     parent='protein_1|residue_11')

    #Set what to do after an action here.
    tv.setAction(event='select', function=foo)
    tv.Select("protein_1")
    
    tv.destroy()




def test_template():
    from mglutil.gui.BasicWidgets.Tk.TreeWidget.tree import TreeView
    tv = TreeView()
    tv.destroy()

def test_NoHistory():
    """ Test the crateion of a TreeView with no history Pane """


    tv = TreeView(nohistory=True)
    # paint canvas red so we see it disappear when distroyed
    tv.canvas.configure(bg='red')

    tv.addNode('protein_1')
    tv.addNode('residue_11',    parent='protein_1')
    tv.addNode('AminoAcid',     parent='protein_1|residue_11')
    tv.addNode('A',             parent='protein_1|residue_11|AminoAcid')
    tv.addNode('H',             parent='protein_1|residue_11|AminoAcid')
    tv.addNode('protein_2')

    tv.Select("protein_1|residue_11|AminoAcid|A")
    
    selection = tv.GetSelected()
    if selection:
        print "Adding", selection.name, "to the history list"
    else:
        print "Nothing is selected"

    # add the current selected node to the list
    tv.AddToHistoryList()
    
    #tv.topFrame.master.update()
    pause(0.2)
    tv.destroy()    



def test_obj2Node():
    """ Test the creation of a TreeView with obj2Node == false
    In that case there is not  1to 1 relation between object and node
    tv.objToNode should stay empty.
    """

    tv = TreeView(obj2Node=False)
    # paint canvas red so we see it disappear when distroyed
    tv.canvas.configure(bg='red')

    tv.addNode('protein_1')
    tv.addNode('residue_11',    parent='protein_1')
    tv.addNode('AminoAcid',     parent='protein_1|residue_11')
    tv.addNode('A',             parent='protein_1|residue_11|AminoAcid')
    tv.addNode('H',             parent='protein_1|residue_11|AminoAcid')
    tv.addNode('protein_2')

    tv.Select("protein_1|residue_11|AminoAcid|A")
    
    assert tv.objToNode == {}
    
    #tv.topFrame.master.update()
    pause(0.2)
    tv.destroy()    


#"""
def test_copyTree():
    """ Test the function to copy a tree into another one"""

    tv = TreeView()
    tv2 = TreeView()
    # paint canvas red so we see it disappear when distroyed
    tv.canvas.configure(bg='red')

    tv.addNode('protein_1')
    tv.addNode('residue_11',    parent='protein_1')
    tv.addNode('AminoAcid',     parent='protein_1|residue_11')
    tv.addNode('A',             parent='protein_1|residue_11|AminoAcid')
    nodetest = tv.addNode('H',parent='protein_1|residue_11|AminoAcid')
    tv.addNode('protein_2')

    tv.copy(tv2)
    node = tv2.findNodeFromName(nodetest.GetFullName())
    assert node.name == 'H'
    Tkinter._default_root.update()
    pause()
    tv.destroy()    
    tv2.destroy()    

def test_delete2TreesWithSameMaster():
    """ Test the function that if a tree does not own its master it does not destgroy it"""

    tv = TreeView()
    tv2 = TreeView()
    Tkinter._default_root.update()
    pause()
    tv.destroy()
    tv2.destroy()

    
def test_moveNode():
    """ Test the function to move a tree node (including subtree) """

    tv = TreeView()
    # paint canvas red so we see it disappear when distroyed
    tv.canvas.configure(bg='red')
    tv.addNode('protein_1')
    tv.addNode('residue_11',parent='protein_1')
    tv.addNode('AminoAcid',parent='protein_1|residue_11')
    tv.addNode('A',parent='protein_1|residue_11|AminoAcid')
    tv.addNode('H',parent='protein_1|residue_11|AminoAcid')
    tv.addNode('protein_2')
    tv.addNode('protein_333')
    tv.addNode('residue_21',parent='protein_2')
    tv.addNode('residue_Root',parent='protein_2')
    tv.addNode('basdfe',parent='protein_2|residue_21')
    tv.addNode('AminoAcidXXX',parent='protein_2|residue_21')
    tv.addNode('etc',parent='protein_1|residue_11')
    tv.addNode('color',parent='protein_1|residue_11|etc')
    tv.addNode('density',parent='protein_1|residue_11|etc')
    tv.addNode('residue_Moving',parent='protein_1')
    tv.addNode('2',parent='protein_2|residue_21')
    tv.addNode('3',parent='protein_2|residue_21')
    tv.addNode('4',parent='protein_2|residue_21')
    tv.addNode('L',parent='protein_2|residue_21|AminoAcidXXX')
    tv.addNode('protein_33xxxxxxx')

    node=tv.roots[0].children[1]
    dest=tv.roots[1].children[1]
    tv.moveNode(node,dest)

    assert node.parent ==dest
    node.GetFullName()=='protein_2|residue_Root|residue_Moving'
    assert dest.children[0] == node
    #tv.topFrame.master.update()
    pause(0.2)
    tv.destroy()


def test_moveNodeUpOrDown():
    """ Test the function to move a tree node up and down in the node list """

    tv = TreeView()
    # paint canvas red so we see it disappear when distroyed

    tv.addNode('protein_1')
    tv.addNode('residue_11',    parent='protein_1')
    tv.addNode('AminoAcid',     parent='protein_1|residue_11')
    tv.addNode('A',             parent='protein_1|residue_11|AminoAcid')
    tv.addNode('H',             parent='protein_1|residue_11|AminoAcid')
    tv.addNode('protein_2')
    tv.addNode('protein_3')
    tv.addNode('residue_21',    parent='protein_2')
    tv.addNode('residue_25',    parent='protein_2')
    tv.addNode('basdfe',        parent='protein_2|residue_21')
    tv.addNode('AminoAcid',     parent='protein_2|residue_21')
    tv.addNode('etc',       parent='protein_1|residue_11')
    tv.addNode('color',     parent='protein_1|residue_11|etc')
    tv.addNode('density',   parent='protein_1|residue_11|etc')
    tv.addNode('residue_12',parent='protein_1')
    tv.addNode('2', parent='protein_2|residue_21')
    tv.addNode('3', parent='protein_2|residue_21')
    tv.addNode('4', parent='protein_2|residue_21')
    tv.addNode('L', parent='protein_2|residue_21|AminoAcid')
    tv.addNode('protein_33xxxxxxx')

    node=tv.roots[0].children[1]
    dest=tv.roots[1].children[1]

    ## NOTE: self.objToNode can also return node instance.

    n2=tv.findNodeFromName('protein_1|residue_12')
    n1=tv.findNodeFromName('protein_1|residue_11')
    p1=tv.findNodeFromName('protein_1')

    n2.moveUp()
    assert p1.children.index(n2)==0    
    n2.moveDown()
    assert p1.children.index(n1)==0
    

    #tv.topFrame.master.update()
    pause(0.2)
    tv.destroy()



if __name__=='__main__':

    test_createTree()
    test_createWithMaster()
    test_hideShow()
    test_addNode()
    test_deleteNode()
    test_Expand_or_Collaps()
    test_SelectNode_deSelectNode()
    test_Add_to_History()
    test_Delete_from_History()
    test_Lock_Unlock_History()
    test_SetAction()
    test_template()
    test_NoHistory()
    test_obj2Node()
    test_copyTree()
    test_delete2TreesWithSameMaster()
    test_moveNode()
    test_moveNodeUpOrDown()
