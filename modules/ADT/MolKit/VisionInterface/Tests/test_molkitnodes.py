#########################################################################
#
# Date: Aug 2004  Author: Daniel Stoffler
#
#       stoffler@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Daniel Stoffler, and TSRI
#
#########################################################################

import sys, os
from time import sleep
from Vision.VPE import VisualProgramingEnvironment

from MolKit.VisionInterface.MolKitNodes import molkitlib, ReadMolecule, NodeSelector, AssignRadii, BondsByDist, AtomsProperty, AtomsAsMSMS

ed = None

###############################
## implement setUp and tearDown
###############################

def setUp():
    global ed
    ed = VisualProgramingEnvironment("test individual molkit nodes",
                                     withShell=0,
                                     visibleWidth=400, visibleHeight=300)

    ed.root.update()
    ed.configure(withThreads=0)
    ed.addLibraryInstance(
        molkitlib, 'MolKit.VisionInterface.MolKitNodes', 'molkitlib')
    ed.root.update()
    

def tearDown():
    ed.exit_cb()
    import gc
    gc.collect()

##########################
## Helper methods
##########################

def pause(sleepTime=0.1):
    ed.master.update()
    sleep(sleepTime)


############################################################################
## Tests
############################################################################


def test_01_ReadMolecule():
    # test the Read Molecule node
    from MolKit.protein import ProteinSet
    net = ed.currentNetwork
    net.runOnNewData.value = True    
    node1 = ReadMolecule(library=molkitlib)
    net.addNode(node1, 20, 20)
    node1.inputPorts[0].widget.set(os.path.abspath("1crn.pdb"))
    pause()
    assert isinstance(node1.outputPorts[0].data, ProteinSet),\
           "Expected %s, got %s"%(
        ProteinSet, node1.outputPorts[0].data.__class__)

    
def test_02_NodeSelector():
    # test the select 'node' node
    net = ed.currentNetwork
    net.runOnNewData.value = True    
    node1 = ReadMolecule(library=molkitlib)
    net.addNode(node1, 20, 20)
    node1.inputPorts[0].widget.set(os.path.abspath("1crn.pdb"))
    node2 = NodeSelector(library=molkitlib)
    net.addNode(node2, 30, 100)
    net.connectNodes(node1, node2, "MolSets", "nodes")
    node2.toggleNodeExpand_cb()
    pause()
    # because we run(), upon connecting we pass the data and the output of
    # node2 should have data
    # default node output is AtomSet
    from MolKit.molecule import AtomSet
    data = node2.outputPorts[0].data
    assert isinstance(data, AtomSet), "Expected %s, got %s"%(
        AtomSet, data.__class__)
    # switch to ResidueSet
    from MolKit.protein import ResidueSet
    node2.inputPorts[1].widget.set("Residue")
    data = node2.outputPorts[0].data
    assert isinstance(data, ResidueSet), "Expected %s, got %s"%(
        ResidueSet, data.__class__)
    # switch to Chain
    from MolKit.protein import ChainSet
    node2.inputPorts[1].widget.set("Chain")
    data = node2.outputPorts[0].data
    assert isinstance(data, ChainSet), "Expected %s, got %s"%(
        ChainSet, data.__class__)
    # switch to Molecule
    from MolKit.molecule import MoleculeSet
    node2.inputPorts[1].widget.set("Molecule")
    data = node2.outputPorts[0].data
    assert isinstance(data, MoleculeSet), "Expected %s, got %s"%(
        MoleculeSet, data.__class__)
    

def test_03_AssignRadii():
    # test the Assign Radii node
    net = ed.currentNetwork
    net.runOnNewData.value = True    
    node1 = ReadMolecule(library=molkitlib)
    net.addNode(node1, 20, 20)
    node1.inputPorts[0].widget.set(os.path.abspath("1crn.pdb"))
    node2 = AssignRadii(library=molkitlib)
    net.addNode(node2, 30, 100)
    net.connectNodes(node1, node2, "MolSets", "molecules")
    node2.toggleNodeExpand_cb()
    pause()
    # both output ports should have data
    # a molecule:
    data = node2.outputPorts[0].data
    from MolKit.protein import ProteinSet
    assert isinstance(data, ProteinSet), "Expected %s, got %s"%(
        ProteinSet, data.__class__)
    # list of radii:
    data = node2.outputPorts[1].data
    assert len(data) == 327, "Expected 327, got %s"%len(data)
    # check if united Radii works
    node2.inputPorts[1].widget.set(1)
    
    
def test_04_BuildBonds():
    # test the Build Bonds By Distance node
    net = ed.currentNetwork
    net.runOnNewData.value = True    
    node1 = ReadMolecule(library=molkitlib)
    net.addNode(node1, 20, 20)
    node1.inputPorts[0].widget.set(os.path.abspath("1crn.pdb"))
    node2 = BondsByDist(library=molkitlib)
    net.addNode(node2, 140, 100)
    net.connectNodes(node1, node2, "MolSets", "molecules")
    node2.toggleNodeExpand_cb()
    pause()
    # node2 should output a molecule:
    data = node2.outputPorts[0].data
    from MolKit.protein import ProteinSet
    assert isinstance(data, ProteinSet), "Expected %s, got %s"%(
        ProteinSet, data.__class__)

    
def test_05_AtomsProperty():
    # test the Extract Atoms property node
    net = ed.currentNetwork
    net.runOnNewData.value = True    
    node1 = ReadMolecule(library=molkitlib)
    net.addNode(node1, 20, 20)
    node1.inputPorts[0].widget.set(os.path.abspath("1crn.pdb"))
    node2 = NodeSelector(library=molkitlib)
    net.addNode(node2, 30, 100)
    net.connectNodes(node1, node2, "MolSets", "nodes")
    node3 = AtomsProperty(library=molkitlib)
    net.addNode(node3, 90, 170)
    net.connectNodes(node2, node3, "nodes", "atoms")
    node3.toggleNodeExpand_cb()
    pause()
    # define a couple of attributes we want to access
    # FIXME: ATTRIBUTE "radius" is not working!!??
    attrs = ["_uniqIndex", "atomicNumber", "bondOrderRadius",
             "covalentRadius", "element", "name", "number", "organic", 
             "segID", "temperatureFactor", "vdwRadius"]
    for a in attrs:
        node3.inputPorts[1].widget.set(a)
        assert len(node3.outputPorts[0].data) == 327,\
               "Expected 327, got %s"%len(node3.outputPorts[0].data)
    

def test_06_MSMS():
    # test the MSMS node
    net = ed.currentNetwork
    net.runOnNewData.value = True    
    node1 = ReadMolecule(library=molkitlib)
    net.addNode(node1, 20, 20)
    node1.inputPorts[0].widget.set(os.path.abspath("1crn.pdb"))
    # assign radii
    node2 = AssignRadii(library=molkitlib)
    net.addNode(node2, 275, 62)
    net.connectNodes(node1, node2, "MolSets", "molecules")
    node2.toggleNodeExpand_cb()
    # select nodes
    node3 = NodeSelector(library=molkitlib)
    net.addNode(node3, 14, 153)
    net.connectNodes(node2, node3, "molecules", "nodes")
    node3.toggleNodeExpand_cb()
    # msms
    node4 = AtomsAsMSMS(library=molkitlib)
    net.addNode(node4, 263, 267)
    net.connectNodes(node3, node4, "nodes", "atoms")
    net.connectNodes(node2, node4, "radii", "radii")
    pause()
    # check data
    from mslib import MSMS
    assert isinstance(node4.outputPorts[0].data, MSMS),\
           "Expected %s, got %s"%(MSMS, node4.outputPorts[0].data.__class__)
    assert len(node4.outputPorts[0].data.coords) == 327,\
           "Expected 327, got %s"%len(node4.outputPorts[0].data.coords)
    # just for the fun of it, click the united radii checkbutton a couple times
    node2.inputPorts[1].widget.set(1)
    node2.inputPorts[1].widget.set(0)
    node2.inputPorts[1].widget.set(1)
    node2.inputPorts[1].widget.set(0)
 

def test_07_MSMSWithAtomSubset():
    # test if we can pass an atom set (C,CA,N) to compute an MSMS
    masterNet = ed.currentNetwork
    ## saving node Read Molecule ##
    node0 = ReadMolecule(library=molkitlib)
    masterNet.addNode(node0,20,20)
    node0.inputPorts[0].widget.set(os.path.abspath("1crn.pdb"),0) # do not run yet
    ## saving node MSMS ##
    node2 = AtomsAsMSMS(library=molkitlib)
    masterNet.addNode(node2,250,275)
    ## saving node Select Nodes ##
    node4 = NodeSelector(library=molkitlib)
    masterNet.addNode(node4,14,153)
    node4.inputPorts[2].widget.set("C,CA,N")
    apply(node4.configure, (), {'expanded': True})
    ## saving node Assign Radii ##
    node5 = AssignRadii(library=molkitlib)
    masterNet.addNode(node5,275,62)
    apply(node5.configure, (), {'expanded': True})
    ## saving node Extract Atom Property ##
    node7 = AtomsProperty(library=molkitlib)
    masterNet.addNode(node7,302,169)
    node7.inputPorts[1].widget.set("radius")
    apply(node7.configure, (), {'expanded': True})
    ## saving connections for network Network 0 ##
    masterNet.connectNodes(node0, node5, "MolSets", "molecules")
    masterNet.connectNodes(node5, node4, "molecules", "nodes")
    masterNet.connectNodes(node4, node2, "nodes", "atoms")
    masterNet.connectNodes(node4, node7, "nodes", "atoms")
    masterNet.connectNodes(node7, node2, "propertyValues", "radii")
    pause()
    masterNet.run()
    # check if the data is correct
    # (compare to above: the entire AtomSet would be 237)
    assert node2.outputPorts[0].data is not None,\
           "Expected data, got %s"%node2.outputPorts[0].data 
    assert len(node2.outputPorts[0].data.coords) == 138,\
           "Expected 128, got %s"%len(node2.outputPorts[0].data.coords)
