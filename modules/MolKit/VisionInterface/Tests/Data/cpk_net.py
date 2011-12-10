from traceback import print_exc

#### Network: cpk ####
#### File written by Vision ####

## loading libraries ##
from MolKit.VisionInterface.MolKitNodes import molkitlib
masterNet.editor.addLibraryInstance(molkitlib,"MolKit.VisionInterface.MolKitNodes", "molkitlib")

from DejaVu.VisionInterface.DejaVuNodes import vizlib
masterNet.editor.addLibraryInstance(vizlib,"DejaVu.VisionInterface.DejaVuNodes", "vizlib")

try:

    ## saving node Read Molecule ##
    from MolKit.VisionInterface.MolKitNodes import ReadMolecule
    node0 = ReadMolecule(constrkw = {}, name='Read Molecule', library=molkitlib)
    masterNet.addNode(node0,94,50)
    widget = node0.inputPorts[0].widget
    widget.set("cv.pdb",0)
except:
    print "WARNING: failed to restore node ReadMolecule called Read Molecule in network masterNet"
    print_exc()
    node0=None
try:

    ## saving node Assign Radii ##
    from MolKit.VisionInterface.MolKitNodes import AssignRadii
    node1 = AssignRadii(constrkw = {}, name='Assign Radii', library=molkitlib)
    masterNet.addNode(node1,86,114)
    widget = node1.inputPorts[1].widget
    widget.set(0,0)
except:
    print "WARNING: failed to restore node AssignRadii called Assign Radii in network masterNet"
    print_exc()
    node1=None
try:

    ## saving node Select Nodes ##
    from MolKit.VisionInterface.MolKitNodes import NodeSelector
    node2 = NodeSelector(constrkw = {}, name='Select Nodes', library=molkitlib)
    masterNet.addNode(node2,89,177)
    widget = node2.inputPorts[1].widget
    widget.set("Atom",0)
    widget = node2.inputPorts[2].widget
    widget.set("",0)
except:
    print "WARNING: failed to restore node NodeSelector called Select Nodes in network masterNet"
    print_exc()
    node2=None
try:

    ## saving node CPK ##
    from MolKit.VisionInterface.MolKitNodes import AtomsAsCPK
    node3 = AtomsAsCPK(constrkw = {}, name='CPK', library=molkitlib)
    masterNet.addNode(node3,85,250)
    widget = node3.inputPorts[4].widget
    widget.set(8,0)
except:
    print "WARNING: failed to restore node AtomsAsCPK called CPK in network masterNet"
    print_exc()
    node3=None
try:

    ## saving node Viewer ##
    from DejaVu.VisionInterface.DejaVuNodes import Viewer
    node4 = Viewer(constrkw = {}, name='Viewer', library=vizlib)
    masterNet.addNode(node4,451,226)
except:
    print "WARNING: failed to restore node Viewer called Viewer in network masterNet"
    print_exc()
    node4=None
try:

    ## saving node Extract Atom Property ##
    from MolKit.VisionInterface.MolKitNodes import AtomsProperty
    node5 = AtomsProperty(constrkw = {}, name='Extract Atom Property', library=vizlib)
    masterNet.addNode(node5,319,52)
    widget = node5.inputPorts[1].widget
    widget.set("number",0)
except:
    print "WARNING: failed to restore node AtomsProperty called Extract Atom Property in network masterNet"
    print_exc()
    node5=None
try:

    ## saving node Color ##
    from DejaVu.VisionInterface.DejaVuNodes import ColorByRamp
    node6 = ColorByRamp(constrkw = {}, name='Color', library=vizlib)
    masterNet.addNode(node6,289,133)
except:
    print "WARNING: failed to restore node ColorByRamp called Color in network masterNet"
    print_exc()
    node6=None

## saving connections for network cpk ##
if node0 is not None and node1 is not None:
    masterNet.connectNodes(node0, node1, "MolSets", "molecules", blocking=True)
if node1 is not None and node2 is not None:
    masterNet.connectNodes(node1, node2, "molecules", "nodes", blocking=True)
if node2 is not None and node3 is not None:
    masterNet.connectNodes(node2, node3, "nodes", "atoms", blocking=True)
if node3 is not None and node4 is not None:
    masterNet.connectNodes(node3, node4, "CPK", "geometries", blocking=True)
if node2 is not None and node5 is not None:
    masterNet.connectNodes(node2, node5, "nodes", "atoms", blocking=True)
if node5 is not None and node6 is not None:
    masterNet.connectNodes(node5, node6, "propertyValues", "values", blocking=True)
if node6 is not None and node3 is not None:
    masterNet.connectNodes(node6, node3, "colors", "colors", blocking=True)
