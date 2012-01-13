## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#
# Last modified on Tue Apr 23 09:20:22 PDT 2002 by lindy
#
# $Header: /opt/cvs/python/packages/share1.5/mglutil/math/statetocoords.py,v 1.11 2008/09/02 22:01:13 gillet Exp $
#

"""statetocoords.py - state to coordinates

The StateToCoords class inherits from Kinematics and Ncoords.
The StateToCoords class handles transformations that apply to
the rootNode of the torTree, changing the coordinates in world
space. The Kinematics class handles those transformations that
apply to the internal nodes to the torTree (ie. torsions) changing
the coordinates in the molecules local coordinate system.
"""
import numpy.oldnumeric as Numeric, math
from mglutil.math.transformation import Transformation
from mglutil.math.kinematics import Kinematics


class StateToCoords(Kinematics):
    def __init__(self, mol, origin, confIndex):
        Kinematics.__init__(self, mol.allAtoms.coords, mol.torTree, tolist=1)
        # this stoc object will leave always deposite it's coords
        # in the given confIndex slot.
        self.confIndex = confIndex
        
        mol.allAtoms.setConformation(confIndex)
        
        def __prepareNode(node, allAtoms, o):
            """Supply each node with atomSet, coords, and atomRange,
            Pre-compute and save the torsionUnitVector, 
            Transform the coords to their local space by subtracting
            the origin
            """
            atomSet = []
            coords = []
            for i in node.atomList:
                atom = allAtoms[i]
                atomSet.append(atom)
                # start with the original coordinates
                c = atom.coords
                # subract the origin
                coords.append((c[0]-o[0], c[1]-o[1], c[2]-o[2], 1.0))
            node.atomSet = atomSet
            node.coords = coords
            node.atomRange = range(len(atomSet))
            if node.bond[0] != None: # skip the root node
                node.a = allAtoms[node.bond[0]]
                node.b = allAtoms[node.bond[1]]

        # add atomSets to each node
        root = mol.torTree.rootNode
        root.pre_traverse(__prepareNode, root, mol.allAtoms, origin)


    def applyState(self, state):
        """
        """
        q = state.quaternion
        t = Numeric.array(state.translation)
        o = Numeric.array(state.origin)

        # construct rootNode transformation matrix
        #mtx = Transformation(t+o, q).getMatrix(transpose=1)
        # Corrected by AG 08/28/2008
        mtx = Transformation(t, q).getMatrix().transpose()

        # apply the torsions
        self.applyAngList(state.torsions, mtx)


    def applyStateOld(self, state):
        """
        """
        q = state.quaternion
        t = Numeric.array(state.translation)
        o = Numeric.array(state.origin)

        # center the coordinates
##          self.resultCoords = (self.resultCoords -
##                               Numeric.array([o[0], o[1], o[2], 0.0]))

        # center the coordinates (node-by-node)
        def __center(node, o): node.coords = node.coords - o
        root = self.torTree.rootNode
        root.pre_traverse(__center, root, Numeric.array([o[0], o[1], o[2], 0.0]))

        # construct rootNode transformation matrix
        mtx = Transformation(t+o, q).getMatrix(transpose=1)

        # apply the torsions
        coords = self.applyAngList(state.torsions, mtx)

        # must "reset" each nodes coords
        def __uncenter(node, o): node.coords = node.coords + o
        root.pre_traverse(__uncenter, root, Numeric.array([o[0], o[1], o[2], 0.0]))

        return coords
    

    def applyOrientation(self, q=(0.,0.,0.,0.), t=(0.,0.,0.), o=(0.,0.,0.)):
        """origin specifies where the local origin is in world coordinates
        (i.e., where is this object's origin in the world)
        """
        # center the coordinates
        self.resultCoords = (self.resultCoords -
                             Numeric.array([o[0], o[1], o[2], 0.0]))
        sum = Numeric.array(t) + Numeric.array(o)
        self.resultCoords = Transformation(sum, q).apply(self.resultCoords)
        return self.getResultCoords()


    def applyQuaternion(self, q, o=(0.,0.,0.)):
        """Apply the given quaterion.
        """
        # center the coordinates
        self.resultCoords = (self.resultCoords -
                             Numeric.array([o[0], o[1], o[2], 0.0]))
        self.resultCoords = Transformation(o, q).apply(self.resultCoords)
        return self.getResultCoords()


    def applyTranslation(self, t=(0., 0., 0.)):
        """Translate by (x, y, z)
        """
        translation = Numeric.array([t[0], t[1], t[2], 0.0])
        self.resultCoords = self.resultCoords + translation
        return self.getResultCoords()

