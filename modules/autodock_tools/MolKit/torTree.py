#
# Last modified on Mon Apr 15 12:16:09 PDT 2002 by lindy
#
#
# $Header: /opt/cvs/python/packages/share1.5/MolKit/torTree.py,v 1.23.4.3 2008/12/12 20:06:56 rhuey Exp $
#


from mglutil.util.tree import TreeNode
from mglutil.math.transformation import Transformation
import string
from molecule import AtomSet, BondSet

global debug
debug = 0



class TorTree:
    """Nodes are mglutil.util.tree.TreeNode instances. In addition to
    the general attributes defined in TreeNode (parent, children, etc.),
    TorTree nodes have (at least) the following specific attributes:
        number   - the torsion number used to __makeTorsionMap
        bond     - a 2-tuple of atom indices
        atomList - a list of atom indices
        angle    - torsion rotation angle (added by setTorsionAngles)

    Additional attributes are added else where.
    """
    def __init__(self, parser, rootAtom=None):
        if rootAtom is not None:
            self.rootNode, allNodes = self.__buildTree(rootAtom)
            self.torsionMap = self.__orderTorsionMap(allNodes)
        else:
            self.rootNode = self.__makeTree(parser.allLines)
            self.torsionMap = self.__makeTorsionMap()


    def __buildTree(self, rootAtom):
        self.tor_number = 0
        bond = rootAtom.bonds[0]
        at2 = bond.atom1
        if at2==rootAtom:
            at2 = bond.atom2
        # _used is necessary to stop recursion 
        rootAtom.top.allAtoms._used = 0
        rootNode, allNodes = self.__buildNode(bond, rootAtom, at2, 1, [])
        delattr(rootAtom.top.allAtoms, '_used')
        return rootNode, allNodes


    def __buildNode(self, bnd, fromAt, startAt, root=0, allNodes=[]):
        # this is called with an activeTors bnd; 
        # always making a newNode
        # first add the fromAtom
        newNode = TreeNode()
        newNode.number = self.tor_number
        allNodes.append(newNode)
        self.tor_number = self.tor_number + 1
        newNode.atoms_to_move = 0
        #fix this: root need special handling
        if root:
            adjAts = [fromAt]
            fromAt.tt_ind = 0
            #atomIndex counts # of atoms put in torTree, 0-based
            self.atomIndex = 0
            atomList = [0, 1]
            ats = [fromAt, startAt]
        else:
            adjAts = [startAt]
            atomList = []
            #ats = []
        #first expand the adjAts to include all atoms linked
        #to startAt by inactive bonds
        for at in adjAts:
            if at._used: continue
            at._used = 1
            for b in at.bonds:
                # recursively add all rootatoms
                # ie atoms connected to rootatoms by non-rotatable bonds
                at2 = b.atom1
                if at2==at:
                    at2 = b.atom2
                if at2._used: continue
                if not b.activeTors:
                    if not hasattr(at2, 'tt_ind'):
                        #and at2!=startAt and at2!=fromAt:
                        self.atomIndex = self.atomIndex + 1
                        at2.tt_ind = self.atomIndex
                    if at2 not in adjAts:
                        adjAts.append(at2)
                        if at2!=startAt and at2 not in atomList:
                            atomList.append(at2.tt_ind)
                            #ats.append(at2)
        #have to redo this after loop to get breadth first:
        for at in adjAts:
            for b in at.bonds:
                at2 = b.atom1
                if at2==at:
                    at2 = b.atom2
                if at2._used: continue
                if not hasattr(at2, 'tt_ind'):
                    #and at2!=startAt and at2!=fromAt:
                    self.atomIndex = self.atomIndex + 1
                    at2.tt_ind = self.atomIndex
                if b.activeTors:
                    nnode, allNodes = self.__buildNode(b, at, at2, 0, allNodes)
                    newNode.children.append(nnode)
                    #keep track of number of atoms in subtree
                    newNode.atoms_to_move = newNode.atoms_to_move + nnode.atoms_to_move
                    #make sure at and at2 are in this node
                    ##test atoms in this bond
                    for a in [at, at2]:
                        if a!=fromAt and a!=startAt:
                            ##if a.tt_ind not in atomList:
                            atomList.append(a.tt_ind)

        #finally, set atoms to move to include atoms in this node
        d = {}
        for a in atomList:
            d[a] = 0
        atList = d.keys()
        atList.sort()
        #newNode.atomList = atomList
        if root:
            newNode.bond = (None, None)
        else:
            newNode.bond = (fromAt.tt_ind, startAt.tt_ind)
        newNode.atoms_to_move = newNode.atoms_to_move + len(atList)
        newNode.atomList = atList
        #newNode.ats = ats
        return newNode, allNodes


    def __orderTorsionMap(self,allNodes):
        """sorts allNodes, a list of TreeNodes in 'torsion order' so that
        torsionMap[tor_number] refers to the corresponding TreeNode.
        """
        def __sortTorsionMap(node1, node2):
            if node1.atoms_to_move < node2.atoms_to_move:
                return -1
            elif node1.atoms_to_move > node2.atoms_to_move:
                return 1
            # equal atoms_to_move, so use tor_number as sort criteria
            elif node1.number < node2.number:
                return -1
            elif node1.number > node2.number:
                return 1
            else:
                # now there's a problem
                raise RuntimeError, "indistinguishable torsion TreeNodes"
                return 0
        allNodes.sort(__sortTorsionMap)
        #don't put rootNode into TorsionMap!!!
        return allNodes[:-1]


    def __makeTree(self, lineList, flexRes=False):
        # initialize 
        nodeStack = []
        atomToParentNode = 0 # first atom after BRANCH goes to parent
        tor_number = 1
        atomIndex = 0
        # process lines/build tree
        for lineStr in lineList:
            if debug: print lineStr
            wordList = string.split(lineStr)
            if not wordList: continue # skip the loop
            #
            # Here lies the main switch for the PDBQ tags
            #
            if wordList[0] == 'HETATM' or wordList[0] == 'ATOM':
                # pdb is one-based; we are zero-based
                if atomToParentNode or flexRes:
                    # The first atom after the BRANCH goes to the parent
                    nodeStack[-1].parent.atomList.append(atomIndex)
                    atomToParentNode = None # unset; set in BRANCH (below)
                    if debug: print "add atom (parent): ", atomIndex, nodeStack[-1].parent
                    flexRes = False
                else:
                    nodeStack[-1].atomList.append(atomIndex)
                    if debug: print "add atom: ", atomIndex, nodeStack[-1]
                atomIndex = atomIndex + 1
            elif (wordList[0] == 'TORS' or wordList[0] == 'BRANCH'):
                atomToParentNode = 1 # set; unset in HETATM (above)
                newNode = TreeNode(parent=nodeStack[-1])
                newNode.number = tor_number
                newNode.bond = (int(wordList[1])-1, int(wordList[2])-1)
                newNode.atomList = []
                tor_number = tor_number + 1;
                nodeStack.append(newNode)
                if debug: print "push node: ", newNode
            elif (wordList[0] == 'ENDTORS' or wordList[0] == 'ENDBRANCH'):
                nodeStack.pop()
            elif wordList[0] == 'ROOT':
                rootNode = TreeNode()
                rootNode.number = 0
                rootNode.bond = (None, None)
                rootNode.atomList = []
                nodeStack.append(rootNode)
                if debug: print "push root: ", rootNode
            elif wordList[0] == 'ENDROOT':
                pass
            else: # ignore it
                pass
        return rootNode


    def __makeTorsionMap(self):
        """Return list of TreeNodes in 'torsion order' so that
        torsionMap[tor_number] refers to the coresponding TreeNode.
        
        Torsions are specified in order of the number of atoms
        to move with lowest first. If two torsions move the same
        number of atoms, then the one with the lower tor_number
        goes first (the one that appears first in the pdbq file).
        The number of atoms to move total number of atoms refered
        to by a node and all of its children.
        """
        global torsionMap
        torsionMap = []

        def __count_atoms(node):
            atom_count = 0
            # sum your children's (if any) atoms
            for child in node.children:
                atom_count = atom_count + child.atoms_to_move
            # add your own atoms
            node.atoms_to_move  = atom_count + len(node.atomList)
            # add this node to the torsionMap
            torsionMap.append(node)
            
        # traverse the torTree counting atoms and building torsionMap
        self.rootNode.post_traverse(__count_atoms, self.rootNode)
        torsionMap.pop() # remove the rootNode from the torsionMap
        
        def __sortTorsionMap(node1, node2):
            if node1.atoms_to_move < node2.atoms_to_move:
                return -1
            elif node1.atoms_to_move > node2.atoms_to_move:
                return 1
            # equal atoms_to_move, so use tor_number as sort criteria
            elif node1.number < node2.number:
                return -1
            elif node1.number > node2.number:
                return 1
            else:
                # now there's a problem
                raise RuntimeError, "indistinguishable torsion TreeNodes"
                return 0
        torsionMap.sort(__sortTorsionMap)
        return torsionMap


    def getTorsionAngles(self):
        """Return the list of torsion angles"""
        torList = []
        for node in self.torsionMap:
            torList.append(node.angle)
        return torList


    def setTorsionAngles(self, angList):
        """Set the torsion angles for the tree.

        This method does not change atom positions"""
        if len(angList) != len(self.torsionMap):
            raise ValueError, "invalid torsion angle list: ", angList
        # @@ should use zip here
        for angle, node in map(None, angList, self.torsionMap):
            node.angle = angle


    def addTorsion(self, atom1, atom2, angle=0.0):
        """Make the bond between atom1 and atom2 rotatable.

        atom1 and atom2 are indeces into mol.allAtoms
        """
        raise NotImplementedError


    def removeTorsion(self, torsion):
        """How should the torsion be specified?
        """
        raise NotImplementedError


    def __printNode(self, node):
        print 'atomList:', node.atomList
        print 'has ', len(node.children),'children\n'
        for c in node.children:
            print 'printing ', c.number, '  child of ', node.number
            self.__printNode(c)


    def printTree(self):
        if not self.rootNode:
            print 'no rootNode'
            return
        print 'printing rootNode '
        self.__printNode(self.rootNode)
                

    def get_bond_from_indicies(self, atoms, indicies):
        d = {}
        for ind in indicies:
            d[atoms[ind]] = 1
        keys = d.keys()
        bnds = atoms.bonds[0].get(lambda x: x.atom1 in keys and x.atom2 in keys)
        assert len(bnds)==1
        return bnds[0]


    def __get_rotatables(self, node, atoms, rotatables, rootNode=False):               
        if node != self.rootNode:
            rotatables.append(self.get_bond_from_indicies(atoms, node.bond))
        for c in node.children:
            self.__get_rotatables(c, atoms, rotatables)
        return rotatables 
        

    def get_rotatable_bonds(self, mol):
        assert mol.allAtoms.bonds[0]
        mol.allAtoms.bonds[0].activeTors = False
        rotatables = self.__get_rotatables(self.rootNode, mol.allAtoms, BondSet(), rootNode=True) 
        return rotatables


    def get_leaf_atoms(self, mol):
        atom_nums = self.get_leaves(self.rootNode, []) 
        atoms = AtomSet()
        if len(atom_nums):
            atoms = mol.allAtoms.get(str(atom_nums[0]))
            for atnum in atom_nums[1:]:
                atoms.append(mol.allAtoms[atnum])
        return atoms


    def get_leaves(self, node, leaves):
        for c in node.children:
            if len(c.children)==0:
                leaves.extend(c.atomList)
                print "added ", c.atomList
            else:
                self.get_leaves(c, leaves)
        return leaves


    def get_depth(self):
        self.rootNode.depth = 0
        return self._depth(self.rootNode)
        

    def _depth(self, node):
        if not hasattr(node, 'depth'):
            node.depth = 1
        num = 0
        for child in node.children:
            new_num =  self._depth(child)
            if new_num>num:
                num = new_num
            #print 'child.bond=', child.bond, ' new_num=', new_num, " num=", num

        node.depth = node.depth + num
        return node.depth


    def printXmlTree(self, allAtoms, index=1, selStr=None):
        """This function is used to generate XML file for FlexTree package"""
        if not self.rootNode: return
        if selStr is None:
            selStr = allAtoms[0].top.name + "::"
        ostr = '<?xml version="1.0" ?>\n'
        ostr = ostr + '\t<root\n\t\tname="Ligand"\n\t\tid="%d"\n\t\tselectionString="%s"\n\t\tconvolve="FTConvolveApplyMatrixToCoords"\n\t\t'%(99,selStr)
        ostr = ostr + 'file="%s">\n\t'%allAtoms[0].top.parser.filename
        ostr = ostr + '\t<node\n\t\t\tname="Core Ligand"\n\t\t\tid="'+ str(index) +'"\n\t\t\t'
        ostr = ostr + 'shapeParams="cutoff: float 1.85"\n\t\t\t'
        sub_ats=allAtoms.get(lambda x: x.number-1 in self.rootNode.atomList)
        selString = sub_ats.full_name()
        ostr = ostr + 'selectionString="%s"\n\t\t\t'%(selString)
        ostr = ostr + 'shape="FTLines"\n\t\t\t'
        ostr = ostr + 'convolve="FTConvolveApplyMatrixToCoords"\n\t\t\t'
        ostr = ostr + '>\n\t\t</node>\n'
        next_index = index + 1
        for c in self.rootNode.children:
            ost, next_index = self.__printXmlNode(c, next_index, index, allAtoms)
            ostr = ostr + ost
        ostr = ostr + "</root>\n\n"
        return ostr


    def __printXmlNode(self, node, next_index, refNode, allAtoms):
        ostr = '\t\t<node\n\t\t\tname="sidechain%d"\n\t\t\tid="%d"\n\t\t\trefNode="%d"\n\t\t\t'%(node.number, next_index, refNode)
        this_nodes_index = next_index
        next_index += 1
        ostr = ostr + 'shapeParams= "cutoff: float 1.85"\n\t\t\t'
        at1 = allAtoms.get(lambda x: x.number-1==node.bond[0])[0]
        at2 = allAtoms.get(lambda x: x.number-1==node.bond[1])[0]
        atmList = node.atomList[:]
        sub_ats=allAtoms.get(lambda x: x.number-1 in atmList)
        ##IS THIS CORRECT??
        ##sub_ats.insert(0, at2)
        #print "len(sub_ats)=", len(sub_ats), " for node number ", node.number
        selectionString = sub_ats.full_name()
        ostr = ostr + 'selectionString="%s"\n\t\t\t'%(selectionString)
        ostr = ostr + 'motion="FTMotion_RotationAboutAxis"\n\t\t\t'
        ostr = ostr + 'shape = "FTLines"\n\t\t\t'
        ostr = ostr + 'convolve="FTConvolveApplyMatrixToCoords"\n\t\t\t'
        mPs = '"'
        ats = [at1, at2]
        for i in [0,1]:
            at = ats[i]
            mPs = mPs + "point%d: list float %f %f %f, "%(i+1, at.coords[0], at.coords[1], at.coords[2])
        mPs = mPs + ' percent: float 1.0, angle: float 0.0, name: str rotatableBond">'
        ostr = ostr + 'motionParams=%s"\n\t\t'%mPs
        ostr = ostr + "</node>\n\n"
        for c in node.children:
            ost, next_index =  self.__printXmlNode(c, next_index, this_nodes_index, allAtoms)
            ostr = ostr + ost
        return ostr, next_index


    def torTree2dot(self, allAtoms, index=1, selStr=None, labelEdges=True,edgeLabelStyle='node numbers', size="8,6",verbose=False):
        """return (dot format) graph specification"""
        if not self.rootNode: return
        offset = allAtoms[0].number-1
        self.labelEdges=labelEdges
        assert edgeLabelStyle in ['node numbers', 'node bond']
        self.edgeLabelStyle = edgeLabelStyle
        gname =  '"' + allAtoms[0].top.name + '"'
        found_pydot=1
        try:
            import pydot
        except:
            found_pydot=0

        if found_pydot:
            dg = pydot.Graph(graph_name=gname, type='digraph',label=gname, size=size)

        rootID =  str(index)
        if verbose: print "1: set rootID to ", rootID
        #change
        atList = self.rootNode.atomList[:]
        # start with c25,c27,c26,n1
        # remove c27,c26, n1 because they're in rotatable bond to c25
        for c in self.rootNode.children:
            next = c.bond[1] - offset
            if next in atList:   #remove atoms connected by rotatable bonds
                atList.remove(next)
                if verbose: print "removed ", next, " from root"
        if verbose: print "atList =", atList
        #sub_ats = AtomSet()
        #for ind in atList:
        #    sub_ats.append(allAtoms[ind-1])
        sub_ats = allAtoms.get(lambda x: x.number-(1+offset) in atList)
        if verbose: print "sub_ats=", sub_ats.full_name()

        #sub_ats=allAtoms.get(lambda x: x.number-1 in self.rootNode.atomList)
        rootLbl = '"' 
        for at in sub_ats:
            rootLbl+="%s,"%at.name
        #remove last ','
        rootLbl = rootLbl[:-1] + '"'
        #rootLbl = sub_ats.full_name()
        #cpos = rootLbl.rfind(':')
        #rootLbl = '"' + rootLbl[cpos+1:] + '"'
        if found_pydot:
            rootNd = pydot.Node(rootID,label=rootLbl,shape="trapezium",style="bold")
            dg.add_node(rootNd)
        else:
            if verbose: print "1: added node %s, label=%s" %(rootID, rootLbl)
            print "would add pydot.Node(%s, label =%s)" %(rootID, rootLbl)
            dg = None

        next_index = index + 1
        for c in self.rootNode.children:
            if verbose: print c.bond, "call self.__torTree2dot(c,%d, %s, %d,dg,[])" %(next_index, rootID, len(allAtoms))
            next_index = self.__torTree2dot(c, next_index, rootID, allAtoms,dg, [], verbose)

        dotstr = "no pydot"
        if found_pydot:
            dotstr = dg.to_string()
        
        return dotstr


    def __torTree2dot(self, ttnode, next_index, parentID, allAtoms, dotGraph, atList, verbose):
        if verbose: print "__tT2d: ttnode.bond=", ttnode.bond, ' next_index=', next_index, 'parentID=', parentID,' atList=', atList
        ndID = str(next_index)
        ndIndex = ttnode.bond[1]
        ndName = allAtoms.get(lambda x: x.number==ttnode.bond[1]+1)[0].name
        ndLbl = '"%s,'%ndName
        if verbose: print "first: ndLbl to ", ndLbl
        #nd = pydot.Node(ndID,label=rootLbl,shape="trapezium",style="bold")
        atmList = ttnode.atomList[:]
        if verbose: print "before: atmList=", atmList
        if ttnode.bond[1] in atmList:
            atmList.remove(ttnode.bond[1])
        if verbose: print "after: atmList=", atmList
        offset = allAtoms[0].number-1
        if verbose: print "using offset=", offset
        for c in ttnode.children:
            if verbose: print "c.bond=", c.bond
            next = c.bond[1]
            if verbose: print "next =", next
            if next-offset in atmList:
                if verbose: print "removing ", next+offset
                index = atmList.index(next-offset)
                if verbose: print 'cutting atmList at ', index
                atmList=atmList[:index]
        #add names of atoms rigidly bonded to ndID
        if verbose: print "finally: atmList=", atmList
        sub_ats=allAtoms.get(lambda x: x.number-(offset+1) in atmList)
        if len(sub_ats) and verbose:
            print "sub_ats=", sub_ats.name, ' w/number ', sub_ats.number
        for i in sub_ats:
            ndLbl += "%s,"%i.name
        if verbose: print "after sub_ats: ndLbl = ", ndLbl 
        ndLbl = ndLbl[:-1] + '"'
        if verbose: print "after cleanup: ndLbl = ", ndLbl
        found_pydot=1
        try:
            import pydot
        except:
            found_pydot=0
        if found_pydot and dotGraph is not None:
            dnode = pydot.Node(ndID,label=ndLbl)
            dotGraph.add_node(dnode)
            bnd0,bnd1 = ttnode.bond
            edgeLbl = "%s-%s" %(parentID,ndID)
            if self.labelEdges:
                if self.edgeLabelStyle=='node numbers':
                    edgeLbl = '"%s-%s"'%(parentID,ndID)
                else:
                    edgeLbl = '"(%d-%d)"' %(bnd0,bnd1)
                edge = pydot.Edge(parentID,ndID, label=edgeLbl)
            else:
                edge = pydot.Edge(parentID,ndID)
            dotGraph.add_edge(edge)
        else:
            print "would add pydot.Node(", ndID, ",label=", ndLbl, ")"
            print "would add pydot.Edge(",parentID, ',', ndID, ")"
        
        currLbl = str(next_index)
        next_index += 1

        atList = ttnode.atomList[:]
        for c in ttnode.children:
            next = c.bond[1]
            if next in atList:
                if verbose: 
                    print "removing ", next
                atList.remove(next)
        if verbose:
            print "END: atList=", atList
        for c in ttnode.children:
            if verbose:
                print "calling __torTree2dot with c.bond=%d,%d and next_index=%d, currLbl=%s" %(c.bond[0], c.bond[1],next_index, currLbl)
            next_index =  self.__torTree2dot(c, next_index, currLbl, allAtoms, dotGraph, atList, verbose)
        return next_index



if __name__ == '__main__':
    import getopt
    import sys
    from MolKit.pdbParser import PdbqParser

    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'f:')
    except getopt.GetoptError, msg:
        print 'torTree.py: %s' % msg
        sys.exit(2)
        
    filename = None
    for o, a in opt_list:
        if o in ('-f', '--f'):
            filename = a

    if filename:
        parser = PdbqParser(filename)
        mol = parser.parse()
        
        # make torTree and print
        tt = TorTree(parser)

        # run as python2.0 -i torTree.py -f <pdbq_file>
        # use interpreter to examine mol and tt

