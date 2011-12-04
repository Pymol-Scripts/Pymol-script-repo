#############################################################################
#
# Author: Michel F. SANNER, Sophie COON
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

#
#$Header: /opt/cvs/python/packages/share1.5/MolKit/data/__init__.py,v 1.1 2002/01/23 20:18:31 rhuey Exp $
#
#$Id: __init__.py,v 1.1 2002/01/23 20:18:31 rhuey Exp $
#

import string

def Read(filename):
    from MolKit.pdbParser import PdbParser, PdbqParser,PdbqsParser, PQRParser
    from MolKit.mol2Parser import Mol2Parser
    ext = string.split(filename, '.')
    if ext[-1]=='pdb':
        parser = PdbParser(filename)

    elif ext[-1]=='pdbq':
        parser = PdbqParser(filename)
    
    elif ext[-1]=='pdbqs':
        parser = PdbqsParser(filename)

    elif ext[-1]=='pqr':
        parser = PQRParser(filename)

    elif ext[-1]=='mol2':
        parser = Mol2Parser(filename)

    else:
        print "File Format unknown can't parse it"
        return []
    molecules = parser.parse()
    return molecules

def WritePDB(filename,node):
    from MolKit.pdbWriter import PdbWriter
    writer = PdbWriter()
    writer.write(filename, node)


##  def getNodesByMolecule(self, nodes, molecules,nodeType=None):
##      """ moleculeSet, [nodeSet, nodeSet] <- getNodesByMolecule(nodes)
##      nodes can be either: a string, a TreeNode or a TreeNodeSet.
##      This method returns a molecule set and for each molecule a TreeNodeSet
##      of the nodes belonging to this molecule.
##      'nodeType' enables a desired type of nodes to be returned for each
##      molecule
##      """

##      # if it is a string, get a bunch of nodes from the string
##      if type(nodes)==types.StringType:
##          nodes = molecules.NodesFromName(nodes)

##      assert issubclass(nodes.__class__, TreeNode) or \
##             issubclass(nodes.__class__, TreeNodeSet)

##      # if nodes is a single TreeNode make it a singleton TreeNodeSet
##      if issubclass(nodes.__class__, TreeNode):
##          nodes = nodes.setClass([nodes])

##      if len(nodes)==0: return MoleculeSet([]), []

##      # catch the case when nodes is already a MoleculeSet
##      if nodes.elementType in [Molecule, Protein]:
##          molecules = nodes
##      else: # get the set of molecules
##          molecules = nodes.top.uniq()

##      # build the set of nodes for each molecule
##      nodeSets = []

##      # find out the type of the nodes we want to return
##      searchType=0
##      if nodeType is None:
##          Klass = nodes.elementType # class of objects in that set
##      else:
##          assert issubclass(nodeType, TreeNode)
##          Klass = nodeType
##          if Klass != nodes.elementType:
##              searchType=1

##      for mol in molecules:
##          # get set of nodes for this molecule
##          mol_nodes = nodes.get(lambda x, mol=mol: x.top==mol)

##          # get the required types of nodes
##          if searchType:
##              if Klass == Atom and hasattr(mol_nodes, 'allAtoms'):
##                  mol_nodes = mol_nodes.allAtoms
##              else:
##                  mol_nodes = mol_nodes.findType( Klass )

##          nodeSets.append( mol_nodes )

##      return molecules, nodeSets

##  from MolKit.protein import ProteinSet, Protein,ResidueSet, Residue, ChainSet, Chain
##  from MolKit.molecule import AtomSet, Atom, MoleculeSet, Molecule

##  def getNodesByLevel(self, nodes, molecules,levelType=Protein, nodeType=None):
##      """ ProteinSet, [nodeSet, nodeSet] <- getNodesByLevel(nodes)
##      nodes can be either: a string, a TreeNode or a TreeNodeSet.
##      This method returns a molecullevel set and for each level a TreeNodeSet
##      of the nodes belonging to this molecule.
##      'nodeType' enables a desired type of nodes to be returned for each
##      molecule
##      """
##      import types
##      # if it is a string, get a bunch of nodes from the string
##      if type(nodes)==types.StringType:
##          nodes = molecules.NodesFromName(nodes)

##      assert issubclass(nodes.__class__, TreeNode) or \
##             issubclass(nodes.__class__, TreeNodeSet)

##      # if nodes is a single TreeNode make it a singleton TreeNodeSet
##      if issubclass(nodes.__class__, TreeNode):
##          nodes = nodes.setClass([nodes])

##      if len(nodes)==0:
##          levelSet = str(levelType)+'Set([])'
##          return eval(levelSet), []

##      # catch the case when nodes is already a MoleculeSet
##      if nodes.elementType == levelType:
##          levelSets = nodes
##      else: # get the set of molecules
##          levelSets = nodes.findType(levelType).uniq()

##      # build the set of nodes for each molecule
##      nodeSets = []

##      # find out the type of the nodes we want to return
##      searchType=0
##      if nodeType is None:
##          Klass = nodes.elementType # class of objects in that set
##      else:
##          assert issubclass(nodeType, TreeNode)
##          Klass = nodeType
##          if Klass != nodes.elementType:
##              searchType=1

##      for level in levelSets:
##          # get set of nodes for this molecule
##          mol_nodes = nodes.get(lambda x, mol=mol: x.top==mol)
##          #level_nodes = nodes.get(lambda x,
##          #nodes.get(lambda x, levelType = levelType, level = level: x.findType(levelType)        # get the required types of nodes
##          if searchType:
##              if Klass == Atom and hasattr(mol_nodes, 'allAtoms'):
##                  mol_nodes = mol_nodes.allAtoms
##              else:
##                  mol_nodes = mol_nodes.findType( Klass )

##          nodeSets.append( mol_nodes )

##      return molecules, nodeSets
