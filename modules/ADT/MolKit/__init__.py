#############################################################################
#
# Author: Michel F. SANNER, Sophie COON
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

#
#$Header: /opt/cvs/python/packages/share1.5/MolKit/__init__.py,v 1.25.4.1 2011/06/25 02:10:58 sanner Exp $
#
#$Id: __init__.py,v 1.25.4.1 2011/06/25 02:10:58 sanner Exp $
#

import string
import os


def Read(filename, modelsAs='molecules'):
    if not os.path.exists(filename):
         raise AssertionError , "%s does't exist" %filename
    from MolKit.pdbParser import PdbParser, PdbqParser,PdbqsParser,\
            PdbqtParser, PQRParser, F2DParser
    
    from MolKit.mol2Parser import Mol2Parser
    from MolKit.mmcifParser import MMCIFParser
    ext = string.split(filename, '.')
    if ext[-1]=='pdb':
        parser = PdbParser(filename, modelsAs=modelsAs)

    elif ext[-1]=='pdbq':
        parser = PdbqParser(filename, modelsAs=modelsAs)
    
    elif ext[-1]=='pdbqt':
        parser = PdbqtParser(filename, modelsAs=modelsAs)

    elif ext[-1]=='pdbqs':
        parser = PdbqsParser(filename, modelsAs=modelsAs)

    elif ext[-1]=='pqr':
        parser = PQRParser(filename, modelsAs=modelsAs)

    elif ext[-1]=='mol2':
        parser = Mol2Parser(filename) #??should modelsAs be available for mol2 format??

    elif ext[-1]=='cif':
        parser = MMCIFParser(filename, modelsAs=modelsAs)

    elif ext[-1]=='f2d':
        parser = F2DParser(filename)

    else:
        print "File Format unknown can't parse it"
        return []
    molecules = parser.parse()
    return molecules


def WritePDB(filename,node):
    from MolKit.pdbWriter import PdbWriter
    writer = PdbWriter()
    writer.write(filename, node)


def makeMoleculeFromAtoms(molname, atomSet):
    """
    create a new molecule from a list of atoms

    mol <- makeMoleculeFromAtoms(molname, atomSet)
"""
    from MolKit.molecule import Atom, AtomSet
    from MolKit.protein import Protein, Chain, Residue


    # create the top object
    mol = Protein(name=molname)

    # find out all residues
    residues = atomSet.parent.uniq()

    # find out all chains
    chains = residues.parent.uniq()

    # create all chains
    chainsd = {}
    for c in chains:
        newchain = Chain(c.id, mol, top=mol)
        chainsd[c] = newchain

    # create all residues
    resd = {}
    for res in residues:
        newres = Residue(res.name[:3], res.name[3:], res.icode,
                         chainsd[res.parent], top=mol)
        resd[res] = newres
        newres.hasCA = 0
        newres.hasO = 0

    # create all the atoms
    newats = []
    for num, at in enumerate(atomSet):
        name = at.name
        res = resd[at.parent]
        if name == 'CA':
            res.hasCA = 1
        if name == 'O' or name == 'OXT' or (len(name)>3 and name[:3]=='OCT'):
            res.hasO = 2
        
        newat = Atom(name, res, at.element, top=mol)
        newats.append(newat)
        # set constructotr attributes
        newat._coords = []
        for coords in at._coords:
            newat._coords.append(coords[:])
        newat.conformation = at.conformation
        newat.chemElem = at.chemElem
        newat.atomicNumber = at.atomicNumber
        newat.bondOrderRadius = at.bondOrderRadius
        newat.covalentRadius = at.covalentRadius
        newat.vdwRadius = at.vdwRadius
        newat.maxBonds = at.maxBonds
        newat.organic = at.organic
        newat.colors = at.colors.copy()
        newat.opacities = at.opacities.copy()
        newat._charges = at._charges.copy()
        newat.chargeSet = at.chargeSet

        # set attributes from PDB parser
        newat.segID = at.segID
        newat.hetatm = at.hetatm
        newat.normalname = at.normalname
        newat.number = num #at.number
        newat.occupancy = at.occupancy
        newat.temperatureFactor = at.temperatureFactor
        newat.altname = at.altname

        # attribute created by PQR parser
        if hasattr(at, 'pqrRadius'):
            newat.pqrRadius = at.pqrRadius

        # attribute created by F2D parser
        if hasattr(at, 'hbstatus'):
            newat.hbstatus = at.hbstatus

        # attribute created by PDBQ parser
        if hasattr(at, 'autodock_element'):
            newat.autodock_element = at.autodock_element

        # attribute created by PDBQT parser
        #if hasattr(at, ''):
        #    newat. = at.

        # attribute created by PDBQS parser
        if hasattr(at, 'AtVol'):
            newat.AtVol = at.AtVol
            newat.AtSolPar = at.AtSolPar

    mol.allAtoms = AtomSet(newats)
    return mol

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

CRITICAL_DEPENDENCIES = ['mglutil', 'numpy']
NONCRITICAL_DEPENDENCIES =['sff', 'PyBabel', 'stride', 'bhtree', 'NetworkEditor', 'DejaVu', 'mslib', 'Vision','Pmv', 'cMolKit', 'symserv', '_xmlplus']
