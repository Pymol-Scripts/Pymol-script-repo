#############################################################################
#
# Author: Sophie I. COON, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################
#
# $Header: /opt/cvs/python/packages/share1.5/MolKit/Tests/test_tree.py,v 1.4 2005/02/04 22:47:20 rhuey Exp $
#
# $Id: test_tree.py,v 1.4 2005/02/04 22:47:20 rhuey Exp $
#


import sys
from MolKit.tree import TreeNodeSet
from MolKit.molecule import AtomSet, MoleculeSet
from MolKit.protein import ResidueSet, ChainSet


def test_NodesFromName():
    # test retrieving nodes using a name.
    # Name is a string as produced by full_name() i.e. node names
    # separated by ':' going from root to leaf
    # comma ',' separated list of names are allowed as well as
    # range specified using the '-' character   

    from MolKit import Read
    mols = Read('Data/1crn.pdb')
    mols += Read('Data/1bsr.pdb')

    # look for non-existing string
    result = mols.NodesFromName("thisisridiculous:::")
    assert len(result)==0
    assert isinstance(result, TreeNodeSet)

    # look for all atoms in a particular protein
    result = mols.NodesFromName("1crn:::")
    assert len(result)==327
    assert isinstance(result, AtomSet)

    # check trailing semicolon, atom level
    result = mols.NodesFromName("1crn:::;")
    assert len(result)==327
    assert isinstance(result, AtomSet)

    # check trailing semicolon, residue level
    result = mols.NodesFromName("1crn::;")
    assert len(result)==46
    assert isinstance(result, ResidueSet)

    # check trailing semicolon, chain level
    result = mols.NodesFromName("1crn:;")
    assert len(result)==1
    assert isinstance(result, ChainSet)

    # check trailing semicolon, molecule level
    result = mols.NodesFromName("1crn;")
    assert len(result)==1
    assert isinstance(result, MoleculeSet)



    # look for all chains in a particular protein
    result = mols.NodesFromName("1bsr::")
    assert len(result)==368
    assert isinstance(result, ResidueSet)

    # look for chain B in a particular protein
    result = mols.NodesFromName("1bsr:B")
    assert len(result)==1
    assert isinstance(result, ChainSet)

    # look for all chains in a particular protein
    result = mols.NodesFromName("1bsr:A,B")
    assert len(result)==2
    assert isinstance(result, ChainSet)

    # look for all chains in all proteins
    result = mols.NodesFromName(":")
    assert len(result)==3
    assert isinstance(result, ChainSet)

    # look for anall atoms in an amino acide in a particular protein
    result = mols.NodesFromName("1crn::THR2:")
    assert len(result)==7
    assert isinstance(result, AtomSet)

    # look for an amino acide in a particular protein
    result = mols.NodesFromName("1crn::THR2")
    assert len(result)==1
    assert isinstance(result, ResidueSet)

    # look for a range of amino acide in a particular protein
    result = mols.NodesFromName("1crn::THR2-CYS16")
    assert len(result)==15
    assert isinstance(result, ResidueSet)


def test_sortOnEmptySet():
    """
    Testing the sort method on a empty set
    """
    from MolKit.tree import TreeNodeSet
    set = TreeNodeSet([])
    set.sort()
    del set


def test_merge():
    # merge mol2 into mol1.
    from MolKit import Read
    mol1 = Read("./Data/protease.pdb")[0]
    mol2 = Read("./Data/indinavir.pdb")[0]
    from MolKit.protein import ProteinSet
    pset = ProteinSet([mol1,mol2])
    chains = pset.chains[:]
    mol1.merge(mol2)
    mol1.chains.sort()
    chains.sort()
    assert mol1.chains == chains



