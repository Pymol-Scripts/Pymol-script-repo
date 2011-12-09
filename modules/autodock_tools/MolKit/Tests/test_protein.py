#############################################################################
#
# Author: Sophie I. COON, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################
#
# $Header: /opt/cvs/python/packages/share1.5/MolKit/Tests/test_protein.py,v 1.2 2003/08/29 17:50:16 sophiec Exp $
#
# $Id: test_protein.py,v 1.2 2003/08/29 17:50:16 sophiec Exp $
#


import sys
from mglutil.regression import testplus
"""
This module implements a set of function to test the functionalities
implemented in the protein module.
"""

def test_findType():
    from MolKit import Read
    mol = Read("./Data/1crn.pdb")[0]
    from MolKit.protein import Residue, Chain, Protein, ResidueSet
    from MolKit.molecule import Atom, AtomSet
    # FindType below from top of the tree
    res = mol.findType(Residue)
    res.sort()
    molres = mol.chains.residues
    molres.sort()
    assert res == molres

    # FindType above from bottom of the tree
    atms = mol.allAtoms
    chains = atms.findType(Chain)
    # Automatically does a uniq when going above
    assert len(chains) == 327
    
    # FindType above uniq from bottom of the tree
    chainsuniq = atms.findType(Chain, uniq=1)
    assert len(chainsuniq)==1

    # FindType above from middle of the tree
    prot = res.findType(Protein, uniq=1)
    assert len(prot)==1 and prot == mol.setClass([mol,])

    # FindType below from middle of the tree
    atoms = chainsuniq.findType(Atom)
    assert len(atoms) == 327

    # FindType on empty set
    emptyres = ResidueSet([])
    atoms = emptyres.findType(Atom)
    assert len(atoms)==0 and isinstance(atoms, AtomSet)
    

    # FindType on same type
    atoms = atms.findType(Atom)
    assert atoms == atms

    # Raise exception
    from MolKit.protein import Helix
    try:
        nodes = atms.findType(Helix)
    except RuntimeError:
        print "passed"

    

harness = testplus.TestHarness( __name__,
                                funs = testplus.testcollect( globals()),
                                )

if __name__ == '__main__':
    print harness
    sys.exit( len( harness))

