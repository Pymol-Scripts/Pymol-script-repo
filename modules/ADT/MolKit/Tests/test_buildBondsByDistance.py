#############################################################################
#
# Author: Sophie I. COON, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################
#
# $Header: /opt/cvs/python/packages/share1.5/MolKit/Tests/test_buildBondsByDistance.py,v 1.2 2003/08/29 17:41:05 sophiec Exp $
#
# $Id: test_buildBondsByDistance.py,v 1.2 2003/08/29 17:41:05 sophiec Exp $
#

"""
This module implements method to test buildbondsbydistance creates
bonds between two atoms, in a residues, in a chain or in molecule.
1- Need to test buildBondsByDistance for a molecule using bhtree
   - molecule with 1 chain, several chains, alternate locations, ....

2- Need to test buildBondsByDistance for a molecule using the old
   buildBondsByDistance.
   - molecule with 1 chain, several chains, with alternate locations.....
Need to test at the Molecule, Chain, Residue and Atoms level.
"""

import sys
from exceptions import AssertionError
from mglutil.regression import testplus
global bhtreeFlag
try:
    import bhtree
    bhtreeFlag = 1
except:
    bhtreeFlag = 0

def test_buildBondsByDistance_1():
    """
    Test the buildBondsByDistance for 1crn.pdb with one chain
    
    """
    from MolKit import Read
    mol = Read('Data/1crn.pdb')[0]
    mol.buildBondsByDistance()
    bonds, nobonds = mol.allAtoms.bonds[0], mol.allAtoms.bonds[1]
    assert len(bonds) == 337 and len(nobonds) == 0


## def test_buildBondsByDistance_2():
##     """
##     Test the molecule.buildBondsByDistance method for a molecule having a lot
##     of alternate locations. Need to make sure that the alternate @A are not
##     bound to the alternates @B.
##     We need to make sure that the right bonds are created as well.
##     """
##     from MolKit import Read
##     mol = Read('Data/1hxblig.pdb')[0]
##     mol.buildBondsByDistance()
##     #1- Make sure that uses
##     # get all the @A alternate
##     altA = mol.allAtoms.get(lambda x: x.name[-2:] == '@A')
##     altB = mol.allAtoms.get(lambda x: x.name[-2:] == '@B')
##     for atmA in altA:
##         for atmB in altB:
##             if atmA.isBonded(atmB):
##                 raise AssertionError(atmA.name+" is bonded to "+atmB.name)
##     del mol

## def test_buildBondsByDistance_3():
##     """
##     Test the molecule.buildBondsByDistance method written to fix a bug with
##     Hydrogen having two many bonds
##     2 Bonds were created using the CONNECT records and 3 Bonds by the
##     buildBondsByDistance.
##     Using bhtree to build those bonds.
##     """
##     from MolKit import Read
##     mol = Read('Data/diversity0744.pdb')[0]
##     mol.buildBondsByDistance()
##     h62 = mol.allAtoms.get(lambda x: x.number == 62)[0]
##     h61 = mol.allAtoms.get(lambda x: x.number == 61)[0]
##     h49 = mol.allAtoms.get(lambda x: x.number == 49)[0]
##     # 1- Assert that those two atoms are only bonded to two other atoms
##     error = ""
##     if len(h62.bonds) != 1:
##         error = error + "H62: "
##         for b in h62.bonds:
##             error = error + b.atom1.name + str(b.atom1.number) \
##                     + "--" + b.atom2.name + str(b.atom2.number) + "; "
##     if len(h61.bonds) != 1:
##         error = error + "H61: " 
##         for b in h61.bonds:
##             error = error + b.atom1.name + str(b.atom1.number) + \
##                     "--" + b.atom2.name + str(b.atom2.number) + "; "
    
##     if len(h49.bonds) != 1:
##         error = error + "H49: " 
##         for b in h49.bonds:
##             error = error + b.atom1.name + str(b.atom1.number) + \
##                     "--" + b.atom2.name + str(b.atom2.number) + "; "
    
##     if error != "":
##         raise AssertionError(error)


## def test_buildBondsByDistance_4():
##     """
##     Test the molecule.buildBondsByDistance method written to fix a bug with
##     Hydrogen having two many bonds same molecule than before but without
##     using the CONNECT records.
##     Using bhtree to build those bonds.
##     """
##     from MolKit import Read
##     mol = Read('Data/diversity0744NoConnect.pdb')[0]
##     mol.buildBondsByDistance()
##     c27 = mol.allAtoms.get(lambda x: x.number == 27)[0]
##     c28 = mol.allAtoms.get(lambda x: x.number == 28)[0]
##     # 1- Assert that those two atoms are only bonded to two other atoms
##     error = ""
##     print len(c27.bonds)
##     if len(c27.bonds) == 4:
##         error = error + "C27: "
##         for b in c27.bonds:
##             error = error + b.atom1.name + str(b.atom1.number) \
##                     + "--" + b.atom2.name + str(b.atom2.number) + "; "
##     print len(c28.bonds)
##     if len(c28.bonds) == 4:
##         error = error + "C28: " 
##         for b in c28.bonds:
##             error = error + b.atom1.name + str(b.atom1.number) + \
##                     "--" + b.atom2.name + str(b.atom2.number) + "; "
    
    
##     if error != "":
##         raise AssertionError(error)

harness = testplus.TestHarness( __name__,
                                funs = testplus.testcollect( globals()),
                                )

if __name__ == '__main__':
    print harness
    sys.exit( len( harness))
