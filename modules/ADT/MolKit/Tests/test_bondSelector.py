#
#
#
#
# $Id: test_bondSelector.py,v 1.12 2009/07/10 22:23:31 rhuey Exp $
#

import unittest
from string import split
from MolKit.bondSelector import BondSelector
from MolKit.bondSelector import BondBranchBondSelector
from MolKit.bondSelector import BondElementBondSelector
from MolKit.bondSelector import LeafBondSelector
from MolKit.bondSelector import HydrogenRotatorBondSelector
from MolKit.bondSelector import AmideBondSelector
from MolKit.bondSelector import PeptideBackBoneBondSelector
from MolKit.bondSelector import CycleBondSelector
from MolKit.bondSelector import BondOrderBondSelector
from MolKit.bondSelector import RotatableBondSelector
from MolKit.bondSelector import GuanidiniumBondSelector
from MolKit.bondSelector import AromaticCycleBondSelector
from MolKit.bondSelector import AromaticCycleBondSelector2
from MolKit.bondSelector import PeptideAromaticCycleBondSelector



class BondSelectorBaseTests(unittest.TestCase):
    def setUp(self):
        from MolKit import Read
        self.mol = Read('Data/1crn_2.pdb')[0]
        self.mol.buildBondsByDistance()
    

    def tearDown(self):
        """
        clean-up
        """
        del(self.mol)


    def test_constructor(self):
        """
        instantiate an BondSelector
        """
        bndSel = BondSelector()
        self.assertEquals(bndSel.__class__, BondSelector)


    def test_constructorOptions(self):
        """
         test possible constructor options
            options = [lambda x: cond1, lambda x: cond2...]
        """
    
        bndSel = BondSelector([lambda x: x])
        self.assertEquals(bndSel.__class__, BondSelector)


    def test_select_with_no_bonds(self):
        """
         test select with no input bonds
        """
        bndSel = BondSelector([lambda x:x])
        bndSel.select()
        #clean-up


    def test_select_with_no_criteria(self):
        """
         test select with no screening criteria
        """
        bndSel = BondSelector()
        ats = self.mol.allAtoms[:10]
        bnds = ats.bonds[0]
        bndSel.select(bnds)


    def test_select_no_bonds(self):
        """
         test select with criteria which selects no bonds
        """
        bndSel = BondSelector([lambda x:x.atom1.number==2000])
        ats = self.mol.allAtoms[:10]
        bnds = ats.bonds[0]
        resultBnds = bndSel.select(bnds)
        self.assertEqual(len(resultBnds) , 0)


    def test_select_one_bond(self):
        """
         test select with criteria which selects one bond
        """
        bndSel = BondSelector([lambda x:x.atom1.number==1])
        ats = self.mol.allAtoms[:10]
        bnds = ats.bonds[0]
        resultBnds = bndSel.select(bnds)
        #self.assertEqual(len(resultBnds) , 1)
        self.assertEqual(len(resultBnds) , 3)


    def test_select_three_bonds(self):
        """
         test select with criteria which selects some bonds
        """
        bndSel = BondSelector([lambda x:x.atom1.number==1 or x.atom1.number==2])
        ats = self.mol.allAtoms[:10]
        bnds = ats.bonds[0]
        resultBnds = bndSel.select(bnds)
        #print "len(resultBnds)=", len(resultBnds)
        #self.assertEqual(len(resultBnds) , 3)
        self.assertEqual(len(resultBnds) , 5)


    def test_select_all_bonds(self):
        """
         test select with criteria which selects all bonds
        """
        bndSel = BondSelector([lambda x:x])
        ats = self.mol.allAtoms[:10]
        bnds = ats.bonds[0]
        resultBnds = bndSel.select(bnds)
        self.assertEqual(len(resultBnds) , len(bnds))


    def test_select_two_criteria(self):
        """
         test select with two criteria, result is sum 
        """
        bndSel = BondSelector([lambda x:x.atom1.number==1, lambda x:\
                    x.atom1.number==2])
        ats = self.mol.allAtoms[:10]
        bnds = ats.bonds[0]
        resultBnds = bndSel.select(bnds)
        #self.assertEqual(len(resultBnds) , 3)
        self.assertEqual(len(resultBnds) , 5)


    def test_select_three_criteria(self):
        """
         test select with three criteria, result is sum 
        """
        bndSel = BondSelector([lambda x:x.atom1.number==1, lambda x:\
                    x.atom1.number==2,lambda x:x.atom1.name=='CA'])
        ats = self.mol.allAtoms[:10]
        bnds = ats.bonds[0]
        resultBnds = bndSel.select(bnds)
        #print 'len(resultBnds)==', len(resultBnds)
        #self.assertEqual(len(resultBnds) , 4)
        self.assertEqual(len(resultBnds) , 5)


    def test_BondBranchBondSelector_constructor(self):
        """
         test BondBranch constructor
        """
        bndLenSel = BondBranchBondSelector(2)
        self.assertEqual( bndLenSel.__class__,BondBranchBondSelector)


    def test_BondBranchBondSelector_select_no_bnds(self):
        """
         test BondBranch select: no bonds
        """
        bndSel = BondBranchBondSelector(2)
        ats = self.mol.allAtoms[0:5]
        bnds = ats.bonds[0]
        resultBnds = bndSel.select(bnds)
        #print 'no_bnds:len(resultBnds)=', len(resultBnds)
        self.assertEqual(len(resultBnds) , 0)


    def test_BondBranchBondSelector_select_some_bnds(self):
        """
         test BondBranch select: some bonds
        """
        bndSel = BondBranchBondSelector(2)
        ats = self.mol.allAtoms[:10]
        bnds = ats.bonds[0]
        resultBnds = bndSel.select(bnds)
        #print 'some_bnds:len(resultBnds)=', len(resultBnds)
        #self.assertEqual(len(resultBnds) , 1)
        self.assertEqual(len(resultBnds) , 2)


    def test_BondBranchBondSelector_select_all_bnds(self):
        """
         test BondBranch select: all bonds
        """
        bndSel = BondBranchBondSelector(2)
        ats = self.mol.allAtoms
        bnds = ats.bonds[0]
        resultBnds = bndSel.select(bnds)
        #print 'all_bnds:len(resultBnds)=', len(resultBnds)
        self.assertEqual(len(resultBnds) , 120)


    def test_BondElementBondSelector_constructor(self):
        """
         test BondElement constructor
        """
        bndLenSel = BondElementBondSelector('C')
        self.assertEqual( bndLenSel.__class__,BondElementBondSelector)


    def test_BondElementBondSelector_select_no_bnds(self):
        """
         test BondElement select: no bonds
        """
        bndSel = BondElementBondSelector('P')
        ats = self.mol.allAtoms
        bnds = ats.bonds[0]
        resultBnds = bndSel.select(bnds)
        #print 'no_bnds:len(resultBnds)=', len(resultBnds)
        self.assertEqual(len(resultBnds) , 0)


    def test_BondElementBondSelector_select_some_bnds(self):
        """
         test BondElement select: some bonds
        """
        bndSel = BondElementBondSelector('C')
        ats = self.mol.allAtoms[:10]
        bnds = ats.bonds[0]
        resultBnds = bndSel.select(bnds)
        #print 'some_bnds:len(resultBnds)=', len(resultBnds)
        #self.assertEqual(len(resultBnds) , 9)
        self.assertEqual(len(resultBnds) , 6)


    def test_BondElementBondSelector_select_all_bnds(self):
        """
         test BondElement select: all bonds
        """
        bndSel = BondElementBondSelector('C')
        ats = self.mol.allAtoms
        bnds = ats.bonds[0]
        resultBnds = bndSel.select(bnds)
        #print 'all_bnds:len(resultBnds)=', len(resultBnds)
        self.assertEqual(len(resultBnds) , 334)


    def test_LeafBondSelector_constructor(self):
        """
         test LeafBond constructor
        """
        leafBndSel = LeafBondSelector()
        self.assertEqual( leafBndSel.__class__, LeafBondSelector)


    def test_LeafBondSelector_select_no_leaves(self):
        """
         test LeafBond select: no bonds
        """
        leafBndSel = LeafBondSelector()
        ats = self.mol.allAtoms[7:10]
        bnds = ats.bonds[0]
        resultBnds = leafBndSel.select(bnds)
        #print 'no_leaves:len(resultBnds)=', len(resultBnds)
        self.assertEqual(len(resultBnds) , 0)


    def test_LeafBondSelector_select_some_leaves(self):
        """
         test LeafBond select: some bonds
        """
        leafBndSel = LeafBondSelector()
        ats = self.mol.allAtoms[:10]
        bnds = ats.bonds[0]
        resultBnds = leafBndSel.select(bnds)
        #print 'some_leaves:len(resultBnds)=', len(resultBnds)
        self.assertEqual(len(resultBnds) , 2)
        #self.assertEqual(len(resultBnds) , 5)


    def test_LeafBondSelector_select_all_leaves(self):
        """
         test LeafBond select: all bonds
        """
        leafBndSel = LeafBondSelector()
        ats = self.mol.allAtoms
        bnds = ats.bonds[0]
        resultBnds = leafBndSel.select(bnds)
        #print 'all_leaves:len(resultBnds)=', len(resultBnds)
        self.assertEqual(len(resultBnds) , 150)


    def test_HydrogenRotatorBondSelector_constructor(self):
        """
         test HydrogenRotator constructor
        """
        hrBndSel = HydrogenRotatorBondSelector()
        self.assertEqual( hrBndSel.__class__, HydrogenRotatorBondSelector)


    def test_HydrogenRotatorBondSelector_select_none(self):
        """
         test HydrogenRotator select: no bonds
        """
        HRSel = HydrogenRotatorBondSelector()
        ats = self.mol.allAtoms[2:6]
        bnds = ats.bonds[0]
        resultBnds = HRSel.select(bnds)
        #print 'no_hydrogenrotators:len(resultBnds)=', len(resultBnds)
        self.assertEqual(len(resultBnds) , 0)


    def test_HydrogenRotatorBondSelector_select_some(self):
        """
         test HydrogenRotator select: some bonds
        """
        HRSel = HydrogenRotatorBondSelector()
        ats = self.mol.chains.residues[:8].atoms
        bnds = ats.bonds[0]
        resultBnds = HRSel.select(bnds)
        #print 'some_hydrogenrotators:len(resultBnds)=', len(resultBnds)
        self.assertEqual(len(resultBnds) , 4)


    def test_HydrogenRotatorBondSelector_select_all(self):
        """
         test HydrogenRotator select: all bonds
        """
        HRSel = HydrogenRotatorBondSelector()
        ats = self.mol.allAtoms
        bnds = ats.bonds[0]
        resultBnds = HRSel.select(bnds)
        #print 'all_hydrogenrotators:len(resultBnds)=', len(resultBnds)
        self.assertEqual(len(resultBnds) , 18)


    def test_AmideBondSelector_constructor(self):
        """
         test AmideBond constructor
        """
        amBndSel = AmideBondSelector()
        self.assertEqual(amBndSel.__class__, AmideBondSelector)


    def test_AmideBondSelector_select_none(self):
        """
         test AmideBond select: no bonds
        """
        from MolKit.bondSelector import AmideBondSelector
        amBndSel = AmideBondSelector()
        ats = self.mol.allAtoms[0:2]
        bnds = ats.bonds[0]
        resultBnds = amBndSel.select(bnds)
        #print 'no_amidebonds:len(resultBnds)=', len(resultBnds)
        self.assertEqual(len(resultBnds) , 0)


    def test_AmideBondSelector_select_some(self):
        """
         test AmideBond select: some bonds
        """
        amBndSel = AmideBondSelector()
        ats = self.mol.allAtoms[:50]
        bnds = ats.bonds[0]
        resultBnds = amBndSel.select(bnds)
        #print 'some_amidebonds:len(resultBnds)=', len(resultBnds)
        #self.assertEqual(len(resultBnds) , 7)
        self.assertEqual(len(resultBnds) , 6)


    def test_AmideBondSelector_select_all(self):
        """
         test AmideBond select: all bonds
        """
        amBndSel = AmideBondSelector()
        ats = self.mol.allAtoms
        bnds = ats.bonds[0]
        resultBnds = amBndSel.select(bnds)
        #print 'all_amidebonds:len(resultBnds)=', len(resultBnds)
        #assert len(resultBnds) == 2
        self.assertEqual(len(resultBnds) , 48)


    def test_PeptideBackBoneBondSelector_constructor(self):
        """
         test PetideBackBone constructor
        """
        ppbbSel = PeptideBackBoneBondSelector()
        self.assertEqual(ppbbSel.__class__, PeptideBackBoneBondSelector)


    def test_PeptideBackBoneBondSelector_select_none(self):
        """
         test PetideBackBone select: no bonds
        """
        ppbbSel = PeptideBackBoneBondSelector()
        ats = self.mol.allAtoms[4:8]
        bnds = ats.bonds[0]
        resultBnds = ppbbSel.select(bnds)
        #print 'no_ppbbbonds:len(resultBnds)=', len(resultBnds)
        self.assertEqual(len(resultBnds) , 0)


    def test_PeptideBackBoneBondSelector_select_some(self):
        """
         test PetideBackBone select: some bonds
        """
        ppbbSel = PeptideBackBoneBondSelector()
        ats = self.mol.allAtoms[:50]
        bnds = ats.bonds[0]
        resultBnds = ppbbSel.select(bnds)
        #print 'some_ppbbbonds:len(resultBnds)=', len(resultBnds)
        #self.assertEqual(len(resultBnds) , 22)
        self.assertEqual(len(resultBnds) , 17)


    def test_PeptideBackBoneBondSelector_select_all(self):
        """
         test PetideBackBone select: all bonds
        """
        ppbbSel = PeptideBackBoneBondSelector()
        ats = self.mol.allAtoms
        bnds = ats.bonds[0]
        resultBnds = ppbbSel.select(bnds)
        #print 'all_ppbbbonds:len(resultBnds)=', len(resultBnds)
        self.assertEqual(len(resultBnds) , 132)


    def test_CycleBondSelector_constructor(self):
        """
         test CycleBondSelector constructor
        """
        cbSel = CycleBondSelector()
        self.assertEqual(cbSel.__class__, CycleBondSelector)


    def test_CycleBondSelector_select_none(self):
        """
         test CycleBondSelector select: no bonds
        """
        cbSel = CycleBondSelector()
        ats = self.mol.allAtoms[4:8]
        bnds = ats.bonds[0]
        resultBnds = cbSel.select(bnds)
        #print 'no_cb_bonds:len(resultBnds)=', len(resultBnds)
        self.assertEqual(len(resultBnds) , 0)


    def test_CycleBondSelector_select_some(self):
        """
         test CycleBondSelector select: some bonds
        """
        cbSel = CycleBondSelector()
        ats = self.mol.allAtoms[:50]
        bnds = ats.bonds[0]


    def test_BondElementBondSelector_select_all_bnds(self):
        """
         test BondElement select: all bonds
        """
        bndSel = BondElementBondSelector('C')
        ats = self.mol.allAtoms
        bnds = ats.bonds[0]
        resultBnds = bndSel.select(bnds)
        #print 'all_bnds:len(resultBnds)=', len(resultBnds)
        self.assertEqual(len(resultBnds) , 334)


    def test_LeafBondSelector_constructor(self):
        """
         test LeafBond constructor
        """
        leafBndSel = LeafBondSelector()
        self.assertEqual( leafBndSel.__class__, LeafBondSelector)


    def test_LeafBondSelector_select_no_leaves(self):
        """
         test LeafBond select: no bonds
        """
        leafBndSel = LeafBondSelector()
        ats = self.mol.allAtoms[7:10]
        bnds = ats.bonds[0]
        resultBnds = leafBndSel.select(bnds)
        #print 'no_leaves:len(resultBnds)=', len(resultBnds)
        self.assertEqual(len(resultBnds) , 0)


    def test_LeafBondSelector_select_some_leaves(self):
        """
         test LeafBond select: some bonds
        """
        leafBndSel = LeafBondSelector()
        ats = self.mol.allAtoms[:10]
        bnds = ats.bonds[0]
        resultBnds = leafBndSel.select(bnds)
        #print 'some_leaves:len(resultBnds)=', len(resultBnds)
        #self.assertEqual(len(resultBnds) , 2)
        self.assertEqual(len(resultBnds) , 5)


    def test_LeafBondSelector_select_all_leaves(self):
        """
         test LeafBond select: all bonds
        """
        leafBndSel = LeafBondSelector()
        ats = self.mol.allAtoms
        bnds = ats.bonds[0]
        resultBnds = leafBndSel.select(bnds)
        #print 'all_leaves:len(resultBnds)=', len(resultBnds)
        self.assertEqual(len(resultBnds) , 150)


    def test_HydrogenRotatorBondSelector_constructor(self):
        """
         test HydrogenRotator constructor
        """
        hrBndSel = HydrogenRotatorBondSelector()
        self.assertEqual( hrBndSel.__class__, HydrogenRotatorBondSelector)


    def test_HydrogenRotatorBondSelector_select_none(self):
        """
         test HydrogenRotator select: no bonds
        """
        HRSel = HydrogenRotatorBondSelector()
        ats = self.mol.allAtoms[2:6]
        bnds = ats.bonds[0]
        resultBnds = HRSel.select(bnds)
        #print 'no_hydrogenrotators:len(resultBnds)=', len(resultBnds)
        self.assertEqual(len(resultBnds) , 0)


    def test_HydrogenRotatorBondSelector_select_some(self):
        """
         test HydrogenRotator select: some bonds
        """
        HRSel = HydrogenRotatorBondSelector()
        ats = self.mol.chains.residues[:8].atoms
        bnds = ats.bonds[0]
        resultBnds = HRSel.select(bnds)
        #print 'some_hydrogenrotators:len(resultBnds)=', len(resultBnds)
        self.assertEqual(len(resultBnds) , 4)


    def test_HydrogenRotatorBondSelector_select_all(self):
        """
         test HydrogenRotator select: all bonds
        """
        HRSel = HydrogenRotatorBondSelector()
        ats = self.mol.allAtoms
        bnds = ats.bonds[0]
        resultBnds = HRSel.select(bnds)
        #print 'all_hydrogenrotators:len(resultBnds)=', len(resultBnds)
        self.assertEqual(len(resultBnds) , 18)


    def test_AmideBondSelector_constructor(self):
        """
         test AmideBond constructor
        """
        amBndSel = AmideBondSelector()
        self.assertEqual(amBndSel.__class__, AmideBondSelector)


    def test_AmideBondSelector_select_none(self):
        """
         test AmideBond select: no bonds
        """
        from MolKit.bondSelector import AmideBondSelector
        amBndSel = AmideBondSelector()
        ats = self.mol.allAtoms[0:2]
        bnds = ats.bonds[0]
        resultBnds = amBndSel.select(bnds)
        #print 'no_amidebonds:len(resultBnds)=', len(resultBnds)
        self.assertEqual(len(resultBnds) , 0)


    def test_AmideBondSelector_select_some(self):
        """
         test AmideBond select: some bonds
        """
        amBndSel = AmideBondSelector()
        ats = self.mol.allAtoms[:50]
        bnds = ats.bonds[0]
        resultBnds = amBndSel.select(bnds)
        #print 'some_amidebonds:len(resultBnds)=', len(resultBnds)
        #self.assertEqual(len(resultBnds) , 7)
        self.assertEqual(len(resultBnds) , 6)


    def test_AmideBondSelector_select_all(self):
        """
         test AmideBond select: all bonds
        """
        amBndSel = AmideBondSelector()
        ats = self.mol.allAtoms
        bnds = ats.bonds[0]
        resultBnds = amBndSel.select(bnds)
        #print 'all_amidebonds:len(resultBnds)=', len(resultBnds)
        #assert len(resultBnds) == 2
        self.assertEqual(len(resultBnds) , 48)


    def test_PeptideBackBoneBondSelector_constructor(self):
        """
         test PetideBackBone constructor
        """
        ppbbSel = PeptideBackBoneBondSelector()
        self.assertEqual(ppbbSel.__class__, PeptideBackBoneBondSelector)


    def test_PeptideBackBoneBondSelector_select_none(self):
        """
         test PetideBackBone select: no bonds
        """
        ppbbSel = PeptideBackBoneBondSelector()
        ats = self.mol.allAtoms[4:8]
        bnds = ats.bonds[0]
        resultBnds = ppbbSel.select(bnds)
        #print 'no_ppbbbonds:len(resultBnds)=', len(resultBnds)
        self.assertEqual(len(resultBnds) , 0)


    def test_PeptideBackBoneBondSelector_select_some(self):
        """
         test PetideBackBone select: some bonds
        """
        ppbbSel = PeptideBackBoneBondSelector()
        ats = self.mol.allAtoms[:50]
        bnds = ats.bonds[0]
        resultBnds = ppbbSel.select(bnds)
        #print 'some_ppbbbonds:len(resultBnds)=', len(resultBnds)
        #self.assertEqual(len(resultBnds) , 22)
        self.assertEqual(len(resultBnds) , 17)


    def test_PeptideBackBoneBondSelector_select_all(self):
        """
         test PetideBackBone select: all bonds
        """
        ppbbSel = PeptideBackBoneBondSelector()
        ats = self.mol.allAtoms
        bnds = ats.bonds[0]
        resultBnds = ppbbSel.select(bnds)
        #print 'all_ppbbbonds:len(resultBnds)=', len(resultBnds)
        self.assertEqual(len(resultBnds) , 132)


    def test_CycleBondSelector_constructor(self):
        """
         test CycleBondSelector constructor
        """
        cbSel = CycleBondSelector()
        self.assertEqual(cbSel.__class__, CycleBondSelector)


    def test_CycleBondSelector_select_none(self):
        """
         test CycleBondSelector select: no bonds
        """
        cbSel = CycleBondSelector()
        ats = self.mol.allAtoms[4:8]
        bnds = ats.bonds[0]
        resultBnds = cbSel.select(bnds)
        #print 'no_cb_bonds:len(resultBnds)=', len(resultBnds)
        self.assertEqual(len(resultBnds) , 0)


    def test_CycleBondSelector_select_some(self):
        """
         test CycleBondSelector select: some bonds
        """
        cbSel = CycleBondSelector()
        ats = self.mol.allAtoms[:50]
        bnds = ats.bonds[0]
        resultBnds = cbSel.select(bnds)
        #print 'some_cb_bonds:len(resultBnds)=', len(resultBnds)
        #self.assertEqual(len(resultBnds) , 5)
        #because of di-sulfide bridges:
        #when increase RingFinder parameter "maxSize"
        self.assertEqual(len(resultBnds) , 42)


    def test_CycleBondSelector_select_all(self):
        """
         test CycleBondSelector select: all bonds
        """
        cbSel = CycleBondSelector()
        ats = self.mol.allAtoms
        bnds = ats.bonds[0]
        resultBnds = cbSel.select(bnds)
        #print 'all_cb_bonds:len(resultBnds)=', len(resultBnds)
        #self.assertEqual(len(resultBnds) , 43)
        #because of di-sulfide bridges:
        #when increase RingFinder parameter "maxSize"
        self.assertEqual(len(resultBnds) , 165)

    
    def test_BondOrderBondSelector(self):
        """
         test BondOrder constructor
        """
        bobSel = BondOrderBondSelector()
        self.assertEqual(bobSel.__class__, BondOrderBondSelector)


    def test_BondOrderBondSelector_select_none(self):
        """
         test BondOrder select: no bonds
        """
        bobSel = BondOrderBondSelector()
        ats = self.mol.allAtoms[4:8]
        bnds = ats.bonds[0]
        resultBnds = bobSel.select(bnds)
        #print 'no_cb_bonds:len(resultBnds)=', len(resultBnds)
        #self.assertEqual(len(resultBnds) , 2)
        self.assertEqual(len(resultBnds) , 3)


    def test_BondOrderBondSelector_select_some(self):
        """
         test BondOrder select: some bonds
        """
        bobSel = BondOrderBondSelector()
        ats = self.mol.allAtoms[:50]
        bnds = ats.bonds[0]
        resultBnds = bobSel.select(bnds)
        #print 'some_cb_bonds:len(resultBnds)=', len(resultBnds)
        #self.assertEqual(len(resultBnds) , 43)
        self.assertEqual(len(resultBnds) , 44)


    def test_BondOrderBondSelector_select_all(self):
        """
         test BondOrder select: all bonds
        """
        bobSel = BondOrderBondSelector()
        ats = self.mol.allAtoms
        bnds = ats.bonds[0]
        resultBnds = bobSel.select(bnds)
        #print 'all_cb_bonds:len(resultBnds)=', len(resultBnds)
        self.assertEqual(len(resultBnds) , 341)


    def test_RotatableBondSelector_constructor(self):
        """
         test RotatableBond constructor
        """
        rotSel = RotatableBondSelector()
        self.assertEqual(rotSel.__class__, RotatableBondSelector)


    def test_RotatableBondSelector_select_none(self):
        """
         test RotatableBond select: no bonds
        """
        rotSel = RotatableBondSelector()
        ats = self.mol.allAtoms[4:8]
        bnds = ats.bonds[0]
        resultBnds = rotSel.select(bnds)
        #print 'no_cb_bonds:len(resultBnds)=', len(resultBnds)
        self.assertEqual(len(resultBnds) , 1)


    def test_RotatableBondSelector_select_some(self):
        """
         test RotatableBond select: some bonds
        """
        rotSel = RotatableBondSelector()
        ats = self.mol.allAtoms[:50]
        bnds = ats.bonds[0]
        resultBnds = rotSel.select(bnds)
        #print 'some_cb_bonds:len(resultBnds)=', len(resultBnds)
        #self.assertEqual(len(resultBnds) , 34)
        #because of disulfide bridges-> large cycles
        #self.assertEqual(len(resultBnds) , 27)
        self.assertEqual(len(resultBnds) , 20)


    def test_RotatableBondSelector_select_all(self):
        """
         test RotatableBond select: all bonds
        """
        rotSel = RotatableBondSelector()
        ats = self.mol.allAtoms
        bnds = ats.bonds[0]
        resultBnds = rotSel.select(bnds)
        #print 'all_cb_bonds:len(resultBnds)=', len(resultBnds)
        #self.assertEqual(len(resultBnds) , 211)
        #because of disulfide bridges-> large cycles
        self.assertEqual(len(resultBnds), 89)


    def test_GuanidiniumBondSelector_select_all(self):
        """
         test GuanidiniumBond select: all bonds
        """
        guanSel = GuanidiniumBondSelector()
        ats = self.mol.allAtoms
        bnds = ats.bonds[0]
        resultBnds = guanSel.select(bnds)
        #print 'all_guanidinium_bonds:len(resultBnds)=', len(resultBnds)
        self.assertEqual(len(resultBnds), 6)
        for b in resultBnds:
            self.assertEqual(b.atom1.parent.type, 'ARG')
            if b.atom1.element=='C':
                self.assertEqual(b.atom1.name, 'CZ')
            if b.atom1.element=='N':
                self.assertEqual(b.atom1.name in ['NE','NH1','NH2'], True)


    def test_PeptideAromaticCycleBondSelector_select_all(self):
        """
         test PeptideAromaticCycleBond select: all aromatic bonds in 1crn
        """
        PABS = PeptideAromaticCycleBondSelector()
        ats = self.mol.allAtoms
        bnds = ats.bonds[0]
        resultBnds = PABS.select(bnds)
        self.assertEqual(len(resultBnds), 18)


    def test_PeptideAromaticCycleBondSelector_select_some(self):
        """
         test PeptideAromaticCycleBond select: 6 bonds
        """
        PABS = PeptideAromaticCycleBondSelector()
        bnds = self.mol.chains.residues[:21].atoms.bonds[0]  #includes 1 PHE
        resultBnds = PABS.select(bnds)
        self.assertEqual(len(resultBnds), 6)


    def test_PeptideAromaticCycleBondSelector_labels_aromatic_bnds(self):
        """
         test PeptideAromaticCycleBond labels aromatic bnds
        """
        PABS = PeptideAromaticCycleBondSelector()
        bnds = self.mol.chains.residues.atoms.bonds[0]  
        resultBnds = PABS.select(bnds)
        for b in resultBnds:
            self.assertEqual(hasattr(b, 'aromatic'), True)



class AromaticCycleBondTests(unittest.TestCase):

    def setUp(self):
        from MolKit import Read
        self.ind = Read('Data/indinavir.pdb')[0]
        self.ind.buildBondsByDistance()
    

    def tearDown(self):
        """
        clean-up
        """
        del(self.ind)



    def test_AromaticCycleBondSelector_select_all(self):
        """
         test AromaticCycleBond select: all aromatic bonds in ind
         This bondSelector uses PyBabel aromatic object 
         #nb doesnot count any in merged ring (?)
        """
        ABS = AromaticCycleBondSelector()
        ats = self.ind.allAtoms
        bnds = ats.bonds[0]
        resultBnds = ABS.select(bnds)
        #print ABS.getAtoms(resultBnds).name
        self.assertEqual(len(resultBnds), 18)


    def test_AromaticCycleBond2Selector_select_all(self):
        """
         test AromaticCycleBond2 select: all aromatic bonds in ind
         This bndSelector uses autotors calculation of normals between
         adjacent atoms in cycles to judge aromaticity 
         #nb does count aromatic bonds in merged ring (?)
        """
        ABS2 = AromaticCycleBondSelector2()
        ats = self.ind.allAtoms
        bnds = ats.bonds[0]
        resultBnds = ABS2.select(bnds)
        #print ABS2.getAtoms(resultBnds).name
        self.assertEqual(len(resultBnds), 18)
        self.assertEqual(len(ABS2.getAtoms(resultBnds)), 18)


    def test_AromaticCycleBond2Selector_select_all2(self):
        """
         This is repeated because repeating it broke when all tests were in
         one class....(???)
        """
        ABS2 = AromaticCycleBondSelector2()
        ats = self.ind.allAtoms
        bnds = ats.bonds[0]
        resultBnds = ABS2.select(bnds)
        #print ABS2.getAtoms(resultBnds).name
        self.assertEqual(len(resultBnds), 18)
        self.assertEqual(len(ABS2.getAtoms(resultBnds)), 18)


    def test_AromaticCycleBond2Selector_select_all3(self):
        """
         This is repeated because repeating it broke when all tests were in
         one class....(???)
        """
        ABS2 = AromaticCycleBondSelector2()
        ats = self.ind.allAtoms
        bnds = ats.bonds[0]
        resultBnds = ABS2.select(bnds)
        #print ABS2.getAtoms(resultBnds).name
        self.assertEqual(len(resultBnds), 18)
        self.assertEqual(len(ABS2.getAtoms(resultBnds)), 18)



if __name__ == '__main__':
    unittest.main()
