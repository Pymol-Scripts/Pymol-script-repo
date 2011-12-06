## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#############################################################################
#
# Authors: Ruth Huey, William M. Lindstrom
#
# Copyright: R. Huey, W. M. Lindstrom TSRI 2003
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/MolKit/bondSelector.py,v 1.25.6.1 2011/01/04 17:45:17 rhuey Exp $
#
# $Id: bondSelector.py,v 1.25.6.1 2011/01/04 17:45:17 rhuey Exp $
#
#
#

"""
This module implements classes which select bonds based on a list of criteria.
They return the selected bonds as a BondSet
"""

from MolKit.molecule import Bond, BondSet, AtomSet
from PyBabel.cycle import RingFinder
from PyBabel.bo import BondOrder
from PyBabel.aromatic import Aromatic
from PyBabel.atomTypes import AtomHybridization
import math, numpy.oldnumeric as Numeric


class BondSelector:
    """ Base class that selects bonds based on a list of criteria.
    """


    def __init__(self, criteria=[None], uniq=1):
        self.criteria = criteria
        self.uniq = uniq


    def select(self, bonds=None):
        """ 
        use self.criteria to select some bnds from input bonds
        """
        #make sure that select is called with some bonds
        if not bonds:
            return BondSet([])
        if not len(bonds):
            return BondSet([])
        assert isinstance(bonds, BondSet)
        returnBonds = BondSet()
        for f in self.criteria:
            newBonds = filter(f, bonds)
            returnBonds.extend(BondSet(filter(f, bonds)))
        return self.makeUniq(returnBonds)


    def getAtoms(self,bnds):
        ats0 = AtomSet()
        for b in bnds:
            ats0.append(b.atom1)
            ats0.append(b.atom2)
        d = {}
        for a in ats0:
            d[a] = 0
        return AtomSet(d.keys())


    def makeUniq(self, bnds):
        if not self.uniq:
            return bnds
        #if called with an empty list
        # an empty list is returned
        d = {}
        for b in bnds:
            d[b] = 0
        uniqBonds = d.keys()
        if len(uniqBonds)==len(bnds):
            return BondSet(bnds)
        else:
            #if there's a duplicate
            #return list with duplicate removed
            return BondSet(uniqBonds)



class BondBranchBondSelector(BondSelector):
    """ class that selects bonds where either atom has specified 
    number of bonds
    eg: if number is 1, the bond is a leaf bond eg -C-H
    """


    def __init__(self, number):
        self.number = number
        criteria = [lambda x, num=number: len(x.atom1.bonds)==num,
                    lambda x, num=number: len(x.atom2.bonds)==num]
        BondSelector.__init__(self, criteria=criteria)



class BondElementBondSelector(BondSelector):
    """ class that selects bonds where either atom is of
        specified element
    eg: if element is S, the bond is a S-X or X-S where X is anything
        including S
    """

    def __init__(self, elem):
        self.elem = elem
        criteria = [lambda x, elem=self.elem: x.atom1.element==elem,
                    lambda x, elem=self.elem: x.atom2.element==elem]
        BondSelector.__init__(self, criteria=criteria)


class LeafBondSelector(BondBranchBondSelector):
    """ class that selects bonds where 1 of the two atoms is a leaf, ie has
    no other bonds
    """

    def __init__(self):
        BondBranchBondSelector.__init__(self, 1)



class HydrogenRotatorBondSelector(BondSelector):
    """ class that selects bonds which only rotate hydrogens: 
    ie: all the other bonds of either atom in the bond are 
        to hydrogen atoms.
    only usefulness of inheriting from BondSelector is reuse of makeUniq
    """


    def __init__(self):
        self.BES = BondElementBondSelector('H')
        BondSelector.__init__(self)
        #can this be fixed??
        #criteria = [lambda x, hsel=self.BES: \
        #    len(x.atom1.bonds)-1==len(hsel(x.bonds)), \
        #    lambda x, hsel=self.BES:\
        #    len(x.atom2.bonds)-1==len(hsel(x.bonds))]


    def select(self, bnds):
        hbnds = self.BES.select(bnds)
        #all atoms in any bond made with a hydrogen, incl hydrogen ats
        all_ats = self.BES.getAtoms(hbnds)
        #atoms in bonds with a hydrogen which are not hydrogens
        not_h = all_ats.get(lambda x:x.element!='H')
        #non_h atoms which only have 1 bond to a non-hydrogen atom
        return not_h.get(lambda x: len(x.findHydrogens())==len(x.bonds)-1)



class AmideBondSelector(BondSelector):
    """ class that selects amide bonds
    only usefulness of inheriting from BondSelector is reuse of makeUniq
    """


    def __init__(self):
        self.OSel = BondElementBondSelector('O')
        self.CSel = BondElementBondSelector('C')
        self.NSel = BondElementBondSelector('N')
        #self.Sel3 = BondBranchBondSelector(3)
        self.lSel = LeafBondSelector()
        BondSelector.__init__(self)
        #can this be fixed??
        #criteria = [lambda x, hsel=self.BES: \
        #    len(x.atom1.bonds)-1==len(hsel(x.bonds)), \
        #    lambda x, hsel=self.BES:\
        #    len(x.atom2.bonds)-1==len(hsel(x.bonds))]


    def select(self, bnds):
        n_bnds = self.NSel.select(bnds)
        #bonds between Nitrogens and Carbons
        nc_bnds = self.CSel.select(n_bnds)
        ansBnds = BondSet()
        for b in nc_bnds:
            a0 = b.atom1
            a2 = b.atom2
            if a0.element=='N':
                a2 = b.atom1
                a0 = b.atom2
            if len(a0.bonds)==3:
                #get bonds from this carbon to oxygen
                o_bnds = self.OSel.select(BondSet(a0.bonds))
                if not len(o_bnds):
                    continue
                #get bonds from this carbon to leaf oxygen
                resBnds = self.lSel.select(o_bnds)
                if len(resBnds):
                    ansBnds.append(b)
        return self.makeUniq(ansBnds)
        


class GuanidiniumBondSelector(BondSelector):
    """ class that selects guanidinium bonds
    """


    def __init__(self):
        self.CSel = BondElementBondSelector('C')
        self.NSel = BondElementBondSelector('N')
        self.Sel3 = BondBranchBondSelector(3)
        BondSelector.__init__(self)
        #can this be fixed??
        #criteria = [lambda x, hsel=self.BES: \
        #    len(x.atom1.bonds)-1==len(hsel(x.bonds)), \
        #    lambda x, hsel=self.BES:\
        #    len(x.atom2.bonds)-1==len(hsel(x.bonds))]


    def select(self, bnds):
        #bonds of carbon bonded to 3 nitrogens
        c_bnds = self.CSel.select(bnds)
        c3_bnds = self.Sel3.select(c_bnds)
        #n_bnds = self.NSel.select(c3_bnds)
        #nb without lambda expression, this selected PRO bonds N-CA, N-CD
        #    and CYS C-N, THR C-N
        #c3_ats = AtomSet(filter(lambda x:x.element=='C', self.getAtoms(c3_bnds)))
        cats = self.getAtoms(c3_bnds)
        c3_ats = cats.get(lambda x: len(x.bonds)==3)
        cycle_bonds = CycleBondSelector().select(bnds)
        ansBnds = BondSet()
        for at in c3_ats:
            incycle_ct = 0  #keep track of bonds in cycles for vina tutorial ligand.pdb
            ok = True
            for b in at.bonds:
                at2 = b.neighborAtom(at)
                if at2.element!='N':
                    ok=0
                if b in cycle_bonds:
                    incycle_ct += 1
            if ok and incycle_ct<2: 
                #if 2 of 3 this carbon's bonds are in cycle, NOT guanidinium
                for b in at.bonds:
                    if b in bnds:
                        ansBnds.append(b)
        return self.makeUniq(ansBnds)



class PeptideBackBoneBondSelector(BondSelector):
    """ class that selects PeptideBackBone bonds
    only usefulness of inheriting from BondSelector is reuse of makeUniq
    """

    std_res_types = ['ILE', 'GLN', 'LYS', 'GLY', 'GLU', 
                     'CYS', 'ASP', 'HSD', 'HSE', 'HSP', 
                     'HID', 'HIE', 'ASN', 'HIP', 'VAL', 
                     'THR', 'HIS', 'TRP', 'SER', 'PHE', 
                     'PRO', 'ALA', 'MET', 'LEU', 'ARG', 
                     'TYR']

    def __init__(self):
        #peptide backbone bonds
        BondSelector.__init__(self)
        self.cycleBondSelector = CycleBondSelector(useMaxSize=0)
        #can this be fixed??
        #criteria = [lambda x, hsel=self.BES: \
        #    len(x.atom1.bonds)-1==len(hsel(x.bonds)), \
        #    lambda x, hsel=self.BES:\
        #    len(x.atom2.bonds)-1==len(hsel(x.bonds))]


    def select(self, bnds):
        ats = self.getAtoms(bnds)
        bbBnds = BondSet()
        bb_ats = ats.get(lambda x: x.parent.type in self.std_res_types and \
                                   x.name in ['N', 'C','CA'])
        if bb_ats:
            bbBnds = bb_ats.bonds[0]
        #remove any (small) cyclebonds first!!
        cycle_bnds = self.cycleBondSelector.select(bbBnds)
        bbBnds = bbBnds - cycle_bnds
        #print "found ", len(bbBnds), " PeptideBackBoneBonds"
        return self.makeUniq(bbBnds)



class CycleBondSelector(BondSelector):
    """ class that selects bonds in cycles
    only usefulness of inheriting from BondSelector is reuse of makeUniq
    detectAll can be used to detect rings not found by PyBabel.RingerFinder
    These are rings formed by bonds between atoms in separate cycles
    """
    def __init__(self, useMaxSize=1):
        BondSelector.__init__(self)
        self.useMaxSize=useMaxSize

    def select(self, bnds, detectAll=False):
        ats = self.getAtoms(bnds)
        rf = RingFinder()
        if self.useMaxSize:
            # increase maxSize of ring 
            rf.findRings2(ats, bnds, maxSize=len(ats))
        else:
            #limit maxSize of ring to default 20 
            rf.findRings2(ats, bnds)
        #check for connected cycles
        bonds= BondSet(rf.allRingBonds)
        if detectAll and len(rf.allRingBonds):
            ats = BondSet(rf.allRingBonds).getAtoms()
            ANS = {}
            for a in ats:
                for b in a.bonds:
                    if b.atom1 in ats and b.atom2 in ats:
                        ANS[b] = 1
            #return BondSet(rf.allRingBonds)
            bonds= BondSet(ANS.keys())
        return bonds
    


class AromaticCycleBondSelector(BondSelector):
    """ class that selects all bonds in aromatic cycles
        (uses PyBabel Aromatic class)
    """

    def select(self, bnds):
        ats = self.getAtoms(bnds)
        rf = RingFinder()
        #nb maxSize of ring is 20 by default
        rf.findRings2(ats, bnds, maxSize=len(ats))
        ats = self.getAtoms(rf.allRingBonds)
        allAts = ats.top.uniq().allAtoms
        atype = AtomHybridization()
        atype.assignHybridization(allAts)
        bo = BondOrder()
        bo.assignBondOrder(allAts, allAts.bonds[0], rf)
        arom = Aromatic(rf)
        arom.find_aromatic_atoms(ats)
        aromatic_bnds = ats.bonds[0].get(lambda x: x.bondOrder=='aromatic')
        aromatic_ats = self.getAtoms(aromatic_bnds)
        aromatic_ats.aromatic = 1
        return aromatic_bnds
        
 

class AromaticCycleBondSelector2(BondSelector):
    """ class that selects aromatic bonds in cycles
        (uses autotors calc of angle between adjacent normals)
    """


    def select(self, bnds, cutoff=7.5):
        cutoffValue = math.cos(cutoff*math.pi/180.)
        #print "cutoffValue=", cutoffValue
        aromaticCs = AtomSet([])
        atD = {}
        ctr = 1
        aromaticBnds = BondSet()

        #taken from autotorsCommands
        ats = self.getAtoms(bnds)
        rf = RingFinder()
        rf.findRings2(ats, bnds)
        cyclecount = rf.ringCount
        aromatic_cycles = []
        ct = 1
        for ring in rf.rings:
            blist = ring['bonds']
            for bnd in blist:
                at = bnd.atom1
                #next find the other bond in this cycle with at as one of the atoms:
                z2 = filter(lambda x:x!=bnd and x.atom1==at or x.atom2==at, blist)
                bnd.nextbond = z2[0]
                neighbor = z2[0].atom1
                if neighbor==at:
                    neighbor = z2[0].atom2
                bnd.next1 = neighbor
                #print "set ", bnd.atom1.name,'-', bnd.atom2.name,".next1 to ", neighbor.name
                #now each bond has 3 atoms specified for it: its own two and the next1 to atom1
                at1 = bnd.next1
                at2 = bnd.atom1
                at3 = bnd.atom2
                pt1 = at1.coords
                pt2 = at2.coords
                pt3 = at3.coords
                a1 = Numeric.subtract(pt2,pt1)
                b1 = Numeric.subtract(pt3,pt2)
                p = [a1[1]*b1[2]-a1[2]*b1[1],a1[2]*b1[0]-a1[0]*b1[2],a1[0]*b1[1]-a1[1]*b1[0]]
                result0 = Numeric.sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])
                bnd.nrmsize = result0
                result1 = p
                bnd.nrms = result1
                #next find the other bond w/atom2:
                z2 = filter(lambda x,bnd=bnd, at2=bnd.atom2, blist=blist:x!=bnd and x.atom1==at2 or x.atom2==at2, blist)
                #finally, return the other atom in this bond
                bnd.nextbond2 = z2[0]
                if bnd.nextbond2==bnd:
                    bnd.nextbond2 = z2[1]
                #neighbor2 = self._getnxtAtom(b.atom2,b.nextbond2)
                neighbor2 = bnd.nextbond2.atom1
                if neighbor2==bnd.atom2:
                    neighbor2 = bnd.nextbond2.atom2
                bnd.next2 = neighbor2
                #next have to check whether the normals are parallel
                #check each pair in each bond, how to keep track??
                #have to get normal at b's atom2: so have to get next2 for this bond:
            #have to loop twice to make sure neighbor has nrms
            for bnd in blist:
                p = bnd.nrms
                psize = bnd.nrmsize
                q = bnd.nextbond2.nrms
                qsize = bnd.nextbond2.nrmsize
                #theta is the critical test for planarity:
                #if angle between 2 nrms is 0, atoms are planar
                #NB>test is comparing theta,cos(angle), w/zero
                bnd.theta = Numeric.dot(p,q)/(psize*qsize)
            #NB: REPEAT STUFF FROM HERE TO THE END IF aromaticCutOff changes
            ct = 0
            for bnd in blist:
                #these are values for the default 7.5degrees:
                #if theta>=0.997 or theta<=-0.997:
                #print bnd.atom1.name, '-', bnd.atom2.name, '->', bnd.theta 
                if bnd.theta>=cutoffValue or bnd.theta<=-cutoffValue:
                    bnd.posAromatic = 1
                    ct = ct + 1
                else:
                    bnd.posAromatic = 0
                #print ' posAromatic=', bnd.posAromatic
            #after checking all the bonds in current cycle, compare #posAromatic w/number
            if ct==len(blist):
                #print ctr," cycle is aromatic"
                #NB: only C=C bonds are considered aromatic 
                #    could change the name and autodock_element here....
                for bnd in blist:
                    # THIS IS WRONG!!!
                    #if bnd.atom1.element=='C' and bnd.atom2.element=='C':
                    #    aromaticBnds.append(bnd)
                    aromaticBnds.append(bnd)
            ctr = ctr + 1
            #print "len(aromaticBnds)=", len(aromaticBnds)
        #for b in aromaticBnds:
        #    print b.atom1.name, '-', b.atom2.name
        return aromaticBnds
        


class PeptideAromaticCycleBondSelector(BondSelector):
    """ class that selects bonds in cycles of peptide residues

    """
    bondD = aromDict = {}
    bondD['TRP'] = ['NE1','CD1','CG','CD2','CE2','CZ2',\
                     'CH2','CZ3','CE3', 'AD1','AG','AD2',\
                     'AE2','AZ2', 'AH2','AZ3','AE3' ]
    bondD['TYR'] = ['CD1','CG','CD2','CE1','CE2','CZ',\
                    'AD1','AG','AD2','AE1','AE2','AZ']
    bondD['PHE'] = ['CD1','CG','CD2','CE1','CE2','CZ',\
                    'AD1','AG','AD2','AE1','AE2','AZ']
    bondD['HIS'] = ['CD2','CE1','CG','ND1','NE2', 
                    'AD2','AE1','AG']
    ##bondD['PRO'] = ['N','CD','CA','CB','CG'] ??why is this included??

    #pep_aromList = ['PHE_CD1', 'PHE_CG', 'PHE_CD2', 'PHE_CE1',\
    #        'PHE_CE2', 'PHE_CZ', 'TYR_CD1', 'TYR_CG', 'TYR_CD2', 'TYR_CE1',\
    #        'TYR_CE2', 'TYR_CZ', 'HIS_CD2', 'HIS_CE1', 'HIS_CG', 'TRP_CD1',\
    #        'TRP_CG', 'TRP_CD2', 'TRP_CE2', 'TRP_CZ2', 'TRP_CH2', 'TRP_CZ3',\
    #        'TRP_AE3', 'PHE_AD1', 'PHE_AG', 'PHE_AD2', 'PHE_AE1',\
    #        'PHE_AE2', 'PHE_AZ', 'TYR_AD1', 'TYR_AG', 'TYR_AD2', 'TYR_AE1',\
    #        'TYR_AE2', 'TYR_AZ', 'HIS_AD2', 'HIS_AE1', 'HIS_AG', 'TRP_AD1',\
    #        'TRP_AG', 'TRP_AD2', 'TRP_AE2', 'TRP_AZ2', 'TRP_AH2', 'TRP_AZ3',\
    #        'TRP_AE3' ]


    def select(self, bnds):
        resSet = self.getAtoms(bnds).parent.uniq()
        aromResSet = resSet.get(lambda x: x.type in self.aromDict.keys())
        if not aromResSet:
            return BondSet()
        bondDict = {}
        pep_arom = BondSet()
        for i in range(len(aromResSet)):
            res = aromResSet[i]
            keys = self.aromDict[res.type]
            res_bnds = bondDict[i+1] = res.atoms.get(lambda x, \
                                keys = keys:x.name in keys).bonds[0]

            for b in res_bnds:
                b.cyclenum = i + 1
                b.aromatic = 1  #???
                if b in bnds:
                    pep_arom.append(b)
        return pep_arom
          


class BondOrderBondSelector(BondSelector):
    """ class that selects bonds with bondOrder>1 

    only usefulness of inheriting from BondSelector is reuse of makeUniq
    """

    def select(self, bnds, bondOrder=1):
        ats = self.getAtoms(bnds)
        mols = ats.top.uniq()
        atype = AtomHybridization()
        for m in mols:
            if not m.chains[0].hasBonds:
                m.buildBondsByDistance()
            allAts = m.chains.residues.atoms
            atype.assignHybridization(allAts)
            rf = RingFinder()
            rf.findRings2(allAts, allAts.bonds[0])
            bo = BondOrder()
            bo.assignBondOrder(allAts, allAts.bonds[0], rf)
        return BondSet(bnds.get(lambda x, ord=bondOrder: x.bondOrder==ord))



class RotatableBondSelector(BondSelector):
    """ class that selects rotatable bonds, according to AutoDock:
    all bonds NOT in cycles and Not leafBonds BUT bondOrder==1
    only usefulness of inheriting from BondSelector is reuse of makeUniq
    """
    def __init__(self):
        BondSelector.__init__(self)


    def select(self, bnds, detectAll=False):
        # first pass: taken from autotorsCommands
        # bondOrder==1
        #THIS IS WRONG! 4/22: guanidiniums possible but not active
        #select the guanidiniums
        ####guanidiniumBonds = GuanidiniumBondSelector().select(bnds)
        #print "guanidiniumBonds=", len(guanidiniumBonds)
        ###possibleBnds = bnds.subtract(guanidiniumBonds)
        ###rotatable = BondOrderBondSelector().select(possibleBnds)
        rotatable = BondOrderBondSelector().select(bnds,1)

        #select the leaves 
        leafBonds = LeafBondSelector().select(rotatable)

        #select the cycles 
        cycleBonds = CycleBondSelector().select(rotatable,detectAll=detectAll)
        #remove the cycle bonds and return the residue
        #remove the bonds which are guanidinium bonds
        ###4/22rotatable = rotatable.subtract(guanidiniumBonds)
        #remove the leaves
        rotatable = rotatable.subtract(leafBonds)
        #remove the cycle bonds and return the residue
        rotatable = rotatable.subtract(cycleBonds)
        
        #print "rotatable=", len(rotatable)
        return rotatable
        
