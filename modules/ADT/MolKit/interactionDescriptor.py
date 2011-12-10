## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

############################################################################
#
# Author:  Ruth Huey
#
# Copyright: M. Sanner TSRI 2007
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/MolKit/interactionDescriptor.py,v 1.22 2010/05/16 18:08:17 rhuey Exp $
#
# $Id: interactionDescriptor.py,v 1.22 2010/05/16 18:08:17 rhuey Exp $
#

"""
This module implements the InteractionDescriptor class which builds sets of
atoms in close contact and hydrogen bonds for two Molecules
"""
from molecule import Atom, AtomSet, HydrogenBond, BondSet
from protein import Residue, ResidueSet
from distanceSelector import CloserThanVDWSelector, DistanceSelector
from hydrogenBondBuilder import HydrogenBondBuilder
from PyBabel.cycle import RingFinder
from PyBabel.aromatic import Aromatic
from bondSelector import AromaticCycleBondSelector
import numpy.oldnumeric as Numeric
from mglutil.math import crossProduct


class InteractionDescriptor:
    """
    object which can detect atoms in close contact and build hydrogen bonds between atoms according
    to their coords and atom type for two sets 
    """

    def __init__(self, lig, macro, percentCutoff=1.0, detect_pi=False, dist_cutoff=6., include_metal_cations=True):
        self.lig_atoms = lig.findType(Atom)
        self.lig = self.lig_atoms[0].top
        self.macro_atoms = macro.findType(Atom)
        self.macro = self.macro_atoms[0].top
        self.percentCutoff = percentCutoff
        self.distanceSelector = CloserThanVDWSelector(return_dist=0)
        self.hydrogen_bond_builder = HydrogenBondBuilder()
        self.distanceSelectorWithCutoff = DistanceSelector()
        self.dist_cutoff=float(dist_cutoff)
        self.include_metal_cations=include_metal_cations
        self.build(detect_pi=detect_pi)



    def build(self, percentCutoff=None, detect_pi=False):
        if not percentCutoff:
            percentCutoff = self.percentCutoff
        # first detect sets of atoms forming hydrogen bonds
        self.buildHydrogenBonds()           #
        # detect sets of atoms in close contact 
        # and detect sets of atoms in close contact not forming hydrogen bonds
        self.buildCloseContactAtoms(percentCutoff)              #
        # detect sequences of >3 contiguous residues which have atoms in close contact
        self.buildContiguousCloseResidueSequences()
        if detect_pi:
            self.detectPiInteractions()


    def buildCloseContactAtoms(self, percentCutoff):
        pairDict = self.distanceSelector.select(self.lig_atoms, 
                        self.macro_atoms, percentCutoff=percentCutoff)
        self.pairDict = pairDict
        #reset here
        lig_close_ats = AtomSet()
        macro_close_ats = AtomSet()
        closeAtoms = AtomSet()  #both sets
        cdict = {}
        for k,v in pairDict.items():
            if len(v):
                cdict[k] = 1
            for at in v:
                if at not in macro_close_ats:
                    cdict[at] = 1
        closeAtoms = AtomSet(cdict.keys())
        
        #macro_close_ats = closeAtoms.get(lambda x: x.top==self.macro)
        #lig_close_ats = closeAtoms.get(lambda x: x.top==self.lig)
        lig_close_ats = closeAtoms.get(lambda x: x in self.lig_atoms)
        macro_close_ats = closeAtoms.get(lambda x: x in self.macro_atoms)
        rdict = self.results
        rdict['lig_close_atoms'] = lig_close_ats
        rdict['lig_close_res'] = lig_close_ats.parent.uniq()
        rdict['lig_close_carbons'] = lig_close_ats.get(lambda x: x.element=='C')
        rdict['lig_close_non_hb'] = lig_close_ats - rdict['lig_hb_atoms']
        rdict['macro_close_atoms'] = macro_close_ats
        rdict['macro_close_res'] = ResidueSet(macro_close_ats.parent.uniq())
        rdict['macro_close_carbons'] = macro_close_ats.get(lambda x: x.element=='C')
        rdict['macro_close_non_hb'] = macro_close_ats - rdict['macro_hb_atoms']
        #deprecate this
        rdict['macro_close_only'] = macro_close_ats - rdict['macro_hb_atoms']


    def buildHydrogenBonds(self):
        self.results = d = {}
        h_pairDict = self.hydrogen_bond_builder.build(self.lig_atoms, self.macro_atoms)
        self.h_pairDict = h_pairDict
        #keys should be from lig, values from macro 
        #sometimes are not...@@check this@@
        h_results = {}
        for k, v in h_pairDict.items():
            h_results[k] = 1
            for at in v:
                h_results[at] = 1
        all_hb_ats = AtomSet(h_results.keys())  #all
        macro_hb_ats = d['macro_hb_atoms'] = all_hb_ats.get(lambda x: x.top==self.macro)
        # process lig
        lig_hb_res = d['lig_hb_res'] = ResidueSet()
        lig_hb_sidechains = d['lig_hb_sidechains'] = AtomSet()
        lig_gly_atoms = AtomSet()
        lig_hb_ats = d['lig_hb_atoms'] = all_hb_ats.get(lambda x: x in self.lig_atoms)
        if len(lig_hb_ats):
            d['lig_hb_res'] = lig_hb_res = lig_hb_ats.parent.uniq()
            d['lig_hb_sidechains'] = lig_hb_sidechains = lig_hb_res.atoms.get('sidechain')
            #to visualize hbonding involving GLY residues which have no side chains, show backbone atoms
            lig_gly_res = d['lig_hb_gly_res'] = lig_hb_res.get(lambda x: x.type=='GLY')
            if len(lig_gly_res):
                lig_gly_atoms = lig_gly_res.atoms
        # build extended set of hbonding_atoms_to_show as lines, just in case
        lig_hbas = AtomSet(lig_hb_sidechains + lig_gly_atoms + lig_hb_ats) #all from lig
        extraAts = AtomSet()
        for at in lig_hbas:
            for b in at.bonds:
                at2 = b.atom1
                if at2==at: 
                    at2 = b.atom2
                #add it to the atomset
                if at2 not in lig_hbas:
                    extraAts.append(at2)
        if len(lig_hbas):
            for at in extraAts:
                lig_hbas.append(at)
        d['lig_hbas'] = lig_hbas
        # process macro
        macro_hb_res =  ResidueSet()
        d['macro_hb_res'] = macro_hb_res
        d['macro_hb_sidechains'] = AtomSet()
        d['macro_hb_gly_res'] = ResidueSet()
        if len(macro_hb_ats):
            macro_hb_res = macro_hb_ats.parent.uniq()
        #4. display sidechains of hbonding residues as sticksNballs
            macro_hb_sidechains = d['macro_hb_sidechains'] = macro_hb_res.atoms.get('sidechain')
            macro_hb_gly_res = d['macro_hb_gly_res'] = macro_hb_res.get(lambda x: x.type=='GLY')
        macro_hb_gly_res = ResidueSet()
        macro_hb_gly_atoms = AtomSet()
        if len(macro_hb_gly_res): 
            macro_hb_gly_atoms = macro_hb_gly_res.atoms
        d['macro_hb_gly_atoms'] = macro_hb_gly_atoms
        # build extended set of hbonding_atoms_to_show as lines
        macro_hbas = d['macro_hbas'] = AtomSet()
        if len(macro_hb_ats):
            macro_hbas = d['macro_hbas'] = AtomSet(macro_hb_sidechains + macro_hb_gly_atoms + macro_hb_ats) #all from macro
        #add atoms bonded to hb atoms to make lines displayed more reasonable
        extraAts = AtomSet()
        for at in macro_hbas:
            for b in at.bonds:
                at2 = b.atom1
                if at2==at: 
                    at2 = b.atom2
                #add it to the atomset
                if at2 not in macro_hbas:
                    extraAts.append(at2)
        if len(macro_hbas):
            for at in extraAts:
                macro_hbas.append(at)
        d['hbas_macro'] = macro_hbas


    def buildContiguousCloseResidueSequences(self):
        #7. attempt to show ribbon for contiguous residues in macromolecule
        rdict = self.results
        res = rdict['macro_close_res']
        chs = res.parent.uniq()
        ss_res = ResidueSet()
        last_ind = 0
        chain1 = 1
        #output = 0
        for c in chs:
            num_res = len(c.residues)
            if num_res <3:
                continue
            rr = res.get(lambda x: x.parent==c)
            rr.sort()
            chain1 = 0
            current_seq = ResidueSet() # contiguous residues
            current_set = ResidueSet() # all pieces in this chain
            skipped_set = ResidueSet() # hole in current contiguous piece
            if len(rr)>3:  #?? min num residues for ss:at least 3??
                #reset all
                first = c.residues.index(rr[0])
                last = c.residues.index(rr[-1])
                for r in c.residues[first:last+1]:
                    if r==c.residues[-1]:
                        if r in rr: 
                            if len(current_seq)>3:
                                current_seq.append(r)
                        if len(current_seq)>4:
                            ss_res.extend(current_seq)
                    if r not in rr:  #process hole
                        skipped_set.append(r)   # one hole ok
                        if len(skipped_set)>=2: # found second hole -> end seq
                            if len(current_seq)>4:
                                if not len(current_set):
                                    current_set = current_seq[:]
                                    current_set.sort()
                                else:
                                    current_set.extend(current_seq)
                                    current_set.sort()
                                if not len(ss_res):
                                    ss_res = current_set[:]
                                else:
                                    ss_res.extend(current_set)
                            skipped_set = ResidueSet()
                            current_seq = ResidueSet()
                    else:
                        #reset skipped_set
                        if len(skipped_set)<=1 and len(current_seq)>=1: #save RR_R 
                            current_seq.extend(skipped_set) #save hole if there is one
                            current_seq.append(r) #save this residue
                        else:  #just save it
                            current_seq.append(r)
                        skipped_set= ResidueSet()
            if len(current_seq)>4:
                for r in current_seq:
                    if r not in ss_res:
                        ss_res.append(r)
            if len(current_set)>4:
                for r in current_set:
                    if r not in ss_res:
                        ss_res.append(r)
        rdict['ss_res'] = ss_res


    def getCations(self, atoms):
        #select atoms in ARG and LYS residues
        arg_cations = atoms.get(lambda x: (x.parent.type=='ARG' and \
                                x.name in ['CZ']))
        lys_cations = atoms.get(lambda x: (x.parent.type=='LYS' and \
                                x.name in ['NZ', 'HZ1', 'HZ2', 'HZ3']))
        #select any positively-charged metal ions... cannot include CA here
        metal_cations = atoms.get(lambda x: x.name in ['Mn','MN', 'Mg',\
                                'MG', 'FE', 'Fe', 'Zn', 'ZN'])
        ca_cations = atoms.get(lambda x: x.name in ['CA', 'Ca'] and x.parent.type=='CA')
        cations = AtomSet() 
        #cations.extend(arg_cations)
        for a in arg_cations:
            cations.append(a)
        #cations.extend(lys_cations)
        for a in lys_cations:
            cations.append(a)
        #cations.extend(metal_cations)
        # including metal_cations and calcium optional
        if self.include_metal_cations:
            for a in metal_cations:
                cations.append(a)
            #cations.extend(ca_cations)
            for a in ca_cations:
                cations.append(a)
        return cations


    def detectPiInteractions(self, tolerance=0.95, debug=False, use_all_cycles=False):
        if debug: print "in detectPiInteractions"
        self.results['pi_pi'] = []        #stacked rings...?
        self.results['t_shaped'] = []     #one ring perpendicular to the other
        self.results['cation_pi'] = []    #
        self.results['pi_cation'] = []    #
        self.results['macro_cations'] = []#
        self.results['lig_cations'] = []  #
        #at this point have self.results
        lig_res = self.results['lig_close_res']
        if not len(lig_res):
            return
        lig_atoms = lig_res.atoms
        macro_res = self.results['macro_close_res']
        if not len(macro_res):
            return
        macro_atoms = macro_res.atoms
        l_rf = RingFinder()
        #Ligand
        l_rf.findRings2(lig_res.atoms, lig_res.atoms.bonds[0])
        #rf.rings is list of dictionaries, one per ring, with keys 'bonds' and 'atoms'
        if debug: print "LIG: len(l_rf.rings)=", len(l_rf.rings)
        if not len(l_rf.rings):
            if debug: print "no lig rings found by l_rf!"
            return
        acbs = AromaticCycleBondSelector()
        lig_rings = []
        for r in l_rf.rings:
            ring_bnds = r['bonds']
            if use_all_cycles: 
                lig_rings.append(ring_bnds)
            else:
                arom_bnds = acbs.select(ring_bnds)
                if len(arom_bnds)>4:
                    lig_rings.append(arom_bnds)
        if debug: print "LIG: len(lig_arom_rings)=", len(lig_rings)
        self.results['lig_rings'] = lig_rings
        self.results['lig_ring_atoms'] = AtomSet()
        #only check for pi-cation if lig_rings exist
        if len(lig_rings):
            macro_cations = self.results['macro_cations'] = self.getCations(macro_atoms)
            lig_ring_atoms = AtomSet()
            u = {}
            for r in lig_rings:
                for a in BondSet(r).getAtoms():
                    u[a] = 1
            if len(u): 
                lig_ring_atoms = AtomSet(u.keys())
                lig_ring_atoms.sort()
                self.results['lig_ring_atoms'] = lig_ring_atoms
            if len(macro_cations):
                if debug: print "check distances from lig_rings to macro_cations here"
                #macro cations->lig rings
                pairDict2 = self.distanceSelector.select(lig_ring_atoms,macro_cations)
                z = {}
                for key,v in pairDict2.items():
                    val = v.tolist()[0]
                    if val in macro_cations:
                        z[val] = [key]
                if len(z):
                    self.results['pi_cation'] = (z.items())
                else:
                    self.results['pi_cation'] = []
        #check the distance between the rings and the macro_cations
        self.results['lig_cations'] = self.getCations(lig_atoms)
        #Macromolecule
        m_rf = RingFinder()
        m_rf.findRings2(macro_res.atoms, macro_res.atoms.bonds[0])
        #rf.rings is list of dictionaries, one per ring, with keys 'bonds' and 'atoms'
        if debug: print "MACRO: len(m_rf.rings)=", len(m_rf.rings)
        if not len(m_rf.rings):
            if debug: print "no macro rings found by m_rf!"
            return
        macro_rings = []
        for r in m_rf.rings:
            ring_bnds = r['bonds']
            if use_all_cycles: 
                macro_rings.append(ring_bnds)
            else:
                arom_bnds = acbs.select(ring_bnds)
                if len(arom_bnds)>4:
                    macro_rings.append(arom_bnds)
        if debug: print "len(macro_arom_rings)=", len(macro_rings)
        self.results['macro_rings'] = macro_rings
        self.results['macro_ring_atoms'] = AtomSet()
        #only check for pi-cation if macro_rings exist
        if len(macro_rings):
            lig_cations = self.results['lig_cations'] = self.getCations(lig_atoms)
            macro_ring_atoms = AtomSet()
            u = {}
            for r in macro_rings:
                for a in BondSet(r).getAtoms(): #new method of bondSets
                    u[a] = 1
            if len(u):
                macro_ring_atoms = AtomSet(u.keys())
                macro_ring_atoms.sort()
                self.results['macro_ring_atoms'] = macro_ring_atoms
            if len(lig_cations):
                if debug: print "check distances from macro_rings to lig_cations here"
                pairDict3 = self.distanceSelector.select(macro_ring_atoms,lig_cations)
                z = {}
                for x in pairDict3.items():
                    #lig cations->macro rings
                    z.setdefault(x[1].tolist()[0], []).append(x[0])
                if len(z):
                    self.results['cation_pi'] = (z.items())
                else:
                    self.results['cation_pi'] = []
                #macro_pi_atoms = AtomSet(pairDict3.keys())
                #l_cations = AtomSet()
                #for v in pairDict3.values():
                #    for x in v:
                #        l_cations.append(x)
                #self.results['cation_pi'] = pairDict3.items()
                #self.results['cation_pi'] = (l_cations, macro_pi_atoms)
        #check for intermol distance <6 Angstrom (J.ComputChem 29:275-279, 2009)
        #compare each lig_ring vs each macro_ring
        for lig_ring_bnds in lig_rings:
            lig_atoms = acbs.getAtoms(lig_ring_bnds)
            lig_atoms.sort()
            if debug: print "len(lig_atoms)=", len(lig_atoms)
            #---------------------------------
            # compute the normal to lig ring
            #---------------------------------
            a1 = Numeric.array(lig_atoms[0].coords)
            a2 = Numeric.array(lig_atoms[2].coords)
            a3 = Numeric.array(lig_atoms[4].coords)
            if debug: print "a1,a2, a3=", a1.tolist(), a2.tolist(), a3.tolist()
            for macro_ring_bnds in macro_rings:
                macro_atoms = acbs.getAtoms(macro_ring_bnds)
                macro_atoms.sort()
                if debug: print "len(macro_atoms)=", len(macro_atoms)
                pD_dist = self.distanceSelectorWithCutoff.select(macro_ring_atoms, lig_atoms, cutoff=self.dist_cutoff)
                if not len(pD_dist[0]):
                    if debug: 
                        print "skipping ligand ring ", lig_rings.index(lig_ring_bnds), " vs ",
                        print "macro ring", macro_rings.index(macro_ring_bnds)
                    continue
                #---------------------------------
                # compute the normal to macro ring
                #---------------------------------
                b1 = Numeric.array(macro_atoms[0].coords)
                b2 = Numeric.array(macro_atoms[2].coords)
                b3 = Numeric.array(macro_atoms[4].coords)
                if debug: print "b1,b2, b3=", b1.tolist(), b2.tolist(), b3.tolist()
                # check for stacking 
                a2_1 = a2-a1
                a3_1 = a3-a1
                b2_1 = b2-b1
                b3_1 = b3-b1
                if debug: print "a2_1 = ", a2-a1
                if debug: print "a3_1 = ", a3-a1
                if debug: print "b2_1 = ", b2-b1
                if debug: print "b3_1 = ", b3-b1
                n1 = crossProduct(a3_1,a2_1) #to get the normal for the first ring
                n2 = crossProduct(b3_1,b2_1) #to get the normal for the second ring
                if debug: print "n1=", n1
                if debug: print "n2=", n2
                n1 = Numeric.array(n1)
                n2 = Numeric.array(n2)
                n1_dot_n2 = Numeric.dot(n1,n2)
                if debug: print "n1_dot_n2", Numeric.dot(n1,n2)
                if abs(n1_dot_n2) >= 1*tolerance: 
                    if debug: print "The rings are stacked vertically" 
                    new_result = (acbs.getAtoms(lig_ring_bnds), acbs.getAtoms(macro_ring_bnds))
                    self.results['pi_pi'].append(new_result)
                if abs(n1_dot_n2) <= 0.01*tolerance: 
                    if debug: print "The rings are stacked perpendicularly" 
                    new_result = (acbs.getAtoms(lig_ring_bnds), acbs.getAtoms(macro_ring_bnds))
                    self.results['t_shaped'].append(new_result)


    def print_ligand_residue_contacts(self, print_ctr=1):
        ctr = 1
        #pairDict is atom-based 
        #for residue-based report:
        # need to build lists of unique parents of keys 
        # and lists of unique parents of corresponding values
        res_d = {}
        for at in self.pairDict.keys():
            if at.parent not in res_d.keys():
                res_d[at.parent] = {}
            for close_at in self.pairDict[at]:
                res_d[at.parent][close_at.parent] = 1
        #print it out
        for lig_res in res_d.keys():
            if print_ctr:
                print ctr, lig_res.parent.name+':'+ lig_res.name + '->',
            else:
                print lig_res.parent.name+':'+ lig_res.name + '->',
            for macro_res in res_d[lig_res]:
                print macro_res.parent.name + ':' + macro_res.name + ',',
            print
            ctr += 1
        return res_d

    def print_macro_residue_contacts(self, print_ctr=1):
        ctr = 1
        #pairDict is atom-based 
        #for residue-based report:
        # need to build lists of unique parents of keys 
        # and lists of unique parents of corresponding values
        res_d = {}
        for at_key, at_list in self.pairDict.items():
            for at in at_list:
                if at.parent not in res_d.keys():
                    res_d[at.parent] = {}
                res_d[at.parent][at_key.parent] = 1
        #print it out
        for macro_res in res_d.keys():
            if print_ctr:
                print ctr, macro_res.parent.name+':'+ macro_res.name + '->',
            else:
                print macro_res.parent.name+':'+ macro_res.name + '->',
            for lig_res in res_d[macro_res]:
                print lig_res.parent.name + ':' + lig_res.name + ',',
            print
            ctr += 1
        return res_d

    def print_report(self, keylist=[]):
        if not len(keylist):
            keylist = [
                'lig_close_atoms',
                'lig_close_res',
                'lig_close_non_hb',
                'lig_close_carbons',
                'lig_hb_atoms',
                'lig_hb_res',
                'lig_hb_sidechains',
                'lig_hbas',
                'macro_close_atoms',
                'macro_close_res',
                'macro_close_non_hb',
                'macro_hb_atoms',
                'macro_hb_res',
                'macro_hb_sidechains',
                'macro_hbas',
                'ss_res',
                    ]
        d = self.results
        for k in keylist:
            print k, ':', len(d[k]), '-', d[k].__class__

    def print_hb_residue(self, print_ctr=1):
        ctr = 1
        #pairDict is atom-based 
        #for residue-based report:
        # need to build lists of unique parents of keys 
        # and lists of unique parents of corresponding values
        res_d = {}
        for at in self.h_pairDict.keys():
            if at.parent not in res_d.keys():
                res_d[at.parent] = {}
            for close_at in self.h_pairDict[at]:
                res_d[at.parent][close_at.parent] = 1
        #print it out
        # Hbond instance are define with donAt - accAtt
        for don_res in res_d.keys():
            if print_ctr:
                print ctr, don_res.top.name+':'+don_res.parent.name+':'+ don_res.name + '->',
            else:
                print don_res.top.name+':'+don_res.parent.name+':'+ don_res.name + '->',
            for acc_res in res_d[don_res]:
                print acc_res.top.name+':'+acc_res.parent.name + ':' + acc_res.name + ',',
            print
            ctr += 1
        return res_d
