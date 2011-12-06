############################################################################
#
# Author:  Ruth Huey
#
# Copyright: M. Sanner TSRI 2004
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/atomTypeTools.py,v 1.23 2010/09/22 21:44:17 rhuey Exp $
#
#

"""
This module implements classes to support AutoDock atomtypes
    including:
        NonpolarHydrogenMerger: a class which detects and then
            merges nonpolarHydrogens with the heavy atoms to 
            which each is bound.
        LonepairMerger: a class which detects and then
            merges lonepairs with the heavy atoms to 
            which each is bound.
        AromaticCarbonManager: a class used for managing 
            names + autodock_element fields of
            aromatic carbons in cycles which it judges
            flat enough according to its 'cutoff' parameter
            to be considered aromatic.
        SolvationParameterizer: class used to assign solvation 
            parameters 'AtSolPar' and 'AtVol' to each atom in 
            receptor molecules for autodock305 experiment

"""

import os
from MolKit.molecule import AtomSet, BondSet
from MolKit.bondSelector import AromaticCycleBondSelector2
from AutoDockTools.sol_par import solvs
from PyBabel.atomTypes import AtomHybridization




class NonpolarHydrogenMerger:
    """
    removes non-polarhydrogens from an atomset;
    WARNING: builds bonds if none have been built previously
"""


    def mergeNPHS(self, atoms, renumber=1):
        if len(atoms.bonds[0])==0:
            print "WARNING atoms have no bonds....BUILDING THEM!!!"
            tops = atoms.top.uniq()
            for t in tops: t.buildBondsByDistance()
        #MAKE sure there are some hydrogens
        hs = atoms.get(lambda x: x.element=='H')
        if hs is None or len(hs)==0:
            return 0
        #Check whether there are any hydrogens with no bonds
        no_bnd_hs = hs.get(lambda x: len(x.bonds)==0)
        if len(no_bnd_hs):
            print "Warning: hydrogens, ", no_bnd_hs.name, " , with no bonds!"
        if len(no_bnd_hs)==len(hs):
            return 0
        bonded_hs = hs.get(lambda x: len(x.bonds)>0)
        #next check is superfluous
        if bonded_hs is None or len(bonded_hs)==0:
            print "Warning: no hydrogens with bonds!"
            return 0
        #Check whether there are any nphs hydrogens
        nphs = bonded_hs.get(lambda x: x.bonds[0].atom1.element=='C' \
                         or x.bonds[0].atom2.element=='C')
        if nphs is None or len(nphs)==0: 
            return 0 #if there aren't any, just return 0
        #if there are nphs, merge the charges + remove nphs
        tops = nphs.top.uniq() 
        for t in tops:
            #t_nphs are the nonpolarHydrogens in this molecule, 't'
            t_nphs = nphs.get(lambda x: x.top==t)

            #clean up t.allAtoms shortcut
            t.allAtoms = t.allAtoms - t_nphs

            #deal with the charges, if there are any
            chList = t_nphs[0]._charges.keys()
            for h in t_nphs:
                chs = h._charges.keys() 
                for c in chList:
                    if c not in chs:
                        chList.remove(c)
            if not len(chList):
                print 'no charges on carbons to increment'
            else:
                for chargeSet in chList:
                    for h in t_nphs:
                        if len(h.bonds)==0:
                            print h.full_name(), ' has no bonds!'
                        else:
                            c_atom = h.bonds[0].atom1
                            if c_atom==h: 
                                c_atom = h.bonds[0].atom2
                            c_atom._charges[chargeSet] = c_atom.charge + h.charge

            #have to do this loop separately because there may not be charges
            len_nphs = len(nphs)
            for h in nphs:
                #b = at.bonds[0]
                #could have hydrogens with >1bond!!!(YIKES)
                for b in h.bonds:
                    c = b.atom1
                    if c==h:
                        c = b.atom2
                    c.bonds.remove(b)
                h.bonds = BondSet()
                h.parent.remove(h)
                del h 
            if renumber:
                lenAts = len(t.chains.residues.atoms)
                assert lenAts==len(t.allAtoms)
                t.allAtoms.number = range(1, lenAts+1) 
        return len_nphs



class LonepairMerger:
    """
    removes lone pair 'atoms' from an atomset;
    WARNING: builds bonds if none have been built previously
This module implements a class which merges lone pair 'atoms' with the heavy
atom to which each is bound.
"""


    def mergeLPS(self, atoms, renumber=1):
        if len(atoms.bonds[0])==0:
            print "WARNING atoms have no bonds....BUILDING THEM!!!"
            tops = atoms.top.uniq()
            for t in tops: t.buildBondsByDistance()
        #lps = atoms.get(lambda x: x.element=='Xx' and \
                        #(x.name[0]=='L' or x.name[1]=='L'))
        lps = atoms.get(lambda x: x.element=='Lp' or x.element=='lp' or \
                        (x.element=='Xx' and (x.name[0]=='L' or x.name[1]=='L')))
        if lps is None or len(lps)==0: 
            return [] #if there aren't any, just return empty list
        tops = lps.top.uniq() 
        for mol in tops:
            #t_lps are the nonpolarHydrogens in this molecule, 'mol'
            mol_lps = lps.get(lambda x: x.top==mol)

            #clean up t.allAtoms shortcut
            mol.allAtoms = mol.allAtoms - mol_lps

            #deal with the charges, if there are any
            chList = mol_lps[0]._charges.keys()
            for lp in mol_lps:
                chs = lp._charges.keys() 
                for c in chList:
                    if c not in chs:
                        chList.remove(c)
            if not len(chList):
                print 'no charges on carbons to increment'
            else:
                for chargeSet in chList:
                    for lp in mol_lps:
                        if len(lp.bonds)==0:
                            print lp.full_name(), ' has no bonds!'
                        else:
                            c_atom = lp.bonds[0].atom1
                            if c_atom==lp: 
                                c_atom = lp.bonds[0].atom2
                            c_atom._charges[chargeSet] = c_atom.charge + lp.charge

            #have to do this loop separately because there may not be charges
            len_lps = len(lps)
            for lp in lps:
                #b = at.bonds[0]
                #could have hydrogens with >1bond!!!(YIKES)
                for b in lp.bonds:
                    c = b.atom1
                    if c==lp:
                        c = b.atom2
                    c.bonds.remove(b)
                lp.bonds = BondSet()
                lp.parent.remove(lp)
                del lp 
            if renumber:
                lenAts = len(mol.chains.residues.atoms)
                assert lenAts==len(mol.allAtoms)
                mol.allAtoms.number = range(1, lenAts+1) 
        return len_lps
            


class AromaticCarbonManager:
    """
    Used on ligands for managing names + autodock_element fields of
    aromatic carbons in cycles
    The parameter 'rename' determines whether aromatic carbons are
    named starting with 'A': 'rename' is True for AutoDock3 ligands 
    and False for AutoDock4 ligands whose written files have added 'type' field.
    """

    def __init__(self, cutoff=7.5, rename=True):
        self.cutoff = 7.5
        self.aromBndSel = AromaticCycleBondSelector2()
        self.rename = rename


    def setAromaticCarbons(self, molecule, cutoff=None, debug=False):
        assert len(molecule.allAtoms.bonds[0])
        if cutoff: 
            if cutoff!=self.cutoff:
                old_aromCs = molecule.allAtoms.get(lambda x:\
                                (x.element=='C' and x.autodock_element=='A'))
                if old_aromCs is not None and len(old_aromCs)!=0:
                    if debug: print "resetting ", len(old_aromCs), " prior aromCs"
                    self.set_carbon_names(old_aromCs, 'C')
            self.cutoff = cutoff
        typed_atoms = molecule.allAtoms.get(lambda x: hasattr(x, 'autodock_element'))
        currentAromCs = AtomSet()
        if typed_atoms is not None and len(typed_atoms)==len(molecule.allAtoms):
            currentAromCs = molecule.allAtoms.get(lambda x: (x.element=='C' and x.autodock_element=='A'))
        if debug:
            if currentAromCs is not None and len(currentAromCs):
                print "now: ", len(currentAromCs)
            else:
                print "now: no aromCs"
        aromBnds = self.aromBndSel.select(molecule.allAtoms.bonds[0],
                                                  self.cutoff)
        aromBndAts = self.aromBndSel.getAtoms(aromBnds)
        result = AtomSet()
        changed = AtomSet()
        if len(aromBndAts):
            aromCs =  AtomSet(filter(lambda x: x.element=='C', \
                                              aromBndAts))           
            aromCs = aromCs.uniq()
            if len(aromCs)>len(result):
                result = aromCs
            if debug: 
                print "len(aromCs)=", len(aromCs)
            changed = self.set_carbon_names(aromCs, 'A')
        #if len(changed)>len(result):
        #    result = changed
        #return  result
        return changed


    def change_aromatic_criteria(self, cutoff):
        """
        cutoff: angle between adjacent normals
        reset Aromatic carbons and redetect with
        new cutoff value. Cutoff is maximum possible value for 
        angle between normals to adjacent atoms in a cycle
        for carbons to considered aromatic. All normals for a
        single cycle must pass for any of its carbons to be
        considered aromatic.  This is a measure of the flatness
        of the cycle.
        """
        self.cutoff = cutoff
        return 


    def set_carbon_names(self, atoms, type):
        #set carbon names explicitly
        if not atoms or not len(atoms):
            return "ERROR: set_carbon_names called with no atoms"
        assert type in ['C','A']
        if not hasattr(atoms, 'autodock_element'):
            atoms.autodock_element = atoms.element
        changed = AtomSet()
        for at in atoms:
            if at.element!='C': continue
            if at.autodock_element!=type:
                if self.rename:
                    if len(at.name)>1:
                        at.name = type + at.name[1:]
                    else:
                        at.name = type
                at.autodock_element = type
                changed.append(at)
        return changed




class SolvationParameterizer:
    """
    adds AtSolPar and AtVol to standard protein atoms
"""
    def __init__(self, bblist=['C','O','N','CA'], 
                       nuc_acid_list= ['  A','  C','  G']):
        self.bblist = bblist
        self.nuc_acid_list = nuc_acid_list
        self.solvsKeys = solvs.keys()


    def addParameters(self, atoms):
        notfound = []
        for at in atoms:
            #use shortcut for peptide backbone atoms
            if at.name in self.bblist:
                atKey = at.name + '---'
            elif (at.name=='C2*' or at.name=='C2\'') \
                    and at.parent.type in self.nuc_acid_list:
                #special treatment for variants of nucleic acid names
                childnames = at.parent.atoms.name
                if 'O2*' in childnames or 'O2\'' in childnames:
                    #treat rna variant even more specially
                    atKey = at.name + ' R' + at.parent.type[2]
                    print 'special rna key', atKey
                else:
                    atKey = at.name + at.parent.type
            else:
                atKey = at.name + at.parent.type

            #if there are solvation parameters for this atom, add them
            if atKey in self.solvsKeys:
                at.AtVol, at.AtSolPar = solvs[atKey]
            else:
                at.AtVol, at.AtSolPar = (0.00, 0.00)
                if at.element!='H':
                    notfound.append(at)
        return notfound


class AutoDock4_AtomTyper:
    """
    sets autodock_element according to AutoDock4 criteria
"""

    def __init__(self, renameAtoms=0, verbose=False, set_aromatic_carbons=True):
        self.renameAtoms = 0 # whether to change first 
                             # char of carbon atom names to 'A' for AD4 type
        self.set_aromatic_carbons = set_aromatic_carbons
        if set_aromatic_carbons:
            #aromatic_carbon_manager:  
            # disable for a huge protein by setting typeAtoms=0
            self.acm = AromaticCarbonManager(rename=False)
        self.verbose = verbose
        self.acceptorList = ['NA','OA','SA']
        self.pep_aromList  = [
             'PHE_CD1', 'PHE_CG', 'PHE_CD2', 'PHE_CE1', 'PHE_CE2', 'PHE_CZ', 
             'PHE_AD1', 'PHE_AG', 'PHE_AD2', 'PHE_AE1', 'PHE_AE2', 'PHE_AZ', 
             'TYR_CD1', 'TYR_CG', 'TYR_CD2', 'TYR_CE1', 'TYR_CE2', 'TYR_CZ', 
             'TYR_AD1', 'TYR_AG', 'TYR_AD2', 'TYR_AE1', 'TYR_AE2', 'TYR_AZ', 
             'TRP_CG', 'TRP_CD1','TRP_CD2', 'TRP_CE2', 'TRP_CE3','TRP_CZ2', 'TRP_CH2', 'TRP_CZ3',
             'TRP_AG', 'TRP_AD1','TRP_AD2', 'TRP_AE2', 'TRP_AE3','TRP_AZ2', 'TRP_AH2', 'TRP_AZ3',
             'HIS_AD2', 'HIS_AE1', 'HIS_AG', 
             'HIS_CD2', 'HIS_CE1', 'HIS_CG'] 
 
     
    def setAutoDockElements(self, mol, typeAtoms=0, reassign=False, splitAcceptors=False):
        if os.path.splitext(mol.parser.filename)[-1]=='.pdbqt' and reassign is False:
            if self.verbose:
                print 'setAutoDockElements unnecessary:\n', mol.name, ' already has AD4 atomtypes: not reassigned!'
            return
        #check that the types have already been set, first
        if self.set_aromatic_carbons:
            #this sets the autodock_element for the planar cyclic carbons to 'A'
            aromatic_carbons = self.acm.setAromaticCarbons(mol)
        d = {}
        ah = AtomHybridization()
        if not len(mol.allAtoms.bonds[0]):
            mol.buildBondsByDistance()
        ats_with_babel_type = filter(lambda x: hasattr(x, 'babel_type'), mol.allAtoms)
        if typeAtoms or len(ats_with_babel_type)!=len(mol.allAtoms):
            if self.verbose:
                print "assigning babel_types"
            ah = AtomHybridization()
            ah.assignHybridization(mol.allAtoms)
        for item in mol.allAtoms:
            if not hasattr(item, 'autodock_element'):
                #item.element = item.element.upper()
                item.autodock_element = item.element
            if item.element=='H':
                if len(item.bonds):
                    n = item.bonds[0].neighborAtom(item)
                    if n.element not in ['C','A']:
                        item.autodock_element = 'HD'
                else:
                    #note: if there are no bonds, assume HD 
                    item.autodock_element = 'HD'
            elif item.element=='N':
                #this has already been done at beginning of for loop
                #if item.babel_type in ['N3+','Ntr', 'Nox']:
                #    item.autodock_element = 'N'
                if item.babel_type in ['N3','N2','N1']:
                    item.autodock_element = item.element + 'A'
                elif item.babel_type in ['Nam', 'Npl', 'Ng+']:
                    if len(item.bonds)==2:
                        item.autodock_element = item.element + 'A'
            elif item.element=='O':
                item.autodock_element = item.element + 'A'
            elif item.element=='S':
                #possible sulphur babel_types->['S3+','S3','S2','Sac','Sox','S']
                if item.babel_type not in ['Sox', 'Sac']:
                    item.autodock_element = item.element + 'A'
            elif item.element=='C':
                if item.parent.type+'_' + item.name in self.pep_aromList:
                    item.autodock_element = 'A'
                    if self.renameAtoms:  # 
                        if len(item.name)==1:
                            item.name = 'A'
                        else:
                            item.name = 'A' + item.name[1:]
                elif hasattr(item, 'autodock_element'):
                    item.autodock_element = item.autodock_element
            if splitAcceptors and item.autodock_element in self.acceptorList:
                #if item has a bond to a hydrogen:
                #   change acceptor type 'XA' to 'XB' for 'B'oth 'A'cceptor and 'D'onor 
                if 'H' in item.bonds.getAtoms().element:
                    #O->OB;N->NB;S->SB
                    item.autodock_element = item.element + 'B'
            d[item.autodock_element] = 1
        type_list = d.keys()
        type_list.sort()
        mol.types = type_list


if __name__ == "__main__":
    from MolKit import Read
    m = Read("/mgl/work4/rhuey/dev23/hsg1_no_conects.pdbq")[0]
    solP = SolvationParameterizer()
    ll = solP.addParameters(m.allAtoms)
    #NB: hydrogens get 0.0 0.0
    print "len(notfound) = ", len(ll)
    for i in ll: 
        print i.full_name(),
    from MolKit.pdbWriter import PdbqsWriter
    writer = PdbqsWriter()
    writer.write("/mgl/work4/rhuey/dev23/test_hsg1.pdbqs", m.allAtoms)

