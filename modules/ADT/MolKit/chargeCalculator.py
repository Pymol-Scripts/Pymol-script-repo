############################################################################
#
# Author:  Ruth Huey
#
# Copyright: M. Sanner TSRI 2004
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/MolKit/chargeCalculator.py,v 1.5 2004/08/11 15:17:53 rhuey Exp $
#
# $Id: chargeCalculator.py,v 1.5 2004/08/11 15:17:53 rhuey Exp $
#

"""
This module implements the ChargeCalculators classes which addCharges of 
various sorts to AtomSets.
"""

from MolKit.molecule import Atom, AtomSet, BondSet
from PyBabel.atomTypes import AtomHybridization
from PyBabel.gasteiger import Gasteiger
#this is for future use
#from AutoDockTools.observer import Subject
#FIX THIS: it shouldn't be in AutoDockTools
from MolKit.data.qkollua import q


class ChargeCalculator:
    """Base Class for computing per atom partial charges, adding entry to 
atom._charge dictionary and setting the atom.chargeSet to self.chtype
    """

    def __init__(self, chtype=''):
        self.chtype = ''


    def addCharges(self, atoms):
        #method to override
        pass



class GasteigerChargeCalculator(ChargeCalculator):
    """Class to compute per atom partial charges, add entry to 
atom._charge dictionary and to set the atom.chargeSet to self.chtype

TWO KNOWN PROBLEMS with gasteiger charges when added to molecules with
hydrogens and 2 c-terminal oxygens:  
1.N-terminus and C-terminus each has +.54 too much
charge. If there are 3 HN, each is set to .257 (down from .275); for pro there
are 2 HN and each is set to .248.  If there are 2 carboxyterminus oxygens (in
lpdb protein set, these are renamed OXT and O1...to fix OCT1 turning into
element Oc), each is set to .619 (from .646)
2.excess charge on single PO4s is +.500.... so -.125 is added to charge of
each of the oxygens bonded to the phosphorus atom.
    """

    def __init__(self):
        ChargeCalculator(self)
        self.chtype = 'gasteiger'


    def addCharges(self, atoms):
        """
        compute gasteiger charges, add them to atoms, return list of new
        charges
        NB: assumes MolKit atoms have bonds
        NB: calls type atoms if all the atoms haven't been babel typed. 
        if babel_types have been set by hand, they will be ok IF every
        atom has a babel type

        """
        #print 'in computeGasteiger'
        #check that atoms have babel_types
        w_babel_types = atoms.get(lambda x: hasattr(x, 'babel_type'))
        if w_babel_types is None or len(w_babel_types)!=len(atoms):
            babel = AtomHybridization()
            #have to have bonds for this
            babel.assignHybridization(atoms)
        Gast = Gasteiger()
        Gast.compute(atoms)
        gastCharges = []
        for c in atoms.gast_charge:
            gastCharges.append(round(c, 4))
        #NB THIS adds entry 'gasteiger' to atom._charges dict
        atoms.addCharges('gasteiger', gastCharges)
        atoms.chargeSet = 'gasteiger'
        #clean-up
        delattr(atoms, 'gast_charge')
        return gastCharges
            



class KollmanChargeCalculator(ChargeCalculator):
    """Class to compute per atom partial charges, add entry to 
atom._charge dictionary and to set the atom.chargeSet to self.chtype
    """

    def __init__(self):
        ChargeCalculator(self)
        self.chtype = 'Kollman'
        self.q = q


    def addCharges(self, atoms):
        #print 'in KollmanChargeCalculator.addCharges'
        for a in atoms:
            if not self.q.has_key(a.parent.type):
                a._charges['Kollman'] = 0.0
            else:
                dict = self.q[a.parent.type]
                if dict.has_key(a.name):
                    a._charges['Kollman'] = dict[a.name]
                else:
                    a._charges['Kollman'] = 0.0
            a.chargeSet = 'Kollman'
        #truncated here... should
        allRes = atoms.parent.uniq()
        #fix histidines
        hisRes = allRes.get(lambda x:x.name[:3]=='HIS')
        if hisRes is not None and len(hisRes):
            for hres in hisRes:
                self.fixHisRes(hres)
        #fix cysteines
        cysRes = allRes.get(lambda x:x.name[:3]=='CYS')
        if cysRes is not None and len(cysRes):
            for cres in cysRes:
                self.fixCysRes(cres)
        #fix termini
        allChains = allRes.parent.uniq()
        for ch in allChains:
            if ch.residues[0] in allRes:
                ##fix the charge on Hs 
                self.fixNterminus(ch.residues[0])
            if ch.residues[-1] in allRes:
                self.fixCterminus(ch.residues[-1])
        #make Kollman the current charge
        self.chargeSet = 'Kollman'
        return atoms.charge


    def fixNterminus(self, nres):
        """newTotal<-fixNterminu(nres)
            nres is the N-terminal residue
        """
        nresSet =  nres.atoms.get('N')
        if nresSet is None or not len(nresSet):
            if self.q.has_key(nres.type):
                for at in nres.atoms:
                    try:
                        at._charges['Kollman'] = self.q[at.parent.type][at.name]
                    except:
                        at._charges['Kollman'] = 0.0
            else:
                print  nres.name, ' not in qkollua; charges set to 0.0'
                for at in nres.atoms:
                    at._charges['Kollman'] = 0.0
            return 
        Natom = nres.atoms.get('N')[0]
        caresSet = nres.atoms.get('CA')
        if caresSet is None or not len(caresSet):
            if self.q.has_key(nres.type):
                for at in nres.atoms:
                    at._charges['Kollman'] = self.q[at.parent.type][at.name]
            else:
                print  nres.name, ' not in qkollua; charges set to 0.0'
                for at in nres.atoms:
                    at._charges['Kollman'] = 0.0
            return 
        CAatom = nres.atoms.get('CA')[0]
        hlist = Natom.findHydrogens()
        #5/5:assert len(hlist), 'polar hydrogens missing from n-terminus'
        if not len(hlist): 
            print 'polar hydrogens missing from n-terminus of chain ' + nres.parent.name
            #self.vf.warningMsg('polar hydrogens missing from n-terminus')
        if nres.type == 'PRO':
            #divide .059 additional charge between CA + CD
            #FIX THIS what if no CD?
            #CDatom = nres.atoms.get('CD')[0]
            CDatom = nres.atoms.get('CD')
            if CDatom is not None:
                CDatom = CDatom[0]
                CDatom._charges['Kollman'] = CDatom._charges['Kollman'] + .029
            else:
                print 'WARNING: no CD atom in ', nres.name
            Natom._charges['Kollman'] = Natom._charges['Kollman'] + .274
            CAatom._charges['Kollman'] = CAatom._charges['Kollman'] + .030
            for ha in hlist:
                ha._charges['Kollman'] = .333
        else:
            Natom._charges['Kollman'] = Natom._charges['Kollman'] + .257
            CAatom._charges['Kollman'] = CAatom._charges['Kollman'] + .055
            for ha in hlist:
                ha._charges['Kollman'] = .312
            
        return 


    def fixCterminus(self, cres):
        """newTotal<-fixCterminu(cres)
            cres is the C-terminal residue
        """
        OXYatomSet = cres.atoms.get('OXT')
        OCTatomSet = cres.atoms.get('OCT')
        if OXYatomSet is not None and len(OXYatomSet):
            OXYatom = OXYatomSet[0]
            OXYatom._charges['Kollman'] = -.706
            #CAUTION!
            CAatom = cres.atoms.get('CA')[0]
            CAatom._charges['Kollman'] = CAatom._charges['Kollman'] - .006
            #CAUTION!
            Catom = cres.atoms.get('C')[0]
            Catom._charges['Kollman'] = Catom._charges['Kollman'] - .082
            #CAUTION!
            Oatom = cres.atoms.get('O')[0]
            Oatom._charges['Kollman'] = Oatom._charges['Kollman'] - .206
            return 
        elif OCTatomSet is not None and len(OCTatomSet):
            OXYatom = OCTatomSet[0]
            OXYatom._charges['Kollman'] = -.706
            #CAUTION!
            CAatom = cres.atoms.get('CA')[0]
            CAatom._charges['Kollman'] = CAatom._charges['Kollman'] - .006
            #CAUTION!
            Catom = cres.atoms.get('C')[0]
            Catom._charges['Kollman'] = Catom._charges['Kollman'] - .082
            #CAUTION!
            Oatom = cres.atoms.get('O')[0]
            Oatom._charges['Kollman'] = Oatom._charges['Kollman'] - .206
            return 
            



    def fixHisRes(self, his):
        """newTotal<-fixHisRes(his)
            his is a HISTIDINE residue
        """
        hisAtomNames = his.atoms.name
        #oldcharge = Numeric.sum(his.atoms._charges['Kollman'])
        oldcharge = 0
        for at in his.atoms:
            oldcharge = oldcharge + at._charges['Kollman']
        assertStr = his.name + ' is lacking polar hydrogens'
        assert 'HD1' or 'HE2' in hisAtomNames, assertStr
        #get the appropriate dictionary
        if 'HD1' in hisAtomNames and 'HE2' in hisAtomNames:
            d = q['HIS+']
        elif 'HD1' in  hisAtomNames:
            d = q['HISD']
        elif 'HE2' in hisAtomNames:
            d = q['HIS']
        else:
            msgStr = his.full_name() + ' missing both hydrogens!'
            print msgStr
            return 

        #assign charges
        for a in his.atoms:
            if d.has_key(a.name):
                a._charges['Kollman'] = d[a.name]
            else:
                a._charges['Kollman'] = 0.0

        #correct total
        #newcharge = Numeric.sum(his.atoms._charges['Kollman'])
        newcharge = 0
        for at in his.atoms:
            newcharge = newcharge + at._charges['Kollman']


    def fixCysRes(self, cys):
        cysAtomNames = cys.atoms.name
        #oldcharge = Numeric.sum(cys.atoms._charges['Kollman'])
        oldcharge = 0
        for at in cys.atoms:
            oldcharge = oldcharge + at._charges['Kollman']
        #get the appropriate dictionary
        if 'HG' in cysAtomNames:
            d = self.q['CYSH']
        else:
            #cystine
            d = self.q['CYS']

        #assign charges
        for a in cys.atoms:
            if d.has_key(a.name):
                a._charges['Kollman'] = d[a.name]
            else:
                a._charges['Kollman'] = 0.0

        #correct total
        #newcharge = Numeric.sum(cys.atoms._charges['Kollman'])
        newcharge = 0
        for at in cys.atoms:
            newcharge = newcharge + at._charges['Kollman']



