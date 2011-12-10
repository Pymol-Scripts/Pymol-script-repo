## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#############################################################################
#
# Author: Ruth HUEY, William  LINDSTROM
#
# Copyright: M. Sanner TSRI 2003
#
#############################################################################


#
#
# $Id: LigandMixin.py,v 1.4 2008/11/21 20:02:01 vareille Exp $
#
#
#
#
#
#
#

"""
This MIXIN Object adds AutoDock Ligand features to a MolKit.Protein
   methods:
       setup
           parameters:
               flags:
                useProteinAromaticList 
                autoMergeNPHS
               maxtors = 32
               aromaticCutOff = 7.5
           creates fields:
               bondClassifier 
               torsionTree = None
               message = ""
               torscount = 0
               possible_torscount = 0
               TORSDOF = 0
               #should this be here?
               aromaticCarbons = AtomSet()
           processes molecule:
               detects whether isPeptide
               checks charges
                   add charges if necessary
               checks hydrogens
                   merge nphs
                   warn if no polar hydrogens
               checks bonds, classifying them
                    via bondClassifier.classify:
               checks Aromatic carbons
       setroot
            self.ROOT
            #fix this: required by write method
            self.ROOT.rnum0 = 0
       autoroot
            self.autoRoot
            plus self.ROOT and self.ROOT.rnum0 = 0
       turnOnTorsion
       turnOffTorsion
       setCarbonName
       setAromaticName
       buildTorsionTree**
       write
       
   bondClassifier:
       classify:
           rotatablebonds
           amidebonds
           peptidebackbonebonds
           aromaticbonds (??)
               ???aromatic carbon detector???
               detect aromatic carbons
           sets torscount
           sets possible_torscount
           sets TORSDOF
           
   torsion tree which has:
       ROOT
       methods:
           set/auto root
           turn_off/turn_on torsions
               ->merge nodes
               ->split nodes
       ??RESIDUES??
       


"""
import numpy.oldnumeric as Numeric, math, types, os
from string import find
from MolKit.pdbWriter import PdbqWriter
from AutoDockTools.AutoDockBondClassifier import AutoDockBondClassifier
from AutoDockTools.observer import Subject
from MolKit.molecule import Atom, AtomSet, BondSet
from MolKit.protein import ProteinSet
from MolKit.torTree import TorTree
from Pmv.qkollua import q
from PyBabel.atomTypes import AtomHybridization
from PyBabel.bo import BondOrder
from PyBabel.cycle import RingFinder
from PyBabel.gasteiger import Gasteiger
from mglutil.util.misc import isInstance

aromDict = {}
aromDict['TRP'] = ['NE1','CD1','CG','CD2','CE2','CZ2',\
                    'CH2','CZ3','CE3']
aromDict['TYR'] = ['CD1','CG','CD2','CE1','CE2','CZ']
aromDict['PHE'] = ['CD1','CG','CD2','CE1','CE2','CZ']
aromDict['HIS'] = ['CD2','CE1','CG','ND1','NE2']
aromDict['PRO'] = ['N','CD','CA','CB','CG']
pep_aromList = ['PHE_CD1', 'PHE_CG', 'PHE_CD2', 'PHE_CE1',\
    'PHE_CE2', 'PHE_CZ', 'TYR_CD1', 'TYR_CG', 'TYR_CD2', 'TYR_CE1',\
    'TYR_CE2', 'TYR_CZ', 'HIS_CD2', 'HIS_CE1', 'HIS_CG', 'TRP_CD1',\
    'TRP_CG', 'TRP_CD2', 'TRP_CE2', 'TRP_CZ2', 'TRP_CH2', 'TRP_CZ3',\
    'TRP_CE3']

#MixIn methods from Chuck Esterbrook, 'Issue 84: Using Mix-ins with Python',
#                   Linux Journal Posted Sunday, April 01, 2001

def MixIn(pyClass, mixInClass, makeAncestor=0):
    if makeAncestor:
        if mixInClass not in pyClass.__bases__:
            pyClass.__bases__ = (mixInClass,) + pyClass.__bases__
    else:
        # Recursively traverse the mix-in ancestor
        # classes in order to support inheritance
        baseClasses = list(mixInClass.__bases__)
        baseClasses.reverse()
        for baseClass in baseClasses:
            MixIn(pyClass, baseClass)
        # Install the mix-in methods into the class
        for name in dir(mixInClass):
            if not name.startswith('__'):
                #skip private members
                member = getattr(mixInClass, name)
                if type(member) is types.MethodType:
                    member = member.im_func
                    setattr(pyClass, name, member)

                    
def SimpleMixIn(pyClass, mixInClass):
    if mixInClass not in pyClass.__bases__:
        pyClass.__bases__ += (mixInClass,)


class LigandMixin(Subject):
    """ lfo adds capability to a protein/molecule to be used as a ligand 
            in an AutoDock docking

    """


    def setup(self, useProteinAromaticList=1,
                    aromaticCutOff=7.5, maxtors=32,
                    autoMergeNPHS=1):
        if self.chains[0].hasBonds==0:
            self.buildBondsByDistance()
        self.useProteinAromaticList = useProteinAromaticList
        self.aromaticCutOff = aromaticCutOff
        self.maxtors = 32
        self.autoMergeNPHS = autoMergeNPHS
        self.bondClassifier = AutoDockBondClassifier()
        self.torsionTree = None
        self.message = ""
        self.torscount = 0
        self.possible_torscount = 0
        self.TORSDOF = 0
        self.resKeys = q.keys()
        self.PdbqWriter = PdbqWriter()
        #should this be here?
        msg = 'initialized ' + self.name +':\n'
        #detect whether isPeptide
        msg += self.checkIsPeptide()
        #process charges
        msg += self.checkCharges()
        ##process hydrogens
        msg += self.checkHydrogens()
        #process bonds
        #this doesn't report anything???
        newmsg, nphs = self.checkBonds()
        msg += newmsg
        #process aromaticCarbons
        msg += self.checkAromatics()
        return msg, nphs


    def checkIsPeptide(self):
        """
        checkIsPeptide
        """
        #check whether each restype is in std list
        #if so, self is a peptide
        resSet = self.chains.residues
        dict = {}
        for r in resSet:
            dict[r.type] = 0
        for t in dict.keys():
            if t not in self.resKeys:
                self.isPeptide = 0
                return '    -it is not a peptide\n'
        #only get to this point if all
        #residue types were found
        self.isPeptide = 1
        return '    -it is a peptide\n'


    def checkCharges(self):
        """
        checkCharges
        """
        msg = '     -already had charges'
        if hasattr(self, 'checked_charges'):
            return msg
        msg = ""
        needToAdd = 0
        chargedAts = self.allAtoms.get(lambda x: hasattr(x, 'charge'))
        if not chargedAts:
            needToAdd = 1
        elif len(chargedAts)!=len(self.allAtoms):
            needToAdd = 1
        elif len(filter(lambda x:x.charge==0, chargedAts))==len(chargedAts):
            #this checks that each atom doesn't have charge=0
            needToAdd = 1
        #to add Kollman need to:
        #   look up each atom's parent in q to get dict=q[at.parent.type]
        #   find the atom's name in dict to get a charge
        #   set atom._charges['Kollman'] to that charge
        #   there are special problems with his, cys and termini
        #   set allAtoms.chargeSet to 'Kollman'
        if needToAdd:     
            if not hasattr(self, 'isPeptide'):
                msg = self.checkIsPeptide()
            if self.isPeptide:
                #self.checked_charges = 'needs Kollman'
                self.addKollman()
                self.checked_charges = 'has Kollman'
                #self.vf.addKollmanCharges(self, topCommand=0, log=0)
                msg = msg + '    -added Kollman charges\n'
            else:
                #self.checked_charges = 'needs gasteiger'
                self.computeGasteiger()
                self.checked_charges = 'has gasteiger'
                #self.vf.computeGasteiger(self, topCommand=0, log=0)
                msg = msg + '    -added gasteiger charges\n'
        else:
            self.checked_charges = 'has charges'
            msg = msg + '   -already had charges\n'
        # adt will have to add charges... AND check them
        #at this point, every atom has a charge and charges are not all 0k
        #have to have unit charges per residue
        #print 'calling checkMolCharges'
        #errCharge, resList = checkMolCharges(self, self.vf)
        #FIX THIS: need to add mechanism to adjust charges
        #gui to let user pick + change
        #or either change on first atom or spread over all
        #self.checked_charges = 1
        #if added charges, msg will say which one
        #if did't add charges, msg will be 'ERROR'
        return msg


    def computeGasteiger(self):
        #to compute Gasteiger need to:
        #   create an AtomHybridization()
        #   call its assignHybridization method on self.allAtoms
        #   create a Gasteiger()
        #   call its compute method on self.allAtoms
        # THEN move gast_charge into _charges with gasteiger key
        #   set allAtoms.chargeSet to 'gasteiger'
        # THEN delattr gast_charge from allAtoms
        allAts = self.allAtoms
        ah = AtomHybridization()
        ah.assignHybridization(allAts)
        Gast = Gasteiger()
        Gast.compute(allAts)
        gastCharges = []
        for c in allAts.gast_charge:
            gastCharges.append(round(c, 3))
        allAts.addCharges('gasteiger', gastCharges)
        del allAts.gast_charge
        allAts.chargeSet = 'gasteiger'


    def addKollman(self):
        #to add Kollman need to:
        #   look up each atom's parent in q to get dict=q[at.parent.type]
        #   find the atom's name in dict to get a charge
        #   set atom._charges['Kollman'] to that charge
        #   there are special problems with his, cys and termini
        #   set allAtoms.chargeSet to 'Kollman'

        for a in self.allAtoms:
            dict = q.get(a.parent.type, {})
            if len(dict):
                a._charges['Kollman'] = dict.get(a.parent.type, 0)
            else:
                a._charges['Kollman'] = 0
            a.chargeSet = 'Kollman'


    def checkHydrogens(self):
        """
        checkHydrogens
        """
        #what if self doesn't have bonds
        if not self.chains[0].hasBonds:
            self.buildBondsByDistance()
        hs = self.allAtoms.get(lambda x: x.element=='H')
        self.nphs = AtomSet()
        if not hs:
            msg = '    -no polar hydrogens found!\n'
        if hs:
            nphs = hs.get(lambda x: x.bonds[0].atom1.element=='C' or \
                    x.bonds[0].atom2.element=='C')
            if nphs:
                self.nphs = nphs
                msg = '    -found ' + str(len(self.nphs)) + ' nonpolar hydrogens\n'
            else:
                msg = '    -no nonpolar hydrogens\n'
        return msg


    def mergeNPHS(self):
        lenNPHS = len(self.nphs)
        if not lenNPHS:
            return
        nphs = self.nphs
        atLen = len(self.allAtoms)
        self.allAtoms = self.allAtoms - nphs
        #first add charge to nph's heavy atom
        chList = nphs[0]._charges.keys()
        if not len(chList):
            print 'no charges present'
            return 'XXX'
        #check that all nphs have a certain kind of charge
        for at in nphs:
            chs = at._charges.keys()
            for c in chList:
                if c not in chs:
                    chList.remove(c)
        if len(chList):
            for chargeSet in chList:
                for h in nphs:
                    b = h.bonds[0]
                    c = b.atom1
                    if c==h: 
                        c = b.atom2
                    c._charges[chargeSet] = c._charges[chargeSet] + h._charges[chargeSet]
        #next delete nphs
        for h in nphs:
            b = h.bonds[0]
            c = b.atom1
            if c==h: 
                c = b.atom2
            c.bonds.remove(b)
            h.bonds = BondSet()
            h.parent.remove(h)
            #nb: 
            #these don't show up in deleteMol.getFreeMemoryInformation
            del h
        self.nphs = AtomSet()
        return '        and merged them\n', nphs


    def checkBonds(self):
        """
        checkBonds
        """
        msg = ""
        #print self.name, ':'
        nphs = AtomSet()
        if len(self.nphs) and self.autoMergeNPHS:
            newmsg, nphs = self.mergeNPHS()
            msg = msg + newmsg
        bc = self.bondClassifier
        bonds = self.chains.residues.atoms.bonds[0]
        for b in bonds:
            b.activeTors = 0
            b.possibleTors = 0
            b.leaf = 0
        #FIX THIS:
        #if isPeptide, don't get cycles this way
        #if self.isPeptide:
            #cycleSelector = bc['cycles']
            #del bc['cycles']
        results = bc.classify(bonds)
        #if self.isPeptide:
            #bc['cycles'] = cycleSelector
            #results['cycles'] = self.getPeptideBondDict()
        #for k,v in results.items():
        #    print k,' = ', len(v)
        self.rotatable = results['rotatable']
        self.torscount = len(self.rotatable) 
        self.possible_torscount = self.torscount
        for b in self.rotatable:
            b.activeTors = 1
            b.possibleTors = 1
        for b in results['leaf']:
            b.leaf = 1
        for b in results['cycle']:
            b.incycle = 1
        hydrogenRotators = results['hydrogenRotators']
        self.hydrogenRotators = hydrogenRotators 
        self.TORSDOF = self.torscount - len(hydrogenRotators)
        ptAts = bc.dict['rotatable'].getAtoms(self.rotatable)
        d = {}
        for a in ptAts:
            d[a] = 0
        self.pTatomset = AtomSet(d.keys())
        #print 'len(pTatomset=', len(self.pTatomset)
        self.leafbonds = results['leaf']
        self.pepbackbonds = results['ppbb']
        self.amidebonds = results['amide']
        self.cyclebonds = results['cycle']
        return msg, nphs


    def checkAromatics(self):
        """
        checkAromatics
        """
        #this depends on userPref useProteinAromaticList
        if not len(self.cyclebonds):
            self.aromaticCs = AtomSet()
            return ""
        if self.isPeptide and self.useProteinAromaticList:
            self.aromaticCs = self.getPeptideAromatics()
            return ""
        if self.isPeptide:
            self.getPeptideBondDict()
        else:
            self.getLigandBondDict()
        counter = 0
        while counter < self.cyclecount:
            counter = counter + 1
            blist = self.bondDict[counter]
            for item in blist:
                at = item.atom1
                self._getAdjAtom(item, blist)
                #now each bond has 3 atoms specified for it: its own two and the next1 to atom1
                result = self._getNormal(item)
                item.nrmsize = result[0]
                item.nrms = result[1]
                #next find the other bond w/atom2:
                z2 = filter(lambda x,item=item, at2=item.atom2, blist=blist:x!=item and x.atom1==at2 or x.atom2==at2, blist)
                #finally, return the other atom in this bond
                item.nextbond2 = z2[0]
                if item.nextbond2==item:
                    item.nextbond2 = z2[1]
                neighbor2 = self._getnxtAtom(item.atom2,item.nextbond2)
                item.next2 = neighbor2
                #next have to check whether the normals are parallel
                #check each pair in each bond, how to keep track??
                #have to get normal at item's atom2: so have to get next2 for this bond:
            #have to loop twice to make sure neighbor has nrms
            for item in blist:
                p = item.nrms
                psize = item.nrmsize
                q = item.nextbond2.nrms
                qsize = item.nextbond2.nrmsize
                #theta is the critical test for planarity:
                #if angle between 2 nrms is 0, atoms are planar
                #NB>test is comparing theta,cos(angle), w/zero
                item.theta = Numeric.dot(p,q)/(psize*qsize)
                for p in ['next1','next2','nextbond','nextbond2']:
                    delattr(item, p)
        self.updateAromatics(self.aromaticCutOff)
        msg = '    -found '+ str(len(self.aromaticCs)) + ' aromatic carbons\n'
        return  msg


    def updateAromatics(self, cutoff):
        self.aromaticCutOff = cutoff
        #cutoff = self.aromaticCutOff
        cutoffValue = math.cos(cutoff*math.pi/180.)
        aromaticCs = AtomSet([])
        atD = {}
        #to keep from overwriting names of aromatic carbons at 
        #junctions of rings use this klug
        for blist in self.bondDict.values():
            for item in blist:
                item.atom1.setThisTime = 0
                item.atom2.setThisTime = 0
        for blist in self.bondDict.values():
            ct = 0
            for item in blist:
                #these are values for the default 7.5degrees:
                #if theta>=0.997 or theta<=-0.997:
                if item.theta>=cutoffValue or item.theta<=-cutoffValue:
                    item.posAromatic = 1
                    ct = ct + 1
                else:
                    item.posAromatic = 0
            #after checking all the bonds in current cycle, compare #posAromatic w/number
            if ct==len(blist):
                #and change the name and autodock_element here....
                for b in blist:
                    at1 = b.atom1
                    at2 = b.atom2
                    if at1.element=='C':
                        at1.name = 'A' + at1.name[1:]
                        at1.autodock_element = 'A'
                        atD[at1] = 0
                        at1.setThisTime = 1
                    if at2.element=='C':
                        at2.name = 'A' + at2.name[1:]
                        at2.autodock_element = 'A'
                        atD[at2] = 0
                        at2.setThisTime = 1
                aromaticCs = AtomSet(atD.keys())
            else:
                #if there were any aromatic carbons which no longer 
                #meet the criteria, change them back
                for b in blist:
                    at1 = b.atom1
                    at2 = b.atom2
                    if at1.name[0]=='A' and not at1.setThisTime:
                        at2.autodock_element = 'C'
                        at1.name = 'C' + at1.name[1:]
                    if at2.name[0]=='A'and not at2.setThisTime:
                        at2.autodock_element = 'C'
                        at2.name = 'C' + at2.name[1:]
        #remove klug
        for blist in self.bondDict.values():
            for item in blist:
                if hasattr(item.atom1, 'setThisTime'):
                    delattr(item.atom1,'setThisTime')
                if hasattr(item.atom2, 'setThisTime'):
                    delattr(item.atom2,'setThisTime')
        for a in aromaticCs:
            a.autodock_element = 'A'
        self.aromaticCs = aromaticCs

            
    def restoreCarbons(self):
        for a in self.aromaticCs:
            if len(a.name)==1:
                a.name = 'C'
            else:
                a.name = 'C' + a.name[1:]
        a.autodock_element = 'C'

   
    def nameAromatics(self):
        for a in self.aromaticCs:
            if len(a.name)==1:
                a.name = 'A'
            else:
                a.name = 'A' + a.name[1:]
        a.autodock_element = 'A'


    def addAromatic(self, at):
        if at.element!='C':
            print at.name, ' is not a carbon'
            return 'ERROR'
        if at in self.aromaticCs:
            print at.name, ' is already in aromatic set'
            return 'ERROR'
        at.autodock_element = 'A'
        self.aromaticCs.append(at)
        print at.name, ' added to aromatic set'
   

    def removeAromatic(self, at):
        if at not in self.aromaticCs:
            print at.name, ' is not in aromatic set'
            return 'ERROR'
        self.aromaticCs.remove(at)
        at.autodock_element = 'C'
        print at.name, ' removed from aromatic set'
   

    def getPeptideAromatics(self):
        atSet = AtomSet()
        allAts = self.allAtoms
        arom_ats = allAts.get(lambda x, l = pep_aromList:\
            x.parent.type+'_'+x.name in l)
        if not arom_ats:
            return AtomSet([])
        for at in arom_ats:
            at.name = 'A' + at.name[1:]
            at.autodock_element = 'A'
        #print 'returning ', len(arom_ats), ' arom_ats'
        return arom_ats

       
    def getPeptideBondDict(self):
        resSet = self.chains.residues.get(lambda x, \
                    d = aromDict.keys():x.type in d)
        if not resSet:
            return BondSet()
        self.cyclecount = numres = len(resSet)
        #NB: cyclenum is 1 based because 0 means not numbered
        for i in range(numres):
            res = resSet[i]
            ats = res.atoms
            #bondDict keys are 1 based too
            keys = aromDict[res.type]
            bnds = bondDict[i+1] = ats.get(lambda x, \
                                keys=keys:x.name in keys).bonds[0]
            for b in bnds:
                bnds.cyclenum = i + 1
        self.bondDict = bondDict


    def getLigandBondDict(self):
        ats = self.allAtoms
        babel = AtomHybridization()
        babel.assignHybridization(self.allAtoms)
        #typeBonds???
        #typeAtoms(ats)
        rf = RingFinder()
        rf.findRings2(ats, ats.bonds[0])
        ct = 1
        bondDict = {}
        for ring in rf.rings:
            bondDict[ct] = ring['bonds']
            for b in ring['bonds']:
                b.cyclenum = ct
            ct = ct + 1
        self.cyclecount = rf.ringCount
        self.bondDict = bondDict

 
    def _numberCycle(self, at):
        for b in at.bonds:
            if hasattr(b,'incycle') and b.cyclenum==0 :
                b.cyclenum = self.cyclecount
                nxtAtom = b.atom1
                if nxtAtom==at:
                    nxtAtom = b.atom2
                self._numberCycle(nxtAtom)


    def _getnxtAtom(self, at,b):
        nxtAtom = b.atom1
        if nxtAtom==at:
            nxtAtom = b.atom2
        return nxtAtom


    def _getAdjAtom(self,b,z):
        #first get the bonds in this cycle
        at = b.atom1
        #next find the other one with at as one of the atoms:
        z2 = filter(lambda x,b=b,at=at,z=z:x!=b and x.atom1==at or x.atom2==at, z)
        #finally, return the other atom in this bond
        b.nextbond = z2[0]
        neighbor = self._getnxtAtom(at,z2[0])
        b.next1 = neighbor


    def _getNormal(self,b):
        at1 = b.next1
        at2 = b.atom1
        at3 = b.atom2
        pt1 = at1.coords
        pt2 = at2.coords
        pt3 = at3.coords
        a = Numeric.subtract(pt2,pt1)
        b = Numeric.subtract(pt3,pt2)
        p = (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
        nrmsize = Numeric.sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])
        return (nrmsize, p)


    def setroot(self, atom):
        """
        setroot to 'C11' or 'hsg1:A:ILE2:CA'
        """
        if type(atom)==types.StringType:
            if find(atom, ':')>-1:
                #have to use full_name list
                #check that it is an atom
                mols = ProteinSet([self])
                nodes = mols.NodesFromName(atom)
                if not nodes:
                    return 'invalid root name'
                if nodes.__class__!=AtomSet:
                    return 'invalid root name: not atom level'
            else:
                nodes = self.allAtoms.get(lambda x, n = atom:x.name==n)
            if not nodes:
                return 'invalid root name'
            if len(nodes)>1:
                return 'root name must be unique'
            atom = nodes[0]
        #elif type(atom)!=types.InstanceType:
        elif isInstance(atom) is False:
            return atom, ' invalid root atom'
        elif atom.__class__!=Atom:
            return atom, ' can only select an Atom instance as root'
        #fix this rnum0 stuff
        #in case toggling back and forth
        if hasattr(self, 'autoRoot') and hasattr(self.autoRoot, 'rnum0'):
            delattr(self.autoRoot, 'rnum0')
        #if there is an old root, remove rnum0
        if hasattr(self, 'ROOT') and hasattr(self.ROOT, 'rnum0'):
            delattr(self.ROOT, 'rnum0')
        self.ROOT = atom
        self.ROOT.rnum0 = 0
        return self.ROOT


    def autoroot(self):
        """
        autoRoot
        """
        #clear old root
        if hasattr(self, 'ROOT') and hasattr(self.ROOT, 'rnum0'):
            delattr(self.ROOT, 'rnum0')
        if hasattr(self, 'autoRoot'):
            self.ROOT = self.autoRoot
            self.ROOT.rnum0 = 0
            return
        if len(self.chains)>1:
            return  "AutoRoot not implemented for molecules with >1 chain"
        self.bestBranch = len(self.allAtoms)
        bestList = []
        for item in self.allAtoms:
            if not hasattr(item, 'bonds'): continue
            if len(item.bonds)==1 and item.bonds[0].leaf: continue
            if hasattr(item,'leaf'): continue
            item.maxbranch = 0
            for b in item.bonds:
                nxtAtom = self._getnxtAtom(item,b)
                if not b.leaf:
                    thistree = self.subTree(item,nxtAtom, self.allAtoms)
                    thisbranch = len(thistree)
                    if thisbranch>item.maxbranch:
                        item.maxbranch = thisbranch
            #bestList holds all current best choices for Root..
            if item.maxbranch<self.bestBranch:
                bestList = []
                bestList.append(item)
                self.bestBranch = item.maxbranch
            if item.maxbranch==self.bestBranch and item not in bestList:
                bestList.append(item)
        if len(bestList)>1:
            foundCycle = 0
            for at in bestList:
                at.cycleatom = 0
                for b in at.bonds:
                    #4/29: peptides won't have aromatic 
                    #but will have incycle
                    #if b.bondOrder=='aromatic':
                    if hasattr(b,'incycle'):
                        at.cycleatom = 1
                        continue
            for at in bestList:
                if at.cycleatom:
                    self.ROOT = at
                    self.autoRoot = at
                    self.ROOT.rnum0 =0
                    foundCycle = 1
                    break
            #if bestList had a cycle atom, it's been set to root..if NOT:
            if not foundCycle:
                self.autoRoot = bestList[0]
                self.ROOT = bestList[0]
                self.ROOT.rnum0 =0
        #no ties for possible root, use bestRoot...
        else:
            self.autoRoot = bestList[0]
            self.ROOT = bestList[0]
            self.ROOT.rnum0 =0
        return self.ROOT


    def buildTorsionTree(self):
        if not hasattr(self, 'ROOT'):
            print 'must set ROOT first!'
            return
        self.torTree = TorTree(self.parser, self.ROOT)


    def turnOnTorsion(self, bnd):
        """
        turnOnTorsion 
        """
        if not bnd in self.rotatable:
            return
        #this is redundant (??)
        if not bnd.possibleTors:
            return
        if not bnd.activeTors:
            self.torscount = self.torscount + 1
            bnd.activeTors = 1
            self.torTree.addTorsion(bnd.atom1.number, bnd.atom2.number)


    def turnOffTorsion(self, bnd):
        """
        turnOffTorsion 
        """
        if not bnd in self.rotatable:
            return
        #this is redundant (??)
        if not bnd.possibleTors:
            return
        if bnd.activeTors:
            self.torscount = self.torscount - 1
            bnd.activeTors = 0
            self.torTree.removeTorsion(bnd.atom1.number, bnd.atom2.number)


    def limitTorsions(self, numTors, type, simpleModel=1):
        """sets number of torsions to specified number by inactivating either those which
move the fewest atoms or those which move the most. if number is > than
current but less than possible, torsions are reactivated
        """

        print 'lT: numTors=', numTors, ' type=', type
        allAts = self.allAtoms
        root = self.ROOT
        torscount = self.torscount
        print 'torscount=', torscount, 
        #NB: torscount is not necessarily the max
        #ie it could have been adjusted already
        at0 = self.allAtoms[0]
        if not hasattr(self, 'torTree') or not hasattr(at0, 'tt_ind'):
            self.torTree = TorTree(self.parser, root)
        torsionMap = self.torTree.torsionMap
        torsionMapNum = len(torsionMap)
        print 'len(tm)=', torsionMapNum
        possibleTors = self.possible_torscount
        if simpleModel:
            self.setTorsions(numTors, type)
            return 

        #FIX THIS: what if want to increase active torsions
        if torscount==numTors:
            msg = 'specified number==number present: no adjustment'
            return 'msg'
        elif torscount<numTors:
            if torscount==possibleTors:
                #if torscount==torsionMapNum:
                msg = 'specified number == number possible: no adjustment'
                return msg
            else:
                #in this case turn on as many as possible
                #if numTors>=torsionMapNum:
                if numTors>=possibleTors:
                    #turn on everything
                    delta = possibleTors - torscount
                    #delta = torsionMapNum - torscount
                else:
                    delta = numTors - torscount
                self.turnOnTorsions(delta, type)
        else:
            #torscount>numTors
            #in this case turn them off 
            delta = torscount - numTors
            self.turnOffTorsions(delta, type)
        msg = str(delta) + ' torsions changed'
        return msg


    def setTorsions(self, numTors, type):
        #this assumes atoms have tt_ind field
        tNum = len(self.torTree.base_torsionMap)
        msg = ""
        if numTors>tNum:
            msg = 'too many torsions specified! '+ str(numTors)+  ' reducing to'+str(tNum)
            numTors = tNum
        if type=='fewest':
            rangeList = range(numTors)
        else:
            rangeList = []
            for k in range(1, numTors+1):
                rangeList.append(-k)
        #turn them all off
        baseTorsionMap = self.torTree.base_torsionMap
        torsionMap = self.torTree.torsionMap
        #turn all active nodes off
        active_nodes = torsionMap[:]
        for node in active_nodes:
            #print 'turning off node ', node.number
            self.torTree.removeTorsion(node.bond[0], node.bond[1])
            b = self.allAtoms.get(lambda x, node=node: x.tt_ind in node.bond).bonds[0][0]
            b.activeTors = 0

        #turn on the right number at correct end
        for i in rangeList:
            node = baseTorsionMap[i]
            #print '2:turning on node ', node.number
            self.torTree.addTorsion(node.bond[0], node.bond[1])
            b = self.allAtoms.get(lambda x, node=node: x.tt_ind in node.bond).bonds[0][0]
            b.activeTors = 1
            #print 'now len(torsionMap)=', len(torsionMap)
        self.torscount = numTors
        return msg


    def turnOnTorsions(self, delta, type = 'fewest'):
        allAts = self.allAtoms
        torsionMap = self.torTree.torsionMap
        baseTorsionMap = self.torTree.base_torsionMap
        torscount = self.torscount
        #turn on delta torsions + adjust torscount in dict 
        if type=='fewest':
            rangeList = range(delta)
        else:
            rangeList = []
            for k in range(1, delta+1):
                rangeList.append(-k)
        #FIX THIS: probably doesn't work
        for i in rangeList:
            node = baseTorsionMap[i]
            b = allAts.get(lambda x, node=node: x.tt_ind in node.bond).bonds[0][0]
            if not b.activeTors:
                self.torTree.addTorsion(node.bond[0], node.bond[1])
                b.activeTors = 1
            else:
                lastInd = rangeList[-1]
                if type=='fewest':
                    rangeList.append(lastInd+1)
                else:
                    rangeList.append(lastInd-1)
        #torscount should be torscount + delta here
        self.torscount = torscount + delta

        
    def turnOffTorsions(self, delta, type = 'fewest'):
        allAts = self.allAtoms
        torsionMap = self.torTree.torsionMap
        baseTorsionMap = self.torTree.base_torsionMap
        torscount = self.torscount
        #turn on delta torsions + adjust torscount in dict 
        if type=='fewest':
            rangeList = range(delta)
        else:
            rangeList = []
            for k in range(1, delta+1):
                rangeList.append(-k)
        for i in rangeList:
            node = baseTorsionMap[i]
            if node.bond==(None,None):
                print 'error in turnOff torsions with ', rangeList
                break
            b = allAts.get(lambda x, node=node: x.tt_ind in node.bond).bonds[0][0]
            if b.activeTors:
                self.torTree.removeTorsion(node.bond[0], node.bond[1])
                b.activeTors = 0
            else:
                lastInd = rangeList[-1]
                if type=='fewest':
                    rangeList.append(lastInd+1)
                else:
                    rangeList.append(lastInd-1)
        self.torscount = torscount - delta


    def setCarbonName(self, at):
        """
        setCarbonName 
        """
        if not at.element=='C':
            return
        if len(at.name)==1:
            at.name = 'C'
        else:
            at.name = 'C' + at.name[1:]



    def setAromaticName(self, at):
        """
        setAromaticName 
        """
        if not at.element=='C':
            return
        if len(at.name)==1:
            at.name = 'A'
        else:
            at.name = 'A' + at.name[1:]


