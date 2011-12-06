## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

############################################################################
#
# Author:  Ruth Huey
#
# Copyright: M. Sanner TSRI 2005
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/MolKit/hydrogenBondBuilder.py,v 1.8.2.1 2011/04/13 18:07:36 rhuey Exp $
#
# $Id: hydrogenBondBuilder.py,v 1.8.2.1 2011/04/13 18:07:36 rhuey Exp $
#

"""
This module implements the HydrogenBondBuilder class which builds hydrogen
bonds between appropriate atoms.

"""
def getHBAtoms(mol):
    """
    Function to identify donor and acceptor atoms for hydrogen bonds

    d2, d3, a2, a3 <- getHBAtoms(mol)

    The function will first build bonds then assign atom hybridization
    and bond orders. Next it will use a HydrogenBondBuilder to identify and'
    report donor and acceptor atoms in sets corresponding to Sp2 and Sp3
    hybridization.
    A .hbstatus is added to all atoms. The value can be:
      NE: NEutral
      D2: Sp2 donor 
      D3: Sp3 donor 
      A2: Sp2 acceptor
      A3: Sp3 acceptor
      B2: Sp2 acceptor and donor
      B3: Sp3 acceptor and donor
    """
    # add bonds
    mol.buildBondsByDistance()
    atoms = mol.allAtoms

    # assign .babel_type
    from PyBabel.atomTypes import AtomHybridization
    atyper = AtomHybridization()
    atyper.assignHybridization(atoms)

    # assign bond order
    from PyBabel.bo import BondOrder
    bond_orderer = BondOrder()
    bond_orderer.assignBondOrder(atoms, atoms.bonds[0])

    # get donors
    from MolKit.hydrogenBondBuilder import HydrogenBondBuilder
    hbbuilder = HydrogenBondBuilder()
    donorTypes = hbbuilder.paramDict['donorTypes']
    donorsp2Ats, donorsp3Ats = hbbuilder.getHBDonors(atoms, donorTypes)

    # get acceptors
    acceptorTypes = hbbuilder.paramDict['acceptorTypes']
    acceptorsp2Ats, acceptorsp3Ats = hbbuilder.getHBAcceptors(atoms, acceptorTypes)
    # create hbstatus attribute
    atoms.hbstatus = 'NE'
    for a in donorsp2Ats:
        a.hbstatus = 'D2'
    for a in donorsp3Ats:
        a.hbstatus = 'D3'
    for a in acceptorsp2Ats:
        a.hbstatus = 'A2'
    for a in acceptorsp3Ats:
        a.hbstatus = 'A3'

    # handle atoms that are both donors and acceptor
    for a in donorsp2Ats.inter(acceptorsp2Ats):
        a.hbstatus = 'B2'
    for a in donorsp3Ats.inter(acceptorsp3Ats):
        a.hbstatus = 'B3'

    return donorsp2Ats, donorsp3Ats, acceptorsp2Ats, acceptorsp3Ats


#need these in order to typeAtoms
from PyBabel.atomTypes import AtomHybridization
from PyBabel.bo import BondOrder

from MolKit.molecule import Atom, AtomSet, HydrogenBond
from MolKit.distanceSelector import DistanceSelector

#needed for some math..
import numpy.oldnumeric as Numeric, math

#set up lists of types for donors/acceptors
sp2Donors = ['Nam', 'Ng+', 'Npl']
sp3Donors = ['N3+', 'O3','S3']
allDonors = ['Nam', 'Ng+', 'Npl', 'N3+', 'O3','S3']
#added Nam because of bdna cytidine N3
#NB: Npl cannot be an acceptor if the sum of its bonds' bondOrder>2
sp2Acceptors = ['O2', 'O-', 'Npl', 'Nam']
sp3Acceptors = ['S3', 'O3']
allAcceptors = ['O2', 'O-', 'Npl', 'Nam', 'S3', 'O3']
#????
#for n in allAcceptors:
#    allAcceptors.append('a'+n)


def dist(c1, c2):
    d = Numeric.array(c2) - Numeric.array(c1) 
    ans = math.sqrt(Numeric.sum(d*d))
    return round(ans, 3)


def getAngle(ac, hat, don ):
    #acCoords = getTransformedCoords(ac)
    #hCoords = getTransformedCoords(hat)
    #dCoords = getTransformedCoords(don)
    acCoords = ac.coords
    hCoords = hat.coords
    dCoords = don.coords
    pt1 = Numeric.array(acCoords, 'f')
    pt2 = Numeric.array(hCoords, 'f')
    pt3 = Numeric.array(dCoords, 'f')
    #pt1 = Numeric.array(ac.coords, 'f')
    #pt2 = Numeric.array(hat.coords, 'f')
    #pt3 = Numeric.array(don.coords, 'f')
    v1 = Numeric.array(pt1 - pt2)
    v2 = Numeric.array(pt3 - pt2)
    dist1 = math.sqrt(Numeric.sum(v1*v1))
    dist2 = math.sqrt(Numeric.sum(v2*v2))
    sca = Numeric.dot(v1, v2)/(dist1*dist2)
    if sca>1.0:
        sca = 1.0
    elif sca<-1.0:
        sca = -1.0
    ang =  math.acos(sca)*180./math.pi
    return round(ang, 5)


def applyTransformation(pt, mat):
    pth = [pt[0], pt[1], pt[2], 1.0]
    return Numeric.dot(mat, pth)[:3]


def getTransformedCoords(atom):
    # when there is no viewer, the geomContainer is None
    if atom.top.geomContainer is None:
        return atom.coords
    g = atom.top.geomContainer.geoms['master']
    c = applyTransformation(atom.coords, g.GetMatrix(g))
    return  c.astype('f')


distCutoff=2.25     #H-Acc Distance cutoff to mean 2.06 std dev 0.19
                    #@@ distCutoff2 increased 4/2011 
                    #   based on data: www.biochem.ucl.ac.uk/bsm/atlas/mc.html
distCutoff2=3.15    #Donor-Acc Distance cutoff @@mean 2.98 std dev 0.17@@
d2min=120           #sp2-donor-hydrogen-acceptor angle 
d2max=180           #
d3min=120           #sp3-donor-hydrogen-acceptor angle
d3max=180           #
a2min=110           #sp2-donor-acceptor-acceptorN angle, 
a2max=150           #    where acceptorN is 'neighbor' atom bonded to acceptor
a3min=100           #sp3-donor-acceptor-acceptorN angle
a3max=150           #    where acceptorN is 'neighbor' atom bonded to acceptor


class HydrogenBondBuilder:
    """
    object which can build hydrogen bonds between atoms according
    to their coords and atom type
    """

    def __init__(self, distCutoff=distCutoff, distCutoff2=distCutoff2,
                        d2min=d2min, d2max=d2max, d3min=d3min, d3max=d3max,
                        a2min=a2min, a2max=a2max, a3min=a3min, a3max=a3max,
                        donorTypes=allDonors, acceptorTypes=allAcceptors,
                        distOnly=False):
        d = self.paramDict = {}
        d['distCutoff'] = distCutoff
        d['distCutoff2'] = distCutoff2
        d['d2min'] = d2min
        d['d2max'] = d2max
        d['d3min'] = d3min
        d['d3max'] = d3max
        d['a2min'] = a2min
        d['a2max'] = a2max
        d['a3min'] = a3min
        d['a3max'] = a3max
        d['donorTypes'] = donorTypes
        d['acceptorTypes'] = acceptorTypes
        d['distOnly'] = distOnly
        self.distSelector = DistanceSelector(return_dist=0)


    def check_babel_types(self, ats):
        num_ats = len(ats)
        num_babel_type = len(ats.get(lambda x: hasattr(x, 'babel_type')))
        num_bnd_type = len(ats.get(lambda x: hasattr(x, 'bnd_type')))
        if (num_babel_type!=num_ats) or (num_bnd_type!=num_ats):
            babel = AtomHybridization()
            bond_orderer = BondOrder()
            tops = ats.top.uniq()
            for mol in tops: 
                babel.assignHybridization(mol.allAtoms)
                bond_orderer.assignBondOrder(mol.allAtoms, mol.allAtoms.bonds[0])
                mol.allAtoms._bndtyped = 1


    def reset(self, ats):
        tops = ats.top.uniq()
        for mol in tops:
            for a in mol.allAtoms:
                if hasattr(a, 'hbonds'):
                    for item in a.hbonds:
                        del item
                    delattr(a, 'hbonds')


    def build(self, group1, group2=None, reset=True, paramDict=None):
        """atDict <- build(group1, group2, reset, paramDict=None, **kw):
            group1: atoms 
            group2: atoms 
            reset: remove all previous hbonds, default True!

            paramDict: a dictionary with these keys and default values
            distCutoff: 2.25  hydrogen--acceptor distance
            distCutoff2: 3.00 donor... acceptor distance
            d2min: 120 <min theta for sp2 hybridized donors>
            d2max: 180 <max theta for sp2 hybridized donors>
            d3min: 120 <min theta for sp3 hybridized donors>
            d3max: 170 <max theta for sp3 hybridized donors>

            a2min: 120 <min phi for sp2 hybridized donors>
            a2max: 150 <max phi for sp2 hybridized donors>
            a3min: 100 <min phi for sp3 hybridized donors>
            a3max: 150 <max phi for sp3 hybridized donors>
            @@FIX THIS: these do not seem to be input here
            donorTypes = allDonors
            acceptorTypes = allAcceptors
           """
        #setup parameter dictionary
        if paramDict is None:
            paramDict = self.paramDict

        #setup group2 
        if group2 is None:
            group2 = group1

        #process each group
        #group1
        if group1.__class__!=Atom:
            group1 = group1.findType(Atom)
            #print "now group1=", group1.full_name()
            if not len(group1):
                return "ERROR"  #@@ OR WHAT?
        self.check_babel_types(group1)

        #group2
        if group2.__class__!=Atom:
            group2 = group2.findType(Atom)
            #print "now group2=", group2.full_name()
            if not len(group2):
                return "ERROR"  #@@ OR WHAT?
        self.check_babel_types(group2)

        if reset:
            #do optional reset: remove all prior hbonds 
            self.reset(group1)
            self.reset(group2)

        #print "group1=", len(group1)
        #print "group2=", len(group2)

        #buildHbonds
        atDict = {}
        dict1 = self.buildD(group1, paramDict)

        #@@what is this doing here???
        #sp2 hybridized atoms
        #dAts2 = ats.get(lambda x, l=sp2: x.babel_type in l)
        if group1==group2:
            #only one step:
            atD1 = self.process(dict1, dict1, paramDict)
            #print 'len(atD1).keys()=', len(atD1.keys())
            atD2 = {}
        else:
            # two steps:
            # 1: group1 donors v group2 acceptors 
            # 2: group1 acceptors vs group2 donors
            dict2 = self.buildD(group2, paramDict)
            atD1 = self.process(dict1, dict2, paramDict)
            atD2 = self.process(dict2, dict1, paramDict)
            #print 'len(atD1).keys()=', len(atD1.keys())
            #print 'len(atD2).keys()=', len(atD2.keys())
        #if called with 1 atom could get tuple of two empty dictionaries
        if type(atD1)==type(atDict):
            if len(atD1):
                atDict.update(atD1)
        if type(atD2)==type(atDict):
            if len(atD2):
                atDict.update(atD2)
        #@@DESCRIBE atDict
        #this dict has atomSets as keys and as values(?check that)
        return atDict


    def checkForPossibleH(self, ats, blen):
        #@@FIX THIS: WHAT IS THE POINT OF THIS???
        #check that if at has all bonds, at least one is to a hydrogen
        # have to do this by element??
        probAts = AtomSet(ats.get(lambda x, blen=blen: len(x.bonds)==blen))
        #probOAts = ats.get(lambda x, blen=blen: len(x.bonds)==blen)
        #probSAts = ats.get(lambda x, blen=blen: len(x.bonds)==blen)
        if probAts:
            rAts = AtomSet([])
            for at in probAts:
                if not len(at.findHydrogens()):
                    rAts.append(at)
            if len(rAts):
                ats =  ats.subtract(rAts)
        return ats
    

    def getHBDonors(self, ats, donorList):
        #getHBDonors 
        sp2 = []
        sp3 = []
        for item in sp2Donors:
            if item in donorList: sp2.append(item)
        for item in sp3Donors:
            if item in donorList: sp3.append(item)
        
        dAts2 = ats.get(lambda x, l=sp2: x.babel_type in l)
        if not dAts2: dAts2=AtomSet([])
        else: 
            dAts2 = self.checkForPossibleH(dAts2, 3)
        
        dAts3 = ats.get(lambda x, l=sp3: x.babel_type in l)
        if not dAts3: dAts3=AtomSet([])
        else: dAts3 = self.checkForPossibleH(dAts3, 4)
        return dAts2, dAts3


    def filterAcceptors(self, accAts):
        ntypes = ['Npl', 'Nam']
        npls = accAts.get(lambda x, ntypes=ntypes: x.babel_type=='Npl')
        nams = accAts.get(lambda x, ntypes=ntypes: x.babel_type=='Nam')
        #nAts = accAts.get(lambda x, ntypes=ntypes: x.babel_type in ntypes)
        restAts = accAts.get(lambda x, ntypes=ntypes: x.babel_type not in ntypes)
        if not restAts: restAts = AtomSet([])
        #if nAts:
        if npls:
            #for at in nAts:
            for at in npls:
                s = 0
                for b in at.bonds:
                    if b.bondOrder=='aromatic':
                        s = s + 2
                    else: s = s + b.bondOrder
                #if s<3:
                #apparently this is wrong
                if s<4:
                    restAts.append(at)
        if nams:
            #for at in nAts:
            for at in nams:
                s = 0
                for b in at.bonds:
                    if b.bondOrder=='aromatic':
                        s = s + 2
                    else: s = s + b.bondOrder
                    #s = s + b.bondOrder
                if s<3:
                    restAts.append(at)
        return restAts
                

    def getHBAcceptors(self, ats, acceptorList):
        #print "getHBAcceptors: acceptorList=", acceptorList
        #getHBAcceptors 
        sp2 = []
        sp3 = []
        for item in sp2Acceptors:
            if item in acceptorList: sp2.append(item)
        for item in sp3Acceptors:
            if item in acceptorList: sp3.append(item)

        dAts2 = AtomSet(ats.get(lambda x, l=sp2: x.babel_type in l))
        if dAts2: 
            dAts2 = self.filterAcceptors(dAts2)
        dAts3 = AtomSet(ats.get(lambda x, l=sp3: x.babel_type in l))
        return dAts2, dAts3
      

    def buildD(self, ats, paramDict=None):
        if paramDict is None:
            paramDict = self.paramDict
        #these are from the __call__ method of vf.buildHBonds
        if not paramDict.has_key('distCutoff'):
            paramDict['distCutoff'] = 2.25
        if not paramDict.has_key('distCutoff2'):
            paramDict['distCutoff2'] = 3.00
        if not paramDict.has_key('d2min'):
            paramDict['d2min'] = 120.
        if not paramDict.has_key('d2max'):
            paramDict['d2max'] = 180.
        if not paramDict.has_key('d3min'):
            paramDict['d3min'] = 120.
        if not paramDict.has_key('d3max'):
            paramDict['d3max'] = 170.
        if not paramDict.has_key('a2min'):
            paramDict['a2min'] = 130.
        if not paramDict.has_key('a2max'):
            paramDict['a2max'] = 170.
        if not paramDict.has_key('a3min'):
            paramDict['a3min'] = 120.
        if not paramDict.has_key('a3max'):
            paramDict['a3max'] = 170.
        if not paramDict.has_key('distOnly'):
            paramDict['distOnly'] = 0
        if not paramDict.has_key('donorTypes'):
            paramDict['donorTypes'] = allDonors
        if not paramDict.has_key('acceptorTypes'):
            paramDict['acceptorTypes'] = allAcceptors

        d = {}
        donorTypes = paramDict['donorTypes']
        donor2Ats, donor3Ats = self.getHBDonors(ats, donorTypes)
        d23 = donor2Ats + donor3Ats
        #hAts = ats.get(lambda x, d23=d23: x.element=='H' \
                    #and x.bonds[0].neighborAtom(x) in d23)
        hydrogen_atoms = ats.get(lambda x: x.element=='H' and len(x.bonds))
        #hAts = AtomSet(ats.get(lambda x, donorTypes=donorTypes: x.element=='H' \
        hAts = AtomSet(hydrogen_atoms.get(lambda x, donorTypes=donorTypes: 
                    x.bonds[0].atom1.babel_type in donorTypes\
                    or x.bonds[0].atom2.babel_type in donorTypes))
        d['hAts'] = hAts
        d['donor2Ats'] = donor2Ats
        d['donor3Ats'] = donor3Ats
    
        acceptorTypes = paramDict['acceptorTypes']
        #print "about to call getHBAcceptors with acceptorTypes=", acceptorTypes
        acceptor2Ats, acceptor3Ats = self.getHBAcceptors(ats, acceptorTypes)
        d['acceptor2Ats'] = acceptor2Ats
        d['acceptor3Ats'] = acceptor3Ats
        if acceptor2Ats:
            acceptorAts = acceptor2Ats
            if acceptor3Ats:
                acceptorAts = acceptorAts + acceptor3Ats
        elif acceptor3Ats:
            acceptorAts = acceptor3Ats
        else:
            #CHECK THIS: should it be None or AtomSet([])
            acceptorAts = None
        d['acceptorAts'] = acceptorAts
        return d
            

    def getMat(self, ats):
        pass
        #tops = ats.top.uniq()
        #if len(tops)>1: 
        #    self.warningMsg('transformation mat=None:>1 mol in atomset!')
        #    return None
        #g = tops[0].geomContainer.geoms['master']
        #return g.GetMatrix(g)

    def process(self, dict1, dict2, paramDict):
        #hAts are keys, aceptorAts are checks
        hAts = dict1['hAts']
        tAts = hAts
        dist = paramDict['distCutoff']
        distOnly = paramDict['distOnly']

        if not hAts:
            #then use donors and a different distance
            tAts = dict1['donor2Ats'] + dict1['donor3Ats']
            dist = paramDict['distCutoff2']
            
        acceptorAts = dict2['acceptorAts']
        #print "acceptorAts=", acceptorAts
        if not acceptorAts or not tAts: #6/14/2004
            return {}, {}

        #call distanceSelector on two groups of atoms with dist
        #keyMat = self.getMat(tAts)
        #checkMat = self.getMat(acceptorAts)
        atDict = self.distSelector.select(tAts, acceptorAts, dist)
        #    keyMat=keyMat, checkMat=checkMat)
        #atDict = self.distSelector.select(tAts, acceptorAts, dist)
        #first remove bonded angles
        atDict = self.removeNeighbors(atDict)

        donor2Ats = dict1['donor2Ats']
        donor3Ats = dict1['donor3Ats']
        acceptor2Ats = dict2['acceptor2Ats']
        acceptor3Ats = dict2['acceptor3Ats']

        if distOnly: 
            #need to build hbonds and return dictionary
            self.makeBonds(atDict, donor2Ats, donor3Ats, \
                    acceptor2Ats, acceptor3Ats, paramDict)
            return atDict

        badAtDict = self.filterBasedOnAngs(atDict, donor2Ats, donor3Ats, \
                    acceptor2Ats, acceptor3Ats, paramDict)
        atDict = self.removeBadAts(atDict, badAtDict)
        if atDict is None:
            atDict = {}
        return atDict


    def makeBonds(self, pD, d2Ats, d3Ats, a2Ats, a3ats, paramDict):
        for k in pD.keys():
            if k.element=='H':
                if hasattr(k, 'hbonds') and len(k.hbonds):
                    continue
                d = k.bonds[0].atom1
                if id(d)==id(k): d = k.bonds[0].atom2
                #d = k.bonds[0].neighborAtom(k)
                h = k
            else: 
                d = k
                h = None
            #pD[k] is a list of close-enough ats
            for ac in pD[k]:
                if ac==d: continue
                dSp2 = d in d2Ats
                aSp2 = ac in a2Ats
                if dSp2:
                    if aSp2: typ = 22
                    else: typ = 23
                elif aSp2: typ = 32
                else: typ = 33
                #THEY could be already bonded
                alreadyBonded = 0
                if hasattr(d, 'hbonds') and hasattr(ac,'hbonds'):
                    for hb in d.hbonds:
                        if hb.donAt==ac or hb.accAt==ac:
                            alreadyBonded = 1
                            
                if not alreadyBonded:
                    newHB = HydrogenBond(d, ac, h, typ=typ)
                    if not hasattr(ac, 'hbonds'):
                        ac.hbonds=[]
                    if not hasattr(d, 'hbonds'):
                        d.hbonds=[]
                    ac.hbonds.append(newHB)
                    d.hbonds.append(newHB)
                    if h is not None:
                        #hydrogens can have only 1 hbond
                        h.hbonds = [newHB]


    def filterBasedOnAngs(self, pD, d2Ats, d3Ats, a2Ats, a3ats, paramDict):
        badAtDict = {}
        d2max = paramDict['d2max']
        d2min = paramDict['d2min']
        d3max = paramDict['d3max']
        d3min = paramDict['d3min']
        #NEED these parameters
        a2max = paramDict['a2max']
        a2min = paramDict['a2min']
        a3max = paramDict['a3max']
        a3min = paramDict['a3min']
        #NB now pD keys could be hydrogens OR donors
        for k in pD.keys():
            if k.element=='H':
                d = k.bonds[0].atom1
                if id(d)==id(k): d = k.bonds[0].atom2
                #d = k.bonds[0].neighborAtom(k)
                h = k
            else: 
                d = k
                h = None
            badAts = AtomSet([])
            ct = 0
            for ac in pD[k]:
                if h is not None:
                    ang = getAngle(ac, h, d)
                else:
                    acN = ac.bonds[0].atom1
                    if id(acN) == id(ac): acN = ac.bonds[0].atom2
                    #acN = ac.bonds[0].neighborAtom(ac)
                    ang = getAngle(d, ac, acN)
                #print 'ang=', ang
                dSp2 = d in d2Ats
                aSp2 = ac in a2Ats
                #these limits could be adjustable
                if h is not None:
                    if dSp2:
                        upperLim = d2max
                        lowerLim = d2min
                        #upperLim = 170
                        #lowerLim = 130
                    else:
                        upperLim = d3max
                        lowerLim = d3min
                        #upperLim = 180
                        #lowerLim = 120
                else:
                    #if there is no hydrogen use d-ac-acN angles
                    if dSp2:
                        upperLim = a2max
                        lowerLim = a2min
                        #upperLim = 150
                        #lowerLim = 110
                    else:
                        upperLim = a3max
                        lowerLim = a3min
                        #upperLim = 150
                        #lowerLim = 100
                if ang>lowerLim and ang <upperLim:
                    #AT THIS MOMENT BUILD HYDROGEN BOND:
                    if dSp2:
                        if aSp2: typ = 22
                        else: typ = 23
                    elif aSp2: typ = 32
                    else: typ = 33
                    #THEY could be already bonded
                    alreadyBonded = 0
                    if hasattr(d, 'hbonds') and hasattr(ac,'hbonds'):
                        for hb in d.hbonds:
                            if hb.donAt==ac or hb.accAt==ac:
                                alreadyBonded = 1
                    if not alreadyBonded:
                        newHB = HydrogenBond(d, ac, h, theta=ang, typ=typ)
                        if not hasattr(ac, 'hbonds'):
                            ac.hbonds=[]
                        if not hasattr(d, 'hbonds'):
                            d.hbonds=[]
                        ac.hbonds.append(newHB)
                        d.hbonds.append(newHB)
                        if h is not None:
                            #hydrogens can have only 1 hbond
                            h.hbonds = [newHB]
                        #    newHB.hlen = dist
                        #else:
                        #    newHB.dlen = dist
                else:
                    badAts.append(ac)
                ct = ct + 1
            badAtDict[k] = badAts
        return badAtDict


    def removeBadAts(self, atDict, badAtDict):
        #clean-up function called after filtering on angles
        badKeys= badAtDict.keys()
        for at in atDict.keys():
            if at not in badKeys:
                continue
            if not len(badAtDict[at]):
                continue
            closeAts = atDict[at]
            badAts = badAtDict[at]
            goodAts = []
            for i in range(len(closeAts)):
                cAt = closeAts[i]
                if cAt not in badAts:
                    goodAts.append(cAt)
            if len(goodAts):
                atDict[at] = goodAts
            else:
                del atDict[at]
        return atDict


    def removeNeighbors(self, atDict):
        #filter out at-itself and at-bondedat up to 1:4
        #NB keys could be hydrogens OR donors
        for at in atDict.keys():
            closeAts = atDict[at]
            bondedAts = AtomSet([])
            for b in at.bonds:
                ###at2 = b.neighborAtom(at)
                at2 = b.atom1
                if id(at2)==id(at): at2 = b.atom2
                bondedAts.append(at2)
                #9/13 remove this:
                ##also remove 1-3
                for b2 in at2.bonds:
                    at3 = b2.atom1
                    if id(at3)==id(at2): at3 = b.atom2
                    #at3 = b2.neighborAtom(at2)
                    if id(at3)!=id(at):
                        bondedAts.append(at3)
                    #for b3 in at3.bonds:
                        #at4 = b2.neighborAtom(at3)
                        #if at4!=at and at4!=at2:
                            #bondedAts.append(at4)
            bondedAts = bondedAts.uniq()
            goodAts = []
            for i in range(len(closeAts)):
                cAt = closeAts[i]
                if cAt not in bondedAts:
                    goodAts.append(cAt)
            if len(goodAts):
                atDict[at] = goodAts
            else:
                del atDict[at]
        return atDict


    def getDonors(self, nodes, paramDict):
        donorList = paramDict['donorTypes']
        #print 'donorList=', donorList
        # currently this is a set of hydrogens
        hats = AtomSet(nodes.get(lambda x: x.element=='H'))
        #hats are optional: if none, process donors
        # if there are hats: dAts are all atoms bonded to all hydrogens
        if hats:
            dAts = AtomSet([])
            for at in hats:
                for b in at.bonds:
                    at2 = b.atom1
                    if id(at2)==id(at): at2 = b.atom2
                    dAts.append(at2)
                    #dAts.append(b.neighborAtom(at))
        else:
            dAts = nodes
        #get the sp2 hybridized possible donors which are all ns
        sp2 = []
        for t in ['Nam', 'Ng+', 'Npl']:
            if t in donorList:
                sp2.append(t)
        #ntypes = ['Nam', 'Ng+', 'Npl']

        sp2DAts = None
        if len(sp2):
            sp2DAts = AtomSet(dAts.get(lambda x, sp2=sp2: x.babel_type in sp2))

        hsp2 = AtomSet([])
        if sp2DAts:
            if hats:
                hsp2 = AtomSet(hats.get(lambda x, sp2DAts=sp2DAts:x.bonds[0].atom1 \
                        in sp2DAts or x.bonds[0].atom2 in sp2DAts))
        if sp2DAts:
            #remove any sp2 N atoms which already have 3 bonds not to hydrogens
            n2Dons = AtomSet(sp2DAts.get(lambda x: x.element=='N'))
            if n2Dons:
                n2Dons.bl=0
                for at in n2Dons:
                    for b in at.bonds:
                        if type(b.bondOrder)==type(2):
                            at.bl = at.bl + b.bondOrder
                        else:
                            at.bl = at.bl + 2
                        #allow that there might already be a hydrogen
                    nH = at.findHydrogens()
                    at.bl = at.bl - len(nH)
                badAts = AtomSet(n2Dons.get(lambda x: x.bl>2))
                if badAts:
                    sp2DAts = sp2DAts - badAts
                delattr(n2Dons,'bl')
        #get the sp3 hybridized possible donors
        sp3 = []
        for t in ['N3+', 'S3', 'O3']:
            if t in donorList:
                sp3.append(t)
        n3DAts = None
        if 'N3+' in sp3:
            n3DAts = AtomSet(dAts.get(lambda x: x.babel_type=='N3+'))
        o3DAts = None
        if 'O3' in sp3:
            o3DAts = AtomSet(dAts.get(lambda x: x.babel_type=='O3'))
        if o3DAts:
            #remove any O3 atoms which already have 2 bonds not to hydrogens
            badO3s = AtomSet([])
            for at in o3DAts:
                if len(at.bonds)<2: continue
                if len(at.findHydrogens()): continue
                else:
                    badO3s.append(at)
            if len(badO3s):
                o3DAts = o3DAts - badO3s
        s3DAts = None
        if 'S3' in sp3:
            s3DAts = AtomSet(dAts.get(lambda x: x.babel_type=='S3'))
        sp3DAts = AtomSet([])
        for item in [n3DAts, o3DAts, s3DAts]:
            if item:
                sp3DAts = sp3DAts + item
        hsp3 = AtomSet([])
        if sp3DAts:
            if hats:
                hsp3 = AtomSet(hats.get(lambda x, sp3DAts=sp3DAts:x.bonds[0].atom1 \
                    in sp3DAts or x.bonds[0].atom2 in sp3DAts))
        hsp = hsp2 + hsp3
        #print 'hsp=', hsp.name
        #print 'sp2DAts=', sp2DAts.name
        #print 'sp3DAts=', sp3DAts.name
        return hsp, sp2DAts, sp3DAts

    def getAcceptors(self, nodes, paramDict):
        acceptorList = paramDict['acceptorTypes']
        #print 'acceptorList=', acceptorList

        sp2 = []
        for t in ['Npl', 'Nam']:
            if t in acceptorList: sp2.append(t)
        n2Accs = None
        if 'Npl' in sp2:
            n2Accs = AtomSet(nodes.get(lambda x: x.babel_type=='Npl'))
        if 'Nam' in sp2:
            n2Accs2 = AtomSet(nodes.get(lambda x: x.babel_type=='Nam'))
            if n2Accs2:
                if n2Accs:
                    n2Accs = n2Accs+n2Accs2
                else:
                    n2Accs = n2Accs2
        if n2Accs is None: 
            n2Accs = AtomSet([])

        o_sp2 = []
        for t in ['O2', 'O-']:
            if t in acceptorList: sp2.append(t)

        o2Accs = None
        if 'O2' in o_sp2:
            o2Accs = AtomSet(nodes.get(lambda x: x.babel_type=='O2'))
        if 'O-' in sp2:
            o2Accs2 = AtomSet(nodes.get(lambda x: x.babel_type=='O-'))
            if o2Accs2:
                if o2Accs:
                    o2Accs = o2Accs+o2Accs2
                else:
                    o2Accs = o2Accs2
        if o2Accs is None: 
            o2Accs = AtomSet([])

        
        o3Accs = None
        if 'O3' in acceptorList:
            o3Accs = AtomSet(nodes.get(lambda x: x.babel_type=='O3'))
        if o3Accs is None: o3Accs = AtomSet([])

        s3Accs = None
        if 'S3' in acceptorList:
            s3Accs = AtomSet(nodes.get(lambda x: x.babel_type=='S3'))
        if s3Accs is None: s3Accs = AtomSet([])

        ret2Ats = AtomSet([])
        for item in [n2Accs, o2Accs]:
            ret2Ats = ret2Ats + item

        ret3Ats = AtomSet([])
        for item in [s3Accs, o3Accs]:
            ret3Ats = ret3Ats + item
        if ret2Ats: print 'ret2Ats=', ret2Ats.name
        else: print 'no ret2Ats'
        if ret3Ats: print 'ret3Ats=', ret3Ats.name
        else: print 'no ret3Ats'
        return ret2Ats, ret3Ats

