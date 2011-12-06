## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

############################################################################
#
# Author: Ruth HUEY, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2001
#
#############################################################################

# $Header: /opt/cvs/python/packages/share1.5/MolKit/amberPrmTop.py,v 1.32 2007/07/24 17:30:40 vareille Exp $
#
# $Id: amberPrmTop.py,v 1.32 2007/07/24 17:30:40 vareille Exp $
#


#from MolKit.molecule import Atom, AtomSet, Bond
from sff.amber import AmberParm

import numpy.oldnumeric as Numeric, types
from math import pi, sqrt, ceil, fabs
from string import split, strip, join
from os.path import basename

from MolKit.data.all_amino94_dat import all_amino94_dat
from MolKit.data.all_aminont94_dat import all_aminont94_dat
from MolKit.data.all_aminoct94_dat import all_aminoct94_dat



class Parm:
    """class to hold parameters for Amber Force Field calcuations
    """
    def __init__(self, allDictList = [all_amino94_dat], ntDictList = [all_aminont94_dat],
                ctDictList = [all_aminoct94_dat]):
        #amber parameter reference dictionaries:
        if len(allDictList)==0:
            allDict = all_amino94_dat
        else:
            allDict = allDictList[0]
            if type(allDict)==types.StringType:
                allDict = self.getDictObj(allDict)
        if len(allDictList)>1:
            for d in allDictList:
                if type(d)== types.StringType:
                    d = self.getDictObj(d)
                allDict.update(d)
                #allDict.extend(d)
        self.allDict = allDict
        if len(ntDictList)==0:
            ntDict = all_aminont94_dat
        else:
            ntDict = ntDictList[0]
            if type(ntDict)==types.StringType:
                ntDict = self.getDictObj(ntDict)
        if len(ntDictList)>1:
            for d in ntDictList:
                if type(d)== types.StringType:
                    d = self.getDictObj(d)
                ntDict.update(d)
                #ntDict.extend(d)
        self.ntDict = ntDict
        if len(ctDictList)==0:
            ctDict = all_aminoct94_dat
        else:
            ctDict = ctDictList[0]
            if type(ctDict)==types.StringType:
                ctDict = self.getDictObj(ctDict)
        if len(ctDictList)>1:
            for d in ctDictList:
                if type(d)== types.StringType:
                    d = self.getDictObj(d)
                ctDict.update(d)
                #ctDict.extend(d)
        self.ctDict = ctDict
        #formatD is used for write method
        formatD = {}
        for k in ['Iac', 'Iblo', 'Cno', 'Ipres', 'ExclAt']:
            formatD[k] = ('%6d', 12, 0)
        for k in ['Charges', 'Masses', 'Rk', 'Req', 'Tk', 'Teq',\
                 'Pk', 'Pn', 'Phase', 'Solty', 'Cn1', 'Cn2']:
            formatD[k] = ('%16.8E', 5, 0)
        for k in ['AtomNames', 'ResNames', 'AtomSym', 'AtomTree']:
            formatD[k] = ('%-4.4s', 20, 0)
        for k in ['allHBnds', 'allBnds']:
            formatD[k] = ('%6d', 12, 3)
        for k in ['allHAngs', 'allAngs']:
            formatD[k] = ('%6d', 12, 4)
        for k in ['allHDihe', 'allDihe']:
            formatD[k] = ('%6d', 12, 5)
        self.formatD = formatD
        #processAtoms results are built in this dictionary
        self.prmDict = {}


    def getDictObj(self, nmstr):
        #mod = __import__('MolKit/data/' + nmstr)
        #dict = eval('mod.'+ nmstr)
        mod = __import__('MolKit')
        b = getattr(mod.data, nmstr)
        dict = getattr(b, nmstr)
        return dict




    def loadFromFile(self, filename):
        """reads a parmtop file"""
        self.prmDict = self.py_read(filename)
        self.createSffCdataStruct(self.prmDict)
        

    def processAtoms(self, atoms, parmDict=None, reorder=1):
        """finds all Amber parameters for the given set of atoms
           parmDict is parm94_dat """


        if atoms:
            self.build(atoms, parmDict, reorder)
        self.createSffCdataStruct(self.prmDict)
        print 'after call to createSffCdataStruct'


    def checkSanity(self):
        d = self.prmDict 
        #length checks:
        Natom = d['Natom']
        assert len(d['Charges']) == Natom
        assert len(d['Masses']) == Natom
        assert len(d['Iac']) == Natom
        assert len(d['Iblo']) == Natom
        assert len(d['AtomRes']) == Natom
        assert len(d['N14pairs']) == Natom
        assert len(d['TreeJoin']) == Natom

        Nres = d['Nres']
        assert len(d['Ipres']) == Nres + 1

        assert len(d['AtomNames']) == Natom * 4 + 81
        assert len(d['AtomSym']) == Natom * 4 + 81
        assert len(d['AtomTree']) == Natom * 4 + 81
        assert len(d['ResNames']) == Nres * 4 + 81

        #Ntypes is number of unique amber_types w/equiv replacement
        Ntypes = d['Ntypes']
        assert len(d['Cno']) == Ntypes**2
        assert len(d['ExclAt']) == d['Nnb']
        assert len(d['Cn1']) == Ntypes*(Ntypes+1)/2.
        assert len(d['Cn2']) == Ntypes*(Ntypes+1)/2.

        #Numbnd is number of bnd types
        Numbnd = d['Numbnd']
        assert len(d['Rk']) == Numbnd
        assert len(d['Req']) == Numbnd 
        #Numang is number of angle types
        Numang = d['Numang']
        assert len(d['Tk']) == Numang
        assert len(d['Teq']) == Numang
        #Numptra is number of dihe types
        Nptra = d['Nptra']
        assert len(d['Pk']) == Nptra
        assert len(d['Pn']) == Nptra
        assert len(d['Phase']) == Nptra

        assert len(d['Solty'])      == d['Natyp']

        #Nbona is number of bonds w/out H
        Nbona = d['Nbona']
        assert len(d['BondAt1']) == Nbona
        assert len(d['BondAt2']) == Nbona
        assert len(d['BondNum']) == Nbona
        #Nbonh is number of bonds w/ H
        Nbonh = d['Nbonh']
        assert len(d['BondHAt1']) == Nbonh
        assert len(d['BondHAt2']) == Nbonh
        assert len(d['BondHNum']) == Nbonh

        #Ntheta is number of angles w/out H
        Ntheta = d['Ntheta']
        assert len(d['AngleAt1']) == Ntheta
        assert len(d['AngleAt2']) == Ntheta
        assert len(d['AngleAt3']) == Ntheta
        assert len(d['AngleNum']) == Ntheta

        #Ntheth is number of angles w/ H
        Ntheth = d['Ntheth']
        assert len(d['AngleHAt1']) == Ntheth
        assert len(d['AngleHAt2']) == Ntheth
        assert len(d['AngleHAt3']) == Ntheth
        assert len(d['AngleHNum']) == Ntheth

        #Nphia is number of dihedrals w/out H
        Nphia = d['Nphia']
        assert len(d['DihAt1']) == Nphia
        assert len(d['DihAt2']) == Nphia
        assert len(d['DihAt3']) == Nphia
        assert len(d['DihAt4']) == Nphia
        assert len(d['DihNum']) == Nphia

        #Nphih is number of dihedrals w/ H
        Nphih = d['Nphih']
        assert len(d['DihHAt1']) == Nphih
        assert len(d['DihHAt2']) == Nphih
        assert len(d['DihHAt3']) == Nphih
        assert len(d['DihHAt4']) == Nphih
        assert len(d['DihHNum']) == Nphih

        ##WHAT ABOUT HB10, HB12, N14pairs, N14pairlist

        #value based on length checks:
        #all values of BondNum and BondHNum in range (1, Numbnd)
        for v in d['BondNum']:
            assert v >0 and v < Numbnd + 1
        for v in d['BondHNum']:
            assert v >0 and v < Numbnd + 1
        #all values of AngleNum and AngleHNum in range (1, Numang)
        for v in d['AngleNum']:
            assert v >0 and v < Numang + 1
        for v in d['AngleHNum']:
            assert v >0 and v < Numang + 1
        #all values of DihNum and DihHNum in range (1, Nptra)
        for v in d['DihNum']:
            assert v >0 and v < Nptra + 1
        for v in d['DihHNum']:
            assert v >0 and v < Nptra + 1


    def createSffCdataStruct(self, dict):
        """Create a C prm data structure"""
        print 'in createSffCdataStruct'
        self.ambPrm = AmberParm('test1', parmdict=dict)
        print 'after call to init'


    def build(self, allAtoms, parmDict, reorder):

        # find out amber special residue name and
        # order the atoms inside a residue to follow the Amber convention 
        self.residues = allAtoms.parent.uniq()

        self.residues.sort()

        self.fixResNamesAndOrderAtoms(reorder)

        # save ordered chains
        self.chains = self.residues.parent.uniq()
        self.chains.sort()

        # save ordered atoms
        self.atoms = self.residues.atoms

        # renumber them
        self.atoms.number = range(1, len(allAtoms)+1)
        print 'after call to checkRes'

        self.getTopology(self.atoms, parmDict)
        print 'after call to getTopology'

        if reorder:
            self.checkSanity()
            print 'passed sanity check'
        else:
            print 'skipping sanity check'
        return
                

    def reorderAtoms(self, res, atList):
        ats = []
        rlen = len(res.atoms)
        if rlen!=len(atList):
            print "atoms missing in residue", res
            print "expected:", atList
            print "found     :", res.atoms.name
        for i in range(rlen):
            a = atList[i]
            for j in range(rlen):
                b = res.atoms[j]
                # DON'T rename HN atom H, HN1->H1, etc...
                # use editCommands instead
                #if b.name=='HN': b.name='H'
                #elif len(b.name)==3 and b.name[:2]=='HN': 
                    #b.name ='H'+b.name[2]
                if b.name==a:
                    ats.append(b)
                    break
        if len(ats)==len(res.atoms):
            res.children.data = ats
            res.atoms.data = ats


    def fixResNamesAndOrderAtoms(self, reorder):
        # level list of atom names used to rename residues
        # check is HIS is HIS, HID, HIP, HIE, etc...

        residues = self.residues
        last = len(residues)-1
        for i in range(len(residues)):
            residue = residues[i]
            chNames = residue.atoms.name

            amberResType = residue.type
            
            if amberResType=='CYS':
                returnVal = 'CYS'
                #3/21:
                if 'HSG' in chNames or 'HG' in chNames:
                    amberResType ='CYS'
                elif 'HN' in chNames:
                    amberResType = 'CYM'
                else:
                    amberResType = 'CYX'
            elif amberResType=='LYS':
                # THIS DOESN'T SUPPORT LYH assigned in all.in
                returnVal = 'LYS'
                if 'HZ1' in chNames or 'HZN1' in chNames:
                    amberResType ='LYS'
                else:
                    amberResType ='LYN'
            elif amberResType=='ASP':
                returnVal = 'ASP'
                #3/21
                if 'HD' in chNames or 'HD2' in chNames:
                    amberResType ='ASH'
                else:
                    amberResType ='ASP'
            elif amberResType=='GLU':
                returnVal = 'GLU'
                #3/21
                if 'HE' in chNames or 'HE2' in chNames:
                    amberResType ='GLH'
                else:
                    amberResType ='GLU'
            elif amberResType=='HIS':
                returnVal = 'HIS'
                hasHD1 = 'HD1' in chNames
                hasHD2 = 'HD2' in chNames
                hasHE1 = 'HE1' in chNames
                hasHE2 = 'HE2' in chNames
                if hasHD1 and hasHE1:
                    if hasHD2 and not hasHE2:
                        amberResType = 'HID'
                    elif hasHD2 and hasHE2:
                        amberResType = 'HIP'
                elif (not hasHD1) and (hasHE1 and hasHD2 and hasHE2):
                    amberResType = 'HIE'
                else:
                    print 'unknown HISTIDINE config'
                    raise ValueError

            residue.amber_type = amberResType
            if residue == residue.parent.residues[0]:
                residue.amber_dict = self.ntDict[amberResType]
            elif residue == residue.parent.residues[-1]:
                residue.amber_dict = self.ctDict[amberResType]
            else:
                residue.amber_dict = self.allDict[amberResType]

            if reorder:
                self.reorderAtoms(residue, residue.amber_dict['atNameList'])


    def processChain(self, residues, parmDict):
        #this should be called with a list of residues representing a chain
        # NOTE: self.parmDict is parm94 which was parsed by Ruth while parmDict is
        #       MolKit.parm94_dat.py
        dict = self.prmDict
        #residues = self.residues
        
        # initialize
        atNames = ''
        atSym = ''
        atTree = ''
        resname = ''
        
        masses = dict['Masses']
        charges = dict['Charges']
        uniqList = []
        uniqTypes = {} # used to build list with equivalent names removed
        atypTypes = {} # used to build list without equivalent names removed
        allTypeList = [] # list of all types
        
        last = len(residues)-1
        dict['Nres'] = dict['Nres'] + last + 1
        atres = dict['AtomRes']
        ipres = dict['Ipres']
        maxResLen = 0
        
        for i in range(last+1):
            res = residues[i]
            atoms = res.atoms

            nbat = len(atoms)
            if nbat > maxResLen: maxResLen = nbat
            ipres.append(ipres[-1]+nbat)
            resname = resname + res.amber_type + ' '

            ad = res.amber_dict
            pdm = parmDict.atomTypes
           
            for a in atoms:
                # get the amber atom type
                name = a.name
                atres.append(i+1)
                atNames = atNames+'%-4s'%name
                atD = ad[name]
                a.amber_type = newtype = '%-2s'%atD['type']
                chg = a._charges['amber'] = atD['charge']*18.2223
                charges.append(chg)
                mas = a.mass = pdm[newtype][0]
                masses.append(mas)
                atTree = atTree+'%-4.4s'%atD['tree']
                allTypeList.append(newtype)
                atSym = atSym+'%-4s'%newtype
                symb = newtype[0]
                if symb in parmDict.AtomEquiv.keys():
                    if newtype in parmDict.AtomEquiv[symb]:
                        newsym = symb + ' '
                        uniqTypes[symb+' '] = 0
                        a.amber_symbol = symb+' '
                        if newsym not in uniqList:
                            uniqList.append(newsym)
                    else:
                        uniqTypes[newtype] = 0
                        a.amber_symbol = newtype
                        if newtype not in uniqList:
                            uniqList.append(newtype)
                else:
                    uniqTypes[newtype] = 0
                    a.amber_symbol = newtype
                    if newtype not in uniqList: uniqList.append(newtype)
                # to get uniq list of all types w/out equiv replacement
                atypTypes[newtype] = 0
                    
        # post processing of some variable
        dict['AtomNames'] = dict['AtomNames'] + atNames
        dict['AtomSym'] = dict['AtomSym'] + atSym
        dict['AtomTree'] = dict['AtomTree'] + atTree
        dict['ResNames'] = dict['ResNames'] + resname

        # save list of unique types for later use
        ###1/10:
        #self.uniqTypeList = uniqList
        uL = self.uniqTypeList
        for t in uniqList:
            if t not in uL:
                uL.append(t)
        #self.uniqTypeList = uniqTypes.keys()
        self.uniqTypeList = uL

        ntypes = len(uL)
        dict['Ntypes'] = ntypes

        aL = self.atypList
        for t in atypTypes.keys():
            if t not in aL:
                aL.append(t)
        self.atypList = aL
        dict['Natyp'] = len(aL)

        dict['Ntype2d'] = ntypes*ntypes
        dict['Nttyp'] = ntypes * (ntypes+1)/2
        if maxResLen > dict['Nmxrs']:
            dict['Nmxrs'] = maxResLen


        newtypelist = []
        for t in residues.atoms.amber_symbol:
            # Iac is 1-based
            newtypelist.append( self.uniqTypeList.index(t) + 1 )
        ###1/10:
        #dict['Iac'] = newtypelist
        dict['Iac'].extend( newtypelist)


    def processBonds(self, bonds, parmDict):
        # NOTE: self,parmDict is parm94 parsed by Ruth while parmDict is
        #       MolKit.parm94_dat.py):
        
        dict = self.prmDict
        bat1 = dict['BondAt1']
        bat2 = dict['BondAt2']
        bnum = dict['BondNum']
        batH1 = dict['BondHAt1']
        batH2 = dict['BondHAt2']
        bHnum = dict['BondHNum']
        rk = dict['Rk']
        req = dict['Req']

        bndTypes = {} # used to build a unique list of bond types
        btDict = parmDict.bondTypes #needed to check for wildcard * in type
        
        for b in bonds:
            a1 = b.atom1
            #t1 = a1.amber_symbol
            t1 = a1.amber_type
            a2 = b.atom2
            #t2 = a2.amber_symbol
            t2 = a2.amber_type
            if t1<t2:
                newtype = '%-2.2s-%-2.2s'%(t1,t2)
            else:
                newtype = '%-2.2s-%-2.2s'%(t2,t1)
            bndTypes[newtype] = 0

            n1 = (a1.number-1)*3
            n2 = (a2.number-1)*3
            if n2<n1: tmp=n1; n1=n2; n2=tmp

            if a1.element=='H' or a2.element=='H':
                bHnum.append(newtype)
                batH1.append(n1)
                batH2.append(n2)
            else:
                bnum.append(newtype)
                bat1.append(n1)
                bat2.append(n2)

        dict['Numbnd'] = len(bndTypes)
        btlist = bndTypes.keys()

        for bt in btlist:
            rk.append( btDict[bt][0] )
            req.append( btDict[bt][1] )

        newbnum = []
        for b in bnum:
            newbnum.append( btlist.index(b) + 1 )
        dict['BondNum'] = newbnum

        newbnum = []
        for b in bHnum:
            newbnum.append( btlist.index(b) + 1 )
        dict['BondHNum'] = newbnum
        
        return


    def processAngles(self, allAtoms, parmDict):
        dict = self.prmDict        
        aa1 = dict['AngleAt1']
        aa2 = dict['AngleAt2']
        aa3 = dict['AngleAt3']
        anum = dict['AngleNum']
        aHa1 = dict['AngleHAt1']
        aHa2 = dict['AngleHAt2']
        aHa3 = dict['AngleHAt3']
        aHnum = dict['AngleHNum']
        tk = dict['Tk']
        teq = dict['Teq']

        angTypes = {}
        atdict = parmDict.bondAngles
        
        for a1 in allAtoms:
            t1 = a1.amber_type
            for b in a1.bonds:
                a2 = b.atom1
                if id(a2)==id(a1): a2=b.atom2
                t2 = a2.amber_type
                for b2 in a2.bonds:
                    a3 = b2.atom1
                    if id(a3)==id(a2): a3=b2.atom2
                    if id(a3)==id(a1): continue
                    if a1.number > a3.number: continue
                    
                    t3 = a3.amber_type

                    nn1 = n1 = (a1.number-1)*3
                    nn2 = n2 = (a2.number-1)*3
                    nn3 = n3 = (a3.number-1)*3
                    if n3<n1:
                        nn3 = n1
                        nn1 = n3

                    rev = 0
                    if (t1==t3 and a1.name > a3.name) or t3 < t1:
                        rev = 1

                    if rev:
                        newtype = '%-2.2s-%-2.2s-%-2.2s'%(t3,t2,t1)
                    else:
                        newtype = '%-2.2s-%-2.2s-%-2.2s'%(t1, t2, t3)
                    #have to check for wildcard *
                    angTypes[newtype] = 0

                    if a1.element=='H' or a2.element=='H' or a3.element=='H':
                        aHa1.append( nn1 )
                        aHa2.append( nn2 )
                        aHa3.append( nn3 )
                        aHnum.append(newtype)
                    else:
                        aa1.append( nn1 )
                        aa2.append( nn2 )
                        aa3.append( nn3 )
                        anum.append(newtype)

        atlist = angTypes.keys()

        torad = pi / 180.0
        atKeys = atdict.keys()
        for t in atlist:
            tk.append( atdict[t][0] )
            teq.append( atdict[t][1]*torad )

        anewlist = []
        for a in anum:
            anewlist.append( atlist.index( a ) + 1 )
        dict['AngleNum'] = anewlist

        anewlist = []
        for a in aHnum:
            anewlist.append( atlist.index( a ) + 1 )
        dict['AngleHNum'] = anewlist

        dict['Numang'] = len(atlist)
        dict['Ntheth'] = len(aHa1)
        dict['Mtheta'] = len(aa1)
        dict['Ntheta'] = len(aa1)
        
        return


    def checkDiheType(self, t, t2, t3, t4, dict):
        #zero X
        newtype = '%-2.2s-%-2.2s-%-2.2s-%-2.2s'%(t,t2,t3,t4)
        if dict.has_key(newtype): return newtype

        newtype = '%-2.2s-%-2.2s-%-2.2s-%-2.2s'%(t4,t3,t2,t)
        if dict.has_key(newtype): return newtype

        #X
        newtype = '%-2.2s-%-2.2s-%-2.2s-%-2.2s'%('X',t2,t3,t4)
        if dict.has_key(newtype): return newtype

        newtype = '%-2.2s-%-2.2s-%-2.2s-%-2.2s'%('X',t3,t2,t)
        if dict.has_key(newtype): return newtype

        #2X
        newtype = '%-2.2s-%-2.2s-%-2.2s-%-2.2s'%('X',t2,t3,'X')
        if dict.has_key(newtype): return newtype

        newtype = '%-2.2s-%-2.2s-%-2.2s-%-2.2s'%('X',t3,t2,'X')
        if dict.has_key(newtype): return newtype

        newtype = '%-2.2s-%-2.2s-%-2.2s-%-2.2s'%('X','X',t3,t4)
        if dict.has_key(newtype): return newtype

        newtype = '%-2.2s-%-2.2s-%-2.2s-%-2.2s'%('X','X',t2,t)
        if dict.has_key(newtype): return newtype

        raise RuntimeError('dihedral type not in dictionary')

## it is slower to check a list if the key is in there than to ask a
## dictionanry if it has this key
##
##          keys = dict.keys()
##          #zero X
##          newtype = '%-2.2s-%-2.2s-%-2.2s-%-2.2s'%(t,t2,t3,t4)
##          if newtype in keys:
##              return newtype
##          newtype2 = '%-2.2s-%-2.2s-%-2.2s-%-2.2s'%(t4,t3,t2,t)
##          if newtype2 in keys:
##              return newtype2
##          #X
##          newtypeX = '%-2.2s-%-2.2s-%-2.2s-%-2.2s'%('X',t2,t3,t4)
##          if newtypeX in keys:
##              return newtypeX
##          newtype2X = '%-2.2s-%-2.2s-%-2.2s-%-2.2s'%('X',t3,t2,t)
##          if newtype2X in keys:
##              return newtype2X
##          #2X
##          newtypeX_X = '%-2.2s-%-2.2s-%-2.2s-%-2.2s'%('X',t2,t3,'X')
##          if newtypeX_X in keys:
##              return newtypeX_X
##          newtype2X_X = '%-2.2s-%-2.2s-%-2.2s-%-2.2s'%('X',t3,t2,'X')
##          if newtype2X_X in keys:
##              return newtype2X_X
##          newtypeXX = '%-2.2s-%-2.2s-%-2.2s-%-2.2s'%('X','X',t3,t4)
##          if newtypeXX in keys:
##              return newtypeXX
##          newtype2XX = '%-2.2s-%-2.2s-%-2.2s-%-2.2s'%('X','X',t2,t)
##          if newtype2XX in keys:
##              return newtype2XX
##          raise RuntimError('dihedral type not in dictionary')


    def processTorsions(self, allAtoms, parmDict):
        # find torsions and also excuded atoms
        dict = self.prmDict
        foundDihedTypes = {}
        ta1 = dict['DihAt1']
        ta2 = dict['DihAt2']
        ta3 = dict['DihAt3']
        ta4 = dict['DihAt4']
        tnum = dict['DihNum']
        taH1 = dict['DihHAt1']
        taH2 = dict['DihHAt2']
        taH3 = dict['DihHAt3']
        taH4 = dict['DihHAt4']
        tHnum = dict['DihHNum']

        nb14 = dict['N14pairs']
        n14list = dict['N14pairlist']

        iblo = dict['Iblo']
        exclAt =  dict['ExclAt']

        dihedTypes = parmDict.dihedTypes
        
        for a1 in allAtoms:
            n14 = []
            excl = []
            t1 = a1.amber_type
            restyp = a1.parent.type
            if restyp in ['PRO', 'TRP', 'HID', 'HIE', 'HIP']:
                ringlist = self.AA5rings[restyp]
            else:
                ringlist = None

            for b in a1.bonds:
                a2 = b.atom1
                if id(a2)==id(a1): a2=b.atom2 
                t2 = a2.amber_type
                if a2.number > a1.number: excl.append(a2.number)
                for b2 in a2.bonds:
                    a3 = b2.atom1
                    if id(a3)==id(a2): a3=b2.atom2
                    if id(a3)==id(a1): continue
                    if a3.number > a1.number: excl.append(a3.number)
                    t3 = a3.amber_type
                    for b3 in a3.bonds:
                        a4 = b3.atom1
                        if id(a4)==id(a3): a4=b3.atom2
                        if id(a4)==id(a2): continue
                        if id(a4)==id(a1): continue
                        if a1.number > a4.number: continue
                        excl.append(a4.number)
                        t4 = a4.amber_type

                        newtype = '%-2.2s-%-2.2s-%-2.2s-%-2.2s'%(t1,t2,t3,t4)
                        dtype = self.checkDiheType(t1,t2,t3,t4,dihedTypes)

                        for i in range(len(dihedTypes[dtype])):
                            tname = dtype+'_'+str(i)
                            foundDihedTypes[tname] = 0
                            sign3 = 1
                            period = dihedTypes[dtype][i][3]
                            if period < 0.0: sign3= -1
                            if a4.parent==a1.parent:
                                if ringlist and a4.name in ringlist \
                                   and a1.name in ringlist:
                                    sign3= -1
                            if a1.element=='H' or a2.element=='H' or \
                               a3.element=='H' or a4.element=='H':
                                taH1.append( (a1.number-1)*3 )
                                taH2.append( (a2.number-1)*3 )
                                taH3.append( sign3*(a3.number-1)*3 )
                                taH4.append( (a4.number-1)*3 )
                                tHnum.append( tname )
                            else:
                                ta1.append( (a1.number-1)*3 )
                                ta2.append( (a2.number-1)*3 )
                                ta3.append( sign3*(a3.number-1)*3 )
                                ta4.append( (a4.number-1)*3 )
                                tnum.append( tname )
                            if sign3>0.0:
                                # this trick work only for 6 rings and
                                # prevents from adding 14 interactions
                                # twice between atoms in the ring
                                # PRO, TRP and HIS and cp. have to be handle
                                # separately
                                num = a4.number-1
                                if num not in n14:
                                    n14.append( num )
                                else: # make 3rd atom in torsion negative
                                   ta3[-1] =  -ta3[-1]
            if len(excl):
                # excl can contain duplicated values (pro tyr phe cycles)
                # we also sort the values (probably only comsetics)
                excl.sort()
                last = excl[0]
                uexcl = [last]
                for i in range(1,len(excl)):
                    if excl[i]!=last:
                        last = excl[i]
                        uexcl.append(last)
                iblo.append(len(uexcl))
                exclAt.extend(uexcl)
            else:
                iblo.append( 1 )
                exclAt.append( 0 )
                
            nb14.append(len(n14))
            ##!##1/28: n14.sort()
            n14list.extend(n14)

        # remember how many proper diehedrals
        lastProper = len(tnum)
        lastHProper = len(tHnum)

        # loop over residues to add improper torsions
        sumAts = 0
        foundImproperDihedTypes = {}
        for res in self.residues:
            foundImproperDihedTypes = self.getImpropTors(
                res, sumAts, foundImproperDihedTypes, parmDict)
            sumAts = sumAts + len(res.atoms)

        #typeDict = foundDihedTypes.copy()
        #typeDict.update(foundImproperDihedTypes)
        #print typeDict.keys()
        
        dict['Nptra'] = len(foundDihedTypes) + len(foundImproperDihedTypes)
        dict['Mphia'] = dict['Nphia'] = len(ta1)
        dict['Nphih'] = len(taH1)
        
        pn = dict['Pn']
        pk = dict['Pk']
        phase = dict['Phase']
        dtlist = foundDihedTypes.keys()
        torad = pi/180.
        for t in dtlist:
            index = int(t[-1])
            val = dihedTypes[t[:-2]][index] # remove the '_x'
            pk.append(val[1]/val[0])
            phase.append(val[2]*torad)
            pn.append(fabs(val[3]))

        dihedTypes = parmDict.improperDihed
        dtlist1 = foundImproperDihedTypes.keys()

        for t in dtlist1:
            val = dihedTypes[t]
            pk.append(val[0])
            phase.append(val[1]*torad)
            pn.append(val[2])

        typenum = []
        dtlist = dtlist + dtlist1
        for t in tnum:
            typenum.append( dtlist.index(t) + 1 ) # types are 1-based
        dict['DihNum'] = typenum

        typenum = []
        for t in tHnum:
            typenum.append( dtlist.index(t) + 1 ) # types are 1-based
        dict['DihHNum'] = typenum

        dict['Nnb'] = len(dict['ExclAt'])
        #print len(tnum), len(dict['DihNum'])

        return

        

    def getImpropTors(self, res, sumAts, foundDihedTypes, parmDict):
        #eg tList:[['CA','+M','C','0'],['-M','CA','N','H']]
        dict = self.prmDict
        offset = sumAts * 3
        nameList = res.atoms.name
        typeList = res.atoms.amber_type

        ta1 = dict['DihAt1']
        ta2 = dict['DihAt2']
        ta3 = dict['DihAt3']
        ta4 = dict['DihAt4']
        tnum = dict['DihNum']
        taH1 = dict['DihHAt1']
        taH2 = dict['DihHAt2']
        taH3 = dict['DihHAt3']
        taH4 = dict['DihHAt4']
        tHnum = dict['DihHNum']

        dihedTypes = parmDict.improperDihed
        atNameList = res.amber_dict['atNameList']
        resat = res.atoms
        for item in res.amber_dict['impropTors']:
            atomNum = []
            atomType = []
            newTors = []
            offset = res.atoms[0].number
            #use hasH to detect 'HZ2' etc
            hasH = 0
            for t in item:
                if t[0]=='H': hasH = 1
                if len(t)==2 and t[1]=='M':
                    if t[0]=='-':
                        atomType.append('C ')
                        atomNum.append(offset - 2)
                    else:
                        atomType.append('N ')
                        atomNum.append(offset + len(res.atoms) )
                else:
                    atIndex = atNameList.index(t)
                    atom = resat[atIndex]
                    atomType.append(atom.amber_type)
                    atomNum.append( atom.number )

            newType = self.checkDiheType(atomType[0], atomType[1],
                                         atomType[2], atomType[3],
                                         dihedTypes)
            foundDihedTypes[newType] = 0

            if hasH:
                taH1.append( (atomNum[0]-1)*3 )
                taH2.append( (atomNum[1]-1)*3 )
                taH3.append(-(atomNum[2]-1)*3 )
                taH4.append(-(atomNum[3]-1)*3 )
                tHnum.append(newType)
            else:
                ta1.append( (atomNum[0]-1)*3 )
                ta2.append( (atomNum[1]-1)*3 )
                ta3.append(-(atomNum[2]-1)*3 )
                ta4.append(-(atomNum[3]-1)*3 )
                tnum.append(newType)

        return foundDihedTypes
    

    def getTopology(self, allAtoms, parmDict):

        dict = self.prmDict
        dict['ititl'] = allAtoms.top.uniq()[0].name + '.prmtop\n'

        natom = dict['Natom'] = len(allAtoms)
        dict['Nat3'] = natom * 3

        dict['AtomNames'] = ''
        dict['AtomSym'] = ''
        dict['AtomTree'] = ''
        dict['Ntypes'] = 0
        dict['Natyp'] = 0
        dict['Ntype2d'] = 0
        dict['Nttyp'] = 0
        dict['Masses'] = []
        dict['Charges'] = []
        dict['Nres'] = 0
        dict['AtomRes'] = []
        dict['ResNames'] = ''
        dict['Ipres'] = [1]
        dict['Nmxrs'] = 0
        ###1/10:
        dict['Iac'] = []
        self.uniqTypeList = []
        #used for construction of Natyp
        self.atypList = []
        
        # fill get all arrays that are of len natom
        # we have to call for each chain
        for ch in self.chains:
            self.processChain( ch.residues, parmDict)

        #PAD AtomNames with 81 spaces
        dict['AtomNames'] = dict['AtomNames'] + 81*' '
        dict['AtomSym'] = dict['AtomSym'] + 81*' '
        dict['AtomTree'] = dict['AtomTree'] + 81*' '
        dict['ResNames'] = dict['ResNames'] + 81*' '

        # create Iac list
        #iac = []
        #tl = self.uniqTypeList
        #for a in allAtoms:
        #    iac.append( tl.index(a.amber_symbol) + 1 )
        #    delattr(a, 'amber_symbol')
        #dict['Iac'] = iac
        
        # to find out the number of bonds with hydrogen we simply count the
        # number of hydrogen atoms
        hlist = allAtoms.get(lambda x: x.element=='H')
        if hlist is not None and len(hlist):
            dict['Nbonh'] = numHs = len(hlist)
        else:
            numHs = 0

        # number of bonds not involving an H atom
        bonds = allAtoms.bonds[0]
        dict['Mbona'] = len(bonds) - numHs

        # since no bonds are constrined, Nbona==Mbona
        dict['Nbona'] = dict['Mbona']

        print 'after call to processChain'

        # new process bond info
        dict['BondAt1'] = []
        dict['BondAt2'] = []
        dict['BondNum'] = []
        dict['BondHAt1'] = []
        dict['BondHAt2'] = []
        dict['BondHNum'] = []
        dict['Rk'] = []
        dict['Req'] = []
        self.processBonds(bonds, parmDict)
        print 'after call to processBonds'

        # now process the angles
        dict['AngleAt1'] = []
        dict['AngleAt2'] = []
        dict['AngleAt3'] = []
        dict['AngleNum'] = []
        dict['AngleHAt1'] = []
        dict['AngleHAt2'] = []
        dict['AngleHAt3'] = []
        dict['AngleHNum'] = []
        dict['Tk'] = []
        dict['Teq'] = []
        self.processAngles(allAtoms, parmDict)
        print 'after call to processAngles'
        
        # now handle the torsions
        dict['Nhparm'] = 0
        dict['Nparm'] = 0
        
        dict['DihAt1'] = []
        dict['DihAt2'] = []
        dict['DihAt3'] = []
        dict['DihAt4'] = []
        dict['DihNum'] = []
        dict['DihHAt1'] = []
        dict['DihHAt2'] = []
        dict['DihHAt3'] = []
        dict['DihHAt4'] = []
        dict['DihHNum'] = []
        dict['Pn'] = []
        dict['Pk'] = []
        dict['Phase'] = []

        dict['Nphih'] = dict['Mphia'] = dict['Nphia'] = dict['Nptra'] = 0
        dict['N14pairs'] = []
        dict['N14pairlist'] = []

        dict['Nnb'] =0 
        dict['Iblo'] = []
        dict['ExclAt'] = []
        # FIXME
        self.AA5rings ={
            'PRO':['N', 'CA', 'CB', 'CG', 'CD'],
            'TRP':['CG', 'CD1', 'CD2', 'NE1', 'CE2'],
            'HID':['CG', 'ND1', 'CE1', 'NE2', 'CD2'],
            'HIE':['CG', 'ND1', 'CE1', 'NE2', 'CD2'],
            'HIP':['CG', 'ND1', 'CE1', 'NE2', 'CD2']
        }
        self.processTorsions(allAtoms, parmDict)
        print 'after call to processTorsions'

        # some unused values
        dict['Nspm'] = 1
        dict['Box'] = [0., 0., 0.]
        dict['Boundary'] = [natom]
        dict['TreeJoin'] = range(natom)
        dict['Nphb'] = 0
        dict['HB12'] = []
        dict['HB10'] = []
        llist = ['Ifpert', 'Nbper','Ngper','Ndper','Mbper', 'Mgper',
            'Mdper','IfBox', 'IfCap', 'Cutcap', 'Xcap', 'Ycap',
            'Zcap', 'Natcap','Ipatm', 'Nspsol','Iptres']
        for item in llist:
            dict[item] = 0

        dict['Cno'] = self.getICO( dict['Ntypes'] )
        dict['Solty'] = self.getSOLTY( dict['Natyp'] )

        dict['Cn1'], dict['Cn2'] = self.getCNList(parmDict)
        
        return


    def getICO(self, ntypes):
        ct = 1
        icoArray = Numeric.zeros((ntypes, ntypes), 'i')
        for i in range(1, ntypes+1):
            for j in range(1, i+1):
                icoArray[i-1,j-1]=ct
                icoArray[j-1,i-1]=ct
                ct = ct+1
        return icoArray.ravel().tolist()


    def getSOLTY(self, ntypes):
        soltyList = []
        for i in range(ntypes):
            soltyList.append(0.)
        return soltyList


    def getCN(self, type1, type2, pow, parmDict, factor=1):
        #pow is 12 or 6 
        #factor is 1 except when pow is 6
        d = parmDict.potParam
        if type1=='N3': type1='N '
        if type2=='N3': type2='N '
        r1, eps1 = d[type1][:2]
        r2, eps2 = d[type2][:2]
        eps = sqrt(eps1*eps2)
        rij = r1 + r2
        newval = factor*eps*rij**pow
        return newval


    def getCNList(self, parmDict):
        ntypes = len(self.uniqTypeList)
        ct = 1
##          size = self.prmDict['Nttyp']
##          cn1List = [0]*size
##          cn2List = [0]*size
##          iac = self.prmDict['Iac']
##          cno = self.prmDict['Cno']
##          for i in range(ntypes):
##              indi = i*ntypes
##              ival = self.uniqTypeList[i]
##              for j in range(i, ntypes):
##                  jval = self.uniqTypeList[j]
##                  ind = cno[indi+j]-1
##                  cn1List[ind] = self.getCN(jval, ival, 12, parmDict)
##                  cn2List[ind] = self.getCN(jval, ival, 6, parmDict, 2)
        nttyp = self.prmDict['Nttyp']
        cn1List = []
        cn2List = []
        for j in range(ntypes):
            jval = self.uniqTypeList[j]
            for i in range(j+1):
                ival = self.uniqTypeList[i]
                cn1List.append(self.getCN(ival, jval, 12, parmDict))
                cn2List.append(self.getCN(ival, jval, 6, parmDict, 2))

        return cn1List, cn2List


    def readSummary(self, allLines, dict):
        #set summary numbers
        ll = split(allLines[1])
        assert len(ll)==12
        #FIX THESE NAMES!
        natom = dict['Natom'] = int(ll[0])
        ntypes = dict['Ntypes'] = int(ll[1])
        nbonh = dict['Nbonh'] = int(ll[2])
        dict['Mbona'] = int(ll[3])
        ntheth = dict['Ntheth'] = int(ll[4])
        dict['Mtheta'] = int(ll[5])
        nphih = dict['Nphih'] = int(ll[6])
        dict['Mphia'] = int(ll[7])
        dict['Nhparm'] = int(ll[8])
        dict['Nparm'] = int(ll[9])
        #called 'next' in some documentation
        #NEXT-> Nnb
        next = dict['Nnb'] = int(ll[10])
        dict['Nres'] = int(ll[11])
        ll = split(allLines[2])
        assert len(ll)==12
        nbona = dict['Nbona'] = int(ll[0])
        ntheta = dict['Ntheta'] = int(ll[1])
        nphia = dict['Nphia'] = int(ll[2])
        numbnd = dict['Numbnd'] = int(ll[3])
        numang = dict['Numang'] = int(ll[4])
        numptra = dict['Nptra'] = int(ll[5])
        natyp = dict['Natyp'] = int(ll[6])
        dict['Nphb'] = int(ll[7])
        dict['Ifpert'] = int(ll[8])
        dict['Nbper'] = int(ll[9])
        dict['Ngper'] = int(ll[10])
        dict['Ndper'] = int(ll[11])
        ll = split(allLines[3])
        assert len(ll)==6
        dict['Mbper'] = int(ll[0])
        dict['Mgper'] = int(ll[1])
        dict['Mdper'] = int(ll[2])
        dict['IfBox'] = int(ll[3])
        dict['Nmxrs'] = int(ll[4])
        dict['IfCap'] = int(ll[5])
        return dict
        

    def readIGRAPH(self, allLines, numIGRAPH, ind=3):
        #the names are not necessarily whitespace delimited
        igraph = []
        for i in range(numIGRAPH):
            ind = ind + 1
            l = allLines[ind]
            for k in range(20):
                it = l[k*4:k*4+4]
                igraph.append(strip(it))
            #igraph.extend(split(l))
        return igraph, ind


    def readCHRG(self, allLines, ind, numCHRG, natom):
        chrg = []
        ct = 0
        for i in range(numCHRG):
            ind = ind + 1
            l = allLines[ind]
            newl = []
            # build 5 charges per line if enough are left
            #otherwise, build the last line's worth
            if natom - ct >=5:
                rct = 5
            else:
                rct = natom - ct
            for q in range(rct):
                lindex = q*16
                item = l[lindex:lindex+16]
                newl.append(float(item))
                ct = ct + 1
            chrg.extend(newl)
        return chrg, ind


    def readNUMEX(self, allLines, ind, numIAC):
        numex = []
        NumexSUM = 0
        for i in range(numIAC):
            ind = ind + 1
            ll = split(allLines[ind])
            newl = []
            for item in ll:
                newent = int(item)
                newl.append(newent)
                NumexSUM = NumexSUM + newent
            numex.extend(newl)
        return numex, ind, NumexSUM


    def readLABRES(self, allLines, ind):
        done = 0
        labres = []
        while not done:
            ind = ind + 1
            ll = split(allLines[ind])
            try:
                int(ll[0])
                done = 1
                break
            except ValueError:
                labres.extend(ll)
        #correct for 1 extra line read here
        ind = ind - 1
        return labres, ind
    

    def readFList(self, allLines, ind, numITEMS):
        v = []
        for i in range(numITEMS):
            ind = ind + 1
            ll = split(allLines[ind])
            newl = []
            for item in ll:
                newl.append(float(item))
            v.extend(newl)
        return v, ind


    def readIList(self, allLines, ind, numITEMS):
        v = []
        for i in range(numITEMS):
            ind = ind + 1
            ll = split(allLines[ind])
            newl = []
            for item in ll:
                newl.append(int(item))
            v.extend(newl)
        return v, ind


    def readILList(self, allLines, ind, numITEMS, n):
        bhlist = []
        for i in range(n):
            bhlist.append([])
        ct = 0
        for i in range(numITEMS):
            ind = ind + 1
            ll = split(allLines[ind])
            for j in range(len(ll)):
                item = ll[j]
                newl = bhlist[ct%n]
                newl.append(int(item))
                ct = ct + 1
        return bhlist, ind


    def py_read(self, filename, **kw):
        #??dict['Iptres'] #dict['Nspsol'] #dict['Ipatm'] #dict['Natcap']
        f = open(filename, 'r')
        allLines = f.readlines()
        f.close()
        dict = {}
        #set title
        dict['ititl'] = allLines[0]

        #get summary numbers
        dict = self.readSummary(allLines, dict)

        #set up convenience fields:
        natom = dict['Natom']
        ntypes = dict['Ntypes']
        dict['Nat3'] = natom * 3
        dict['Ntype2d'] = ntypes ** 2
        nttyp = dict['Nttyp'] = ntypes * (ntypes+1)/2

        # read IGRAPH->AtomNames
        numIGRAPH = int(ceil((natom*1.)/20.))
        anames, ind = self.readIGRAPH(allLines, numIGRAPH)
        dict['AtomNames'] = join(anames)

        # read CHRG->Charges
        numCHRG = int(ceil((natom*1.)/5.))
        dict['Charges'], ind = self.readCHRG(allLines, ind, numCHRG, natom)

        # read AMASS **same number of lines as charges->Masses
        dict['Masses'], ind = self.readFList(allLines, ind, numCHRG)

        # read IAC **NOT same number of lines as IGRAPH 12!!
        numIAC = int(ceil((natom*1.)/12.))
        dict['Iac'], ind  = self.readIList(allLines, ind, numIAC)

        # read NUMEX **same number of lines as IAC
        dict['Iblo'], ind, NumexSUM = self.readNUMEX(allLines, ind, numIAC)

        # read ICO *Ntype2d/12
        numICO = int(ceil((ntypes**2*1.0)/12.))
        dict['Cno'], ind = self.readIList(allLines, ind, numICO)

        ##NB this should be half of a matrix
        # read LABRES....no way to know how many
        dict['ResNames'], ind = self.readLABRES(allLines, ind)
        labres = dict['ResNames']

        # read IPRES....depends on len of LABRES
        numIPRES = int(ceil((len(labres)*1.)/20.))
        dict['Ipres'], ind = self.readIList(allLines, ind, numIPRES)

        # read RK + REQ-> depend on numbnd
        numbnd = dict['Numbnd']
        numRK = int(ceil((numbnd*1.)/5.))
        dict['Rk'], ind = self.readFList(allLines, ind, numRK)
        dict['Req'], ind = self.readFList(allLines, ind, numRK)

        # read TK + TEQ-> depend on numang
        numang = dict['Numang']
        numTK = int(ceil((numang*1.)/5.))
        dict['Tk'], ind = self.readFList(allLines, ind, numTK)
        dict['Teq'], ind = self.readFList(allLines, ind, numTK)

        # read PK, PN + PHASE-> depend on numptra
        nptra = dict['Nptra']
        numPK = int(ceil((nptra*1.)/5.))
        dict['Pk'], ind = self.readFList(allLines, ind, numPK)
        dict['Pn'], ind = self.readFList(allLines, ind, numPK)
        dict['Phase'], ind = self.readFList(allLines, ind, numPK)

        # read SOLTY
        natyp = dict['Natyp']
        numSOLTY = int(ceil((natyp*1.)/5.))
        dict['Solty'], ind = self.readFList(allLines, ind, numSOLTY)

        # read CN1 and CN2
        numCN = int(ceil((nttyp*1.)/5.))
        dict['Cn1'], ind = self.readFList(allLines, ind, numCN)
        dict['Cn2'], ind = self.readFList(allLines, ind, numCN)

        # read IBH, JBH, ICBH 12
        nbonh = dict['Nbonh']
        numIBH = int(ceil((nbonh*3.0)/12.))
        [dict['BondHAt1'], dict['BondHAt2'], dict['BondHNum']], ind = \
            self.readILList(allLines, ind, numIBH, 3)
        # read IB, JB, ICB 12
        nbona = dict['Nbona']
        numIB = int(ceil((nbona*3.0)/12.))
        [dict['BondAt1'], dict['BondAt2'], dict['BondNum']], ind = \
            self.readILList(allLines, ind, numIB, 3)

        # read ITH, JTH, KTH, ICTH 12
        ntheth = dict['Ntheth']
        numITH = int(ceil((ntheth*4.0)/12.))
        [dict['AngleHAt1'], dict['AngleHAt2'], dict['AngleHAt3'],\
            dict['AngleHNum']], ind = self.readILList(allLines, ind, numITH, 4)
        # read IT, JT, KT, ICT 12
        ntheta = dict['Ntheta']
        numIT = int(ceil((ntheta*4.0)/12.))
        [dict['AngleAt1'], dict['AngleAt2'], dict['AngleAt3'],\
            dict['AngleNum']], ind = self.readILList(allLines, ind, numIT, 4)

        # read IPH, JPH, KPH, LPH, ICPH 12
        nphih = dict['Nphih']
        numIPH = int(ceil((nphih*5.0)/12.))
        [dict['DihHAt1'], dict['DihHAt2'], dict['DihHAt3'], dict['DihHAt4'],\
            dict['DihHNum']], ind = self.readILList(allLines, ind, numIPH, 5)
        # read IP, JP, KP, LP, ICP 12
        nphia = dict['Nphia']
        numIP = int(ceil((nphia*5.0)/12.))
        [dict['DihAt1'], dict['DihAt2'], dict['DihAt3'], dict['DihAt4'],\
            dict['DihNum']], ind = self.readILList(allLines, ind, numIP, 5)

        # read NATEX 12
        #FIX THIS: has to be the sum of previous entries
        numATEX = int(ceil((NumexSUM*1.0)/12.))
        dict['ExclAt'], ind = self.readIList(allLines, ind, numATEX)

        # read CN1 and CN2
        # skip ASOL
        # skip BSOL
        # skip HBCUT
        ind = ind + 3
        # read ISYMBL 20
        asym, ind = self.readIGRAPH(allLines, numIGRAPH, ind)
        dict['AtomSym'] = join(asym)

        # read ITREE 20
        atree, ind = self.readIGRAPH(allLines, numIGRAPH, ind)
        dict['AtomTree'] = join(atree)
        return dict


    def makeList(self, llist, num):
        newL = []
        for i in range(len(llist[0])):
            ni = []
            for j in range(num):
                ni.append(llist[j][i])
            newL.append(ni)
        return newL

            
    #functions to write self
    def write(self, filename, **kw):
        fptr = open(filename, 'w')
        dict = self.prmDict
        self.writeItitl(fptr, dict['ititl'])
        self.writeSummary(fptr)
        #WHAT ABOUT SOLTY???
        self.writeString(fptr,dict['AtomNames'][:-81])
        for k in ['Charges', 'Masses', 'Iac','Iblo','Cno']:
            item = dict[k]
            f = self.formatD[k]
            if f[2]:
                self.writeTupleList(fptr, item, f[0], f[1], f[2])
            else:
                self.writeList(fptr, item, f[0], f[1])
        self.writeString(fptr,dict['ResNames'][:-81])
        self.writeList(fptr, dict['Ipres'][:-1], '%6d', 12 )
        for k in ['Rk', 'Req', 'Tk', 'Teq',
            'Pk', 'Pn', 'Phase', 'Solty', 'Cn1','Cn2']:
            item = dict[k]
            f = self.formatD[k]
            if f[2]:
                self.writeTupleList(fptr, item, f[0], f[1], f[2])
            else:
                self.writeList(fptr, item, f[0], f[1])
        #next write bnds angs and dihe
        allHBnds = zip(dict['BondHAt1'], dict['BondHAt2'], 
                dict['BondHNum'])
        self.writeTupleList(fptr, allHBnds, "%6d", 12, 3)
        allBnds = zip(dict['BondAt1'], dict['BondAt2'], 
                dict['BondNum'])
        self.writeTupleList(fptr, allBnds, "%6d", 12, 3)

        allHAngs = zip(dict['AngleHAt1'], dict['AngleHAt2'],
                dict['AngleHAt3'], dict['AngleHNum'])
        self.writeTupleList(fptr, allHAngs, "%6d", 12,4)
        allAngs = zip(dict['AngleAt1'], dict['AngleAt2'],
                dict['AngleAt3'], dict['AngleNum'])
        self.writeTupleList(fptr, allAngs, "%6d", 12, 4)

        allHDiHe = zip(dict['DihHAt1'], dict['DihHAt2'],
                dict['DihHAt3'], dict['DihHAt4'], dict['DihHNum'])
        self.writeTupleList(fptr, allHDiHe, "%6d", 12,5)
        allDiHe = zip(dict['DihAt1'], dict['DihAt2'],
                dict['DihAt3'], dict['DihAt4'], dict['DihNum'])
        self.writeTupleList(fptr, allDiHe, "%6d", 12, 5)
        self.writeList(fptr, dict['ExclAt'], '%6d', 12)
        fptr.write('\n')
        fptr.write('\n')
        fptr.write('\n')
        for k in ['AtomSym', 'AtomTree']:
            item = dict[k][:-81]
            self.writeString(fptr, item)
        zList = []
        for i in range(dict['Natom']):
            zList.append(0)
        self.writeList(fptr, zList, "%6d", 12)
        self.writeList(fptr, zList, "%6d", 12)
        fptr.close()


    def writeString(self, fptr, item):
        n = int(ceil(len(item)/80.))
        for p in range(n):
            if p!=n-1:
                fptr.write(item[p*80:(p+1)*80]+'\n')
            else:
                #write to the end, whereever it is
                fptr.write(item[p*80:]+'\n')


    def writeList(self, fptr, outList, formatStr="%4.4s", lineCt=12):
        ct = 0
        s = ""
        nlformatStr = formatStr+'\n'
        lenList = len(outList)
        for i in range(lenList):
            #do something with outList[i]
            s = s + formatStr%outList[i]
            #ct is how many item are in s
            ct = ct + 1
            #if line is full, write it and reset s and ct
            if ct%lineCt==0:
                s = s +  '\n'
                fptr.write(s)
                s = ""
                ct = 0
            #if last entry write it and exit
            elif i == lenList-1:
                s = s +  '\n'
                fptr.write(s)
                break


    def writeTupleList(self, fptr, outList, formatStr="%4.4s", lineCt=12, ll=2):
        ct = 0
        s = ""
        nlformatStr = formatStr+'\n'
        for i in range(len(outList)):
            if i==len(outList)-1:
                for k in range(ll):
                    s = s + formatStr%outList[i][k]
                    ct = ct + 1
                    if ct%lineCt==0:
                        s = s + '\n'
                        fptr.write(s)
                        s = ""
                        ct = 0
                #after adding last entry, if anything left, print it
                if ct!=0:
                    s = s + '\n'
                    fptr.write(s)
            else:
                for k in range(ll):
                    s = s + formatStr%outList[i][k]
                    ct = ct + 1
                    if ct%lineCt==0:
                        s = s + '\n'
                        fptr.write(s)
                        s = ""
                        ct = 0

    def writeItitl(self, fptr, ititl):
        fptr.write(ititl)


    def writeSummary(self, fptr):
        #SUMMARY
        #fptr.write('SUMMARY\n')
        ##FIX THESE NAMES!!!
        kL1 = ['Natom','Ntypes','Nbonh','Mbona',\
            'Ntheth','Mtheta','Nphih','Mphia','Nhparm',\
            'Nparm','Nnb','Nres']
        kL2 = ['Nbona','Ntheta','Nphia','Numbnd',\
            'Numang','Nptra','Natyp','Nphb','Ifpert',\
            'Nbper','Ngper','Ndper']
        kL3 = ['Mbper','Mgper','Mdper','IfBox','Nmxrs',\
            'IfCap']

        for l in [kL1, kL2, kL3]:
            newL = []
            for item in l:
                newL.append(self.prmDict[item])
            #print 'newL=', newL
            self.writeList(fptr, newL, "%6d", 12)


if __name__ == '__main__':
    # load a protein and build bonds
    from MolKit import Read
    p = Read('sff/testdir/p1H.pdb')
    p[0].buildBondsByDistance()

    # build an Amber parameter description objects
    from MolKit.amberPrmTop import ParameterDict
    pd = ParameterDict()

    from MolKit.amberPrmTop import Parm
    prm = Parm()
    prm.processAtoms(p.chains.residues.atoms)
    
