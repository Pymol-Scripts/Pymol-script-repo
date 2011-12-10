#############################################################################
#
# Author: Sophie COON, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/MolKit/getsecondarystructure.py,v 1.28.2.2 2011/06/22 19:54:14 sargis Exp $
#
# $Id: getsecondarystructure.py,v 1.28.2.2 2011/06/22 19:54:14 sargis Exp $
#

"""
This module implements the classes GetSecondaryStructure,
GetSecondaryStructureFromFile, GetSecondaryStructureWeb,
These classes are expected to be specialized in getting  the informations on
the secondary Structure of a protein and create a secondarystructureset for a
chain
"""

from MolKit.protein import Helix, Strand, Turn, Coil,\
     SecondaryStructureSet, Residue, ResidueSet
from MolKit.molecule import Atom
import string, types


class GetSecondaryStructure:
    """
    Base class to get the datas of secondary structure of a molecule..
    This creates a dictionary called self.ssDataForMol with key
    chain.id and value is a list of [helStartEndForChain,
    strandStartEndForChain, turnStartEndForChain, None]
"""


    def Compare(self, struct1, struct2):
        """ compares 2 structures looking at the start residues index."""
        assert(struct1.chain.id == struct2.chain.id)
        res = struct1.chain.residues
        x = cmp(res.index(struct1.start), res.index(struct2.start))
        return x

          
    def createSSstructures(self, chain, heldatas, strdatas, turndatas,
                           coildatas):
        # create a secondary structure set
        #chainStructureSet = SecondaryStructureSet()
        chain.secondarystructureset = SecondaryStructureSet()
        # helices objects
        self.buildSecondaryStructureObjects(heldatas,chain)
        # sheets objects
        self.buildSecondaryStructureObjects(strdatas,chain)
        # turns objects
        if turndatas :
            self.buildSecondaryStructureObjects(turndatas,chain)

        # coils
        self.getCoils(chain)
        if hasattr(chain,'residueInSS'):
            chain.residuesInSS = chain.residues.get(lambda x: hasattr(x, 'secondarystructure'))
##         # coils
##         if coildatas:
##             self.buildSecondaryStructureObjects(coildatas,chain)
##             chain.secondarystructureset.sort(self.Compare)

##         else:
##             self.getCoils(chain)
##             # When FromFile the set of residues belonging to a ss is:
##             chain.residuesInSS = chain.residues.get(lambda x: hasattr(x, 'secondarystructure'))
    

    

    def buildSecondaryStructureObjects(self, datas, chain):
        ssClass = datas[0]
        chainStructureSet = chain.secondarystructureset
        for i in range(1,len(datas)):
            if type(datas[i]) is types.DictionaryType:
                start = datas[i]['start']
                end = datas[i]['end']
                kw = datas[i]
            elif isinstance(datas[i], ResidueSet):
                start = datas[i][0]
                end = datas[i][1]
                kw = {}
                kw['start']=start
                kw['end']=end
            elif isinstance(datas[i], list):
                start = datas[i][0]
                end = datas[i][1]
                kw = {}
                kw['start']=start
                kw['end']=end
            # check that the parser has found the correct start and end
            if not start or not end:
                return None   # error: no secondary structure set created

            if chainStructureSet:
                # Check if the structure has not been already built
                if start in chainStructureSet.start\
                   and end in chainStructureSet.end:
                    continue
                # Check if the entire new structure doesn't belong to a
                # previous one..
                # Check if the structure has already beein built.
                # This is for strands belonging to two different sheets.

                if start in chainStructureSet.residues\
                   and end in chainStructureSet.residues:
                    print "In %s, %s and %s already belongs to %s"\
                          %(start.parent.id, start.name, end.name,
                            start.secondarystructure.name)
                            
                    continue

                if start in chainStructureSet.residues:
                    print "In %s, %s is  already the end of %s"\
                          %(start.parent.id, start.name,
                            start.secondarystructure.name)
                    #continue
                if end in chainStructureSet.residues:
                    print "In %s, %s is  already the start of %s"\
                          %(end.parent.id, end.name,
                            end.secondarystructure.name)
                    continue
            #kw = datas[i]
            kw['index'] = i
            kw['chain'] = chain
            secondarystructure = apply(ssClass, (), kw)
            chainStructureSet.append(secondarystructure)
            secondarystructure.parent = chain

    def close(self, ca1, ca2):
        x1,y1,z1 = ca1.coords
        x2,y2,z2 = ca2.coords
        d2 = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1)
        #print d2
        return d2 < 16. # 4 Angs **2 cutoff


    def getCoils(self, chain):
        """
        gets the coils of a proteic chain by finding which AA are not yet
        in a secondary structure element
        """
        chain.ribbonType()
        if chain.ribbonType()=='NA':
            coil = Coil(chain=chain, index='Nucleic', 
                        start=chain.DNARes[0], end=chain.DNARes[-1])
            chain.secondarystructureset.append(coil)
            coil.parent = chain
            coil.isNA = True

        else:
            res = chain.AARes
            chainStructureSet = chain.secondarystructureset

            j = 1 # coil counter
            
            # loop over residues and build coils for continuous segments of
            # unassigned residues 
            endedWithGap = False
            #print 'FOR CHAIN ', chain.id

            residuesCoil = [] # list of residues in current Coil
            # look over all residues to build list of residues in current Coil

            for r in res:
                #print r
                # if we hit a residue
                #cCAatom = None
                if not r.hasCA: # we hit residue that breaks coil
                    if len(residuesCoil):
                        coil = self.buildCoil(residuesCoil, j, chain)
                        coil.gapAfter = True
                        if endedWithGap:
                            coil.gapBefore = True
                        chainStructureSet.append(coil)
                        coil.parent = chain
                        r1 = residuesCoil[0]
                        r2 = residuesCoil[-1]
                        #print '0', coil, r1, r2,coil.gapBefore, coil.gapAfter
                        endedWithGap = True
                        residuesCoil = [r]
                        j = j+1
                        
                else:
                    anames = r.atoms.name
                    # CA atom of current res
                    # using get('CA*')[0] handle alternate location as in
                    # THR 4 in 1KZK.pdb
                    cCAatom = r.atoms.get('CA*')[0]
                    #cNatom = r.atoms[anames.index('N')] # N atom of current res
                    #cCatom = r.atoms[anames.index('C')] # C atom of current res

                    if hasattr(r, 'secondarystructure'):
                        if len(residuesCoil):
                            coil = self.buildCoil(residuesCoil, j, chain)
                            if endedWithGap:
                                coil.gapBefore = True
                            chainStructureSet.append(coil)
                            coil.parent = chain
                            r1 = residuesCoil[0]
                            r2 = residuesCoil[-1]
                            #print '1', coil, r1, r2, coil.gapBefore, coil.gapAfter
                            endedWithGap = False
                            residuesCoil = []
                            j = j+1

                    else:
                        if len(residuesCoil)==0:
                            residuesCoil = [r]

                        else:
                            if self.close(pCAatom, cCAatom):
                                residuesCoil.append(r)
                                if r==res[-1]: #last seg is coil
                                    coil = self.buildCoil(residuesCoil, j, chain)
                                    chainStructureSet.append(coil)
                                    if endedWithGap:
                                        coil.gapBefore = True
                                    coil.parent = chain
                                   
                            elif len(residuesCoil): # distance gap
                                coil = self.buildCoil(residuesCoil, j, chain)
                                chainStructureSet.append(coil)
                                coil.gapAfter = True
                                if endedWithGap:
                                    coil.gapBefore = True
                                coil.parent = chain
                                r1 = residuesCoil[0]
                                r2 = residuesCoil[-1]
                                #print '2', coil, r1, r2,coil.gapBefore, coil.gapAfter
                                endedWithGap = True
                                residuesCoil = [r]
                                j = j+1
                
                    pCAatom = cCAatom
                #pNatom = cNatom
                #pCatom = cCatom
                        
            ## res = chain.AARes

            ## chainStructureSet = chain.secondarystructureset
            ## resNum = None
            ## try:
            ##     resNum = int(res[0].number)  
            ## except:
            ##     pass

            ## j = 1 # coil counter

            ## for rn, r in enumerate(res):
            ##     print 'AAA1', r, resNum, r.number
            ##     if resNum:
            ##         if resNum != int(r.number):
            ##             if residuesCoil:
            ##                 print 'AAA2', j, len(residuesCoil)
            ##                 coil = self.buildCoil(residuesCoil, j, chain)
            ##                 chainStructureSet.append(coil)
            ##                 coil.parent = chain
            ##                 coil.gap =int(r.number) - resNum #used in secondaryStructureCommands to fix Bug 1120
            ##                 residuesCoil = [r]
            ##                 r.gap = True
            ##                 j = j+1
            ##         try:
            ##             if res[rn+1].number != r.number: #icode is not '':
            ##                 resNum = int(r.number) +1 
            ##         except IndexError:
            ##             pass

            ##     if not hasattr(r, 'secondarystructure') and r.hasCA and r.hasO:
            ##         print 'AAA3', r, len(residuesCoil), r == res[-1]
            ##         residuesCoil.append(r)
            ##         #test if this residue is the last one if yes it build the
            ##         #information on the COIL, if not i is incremented.

            ##         if r == res[-1] and residuesCoil != []:
            ##             coil = self.buildCoil(residuesCoil, j, chain)
            ##             chainStructureSet.append(coil)
            ##             coil.parent = chain
            ##             if hasattr(residuesCoil[0],'gap'):
            ##                 coil.gap = True 
            ##             residuesCoil = []
            ##             j = j+1
            ##     else:
            ##         print 'AAA4', len(residuesCoil)
            ##         if residuesCoil != []:
            ##             coil = self.buildCoil(residuesCoil, j, chain)
            ##             chainStructureSet.append(coil)
            ##             coil.parent = chain
            ##             if hasattr(residuesCoil[0],'gap'):
            ##                 coil.gap = True 
            ##             residuesCoil = []
            ##             j = j+1

        chain.secondarystructureset.sort(self.Compare)


    def buildCoil(self, residuesCoil, j, chain):
        start = residuesCoil[0]
        end = residuesCoil[-1]
        coil = Coil(chain=chain, index=j, start=start, end=end)
        coil.isNA = False
        return coil

  
    def createSSNodesForChain(self, chain):
        """Create the secondarystructureset of a chain. It gets the proper
        informations in processing the datas returned by the functions
        parseSSData. ."""
        # Add a test to check if all the residues of the chain only have a CA.
##         if chain.findType(Atom).name == len(chain.residues)*['CA']:
##             print " Can't get the Secondary Structure because only the CA of the residues of %s in %s are described in the pdb File."%(chain.id,chain.parent.name)
##             return None
        #import traceback
        #traceback.print_stack()
        if not self.ssDataForMol.has_key(chain.id): return
        helData, strandData, turnData, coilData= self.ssDataForMol[chain.id]
        self.createSSstructures(chain, helData, strandData,
                                turnData, coilData)


    
class GetSecondaryStructureFromFile(GetSecondaryStructure):
    """class to extend GetSecondaryStructure with specific methods to
    get informations on the secondary structure from files (PDB file,
    MOL2) and build the equivalent nodes for the molecule.. """
    def __init__(self, mol):
        """ Get the information describing the secondary structure for
        the molecule from the molecule."""
        self.ssDataForMol = mol.parser.parseSSData(mol)
        


class GetSecondaryStructureFromPross(GetSecondaryStructure):
    """class to extend GetSecondaryStructure with specific methods to
    get informations on the secondary structure from PROSS This module 
    uses much of the code from the original BIOMOL collection of utilities 
    written by Raj Srinivasani with enhancements by Nick Fitzkee.

    The script was put together by Pat Fleming so that a user would not need
    to have the BIOMOL distribution installed to run PROSS.

    Note that since Raj's time the definitions of mesostates has been superceded
    by the fine grained 30 deg x 30 deg grid for most purposes. Either mesostate
    grid will work for PROSS. Give your choice as an argument (see USAGE below).

    Date: September 2004
    Author: Pat Fleming, pat.fleming@jhu.edu
    http://www.roselab.jhu.edu/utils/pross.html"""
    def __init__(self, mol, default='fgmeso'):
        """ # Setup Fine Grain Mesostate Bins (ala Pat Fleming)"""
        """ see PROSS.py """
        from PROSS import MSDEFS
        self.mode=default
        self.MSDEFS=MSDEFS[default]
        self.mol = mol
        #self.updatePhiPsi() # this called by rc_ss
        self.ssDataForMol = {}
        self.sst={}
        for chain in self.mol.chains:
        #    self.ssDataForMol[chain.id] = []
            self.ssDataForMol[chain.id] = self.rc_ss(chain)
    #print 'SSABC', self.ssDataForMol
        
    def updatePhiPsi(self, chain):
        aa = [res for res in chain.residues if res.hasCA]
        map(lambda x: x.getPhi(), aa)
        map(lambda x: x.getPsi(), aa)

        
    def res_rc(self,r1, r2, r3=180):
        """res_rc(r1, r2, r3) - get mesostate code for a residue

        Given a phi (r1), psi (r2), and omega (r3) torsion angle, calculate
        the mesostate code that describes that residue. 
        A mesostate will be returned in
        all but two cases:  if omega deviates from planarity by more than 90
        degrees, '*' is returned.  Also, if any torsions are greater than
        180.0 (biomol assignes 999.0 degrees to angles that that are
        indeterminate), prosst.INVALID is returned.  Here, r3 defaults to
        180.0.
        """

        ms = self.MSDEFS
        OMEGA   = ms['OMEGA']
        INVALID = ms['INVALID']
        PHI_OFF = ms['PHI_OFF']
        PSI_OFF = ms['PSI_OFF']
        DELTA   = ms['DELTA']
        RC_DICT = ms['RC_DICT']
        if r1 == None or r2 == None : return INVALID 
        if (abs(r3) <= 90.0):
            return OMEGA
        elif r1>180.0 or r2>180.0 or r3>180.0:
            return INVALID

        ir1 = -int(PHI_OFF) + int(round((r1+PHI_OFF)/DELTA )) * int(DELTA)
        ir2 = -int(PSI_OFF) + int(round((r2+PSI_OFF)/DELTA )) * int(DELTA)

        while ir1 <= -180: ir1 = ir1 + 360
        while ir1 >   180: ir1 = ir1 - 360
        while ir2 <= -180: ir2 = ir2 + 360
        while ir2 >   180: ir2 = ir2 - 360

        return RC_DICT[(ir1,ir2)]

    def rc_codes(self,chain, phi=None, psi=None, ome=None):
        """rc_codes(chain, phi, psi, ome) - return rotamer codes

        Given a protein chain (and optionally phi, psi, omega), this
        function will return a list of mesostate codes that
        applies to the chain, as determined by res_rc.
        """
        n = range(len(chain.residues))
        if phi is None: phi = chain.residues.phi
        if psi is None: psi = chain.residues.psi
        #if ome is None: ome = map(chain.omega, n)
        return map(lambda x, y: self.res_rc(x, y), phi, psi)


    def rc_ss(self, chain, phi=None, psi=None, ome=None):
        """rc_ss(chain, phi, psi, ome) - calculate secondary structure

        This function calculates the secondary structure using the PROSS method
        with rotamer codes.  Given a chain, and optionally a list of phi,
        psi, and omega, it calculates the backbone secondary structure of
        the chain.  The return value is (phi, psi, ome, sst), where
        phi, psi, and ome are calculated if not specified, and sst is the
        secondary structure codes: H = helix, E = strand, P = PII, C = coil.
        """
        
        ms = self.MSDEFS
        PII    = ms['PII']
        TURNS  = ms['TURNS']
        HELIX  = ms['HELIX']
        STRAND = ms['STRAND']

        nres = len(chain.residues)
        self.updatePhiPsi(chain)
        if phi is None:
            #chain.gaps()
            phi = chain.residues.phi
        if psi is None: psi = chain.residues.psi
        #if ome is None: ome = None #map(chain.residues.omega, xrange(nres))
        #print psi
        codes = self.rc_codes(chain, phi, psi)

        #chain.gaps()
        
        sst = ['C']*nres
        hData = [Helix]
        sData = [Strand]
        tData = [Turn]
        cData = [Coil]

        is_PII = PII.has_key

        for i in xrange(nres-1):
            code = codes[i]
            if is_PII(code):
                sst[i] = 'P'

        is_turn = TURNS.has_key

        for i in xrange(nres-1):
            code = codes[i] + codes[i+1]
            if is_turn(code):        
                sst[i] = sst[i+1] = 'T'

        helices = self._rc_find(codes, HELIX)
        strands = self._rc_find(codes, STRAND)
        #turn = self._rc_find(codes, TURNS)
        #coil = self._rc_find(codes, STRAND)

        for helix in helices:
            i, j = helix
            #hData.append([chain.residues[i],chain.residues[j]])
            for k in range(i, j):
                sst[k] = 'H'

        for strand in strands:
            i, j = strand
            #sData.append([chain.residues[i],chain.residues[j]])
            for k in range(i, j):
                if sst[k] in ('C', 'P'): sst[k] = 'E'
        #self.sst = sst
        
        turns = self.findSS(sst,('T'))
        #print 'FOGO TURNS', turns
        for turn in turns:
            i, j = turn
            tData.append([chain.residues[i],chain.residues[j]])
        strands = self.findSS(sst,('E'))
        for strand in strands:
            i, j = strand
            sData.append([chain.residues[i],chain.residues[j]])
        helices = self.findSS(sst,('H'))
        for helix in helices:
            i, j = helix
            hData.append([chain.residues[i],chain.residues[j]])
        cData = None#self.findSS(self.sst,('C','P'))
        self.sst[chain.id]=sst            
        return (hData, sData, tData, cData)


    def findSS(self,sst,ssType):
        start=0
        end=0
        prev=None
        next=None
        res=[]            
        for i,ss in enumerate(sst):    
            if ss == ssType :
                if ss != prev :
                    start=i
                else : 
                    if i == len(sst)-1 or ss != sst[i+1] :
                        end = i
                        res.append([start,end])
            prev = ss
        return res


    def _rc_find(self,codes, pattern):
        """_rc_find(codes, pat_obj) - find a endpoints of a regexp

        Given a list of mesostate codes, this function identifies a endpoints
        of a match  <pattern>.  <pat_obj> is a compiled regular expression
        pattern whose matches will be returned as pairs indicated start,
        end in <codes>
        """

        CODE_LENGTH = self.MSDEFS['CODE_LENGTH']

        if not type(codes) == type(''):
            codes = string.join(codes, '')

        matches = []
        it = pattern.finditer(codes)

        try:
            while 1:
                mat = it.next()
                matches.append((mat.start()/CODE_LENGTH, mat.end()/CODE_LENGTH))
        except StopIteration:
            pass

        return matches



class GetSecondaryStructureFromStride(GetSecondaryStructure):
    """class to extend GetSecondaryStructure with specific methods to
    get informations on the secondary structure from the output of STRIDE
    """
    def __init__(self, mol):
        import traceback
        traceback.print_stack()
        
        import stride
        from MolKit.pdbParser import PdbParser
        s = stride.STRIDE()
        if not isinstance(mol.parser, PdbParser):
            print "cannot use stride to get the secondary structure for the %s"%mol.name
            return None
        
        if not 'ATOM' in mol.parser.keys:
            print "cannot use stride to get the secondary structure for the %s, the file doens't have any ATOM record"%mol.name
            return None

        if not s.getPDBRecords( mol.parser.allLines, len(mol.parser.allLines)):
            print "STRIDE has failed to parse ", mol
            self.ssDataForMol = {}
            return
        #s.run( report = 0 )
        s.run( )

        self.ssDataForMol = self.parseSSData( s, mol )


    def getAsn(self, chain):
        asn = []
        chain.residuesInSS = ResidueSet()
        from stride.stride import RESIDUE
        chain.residuesInSS.elementType = RESIDUE
        for j in range(chain.NRes):
            res = chain.getResidue(j)
            chain.residuesInSS.append(res)
            asn.append( (res.ResType+res.PDB_ResNumb, res.Prop.Asn) )
        return asn

    def parseSSData(self, s, mol):
        ssDataForMol = {}
        for j in xrange(s.NChain):
            k = 0
            # k is the constant representing the difference in indices between
            # the line index in the stride report and the residue number.
            # it is equals to 1 because the line number is 0 based and the
            # residue number is 1 based.
            chain = s.getChain(j)
            if chain is None: continue
            asnForChain = self.getAsn(chain)
            if hasattr(chain, 'Id'):
                cId = chain.Id
            else:
                cId = chain.__getmethods__['Id'](chain)
            if mol.chains.id == ['UNK']: #this is needed for pqr files
                c = mol.chains[0]
            else:
                c = mol.chains.get(lambda x: x.id == cId)[0]
            #c = mol.chains[j]
            hData = [Helix]
            sData = [Strand]
            tData = [Turn]
            cData = [Coil]
            ssType = None
            lenResidues = len(c.residues)
            for i in xrange(len(asnForChain)):
                l = asnForChain[i]
                if l[1] == 'B': continue # Don't know what B stands for...
                elif (i+k) >=lenResidues: continue # Bugfix for  #1033
                elif l[0] != c.residues[i+k].name:
                    # check if there is a residue in the chain but not in the
                    # stride output file.
                    k = k+1
                    discontinue = 1
                else:
                    discontinue = 0
                        
                if ssType != l[1] or (ssType == l[1] and discontinue==1):
                    if ssType is None:
                        start = c.get(l[0])
                        if start is None:
                            raise
                        ssType = l[1]
                    else:
                        end = c.get(asnForChain[i-1][0])
                        if end is None:
                            raise
                        startend = start+end
                        hData, sData, tData, cData  = \
                               self.buildSSData(ssType,startend,hData, sData,
                                                tData, cData)
                        start = c.get(l[0])

                        if start is None:
                            raise
                        ssType = l[1]
                        
                        if i == len(asnForChain)-1:
                            end = c.get(l[0])
                            if end is None :
                                raise
                            startend = start+end
                            hData, sData, tData, cData  = \
                                   self.buildSSData(ssType,startend,
                                                    hData, sData,
                                                    tData, cData)
            ssDataForMol[c.id] = [hData, sData, tData, cData]
        return ssDataForMol            

    def buildSSData(self, ssType,startend, hData, sData, tData, cData):
        if ssType in ['C']:
            cData.append({'start':startend[0], 'end':startend[1]})
            #cData.append( startend )
        elif ssType in ['H', 'G']:
            hData.append({'start':startend[0], 'end':startend[1], 'helClass':ssType})
            #hData.append( startend )
        elif ssType in ['T']:
            tData.append({'start':startend[0], 'end':startend[1]})
            #tData.append( startend )
        elif ssType in ['E']:
            sData.append({'start':startend[0], 'end':startend[1]})
            #sData.append( startend )
        return hData, sData, tData, cData
