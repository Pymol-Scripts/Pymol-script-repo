#############################################################################
#
# Author: Sophie COON, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/MolKit/getsecondarystructure.py,v 1.18.8.2 2009/06/11 20:10:52 sargis Exp $
#
# $Id: getsecondarystructure.py,v 1.18.8.2 2009/06/11 20:10:52 sargis Exp $
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
    """ Base class to get the datas of secondary structure of a molecule.."""


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
                    #import pdb; pdb.set_trace()
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
           

    def getCoils(self, chain):
        """ gets the coils of a proteic chain by finding which AA are not yet
        in a secondary structure element"""
        residuesCoil = []
        j = 1
        if chain.isDna():
            coil = Coil(chain=chain, index='Nucleic', 
                                start=chain.residues[0], end=chain.residues[-1])
            chain.secondarystructureset.append(coil)
            coil.parent = chain
        if not hasattr(chain, 'AARes') and not chain.isProteic():
            protRes = chain.AARes
            if not len(protRes): 
                res = chain.residues
            else:
                res = protRes
        else:
            res = chain.residues
        chainStructureSet = chain.secondarystructureset
        for r in res:
            # test if the residues we want in the coil has a CA and O.
##              atmCA = r.atoms.get(lambda x:
##                                  string.split(x.name, '@')[0] == 'CA')
##              atmO = r.atoms.get(lambda x:
##                                  string.split(x.name, '@')[0] == 'O')
##              if not atmO:
##                  atmO = r.atoms.get('OXT')
##              if not atmO  or not atmCA:
##             if  r.hasCA | r.hasO != 3:
##                 # if the residue doesn't have a CA nor a O
##                 # it is not condidered !
##                 if residuesCoil != []:
##                     coil = self.buildCoil(residuesCoil, j, chain)
##                     chainStructureSet.append(coil)
##                     coil.parent = chain
##                     residuesCoil = []
##                     j = j+1
##                 break
            if not hasattr(r, 'secondarystructure'):
                residuesCoil.append(r)
                #test if this residue is the last one if yes it build the
                #information on the COIL, if not i is incremented.
                if r == res[-1] and residuesCoil != []:
                    coil = self.buildCoil(residuesCoil, j, chain)
                    chainStructureSet.append(coil)
                    coil.parent = chain
                    residuesCoil = []
                    j = j+1
            else:
                if residuesCoil != []:
                    coil = self.buildCoil(residuesCoil, j, chain)
                    chainStructureSet.append(coil)
                    coil.parent = chain
                    residuesCoil = []
                    j = j+1
        chainStructureSet.sort(self.Compare)

    def buildCoil(self, residuesCoil, j, chain):
        start = residuesCoil[0]
        end = residuesCoil[-1]
        coil = Coil(chain=chain, index=j, start=start, end=end)
        
        return coil
    
    def createSSNodesForChain(self, chain):
        """Create the secondarystructureset of a chain. It gets the proper
        informations in processing the datas returned by the functions
        parseSSData. ."""
        # Add a test to check if all the residues of the chain only have a CA.
##         if chain.findType(Atom).name == len(chain.residues)*['CA']:
##             print " Can't get the Secondary Structure because only the CA of the residues of %s in %s are described in the pdb File."%(chain.id,chain.parent.name)
##             return None
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
        

class GetSecondaryStructureFromStride(GetSecondaryStructure):
    """class to extend GetSecondaryStructure with specific methods to
    get informations on the secondary structure from the output of STRIDE
    """
    def __init__(self, mol):
        import stride
        from MolKit.pdbParser import PdbParser
        s = stride.STRIDE()
        if not isinstance(mol.parser, PdbParser):
            print "cannot use stride to get the secondary structure for the %s"%mol.name
            return None
        
        if not 'ATOM' in mol.parser.keys:
            print "cannot use stride to get the secondary structure for the %s, the file doens't have any ATOM record"%mol.name
            return None
        
        if not s.getPDBRecords( mol.parser.allLines, len(mol.parser.allLines) ):
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
