## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/MolKit/protein.py,v 1.61.2.2 2011/06/24 22:16:45 sanner Exp $
#
# $Id: protein.py,v 1.61.2.2 2011/06/24 22:16:45 sanner Exp $
#

"""
This Module implements the classes Residue, Chain and Protein.

Residue and Chain are TreeNode objects similarly to the Atom object.

Protein is a specialization of the Molecule class to represent a protein
in a 4 level tree (from leafs to root: atoms, residues, chains, molecule)
"""

from MolKit.tree import TreeNode, TreeNodeSet, TreeNodeSetSelector
from MolKit.molecule import Molecule, MoleculeSet, Atom, AtomSet, Bond
import re
from types import NoneType, StringType, IntType, TupleType, FloatType, ListType
from string import split, upper, find
from numpy.oldnumeric import sum
global bhtreeFlag
try:
    import bhtree
    bhtreeFlag = 1
except:
    bhtreeFlag = 0

# FIXME had to overwrite isBelow because of the Protein vs Molecule pb
class ProteinMolecule(Molecule):
##     def __del__(self):
##         self.__dict__.clear()
##         Molecule.__del__(self)

    def findType(self, _type, uniq=0):
        if self.setClass is MoleculeSet and _type is Protein:
            n = ProteinSet([self], stringRepr=self.name)
        elif self.setClass is ProteinSet and _type is Molecule:
            n = MoleculeSet([self], stringRepr=self.name)
        else:
            n = self.setClass([self], stringRepr=self.name)
        if n.elementType == _type:
            return n
        result = n.findType(_type, uniq=uniq)
        return result
        
    def isBelow(self, Klass):
        if Klass is Molecule: Klass=Protein
        return Molecule.isBelow(self, Klass)

#######################################################################
##  RESIDUES                                                         ##
#######################################################################

class ResidueSet(TreeNodeSet):
    """Class to extend a TreeNodeSet with residue specific methods"""

    def __init__(self, objects=None, stringRepr=None, comments="", keywords=[]):
        TreeNodeSet.__init__(self, objects, Residue, stringRepr, comments, keywords)
        if stringRepr is None:
            ## build string repr
            # find all molecules
            strr = ''
            if objects is not None and len(objects):
                mols = [x.top for x in objects]
                molDict = {}.fromkeys(mols, [])

                # find residues in each molecule
                for a in objects:
                    molDict[a.top].append(a)

                for k,v in molDict.items():
                    ## this line was causing an endless loop when we select
                    ## Chain A in protease in Dashboard and then clcik on
                    ## display line column for chain A
                    ## I (MS) replaced len(k.chains.residues) by the for loop
                    ## which avoids calling TreeNodeSet.__getattr__ which would
                    ## try to create a ResidueSet, thus creating the loop
                    #if len(v)==len(k.chains.residues): # special case all res
                    nbres = 0
                    for c in k.chains:
                        nbres += len(c.residues)
                    if len(v)==nbres: # special case all res
                        strr += k.name+'::;'
                    else:
                        strr += k.name+"::"
                        for a in v:
                            strr += a.name+','
                        strr += ';'
            
            self.stringRepr = strr
        

#    def get(self, selectionString, selector=None, sets=None,
#                    caseSensitive=True, escapeCharacters=False):
#        if selector is None: 
#            selector = ResidueSetSelector()
#            selector.caseSensitive = caseSensitive
#            selector.escapeCharacters = escapeCharacters
#        #results, msg = selector.select(self, selectionString, sets, 
#        selectionStringRepr = str(selectionString)
#        results = selector.processListItem(self, selectionString, sets) 
#        selectionStringRepr = '&' + selectionStringRepr
#        results.setStringRepr(selectionStringRepr)
#        return results


    def getSelector(self):
        if self.selector==None:
            self.selector = ResidueSetSelector()
        return self.selector

##     def get(self, selectionString, selector=None, sets=None,
##                     caseSensitive=True, escapeCharacters=False,
##                     returnMsg=False):
##         """
##         selects from self.data, objects of self.elementType specified by selectionString
##         """
##         #print "RS: get: selectionString=", selectionString

##         if selector is None: 
##             selector = ResidueSetSelector()
##             selector.caseSensitive = caseSensitive
##             selector.escapeCharacters = escapeCharacters
##         selectionStringRepr = '&' + str(selectionString)
##         if type(selectionString)==StringType:
##             result, msg = selector.select(self, selectionString)
##             result = self.ReturnType(result)
##             result.setStringRepr(selectionStringRepr)
##             if returnMsg: result = (result, msg)
##             return result
##             #if selector.select(self, selectionString)[0] is None:
##             #    return self.ReturnType([])
##             #else:
##             #    return selector.select(self, selectionString)[0]
##         elif callable(selectionString):
##             result = filter(selectionString, self.data)
##             if len(result)==len(self.data):
##                 return self
##             else:
##                 result = self.ReturnType(result)
##                 result.setStringRepr(selectionStringRepr)
##                 return result
##         else:
##             raise RuntimeError("argument has to be a function or a string")

        
from mglutil.math.torsion import torsion

class Residue(ProteinMolecule):
    """Class to represent an amino acide. Inherits from tree element"""

    _bbname = ['N', 'CA', 'C', 'O']
    _bbtype = ['N', 'C', 'C', 'O']



##     def __del__(self):
##         self.__dict__.clear()
##         ProteinMolecule.__del__(self)


    def __init__(self, type='UNK', number=-1, icode='', parent=None,
                 elementType=Atom, list=None, childrenName='atoms',
                 setClass=ResidueSet,
                 childrenSetClass=AtomSet, top=None, childIndex=None,
                 assignUniqIndex=1):
        """Residue constructor.
        Arguments:
        type (string)
        number (integer or string)
        icode (1character) insertion code
        optional parent (instance of a TreeNode)
        optional elementType (instance of class inheriting from TreeNode)"""
        name = type+str(number)+icode
        Molecule.__init__(self, name, parent, elementType,
                          list, childrenName, setClass, childrenSetClass, top,
                          childIndex, assignUniqIndex)
        self.type = type
        self.number = number
        self.icode = icode
        self.psi  = None # not calculated
        self.phi  = None # not calculated
        
    def getPsi(self):
        """  compute PSI N(i),CA(i),C(i),N(i+1) """
        nextResidue = self.getNext()
        if self.getNext() is not None and nextResidue.hasCA:
            try:
                names = [name.split("@")[0] for name in self.atoms.name]
                idx=names.index('N') ; at1 = self.atoms[idx]
                idx=names.index('CA'); at2 = self.atoms[idx]
                idx=names.index('C') ; at3 = self.atoms[idx]
                nextResidueAtoms = nextResidue.atoms
                
                names = [name.split("@")[0] for name in nextResidueAtoms.name]
                idx= names.index('N')
                at4 = nextResidueAtoms[idx]            
    
                self.psi = torsion(at1.coords, at2.coords, at3.coords, at4.coords)
            except:
                self.psi = 0
        return self.psi

    def getPhi(self):
        """  compute PHI C(i-1),N(i),CA(i),c(i)  """
        if self.getPrevious() is not None:
            from mglutil.math.torsion import torsion
            prevResidue = self.getPrevious()
            if prevResidue is None or not prevResidue.hasCA:
                self.phi = 0.
            else:
                try:
                    prevResidueAtoms = self.getPrevious().atoms
                    names = [name.split("@")[0] for name in prevResidueAtoms.name]
                    idx= names.index('C')
                    at1 = prevResidueAtoms[idx]            
    
                    names = [name.split("@")[0] for name in self.atoms.name]
                    idx=names.index('N') ; at2 = self.atoms[idx]
                    idx=names.index('CA'); at3 = self.atoms[idx]
                    idx=names.index('C') ; at4 = self.atoms[idx]
    
                    self.phi = torsion(at1.coords, at2.coords,
                                       at3.coords, at4.coords)
                except:
                    self.phi = 0
        return self.phi


    def buildBondsByDistance(self):
        """Build bonds between atoms inside this residue, based on distance
        WARNING this is a n^2 process, should only be used for small
        molecules like residues"""

        if self.hasBonds: return
#        atoms = findType(Atom)
        atoms = self.children
        if not self.hasBonds:
            self.buildBondsByDistanceOnAtoms(atoms)
            self.hasBonds = 1
        return len(atoms)


    def hetatm(self):
        """Return atomset containing hetatm atoms"""
        
        heta = self.atoms.get(lambda x: x.hetatm==1)
        if heta is None: return AtomSet( []  )
        else: return heta


    def backbone(self):
        """Return atomset containing backbone atoms"""
        
        bb = self.get(lambda x, name=self._bbname,type=self._bbtype:
                      x.name in name and x.element in type)
        if bb is None: return AtomSet( []  )
        else: return bb


    def sidechain(self):
        """Return atomset containing sidechain atoms"""

        return self.atoms - self.backbone()

    def getAtmsAndCoords(self, atmNames):
        """
        Function returning the coords of all the atoms of the given atom name
        or None if one is not in the atoms residues
        """
        childNames = self.childByName.keys()
        check = map(lambda x, cn = childNames: x in cn, atmNames)
        # all the given names should be atoms of the residue
        if 0 in check:
            return 0, None
        coords = []
        # short cut to coords.append
        coordsappend = coords.append
        cName = self.childByName
        for name in atmNames:
            atm = cName[name]
            if atm.alternate: coordsappend(atm.getAverageCoords())
            else:
                coordsappend(atm.coords)
        return 1, coords

            

class ResidueSetSelector(TreeNodeSetSelector):

    residueList = {}
    residueList['std']={
        'ala':True, 'ALA':True,
        'arg':True, 'ARG':True,
        'asn':True, 'ASN':True,
        'asp':True, 'ASP':True,
        'cys':True, 'CYS':True,
        'gln':True, 'GLN':True,
        'glu':True, 'GLU':True,
        'gly':True, 'GLY':True,
        'his':True, 'HIS':True,
        'ile':True, 'ILE':True,
        'leu':True, 'LEU':True,
        'lys':True, 'LYS':True,
        'met':True, 'MET':True,
        'phe':True, 'PHE':True,
        'pro':True, 'PRO':True,
        'ser':True, 'SER':True,
        'thr':True, 'THR':True,
        'trp':True, 'TRP':True,
        'tyr':True, 'TYR':True,
        'val':True, 'VAL':True,
        }
    residueList['acidic']={
        'ASP':True, 'asp':True,
        'glu':True, 'GLU':True,
        }
    residueList['acyclic']={
        'ala':True, 'ALA':True,
        'arg':True, 'ARG':True,
        'asn':True, 'ASN':True,
        'asp':True, 'ASP':True,
        'cys':True, 'CYS':True,
        'glu':True, 'GLU':True,
        'gln':True, 'GLN':True,
        'gly':True, 'GLY':True,
        'ile':True, 'ILE':True,
        'leu':True, 'LEU':True,
        'lys':True, 'LYS':True,
        'met':True, 'MET':True,
        'ser':True, 'SER':True,
        'thr':True, 'THR':True,
        'val':True, 'VAL':True,
        }
    residueList['aliphatic']={
        'ala':True, 'ALA':True,
        'gly':True, 'GLY':True,
        'ile':True, 'ILE':True,
        'leu':True, 'LEU':True,
        'val':True, 'VAL':True,
        }
    residueList['aromatic']={
        'his':True, 'HIS':True,
        'phe':True, 'PHE':True,
        'trp':True, 'TRP':True,
        'tyr':True, 'TYR':True,
        }
    residueList['basic']={
        'arg':True, 'ARG':True,
        'his':True, 'HIS':True,
        'lys':True, 'LYS':True,
        }
    residueList['buried']={
        'ala':True, 'ALA':True,
        'cys':True, 'CYS':True,
        'ile':True, 'ILE':True,
        'leu':True, 'LEU':True,
        'met':True, 'MET':True,
        'phe':True, 'PHE':True,
        'trp':True, 'TRP':True,
        'val':True, 'VAL':True,
        }
    residueList['charged']={
        'arg':True, 'ARG':True,
        'asp':True, 'ASP':True,
        'glu':True, 'GLU':True,
        'his':True, 'HIS':True,
        'lys':True, 'LYS':True,
        }
    residueList['cyclic']={
        'his':True, 'HIS':True,
        'phe':True, 'PHE':True,
        'pro':True, 'PRO':True,
        'trp':True, 'TRP':True,
        'tyr':True, 'TYR':True,
        }
    residueList['hydrophobic']={
        'ala':True, 'ALA':True,
        'gly':True, 'GLY':True,
        'ile':True, 'ILE':True,
        'leu':True, 'LEU':True,
        'met':True, 'MET':True,
        'phe':True, 'PHE':True,
        'pro':True, 'PRO':True,
        'trp':True, 'TRP':True,
        'tyr':True, 'TYR':True,
        'val':True, 'VAL':True,
        }
    residueList['large']={
        'arg':True, 'ARG':True,
        'glu':True, 'GLU':True,
        'gln':True, 'GLN':True,
        'his':True, 'HIS':True,
        'ile':True, 'ILE':True,
        'leu':True, 'LEU':True,
        'lys':True, 'LYS':True,
        'met':True, 'MET':True,
        'phe':True, 'PHE':True,
        'trp':True, 'TRP':True,
        'tyr':True, 'TYR':True,
        }
    residueList['medium']={
        'asn':True, 'ASN':True,
        'asp':True, 'ASP':True,
        'cys':True, 'CYS':True,
        'pro':True, 'PRO':True,
        'thr':True, 'THR':True,
        'val':True, 'VAL':True,
        }
    residueList['negative']={
        'asp':True, 'ASP':True,
        'glu':True, 'GLU':True,
        }
    residueList['neutral']={
        'ala':True, 'ALA':True,
        'asn':True, 'ASN':True,
        'cys':True, 'CYS':True,
        'gln':True, 'GLN':True,
        'gly':True, 'GLY':True,
        'his':True, 'HIS':True,
        'ile':True, 'ILE':True,
        'leu':True, 'LEU':True,
        'met':True, 'MET':True,
        'phe':True, 'PHE':True,
        'pro':True, 'PRO':True,
        'ser':True, 'SER':True,
        'thr':True, 'THR':True,
        'trp':True, 'TRP':True,
        'tyr':True, 'TYR':True,
        'val':True, 'VAL':True,
        }
    residueList['polar']={
        'arg':True, 'ARG':True,
        'asn':True, 'ASN':True,
        'asp':True, 'ASP':True,
        'cys':True, 'CYS':True,
        'glu':True, 'GLU':True,
        'gln':True, 'GLN':True,
        'his':True, 'HIS':True,
        'lys':True, 'LYS':True,
        'ser':True, 'SER':True,
        'thr':True, 'THR':True,
        }
    residueList['positive']={
        'arg':True, 'ARG':True,
        'his':True, 'HIS':True,
        'lys':True, 'LYS':True,
        }
    residueList['small']={
        'ala':True, 'ALA':True,
        'gly':True, 'GLY':True,
        'ser':True, 'SER':True,
        }
    residueList['surface']={
        'arg':True, 'ARG':True,
        'asn':True, 'ASN':True,
        'asp':True, 'ASP':True,
        'glu':True, 'GLU':True,
        'gln':True, 'GLN':True,
        'gly':True, 'GLY':True,
        'his':True, 'HIS':True,
        'lys':True, 'LYS':True,
        'pro':True, 'PRO':True,
        'ser':True, 'SER':True,
        'thr':True, 'THR':True,
        'tyr':True, 'TYR':True,
        }
    
    from MolKit.PDBresidueNames import DNAnames, RNAnames, Nucleotides, \
         AAnames, ionNames, waterNames, allResidueNames
    residueList['ions']=ionNames
    residueList['water']=waterNames

    residueList['dna']=DNAnames
    residueList['rna']=RNAnames
    residueList['nucleotides']=Nucleotides
    residueList['aminoacids']=AAnames
    residueList['all']=allResidueNames
    residueList['ligand']={}
    
    #dictionary mapping from residue types to 1 char identifier
    r_keyD = {}
    r_keyD['ALA']='A'
    r_keyD['ARG']='R'
    r_keyD['ASN']='N'
    r_keyD['ASP']='D'
    r_keyD['ASX']='B'
    r_keyD['CYS']='C'
    r_keyD['GLU']='E'
    r_keyD['GLN']='Q'
    r_keyD['GLX']='Z'
    r_keyD['GLY']='G'
    r_keyD['HIS']='H'
    r_keyD['ILE']='I'
    r_keyD['LEU']='L'
    r_keyD['LYS']='K'
    r_keyD['MET']='M'
    r_keyD['PHE']='F'
    r_keyD['PRO']='P'
    r_keyD['SER']='S'
    r_keyD['THR']='T'
    r_keyD['TRP']='W'
    r_keyD['TYR']='Y'
    r_keyD['VAL']='V'
    #FIX THIS HACK
    r_keyD['XAA']='X' #from www.ensembl.org/Docs/Pdoc/bioperl-live/Bio/SeqUtils.html
    r_keyD['SEL']='U'
    r_keys = r_keyD.values()


    def __init__(self):
        TreeNodeSetSelector.__init__(self)
        self.level = ResidueSet


    def matchSequence(self, nodes, item):
        result= ResidueSet()
        #print 'seeking ', item, ':'
        parents = nodes.parent.uniq()
        for ch in parents:
            res_nodes = ch.residues
            n_str = ""
            for r in res_nodes:
                #J is not currently used for any res type so
                #use it for residues whose type is outside the dict
                n_str = n_str + self.r_keyD.get(r.type,'J')
                #n_str = n_str + r_keyD[r.type]
            #print 'in chain ', ch.name, '->', n_str
            ind1 = find(n_str, item)
            if ind1==-1:
                continue
                #return None
            result.extend(res_nodes[ind1:ind1+len(item)])
            for i in range(ind1+1, len(res_nodes)):
                ind1 = find(n_str[i], item)
                if ind1>-1:
                    result.extend(res_nodes[ind1:ind1+len(item)])
        if not len(result):
            return None
        return result


    #FOR RESIDUES:
    def getRange(self, nodes, item):
        if len(nodes)<2:
            return None
        levItList=split(item, '-')
        if len(levItList)!=2: return None
        #if levItList[0][0]=='#' or levItList[1][0]=='#':
        if levItList[0][0]=='#' and levItList[1][0]=='#':
            return self.getResidueRelRange(nodes, item)
        firstNodes = self.processListItem(nodes, levItList[0])
        lastNodes = self.processListItem(nodes, levItList[1])
        if firstNodes and lastNodes:
            return self.rangeMatch(nodes, firstNodes[0],lastNodes[-1])
        else:
            return None

    
    def getResidueRelRange(self, nodes, item):
        #this needs to be done on a PER CHAIN basis:
        levItList = split(item, '-')
        selNodes = None
        parentNodes = ChainSet(nodes.parent.uniq())
        for par in parentNodes:
            nds = ResidueSet(filter(lambda x, par=par: x.parent==par, nodes))
            if len(nds)<2: continue
            firstNodes = self.processListItem(nds, levItList[0])
            lastNodes = self.processListItem(nds, levItList[1])
            if firstNodes and lastNodes: 
                newNodes = self.rangeMatch(nds,firstNodes[0],lastNodes[-1])
                if newNodes:
                    if selNodes: selNodes = selNodes + newNodes
                    else: selNodes = newNodes
            else: continue
        return selNodes


    def processListItem(self, nodes, item, sets=None):

        # check for pre-defined filtering lists
        if item.lower() in self.residueList.keys():
            item = item.lower()

            # lists might have been extended
            from MolKit.PDBresidueNames import DNAnames, RNAnames, \
                 Nucleotides, AAnames, ionNames, waterNames, allResidueNames
            self.residueList['ions']=ionNames
            self.residueList['water']=waterNames
            self.residueList['DNA']=DNAnames
            self.residueList['RNA']=RNAnames
            self.residueList['Nucleotides']=Nucleotides
            self.residueList['AminoAcids']=AAnames
            self.residueList['all']=allResidueNames

            if item=='ligand':
                d = self.residueList['all']
                newNodes = [x for x in nodes if not d.has_key(x.type.strip().upper())]
            else:
                d = self.residueList[item]
                newNodes = [x for x in nodes if d.has_key(x.type.strip().upper())]

            #names = self.residueList[item]
            #newNodes = filter(lambda x, names=names, nodes=nodes:
            #                  x.type in names, nodes)
            return self.level(newNodes)

        elif self.testSequence(item):
            # detect ambiguous residue sequences such as ASP vs Alanine-Serine-Proline
            # and all char in 'ARNDCEQGHILKMFPSTWYV'
            #print "testSequence: item=", item 
            newNodes = self.matchSequence(nodes, item)
            #newNodes = FD['resSeq'](item, nodes)
            return newNodes

        else:
            return TreeNodeSetSelector.processListItem(self, nodes, item, sets=sets)


    def testR(self,c):
        #t = 'ARNDCEQGHILKMFPSTWYV'
        return c in self.r_keys
        #return c in r_keyD.values()
        #return c in t


    def testSequence(self, item):
        import numpy.oldnumeric as Numeric
        try:
            ans = Numeric.add.reduce(map(self.testR,item))==len(item)
        except:
            ans = 0
        return ans


    def regexp(self, nodes, item):
        #use self.procFunction to build a regexp
        #Residue strings match to 'type', so use 'objectsFromStringField'
        if not self.caseSensitive:
            if self.escapeCharacters:
                item = self.processStringcIWEC(item)
                #newNodes = self.processStringcIWEC(item)
            else:
                item = self.processStringcI(item)
                #newNodes = self.processStringcI(item)
        else:
            item = self.processStringcS(item)
            #newNodes = self.processStringcS(item)
        try:
            #to match 13 to residue whose number is '13'
            #t = int(item[0])
            t = int(item)
            return self.level(nodes.objectsFromString(item, 'number'))
            #return self.level(nodes.objectsFromStringField(item, 'number'))
        except ValueError:
            #if not a number,  
            #do something like match 'THR1' to residue whose name is 'THR1'
            return self.level(nodes.objectsFromString(item))
            #return nodes.ReturnType(nodes.objectsFromStringField(reItem, 'type'))


#######################################################################
##  CHAINS                                                           ##
#######################################################################
  
class ChainSet(TreeNodeSet):
    """Class to extend a TreeNodeSet with Chain specific methods"""

    def __init__(self, objects=None, stringRepr=None, comments="", keywords=[]):
        TreeNodeSet.__init__(self, objects, Chain, stringRepr, comments, keywords)

        if stringRepr is None:
            ## build string repr
            # find all molecules
            strr = ''
            if objects is not None and len(objects):
                mols = [x.top for x in objects]
                molDict = {}.fromkeys(mols, [])

                # find chains in each molecule
                for a in objects:
                    molDict[a.top].append(a)

                for k,v in molDict.items():
                    if len(v)==len(k.chains): # special case all chains
                        strr += k.name+':;'
                    else:
                        strr += k.name+":"
                        for a in v:
                            strr += a.name+','
                        strr += ';'
            
            self.stringRepr = strr

#    def get(self, selectionString, selector=None, sets=None,
#                    caseSensitive=True, escapeCharacters=False):
#        if selector is None: 
#            selector = ChainSetSelector()
#            selector.caseSensitive = caseSensitive
#            selector.escapeCharacters = escapeCharacters
#        #results, msg = selector.select(self, selectionString, sets, 
#        selectionStringRepr = str(selectionString)
#        results = selector.processListItem(self, selectionString, sets) 
#        selectionStringRepr = '&' + selectionStringRepr
#        results.setStringRepr(selectionStringRepr)
#        return results

    def getSelector(self):
        if self.selector==None:
            self.selector = ChainSetSelector()
        return self.selector

##     def get(self, selectionString, selector=None, sets=None,
##                     caseSensitive=True, escapeCharacters=False,
##                     returnMsg=False):
##         """
##         selects from self.data, objects of self.elementType specified by selectionString
##         """
##         #print "CS: get: selectionString=", selectionString
##         if selector is None: 
##             selector = ChainSetSelector()
##             selector.caseSensitive = caseSensitive
##             selector.escapeCharacters = escapeCharacters
##         selectionStringRepr = '&' + str(selectionString)
##         if type(selectionString)==StringType:
##             result, msg = selector.select(self, selectionString)
##             result = self.ReturnType(result)
##             result.setStringRepr(selectionStringRepr)
##             if returnMsg: result = (result, msg)
##             return result
##         elif callable(selectionString):
##             result = filter(selectionString, self.data)
##             if len(result)==len(self.data):
##                 return self
##             else:
##                 result = self.ReturnType(result)
##                 result.setStringRepr(selectionStringRepr)
##                 return result
##         else:
##             raise RuntimeError("argument has to be a function or a string")

        

class Chain(ProteinMolecule):
    """Class to represent chains or residues. Inherits from tree element"""


##     def __del__(self):
##         self.__dict__.clear()
##         ProteinMolecule.__del__(self)



    def __init__(self, id=None, parent=None, elementType=Residue, list=None,
                 childrenName='residues', setClass=ChainSet,
                 childrenSetClass=ResidueSet, top=None, childIndex=None,
                 assignUniqIndex=1):
        """Chain constructor.
        Arguments:
        id (string)
        optional parent (instance of a TreeNode)
        optional elementType (instance of class inheriting from TreeNode)"""

        Molecule.__init__(self, str(id), parent, elementType, list,
                          childrenName, setClass, childrenSetClass, top,
                          childIndex, assignUniqIndex)
        self.id = id
        self.hasBonds = 0 # 1 for bondsByDistance is supported
        self.gaps = []    # list to store tuple (Res, Res) define the edge of a gap

    def shortestDist(self, atoms1, atoms2):
        min = 99999.9
        for a1 in atoms1:
            c1 = a1.coords
            for a2 in atoms2:
                c2 = a2.coords
                diff = (c1[0]-c2[0], c1[1]-c2[1], c1[2]-c2[2])
                d2 = (diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2])
                if d2 < min:
                    min = d2
                    as1 = a1
                    as2 = a2
        return as1, as2, min

    def connectResidues(self, res1, res2, cut_off = 1.85):
        """ Connect residues  based on distance"""
        self.bbDisconnectedAfter = []
        #cut_off2 = cut_off*cut_off
        #diff = [0,0,0]
        p = 0
        # Get the C atom of the first residue res1 :
        c = res1.atoms.get(lambda x: split(x.name)[0]=='C')
        
        if c is None or len(c) == 0: c = res1.atoms.get(lambda x:
                                        split(x.name)[0]=='O3*')
        # Get the N atom of the second Residue only if the first residue
        # has a C atom.
        if not c is None and len(c) != 0:
            # get c coords
            cx, cy, cz = c[0].coords
            cov_radc = c[0].bondOrderRadius

            n = res2.atoms.get(lambda x: split(x.name)[0]=='N')
            if n is None or len(n) == 0:
                n = res2.atoms.get(lambda x:split(x.name)[0]=='P')
            if n is not None and len(n)!=0:
                nx, ny,nz = n[0].coords
                cov_radsum = (cov_radc + n[0].bondOrderRadius)*1.1
                diffx = cx-nx
                diffy = cy-ny
                diffz = cz-nz
                d = (diffx*diffx) + (diffy*diffy) + (diffz*diffz) 
                if d<cov_radsum*cov_radsum:
                    resConnectionFound = 1
                    atom1 = c[0]
                    atom2 = n[0]
                    if not atom1.isBonded(atom2):
                        Bond( atom1, atom2, origin='BuiltByDistance', check=0 )
                else:
                    resConnectionFound = 0
                    self.bbDisconnectedAfter.append(res1)
            else: # N not found
                resConnectionFound = 0
        else: # C not found
            resConnectionFound = 0

        if not resConnectionFound:
            at1,at2,dist2 = self.shortestDist(res1.atoms,
                                              res2.atoms)
            cov_radsum = (at1.bondOrderRadius + at2.bondOrderRadius) *1.1
            if dist2<cov_radsum *cov_radsum:
                if not at1.isBonded(at2):
                    Bond( at1, at2, origin='BuiltByDistance', check=0 )
            else:
                self.bbDisconnectedAfter.append(res1)
        
    def buildBondsByDistance(self, cut_off=1.85):
        """Build bonds between atoms inside this chain, based on distance"""
        if self.hasBonds: return
        for i in xrange(len(self.residues)):
            res = self.residues[i]
            length = res.buildBondsByDistance()
            if upper(res.type) in ['HOH', 'DOD']: continue
            # Now we try to connect residues together.
            if i < len(self.residues)-1:
                self.connectResidues(res, self.residues[i+1], cut_off)
        self.hasBonds = 1


    def ribbonType(self, noCache=False):
        """
        type <- chain.ribbonType(noCache=False)

        this function compares the number of amino acids and nucleotides in
        the chain. If there are no amino acids and no nucleotides it will set
        chain._ribbonType to None else it will set this attribute to 'NA' if
        there are more nucleotides than amino acids or 'AA' if there are
        more amino acids than nucleic acids.
        the _ribbonType attribute is returned
        if the attribute self._ribbonType is found we return it  unless
        noCache is True.
        The list of nucleic acides is saves in self.DNARes and the list of
        amino acids in self.AARes
        """

        if hasattr(self, '_ribbonType') and noCache==False:
            return self._ribbonType
            
        # our implementation allows only for one Sheet2D per chain
        # this means we cannot have a ribbon for amino acids and one for
        # nucleic acids in the same chain. hence, we decide make ribbons
        # for the type that is mostly present in the chain
        from MolKit.PDBresidueNames import Nucleotides, AAnames

        self.DNARes = ResidueSet([x for x in self.residues if \
                                  Nucleotides.has_key(x.type.strip().upper())])

        self.AARes = ResidueSet([x for x in self.residues if \
                                 AAnames.has_key(x.type.strip().upper())])

        if len(self.DNARes)==0 and len(self.AARes)==0:
            self._ribbonType = None

        elif len(self.DNARes)>len(self.AARes):
            self._ribbonType = 'NA'
        else:
            self._ribbonType = 'AA'

        return self._ribbonType

    
    def isDna(self):
        """
        checks if the chain is DNA or not.
        """
        from MolKit.PDBresidueNames import Nucleotides

        dnaRes = [x for x in self.residues if \
                  Nucleotides.has_key(x.type.strip())]

        water = [x for x in self.residues if x.type in ['HOH', 'WAT']]

        if len(dnaRes) and len(dnaRes)+len(water) == len(self.residues):
            self.isDNA = True
            return True
        else:
            self.isDNA = False
            return False


    def isProteic(self):
        """ checks if the chain is proteic or not."""
        from MolKit.PDBresidueNames import AAnames

        self.AARes = [x for x in self.residues if AAnames.has_key(x.type)]

        water = [x for x in self.residues if x.type in ['HOH', 'WAT']]

        if len(self.AARes) and len(self.AARes)+len(water) == len(self.residues):
            return True
        else:
            return False


    def isHetatmChain(self):
        """ checks if is whole chain of hetatms """
        n = filter(lambda x: not x.hetatm, self.residues.atoms)
        if n: return 0
        else: return 1


    def secondaryStructure(self, ssBuilder):
        """ create a secondarystructureset. If secondarystructureset can't be
        obtained, none is created and a warning is printed."""
        from MolKit.getsecondarystructure import GetSecondaryStructure
        assert isinstance(ssBuilder, GetSecondaryStructure)
##          ss = ssBuilder.createSSNodesForChain(self)
        ssBuilder.createSSNodesForChain(self)
        if not self.secondarystructureset:
            delattr(self, 'secondarystructureset')


class ChainSetSelector(TreeNodeSetSelector):
    #what else to do here?
    def __init__(self):            
        TreeNodeSetSelector.__init__(self)
        self.level = ChainSet


    def processListItem(self, nodes, item, sets=None):
        # special treatment for items 'proteic'  and 'dna'
        #chains = nodes.findType(Chain).uniq()
        #for chain in chains:
        #    chain.ribbonType()
        if item =='proteic':
            newNodes = filter(lambda x, nodes=nodes:
                              x.isProteic(), nodes)
            return self.level(newNodes)
        elif item == 'dna':
            newNodes = filter(lambda x, nodes=nodes:
                              x.isDna(), nodes)
            return self.level(newNodes)

        else:
            return TreeNodeSetSelector.processListItem(self, nodes, item, 
                                            sets=sets)


    def regexp(self, nodes, item):
        #use self.procFunction to build a regexp
        #Chain strings match to 'id', so use 'objectsFromStringField'
        if not self.caseSensitive:
            if self.escapeCharacters:
                item = self.processStringcIWEC(item)
                #newNodes = self.processStringcIWEC(item)
            else:
                item = self.processStringcI(item)
                #newNodes = self.processStringcI(item)
        else:
            item = self.processStringcS(item)
        return self.level(nodes.objectsFromString(item, 'id'))
        #return self.level(nodes.objectsFromStringField(item, 'id'))
    
   



#######################################################################
##  PROTEIN                                                          ##
#######################################################################
from MolKit.molecule import MoleculeSet

class ProteinSet(MoleculeSet):
    """Class to extend a TreeNodeSet with molecule specific methods"""

    def __init__(self, objects=None, stringRepr=None, comments="", keywords=[]):
        MoleculeSet.__init__(self, objects, stringRepr, comments=comments,
                                keywords=keywords)
        self.elementType = Protein


#    def get(self, selectionString, selector=None, sets=None,
#                    caseSensitive=True, escapeCharacters=False):
#        if selector is None: 
#            selector = ProteinSetSelector()
#            selector.caseSensitive = caseSensitive
#            selector.escapeCharacters = escapeCharacters
#        #results, msg = selector.select(self, selectionString, sets, 
#        selectionStringRepr = str(selectionString)
#        results = selector.processListItem(self, selectionString, sets) 
#        selectionStringRepr = '&' + selectionStringRepr
#        results.setStringRepr(selectionStringRepr)
#        return results

    def getSelector(self):
        if self.selector==None:
            self.selector = ProteinSetSelector()
        return self.selector

        
##     def get(self, selectionString, selector=None, sets=None,
##                     caseSensitive=True, escapeCharacters=False,
##                     returnMsg=False):
##         """
##         returns from self.data, objects of self.elementType specified by selectionString
##         """
##         #print "PS: get: selectionString=", selectionString
##         selectionStringRepr = '&' + str(selectionString)
##         if selector is None: 
##             selector = ProteinSetSelector()
##             selector.caseSensitive = caseSensitive
##             selector.escapeCharacters = escapeCharacters
##         if type(selectionString)==StringType:
##             result, msg = selector.select(self, selectionString)
##             result = self.ReturnType(result)
##             result.setStringRepr(selectionStringRepr)
##             if returnMsg: result = (result, msg)
##             return result
##             #if selector.select(self, selectionString)[0] is None:
##             #    return self.ReturnType([])
##             #else:
##             #    return selector.select(self, selectionString)[0]
##         elif callable(selectionString):
##             result = filter(selectionString, self.data)
##             if len(result)==len(self.data):
##                 return self
##             else:
##                 result = self.ReturnType(result)
##                 result.setStringRepr(selectionStringRepr)
##                 return result
##         else:
##             raise RuntimeError("argument has to be a function or a string")


class Protein(ProteinMolecule):
    """Class to represent a protein.
     A protein is a hierarchical structure made of chains, residues and atoms.
     By definition a Protein is a list of chains (inheritence from TreeNode)
     For efficiency reasons the protein also stores a list of residues
     and atoms

     Read methods are provided to handle various PDB file format flavors 
     """

    def __init__(self, name='NoName', parent=None, elementType=Chain,
                 list=None, childrenName='chains', setClass=ProteinSet,
                 childrenSetClass=ChainSet, top=None, childIndex=None,
                 assignUniqIndex=1):
        """Protein constructor.
        Arguments:
        name (string)
        optional parent (instance of a TreeNode)
        optional elementType (instance of class inheriting from TreeNode,
        defaults to Chain)"""

        Molecule.__init__(self, name=name, parent=parent,
                          elementType=elementType, list=list,
                          childrenName=childrenName,
                          setClass=setClass,
                          childrenSetClass=childrenSetClass, top=top,
                          childIndex=childIndex, assignUniqIndex=assignUniqIndex)
        self.bondsflag = 0
        self.hasSS = []
        self.hasBonds = 0 # 1 for bondsByDistance is supported

    def copy(self, newname=None):
        """copy makes a new Protein instance with 'newname' and 
        other protein level parameters from self. Next,self.allAtoms is copied
        atom by atom. First: '_fit_atom_into_tree', which uses the same
        logic as pdbParser, builds up new instances of residues and chains
        as necessary.  Then: _copy_atom_attr copies the remaining
        String, Int, Float, None, List and Tuple attributes into new atom
        instances. The new molecule is returned by copy. 
        NB: subsequently the two copies can be visualized: 
        copy2=mv.Mols[0].copy()
        mv.addMolecule(copy2)
        mv.GUI.VIEWER.TransformRootOnly( yesno=0)
        mv.GUI.VIEWER.currentObject=copy2.geomContainer.geoms['master']
        then mouse movements would move only copy2, the new object """

        if not newname: newname = self.name + "_copy"
        newmol=Protein(name=newname, parent=self.parent,
            elementType=self.elementType, childrenName=self.childrenName,
            setClass=self.setClass, childrenSetClass=self.childrenSetClass,
            top=self.top)
        newmol.curChain=Chain()
        newmol.curRes=Residue()
        newmol.allAtoms= AtomSet()
        newmol.parser = self.parser
        for at in self.allAtoms:
            self._fit_atom_into_tree(newmol, at)
        newmol.buildBondsByDistance()
        return newmol

    def _fit_atom_into_tree(self,newmol,at):
        chainID = at.parent.parent.id
        if chainID != newmol.curChain.id:
            newmol.curChain = Chain(chainID, newmol, top=newmol)
        resName = at.parent.name
        resNum =at.parent.number
        if resName != newmol.curRes.name or resNum !=newmol.curRes.number:
            newmol.curRes=Residue(resName[:3], resNum, newmol.curChain,
                                  top=newmol)
        newat=Atom(at.name, newmol.curRes, at.element, top=newmol)
        self._copy_atom_attr(newat,at)
        newmol.allAtoms.append(newat)

    def _copy_atom_attr(self, newat, at):
        for item in at.__dict__.items():
            if type(item[1]) in [NoneType,StringType,
                                 IntType,FloatType]:
                exec('newat.%s=item[1]' %item[0])
            if type(item[1]) in [ListType, TupleType]:
                exec('newat.%s=item[1][:]' %item[0]) 


    def buildBondsByDistance(self, cut_off=1.85):
        """Build bonds between atoms inside this chain, based on distance"""
        if self.hasBonds: return

        if bhtreeFlag:
            for c in self.chains:
                c.buildBondsBhtree()
                c.hasBonds = 1
                c.residues.hasBonds = 1
            self.hasBonds=1
        else:
            for c in self.chains:
                c.buildBondsByDistance(cut_off)
                c.hasBonds = 1
                c.residues.hasBonds = 1
            self.hasBonds = 1

        self.bondsflag = 1    


    def secondaryStructure(self,ssBuilder) :

        from MolKit.getsecondarystructure import GetSecondaryStructure
        assert isinstance(ssBuilder, GetSecondaryStructure)

        for c in self.chains:
            #if c.isDna(): continue
            #if c.isHetatmChain(): continue
            if not ssBuilder.ssDataForMol.has_key(c.id): continue
            else:
                c.secondaryStructure(ssBuilder)     


    def secondaryStructureFromFile(self):
        """Function which a crate an instance of
        GetSecondaryStructureFromFile to add the secondarystructurelevel
        to the molecule hierarchy"""
        from MolKit.getsecondarystructure import GetSecondaryStructureFromFile
        ssBuilder = GetSecondaryStructureFromFile(self)
        self.builder = ssBuilder
        self.secondaryStructure(ssBuilder)
        self.hasSS = ['From File']       


    def secondaryStructureFromStride(self):
        """Function which a creat an instance of
        GetSecondaryStructureFromStride to add the secondarystructurelevel
        to the molecule hierarchy."""
        from MolKit.getsecondarystructure import GetSecondaryStructureFromStride
        ssBuilder = GetSecondaryStructureFromStride(self)
        self.secondaryStructure(ssBuilder)
        self.hasSS = ['From Stride']


    def secondaryStructureFromPross(self):
        """Function which create an instance of
        GetSecondaryStructureFromPross to add the secondarystructurelevel
        to the molecule hierarchy and make an attribute to the builder"""
        from MolKit.getsecondarystructure import GetSecondaryStructureFromPross
        ssBuilder = GetSecondaryStructureFromPross(self)
        self.secondaryStructure(ssBuilder)
        self.hasSS = ['From Pross']       
        self.builder = ssBuilder



class ProteinSetSelector(TreeNodeSetSelector):
    #what else to do here?
    def __init__(self):
        TreeNodeSetSelector.__init__(self)
        self.level = ProteinSet




#############################################################################
#
#  section defining objects for protein's secondary structure
#
#############################################################################

class SecondaryStructureSet(TreeNodeSet):
    """class to represent a set of secondary structure elements
    typically for a protein's chain"""

    def __init__(self, objects=None, stringRepr=None, comments='', keywords=[]):
        TreeNodeSet.__init__(self, objects, SecondaryStructure, stringRepr, 
                                comments, keywords)

    def __repr__(self):
        if len(self.data):
            ob = self.data[0]
            return '<%s instance> holding %d %s' %(
                self.__class__.__name__, len(self.data),
                ob.__class__.__bases__[0].__name__)
        else:
            return '<%s instance> empty' %(self.__class__.__name__)


class SecondaryStructure(TreeNode):
    """Base class to represent a Secondary Structure element such as Helix,
    Sheet, etc..."""
##      def __del__(self):
##          print 'free %s'%self.__class__, self.name
##          Tree.__del__()
##          print self.__class__._numberOfDeletedNodes
    def __init__(self, chain=None, structureType=None, index=None,
                 start=None, end=None, parent=None,
                 elementType=Residue, list=None,
                 childrenName='residues',
                 setClass=SecondaryStructureSet,
                 childrenSetClass=ResidueSet, top=None,
                 childIndex=None, assignUniqIndex=1, createNewLevel=1):
    
        TreeNode.__init__(self,structureType+str(index), parent,
                          elementType, list, childrenName, setClass,
                          childrenSetClass, top, childIndex, assignUniqIndex)
        self.index = index
        self.structureType = structureType
        self.start = start
        self.end = end
        self.chain = chain
        self.children = self.chain.residues.get(start.name+'-'+end.name)
        if createNewLevel:
            self.createNewLevel()
        self.residues = self.children

    def createNewLevel(self):
        
        if hasattr(self.chain,'secondarystructureset'):
            for i in self.children:
                i.secondarystructure = self
    
class Helix(SecondaryStructure):
    """Class to represent an helix inherits from SecondaryStructure."""

    def __init__(self, chain=None, index=None,
                 start=None, end=None, parent=None,
                 list=None, top=None, createNewLevel=True,
                 helClass=1, comment=None ):

        """
        optional argument:
        chain          -- Chain instance to which the secondary structure
                          belongs to
        index          -- Helix index
        start          -- N-terminal residue of the helix
        end            -- last residue of the helix
        parent 
        list
        top
        createNewLevel -- Boolean flag to specify whether or not to create a
                          new level in the tree representation of the
                          molecule for the SS.
        helClass       -- Helix class number (PDB record 39-40)
                          1 (default)   Right-handed alpha
                          2             Right-handed omega
                          3             Right-handed pi
                          4             Right-handed gamma
                          5             Right-handed 310
                          6             Left-handed alpha
                          7             Left-handed omega
                          8             Left-handed gamma
                          9             27 ribbon/helix
                          10            Polyproline
        
        comment         -- String describing the helix (PDB record 41-70
        """
        self.comment = comment
        # self.helDescr:
        # key: PDB helClass or Stride HelClass
        # value: {'helType':type of Helix, 'helDir': direction of the helix}
        # the direction of the helix can either be Right-handed, Left-handed
        # or None for unknown
        self.helDescr = {1:{'helType':'alpha',
                            'helDir':'Right-handed'},
                         2:{'helType':'omega',
                            'helDir':'Right-handed'},
                         3:{'helType':'pi',
                            'helDir':'Right-handed pi'},
                         4:{'helType':'gamma',
                            'helDir':'Right-handed'},
                         5:{'helType':'310',
                            'helDir':'Right-handed'},
                         6:{'helType':'alpha',
                            'helDir':'Left-handed'},
                         7:{'helType':'omega',
                            'helDir':'Left-handed'},
                         8:{'helType':'gamma',
                            'helDir':'Left-handed'},
                         9:{'helType':'27 ribbon.helix',
                            'helDir':None},
                         10:{'helType':'Polyproline',
                             'helDir':None},
                         'H':{'helType':'alpha',
                              'helDir':None},
                         'G':{'helType':'310',
                              'helDir':None},
                         'I':{'helType':'pi',
                              'helDir':None}}

        if helClass is None or not self.helDescr.has_key(helClass):
            self.helType = None
            self.helDir = None
        else:
            self.helType = self.helDescr[helClass]['helType']
            self.helDir = self.helDescr[helClass]['helDir']

        # This is needed to be able to write out the Helix information
        # in a PDB file
        if helClass in xrange(1,10): self.helClass = helClass
        elif helClass == 'H': self.helClass = 1
        elif helClass == 'G': self.helClass = 5
        elif helClass == 'I': self.helClass = 3
        else: self.helClass = 1
            
    
        SecondaryStructure.__init__(self, chain=chain,
                                    structureType='Helix', index=index,
                                    start=start,end=end, parent=parent,
                                    elementType=Residue, list=list,
                                    childrenName='residues',
                                    setClass=SecondaryStructureSet,
                                    childrenSetClass=ResidueSet,
                                    top=top, createNewLevel=createNewLevel)

class Strand(SecondaryStructure):
    """Class to represent a sheet inherits from SecondaryStructure."""

    def __init__(self, chain=None, index=None,
                 start=None, end=None, parent=None,
                 list=None,top=None,createNewLevel=1, 
                 nbStrand=None, sense=None):
        """
        optional argument:
        chain          -- Chain instance to which the secondary structure
                          belongs to
        index          -- Helix index
        start          -- N-terminal residue of the strand
        end            -- last residue of the strand
        parent 
        list
        top
        createNewLevel -- Boolean flag to specify whether or not to create a
                          new level in the tree representation of the molecule
                          for the SS.
        nbStrand       -- Number of strand in the sheet None if not known
        sense          -- Sense of strand with respect to previous strand in
                          the sheet
                           0 if first strand
                          -1 if anti-parallel
                           1 if parallel
                           None if not known (from stride or MOL2 files
                           
        """
        self.nbStrand = nbStrand
        self.sense = sense

        SecondaryStructure.__init__(self, chain=chain,
                                    structureType='Strand', index=index,
                                    start=start,end=end, parent=parent,
                                    elementType=Residue, list=list,
                                    childrenName='residues',
                                    setClass=SecondaryStructureSet,
                                    childrenSetClass=ResidueSet,
                                    top=top, createNewLevel=createNewLevel)


class Turn(SecondaryStructure):
    """Class to represent a turn inherits from SecondaryStructure."""
    def __init__(self, chain=None, index=None,
                 start=None, end=None, parent=None,
                 list=None, top=None, createNewLevel=1, 
                 comment=None):

        """
        optional argument:
        chain          -- Chain instance to which the secondary structure
                          belongs to
        index          -- Helix index
        start          -- N-terminal residue of the strand
        end            -- last residue of the strand
        parent 
        list
        top
        createNewLevel -- Boolean flag to specify whether or not to create
                          a new level in the tree representation of the
                          molecule for the SS.
        comment        -- String describing the turn
        """
        self.comment = comment
        SecondaryStructure.__init__(self, chain=chain,
                                    structureType='Turn', index=index,
                                    start=start,end=end, parent=parent,
                                    elementType=Residue, list=list,
                                    childrenName='residues',
                                    setClass=SecondaryStructureSet,
                                    childrenSetClass=ResidueSet,
                                    top=top, createNewLevel=createNewLevel)


class Coil(SecondaryStructure):
    """Class to represent a coil inherits from SecondaryStructure."""
    
    def __init__(self, chain=None, index=None,
                 start=None, end=None, parent=None,
                 list=None,
                 structureType='Coil',
                 top=None,
                 createNewLevel=1):
        SecondaryStructure.__init__(self, chain=chain,
                                    structureType=structureType, index=index,
                                    start=start,end=end, parent=parent,
                                    elementType=Residue, list=list,
                                    childrenName='residues',
                                    setClass=SecondaryStructureSet,
                                    childrenSetClass=ResidueSet,
                                    top=top, createNewLevel=createNewLevel)
        
        self.gapBefore = False
        self.gapAfter=  False

        
def test_secondaryStructure():
    from MolKit.pdbParser import PdbParser
    from MolKit.protein import Protein
    print 'create an object Protein crn'
    crn = Protein()
    print 'read the pdb file'
    crn.read('/tsri/pdb/struct/1crn.pdb', PdbParser())
    print 'create an object secondarystructureSet for each chain of crn'
    crn.getSS()
    print 'create the geometries for each structures of crn'
    extrudestructure = []
    for c in range(len(crn.chains)):
        for i in range(len(crn.chains[c].secondarystructureset)):
            extrudestructure.append(crn.chains[c].secondarystructureset[i].extrudeSS())
