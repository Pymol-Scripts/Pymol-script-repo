## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/MolKit/molecule.py,v 1.77.4.1 2008/11/24 21:01:13 rhuey Exp $
#
# $Id: molecule.py,v 1.77.4.1 2008/11/24 21:01:13 rhuey Exp $
#

"""
This module implements the classes Atom, AtomSet, Bond, BondSet and Molecule.
The class Molecule is expected to be specialized to represent molecular
structures with higher level of hierarchies (see protein.py)
"""

import string, types, re
#import misc
from mglutil.util import misc
from numpy.oldnumeric import sum, array, less_equal, take, nonzero, argsort
from MolKit.tree import TreeNode, TreeNodeSet, TreeNodeSetSelector
from math import sqrt
# from hbtree import bhtreelib
import numpy.oldnumeric as Numeric

global bhtreeFlag
try:
    import bhtree
    bhtreeFlag = 1
except:
    bhtreeFlag = 0
# implements a possible userPreference for what happens when
#keys are missing from some dictionaries in AtomSet._getattribute:
#if 1, an exception is raised, otherwise 'None' entry created 
raiseExceptionForMissingKey = 1

#######################################################################
##  ATOMS                                                            ##
#######################################################################
class AtomSet(TreeNodeSet):
    """Class to extend a TreeNodeSet with atom specific methods"""


    def __init__(self, objects=None, stringRepr=None, comments='', keywords=[]):
        TreeNodeSet.__init__(self, objects, Atom, stringRepr,
                             comments=comments, keywords=keywords)
        

#    def get(self, selectionString, selector=None, sets=None,
#                    caseSensitive=True, escapeCharacters=False):
#        if selector is None: 
#            selector = AtomSetSelector()
#            selector.caseSensitive = caseSensitive
#            selector.escapeCharacters = escapeCharacters
#        #results, msg = selector.select(self, selectionString, sets, 
#        selectionStringRepr = str(selectionString)
#        results = selector.processListItem(self, selectionString, sets) 
#        selectionStringRepr = '&' + selectionStringRepr
#        results.setStringRepr(selectionStringRepr)
#        return results
#        
#
    def getSelector(self):
        if self.selector==None:
            self.selector = AtomSetSelector()
        return self.selector
    
##     def get(self, selectionString, selector=None, sets=None,
##             caseSensitive=True, escapeCharacters=False,
##             returnMsg=False):
##         """
##         selects from self.data, objects of self.elementType specified by selectionString
##         """
##         #print "AS: get: selectionString=", selectionString
##         if selector is None: 
##             selector = AtomSetSelector()
##             selector.caseSensitive = caseSensitive
##             selector.escapeCharacters = escapeCharacters
##         selectionStringRepr = '%s/s/%s'%(self.stringRepr, str(selectionString))
##         if type(selectionString)==types.StringType:
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
##                 return self.ReturnType(result)
##         else:
##             raise RuntimeError("argument has to be a function or a string")



                
    def __getattr__(self, member):
        if member=='coords':# added because calling coords method is expensive
            res = []
            for a in self.data:
                res.append(a._coords[a.conformation])
            return res

        elif member=='charge':# added because multiple charges possible
            res = []
            for a in self.data:
                res.append(a._charges[a.chargeSet])
            return res

        result = TreeNodeSet.__getattr__(self, member)
        if result is None:
            return []
        if member=='bonds':
            # return bonds within bnd set and AtomSet of atoms with 0 bonds
            flat = result.uniq()
            # make the set of atoms uniq (not sure we win here)
            # uni = misc.uniq(self.data)

            # WARNING this assumes no duplicate atoms in the set
            uni = AtomSet(self.data)
            uni._bndIndex_ = range(len(uni))
            uni._inset=1
            uni._inbnd=0

            l = lambda x: hasattr(x.atom1, '_inset') and \
                hasattr(x.atom2, '_inset')
            # need to filter out the bonds with atom which don't belong to
            # the atom set.
            res = []
            for b in flat:
                if l(b):
                    res.append( b )
                    b.atom1._inbnd = 1
                    b.atom2._inbnd = 1
            result = BondSet( res )
            atnobnd = AtomSet( filter( lambda x: x._inbnd == 0, uni ) )
            for a in uni:
                if hasattr(a, '_inset'):
                    del a._inset
                if hasattr(a, '_inbnd'):
                    del a._inbnd
            return result, atnobnd

        elif member=='hbonds':
            # return hbonds 
            flat = []
            for bl in result:
                n = len(flat)
                flat[n:len(bl)] = bl
            # make uniq set of bonds (since each bond would be there thrice)
            return misc.uniq(flat)

        elif len(result) and type(result[0]) == types.DictType:
            # result is a list of dictionaries.
            # keys is a list of all keys in any dictionary
            # we build and return a dictionary having the same keys
            # and lists of values for each key.
            # Depending on 'raiseExceptionForMissingKey',
            # either an exception is raised if a dict in result is 
            # missing a key OR None is added to the list for that dict

            l = []
            for d in result:
                l = l + d.keys()
            keys = misc.uniq(l)

            if raiseExceptionForMissingKey:
                try:
                    res = {}
                    for k in keys:
                        res[k] = map(lambda x, k=k:x[k], result)
                    result = res
                except:
                    raise RuntimeError('inconsistent set of keys across nodes')
            else:
                res = {}
                for k in keys:
                    res[k] = map(lambda x, k=k: x.get(k), result)
                result = res


        return result


    def addConformation(self, coords):
        """Add a set of coordinates for alternate conformation"""
        assert len(coords)==len(self.data)
        for i in range(len(coords)):
            assert(len(coords[i])==3)
            self.data[i]._coords.append( list(coords[i]) )
            
  
    def setConformation(self, confNum):
        """Set the conformation attribute for each atom in the set"""
        assert type(confNum)==type(0)
        for a in self.data:
            assert confNum < len(a._coords)
            a.conformation = confNum


    def updateCoords(self, coords, ind=None):
        """Set ats[i]._coords[ind] to coords[i]  for each atom in the set"""
        if ind:
            self.setConformation(ind)
        else:
            #by default, change the current coords
            ind = self.data[0].conformation
        #check that the appropriate number of coords exist
        assert len(coords)==len(self.data)
        for i in range(len(self.data)):
            self.data[i]._coords[ind] = coords[i][:]


    def addCharges(self, chargeSet, newCharges):
        """Add a set of charges with key 'chargeSet'"""
        assert len(newCharges) == len(self.data)
        for a, c in map(None, self.data, newCharges):
            a._charges[chargeSet] = c
  

    def setCharges(self, chargeSet):
        """Set the charge attribute for each atom in the set"""
        assert type(chargeSet)==type('abc')
        assert chargeSet in self.data[0]._charges.keys()
        #previously this made no sense:
        #assert chargeSet not in self.data[0]._charges.keys()
        for a in self.data:
            a.chargeSet = chargeSet
   

    def delCharges(self, chargeSet):
        """Delete the _charges entry for each atom in the set"""
        assert type(chargeSet)==type('abc')
        assert chargeSet in self.data[0]._charges.keys()
        for a in self.data:
            del a._charges[chargeSet]
   

from PyBabel.babelElements import babel_elements

class Atom(TreeNode):
    """Class to represent an Atom. Inherits from tree element"""

##     def __del__(self):
##         self.__dict__.clear()
##         TreeNode.__del__(self)
##         #print 'Atom del on ', id(self)

# FIXME had to overwrite isBelow because of the Protein vs Molecule pb
    def __init__(self, name='NoName', parent=None, chemicalElement=None,
                 elementType=None,
                 list=None, childrenName=None, setClass=AtomSet,
                 childrenSetClass=TreeNodeSet, top=None,
                 childIndex=None, assignUniqIndex=1):

        """Atom constructor.
        Arguments:
        optional name (string)
        optional parent (instance of a TreeNode)
        optional chemicalElement (string) if omitted we use the atom name
                 the chemical type is used to find radii, cov_rad etc...
        optional elementType (instance of class inheriting from TreeNode)"""

        TreeNode.__init__(self, name, parent, elementType, list,
                          childrenName, setClass, childrenSetClass, top,
                          childIndex, assignUniqIndex)

        self.conformation = 0
        self.element = chemicalElement
        
        ok = 0
        if chemicalElement!=None:
            try:
                babel_elements[chemicalElement]
                ok = 1
                self.chemElem = chemicalElement
            except:
                pass

        if not ok:
            l = len(name)
            if l<5: # name is not "NoName"
                celem = string.upper(name[0])
                if l>1:
                    chemElem = celem+string.lower(name[1])
                    try:
                        babel_elements[chemElem]
                        self.chemElem=chemElem
                    except:
                        try:
                            babel_elements[celem]
                            self.chemElem = celem
                        except:
                            self.chemElem = 'Xx'

                else:
                    try:
                        babel_elements[celem]
                        self.chemElem = celem
                    except:
                        self.chemElem = 'Xx'
            else: 
                self.chemElem = 'Xx'

        data = babel_elements[self.chemElem]
        self.atomicNumber = data['num']
        self.bondOrderRadius = data['bond_ord_rad']
        self.covalentRadius = data['cov_rad']
        self.vdwRadius = data['vdw_rad']
        self.maxBonds = data['max_bonds']

        if self.chemElem[0] in [ 'C', 'H', 'O', 'N', 'S', 'P' ]:
            self.organic=1
        else: self.organic = 0

        self.bonds = BondSet([])
        self.alternate = AtomSet()
        self.altname = None
        # the colors attribute is a dictionary the ('geom':color)
        self.colors = {}
        self.opacities = {}
        
        # the _charges attribute is a dictionary, eg  {'gasteiger': 0.234}
        # self.chargeSet, a string, is key into dictionary for current charge
        # AtomSet has methods: setCharges and addCharges to change current
        # charge and to add a new set 
        self._charges = {}
        self.chargeSet = None
       
        
    def isBelow(self, Klass):
        from protein import Protein
        if Klass in (Molecule, Protein): # try both Molecule and Protein
            l = TreeNode.isBelow(self, Molecule)
            if l>0: # not below
                return l
            return TreeNode.isBelow(self, Protein)
        else:
            return TreeNode.isBelow(self, Klass)


    def isBonded(self, at2):
        l1 = len(self.bonds)
        if l1 ==0: return 0
        l2 = len(at2.bonds)
        if l2 ==0: return 0
        if l1 < l2:
            for b in self.bonds:
                a2 = b.atom1
                if a2==self: a2 = b.atom2
                if a2==at2: return 1# already bonded
            return 0
        else:
            for b in at2.bonds:
                a2 = b.atom1
                if a2==at2: a2 = b.atom2
                if a2==self: return 1# already bonded
            return 0
            

    def findHydrogens(self):
        """hlist<-self.findHydrogens()
           hlist: list of hydrogen atoms bonded to self
        """
        hlist = []
        for b in self.bonds:
            at2 = b.atom1
            if at2 == self:
                at2 = b.atom2
            if at2.element == 'H':
                hlist.append(at2)
        return hlist

    def __getattr__(self, member):
        if member=='coords':# added because calling coords method is expensive
            return self._coords[self.conformation]
        elif member=='charge':# added because calling coords method is expensive
            return self._charges[self.chargeSet]
        if not hasattr(self.__dict__, member):
            # DO NOT USE repr() in here else object becomes unpicklable
            # raise AttributeError('member %s not found in %s'%(member, repr(self)))
            raise AttributeError('member %s not found'%member)
        else:
            return self.__dict__[member]


    def getAverageCoords(self):
        """ Compute the average coords of the alternate locations of
        an atom. """
        if self.alternate != []:
            # get the alternate coords.
            coords = [self.coords,] + self.alternate.coords
            # compute the average as a list of coords.
            avCoords = (sum(coords)/len(coords)).tolist()
        else:
            # only one coords.
            avCoords = self.coords
        return avCoords

    def addBond(self, bond):
        """Add a bond to this atom"""

        assert isinstance(bond, Bond)
        self.bonds.append(bond)


    def addConformation(self, coords):
        """Add a set of coordinates for alternate conformation"""
        assert(len(coords)==3)
        self._coords.append( coords.tolist() )


class AtomSetSelector(TreeNodeSetSelector):

    atomList = {}
    atomList['backbone']=['C','CA','N','O']
    atomList['backbone+h']=['C','CA','N','O','HN','HN1','HN2', 'HA']
    atomList['sidechain']=['SG', 'SD', 'CB','CG','CD','CD1', 'CD2', 'CE', 
    'CE1', 'CE2', 'CE3', 'CG1', 'CG2', 'CZ','CZ2','CZ3','CH2',
    'ND1', 'ND2', 'NE', 'NE1', 'NE2','NH1','NH2','NZ',
    'OD1', 'OD2', 'OG', 'OG1',  'OE1', 'OE2', 'OH',
    'HD1', 'HE1', 'HE2', 'HE', 'HE21', 'HE22','HH22','HH11','HH12', 
    'HG', 'HG1', 'HH21', 'HD22', 'HD21', 'HZ1','HZ2','HZ3','HH','HB1',
    'HB2','HB3', 'HG2', 'HD2']
    atomList['hetatm'] = []

    std_residues_types = ['ala', 'ALA','arg', 'ARG','asn', 'ASN', 'asp', 'ASP',
                    'cys', 'CYS','gln', 'GLN','glu', 'GLU', 'gly', 'GLY',
                    'his', 'HIS','ile', 'ILE', 'leu', 'LEU','lys', 'LYS',
                    'met', 'MET', 'phe', 'PHE','pro', 'PRO','ser', 'SER',
                    'thr', 'THR','trp', 'TRP','tyr', 'TYR','val', 'VAL']

    def __init__(self):
        TreeNodeSetSelector.__init__(self)
        self.level = AtomSet


    def processListItem(self, nodes, item, sets=None):
        # check for pre-defined filtering lists
        if item=='hetatm':
            return nodes.get(lambda x: x.hetatm==1)
        elif item in self.atomList.keys():
            #cannot just lists of names for filtering
            #want to restrict to atoms in std residues
            return self.getNamedAtomSet( nodes, item)
        else:
            return TreeNodeSetSelector.processListItem(self, nodes, item, sets)


    #FOR ATOMS:
    def getRange(self, nodes, item):
        if len(nodes)<2:
            return None
        levItList=string.split(item, '-')
        if len(levItList)!=2: return None
        #if levItList[0][0]=='#' or levItList[1][0]=='#':
        if levItList[0][0]=='#' and levItList[1][0]=='#':
            return self.getAtomRelRange(nodes, item)
        #FIX THIS: why was chainFD used here?
        #firstNodes = self.processListItem(nodes, levItList[0], self.chainFD)
        #lastNodes = self.processListItem(nodes, levItList[1], self.chainFD)
        firstNodes = self.processListItem(nodes, levItList[0])
        lastNodes = self.processListItem(nodes, levItList[1])
        if firstNodes and lastNodes:
            return self.rangeMatch(nodes, firstNodes[0],lastNodes[-1])
        else:
            return None


    def getAtomRelRange(self, nodes, item):
        levItList=string.split(item, '-')
        #now the hard part: need to call pLI w/ each set of parent nodes
        selNodes = None
        parentNodes = nodes.parent.uniq()
        for par in parentNodes:
            nds = AtomSet(filter(lambda x, par=par: x.parent==par, nodes))
            firstNodes = self.processListItem(nds, levItList[0])
            lastNodes = self.processListItem(nds, levItList[-1])
            newNodes = None
            if firstNodes and lastNodes:
                newNodes = self.rangeMatch(nds,firstNodes[0],lastNodes[-1])
            if newNodes:
                if selNodes: selNodes=selNodes + newNodes
                else: selNodes = newNodes
        return selNodes


    def getNamedAtomSet(self, nodes, item):
        #here get all atoms w/ name in  atomList[item]
        #AND 8/2004: whose parents are std residues...
        alist = self.atomList[item]
        #only get atoms in standard residues
        reslist = self.std_residues_types
        res_atoms = filter(lambda x, nodes=nodes: x.parent.type in reslist, nodes)
        ans = filter(lambda x, alist=alist, res_atoms=res_atoms: x.name in alist, res_atoms)
        #previously:
        #ans = filter(lambda x, alist=alist, nodes=nodes: x.name in alist, nodes)
        return AtomSet(ans)




#######################################################################
##  BONDS                                                            ##
#######################################################################
class BondSetSelector(TreeNodeSetSelector):
    """Class for selecting bnds, primarily based on strings"""
    def __init__(self):
        TreeNodeSetSelector.__init__(self)
        self.level = BondSet

    
class BondSet(TreeNodeSet):
    """Class to extend a TreeNodeSet with bond specific methods"""

    def __init__(self, objects=None, stringRepr=None, comments="", keywords=[]):
        TreeNodeSet.__init__(self, objects, Bond, stringRepr, comments=comments, keywords=keywords)


    def full_name(self):
        if len(self) == 0: return ''
        name  = ""
        for d in self.data:
            name = name + d.name + ","
        return name


    def getAtoms(self):
        """ return atoms in bonds in this BondSet"""
        u = {} 
        for d in self.data:
            u[d.atom1] = 1
            u[d.atom2] = 1
        return AtomSet(u.keys())


    def getSelector(self):
        if self.selector==None:
            self.selector = BondSetSelector()
        return self.selector


class Bond:
    """Class to represent an Atomic bond. """

    _numberOfDeletedBonds = 0

    def full_name(self):
        return self.name

    def __init__(self, atom1, atom2, bondOrder=None, origin='File',
                 name='NoName', parent=None, elementType=None, list=None,
                 setClass=BondSet, check=1, addBond=1):
        """Bond constructor.
        Arguments:
        if check: 
            atom1.isBonded returns 1 if bond to atom2 already exists 
                 which raises an exception
        """
        assert isinstance(atom1, Atom)
        assert isinstance(atom2, Atom)
        self.setClass = setClass
        if check:
            if atom1.isBonded(atom2):
                msg = atom1.full_name() + ' is already bonded to '+ \
                      atom2.full_name()
                raise RuntimeError(msg)
                
        self.atom1 = atom1
        self.atom2 = atom2
        self.name = '%s == %s'%(self.atom1.full_name(), self.atom2.full_name())
        self.children = []
        
        self.bondOrder = bondOrder
        self.origin = origin
        # shared is a flag set to 1 when a bond is shared. example in aromatic
        self.shared = 0

        if addBond:
            self.atom1.addBond(self)
            self.atom2.addBond(self)

    def neighborAtom(self, atom):
        """atom <- neighborAtom(atom)  return the 'other' atom in the bond"""
        if self.atom1==atom: return self.atom2
        elif self.atom2==atom: return self.atom1
        else:
            raise RuntimeError("Atom %s not in bond %s" % (atom.full_name(),
                                                           self))
        

##     def __del__(self):
##         self.__class__._numberOfDeletedBonds =self.__class__._numberOfDeletedBonds +1

      
        
    def computeDispVectors(self):
        #if not self.atom1
        if self.bondOrder in [None, 1]\
           or (hasattr(self.atom1, 'dispVec') and \
               hasattr(self.atom2, 'dispVec')):
            return

        # Need to compute the displacement vector for each atom of the
        # bond.
        # Start with the atom1vect1 and atom2vect1 disp Vector:

        # 1- Atom1 is the origin atom and atom2 is the end atom
        atm1 = self.atom1
        #inring = hasattr(self.atom1, 'ring')
        atm1Coords = array(atm1.coords)

        #1- First vector is the vector that linked self.atom1 -> self.atom2
        # atm2 = CZ
        atm2 = self.atom2
        atm2Coords = array(atm2.coords)
        # Vect1 = atm1 -> atm2 is defined selfy the selfond you are trying
        # to draw.
        a1vect1 = atm2Coords-atm1Coords
        a2vect1 = atm1Coords - atm2Coords

        
        #2- Second vector:
        # 
        # The second vector is atm1-> atom bound to the smallest chain or
        # in a non aromatic cycle.
        mol = atm1.top
        if hasattr(mol, 'rings') and not mol.rings is None \
               and mol.rings.bondRings.has_key(self):
            indices = mol.rings.bondRings[self]
            aromCycles = filter(lambda x: mol.rings.rings[x]['aromatic'],
                                indices)
            if not aromCycles:
                # chose the first ring...
                ring = mol.rings.rings[indices[0]]
                atmsInRing = ring['atoms']

                atm1ind = atmsInRing.index(atm1)
                atm2ind = atmsInRing.index(atm2)

                if atm1ind < atm2ind :
                    atm0 = atmsInRing[atm1ind-1]
                    if atm2ind+1>=len(atmsInRing):
                        atm3 = atmsInRing[0]
                    else:
                        atm3 = atmsInRing[atm2ind+1]
                else:
                    if atm1ind+1>=len(atmsInRing):
                        atm0 = atmsInRing[0]
                    else:
                        atm0 = atmsInRing[atm1ind+1]

                    atm3 = atmsInRing[atm2ind-1]
                atm0Coords = array(atm0.coords)
                atm3Coords = array(atm3.coords)
                a1vect2 = atm0Coords-atm1Coords
                a2vect2 = atm3Coords-atm2Coords
            else:
                a1vect2 = self.getSecondVector(atm1, atm2)
                a2vect2 = self.getSecondVector(atm2, atm1)
        else:
            a1vect2 = self.getSecondVector(atm1, atm2)
            a2vect2 = self.getSecondVector(atm2, atm1)


        if a1vect2 is None:
            a1vect2 = a2vect2
            a1sum = a2vect2
        else:
            # SumVector:
            a1sum = a1vect1 + a1vect2

        # Normalize this vector,
        a1len = sqrt(sum(a1sum*a1sum))
        a1norm = a1sum/a1len
        v1 = atm1.dispVec = 0.20*a1norm

        if a2vect2 is None:
            a2vect2 = a1vect2
            a2sum = a2vect2
        else:
            # SumVector:
            a2sum = a2vect1 + a2vect2
        # Normalize this vector,
        a2len = sqrt(sum(a2sum*a2sum))
        a2norm = a2sum/a1len
        v2 = atm2.dispVec = 0.20 * a2norm

        # make sure they are on the same side
        # does not seem to be necessary
        scalar = v1[0]*v2[0] + v1[1]*v2[1] +v1[2]*v2[2]
        if scalar < 0.0:
            if a2vect2 is None:
                a2vect2 = a1vect2
                a2sum = -a2vect2
            else:
                # SumVector:
                a2sum = a2vect1 - a2vect2
            # Normalize this vector,
            a2len = sqrt(sum(a2sum*a2sum))
            a2norm = a2sum/a1len
            v2 = atm2.dispVec = 0.20 * a2norm

    def getSecondVector(self, oriAtm, endAtm):
        branches = oriAtm.bonds
        oriAtmCoords = array(oriAtm.coords)
        if len(branches) == 1:
            # No other bonds than the one we are looking at.
            # have to take the other atom of the bond vect3.
            # Need to wait to compute atom2 displacement vector
            a1vect2 = None
            #atm3 = oriAtm
            
        elif len(branches) == 2:
            # Easy only one choice for the second vector
            if branches[0] == self: branch = branches[1]
            else: branch = branches[0]
            
            if branch.atom1 == oriAtm:
                atm3 = branch.atom2
            else:
                atm3 = branch.atom1
            atm3Coords = array(atm3.coords)
            a1vect2= atm3Coords - oriAtmCoords
        else:
            results = []
            # need to look for the shortest chain or if in a cycle get
            # the atom in a the cycle.
            for branch in branches:
                if branch == self:
                    # Do not consider the double bond
                    continue
                nbBonds = branch.getShortestBranch(oriAtm, endAtm)
                results.append(nbBonds)
            distances = map(lambda x: x[1], results)

            # sorting is too expansif, when only minimum needed
            vect2 = results[ distances.index( min(distances) ) ][0]
            if vect2.atom1 == oriAtm :
                atm3 = vect2.atom2
            else:
                atm3 = vect2.atom1
            atm3Coords = array(atm3.coords)
            a1vect2 = atm3Coords - oriAtmCoords
        return a1vect2

    def getShortestBranch(self, oriAtm, endAtm):
        # Looking for the shortest branch to determine which way to translate
        # the dispVector.
        ### FIXME .. needs description
        bondStack = [self,]
        bondStackIndex = 0
        currentLength = 0
        lenStack = [0,]
        oriStack = [oriAtm,]
        while (1):
            # Get Bond to look atL
            bond = bondStack[bondStackIndex]
            # If the bond is inter Residue stop

            # if length get longer than 10 bonds we stop
            # (we are probably in main chain)
            if len(bondStack) > 10:
                return (bondStack[0], lenStack[bondStackIndex])
            
            elif self.atom1 == endAtm or self.atom2 == endAtm:
                # the bond in the cycle wins over the bond
                # outside the cycle.
                return (bondStack[0],-1)
            
            elif hasattr(self.atom1, 'babel_type') and \
                 self.atom1.babel_type == 'Car' and \
                 hasattr(self.atom2, 'babel_type') and \
                 self.atom2.babel_type == 'Car':
                return ( bondStack[0], -1)
            
            # Get the branches coming from that branch
            if self.atom1 == oriStack[bondStackIndex]:
                bonds = filter(lambda x, b= bond: x!=b, self.atom2.bonds)
                oriAtm = self.atom2
            else:
                bonds = filter(lambda x, b= bond: x!=b, self.atom1.bonds)
                oriAtm = self.atom1

            if len(bonds) == 0:
                # Done, no other branch from that branch
                return (bondStack[0], lenStack[bondStackIndex])
            else:
                for b in bonds:
                    bondStack.append(b)
                    oriStack.append(oriAtm)
                    lenStack.append(currentLength + 1)
                currentLength = currentLength + 1
                bondStackIndex = bondStackIndex + 1
        
        
    def __repr__(self):
        if hasattr(self, 'bondOrder'):
            return "<%s instance '%s','%s'>" % (self.__class__.__name__,
                                                self.name,
                                                str(self.bondOrder))
        else:
             return "<%s instance '%s'>" % (self.__class__.__name__,
                                                self.name)
           
class HydrogenBond:
    """Class to represent a Hydrogen bond. Inherits from Bond"""

      
    def __init__(self, donAt,  accAt, hAt=None, 
                hlen=None, dlen=None, theta=None, phi=None, 
                torsion=None, energy=None, origin='Hydbnd', 
                name='NoName', typ=None, check=1):
        """
        donAt: heavy atom which 'donates' hydrogen to bond
        accAt: heavy atom which 'accepts' hydrogen from bond
        hAt: 'donated' hydrogen, optional
        hlen: distance between hAt + accAt
        dlen: distance between donAt + accAt
        theta: donAt, hAt, accAt angle
        phi: hAt, accAt angle, accAtNeighbor
        torsion: angle between normals to d-h-a and h-a-aa
        energy: calc from j.Mol.Graphics Modell., 2001, Vol19, No.1, p62
            energy=Vo{5(do/d)**12-6(do/d)**10}F(theta,phi,gamma)
                different forms of F are used depending on whether
                donAt +  accAt are sp2 or sp3 hybridized 
                gamma is torsion defined by donAt-hAt-accAt-accAtNeighbor
        if check: 
            atom1.isBonded returns 1 if hbond to atom2 already exists 
                 which raises an exception
        """
        assert isinstance(donAt, Atom)
        assert isinstance(accAt, Atom)
        #assert isinstance(hAt, Atom)
        if check: 
            bonded = self.checkIsBonded(donAt, accAt)
            if bonded:
                msg = donAt.full_name() + ' is already hbonded to '+ accAt.full_name()
                raise RuntimeError(msg)
        self.donAt = donAt
        self.accAt = accAt
        self.hAt = hAt
        self.hlen = hlen
        self.dlen = dlen
        self.theta = theta
        self.phi = phi
        self.torsion = torsion
        self.energy = energy
        self.type = typ
        self.origin = 'Hydbnd'

        ##NB: when addBond not set to zero, this buils bond between
        ##hat and accAt...
        #Bond.__init__(self, donAt, accAt, origin='Hydbnd', name=name,
                        #check=0, addBond=0)
            
        

    def checkIsBonded(self, donAt, accAt):
        if not hasattr(donAt, 'hbonds'): return 0
        if not hasattr(accAt, 'hbonds'): return 0
        for b in donAt.hbonds:
            if id(b.donAt)==id(donAt) and id(b.accAt)==id(accAt):
                return 1
            if id(b.accAt)==id(donAt) and id(b.donAt)==id(accAt):
                return 1
        return 0


    def __repr__(self):
        if not self.hAt is None:
            return "<%s instance '%s' - '%s'....'%s'>" % \
                   ( self.__class__.__name__, self.donAt.full_name(),
                     self.hAt.full_name(), self.accAt.full_name())
        else:
            return "<%s instance '%s' - 'None'....'%s'>" % \
                   ( self.__class__.__name__, self.donAt.full_name(),
                      self.accAt.full_name())

    
        
class HydrogenBondSet(TreeNodeSet):
    """Class to extend a TreeNodeSet with hbond specific methods"""

    def __init__(self, objects=None, stringRepr=None, comments="", keywords=[]):
        TreeNodeSet.__init__(self, objects, HydrogenBond, stringRepr=None,
                            comments="", keywords=[])


    def __getattr__(self, member):
        result = TreeNodeSet.__getattr__(self, member)
        return self.ReturnType(result)
    

#######################################################################
##  MOLECULE                                                         ##
#######################################################################

class MoleculeSet(TreeNodeSet):
    """Class to extend a TreeNodeSet with molecule specific methods"""

    def __init__(self, objects=None, stringRepr=None, comments="", keywords=[]):
        TreeNodeSet.__init__(self, objects, Molecule, stringRepr,
                             comments, keywords)
        

    def getSelector(self):
        if self.selector==None:
            self.selector = MoleculeSetSelector()
        return self.selector

#    def get(self, selectionString, selector=None, sets=None,
#                    caseSensitive=True, escapeCharacters=False):
#        if selector is None: 
#            selector = MoleculeSetSelector()
#            selector.caseSensitive = caseSensitive
#            selector.escapeCharacters = escapeCharacters
#        #results, msg = selector.select(self, selectionString, sets, 
#        selectionStringRepr = str(selectionString)
#        results = selector.processListItem(self, selectionString, sets) 
#        selectionStringRepr = '&' + selectionStringRepr
#        results.setStringRepr(selectionStringRepr)
#        return results

##     def get(self, selectionString, selector=None, sets=None,
##                     caseSensitive=True, escapeCharacters=False,
##                     returnMsg=False):
##         """
##         selects from self.data, objects of self.elementType specified by selectionString
##         """
##         #print "MS: get: selectionString=", selectionString
##         selectionStringRepr = '&' + str(selectionString)
##         if selector is None: 
##             selector = MoleculeSetSelector()
##             selector.caseSensitive = caseSensitive
##             selector.escapeCharacters = escapeCharacters
##         if type(selectionString)==types.StringType:
##             result, msg = selector.select(self, selectionString)
##             # should this be only if len(result)?
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



    def findType(self, what, uniq=0):
        from MolKit.protein import Protein
        if self.elementType in [Molecule, Protein] and \
           what in [Molecule, Protein]: return self

        elif len(self) == 0:
            try:
                return what().setClass([])
            except:
                raise RuntimeError ("could not find level %s"%what)

        else:
            levelBelow = self[0].isAbove(what)
            if levelBelow>1:
                exec("result=self"+".children"*levelBelow)
                if uniq:
                    result = result.uniq()
                return result

            elif levelBelow==1:
                result = self.children
                if uniq:
                    result=result.uniq()
                return result
            else:
                levelAbove = self[0].isBelow(what)
                if levelAbove>1:
                    exec("result=self"+".parent"*levelAbove)
                    if uniq:
                        result = result.uniq()
                    return result

                elif levelAbove==1:
                    result = self.parent
                    if uniq:
                        result=result.uniq()
                    return result

                else:
                    raise RuntimeError ("could not find level of type %s"%type)
            
            
def dist(a1, a2):
    c1 = array(a1.coords)
    c2 = array(a2.coords)
    d = c2-c1
    return sqrt(sum(d*d))
                  
            
class Molecule(TreeNode):
    """
    Class to represent an Molecule.
    """
    from MolKit.radii_patterns import AtomRadiiPatterns, AAradii
    compiled_patterns = []

    for p in AtomRadiiPatterns:
        compiled_patterns.append( (re.compile(p[0]),re.compile(p[1]),p[2]) )

##     def __del__(self):
##         self.__dict__.clear()
##         TreeNode.__del__(self)
##         #print 'Molecule del'


    def __init__(self, name='NoName', parent=None, elementType=None,
                 list=None, childrenName='atoms', setClass=MoleculeSet,
                 childrenSetClass=AtomSet, top=None, childIndex=None,
                 assignUniqIndex=1):
        """Molecule constructor.
        Arguments:
        name (string)
        optional parent (instance of a TreeNode)
        optional elementType (instance of class inheriting from TreeNode)"""

        TreeNode.__init__(self, name, parent, elementType, list,
                          childrenName, setClass, childrenSetClass, top,
                          childIndex, assignUniqIndex)
        self.bonds = BondSet([])
        self.parser = None  # parser object to read files
        self.hasBonds = 0 # 1 for bondsByDistance is supported
        # dictionary used to parse CONECT records in pdbParsers.
        if parent==None:
            self.atmNum = {}
            self.allAtoms = AtomSet([])


    def mergeNPH(self, ats=None):
        """Merge non polar atoms in atomSet
"""
        if ats is None:
            ats = self.allAtoms

        mol = self
        mol.buildBondsByDistance()
	mol.buildBrokenBonds()
        npHSet = ats.get(lambda x:x.element=='H'\
             and (x.bonds[0].atom1.element=='C' or \
             x.bonds[0].atom2.element=='C'))
        if npHSet is None or len(npHSet)==0:
            return None
        #process each field in npHSet._charges.keys()
        #to get only the ones every npH has:
        #for a in npHSet: print a.name, ':', len(a.bonds)
        problem_npHs = npHSet.get(lambda x: len(x.bonds)>1)
        if problem_npHs is not None and len(problem_npHs)>0:
            print "The following nonpolar hydrogens had more than one bond"
            for a in problem_npHs:
                print a.full_name() + "-" + str(len(a.bonds)) + " bonds"
        mol.allAtoms = mol.allAtoms - npHSet

        chList = npHSet[0]._charges.keys()
        for at in npHSet:
            chs = at._charges.keys()
            for c in chList:
                if c not in chs: 
                    chList.remove(c)
        if not len(chList):
            s = 'charges on carbons unchanged'
            print s
        else:
            for chargeSet in chList:
                for h in npHSet:
                    Catom = h.bonds[0].atom1
                    if Catom==h: 
                        Catom = h.bonds[0].atom2
                    Catom._charges[chargeSet] = Catom.charge + h.charge
        #have to do this loop separately because there may not be charges
        for at in npHSet:
            #b = at.bonds[0]
            #could have hydrogens with >1bond!!!(YIKES)
            for b in at.bonds:
                c = b.atom1
                if c==at:
                    c = b.atom2
                c.bonds.remove(b)
            at.bonds = BondSet()
            at.parent.remove(at)
            del at 
        return npHSet


    def removeAllHydrogens(self, ats=None):
        """Merge non polar atoms in atomSet
"""
        if ats is None:
            ats = self.allAtoms

        mol = self
        mol.buildBondsByDistance()
	mol.buildBrokenBonds()
        allHSet = ats.get(lambda x:x.element=='H')
        if allHSet is None or len(allHSet)==0:
            return allHSet
        problem_allHs = allHSet.get(lambda x: len(x.bonds)>1)
        if problem_allHs is not None and len(problem_allHs)>0:
            print "The following hydrogens had more than one bond"
            for a in problem_allHs:
                print a.full_name() + "-" + str(len(a.bonds)) + " bonds"
        mol.allAtoms = mol.allAtoms - allHSet

        for at in allHSet:
            for b in at.bonds:
                c = b.atom1
                if c==at:
                    c = b.atom2
                c.bonds.remove(b)
            at.bonds = BondSet()
            at.parent.remove(at)
            del at 
        return allHSet

    

    
    def buildBrokenBonds(self):
        """Merge non polar atoms in atomSet
"""
	# Initialization
        mol = self
        mol.buildBondsByDistance()
	ats = mol.allAtoms
	atomCoords = Numeric.array( ats.coords, 'f' )
	bht = bhtree.bhtreelib.BHtree( atomCoords, None, 10 )
	closePD = bht.closePointsDist2
	cutoff  = 2.0
	indices = Numeric.zeros( (len(atomCoords),)).astype('i')
	dist2   = Numeric.zeros( (len(atomCoords),)).astype('f')

	# Search atom without bond and build one bond for it with its closest atom
	for a in ats:
	    if len(a.bonds) != 0:
	        continue
	    nb = closePD(tuple(a.coords), cutoff, indices, dist2)
	    mind = 9999999.9
	    j = 0
            for d in dist2[0:nb]:
                if d < mind:
                    mini = indices[j]
                    mind = d
		j += 1
	    atomB = ats[mini]
	    Bond( a, atomB )
        return mol


	
    def addBond(self, atom1, atom2, bondOrder=1, origin='UserDefined'):
        """
        Method that creates a bond between the two given atoms.
        """
        # The two atoms are already bonded.
        if atom1.isBonded(atom2): return
        Bond(atom1, atom2, bondOrder=bondOrder, origin=origin)

    def removeBond(self, atom1, atom2):
        """
        Method to remove the bond between the two given atoms.
        """
        if not atom1.isBonded(atom2): return
        theBond = None
        for b in atom1.bonds:
            at2 = b.atom1
            if at2 == atom1:
                at2 = b.atom2
            if at2 == atom2:
                theBond = b
                #remove this bond from the atom2.bonds list
                atom2.bonds.remove(theBond)
                #remove this bond from the atom1.bonds list
                atom1.bonds.remove(theBond)
                atom1.parent.hasBonds=0
                if atom2.parent!=atom1.parent:
                    atom2.parent.hasBonds = 0
                    if atom1.parent.parent:
                        atom1.parent.parent.hasBonds = 0
                    if atom2.parent.parent:
                        atom2.parent.parent.hasBonds = 0
                del theBond
                return
        if not theBond: return 'ERROR'
        
    def closerThan(self, point, atomList, cut_off, sort=0):
        """ AtomSet <- closerThan(point, atomList, cut_off, sort=0)
        Return all atoms of atomList withing cut_off of atom.
        if sort==1 the returned atoms are sorted from closest to furthest
        """
        c1 = array(point)
        c2 = array(atomList.coords)
        diff = c2 - c1
        dist = sum(diff*diff, 1)
        close = less_equal(dist, cut_off*cut_off)
        closeAtoms = take( atomList, nonzero(close) )
        if sort:
            d = take(dist, nonzero(close) )
            closeAtoms = take( closeAtoms, argsort(d) )
        return AtomSet(list(closeAtoms))


    def _atomRadius(self, atom, united=1 ):
        """Assign an atomic radius based on the residue type and atom name

        arguments:
        atom\t\t:atom for which the radius will be assigned
        united='yes'\t:use larger atomic radii when there are no hydrogens
        """
        for p in self.compiled_patterns:
            match = p[1].search( atom.name )
            if not match:
                #match = p[1].search( string.upper(atom.element) )
                match = p[1].search( atom.element )
            if match:  # atom name matches
                if hasattr(atom.parent, 'type') and \
                   type(atom.parent.type)==types.StringType:
                    resname = atom.parent.type
                else:
                    resname = 'XXX'

                match = p[0].search( resname )
                if match:  # residue name matches
                    n = p[2]
                    if united==1: r = self.AAradii[n][2]
                    else: r = self.AAradii[n][1]
                    atom.radius = r
                    return r

        print 'no radius found for atom %s in residue %s' % \
              (atom.name, atom.parent.name )
        atom.radius = 0.2
        return 0.2
    

    def defaultRadii(self, atomSet=None, united=None, overwrite=0):
        """Assign atom radii to the selected atoms

        radiiLst <- defaultRadii(self, atomSet=None, united=None, overwrite=0)

if atomSet is None, all atoms i nthe molecule are used
if overwrite is true pre-exiting atom.radius attribute will be overwritten
if united is true large atomic radii are used for atoms which have no hydrogens
"""
        if atomSet==None:
            atomSet = self.findType(Atom)

        r2 = []
        strType = types.StringType
        curRes = None
        mols = atomSet.top.uniq()
        if len(mols)==1: labeltext = 'Assign Atoms Radii for %s'%mols[0].name
        else:
            labeltext='Assign Atoms Radii'
            
            
        self.configureProgressBar(init=1, mode='increment',
                                  max=len(atomSet),
                                  labeltext=labeltext)
        for a in atomSet:
            withH = 0
            if hasattr(a, 'radius') and not overwrite:
                r2.append(a.radius)
                continue
            if united:
                for b in a.bonds:
                    if b.atom1.chemElem=='H' or b.atom2.chemElem=='H':
                        withH = 1
                        break
            for resPat, atomPat, pnum in self.compiled_patterns:
                if a.parent != curRes:
                    curRes = a.parent
                    if hasattr(curRes, 'type') and \
                       type(curRes.type)==strType:
                        resname = curRes.type
                    else:
                        resname = 'XXX'
                    
                if atomPat.search( a.name ) and resPat.search( resname ):
                    break

            if pnum==len(self.compiled_patterns): # not found:
                print "WARNING: no pattern match for ",a.full_name()
                r = 0.2
            else:
                if withH and united:
                    r = self.AAradii[pnum][1]
                else:
                    r = self.AAradii[pnum][2]
            a.radius = r
            r2.append(r)
            self.updateProgressBar()
            
        return r2
    


    def read(self, filename, parser):
        """parse the given file using the given  parser"""

        from MolKit.moleculeParser import MoleculeParser
        assert isinstance(parser, MoleculeParser)
        self.parser = parser
        self.parser.parse(self, filename)
        if self.name == 'NoName':
            import os
            self.name = os.path.basename(os.path.splitext(filename)[0])

    def write(self, filename, nodes, writer, sort=False, sortFunc=None,
              records=['ATOM', 'CONECT'], bondOrigin=('File',),
              ssOrigin='File'):
        """ write the data from the molecule to a file in the given writer format
        """
        from MolKit.moleculeWriter import MoleculeWriter
        assert isinstance(writer, MoleculeWriter)
        self.writer = writer
        self.writer.write(filename, nodes, sort=sort, sortFunc=sortFunc,
                          records=records, bondOrigin=bondOrigin,
                          ssOrigin=ssOrigin)


    def getAtoms(self, root=None):
        """Return an atom set of all atoms in the subtree with root"""
        n=root
        if n==None: n=self
        while len(n.children) > 0 and not isinstance(n.children[0], Atom):
            n = n.children
        if len(n.children)==0:
            raise RuntimeError(" no Atom objects found\n")
        return n.children

    def buildBondsByDistanceOnAtoms(self, atoms):
        # bond won't be built between two atoms a1 and a2 if the 2 atoms are
        # alternate of each other or if they don't have the same alternate
        # name.
        la = len(atoms)
        for i in xrange(la-1):
            a1 = atoms[i]
            c1 = array(a1.coords)
            cov_rad1 = a1.bondOrderRadius
            for j in range(i+1,la):
                a2 = atoms[j]
                if not a1.altname is None and not a2.altname is None:
                    if a1 in a2.alternate or a2 in a1.alternate:
                        continue
                    else:
                        if a1.altname != a2.altname:
                            continue
                c2 = array(a2.coords)
                cov_radsum = (cov_rad1 + a2.bondOrderRadius)*1.1
                diff = c1-c2
                d = sum(diff*diff)
                if d<cov_radsum*cov_radsum:
                    if not a1.isBonded(a2):
                        Bond( a1, a2, origin='BuiltByDistance', check=0)


    def configureProgressBar(self, **kw):
        # this method is to be implemented by the user from outside
        pass


    def updateProgressBar(self, progress=None):
        # this method is to be implemented by the user from outside
        #print 'Prorgess: ' + `progress` + '%'
        pass


    def buildBondsBhtree(self):
        from bhtree import bhtreelib
        # Create a bhtree for the molecule using a granility of 10
        atoms = self.findType(Atom)
        bht = bhtreelib.BHtree( atoms.coords, atoms.bondOrderRadius, 10)
        # find all pairs of atoms for which the distance is less than 1.1
        # times the sum of the radii
        pairs = bht.closePointsPairsInTree(1.1)
        # pairs is list of tuple of atom indices.

        if len(pairs) == 0:
            return

        # set up progress bar
        self.configureProgressBar(init=1, mode='increment', max=len(pairs),
                                  labeltext='build bonds by distance')

        for pair in pairs:
            # 1- Get the atoms corresponding of the indices of the pair
            atm1 = atoms[int(pair[0])]
            if len(atm1.bonds) > atm1.maxBonds:
                continue
            atm2 = atoms[int(pair[1])]
            # 2- Need to make sure that the two atoms are not already
            # bonded.
            if atm1.isBonded(atm2): continue

            if len(atm2.bonds) > atm2.maxBonds:
                continue

            #elif atm1.name == 'H' and atm2.name == 'H':
            #    continue
            #elif atm1.name == 'H' and len(atm1.bonds) == 1:
            #    continue
            #elif atm2.name == 'H' and len(atm2.bonds) == 1:
            #    continue

            # 3- Need to make sure that the two atoms of this pairs are not
            # alternate location.
            if not atm1.altname is None and not atm2.altname is None:
                # atm1 and atm2 are alternate locations
                if atm1.altname == atm2.altname :
                    # Can only be bound if they have the same altname
                    bond = Bond(atm1, atm2, origin='BuiltByDistance', check=0)
            # 4- Create a bonds between the two atoms:
            else:
                #bond = Bond(atm1, atm2, origin='Bhtree', check=0)
                bond = Bond(atm1, atm2, origin='BuiltByDistance', check=0)

            # set progress bar
            self.updateProgressBar()

            
    def buildBondsByDistance(self):
        """Build bonds between atoms inside this residue, based on distance
        WARNING this is a n^2 process, should only be used for small
        molecules like residues"""
        # Build the bonds in self.
        if self.hasBonds: return
        # Uses bhtree if available
        if bhtreeFlag:
            # Build the bonds between the atoms in a chain
            self.buildBondsBhtree()
            self.hasBonds = 1
            self.chains.hasbonds = 1
            self.chains.residues.hasbonds = 1
        else:
            # use the buildBondsByDistance of all the atoms.
            atoms = self.findType(Atom)
            self.buildBondsByDistanceOnAtoms(atoms)
            self.hasBonds = 1
        return len(atoms)
                        

    def _MarkTree(self,startatom):
        """Tag all atoms in a subtree. called by atomsSubTree"""
        startatom._subtree_selector=1
        for i in startatom.bonds:
                nxtAtom = i.atom1
                if nxtAtom==startatom: nxtAtom=i.atom2
                if nxtAtom._subtree_selector==0:
                        self._MarkTree(nxtAtom)
        return


    def subTree(self, atom1, atom2, nodes):
        """ Returns the list of atoms found beyond atom2 starting at atom1
        in the subtree defined by chemical bonds
        """
        #FIXME .. should first check that bonds are available
        atoms = nodes.findType( Atom )
        atoms._subtree_selector = 0
        atom1._subtree_selector = atom2._subtree_selector = 1
        self._MarkTree(atom2)
        atom1._subtree_selector = atom2._subtree_selector = 0
        result = atoms.get(lambda x: x._subtree_selector==1)
        del atoms._subtree_selector
        return AtomSet(result)


    def getCenter(self):
        """sets self.center<-getCenter(self)"""
        coords = self.allAtoms.coords
        self.center = sum(coords)/(len(coords)*1.0)
        self.center = list(self.center)
        for i in range(3):
            self.center[i] = round(self.center[i], 4)
        #print "center =", self.center
        return self.center


    def getNbrAtoms(self, centers, radius=5.0):
        """get atoms within <radius> of <centers> """
        centers=centers.findType(Atom)
        centerList=centers.coords
        coords=self.allAtoms.coords
        atSet=AtomSet([])
        ats = AtomSet([])
        atar = Numeric.array(self.allAtoms.data)
        for center in centerList:
            d = coords - Numeric.array(center, 'f')
            d = d * d
            d = Numeric.sum(d,1)
            atindex = Numeric.nonzero(Numeric.less(d, radius * radius)) 
            newats = Numeric.take(atar, atindex)        
            if len(newats)>0: 
                ats = ats + AtomSet(newats)

        if len(ats)>0:
            atSet = AtomSet(list(ats))
        return atSet.uniq()-centers
        


    def getTypeList(self, sortKeyList = None):
        """sets self.typeList<-getTypesList(self, sortKeyList = None)
           sortKeyList is optional List to reorder types"""
        elementList = self.allAtoms.element
        elementList = misc.uniq(elementList)

        if sortKeyList:
            sorted_List = []
            for item in sortKeyList:
                if item in elementList:
                    sorted_List.append(item)
            self.typeList = sorted_List
        else:
            self.typeList =  elementList


    def attach_nonbonded_fragments(self, attach_singletons=False):
        """ 
        detect bonded fragments and try to attach them via a bond
        between the closest two atoms
        detect atoms not bonded to anything and attach them to 
        closest atom
        """
        atoms = self.allAtoms
        ct = 0
        #first: build a list of bonded fragments:
        l = []   # the list to hold dicts of bonded fragments
        fb = []  # keep track of which bonds are found
        #start with unique list of bonds 'bnds'
        bond_dict = {}
        for b in atoms.bonds[0]:
            bond_dict[b] = 1  
        bnds = bond_dict.keys()
        #first pass: try to connect pieces
        for b in bnds:
            ind1 = b.atom1
            ind2 = b.atom2
            found = b in fb
            if not found:
                for d in l:
                    d_keys = d.keys()
                    if ind1 in d_keys:
                        d[ind2] = 1
                        fb.append(b)
                    elif ind2 in d_keys:
                        d[ind1] = 1
                        fb.append(b)
                if not b in fb:
                    #start new fragment 
                    l.append({ind1:1, ind2:1})
                    fb.append(b)
        #now to try to merge fragments, use lists of keys
        key_list = []
        for dict in l:
            key_list.append(dict.keys())
        #check each list of keys against following lists
        #if there are any duplications, merge current
        #into the following one..
        for i in range(len(key_list)):
            #check this list 
            kl = key_list[i]
            found = 0
            #...against the each of the subsequent ones
            for j in range(i+1, len(key_list)):
                jl = key_list[j]
                #....check each entry in subsequent one
                #.........against this list 'kl'
                for entry in jl:
                    if entry in kl:
                        #if a match
                        #...merge this dict into jth one..
                        l[j].update(l[i])
                        #.....reset this dict to {}
                        l[i] = {}
                        #.......set found flag
                        found = 1
                        #..........update jth list of keys
                        key_list[j]= l[j].keys()
                        #............skip rest of jl
                        break
                if found:
                    #.................and skip rest of key_list
                    break
        #now build a list of non-zero-length dicts
        cl = []
        for dict in l:
            len_d = len(dict)
            if len_d >0:
                cl.append(dict)
        #print "len(cl)=", len(cl)
        #now connect the closest pieces
        cdist = 10000
        closest_pair = []
        for i in range(len(cl)-1):
            cdist = 10000
            for a in cl[i].keys():
                #only allow one bond for hydrogen atoms
                if a.element=='H':
                    continue
                #need to check i-th dict vs all the other dicts in list,not just the next one
                for next_dict in cl:
                    if next_dict==cl[i]:
                        continue
                    #for next_dict in cl[i+1:]:
                    for a2 in next_dict.keys():
                        #only allow one bond for hydrogen atoms
                        if a2.element=='H':
                            continue
                        d = dist(a, a2)
                        if d<cdist:
                            cdist = d
                            #print "new closest_pair=", a.name, ':', a2.name,
                            #print "new cdist =", cdist
                            closest_pair = [a, a2]
            #AT THIS POINT MAKE A BOND between a, a2            
            print "building bond between :", closest_pair[0].name, ' and ', closest_pair[1].name
            if len(AtomSet(closest_pair).bonds[0]):
                print "SKIPPING!"
            else:
                Bond(closest_pair[0], closest_pair[1])
            ct = ct + 1
        if attach_singletons:
            #detect and deal with any orphans
            orphans = []
            for a in atoms:
                a.found = 0
                if len(a.bonds): 
                    a.found = 1
                    continue
                for d in cl:
                    if a in dict.keys():
                        a.found = 1
                        break
                if not a.found and len(a.bonds)==0:
                    orphans.append(a)

            if len(orphans):
                #print "len(orphans)=", len(orphans)
                for orph in orphans:
                    #print "processing ", orph.name
                    if orph.element=="H" and len(orph.bonds):
                        b = orph.bonds[0]
                        a2 = b.atom1
                        if a2==orph:
                            a2 = b.atom2
                        #instead of trying to attach this hydrogen
                        #attach the atom already bonded to it
                        #print "processing ", a2.name, ' instead!'
                        orph = a2
                    odist = 100000
                    o_closest_pair = []
                    attached_to_fragment = False
                    for b in orph.bonds:
                        at2 = b.atom1
                        if at2==orph:
                            at2 = b.atom2
                        if at2 not in orphans:
                            attached_to_fragment = True
                    if attached_to_fragment:
                        #print "skipping loop over atoms for ", orph.name
                        continue
                    else:
                        #find the atom closest overall 
                        for a in atoms:
                            if a==orph:
                                continue
                            #do not build any H-H bonds
                            if a.element=='H' and orph.element=='H':
                                continue
                            d = dist(orph, a)
                            if d<odist:
                                #check whether there is already a bond orph-a
                                already_bonded = False
                                for b in orph.bonds:
                                    a2 = b.atom1
                                    if a2==orph: a2=b.atom2
                                    if a2==a:
                                        already_bonded = True
                                if not already_bonded:        
                                    odist = d
                                    o_closest_pair = [orph, a]
                        #at this point, build the bond 
                        Bond(o_closest_pair[0], o_closest_pair[1])
                        ct = ct + 1

        return ct



class MoleculeSetSelector(TreeNodeSetSelector):

    def __init__(self):
        TreeNodeSetSelector.__init__(self)
        self.level = MoleculeSet



class MolecularSystem(TreeNode):
    """
    Class to represent an multiple Molecules.
    """

    def __init__(self, name='NoName', parent=None, elementType=None,
                 list=None, childrenName='molecules', setClass=None,
                 childrenSetClass=MoleculeSet, top=None, childIndex=None,
                 assignUniqIndex=1):
        """MolecularSystem constructor.
        Arguments:
        name (string)
        optional parent (instance of a TreeNode)
        optional elementType (instance of class inheriting from TreeNode)"""

        TreeNode.__init__(self, name, parent, elementType, list,
                          childrenName, setClass, childrenSetClass, top,
                          childIndex, assignUniqIndex)
        self.bonds = BondSet([])
