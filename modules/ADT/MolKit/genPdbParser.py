#############################################################################
#
# Author: Kevin Chan, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

#
#$Header: /opt/cvs/python/packages/share1.5/MolKit/genPdbParser.py,v 1.2 2007/12/04 17:53:42 annao Exp $
#
#$Id: genPdbParser.py,v 1.2 2007/12/04 17:53:42 annao Exp $
#

from MolKit.pdbParser import PdbParser
import string, sys, copy, os, types
from MolKit.moleculeParser import MoleculeParser
from MolKit.protein import Protein, Chain, ChainSet, Residue, ResidueSet, ProteinSet
from MolKit.molecule import Atom, AtomSet, Bond, BondSet


class GenPdbParser(PdbParser):

    def __init__(self, filename, recSpecs):
        PdbParser.__init__(self, filename)
        assert isinstance(recSpecs, AtomRecFieldSpecs)
        self.recSpecs = recSpecs
        # what type of parsing, by column or by index
        self.specType = recSpecs.specType
        self.atomCounter = 0  # what about model??? i think it's ok
        self.residueCounter = 0
        self.HOHCounter = 0
        self.chaincounter = 0

        
    def get_Field_Value(self, rec, fieldName):
        """ gets the value for fieldName by looking at rec according to the
        specification of self.recSpecs.  If no specification, returns ''. If
        type wrong, ..."""
        if self.specType=='c':
            beg = self.recSpecs.get(fieldName, 'from')
            end = self.recSpecs.get(fieldName, 'to')
            if beg is not None and end is not None:
                beg = beg-1    # because ColumnSpecs columns start at 1 not 0
                val = rec[beg:end]
            else: val = ''
        elif self.specType=='i':
            index = self.recSpecs.get(fieldName, 'index')
            if index is not None: val = rec[index]
            else: val = ''
        else: val = ''
        ### how to check val according to var_type ?
        # if type wrong, WARNING or Error?
        type = self.recSpecs.get(fieldName, 'var_type')
        ok = 1
        if type=='int':
            for c in val:
                if c not in (string.digits + ' -'): ok = 0
        elif type=='character':
            if len(val) != 1 or val not in (string.letters + ' '): ok = 0
        elif type=='float':
            for c in val:
                if c not in (string.digits + '. -'): ok = 0
        elif type=='alphabetic':
            for c in val:
                if c not in (string.letters + ' '): ok = 0
        elif type=='string': pass
        if not ok:
            #self.f.write('WARNING: '+ fieldName+' value not of type '+type+ ' in record:\n'+rec)
            print 'WARNING:', fieldName, 'value not of type', type, 'in record:\n', rec
            val = ''            
        return val

    def parse_PDB_atoms(self, atoms):

        #name = '/home/ktchan/pmvproject/debug.txt'
        #self.f = open(name, 'w')
        self.mol.curChain = Chain()
	self.mol.curRes = Residue()
        self.mol.levels = [Protein, Chain, Residue, Atom]
        cid = self.get_Field_Value(atoms[0], 'chainID')
        if cid:   # we know the chains
            atm = filter(lambda x, cid = cid, f = self.get_Field_Value:
                         f(x, 'chainID')==cid, atoms)
        else: atm = atoms  # assume they are all of one chain
        map( self.parse_PDB_ATOM_record,atm)
        if cid:   # only if we know the chains
            atoms = filter(lambda x, cid = cid, f = self.get_Field_Value:
                           f(x, 'chainID')!=cid, atoms)
            if len(atoms)>0: self.parse_PDB_atoms(atoms)
        #self.f.close()

    def getPDBAtomName(self, name, element):
        """Figure out the real atom name as well as the atom element, if name
        string is in pdb format"""
        # To be DONE
##          bad = 0
##          for c in element:
##              if c not in string.letters: bad = 1
        if element == '  ' or element == '':# or bad:
##              if name[0] != ' ':
##                  element = string.strip(name[0:1])
##                  # Won't work if the atomType is FE !!!!
##              else:
            element = string.strip(name[0:2])
        else:
            element = string.strip(element)
            # should I still check that it really is an element ???
        if name[1]=='H':
            if name[0] in ('1','2','3'):
                name = string.strip(name[1:])+name[0]
                element = 'H'
        if len(element)==2:
            element = element[0]+string.lower(element[1])
        elem = string.lower(element)
        if elem =='lp' or elem =='ld':
            element = 'Xx'
        return string.strip(name), element

    # if there is no specs, I'm not sure about the spaces for each variable
    def parse_PDB_ATOM_record(self, rec):
        """Parse PDB ATOM records using the pdb columns specifications"""
        self.atomCounter = self.atomCounter + 1  # not sure about altLoc
        if self.specType=='i': rec = string.split(rec)
        # Handle the alternate location using a flag.
        altLoc = self.get_Field_Value(rec, 'altLoc')
        if altLoc!= ' ': self.altLoc = altLoc
        else: self.altLoc = ''   # changed from None to ''

        # check for chains break
        #self.modlflag = modlflag
        #chainID = rec[21]+ modlflag
        hascid = 1
        chainID = self.get_Field_Value(rec, 'chainID')
        if not chainID:
            hascid = 0
            chainID = str(self.chaincounter)  ## should be unk???
        if chainID != self.mol.curChain.id :
            # has to check if the chain exists already or not !!!
            if not self.mol.chains.id or not chainID in self.mol.chains.id or \
               hascid==0:
                self.chaincounter = self.chaincounter + 1
                if hascid==0: chainID = str(self.chaincounter)
                self.mol.curChain = Chain(chainID, self.mol, top=self.mol)
                self.residueCounter = 0
            else:
                self.mol.curChain = self.mol.chains.get(chainID)[0]

        # check for residue break
        resName = self.get_Field_Value(rec, 'resName')
        resSeq = string.strip(self.get_Field_Value(rec, 'resSeq'))
        #WARNING reSeq is a STRING
        noresSeq = 0
        if not resSeq and resName==self.mol.curRes.type and resName!='HOH':
            noresSeq = 1
            resSeq = self.mol.curRes.number
        if resSeq != self.mol.curRes.number or \
           resName != self.mol.curRes.type:
            # check if this residue already exists
            na = string.strip(resName) + string.strip(resSeq)
            res = self.mol.curChain.get( na )
            if res:
                self.mol.curRes = res[0]
            else:
                self.residueCounter = self.residueCounter + 1
                if resName=='HOH': self.HOHCounter = self.HOHCounter + 1
                if not resSeq: 
                    if resName=='HOH': resSeq = self.HOHCounter
                    else: resSeq = self.residueCounter
                ## FIXME icodes are ignored
                self.mol.curRes = Residue(resName, resSeq, '',
                                          self.mol.curChain,
                                          top=self.mol)
        icode = self.get_Field_Value(rec, 'iCode')
        if not icode: pass
        elif icode != ' ': self.mol.curRes.icode = icode

        # parse atom info

        # handle atom names (calcium, hydrogen) and find element type
        # check validity of chemical element column and charge column
        ## only works if 'name' is in the pdb format!  FIX!
        n = self.get_Field_Value(rec, 'name')
        el = self.get_Field_Value(rec, 'element')
        if n:
            name, element = self.getPDBAtomName(n, el)
            # if there is not resSeq spec, use first N to figure out new res
            if noresSeq and name=='N':
                At = self.mol.curRes.get('N')
                if At:
                    self.residueCounter = self.residueCounter + 1
                    resSeq = self.residueCounter
                    self.mol.curRes = Residue(resName, resSeq,
                                              self.mol.curChain,
                                              top=self.mol)                    
            atom = Atom(name, self.mol.curRes, element, top=self.mol)
        else:
            element = el
            if element: atom = Atom(parent = self.mol.curRes,
                                    chemicalElement = element, top=self.mol)
            else: atom = Atom(parent = self.mol.curRes, top=self.mol)
##          elem = string.lower(element)          # moved to getPDBAtomName
##          if elem =='lp' or elem =='ld':
##              element = 'Xx'
        atom.charge = self.get_Field_Value(rec, 'charge')
        #should have atom.charge if no charge?
        # coords are required; where to set default or check?
        xcoord = self.get_Field_Value(rec, 'x')
        ycoord = self.get_Field_Value(rec, 'y')
        zcoord = self.get_Field_Value(rec, 'z')
        assert xcoord and ycoord and zcoord
        atom._coords = [ [ float(xcoord), float(ycoord), float(zcoord) ] ]
        atom.segID = string.strip(self.get_Field_Value(rec, 'segID'))
        if rec[:4]=='ATOM' or rec[0]=='ATOM': atom.hetatm = 0
        else: atom.hetatm = 1
        #atom.alternate = []
        atom.element = element
        num = self.get_Field_Value(rec, 'serial')
        if num: atom.number = int(num)
        else: atom.number = self.atomCounter
        occupancy = self.get_Field_Value(rec, 'occupancy')
        if occupancy: atom.occupancy = float(occupancy)
        # must check that it is a number
        atom.conformation = 0
        tempFactor = self.get_Field_Value(rec, 'tempFactor')
        if tempFactor: atom.temperatureFactor = float(tempFactor)

        # add in user defined fields to atom attributes
        for field_name in self.recSpecs.UserFieldsDict.keys():
            value = self.get_Field_Value(rec, field_name)
            type = self.recSpecs.get(field_name, 'var_type')
            if value:
                if type=='int': atom.__setattr__(field_name, int(value))
                elif type=='float': atom.__setattr__(field_name, float(value))
                else: atom.__setattr__(field_name, value)
            else: atom.__setattr__(field_name, value)
        
        if self.altLoc :
            # check if the name of the atom is the same than the
            #name of the previous atom .
            name = name + '@'+self.altLoc
            atom.name = name
            if len(self.mol.curRes.atoms)>1:
                # the new atom has been add to the current residue
                # You have to go to the one before.
                lastAtom = self.mol.curRes.atoms[-2]
                altname = string.split(lastAtom.name, '@')[0]
                if string.split(name, '@')[0] == altname:
                    # Add the new alternate atom to the LastAtom.alternate and
                    # add the lastAtom to the atom.alternate.
                    lastAtom.alternate.append(atom)
                    atom.alternate.append(lastAtom)
                    for l in lastAtom.alternate:
                        if atom.name != l.name:
                            atom.alternate.append(l)
                            l.alternate.append(atom)
        return atom


Default_FieldNames = ['serial', 'name', 'altLoc', 'resName', 'chainID',
                      'resSeq', 'iCode', 'x', 'y', 'z', 'occupancy',
                      'tempFactor', 'segID', 'element', 'charge']

# possible entries for a dict: 'field_name', 'var_type', 'from', 'to', 'index'
# will only contain the dicts the user defines
# check here that columns don't overlap and indices are not the same and
# values are integers
# if 'index', should not be 0
variableTypes = ['int', 'character', 'float', 'alphabetic', 'string']

class AtomRecFieldSpecs:

    def __init__(self):
        self.DefFieldsDict = {}
        self.UserFieldsDict = {}
        self.specType = None

    def checkEntry(self, FDict):
        pass
    
    def define_field(self, Fieldprops):
        """ Fieldprops is a dictionary describing a field specification """
        ok = self.checkEntry(Fieldprops)
        if not ok: return
        fieldName = Fieldprops['field_name']
        # if fieldName == ... (something known but called different)
        if fieldName in Default_FieldNames:
            self.DefFieldsDict[fieldName] = Fieldprops
        else: self.UserFieldsDict[fieldName] = Fieldprops

    def remove_field(self, Fieldprops):
        fieldName = Fieldprops['field_name']
        if self.DefFieldsDict.has_key(fieldName):
            del self.DefFieldsDict[fieldName]
        elif self.UserFieldsDict.has_key(fieldName):
            del self.UserFieldsDict[fieldName]
        else: pass 

    def get(self, fieldName, what):
        """ what is 'to', 'from', 'index', or 'var_type' """
        if self.DefFieldsDict.has_key(fieldName):
            if self.DefFieldsDict[fieldName].has_key(what):
                return self.DefFieldsDict[fieldName][what]
            else: return None
        elif self.UserFieldsDict.has_key(fieldName):
            if self.UserFieldsDict[fieldName].has_key(what):
                return self.UserFieldsDict[fieldName][what]
            else: return None
        else: return None

class ColumnSpecs(AtomRecFieldSpecs):

    def __init__(self):
        AtomRecFieldSpecs.__init__(self)
        self.specType = 'c'

    def checkEntry(self, FDict):
        if not ( FDict.has_key('from') and FDict.has_key('to') ):
            print 'column specifications missing'
            return 0
        if not FDict.has_key('var_type'):
            print 'var_type missing'
            return 0
        if not FDict.has_key('field_name'):
            print 'field_name missing'
            return 0
        f, t = FDict['from'], FDict['to']
        # check that specs are integers
        if type(f) != types.IntType or type(t) != types.IntType or \
           f < 1 or t < 1 or f > 80 or t > 80:
            print 'ERROR: column specification must be an integer between 1 and 80'
            return 0
        if f > t: print "ERROR: 'from' is greater than 'to'"
        if FDict['var_type'] not in variableTypes:
            print 'ERROR: unkown variable type'
            return 0 
        # check that column specs don't overlap
        for fName in (self.DefFieldsDict.keys() + self.UserFieldsDict.keys()):
            if fName != FDict['field_name']:
                beg, end = self.get(fName, 'from'), self.get(fName, 'to')
                if (f in range(beg, end+1)) or (t in range(beg, end+1)) or \
                   (beg in range(f, t+1)) or (end in range(f, t+1)):
                    print 'ERROR: overlapping column specifications'
                    return 0
        return 1          

##      def get(self, fieldName):
##          if self.FieldDict.has_key(fieldName):
##              return self.FieldDict[fieldName]['from'],\
##                     self.FieldDict[fieldName]['to'],\
##                     self.fieldDict[fieldName]['var_type']
##          else: return None

class IndexSpecs(AtomRecFieldSpecs):

    def __init__(self):
        AtomRecFieldSpecs.__init__(self)
        self.specType = 'i'

    def checkEntry(self, FDict):
        if not FDict.has_key('index'):
            print 'index missing'
            return 0
        if not FDict.has_key('var_type'):
            print 'var_type missing'
            return 0
        if not FDict.has_key('field_name'):
            print 'field_name missing'
            return 0
        i = FDict['index']
        # check that specs are integers
        if type(i) != types.IntType or i < 1:
            print 'ERROR: index specification must be a positive integer'
            return 0 
        if FDict['var_type'] not in variableTypes:
            print 'ERROR: unkown variable type'
            return 0 
        # check that index specs don't overlap
        for fName in (self.DefFieldsDict.keys() + self.UserFieldsDict.keys()):
            if fName != FDict['field_name']:
                index = self.get(fName, 'index')
                if i==index:
                    print 'ERROR: index specification already defined'
                    return 0 
        return 1
          
##      def get(self, fieldName):
##          if self.FieldDict.has_key(fieldName):
##              return self.FieldDict[fieldName]['index'],\
##                     self.FieldDict[fieldName]['var_type']
##          else: return None
            
