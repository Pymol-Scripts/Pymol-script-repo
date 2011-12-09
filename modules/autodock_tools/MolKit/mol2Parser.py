## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#############################################################################
#
# Author: Sophiec COON, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/MolKit/mol2Parser.py,v 1.14 2010/08/12 17:33:45 rhuey Exp $
#
# $Id: mol2Parser.py,v 1.14 2010/08/12 17:33:45 rhuey Exp $
#

import string, sys, types, os

from MolKit.protein import Protein, Chain, ChainSet, Residue, ResidueSet, ProteinSet
from MolKit.protein import Helix, Strand, Turn, Coil
from MolKit.molecule import Molecule, Atom, AtomSet, Bond, BondSet

from MolKit.moleculeParser import MoleculeParser

class Mol2Parser(MoleculeParser):

    Mol2Tags = ["@<TRIPOS>AlT_TYPE","@<TRIPOS>ANCHOR_ATOM",
                "@<TRIPOS>ASSOCIATED_ANNOTATION", "@<TRIPOS>ATOM",
                "@<TRIPOS>BOND","@<TRIPOS>CENTER_OF_MASS",
                "@<TRIPOS>CENTROID", "@<TRIPOS>COMMENT",
                "@<TRIPOS>CRYSIN", "@<TRIPOS>CURR_POS", "@<TRIPOS>DICT",
                "@<TRIPOS>DATA_FILE", "@<TRIPOS>EXTENSION_POINT",
                "@<TRIPOS>FF_PBC", "@<TRIPOS>FFCON_ANGLE",
                "@<TRIPOS>FFCON_DIST","@<TRIPOS>FFCON_RANGE",
                "@<TRIPOS>FFCON_TORSION", "@<TRIPOS>LINE",
                "@<tripos>LSPLANE",
                "@<TRIPOS>MOLECULE", "@<TRIPOS>NORMAL",
                "@<TRIPOS>POLYBUILD_HIST", "@<TRIPOS>QSAR_ALIGN_RULE",
                "@<TRIPOS>RING_CLOSURE", "@<TRIPOS>ROTABLE_BOND",
                "@<TRIPOS>SEARCH_DIST", "@<TRIPOS>SEARCH_OPTIONS",
                "@<TRIPOS>SUBSTRUCTURE", "@<TRIPOS>U_FEAT"]
    
               
    def __init__(self, filename):
        MoleculeParser.__init__(self, filename)
        self.mol2RecordParser = {}
        self.defaultReadOptions = ['@<TRIPOS>ATOM','@<TRIPOS>BOND',
                                   '@<TRIPOS>MOLECULE',
                                   '@<TRIPOS>SET','@<TRIPOS>SUBSTRUCTURE',
                                   '@<TRIPOS>DICT'] 
        self.keysAndLinesIndices = {} # stores all Mol2 keys .
        self.counter = 0
        self.setsDatas = []
        #self.molList = []

    def getKeysAndLinesIndices(self):
        """ Function to build a dictionary where the keys will be the
        records name of the mol2 files (@<TRIPOS>ATOM, @<TRIPOS>BOND...) and
        the value will be the index of the starting line of that record.
        """
        #this removes all comment and blank lines to fix bug #846
        for i,line in enumerate(self.allLines):
            if not line:
                self.allLines.pop(i)
            elif line[0] == '#':
                self.allLines.pop(i)
                
        i = 0
        record = None
        while  i != len(self.allLines):
            if self.allLines[i][:9] == '@<TRIPOS>':
                if self.keysAndLinesIndices:
                    self.keysAndLinesIndices[record].append(i)
                record = string.strip(self.allLines[i])
                self.keysAndLinesIndices[record] = [i+1]
                i = i+1
            else:
                i = i+1
        if record:
            self.keysAndLinesIndices[record].append(i)
        else:
            print " the file %s doesn't contain any mol2 records"%self.filename
            

        
    def parse(self):
        """ This function read a file and create the corresponding
        data hierarchy. """
        self.readFile()
        #molList = []
        molList = ProteinSet()
        if self.allLines is None:
            return
        elif len(self.allLines)!=0:
            self.getKeysAndLinesIndices()
        else:
            print "The file %s is empty"%self.filename
            return molList

        if not self.keysAndLinesIndices.has_key("@<TRIPOS>ATOM"):
            print "The file %s doesn't have Atom records, molecules can't be built"%self.filename
            return molList
        if self.keysAndLinesIndices.has_key('@<TRIPOS>SUBSTRUCTURE'):
            
            self.parse_MOL2_Substructure(self.allLines
                                         [self.keysAndLinesIndices
                                          ['@<TRIPOS>SUBSTRUCTURE'][0]:
                                          self.keysAndLinesIndices
                                          ['@<TRIPOS>SUBSTRUCTURE'][1]])
            molList.append(self.mol)

        else:
            atmlines = map(string.split, self.allLines
                           [self.keysAndLinesIndices
                            ['@<TRIPOS>ATOM'][0]:
                            self.keysAndLinesIndices
                            ['@<TRIPOS>ATOM'][1]])
            self.build4LevelsTree({},atmlines)
##              self.build2LevelsTree(map(string.split, self.allLines
##                                        [self.keysAndLinesIndices
##                                         ['@<TRIPOS>ATOM'][0]:
##                                         self.keysAndLinesIndices
##                                         ['@<TRIPOS>ATOM'][1]]))
            molList.append(self.mol)

        if self.keysAndLinesIndices.has_key('@<TRIPOS>BOND'):
            self.parse_MOL2_Bonds(self.allLines
                                  [self.keysAndLinesIndices
                                   ['@<TRIPOS>BOND'][0]:
                                   self.keysAndLinesIndices['@<TRIPOS>BOND']
                                   [1]])

        if self.keysAndLinesIndices.has_key('@<TRIPOS>SET'):
            self.parse_MOL2_Sets(self.keysAndLinesIndices['@<TRIPOS>SET'])

        return molList

    def parse_MOL2_Substructure(self, substlines):
        """build a dictionary with the chain id as keys and the
        list of residues belonging to that chain as values. If
        the id of the chain is not here then the keys is '', if
        two residues with the same name but belonging to two different
        chains then the value corresponding to that key is a list of
        the chains ID."""
        atomlines = map(string.split, self.allLines
                                 [self.keysAndLinesIndices
                                  ['@<TRIPOS>ATOM'][0]:
                                  self.keysAndLinesIndices
                                  ['@<TRIPOS>ATOM'][1]])
        subst_chain = {}
        if len(substlines) == 0:
            # case 1: no substructures are defined --> 2 levels tree.
            #self.build2LevelsTree(atomlines)
            #subst_chain = {'t
            self.build4LevelsTree(subst_chain,atomlines)

##          else:
##              substlines = map(string.split, substlines)
##              lines = filter(lambda x: len(x)>5, substlines)
##              if lines == [] or (lines != [] and \
##                                 filter(lambda x: x[5] != '****', lines)==[]):
                
##                  #self.build3LevelsTree(atomlines)

##              else:
##                  # case 3: at least 1 substructure and 1 chain --> 4 levels tree.
##                  #subst_chain = {}
                
##                  for line in substlines:
##                      try:
##                          if line[1] in subst_chain.keys():
##                              subst_chain[line[1]] = [subst_chain[line[1]]]
##                              subst_chain[line[1]].append(line[5])
##                          else:
##                              subst_chain[line[1]] = line[5]
##                      except:
##                          if line[1] in subst_chain.keys():
##                              list(subst_chain[line[1]]).append('')
##                          else:
##                              subst_chain[line[1]] = ''
##                  self.subst_chain = subst_chain
##                  self.build4LevelsTree(subst_chain,atomlines)
        else:
            # case 3: at least 1 substructure and 1 chain --> 4 levels tree.
            #subst_chain = {}
            substlines = map(string.split, substlines)
            for line in substlines:
                if len(line)<6 or line[5] == '****':
                    continue
                else:
                    try:
                        if line[1] in subst_chain.keys():
                            subst_chain[line[1]] = [subst_chain[line[1]]]
                            subst_chain[line[1]].append(line[5])
                        else:
                            subst_chain[line[1]] = line[5]
                    except:
                        if line[1] in subst_chain.keys():
                            list(subst_chain[line[1]]).append('')
                        else:
                            subst_chain[line[1]] = ''
            self.subst_chain = subst_chain
            self.build4LevelsTree(subst_chain,atomlines)


    def build2LevelsTree (self, atomlines):
        """
        Function to build a two level tree. 
        """
        print 'try to build a 2 level tree'
        self.mol= Molecule()
        self.mol.allAtoms = AtomSet()
        self.mol.atmNum = {}
        self.mol.parser = self
        if self.mol.name == 'NoName':
            self.mol.name = os.path.basename(os.path.splitext
                                             (self.filename)[0])

        self.mol.children = AtomSet([])
        self.mol.childrenName = 'atoms'
        self.mol.childrenSetClass = AtomSet
        self.mol.elementType = Atom
        self.mol.levels = [Molecule, Atom]
        ##1/18:self.mol.levels = [Protein, Atom]
        for atmline in atomlines:
            atom = Atom(atmline[1], self.mol,
                        chemicalElement = string.split(atmline[5], '.')[0],
            top = self.mol)
            #atom.element = atmline[5][0]
            atom.element = atom.chemElem
            atom.number = int(atmline[0])
            self.mol.atmNum[atom.number] = atom
            atom._coords = [ [float(atmline[2]), float(atmline[3]),
                                  float(atmline[4]) ] ]
            if len(atmline)>=9:
                atom._charges['mol2'] = float(atmline[8])
                atom.chargeSet = 'mol2'
#            atom.conformation = 0
            atom.hetatm = 0
            #add altname so buildBondsByDist doesn't croak
            atom.altname = None
            self.mol.allAtoms.append(atom)
        self.mol.atoms = self.mol.children


    def build3LevelsTree(self,atomlines):
        """ Function to build a 3 levels hierarchy Molecule-substructure-atoms."""

        self.mol= Protein()
        self.mol.allAtoms = AtomSet()
        self.mol.atmNum = {}
        self.mol.parser = self
        if self.mol.name == 'NoName':
            self.mol.name = os.path.basename(os.path.splitext
                                             (self.filename)[0])
        self.mol.children = ResidueSet([])
        self.mol.childrenName = 'residues'
        self.mol.childrenSetClass = ResidueSet
        self.mol.elementType = Residue
        self.mol.curRes = Residue()
        self.mol.curRes.hasCA = 0
        self.mol.curRes.hasO = 0
        
        self.mol.levels = [Protein, Residue, Atom]
        for atmline in atomlines:
            if len(atmline)>= 10:
                status = string.split(atmline[9], '|')
            else:
                status = None
            resName = atmline[7][:3]
            resSeq = atmline[7][3:]
            if resSeq != self.mol.curRes.number or \
               resName != self.mol.curRes.type:
                # check if this residue already exists
                na = string.strip(resName) + string.strip(resSeq)
                res = self.mol.get(na)
                if res:
                    self.mol.curRes = res[0]
                else:
                    self.mol.curRes = Residue(resName, resSeq, '',
                                              self.mol,
                                              top = self.mol)
            name = atmline[1]
            if name == 'CA': self.mol.curRes.hasCA = 1
            if name == 'O' : self.mol.curRes.hasO = 2
            atom = Atom(name, self.mol.curRes, top = self.mol,
            chemicalElement = string.split(atmline[5], '.')[0])
            #atom.element = atmline[5][0]
            atom.element = atom.chemElem
            atom.number = int(atmline[0])
            self.mol.atmNum[atom.number] = atom
            atom._coords = [ [float(atmline[2]), float(atmline[3]),
                              float(atmline[4]) ] ]
            atom._charges['mol2'] = float(atmline[8])
            atom.chargeSet = mol2
#            atom.conformation = 0
            atom.hetatm = 0
            #Add a data member containing a list of string describing
            # the Sybyl status bis of the atoms.
            atom.status = status 
            #add altname so buildBondsByDist doesn't croak
            atom.altname = None
            self.mol.allAtoms.append(atom)
            
        self.mol.residues = self.mol.children
        assert hasattr(self.mol, 'chains')
        delattr(self.mol, 'chains')
        delattr(self.mol, 'curRes')

    def build4LevelsTree(self, subst_chain, atomlines):
        """
        Function to build a 4 level hierarchy Protein-Chain-Residue-Atom.

        """
        self.mol= Protein()
        self.mol.allAtoms = AtomSet()
        self.mol.atmNum = {}
        self.mol.parser = self
        if self.mol.name == 'NoName':
            self.mol.name = os.path.basename(os.path.splitext
                                             (self.filename)[0])
        self.mol.curChain = Chain()
        self.mol.curRes = Residue()
        self.mol.levels = [Protein, Chain, Residue, Atom]
        i = 1
        for atmline in atomlines:
            if len(atmline)>= 10:
                status = string.split(atmline[9], '|')
            else: status = None
            if len(atmline) == 8:
                tmp = [atmline[5][:5], atmline[5][5:]]
                atmline[5] = tmp[0]
                atmline.insert(6, tmp[1])

            if status and status[0]=='WATER':
                chainID = 'W'
                atmline[7] = 'HOH'+str(i)
                subst_chain[atmline[7]] = chainID
                i = i+1

            if subst_chain == {}:
                chainID = 'default'

            elif not subst_chain.has_key(atmline[7]):
                if subst_chain.has_key('****'):
                    try:
                        chainID = subst_chain[atmline[7]]
                    except:
                        chainID = 'default'
                else:
                    chainID = 'default'

            elif type(subst_chain[atmline[7]]) is types.StringType:
                # that is to say that only chains has this substructure name.
                chainID = subst_chain[atmline[7]]

            elif type(subst_chain[atmline[7]]) is types.ListType:
                # That is to say that several chains have the same substructure.
                 chainID = subst_chain[atmline[7]][0]
                 subst_chain[atmline[7]] = subst_chain[atmline[7]].remove(chainID)
                 
            if chainID != self.mol.curChain.id:
                if not self.mol.chains.id or not chainID in self.mol.chains.id:
                    self.mol.curChain = Chain(chainID, self.mol,
                                          top = self.mol)
                else:
                    self.mol.curChain = self.mol.chains.get(chainID)[0]

            if len(atmline)<7:
                # test if the atmline has a res name and resseq:
                resName = 'RES'
                resSeq = '1'
            else:
                resName = atmline[7][:3]
                resSeq = atmline[7][3:]

            if resSeq != self.mol.curRes.number or \
               resName != self.mol.curRes.type:
                # check if this residue already exists
                na = string.strip(resName) + string.strip(resSeq)
                res = self.mol.curChain.get( na )
                if res:
                    self.mol.curRes = res[0]
                else:
                    self.mol.curRes = Residue(resName, resSeq, '',
                                              self.mol.curChain,
                                              top = self.mol)
            name = atmline[1]
            if name == 'CA': self.mol.curRes.hasCA = 1
            if name == 'O' : self.mol.curRes.hasO = 2
            atom = Atom(name, self.mol.curRes, top = self.mol,
            chemicalElement = string.split(atmline[5], '.')[0])
            #atom.element = atmline[5][0]
            atom.element = atom.chemElem
            atom.number = int(atmline[0])
            self.mol.atmNum[atom.number] = atom
            atom._coords = [ [float(atmline[2]), float(atmline[3]),
                              float(atmline[4]) ] ]
            if len(atmline)>=9:                
                atom._charges['mol2'] = float(atmline[8])
                atom.chargeSet = 'mol2'
#            atom.conformation = 0
            atom.hetatm = 0
            #Add a data member containing a list of string describing
            # the Sybyl status bis of the atoms.
            atom.status = status
            #add altname so buildBondsByDist doesn't croak
            atom.altname = None
            self.mol.allAtoms.append(atom)
        delattr(self.mol, 'curRes')
        delattr(self.mol, 'curChain')
        
    def parse_MOL2_Molecule(self, mollines):
        """Function to parse the Molecule records"""
        mollines = map(string.split, mollines)
        return mollines

    def parse_MOL2_Bonds(self, bondlines):
        """ Function to build the bonds object using the bond record of
        the mol2 file."""
        bondlines = map(string.split, bondlines)
        for bd in bondlines:
            at1 = self.mol.atmNum[int(bd[1])]
            at2 = self.mol.atmNum[int(bd[2])]

            if at1.isBonded(at2): continue
            bond = Bond(at1, at2, check=0)
            bond.type = bd[3]
            try:
                bond.bondOrder = int(bd[3])
            except:
                if bd[3]=='ar':
                    bond.bondOrder = 'aromatic'
                elif bd[3]=='am':
                    bond.bondOrder = 'amide'
                else:
                    bond.bondOrder = bd[3]
        self.mol.bondsflag = 1
        self.mol.hasBonds = 1
        
    def parse_MOL2_Sets(self, setRecords):
        """ Function to parse the Sets records"""
        setRecords = map(string.split, self.allLines[setRecords[0]:
                                                     setRecords[1]])
        i = 0
        while i!=len(setRecords):
            rec = []
            if len(setRecords[i]) <= 5:
                comments = None
                for j in xrange(len(setRecords[i])):
                    rec.append(setRecords[i][j])
                rec.append(comments)
            else :
                for j in xrange(len(setRecords[i][:5])):
                    rec.append(setRecords[i][j])
                
                comments = setRecords[i][5]
                for j in xrange(6, len(setRecords[i])):
                    comments = comments+' '+setRecords[i][j]
                rec.append(comments)
            number = []
            
            self.setsDatas.append(rec)

##              self.setsDatas.append([setRecords[i][0], setRecords[i][1],
##                                     setRecords[i][2], setRecords[i][3],
##                                     setRecords[i][4],comments])
            while len(setRecords[i+1])!=0 and setRecords[i+1][-1] == '\\':
                number = number+(map(lambda x: int(x), setRecords[i+1][:-1]))
                i = i+1

            number = number+map(lambda x: int(x),setRecords[i+1])
            self.setsDatas[-1].append(number)
            i = i+2

    def hasSsDataInFile(self):
        """ Function to extract the data on the secondarystructure and
        that replace the root atom number by the residue instance
        corresonding. """
        hData = filter(lambda x: x[0][:4] == 'HELI',self.setsDatas)
        sData = filter(lambda x: x[0][:4] == 'SHEE',self.setsDatas)
        tData = filter(lambda x: x[0][:4] == 'TURN',self.setsDatas)
        self.processSSEltData(sData, self.mol)
        self.processSSEltData(hData, self.mol )
        self.processSSEltData(tData, self.mol)
        self.ssData = [hData, sData, tData]

        if self.ssData == []:
            return 0
        else:
            return 1
        
    def parseSSData(self, mol):
        """
        Function to parse the info and return a list containing,
        the record name, and then the first and last residues for each
        secondary structure .
        """
        if not hasattr(self, 'ssData'):
            self.hasSsDataInFile()

        # Step 1: Create a list containing the information describing the
        # the secondary structures organized the following way:
        # [ ['chain1ID', [Helix, [startHel1, endHel1],[startHel2, endHel2]],
        # [Strand, [startSheet1, endSheet1]] ], ['chain2ID', [Helix .....]] ]
        ssDataForMol = {}
        for c in mol.chains:
            helStartEndForChain = self.processSSData(self.ssData[0], c)
            helStartEndForChain.insert(0, Helix)

            strandData = self.findStrands(self.ssData[1])
            strandStartEndForChain = self.processSSData(strandData, c)
            strandStartEndForChain.insert(0, Strand)
        
            turnStartEndForChain = self.processSSData(self.ssData[2], c)
            turnStartEndForChain.insert(0, Turn)

            ssDataForMol[c.id] = [ helStartEndForChain,strandStartEndForChain,
                                   turnStartEndForChain, None]

        return ssDataForMol

    def findStrands(self, data):
        """ Function to separate each strands of a sheet."""
        if len(data) == 0: return data
        else:
            for sheet in data:
                strandsBreak = []
                strandData = []
                for i in xrange(1,len(sheet[6])):
                    if i != 1 and \
                       int(sheet[6][i].number) - int(sheet[6][i-1].number)!=1:
                        strandsBreak.append(i)

                if len(strandsBreak) == 0:
                    strandData = sheet
                else:
                    i = 0
                    strandData.append(sheet[0],sheet[1],sheet[2],
                                       sheet[3],sheet[4],sheet[5],
                                       sheet[6][:strandsBreak[i]])
                    i = i+1
                    while i!= len(strandsBreak):
                        strandData.append(sheet[0],sheet[1],sheet[2],
                                       sheet[3],sheet[4],sheet[5],
                                       sheet[6][strandsBreak[i-1]:
                                                strandsBreak[i]])
                        i = i+1

                    strandData.append(sheet[0],sheet[1],sheet[2],
                                       sheet[3],sheet[4],sheet[5],
                                       sheet[6][strandsBreak[i-1]:])
            return strandData

    def processSSData(self, data, chain):
        """
        Function returning the information on the secondary structure of
        a given chain as a list which format is the following:
        - the first element of the list is the name of the secondary structure
        type ('Helix', 'Sheet', 'Turn')
        - the other are tuple containing the first residue of the structure,
        and the last one.
        This information is used by the class GetSecondarySTructureFromFile.
        """
        dataByChainID = filter(lambda x, id = chain.id:
                                x[-1][1].parent.id == id,
                                data)
        startEnd = map(lambda x: (x[-1][1],x[-1][-1]), dataByChainID)
    
        return startEnd

    def processSSEltData(self, ssData, mol):
        """
        Function to get the residue corresponding to the root atom number.
        """
        atoms = mol.chains.residues.atoms
        for data in ssData:
            for i in xrange(1,len(data[6])):
                if isinstance(data[6][i], types.IntType):
                    data[6][i] = atoms[data[6][i]-1].parent
                else:
                    return
            
    def getMoleculeInformation(self):
        """ Function to get the information on a molecule"""
        molStr = self.parse_MOL2_Molecule(self.allLines
                                          [self.keysAndLinesIndices
                                           ['@<TRIPOS>MOLECULE'][0]:
                                           self.keysAndLinesIndices
                                           ['@<TRIPOS>MOLECULE'][1]])
        chemical_formula = None
        if molStr != []:
            try:
                chemical_formula = molStr[-1][0]
            except:
                pass
            molStr = molStr[0][0]
        else:
            molStr = ''
        if chemical_formula in ["USER_CHARGES","NO_CHARGES"]:
            return molStr
        elif chemical_formula is not None:
            return "%s %s" %(molStr, chemical_formula)
        return molStr

