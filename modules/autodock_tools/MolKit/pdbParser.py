#############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/MolKit/pdbParser.py,v 1.97.4.2 2009/08/20 23:34:37 rhuey Exp $
#
# $Id: pdbParser.py,v 1.97.4.2 2009/08/20 23:34:37 rhuey Exp $
#

"""Module pdbParser.

Implements PdbParser and PQRParser. The objects are PDB file readers.
parsers can be registered fro every PDB card. By default, ATOM and HETATM
are parsed. Such objects are typically passed as an argument to the read
method of a molecule object. They build a 4 level tree mol:chain:residue:atom

Example:
    
    from MolKit.protein import Protein
    from MolKit.PdbParser import PdbParser
    mol = Protein()
    mol.read( '/tsri/pdb/struct/1crn.pdb', PdbParser() )

Atom parsing conventions:

    - we split the file into blocks delineated by MODEL/ENDMDL cards
    - for each block we split again into blocks starting at the first ATOM
      record and ending before the first ATOM record immediately follwing
      the first TER record, i.e. HETATMs placed behind a TER record will
      belong to the chain ended by the this TER
    - The order of ATOM and HETATM records is observed
    - Multiple models can create a chain for each model (default), or a unique
      chain with alternate conformations for each atom. The latter option
      requires that all models have the same number of atoms. This behavior is
      selected using the keyword 'modelsAs' which can assume the following
      values: 'molecules', 'chains', 'conformations'. The default is 'molecules'.
      be either 0 or 1. The options 'chains' is not yet implemented.
    - in PdbParser All atoms fields are parsed according to PDB specifications
    - in PQRParser we assume the following content for ATOM records:
         *ATOM num name resType resNum X Y Z charge radius*
"""

from os.path import splitext, basename
from string import split, strip, digits, lower, find
from MolKit.moleculeParser import MoleculeParser
from MolKit.protein import Protein, Chain, ChainSet, Residue, ResidueSet, ProteinSet
from MolKit.molecule import Atom, AtomSet, Bond, BondSet, HydrogenBond
from warnings import warn


class PdbParser(MoleculeParser):

    PDBtags = [
        'ANISOU', 'ATOM', 'AUTHOR', 'CAVEAT', 'CISPEP', 'COMPND', 'CONECT',
        'CRYST1', 'DBREF', 'END', 'ENDMDL', 'EXPDTA', 'FORMUL', 'HEADER',
        'HELIX', 'HET', 'HETATM', 'HETNAM', 'HETSYN', 'HYDBND', 'JRNL',
        'KEYWDS', 'LINK', 'MASTER', 'MODEL', 'MODRES', 'MTRIX1', 'MTRIX2',
        'MTRIX3', 'OBSLTE', 'ORIGX1', 'ORIGX2', 'ORIGX3', 'REMARK',
        'REVDAT', 'SCALE1', 'SCALE2', 'SCALE3', 'SEQADV', 'SEQRES',
        'SHEET', 'SIGATM', 'SIGUIJ', 'SITE', 'SLTBRG', 'SOURCE',
        'SPRSDE', 'SSBOND', 'TER', 'TITLE', 'TURN', 'TVECT'
        ]

##      def __del__(self):
##          print 'FreeParser'


    def __init__(self, filename=None, allLines=None, modelsAs='molecules'):
        """
        build a PDB parser object

        obj <- PDBParse(filename=None, allLines=None, modelsAs='molecules')

        filename: PDB file to be parsed
        allLines: ??? # FIXME
        modelsAs: can be set to 'molecules', 'chains', 'conformations'
                  'chains' is not yet implemented (7/09)
        """
        MoleculeParser.__init__(self, filename, allLines)
        self.pdbRecordParser = {}
        self.defaultReadOptions()
        self.keys = [] # stores all PDB keys from input file
        self.altLoc = None #Flag to handle alternate locations.
        self.model = False #Flag to indicate if Models or not
        self.useAbove73 = 0
        self.recordInfo = {} # stores the info when a record is parsed
        self.modelsAs = modelsAs
        
    def getKeys(self):
        if self.allLines:
            #l = map(split, self.allLines)
            #self.keys = map(lambda x: x[0], l)
            self.keys = map(lambda x: strip(x[:6]), self.allLines)
                
            

    def indexOfFirst(self, name, beginAt=0, endAt=None):
        """return index of first pdb card NAME starting at bginAt.
        self.allLines[index] is the card
        """
        if endAt is None: endAt=len(self.allLines)
        if name in self.keys[beginAt:endAt]:
            n = self.keys[beginAt:].index(name)
            return beginAt+n
        else: return -1

    
    def findBegEndFlagModels(self):
        """returns a list of ENDMDL cards indices"""
        mod = []
        end = self.indexOfFirst('ENDMDL')
        beg = self.indexOfFirst('MODEL')
        while end != -1 and beg != -1:
            flag = self.allLines[beg][11:15]
            mod.append( [beg,end,flag] )
            end = self.indexOfFirst('ENDMDL', end+1 )
            beg = self.indexOfFirst('MODEL', beg+1)
        return mod


    def findEndChains(self, beg=0, end=None):
        """Find what part of the file will be parsed as a chain.
        We search for the TER record and then go to the first ATOM.
        This means that HETATMS following will be in the block of atoms
        preceding the TER. If they have different chainID they will still
        create their own chain.
        """

        if end is None:
            end = len(self.allLines)
        ter = []
        t = self.indexOfFirst('TER', beg, end)
        while t != -1:
            trec = self.allLines[t]
            if trec.strip()=='TER':
                chainID = self.allLines[t-1][21]
            else:
                if len(trec) > 21:
                    chainID = trec[21]
                else: 
                    chainID = ' '
            ter.append((t, chainID))
            beg = t+1
            t = self.indexOfFirst('TER', beg)
            
        if t == -1:
            a = self.indexOfFirst('ATOM', beg)
            h = self.indexOfFirst('HETATM', beg)
            if a != -1 or h!=-1:
                ter.append((end, None))
#
        if len(ter)==0:
            ter.append((end,None))
        return ter

        
    def getRecords(self, reclist, key):
        """get all records starting with a given PDB tag"""
        return filter(lambda x: x[:len(key)]==key, reclist)


    def checkForRemark4(self):
        """check wether the file contains a 'REMARK   4' line. If we find it we
        will try to interpret colums 73-80, else we will ignore them
        """
        rem4 = filter(lambda x: x[:10] == 'REMARK   4', self.allLines)
        rem4c = filter(lambda x: x[16:24]=='COMPLIES', rem4)
        self.useAbove73 = len(rem4c)


    def defaultReadOptions(self):
        """define which PDB records get parsed and which function parse a
        specific record. By default only ATOM records are read"""
        self.pdbRecordParser['ATOM'] = self.parse_PDB_atoms
        self.pdbRecordParser['CONECT']= self.parse_PDB_CONECT
        #ADDED 8/10
        self.pdbRecordParser['HYDBND'] = self.parse_PDB_HYDBND


    def SetReadOptions(self, key, parser='default'):
        """Customize PDB reader. Arguments are a 'key' and an optional parser
The key has to be in the list of PDB tags see:
(http://www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html)
The optional parser argument can be:
    None: to prevent a record from being parsed
    'default': to use the default parser for that record. (default values)
    an executable object.

NOTE: The list currently registered parsers is in
      self.pdbRecordParser[key] = parser
"""

        assert key in self.PDBtags
        if callable(parser):
            self.pdbRecordParser[key] = parser
        elif parser is None:
            if key in self.pdbRecordParser.keys():
                del self.pdbRecordParser[key]
        else:
            if key   == 'ATOM'  : parser = self.parse_PDB_atoms
            elif key == 'HETATM': parser = self.parse_PDB_atoms
            elif key == 'CRYST1': parser = self.parse_PDB_CRYST1
            elif key == 'CONECT': parser = self.parse_PDB_CONECT
            elif key == 'SSBOND': parser = self.parse_PDB_SSBOND
            elif key == 'HELIX' : parser = self.parse_PDB_Helix
            elif key == 'TURN'  : parser = self.parse_PDB_Turn
            elif key == 'SHEET' : parser = self.parse_PDB_Strand
            elif key == 'HYDBND': parser = self.parse_PDB_HYDBND
            else:
                raise ValueError('No default parser for %s record, \
You have to provide one as the second argument' % key);
            self.pdbRecordParser[key] = parser


    def parse_PDB_CRYST1(self, lines):
        """Parse the CRYST1 record. Create the following members:
        cellLength, cellAngles, spaceGroup, Zvalue"""
        
        if len(lines)==0:
            raise RuntimeError( 'No CRYST1 record available in %s' %
                  self.filename )
        a = float(lines[0][6:15])
        b = float(lines[0][15:24])
        c = float(lines[0][24:33])
        alpha = float(lines[0][33:40])
        beta = float(lines[0][40:47])
        gamma = float(lines[0][47:54])
        self.mol.cellLength = [ a, b, c ]
        self.mol.cellAngles = [ alpha, beta, gamma ]
        self.mol.spaceGroup = strip(lines[0][55:66])
        try:
            self.mol.Zvalue = int(lines[0][66-70])
        except:
            self.mol.Zvalue = None
            

    def parse_PDB_CONECT(self, lines):
        """Parse the CONECT record. Create the bonds described in
        that record."""
        
        if len(lines)==0:return
        for c in lines:
            indices = [c[6:11], c[11:16], c[16:21], c[21:26], c[26:31],
                       c[31:36],c[36:41], c[41:46], c[46:51], c[51:56],
                       c[56:61]]
            indices = map(strip, indices)
            if self.mol.atmNum.has_key(int(indices[0])):
                a1 = self.mol.atmNum[int(indices[0])]
            else:
                continue

            for i in xrange(len(indices)-1):
                if indices[i+1] in [ '     ', '\n', '', '   \012']:
                    break

                if self.mol.atmNum.has_key(int(indices[i+1])):
                    a2 = self.mol.atmNum[int(indices[i+1])]
                else:
                    continue

                atms = []
                for b in a1.bonds:
                    atms.append(b.atom1)
                    atms.append(b.atom2)

                if not a2 in atms  and not a1.isBonded(a2):
                    # check if the bond a1-a2 has not been created yet.
                    bond = Bond(a1, a2, check=0)


    def parse_PDB_HYDBND(self, lines):
        """Parse the HYDBND record. Create the hbond described in
        that record by finding dAt, hAt and aAt, the donor, hatom and
        acceptorAtoms respectively."""
        if len(lines)==0:return
        #6:11 corresponds to documentation's 7-11
        for c in lines:
            dAtName = strip(c[12:16])
            dAtPType = c[17:20]
            dAtPNum = int(c[22:27])
            dAtPName = dAtPType + str(dAtPNum)
            dAtPIcode = c[27]
            if not dAtPIcode==' ':
                dAtPName = dAtPName + dAtPIcode
            dAtChId = c[21]
            dname = self.mol.name+':'+dAtChId+':'+dAtPName+':'+dAtName

            hAtName = strip(c[29:33])
            if len(hAtName):
                if len(hAtName)==4:
                    hAtName = hAtName[1:]+hAtName[0]
                elif len(hAtName) > 1:
                    if hAtName[1]=='H':
                        if hAtName[0] in ('1','2','3'):
                            hAtName = strip(hAtName[1:])+hAtName[0]
                #construct the full name of this hydrogen
                hAtChId = c[33]
                hAtPNum = c[36:41]
                h_full_name = self.mol.name+':'+dAtChId+':'+dAtPName+':'+hAtName

            aAtName = strip(c[43:47])
            aAtPType = c[48:51]
            aAtPNum = int(c[53:58])
            aAtPName = aAtPType + str(aAtPNum)
            aAtPIcode = c[58]
            if not aAtPIcode == ' ':
                aAtPName = aAtPName + aAtPIcode
            aAtChId = c[52]
            aname = self.mol.name+':'+aAtChId+':'+aAtPName+':'+aAtName
            try:               
                dAt = self.mol.allAtoms.get(lambda x: x.full_name()==dname)[0]
                aAt = self.mol.allAtoms.get(lambda x: x.full_name()==aname)[0]
                #dAt = self.mol.NodesFromName(dname)[0]
                #aAt = self.mol.NodesFromName(aname)[0]
                # hydrogen atom is optional
                if len(hAtName):
                    hAt = self.mol.allAtoms.get(lambda x: x.full_name()==h_full_name)[0]
                    #hAt = self.mol.NodesFromName(hname)[0]
                else:
                    hAt = None
                hbond = HydrogenBond(dAt, aAt, hAt, check=0)
                for item in [dAt, aAt]:
                    if not hasattr(item, 'hbonds'):
                        item.hbonds = [hbond]
                    else:
                        item.hbonds.append(hbond)
                if hAt is not None:
                    hAt.hbonds = [hbond]
            except:
                import sys
                print >>sys.stderr,"Unable to parse Hydrogen Bond Record in",\
                self.filename
                print >>sys.stderr,c

    def parse_PDB_atoms(self, atoms, mol):
        """find and parse the PDB compliant ATOM.
        The tree structure for chains residues and atoms is built here"""
        mol.curChain = Chain()
        mol.curRes = Residue()
        mol.levels = [Protein, Chain, Residue, Atom]
        
        # create a allAtom set which contains all the atoms of the molecule
        ats = []
        for atRec in atoms:
            if atRec[:3]=='TER':
                self.curChain = Chain()
                continue

            if atRec[:4]!='ATOM' and atRec[:6]!='HETATM':
                continue

            if atRec.startswith('ATOM   5050'):
                pass
            try:
                at = self.parse_PDB_ATOM_record(atRec, mol)
                if at is not None:
                    ats.append(at)
            except Exception, inst:
                print inst
                print "Error parsing the following line in pdb:\n"+atRec

        mol.allAtoms = mol.allAtoms + AtomSet(ats)
        
        #map( self.parse_PDB_ATOM_record,atm, [modlflag,]*len(atm) )
        #atoms = filter(lambda x, cid=atoms[0][21]: x[21]!=cid, atoms)
        #print len(atoms)
        #if len(atoms)>0: self.parse_PDB_atoms(atoms, mol)
    

    def parse_PDB_ATOM_record(self, rec, mol):
        """Parse PDB ATOM records using the pdb columns specifications"""
        # Handle the alternate location using a flag.
        rec = strip(rec)

        if rec[16]!= ' ': self.altLoc = rec[16]
        else: self.altLoc = None
        
        # check for chains break
        chainID = rec[21]

        if chainID != mol.curChain.id:
            # If not creating a new chain when new ID
            # 1- Try to get the old chain with same ID
            if not mol.chains.get(lambda x: x.id == chainID):
                # create a new chain
                mol.curChain = Chain(chainID, mol, top=mol)
            else:
                # reuse the first chain with the same ID
                chains = mol.chains.get(lambda x: x.id == chainID)
                mol.curChain = chains[0]
        
        # check for residue break
        resName = rec[17:20]
        resSeq = strip(rec[22:26]) #WARNING reSeq is a STRING
        if rec[26] != ' ': icode = rec[26]
        else: icode = ''

        curRes = mol.curRes
        if curRes.parent is None or curRes.parent.name!=chainID or \
                resSeq != curRes.number or  resName != curRes.type or \
                icode != curRes.icode:
            #if resSeq != curRes.number or resName != curRes.type or \
            #   icode != curRes.icode:
            # check if this residue already exists
            na = strip(resName) + resSeq + icode
            try:
                mol.curRes = mol.curChain.childByName[na]
            except KeyError: # child not found
                # create a  new residues
                mol.curRes = Residue(resName, resSeq, icode,
                                     mol.curChain, top=mol)
                # Add a hasCA and hasO flags to the residue and set it to 0
                mol.curRes.hasCA = 0
                mol.curRes.hasO = 0
        # parse atom info

        # handle atom names (calcium, hydrogen) and find element type
        # check validity of chemical element column and charge column
        name, element, charge = self.getPDBAtomName(rec[12:16], rec[76:78],
                                                    rec[78:80])
        
        if name == 'CA':
            # Set the flag hasCA to 1 if the atom name is CA
            mol.curRes.hasCA = 1

        if name == 'O' or name == 'OXT' or (len(name)>3 and name[:3]=='OCT'):
            # Set the hasO flag to 2 if atom name is O or OXT
            mol.curRes.hasO = 2

        # WHY lower ??? why not test elem=='L'
        elem = lower(element)
        if elem =='l':
            element = 'Xx'
        atom = Atom(name, mol.curRes, element, top=mol)
        atom._coords = [ [ float(rec[30:38]), float(rec[38:46]),
                           float(rec[46:54]) ] ]
        atom.segID = strip(rec[72:76])
        if rec[:4]=='ATOM': atom.hetatm = 0
        else: atom.hetatm = 1
        #atom.alternate = []
        #atom.element = element  # done in Atom constructor
        atom.normalname = name
        atom.number = int(rec[6:11])
        mol.atmNum[atom.number] = atom
        if rec[54:60] == "      " or len(rec[54:60]) == 0:
            atom.occupancy = 0.0
        else:
            atom.occupancy = float(rec[54:60])
        # conformation is now assigned in constructor
        #        atom.conformation = 0
        if len(rec) > 61:
            if rec[60:66]=="      " or len(rec[60:66]) == 0:
                atom.temperatureFactor = 0.0
            else:
                atom.temperatureFactor=float(rec[60:66])
        else: atom.temperatureFactor= 0.0
        
        atom.altname = None
        if self.altLoc :
            # check if the name of the atom is the same than the
            #name of the previous atom .
            atom.normalname = name
            name = name + '@'+self.altLoc
            atom.name = name
            atom.altname = self.altLoc
            if len(mol.curRes.atoms)>1:
                # the new atom has been added to the current residue
                # You have to go to the one before.
                lastAtom = mol.curRes.atoms[-2]
                if atom.normalname == lastAtom.normalname:
                    # Add the new alternate atom to the LastAtom.alternate and
                    # add the lastAtom to the atom.alternate.
                    lastAtom.alternate.append(atom)
                    atom.alternate.append(lastAtom)
                    for l in lastAtom.alternate:
                        if atom.name != l.name:
                            atom.alternate.append(l)
                            l.alternate.append(atom)

        # call the progress bar
        self.updateProgressBar()
                
        return atom


    def getPDBAtomName(self, name, element, charge):
        """Figure out the real atom name as well as the atom element"""
        # To be DONE
        orig_element = element.strip()
        if not self.useAbove73 or element == '  ' or element == '':
            
            element = strip(name[0:2])
            if not element:     
                warn("Chemical element type is missing for %s!, Please corrent the pdb file"%name)
                
        if element and element[0] in digits:
            if orig_element:
                element = orig_element
            else:
                element = name[1]
            charge = None        
                
        if element and len(element) > 1:
            if element[1] in digits:
                if orig_element:
                    element = orig_element
                else:
                    element = name[0]
                charge = None                    
        if element:
            element = strip(element)
            
        
        if len(name)>1 and name[1]=='H':
            if name[0] in ('1','2','3'):
                name = strip(name[1:])+name[0]
                element = 'H'
        #old style names: 'HE12' in columns 12-15
        if name[0]=='H' and len(name)==4 and name[2] in ('1','2','3') \
                and name[3] in ('1','2','3'):
            name = name 
            element = 'H'
        if len(element)==2:
            element = element[0]+lower(element[1])
        return strip(name), element, charge


    def configureProgressBar(self, **kw):
        # this method is to be implemented by the user from outside
        pass


    def updateProgressBar(self, progress=None):
        # this method is to be implemented by the user from outside
        #print 'Progess: ' + `progress` + '%'
        pass
    

    def addModel(self, atoms):
        """Add a conformation to the first chains of the molecule"""
        l = lambda x: [float(x[30:38]), float(x[38:46]), float(x[46:54])]
        coords = map(l, atoms)
        map( Atom.addConformation, self.mol.chains[0].residues.atoms, coords )
        
        
    def parse(self, objClass=Protein):
        """Reads a PDB that is PDB specs compliants"""

        if self.allLines is None and self.filename:
            self.readFile()
            if self.allLines is None or len(self.allLines)==0: return

        self.getKeys()
        # check if the file describe a molecule in a pdb format
        pdbKeys = filter(lambda x: x in self.PDBtags, self.keys)
        if len(pdbKeys)==0:
            return 

        self.checkForRemark4()
        totalAtomNumber = len(self.getAtomsLines(-2, 0))
        # configure the progress bar, initialize it, set mode to 'increment'
        self.configureProgressBar(init=1, mode='increment',
                                  labeltext='parse atoms', max=totalAtomNumber)
        
        if self.pdbRecordParser.has_key('ATOM'):
            # deal with ATOM first to make sure self.mol,
            # an objClass instance, is created
            record='ATOM'
            modelEnd = self.findBegEndFlagModels()
            if not modelEnd: # No MODEL records.
                self.mol = objClass()
                self.mol.allAtoms = AtomSet([])
                self.mol.parser = self
                if self.mol.name == 'NoName':
                    if not self.filename is None:
                        self.mol.name = basename(splitext(self.filename)[0])

                molList = self.mol.setClass()
                molList.append(self.mol)

                self.pdbRecordParser[record](self.allLines, self.mol)

## MS Not sure what .gap is used for ? might affect ribbons ?
##
##                 ter = self.findEndChains()
##                 beg = 0
## ##                 #if __debug__:
## ##                     #print 'End of chains lines: ',ter
##                 for t, chainID in ter:
##                     #if self.pdbRecordParser.has_key('ATOM'):
##                     # getAtomsLines already returns atom + hetatm
##                     atoms = self.getAtomsLines(t, beg)
##                     parseFunc = self.pdbRecordParser['ATOM']
                        
##                     if self.pdbRecordParser.has_key('HETATM'):
##                         # if 'ATOM' and 'HETATM' in pdbRecordParser.keys
##                         # redondant...
##                         atoms = atoms + self.getRecords(self.allLines,
##                                                         'HETATM')
##                         parseFunc = self.pdbRecordParser['HETATM']
                        
##                     if chainID is None:
##                         chainID = atoms[0][21]
                        
##                     if self.mol.chains.id and chainID in self.mol.chains.id:
##                         self.mol.curChain = self.mol.chains.get(lambda x: x.id==chainID)[0]
##                         gap=True
##                     else:
##                         self.mol.curChain = Chain(id = chainID,
##                                                   parent=self.mol,
##                                                   top=self.mol)
##                         gap=False

##                     if gap:
##                         curResName = self.mol.curRes.name
##                         rec = atoms[0]
##                         nextResName= rec[17:20] + strip(rec[22:26])
##                         self.mol.curChain.gaps.append((curResName,
##                                                        nextResName))

##                     parseFunc(atoms)
                    
##                     beg = t+1

            else: 
                self.model = True
                self.mol = objClass()
                #self.mol.model = ProteinSet()
                self.mol.parser = self
                molList = self.mol.setClass()
                name = basename(splitext(self.filename)[0])
                # create a molecule using the first model
                mBeg, mEnd, mNum = modelEnd[0]
                
                if mBeg:
                    #needed for Complex visualisation. See http://mgl.scripps.edu/forum/viewtopic.php?f=11&t=280
                    mol = objClass()
                    mol.allAtoms = AtomSet([])
                    mol.parser = self
                    if mol.name == 'NoName':
                        if not self.filename is None:
                            mol.name = basename(splitext(self.filename)[0])
                    self.pdbRecordParser[record](self.allLines[:mBeg], mol)
                    if len(mol.allAtoms):
                        molList.append(mol)
                                                    
                noNumbers = False
                if mNum=="": 
                    noNumbers = True
                    mNum="1"
                    mnum=1
                else:
                    mnum = int(mNum.strip())
                self.configureProgressBar(labeltext='parse atoms Model '+\
                                          mNum.strip()+'/%d'%len(modelEnd))
                self.mol = objClass()
                self.mol.parser = self
                if self.mol.name == 'NoName':
                    self.mol.name = name
                    if self.modelsAs == 'molecules':
                        self.mol.name = name+'_model%d'%(mnum)
                    #self.mol.name = name+'_model'+mNum.strip()
                    firstName = self.mol.name
                molList.append(self.mol)
                self.pdbRecordParser[record](self.allLines[mBeg+1:mEnd],
                                             self.mol)
                if  self.pdbRecordParser.has_key('CONECT'):
                    records = self.getRecords(self.allLines, 'CONECT' )
                    self.pdbRecordParser['CONECT'](records)    
                 
                if self.modelsAs=='molecules': # create a new molecule for each model
                    # loop over block of lines for each model
                    for mBeg, mEnd, mNum in modelEnd[1:]:
                        if mNum=="": mNum = str(mnum)
                        labeltext='parse atoms Model%d  of %d'%(mnum, len(modelEnd))
                        self.configureProgressBar(labeltext='parse atoms Model%d  of %d'%(mnum, len(modelEnd)))
                                                  #mNum.strip()+'/%d'%len(modelEnd))

                        newmol = objClass()
                        newmol.parser = self
                        if self.mol.name == firstName:
                            mnum = mnum + 1
                            newmol.name = name+'_model%d'%(mnum)
                        molList.append(newmol)
                        self.pdbRecordParser[record](self.allLines[mBeg+1:mEnd],
                                                     newmol)
                        if  self.pdbRecordParser.has_key('CONECT'):
                            records = self.getRecords(self.allLines, 'CONECT' )
                            self.pdbRecordParser['CONECT'](records)    

                elif self.modelsAs=='chains': # create a new chain for each model
                    print 'NOT YET IMPLEMENTED'
                    raise RuntimeError

                elif self.modelsAs=='conformations': # add a conformation to molecule
                    atms = self.mol.allAtoms
                    for mBeg, mEnd, mNum in modelEnd[1:]:
                        self.configureProgressBar(
                            labeltext='parse atoms Model '+\
                            mNum.strip()+'/%d'%len(modelEnd))

                        mol = objClass()
                        self.pdbRecordParser[record](self.allLines[mBeg+1:mEnd], mol)
                        atms.addConformation(mol.allAtoms.coords)
                        
                else:
                    raise ValueError("bad value modelsAs, expected 'molecules', 'chains' or 'conformations', got" + self.modelsAs)
                
##                 ter = self.findEndChains(modelEnd[0][0],
##                                          modelEnd[0][1])
##                 begChain = 0
##                 assert (self.mol.setClass == molList.__class__)
##                 molList.append(self.mol)
##                 for t, chainID in ter:
##                     atoms = self.getAtomsLines(t, begChain)
##                     self.pdbRecordParser[record](atoms)
##                     begChain = t+1
##                 nbAtoms = len(atoms)
##                 print 'MODEL1', nbAtoms
##                 raise
##                 for i in range(1,len(modelEnd)):
##                     # Create a new molecule for each model !
##                     self.mol = objClass()
##                     self.mol.parser = self
##                     if self.mol.name == 'NoName':
##                         self.mol.name = name+modelEnd[i][2]
##                     print 'Model', i 
##                     assert (molList.__class__ == self.mol.setClass)
##                     molList.append(self.mol)
##                     ter = self.findEndChains(modelEnd[i][0],
##                                              modelEnd[i][1])
##                     begChain = modelEnd[i][0]
##                     for t, chainID in ter:
##                         atoms = self.getAtomsLines(t, begChain)
##                         print 'FAGA', t, begChain, len(atoms)
##                         self.pdbRecordParser[record](atoms)
##                         #self.addModelToMolecules(molList)
##                         begChain = t+1
            delattr(self.mol, 'curChain')
            delattr(self.mol,'curRes')
            
        else:
            self.mol = objClass()
            molLis = self.mol.setClass([self.mol])

        # MS commented out because progress should go over all atoms in MODELS
        # set progress bar mode to 'increment' and initialize it
        #lenRecs = len(self.pdbRecordParser.keys())
        #self.configureProgressBar(init=1, mode='increment', max=lenRecs)
        for record in self.pdbRecordParser.keys():
            #self.configureProgressBar(labeltext='parse '+record)
            #self.updateProgressBar()
            if record=='HETATM': continue # they get parsed with ATOM
            elif record=='ATOM': continue
            else:
##                 data = self.getParsedRecordData(record,
##                                                 self.pdbRecordParser[record])
                if record=='CONECT' and self.model==True: 
                    continue
                else:
                    records = self.getRecords(self.allLines, record )
                    self.pdbRecordParser[record](records)

        for mol in molList:
            self.set_stringReprs(mol)
        
        #also set the stringReprs for this MoleculeSet
        name = ''
        for n in molList.name:
            name = n + ','
        

            
        #remove the terminal comma
        name = name[:-1]
        molList.setStringRepr(name)
        strRpr = name + ':::'
        molList.allAtoms.setStringRepr(strRpr)
        
        #have to assign vina attributes here
        if hasattr(molList[0],'vina_results'):
            results = molList[0].vina_results[:]
            if len(molList)==len(results):
                for m, r in zip(molList, results):
                    m.vina_energy, m.vina_rmsd_lb, m.vina_rmsd_ub = map(float, r)
                    print m.name, ':',  m.vina_energy, ',', m.vina_rmsd_lb,',', m.vina_rmsd_ub
            elif len(results):
                m = molList[0]
                m.vina_results = results
                m.vina_energy, m.vina_rmsd_lb, m.vina_rmsd_ub = map(float, results[0])

        return molList


    def set_stringReprs(self, mol):
        # after the hierarchy has been built go set all stringRepr
        # their concise values
        mol.allAtoms = mol.chains.residues.atoms
        mname = mol.name
        strRpr = mname + ':::'
        mol.allAtoms.setStringRepr(strRpr)

        strRpr = mname + ':'
        mol.chains.setStringRepr(strRpr)

        for c in mol.chains:
            cname = c.id
            strRpr = mname + ':' + cname + ':'
            c.residues.setStringRepr(strRpr)
            for r in c.residues:
                rname = r.name
                strRpr = mname + ':' + cname + ':' + rname + ':'
                r.atoms.setStringRepr(strRpr)


    def addModelToMolecules(self, listOfMol):
        length = len(listOfMol)
        for i in xrange(length):
            listOfMol[i].model = ProteinSet()
            for j in xrange(length):
                if listOfMol[i]!= listOfMol[j]:
                    listOfMol[i].model.append(listOfMol[j])


    def getAtomsLines(self, ter, begChain):
        l = lambda x: x[:4]=='ATOM' or x[:6]=='HETATM'
        atoms = filter(l, self.allLines[begChain:ter+1])
        return atoms


    def parse_PDB_Helix(self, lines):
        """
        Parse the Helix records and stores the list of tuple containing all
        the datas in self.recordInfo['HELIX'].
        hdatas will contain the following:
        0 - HELIX
        1 - int serial number of the helix
        2 - Helix Identifier
        3 - Name of Nterminus Res (start)
        4 - Chain ID
        5 - Seq number of the NTerm Res
        6 - Insertion Code of the NTerm Res
        7 - Name of CTerm Residue (end)
        8 - Chain ID
        9 - Seq number of the CTerm Res
        10- Insertion Code of the CTerm Res
        11- Helix class (int)
        12- comment
        13- Length of Helix
        """
        hdatas = []
        for line in lines:
            try:
##                 hdatas.append( (line[0:6], int(line[7:10]), line[11:14],
##                                 line[15:18], line[19], line[21:26],
##                                 line[27:30], line[31], line[33:38],
##                                 int(line[38:40]), line[40:70],
##                                 line[71:76] ))
                hdatas.append( (line[0:6],        # RECORD
                                int(line[7:10]),  # SER NUM
                                line[11:14],      # HEL ID
                                line[15:18],      # RES NAME
                                line[19],         # CHAIN ID
                                line[21:25],      # RES SEQNUMBER
                                line[25],         # RES INSERTION
                                line[27:30],      # RES NAME
                                line[31],         # CHAIN ID
                                line[33:37],      # RES SEQNUMBER
                                line[37],         # RES INSERTION
                                int(line[38:40]), # HEL CLASS
                                line[40:70],      # COMMENT
                                line[71:76]       # HEL LEN
                                ))
            except ValueError:
                hdatas.append( (line[0:6],        # RECORD
                                int(line[7:10]),  # SER NUM
                                line[11:14],      # HEL ID
                                line[15:18],      # RES NAME
                                line[19],         # CHAIN ID
                                line[21:25],      # RES SEQNUMBER
                                line[25],         # RES INSERTION
                                line[27:30],      # RES NAME
                                line[31],         # CHAIN ID
                                line[33:37],      # RES SEQNUMBER
                                line[37],         # RES INSERTION
                                int(line[38:40]), # HEL CLASS
                                line[40:70],      # COMMENT
                                line[71:76]       # HEL LEN
                                ))
        
        if not self.recordInfo.has_key('HELIX'):
            self.recordInfo['HELIX']=hdatas
        else:
            self.recordInfo['HELIX']=hdatas
        #return hdatas

    def parse_PDB_Strand(self, lines):
        """
        Parse the Strand records and stores the list of tuple containing
        all the datas in self.recordInfo['SHEET'].
        sdatas will contain the following information
        
        """
        sdatas = []
        for line in lines:
            try:
                sdatas.append( ( line[0:6],       # 0 RECORD
                                 int(line[7:10]), # 1 SER NUM
                                 line[11:14],     # 2 SHEET ID
                                 int(line[14:16]),# 3 STR NB 
                                 line[17:20],     # 4 RES NAME
                                 line[21],        # 5 CHAIN ID
                                 line[22:26],     # 6 RES SEQNUMBER
                                 line[27],        # 7 RES INSERTION
                                 line[28:31],     # 8 RES NAME
                                 line[32],        # 9 CHAIN ID
                                 line[33:37],     # 10 RES SEQNUMBER
                                 line[37],        # 11 RES INSERTION
                                 int(line[38:40]),# 12 STR SENSE
                                 line[41:45],     
                                 line[45:48],
                                 line[49],
                                 line[50:55],
                                 line[56:60],
                                 line[60:63],
                                 line[64],
                                 line[65:69],
                                 line[69]))
            except ValueError:
                sdatas.append( ( line[0:6],       # RECORD
                                 int(line[7:10]), # SER NUM
                                 line[11:14],     # SHEET ID
                                 int(line[14:16]),# STR NB 
                                 line[17:20],     # RES NAME
                                 line[21],        # CHAIN ID
                                 line[22:26],     # RES SEQNUMBER
                                 line[27],        # RES INSERTION
                                 line[28:31],     # RES NAME
                                 line[32],        # CHAIN ID
                                 line[33:37],     # RES SEQNUMBER
                                 line[37],        # RES INSERTION
                                 int(line[38:40]),# STR SENSE
                                 line[41:45],     
                                 line[45:48],
                                 line[49],
                                 line[50:55],
                                 line[56:60],
                                 line[60:63],
                                 line[64],
                                 line[65:69],
                                 line[69]))


        if not self.recordInfo.has_key('SHEET'):
            self.recordInfo['SHEET']=sdatas
        else:
            self.recordInfo['SHEET']=sdatas


    def parse_PDB_Turn(self, lines):
        """
        Parse the Turn records and stores the list of tuple containing
        all the datas in the self.recordInfo['TURN'].
        """
        tdatas = []
        for line in lines:
            tdatas.append(( line[0:6],       # 0 RECORD
                            int(line[7:10]), # 1 SER NUM
                            line[11:14],     # 2 TURN ID
                            line[15:18],     # 3 RES NAME
                            line[19],        # 4 CHAIN ID
                            line[20:24],     # 5 RES SEQNUMBER
                            line[24],        # 6 RES INSERTION
                            line[26:29],     # 7 RES NAME
                            line[30],        # 8 CHAIN ID
                            line[31:35],     # 9 RES SEQNUMBER
                            line[35],        # 10 RES INSERTION
                            line[40:70]      # 11 COMMENT
                            ))

        if not self.recordInfo.has_key('TURN'):
            self.recordInfo['TURN']=tdatas
        else:
            self.recordInfo['TURN']=tdatas


    def hasSsDataInFile(self):
        """ Function testing if the informations on the secondary structure
        are in the file"""
        test = filter(lambda x: x in self.keys,['HELIX','SHEET', 'TURN'])
        if test: return 1
        else: return 0 


    def parseSSData(self, mol):
        """
        Function to parse the information describing the secondary structure
        of the protein grouped as chain ID, the information is provided
        as a list of the following structure:
        [ ['chainID',[ Helix,[1stResHelix1,lastResHelix1], ...],
        [ Strand, [1stResSheet1,lastResSheet1] ],...],.... ]
        """
        from MolKit.protein import Helix, Strand, Turn

        # Step 1: Parse the information describing each type of secondary
        # structure for the whole molecule.
        if not self.recordInfo.has_key('HELIX'):
            hRecords = self.getRecords(self.allLines,'HELIX')
            self.parse_PDB_Helix(hRecords)
        helDataForMol = self.recordInfo['HELIX']
        
        if not self.recordInfo.has_key('SHEET'):
            sRecords = self.getRecords(self.allLines,'SHEET')
            self.parse_PDB_Strand(sRecords)
        strandDataForMol = self.recordInfo['SHEET']
        if not self.recordInfo.has_key('TURN'):
            tRecords = self.getRecords(self.allLines,'TURN')
            self.parse_PDB_Turn(tRecords)
        turnDataForMol=self.recordInfo['TURN']

        # Step 2: Create a list containing the information describing the
        # the secondary structures organized the following way:
        # [ ['chain1ID', [Helix, [startHel1, endHel1],[startHel2, endHel2]],
        # [Strand, [startSheet1, endSheet1]] ], ['chain2ID', [Helix .....]] ]
        ssDataForMol = {}
        for c in mol.chains:
            helStartEndForChain  = self.processHelData(helDataForMol,
                                                       (4,3,5,6,7,9,10,11,12),
                                                       c)
            # fieldIndices
            # 4: chain ID, 3: RESNAME, 5: RES SEQNB, 6: RES INSER,
            # 7: RESNAME, 9: RES SEQNB, 10: RES INSER,
            # 11: HEL CLASS, 12: COMMENT
            helStartEndForChain.insert(0, Helix)
            
            # fieldIndices
            # 5: chain ID, 4: RESNAME, 6: RES SEQNB, 7: RES INSER,
            # 8: RESNAME, 10: RES SEQNB, 11: RES INSER, 3: NB STRAND,
            # 12: STRAND SENSE
            strandStartEndForChain = self.processStrData(strandDataForMol,
                                                         (5,4,6,7,8,10,11,
                                                          3,12),
                                                         c)
            strandStartEndForChain.insert(0,Strand)

            # fieldIndices
            # 4: chain ID, 3: RESNAME, 5: RES SEQNB, 6: RES INSER,
            # 7: RESNAME, 9: RES SEQNB, 10: RES INSER,
            # 11: COMMENT
            turnStartEndForChain = self.processTurnData(turnDataForMol,
                                                        (4,3,5,6,7,9,10,11),
                                                        c)
            turnStartEndForChain.insert(0,Turn)
            
            ssDataForMol[c.id] = [helStartEndForChain,strandStartEndForChain,
                                   turnStartEndForChain, None]
            
        return ssDataForMol

    def processHelData(self, recordDataForMol, fieldIndices, chain):
        # fieldIndices[0] = chain ID
        recordDataForChain = filter(lambda x: x[fieldIndices[0]]==chain.id,
                                    recordDataForMol)
        helStartEndData = []
        for rec in recordDataForChain:
            # ResName+ResNumber+ResInsertion
            startData = rec[fieldIndices[1]]+rec[fieldIndices[2]].strip()+\
                        rec[fieldIndices[3]]
            startData = startData.strip()
            startRes = chain.residues.get(startData)
            if len(startRes)!= 1:
                raise ValueError("ERROR: When parsing the helix information %s \
was not found in chain %s"%(startData, chain.id))

            endData = rec[fieldIndices[4]]+rec[fieldIndices[5]].strip()+\
                      rec[fieldIndices[6]]
            endData = endData.strip()
            endRes = chain.residues.get(endData)
            if len(endRes)!= 1:
                raise ValueError("ERROR: When parsing the helix information %s \
was not found in chain %s"%(endData, chain.id))

            helClass = rec[fieldIndices[7]]
            comment = rec[fieldIndices[8]]
            helStartEndData.append({'start':startRes[0], 'end':endRes[0],
                                    'helClass':helClass, 'comment':comment})
        return helStartEndData
                                                 
    def processStrData(self, recordDataForMol, fieldIndices, chain):
        # Chain id
        recordDataForChain = filter(lambda x: x[fieldIndices[0]] == chain.id,
                                    recordDataForMol)
        strStartEndData = []
        for rec in recordDataForChain:
            # Res name + Res seqnumber + Res insertion
            startData = rec[fieldIndices[1]]+rec[fieldIndices[2]].strip()+\
                        rec[fieldIndices[3]]
            startData = startData.strip()
            startRes = chain.residues.get(startData)
            if len(startRes)!= 1:
                raise ValueError("ERROR: When parsing the helix information %s \
was not found in chain %s"%(startData, chain.id))
            
            # Res name + Res seqnumber + Res insertion
            endData = rec[fieldIndices[4]]+rec[fieldIndices[5]].strip()+\
                      rec[fieldIndices[6]]
            endData = endData.strip()
            endRes = chain.residues.get(endData)
            if len(endRes)!= 1:
                raise ValueError("ERROR: When parsing the helix information %s \
was not found in chain %s"%(endData, chain.id))
            nbStrand = rec[fieldIndices[7]]
            sense = rec[fieldIndices[8]]
            strStartEndData.append({'start':startRes[0], 'end':endRes[0],
                                    'nbStrand':nbStrand,
                                    'sense':sense})
        return strStartEndData

    def processTurnData(self, recordDataForMol, fieldIndices, chain):
        # chain ID
        recordDataForChain = filter(lambda x: x[fieldIndices[0]] == chain.id,
                                    recordDataForMol)
        turnStartEndData = []
        for rec in recordDataForChain:
            startData = rec[fieldIndices[1]]+rec[fieldIndices[2]].strip()+\
                        rec[fieldIndices[3]]
            startData = startData.strip()
            startRes = chain.residues.get(startData)
            if len(startRes)!= 1:
                raise ValueError("ERROR: When parsing the helix information %s\
was not found in chain %s"%(startData, chain.id))

            endData = rec[fieldIndices[4]]+rec[fieldIndices[5]].strip()+\
                      rec[fieldIndices[6]]
            endData = endData.strip()
            endRes = chain.residues.get(endData)
            if len(endRes)!= 1:
                raise ValueError("ERROR: When parsing the helix information %s\
was not found in chain %s"%(endData, chain.id))

            comment = rec[fieldIndices[7]]
            turnStartEndData.append({'start':startRes[0], 'end':endRes[0],
                                    'comment':comment})
        return turnStartEndData

    def getMoleculeInformation(self):
        """ Function to retrieve the general informations on the molecule.
        They can be contained in the HEADER records or in the COMPND records.
        This information is used by the molecule chooser to provide
        informations on the molecule selected.
        """
        molStr = self.getRecords(self.allLines, 'HEADER')
        if molStr == []:
            molStr = self.getRecords(self.allLines, 'COMPND')
        if molStr != []: molStr = molStr[0][10:]
        else:  molStr = ''
        return molStr


    def buildMolFromAtomLines(self, atLines, nameStr='NoName',
                              filename=None, objClass=Protein):
        """
        Function to build a molecule from a list of ATOM lines
        """
        mol = objClass()
        self.mol = mol
        mol.parser = self
        self.allLines = atLines
        self.getKeys()
        self.parse_PDB_atoms(atLines, mol)
        self.filename = filename
        mol.name = nameStr
        return mol


    def write_with_new_coords(self, coords, filename=None):
        if len(coords)!=len(self.mol.allAtoms):
            print "length of new coordinates %d does not match %d " %(len(coords), len(self.mol.allAtoms))
        has_filename = filename is not None
        newLines = []
        if has_filename:
            fptr = open(filename, 'w')
        ct = 0
        for l in self.allLines:
            if l.find("ATOM")==0 or l.find("HETATM")==0:
                cc = coords[ct]
                ct = ct + 1
                new_l = l[:30]+"%8.3f%8.3f%8.3f" %(cc[0],cc[1],cc[2]) + l[54:]
            else:
                new_l = l
            newLines.append(new_l)
        if has_filename:
            for new_l in newLines:
                fptr.write(new_l)
        else:
            return newLines
        


class PQRParser(PdbParser):
    """Parser object for the MEAD file format"""
    
    def __init__(self, filename=None, allLines=None, modelsAs='molecules'):
        """Constructor"""
        PdbParser.__init__(self, filename, allLines, modelsAs=modelsAs)

    def findEndChains(self, beg=0, end=None):
        end = len(self.allLines)
        return [(end, 'UNK'),]

    def parse_PDB_ATOM_record(self, rec, mol):
        """Parse PDB ATOM records using white space separated fields.
        We assume 'ATOM' num name resType resNum X Y Z charge radius"""
        orig_rec = rec
        rec = split(rec)
        if len(rec) == 10:
            use_split = True
        else:
            use_split = False
        chainID = 'UNK'

        if chainID != mol.curChain.id:
            if not mol.childByName.has_key(chainID):
                # create a new chain
                mol.curChain = Chain(chainID, mol, top=mol)
            else:
                # reuse the chain with the same name
                mol.curChain = mol.childByName[chainID]
            

        if use_split:
            resName = rec[3]
            resSeq = rec[4]
        else:            
            resName = orig_rec[17:20]
            resSeq = strip(orig_rec[22:26])
            
#       if resSeq != mol.curRes.number or resName != mol.curRes.type:
        if resSeq != mol.curRes.number:
            # check if this residue already exists
            na = strip(resName) + strip(resSeq)
            res = mol.curChain.get( na )
            if res:
                mol.curRes = res[0]
            else:
                mol.curRes = Residue(type=resName, number=resSeq,
                                          parent=mol.curChain,
                                          top=mol)
                # Add a hasCA and hasO flags to the residue and set it to 0
                mol.curRes.hasCA = 0
                mol.curRes.hasO = 0

        if use_split:
            name = rec[2]
        else:
            name = strip(orig_rec[12:17])
        if name == 'CA':
            # Set the flag hasCA to 1 if the atom name is CA
            mol.curRes.hasCA = 1

        if name == 'O' or name == 'OXT':
            # Set the hasO flag to 2 if atom name is O or OXT
            mol.curRes.hasO = 2

        if use_split:
            element = rec[2][0]   
        else:
            element = name[0]
        elem = lower(element)
        if elem =='l':
            element = 'Xx'
        atom = Atom(name, mol.curRes, element, top=mol)

        if use_split:
            if rec[0]=='ATOM': atom.hetatm = 0
            else: atom.hetatm = 1
        else:
            if orig_rec[0:4]=='ATOM': atom.hetatm = 0
            else: atom.hetatm = 1
        #atom.element = element   # done in Atom constructor
        if use_split:
            atom.number = int(rec[1])        
        else:
            atom.number = int(strip(orig_rec[7:12]))

        mol.atmNum[atom.number] = atom

        if use_split:
            atom._coords = [ [ float(rec[5]), float(rec[6]), float(rec[7]) ] ]
            atom._charges['pqr'] = float(rec[8])
            atom.radius = float(rec[9])
        else:
            atom._coords = [ [ float(strip(orig_rec[30:38])), float(strip(orig_rec[38:46])),
                               float(strip(orig_rec[46:54])) ] ]
            atom._charges['pqr'] = float(strip(orig_rec[55:62]))
            atom.radius = float(orig_rec[63:])
        
        atom.chargeSet = 'pqr'
        atom.pqrRadius = atom.radius
        atom.altname = None
        #atom.alternate = []
        return atom
        

class PdbqParser(PdbParser):
    """Class pdbqParser.

    Implements PdbqParser. The objects are PDBQ file readers.
        parsers can be registered for every PDB card. By default, ATOM and HETATM
    are parsed. Such objects are typically passed as an argument to the read
    method of a molecule object. They build a 4 level tree mol:chain:residue:atom

    Example:

        from MolKit.protein import Protein
        from MolKit.pdbqParser import PdbqParser
        mol = Protein()
        mol.read( '../btn-1.pdbq', PdbqParser() )

    NB: - in PdbqParser All atoms fields are parsed according to PDB specifications
      PLUS added field for charge
    """

    def set_stringReprs(self, mol):
        #do not reset allAtoms if there is a torTree which has a ROOT
        if not self.pdbRecordParser.has_key('ROOT'):
            mol.allAtoms = mol.chains.residues.atoms
            # after the hierarchy has been built go set all stringRepr
            # their concise values
            mname = mol.name
            strRpr = mname + ':::'
            mol.allAtoms.setStringRepr(strRpr)

            strRpr = mname + ':'
            mol.chains.setStringRepr(strRpr)

            for c in mol.chains:
                cname = c.id
                strRpr = mname + ':' + cname + ':'
                c.residues.setStringRepr(strRpr)
                for r in c.residues:
                    rname = r.name
                    strRpr = mname + ':' + cname + ':' + rname + ':'
                    r.atoms.setStringRepr(strRpr)


    def defaultReadOptions(self, readOptions=None):
        self.pdbRecordParser['ATOM']  = self.parse_PDB_atoms
        self.pdbRecordParser['HYDBND'] = self.parse_PDB_HYDBND
        self.pdbRecordParser['ROOT'] = self.parse_PDB_ROOT
        self.pdbRecordParser['TORSDOF'] = self.parse_PDB_TORSDOF
        self.pdbRecordParser['REMARK'] = self.parse_PDB_REMARK
        self.pdbRecordParser['CONECT']=self.parse_PDB_CONECT
        ## ADDED 5/12/03, DST
        #if readOptions is not None and len(readOptions) != 0:
        #    for key in readOptions:
        #        if key in self.PDBtags:
        #            self.SetReadOptions(key)



    def parse_PDB_ROOT(self, rec):
        from MolKit.torTree import TorTree
        if len(rec): 
            self.mol.torTree = TorTree(self)
            self.mol.ROOT = self.mol.chains.residues.atoms.get(lambda x: x._uniqIndex==0)[0]
        

    def parse_PDB_TORSDOF(self, rec):
        if not len(rec): return
        self.mol.TORSDOF = int(split(rec[0])[1])


    def parse_PDB_REMARK(self, rec):
        if not len(rec): 
            return
        found = 0
        for item in rec:
            if find(item, 'active torsions')>-1:
                l = item
                found = 1
                llist = split(l)
                self.mol.ndihe = int(llist[1])
                break


    ### Overwrite the pdbParser parse_Atom_Record function.
    def parse_PDB_ATOM_record(self, rec, mol):
        """Parse PDB ATOM records using the pdb columns specifications"""
        # Handle the alternate location using a flag.
        rec = strip(rec)
        if rec[16]!= ' ': self.altLoc = rec[16]
        else: self.altLoc = None

        chainID = rec[21]
        # check for chains break
        if chainID != mol.curChain.id:
        # parse atom info
            #if not mol.childByName.has_key(chainID):
            if not mol.chains.get(lambda x: x.id == chainID):
                # create a new chain
                mol.curChain = Chain(chainID, mol, top=mol)
            else:
                # reuse the chain with the same name
                #mol.curChain = mol.childByName[chainID]
                chains = mol.chains.get(lambda x: x.id == chainID)
                mol.curChain = chains[0]

        # check for residue break
        resName = rec[17:20]
        resSeq = strip(rec[22:26]) #WARNING resSeq is a STRING
        if rec[26] != ' ': icode = rec[26]
        else: icode = ''

        curRes = mol.curRes
        if curRes.parent is None or curRes.parent.name!=chainID or \
                resSeq != curRes.number or  resName != curRes.type or \
                icode != curRes.icode:
            #if resSeq != curRes.number or  resName != curRes.type or \
            #        icode != curRes.icode:
            # check if this residue already exists
            na = strip(resName) + resSeq + icode
            try:
                mol.curRes = mol.curChain.childByName[na]
            except KeyError:   #child not found
                mol.curRes = Residue(resName, resSeq, icode, mol.curChain,
                                     top=mol)
            #res = mol.curChain.get( na )
            #if res:
               # mol.curRes = res[0]
            #else:
                #mol.curRes = Residue(resName, resSeq, icode,
                                          #mol.curChain,
                                          #top=mol)
                # Add a hasCA and hasO flags to the residue and set it to 0
                mol.curRes.hasCA = 0
                mol.curRes.hasO = 0
            

        # handle atom names (calcium, hydrogen) and find element type
        # check validity of chemical element column and charge column
        name, element, charge = self.getPDBAtomName(rec[12:16],"  ",None)
        if name == 'CA':
            # Set the flag hasCA to 1 if the atom name is CA
            mol.curRes.hasCA = 1

        if name == 'O' or name == 'OXT':
            # Set the hasO flag to 2 if atom name is O or OXT
            mol.curRes.hasO = 2

        elem = lower(element)
        if elem =='l':
            autodock_element='element'
            element = 'Xx'
        elif element=='A':
            autodock_element='A'
            element='C'
        elif element == 'n':
            autodock_element='n'
            element='N'
        elif element == 'f':
            autodock_element='f'
            element='Fe'
        elif element == 'Fe':
            autodock_element='f'
            element='Fe'
        elif element == 'Cl':
            autodock_element='c'
            element='Cl'
        elif element == 'c':
            autodock_element='c'
            element='Cl'
        elif element == 'Br':
            autodock_element='b'
            element='Br'
        elif element == 'b':
            autodock_element='b'
            element='Br'
        elif element == 'Zn':
            autodock_element='Zn'
            element='Zn'
        elif element[0] == 'Z':
            autodock_element='Z'
            if len(element)>1:
                element=element[1]
            else:
                element='C'
        else:
            autodock_element=element

        atom = Atom(name, mol.curRes, element, top=mol)
        atom._coords = [ [ float(rec[30:38]), float(rec[38:46]),
                     float(rec[46:54]) ] ]
##          atom.segID = strip(rec[72:76])
        if rec[:4]=='ATOM': atom.hetatm = 0
        else: atom.hetatm = 1
        #atom.alternate = []
        #atom.element = element   # done in Atom constructor
        atom.autodock_element = autodock_element
        atom.number = int(rec[6:11])
        mol.atmNum[atom.number] = atom
#        atom.conformation = 0
        if rec[54:60] == "      " or len(rec[54:60]) == 0:

            atom.occupancy = 0.0
        else:
            atom.occupancy = float(rec[54:60])

        if len(rec) > 61:
            if rec[60:66] == "      ":
                atom.temperatureFactor = 0.0
            else:
                atom.temperatureFactor = float(rec[60:66])
        else: atom.temperatureFactor= 0.0
            
        ##THIS IS THE KEY PDBQ DIFFERENCE:
        ##NB:in older versions of pdbq charge is in columns 55-60
        if len(rec) >=76:
            charge = float(strip(rec[70:76]))
        else:
            charge = float(strip(rec[55:61]))
        if charge is not None: 
            atom._charges['pdbq'] = charge
            atom.chargeSet = 'pdbq'

        atom.altname = None
        if self.altLoc :
            # check if the name of the atom is the same than the
            #name of the previous atom .
            name = name + '@'+self.altLoc
            atom.name = name
            if len(mol.curRes.atoms)>1:
                # the new atom has been add to the current residue
                # You have to go to the one before.
                lastAtom = mol.curRes.atoms[-2]
                atom.altname = split(lastAtom.name, '@')[0]
                if split(name, '@')[0] == atom.altname:
                    # Add the new alternate atom to the LastAtom.alternate and
                    # add the lastAtom to the atom.alternate.
                    lastAtom.alternate.append(atom)
                    atom.alternate.append(lastAtom)
                    for l in lastAtom.alternate:
                        if atom.name != l.name:
                            atom.alternate.append(l)
                            l.alternate.append(atom)
        return atom




"""Class pdbqtParser.

Implements PdbqtParser. The objects are PDBQT file readers.
parsers can be registered for every PDB card. By default, ATOM and HETATM
are parsed. Such objects are typically passed as an argument to the read
method of a molecule object. They build a 4 level tree mol:chain:residue:atom

Example:
    
    from MolKit.protein import Protein
    from MolKit.pdbqtParser import PdbqtParser
    mol = Protein()
    mol.read( '../parts.pdbqt', PdbqtParser() )

NB: - in PdbqtParser All atoms fields are parsed according to PDB specifications
      PLUS added fields for charge (q) and autodock type (t). Currently
      autodock types are elementname+flag for hydrogen bond donors ('D') or
      hydrogen bond acceptors ('A'). For example the autodock_type of
      a nitrogen atom able to hydrogen bond is written as 'NA' with 
      the N in column 77(zero-based) and the autodock_type in column 78.
      ALSO 'A' for aromatic carbons....this causes problems
"""

class PdbqtParser(PdbqParser):


    def __init__(self, filename=None, allLines=None, modelsAs='molecules'):
        PdbqParser.__init__(self, filename, allLines, modelsAs=modelsAs)
        self.useAbove73 = 1


    def parse_PDB_REMARK(self, rec):
        if not len(rec): 
            return
        found = 0
        vina_results = []
        for item in rec:
            if find(item, 'active torsions')>-1:
                l = item
                found = 1
                llist = split(l)
                self.mol.ndihe = int(llist[1])
            if find(item, 'VINA RESULT')>-1:
                #print "processing vr=", item
                l = item
                llist = split(l)
                #vina_results.append("VINA RESULT: ")
                vina_results.append(llist[3:])
        self.mol.vina_results = vina_results
        #if not found: return
        #llist = split(l)
        #self.mol.ndihe = int(llist[1])


    def defaultReadOptions(self, readOptions=None):
        self.pdbRecordParser['ATOM']  =self.parse_PDB_atoms
        self.pdbRecordParser['HYDBND']=self.parse_PDB_HYDBND
        self.pdbRecordParser['ROOT']  =self.parse_PDB_ROOT
        self.pdbRecordParser['TORSDOF']=self.parse_PDB_TORSDOF
        self.pdbRecordParser['REMARK']=self.parse_PDB_REMARK
        self.pdbRecordParser['CONECT']=self.parse_PDB_CONECT
        self.pdbRecordParser['BEGIN_']=self.parse_PDB_BEGIN_RES
        ## ADDED 5/12/03, DST
        #if readOptions is not None and len(readOptions) != 0:
        #    for key in readOptions:
        #        if key in self.PDBtags:
        #            self.SetReadOptions(key)


    def parse_PDB_BEGIN_RES(self, rec):
        #print 'in parse_PDB_BEGIN_RES with rec=', rec
        if len(rec): 
            if not hasattr(self.mol, 'flex_res_list'): 
                self.mol.flex_res_list = []
            for s in rec:
                self.mol.flex_res_list.append(s[10:-1])
        

    def parse_PDB_ROOT(self, rec):
        from MolKit.torTree import TorTree
        if not len(rec): return
        #print "in parse_PDB_ROOT: len(rec)=", len(rec)
        if hasattr(self.mol, 'torTree'):
            print "skipping parse_PDB_ROOT"
            return
        if len(rec)==1:
            self.mol.torTree = TorTree(self)
            self.mol.ROOT = self.mol.chains.residues.atoms.get(lambda x: x._uniqIndex==0)[0]
            return
        if len(rec)>1:
            #build a complicated molecule with many torTrees
            num_res = len(rec)-1
            allLines = self.allLines
            torsdof_line = None
            for l in allLines:
                if l.find("TORSDOF")==0:
                    torsdof_line = l
                    #print "found TORSDOF so end of true_ligand_atoms"
                    break
            if torsdof_line:
                end_lig = allLines.index(torsdof_line)+1
                self.allLines = allLines[:end_lig]
            #print "now len(self.allLines) is ", len(self.allLines)
            # build torTree for true_ligand_atoms
            self.mol.torTree = TorTree(self)
            self.mol.ROOT = self.mol.chains.residues.atoms.get(lambda x: x._uniqIndex==0)[0]
            #print "built tortree for true_ligand_atoms"
            #now try to build a 'torTree' for each flexible residues
            ctr = 0
            all_res_lines = []
            res_lines = []
            if torsdof_line:
                for l in allLines[end_lig:]:
                    if l.find("BEGIN_RES")==0:
                        #start a new section
                        res_lines = [l]
                        ctr += 1
                        #print "found another residue! now there are ", ctr
                    elif l.find("END_RES")==0:
                        # store this section, possibly
                        res_lines.append(l)
                        all_res_lines.append(res_lines)
                        #print "found end of a residue which had ", len(res_lines), " lines"
                    else:
                        res_lines.append(l)
            flex_res = self.flex_res = []
            for i in range(ctr):
                res_lines = all_res_lines[i]
                self.allLines = res_lines
                keys = map(lambda x: strip(x[:6]), res_lines)
                ntors = keys.count("BRANCH")
                #index of this chain is -ctr + i
                ind = ctr - i
                #build a torTree and assign it to this residue, 
                # which is the last residue in its chain (?)
                #find the chain from the name
                #1/11/2007 HACK!!
                #res = self.mol.chains[-ind].residues[-1]
                res = self.mol.chains.residues[-ind]
                res.torTree = TorTree(self)
                res.ROOT = res.atoms.get(lambda x: x._uniqIndex==0)[0]
                res.ntors = ntors
                #print 'set ntors to ', ntors
                flex_res.append(res)
            self.allLines = allLines
            self.mol.hasFlexRes = ctr


    def parse_PDB_TORSDOF(self, rec):
        if not len(rec): return
        #print "TORSDOF:rec=", rec
        for rr in rec:
            if not hasattr(self.mol, 'TORSDOF'):
                self.mol.TORSDOF = int(split(rr)[1])
            else:
                self.mol.TORSDOF += int(split(rr)[1])


    def parse_PDB_ATOM_record(self, rec, mol):
        """Parse PDB ATOM records using the pdb columns specifications"""
        # Handle the alternate location using a flag.
        self.useAbove73 = 1
        rec = strip(rec)
        if rec[16]!= ' ': self.altLoc = rec[16]
        else: self.altLoc = None

        # check for chains break
        chainID = rec[21]

        if chainID != mol.curChain.id:
            #if not mol.childByName.has_key(chainID):
            if not mol.chains.get(lambda x: x.id == chainID):
                # create a new chain
                mol.curChain = Chain(chainID, mol, top=mol)
            else:
                # reuse the chain with the same name
                #mol.curChain = mol.childByName[chainID]
                chains = mol.chains.get(lambda x: x.id == chainID)
                mol.curChain = chains[0]

        # check for residue break
        resName = rec[17:20]
        resSeq = strip(rec[22:26]) #WARNING resSeq is a STRING
        if rec[26] != ' ': icode = rec[26]
        else: icode = ''

        curRes = mol.curRes
        #if resSeq != curRes.number or  resName != curRes.type or \
        if curRes.parent is None or curRes.parent.name!=chainID or \
                resSeq != curRes.number or  resName != curRes.type or \
                icode != curRes.icode:
            # check if this residue already exists
            na = strip(resName) + resSeq + icode
            try:
                mol.curRes = mol.curChain.childByName[na]
            except KeyError:   #child not found
                mol.curRes = Residue(resName, resSeq, icode,
                                          mol.curChain,
                                          top=mol)
            #res = mol.curChain.get( na )
            #if res:
               # mol.curRes = res[0]
            #else:
                #mol.curRes = Residue(resName, resSeq, icode,
                                          #mol.curChain,
                                          #top=mol)
                # Add a hasCA and hasO flags to the residue and set it to 0
                mol.curRes.hasCA = 0
                mol.curRes.hasO = 0
            
            
        if rec[26] != ' ': mol.curRes.icode = rec[26]

        # parse atom info

        # handle atom names (calcium, hydrogen) and find element type
        # check validity of chemical element column and charge column
        #name, element, charge = self.getPDBAtomName(rec[12:16],"  ",None)
        name, z, y = self.getPDBAtomName(rec[12:16],rec[76:],None)
        if name == 'CA':
            # Set the flag hasCA to 1 if the atom name is CA
            mol.curRes.hasCA = 1

        if name == 'O' or name == 'OXT':
            # Set the hasO flag to 2 if atom name is O or OXT
            mol.curRes.hasO = 2

        autodock_element = strip(rec[76:])

        if autodock_element=='A':
            element = 'C'
        elif autodock_element=='Z':
            #AD4 covalent map type: 
            #  try set element according to name
            if len(name)>1:
                element = name[1]
            else:
                element = 'C'
        elif len(autodock_element)==1:
            element = autodock_element
        elif autodock_element in ['NA','SA','OA','HD']:
            element = autodock_element[0]
        else:
            element = autodock_element
            

        atom = Atom(name, mol.curRes, element, top=mol)
        atom._coords = [ [ float(rec[30:38]), float(rec[38:46]),
                     float(rec[46:54]) ] ]

        atom.autodock_element = autodock_element

        if rec[:4]=='ATOM': atom.hetatm = 0
        else: atom.hetatm = 1
        #atom.element = element   # done in Atom constructor
        #NEW FORMAT: autodock_element written in element columns 76+77
        ##10/7: in columns 77 and 78
        #atom.autodock_element = strip(rec[76:78])
        #5/19:CAUTION:this may break!
        #atom.element = strip(rec[76:78])
        #atom.autodock_element = strip(rec[77:])
        #repair oddities caused by autodock element in 76+77
        #if len(atom.element)>1 and atom.element[0] in ['H','N','O','S']:
        #    atom.element = atom.element[0]
        #need to make aromatic carbons, carbons
        #if len(atom.element)==1 and atom.element[0]=='A':
        #    atom.element = 'C'
        #    atom.autodock_element= 'A'
        atom.number = int(rec[6:11])
        mol.atmNum[atom.number] = atom
        if rec[54:60]=="      " or len(rec[54:60])==0:
            atom.occupancy = 0.0
        else:
            atom.occupancy = float(rec[54:60])
        if len(rec) > 61:
            if rec[60:66]=="      " or len(rec[60:66])==0:
                atom.temperatureFactor = 0.0
            else:
                atom.temperatureFactor=float(rec[60:66])
        else: atom.temperatureFactor= 0.0
        ##THIS IS THE KEY PDBQ DIFFERENCE:
##          atom.segID = strip(rec[72:76])
        charge = float(strip(rec[70:76]))
        if charge is not None: 
            atom._charges['pdbqt'] = charge
            atom.chargeSet = 'pdbqt'

        #5/19:CAUTION this may break!
        ##atom.autodock_element = atom.element
        #For FUTURE use for aromatic vs aliphatic carbon distinction
        #if atom.element=='A':
            #atom.element = 'C'
        ##don_acc = strip(rec[78:80])
        ##if len(don_acc):
            ##atom.autodock_element = atom.element + don_acc
            #eg HD or OA or NA
            #atom.autodock_element = atom.element + 'H'

#        ##THIS IS THE KEY PDBQS DIFFERENCE:
#        #correct 5/15: these were switched
#        #atom.AtSolPar=float(strip(rec[78:84]))
#        #atom.AtVol=float(strip(rec[86:92]))
#        atom.AtVol = float(strip(rec[78:84]))
#        atom.AtSolPar = float(strip(rec[86:92]))

#        atom.altname = None
#        if self.altLoc :
#            # check if the name of the atom is the same than the
#            #name of the previous atom .
#            name = name + '@'+self.altLoc
#            atom.name = name
#            if len(mol.curRes.atoms)>1:
#                # the new atom has been add to the current residue
#                # You have to go to the one before.
#                lastAtom = mol.curRes.atoms[-2]
#                atom.altname = split(lastAtom.name, '@')[0]
#                if split(name, '@')[0] == atom.altname:
#                    # Add the new alternate atom to the LastAtom.alternate and
#                    # add the lastAtom to the atom.alternate.
#                    lastAtom.alternate.append(atom)
#                    atom.alternate.append(lastAtom)
#                    for l in lastAtom.alternate:
#                        if atom.name != l.name:
#                            atom.alternate.append(l)
#                            l.alternate.append(atom)
        return atom



"""Class pdbqsParser.

Implements PdbqsParser. The objects are PDBQ file readers.
parsers can be registered for every PDB card. By default, ATOM and HETATM
are parsed. Such objects are typically passed as an argument to the read
method of a molecule object. They build a 4 level tree mol:chain:residue:atom

Example:
    
    from MolKit.protein import Protein
    from MolKit.pdbqsParser import PdbqsParser
    mol = Protein()
    mol.read( '../parts.pdbqs', PdbqsParser() )

NB: - in PdbqsParser All atoms fields are parsed according to PDB specifications
      PLUS added fields for charge, AtSolPar-atomic solvent parameter and AtVol
      atomic volume (which are used by AutoDock in force field calculations)
"""

class PdbqsParser(PdbParser):


    def parse_PDB_ATOM_record(self, rec, mol):
        """Parse PDB ATOM records using the pdb columns specifications"""
        # Handle the alternate location using a flag.
        rec = strip(rec)
        if rec[16]!= ' ': self.altLoc = rec[16]
        else: self.altLoc = None

        # check for chains break
        chainID = rec[21]

        if chainID != mol.curChain.id:
            #if not mol.childByName.has_key(chainID

            if not mol.chains.get(lambda x: x.id == chainID):
                # create a new chain
                mol.curChain = Chain(chainID, mol, top=mol)
            else:
                # reuse the chain with the same name
                #mol.curChain = mol.childByName[chainID]
                chains = mol.chains.get(lambda x: x.id == chainID)
                mol.curChain = chains[0]

        # check for residue break
        resName = rec[17:20]
        resSeq = strip(rec[22:26]) #WARNING resSeq is a STRING
        if rec[26] != ' ': icode = rec[26]
        else: icode = ''

        curRes = mol.curRes
        if curRes.parent is None or curRes.parent.name!=chainID or \
                resSeq != curRes.number or  resName != curRes.type or \
                icode != curRes.icode:
            #   if resSeq != curRes.number or  resName != curRes.type or \
            #        icode != curRes.icode:
            # check if this residue already exists
            na = strip(resName) + resSeq + icode
            try:
                mol.curRes = mol.curChain.childByName[na]
            except KeyError:   #child not found
                mol.curRes = Residue(resName, resSeq, icode,
                                          mol.curChain,
                                          top=mol)
            #res = mol.curChain.get( na )
            #if res:
               # mol.curRes = res[0]
            #else:
                #mol.curRes = Residue(resName, resSeq, icode,
                                          #mol.curChain,
                                          #top=mol)
                # Add a hasCA and hasO flags to the residue and set it to 0
                mol.curRes.hasCA = 0
                mol.curRes.hasO = 0
            
            
        if rec[26] != ' ': mol.curRes.icode = rec[26]

        # parse atom info

        # handle atom names (calcium, hydrogen) and find element type
        # check validity of chemical element column and charge column
        name, element, charge = self.getPDBAtomName(rec[12:16],"  ",None)
        if name == 'CA':
            # Set the flag hasCA to 1 if the atom name is CA
            mol.curRes.hasCA = 1

        if name == 'O' or name == 'OXT':
            # Set the hasO flag to 2 if atom name is O or OXT
            mol.curRes.hasO = 2

        elem = lower(element)
        if elem =='l':
            autodock_element='element'
            element = 'Xx'
        elif element=='A': 
            autodock_element='A'
            element='C'
        elif element == 'n':
            autodock_element='n'
            element='N'
        elif element == 'f':
            autodock_element='f'
            element='Fe'
        elif element == 'Fe':
            autodock_element='f'
            element='Fe'
        elif element == 'Cl':
            autodock_element='c'
            element='Cl'
        elif element == 'c':
            autodock_element='c'
            element='Cl'
        elif element == 'Br':
            autodock_element='b'
            element='Br'
        elif element == 'b':
            autodock_element='b'
            element='Br'
        elif element == 'Zn':
            autodock_element='Zn'
            element='Zn'
        elif element[0] == 'Z':
            autodock_element='Z'
            if len(element)>1:
                element=element[1]
            else:
                element='C'
        else:
            autodock_element=element

        atom = Atom(name, mol.curRes, element, top=mol)
        atom._coords = [ [ float(rec[30:38]), float(rec[38:46]),
                     float(rec[46:54]) ] ]
        if rec[:4]=='ATOM': atom.hetatm = 0
        else: atom.hetatm = 1
        #atom.alternate = []
        #atom.element = element   # done in Atom constructor
        atom.autodock_element = autodock_element
        atom.number = int(rec[6:11])
        mol.atmNum[atom.number] = atom
        if rec[54:60]=="      " or len(rec[54:60])==0:
            atom.occupancy = 0.0
        else:
            atom.occupancy = float(rec[54:60])
        if len(rec) > 61:
            if rec[60:66]=="      " or len(rec[60:66])==0:
                atom.temperatureFactor = 0.0
            else:
                atom.temperatureFactor=float(rec[60:66])
        else: atom.temperatureFactor= 0.0

        ##THIS IS THE KEY PDBQ DIFFERENCE:
##          atom.segID = strip(rec[72:76])
        charge = float(strip(rec[70:76]))
        if charge is not None: 
            atom._charges['pdbqs'] = charge
            atom.chargeSet = 'pdbqs'

        ##THIS IS THE KEY PDBQS DIFFERENCE:
        #correct 5/15: these were switched
        #atom.AtSolPar=float(strip(rec[78:84]))
        #atom.AtVol=float(strip(rec[86:92]))
        atom.AtVol = float(strip(rec[78:84]))
        atom.AtSolPar = float(strip(rec[86:92]))

        atom.altname = None
        if self.altLoc :
            # check if the name of the atom is the same than the
            #name of the previous atom .
            name = name + '@'+self.altLoc
            atom.name = name
            if len(mol.curRes.atoms)>1:
                # the new atom has been add to the current residue
                # You have to go to the one before.
                lastAtom = mol.curRes.atoms[-2]
                atom.altname = split(lastAtom.name, '@')[0]
                if split(name, '@')[0] == atom.altname:
                    # Add the new alternate atom to the LastAtom.alternate and
                    # add the lastAtom to the atom.alternate.
                    lastAtom.alternate.append(atom)
                    atom.alternate.append(lastAtom)
                    for l in lastAtom.alternate:
                        if atom.name != l.name:
                            atom.alternate.append(l)
                            l.alternate.append(atom)
        return atom


if __name__ == '__main__':

    import sys, pdb
    from MolKit.protein import Protein
    print "reading molecule"
    mol = Protein()
    mol.read("/tsri/pdb/struct/%s.pdb"%sys.argv[1], PdbParser())
    print "Done"

# bond stuff
    print "building bonds"
    mol.parser.parse_PDB_CONECT(mol.parser.getRecords(mol.parser.parser.allLines, 'CONECT'))
    mol.buildBondsByDistance()
    bonds = mol.chains.residues.atoms.bonds
    print "Done"

    #mol1 = Protein()
    #mol1.read("xaa.aypyd.pqr", PQRParser())

    # parse multi model file into different conformations
    #mol2 = Protein()
    #mol2.read("/tsri/pdb/struct/1b1g.pdb", PdbParser(0))
    
    #from pdbParser import PdbParser
    #mol.read("/tsri/pdb/struct/2gls.pdb", PdbParser())
    #sys.exit()
    
##      print mol
##      print mol.chains
##      print mol.chains[0]
##      print mol.chains[0].residues
##      print mol.chains[0].residues[6]
##      print mol.chains[0].residues[6].atoms
##      print mol.chains[0].residues[8:16]
    
# atoms selection
    s1 = mol.chains.residues.atoms.get(lambda x:x.name[0]=='C')
    s2 = mol.chains.residues.atoms.get(lambda x:x.name=='CA' and x.parent.type=='CYS')

#residues selection
    polarNegativeCharged = [ 'ASP', 'GLN' ]
    r1 = mol.chains.residues.get(lambda x, rl=polarNegativeCharged: x.type in rl)

# get TreeNodeSet attributs
    s1.name
    s1.parent.name
    s1.parent.uniq().name

# set TreeNodeSet attributs
    mol.chains.residues.atoms.tag = 0
    s1.tag = 1
    mol.chains.residues.atoms.get(lambda x: x.tag==1)
    mol.chains.residues.atoms.color = ((1.,0., 0.),)

# inspect content of TreeNode object
    mol.dump()
    mol.chains[0].dump()
    mol.chains[0].residues[3].dump()
    mol.chains[0].residues[3].atoms[2].dump()

#subtree
    ca = mol.chains[0].residues[0].get(lambda x:x.name=='CA')
    cb = mol.chains[0].residues[0].get(lambda x:x.name=='CB')
    at = mol.subTree(ca[0], cb[0], mol.chains[0])

    a = mol.chains.residues.atoms[0]
    a.getRoot()
    a.parent.getRoot()
    a.parent.parent.getRoot()
    mol.getRoot()
    from MolKit.molecule import AtomSet

# find atoms with laternate lcoations
    mol.chains.residues.atoms.get(lambda x: hasattr(x, 'altLoc'))
    
#compute molecular surfaces for all residues
#    import msms
#    srf = []
#    mol.defaultRadii()
#    for res in mol.chains.residues:
#        s = msms.MSMS( coords = res.atoms.coord, radii=res.atoms.radius )
#        s.compute()
#        s.display()


# test secondary structure stuff
    mol.getSSFromFile()

# test tree spliting and merging

# first test split the fisrt chain in 2
    nbc = len(mol.chains)
    tot = len(mol.chains.residues)
    c2 = mol.chains[0].split( mol.chains[0].residues[20:] )
    assert len(mol.chains) == nbc+1
    assert len(mol.chains.residues)==tot
    
# split all residues in a chain
    resc = mol.chains[0].split( )
    
