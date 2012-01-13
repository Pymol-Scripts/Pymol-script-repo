#############################################################################
#
# Author: Kevin Chan, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/MolKit/pdbWriter.py,v 1.25.2.2 2011/06/09 19:45:43 sanner Exp $
#
# $Id: pdbWriter.py,v 1.25.2.2 2011/06/09 19:45:43 sanner Exp $
#

from MolKit.moleculeWriter import MoleculeWriter
from MolKit.pdbParser import PdbParser
from MolKit.protein import Protein, Chain, ChainSet
from MolKit.protein import Helix, Turn, Strand, SecondaryStructureSet
from MolKit.molecule import Molecule, AtomSet, Atom
from MolKit.tree import TreeNode, TreeNodeSet
import string, os, types
from MolKit.PDBdict import PDBformat, PDBFormatConstr

###############################################################################

class PdbWriter(MoleculeWriter):
    """Class to write data records from a molecule tree to a pdb file.
    Has methods for the user to add own records and to write the record."""
    
##     numRemark = 0
##     numHet = 0
##     numSite = 0
##     numXform = 0
##     numCoord = 0    numTer = 0
##     numConect = 0
##     numSeq = 0

    def __init__(self):
        """Constructor:
        userRecords contains input from user to be written to file.
        missingRecords contains types of mandatory records that are missing
        from userReconds and PdbParser records"""
        self.recordsToWrite = {}
        # Need to add a USER record saying that this PDB file has been generated using PMV vers...
        # by user ... at date.
        # All the PDB records in order
        self.PDBRECORDS = ['HEADER', 'OBSLTE','TITLE', 'CAVEAT',
                           'COMPND', 'SOURCE','KEYWDS', 'EXPDTA',
                           'AUTHOR', 'REVDAT', 'SPRSDE', 'JRNL', 'USER',
                           'REMARK', 'DBREF', 'SEQADV', 'SEQRES',
                           'MODRES', 'HET','HETNAM', 'HETSYN', 'FORMUL',
                           'HELIX', 'SHEET', 'TURN','SSBOND', 'LINK',
                           'HYDBND', 'SLTBRG', 'CISPEP', 'SITE',
                           'CRYST1', 'ORIGX1', 'ORIGX2', 'ORIGX3',
                           'SCALE1', 'SCALE2','SCALE3', 'MTRIX1', 'MTRIX2',
                           'MTRIX3', 'TVECT', 'MODEL','ATOM','SIGATM',
                           'ANISOU','SIGUIJ','TER','HETATM','ENDMDL',
                           'CONECT','MASTER','END'
                           ]
        # PDB records that will be created from the data structure in the right order
        self.FROMDATASTRUCT = ['HELIX','SHEET', 'TURN', 'HYDBND', 'ATOM','TER',
                               'HETATM','CONECT']

    def write(self, filename, nodes, sort=False, sortFunc=None,
              records=['ATOM', 'CONECT'], bondOrigin=('File','UserDefined'),
              ssOrigin='File'):
        """
        required argument:
        filename  -- path to the new file, a .pdb extension will be added
                     when missing.
        nodes     -- TreeNode, TreeNodeSet instance to save as PDB

        optional arguments:
        sort  -- (False) Boolean flag to specify whether or not to sort the
                 given nodes
        sortFunc  -- (None) sort function that will be used to sort the nodes
                     when specified.This function has to return (-1, 0 or 1).
        records -- list of PDB record to write out
        bondOrigin -- (('File', 'UserDefined')) This will be used if the CONECT records are
                      written out. Can be any combination of 'File',
                      'BuiltByDistance' and 'UserDefined'.
        ssOrigin -- 'File' This will be used if the TURN, HELIX and SHEET
                     records are written. Can be either from the originating
                     PDB file or from the data structure.
        """
        self.records = records
        # If the filename doesn't have a pdb extension add it
        fileExt = os.path.splitext(filename)[1]
        if fileExt=='':
            filename = '%s.pdb' %filename
        # Nodes need to be either a TreeNode or a TreeNodeSet instance.
        assert isinstance(nodes, TreeNode) or isinstance(nodes, TreeNodeSet)
        
        # sort the nodes
        if sort and hasattr(nodes, 'sort'):
            nodes.sort(sortFunc)

        if isinstance(nodes, TreeNode): mol = nodes.top   
        elif isinstance(nodes, TreeNodeSet): mol = nodes.top.uniq()[0]
        # get a handle on the molecule parser
        parser = mol.parser

        # get the atoms in nodes
        atmInNode = nodes.findType(Atom)
        
        # Create all the records to be written out.
        # Get the from file records this is possible only if the file
        # comes from a PDB parser.
        if isinstance(parser, PdbParser):
            fileRec = filter(lambda x: not x in self.FROMDATASTRUCT, records)
            for rec in fileRec:
                self.recordsToWrite[rec] = parser.getRecords(parser.allLines,rec)

        # Create the records from the data structure:
        # secondary structure 'HELIX', 'SHEET', 'TURN', 'REMARK 650 HELIX',
        # 'REMARK 700 SHEET', 'REMARK 750 TURN'
        ssRec = filter(lambda x: x in ['HELIX', 'SHEET', 'TURN'], records)
        if len(ssRec):
            self.defineSecondaryStructureSection(mol, origin=ssOrigin)

        # Atom records (ATOM, TER and HETATM)
        atmRec = 'ATOM' in records
        hetRec = 'HETATM' in records
        if atmRec is True or hetRec is True : 
            #self.defineCoordsSection(nodes, sort=sort, sortFunc=sortFunc, atmRec=atmRec, hetRec=hetRec)
            self.defineCoordsSection(atmInNode, sort=sort, sortFunc=sortFunc, atmRec=atmRec, hetRec=hetRec)
            
        # Hydrogen bonds (HYDBND)
        if 'HYDBND' in records:
            self.defineHYDBNDRecords(atmInNode)

        # CONECT records
        #build excluded Bonds from recType here
        #possible recTypes are:
        #all, fileUser, fileDist, file, user, userDist, dist, none
        if 'CONECT' in records:
            self.defineConnectSection(atmInNode, bondOrigin)

        file = open(filename, 'w')
        if not ('  ' or ' ' or '') in [x.chemElem for x in mol.allAtoms.data]:
            file.write('REMARK   4 XXXX COMPLIES WITH FORMAT V. 2.0\n')       

        for rec in self.PDBRECORDS:
            if self.recordsToWrite.has_key(rec):
                recLine = self.recordsToWrite[rec]
                if type(recLine) is types.ListType:
                    for line in recLine:
                        file.write(line)
                else:
                    file.write(line)
        file.close()
        
        #self.recordsToWrite = {}

    def write_atom(self, f, atm):
        """
        Takes a file object and an Atom instance.
        Writes the atom record to the given file."""
        self.recordsToWrite['ATOM'] = []
        # 1- Need to define the ATOM record of the given atom
        self.recordsToWrite['ATOM'].append(self.defineATOM_HETATMRecord(atm))
        # 2- Need to write the ATOM record to the given file.
        f.write(self.recordsToWrite['ATOM'][0])
        
    def formatName(self, at):
        altLoc = None
        name = at.name
        #if len(name)>2 and name[-2]=='@':
        if string.find(name,'@')>-1:
            ind = string.index(name,'@')
            altLoc = name[ind+1:]
            #altLoc = name[-1]
            name = name[:ind]
        if len(name)==4:
            if at.element=='H':
                name = name[-1] + name[:-1]
                nameStr = '%-4.4s' % name
            else:
                nameStr = '%4.4s' % name
        elif len(at.element)==2:
            nameStr = '%-4.4s' % at.element
        else:
            nameStr = ' %-3s' % name
        return nameStr, altLoc
                
    def addRecord(self, key, record=[]):
        """Allows user to enter own record for the record type given by
        key. record should be a list of tuples with a tuple for each line
        in the record for that type. The method checks that
        the user's record fits the PDB format, else there is an assertion
        error.  If if a type is entered but no record, when write() is
        called, it writes the record given by the molecule's parser
        records."""

        self.userRecords[key] = []
        constraints=0
        for line in record:
            assert len(line)==len(PDBFormatConstr[key])
            for num in PDBFormatConstr[key]:
                if num != None:
                    constraints = 1
            if constraints==1:
                i = 0
                for value in line:
                    if PDBFormatConstr[key][i]==None:
                        pass
                    else:
                        assert len(str(value)) <= \
                        PDBFormatConstr[key][i]
                    i = i + 1
            self.userRecords[key].append(PDBformat[key] % line)

    def defineHYDBNDRecords(self, atoms):
        self.userRecords['HYDBND'] = []
        for a in atoms:
            if not hasattr(a, 'hbonds'): continue 
            for b in a.hbonds:
                #only write a record if a is the donor 
                if b.donAt!=a: continue
                #only write a record if all atoms are in same molecule
                if b.donAt.top!=b.accAt.top: continue
                if not b.hAt is None and b.donAt.top!=b.hAt.top: continue
                #columns 1-6 + spaces for columns 7-12
                s = 'HYDBND      '
                #columns 13-16 donor name
                #strip off altloc + save it
                nameStr, altLoc = self.formatName(a)
                s = s + nameStr
                #column 17 donor altLoc indicator
                if altLoc: s = s + altLoc
                else: s = s + ' '
                #columns 18-20 donor parent name (residue type)
                #column 21 space + column 22 chain id
                s = s + a.parent.type + ' ' + a.parent.parent.id 
                #columns 23-27 donor parent number 
                #column 27 insertion code
                #column 28 space
                if not a.parent.icode:
                    s = s + '%5d' %int(a.parent.number) + '  '
                else:
                    s = s + '%5d' %int(a.parent.number) + a.parent.icode + ' '

                # write the OPTIONAL hydrogen atom info
                hAt = b.hAt
                if not hAt is None:
                    s = s + "              "
                else:
                #columns 30-33 hydrogen atom name
                    nameStr, altLoc = self.formatName(hAt)
                    s = s + nameStr
                    #column 34 hydrogen atom  altLoc indicator
                    if altLoc: s = s + altLoc
                    else: s = s + ' '
                    #column 35 space + column 36 chain id
                    s = s + ' ' + hAt.parent.parent.id
                    #columns 37-41 hydrogen atom parent number 
                    # nb: 42 would be insertion code then 43 a space
                    if not hAt.parent.icode:
                        s = s + '%5d' %int(hAt.parent.number) + '  '
                    else:
                        s = s + '%5d' %int(hAt.parent.number)+hAt.parent.icode+' '
                

                # write the acceptor atom info
                acc = b.accAt 
                #columns 44-47 acceptor name
                nameStr, altLoc = self.formatName(acc)
                s = s + nameStr
                #column 48 acceptor altLoc indicator
                if altLoc: s = s + altLoc
                else: s = s + ' '
                #columns 49-51 acceptor parent name (residue type)
                #column 52 space + column 53 chain id
                s = s + acc.parent.type + ' ' + acc.parent.parent.id
                #columns 54-58 acceptor parent number 
                # nb: 59 would be insertion code 
                # ???? non implemented->
                # 60-65 symmetry operator for 1st non-hyd. atom
                if not acc.parent.icode:
                    s = s + '%5d' %int(acc.parent.number) +'  \n'
                else:
                    s = s + '%5d' %int(acc.parent.number) + acc.parent.icode+' \n'
                #s = s + '%5d' %int(acc.parent.number) +'  \n'
                self.recordsToWrite['HYDBND'].append(s)
        

    
    def defineHELIXRecords(self, helix):
        """
        Takes a list of Helix objects and define the corresponding HELIX records
        """
        for ss in helix:
            # 1- Get information from the sheet object
            ssNumber = int(ss.name.split("Helix")[1])
            startRes = ss.start
            endRes = ss.end
            ssChain = ss.chain
            pdbHelClass = ss.helClass
            if ss.comment is None:
                comment = ' '
            else:
                cxomment = ss.comment
            # 2- Create the record with this information
            # [0:5] Record name
            rec = "HELIX "
            # Empty space
            rec = rec + "%1.1s"%''
            # [7:9] Helix number
            rec = rec + "%3d"%ssNumber
            # 10 empty space
            rec = rec + "%1.1s"%''
            # 11-13 Helix identifier
            hID = "H%d"%ssNumber
            rec = rec + "%3.3s"%hID
            # 14 empty space
            rec = rec + "%1.1s"%''
            # 15-17 start residue name
            rec = rec + "%3.3s"%startRes.type
            # 18 empty space
            rec = rec + "%1.1s"%''
            # 19 Chain ID
            rec = rec + "%1.1s"%ssChain.id
            # 20 empty space
            rec = rec + "%1.1s"%''
            # 21-24 res number
            rec = rec + "%4.4s"%startRes.number
            # 25 start residue icode
            rec = rec + "%1.1s"%startRes.icode
            # 26 empty space
            rec = rec + "%1.1s"%''
            # 27-29 end residue type
            rec = rec + "%3.3s"%endRes.type
            # 30 empty space
            rec = rec + "%1.1s"%''
            # 31 Chain ID
            rec = rec + "%1.1s"%ssChain.id
            # 32 empty space
            rec = rec + "%1.1s"%''
            # 33-36 res number
            rec = rec + "%4.4s"%endRes.number
            # 37 start residue icode
            rec = rec + "%1.1s"%endRes.icode
            # 38-39 helClass
            rec = rec + "%2d"%pdbHelClass

            # 40-69 comment
            rec = rec + "%29s"%comment
            rec = rec +"\n"
            self.recordsToWrite['HELIX'].append(rec)

    def defineSHEETRecords(self, sheet):
        """
        Takes a set of Strand objects and define the SHEET records
        """
        for ss in sheet:
            # 1- Get information from the sheet object
            ssNumber = int(ss.name.split("Strand")[1])
            startRes = ss.start
            endRes = ss.end
            ssChain = ss.chain
            nbStrand = ss.nbStrand
            if nbStrand is None: nbStrand=1
            sense = ss.sense
            if sense is None: sense = 0
            # 2- Create the record with this information
            # [0:5] Record name
            rec = "SHEET "
            # Empty space
            rec = rec + "%1.1s"%''
            # [7:9] Strand number
            rec = rec + "%3d"%ssNumber
            # 10 empty space
            rec = rec + "%1.1s"%''
            # 11-13 Sheet identifier
            sID = "S%d"%ssNumber
            rec = rec + "%3.3s"%sID
            # 14-15 Number of strand in the sheet
            rec = rec + "%2d"%nbStrand
            # 16 empty space
            rec = rec + "%1.1s"%''
            # 17-19 start residue name
            rec = rec + "%3.3s"%startRes.type
            # 20 empty space
            rec = rec + "%1.1s"%''
            # 21 Chain ID
            rec = rec + "%1.1s"%ssChain.id
            # 22-25 res number
            rec = rec + "%4.4s"%startRes.number
            # 26 start residue icode
##             if startRes.icode == '':
##                 rec = rec + ' '
##             else:
##                 rec = rec + startRes.icode
            rec = rec + "%1.1s"%startRes.icode
            # 27 empty space
            rec = rec + "%1.1s"%""
            # 28-30 end residue type
            rec = rec + "%3.3s"%endRes.type
            # 31 empty space
            rec = rec + "%1.1s"%""
            # 32 Chain ID
            rec = rec + "%1.1s"%ssChain.id
            # 33-36 res number
            rec = rec + "%4.4s"%endRes.number
            # 37 start residue icode
##             if endRes.icode == '':
##                 rec = rec + ' '
##             else:
##                 rec = rec + endRes.icode
            rec = rec + "%1.1s"%endRes.icode
            # 38-39 strand sense 0, -1 or 1
            rec = rec + "%2d"%sense

            # 40-69 not parsed information
            rec = rec + "%29s"%' '
            rec = rec +"\n"
            self.recordsToWrite['SHEET'].append(rec)
        
    def defineTURNRecords(self, turn):
        """
        Takes a set of  Turn objects and define  the TURN
        records
        """
        for ss in turn:
            # 1- Get information from the sheet object
            ssNumber = int(ss.name.split("Turn")[1])
            startRes = ss.start
            endRes = ss.end
            ssChain = ss.chain
            comment = ss.comment
            if comment is None:
                comment = " "

            # 2- Create the record with this information
            # [0:5] Record name
            rec = "TURN  "
            # Empty space
            rec = rec + "%1.1s"%''
            # [7:9] Turn number
            rec = rec + "%3d"%ssNumber
            # 10 empty space
            rec = rec + "%1.1s"%''
            # 11-13 Turn identifier
            tID = "T%d"%ssNumber
            rec = rec + "%3.3s"%tID
            # 14 empty space
            rec = rec + "%1.1s"%''
            # 15:17 start res type
            rec = rec + "%3.3s"%startRes.type
            # 18 empty space
            rec = rec + "%1.1s"%''
            # 19 Chain ID
            rec = rec + "%1.1s"%ssChain.id
            # 20:23 start res number
            rec = rec + "%4.4s"%startRes.number
            # 24 start residue icode
##             if startRes.icode == '':
##                 rec = rec + ' '
##             else:
##                 rec = rec + startRes.icode
            rec = rec + "%1.1s"%startRes.icode
            # 25 empty space
            rec = rec + "%1.1s"%''
            # 26-28 end residue type
            rec = rec + "%3.3s"%endRes.type
            # 29 empty space
            rec = rec + "%1.1s"%''
            # 30 Chain ID
            rec = rec + "%1.1s"%ssChain.id
            # 31:34 end res number
            rec = rec + "%4.4s"%endRes.number
            # 35 end residue icode
##             if endRes.icode == '':
##                 rec = rec + ' '
##             else:
##                 rec = rec + endRes.icode
            rec = rec + "%1.1s"%endRes.icode

            # 36:39 empty spaces
            rec = rec +"%4s"%''

            # 40:69 comment
            rec = rec + "%29s"%comment
            rec = rec +"\n"
            self.recordsToWrite['TURN'].append(rec)


    def defineATOM_HETATMRecord(self, atm):
        """
        Define the ATOM or HETATM rec for the given atm
        """
        # Here the column are 0 based while in the PDB description it is 1 based.
        # [0:5] Record name
        if atm.hetatm==0: rec = "ATOM  "
        else: rec = "HETATM"
        # [6:11] Atom serial number + a space
        rec = rec+ '%5i '%atm.number
        # [12:15] Atom atmName
        # get the atmName
        atmName = atm.name
        # Separate the alternate location indicator from the name
        if '@' in atmName:
            atmName, altLoc = atmName.split('@')
        else:
            altLoc = ' '
        # Creating the atom name field
        # If the atom name is 4 characters long
        #if len(atmName)==4:
        if len(atmName)>=4:
            if len(atmName)>4:
                atmName = atmName[:4]
            # do something if it is an H
            if atm.element == 'H':
                atmName = atmName[-1]+atmName[:-1]
                rec = rec + "%-4.4s"%atmName
            #else:
            #    rec = rec + "%4.4s"%atmName
            elif len(atm.element)==2:
                rec = rec + "%-4.4s"%atmName
            else:
                rec = rec + ' %-3s'%atmName[:-1]
        elif len(atm.element)==2: # the Ca i.e else of if(atm.name)==4
            rec = rec + '%-4.4s'%atm.element
        else:
            rec = rec + ' %-3s'%atmName

        # [16] Alternate location indicator
        rec = rec + altLoc

        # get the res name, res sequence number and the chain id.
        resName = ''
        resSeq = ''
        chainID = ''
        resIcode = ''
        if hasattr(atm, 'parent'):
            if hasattr(atm.parent, 'type'):
                resName = atm.parent.type
            if hasattr(atm.parent, 'number'):
                resSeq = atm.parent.number
            if hasattr(atm.parent, 'icode'):
                resIcode = atm.parent.icode
            if hasattr(atm.parent, 'parent') and hasattr(atm.parent.parent, 'id'):
                chainID = atm.parent.parent.id
                
        # [17:19] Residue name:
        rec = rec + "%3.3s"%resName
        # [20] space
        rec = rec + ' '
        # [21] Chain identifier.
        rec = rec + "%1.1s"%chainID
        # [22:25] residue sequence number
        rec = rec + '%4.4s'%resSeq
        # [26] res icode
        rec = rec + '%1.1s'%resIcode
        # [27:29] 3 empty spaces
        rec = rec + "%3.3s"%''
        
        coords = atm.coords
        # [30:37] x coordinate
        rec = rec + "%8.3f"%coords[0]
        # [38:45} y coordinate
        rec = rec + "%8.3f"%coords[1]
        # [46:53] z coordinate
        rec = rec + "%8.3f"%coords[2]
        # [54:59] occupancy
        if hasattr(atm, 'occupancy'):
            occupancy = atm.occupancy
        else:
            occupancy = 0.0
        rec = rec + '%6.2f'%occupancy
        # [60:65] Temperature factor
        if hasattr(atm, 'temperatureFactor'):
            tf = atm.temperatureFactor
        else:
            tf = 0.0
        rec = rec + '%6.2f'%tf
        # [66:71] 6 empty spaces
        rec = rec + '%6.6s'% ''
        # [72:75] Segment identifier, left-justified.
        if hasattr(atm, 'segID'):
            rec = rec + '%-4.4s'%atm.segID
        else:
            rec = rec + '%-4.4s'%atm.top.name
        # [76:77] Element symbol, right justified.
        rec = rec + '%2.2s'%atm.element
        # [78:79] Charge of the atom
        if hasattr(atm, 'PDBcharge'):
            rec = rec + '%2.2s'%atm.PDBcharge
        else: 
            rec = rec + '%2.2s'%''
        rec = rec + '\n'
        return rec


    def defineTERRecord(self, atom):
        """
        TER record to the file.
        """
        rec = 'TER   '
        number = atom.number + 1
        rec = rec + '%5i      ' % number
        rec = rec + '%3s %1s'%(atom.parent.type, atom.parent.parent.id)
        rec = rec + '%4s' % atom.parent.number
        if hasattr(atom.parent, 'icode'):
            rec = rec + '%1.1s' % atom.parent.icode
        else: 
            rec = rec + ' '
        rec = rec + '\n'
        return rec

    def defineSecondaryStructureSection(self, mol, origin='File'):
        """
        The Secondary structure section contains the following records:
        HELIX, SHEET, TURN
        Information will taken from either the file or the data structure
        required argument:
        origin can either be '', File or Stride
        
        """
        # 1- Create the SS records from the records contained in the File
        self.recordsToWrite['HELIX'] = []   
        self.recordsToWrite['SHEET'] = []   
        self.recordsToWrite['TURN'] = []   
        if origin == 'File':
            # Directly get the records from the file if molecule comes
            # from a PDB file and the file contains the info
            parser = mol.parser
            if isinstance(parser, PdbParser):
                if not parser.hasSsDataInFile(): return
                self.recordsToWrite['HELIX'] = parser.getRecords(parser.allLines,
                                                            'HELIX')
                self.recordsToWrite['SHEET'] = parser.getRecords(parser.allLines,
                                                            'SHEET')
                self.recordsToWrite['TURN'] = parser.getRecords(parser.allLines,
                                                           'TURN')

            # Else get it from the data structure if information has been
            # parsed
            elif mol.hasSS == ['From File']:
                allstrands = SecondaryStructureSet([])
                allhelices = SecondaryStructureSet([])
                allturns = SecondaryStructureSet([])
                for chain in mol.chains:
                    if not hasattr(chain,'secondarystructureset'): continue
                    sSet = chain.secondarystructureset
                    helices = sSet.get(lambda x: x.name.startswith('Helix'))
                    if not helices is None: allhelices = allhelices + helices
                    strands = sSet.get(lambda x: x.name.startswith('Strand'))
                    if not strands is None: allstrands = allstrands + strands
                    turns = sSet.get(lambda x: x.name.startswith("Turn"))
                    if not turns is None: allturns = allturns + turns
                self.defineHELIXRecords(allhelices)
                self.defineSHEETRecords(allstrands)
                self.defineTURNRecords(allturns)
            else: return

        # 2- Create the SS record from the information generated by stride
        else:
            # only if the data structure has been created
            if not mol.hasSS in ['From Stride', "From PROSS"]: return
            else:
                allstrands = SecondaryStructureSet([])
                allhelices = SecondaryStructureSet([])
                allturns = SecondaryStructureSet([])
                for chain in mol.chains:
                    if not hasattr(chain,'secondarystructureset'): continue
                    sSet = chain.secondarystructureset
                    helices = sSet.get(lambda x: x.name.startswith('Helix'))
                    if not helices is None: allhelices = allhelices + helices
                    strands = sSet.get(lambda x: x.name.startswith('Strand'))
                    if not strands is None: allstrands = allstrands + strands
                    turns = sSet.get(lambda x: x.name.startswith("Turn"))
                    if not turns is None: allturns = allturns + turns
                    
                self.defineHELIXRecords(allhelices)
                self.defineSHEETRecords(allstrands)
                self.defineTURNRecords(allturns)        

    def defineCoordsSection(self, nodes, sort=False, sortFunc=None, atmRec=True, hetRec=True):
        """
        The coordinate section should contain the following records:
        MODEL, ATOM, SIGATM, ANISOU, SIGUIJ, TER, HETATM, ENDMDL
        Here we only save the current conformation and only
        the ATOM, TER, and HETATM records.
        """
        # Only save the current conformation.
        # get the nodes by chains.
        # define the ATOM TER if atmRec is True and HETATM if hetRec is Ture
        # for the given nodes.
        
        # Get all the chains in the nodes.
        chains = nodes.findType(Chain, uniq=1)
        #selAtoms = nodes.findType(Atom)
        self.recordsToWrite['ATOM']=[]
        for c in chains:
            allAtoms = filter(lambda x: x.parent.parent == c, nodes)
            if atmRec:
                atoms = filter(lambda x: x.hetatm==0, allAtoms)
                for atm in atoms:
                    self.recordsToWrite['ATOM'].append(self.defineATOM_HETATMRecord(atm))
                if len(atoms):
                    self.recordsToWrite['ATOM'].append(self.defineTERRecord(atoms[-1]))
            if hetRec:
                hetatm = filter(lambda x: x.hetatm==1, allAtoms)
                for atm in hetatm:
                    self.recordsToWrite['ATOM'].append(self.defineATOM_HETATMRecord(atm))


    def defineConnectSection(self, atms, bondOrigin=('File','UserDefined')):
        """
        The Connectivity section contains the following records:
        CONECT
        bondOrigin -- either a string 'all' or a tuple of string describing the
                      origin of the bonds:
                      'File' : CONECT records of the originating file describing the molecule
                      'BuiltByDistance': Bonds created by distance.
                      'UserDefined' : Bonds added by the user.
        """
        # only write the regular bonds in the CONECT records.
        # the Hbond are defined by the HYBND section
        if bondOrigin == 'all':
            bondOrigin = ('File', 'BuiltByDistance', 'UserDefined')
        elif type(bondOrigin) is types.StringType and \
             bondOrigin in ['File', 'BuiltByDistance', 'UserDefined']:
            bondOrigin = (bondOrigin,)
        elif type(bondOrigin) is types.ListType:
            bondOrigin = tuple(bondOrigin)
        self.recordsToWrite['CONECT'] = []

        ##
        ## MS June 2011
        ## replace thsi code with a loop over bonds whihc is faster but creates
        ## a CONECT record for each bond rather than one for all bonds for a
        ## given atom
        ##
##         from time import time
##         t1 = time()
##         for atm in atms:
##             rec = 'CONECT%5i'%atm.number
##             alLeast1Bond = 0
##             for b in atm.bonds:
##                 if not b.origin in bondOrigin: continue
##                 a2 = b.atom1
##                 if a2==atm: a2 = b.atom2
##                 if not a2 in atms: continue
##                 alLeast1Bond = 1
##                 rec = rec + '%5i'%a2.number
##             if alLeast1Bond:
##                 self.recordsToWrite['CONECT'].append(rec+'\n')
##         print 'Loop1', time()-t1

##         t1 = time()
        allBonds = atms.bonds[0]
        bl = self.recordsToWrite['CONECT']
        for b in allBonds:
           if not b.origin in bondOrigin: continue
           rec = 'CONECT%5i%5i'%(b.atom1.number, b.atom2.number)
           bl.append(rec+'\n')
##         print 'Loop2', time()-t1
        
class PdbqWriter(PdbWriter):
    """Class to write data records from a molecule tree to a pdbq file.
    Has methods for the user to add own records and to write the record."""


    def __init__(self):
        """Constructor:
        """
        PdbWriter.__init__(self)


    def write_records(self, file, molecule):
        """Writes the record types up to ATOM to the pdb file.  For each
        record type, write_records first looks at userRecords to get the
        record info; if they record type is not there, the method looks
        in the parser records, but only for the 'mandatory' or 'required'
        record types, or if the user_Records has the keyword but no record.
        If there is no record for a 'mandatory' type, a warning is printed.
        'required' types are those written if the record is specified by
        the user or if it is in the parser records, but no waring is
        printed.  All record types not mandatory or required are
        optional and must be specified by the user to be written."""

        self.missingRecords = {}
        Parser = molecule.parser
        isPdbParser = isinstance(molecule.parser, PdbParser)
        tags = Parser.PDBtags[0:42]    # or only mandatory types
        mandatory = []
        # added SEQRES to mandatory record types
        # 'REMARK' no longer mandatory
        required = []
        # required: write only if it is specified by user or in old file
        for key in tags:
            if self.userRecords.has_key(key) and \
                    len(self.userRecords[key])!=0:
                RecordList = self.userRecords[key]
            elif isPdbParser and len(Parser.getRecords(Parser.allLines, '%s' % key))!=0 \
                    and ((key in mandatory) or (key in required) or \
                    self.userRecords.has_key(key)):
                RecordList = Parser.getRecords(Parser.allLines, '%s' % key)
            else:
                RecordList = []
                #if key in mandatory:
                #    print 'warning: mandatory record %s missing' % key
                #    self.missingRecords[key] = []
        for line in RecordList:
            file.write('%s' % line)
            if key=='REMARK':
                self.numRemark = self.numRemark + 1
            elif key=='HET   ':
                self.numHet = self.numHet + 1
            elif key=='SITE  ':
                self.numSite = self.numSite + 1
            elif key in ['ORIGX1', 'ORIGX2', 'ORIGX3', 'SCALE1',\
                    'SCALE2', 'SCALE3', 'MTRIX1', 'MTRIX2', 'MTRIX3']:
                self.numXform = self.numXform + 1
            elif key=='SEQRES':
                self.numSeq = self.numSeq + 1
            else: pass 


    def defineATOM_HETATMRecord(self, atm):
        """
        Define the ATOM or HETATM rec for the given atm
        """
        # Here the column are 0 based while in the PDB description it is 1 based.
        # [0:5] Record name
        if atm.hetatm==0: rec = "ATOM  "
        else: rec = "HETATM"
        # [6:11] Atom serial number + a space
        rec = rec+ '%5i '%atm.number
        # [12:15] Atom atmName
        # get the atmName
        atmName = atm.name
        # Separate the alternate location indicator from the name
        if '@' in atmName:
            atmName, altLoc = atmName.split('@')
        else:
            altLoc = ' '
        # Creating the atom name field
        # If the atom name is 4 characters long
        #if len(atmName)==4:
        if len(atmName)>=4:
            if len(atmName)>4:
                atmName = atmName[:4]
            # do something if it is an H
            if atm.element == 'H':
                atmName = atmName[-1]+atmName[:-1]
                rec = rec + "%-4.4s"%atmName
            #else:
            #    rec = rec + "%4.4s"%atmName
            elif len(atm.element)==2:
                rec = rec + "%-4.4s"%atmName
            else:
                rec = rec + ' %-3s'%atmName[:-1]
        elif len(atm.element)==2: # the Ca i.e else of if(atm.name)==4
            rec = rec + '%-4.4s'%atm.element
        else:
            rec = rec + ' %-3s'%atmName
        # [16] Alternate location indicator
        rec = rec + altLoc
        # get the res name, res sequence number and the chain id.
        resName = ''
        resSeq = ''
        chainID = ''
        resIcode = ''
        if hasattr(atm, 'parent'):
            if hasattr(atm.parent, 'type'):
                resName = atm.parent.type
            if hasattr(atm.parent, 'number'):
                resSeq = atm.parent.number
            if hasattr(atm.parent, 'icode'):
                resIcode = atm.parent.icode
            if hasattr(atm.parent, 'parent') and hasattr(atm.parent.parent, 'id'):
                chainID = atm.parent.parent.id
        # [17:19] Residue name:
        rec = rec + "%3.3s"%resName
        # [20] space
        rec = rec + ' '
        # [21] Chain identifier.
        rec = rec + "%1.1s"%chainID
        # [22:25] residue sequence number
        rec = rec + '%4.4s'%resSeq
        # [26] res icode
        rec = rec + '%1.1s'%resIcode
        # [27:29] 3 empty spaces
        rec = rec + "%3.3s"%''

        coords = atm.coords
        # [30:37] x coordinate
        rec = rec + "%8.3f"%coords[0]
        # [38:45} y coordinate
        rec = rec + "%8.3f"%coords[1]
        # [46:53] z coordinate
        rec = rec + "%8.3f"%coords[2]
        # [54:59] occupancy
        if hasattr(atm, 'occupancy'):
            occupancy = atm.occupancy
        elif hasattr(atm, '_charges') and atm._charges.has_key('pqr'):
            occupancy = atm._charges['pqr']
        else:
            occupancy = 0.0
        rec = rec + '%6.2f'%occupancy
        # [60:65] Temperature factor
        if hasattr(atm, 'temperatureFactor'):
            tf = atm.temperatureFactor
        elif hasattr(atm, 'pqrRadius'):
            tf = atm.pqrRadius
        else:
            tf = 0.0
        rec = rec + '%6.2f'%tf
        # [66:70] 4 empty spaces
        rec = rec + '%4.4s'%''
        
        # [71:77] charges
        if hasattr(atm, 'charge'):
            rec = rec + '%6.3f'%atm.charge
        elif hasattr(atm, 'gast_charge'):
            rec = rec + '%6.3f'%atm.gast_charge
        else:
            rec = rec + '%6s'%''
        rec = rec + '\n'
        return rec


class PdbqsWriter(PdbqWriter):
    """Class to write data records from a molecule tree to a pdbq file.
    Has methods for the user to add own records and to write the record."""

    def __init__(self):
        """Constructor:
        """
        PdbWriter.__init__(self)


    def defineATOM_HETATMRecord(self, atm):
        """
        Define the ATOM or HETATM rec for the given atm
        """
        # Here the column are 0 based while in the PDB description it is 1 based.
        # [0:5] Record name
        if atm.hetatm==0: rec = "ATOM  "
        else: rec = "HETATM"
        # [6:11] Atom serial number + a space
        rec = rec+ '%5i '%atm.number
        # [12:15] Atom atmName
        # get the atmName
        atmName = atm.name
        # Separate the alternate location indicator from the name
        if '@' in atmName:
            atmName, altLoc = atmName.split('@')
        else:
            altLoc = ' '
        # Creating the atom name field
        # If the atom name is 4 characters long
        #if len(atmName)==4:
        if len(atmName)>=4:
            if len(atmName)>4:
                atmName = atmName[:4]
            # do something if it is an H
            if atm.element == 'H':
                atmName = atmName[-1]+atmName[:-1]
                rec = rec + "%-4.4s"%atmName
            #else:
            #    rec = rec + "%4.4s"%atmName
            elif len(atm.element)==2:
                rec = rec + "%-4.4s"%atmName
            else:
                rec = rec + ' %-3s'%atmName[:-1]
        elif len(atm.element)==2: # the Ca i.e else of if(atm.name)==4
            rec = rec + '%-4.4s'%atm.name   #7/8/05. why was this changed here?
            #rec = rec + '%-4.4s'%atm.element
        else:
            rec = rec + ' %-3s'%atmName
        # [16] Alternate location indicator
        rec = rec + altLoc
        # get the res name, res sequence number and the chain id.
        resName = ''
        resSeq = ''
        chainID = ''
        resIcode = ''
        if hasattr(atm, 'parent'):
            if hasattr(atm.parent, 'type'):
                resName = atm.parent.type
            if hasattr(atm.parent, 'number'):
                resSeq = atm.parent.number
            if hasattr(atm.parent, 'icode'):
                resIcode = atm.parent.icode
            if hasattr(atm.parent, 'parent') and hasattr(atm.parent.parent, 'id'):
                chainID = atm.parent.parent.id
        # [17:19] Residue name:
        rec = rec + "%3.3s"%resName
        # [20] space
        rec = rec + ' '
        # [21] Chain identifier.
        rec = rec + "%1.1s"%chainID
        # [22:25] residue sequence number
        rec = rec + '%4.4s'%resSeq
        # [26] res icode
        rec = rec + '%1.1s'%resIcode
        # [27:29] 3 empty spaces
        rec = rec + "%3.3s"%''

        coords = atm.coords
        # [30:37] x coordinate
        rec = rec + "%8.3f"%coords[0]
        # [38:45} y coordinate
        rec = rec + "%8.3f"%coords[1]
        # [46:53] z coordinate
        rec = rec + "%8.3f"%coords[2]
        # [54:59] occupancy
        if hasattr(atm, 'occupancy'):
            occupancy = atm.occupancy
        else:
            occupancy = 0.0
        rec = rec + '%6.2f'%occupancy
        # [60:65] Temperature factor
        if hasattr(atm, 'temperatureFactor'):
            tf = atm.temperatureFactor
        else:
            tf = 0.0
        rec = rec + '%6.2f'%tf
        # [66:70] 4 empty spaces
        rec = rec + '%4.4s'%''
        
        # charges
        if hasattr(atm, 'charge'):
            rec = rec + '%6.3f'%atm.charge
        elif hasattr(atm, 'gast_charge'):
            rec = rec + '%6.3f'%atm.gast_charge
        else:
            rec = rec + '%6s'%''

        #5/15:corrected the order 
        rec = rec + '  % 6.2f' % atm.AtVol
        rec = rec + '  % 6.2f' % atm.AtSolPar
        rec = rec + '\n'
        return rec



class PdbqtWriter(PdbqWriter):
    """Class to write data records from a molecule tree to a pdbq file.
    Has methods for the user to add own records and to write the record."""


    def __init__(self):
        """Constructor:
        """
        PdbWriter.__init__(self)


    def defineATOM_HETATMRecord(self, atm):
        """
        Define the ATOM or HETATM rec for the given atm
        """
        # Here the column are 0 based while in the PDB description it is 1 based.
        # [0:5] Record name
        if atm.hetatm==0: rec = "ATOM  "
        else: rec = "HETATM"
        # [6:11] Atom serial number + a space
        rec = rec+ '%5i '%atm.number
        # [12:15] Atom atmName
        # get the atmName
        atmName = atm.name
        # Separate the alternate location indicator from the name
        if '@' in atmName:
            atmName, altLoc = atmName.split('@')
        else:
            altLoc = ' '
        # Creating the atom name field
        # If the atom name is 4 characters long
        #if len(atmName)==4:
        if len(atmName)>=4:
            if len(atmName)>4:
                atmName = atmName[:4]
            # do something if it is an H
            if atm.element == 'H':
                atmName = atmName[-1]+atmName[:-1]
                rec = rec + "%-4.4s"%atmName
            elif len(atm.element)==2:
                rec = rec + "%-4.4s"%atmName
            else:
                rec = rec + ' %-3s'%atmName[:-1]
        elif len(atm.element)==2: # the Ca i.e else of if(atm.name)==4
            #rec = rec + '%-4.4s'%atm.element
            rec = rec + '%-4.4s'%atm.name   #7/8/05. why was this changed here?
        else:
            rec = rec + ' %-3s'%atmName
        # [16] Alternate location indicator
        rec = rec + altLoc
        # get the res name, res sequence number and the chain id.
        resName = ''
        resSeq = ''
        chainID = ''
        resIcode = ''
        if hasattr(atm, 'parent'):
            if hasattr(atm.parent, 'type'):
                resName = atm.parent.type
            if hasattr(atm.parent, 'number'):
                resSeq = atm.parent.number
            if hasattr(atm.parent, 'icode'):
                resIcode = atm.parent.icode
            if hasattr(atm.parent, 'parent') and hasattr(atm.parent.parent, 'id'):
                chainID = atm.parent.parent.id
        # [17:19] Residue name:
        rec = rec + "%3.3s"%resName
        # [20] space
        rec = rec + ' '
        # [21] Chain identifier.
        rec = rec + "%1.1s"%chainID
        # [22:25] residue sequence number
        rec = rec + '%4.4s'%resSeq
        # [26] res icode
        rec = rec + '%1.1s'%resIcode
        # [27:29] 3 empty spaces
        rec = rec + "%3.3s"%''

        coords = atm.coords
        # [30:37] x coordinate
        rec = rec + "%8.3f"%coords[0]
        # [38:45} y coordinate
        rec = rec + "%8.3f"%coords[1]
        # [46:53] z coordinate
        rec = rec + "%8.3f"%coords[2]
        # [54:59] occupancy
        if hasattr(atm, 'occupancy'):
            occupancy = atm.occupancy
        elif hasattr(atm, '_charges') and atm._charges.has_key('pqr'):
            occupancy = atm._charges['pqr']
        else:
            occupancy = 0.0
        rec = rec + '%6.2f'%occupancy
        # [60:65] Temperature factor
        if hasattr(atm, 'temperatureFactor'):
            tf = atm.temperatureFactor
        elif hasattr(atm, 'pqrRadius'):
            tf = atm.pqrRadius
        else:
            tf = 0.0
        rec = rec + '%6.2f'%tf
        # [66:70] 4 empty spaces
        rec = rec + '%4.4s'%''
        
        # charges
        if hasattr(atm, 'charge'):
            rec = rec + '%6.3f'%atm.charge
        elif hasattr(atm, 'gast_charge'):
            rec = rec + '%6.3f'%atm.gast_charge
        else:
            rec = rec + '%6s'%''

        #5/19:
        rec = rec + ' %-2.2s'%atm.autodock_element
        #rec = rec + ' %-2.2s'%atm.autodock_element.upper()
##         #NB: write 'A' in element slot for aromatic carbons
##         if atm.autodock_element=='A':
##             #in this case, columns 78+79 are blanks
##             rec = rec + 'A  '
##         else:
##             #rec = rec + '%2.2s'%atm.element
##             #5/19:
##             #columns 78+79: autodock_element
##             rec = rec + '%s '%atm.autodock_element
##             #if atm.element!=atm.autodock_element:
##             #    #eg HD or NA or SA or OA, always 2 chars
##             #    rec = rec + '%s '%atm.autodock_element[1]
##             #else:
##             #    rec = rec + '  '
        rec = rec + '\n'
        return rec




if __name__=='__main__':
    from MolKit.protein import Protein
    from MolKit.pdbParser import PdbParser
    mol = Protein()
    mol.read('/tsri/pdb/struct/4tpi.pdb', PdbParser())
    writer = PdbWriter()
    writer.add_userRecord('REMARK', )
    writer.add_userRecord('TITLE ', [('', 'This is the title record\n')])
    writer.write('/home/ktchan/jumble.pdb', mol)
