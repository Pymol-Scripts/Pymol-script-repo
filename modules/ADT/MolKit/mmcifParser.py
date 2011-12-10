# $Header: /opt/cvs/python/packages/share1.5/MolKit/mmcifParser.py,v 1.20.6.2 2011/05/19 19:20:56 sargis Exp $
#
# $Id: mmcifParser.py,v 1.20.6.2 2011/05/19 19:20:56 sargis Exp $

""" Module mmCIF_Parser. 

Implements MMCIFParser that builds dictionary (mmCIF_dict)
from cif files, and creates MolKit.molecule object from
this dictionary.

Example:
    #read and build dictionary
    parser = MMCIFParser( filename='Tests/Data/1A56.cif' )
    #create MolKil Protein object
    mol = parser.parse()
    #parse secondary structure information (optional)
    parser.parseSSData( mol )
"""
    
from os.path import splitext, basename
import types, sys
from string import split, strip, digits, lower, find
from MolKit.moleculeParser import MoleculeParser
from MolKit.protein import Protein, Chain, ChainSet, Residue, ResidueSet, ProteinSet
from MolKit.molecule import Atom, AtomSet, Bond, BondSet, HydrogenBond

class MMCIFParser(MoleculeParser):

    def __init__(self, filename=None, allLines=None):
        """Constructor for mmCIFParser: adopted form PdbParser"""
        #self.filename = filename
        #self.allLines = allLines
        MoleculeParser.__init__(self, filename, allLines)
#       self.altLoc = None #Flag to handle alternate locations.
#        self.model = 0 #Flag to indicate if Models or not

        self.mmCIF_dict = {} # dictionary used to store mmcif data structure
        #if filename:
        #    fptr = open( filename )
        #    self.allLines = fptr.readlines()
        #    fptr.close()
        #if not len( self.allLines ):
        #    print "ERROR in mmcifParser.py: no lines to parse!"
        #    return "ERROR"
        #self.mmCIF2Dict()
        

    def mmCIF2Dict(self):
        """Converts .cif file into dictionary (self.mmCIF_dict)  """
        allLines = self.allLines #
        len_allLines = len(allLines)
        i_line = 0 #i_line counts lines in allLines
        mmCIF_dict = self.mmCIF_dict
        mmCIF_dict['data_'] = []
        while i_line < len_allLines: 
            line = allLines[i_line]
            if line[:5] == 'data_': # check for the data block
                i_line += 1 
                mmCIF_dict['data_'].append(line[5:])
                continue
            elif line[0] == "_": # check for the data name tag
                tmp_list = line.split()
                if len( tmp_list ) == 1:
                    i_line += 1
                    tmp_string =  allLines[i_line]
                    if tmp_string[0] == ';': #multiline data value tag
                        i_line += 1          #collect all the lines untill
                        while not allLines[i_line][0] == ';': #hittinh another ";"
                            tmp_string += allLines[i_line].strip()
                            i_line += 1
                        tmp_string = tmp_string.split(';')[1]
                    mmCIF_dict[tmp_list[0]] = tmp_string
                else:
                    mmCIF_dict[tmp_list[0]] = ' '.join(tmp_list[1:])
                i_line += 1
                if i_line >= len_allLines: break
                continue
            elif line[0] == '#':
                i_line += 1
                continue #comments are ignored
            elif line[:5] == 'loop_': #parsing data structure in the loop
                i_line += 1
                tmp_dict = []
                for Line_ in allLines[i_line:]: #this loop collects data names and data values
                    Line_ = Line_.strip()       #in the temporary dictionary(tmp_dict) 
                    if Line_[0] == '_':  #First get data name tag which start with '_'
                        tmp_dict.append((Line_, []))
                        i_line += 1
                    else:  #now parse the rest of the loop which the data part
                        key_i = 0
                        while allLines[i_line][0] != '_' and allLines[i_line][:5]!= 'loop_':
                                tmp_string = allLines[i_line].strip()
                                if tmp_string == "":
                                    i_line += 1
                                if i_line >= len_allLines: break    
                                if tmp_string[0] == ';': #multiline data-value 
                                    i_line += 1
                                    while not allLines[i_line][0] == ';':
                                        tmp_string += allLines[i_line].strip()
                                        i_line += 1
                                    tmp_string = tmp_string.split(';')[1]
                                    tmp_dict[key_i][1].append(tmp_string)
                                    key_i += 1
                                    i_line += 1 
                                    if i_line >= len_allLines: break
                                    if key_i == len(tmp_dict):
                                        key_i = 0
                                    continue                                   
                                elif tmp_string[0] == '#': #comments are ignored
                                    i_line += 1
                                    if i_line >= len_allLines: break
                                    continue 
                                tmp_list = tmp_string.split();
                                if len(tmp_list) > len(tmp_dict):
                                    print "WARNING!!! in mmcifParser.py" 
                                    print "More data-values was provided than data-name"
                                for key_ii in range(len(tmp_list)):
                                    tmp_dict[key_i][1].append(tmp_list[key_ii])
                                    key_i += 1
                                    if key_i == len(tmp_dict):
                                        key_i = 0
                                        break
                                i_line += 1
                                if i_line >= len_allLines: break
                        new_dict = {}
                        for key, value in tmp_dict:
                            new_dict[key] = value
                        mmCIF_dict.update(new_dict)    
                        break
            else:
                i_line += 1 
        
        
    def parse(self, objClass=Protein):
        """Parses mmCIF dictionary (self.mmCIF_dict) into MolKit object"""
        if self.allLines is None and self.filename:
            self.readFile()
            if self.allLines is None or len(self.allLines)==0:
                return
            self.mmCIF2Dict()
        type_symbol = None
        B_iso_or_equiv = None
        mmCIF_dict = self.mmCIF_dict
        if mmCIF_dict.has_key('_atom_site.id'):
            #The description of the data names can be found in the following link
            #http://mmcif.pdb.org/dictionaries/mmcif_pdbx.dic/Items   
            ids = mmCIF_dict['_atom_site.id'] #1 number
            group_PDB = mmCIF_dict['_atom_site.group_PDB']          #2 atom/hetatm
            
            atom_id = mmCIF_dict['_atom_site.label_atom_id']  #3 name

            comp_id = mmCIF_dict['_atom_site.label_comp_id']  #4 residue type
            label_asym_id = mmCIF_dict['_atom_site.label_asym_id']  #5 chain 
            #Note: chain ID from mmCIF file might be different from PDB file
            seq_id = mmCIF_dict['_atom_site.label_seq_id']    #6 residue number
            x_coords = mmCIF_dict['_atom_site.Cartn_x']             #7 xcoord
            y_coords = mmCIF_dict['_atom_site.Cartn_y']             #8 ycoord
            z_coords = mmCIF_dict['_atom_site.Cartn_z']             #9 zcoord
            occupancy = mmCIF_dict['_atom_site.occupancy']          #10    
            B_iso_or_equiv = mmCIF_dict['_atom_site.B_iso_or_equiv']#11
            type_symbol = mmCIF_dict['_atom_site.type_symbol']
            molName = mmCIF_dict['_entry.id']
                
        elif mmCIF_dict.has_key('_atom_site_label'):
            #ftp://ftp.iucr.org/pub/cif_core.dic
            atom_id = mmCIF_dict['_atom_site_label']
            len_atoms = len(atom_id)
            ids = range(len_atoms)
            
            group_PDB = len_atoms*['HETATM']
            comp_id = len_atoms*["CIF"]
            label_asym_id = len_atoms*['1']
            seq_id = len_atoms*[1]
            
            from mglutil.math.crystal import Crystal
            a = mmCIF_dict['_cell.length_a'] = float(mmCIF_dict['_cell_length_a'].split('(')[0])
            b = mmCIF_dict['_cell.length_b'] = float(mmCIF_dict['_cell_length_b'].split('(')[0])
            c = mmCIF_dict['_cell.length_c'] = float(mmCIF_dict['_cell_length_c'].split('(')[0])
            alpha = mmCIF_dict['_cell.angle_alpha'] = float(mmCIF_dict['_cell_angle_alpha'].split('(')[0])
            beta = mmCIF_dict['_cell.angle_beta'] = float(mmCIF_dict['_cell_angle_beta'].split('(')[0])
            gamma = mmCIF_dict['_cell.angle_gamma'] = float(mmCIF_dict['_cell_angle_gamma'].split('(')[0])
            cryst = Crystal((a, b, c), (alpha, beta, gamma))
            x = []
            for item in mmCIF_dict['_atom_site_fract_x']:
                x.append(float(item.split('(')[0]))
            y = []
            for item in mmCIF_dict['_atom_site_fract_y']:
                y.append(float(item.split('(')[0]))
            z = []
            for item in mmCIF_dict['_atom_site_fract_z']:
                z.append(float(item.split('(')[0]))
                
            x_coords = []
            y_coords = []
            z_coords = []
            B_iso_or_equiv = []
            for i in ids:
                trans = cryst.toCartesian([x[i], y[i], z[i]])
                
                x_coords.append(trans[0]) 
                y_coords.append(trans[1])
                z_coords.append(trans[2])
                if mmCIF_dict.has_key('_atom_site_U_iso_or_equiv'):
                    B_iso_or_equiv.append(mmCIF_dict['_atom_site_U_iso_or_equiv'][i].split('(')[0])
            if mmCIF_dict.has_key('_atom_site_type_symbol'):
                type_symbol = mmCIF_dict['_atom_site_type_symbol']
            if mmCIF_dict.has_key('_atom_site_occupancy'):
                occupancy = mmCIF_dict['_atom_site_occupancy']
            if mmCIF_dict.has_key('_chemical_name_common'):   
                molName = mmCIF_dict['_chemical_name_common']
            elif mmCIF_dict.has_key('_chemical_name_mineral'):
                molName = mmCIF_dict['_chemical_name_mineral']
                                
            if mmCIF_dict.has_key('_symmetry_space_group_name_H-M'):   
                mmCIF_dict['_symmetry.space_group_name_H-M'] = mmCIF_dict['_symmetry_space_group_name_H-M']
        else:
            print 'No _atom_site.id or _atom_site_label record is available in %s' % self.filename
            return  None  
        
        mol = Protein()
        self.mol = mol
        self.mol.allAtoms = AtomSet([])
        molList = mol.setClass()
        molList.append( mol )
        current_chain_id = None
        current_residue_number = None
        current_chain = None
        current_residue = None
        
        number_of_atoms = len(ids)

        self.configureProgressBar(init=1, mode='increment', 
                                  authtext='parse atoms', max=number_of_atoms)
        for index in range(number_of_atoms):              
            #make a new atom for the current index
            chain_id = label_asym_id[index]
            if chain_id != current_chain_id:         #make a new chain
                #molecule should adopt the current chain if there is one
                current_chain = Chain(id=chain_id)
                # FIXME: current_chain should not have allAtoms attribute
                delattr(current_chain, "allAtoms")
                current_chain_id = chain_id
                
                if current_chain is not None:    #REMEMBER TO ADOPT THE LAST ONE!!!
                    mol.adopt(current_chain, setChildrenTop=1)                    
            residue_number = seq_id[index]   

            if residue_number != current_residue_number or chain_id != label_asym_id[index-1]:         #make a new chain:
                #current_chain should adopt the current residue if there is one
                #create new residue
                residue_type = comp_id[index]
                current_residue = Residue(type=residue_type, number=residue_number)
                current_residue_number = residue_number
                if current_residue is not None:    #REMEMBER TO ADOPT THE LAST ONE!!!
                    current_chain.adopt(current_residue, setChildrenTop=1)
                
            
            name = atom_id[index]
            if type_symbol:
                element = type_symbol[index]
            else:
                element = None
            atom = Atom( name, current_residue, element, top=mol )
            atom._coords = [[float(x_coords[index]), float(y_coords[index]), float(z_coords[index])]]
            atom._charges = {}
            atom.segID =  mol.name   
            atom.normalname = name
            atom.number = int(ids[index])
            mol.atmNum[atom.number] = atom
            atom.occupancy = float(occupancy[index])
            if B_iso_or_equiv:
                atom.temperatureFactor = float(B_iso_or_equiv[index])
            atom.altname = None    
            atom.hetatm = 0
            if group_PDB[index]=='HETATM':
                atom.hetatm = 1
            self.updateProgressBar()
                           
        self.parse_MMCIF_CELL()
        try:
            self.parse_MMCIF_HYDBND()       
        except:
             print >>sys.stderr,"Parsing Hydrogen Bond Record Failed in",self.filename
               
        mol.name = molName
        mol.allAtoms = mol.chains.residues.atoms
        
        mol.parser = self
        mol.levels = [Protein, Chain, Residue, Atom]
        name = ''
        for n in molList.name:
            name = n + ','
        name = name[:-1]
        molList.setStringRepr(name)
        strRpr = name + ':::'
        molList.allAtoms.setStringRepr(strRpr)
        for m in molList:
            mname = m.name
            strRpr = mname + ':::'
            m.allAtoms.setStringRepr(strRpr)
            strRpr = mname + ':'
            m.chains.setStringRepr(strRpr)
            for c in m.chains:
                cname = c.id
                strRpr = mname + ':' + cname + ':'
                c.residues.setStringRepr(strRpr)
                for r in c.residues:
                    rname = r.name
                    strRpr = mname + ':' + cname + ':' + rname + ':'
                    r.atoms.setStringRepr(strRpr)                            
        return molList
        

    def configureProgressBar(self, **kw):
        # this method is to be implemented by the user from outside
        pass
        
    def getMoleculeInformation(self):
        """Function to retrieve the general informations on the molecule.
        FIXME: Needs to be modified"""
        
        return self.mmCIF_dict['_entry.id']
        
    def updateProgressBar(self, progress=None):
        # this method is to be implemented by the user from outside
        #print 'Prorgess: ' + `progress` + '%'
        pass
   
    def parseSSData(self, mol):
        """
        Function to parse the information describing the secondary structure
        of the protein grouped as chain ID, the information is provided
        as a list of the following structure:
        [ ['chainID',[ Helix,[1stResHelix1,lastResHelix1], ...],
        [ Strand, [1stResSheet1,lastResSheet1] ],...],.... ]
        """
        from MolKit.protein import Helix, Strand, Turn

        # Step 1: Create a list containing the information describing the
        # the secondary structures organized the following way:
        # [ ['chain1ID', [Helix, [startHel1, endHel1],[startHel2, endHel2]],
        # [Strand, [startSheet1, endSheet1]] ], ['chain2ID', [Helix .....]] ]
        ssDataForMol = {}
        for chain in mol.chains:
            helStartEndForChain = self.processHelData(chain)
            if helStartEndForChain:
                helStartEndForChain.insert(0, Helix)
            else:
                helStartEndForChain = [0] 
            
            # fieldIndices
            # 5: chain ID, 4: RESNAME, 6: RES SEQNB, 7: RES INSER,
            # 8: RESNAME, 10: RES SEQNB, 11: RES INSER, 3: NB STRAND,
            # 12: STRAND SENSE
            strandStartEndForChain = self.processStrData(chain)
            if strandStartEndForChain:
                strandStartEndForChain.insert(0, Strand)
            else: 
                strandStartEndForChain = [0]                    

            # fieldIndices
            # 4: chain ID, 3: RESNAME, 5: RES SEQNB, 6: RES INSER,
            # 7: RESNAME, 9: RES SEQNB, 10: RES INSER,
            # 11: COMMENT
            turnStartEndForChain = self.processTurnData(chain)
            if turnStartEndForChain:
                turnStartEndForChain.insert(0, Turn)
            else:
                turnStartEndForChain = [0] 
                
            ssDataForMol[chain.id] = [helStartEndForChain, strandStartEndForChain, 
                                   turnStartEndForChain, None]
            
        return ssDataForMol

    def processHelData(self, chain):
        """ Processes Helix data record"""
        chain_id = chain.id
        mmCIF_dict = self.mmCIF_dict            
        #The description of the data names can be found in the following link
        #http://mmcif.pdb.org/dictionaries/mmcif_pdbx.dic/Items/_data_name.html
        try:
            conf_id = mmCIF_dict['_struct_conf.conf_type_id']
        except KeyError:      
            print  'No STRUCT_CONF category record is available in %s' % self.filename
            return  None

        beg_comp_id  = mmCIF_dict['_struct_conf.beg_label_comp_id']
        beg_asym_id  = mmCIF_dict['_struct_conf.beg_label_asym_id']
        beg_seq_id   = mmCIF_dict['_struct_conf.beg_label_seq_id']
        beg_PDB_ins_code  = mmCIF_dict['_struct_conf.pdbx_beg_PDB_ins_code'] 
        end_comp_id  = mmCIF_dict['_struct_conf.end_label_comp_id']
        end_asym_id  = mmCIF_dict['_struct_conf.end_label_asym_id']
        end_seq_id   = mmCIF_dict['_struct_conf.end_label_seq_id']
        end_PDB_ins_code  = mmCIF_dict['_struct_conf.pdbx_end_PDB_ins_code']
        
        helix_class       = mmCIF_dict['_struct_conf.pdbx_PDB_helix_class']
        conf_details      = mmCIF_dict['_struct_conf.details']
        number_of_records = len(conf_id)
        
        helStartEndData = []
        for index in range(number_of_records):
            if beg_asym_id[index] == chain.id and conf_id[index].find( 'HELX' ) == 0:
                startData = beg_comp_id[index]+ beg_PDB_ins_code[index].strip( '?' )\
                            +  beg_seq_id[index]
                startData = startData.strip()
                startRes = chain.residues.get(startData)

                if len(startRes) != 1:
                    print ("ERROR: When parsing the HELIX information found \
%d %s in chain %s"%(len(startRes), startData, chain.id))
                    continue
                endData = end_comp_id[index]+ end_PDB_ins_code[index].strip('?')\
                            +  end_seq_id[index]
                endData = endData.strip()
                endRes = chain.residues.get(endData)

                if len(endRes) != 1:
                    print ("ERROR: When parsing the HELIX information found \
%d %s in chain %s"%(len(endRes), endData, chain.id))
                    continue
                helClass = helix_class[index]
                comment = conf_details[index]

                
                helStartEndData.append({'start':startRes[0], 'end':endRes[0], 
                                    'helClass':helClass, 'comment':comment})
        
        return helStartEndData
                                                 
    def processStrData(self, chain):
        """ Processes Sheet data record"""
        chain_id = chain.id
        mmCIF_dict = self.mmCIF_dict
        strStartEndData = []
        try:
            sheet_range_id = mmCIF_dict['_struct_sheet_range.sheet_id']
        except KeyError:      
            print  'No STRUCT_SHEET category record is available in %s' % self.filename
            return  None

        range_id = mmCIF_dict['_struct_sheet_range.id']
        beg_comp_id = mmCIF_dict['_struct_sheet_range.beg_label_comp_id']
        beg_asym_id = mmCIF_dict['_struct_sheet_range.beg_label_asym_id']
        beg_seq_id  = mmCIF_dict['_struct_sheet_range.beg_label_seq_id']
        beg_PDB_ins_code  = mmCIF_dict['_struct_sheet_range.pdbx_beg_PDB_ins_code'] 
        end_comp_id = mmCIF_dict['_struct_sheet_range.end_label_comp_id']
        end_asym_id = mmCIF_dict['_struct_sheet_range.end_label_asym_id']
        end_seq_id  = mmCIF_dict['_struct_sheet_range.end_label_seq_id']
        end_PDB_ins_code  = mmCIF_dict['_struct_sheet_range.pdbx_end_PDB_ins_code']
        # The following 4 lists are needed to get the sense of the strand
        order_sheet_id    = mmCIF_dict['_struct_sheet_order.sheet_id']
        order_range_id_1  = mmCIF_dict['_struct_sheet_order.range_id_1']
        order_range_id_2  = mmCIF_dict['_struct_sheet_order.range_id_2']
        order_sense       = mmCIF_dict['_struct_sheet_order.sense']
        
        number_of_records = len(sheet_range_id)
 
 
        #This is the main loop where we go through each record and extract
        #strStartEndData: similar to what has been done in PdbParser
        for index in range(number_of_records):
            if beg_asym_id[index] == chain.id:
                startData = beg_comp_id[index] + beg_PDB_ins_code[index].strip('?')\
                            +  beg_seq_id[index]
                startData = startData.strip()
                startRes = chain.residues.get(startData)

                if len(startRes) != 1:
                    print ("ERROR: When parsing the SHEET information found \
%d %s in chain %s"%(len(startRes), startData, chain.id))
                    continue

                endData = end_comp_id[index]+ end_PDB_ins_code[index].strip('?')\
                            +  end_seq_id[index]
                endData = endData.strip()
                endRes = chain.residues.get(endData)

                if len(endRes) != 1:
                    print ("ERROR: When parsing the SHEET information found \
%d %s in chain %s"%(len(endRes), endData, chain.id))
                    continue
                
                nbStrand = len(filter(lambda x: x[0] == sheet_range_id[index], \
                               sheet_range_id))
                if type( order_sheet_id ) == types.ListType:
                    for tmp_index in range(len(order_sense)):
                        if sheet_range_id[index] == order_sheet_id[tmp_index]:
                            if order_range_id_1[tmp_index] == range_id[index] or \
                               order_range_id_2[tmp_index] == range_id[index]:
                                   if order_range_id_1[tmp_index] == '1':
                                       sense = 0 #NOTE: 0 if first strand in PDB
                                   if order_sense[tmp_index] == 'parallel':
                                       sense = 1
                                   elif order_sense[tmp_index] == 'anti-parallel':
                                       sense = -1
                else:
                    if order_range_id_1 == '1':
                        sense = 0 #NOTE: 0 if first strand in PDB
                    if order_sense == 'parallel':
                        sense = 1
                    elif order_sense == 'anti-parallel':
                        sense = -1
                                           
                strStartEndData.append({'start':startRes[0], 'end':endRes[0], 
                                    'nbStrand':nbStrand, 'sense':sense})
        return strStartEndData

    def processTurnData(self, chain):
        """ Processes Turn data record"""
        chain_id = chain.id
        mmCIF_dict = self.mmCIF_dict            
        try:
            conf_id = mmCIF_dict['_struct_conf.id']
        except KeyError:      
            #print  'No STRUCT_CONF category record is available in %s' % self.filename
            return None
        conf_type_id = mmCIF_dict['_struct_conf.conf_type_id']
        beg_comp_id = mmCIF_dict['_struct_conf.beg_label_comp_id']
        beg_asym_id = mmCIF_dict['_struct_conf.beg_label_asym_id']
        beg_seq_id  = mmCIF_dict['_struct_conf.beg_label_seq_id']
        beg_PDB_ins_code  = mmCIF_dict['_struct_conf.pdbx_beg_PDB_ins_code'] 
        end_comp_id = mmCIF_dict['_struct_conf.end_label_comp_id']
        end_asym_id = mmCIF_dict['_struct_conf.end_label_asym_id']
        end_seq_id  = mmCIF_dict['_struct_conf.end_label_seq_id']
        end_PDB_ins_code  = mmCIF_dict['_struct_conf.pdbx_end_PDB_ins_code']
        helix_class       = mmCIF_dict['_struct_conf.pdbx_PDB_helix_class']
        conf_details      = mmCIF_dict['_struct_conf.details']
        number_of_records = len(conf_id)
        
        turnStartEndData = []
        #This is the main loop where we go through each record and extract
        #turnStartEndData: similar to what has been done in PdbParser
        for index in range(number_of_records):
            if beg_asym_id[index] == chain.id and conf_id[index].find('TURN') == 0:
                startData = beg_comp_id[index] + beg_PDB_ins_code[index].strip('?')\
                            +  beg_seq_id[index]
                startData = startData.strip()
                startRes = chain.residues.get(startData)
                if len( startRes ) != 1:
                    print ("ERROR: When parsing the TURN information found \
%d %s in chain %s"%(len(startRes), startData, chain.id))
                    continue
                endData = end_comp_id[index] + end_PDB_ins_code[index].strip('?')\
                            +  end_seq_id[index]
                endData = endData.strip()
                endRes = chain.residues.get(endData)
                if len(endRes) != 1:
                    print ("ERROR: When parsing the TURN information found \
%d %s in chain %s"%(len(endRes), endData, chain.id))
                    continue

                comment = conf_details[index]
                turnStartEndData.append({'start':startRes[0], 'end':endRes[0], 
                                    'comment':comment})
        return turnStartEndData

    def addModelToMolecules(self, listOfMol):
        length = len(listOfMol)
        for i in xrange(length):
            listOfMol[i].model = ProteinSet()
            for j in xrange(length):
                if listOfMol[i]!= listOfMol[j]:
                    listOfMol[i].model.append(listOfMol[j])
    def hasSsDataInFile(self):
        """ Function testing if the informations on the secondary structure
        are in the file"""
        test = filter(lambda x: x in self.mmCIF_dict, ['_struct_conf.id', '_struct_sheet_range.id'])
        if test: return 1
        else: return 0 

    def parse_MMCIF_CELL(self):
        """Parse the CELL category record. Create the following members:
        cellLength, cellAngles, spaceGroup, Zvalue"""
        mmCIF_dict = self.mmCIF_dict            
        try:
            mmCIF_dict['_cell.length_a']
        except KeyError:      
            print  'No CELL category record is available in %s' % self.filename
            return
        a = mmCIF_dict['_cell.length_a']
        b = mmCIF_dict['_cell.length_b']
        c = mmCIF_dict['_cell.length_c']
        alpha = mmCIF_dict['_cell.angle_alpha']
        beta = mmCIF_dict['_cell.angle_beta']
        gamma = mmCIF_dict['_cell.angle_gamma']
        self.mol.cellLength = [ a, b, c ]
        self.mol.cellAngles = [ alpha, beta, gamma ]
        self.mol.spaceGroup = mmCIF_dict['_symmetry.space_group_name_H-M'][1:-1]
        try:
            self.mol.Zvalue = mmCIF_dict['_cell.Z_PDB']
        except:
            self.mol.Zvalue = None
            
    def parse_MMCIF_HYDBND(self):
        """Parse the HYDBND record. Create the hbond described in
        that record by finding dAt, hAt and aAt, the donor, hatom and
        acceptorAtoms respectively."""
        mmCIF_dict = self.mmCIF_dict
        try:
            struct_conn_id = mmCIF_dict['_struct_conn.conn_type_id']
        except KeyError:      
            print  'No STRUCT_CONN category record is available in %s' % self.filename
            return
        ptnr1_asym_id = mmCIF_dict['_struct_conn.ptnr1_label_asym_id']
        ptnr1_atom_id = mmCIF_dict['_struct_conn.ptnr1_label_atom_id']
        ptnr1_alt_id  = mmCIF_dict['_struct_conn.pdbx_ptnr1_label_alt_id']
        ptnr1_comp_id = mmCIF_dict['_struct_conn.ptnr1_label_comp_id']
        ptnr1_seq_id  = mmCIF_dict['_struct_conn.ptnr1_label_seq_id']
        ptnr2_asym_id = mmCIF_dict['_struct_conn.ptnr2_label_asym_id']
        ptnr2_atom_id = mmCIF_dict['_struct_conn.ptnr2_label_atom_id']        
        ptnr2_alt_id  = mmCIF_dict['_struct_conn.pdbx_ptnr2_label_alt_id']
        ptnr2_comp_id = mmCIF_dict['_struct_conn.ptnr2_label_comp_id']
        ptnr2_seq_id  = mmCIF_dict['_struct_conn.ptnr2_label_seq_id']
        ptnr3_atom_id = mmCIF_dict['_struct_conn.pdbx_ptnr3_label_atom_id']
        ptnr3_seq_id  = mmCIF_dict['_struct_conn.pdbx_ptnr3_label_seq_id']
        ptnr3_asym_id = mmCIF_dict['_struct_conn.pdbx_ptnr3_label_asym_id']        
        
        number_of_records = len(struct_conn_id)

        for index in range (number_of_records):
            if struct_conn_id[index].find('hydrog') == 0:
                dAtName   = ptnr1_atom_id[index]
                dAtPType  = ptnr1_comp_id[index]
                dAtPNum   = ptnr1_seq_id[index]
                dAtPName  = dAtPType + str(dAtPNum)
                dAtPIcode = ptnr1_alt_id[index]
                if not dAtPIcode == '?':
                    dAtPName = dAtPName + dAtPIcode
                dAtChId = ptnr1_asym_id[index]
                dname = self.mol.name+':'+dAtChId+':'+dAtPName+':'+dAtName

                hAtName = ptnr3_atom_id[index]
                if len(hAtName):
                    if len(hAtName) == 4:
                        hAtName = hAtName[1:]+hAtName[0]
                    hAtChId = ptnr3_asym_id[index]
                    hAtPNum = ptnr3_seq_id[index]
                    hname = self.mol.name+':'+dAtChId+':'+dAtPName+':'+hAtName

                aAtName  = ptnr2_atom_id[index]
                aAtPType = ptnr2_comp_id[index]
                aAtPNum  = ptnr2_seq_id[index]
                aAtPName = aAtPType + str(aAtPNum)
                aAtPIcode = ptnr2_alt_id[index]
                if not aAtPIcode == '?':
                    aAtPName = aAtPName + aAtPIcode
                aAtChId = ptnr2_asym_id[index]
                aname = self.mol.name+':'+aAtChId+':'+aAtPName+':'+aAtName
                
                dAt = self.mol.allAtoms.get(lambda x: x.full_name()==dname)[0]
                aAt = self.mol.allAtoms.get(lambda x: x.full_name()==aname)[0]
                if len(hAtName):
                   hAt = self.mol.allAtoms.get(lambda x: x.full_name()==hname)[0]
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
if __name__ == '__main__':
    parser = MMCIFParser( filename='Tests/Data/1CRN.cif' )
    print "Reading molecule"
    mol = parser.parse()
    print "Done parsing"
    SS_Data  = parser.parseSSData( mol )
    print "Done parsing secondary structure"
    print "Done"
