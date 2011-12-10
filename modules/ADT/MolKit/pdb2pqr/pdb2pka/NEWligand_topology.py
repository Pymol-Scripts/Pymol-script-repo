#
# $Id: NEWligand_topology.py,v 1.2 2007/10/10 22:15:40 vareille Exp $
# PC 2005/09/23
# Get ligand topologies
#
import numpy.oldnumeric as Numeric
    
from sets import Set
from ligandclean.trial_templates import *
from ligandclean.lookuptable import *
from substruct import Algorithms
from types import *

def length(vector):
    # This function returns the length of vector
    import math
    sum=0.0
    for value in vector:
        sum=sum+math.pow(value,2)
    return math.sqrt(sum)


class get_ligand_topology:
    ### PC
    #
    # here we need to check if we have MOL2 file, then guess_atom_types
    # are not necessary
    #
    def __init__(self,lines,MOL2FLAG):
        #
        # Given the atoms, this routine tells you everything you want to know about
        # the ligand
        #
        #
        # Store the atoms
        #
        if MOL2FLAG == False:
            self.atoms={}
            import string
            for line in lines:
                split=string.split(line)
                name=split[0]
                self.atoms[name]={'coords':Numeric.array([float(split[1]),float(split[2]),float(split[3])])}
                self.atoms[name]['bonds']=[]
            #
            # Get the likely types from names
            #
            trivial_types=['N','O','C','H']
            for atom in self.atoms.keys():
                if atom[0] in trivial_types:
                    if atom[0]!='H': # Get rid of all the hydrogens
                        self.atoms[atom]['type']=atom[0]
                        self.atoms[atom]['sybylType']='Unknown'
            #
            # Get the bonds
            # First approximation: Anything closer than 2.0 A is bonded
            #
            self.dists={}
            for atom1 in self.atoms.keys():
                self.dists[atom1]={}
                for atom2 in self.atoms.keys():
                    if atom1==atom2:
                        continue
                    #
                    # Calculate the distance
                    #
                    self.dists[atom1][atom2]=length(self.atoms[atom1]['coords']-self.atoms[atom2]['coords'])
                    if self.dists[atom1][atom2]<2.0:
                        self.atoms[atom1]['bonds'].append(atom2)
            #
            # Count number of bonds to non-H atoms and guess atom type
            #
            bond_lengths={'C-C':[1.5,0.2]}
            #
            # Get the torsion angles
            #
            atoms=self.atoms.keys()
            atoms.sort()
            for atom in atoms:
                self.atoms[atom]['torsions']=self.get_torsions(atom)
            #
            # Produce the definition lines
            #
            self.lines=self.create_deflines()
            #
            # Now we try to guess the atom types
            #
            self.guess_atom_types()
        else:
            #
            # We have a mol2 file
            #
            LIG = lines
            self.atoms={}
            for line in lines:
                name = line.name
                #            self.atoms[name] = name
                self.atoms[name] = {'coords': Numeric.array([float(line.x),float(line.y),float(line.z)])}
                self.atoms[name]['sybylType'] = line.sybylType
                #
                # we don't have this information when coming from PDB!
                self.atoms[name]['lBondedAtoms'] = line.lBondedAtoms
                self.atoms[name]['lBonds'] = line.lBonds
                ###PC
                # one bond is lost!
                self.atoms[name]['bonds']=[]
                for BBonds in line.lBondedAtoms: #line.lBonds:
                    self.atoms[name]['bonds'].append(BBonds.name)
                ###PC
                # save the atomname & id
                self.atoms[name]['atomname'] = name
                self.atoms[name]['serial'] = line.serial
                #
                # USEFUL information
                # bonded heavy atoms
                self.atoms[name]['nbhvy'] = len([x for x in self.atoms[name]['lBondedAtoms'] if x.sybylType != "H"])
                # number of bonds (including hydrogens)
                self.atoms[name]['nbds'] = len(self.atoms[name]['lBonds'])
                # number of bonded hydrogens
                self.atoms[name]['nbhyd'] = self.atoms[name]['nbds']- self.atoms[name]['nbhvy']
                # element
                self.atoms[name]['ele'] = self.atoms[name]['sybylType'].split('.')[0]
            #
            # Get the torsion angles
            #
            atoms=self.atoms.keys()
            atoms.sort()
            for atom in atoms:
                self.atoms[atom]['torsions']=self.get_torsions(atom)
            self.lines=self.create_deflines()
        #
        # we have the sybylType in self.atoms!!!
        return


    #
    # -------
    #

    def guess_atom_types(self):
        #
        # Phase I
        # Loop over all atoms and count number of bonds
        # + determine their likely order (e.g. single, double, or triple)
        #
        ambs={}
        atoms=self.atoms.keys()
        for atom_name in atoms:
            bonds=self.atoms[atom_name]['bonds']
            atype=self.atoms[atom_name]['type']
            numbonds=0
            aromatic=None
            for bonded_atom in bonds:
                #
                # Get the bond order from the distance
                #
                bond_order=self.get_bond_order(atom_name,bonded_atom)
                if bond_order<4:
                    numbonds=numbonds+bond_order
                else:
                    aromatic=1
            self.atoms[atom_name]['sum_bondorder']=numbonds
            self.atoms[atom_name]['aromatic']=aromatic
        #
        # ok, now we have info on all atoms.
        # Solve ambiguities starting with the simplest case
        #
        valences={'C':4,'O':2,'N':3}
        for atom in self.atoms:
            print atom, self.atoms[atom]
        #
        # ok, now it gets hairy
        #
        print
        print 'Guessing sybyl atom types'
        for atom in self.atoms.keys():
            stype=None
            at=self.atoms[atom]
            #
            # Get the precalculated characteristics
            #
            number_of_bonds=len(at['bonds'])
            sum_of_bondorder=at['sum_bondorder']
            if at['type']=='C':
                #
                # Carbon
                #
                #
                # Table or sum of bond-order, number of atoms bound to
                #
                C_types={5:{3:'C.2'}} # Carboxylic acid typically
                C_types[4]={4:'C.3',3:'C.2',2:'C.2'} # sum of bondorder is 4
                C_types[3]={3:'C.3',2:'C.2'}
                C_types[2]={2:'C.3',1:'C.2'}
                C_types[1]={1:'C.3'}
                #
                # Get the bond order
                #
                stype=C_types[sum_of_bondorder][number_of_bonds]
            elif at['type']=='O':
                #
                # OXYGEN
                #
                # Table or sum of bond-order, number of atoms bound to
                #
                O_types={2:{2:'O.3',1:'O.2'},1:{1:'O.3'}}
                stype=O_types[sum_of_bondorder][number_of_bonds]
            else:
                pass
            self.atoms[atom]['sybylType']=stype 
            #print atom,stype
        #
        # Do some postchecks
        # - right now only for Carboxylic acids
        #
        for atom in self.atoms.keys():
            at=self.atoms[atom]
            if at['sybylType']=='C.2':
                #
                # See if we have two oxygens bound + an extra bond
                #
                if len(at['bonds'])==3:
                    Os=[]
                    for bound_atom_name in at['bonds']:
                        bound_atom=self.atoms[bound_atom_name]
                        if bound_atom['sybylType']=='O.2':
                            Os.append(bound_atom_name)
                    #
                    # If we had two O.2s, then change their hybridisation to O.3
                    #
                    if len(Os)==2:
                        for O in Os:
                            self.atoms[O]['sybylType']='O.3'
        #
        # All Done
        #
        atoms=self.atoms.keys()
        atoms.sort()
        print '\nFinal Sybyl type results'
        for atom in atoms:
            print atom,self.atoms[atom]['sybylType']
        return

    #
    # --------
    #

    def get_bond_order(self,atom1,atom2):
        #
        # Get the bond order
        #
        # Bond lengths from
        # http://www.chem.swin.edu.au/modules/mod2/bondlen.html
        # We should get a better reference
        #
        # Returns:
        # 1: single bond, 2: double bond, 3: triple bond, 4: aromatic
        #
        #
        at1=self.atoms[atom1]
        at2=self.atoms[atom2]
        bond_props={'C-C':[1.54,1],'C=C':[1.34,2],'CtC':[1.20,3],'CaC':[1.40,4],
                    'C-O':[1.43,1],'C=O':[1.21,2],
                    'C-N':[1.47,1],'C=N':[1.25,2],'CtN':[1.16,3],'CaN':[1.34,4],
                    'NaN':[1.35,4]}
        dist=length(at1['coords']-at2['coords'])
        tps=[at1['type'],at2['type']]
        tps.sort() # To agree with dictionary layout
        best_fit=2.00
        best_type=None
        for bond in bond_props.keys():
            if bond[0]==tps[0] and bond[-1]==tps[1]:
                if abs(dist-bond_props[bond][0])<best_fit:
                    best_fit=abs(dist-bond_props[bond][0])
                    best_type=bond
        #
        # convert to bond order
        #
        return bond_props[best_type][1]

    #
    # ----
    #

    def get_torsions(self,start_atom):
        #
        # Get the torsion angles that start with this atom
        #
        #print '---------------------'
        #print 'Starting atom',start_atom
        possible_torsions=[]
        for bonded1 in self.atoms[start_atom]['bonds']:
            for bonded2 in self.atoms[bonded1]['bonds']:
                if bonded2==start_atom:
                    continue
                for end_atom in self.atoms[bonded2]['bonds']:
                    if end_atom==bonded1:
                        continue
                    #
                    # Add the torsion
                    #
                    possible_torsions.append([start_atom,bonded1,bonded2,end_atom])
        #
        # Filter the torsions
        #
        #
        # commented out by Paul 010206
#        print 'Jens has to write the stuff for filtering torsions..\n'
        return possible_torsions
                        
    #
    # ---------
    #

    def create_deflines(self):
        #
        # Make the lines for the pdb2pqr definition
        #
        self.numbers={}
        atoms=self.atoms.keys()
        atoms.sort()
        number=0
        for atom in atoms:
            number=number+1
            self.numbers[atom]=number
        #
        # Produce the lines
        #
        lines=[]
        for atom in atoms:
            lines.append('%s     %.2f  %.2f  %.2f' %(atom,self.atoms[atom]['coords'][0],
                                                     self.atoms[atom]['coords'][1],
                                                     self.atoms[atom]['coords'][2]))
        #
        # Bonds
        #
        bonds=''
        for atom in atoms:
            for bond in self.atoms[atom]['bonds']:
                start_num=self.numbers[atom]
                end_num=self.numbers[bond]
                #
                # Only write bonds one way (small number -> big number)
                #
                if end_num>start_num:
                    bonds=bonds+'%d %d ' %(start_num,end_num)
        lines.append(bonds)
        #
        # Torsions
        #
        tors=''
        written=0
        for atom in atoms:
            for torsion in self.atoms[atom]['torsions']:
                atom1=self.numbers[torsion[0]]
                atom2=self.numbers[torsion[1]]
                atom3=self.numbers[torsion[2]]
                atom4=self.numbers[torsion[3]]
                if atom1<atom2 and atom2<atom3 and atom4<atom4:
                    tors=tors+'%d %d %d %d ' %(atom1,atom2,atom3,atom4)
                    written=written+1
        lines.append('%d %s' %(written,tors))
        return lines
                    

    #
    # ---------
    #
    def ring_detection(self,start_atom,already_visited=[],level=0):
        if start_atom in already_visited and len(already_visited) >= 2:
            if start_atom==already_visited[-2]:
                return []
            if already_visited[0]!=start_atom:
                return []
            return already_visited+[start_atom]
        #
        return_lists=[]
        for bonded_atom in self.atoms[start_atom]['bonds']:
            this_list=already_visited[:]+[start_atom]
            this_list=self.ring_detection(bonded_atom,this_list,level+1)
            if this_list!=[]:# and len(this_list)>1:
                return_lists.append(this_list)
        return return_lists

        #
    # ------
    #

    def get_items(self,item):
        #
        # Reformat the lists of lists of lists of ... that we get from the ring detection
        #

        if type(item) is ListType:
            real_list=[]
            for sub_item in item:
                if not type(sub_item) is ListType:
                    real_list.append(sub_item)
                else:
                    self.get_items(sub_item)
            #
            # If we got something in the real_list then add it to the biglist
            #
            if real_list!=[]:
                self.biglist.append(real_list)
            return 
        else:
            raise 'this should not happen'

    #
    # -----
    #

    def assignRingAttribute(self,ring,atoms,current_atom):
        if self.atoms[current_atom]['in_ring'] != 0:
           self.atoms[current_atom]['in_ring'] += +1
        else:
            self.atoms[current_atom]['in_ring'] = 1
        return

                 

    def find_titratable_groups(self):
        #
        # Look for simple substructures that would be titratable groups in the ligand
        #
        atoms=self.atoms.keys()
        #
        # ring detection (including deleting redundancies & sorting issues)
        ring_list = []
        tmp=[]
        for atom in self.atoms.keys():
            temp_ring_list = []
            tmp.append(self.ring_detection(atom))
        #
        # Get just a single list of lists
        self.biglist=[]
        self.get_items(tmp)
        # bigList = self.biglist
        ring_list=self.biglist
        sorted_ring_list = []
        for rring in ring_list:
            rring = rring[:-1]
            rring.sort()
            sorted_ring_list.append(rring)
        sorted_ring_list.sort()   
        #
        # delete ring redundancies - only if ring present
        if len(sorted_ring_list) > 0:
            last = (sorted_ring_list)[-1]
            for i  in range(len(sorted_ring_list)-2,-1,-1):
                if last == sorted_ring_list[i]:
                    del sorted_ring_list[i]
                else:
                    last = sorted_ring_list[i]
        #print "# overall rings (including potentially fused rings) :", len(sorted_ring_list)
        #
        #
        # assigning ring attribute for every ring atom
        for at in self.atoms.keys():
            self.atoms[at]['in_ring'] = 0
        for rring in sorted_ring_list:
            for current_atom in rring:
                self.assignRingAttribute(rring,atoms,current_atom)
        # new attribute for each ring atom: appending the complete ring
        # atom names to which the atom belongs
        for atom in atoms:
            at = self.atoms[atom]
            at['ring_list'] = []
        non_fused_counter = 0
        for rring in sorted_ring_list:
            already_detected_false = False
            for atom in rring:
                at = self.atoms[atom]
                if already_detected_false == False and at['in_ring'] == 1:
                    non_fused_counter += 1
                already_detected_false = True
                at = self.atoms[atom]
                if at['ring_list'] == []:
                      at['ring_list'] = [rring]
                elif rring not in at['ring_list']:
                      at['ring_list'].append(rring)
        #print  "# non-fused rings                                   :", non_fused_counter

        ## LET'S START THE MATCHING...
        def match(t,l):
            class Node:
                def __init__(self, idx1, idx2):
                    self.at_idx1 = idx1
                    self.at_idx2 = idx2
            nodes = []
            for ligand_atoms in l:
                self.atoms[ligand_atoms]['alreadyvisited'] = False
            for ligand_atoms in l:
                for template_atoms in t:
                    if templates[current_template][template_atoms]['sybylType'] == self.atoms[ligand_atoms]['sybylType']:
                        matched_type = Node(ligand_atoms,template_atoms)
                        nodes.append(matched_type)
                        #print "matching ", self.atoms[ligand_atoms]['sybylType'], self.atoms[ligand_atoms]['atomname']
                    else:
                        pass
            if len(nodes) != 0:
                AtomNameListList = []
                for i in nodes:
                    for j in nodes:
                        AtomNameList = []
                        # the O.co2 atoms of the ligand can only match ONCE on EACH template O.co2 template atom
                        if i.at_idx1 == j.at_idx1 and (i.at_idx2 == j.at_idx2) and (self.atoms[i.at_idx1]['sybylType'] != "O.co2"):
                            AtomNameList.append(i.at_idx1)
                            AtomNameList.append(i.at_idx2)
                            AtomNameListList.append(AtomNameList)
                        elif i.at_idx1 == j.at_idx1 and (i.at_idx2 == j.at_idx2) \
                                 and (self.atoms[i.at_idx1]['sybylType'] == "O.co2")\
                                 and (self.atoms[i.at_idx1]['alreadyvisited'] != True)\
                                 and (templates[current_template][i.at_idx2]['alreadyvisited'] != True):
                            AtomNameList.append(i.at_idx1)
                            AtomNameList.append(i.at_idx2)
                            self.atoms[i.at_idx1]['alreadyvisited'] = True
                            templates[current_template][i.at_idx2]['alreadyvisited'] = True
                            AtomNameListList.append(AtomNameList)
                AGL = Numeric.zeros((len(AtomNameListList),len(AtomNameListList)))
                daseinecounter = 0
                dasanderecounter = 0
                # what again are these counters counting?
                for daseine in AtomNameListList:
                    dasanderecounter = 0
                    for dasandere in AtomNameListList:
                        if (daseine[0] == dasandere[0]) and (daseine[1] == dasandere[1]):
                            AGL[daseinecounter][dasanderecounter] = 1
                            # bonded
                        elif dasandere[0] in self.atoms[daseine[0]]['bonds'] and\
                           dasandere[1] in templates[current_template][daseine[1]]['neighbours']:
                            AGL[daseinecounter][dasanderecounter] = 1
                        # non-bonded atoms in template must also be non-bonded in the ligand
                        elif dasandere[0] not in self.atoms[daseine[0]]['bonds'] and\
                           dasandere[1] not in templates[current_template][daseine[1]]['neighbours']:
                            AGL[daseinecounter][dasanderecounter] = 1
                        dasanderecounter += 1
                    daseinecounter += 1
                allcliques = Algorithms.find_cliques(AGL)
                return allcliques,AtomNameListList,current_template
            else:
                return [],[],None

        # Looping over the entries of the template entries
        AllCliqueList = []
        dictCounter = 0
        dict_of_matched_lig_fragments = {}
        matched_lig_template = []
        for current_template in templates.keys():
            #print "MATCHING", current_template
            #print "######## ########"
            output = match(templates[current_template],atoms)
            for ee in output[0]:
                templateDoubleCounter = 0
                templateList = []
                ligList = []
                goodList = []
                for cliq_ee in ee:
                    templateList.append(output[1][cliq_ee][1])
                    ligList.append(output[1][cliq_ee][0])
                temp_templateList = list(Set(templateList))
                matchedLigAtoms =[]
                temptemp = []
                if len(temp_templateList) == len(ligList) and (len(ee) == len(templates[current_template].keys())):
                    for xxee in ee:
                        # output[1][xxee] => maching pair of ligand and template atoms
                        # temporary depostion
                        temptemp.append(output[1][xxee])
                        matchedLigAtoms.append(output[1][xxee][0])
                    AllCliqueList.append(matchedLigAtoms)
                    dict_of_matched_lig_fragments[dictCounter] = {}
                    dict_of_matched_lig_fragments[dictCounter]['matchedligatoms'] = matchedLigAtoms
                    dict_of_matched_lig_fragments[dictCounter]['templatename'] = current_template
                    dict_of_matched_lig_fragments[dictCounter]['type'] = templates_attributes[output[2]]['type']
                    dict_of_matched_lig_fragments[dictCounter]['modelpka'] = templates_attributes[output[2]]['modelpka']
                    #print "ligList",ligList
                    #print "dict",templates_attributes[output[2]]['titratableatoms']
                    #
                    # detection of titratable atoms by comparison with template_dict
                    deposit_list = []
                    # that's for the hydrogen injection
                    # in the template library, STANDARD are being used
                    deposit_template_list = []
                    dict={}
                    for looping_template in templates_attributes[output[2]]['titratableatoms']:
                        for looping_lig in output[1]:
                            #print looping_lig
                            if looping_lig[1] == looping_template and (looping_lig in temptemp) :
                                deposit_list.append(looping_lig[0])
                                deposit_template_list.append(looping_lig[1])
                                dict[looping_lig[1]]=looping_lig[0]
                                #print "tt", looping_lig[1],"  ll  ", looping_lig[0],looping_lig
                    dict_of_matched_lig_fragments[dictCounter]['titratableatoms']=deposit_list
                    dict_of_matched_lig_fragments[dictCounter]['matching_atoms']=dict
                    dictCounter += 1
        print
        print
        
        #if len(AllCliqueList) == 0:
        #    print
        #    print
        #    print "NOTHING CLIQUED."
        #    print
        #    print
        #
        # AVOID REDUNDANCIES FROM HERE ON
        counter = 0
        NonRedundantCliques = []
        # looping over all entries of maximum_cliques
        if len(AllCliqueList) > 1:
            for xxxx in AllCliqueList:
                if counter < len(AllCliqueList):
                    internalCounter = counter
                    # comparison of an individual clique with all other cliques
                    for comparing in range(internalCounter,len(AllCliqueList)):
                        # do not compare a clique with itself
                        if xxxx != AllCliqueList[internalCounter]:
                            if len(Set(xxxx).intersection(Set(AllCliqueList[internalCounter]))) > 0:
                                if len(xxxx) >= len(AllCliqueList[internalCounter]) and xxxx not in NonRedundantCliques:  # >= INSTEAD of >
                                                                                                                          # piperidine case!
                                    # avoid that subset is added to this list
                                    if len(NonRedundantCliques) != 0:
                                        # loop over all entries
                                        for possiblyredundantentries in NonRedundantCliques:
                                            if Set(possiblyredundantentries).issubset(Set(xxxx)):
                                                print NonRedundantCliques
                                                NonRedundantCliques.remove(possiblyredundantentries)
                                                NonRedundantCliques.append(xxxx)
                                                print NonRedundantCliques
                                            elif Set(xxxx).issubset(Set(possiblyredundantentries)):
                                                #print "found subset which is not added to the list"
                                                pass
                                    else:
                                        NonRedundantCliques.append(xxxx)
                            # the other way around
                                elif len(AllCliqueList[internalCounter]) >=  len(xxxx) and AllCliqueList[internalCounter] not in NonRedundantCliques:
                                    if len(NonRedundantCliques) != 0:
                                        for possiblyredundantentries in NonRedundantCliques:
                                            if Set(possiblyredundantentries).issubset(Set(AllCliqueList[internalCounter])):
                                                NonRedundantCliques.remove(possiblyredundantentries)
                                                NonRedundantCliques.append(AllCliqueList[internalCounter])
                                            elif Set(AllCliqueList[internalCounter]).issubset(Set(possiblyredundantentries)):
                                                #print "found subset which is not added to the list"
                                                pass
                                    else:
                                        NonRedundantCliques.append(AllCliqueList[internalCounter])
                        #
                        # how to adoid, that the actual matching (pair of ligand and template atoms) is not lost?
                        internalCounter += 1
                counter += 1
        # if we only find one titratable group
        elif len(AllCliqueList) == 1:
            NonRedundantCliques = AllCliqueList
        # redundancies should be removed now...
        #
       
         
        for allCl in dict_of_matched_lig_fragments:
            if dict_of_matched_lig_fragments[allCl]['matchedligatoms'] == NonRedundantCliques[0]:
                print "WE MATCHED", dict_of_matched_lig_fragments[allCl]['templatename']
                print "matchedligatoms            : ", dict_of_matched_lig_fragments[allCl]['matchedligatoms']
                print "type                       : ", dict_of_matched_lig_fragments[allCl]['type']
                print "modelpka                   : ", dict_of_matched_lig_fragments[allCl]['modelpka']
                print "titratableatoms            : ", dict_of_matched_lig_fragments[allCl]['titratableatoms']
                print "matching atoms             : ", dict_of_matched_lig_fragments[allCl]['matching_atoms']
        # re-run matching to get mutiple titratable sites?
        
        # TJD: This is to resolve the bug fix when allCl is None
        if dict_of_matched_lig_fragments != {}:
            print dict_of_matched_lig_fragments[allCl]
            return dict_of_matched_lig_fragments[allCl]
        else:
            return {}
