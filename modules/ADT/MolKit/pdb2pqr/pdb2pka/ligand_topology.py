#
# $Id: ligand_topology.py,v 1.2 2007/10/10 22:15:40 vareille Exp $
# PC 2005/09/23
# Get ligand topologies
#
import numpy.oldnumeric as Numeric
    
from sets import Set
from ligandclean.trial_templates import *

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
        print 'Jens has to write the stuff for filtering torsions..\n'
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
        if start_atom in already_visited:
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
        from types import *
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

    #
    # ------
    #
        
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
        print "# overall rings (including potentially fused rings) :", len(sorted_ring_list)
        stop ## PC 03.01.06
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
        print  "# non-fused rings                                   :", non_fused_counter

 
        


    def matched_atom_types(atom2match,t):
        match_list=[]
            #
        # match ligand atom type with atom type from template
        for at in t.keys():
            if t[at]['sybylType'] == self.atoms[atom2match]['sybylType']:
                match_list.append(at)
            if len(match_list) != 0:
                return match_list,t[at]['sybyl_neighbours']
        if len(match_list) == 0:
            return None,None

        def match(t,l,already_visited=[],type_matches=[]):
            for counter in range(len(atoms)):
                at_lig = atoms[counter]
                already_visited.append(at_lig)
                # 1st matching: based on atom types
                matched_atom_in_template, nbs_in_template = matched_atom_types(at_lig,t)
                if matched_atom_in_template != None:
                    for entries in  matched_atom_in_template:
                        ligand_list = []
                        # Create sybyl_neighbors on-the-fly for ligand
                        for sybyl_bonded_at in self.atoms[at_lig]['lBondedAtoms']:
                            ligand_list.append(sybyl_bonded_at.sybylType)
                            ligand_set = Set(ligand_list)
                            template_set = Set(nbs_in_template)
                            # Now match simultaneously atom_type and neighbouring atom_types for ligand AND template
                            if len(ligand_set.difference(template_set)) == 0 and len(ligand_list) == len(nbs_in_template):
                                for entry in matched_atom_in_template:
                                    print "%3d"%(counter),"  Ligand %4s %5s %28s " \
                                          %(at_lig,self.atoms[at_lig]['sybylType'],ligand_list),\
                                          "template %s %s %s %s" \
                                          %(matched_atom_in_template,t[entry]['sybylType'],nbs_in_template,t[entry]['neighbours'])
                                    for neighboured_template_atoms in t[entry]['neighbours']:
                                        print neighboured_template_atoms,t[neighboured_template_atoms]['sybylType'],t[neighboured_template_atoms]['sybyl_neighbours']
                                    for neighboured_ligand_atoms in self.atoms[at_lig]['lBondedAtoms']:
                                        print neighboured_ligand_atoms.name, neighboured_ligand_atoms.sybylType,neighboured_ligand_atoms.lBondedAtoms
                                    stop
                counter += 1

        def matched_atom_types2(atom2match,t,stored_nbs_of_atom2match=[],already_visited=[],matching_template={}):
            #
            # match ligand atom type with atom type from template
            if atom2match == "F14":
                print "YYY_atom2match_YYY", atom2match
#            print "alrvis",len(already_visited),already_visited
            if matching_template == {}:
                matching_template['MatchedFragments'] = {}
            if len(stored_nbs_of_atom2match) != 0 and stored_nbs_of_atom2match[-1] == atom2match:
                print "bis zum erbrechen schreien!!!!", self.atoms[atom2match]['bonds']
                for e in  self.atoms[atom2match]['bonds']:
                    atom2match = e
            for at in t.keys():
                # TODO:matching ALL atom types in template => gives a match_list
                if t[at]['sybylType'] == self.atoms[atom2match]['sybylType'] \
                   and atom2match not in already_visited:
                    already_visited.append(self.atoms[atom2match]['atomname'])
                    Lig_nbs_SybylList = []
                    Lig_nbs_AtomnameList = []
                    # Create sybyl_neighbors on-the-fly for ligand
                    for att in self.atoms[atom2match]['lBondedAtoms']:
                        Lig_nbs_SybylList.append(att.sybylType)
                        Lig_nbs_AtomnameList.append(att.name)
                    ligand_set = Set(Lig_nbs_SybylList)
                    template_set = Set(t[at]['sybyl_neighbours'])
                    diff = ligand_set.difference(template_set)
                    if len(diff) == 0:
                        stored_nbs_of_atom2match = Lig_nbs_AtomnameList
                        matching_template['MatchedFragments'][atom2match] = {}
                        matching_template['MatchedFragments'][atom2match]['sybyl_neighbours'] = Lig_nbs_SybylList
                        # go through all bonded atoms
                        if len(stored_nbs_of_atom2match) != 0:
                            for bb in stored_nbs_of_atom2match:
                                if bb not in already_visited:
                                    already_visited.append(bb)
                                    matching_template['MatchedFragments'][bb] = {}
                                    bb_list = []
                                    for bat in self.atoms[bb]['lBondedAtoms']:
                                        bb_list.append(bat.sybylType)
                                    matching_template['MatchedFragments'][bb]['sybyl_neighbours'] = bb_list
                                    #
                                    # here we call the routine by itself
                                    matched_atom_types2(bb,t,stored_nbs_of_atom2match)
                    else: # NO MATCH
                        for nbat in self.atoms[atom2match]['bonds']:
                            if nbat in already_visited:
                                start_id = 1
                            for id in range(len(self.atoms[atom2match]['bonds'])-1):
                                next_nbat_id  = id+1
                                next_nbat_at = self.atoms[atom2match]['bonds'][next_nbat_id]
                                if next_nbat_at not in already_visited:
                                    #not 100% sure, if if append the correct atom
#                                    already_visited.append(next_nbat_at)
                                    already_visited.append(self.atoms[atom2match]['bonds'][next_nbat_id])
                                    matched_atom_types2(self.atoms[atom2match]['bonds'][next_nbat_id],t)
                                else:
                                    pass

                        else:
                            matched_atom_types2(nbat,t)
                else:
                    print "sybylType s don't match", atom2match
            # 2nd loop to go over to the neighboured atoms
            for at in t.keys():
                if atom2match in already_visited:
                    for nbat in self.atoms[atom2match]['bonds']:
                        if nbat in already_visited:
                            #
                            # TODO: Do not go beyond the last list member
                            start_id = 1
                            for id in range(len(self.atoms[atom2match]['bonds'])-1):
                                next_nbat_id  = id+1
                                next_nbat_at = self.atoms[atom2match]['bonds'][next_nbat_id]
                                if next_nbat_at not in already_visited:
                                    already_visited.append(next_nbat_at)
                                    matched_atom_types2(self.atoms[atom2match]['bonds'][next_nbat_id],t)
                                else:
                                    pass

                        else:
                            matched_atom_types2(nbat,t)
            print "\t\t\tlen alrvis %3d" % (len(already_visited))

        def createsybyllistonthefly(lig_atom):
            # look in matched_atom_types2 - line 656
            sybyllist = []
            for att in self.atoms[lig_atom]['lBondedAtoms']:
                sybyllist.append(att.sybylType)
            return sybyllist

        def gothroughallnbsofmatchlistatom(stored_nbs_of_atom2match,t,already_visited,hit_list):
            putative_next_a2m_list = []
            for ent_lig in stored_nbs_of_atom2match:
                matchlist = []
                if ent_lig not in already_visited:
                    already_visited.append(ent_lig)
                    # look for matching neighbours
                    for at in t.keys():
                        if t[at]['sybylType'] == self.atoms[ent_lig]['sybylType'] and ent_lig not in matchlist:
                            matchlist.append(at)
                    for matches in matchlist:
                        if len(Set(t[matches]['sybyl_neighbours']).difference(Set(createsybyllistonthefly(ent_lig)))) == 0:
                            if ent_lig not in hit_list:
                                hit_list.append(ent_lig)
                            for putative_next_atom2match in self.atoms[ent_lig]['bonds']:
                                if putative_next_atom2match not in putative_next_a2m_list:
                                    putative_next_a2m_list.append(putative_next_atom2match)
                        else:
                             print "sybyl neighbours don't match"
                else:
                    # what's here?
                    pass
            # delete the stored nbs
            stored_nbs_of_atom2match = []
            # next atom2match???
            return already_visited,stored_nbs_of_atom2match,putative_next_a2m_list,hit_list

        def matchatomtypeintemplateandgetliglist(atom2match,t,stored_nbs_of_atom2match=[],been_here_flag=False,\
                                                 already_visited=[],hit_list=[]):
            print atom2match,"hit_list",hit_list,been_here_flag
            putative_next_a2m_list = []
            # we don't want to miss the nbs of a matched atom (see [1])
            if atom2match in stored_nbs_of_atom2match:
                already_visited,stored_nbs_of_atom2match,putative_next_a2m_list,hit_list = \
                 gothroughallnbsofmatchlistatom(stored_nbs_of_atom2match,t,already_visited,hit_list)
            # does this really work? - to which position of the routine do we go now?
            if been_here_flag == True:
                print "it's true...", putative_next_a2m_list
                for next_at in putative_next_a2m_list:
                    print "TRUE (been_here_flag)", putative_next_a2m_list
                    matchatomtypeintemplateandgetliglist(next_at,t,been_here_flag=True)
            matchlist = []
            for at in t.keys():
                if t[at]['sybylType'] == self.atoms[atom2match]['sybylType'] and atom2match not in already_visited:
                    already_visited.append(atom2match)
                    print "we found a match for %4s " %(atom2match)
                    matchlist.append(at)
            # look for sybylnbs of all stored entries in matchlist
            for entries in matchlist:
                # TODO: Set deletes redundancies in 'sybyl_neighbours': AVOID THIS!
                if len(Set(t[entries]['sybyl_neighbours']).difference(Set(createsybyllistonthefly(atom2match)))) == 0:
                    hit_list.append(atom2match)
                    stored_nbs_of_atom2match = self.atoms[atom2match]['bonds']
                    print "nbs %s of hit %s" %(stored_nbs_of_atom2match,atom2match)
                    for nbs in stored_nbs_of_atom2match:
                        # call itself!
                        #
                        matchatomtypeintemplateandgetliglist(nbs,t,stored_nbs_of_atom2match)
                    # when passing stored_nbs_of_atom2match - always control atom2match with list entries! [1]
                    # (we have to do this at the beginning of our routine)
                    # ELSE case:

        def match2(t,l,start_atom,already_visited=[],type_matches=[]):
            matchatomtypeintemplateandgetliglist(start_atom,t)
#            matched_atom_types2(start_atom,t)

        for current_template in templates.keys():
            match2(templates[current_template],atoms,start_atom=atoms[4])  # start_atom should be 0
#            match(templates[current_template],atoms)
        
