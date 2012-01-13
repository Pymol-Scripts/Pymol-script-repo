#!/usr/bin/python
#
# pKa calculations with APBS
#
# Jens Erik Nielsen & Todd Dolinsky
#
# University College Dublin & Washington University St. Louis 2004-2006
#
__date__="16 August 2005"
__author__="Jens Erik Nielsen, Todd Dolinsky"

debug=None
import getopt
import sys, os
from pKa_base import *

if debug:
    from Tkinter import *

    class charge_mon(Frame):

        def __init__(self):
            Frame.__init__(self)
            self.master.title('Charge monitor')
            height=800
            width=1200
            self.cv=Canvas(self.master,bd=5,bg='white',
                           width=width,
                           height=height,
                           scrollregion=(0,0,width,height))
            self.cv.grid(row=0,column=0)
            self.calc=0
            self.seqstart=200
            self.text=''
            return

        def init_protein(self,pkaroutines):

            self.cv.create_text(0,self.calc,text='Setup',anchor='nw')
            x_count=self.seqstart
            self.res_pos={}
            for chain in pkaroutines.protein.chains:
                print chain.chainID
                for residue in chain.residues:
                    print residue.name, residue.resSeq
                    if residue.name=='ASP' or residue.name=='GLU':
                        fill='red'
                    elif residue.name=='LYS' or residue.name=='ARG':
                        fill='blue'
                    else:
                        fill='black'
                    self.cv.create_text(x_count,self.calc,text='%3d' %residue.resSeq,anchor='nw',fill=fill)
                    self.res_pos[residue.resSeq]=x_count
                    x_count=x_count+24
            self.calc=self.calc+15
            self.master.update()
            return

        def set_calc(self,text):
            self.text=text
            return

        def display_charges(self,charge_list):
            #
            # Print the calc
            #
            self.cv.create_text(0,self.calc,text=self.text,anchor='nw')
            charges={}
            for resnum,atomname,charge in charge_list:
                if not charges.has_key(resnum):
                    charges[resnum]=[]
                charges[resnum].append(charge)
            #
            # Sum all charges
            #
            for res in charges.keys():
                non_zero=None
                sum=0.0
                for crg in charges[res]:
                    sum=sum+crg
                    if crg!=0.0:
                        non_zero=1
                if non_zero:
                    charges[res]=sum
                else:
                    charges[res]=None
            #
            #
            #
            later=[]
            for resid in charges.keys():
                x_count=self.res_pos[resid]
                if charges[resid] is None:
                    fill='white'
                elif abs(charges[resid])<0.001:
                    fill='grey'
                elif abs(charges[resid]-1.0)<0.001:
                    fill='blue'
                elif abs(charges[resid]+1.0)<0.001:
                    fill='red'
                else:
                    fill='yellow'

                self.cv.create_rectangle(x_count,self.calc,x_count+24,self.calc+10,fill=fill)
                if fill=='yellow':
                    later.append([x_count,'%4.2f' %charges[resid]])
            #
            # Print all the wrong charges
            #
            for x_count,text in later:
                self.cv.create_text(x_count,self.calc,text=text,anchor='nw',fill='black')
            #
            # Update and increment row
            #
            self.master.update()
            self.calc=self.calc+15
            return
        
if debug:
    CM=charge_mon()
else:
    CM=None
#
# find the path to the script, and to pdb2pqr
#
scriptpath=os.path.split(sys.argv[0])[0]
if scriptpath=='.':
    scriptpath=os.getcwd()
pdb2pqr_path=os.path.split(scriptpath)[0]
sys.path.append(pdb2pqr_path)
#
# Set the phidir - where results of apbscalcs are stored
#
phidir=os.path.join(os.getcwd(),'phidir')
if not os.path.isdir(phidir):
    os.mkdir(phidir)
#
# --
#
import string
import math
import string
import getopt
import time
from src import pdb
from src import utilities
from src import structures
from src import routines
from src import protein
from src import server
from src.pdb import *
from src.utilities import *
from src.structures import *
from src.definitions import *
from src.forcefield import *  
from src.routines import *
from src.protein import *
from src.server import *
from StringIO import *
from src.hydrogens import *

import ligandclean.ligff
from apbs import *

##
# ---
# Function for comparing two atoms
#

def is_sameatom(atom1,atom2):
    #
    # Compare atom1 and atom2
    #
    properties=['name','resSeq','chainID']
    for attr in properties:
        a1_prop=getattr(atom1,attr,None)
        a2_prop=getattr(atom2,attr,None)
        if (attr!='chainID' and (not a1_prop or not a2_prop)) or a1_prop!=a2_prop:
            #print 'Failing because of %s "%s" "%s" \n' %(attr,str(a1_prop),str(a2_prop))
            return None
    return 1

#
# ----
#

TITRATIONFILE = os.path.join(scriptpath,"TITRATION.DAT")

class pKaRoutines:
    """
        Class for running all pKa related functions
    """
    def __init__(self, protein, routines, forcefield,apbs_setup):
        """
            Initialize the class using needed objects

            Parameters
                protein:    The PDB2PQR protein object
                routines:   The PDB2PQR routines object
                forcefield: The PDB2PQR forcefield object
        """
        self.protein = protein
        self.routines = routines
        self.forcefield = forcefield
        self.apbs_setup=apbs_setup
        self.pKagroups = self.readTitrationDefinition()
        
        self.pKas = []

        myHydrogenRoutines = hydrogenRoutines(routines)
        myHydrogenRoutines.readHydrogenDefinition()
        self.hydrogenRoutines = myHydrogenRoutines
        #
        # Not sure this is the best place for the interaction energies...
        #
        self.matrix={}
        return
    
    #
    # -----------------------------------------
    #
    

    def insert_new_titratable_group(self,ligand_titratable_groups):
        """Insert the new titratable groups in to self.pkagroups"""
        group_type=ligand_titratable_groups['type']
        if self.pKagroups.has_key(group_type):
            #
            # Now modify the group so that it will correspond to the group
            # we have in the ligand
            #
            ligand_name='LIG' #Note: we have to implement automatic determination of ligand name
            import copy
            new_group=copy.deepcopy(self.pKagroups[group_type])
            new_group.DefTitrations[0].modelpKa=ligand_titratable_groups['modelpka']
            new_group.name='LIG'
            new_group.resname='LIG'
            #print new_group.Residue
            
            self.pKagroups['LIG']=copy.deepcopy(new_group)
            atom_map=ligand_titratable_groups['matching_atoms']
            #
            # Insert definition into HYDROGEN arrays
            #
            for hdef in self.hydrogenRoutines.hydrodefs:
                if hdef.name==group_type:
                    newdef=copy.deepcopy(hdef)
                    #print newdef
                    newdef.name=ligand_name
                   
                    #
                    # Change the names in each of the conformatinos
                    #
                    # The name of the H is not changed!
                    #
                    for conformation in newdef.conformations:
                        #
                        # Change the name of the atom that the H is bound to
                        #
                        if atom_map.has_key(conformation.boundatom):
                            conformation.boundatom=atom_map[conformation.boundatom]
                        #
                        # Change the name of the hydrogen
                        #
                        oldhname=conformation.hname
                        conformation.hname='H'+conformation.boundatom
                        #
                        # And then for the individual atom names
                        #
                        for atom in conformation.atoms:
                            if atom_map.has_key(atom.name):
                                atom.name=atom_map[atom.name]
                            elif atom.name==oldhname:
                                atom.name=conformation.hname
                    self.hydrogenRoutines.hydrodefs.append(copy.deepcopy(newdef))
        #stop
        return

    #
    # -----------------------------------------
    #
        
    def runpKa(self):
        #
        #    Main driver for running pKa calculations
        #
        self.pKas = self.findTitratableGroups()
        
        self.calculateIntrinsicpKa()

        """ Calculate Pairwise Interactions """
        self.calculatePairwiseInteractions()
        
        """ Calculate Full pKa Value """
        self.calculatepKaValue()

    #
    # -----------------------------------
    #

    def calculatePairwiseInteractions(self):
        #
        # Calculate the pairwise interaction energies
        #
        for pKa in self.pKas:
            residue = pKa.residue
            pKaGroup = pKa.pKaGroup
            ambiguity = pKa.amb

            #
            # TODD: What does this one do?
            #
            #fixedstates = self.hydrogenRoutines.getstates(ambiguity)
            #
            # Loop over each titration
            #
            if not self.matrix.has_key(pKa):
                self.matrix[pKa]={}
            #
            for titration in pKaGroup.DefTitrations:
                if not self.matrix[pKa].has_key(titration):
                    self.matrix[pKa][titration]={}
                #
                # Get the atomnames
                #
                atomnames = self.getAtomsForPotential(pKa,titration)
                #
                # Get all states
                #
                possiblestates = titration.allstates
                for state in possiblestates:
                    #
                    # Here we switch the center group to a particular state
                    #
                    self.hydrogenRoutines.switchstate('pKa', ambiguity, state) 
                    self.zeroAllRadiiCharges()
                    self.setAllRadii()
                    self.setCharges(residue, atomnames)
                    #
                    # get_interaction_energies get the potential at all titratable groups due the charges
                    # this group
                    #
                    self.matrix[pKa][titration][state]=self.get_interaction_energies(pKa,titration,state)
        return

    #
    # ----
    #

    def get_interaction_energies(self,pKa_center,titration,state):
        """Get the potentials and charges at all titratable groups"""
        print 'In get_interaction_energies. center: group: %s, titration: %s, state: %s' %(pKa_center,titration,state)
        #
        # Set the calc type and center
        #
        self.apbs_setup.set_type('intene')
        #
        # Run APBS and get the interaction with all other states
        #
        if debug:
            CM.set_calc('IE %s %s' %(pKa_center.residue.resSeq,state))
        potentials=self.getAPBSPotentials(pKa_center,titration,state)
        #
        # construct this side
        #
        energies={}
        #
        # Loop over all groups
        #
        for pKa in self.pKas:
            residue = pKa.residue
            pKaGroup = pKa.pKaGroup
            ambiguity = pKa.amb
            #
            # Loop over each titration
            #
            if not energies.has_key(pKa):
                energies[pKa]={}
            #
            for titration in pKaGroup.DefTitrations:
                if not energies[pKa].has_key(titration):
                    energies[pKa][titration]={}
                #
                # Get all states
                #
                possiblestates = titration.allstates
                for state in possiblestates:
                    #
                    # Switch to the particular state
                    #
                    atomnames = self.getAtomsForPotential(pKa,titration)
                    self.hydrogenRoutines.switchstate('pKa', ambiguity, state) 
                    self.zeroAllRadiiCharges()
                    self.setAllRadii()
                    self.setCharges(residue, atomnames)
                    #
                    # Get atoms for potential
                    #
                    #print 'Atoms for measuring potential',atomnames
                    atomlist=[]
                    for atomname in atomnames:
                        atomlist.append(residue.getAtom(atomname))
                    energy=0.0
                    count=0
                    for atom in protein.getAtoms():
                        for atom2 in atomlist:
                            if is_sameatom(atom,atom2):
                                energy=energy+potentials[1][count]*atom.get("ffcharge")
                                #print 'Getting potential',residue.get('name'),atom.name,atom.get('ffcharge')
                        count=count+1
                    #
                    # We set all energies with self to zero
                    #
                    if pKa==pKa_center:
                        energies[pKa][titration][state]=0.0
                    else:
                        energies[pKa][titration][state]=energy
        return energies

    #
    # ----------------------------------
    #

    def calculatepKaValue(self):
        #
        #  Calculate the pKa Value
        #
        # We use a c++ class for the MC steps since it's a lot faster..
        #
        import pMC_mult
        #
        # Matrix needs to be linear for transport to c++
        # TODD: Do you know how to pass a full dictionary or list of lists to c++?
        #
        linear=[]
        intpkas=[]
        acidbase=[]
        state_counter=[]
        is_charged_state=[]
        print
        print
        print 'Residue\tIntrinsic pKa'
        for pKa in self.pKas:
            pKaGroup = pKa.pKaGroup
            for titration in pKaGroup.DefTitrations:
                #
                # Acid/Base
                #
                if pKaGroup.type=='acid':
                    acidbase.append(-1)
                else:
                    acidbase.append(1)
                #
                # 
                #
                print pKa.residue.get('name'),pKa.intrinsic_pKa
                possiblestates = titration.allstates
                #
                # Record the number of states for each titratable group
                #
                state_counter.append(len(possiblestates))

                pos_states=possiblestates
                pos_states.sort()
                for state in pos_states:
                    #
                    # Is this a charged state?
                    #
                    crg=self.is_charged(pKa,titration,state)
                    is_charged_state.append(crg)
                    print 'State: %s, is_charged  is %d' %(state,crg)
                    #if crg!=0:
                    intpkas.append(pKa.intrinsic_pKa[state])
                    #else:
                    #    intpkas=intpkas+[0.0]
                    #
                    # Now for the states that this state interacts with
                    #
                    for pKa2 in self.pKas:
                        pKaGroup2 = pKa2.pKaGroup
                        for titration2 in pKaGroup2.DefTitrations:
                            states2=titration2.allstates
                            states2.sort()
                            for state2 in states2:
                                linear.append(self.matrix[pKa][titration][state][pKa2][titration2][state2])
        #print
        print 'The charged_state array'
        print is_charged_state
        print 'These are the old intpkas'
        print intpkas

        #intpkas=[0.0,5.0,5.0,5.0,5.0]
        print 'Linearized matrix',linear
        print 'Length of linearised matrix',len(linear)
        mcsteps=50000
        phstart=2.0
        phend=12.0
        phstep=0.1
        #
        # Call our little C++ module
        #
        print
        print 'state_counter',state_counter
        print '-----------------------------------------------'
        print is_charged_state
        print
        #vector<double> intpKas, vector<double> lin_matrix,vector<double> ab,vector<int> state_counter,vector<int> charged_state
        FAST=pMC_mult.MC(intpkas,linear,acidbase,state_counter,is_charged_state)
        print
        print 'Setting MC steps'
        FAST.set_MCsteps(int(mcsteps))
        print
        print 'Starting to calculate pKa values'
        pKavals=FAST.calc_pKas(phstart,phend,phstep)
        count=0
        print
        print
        print 'Final pKa values'
        print
        for pKa in self.pKas:
            pKaGroup = pKa.pKaGroup
            for titration in pKaGroup.DefTitrations:
                titration.pKa=pKavals[count]
                count=count+1
                print '%s %s %5.2f' %(pKa.residue.get('name'),str(pKa.residue.get('resSeq')),titration.pKa)
        return

    #
    # ----------------------------------
    # 

    def is_charged(self,pKa,titration,state):
        #
        # Check if this state is a charged state
        #
        ambiguity = pKa.amb
        #states = self.hydrogenRoutines.getstates(ambiguity)
        self.hydrogenRoutines.switchstate('pKa', ambiguity, state)
        residue = pKa.residue
        pKaGroup = pKa.pKaGroup
        sum=0.0
        for atom in residue.atoms:
            atomname = atom.get("name")
            if atomname.find('FLIP')!=-1:
                continue
            charge, radius = self.forcefield.getParams(residue, atomname)
            sum=sum+charge
        print 'Total charge for %s is %.2f' %(residue.name,sum)
        if abs(sum)>0.05:
            return 1
        return 0

    #
    # ----------------------------------
    #
        
    def calculateIntrinsicpKa(self):
        #
        #  Calculate the intrinsic pKa
        #
        self.calculateDesolvation()
        self.calculateBackground()  
        #
        # Calculate the intrinsic pKas
        #
        # Print what we got
        #
        for pKa in self.pKas:
            print 'State\tModel pKa\tDesolvation\tBackground'
            for titration in pKa.pKaGroup.DefTitrations:
                for state in titration.allstates:
                    print '%10s\t%5.3f\t\t%5.3f\t\t%5.3f' %(state,titration.modelpKa,pKa.desolvation[state],pKa.background[state])
        print
        print
        #
        # We calculate an intrinsic pKa for every possible <startstate> -> <endstate> transition
        #
        import math
        ln10=math.log(10)
        for pKa in self.pKas:
            pKaGroup=pKa.pKaGroup
            #
            # We measure intrinsic pKa values against a single reference state
            #
            for titration in pKaGroup.DefTitrations:
                #
                # Find an uncharged reference state
                #
                ref_state=None
                for state in titration.allstates:
                    crg=self.is_charged(pKa,titration,state)
                    if crg==0.0:
                        ref_state=state
                if ref_state:
                    print 'State %s is the reference state' %ref_state
                else:
                    print 'No uncharged state for group',pKa
                    stop
                #
                all_states=titration.allstates
                all_states.sort()
                for state in all_states:
                    if self.is_charged(pKa,titration,state)==1:
                        dpKa_desolv=(pKa.desolvation[state]-pKa.desolvation[ref_state])/ln10
                        dpKa_backgr=(pKa.background[state]-pKa.background[ref_state])/ln10
                        intpKa=titration.modelpKa+dpKa_desolv+dpKa_backgr
                        print 'Energy difference for %s -> %s [reference state] is %5.2f pKa units' %(state,ref_state,intpKa)
                        pKa.intrinsic_pKa[state]=intpKa
                    else:
                        dpKa_desolv=(pKa.desolvation[state]-pKa.desolvation[ref_state])/ln10
                        dpKa_backgr=(pKa.background[state]-pKa.background[ref_state])/ln10
                        dpKa=dpKa_desolv+dpKa_backgr
                        print 'Energy difference for %s -> %s [reference state] is %5.2f kT' %(state,ref_state,dpKa)
                        pKa.intrinsic_pKa[state]=dpKa
            #
            # Loop over all titrations (the diffrence between titrations and transitions needs to become clearer)
            #
           ##  for titration in pKaGroup.DefTitrations:
##                 for transition_num in titration.transitions.keys():
##                     transition=titration.transitions[transition_num]
##                     #
##                     # Get the start and end states for this transition
##                     #
##                     start_s=transition['start']
##                     end_s=transition['end']
##                     print 'Transition: %3d, start: %3d, end: %3d' %(transition_num,start_s,end_s)
##                     dpKa_desolv=(pKa.desolvation[end_s]-pKa.desolvation[start_s])/ln10
##                     dpKa_backgr=(pKa.background[end_s]-pKa.background[start_s])/ln10
##                     intpKa=titration.modelpKa+dpKa_desolv+dpKa_backgr
##                     print 'dpKa_desolv: %.2f, dpKa_backgr: %.2f = Intrinsic pKa: %5.3f' %(dpKa_desolv,dpKa_backgr,intpKa)
##                     titration.intrinsic_pKa.append(intpKa)
        return

    #
    # --------------------------------
    #

    def calculateBackground(self):
        #
        #    Calculate background interaction energies
        #
        for pKa in self.pKas:
            residue = pKa.residue
            pKaGroup = pKa.pKaGroup
            ambiguity = pKa.amb
            
            print "Finding Background Interaction Energy for %s %s" %(residue.name, residue.resSeq)
            #
            # Loop over all titrations in this group
            #
            for titration in pKaGroup.DefTitrations:
                #
                # Get all the states for this titration
                #
                possiblestates = titration.startstates + titration.endstates
                print 'Possible states',possiblestates
                #fixedstates = self.hydrogenRoutines.getstates(ambiguity)
                #
                # Loop over all states and calculate the Backgound Interaction energy for each
                #
                for state in possiblestates:
                    #
                    # Get the atoms where we will measure the potential
                    #
                    firststate = possiblestates[0]
                    atomnames = self.getAtomsForPotential(pKa,titration)
                    print 'Atoms for measuring potential',atomnames
                    atomlist=[]
                    for atomname in atomnames:
                        atomlist.append(residue.getAtom(atomname))
                    #
                    # Switch the states of all other titratable groups to a neutral state
                    #
                    for other_pKa in self.pKas:
                        if pKa==other_pKa:
                            continue
                        for other_titration in other_pKa.pKaGroup.DefTitrations:
                            self.hydrogenRoutines.switchstate('pKa',other_pKa.amb,self.neutral_ref_state[other_pKa])
                    #
                    # Switch the state for the group in question
                    #
                    print "Calculating Backgound for state %s" % state
                    self.hydrogenRoutines.switchstate('pKa', ambiguity, state) 
                    self.zeroAllRadiiCharges()
                    self.setAllRadii()
                    #
                    # Set charges on all other residues
                    #
                    for chain in self.protein.getChains():
                        for otherresidue in chain.get("residues"):
                            if residue == otherresidue:
                                continue
                            otherlist = []
                            for atom in otherresidue.atoms:
                                if not atom:
                                    continue
                                if atom.get('name').find('FLIP')==-1:
                                    otherlist.append(atom.get("name"))
                            self.setCharges(otherresidue, otherlist)
                   
                    #
                    # Center the map on our residue
                    #
                    center=self.get_atoms_center(atomlist)
                    self.apbs_setup.setfineCenter(center)
                    self.apbs_setup.set_type('background')
                    #
                    # Run APBS
                    #
                    if debug:
                        CM.set_calc('background %s %s' %(pKa.residue.resSeq,state))
                    potentials=self.getAPBSPotentials(pKa,titration,state)
                    #
                    # Assign charges to our residue
                    #
                    self.setCharges(residue, atomnames)
                    #
                    # Return the potentials - same order as in atomnames
                    #
                    energy=self.get_elec_energy(potentials,atomlist)
                    #
                    # Done with Background calc for this state
                    #
                    pKa.background[state] = energy
        return


    #
    # --------------------------------
    #
                    
    def calculateDesolvation(self):
        #
        #   Calculate the Desolvation Energies
        #
        for pKa in self.pKas:
            residue = pKa.residue
            pKaGroup = pKa.pKaGroup
            ambiguity = pKa.amb
            print "Calculating Desolvation Energy for %s %s" %(residue.name, residue.resSeq)
            for titration in pKaGroup.DefTitrations:
                possiblestates = titration.allstates
                #
                # Get atoms for potential
                #
                atomnames = self.getAtomsForPotential(pKa,titration)
                print "atomnames", atomnames ### PC
                print 'Atoms for measuring potential',atomnames
                atomlist=[]
                for atomname in atomnames:
                    atomlist.append(residue.getAtom(atomname))
                #
                # Calculate the self energy for each state
                #        
                #fixedstates = self.hydrogenRoutines.getstates(ambiguity)
                for state in possiblestates:
                    print "Calculating desolvation energy for residue %s state %s in solvent" %(residue.name,state)
                    #
                    # Center the map on our set of atoms
                    #
                    center=self.get_atoms_center(atomlist)
                    self.apbs_setup.setfineCenter(center)
                    self.apbs_setup.set_type('desolv')
                    #
                    # Switch to the state
                    # Assign, radii, charges
                    self.hydrogenRoutines.switchstate('pKa', ambiguity, state) 
                    self.zeroAllRadiiCharges()
                    print 'Setting charges for residue',residue.resSeq,residue.name,atomnames
                    self.setCharges(residue, atomnames)
                    self.setRadii(residue, atomnames)
                    #
                    # Run APBS first time for the state in solvent
                    #
                    if debug:
                        CM.set_calc('Desolv solv %s %s' %(pKa.residue.resSeq,state))
                    solutionEnergy=self.get_elec_energy(self.getAPBSPotentials(pKa,titration,state),atomlist)
                    print "Calculating self energy for state %s in solvent" % state
                    #
                    # Now we set all radii (= in protein)
                    #
                    self.setAllRadii()
                    #
                    # Run APBS again, - this time for the state in the protein
                    #
                    print
                    print 'Calculating self energy for residue %s state %s in the protein' %(residue.name,state)
                    if debug:
                        CM.set_calc('Desolv prot %s %s' %(pKa.residue.resSeq,state))
                    proteinEnergy = self.get_elec_energy(self.getAPBSPotentials(pKa,titration,state),atomlist)
                    #
                    # Calculate the difference in self energy for this state
                    #
                    desolvation = proteinEnergy - solutionEnergy
                    print 'Desolvation for %s %d in state %s is %5.3f'  \
                          %(residue.name,residue.resSeq,state,desolvation)
                    print
                    print '======================================='
                    pKa.desolvation[state] = desolvation
        return

    #
    # ---------
    #

    def get_atoms_center(self,atomlist):
        #
        # Get the centre of a list of atoms
        #
        minmax={'x':[999.9,-999.9],'y':[999.9,-999.9],'z':[999.9,-999.9]}
        for atom in atomlist:
            if atom:
                for axis in ['x','y','z']:
                    coord=getattr(atom,axis)
                    if coord<minmax[axis][0]:
                        minmax[axis][0]=coord
                    if coord>minmax[axis][1]:
                        minmax[axis][1]=coord
        #
        # Calc the geometric center and extent
        #
        center={}
        extent={}
        for axis in minmax.keys():
            extent[axis]=minmax[axis][1]-minmax[axis][0]
            center[axis]=extent[axis]/2.0+minmax[axis][0]
        return [center['x'],center['y'],center['z']]
        

    #
    # -------------------------------
    #

    def get_elec_energy(self,potentials,atomlist):
        #
        # Given the electrostatic potential from getAPBSPotentials and a list
        # of atoms, this routine returns the energy in kT
        #
        # This function could be made a lot smarter!! (JN)
        # 
        energy=0.0
        count=0
        totphi=0.0
        totcrg=0.0
        found=0
        for atom in protein.getAtoms():
            if not atom:
                continue
            for atom_2 in atomlist:
                if not atom_2:
                    continue
                if is_sameatom(atom,atom_2):
                    totcrg=totcrg+abs(atom.get("ffcharge"))
                    totphi=totphi+abs(potentials[-1][count])
                    energy=energy+(potentials[-1][count])*atom.get("ffcharge")
                    print atom.name,atom.get('ffcharge'),(potentials[-1][count])
                    found=found+1
                    break
            #
            # This counter is outside the atom_2 loop!!
            #
            count=count+1
            if found==len(atomlist):
                break
        print '---------------'
        print 'Total energy',energy
        if abs(totphi)<0.01 or abs(totcrg)<0.01:
            print 'total abs phi',totphi
            print 'total abs crg',totcrg
            raise 'Something is rotten'
        print
        print
        return energy

    #
    # ----------------------------------
    #

    def getAPBSPotentials(self,group,titration,state):
        """
        #    Run APBS and get the potentials 
        #
        #    Returns
        #        list of potentials (list of floats)
        #"""
        #
        # Do we have results for this calculation?
        #
        result_file=os.path.join(phidir,'%s_%s_%s.potentials' %(self.apbs_setup.type,group.uniqueid,state))
        import cPickle
        loaded=None
        if os.path.isfile(result_file):
            #
            # Yes!
            #
            fd=open(result_file,'rb')
            try:
                potentials=cPickle.load(fd)
                fd.close()
                loaded=1
            except:
                fd.close()
        loaded=None
        #
        # Run calc again if needed
        #
        if not loaded:
            print 'Running APBS calculation'
            apbs_inputfile=self.apbs_setup.printInput()
            potentials = runAPBS(self.protein, apbs_inputfile, CM)
            fd=open(result_file,'wb')
            cPickle.dump(potentials,fd)
            fd.close()
        return potentials



    #
    # ----------------------
    #
                
    def setRadii(self, residue, atomlist):
        """
            Set the radii for specific atoms in a residue

            Parameters
                residue:  The residue to set (residue)
                atomlist: A list of atomnames (list)
        """
        for atom in residue.get("atoms"):
            atomname = atom.get("name")
            if atomname not in atomlist: continue
            charge, radius = self.forcefield.getParams(residue, atomname)
            if radius != None:
                atom.set("radius", radius)
            else:
                text = "Could not find radius for atom %s" % atomname
                text += " in residue %s %i" % (residue.name, residue.resSeq)
                text += " while attempting to set radius!"
                raise ValueError, text

    #
    # ------------------------------------
    #
            
    def setCharges(self, residue, atomlist):
        """
            Set the charges for specific atoms in a residue
            
            Parameters
                residue:  The residue to set (residue)
                atomlist: A list of atomnames (list)
        """
        for atom in residue.get("atoms"):
            atomname = atom.get("name")
            if atomname not in atomlist:
                continue
            charge, radius = self.forcefield.getParams(residue, atomname)
            if charge != None:
                atom.set("ffcharge", charge)
            else:
                text = "Could not find charge for atom %s" % atomname
                text += " in residue %s %i" % (residue.name, residue.resSeq)
                text += " while attempting to set charge!"
                raise ValueError, text
        return
    #
    # ----------------------------
    #

    def setAllRadii(self):
        """
            Set all radii for the entire protein
        """
        for chain in self.protein.getChains():
            for residue in chain.get("residues"):
                for atom in residue.get("atoms"):
                    atomname = atom.get("name")
                    if atomname.find('FLIP')!=-1:
                        continue
                    if atomname == "HD1": ###PC
                        charge = 0.44
                        radiues = 1.05
                    else:
                     charge, radius = self.forcefield.getParams(residue, atomname)
                    ###PC
                    if radius != None:
                        atom.set("radius", radius)
                    else:
                        if residue.type != 2:
                            text = "Could not find radius for atom %s " % atomname
                            text +="in residue %s %i" % (residue.name, residue.resSeq)
                            text += " while attempting to set all radii!"
                            raise ValueError, text
    #
    # -------------------------------
    #
                        
    def zeroAllRadiiCharges(self):
        """
            Set all charges and radii for the protein to zero
        """
        for chain in self.protein.getChains():
            for residue in chain.get("residues"):
                for atom in residue.get("atoms"):
                    atom.set("ffcharge",0.0)
                    atom.set("radius",0.0)

    #
    # --------------------------------
    #
                
    def getAtomsForPotential(self, pKa,titration, get_neutral_state=None):
        #
        #    Find the atoms that are needed for measuring the potential,
        #    only selecting atoms where the charge changes.
        #    Parameters
        #        pKa:  The pKa object (pKa)
        #    Returns:
        #        atomnames:  A list of atomnames to measure (list)
        neutral_state=None
        atomnames = []
        initialmap = {}
        residue = pKa.residue
        pKaGroup = pKa.pKaGroup
        ambiguity = pKa.amb
        #states = self.hydrogenRoutines.getstates(ambiguity)
        #
        # Change to the start state
        #
        start_state=titration.startstates[0]
        self.hydrogenRoutines.switchstate('pKa', ambiguity, start_state)
        sum=0.0
        for atom in residue.atoms:
            atomname = atom.get("name")
            if atomname.find('FLIP')!=-1:
                continue
            charge, radius = self.forcefield.getParams(residue, atomname) 
            initialmap[atomname] = charge
            sum=sum+charge
        if abs(sum)<0.001:
            neutral_state=start_state
        #
        # Check if charges change in all other states
        #
        for state in titration.endstates+titration.startstates[1:]:
            self.hydrogenRoutines.switchstate('pKa', ambiguity, state)
            sum=0.0
            for atom in residue.atoms:
                atomname = atom.get("name")
                if atomname.find('FLIP')!=-1:
                    continue
                charge, radius = self.forcefield.getParams(residue, atomname)
                sum=sum+charge
                try:
                    initcharge = initialmap[atomname]
                    if charge != initcharge:
                        if not atomname in atomnames:
                            atomnames.append(atomname)
                except KeyError:
                    if not atomname in atomnames:
                        atomnames.append(atomname)
            #
            # Is this the first neutral state?
            #
            if abs(sum)<0.0001:
                if not neutral_state:
                    neutral_state=state
        #
        # Did we just want a neutral state identification?
        #
        if get_neutral_state:
            if not neutral_state:
                raise "no neutral state for",residue.resSeq
            return neutral_state
        #
        # Make sure that the charges add up to integers by adding extra atoms
        #
        sum=0.01
        while sum>0.001:
            sum=0.0
            added=None
            #
            # Loop over all states to find atoms to add
            #
            for state in titration.endstates+titration.startstates:
                self.hydrogenRoutines.switchstate('pKa', ambiguity, state)
                #
                # Sum this state
                #
                this_sum=0.0
                for atom in residue.atoms:
                    atomname = atom.get("name")
                    if atomname.find('FLIP')!=-1:
                        continue
                    if not atomname in atomnames:
                        continue
                    charge, radius = self.forcefield.getParams(residue, atomname)
                    this_sum=this_sum+charge
                #
                # Is this an integer charge?
                #
                diff=abs(1000.0*this_sum)-abs(1000.0*int(this_sum))
                sum=sum+diff
                if diff>0.001:
                    #
                    # Find all atoms one bond away
                    #
                    for atom in residue.atoms:
                        atomname=atom.get('name')
                        if atomname.find('FLIP')!=-1:
                            continue
                        if not atomname in atomnames:
                            continue
                        #
                        # Add all atoms that are not already added
                        #
                        for bound_atom in atom.intrabonds:
                            if not bound_atom in atomnames:
                                atomnames.append(bound_atom)
                                added=1
                #
                # Next state
                #
                pass
            #
            # Did we add anything?
            #
            if added is None and sum>0.001:
                print sum
                print atomnames
                raise 'Could not find integer charge state'
        return atomnames

    #
    # -------------------------------
    #


    def findTitratableGroups(self):
        """
            Find all titratable groups in the protein based on the definition

            Returns
                pKalist:  A list of pKa objects (list)
        """
        pKalist = []
        
        print "Finding Titratable residues:"
        for chain in self.protein.getChains():
            for residue in chain.get("residues"):
                resname = residue.get("name")
                for group in self.pKagroups:
                    if resname == group:
                        amb = None
                        for hydrodef in self.hydrogenRoutines.hydrodefs:
                            hydname = hydrodef.name
                            if hydname == group:
                                amb = hydrogenAmbiguity(residue, hydrodef,self.routines)
                        if amb == None:
                            text = "Could not find hydrogen ambiguity "
                            text += "for titratable group %s!" % group
                            raise ValueError, text
                             
                        thispKa = pKa(residue, self.pKagroups[group], amb)
                        pKalist.append(thispKa)
                        print "\t%s %s" % (resname, residue.resSeq)
 
        #
        # Print the residues that we have selected
        #
        print
        print
        print 'Titratable residues'
        for pKa_v in pKalist:
            print pKa_v.residue.name,pKa_v.residue.resSeq
        print
        print
        #
        # Find a neutral state for each group
        #
        self.neutral_ref_state={}
        for this_pka in pKalist:
            residue = this_pka.residue
            pKaGroup = this_pka.pKaGroup
            ambiguity = this_pka.amb
            for titration in pKaGroup.DefTitrations:
                neutral_state = self.getAtomsForPotential(this_pka,titration,get_neutral_state=1)
                self.neutral_ref_state[this_pka]=neutral_state
        return pKalist

    #
    # ----------------------------------
    #
        
    def readTitrationDefinition(self):
        """
            Read the Titration Definition

            Returns:
               mygroups: A dictionary of pKaGroups
        """
        mygroups = {}
        filename = TITRATIONFILE
        if not os.path.isfile(TITRATIONFILE):
            raise ValueError, "Could not find TITRATION.DAT!"
        file = open(filename)
        
        while 1:
            line=file.readline()
            if line.startswith("//"): pass
            elif line == '': break
            elif line[0]=='*':
                name = ""
                resname = ""
                type = ""
                titrations = []

                name = string.strip(line[1:])
                line = file.readline()
                if line[:8] != 'Residue:':
                    text = "Wrong line found when looking for 'Residue'"
                    raise ValueError, "%s: %s" % (text, line)
                
                resname = string.strip(string.split(line)[1])
            
                line = file.readline()
                if line[:10] != 'Grouptype:':
                    text = "Wrong line found when looking for 'Grouptype'"
                    raise ValueError, "%s: %s" % (text, line)
                
                type = string.lower(string.strip(string.split(line)[1]))
                if type != 'acid' and type != 'base':
                    raise ValueError, 'Group type must be acid or base!'

                line = file.readline()
                while 1:  
                    """ Find next transition """
                    #
                    # Skip comments
                    #
                    while line[:2]=='//':
                        line=file.readline()
                    
                    startstates = []
                    endstates = []
                    modelpKa = None

                    if line[:11] != 'Transition:':
                        text = "Wrong line found when looking for 'Transition:'"
                        raise ValueError, "%s: %s" % (text, line)
                    
                    split=string.split(line[11:],'->')
                    for number in string.split(split[0], ','):
                        startstates.append(string.strip(number))                        
                    for number in string.split(split[1], ','):
                        endstates.append(string.strip(number))

                    line = file.readline()
                    #
                    # Skip comments
                    #
                    while line[:2]=='//':
                        line=file.readline()
                    #
                    # Must be the model pKa line
                    #
                    if line[:10]!='Model_pKa:':
                        text = "Wrong line found when looking for 'Model_pKa'"
                        raise ValueError, "%s: %s" % (text, line)
                       
                    modelpKa = float(string.split(line)[1])

                    thisTitration = DefTitration(startstates, endstates,modelpKa)
                    titrations.append(thisTitration)
                    
                    line = file.readline()
                    if string.strip(line) == 'END': break

                thisGroup = pKaGroup(name, resname, type, titrations)
                mygroups[name] = thisGroup

                line = file.readline()
                if string.strip(line) == 'END OF FILE': break
        
        return mygroups


#
# -------------
#

def usage(x):
    #
    # Print usage guidelines
    #
    print 'Usage: pka.py --ff <forcefield> --lig <ligand in MOL2> <pdbfile>'
    print 'Force field can be amber, charmm and parse'
    print
    return

#
# --------------------------------------------------
#

def startpKa():
    """
        Function for starting pKa script from the command line.

        Returns
            protein:    The protein object as generated by PDB2PQR
            routines:   The routines object as generated by PDB2PQR
            forcefield: The forcefield object as generated by PDB2PQR
    """
    print
    print 'PDB2PQR pKa calculations'
    print
    shortOptlist = "h,v"
    longOptlist = ["help","verbose","ff=",'lig=']

    try:
        opts, args = getopt.getopt(sys.argv[1:], shortOptlist, longOptlist)
    except getopt.GetoptError, details:
        sys.stderr.write("GetoptError:  %s\n" % details)
        usage(2)
        sys.exit(0)
    #
    #
    #
    if len(args) < 1 or len(args) > 3:
        sys.stderr.write("Incorrect number (%d) of arguments!\n" % len(args))
        usage(2)
        sys.exit(0)

    verbose = 0
    ligfilename=None
    ff = None
    for o,a in opts:
        if o in ("-v","--verbose"):
            verbose = 1
        elif o in ("-h","--help"):
            usage(2)
            sys.exit()
        elif o == "--ff":
            if a in ["amber","AMBER","charmm","CHARMM","parse","PARSE"]:
                ff = string.lower(a)
            else:
                raise ValueError, "Invalid forcefield %s!" % a
        elif o == "--lig":
            ligfilename=a

    #
    # No forcefield?
    #
    if ff == None:
        raise ValueError, "Forcefield not specified!"
    #
    # Find the PDB file
    #
    path = args[0]
    #
    # Call the pre_init function
    #
    return pre_init(pdbfilename=path,ff=ff,ligand=ligfilename)

#
# ----
#
    

def pre_init(pdbfilename=None,ff=None,ligand=None,verbose=1):
    #
    # Get the PDBfile
    #
    pdbfile = getFile(pdbfilename)
    pdblist, errlist = readPDB(pdbfile)
    
    if len(pdblist) == 0 and len(errlist) == 0:
        print "Unable to find file %s!\n" % path
        os.remove(path)
        sys.exit(2)

    if len(errlist) != 0 and verbose:
        print "Warning: %s is a non-standard PDB file.\n" %pdbfilename
        print errlist

    if verbose:
        print "Beginning PDB2PQR...\n"
    #
    # Read the definition file
    #
    myDefinition = Definition()
    ligand_titratable_groups=None
    #
    #
    # Choose whether to include the ligand or not
    #
    # Add the ligand to the pdb2pqr arrays
    #
    Lig=None
    MOL2FLAG = False 
    if not ligand:
        myProtein = Protein(pdblist)
    else:
        #
        # Read the ligand into Paul's code
        #
        Lig = ligandclean.ligff.ligand_charge_handler()
        Lig.read(ligand)
        #
        # Create the ligand definition from the mol2 data
        #
        import NEWligand_topology
        MOL2FLAG = True
        #
        X=NEWligand_topology.get_ligand_topology(Lig.lAtoms,MOL2FLAG)
        #
        # Add it to the 'official' definition
        #
        ligresidue=myDefinition.parseDefinition(X.lines, 'LIG', 2)
        myDefinition.AAdef.addResidue(ligresidue)
        #
        # Look for titratable groups in the ligand
        #
        print '==============================\n================================\n======================='
        ligand_titratable_groups=X.find_titratable_groups()
        print '==============================\n================================\n======================='
        print "ligand_titratable_groups", ligand_titratable_groups
        #
        # ------------------------------------------------------
        # Creation of ligand definition and identification of ligand titgrps done
        # Start loading everything into PDB2PQR
        #
        # Append the ligand data to the end of the PDB data
        #
        newpdblist=[]
        # First the protein
        for line in pdblist:
            if isinstance(line, END) or isinstance(line,MASTER):
                continue
            newpdblist.append(line)
        # Now the ligand
        for e in Lig.lAtoms:
            newpdblist.append(e)
        #
        # Add a TER and an END record for good measure
        #
        newpdblist.append(TER)
        newpdblist.append(END)
        #
        # Let PDB2PQR parse the entire file
        #
        myProtein = Protein(newpdblist)

        #
        # Post-Processing for adding sybylTypes to lig-atoms in myProtein
        # Jens: that's the quick and easy solution
        for rrres in  myProtein.chainmap['L'].residues:
            for aaat in rrres.atoms:
                for ligatoms in Lig.lAtoms:
                    if ligatoms.name == aaat.name:
                        aaat.sybylType = ligatoms.sybylType
                        #
                        # setting the formal charges
                        if ligatoms.sybylType == "O.co2":
                            aaat.formalcharge = -0.5
                        else: aaat.formalcharge = 0.0
                        xxxlll = []
                        for xxx in ligatoms.lBondedAtoms:
                            xxxlll.append(xxx.name)
                        aaat.intrabonds = xxxlll
                        #
                        # charge initialisation must happen somewhere else
                        # but i don't know where...
                        aaat.charge = 0.0
        #
    #
    # =======================================================================
    #
    # We have identified the structural elements, now contiue with the setup
    #
    # Print something for some reason?
    #
    if verbose:
        print "Created protein object -"
        print "\tNumber of residues in protein: %s" % myProtein.numResidues()
        print "\tNumber of atoms in protein   : %s" % myProtein.numAtoms()
    #
    # Set up all other routines
    #
    myRoutines = Routines(myProtein, verbose, myDefinition)
    myRoutines.updateResidueTypes()
    myRoutines.updateSSbridges()
    myRoutines.updateIntraBonds()
    myRoutines.updateExtraBonds()
    myRoutines.correctNames()
    myRoutines.findMissingHeavy()
    myRoutines.addHydrogens()
    #myRoutines.randomizeWaters()
    myProtein.reSerialize()
    #
    # Inject the information on hydrogen conformations in the HYDROGENS.DAT arrays
    # We get this information from ligand_titratable_groups
    #
    from src.hydrogens import hydrogenRoutines
    myRoutines.updateIntraBonds()
    myRoutines.calculateChiangles()
    myhydRoutines = hydrogenRoutines(myRoutines)
    myhydRoutines.readHydrogenDefinition()
    #
    # Here we should inject the info!!
    #
    optflag=1 # 1 for optimising everything, otherwise optimize only waters..
    myhydRoutines.optimizeHydrogens(optflag)
    #
    # Choose the correct forcefield
    #
    if not Lig:
        myForcefield = Forcefield(ff)
    else:
        # Here I am now passing the instance of Lig that's already been create - I have not
        # yet changed it in ligforcefield.... It will give an error for now
        #
        myForcefield = ligandclean.ligff.ligforcefield(ff,Lig)
        print "That's Lig: ", Lig
        print "That's X:   " , X
    
    if verbose:
        print "Created protein object (after processing myRoutines) -"
        print "\tNumber of residues in protein: %s" % myProtein.numResidues()
        print "\tNumber of atoms in protein   : %s" % myProtein.numAtoms()
    #
    # Create the APBS input file
    #
    import src.psize
    size=src.psize.Psize()

    method=""
    async=0
    split=0
    import inputgen_pKa
    igen = inputgen_pKa.inputGen(pdbfilename)
    #
    # For convenience
    #
    igen.pdie=8
    igen.sdie=80
    #
    # test for Paul
    #
    #for chain in myProtein.getChains():
    #    for residue in chain.get("residues"):
    #        resname = residue.get("name")
    #        for atom in residue.get("atoms"):
    #            print atom.name,atom.intrabonds
    #
    # Return all we need
    #
    return myProtein, myRoutines, myForcefield,igen, ligand_titratable_groups

#
# -----------------------------------------------
#

def test_interface():
    """Test the interface with pKaTool"""
    import pKaTool.pKa_calc
    X=pKaTool.pKa_calc.Monte_Carlo_Mult_CPP()
    
    X.intrinsic_pKa={':0001:ASP':[0.0,4.0,5.0]}
    X.charged_state={':0001:ASP':[0,1,1]}
    X.acid_base={':0001:ASP':-1}
    X.intene_mult={':0001:ASP':{':0001:ASP':[[0,0,0],[0,0,0],[0,0,0]]}}
    X._calc_pKas(0.0,10.0,0.5)
    return
 
            
if __name__ == "__main__":
    #test_interface()
    #stop
    state=1
    if state==0:
        import pMC_mult
        groups=1
        states=3
        ip=[6.0,3.0,4.0,6.0]
        intpkas=[5.0,2.0,2.0]
        state_counter=[]
        is_charged_state=[1,0,0]
        linear=[]
        acidbase=[]
        names=[]
        for group in range(groups):
            #
            # Names
            #
            names.append('Group %d' %group)
            #
            # Add intrinsic pKas and is_charged_state
            #
            #intpkas.apend(ip[group])
            #is_charged_state.append(1)
            #for state in range(states-1):
                #intpkas.append(0.0)
            #    is_charged_state.append(0)
            #
            # Add state_counter
            #
            state_counter.append(states)
            #
            # Acid_base
            #
            acidbase.append(-1)
            #
            # Matrix
            #
            for group2 in range(groups):
                for state1 in range(states):
                    for state2 in range(states):
                        linear.append(0.0)
            
        #names=['group1','group2']
        #intpkas=[6.0,0.0,0.0,3.0,0.0,0.0]
        #linear=[0.0,0.0,0.0,0.0,
        #0.0,0.0,0.0,0.0,
        #0.0,0.0,0.0,0.0,
        #0.0,0.0,0.0,0.0]
        #state_counter=[3,3]
        #acidbase=[-1,-1]
        #is_charged_state=[1,0,0,1,0,0]
        print 'Linearized matrix',linear
        print 'Length of linmatrix',len(linear)
        mcsteps=50000
        phstart=2.0
        phend=12.0
        phstep=0.1
        #
        # Call our little C++ module
        #
        print
        print 'state_counter',state_counter
        print '-----------------------------------------------'
        print is_charged_state
        print
        #vector<double> intpKas, vector<double> lin_matrix,vector<double> ab,vector<int> state_counter,vector<int> charged_state
        FAST=pMC_mult.MC(intpkas,linear,acidbase,state_counter,is_charged_state)
        print
        print 'Setting MC steps'
        FAST.set_MCsteps(int(mcsteps))
        print
        print 'Starting to calculate pKa values'
        pKavals=FAST.calc_pKas(phstart,phend,phstep)
        count=0
        print
        print 'Final pKa values'
        for name in names:
            print '%s pKa: %5.2f ' %(name,pKavals[count])
            count=count+1
        stop
    elif state==1:
        #
        # Do a real pKa calculation
        #
        #try:
        protein, routines, forcefield,apbs_setup, ligand_titratable_groups = startpKa()
        print "chains of current protein ", protein.chainmap#     print "currently used forefield  ", forcefield
        mypkaRoutines = pKaRoutines(protein, routines, forcefield,apbs_setup)
        #
        # Debugging
        #
        if debug:
            CM.init_protein(mypkaRoutines)
        #
        # Deal with ligands
        #
        if ligand_titratable_groups:
            print "lig (before mypKaRoutines) ", ligand_titratable_groups['titratableatoms']
            mypkaRoutines.insert_new_titratable_group(ligand_titratable_groups)
        mypkaRoutines.runpKa()
        #finally:
        #    
        #    CM.mainloop()
    elif state==2:
        #
        # Just assign charges
        #
        protein, routines, forcefield,apbs_setup, ligand_titratable_groups = startpKa()
        for chain in protein.getChains():
            for residue in chain.get("residues"):
                for atom in residue.get("atoms"):
                    atomname = atom.get("name")
                    charge, radius = forcefield.getParams(residue, atomname)
                    print '%2s %4s %3d %4s %5.2f %5.2f' %(chain.chainID,residue.name,residue.resSeq,atomname,charge,radius)
