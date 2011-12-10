"""
    Hydrogen optimization routines for PDB2PQR

    This module contains the hydrogen optimization routines and classes for
    PDB2PQR.  It is (optionally) used to check protonation states and
    improve hydrogen networks within a protein.

    Based on C code from Jens Erik Nielsen
    UCSD/HHMI
    
    Ported to Python by Todd Dolinsky (todd@ccb.wustl.edu)
    Washington University in St. Louis
"""

__date__ = "8 March 2004"
__author__ = "Jens Erik Nielsen, Todd Dolinsky"

HYDROGENFILE = "HYDROGENS.DAT"
HACCEPTOR=[0,1,1,0,1,1,0,1,0,1,0,0,1,0]
HDONOR   =[0,0,1,1,1,0,1,0,1,0,0,0,1,1]
HYDROGEN_DIST = 6.0
WATER_DIST=4.0

import os
import string

from definitions import *
from utilities import *
from math import *
from quatfit import *
from random import *
from time import *

class hydrogenAmbiguity:
    """
        A class containing information about the ambiguity
    """
    def __init__(self, residue, hdef):
        """
            Initialize the class

            Parameters
                residue:  The residue in question (residue)
                hdef:     The hydrogen definition matching the residue
        """
        self.residue = residue
        self.hdef = hdef
        self.nearatoms = []

    def setNearatoms(self, allatoms):
        """
            Set the nearby atoms to this residue.  The only donors/acceptors
            that will be changing positions are the flips, with maximum change
            of 2.5 A.

            Parameters
                allatoms:  A list of all donors/acceptors (list)
        """
        nearatoms = []
        resname = self.residue.get("name")
        for atom in self.residue.get("atoms"):
            if (atom.get("hacceptor") or atom.get("hdonor")):
                for nearatom in allatoms:
                    if nearatom.get("residue") == self.residue: continue
                    nearres = atom.get("residue")
                    nearname = nearres.get("name")
                    dist = distance(atom.getCoords(), nearatom.getCoords())
                    compdist = HYDROGEN_DIST
                    if resname in ["HIS","ASN","GLN"]: compdist += 2.5
                    if nearname in ["HIS","ASN","GLN"]: compdist += 2.5
                    if dist < compdist and nearatom not in nearatoms:
                        nearatoms.append(nearatom)     
        self.nearatoms = nearatoms
        
class hydrogenRoutines:
    """
        The main class of routines in the hydrogen optimization process.
    """
    
    def __init__(self, routines):
        """
            Initialize the routines and run the hydrogen optimization

            Parameters
                routines:  The parent routines object (Routines)
        """
        self.hdebug = 0
        self.routines = routines
        self.protein = routines.protein
        self.hydrodefs = []
        self.groups = []
        self.watermap = {}
        self.accmap = {}
    
    def debug(self, text):
        """
            Print text to stdout for debugging purposes.

            Parameters
                text:  The text to output (string)
        """
        if self.hdebug:
            print text

    def getstates(self, amb):
        """
            Get all possible states for a conformation/protonation
            ambiguity and store them in a list. Each.
            hydrogen type must be explicitly defined.

            Parameters
                amb   : The ambiguity to get the states of (tuple)
            Returns
                states: A list of states, where each state
                        is a list of conformations of the atom. (list)
        """
        states = []
        residue = getattr(amb,"residue")
        hdef = getattr(amb,"hdef")
        type = hdef.type

        confs = hdef.conformations

        # Now make the states

        if type == 1: # CTR/ASP/GLU
            states.append([()])
            states.append([confs[0]]) # H on O1 cis
            states.append([confs[1]]) # H on O1 trans
            states.append([confs[2]]) # H on O2 cis
            states.append([confs[3]]) # H on O2 trans
        elif type == 3: # NTR
            states.append([()]) # No Hs
            states.append([confs[0]]) # H3 only
            states.append([confs[1]]) # H2 only
            states.append([confs[0], confs[1]]) # H3 and H2
        elif type == 4: # HIS 
            states.append([confs[0]]) # HSD
            states.append([confs[1]]) # HSE
            states.append([()]) # Negative HIS
            states.append([confs[0], confs[1]]) #HSP
        elif type == 10: #PNTR
            states.append([()])
            states.append([confs[0]])
            states.append([confs[1]])
            states.append([confs[0], confs[1]])
        elif type == 13: #GLH
            states.append([confs[0]]) # H on O1 cis
            states.append([confs[1]]) # H on O1 trans
            states.append([confs[2]]) # H on O2 cis
            states.append([confs[3]]) # H on O2 trans
        return states          
                
    def switchstate(self, states, amb, id):
        """
            Switch a residue to a new state by first removing all
            hydrogens.

            Parameters
                states: The list of states (list)
                amb   : The amibiguity to switch (tuple)
                id    : The state id to switch to (int)
        """
        if id > len(states):
            raise ValueError, "Invalid State ID!"
        
        # First Remove all Hs
        residue = getattr(amb,"residue")
        hdef = getattr(amb,"hdef")
        type = hdef.type
        for conf in hdef.conformations:
            hname = conf.hname
            boundname = conf.boundatom
            if residue.getAtom(hname) != None:
                residue.removeAtom(hname)
            residue.getAtom(boundname).hacceptor = 1
            residue.getAtom(boundname).hdonor = 0

        # Update the IntraBonds
        name = residue.get("name")
        defresidue = self.routines.aadef.getResidue(name)
        residue.updateIntraBonds(defresidue)
                              
        # Now build appropriate atoms
        state = states[id]
        for conf in state:
            refcoords = []
            defcoords = []
            defatomcoords = []
            if conf == (): continue # Nothing to add
            hname = conf.hname
            for atom in conf.atoms:
                atomname = atom.get("name")
                resatom = residue.getAtom(atomname)
                if atomname == hname:
                    defatomcoords = atom.getCoords()
                elif resatom != None:
                    refcoords.append(resatom.getCoords())
                    defcoords.append(atom.getCoords())
                else:
                    raise ValueError, "Could not find necessary atom!"
        
            newcoords = findCoordinates(3, refcoords, defcoords, defatomcoords)
            boundname = conf.boundatom
            residue.createAtom(hname, newcoords, "ATOM")
            residue.addDebumpAtom(residue.getAtom(hname))
            residue.getAtom(boundname).addIntraBond(hname)    
            residue.getAtom(boundname).hacceptor = 0
            residue.getAtom(boundname).hdonor = 1

    def optimizeSingle(self,amb):
        """
            Use brute force optimization for a single ambiguity - try all
            energy configurations and pick the best.

            Parameters
                amb:  The ambiguity object (tuple)
        """

        residue = getattr(amb,"residue")
        hdef = getattr(amb,"hdef")
        type = hdef.type
        self.debug("Brute Force Optimization for residue %s %i - type %i" %\
                   (residue.get("name"), residue.get("resSeq"), type))
        
        best = 0
        energy = None
        bestenergy = 1000.0
        
        # Brute force for fixed states
        
        if type in [1,4,3,10,13]:
            if type == 4:
                raise ValueError, "We shouldn't have a brute force HIS without the FLIP!"
            states = self.getstates(amb)
            for i in range(len(states)): 
                self.switchstate(states, amb, i)
                energy = self.getHbondEnergy(amb)
                if energy < bestenergy:
                    bestenergy = energy
                    best = i        
            self.switchstate(states, amb, best)
            
        # Brute force for chi angle
             
        elif type in [2]:
            name = residue.get("name")
            defresidue = self.routines.aadef.getResidue(name)
            chinum = hdef.chiangle - 1
            for angle in range(-180, 180, 5):
                self.routines.setChiangle(residue, chinum, angle, defresidue)
                energy = self.getHbondEnergy(amb)    
                if energy < bestenergy:
                    bestenergy = energy
                    best = angle
            self.routines.setChiangle(residue, chinum, best, defresidue)

        # Brute force for flips

        elif type in [11]:
            bestenergy = self.getHbondEnergy(amb)
            name = residue.get("name")
            defresidue = self.routines.aadef.getResidue(name)
            chinum = hdef.chiangle - 1
            oldangle = residue.get("chiangles")[chinum]
            angle = 180.0 + oldangle
            self.routines.setChiangle(residue, chinum, angle, defresidue)
            energy = self.getHbondEnergy(amb)
            if energy >= bestenergy: # switch back!
                self.routines.setChiangle(residue, chinum, oldangle, defresidue)
        else:
            raise ValueError, "Invalid Hydrogen type %i in %s %i!" % \
                  (type, residue.get("name"), residue.get("resSeq"))

    def optimizeHydrogens(self):
        """
            Optimize hydrogens according to HYDROGENS.DAT.  This
            function serves as the main driver for the optimizing
            script.
        """
        starttime = time()
        allatoms = self.findAmbiguities(0)
        self.printAmbiguities()
        networks = self.findNetworks(HYDROGEN_DIST)
        for cluster in networks:
            #clusteratoms, compatoms = self.initHbondEnergy(cluster, allatoms)
            if len(cluster) == 1:
                amb = self.groups[cluster[0]]
                residue = getattr(amb, "residue")
                amb.setNearatoms(allatoms)
                self.optimizeSingle(amb)

            else:     
                # Use Monte Carlo algorithm to optimize
                
                steps = 0
                if len(cluster) == 2: steps = 10
                elif len(cluster) == 3: steps = 15
                elif len(cluster) >= 4 and len(cluster) < 10: steps = pow(2,len(cluster))
                if steps > 200 or len(cluster) >= 10: steps = 200
            
                # Initialize

                statemap = {}
                curmap = {}
                bestmap = {}
                energy = 0.0
                
                for id in range(len(cluster)):
                    amb = self.groups[cluster[id]] 
                    residue = getattr(amb,"residue")
                    hdef = getattr(amb,"hdef")
                    type = hdef.type
                    if type in [1,4,3,10,13]:
                        states = self.getstates(amb)
                        statemap[id] = states
                        self.switchstate(states, amb, 0)
                        curmap[id] = 0
                        bestmap[id] = 0                       
                    elif type in [2,11]:
                        chinum = hdef.chiangle - 1
                        oldangle = residue.get("chiangles")[chinum]
                        curmap[id] = oldangle
                        bestmap[id] = oldangle

                    if getattr(amb,"nearatoms") == []:
                        amb.setNearatoms(allatoms)
                        energy += self.getHbondEnergy(amb)
                  
                maxenergy = energy + 1000
                bestenergy = energy
                newenergy = energy

                self.debug("Initial Best energy: %.2f" % bestenergy)
                self.debug("Initial Cur energy: %.2f" % energy)
                self.debug("Initial Max energy: %.2f" % maxenergy)
                self.debug("Initial Current map: %s" % curmap)

                # Now go through the steps

                for step in range(steps):
                    #self.debug("\n************************")
                    #self.debug("Current map: %s" %  curmap)
                    #self.debug("Current energy: %.2f" % energy)
                    #self.debug("Maximum energy: %.2f" % maxenergy)
                    #self.debug("Trying step %i out of %i" % (step, steps))
                    id = randint(0, len(cluster) - 1)
                    amb = self.groups[cluster[id]]
                    residue = getattr(amb,"residue")
                    hdef = getattr(amb,"hdef")
                    type = hdef.type
                    states = []

                    #oldenergy = self.getHbondEnergy(clusteratoms, compatoms, residue)
                    oldenergy = self.getHbondEnergy(amb)
                    
                    newstate = None
                    if type in [1,4,3,10,13]:
                        states = statemap[id]
                        newstate = randint(0, len(states) - 1)
                        while newstate == curmap[id] and type != 3: #Don't waste a step switching to same state
                            newstate = randint(0, len(states) - 1)
                        self.switchstate(states, amb, newstate)
                    elif type in [2]:                      
                        name = residue.get("name")
                        defresidue = self.routines.aadef.getResidue(name)
                        chinum = hdef.chiangle - 1
                        newstate = randint(0,71)*5.0 - 180 
                        self.routines.setChiangle(residue, chinum, newstate, defresidue)
                    elif type in [11]:
                        name = residue.get("name")
                        defresidue = self.routines.aadef.getResidue(name)
                        chinum = hdef.chiangle - 1
                        oldangle = residue.get("chiangles")[chinum]
                        newstate = 180.0 + oldangle
                        self.routines.setChiangle(residue, chinum, newstate, defresidue)
              
                    #self.debug("Trying to change amb %i to new state %.2f" % (cluster[id], newstate))

                    # Evaluate the change

                    #newenergy = energy - oldenergy + self.getHbondEnergy(clusteratoms, compatoms, residue)
                    newenergy = energy - oldenergy + self.getHbondEnergy(amb)
                    ediff = newenergy - energy

                    rejected = 0
                    if ediff > 0.0 and ediff < 50.0:
                        exp = math.exp(-1 * ediff)
                        if random() >= exp or newenergy > maxenergy:
                            rejected = 1
                    elif ediff >= 50.0:
                        rejected = 1

                    if rejected == 1: # Switch Back
                        #self.debug("Rejected!")
                        if type in [1,4,3,10,13]:
                            self.switchstate(states, amb, curmap[id])
                        elif type in [2,11]:
                            self.routines.setChiangle(residue, chinum, curmap[id], defresidue)
                    else:
                        #self.debug("Accepted!")
                        energy = newenergy
                        curmap[id] = newstate
                        if energy < bestenergy:
                            bestenergy = energy
                            for cid in range(len(cluster)):
                                bestmap[cid] = curmap[cid]
                                
                    # Cool the system
           
                    if step > steps/2 and energy < (bestenergy + 5.0):
                        maxenergy = bestenergy + 5

                # Now switch to best energy
                self.debug("Final state map: %s" % curmap)
                for id in range(len(cluster)):
                    if bestmap[id] == curmap[id]: continue
                    amb = self.groups[cluster[id]]
                    residue = getattr(amb,"residue")
                    hdef = getattr(amb,"hdef")
                    type = hdef.type
                    if type in [1,4,3,10,13]:
                        newstate = bestmap[id]
                        states = statemap[id]
                        self.switchstate(states, amb, newstate)
                    elif type in [2,11]:
                        name = residue.get("name")
                        defresidue = self.routines.aadef.getResidue(name)
                        chinum = hdef.chiangle - 1
                        self.routines.setChiangle(residue, chinum, bestmap[id], defresidue)
                    
                #finalenergy = self.getHbondEnergy(clusteratoms, compatoms)        
                self.debug("Final Best energy: %.2f" % bestenergy)
                #self.debug("Final Actual energy: %.2f" % finalenergy)
                self.debug("Best state map: %s" % bestmap)
                self.debug("*******************\n")

        self.debug("Total time %.2f" % (time() - starttime))
        self.liststates()

    def liststates(self):
        """
            List the final results of all conformation/protonation
            ambiguities to stdout.
        """
        
        for i in range(len(self.groups)):
            state = ""
            amb = self.groups[i]
            residue = getattr(amb,"residue")
            hdef = getattr(amb,"hdef")
            resname = residue.get("name")
            if residue.get("isNterm"):
                H2 = residue.getAtom("H2")
                H3 = residue.getAtom("H3")
                if H2 != None and H3 != None:
                    state = "Positively charged N-terminus"
                elif H2 != None or H3 != None:
                    state = "Neutral N-terminus"
                else:
                    state = "N TERMINUS ERROR!"

            elif residue.get("isCterm") and hdef.name == "CTR":
                HO = residue.getAtom("HO")
                if HO != None:
                    state = "Neutral C-Terminus"
                else:
                    state = "Negative C-Terminus"

            elif (resname == "HIS" and hdef.name != "HISFLIP") or \
                  resname in ["HSP","HSE","HSD","HID","HIE","HIP"]:
                HD1 = residue.getAtom("HD1")
                HE2 = residue.getAtom("HE2")
                if HD1 != None and HE2 != None:
                    state = "Positive HIS"
                elif HD1 != None:
                    state = "HIS HD1"
                elif HE2 != None:
                    state = "HIS HE2"
                else:
                    state = "Negative HIS"

            elif resname == "ASP":
                HD1 = residue.getAtom("HD1")
                HD2 = residue.getAtom("HD2")
                if HD1 != None or HD2 != None:
                    state = "Neutral ASP"
                else:
                    state = "Negative ASP"

            elif resname == "GLU":
                HE1 = residue.getAtom("HE1")
                HE2 = residue.getAtom("HE2")
                if HE1 != None or HE2 != None:
                    state = "Neutral GLU"
                else:
                    state = "Negative GLU"

            elif resname == "GLH":
                HE1 = residue.getAtom("HE1")
                HE2 = residue.getAtom("HE2")
                if HE1 != None: state = "Neutral GLU (HE1)"
                elif HE2 != None: state = "Neutral GLU (HE2)"
                else: raise ValueError, "GLH should always be neutral!"

            if state != "":
                self.routines.write("Ambiguity #: %i, chain: %s, residue: %i %s - %s\n" \
                                    % (i, residue.chainID, residue.resSeq, residue.name, state),1)

    def getHbondEnergy2(self,clusteratoms, compatoms, res=None):
        """
            Get the hydrogen bond energy for a cluster of donors
            and acceptors. If res is not present, compare each atom
            in clusteratoms to all nearby atoms in compatoms.  If res is
            present, we are trying to find the new energy of the residue
            res that has just switched states, and thus need to include
            all comparisons where atoms within the residue is an
            acceptor.

            Parameters
                clusteratoms: A list of hydrogen donor/acceptor atoms in
                              the cluster(list)
                compatoms:    A list of all hydrogen donor/acceptor
                              atoms within a given distance of the
                              cluster (list)
                res:          (Optional) Residue to get the energy of
                              (Residue)
            Returns
                energy:       The energy of this hydrogen bond network
                              (float)
        """
        energy = 0.0
        if res != None:
            for atom2 in compatoms:
                res2 = atom2.get("residue")
                if res2 != res: continue
                elif not atom2.get("hacceptor"): continue
                for atom1 in clusteratoms:        
                    if atom1.get("residue") == res2: continue
                    elif not atom1.get("hdonor"): continue
                    elif distance(atom1.getCoords(), atom2.getCoords()) < HYDROGEN_DIST:
                        pair = self.getPairEnergy(atom1, atom2)
                        energy = energy + pair
                        #print "\tComparing %s %i to %s %i %.2f" % (atom1.name, atom1.residue.resSeq, atom2.name, atom2.residue.resSeq, pair)
        
        for atom1 in clusteratoms:
            energy = energy + self.getPenalty(atom1)
            res1 = atom1.get("residue")
            if res != None and res1 != res: continue
            elif not res1.get("isNterm") and atom1.get("name") == "N": continue
            elif not atom1.get("hdonor"): continue
            for atom2 in compatoms:
                if res1 == atom2.get("residue"): continue
                elif not atom2.get("hacceptor"): continue
                elif distance(atom1.getCoords(), atom2.getCoords()) < HYDROGEN_DIST:
                    pair = self.getPairEnergy(atom1, atom2)
                    energy = energy + pair
                    #print "\tComparing %s %i to %s %i %.2f" % (atom1.name, atom1.residue.resSeq, atom2.name, atom2.residue.resSeq, pair)
        return energy

    def getHbondEnergy(self, amb):
        """
        """
        energy = 0.0
        residue = getattr(amb,"residue")
        nearatoms = getattr(amb,"nearatoms")
        energy = energy + self.getPenalty(residue)
        for nearatom in nearatoms:
            for atom in residue.get("atoms"):
                if atom.get("hdonor"):
                    if not nearatom.get("hacceptor"): continue
                    elif atom.get("name") == "N" and not residue.get("isNterm"): continue
                    elif nearatom.get("name") == "O" and not nearatom.residue.get("name") == "WAT": continue
                    if distance(atom.getCoords(), nearatom.getCoords()) < HYDROGEN_DIST:                
                        pair = self.getPairEnergy(atom, nearatom)
                        energy = energy + pair
                if atom.get("hacceptor"):
                    if not nearatom.get("hdonor"): continue
                    elif atom.get("name") == "O" and not residue.get("name") == "WAT": continue
                    elif nearatom.get("name") == "N" and not nearatom.residue.get("isNterm"): continue
                    if distance(atom.getCoords(), nearatom.getCoords()) < HYDROGEN_DIST:             
                        pair = self.getPairEnergy(nearatom, atom)
                        energy = energy + pair
              
        return energy
                

    def getPairEnergy(self, donor, acceptor):
        """
            Get the energy between two atoms

            Parameters
                donor:    The first atom in the pair (Atom)
                acceptor: The second atom in the pair (Atom)
            Returns
                energy:   The energy of the pair (float)
        """
        max_hbond_energy = -10.0
        max_ele_energy = -1.0
        maxangle = 20.0
        max_dha_dist = 3.3
        max_ele_dist = 5.0
        energy = 0.0
        donorh = None
        acceptorh = None
        donorhs = []
        acceptorhs = []
        
        if not (donor.get("hdonor") and acceptor.get("hacceptor")):
            return energy

        # See if hydrogens are presently bonded to the acceptor and donor

        for bond in donor.get("intrabonds"):
            if bond[0] == "H":
                donorhs.append(bond)
        for bond in acceptor.get("intrabonds"):
            if bond[0] == "H":
                acceptorhs.append(bond)
            
        if donorhs == []: return energy

        # Case 1: Both donor and acceptor hydrogens are present

        if acceptorhs != []:
            for donorh in donorhs:
                for acceptorh in acceptorhs:
                    donorhatom = donor.get("residue").getAtom(donorh)
                    acceptorhatom = acceptor.get("residue").getAtom(acceptorh)
                    if donorhatom == None or acceptorhatom == None:
                        text = "Couldn't find bonded hydrogen even though "
                        text = text + "it is present in intrabonds!"
                        raise ValueError, text
            
                    dist = distance(donorhatom.getCoords(), acceptor.getCoords())
                    if dist > max_dha_dist and dist < max_ele_dist: # Are the Hs too far?
                        energy += max_ele_energy/(dist*dist)
                        continue
                    
                    hdist = distance(donorhatom.getCoords(), acceptorhatom.getCoords())
                    if hdist < 1.5: # Are the Hs too close?
                        energy += -1 * max_hbond_energy
                        continue
                    
                    angle1 = self.getHbondangle(acceptor, donor, donorhatom)
                    if angle1 <= maxangle:
                        angleterm = (maxangle - angle1)/maxangle
                        angle2 = self.getHbondangle(donorhatom, acceptorhatom, acceptor)
                        if angle2 < 110.0: angle2 = 1.0
                        else: angle2=-1.0/110.0*(angle2-110.0)
                        energy+=max_hbond_energy/pow(dist,3)*angleterm*angle2
                    
            return energy

        # Case 2: Only the donor hydrogen is present

        elif acceptorhs == []:
            for donorh in donorhs:
                donorhatom = donor.get("residue").getAtom(donorh)
                if donorhatom == None:
                    text = "Couldn't find bonded hydrogen even though "
                    text = text + "it is present in intrabonds!"
                    raise ValueError, text
                
                dist = distance(donorhatom.getCoords(), acceptor.getCoords())
                if dist > max_dha_dist and dist < max_ele_dist: # Or too far?
                    energy += max_ele_energy/(dist*dist)               
                    continue
              
                angle1 = self.getHbondangle(acceptor, donor, donorhatom)
                if angle1 <= maxangle:
                    angleterm = (maxangle - angle1)/maxangle
                    energy += max_hbond_energy/pow(dist,2)*angleterm
                  
            return energy

    def getHbondangle(self, atom1, atom2, atom3):
        """
            Get the angle between three atoms

            Parameters
                atom1:  The first atom (atom)
                atom2:  The second (vertex) atom (atom)
                atom3:  The third atom (atom)
            Returns
                angle:  The angle between the atoms (float)
        """
        angle = 0.0
        atom2Coords = atom2.getCoords()
        coords1 = subtract(atom3.getCoords(), atom2Coords)
        coords2 = subtract(atom1.getCoords(), atom2Coords)
        norm1 = normalize(coords1)
        norm2 = normalize(coords2)
        dotted = dot(norm1, norm2)
        if dotted > 1.0: # If normalized, this is due to rounding error
            dotted = 1.0
        rad = abs(acos(dotted))
        angle = rad*180.0/pi
        if angle > 180.0:
            angle = 360.0 - angle
        return angle

    def getPenalty(self, residue):
        """
            Add penalties for unusual protonation states.

            Parameters
                atom:    The residue to examine (Atom)
            Returns
                penalty: The amount of the penalty (float)
        """
        acidpenalty = 25.0
        hispos = 0.1
        hisminus = 10.0
        nterm = 5.0
        penalty = 0.0

        resname = residue.get("name")

        if residue.get("isNterm"):
            charge = 1
            if residue.getAtom("H2") == None: charge = charge - 1
            if residue.getAtom("H3") == None: charge = charge - 1
            penalty = penalty + (1- charge)*nterm

        if resname == "HIS":
            hd1 = residue.getAtom("HD1")
            he2 = residue.getAtom("HE2")
            if hd1 == None and he2 == None: 
                penalty = penalty + hisminus
            elif hd1 != None and he2 != None:
                penalty = penalty + hispos

        if resname != "WAT":
            for atom in residue.get("atoms"):
                atomname = atom.get("name")
                if atomname in ["OD1","OD2","OE1","OE2","O","OXT"] and atom.get("hdonor"):
                    penalty = penalty + acidpenalty
                    break

        return penalty
                
                
    def initHbondEnergy(self, cluster, allatoms):
        """
            Create a list of hydrogen donors/acceptors within this cluster
            and another list of donors/acceptors throughout the
            entire protein.

            Parameters
                cluster:      A list of group ids that are networked (list)
            Returns
                clusteratoms: A list of hydrogen donor/acceptor atoms in
                              the cluster (list)
        """
        clusteratoms = []
        compatoms = []
        for id in cluster:
            residue = self.groups[id][0]
            for atom in residue.get("atoms"):
                if (atom.get("hacceptor") or atom.get("hdonor")) \
                       and atom not in clusteratoms:
                    clusteratoms.append(atom)
                    for atom2 in allatoms:
                        dist = distance(atom.getCoords(), atom2.getCoords())
                        if dist < HYDROGEN_DIST and atom2 not in compatoms:
                            compatoms.append(atom2)
        return clusteratoms, compatoms
            
    def findNetworks(self,limit):
        """
            Find hydrogen networks that should be optimized together.

            Parameters
                limit:    The limit to see how close two boundatoms are

            Returns
                networks: A list of group ID networks (list)
        """
        map = {}
        networks = []
        done = []
        groups = self.groups
        for i in range(len(groups)):
            amb1 = groups[i]
            residue1 = getattr(amb1, "residue")
            hydrodef1 = getattr(amb1,"hdef")
            for conf1 in hydrodef1.conformations:
                boundatom1 = residue1.getAtom(conf1.boundatom)
                for j in range(i+1, len(groups)):
                    amb2 = groups[j]
                    residue2 = getattr(amb2, "residue")
                    hydrodef2 = getattr(amb2,"hdef")
                    for conf2 in hydrodef2.conformations:
                        boundatom2 = residue2.getAtom(conf2.boundatom)
                        if distance(boundatom1.getCoords(), boundatom2.getCoords()) < limit:
                            if i not in map:
                                map[i] = [j]
                            elif j not in map[i]:
                                map[i].append(j)
                            break

        for i in range(len(groups)):
            if i in map and i not in done:
                list = analyzeMap(map,i,[])
                for item in list:
                    done.append(item)
                networks.append(list)
            elif i not in map and i not in done:
                networks.append([i])

        self.debug(networks)
        return networks

    def randomizeWaters(self):
        """
            Randomize the waters found in a protein. Mimics the
            optimizeWaters function, but instead of going through
            all possible 5 degree increments, simply choose a random
            angle.
        """
        allatoms = self.findAmbiguities(1)
        closeatoms = {}
        overallenergy = 0.0

        # Step 1: Get list of water residues

        waters = []
        for group in self.groups:
            residue = getattr(group,"residue")
            waters.append(residue)

        # Step 2: Shuffle waters

        shuffle(waters)
        
        # Step 3: Satisfy all protein donors
    
        for residue in waters:
            self.watermap[residue] = None
            closeatoms[residue] = []
            bestdon = None       
            bestdondist = 999.99
            bestnh = None
            oxygen = residue.getAtom("O")
            for atom in allatoms:
                closedist = distance(oxygen.getCoords(), atom.getCoords())
                if oxygen != atom and closedist < WATER_DIST: closeatoms[residue].append(atom) 
                if atom.get("residue").name == "WAT": continue  
                if atom.get("hdonor"):
                    for bond in atom.get("intrabonds"):
                        if bond[0] != "H": continue
                        bondatom = atom.get("residue").getAtom(bond)
                        dist = distance(oxygen.getCoords(), bondatom.getCoords())
                        angle = self.getHbondangle(atom, oxygen, bondatom)
                        if dist < bestdondist and angle <= 20.0 : 
                            bestdon = bondatom
                            bestdondist = dist
                elif atom.get("name")[0:2] == "NH":
                    arg = atom.get("residue")
                    if atom.get("name") == "NH1": other = "NH2"
                    else: other = "NH1"
                    if arg.getAtom(other).get("hdonor") != 1:
                        for bond in atom.get("intrabonds"):
                            if bond[0] != "H": continue
                            bondatom = arg.getAtom(bond)
                            dist = distance(oxygen.getCoords(), bondatom.getCoords())
                            angle = self.getHbondangle(atom, oxygen, bondatom)
                            if dist < bestdondist and angle <= 20.0:
                                bestdon = bondatom
                                bestdondist = dist
                                bestnh = atom

            if bestnh != None:
                if bestdon.get("name") in bestnh.get("intrabonds"):
                    bestnh.set("hdonor",1)
                    bestnh.set("hacceptor",0)
          
            if bestdondist < 3.3:
                #print "Optimizing WAT %i" % residue.resSeq
                #print "\tBest Donor: ", bestdon.residue.name, bestdon.residue.resSeq, bestdon.name, bestdondist
                R = bestdondist
                refcoords, defcoords = [], []
          
                defcoords.append([0,0,.3333])   # Oxygen
                defcoords.append([0,0,.3333+R]) # Donor
                refcoords.append(oxygen.getCoords())
                refcoords.append(bestdon.getCoords())

                defatomcoords = [.9428,0,0] # Location 1
                newcoords = findCoordinates(2, refcoords, defcoords, defatomcoords)
                residue.createAtom("H1",newcoords,"HETATM")

                defatomcoords = [-.4714,.8165,0] # Location 2
                newcoords = findCoordinates(2, refcoords, defcoords, defatomcoords)
                residue.createAtom("H2",newcoords,"HETATM")

                oxygen.intrabonds = ["H1","H2"]
                self.watermap[residue] = bestdon

                # Now randomize

                randangle = randint(0,71)*5.0 - 180
                self.setWaterHydrogens(residue, randangle)

            else: closeatoms[residue] = []

        # Step 4: Place and orient hydrogens

        for residue in waters:
            #print "Optimizing WAT %i" % residue.resSeq
            oxygen = residue.getAtom("O")
            if self.watermap[residue] != None: continue
            bestdon = None
            bestacc = None        
            bestdondist = 999.99
            bestaccdist = 999.99
            for atom in allatoms:
                if atom == oxygen: continue
                closedist = distance(oxygen.getCoords(), atom.getCoords())
                if closedist < WATER_DIST: closeatoms[residue].append(atom)
                if atom.get("hacceptor"):
                    dist = closedist
                    if dist < bestaccdist:
                        bestacc = atom
                        bestaccdist = dist
                if atom.get("hdonor"):
                    for bond in atom.get("intrabonds"):
                        if bond[0] != "H": continue
                        bondatom = atom.get("residue").getAtom(bond)
                        dist = distance(oxygen.getCoords(), bondatom.getCoords())      
                        if dist < bestdondist:
                            bestdon = bondatom
                            bestdondist = dist

            #print "\tBest Donor: ", bestdon.residue.name, bestdon.residue.resSeq, bestdon.name, bestdondist
            #print "\tBest Acc: ", bestacc.residue.name, bestacc.residue.resSeq, bestacc.name, bestaccdist
    
            if bestdondist < bestaccdist: # Donor is closest
                R = bestdondist
                refcoords, defcoords = [], []
                
                defcoords.append([0,0,.3333])   # Oxygen
                defcoords.append([0,0,.3333+R]) # Donor
                refcoords.append(oxygen.getCoords())
                refcoords.append(bestdon.getCoords())
                
                defatomcoords = [.9428,0,0] # Location 1
                newcoords = findCoordinates(2, refcoords, defcoords, defatomcoords)
                residue.createAtom("H1",newcoords,"HETATM")
                
                defatomcoords = [-.4714,.8165,0] # Location 2
                newcoords = findCoordinates(2, refcoords, defcoords, defatomcoords)
                residue.createAtom("H2",newcoords,"HETATM")
                
                oxygen.intrabonds = ["H1","H2"]
                self.watermap[residue] = bestdon

            else: # Acceptor is closest
                b = 1.0  # Oxygen (Donor) - H1 dist
                r = 1.5  # Acceptor - H1 dist
                if bestaccdist >= (b + r):
                    # Then H1 is perfectly placed on the vector
                    vec = subtract(bestacc.getCoords(), oxygen.getCoords())
                    x = oxygen.get("x") + b/bestaccdist * vec[0]
                    y = oxygen.get("y") + b/bestaccdist * vec[1]
                    z = oxygen.get("z") + b/bestaccdist * vec[2]
                    newcoords = [x,y,z]
                    residue.createAtom("H1",newcoords,"HETATM")
                    
                else:
                    # Minimize the H-O-A angle 
                    defcoords, refcoords = [], []
                    R = distance(oxygen.getCoords(), bestacc.getCoords())
                    psi = acos((b*b + R*R - r*r)/(2*b*R))
                    
                    defcoords.append([0,0,0])
                    refcoords.append(oxygen.getCoords())
                    defcoords.append([R,0,0])
                    refcoords.append(bestacc.getCoords())
                    
                    y = random()
                    while y > sin(psi):
                        y = y/2
                    z = sqrt(sin(psi)*sin(psi) - y*y)
                        
                    defatomcoords = [cos(psi), y, z]
                    newcoords = findCoordinates(2, refcoords, defcoords, defatomcoords)
                    residue.createAtom("H1",newcoords,"HETATM")

                defcoords, refcoords = [], []
                defcoords.append([0,0,.3333])   # Oxygen
                defcoords.append([.9428,0,0])   # H1
                refcoords.append(oxygen.getCoords())
                refcoords.append(newcoords)
                defatomcoords = [-.4714,.8165,0] # Location 2
                h2coords = findCoordinates(2, refcoords, defcoords, defatomcoords)
                residue.createAtom("H2",h2coords,"HETATM")
                oxygen.intrabonds = ["H1","H2"]
                self.watermap[residue] = residue.getAtom("H1")

            # Now randomize
                    
            randangle = randint(0,71)*5.0 - 180
            self.setWaterHydrogens(residue, randangle)
            
    def setWaterHydrogens(self, residue, newangle):
        """
            Optimize a Water molecule

            Parameters
                residue:  The water residue
                newangle: The new chi angle (float)
        """
        movenames = []
        movecoords = []

        hydrogen = self.watermap[residue]
        oxygen = residue.getAtom("O")
        
        initcoords = subtract(oxygen.getCoords(), hydrogen.getCoords())

        # Determine which atoms to rotate

        if hydrogen not in residue.get("atoms"):
            movenames = ["H1","H2"]
        elif hydrogen.get("name") == "H1":
            movenames = ["H2"]
        else:
            raise ValueError, "Got improperly mapped water hydrogen!"
            
        for name in movenames:
            atom = residue.getAtom(name)
            movecoords.append(subtract(atom.getCoords(), hydrogen.getCoords()))

        newcoords = qchichange(initcoords, movecoords, newangle)
    
        for i in range(len(movenames)):
            name = movenames[i]
            atom = residue.getAtom(name)
            self.routines.removeCell(atom)
            x = (newcoords[i][0] + hydrogen.get("x"))
            y = (newcoords[i][1] + hydrogen.get("y"))
            z = (newcoords[i][2] + hydrogen.get("z"))
            atom.set("x", x)
            atom.set("y", y)
            atom.set("z", z)
            self.routines.addCell(atom)
            
    def optimizeWaters(self):
        """
            Optimize the waters found in a protein
        """
        allatoms = self.findAmbiguities(1)
        closeatoms = {}
        overallenergy = 0.0

        # Step 1: Get list of water residues

        waters = []
        for group in self.groups:
            residue = getattr(group,"residue")
            waters.append(residue)

        # Step 2: Shuffle waters

        shuffle(waters)
        
        # Step 3: Satisfy all protein donors
    
        for residue in waters:
            self.watermap[residue] = None
            closeatoms[residue] = []
            bestdon = None       
            bestdondist = 999.99
            bestnh = None
            oxygen = residue.getAtom("O")
            for atom in allatoms:
                closedist = distance(oxygen.getCoords(), atom.getCoords())
                if oxygen != atom and closedist < WATER_DIST: closeatoms[residue].append(atom) 
                if atom.get("residue").name == "WAT": continue  
                if atom.get("hdonor"):
                    for bond in atom.get("intrabonds"):
                        if bond[0] != "H": continue
                        bondatom = atom.get("residue").getAtom(bond)
                        dist = distance(oxygen.getCoords(), bondatom.getCoords())
                        angle = self.getHbondangle(atom, oxygen, bondatom)
                        if dist < bestdondist and angle <= 20.0 : 
                            bestdon = bondatom
                            bestdondist = dist
                elif atom.get("name")[0:2] == "NH":
                    arg = atom.get("residue")
                    if atom.get("name") == "NH1": other = "NH2"
                    else: other = "NH1"
                    if arg.getAtom(other).get("hdonor") != 1:
                        for bond in atom.get("intrabonds"):
                            if bond[0] != "H": continue
                            bondatom = arg.getAtom(bond)
                            dist = distance(oxygen.getCoords(), bondatom.getCoords())
                            angle = self.getHbondangle(atom, oxygen, bondatom)
                            if dist < bestdondist and angle <= 20.0:
                                bestdon = bondatom
                                bestdondist = dist
                                bestnh = atom

            if bestnh != None:
                if bestdon.get("name") in bestnh.get("intrabonds"):
                    bestnh.set("hdonor",1)
                    bestnh.set("hacceptor",0)
          
            if bestdondist < 3.3:
                #print "Optimizing WAT %i" % residue.resSeq
                #print "\tBest Donor: ", bestdon.residue.name, bestdon.residue.resSeq, bestdon.name, bestdondist
                R = bestdondist
                refcoords, defcoords = [], []
          
                defcoords.append([0,0,.3333])   # Oxygen
                defcoords.append([0,0,.3333+R]) # Donor
                refcoords.append(oxygen.getCoords())
                refcoords.append(bestdon.getCoords())

                defatomcoords = [.9428,0,0] # Location 1 on tetrahedral
                newcoords = findCoordinates(2, refcoords, defcoords, defatomcoords)
                residue.createAtom("H1",newcoords,"HETATM")

                defatomcoords = [-.4714,.8165,0] # Location 2 on tetrahedral
                newcoords = findCoordinates(2, refcoords, defcoords, defatomcoords)
                residue.createAtom("H2",newcoords,"HETATM")

                oxygen.intrabonds = ["H1","H2"]
                self.watermap[residue] = bestdon

                # Now optimize

                amb = hydrogenAmbiguity(residue,None)
                setattr(amb,"nearatoms",closeatoms[residue])
                #energy = self.getHbondEnergy([oxygen], closeatoms[residue])
                energy = self.getHbondEnergy(amb)
                bestangle = None
                bestenergy = energy
                
                for angle in range(-180,180,5):
                    self.setWaterHydrogens(residue, angle)
                    #energy = self.getHbondEnergy([oxygen], closeatoms[residue])
                    energy = self.getHbondEnergy(amb)
                    if energy < bestenergy:
                        bestenergy = energy
                        bestangle = angle
             
                if bestangle == None:
                    bestangle = randint(0,71)*5.0 - 180
                self.setWaterHydrogens(residue, bestangle)
                overallenergy += bestenergy

            else: closeatoms[residue] = []

        # Step 4: Place and orient hydrogens

        for residue in waters:
            #print "Optimizing WAT %i" % residue.resSeq
            oxygen = residue.getAtom("O")
            if self.watermap[residue] != None: continue
            bestdon = None
            bestacc = None        
            bestdondist = 999.99
            bestaccdist = 999.99
            for atom in allatoms:
                if atom == oxygen: continue
                closedist = distance(oxygen.getCoords(), atom.getCoords())
                if closedist < WATER_DIST: closeatoms[residue].append(atom)
                if atom.get("hacceptor"):
                    dist = closedist
                    if dist < bestaccdist:
                        bestacc = atom
                        bestaccdist = dist
                if atom.get("hdonor"):
                    for bond in atom.get("intrabonds"):
                        if bond[0] != "H": continue
                        bondatom = atom.get("residue").getAtom(bond)
                        dist = distance(oxygen.getCoords(), bondatom.getCoords())      
                        if dist < bestdondist:
                            bestdon = bondatom
                            bestdondist = dist

            #print "\tBest Donor: ", bestdon.residue.name, bestdon.residue.resSeq, bestdon.name, bestdondist
            #print "\tBest Acc: ", bestacc.residue.name, bestacc.residue.resSeq, bestacc.name, bestaccdist
    
            if bestdondist < bestaccdist: # Donor is closest
                R = bestdondist
                refcoords, defcoords = [], []
                
                defcoords.append([0,0,.3333])   # Oxygen
                defcoords.append([0,0,.3333+R]) # Donor
                refcoords.append(oxygen.getCoords())
                refcoords.append(bestdon.getCoords())
                
                defatomcoords = [.9428,0,0] # Location 1 on tetrahedral
                newcoords = findCoordinates(2, refcoords, defcoords, defatomcoords)
                residue.createAtom("H1",newcoords,"HETATM")
                
                defatomcoords = [-.4714,.8165,0] # Location 2 on tetrahedral
                newcoords = findCoordinates(2, refcoords, defcoords, defatomcoords)
                residue.createAtom("H2",newcoords,"HETATM")
                
                oxygen.intrabonds = ["H1","H2"]
                self.watermap[residue] = bestdon

            else: # Acceptor is closest
                b = 1.0  # Oxygen (Donor) - H1 dist
                r = 1.5  # Acceptor - H1 dist
                if bestaccdist >= (b + r):
                    # Then H1 is perfectly placed on the vector
                    vec = subtract(bestacc.getCoords(), oxygen.getCoords())
                    x = oxygen.get("x") + b/bestaccdist * vec[0]
                    y = oxygen.get("y") + b/bestaccdist * vec[1]
                    z = oxygen.get("z") + b/bestaccdist * vec[2]
                    newcoords = [x,y,z]
                    residue.createAtom("H1",newcoords,"HETATM")
                    
                else:
                    # Minimize the H-O-A angle 
                    defcoords, refcoords = [], []
                    R = distance(oxygen.getCoords(), bestacc.getCoords())
                    psi = acos((b*b + R*R - r*r)/(2*b*R))
                    
                    defcoords.append([0,0,0])
                    refcoords.append(oxygen.getCoords())
                    defcoords.append([R,0,0])
                    refcoords.append(bestacc.getCoords())
                    
                    y = random()
                    while y > sin(psi):
                        y = y/2
                    z = sqrt(sin(psi)*sin(psi) - y*y)
                        
                    defatomcoords = [cos(psi), y, z]
                    newcoords = findCoordinates(2, refcoords, defcoords, defatomcoords)
                    residue.createAtom("H1",newcoords,"HETATM")

                defcoords, refcoords = [], []
                defcoords.append([0,0,.3333])   # Oxygen
                defcoords.append([.9428,0,0])   # H1
                refcoords.append(oxygen.getCoords())
                refcoords.append(newcoords)
                defatomcoords = [-.4714,.8165,0] # Location 2
                h2coords = findCoordinates(2, refcoords, defcoords, defatomcoords)
                residue.createAtom("H2",h2coords,"HETATM")
                oxygen.intrabonds = ["H1","H2"]
                self.watermap[residue] = residue.getAtom("H1")

            # Now optimize

            amb = hydrogenAmbiguity(residue,None)
            setattr(amb,"nearatoms",closeatoms[residue])

            energy = self.getHbondEnergy(amb)
            #energy = self.getHbondEnergy([oxygen], closeatoms[residue])
            bestangle = None
            bestenergy = energy

            for angle in range(-180,180,5):
                self.setWaterHydrogens(residue, angle)

                energy = self.getHbondEnergy(amb)
                #energy = self.getHbondEnergy([oxygen], closeatoms[residue])
    
                if energy < bestenergy:
                    bestenergy = energy
                    bestangle = angle
             
            if bestangle == None:
                bestangle = randint(0,71)*5.0 - 180
            self.setWaterHydrogens(residue, bestangle)
            overallenergy += bestenergy
       
        #print "Overall energy: %.2f" % overallenergy
             
    def findViableAngles(self, residue, nearatoms):
        """
            Find the viable angles that a water molecule can be rotated to.
            If there are no donor/acceptor atoms within 

            Parameters
                residue:   The water residue to examine
                nearatoms: A list of nearby donor/acceptors
            Returns
                angle:   A list of viable angles 
        """
        angles = []
        bestmap = {} # Store best values by tuple (nearatom, hydrogen)
        
        hydrogen = self.watermap[residue]
        if hydrogen not in residue.get("atoms"):
            moveatoms = [residue.getAtom("H1"),residue.getAtom("H2")]
        elif hydrogen.get("name") == "H1":
            moveatoms = [residue.getAtom("H2")]

        for atom in nearatoms:
            for moveableatom in moveatoms:
                bestmap[(atom, moveableatom)] = (999.99, None)
    
        for angle in range(-180,180,5):
            self.optimizeWaterHydrogens(residue, angle)                
            for atom in nearatoms:
                for moveableatom in moveatoms:
                    dist = distance(moveableatom.getCoords(), atom.getCoords())
                    bestdist = bestmap[(atom, moveableatom)][0]
                    if dist < bestdist:
                        bestmap[(atom, moveableatom)] = (dist, angle)

        for atom in nearatoms:
            for moveableatom in moveatoms:
                angle = bestmap[(atom,moveableatom)][1]
                if angle not in angles: angles.append(angle)
                
        return angles
    
    def findAmbiguities(self, water):
        """
            Find the amibiguities within a protein according to the
            DAT file, and set all boundatoms to their hydrogen donor/
            acceptor state.  Store the ambiguities as (residue, hydrodef)
            tuples in self.groups.

            Returns
                allatoms:  A list of all donors and acceptors in the
                           protein (list)
                water:     If 1, only put waters in groups, but fill allatoms
                           appropriately
        """
        allatoms = []
        hydrodefs = self.hydrodefs
        for chain in self.protein.getChains():
            for residue in chain.get("residues"):
                resname = residue.get("name")
                nter = residue.get("isNterm")
                cter = residue.get("isCterm")
                type = residue.get("type")
                if type == 2: continue
                for group in hydrodefs:
                    groupname = group.name
                    htype = group.type
                    if resname == groupname or \
                       (groupname == "APR") or \
                       (groupname == "APP" and resname != "PRO") or \
                       (groupname.endswith("FLIP") and resname == groupname[0:3]) or \
                       (groupname == "HISFLIP" and resname in ["HIP","HID","HIE","HSP","HSE","HSD"]) or \
                       (groupname == "NTR" and nter and resname != "PRO") or \
                       (groupname == "PNTR" and nter and resname == "PRO") or \
                       (groupname == "CTR" and cter):

                        if group.method != 0:
                            if not water and type == 1:
                                amb = hydrogenAmbiguity(residue, group)
                                self.groups.append(amb)
                            if water and type == 3:
                                amb = hydrogenAmbiguity(residue, group)
                                self.groups.append(amb)

                        for conf in group.conformations:
                            boundatom = conf.boundatom
                            atom = residue.getAtom(boundatom)
                            if atom != None:
                                atom.set("hacceptor",HACCEPTOR[htype])
                                atom.set("hdonor",HDONOR[htype])   
                                if atom not in allatoms:
                                    allatoms.append(atom)
                                    
        return allatoms
                                
    def printAmbiguities(self):
        """
            Print the list of ambiguities to stdout
        """
        if self.hdebug == 0: return
        i = 0
        for amb in self.groups:
            residue = getattr(amb,"residue")
            hydrodef = getattr(amb,"hdef")
            conf = hydrodef.conformations[0]
            self.routines.write("Ambiguity #: %i, chain: %s, residue: %i %s, hyd_type: %i state: %i Grp_name: %s Hname: %s Boundatom: %s\n" % (i, residue.chainID, residue.resSeq, residue.name, hydrodef.type, hydrodef.standardconf, hydrodef.name, conf.hname, conf.boundatom))
            i += 1
            
    def parseHydrogen(self, lines):
        """
            Parse a list of lines in order to make a hydrogen
            definition

            Parameters
                lines:  The lines to parse (list)
            Returns
                mydef:  The hydrogen definition object (HydrogenDefinition)
        """
        maininfo = string.split(lines[0])
        name = maininfo[0]
        group = maininfo[1]
        numhydrogens = int(maininfo[2])
        standardconf = int(maininfo[4])
        type = int(maininfo[5])
        chiangle = int(maininfo[6])
        method = int(maininfo[7])

        mydef = HydrogenDefinition(name, group, numhydrogens, standardconf, \
                                   type, chiangle, method)

        conf = []
        for newline in lines[1:]:
            if newline.startswith(">"):
                if conf == []: continue
                confinfo = string.split(conf[0])
                hname = confinfo[0]
                boundatom = confinfo[1]
                bondlength = float(confinfo[2])
                myconf = HydrogenConformation(hname, boundatom, bondlength)
                count = 0
                for line in conf[1:]:
                    textatom = string.split(line)
                    name = textatom[0]
                    x = float(textatom[1])
                    y = float(textatom[2])
                    z = float(textatom[3])
                    atom = DefinitionAtom(count, name, "", x,y,z)
                    myconf.addAtom(atom)

                mydef.addConf(myconf)
                conf = []
            else:
                conf.append(newline)
                
        if conf != []:
            confinfo = string.split(conf[0])
            hname = confinfo[0]
            boundatom = confinfo[1]
            bondlength = float(confinfo[2])
            myconf = HydrogenConformation(hname, boundatom, bondlength)
            count = 0
            for line in conf[1:]:
                textatom = string.split(line)
                name = textatom[0]
                x = float(textatom[1])
                y = float(textatom[2])
                z = float(textatom[3])
                atom = DefinitionAtom(count, name, "", x,y,z)
                myconf.addAtom(atom)
            mydef.addConf(myconf)

        if lines[1:] == []:  # FLIPS
            myconf = HydrogenConformation(None, "CA", 0.0)
            mydef.addConf(myconf)
            
        return mydef

    def readHydrogenDefinition(self):
        """
            Read the Hydrogen Definition file

            Returns
                hydrodef:  The hydrogen definition ()
        """
        defpath = HYDROGENFILE
        if not os.path.isfile(defpath):
            for path in sys.path:
                testpath = "%s/%s" % (path, defpath)
                if os.path.isfile(testpath):
                    defpath = testpath
                    break
        if not os.path.isfile(defpath):
            raise ValueError, "%s not found!" % defpath

       
        file = open(defpath)
        lines = file.readlines()
        file.close()
        info = []

        for line in lines:
            if line.startswith("//"): pass
            elif line.startswith("*") or line.startswith("!"):
                if info == []: continue
                mydef = self.parseHydrogen(info)
                self.hydrodefs.append(mydef)
                info = []
            else:
                info.append(string.strip(line))

class HydrogenDefinition:
    """
        HydrogenDefinition class

        The HydrogenDefinition class provides information on possible
        ambiguities in amino acid hydrogens.  It is essentially the hydrogen
        definition file in object form.
    """
    
    def __init__(self, name, group, numhydrogens, standardconf, type, \
                 chiangle, method):
        """
            Initialize the object with information from the definition file

            Parameters:
                name:          The name of the grouping (string)
                group:         The group of the definition
                               (acid/base/none, string)
                numhydrogens:  The number of hydrogens that can be added (int)
                standardconf:  The number of standard conformations (int)
                type        :  Type of Hydrogen (int)
                chiangle    :  The chiangle to be changed (int)
                method      :  The standard optimization method (int)

                See HYDROGENS.DAT for more information
        """
        self.name = name
        self.group = group
        self.numhydrogens = numhydrogens
        self.standardconf = standardconf
        self.type = type
        self.chiangle = chiangle
        self.method = method
        self.conformations = []

    def __str__(self):
        """
            Used for debugging purposes

            Returns
                output:  The information about this definition (string)
        """
        output =  "Name:                  %s\n" % self.name
        output += "Group:                 %s\n" % self.group
        output += "# of Hydrogens:        %i\n" % self.numhydrogens
        output += "# of Conformations:    %i\n" % len(self.conformations)
        output += "Standard Conformation: %i\n" % self.standardconf
        output += "Type of Hydrogen:      %i\n" % self.type
        output += "Chiangle to change:    %i\n" % self.chiangle
        output += "Optimization method:   %i\n" % self.method
        output += "Conformations:\n"
        for conf in self.conformations:
            output += "\n%s" % conf
        output += "*****************************************\n"
        return output

    def addConf(self, conf):
        """
            Add a HydrogenConformation to the list of conformations

            Parameters
                conf:  The conformation to be added (HydrogenConformation)
        """
        self.conformations.append(conf)

class HydrogenConformation:
    """
        HydrogenConformation class

        The HydrogenConformation class contains data about possible
        hydrogen conformations as specified in the hydrogen data file.
    """

    def __init__(self, hname, boundatom, bondlength):
        """
           Initialize the object

           Parameters
               hname      : The hydrogen name (string)
               boundatom  : The atom the hydrogen is bound to (string)
               bondlength : The bond length (float)
        """
        self.hname = hname
        self.boundatom = boundatom
        self.bondlength = bondlength
        self.atoms = []

    def __str__(self):
        """
            Used for debugging purposes

            Returns
                output:  Information about this conformation (string)
        """
        output  = "Hydrogen Name: %s\n" % self.hname
        output += "Bound Atom:    %s\n" % self.boundatom
        output += "Bond Length:   %.2f\n" % self.bondlength
        for atom in self.atoms:
            output += "\t%s\n" % atom
        return output
    
    def addAtom(self, atom):
        """
            Add an atom to the list of atoms

            Parameters
                atom: The atom to be added (DefinitionAtom)
        """
        self.atoms.append(atom)
        
