## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#
# Last modified on Thu Aug  9 16:32:21 PDT 2001 by lindy
#
# This is the python translation of cluster_analysis.cc from
# the AutoDock code.

#import numpy.oldnumeric as N
import numpy.oldnumeric as Numeric
N = Numeric
import UserList
from mglutil.math.rmsd import RMSDCalculator
from Pmv.setangleCommands import SetRelativeTorsion, SetQuaternion

class Clust(UserList.UserList):
    def __init__(self, mol, seed, coords):
        UserList.UserList.__init__(self)
        self.seed = seed
        self.mol = mol
        self.seed_coords = coords 
        self.ruler = RMSDCalculator(coords)
        self.min_energy = self.max_energy = mol.autodock_states[seed].e_binding ##### get seed energy HERE
        self.append(seed)

    def append(self, item):
        # override append to maintain min, max
        if self.mol.autodock_states[item].e_binding < self.min_energy:
            self.min_energy = self.mol.autodock_states[item].e_binding
        if self.mol.autodock_states[item].e_binding > self.max_energy:
            self.max_energy = self.mol.autodock_states[item].e_binding
        self.data.append(item)
   
class Cluster_AD:
    """
    """
    def __init__(self, mol):
        self.mol = mol
        energies = []
        for state in mol.autodock_states:
            energies.append(state.e_binding)
        self.energies = N.array(energies)
        self.argsort = N.argsort(self.energies)
        self.test_coords = None
        self.rmsdarray=N.zeros([len(self.mol.autodock_states),
                                len(self.mol.autodock_states)])-1.
        
    def set_reference(self, reference=None):
        """reference should be an index into the data array.
        If it's a problem to have the reference be part of the data
        array, reference could be a separate instance of data.
        """
        if not reference:
            self.reference = self.argsort[0]  # reference lowest energy conformation
        else:
            self.reference = reference
        self.get_coords(self.mol, self.reference)
       
    def set_tolerance(self, tolerance):
        self.tolerance = tolerance

    def get_dist(self, clust, probeidx):
        if self.rmsdarray[clust.seed, probeidx] == -1:  
            self.get_coords(self.mol, probeidx) # sets self.test_coords
            clust.ruler.computeRMSD(self.test_coords)
            self.rmsdarray[clust.seed, probeidx] = clust.ruler.rmsd
        return self.rmsdarray[clust.seed, probeidx]

    def get_coords(self, mol, val):
        SetTors = SetRelativeTorsion()
        SetQuat = SetQuaternion()
        mol.torTree.fillTorAtoms(mol)
        #BEFORE CALCULATING NEW COORDS: RESET TO initial ones
        if len(mol.allAtoms[0]._coords)==1:
            mol.allAtoms.addConformation(mol.allAtoms.coords)
            mol.allAtoms.setConformation(1)
        else:
            for i in range(len(mol.allAtoms)):
                at = mol.allAtoms[i]
                at._coords[at.conformation] = at._coords[0]
        #do something special for 0
        if val!=-1:
            state = mol.autodock_states[val]
            tList = mol.torTree.tList
            angList = state.torsions
            allAtoms = mol.allAtoms
            for k in range(state.ntorsions):
                thisT = tList[k]
                at1 = allAtoms[thisT.atm1_ix]
                at2 = allAtoms[thisT.atm2_ix]
                SetTors.doit(at1, at2, angList[k], 
                             thisT.mov_atoms)
            SetQuat.doit(allAtoms, state.quaternion, 
                         state.origin, state.translation)
        self.test_coords = mol.allAtoms.coords

    def make_clusters(self, tolerance=None, ref=None):
        if tolerance:
            self.set_tolerance(tolerance)
        self.set_reference(ref) # set self.test_coords
        clusters = []
        clusters.append(Clust(self.mol, self.reference, self.test_coords))
        self.test_coords = None # reset self.test_coords
        # Go through all conformations
        for i in range(len(self.argsort)):
            probe_conf_ix = self.argsort[i]
            if probe_conf_ix == self.reference: continue
            #Assume this is a new conformation until proven otherwise...
            new_conf = 1
            for c in clusters:
                #get dist to lowest energy in cluster (seed)
                dist = self.get_dist(c, probe_conf_ix)
                #Check rms; if greater than tolerance
                if dist > self.tolerance:
                    continue
                else: # add this conformation to the current cluster
                    c.append(probe_conf_ix)
                    new_conf = 0
                    break
            if new_conf: # start a new cluster
                if not self.test_coords:
                    self.get_coords(self.mol, probe_conf_ix)
                clusters.append(Clust(self.mol, probe_conf_ix, self.test_coords))
            self.test_coords = None
        self.clusters = clusters
        return self.clusters













