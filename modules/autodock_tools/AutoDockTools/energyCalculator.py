## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#############################################################################
#
# Author: Ruth HUEY, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################


# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/energyCalculator.py,v 1.5 2007/07/24 17:30:44 vareille Exp $
#
# $Id: energyCalculator.py,v 1.5 2007/07/24 17:30:44 vareille Exp $
#
#

"""
This class calculates the AutoDock internal energy for a set of atoms.
"""

import numpy.oldnumeric as Numeric
from math import sqrt
from AutoDockTools.energyConstants import Rij, epsij
from MolKit.molecule import AtomSet


        
class EnergyCalculator:
    """This class knows how to calculate the energy of its set of atoms
    the Autodock way.  

    self.data: the list of atoms for the calculation
    self.reset: if 1, recalculates all the energies each time
    self.lowLim: number to be flag that calculation hasn't been done
    self.energy_type: internal, some day we may have others...
    
    """
    def __init__(self, data, reset=0, lowLim=-999999.0, energy_type='internal'):
        """data is a set of atoms
        """
        self.data = data
        self.data_ids = []
        for a in data:
            self.data_ids.append(id(a))
        #have to make sure that atom.number is in order
        self.data.sort()
        # save the pair-wise energies for reuse by get_energy
        # code: -1 if not set, 0 if within 1-4 distance, otherwise pairwise
        # energy
        self.energy_matrix = Numeric.zeros([len(data), len(data)]) + lowLim
        self.lowLim = lowLim
        self.setRotatable()
        self.setNeighbors()
        self.weedbonds()
        self.reset = reset
        self.NBCUTOFF = 8.0


    def setRotatable(self):
        #find the bonds which are in the torsion Map and label rotatable
        mols = self.data.top.uniq()
        for m in mols:
            if not hasattr(m, 'torTree'):
                continue
            l = []
            m.allAtoms.bonds[0].rotatable = 0
            for n in m.torTree.torsionMap:
                ats = m.allAtoms.get(lambda x, n=n: x.number-1 in n.bond)
                #ats = AtomSet([n.a,n.b])
                b = ats.bonds[0][0]
                b.rotatable = 1


    def weedbonds(self):
        mols = self.data.top.uniq()
        for m in mols:
            if not hasattr(m, 'torTree'):
                continue
            l = []
            atL = m.torTree.rootNode.atomList
            #print 'setting 0 for rootNode atomList:', atL
            for i in atL:
                for j in atL:
                    self.energy_matrix[i][j] = 0
            for n in m.torTree.torsionMap:
                atL = n.atomList
                atL.extend(n.bond)
                #print 'setting 0 for atomList:', atL
                for i in atL:
                    for j in atL:
                        self.energy_matrix[i][j] = 0


    def setNeighbors(self):
        """
        detect 1-2, 1-3 and 1-4 neighbors and set their entries in
        self.energy_matrix to 0
        FIX THIS: if bond is rotatable, don't reset the 1-4 pair
        """
        
        self.nb = {}
        for a1 in self.data:
            ind1 = self.data.index(a1)
            nbd =  {}
            nbd[ind1] = 0
            #nblist = self.nb[ind1] = [ind1]
            self.energy_matrix[ind1, ind1] = 0
            #set 1-2 interactions to 0
            for b in a1.bonds:
                a2 = b.atom1
                if id(a2)==id(a1): a2 = b.atom2
                if a2.number > a1.number:
                    #then do something here
                    ind2 = self.data.index(a2)
                    nbd[ind2] = 0
                    self.energy_matrix[ind1, ind2] = 0
                    self.energy_matrix[ind2, ind1] = 0
                #set 1-3 interactions to 0
                for b2 in a2.bonds:
                    a3 = b2.atom1
                    if id(a3)==id(a2): a3 = b2.atom2
                    if id(a3)==id(a1): continue
                    if a3.number > a1.number:
                        #then do something here
                        ind3 = self.data.index(a3)
                        nbd[ind3] = 0
                        self.energy_matrix[ind1, ind3] = 0
                        self.energy_matrix[ind3, ind1] = 0
                    #set 1-4 interactions to 0
                    #if rotatable, skip 1-4s for b2
                    if b2.rotatable:
                        #print 'skipping ', b2.atom1.name, '-', b2.atom2.name
                        continue
                    for b3 in a3.bonds:
                        #if hasattr(b3, 'activeTors') and b3.activeTors: continue
                        a4 = b3.atom1
                        if id(a4)==id(a3): a4 = b3.atom2
                        if id(a4)==id(a2): continue
                        if id(a4)==id(a1): continue
                        if a4.number > a1.number:
                            #then do something here
                            ind4 = self.data.index(a4)
                            nbd[ind4] = 0
                            self.energy_matrix[ind1, ind4] = 0
                            self.energy_matrix[ind4, ind1] = 0
            self.nb[ind1] = nbd.keys()


    def dist(self, a, b):
        """return distance between two atoms, a and b.
        """
        d = Numeric.array(b.coords) - Numeric.array(a.coords)
        return sqrt(Numeric.sum(d*d))

                        

    #def _get_energy_custom(self, a, b):
    #    """return energy between two atoms, a and b.
    #    """
    #    ax = self.data.index(a)
    #    bx = self.data.index(b)
#
#        if self.energy_matrix[ax][bx] >= 0.0:
#            # return previously saved energy
#            return self.energy_matrix[ax][bx]
#        else:
#            # compute, save, and return energy
#            #FIX THIS
#            dist = a.getRMSD_custom(b.getCoords())
#            self.energy_matrix[ax][bx] = self.energy_matrix[bx][ax] = dist
#            return dist


    def get_energy(self, a, b):
        """return energy between two atoms, a and b.
        """
        ax = self.data.index(a)
        bx = self.data.index(b)

        e = self.energy_matrix[ax][bx]
        if e > self.lowLim and not self.reset:
            # return previously saved energy
            return e
        elif e==0:
            # if getNeighbors set this to 0
            return e
        else:
            # compute, save, and return  energy
            energy = self.nbenergy(a,b)
            self.energy_matrix[ax][bx] = self.energy_matrix[bx][ax] = energy
            return energy


    def nbenergy(self, a, b):
        r = self.dist(a,b)
        if r>=self.NBCUTOFF:
            #print 'skipping ', a.number, '-', b.number
            return 0.0
        key = 'lj'+ a.element + b.element
        r_eq = Rij[key]
        ep = epsij[key]
        n = 12
        m = 6
        energy = (m/(n-m))*ep*r_eq**n/(r**n) - (n/(n-m))*ep*r_eq**m/r**m
        #energy = m*epsilon*r_eq**n/(r**n)/(n-m) - n*epsilon*r_eq**m/r**m/(n-m)
        return  energy


    def get_internal_energy(self):
        intEnergy = 0
        datalen = len(self.data)
        for i in range(datalen):
            a = self.data[i]
            for j in range(i+1, datalen):
                b = self.data[j]
                eij = self.get_energy(a,b)
                intEnergy = intEnergy + eij
        return intEnergy
                
    
    def get_subset_internal_energy(self, ats):
        subset_intEnergy = 0
        for a1 in ats:
            if id(a1) not in self.data_ids:
                print 'warning', a1.name, ' not in self.data'
            for a2 in ats:
                if id(a2) not in self.data_ids:
                    print 'warning', a2.name, ' not in self.data'
                eij = self.get_energy(a1, a2)
                subset_intEnergy = subset_intEnergy + eij
        return subset_intEnergy


    def reset_matrix(self):
        self.energy_matrix = Numeric.zeros([len(self.data), len(self.data)]) + self.lowLim
        self.weedbonds()
        self.setNeighbors()



    def print_nb_matrix(self):
        n = len(self.data)
        s = 'Atom: ID:'
        for i in range(1,n+1):
            s = s + '%2d'%i
        print s
        s = '_'*n
        print s
        em = self.energy_matrix
        for i in range(n):
            at = self.data[i]
            s = '  %2s     %2d'%(at.name,at.number)
            for j in range(n):
                if em[i][j]!=0:
                    s = s + '|x'
                else:
                    s = s + '|_'
            print s


    def print_half_nb_matrix(self):
        n = len(self.data)
        em = self.energy_matrix
        for i in range(n):
            s = '' + i*'_'
            for j in range(i, n):
                if em[i][j]!=0:
                    s = s + 'x'
                else:
                    s = s + '_'
            print s

