############################################################################
#
# Author: Michel F. SANNER
# Reimplemented from Babel v1.6 from Pat Walters and Math Stahl
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################
#
# $Id: cycle.py,v 1.8 2003/09/04 23:53:56 lindy Exp $
#

"""
This file implements the RingFinder class that can be used to identify
rings in molecules. When rings are nested the smalest rings are reported.
The algorithms in here are different from thoses in Babel1.6. They might
fail to report all rings for rings in which any atom of a ring belongs
also to more than 1 ring.
findRings2 is more robust and should reports all the smallest cycles. It is
more expansive too. Should be used for small molecules

example:
    
    >>> r = RingFinder()
    >>> r.findRings2(atoms, bonds)
    >>> r.printRings()

      atoms has to be a list of Atom objects
      Atom:
          a.coords : 3-sequence of floats
          a.bonds : list of Bond objects
      Bond:
          b.atom1 : instance of Atom
          b.atom2 : instance of Atom

      after completion the RingFinder object has the following members:
          ringCount: number of rings
          rings : a list of rings. A ring is a dictionary with 2 keys
                  'atoms' and 'bonds'. The atoms in rings['atoms'] are
                  ordered along the cycle.
          allRingAtoms: list of atoms that are in rings(atoms may appear twice)
          allRingBonds: list of bonds that are in rings(bonds may appear twice)

      In addition:
          atoms involved in rings have a member 'rings' that is a list of
          rings they belong to (0-based list of integers)

Michel Sanner April 2000
"""




class RingFinder:
    """ """
    def __init__(self):
        """ """
        self.num = 0
        self.rings = []


    def tagOneAtomBond(self, atom, inbond):
        # tag all bond that cannot be in a ring by adding them as keys to
        # self.bondInCycle[bond], value is set to 2 (no particular meaning)

        # since findRings2 will oly consider atoms for which at least one
        # bond is not yet in a cycle, pretending that the only bond to an
        # atom is already in a cycle will pevent it from being considered
        if len(atom.bonds)==1:
            bond = atom.bonds[0]
            self.bondInCycle[bond] = 2
            #print 'REMOVE ',atom.name,bond.atom1.name,'-',bond.atom2.name
            if bond!=inbond:
                atom2 = bond.atom1
                if atom2==atom:
                    atom2=bond.atom2
                self.tagOneAtomBond(atom2, bond)
        else:
            num = 0
            for b in atom.bonds:
                if not self.bondInCycle.has_key(b):
                    num = num + 1
                    bond = b
            if num == 1:
                #print 'REMOVE ',atom.name,bond.atom1.name,'-',bond.atom2.name
                self.bondInCycle[bond] = 2
                atom2 = bond.atom1
                if atom2==atom:
                    atom2=bond.atom2
                self.tagOneAtomBond(atom2, bond)

                
    def findRings2(self, atoms, bonds, maxSize=20):
        # EXPANSIVE but robust
        # for each atom for which at least one bond is not yet in a cycle
        # find the smallest cycle. Then check if the cycle is already known
        # if maxSize is specified, only ring of size maxSize will be found
        
        self.bondInCycle = {} # key is bond, is created when bond has been
                              # traversed (persistant across search over atoms
        self._rings = []      # dict. with atoms as keys per cycle
        self.allRingAtoms = {} # dict with keys being all atoms in cycles
        self.allRingBonds = {} # dict with keys being all bonds in cycles
        self.ringCount = 0     # number of cycles
        maxLevel = maxSize/2 + 1

        # tag all bonds that cannot be in cycles, i.e. they connect an atom
        # having only 1 neighbor (recursively)
        #for a in atoms:
        #    self.tagOneAtomBond(a, None)
            
        # loop over all atoms
        for a in atoms:
            # being in a cycle requires at least 2 neighbors
            if len(a.bonds)==1:
                continue
            # count number of bond already in cycles or connecting
            # to an leaf atom (i.e. atom with only 1 neighbor)
            num = 0
            for b in a.bonds:
                atom2 = b.atom1
                if atom2==a: atom2=b.atom2
                if len(atom2.bonds)==1 or self.bondInCycle.has_key(b):
                    num = num + 1

            # if at least on bond not in a cycle find smallest cycle for a
            if num < len(a.bonds):
                #print 'Find smallest Ring for', a.name, num, len(a.bonds)
                # ra is a list of atoms in cycle
                # rb is a list of bonds in cycle
                #print 'Find smallest cycle for', a.name
                ra, rb = self.findSmallestRing(a, maxLevel)
                #print 'Found cycle:', ra
                # check is that cycle has been found before
                if len(ra):
                    same = 0
                    for r in self._rings:
                        if len(ra)==len(r):
                            same = 1
                            for a in ra:
                                if not r.has_key(a):
                                    same = 0
                                    break
                            if same==1:
                                break

                    # if it hasn't be found before, add it
                    if not same:
                        # atoms in ra have to be sorted along the cycle
                        ras = [ra[0]]
                        # Take all odds first going trough list toward end
                        for i in range(1, len(ra), 2):
                            ras.append(ra[i])
                        # take all even going through list backwards and
                        # starting a last of cycles woth even number of atoms
                        # and one before last of rinfs with odd number of at.
                        l = len(ra)
                        if (l/2*2) != l:
                            end = l-1 # odd -> skip last
                        else:
                            end = l-2 # end, use last
                        for i in range(end, 1, -2):
                            ras.append(ra[i])
                        ra = ras
                        
                        # build the dict of atoms for that cycle
                        d = {}
                        ringnum = len(self._rings)
                        for a in ra:
                            d[a] = 1
                            if not hasattr(a, 'rings'):
                                a.rings = [ringnum]
                            else:
                                a.rings.append(ringnum)
                        # add it to _rings
                        self._rings.append(d)
                        # add the cycle to self.rings
                        self.rings.append( {'atoms':ra, 'bonds':rb } )
                        # update dict of all atoms in cycles
                        self.allRingAtoms.update(d)
                        # update dict of all bonds in cycles
                        for b in rb:
                            self.allRingBonds[b] = 1
                        # increment cycl counter
                        self.ringCount = self.ringCount + 1

##                          print 'RING ======================'
##                          for a in ra:
##                              print a.name
##                          for b in rb:
##                              print b.atom1.name,'-',b.atom2.name
##                          print 'END RING ======================'
        self.allRingAtoms = self.allRingAtoms.keys()
        self.allRingBonds = self.allRingBonds.keys()


    def findSmallestRing(self, root, maxLevel):
        # each new generation adds a level in the width first traversal
        level = 0
        # each level has a stack of atoms
        stacks = []
        # each level has a dict of atom:bond telling through which bond we
        # came to this atom
        bndDicts = []
        # dict of atoms used to check if an atom has been seen before
        atinstack = {}
        # dict of bonds used to check if a bond has been seen before
        bondseen = {}

        # bond through which we came (None for root)
        bond = None
        atinstack[root] = 1
        # list of atoms a this level
        levelstack = []
        stacks.append(levelstack)
        # dict of atom:bond for this level
        levelDict = {}
        bndDicts.append(levelDict)

        # add children of root to this level's stack and dict
        for b in root.bonds:
            bondseen[bond] = 1
            atom2 = b.atom1
            if atom2==root:
                atom2 = b.atom2
            if len(atom2.bonds)>1:
                #print 'adding to stack 0', atom2.name
                levelDict[atom2] = b
                levelstack.append(atom2)
                atinstack[atom2] = 1

        maxLen = len(levelstack) # number of bonds with level 0
        stackPtr = 0

        # width first traversal, i.e. loop over stack adding levels
        while stackPtr < maxLen:
            # add a level
            levelstack = []
            stacks.append(levelstack)
            levelDict = {}
            bndDicts.append(levelDict)
            #print "Looping over LEVEL: ", level
            # loop over atoms at this level
            for levelroot in stacks[level]:
##                  if len(levelroot.bonds)==1:
##                      continue
                # find bond through with we came to levelroot atom
                bond = bndDicts[level][levelroot]
                # if already seen get next level root

##                  if bondseen.has_key(bond):
##                      continue

##                  if len(bond.atom1.bonds)==1 or len(bond.atom2.bonds)==1:
##                      continue # this bond cannot be in a cycle
##                  if self.bondInCycle.has_key(bond) and \
##                     self.bondInCycle[bond]==2:
##                      continue
                # else make this bond as seen
                bondseen[bond] = 1
                #print 'level ROOT', levelroot.name, bond, atinstack.has_key(levelroot)
                # add children of levelroot to the current level
                for b in levelroot.bonds:
                    # except for the parent of levelroot
                    if b is bond:
                        continue
                    atom2 = b.atom1
                    if atom2==levelroot:
                        atom2 = b.atom2
                    if len(atom2.bonds)==1:
                        continue
                    
                    #print levelroot.name, atom2.name
                    # if the child of levelroot is instack we found a cycle
                    # so we start back tracking from both sides through the
                    # levels until we reach a common atom. If that atom is
                    # root we have the smallest cycle for root, else we
                    # continue

                    # 2 cases are possible: either we have a even or an odd
                    # number of atoms in the cycle. For even numbers atom2
                    # is found in atinstack which is associated with this
                    # level. else atom2 should be found in
                    # bndDicts[level]

                    if atinstack.has_key(atom2):
                        #print 'CYCLE', atom2.name
                        # even number of atoms
                        #print 'instack'
                        #for a in atinstack:
                        #    print a.name
                        if levelDict.has_key(atom2):
                            #print 'EVEN ******************'
                            # cycle with even number of atoms
                            # b1 is the bond from which we arrived at atom2
                            # previousely at this level
                            b1 = levelDict[atom2]
                            at1 = levelroot
                            at2 = b1.atom1
                            if at2==atom2: at2 = b1.atom2

                            #print 'ringAtoms1', atom2.name, at1.name, at2.name
                            ringAtoms = [atom2, at1, at2]
                            ringBonds = [b, b1]
                            # backtrack through level
                        elif bndDicts[level].has_key(atom2):
                            # odd number of atoms in cycle
                            # b1 is the bond from which we arrived at atom2
                            # previousely at this level
                            #print 'ODD ******************'
                            at1 = levelroot
                            at2 = atom2
                            #print 'ringAtoms2', at1.name, at2.name

                            ringAtoms = [at1, at2]
                            ringBonds = [b]
                        else:
                            continue

                        # backtrack
                        for i in range(level,-1,-1):
                            #print 'level:', i
                            b1 = bndDicts[i][at1]
                            other1 = b1.atom1
                            if other1==at1: other1 = b1.atom2
                            #print other1.name
                            b2 = bndDicts[i][at2]
                            other2 = b2.atom1
                            if other2==at2: other2 = b2.atom2
                            #print other2.name
                            at1 = other1
                            at2 = other2
                            ringBonds.append(b1)
                            ringBonds.append(b2)
                            ringAtoms.append(other1)
                            if other1!=other2:
                                ringAtoms.append(other2)
                            else:
                                break
                        if other1==root or other2==root:
                            for b in ringBonds:
                                self.bondInCycle[b] = 1
                            return (ringAtoms, ringBonds)

                    else:
                        #print 'adding to stack', atom2.name, len(levelstack)
                        levelstack.append(atom2)
                        levelDict[atom2] = b
                        atinstack[atom2] = 1

            level = level + 1
            if level==maxLevel:
                return [],[]
            maxLen = maxLen + len(levelstack)
            stackPtr = stackPtr + 1
        return [],[]
    

    def backtrack(self, atom1, atom2, b):
        """go up ancestor tree until first common parent is found"""
        ringAtoms = [ atom2 ]
        ringbonds = [ b ]
        b._ring_seen = 1
        while atom1!=atom2:
            ringAtoms.append( atom1 )
            ringbonds.append( atom1._ring_ancestor_bond )
            atom1 = atom1._ring_ancestor
        return { 'atoms':ringAtoms, 'bonds':ringbonds }


    def tag_neighbors(self, atom1, bond):
        """ """
        atom2 = bond.atom1
        if atom2==atom1: atom2 = bond.atom2
        #print 'atom2', atom2.name, bond
        if hasattr(atom2,'_ring_ancestor'):
            self._rings.append( self.backtrack(atom1, atom2, bond) )
            #print '\n********************* ring found'
            #for a in self._rings[-1]['atoms']:
                #print a.name
            return

        atom2._ring_ancestor = atom1
        atom2._ring_ancestor_bond = bond
        bond._ring_seen = 1
        for b in atom2.bonds:
            #FIXME: we are supposed to skip di-sulfite bridges here
            if hasattr(b, '_ring_seen') and b._ring_seen:
                continue
            self.tag_neighbors( atom2, b )
        #print 'step back from', atom2.name, bond
        
                
    def findRings(self, atoms, bonds):
        """method to find cycles in molecules, first we simply tag all atoms
        and bonds in rings using a depth first traversal then we call
        checkRings for identifying the smallest cycles when fused rings are
        present
        """

        i = 0
        for b in bonds:
            b._ring_seen = 0
            
        self._rings = []

        # first try with leaf atoms (i.e only 1 neighbor)
        done = 0
        for a in atoms:
            if len(a.bonds)==1:
                a._ring_ancestor = None
                a._ring_ancestor_bond = None
                #print 'start at ',a, a.bonds[0]
                self.tag_neighbors(a, a.bonds[0])
                done = 1
                break

        #print 'AAAA', done
        
        # for molecules with all atoms in cycles we did not find a leaf
        if not done:
            for a in atoms:
                if not hasattr(a, '_ring_ancestor'):
                    l = len(a.bonds)
                    if l==2:
                        atom_before = a.bonds[0].atom1
                        if atom_before==a: atom_before=a.bonds[0].atom2
                        a._ring_ancestor = atom_before
                        a._ring_ancestor_bond = a.bonds[0]
                        self.tag_neighbors(a, a.bonds[1])

                        if a.bonds[0]._ring_seen:
                            continue
                        atom_before = a.bonds[1].atom1
                        if atom_before==a: atom_before=a.bonds[1].atom2
                        a._ring_ancestor = atom_before
                        a._ring_ancestor_bond = a.bonds[0]
                
##          for a in atoms:
##              if not hasattr(a, '_ring_ancestor'):
##                  l = len(a.bonds)
##                  if l==0 or l > 2: continue
##                  if l==2:
##                      atom_before = a.bonds[0].atom1
##                      if atom_before==a: atom_before=a.bonds[0].atom2
##                      a._ring_ancestor = atom_before
##                      a._ring_ancestor_bond = a.bonds[0]
##                      self.tag_neighbors(a, a.bonds[1])

##                      if a.bonds[0]._ring_seen:
##                          continue
##                      atom_before = a.bonds[1].atom1
##                      if atom_before==a: atom_before=a.bonds[1].atom2
##                      a._ring_ancestor = atom_before
##                      a._ring_ancestor_bond = a.bonds[0]
##                  else:
##                      a._ring_ancestor = None
##                      a._ring_ancestor_bond = None
##                      self.tag_neighbors(a, a.bonds[0])

        for a in atoms:
            if hasattr(a,'_ring_ancestor'):
                delattr(a, '_ring_ancestor')
            if hasattr(a,'_ring_ancestor_bond'):
                delattr(a, '_ring_ancestor_bond')
        for b in bonds:
            delattr(b, '_ring_seen')

        self.checkRings()
        delattr(self, '_rings')

        
    def smallestCycle(self, atom):
        """
        find smalest cycle containing starting at atom and traversing
        tree in breadth first order. When an atom with _ancestor is found
        we backtrack and both sides to build listst of atoms and bonds.
        It is possible that the first cycle does not contain the initial
        atom.
        """
        if self._result:
            return self._result
        
        for b in atom.bonds:
            if self._result:
                return self._result
            if not hasattr(b, '_ring'): continue
            if b._seen: continue
            atom2 = b.atom1
            if atom2==atom: atom2 = b.atom2
            if atom2._ancestor:
                l1 = [atom2]
                atom2.rings.append( self.ringCount )
                l2 = [atom]
                atom.rings.append( self.ringCount )
                b1 = [b]
                b2 = []
                a2 = atom2
                a1 = atom
                done = 0
                while not done:
                    b1.append(a2._ancestor_bond)
                    a = a2._ancestor
                    if a != l2[-1]:
                        l1.append(a)
                        a.rings.append( self.ringCount )
                        a2 = a
                    else:
                        break

                    b2.append(a1._ancestor_bond)
                    a = a1._ancestor
                    if a != l1[-1]:
                        l2.append(a)
                        a.rings.append( self.ringCount )
                        a1 = a
                    else:
                        break
                    
                l2.reverse()
                b2.reverse()
                self._result = (l2+l1, b2+b1)
                return
            else:
                atom2._ancestor = atom
                atom2._ancestor_bond = b
                self.stack.append(atom2)
                b._seen = 1
           
        if len(self.stack):
            a = self.stack[0]
            self.stack.remove(a)
            self.smallestCycle(a)

        return self._result


    def checkRings(self):
        """
        this functions uses the rings found by findRings to identify smalest
        cycles in structure.

        After this method was called this object as the following new members:
          rings: a list of dictionnaries {'atoms':list of atoms, 'bonds': list}
          allRingAtoms: list of all atoms in rings
          allRingBonds: list of all bonds in rings
          ringCount: number of rings
        """
        self.rings = []
        self.allRingAtoms = []
        self.allRingBonds = []
        self.ringCount = 0
        for ring in self._rings:
            self.allRingAtoms = self.allRingAtoms + ring['atoms']
            self.allRingBonds = self.allRingBonds + ring['bonds']

        for a in self.allRingAtoms: a.rings = []
        for b in self.allRingBonds: b._ring = 1

        # breadth first walking over atoms and bonds in rings
        # tag atoms in shortest cycle
        for a in self.allRingAtoms:
            if len(a.rings)>0: continue

            self.stack = []
            for b in self.allRingBonds: b._seen = 0
            for at in self.allRingAtoms: at._ancestor = None
            a._ancestor_bond = None
            self._result = None

            ratoms, rbonds = self.smallestCycle(a)

            if a in ratoms:
                self.ringCount = self.ringCount+1
                self.rings.append( {'atoms':ratoms, 'bonds':rbonds} )
            else:
                for at in ratoms:
                    at.rings = at.rings[:-1]

        # clean up
        for a in self.allRingAtoms:
            if hasattr(a, '_ancestor'): delattr(a, '_ancestor')
            if hasattr(a, '_ancestor_bond'): delattr(a, '_ancestor_bond')

        for b in self.allRingBonds:
            if hasattr(b, '_ring'): delattr(b, '_ring')
            if hasattr(b, '_seen'): delattr(b, '_seen')

        if hasattr(self, '_result'): delattr(self, '_result')

        
    def printRings(self):
        """ """
        if not hasattr(self, 'rings'):
            return
        i = 0
        for r in self.rings:
            print 'RING ',i
            for j in range(len(r['atoms'])):
                a = r['atoms'][j]
                b = r['bonds'][j]
                print '%10s %4d %s'%(a.name, a.number, repr(b))
            i = i + 1

            
if __name__ == '__main__':
    import pdb, sys
    from MolKit.pdbParser import NewPdbParser
    parser = NewPdbParser("/tsri/pdb/struct/%s.pdb"%sys.argv[1])
    mols = parser.parse()
    mol = mols[0]
    mol.buildBondsByDistance()
    allAtoms = mol.chains.residues.atoms

    print 'Looking for rings'
    r = RingFinder()
    bonds = (allAtoms.bonds)[0]
    r.findRings(allAtoms, bonds)
    r.printRings()


    from MolKit.pdbParser import NewPdbqParser
    parser = NewPdbqParser("./txp.pdbq")
    mols = parser.parse()
    mol = mols[0]
    mol.buildBondsByDistance()
    allAtoms = mol.chains.residues.atoms


    print 'Looking for rings ...'
    r = RingFinder()
    print "Done"
    bonds = (allAtoms.bonds)[0]
    r.findRings(allAtoms, bonds)

    r.printRings()
