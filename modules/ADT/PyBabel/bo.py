#############################################################################
#
# Author: Michel F. SANNER
# Reimplemented from Babel v1.6 from Pat Walters and Math Stahl
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################
#
# $Header: /opt/cvs/python/packages/share1.5/PyBabel/bo.py,v 1.4 2002/10/22 23:11:38 sanner Exp $
#
# $Id: bo.py,v 1.4 2002/10/22 23:11:38 sanner Exp $
#

"""
This file implements the BondOrder class that can be used to compute
bond order.

Before a BondOrder object can be used, atoms must have been assigned
a type see (AtomHybridization in types.py).

Bond order can be calculated using 2 different methods depending on whether
rings have been identified previously or not. Babel decides to use the first
method for molecules with more than 200 atoms and the second one else.
    
example:
      >>> from PyBabel.atomTypes import AtomHybridization
      >>> from PyBabel.cycle import RingFinder
      >>> atype = AtomHybridization()
      >>> atype.assignHybridization(atoms)
      >>> bo = BondOrder()
      >>> bo.assignBondOrder( atoms, bonds )

      or

      >>> atype = AtomHybridization()
      >>> atype.assignHybridization(atoms)
      >>> rings = RingFinder()
      >>> rings.findRings(allAtoms, bonds)
      >>> bo = BondOrder()
      >>> bo.assignBondOrder( atoms, bonds, rings )

      atoms has to be a list of atom objects
      Atom:
          a.coords : 3-sequence of floats
          a.bonds : list of Bond objects
          babel_type: string
          babel_atomic_number: int

      Bond:
          b.atom1 : instance of Atom
          b.atom2 : instance of Atom

      after completion each bond has a 'bondOrder' attribute (integer)

reimplmentation of Babel1.6 in Python by Michel Sanner April 2000
Original code by W. Patrick Walters and Matthew T. Stahl 
"""



import string, math

from babelAtomTypes import babel_types
from babelElements import babel_elements
from util import *
from atomTypes import TypeConverter


SINGLE_DOUBLE_CUTOFF = 0.95
#SINGLE_DOUBLE_CUTOFF = 0.955
DOUBLE_TRIPLE_CUTOFF = 0.81

class BondOrder:
    """ """

    def assignBondOrder(self, atoms, bonds, rings=None):
        """ """
        if not rings:
            self.assign_bond_order1(atoms, bonds)
        else:
            self.assign_bond_order2(atoms, bonds, rings)


    def assign_bond_order1(self, atoms, bonds):
        """ """

        hyb_val = [0,3,2,1]
        converter = TypeConverter("HYB")
        for a in atoms:
            hyb = converter.convert(a.babel_type, 'dummy')
            a._redo = hyb_val[int(hyb)]
            #print a.full_name(), a.babel_type, hyb, a._redo

        for b in bonds:
            # initialize bondOrder attribute
            if b.bondOrder is None:
                b.bondOrder = 1
            sum_code = b.atom1._redo + b.atom2._redo
            #print b, sum_code
            if sum_code == 6:
                b.bondOrder = 3
            elif sum_code == 4:
                b.bondOrder = 2
            else:
                b.bondOrder = 1
            
            if self.is_carboxyl(b):
                b.bondOrder = 2

            if b.bondOrder < 1 or b.bondOrder > 3:
                print "Bond %s is wierd - Bond order is %d\n" % \
                      (b, b.bondOrder)

        self.check_for_conjugation(atoms)

        # cleanup
        for a in atoms:
            delattr(a, '_redo')


    def is_carboxyl(self, bond):
        """ """
        c_end = 0
        o_end = 0
        check = 0
        
        if bond.atom1.babel_type == "Cac" and bond.atom2.babel_type == 'O':
            c_end = bond.atom1
            o_end = bond.atom2
            check = 1

        if bond.atom2.babel_type == "Cac" and bond.atom1.babel_type == 'O':
            check = 1
            c_end = bond.atom2
            o_end = bond.atom1

        if check and len(o_end.bonds) == 1: return 1
        else: return 0


    def check_for_conjugation(self, atoms):
        """ """

        for a in atoms:
            #if a.full_name()=='1crn: :ASN14:CG': raise
            for b1 in a.bonds:
                if b1.bondOrder<=1: continue
                for b2 in a.bonds:
                    if b1==b2: continue
                    if b2.bondOrder<=1: continue
                    if len(b2.atom1.bonds) > 1 and len(b2.atom2.bonds) > 1:
                        b2.bondOrder = 1


    def check_for_carbonyl(self, atom):
        """ """
        for b in atom.bonds:
            bonded_atom = b.atom1
            if bonded_atom==atom:
                bonded_atom = b.atom2
            if bonded_atom.babel_type=="O2" or bonded_atom.babel_type=="S2":
                return 3
            return 2


    def assign_bond_order2(self, atoms, bonds, rings):
        """ """

        self.rings = rings

        for a in atoms:
            if hasattr(a, 'rings'):
                a._redo = 1
            else:
                a._redo = 0
            a._dot = 0
            a._dbatom = 0
            
        self.assign_hybrid_radii(atoms)
        self.estimate_bond_order2(bonds)

        for ring in self.rings.rings:
            if len(ring['atoms'])==5:
                self.process_5_ring(ring)

        for b in bonds:
            # initialize bondOrder attribute
            if b.bondOrder is None:
                b.bondOrder = 1
            b._dbbond = 0
            if b.bondOrder == 2:
                if len(b.atom1.bonds) > 1 and b.atom1.babel_type[0] == 'O':
                    b.bondOrder = 1
                elif len(b.atom2.bonds) > 1 and b.atom2.babel_type[0] == 'O':
                    b.bondOrder = 1

        for b in bonds:
            if b.bondOrder == 2:
                if len(b.atom1.bonds) == 1 and b.atom1.babel_type[0] == 'N':
                    b.bondOrder = 1
                elif len(b.atom2.bonds) == 1 and b.atom2.babel_type[0] == 'N':
                    b.bondOrder = 1

##          for b in bonds:
##              if b.bondOrder > 1:
##                  print "%3d %3d"%(atoms.index(b.atom1)+1, atoms.index(b.atom2)+1)
        for b in bonds:
            if b.bondOrder > 1:
                a1 = b.atom1
                if a1._redo and self.check_for_carbonyl(a1) != 3:
                    a1._dot = 1
                a2 = b.atom2
                if a2._redo and self.check_for_carbonyl(a2) != 3:
                    a2._dot = 1
                if len(a1.bonds) == 1 or len(a2.bonds) == 1:
                    a1._dot = 0
                    a2._dot = 0

##          for a in atoms:
##              if a._dot==1:
##                  print atoms.index(a)+1, a._dot

        for a in atoms:
            if a.babel_type=="Npl" and len(a.bonds) == 3:
                a._dot = 0

##          for a in atoms:
##              if a._dot==1:
##                  print atoms.index(a)+1, a._dot

        self.atoms = atoms
        self.cycles = 0
        self.bondStack = []
        self.bondChoice = []

        # for PYTHON interpreters after 1.5 the recursion depth is limited
        # this method can exceed the default recursion depth
        import sys
        if float(sys.version[:3]) > 1.5:
            sys.setrecursionlimit(20000)

        self.connect_the_dots(0,0)
        
        for b in self.bondStack:
            if b.bondOrder > 1:
                b._dbbond = 1
                b.atom1._dbatom = 1
                b.atom2._dbatom = 1

        for b in bonds:
            if b.atom1.babel_type=="O2" or b.atom2.babel_type=="O2":
                b._dbbond = 1
                b.atom1._dbatom = 1
                b.atom2._dbatom = 1
            elif b.atom1.babel_type=="O-" and len(b.atom1.bonds) == 1:
                b._dbbond = 1
                b.atom1._dbatom = 1
                b.atom2._dbatom = 1

            elif b.atom2.babel_type=="O-" and len(b.atom2.bonds) == 1:
                b._dbbond = 1
                b.atom1._dbatom = 1
                b.atom2._dbatom = 1

        for a in atoms:
            a._dot = 0

        for b in bonds:
            if b.bondOrder > 1:
                a1 = b.atom1
                a2 = b.atom2
                if a1._dbatom==0 and a2._dbatom==0:
                    a1._dot = 1
                    a2._dot = 1

        self.bondStack = []
        self.bondChoice = []
        self.connect_the_dots(0,0)

        for b in self.bondStack:
            if b.bondOrder > 1:
                b._dbbond = 1

##          for b in bo.bondStack:
##              a1 = b.atom1
##              a2 = b.atom2
##              print a1.number, a2.number, a1._dot, a2._dot, b.bondOrder


        for b in self.bondStack:
            b._dbbond = 1

        for b in bonds:
            if not b._dbbond:
                b.bondOrder = 1

        for a in atoms:
            a._redo = 0
            for b in a.bonds:
                a._redo = a._redo + b.bondOrder

            if (a.babel_atomic_number == 6 or a.babel_atomic_number) == 7 and \
               a._redo > 4:

                for b in a.bonds:
                    conn=None
                    if b.bondOrder == 2:
                        b.bondOrder = 1

        #cleanup
        for a in atoms:
            delattr(a, '_dot')
            delattr(a, '_dbatom')
            delattr(a, '_redo')
            
        for b in bonds:
            delattr(b, '_dbbond')
            
        delattr(self, 'atoms')
        delattr(self, 'cycles')
        delattr(self, 'rings')

##          for b in bonds:
##              n = b.atom1.number
##              print "%4i %4i %2i"%(n, b.atom2.number, b.bondOrder)


    def connect_the_dots(self, atom, start):
        """ """
        a = self.atoms[atom]
        if start == len(a.bonds): return

        #print 'AAA', atom+1

        if a._dot:
            done = 0
            i = start
            for b in a.bonds[start:]:
                con = b.atom1
                if con==a: con=b.atom2
                if con._dot:
                    self.bondStack.append(b)
                    self.bondChoice.append(0)
                    if a==b.atom1:
                        self.bondChoice[-1] = i+1
                    else:
                        self.bondChoice[-1] = -i-1
                    a._dot = a._dot-1
                    con._dot = con._dot-1
                    done = 1
                    break
                i = i + 1
                
            if not done and len(self.bondStack):
                b = self.bondStack[-1]
                if self.bondChoice[-1] > 0: new_atm = b.atom1
                else: new_atm = b.atom2
                choice_bnd = abs(self.bondChoice[-1])
                self.bondChoice = self.bondChoice[:-1]
                self.bondStack = self.bondStack[:-1]
                b.atom1._dot = b.atom1._dot + 1
                b.atom2._dot = b.atom2._dot + 1
                #print 'BBB', self.atoms.index(new_atm)+1
                self.connect_the_dots(self.atoms.index(new_atm), choice_bnd )

        if self.cycles > 10000:
           #print 'EEE'
           return

        if atom+1 == len(self.atoms):
           #print 'DDD'
           return
        else:
            self.cycles = self.cycles+1
            #print 'CCC', atom+2
            self.connect_the_dots(atom+1,0)


    def get_bond(self, bonds, a1, a2):
        """ """
        for b in bonds:
            if (b.atom1==a1 and b.atom2==a2) or (b.atom1==a2 and b.atom2==a1):
                return b


    def process_5_ring(self, ring):
        """ """
        atoms = ring['atoms']
        t1 = torsion_angle(atoms[4].coords,atoms[0].coords,
                           atoms[1].coords,atoms[2].coords)
        t2 = torsion_angle(atoms[0].coords,atoms[1].coords,
                           atoms[2].coords,atoms[3].coords)
        t3 = torsion_angle(atoms[1].coords,atoms[2].coords,
                           atoms[3].coords,atoms[4].coords)
        t4 = torsion_angle(atoms[2].coords,atoms[3].coords,
                           atoms[4].coords,atoms[0].coords)
        t5 = torsion_angle(atoms[3].coords,atoms[4].coords,
                           atoms[0].coords,atoms[1].coords)

        if math.fabs(t1) < 7.0:
            a1 = atoms[0]
            a2 = atoms[1]
            bond = self.get_bond(ring['bonds'], a1, a2 )
            bond.bondOrder = 1
            dist = distance(a1.coords, a2.coords)
            cov_sum = a1.babel_bond_ord_rad + a2.babel_bond_ord_rad
            ratio = dist/cov_sum
            if ratio < SINGLE_DOUBLE_CUTOFF:
                bond.bondOrder = 2

        if math.fabs(t2) < 7.0:
            a1 = atoms[1]
            a2 = atoms[2]
            bond = self.get_bond(ring['bonds'], a1, a2 )
            bond.bondOrder = 1
            dist = distance(a1.coords, a2.coords)
            cov_sum = a1.babel_bond_ord_rad + a2.babel_bond_ord_rad
            ratio = dist/cov_sum
            if ratio < SINGLE_DOUBLE_CUTOFF:
                bond.bondOrder = 2

        if math.fabs(t3) < 7.0:
            a1 = atoms[2]
            a2 = atoms[3]
            bond = self.get_bond(ring['bonds'], a1, a2 )
            bond.bondOrder = 1
            dist = distance(a1.coords, a2.coords)
            cov_sum = a1.babel_bond_ord_rad + a2.babel_bond_ord_rad
            ratio = dist/cov_sum
            if ratio < SINGLE_DOUBLE_CUTOFF:
                bond.bondOrder = 2

        if math.fabs(t4) < 7.0:
            a1 = atoms[3]
            a2 = atoms[4]
            bond = self.get_bond(ring['bonds'], a1, a2 )
            bond.bondOrder = 1
            dist = distance(a1.coords, a2.coords)
            cov_sum = a1.babel_bond_ord_rad + a2.babel_bond_ord_rad
            ratio = dist/cov_sum
            if ratio < SINGLE_DOUBLE_CUTOFF:
                bond.bondOrder = 2

        if math.fabs(t5) < 7.0:
            a1 = atoms[4]
            a2 = atoms[0]
            bond = self.get_bond(ring['bonds'], a1, a2 )
            bond.bondOrder = 1
            dist = distance(a1.coords, a2.coords)
            cov_sum = a1.babel_bond_ord_rad + a2.babel_bond_ord_rad
            ratio = dist/cov_sum
            if ratio < SINGLE_DOUBLE_CUTOFF:
                bond.bondOrder = 2


    def estimate_bond_order2(self, bonds):
        """ """
        converter = TypeConverter("HYB")
        for b in bonds:
            bo = 1
            a1 = b.atom1
            a2 = b.atom2
            dist = distance(a1.coords, a2.coords)
            cov_sum = a1.babel_bond_ord_rad + a2.babel_bond_ord_rad
            ratio = dist/cov_sum
            start_type = converter.convert(a1.babel_type, "all_caps")
            end_type = converter.convert(a2.babel_type, "all_caps")

            if ratio <= DOUBLE_TRIPLE_CUTOFF:
                if start_type[0] == '1' and end_type[0] == '1':
                    bo = 3
            elif ratio <= SINGLE_DOUBLE_CUTOFF:
                if start_type[0] == '2' and end_type[0] == '2':
                    bo = 2
            b.bondOrder = bo
       
        
    def assign_hybrid_radii(self, atoms):
        """ """
        converter = TypeConverter("XYZ")
        for a in atoms:
            atm_type = converter.convert(a.babel_type, 'zero')
            if atm_type == 0:
                atm_type = a.babel_type
            atm_type = converter.clean_atom_type(atm_type)

            a.babel_cov_rad = babel_elements[atm_type]['cov_rad']
            a.babel_bond_ord_rad = babel_elements[atm_type]['bond_ord_rad']
            a.babel_max_bonds = babel_elements[atm_type]['max_bonds']


if __name__ == '__main__':
    import pdb, sys
    from atomTypes import AtomHybridization
    from cycle import RingFinder
    
    from MolKit.pdbParser import NewPdbParser
    parser = NewPdbParser("/tsri/pdb/struct/%s.pdb"%sys.argv[1])
    mols = parser.parse()
    mol = mols[0]
    mol.buildBondsByDistance()
    allAtoms = mol.chains.residues.atoms
    bonds = allAtoms.bonds[0]

    print "assigning atom types"
    babel = AtomHybridization()
    babel.assignHybridization(allAtoms)

    print "looking for rings"
    rings = RingFinder()
    rings.findRings(allAtoms, bonds)

    print "assigning bond order"
    bo = BondOrder()
    #pdb.run("bo.assignBondOrder(allAtoms, bonds)")
    bo.assignBondOrder(allAtoms, bonds)

    for b in bonds:
        if b.bondOrder > 1:
            a1 = b.atom1
            a2 = b.atom2
            print '%-20s %-20s %d'% ( a1.full_name(), a2.full_name(),
                                      b.bondOrder )
