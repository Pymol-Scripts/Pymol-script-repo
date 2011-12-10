#############################################################################
#
# Author: Michel F. SANNER
# Reimplemented from Babel v1.6 from Pat Walters and Math Stahl
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################
#
# $Header: /opt/cvs/python/packages/share1.5/PyBabel/addh.py,v 1.6 2007/10/11 17:43:48 sargis Exp $
#
# $Id: addh.py,v 1.6 2007/10/11 17:43:48 sargis Exp $
#

"""
This file implements the AddHydrogens class.

Before this AddHydrogens object can be used, atoms must have been assigned
a type see (AtomHybridization in types.py).

Hydrogen atoms can be added using 2 different methods. The first one requires
bondOrders to have been calculated previousely.

example:
    
      >>> atype = AtomHybridization()
      >>> atype.assignHybridization(atoms)
      >>> addh = AddHydrogens()
      >>> hat = addh.addHydrogens(atoms)

      atoms has to be a list of atom objects
      Atom:
          a.coords : 3-sequence of floats
          a.bonds : list of Bond objects
          babel_type: string
          babel_atomic_number: int

      Bond:
          b.atom1 : instance of Atom
          b.atom2 : instance of Atom

      or
      
      >>> addh = AddHydrogens()
      >>> hat = addh.addHydrogens(atoms, method='withBondOrder')

      atoms has to be a list of atom objects as above and
      Bond:
          b.atom1 : instance of Atom
          b.atom2 : instance of Atom
          b.bondOrder : integer
          
      these calls return a list 'hat' containing a tuple for each Hydrogen
      atom to be added. This tuple provides:
          coordsH       : 3-float coordinates of the H atom
          atom          : the atom to which the H atom is to be connected
          atomic_number : the babel_atomic_number of the H atom
          type          : tyhe babel_type of the H atom

reimplmentation of Babel1.6 in Python by Michel Sanner April 2000
Original code by W. Patrick Walters and Matthew T. Stahl 
"""


import math

from atomTypes import TypeConverter
from util import *

ONE_OVER_SQRT3 = 0.577350269
SQRT_TWO_THIRDS = 0.816496581

SP3_C_H_DIST = 1.115
SP2_C_H_DIST = 1.103
SP_C_H_DIST = 1.090

SP3_N_H_DIST = 1.020
SP2_N_H_DIST = 1.020

SP3_O_H_DIST = 0.950

class AddHydrogens:
    """ """

    def addHydrogens(self, atoms, method='noBondOrder'):
        """ """
	Hatoms = []
#        method = 'noBondOrder'
        if method=='noBondOrder':
            num_H_to_add = self.count_missing_hydrogens(atoms)
            if num_H_to_add:
                Hatoms = self.place_hydrogens1(atoms, num_H_to_add)

        else:
            num_H_to_add = self.count_missing_bo_hydrogens(atoms)
            if num_H_to_add:
                Hatoms = self.place_hydrogens2(atoms, num_H_to_add)

            # cleanup
            for a in atoms:
                delattr(a, '_redo')

        return Hatoms


    def place_hydrogens1(self, atoms, num_H_to_add):
        """ """

        Hat = []
        for a in atoms:

            if a.babel_type == "C3":
                val = len(a.bonds)
                if val == 3:
                    Hat = Hat + self.add_tertiary_hydrogen(a, SP3_C_H_DIST)
                elif val == 2:
                    Hat = Hat + self.add_methylene_hydrogens(a ,SP3_C_H_DIST)
                elif val==1:
                    Hat = Hat + self.add_methyl_hydrogen(a, SP3_C_H_DIST)
                    Hat = Hat + self.add_methylene_hydrogens(a, SP3_C_H_DIST,
                                                             Hat[-1])
                    
            elif a.babel_type == "N3+": 
                val = len(a.bonds)
                if val == 2:
                    Hat = Hat + self.add_methylene_hydrogens(a ,SP3_N_H_DIST)
                elif val==1:
                    Hat = Hat + self.add_methyl_hydrogen(a, SP3_N_H_DIST)
                    Hat = Hat + self.add_methylene_hydrogens(a ,SP3_N_H_DIST,
                                                             Hat[-1])

            elif a.babel_type == "C2" or a.babel_type == "Car":
                val = len(a.bonds)
                if val == 2:
                    Hat = Hat + self.add_sp2_hydrogen(a ,SP2_C_H_DIST)
                
            elif a.babel_type == "Npl" or a.babel_type == "Nam" or \
                 a.babel_type == "Ng+":
                val = len(a.bonds)
                if val == 2:
                    Hat = Hat + self.add_sp2_hydrogen(a ,SP2_N_H_DIST)

            elif a.babel_type == "C1":
                if len(a.bonds) == 1:
                    Hat = Hat + self.add_sp_hydrogen(a ,SP_C_H_DIST)

            elif a.babel_type == "O3":
                if len(a.bonds) == 1:
                    Hat = Hat + self.add_methyl_hydrogen(a, SP3_O_H_DIST)


        for a in atoms:

            if a.babel_type == "C2":
                if len(a.bonds) == 1:
                    Hat = Hat + self.add_vinyl_hydrogens(a ,SP2_C_H_DIST)

            elif a.babel_type == "Npl" or a.babel_type == "Nam" or \
                 a.babel_type == "Ng+":
                if len(a.bonds) == 1:
                    # FIXME babel C code says SP2_C_H_DIST here ???
                    Hat = Hat + self.add_vinyl_hydrogens(a ,SP2_N_H_DIST)

        return Hat


    def place_hydrogens2(self, atoms, num_H_to_add):
        """ """

        Hat = []
        converter = TypeConverter("HYB")
        
        for a in atoms:
            type_name = converter.convert(a.babel_type, 'zero')
            hyb = int(type_name)
            code = a.babel_atomic_number * 10 + hyb
            to_add = a._redo
            
            if code == 63: # sp3 carbon
                if to_add==1:
                    Hat = Hat + self.add_tertiary_hydrogen(a, SP3_C_H_DIST)
                elif to_add==2:
                    Hat = Hat + self.add_methylene_hydrogens(a ,SP3_C_H_DIST)
                elif to_add==3:
                    Hat = Hat + self.add_methyl_hydrogen(a, SP3_C_H_DIST)
                    Hat = Hat + self.add_methylene_hydrogens(a, SP3_C_H_DIST,
                                                             Hat[-1])
            elif code == 73: # sp3 nitrogen
                if to_add==1:
                    if a.babel_type=="N3+":
                        Hat = Hat + self.add_tertiary_hydrogen(a, SP3_N_H_DIST)
                    else:
                        Hat = Hat + self.add_sp3_N_hydrogen(a, SP3_N_H_DIST)
                elif to_add==2:
                    Hat = Hat + self.add_methylene_hydrogens(a ,SP3_N_H_DIST)
                elif to_add==3:
                    Hat = Hat + self.add_methyl_hydrogen(a, SP3_N_H_DIST)
                    Hat = Hat + self.add_methylene_hydrogens(a ,SP3_N_H_DIST,
                                                             Hat[-1])
                    
            elif code == 62: # sp2 carbon
                if to_add==1:
                    Hat = Hat + self.add_sp2_hydrogen(a ,SP2_C_H_DIST)

            elif code == 72: # sp2 nitrogen
                if to_add==1:
                    Hat = Hat + self.add_sp2_hydrogen(a ,SP2_N_H_DIST)

            elif code == 61: # sp carbon
                if to_add==1:
                    Hat = Hat + self.add_sp_hydrogen(a ,SP_C_H_DIST)

            elif code == 83: # sp3 oxygen
                if to_add==1:
                    Hat = Hat + self.add_methyl_hydrogen(a, SP3_O_H_DIST)

            # save vinyl and amide protons for last,
            # this way we know where to put them

        for a in atoms:

            type_name = converter.convert(a.babel_type, 'zero')
            hyb = int(type_name)
            code = a.babel_atomic_number * 10 + hyb
            to_add = a._redo

            if code == 62: # sp2 carbon
                if to_add==2:
                    Hat = Hat + self.add_vinyl_hydrogens(a ,SP2_C_H_DIST)

            elif code == 72: # sp2 nitrogens
                if to_add==2:
                    Hat = Hat + self.add_vinyl_hydrogens(a ,SP2_N_H_DIST)

        return Hat


    def type_added_hydrogen(self, atom):
        """ return babel_atomic_number and babel_type for adde H atom"""
  
        atomic_number = 1
        if atom.babel_atomic_number == 6:
            htype = "HC"
        else:
            htype = "H"
            
        return atomic_number, htype


    def add_sp3_N_hydrogen(self, atom, b_length):
        """ """
        c2 = atom.bonds[0].atom1
        if c2==atom: c2 = atom.bonds[0].atom2

        if not addedH:
            c3 = atom.bonds[1].atom1
            if c3==atom: c3 = atom.bonds[1].atom2
            c3 = c3.coords
        else:
            c3 = addedH[0]

        c = atom.coords
        v = vec3(c2.coords, c)
        v1 = vec3(c3, c)
        s = [ v[0]+v1[0], v[1]+v1[1], v[2]+v1[2] ]

        n = [ v[1]*v1[2] - v1[1]*v[2],
              v1[0]*v[2] - v[0]*v1[2],
              v[0]*v1[1] - v1[0]*v[1] ]
       
        s = [ s[0]*ONE_OVER_SQRT3, s[1]*ONE_OVER_SQRT3, s[2]*ONE_OVER_SQRT3 ]
        n = [ n[0]*SQRT_TWO_THIRDS, n[1]*SQRT_TWO_THIRDS,n[2]*SQRT_TWO_THIRDS ]

        h1 = [ s[0]+n[0], s[1]+n[1], s[2]+n[2] ]
        mag= b_length / math.sqrt( h1[0]*h1[0] + h1[1]*h1[1] + h1[2]*h1[2])
        cH1 = [c[0] + (h1[0]*mag), c[1] + (h1[1]*mag), c[2] + (h1[2]*mag)]

        atomic_number, type = self.type_added_hydrogen(atom)

        return [(cH1, atom, atomic_number, type)]

            
    def add_methyl_hydrogen(self, atom, b_length):
        """ """

        b = atom.bonds[0]
        c = atom.coords

        # find atom2 such that (1) atom2 is bonded to atom, and
        #                      (2) atom2 != atom
        atom2 = b.atom1
        c2 = b.atom1.coords
        if b.atom1 == atom:
            atom2 = b.atom2
            c2 = b.atom2.coords

        # find c3 such that (1) atom3 is bonded to atom2, and
        #                   (2) atom3 != atom2, and
        #                   (3) atom3 != atom
        if atom2.bonds[0].atom1 == atom:
            # atom has been located among the bonds of atom2, so
            # move onto the next bond
            c3 = atom2.bonds[1].atom1.coords
            if atom2.bonds[1].atom1==atom2:
                # since atom was part of bonds[0], the other atom of
                # this bond is safe to use (ie. not atom or atom2)
                c3 = atom2.bonds[1].atom2.coords
        else:
            c3 = atom2.bonds[0].atom1.coords
            if atom2.bonds[0].atom1==atom2:
                c3 = atom2.bonds[0].atom2.coords
                # but atom2.bonds[0].atom2 might be atom, so check
                if atom2.bonds[0].atom2==atom:
                    # same logic as above...
                    c3 = atom2.bonds[1].atom1.coords
                    if atom2.bonds[1].atom1==atom2:
                        c3 = atom2.bonds[1].atom2.coords

        # FIX ME: this only works if atom2 has reasonable sp3 geometry!!!
        v = vec3(c3, c2, b_length)
        coordsH = [ c[0]+v[0], c[1]+v[1], c[2]+v[2] ]
        
        atomic_number, type = self.type_added_hydrogen(atom)
        return [(coordsH, atom, atomic_number, type)]


    def add_tertiary_hydrogen(self, atom, b_length):
        """ """
        c2 = atom.bonds[0].atom1
        if c2==atom: c2 = atom.bonds[0].atom2

        c3 = atom.bonds[1].atom1
        if c3==atom: c3 = atom.bonds[1].atom2

        c4 = atom.bonds[2].atom1
        if c4==atom: c4 = atom.bonds[2].atom2

        c = atom.coords
        v1 = vec3(c2.coords, c)
        v2 = vec3(c3.coords, c)
        v3 = vec3(c4.coords, c)

        m = invert_3x3( (v1,v2,v3) )

        s = [ m[0][0] + m[0][1] + m[0][2],
              m[1][0] + m[1][1] + m[1][2],
              m[2][0] + m[2][1] + m[2][2] ]

        mag= b_length / math.sqrt( s[0]*s[0] + s[1]*s[1] + s[2]*s[2])
        cH = [c[0]+ (s[0]*mag), c[1] + (s[1]*mag), c[2] + (s[2]*mag)]

        atomic_number, type = self.type_added_hydrogen(atom)
        return [(cH, atom, atomic_number, type)]


    def add_methylene_hydrogens(self, atom, b_length, addedH=None):
        """ """
        c2 = atom.bonds[0].atom1
        if c2==atom: c2 = atom.bonds[0].atom2

        if not addedH:
            c3 = atom.bonds[1].atom1
            if c3==atom: c3 = atom.bonds[1].atom2
            c3 = c3.coords
        else:
            c3 = addedH[0]
            
        c = atom.coords
        v = vec3(c2.coords, c)
        v1 = vec3(c3, c)
        s = [ v[0]+v1[0], v[1]+v1[1], v[2]+v1[2] ]
        
        n = [ v[1]*v1[2] - v1[1]*v[2],
              v1[0]*v[2] - v[0]*v1[2],
              v[0]*v1[1] - v1[0]*v[1] ]

        s = [ s[0]*ONE_OVER_SQRT3, s[1]*ONE_OVER_SQRT3, s[2]*ONE_OVER_SQRT3 ]
        n = [ n[0]*SQRT_TWO_THIRDS, n[1]*SQRT_TWO_THIRDS,n[2]*SQRT_TWO_THIRDS ]

        h1 = [ s[0]+n[0], s[1]+n[1], s[2]+n[2] ]
        mag= b_length / math.sqrt( h1[0]*h1[0] + h1[1]*h1[1] + h1[2]*h1[2])
        cH1 = [c[0] + (h1[0]*mag), c[1] + (h1[1]*mag), c[2] + (h1[2]*mag)]

        h2 = [ s[0]-n[0], s[1]-n[1], s[2]-n[2] ]
        mag= b_length / math.sqrt( h2[0]*h2[0] + h2[1]*h2[1] + h2[2]*h2[2])
        cH2 = [c[0] + (h2[0]*mag), c[1] + (h2[1]*mag), c[2] + (h2[2]*mag)]

        atomic_number, type = self.type_added_hydrogen(atom)
 
        return [(cH1, atom, atomic_number, type),
                (cH2, atom, atomic_number, type)]


    def add_sp2_hydrogen(self, atom, b_length):
        """ """

        c2 = atom.bonds[0].atom1
        if c2==atom: c2 = atom.bonds[0].atom2

        c3 = atom.bonds[1].atom1
        if c3==atom: c3 = atom.bonds[1].atom2

        c = atom.coords
        v = vec3(c2.coords, c )
        v1 = vec3(c3.coords, c )
        s = [ v[0]+v1[0], v[1]+v1[1], v[2]+v1[2] ]
        mag= b_length / math.sqrt( s[0]*s[0] + s[1]*s[1] + s[2]*s[2])

        coordsH = [c[0] + (s[0]*mag), c[1] + (s[1]*mag), c[2] + (s[2]*mag)]

        atomic_number, type = self.type_added_hydrogen(atom)
        return [(coordsH, atom, atomic_number, type)]


    def add_sp_hydrogen(self, atom, b_length):
        """ """
        b = atom.bonds[0]
        c = atom.coords
        if b.atom1 == atom:
            c2 = b.atom2.coords
        else:
            c2 = b.atom1.coords

        v = vec3(c2, c, b_length)
        coordsH = [ c[0]+v[0], c[1]+v[1], c[2]+v[2] ]
        
        atomic_number, type = self.type_added_hydrogen(atom)
        return [(coordsH, atom, atomic_number, type)]


    def add_vinyl_hydrogens(self, atom, b_length):
        """ """

        # c2 is neigbour of c1
        c2 = atom.bonds[0].atom1
        if c2==atom:
            c2= atom.bonds[0].atom2

        # c3 is neigbour of c2 different from c1
        c3 = c2.bonds[0].atom1
        if c3==c2: c3 = c2.bonds[0].atom2

        if c3==atom:
            c3 = c2.bonds[1].atom1
            if c3 == c2:
                c3 = c2.bonds[1].atom2

        # c4 is neighbor of c2 different from c1 and c3
        c4 = c2.bonds[0].atom1
        if c4==c2:
            c4 = c2.bonds[0].atom2
        if c4 == atom or c4 == c3:
            c4 = c2.bonds[1].atom1
            if c4 == c2:
                c4 = c2.bonds[1].atom2
            if c4==atom or c4 == c3:
                if len(c2.bonds)>2:
                    c4 = c2.bonds[2].atom1
                    if c4==c2:
                        c4 = c2.bonds[2].atom2
                else:
                    c4 = None
        # here we diverge from Babel1.6 because they use H atoms
        # that might have laready been added. We will only add them
        # later. Therefore c2 might not have 3 neighbors
        c = atom.coords
        v = vec3(c3.coords, c2.coords, b_length)
        cH1 = [ c[0]+v[0], c[1]+v[1], c[2]+v[2] ]
        if c4 is not None and len(c2.bonds)==3:
            v1 = vec3(c4.coords, c2.coords, b_length)
            cH2 = [ c[0]+v1[0], c[1]+v1[1], c[2]+v1[2]]
        else:
            # we build second H like a sp2_hydrogen using c2,c,H
            v = vec3(c2.coords, c)
            v1 = vec3(cH1, c, b_length )
            s = [ v[0]+v1[0], v[1]+v1[1], v[2]+v1[2] ]
            mag= b_length / math.sqrt( s[0]*s[0] + s[1]*s[1] + s[2]*s[2])
            cH2 = [c[0] + (s[0]*mag), c[1] + (s[1]*mag), c[2] + (s[2]*mag)]

        atomic_number, htype = self.type_added_hydrogen(atom)

        return [ (cH1, atom, atomic_number, htype),
                 (cH2, atom, atomic_number, htype) ]

    
    def count_missing_hydrogens(self, atoms):
        """ """
        missing = 0
        converter = TypeConverter("HAD")
        
        for a in atoms:
            temp_type = converter.convert(a.babel_type, 'zero')

            if temp_type == 0:
                print "Unable to assign valence to atom %s type = %s" % \
                      (a.full_name(), a.babel_type)

            type_valence = int(temp_type)
            val = len(a.bonds)
            if val < type_valence and val > 0:
                missing = missing + type_valence - val

        return missing


    def count_missing_bo_hydrogens(self, atoms):
        """ """
        missing = 0
        for a in atoms:
            type_valence = 0

            if a.babel_atomic_number==6:
                type_valence = 4

            elif a.babel_atomic_number==7:
                if a.babel_type=='N2' and len(a.bonds)==1:
                    type_valence = 2
                elif a.babel_type=='N3+':
                    type_valence = 4
                else:
                    type_valence = 3

            elif a.babel_atomic_number==8:
                if a.babel_type=="O-" or a.babel_type=="O2":
                    type_valence = 1
                else:
                    type_valence = 2
  
            attached = self.count_attached_bonds(a)
            
            a._redo = 0
            to_add = 0
            if len(a.bonds) < type_valence and len(a.bonds) > 0:
                to_add = type_valence - attached
                if to_add > 0:
                    a._redo = to_add
                    missing = missing + to_add

        return missing


    def count_attached_bonds(self, atom):
        """ """
        bonds = 0.0;
        for b in atom.bonds:
            if b.bondOrder == 5:
                bonds = bonds + 1.5
            elif b.bondOrder == 'aromatic':
                bonds = bonds + 1.5
            elif b.bondOrder == 'amide':
                bonds = bonds + 1.5
            else:
                bonds = bonds + b.bondOrder
        return int(bonds)


if __name__ == '__main__':
    import pdb, sys
    from cycle import RingFinder
    from bo import BondOrder
    from atomTypes import AtomHybridization
    from aromatic import Aromatic
    
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

    print "looking for aromatic rings"
    arom = Aromatic(rings)
    arom.find_aromatic_atoms(allAtoms)
    
    print "done"
##      db = filter(lambda x:x.bondOrder==2, bonds)
##      for b in db:
##          print b
    
    addh = AddHydrogens()
    #pdb.run("hat = addh.addHydrogens(allAtoms)")
    hat = addh.addHydrogens(allAtoms)
    
    from MolKit.molecule import Atom, Bond
    for a in hat:
        atom = Atom('H', a[1].parent, top=a[1].top)
        atom._coords = [ a[0] ]
        atom.segID = a[1].segID
        atom.hetatm = 0
        atom.alternate = []
        atom.element = 'H'
        atom.number = -1
        atom.occupancy = 1.0
        atom.conformation = 0
        atom.temperatureFactor = 0.0
        atom.babel_atomic_number = a[2]
        atom.babel_type = a[3]
        atom.babel_organic=1
        bond = Bond( a[1], atom )

    from Pmv.moleculeViewer import MoleculeViewer
    mv = MoleculeViewer()
    mv.addMolecule(mol)
    mv.lines( mol )
