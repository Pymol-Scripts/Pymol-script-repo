#############################################################################
#
# Author: Michel F. SANNER
# Reimplemented from Babel v1.6 from Pat Walters and Math Stahl
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################
#
# $Header: /opt/cvs/python/packages/share1.5/PyBabel/aromatic.py,v 1.3 2008/11/15 00:14:40 sargis Exp $
#
# $Id: aromatic.py,v 1.3 2008/11/15 00:14:40 sargis Exp $
#

"""
This file implements the Aromatic class.

Before this Aromatic object can be used to identify aromatic rings,
 1 - atoms must have been assigned a type see (AtomHybridization in types.py)
 2 - bond orders must have been defined
 3 - rings must have been identified (see RingFinder in cycle.py)

example:
    
      >>> atype = AtomHybridization()
      >>> atype.assignHybridization(atoms)
      >>> rings = RingFinder()
      >>> rings.findRings(atoms, bonds)
      >>> bo = BondOrder()
      >>> bo.assignBondOrder(atoms, bonds, rings)
      >>> arom = Aromatic(rings)
      >>> arom.find_aromatic_atoms(atoms)

      atoms has to be a list of atom objects
      Atom:
          a.coords : 3-sequence of floats
          a.bonds : list of Bond objects
          babel_type: string
          babel_atomic_number: int

      Bond:
          b.atom1 : instance of Atom
          b.atom2 : instance of Atom

      After completion
      
reimplmentation of Babel1.6 in Python by Michel Sanner April 2000
Original code by W. Patrick Walters and Matthew T. Stahl 
"""


import math

from atomTypes import TypeConverter
from util import distance, bond_angle

class Aromatic:
    """ """

    def __init__(self, rings):
        """ """
        self.rings = rings
        
    
    def find_aromatic_atoms(self, atoms):
        """ """
        info = self.setup_ring_info( atoms )

        for r in self.rings.rings:
            for a in r['atoms']:
                if a.arom_atm: # atom is aromatic
                    if a.babel_atomic_number == 6:
                        a.babel_type = "Car"
                    elif a.babel_atomic_number == 7:
                        a.babel_type = "Nar"

            for b in r['bonds']:
                if b.atom1.arom_atm and b.atom2.arom_atm:
                    b.bondOrder = 'aromatic'

        for r in self.rings.rings:
            for a in r['atoms']:
                if hasattr(a, 'arom_atm'): delattr(a, 'arom_atm')

                
    def find_aromatic_rings(self):
        """ """
        converter = TypeConverter("HYB")
       
        for r in self.rings.rings:
            is_aromatic = 0
            if len(r['atoms']) == 5 or len(r['atoms']) == 6:
                is_aromatic = 1
                for a in r['atoms']:
                    # check to see if the atom is sp2 
	            # if not, the ring isn't aromatic
                    hyb_str = converter.convert(a.babel_type, 'all_caps')
                    hyb = int(hyb_str)
                    if hyb != 2:
                        is_aromatic = 0
                        break

                    # check to see if the atom has 3 bonds to
                    # heavy atoms in the same ring.
                    # If not, the ring isn't aromatic

                    if self.count_arom_atm_bonds(a, r) != 3:
                        is_aromatic = 0
                        break
    
            if is_aromatic:
                r['aromatic'] = 1
                for a in r['atoms']: a.arom_atm = 1
                

    def find_aromatic_rings2(self):
        """ """
        converter = TypeConverter("HYB")
       
        for r in self.rings.rings:
            is_aromatic = 0
            if len(r['atoms']) == 5 or len(r['atoms']) == 6:
                is_aromatic = 1
                for a in r['atoms']:
                    # check to see if the atom is sp2 
                    # if not, the ring isn't aromatic

                    hyb_str = converter.convert(a.babel_type, 'all_caps')
                    hyb = int(hyb_str)
                    if hyb != 2:
                        is_aromatic = 0
                        break

                    # check to see if the atom has 3 bonds to
                    # heavy atoms in the same ring.
                    # If not, the ring isn't aromatic

#                    if not a.arom_atm:
#                        if self.count_arom_atm_bonds(a, r) != 3:
#                            is_aromatic = 0
#                            break
    
            if is_aromatic:
                r['aromatic'] = 1
                for a in r['atoms']: a.arom_atm = 1


    def setup_ring_info(self, atoms):
        """ """        
        for r in self.rings.rings:
            r['aromatic'] = 0
            for a in r['atoms']: a.arom_atm = 0
            
        self.find_aromatic_rings()
        self.find_aromatic_rings2()

        for r in self.rings.rings:
            # compute cycle center
            center = [0,0,0]
            sca = 1.0/len(r['atoms'])
            for a in r['atoms']:
                c = a.coords
                center[0] = center[0] + c[0]*sca
                center[1] = center[1] + c[1]*sca
                center[2] = center[2] + c[2]*sca
            r['center'] = center

            # compute cycle radius
            b = r['bonds'][0]
            c1 = b.atom1.coords
            c2 = b.atom2.coords
            mid = ( (c1[0]+c2[0])*0.5, (c1[1]+c2[1])*0.5, (c1[2]+c2[2])*0.5)
            r['radius'] = 0.9 * distance(center, mid)
            
            # compute cycle normal vector
            n_avg = [0,0,0]
            for i in range(len(r['bonds'])-1):
                b = r['bonds'][i]
                c1 = b.atom1.coords 
                c2 = b.atom2.coords 
                a3 = r['bonds'][i+1].atom1
                if a3==b.atom1 or a3==b.atom2: a3 = r['bonds'][i+1].atom2
                c3 = a3.coords
                v = (c2[0]-c1[0], c2[1]-c1[1], c2[2]-c1[2])
                v1 = (c3[0]-c1[0], c3[1]-c1[1], c3[2]-c1[2])
                n = [ v[1]*v1[2] - v1[1]*v[2],
                      v1[0]*v[2] - v[0]*v1[2],
                      v[0]*v1[1] - v1[0]*v[1] ]
                sca = 1.0 / math.sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2])
                n_avg = (n[0]*sca, n[1]*sca, n[2]*sca)

            r['normal'] = n_avg

            
    def count_arom_atm_bonds(self, atom, ring):
        """ """
        hvy_bonds = 0
        for b in atom.bonds:
            conn = b.atom1
            if conn==atom: conn=b.atom2
            if conn in ring['atoms']:
            
                if conn.babel_atomic_number != 1:
                    # add bond order of bonds to non hydrogen atoms
                    hvy_bonds = hvy_bonds+b.bondOrder

        return hvy_bonds


if __name__ == '__main__':
    import pdb, sys
    from atomTypes import AtomHybridization
    from cycle import RingFinder
    from bo import BondOrder
    from aromatic import Aromatic
    from addh import AddHydrogens
    
    from MolKit.pdbParser import PdbqParser
    #parser = NewPdbqParser("./txp_noh.pdbq") # doesn't work
    #parser = NewPdbqParser("./txp.pdbq")     # OK
    #parser = NewPdbqParser("./gmmtest.pdbq")  # 
    from MolKit.mol2Parser import NewMol2Parser
    #mols = parser.parse()
    #mol = mols[0]
    #mol.buildBondsByDistance()

    parser = NewMol2Parser("./mtx.mol2") # OK
    mols = parser.parse()
    mol = mols[0]

    allAtoms = mol.chains.residues.atoms

    print "assigning atom types"
    babel = AtomHybridization()
    babel.assignHybridization(allAtoms)

    print "looking for rings"
    bonds = allAtoms.bonds[0]
    rings = RingFinder()
    rings.findRings(allAtoms, bonds)

    print "assigning bond order"
    bo = BondOrder()
    #pdb.run("bo.assignBondOrder(allAtoms, bonds)")
    bo.assignBondOrder(allAtoms, bonds)

    def countMissingBonds(atoms):
        converter = TypeConverter('HYB')
        needDBatoms = []
        canHaveDBatoms = []
        for a in atoms:
            a._needsDB = 0
            a._canHaveDB = 0
            hyb = int(converter.convert(a.babel_type, 'all_caps'))
            if hyb != 2: continue         
            if a.babel_type[0] == 'C':
                sum = 0
                for b in a.bonds:
                    sum = sum + b.bondOrder

                if sum == 3:
                    if len(a.bonds)==3:
                        a._needsDB = 1
                        needDBatoms.append(a)
                    elif len(a.bonds)==2: # ==C-- or ==C/
                        a1 = a.bonds[0].atom1
                        if a1==a: a1 = a.bonds[0].atom2
                        a3 = a.bonds[1].atom1
                        if a3==a: a3 = a.bonds[1].atom2
                        ang = bond_angle( a1.coords, a.coords, a3.coords )
                        if math.fabs(ang-180) < 10: # ==C--
                            a._needDB = 1
                            needDBatoms.append(a)
                            
                elif sum == 2:
                    if len(a.bonds)==2:
                        a1 = a.bonds[0].atom1
                        if a1==a: a1 = a.bonds[0].atom2
                        a3 = a.bonds[1].atom1
                        if a3==a: a3 = a.bonds[1].atom2
                        ang = bond_angle( a1.coords, a.coords, a3.coords )
                        if math.fabs(ang-180) < 10:
                            a._needsDB = 2
                            needDBatoms.append(a)
                        else:
                            a._needsDB = 1
                            needDBatoms.append(a)

                elif sum == 1:
                    a._needsDB = 1
                    needDBatoms.append(a)
                    
            elif a.babel_type[0] == 'N':
                sum = 0
                for b in a.bonds:
                    sum = sum + b.bondOrder
                if sum == 2 and len(a.bonds)==2:
                    canHaveDBatoms.append(a)
                    a._canHaveDB = 1

            elif a.babel_type[0] == 'O':
                sum = 0
                for b in a.bonds:
                    sum = sum + b.bondOrder
                if sum == 1:
                    canHaveDBatoms.append(a)
                    a._canHaveDB = 1

        return needDBatoms, canHaveDBatoms


    def fixMissingBonds(pbatoms):
        i = 0
        while (len(pbatoms) and i <len(pbatoms)):
            atom = pbatoms[i]
            possibleNeighbors = []
            for b in atom.bonds:
                atom2 = b.atom1
                if atom2==atom: atom2=b.atom2
                if atom2._needsDB: possibleNeighbors.append( (atom2, b) )
            if len(possibleNeighbors)==atom._needsDB:
                for entry in possibleNeighbors:
                    b = entry[1]
                    a2 = entry[0]
                    b.bondOrder = b.bondOrder+1
                    a2._needsDB = a2._needsDB -1
                    atom._needsDB = atom._needsDB -1
                    if a2._needsDB==0: pbatoms.remove(a2)
                pbatoms.remove(atom)
            else:
                i = i +1
        return pbatoms


    def fixMissingBondOrders():
        needDB, canHave = countMissingBonds(allAtoms)
        while (1):
            l = len(needDB)
            needDB = fixMissingBonds(needDB)
            if len(needDB)==0: break
            if l==len(needDB): break

        print 'NEED'
        for a in needDB:
            print a
        print 'canHave'
        for a in canHave:
            print a
        all = []
##          while (1):
##              l = len(all)
##              all = fixMissingBonds(all)
##              if len(all)==0: break
##              if l==len(all): break
        return all
    
    pbatoms = fixMissingBondOrders()
    print "WARNING couldn't fix all bond orders", pbatoms

    print "looking for aromatic rings"
    arom = Aromatic(rings)
    arom.find_aromatic_atoms(allAtoms)
    print "done"

    import math

    def addDispVector():
        for b in bonds:                
            if b.bondOrder != 1:
                inring = hasattr(b.atom1, 'rings')
                c1 = b.atom1.coords
                c2 = b.atom2.coords
                if inring:
                    found = 0
                    for nb in b.atom1.bonds:
                        if found: break
                        if nb==b: continue
                        a3 = nb.atom1
                        if a3==b.atom1: a3 = nb.atom2
                        if hasattr(a3, 'rings'):
                            for ringNum in a3.rings:
                                if hasattr(b.atom1, 'rings'): r1=b.atom1.rings
                                else: r1 = []
                                if hasattr(b.atom2, 'rings'): r2=b.atom2.rings
                                else: r2 = []
                                if ringNum in r1 and ringNum in r2:
                                    found=1
                                    break
                else:
                    nb = b.atom1.bonds[0]
                    if nb==b: nb = b.atom2.bonds[0]
                    if nb==b:
                        if len(b.atom1.bonds)>1:
                            nb = b.atom1.bonds[1]
                        else:
                            nb = b.atom2.bonds[1]
                    a3 = nb.atom1
                    if a3==b.atom1: a3 = nb.atom2
                    
                c3 = a3.coords
                v = (c2[0]-c1[0], c2[1]-c1[1], c2[2]-c1[2])
                v1 = (c3[0]-c1[0], c3[1]-c1[1], c3[2]-c1[2])
                n = [ v[1]*v1[2] - v1[1]*v[2],
                      v1[0]*v[2] - v[0]*v1[2],
                      v[0]*v1[1] - v1[0]*v[1] ]
                nn = [ n[1]*v[2] - v[1]*n[2],
                       v[0]*n[2] - n[0]*v[2],
                       n[0]*v[1] - v[0]*n[1] ]
                sca = 0.15 / math.sqrt(nn[0]*nn[0] + nn[1]*nn[1] + nn[2]*nn[2])
                b.dispVec = (nn[0]*sca, nn[1]*sca, nn[2]*sca)

    addDispVector()
    
##      addh = AddHydrogens()
##      hat = addh.addHydrogens(allAtoms)
    
##      from MolKit.molecule import Atom, Bond
##      for a in hat:
##          atom = Atom('H', a[1].parent, top=a[1].top)
##          atom._coords = [ a[0] ]
##          if hasattr(atom, 'segID'):
##              atom.segID = a[1].segID
##          atom.hetatm = 0
##          atom.alternate = []
##          atom.element = 'H'
##          atom.number = -1
##          atom.occupancy = 1.0
##          atom.conformation = 0
##          atom.temperatureFactor = 0.0
##          atom.babel_atomic_number = a[2]
##          atom.babel_type = a[3]
##          atom.babel_organic=1
##          bond = Bond( a[1], atom )

    from Pmv.moleculeViewer import MoleculeViewer
    mv = MoleculeViewer()
    mv.addMolecule(mol)
    mol.bondsflag = 1
    mv.lines( mol )

    v = []
    l = []
    g = mv.Mols[0].geomContainer.geoms['lines']
    vc = len(g.vertexSet)
    scale = [0, 1, -1, 2, -2, 3, -3, 4, -4]
    add = -1

##      for b in bonds:
##          if b.bondOrder>1:
##              b.bondOrder = 5
##              break

    for b in bonds:
        n = b.bondOrder
        if n=='aromatic': continue
        for i in range(1, n):
            print i, scale[i]
            l.append( (vc,vc+1) )
            c1 = b.atom1.coords
            v.append( c1[0]+(scale[i]*b.dispVec[0]),
                      c1[1]+(scale[i]*b.dispVec[1]),
                      c1[2]+(scale[i]*b.dispVec[2]) )
            vc = vc+1
            c1 = b.atom2.coords
            v.append( c1[0]+(scale[i]*b.dispVec[0]),
                      c1[1]+(scale[i]*b.dispVec[1]),
                      c1[2]+(scale[i]*b.dispVec[2]) )
            vc = vc+1

    g = mv.Mols[0].geomContainer.geoms['lines']
    g.Add(vertices = v, faces=l)

    from DejaVu.Arcs3D import Arcs3D
    aromC = []
    aromR = []
    aromN = []
    for r in rings.rings:
        if r['aromatic']:
            aromC.append(r['center'])
            aromR.append(r['radius'])
            aromN.append(r['normal'])

    c = Arcs3D('aromCircles', vertices = aromC, vnormals = aromN,
                materials = ( (0,1,0), ), radius = aromR,
                angles = ( 360, ), stippleLines = 1, lineWidth=2)

    mv.Mols[0].geomContainer.geoms['aromCircles'] = c
    mv.GUI.VIEWER.AddObject(c, parent = mv.Mols[0].geomContainer.masterGeom)
