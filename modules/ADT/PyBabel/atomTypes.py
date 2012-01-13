#############################################################################
#
# Author: Michel F. SANNER
# Reimplemented from Babel v1.6 from Pat Walters and Math Stahl
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################
#
# $Header: /opt/cvs/python/packages/share1.5/PyBabel/atomTypes.py,v 1.5 2009/04/06 21:35:17 sargis Exp $
#
# $Id: atomTypes.py,v 1.5 2009/04/06 21:35:17 sargis Exp $
#

"""
This file implements the AtomHybridization class that can be used to assign
atom types.

example:
    
      >>> atype = AtomHybridization()
      >>> atype.assignHybridization(atoms)

      atoms has to be a list of atom objects
      Atom:
          a.element : atom's chemical element symbol (string)
          a.coords : 3-sequence of floats
          a.bonds : list of Bond objects
      Bond:
          b.atom1 : instance of Atom
          b.atom2 : instance of Atom

      after completion each atom has the following additional members:
          babel_type: string
          babel_atomic_number: int
          babel_organic
      
reimplmentation of Babel1.6 in Python by Michel Sanner April 2000
Original code by W. Patrick Walters and Matthew T. Stahl 
"""



import string

from babelAtomTypes import babel_types
from babelElements import babel_elements
from util import *


# for valence_three
SP3_MAX     = 114.8
MAY_BE_SP2  = 122.0
SP_MIN      = 160.0

# for valence_one
V1_C1_C1_CUTOFF = 1.22
V1_C2_C_CUTOFF  = 1.41
V1_C2_N_CUTOFF  = 1.37

V1_N1_C1_CUTOFF = 1.20 
V1_N3_C_CUTOFF  = 1.38
V1_N3_N3_CUTOFF = 1.43
V1_N3_N2_CUTOFF = 1.41

V1_O2_C2_CUTOFF = 1.30
V1_O2_AS_CUTOFF = 1.685

V1_S2_C2_CUTOFF = 1.76
V1_S2_AS_CUTOFF = 2.11

V2_C3_C_CUTOFF  = 1.53
V2_C3_N_CUTOFF  = 1.46
V2_C3_O_CUTOFF  = 1.44

V2_N2_C_CUTOFF  = 1.38
V2_N2_N_CUTOFF  = 1.32

V2_C2_C_CUTOFF  = 1.42
V2_C2_N_CUTOFF  = 1.41
GEN_C3_C_CUTOFF = 1.45


class AtomHybridization:


    def __init__(self):
        """constructor"""
        self.atoms = None
        
        
    def get_atomic_number(self, name):
        """return the element number for a given name or raises a
        ValueError exception if the element is not known"""
        _name = string.upper(name[0])
        if len(name)>1:
            if not name[1] in string.digits:
                _name = _name + string.lower(name[1])
        if _name in babel_elements.keys():
            return babel_elements[_name]['num'] 
        else:
            raise ValueError( "Could not find atomic number for %s %s"% \
                              (name,_name) )


    def assignHybridization(self, atoms):
        """atoms is a list of objects of type Atom having the following
        members:
        Atom:
            a.element : atom's chemical element symbol (string)
            a.coords : 3-sequence of floats
            a.bonds : list of Bond objects
        Bond:
            b.atom1 : instance of Atom
            b.atom2 : instance of Atom

        after completion each atom has the following additional members:
        babel_type: string
        babel_atomic_number: int
        babel_organic
        """
        
        self.atoms = atoms
        for a in self.atoms:
            a.babel_type = a.element
            a._babel_redo = 0
            a.babel_atomic_number = self.get_atomic_number(a.babel_type)

            if a.babel_type[0] in [ 'C', 'H', 'O', 'N', 'S', 'P' ]:
                a.babel_organic=1
            else: a.babel_organic = 0
            if a.element=='Ca': a.babel_organic = 0
        
        self.phase1()
        self.valence_four()
        self.valence_three()
        self.valence_two()
        self.valence_one()

        self.phase4()
        self.phase5() 
        self.phase6()
        self.check_for_amides()

        # cleanup
        for a in self.atoms:
            delattr(a,'_babel_redo')

        delattr(self,'atoms')


    def count_heavy_atoms(self, atom):
        count = 0
        for b in atom.bonds:
            bonded_atom = b.atom1
            if bonded_atom==atom: bonded_atom=b.atom2
            if bonded_atom.babel_type[0] == 'H': count = count + 1
        return len(atom.bonds) - count


    def count_free_ox(self, atom):
        free_O_count=0
        for b in atom.bonds:
            bonded_atom = b.atom1
            if bonded_atom==atom: bonded_atom=b.atom2
            if bonded_atom.babel_type[0] == 'O' and \
               self.count_heavy_atoms(bonded_atom) == 1:
                free_O_count = free_O_count+1
        return free_O_count

        
    def phase1(self):
        for a in self.atoms:
            if a.babel_type[0] == 'H':
                a.babel_type = 'H'

                if len(a.bonds):
                    k = a.bonds[0].atom1
                    if k==a: k = a.bonds[0].atom2
                    if k.babel_type[0] == 'C' and k.element == 'C':
                        a.babel_type = 'HC'


    def valence_four(self):

      for a in self.atoms:
        if len(a.bonds) == 4 and a.babel_organic:

            if a.babel_type[0] == 'C' and a.element == 'C':
                if a.babel_type=='C': a.babel_type = "C3"

            elif a.babel_type[0] == 'N':
                if self.count_free_ox(a) >= 1: a.babel_type = "Nox"
                else: a.babel_type = "N3+"

            elif a.babel_type[0] == 'P':
                if len(a.babel_type) == 1:
                    count = self.count_free_ox(a)
                    if count >= 2: a.babel_type = "Pac"
                    elif count == 1: a.babel_type = "Pox"
                    else: a.babel_type = "P3+"

            elif a.babel_type[0] == 'S':
                if a.babel_type=='S':
                    count = self.count_free_ox(a)
                    if count >= 3: a.babel_type = "Sac"
                    elif count >= 1: a.babel_type = "Sox"
                    else: a.babel_type = "S"

            elif a.babel_type[0] == 'B':
                count = self.count_free_ox(a)
                if count >= 3: a.babel_type = "Bac"
                if count >= 1: a.babel_type = "Box"
                else: a.babel_type = "B"


    def valence_three(self):

      for a in self.atoms:
          if len(a.bonds) == 3 and a.babel_organic:

              k = a.bonds[0].atom1
              if k==a: k = a.bonds[0].atom2
              l = a.bonds[1].atom1
              if l==a: l = a.bonds[1].atom2
              m = a.bonds[2].atom1
              if m==a: m = a.bonds[2].atom2

              angle1 = bond_angle(k.coords, a.coords, l.coords)
              angle2 = bond_angle(k.coords, a.coords, m.coords)
              angle3 = bond_angle(l.coords, a.coords, m.coords)
              avg_angle = (angle1 + angle2 + angle3)/3

              if a.babel_type[0] =='C' and a.element == 'C':
                  if avg_angle < SP3_MAX: a.babel_type="C3"
                  elif self.count_free_ox(a) >= 2: a.babel_type="Cac"
                  else: a.babel_type = "C2"

              elif a.babel_type[0] =='N':
                  if avg_angle < SP3_MAX: a.babel_type = "N3"
                  elif self.count_free_ox(a) >= 2: a.babel_type = "Ntr"
                  else: a.babel_type = "Npl"

              elif a.babel_type[0] =='B':
                  if self.count_free_ox(a) >= 1: a.babel_type = "Box"
                  else: a.babel_type = "B"

              elif a.babel_type =='S':
                  if self.count_free_ox(a) >= 1: a.babel_type = "Sox"
                  else: a.babel_type = "S3+"      


    def valence_two(self):

      for a in self.atoms:
          if len(a.bonds) == 2 and a.babel_organic:

              k = a.bonds[0].atom1
              if k==a: k = a.bonds[0].atom2
              l = a.bonds[1].atom1
              if l==a: l = a.bonds[1].atom2
              
              if a.coords == l.coords:
                  print  a.full_name() +" and " +l.full_name() +" have the same coordinates"
                  
              angle1 = bond_angle(k.coords, a.coords, l.coords)

              if a.babel_type[0] == 'C' and a.element == 'C':
                  if a.babel_type =="C":
                      if angle1 < SP3_MAX:
                          a.babel_type = "C3"
                          a._babel_redo = 1
                      elif angle1 < SP_MIN:
                          a.babel_type = "C2"
                          if angle1 < MAY_BE_SP2:
                              a._babel_redo = 3
                      else: a.babel_type = "C1"

              elif a.babel_type[0] == 'N':
                  if angle1 <= SP3_MAX:
                      a.babel_type = "N3"
                      a._babel_redo = 2
                  elif angle1 <= SP_MIN: a.babel_type = "Npl"
                  else: a.babel_type = "N1"

              elif a.babel_type[0] == 'O':
                a.babel_type = "O3"

              elif a.babel_type[0] == 'S':
                  if a.babel_type =="S": a.babel_type = "S3"


    def valence_one(self):

      for a in self.atoms:

          if len(a.bonds) == 1 and a.babel_organic:
              k = a.bonds[0].atom1
              if k==a: k = a.bonds[0].atom2
              bond_length = distance(a.coords, k.coords)

              if a.babel_type[0] == 'C' and a.element == 'C':
                  if a.babel_type=="C":
                      if k.babel_type[:2]=='C1' and \
                         bond_length <= V1_C1_C1_CUTOFF:
                          a.babel_type = "C1"
                      elif k.babel_type[0] == "C" and \
                           bond_length <= V1_C2_C_CUTOFF:
                          a.babel_type = "C2"
                      else: a.babel_type = "C3"

                  if k.babel_type[0]=="N":
                      if bond_length <= V1_C2_N_CUTOFF: a.babel_type = "C2"
                      else: a.babel_type = "C3"

              if a.babel_type[0] == 'N':
                   if a.babel_type=="N":
                      if k.babel_type[:2]=='C1' and \
                         bond_length <= V1_N1_C1_CUTOFF:
                          a.babel_type = "N1"
                      elif (k.babel_type[:2] == "C2" or \
                            k.babel_type[:2] == "C3") \
                           and bond_length > V1_N3_C_CUTOFF:
                          a.babel_type = "N3"
                      elif a.babel_type[:2]== "N3" and \
                           bond_length > V1_N3_N3_CUTOFF:
                          a.babel_type = "N3"
                      elif a.babel_type[:2]== "Npl" and \
                           bond_length > V1_N3_N2_CUTOFF:
                          a.babel_type = "N3"
                      else:
                          a.babel_type = "Npl"

              if a.babel_type[0] == 'O':
                  if a.babel_type=="O":
                      if k.babel_type in ["Cac", "Pac", "Sac", "Ntr"]:
                          a.babel_type = "O-"
                      elif k.babel_type in ["Nox", "Pox", "Sox"]:
                          a.babel_type = "O2"
                      elif k.babel_type[0] =="C" and \
                           bond_length <= V1_O2_C2_CUTOFF:
                          a.babel_type = "O2"
                          k.babel_type = "C2"
                          k._babel_redo = 0
                      elif k.babel_type=="As" and \
                           bond_length <= V1_O2_AS_CUTOFF:
                          a.babel_type = "O2"
                      else: a.babel_type = "O3"

              if a.babel_type[0] == 'S':
                  if a.babel_type=="S":
                      if k.babel_type[0] =="P": a.babel_type = "S2"
                      elif k.babel_type[0]=="C" and \
                           bond_length <= V1_S2_C2_CUTOFF:
                          a.babel_type = "S2"
                          k.babel_type = "C2"
                          a._babel_redo = 0
                      elif k.babel_type=="As" and \
                           bond_length <= V1_S2_AS_CUTOFF:
                          a.babel_type = "S2"
                      else: a.babel_type = "S3"


    def phase4(self):

      for a in self.atoms:
          if a._babel_redo==1:
              for b in a.bonds:
                  j = b.atom1
                  if j==a: j = b.atom2
                  bond_length = distance(a.coords, j.coords)
                  if bond_length <= V2_C2_C_CUTOFF and j.babel_type[0] == 'C' and j.element == 'C':
                      a.babel_type = "C2"
                  elif bond_length <= V2_C2_N_CUTOFF and j.babel_type[0] =='N':
                      a.babel_type = "C2"

              for b in a.bonds:
                  j = b.atom1
                  if j==a: j = b.atom2
                  bond_length = distance(a.coords, j.coords)
                  if bond_length > V2_C3_C_CUTOFF and j.babel_type[0] == 'C' and j.element == 'C':
                      a.babel_type = "C3"
                  elif bond_length > V2_C3_N_CUTOFF and j.babel_type[0] == 'N':
                      a.babel_type = "C3"
                  elif bond_length > V2_C3_O_CUTOFF and j.babel_type[0] == 'O':
                      a.babel_type = "C3"

          elif a._babel_redo==2:
              for b in a.bonds:
                  j = b.atom1
                  if j==a: j = b.atom2
                  bond_length = distance(a.coords, j.coords)
                  if bond_length <= V2_N2_C_CUTOFF and j.babel_type[0] == 'C' and j.element == 'C':
                      a.babel_type = "Npl"
                  elif bond_length <= V2_N2_N_CUTOFF and j.babel_type[0] =='N':
                      a.babel_type = "Npl"

          elif a._babel_redo==3:
              flag = 0
              for b in a.bonds:
                  j = b.atom1
                  if j==a: j = b.atom2
                  bond_length = distance(a.coords, j.coords)

                  if bond_length <= V2_C2_C_CUTOFF and j.babel_type[0] == 'C' and j.element == 'C':
                      a.babel_type = "C2"
                      flag = 1
                  elif bond_length <= V2_C2_N_CUTOFF and j.babel_type[0] =='N':
                      a.babel_type = "C2"
                      flag = 1

              if flag == 0:
                  for b in a.bonds:
                      j = b.atom1
                      if j==a: j = b.atom2
                      bond_length = distance(a.coords, j.coords)
                      if bond_length > V2_C3_C_CUTOFF and j.babel_type[0]=='C' and j.element == 'C':
                          a.babel_type = "C3"
                          flag = 1
                      elif bond_length>V2_C3_N_CUTOFF and j.babel_type[0]=='N':
                          a.babel_type = "C3"
                          flag = 1
                      elif bond_length>V2_C3_O_CUTOFF and j.babel_type[0]=='O':
                          a.babel_type = "C3"
                          flag = 1
                      elif flag == 0:
                          if bond_length > GEN_C3_C_CUTOFF and \
                             j.babel_type[0] == 'C' and j.element == 'C':
                              a.babel_type = "C3"
                              flag = 1

    def phase5(self):

        for a in self.atoms:
            if a.babel_type == "C2":
                flag = 0;
                for b in a.bonds:
                    j = b.atom1
                    if j==a: j = b.atom2
                    if not j.babel_type in ["C3", "DC", "HC", "N3", "N3+", "O3"] and\
                       not j.babel_type in ["Pac", "Sac", "Sox", "C1", "S3", "Cac" ]:
                        flag = 1
                if flag == 0:
                    a.babel_type = "C3"


    def phase6(self):

        for a in self.atoms:
            no_plus = 1
            protonated = 1
            if a.babel_type=="N3":
                for b in a.bonds:
                    conn = b.atom1
                    if conn==a: conn = b.atom2
                    # If an unsaturated atom is attached to this nitrogen then
                    # it should be Npl
                    if len(a.bonds) == 2 and \
                       conn.babel_type in ["Car","C2","Sox","Sac","Pac","So2"]:
                        protonated = 0
                        a.babel_type = "Npl"

                    # If the attached atom is something other that C3,
                    # H or D the nitrogen is not positively charged
                    if conn.babel_type != "C3" and conn.babel_atomic_number != 1:
                        protonated = 0
                if protonated: a.babel_type = "N3+"

            # look for guanadinium nitrogens
            elif a.babel_type == "C2":
                # First see if we have an sp2 carbon surrounded by 3 sp2
                # nitrogens

                m = 0;
                for b in a.bonds:
                    k = b.atom1
                    if k==a: k = b.atom2
                    if k.babel_type=="Npl" or k.babel_type=="N2" or k.babel_type=="Ng+": m=m+1

                if m == 3:
                    a.babel_type = "C+"
                    for b in a.bonds:
                        k = b.atom1
                        if k==a: k = b.atom2
                        k.babel_type = "Ng+"

            elif a.babel_type == "Cac":
                for b in a.bonds:
                    k = b.atom1
                    if k==a: k = b.atom2
                    if k.babel_type[0]=="O" and self.count_heavy_atoms(k) == 1:
                        k.babel_type = "O-"


    def check_for_carbonyl(self, atom):
        for b in atom.bonds:
            bonded_atom = b.atom1
            if bonded_atom==atom: bonded_atom = b.atom2
            if bonded_atom.babel_type=="O2" or bonded_atom.babel_type=="S2":
                return 3
        return 2


    def check_for_amides(self):
  
        for a in self.atoms:
            if a.babel_type=="Npl":
                for b in a.bonds:
                    conn = b.atom1
                    if conn==a: conn = b.atom2
                    if conn.babel_type=="Cac" or conn.babel_type=="Sox" or \
                       conn.babel_type=="So2":
                        a.babel_type = "Nam"
                        break

                    if conn.babel_type=="C2":
                        if self.check_for_carbonyl(conn) == 3:
                            a.babel_type = "Nam"
                            break


    def getProperty(self, property, elements):
        """list <- getProperty(property, elements)

        property has to be 'num', 'cov_rad', 'bond_ord_rad', 'vdw_rad',
                           'bs_rad', 'max_bonds' or 'rgb'
        elements is a list of 1 or 2 character(s) strings
        """
        if property not in babel_elements["C"].keys():
            raise RuntimeError("Invalid property %s, has to be in %s\n" % \
                               (property, babel_elements["C"].keys()))
        prop = []
        for el in elements:
            prop.append(babel_elements[el][property])
        return prop

    
                               
class TypeConverter:

    def __init__(self, outputType):
        if outputType not in babel_types.keys():
            raise RuntimeError("Invalid format %s\n"%outputType)
        self.outputType = outputType


    def convert(self, input, mode='all_caps'):

        try:
            i = babel_types['INT'].index(input)
            return babel_types[self.outputType][i]
        except:
            print "Unable to assign %s type to atom %s"%(self.outputType,input)
            if mode=='zero':
                return 0
            elif mode=='dummy':
                i = babel_types['INT'].index("X")
                return babel_types[self.outputType][i]
            elif mode=='all_caps':
                return string.upper(input)
            else: return input

    def clean_atom_type(self, type_name):
       name = string.upper(type_name[0])
       if len(type_name) > 1:
           name = name + string.lower(type_name[1])
           if name[1] not in string.letters:
               return name[0]
       return name


if __name__ == '__main__':
    import pdb, sys
    from cycle import RingFinder
    from bo import BondOrder
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

