#############################################################################
#
# Author: Michel F. SANNER
# Reimplemented from Babel v1.6 from Pat Walters and Math Stahl
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################
#
# $Header: /opt/cvs/python/packages/share1.5/PyBabel/__init__.py,v 1.10 2006/03/15 00:45:52 annao Exp $
#
# $Id: __init__.py,v 1.10 2006/03/15 00:45:52 annao Exp $
#
#

"""

PyBabel v0.1alpha

    This a Python re-implmentation of some of Babel v1.6 functionalities, including:

    - Atom hybridization assignment (from geometry)
    - Gasteiger charges calculations
    - atom type conversions
    - properties from atom element name
    - finding rings
    - assigning bond orders
    - identifying aromatic rings
    - adding hydrogen atoms




I - Atom hybridization assignment:

    The AtomHybridization class allows to assign atom types.
  
  example:
    
    >>> #create an instance of AtomHybridization
    >>> from PyBabel.atomTypes import AtomHybridization
    >>> babel = AtomHybridization()
    >>> babel.assignHybridization(allAtoms)

    allAtoms is a list of 'Atom' objects having the following members:

        - a.element : atom's chemical element symbol (string,
                      2 character max, second character lower case)

        - a.coords  : 3-sequence of floats (tuple, float or array)
  
        - a.bonds   : list of Bond objects

    Bond are objects having the following members:

        - b.atom1 : instance of Atom

        - b.atom2 : instance of Atom


  After completion each atom has the following additional members:

      - babel_type          : string
      - babel_atomic_number : int
      - babel_organic       : 0/1




II - Gasteiger charges calculations:

  Before Gasteiger charges can be calculated, babel AtomHybridization as to
  be assigned (see I).

  example:
    
      >>> #create an instance of AtomHybridization
      >>> from PyBabel.atomTypes import AtomHybridization
      >>> babel = AtomHybridization()
      >>> babel.assignHybridization(allAtoms)
      >>>
      >>> from PyBabel.gasteiger import Gasteiger
      >>> Gast = Gasteiger()
      >>> Gast.compute(allAtoms)


  allAtoms is a list of 'Atom' objects having the following members:

      * Atom:

          - a.babel_type

          - a.babel_atomic_number

          - a.bonds

  After completion each atom has a gast_charge member




III - types conversion

  Before atom types can be obtained, babel's internal types have to be assigned
  (see I)

  example:
    
      >>> #create an instance of AtomHybridization
      >>> from PyBabel.atomTypes import AtomHybridization
      >>> babel = AtomHybridization()
      >>> babel.assignHybridization(allAtoms)
      >>>
      >>> from PyBabel.atomTypes import AtomHybridization
      >>> converter = TypeConverter('SYB')
      >>> newTypes = map( converter.convert, allAtoms.babel_type)




IV - properties from atom element name

  Before atom properties can be obtained, babel's internal types have to 
  be assigned (see I)

  example:
    
      >>> #create an instance of AtomHybridization
      >>> from PyBabel.atomTypes import AtomHybridization
      >>> babel = AtomHybridization()
      >>> babel.assignHybridization(allAtoms)
      >>>
      >>> rad = babel.getProperty('vdw_rad', allAtoms.element)

  valid property names are:

      - 'num'

      - 'cov_rad'

      - 'bond_ord_rad'

      - 'vdw_rad'

      - 'bs_rad'

      - 'max_bonds'

      - 'rgb'




V - finding rings:

  This algorithm is different from the one used in Babel to avoid the N^2
  complexity in number of bonds.

  When rings are nested the smalest rings are reported.

  The algorithms might fail to report all rings if all the atoms in some
  rings belong to multiple rings.

  example:
    
      >>> #create an instance of RingFinder object
      >>> from PyBabel.cycle import RingFinder
      >>> r = RingFinder()
      >>> r.findRings(atoms, bonds)
      >>> r.printRings()

  atoms has to be a list of Atom objects

      * Atom:
          
          - a.coords : 3-sequence of floats
          
          - a.bonds : list of Bond objects
          
      * Bond:
          
          - b.atom1 : instance of Atom
          
          - b.atom2 : instance of Atom
          

  After completion the RingFinder object has the following members:
      
      - ringCount: number of rings
      
      - rings : a list of rings. A ring is a dictionary with 2 keys
              'atoms' and 'bonds'
              
      - allRingAtoms: list of atoms that are in rings(atoms may appear twice)

      - allRingBonds: list of bonds that are in rings(bonds may appear twice)


  In addition:

      atoms involved in rings have a member 'ring' that is a list of
      rings they belong to (0-based list of integers)




VI - assigning bond orders:

  The BondOrder class that can be used to compute bond order.

  Before a BondOrder object can be used, atoms must have been assigned
  a type (see I)

  Bond order can be calculated using 2 different methods depending on whether
  rings have been identified previously or not. Babel decides to use the first
  method for molecules with more than 200 atoms and the second one else.

  example:
    
      >>> #create an instance of AtomHybridization
      >>> from PyBabel.atomTypes import AtomHybridization
      >>> atype = AtomHybridization()
      >>> atype.assignHybridization(atoms)

      >>> #create an instance of BondOrder object
      >>> from PyBabel.bo import BondOrder
      >>> bo = BondOrder()
      >>> bo.assignBondOrder( atoms, bonds )

      or

      >>> #create an instance of AtomHybridization
      >>> from PyBabel.atomTypes import AtomHybridization
      >>> atype = AtomHybridization()
      >>> atype.assignHybridization(atoms)

      >>> from PyBabel.cycle import RingFinder
      >>> rings = RingFinder()
      >>> rings.findRings(allAtoms, bonds)

      >>> from PyBabel.bo import BondOrder
      >>> bo = BondOrder()
      >>> bo.assignBondOrder( atoms, bonds, rings )

  
  atoms has to be a list of atom objects

      * Atom:
          
          - a.coords : 3-sequence of floats
          
          - a.bonds : list of Bond objects
          
          - babel_type: string
          
          - babel_atomic_number: int

      * Bond:
          
          - b.atom1 : instance of Atom
          
          - b.atom2 : instance of Atom
          

  After completion each bond has a 'bondOrder' attribute (integer)




VII - identifying aromatic rings

  Before this Aromatic object can be used to identify aromatic rings,
  
        1. atoms must have been assigned a type (see I)

	2. rings must have been identified (see V)

	3. bond orders must have been defined (see VI)

  example:
    
      >>> #create an instance of AtomHybridization
      >>> from PyBabel.atomTypes import AtomHybridization
      >>> atype = AtomHybridization()
      >>> atype.assignHybridization(atoms)

      >>> from PyBabel.cycle import RingFinder
      >>> rings = RingFinder()
      >>> rings.findRings(atoms, bonds)

      >>> from PyBabel.bo import BondOrder
      >>> bo = BondOrder()
      >>> bo.assignBondOrder(atoms, bonds, rings)

      >>> from PyBabel.aromatic import Aromatic
      >>> arom = Aromatic(rings)
      >>> arom.find_aromatic_atoms(atoms)

  atoms has to be a list of atom objects

      * Atom:

          - a.coords : 3-sequence of floats

          - a.bonds : list of Bond objects

          - babel_type: string

          - babel_atomic_number: int

      * Bond:

          - b.atom1 : instance of Atom

          - b.atom2 : instance of Atom

  After completion the bondOrder member of bonds in aromatic bonds is
  is 'aromatic'. 

  The aromatic rings has a new key called 'aromatic' set to 1 or 0
	



VIII - adding hydrogen atoms

  Before this AddHydrogens object can be used, atoms must have been assigned
  a type (see I)

  Hydrogen atoms can be added using 2 different methods. The first one 
  requires bondOrders to have been calculated previousely.

  example:
    
      >>> #create an instance of AtomHybridization
      >>> from PyBabel.atomTypes import AtomHybridization
      >>> atype = AtomHybridization()
      >>> atype.assignHybridization(atoms)

      >>> from PyBabel.addh import AddHydrogens
      >>> addh = AddHydrogens()
      >>> hat = addh.addHydrogens(atoms)

  atoms has to be a list of atom objects:

      * Atom:
          
          - a.coords : 3-sequence of floats

          - a.bonds : list of Bond objects

          - babel_type: string

          - babel_atomic_number: int

      * Bond:
          
          - b.atom1 : instance of Atom

          - b.atom2 : instance of Atom

  example:
  
      >>> from PyBabel.addh import AddHydrogens
      >>> addh = AddHydrogens()
      >>> hat = addh.addHydrogens(atoms, method='withBondOrder')

  atoms has to be a list of atom objects as above and

       * Bond:
           
           - b.atom1 : instance of Atom

           - b.atom2 : instance of Atom

           - b.bondOrder : integer

          
  These calls return a list 'hat' containing a tuple for each Hydrogen
  atom to be added. This tuple provides:

      - coordsH       : 3-float coordinates of the H atom
      
      - atom          : the atom to which the H atom is to be connected
      
      - atomic_number : the babel_atomic_number of the H atom
      
      - type          : tyhe babel_type of the H atom





Michel Sanner April 2000
sanner@scripps.edu
"""

CRITICAL_DEPENDENCIES = ['MolKit', 'mglutil']
NONCRITICAL_DEPENDENCIES = ['Pmv', 'DejaVu']
