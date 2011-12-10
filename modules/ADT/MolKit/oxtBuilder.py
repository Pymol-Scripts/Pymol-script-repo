############################################################################
#
# Author:  Ruth Huey
#
# Copyright: M. Sanner TSRI 2004
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/MolKit/oxtBuilder.py,v 1.1 2005/09/15 16:02:04 rhuey Exp $
#
# $Id: oxtBuilder.py,v 1.1 2005/09/15 16:02:04 rhuey Exp $
#

"""
This module implements the OxtBuilder class which adds an oxygen atom to
a specific carbon atom, presumably in the c-terminus residue. This class
assumes that bonds have been previously built in the molecule.

"""

from MolKit.molecule import Atom, AtomSet, Bond
from PyBabel.atomTypes import AtomHybridization
from PyBabel.addh import AddHydrogens


class OxtBuilder:
    """Base Class for adding an oxygen atom to a carbon atom.
       NOTE: molecule must have bonds built between atoms
    """


    def __init__(self):
        self.addh = AddHydrogens()


    def add_oxt(self, catom):       
        if catom.element!='C':
            return
        mol = catom.top    
        ##check for bonds 
        #if len(mol.allAtoms.bonds[0])==0:
        #    mol.buildBondsByDistance()
        #check whether residue already has OXT
        res = catom.parent
        if 'OXT' in res.atoms.name:
            print 'not adding OXT to ',res.full_name(),'\n', 'it already has an OXT atom'
            return
        #check whether catom has a hydrogen to delete
        hatoms = catom.parent.atoms.get(lambda x: x.name=='HC')
        if len(hatoms):
            hatom = hatoms[0]
            #check for hbonds
            if hasattr(hatom, 'hbonds'):
                #remove hbonds
                for b in hatom.hbonds:
                    atList = [b.donAt, b.accAt]
                    if b.hAt is not None:
                        atList.append(b.hAt)
                    for at in atList:
                        #hbonds might already be gone
                        if not hasattr(at, 'hbonds'):
                            continue
                        okhbnds = []
                        for hb in at.hbonds:
                            if hb!=b:
                                okhbnds.append(hb)
                        if len(okhbnds):
                            at.hbonds = okhbnds
                        else:
                            delattr(at, 'hbonds')
            #remove covalent bonds                
            for b in hatom.bonds:
                at2 = b.atom1
                if at2 == hatom: at2 = b.atom2
                at2.bonds.remove(b)
            hatom.parent.remove(hatom, cleanup=1)

        #have to type atoms before call to add_sp2_hydrogen:
        if not hasattr(catom,'babel_type'):
            print 'catom has no babel_type: calling typeAtoms'
            #self.warningMsg(msg)
            #typeAtoms does whole molecule
            babel = AtomHybridization()
            babel.assignHybridization(mol.allAtoms)

        #NB: bond_length 1.28 measured from OXT-C bond in 1crn
        tup1 = self.addh.add_sp2_hydrogen(catom, 1.28)
        res = catom.parent

        # find where to insert H atom
        childIndex = res.children.index(catom)+1
        name = 'OXT'

        # create the OXT atom object
        atom = Atom(name, res, top=mol,
                    childIndex=childIndex, assignUniqIndex=0)

        # set atoms attributes
        atom._coords = [ tup1[0][0] ]
        if hasattr(catom, 'segID'): atom.segID = catom.segID
        atom.hetatm = 0
        atom.alternate = []
        atom.element = 'O'
        atom.occupancy = 1.0
        atom.conformation = 0
        atom.temperatureFactor = 0.0
        atom.babel_atomic_number = 8
        atom.babel_type = 'O-'
        atom.babel_organic = 1

        # create the Bond object bonding Hatom to heavyAtom
        bond = Bond( catom, atom, bondOrder=2)

        # create the color entries for all geometries
        # available for the other oxygen atom attached to 'C'
        oatom = res.atoms.get(lambda x: x.name=='O')[0]
        if oatom is not None:
            for key, value in oatom.colors.items():
                atom.colors[key] = value
                #atom.opacities[key] = oatom.opacities[key]
                
                # update the allAtoms set in the molecule
        mol.allAtoms = mol.chains.residues.atoms

        # update numbers of allAtoms
        fst = mol.allAtoms[0].number
        mol.allAtoms.number = range(fst, len(mol.allAtoms)+fst)

        # update _uniqIndex of this residues atoms
        res.assignUniqIndex()
        #return AtomSet([atom])
        return atom

