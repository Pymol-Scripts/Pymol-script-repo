#
# Last modified on Thu Jul 17 10:27:13 PDT 2003 by lindy
#
# $Id: MolecularDescriptor.py,v 1.3 2003/07/17 17:27:07 lindy Exp $
#



class MolecularDescriptorMixIn:
    """Compute descriptors for mixed in class.

    @@ How can this class be made less dependent on the class it's
       mixed-in with ?? @@

    Assumptions about instances of the main class:
    1. attribute allAtoms is a list of atoms with
       the following attributes:
       1.1. element -  used as key to atomic_weight dictionary
       1.2. bonds - used by Hbond_donors
       1.2.1 calls molecule.buildBondsByDistance() if necessary
       1.2.2 bonds have neighborAtom() method

    Descriptors:
        atomic_formula
        molecular_weight:
        Hbond_donors:
        Hbond_acceptors:
        rotabable_bonds:
        rigid_rings:
        logP:                    NotImplemented
    """

    def atomic_formula(self):
        """Return a {'element': count} dictionary
        """
        formula = {}
        for atom in self.allAtoms:
            try:
                formula[atom.element] += 1
            except KeyError:
                formula[atom.element] = 1
        return formula

    
    # Coplen, TB (1999) Atomic Weights of the Elements 1999.
    # US Geological Survey, Reston, Virginia, USA
    # http://physics.nist.gov/PhysRefData/Compositions/index.html
    atomic_weight = {
        'H'  :  1.00794,
        'He' :  4.002602,
        'Li' :  6.941,
        'Be' :  9.012182,
        'B'  :  10.811,
        'C'  :  12.0107,
        'N'  :  14.0067,
        'O'  :  15.9994,
        'F'  :  18.9984032,
        'Ne' :  20.1797,
        'Na' :  22.989770,
        'Mg' :  24.3050,
        'Al' :  26.981538,
        'Si' :  28.0855,
        'P'  :  30.973761,
        'S'  :  32.065,
        'Cl' :  35.453,
        'Ar' :  39.948,
        'K'  :  39.0983,
        'Ca' :  40.078,
        'Sc' :  44.955910,
        'Ti' :  47.867,
        'V'  :  50.9415,
        'Cr' :  51.9961,
        'Mn' :  54.938049,
        'Fe' :  55.845,
        'Co' :  58.933200,
        'Ni' :  58.6934,
        'Cu' :  63.546,
        'Zn' :  65.39,
        'Ga' :  69.723,
        'Ge' :  72.64,
        'As' :  74.92160,
        'Se' :  78.96,
        'Br' :  79.904,
        'Kr' :  83.80,
        'Rb' :  85.4678,
        'Sr' :  87.62,
        'Y'  :  88.90585,
        'Zr' :  91.224,
        'Nb' :  92.90638,
        'Mo' :  95.94,
        'Ru' :  101.07,
        'Rh' :  102.90550,
        'Pd' :  106.42,
        'Ag' :  107.8682,
        'Cd' :  112.411,
        'In' :  114.818,
        'Sn' :  118.710,
        'Sb' :  121.760,
        'Te' :  127.60,
        'I'  :  126.90447,
        'Xe' :  131.293,
        'Cs' :  132.90545,
        'Ba' :  137.327,
        'La' :  138.9055,
        'Ce' :  140.116,
        'Pr' :  140.90765,
        'Nd' :  144.24,
        'Sm' :  150.36,
        'Eu' :  151.964,
        'Gd' :  157.25,
        'Tb' :  158.92534,
        'Dy' :  162.50,
        'Ho' :  164.93032,
        'Er' :  167.259,
        'Tm' :  168.93421,
        'Yb' :  173.04,
        'Lu' :  174.967,
        'Hf' :  178.49,
        'Ta' :  180.9479,
        'W'  :  183.84,
        'Re' :  186.207,
        'Os' :  190.23,
        'Ir' :  192.217,
        'Pt' :  195.078,
        'Au' :  196.96655,
        'Hg' :  200.59,
        'Tl' :  204.3833,
        'Pb' :  207.2,
        'Bi' :  208.98038
        } # atomic_weight


    def molecular_weight(self):
        """Return the molecular weight of the molecule.
        """
        mw = 0.0
        for a in self.allAtoms:
            try:
                mw += self.atomic_weight[a.element]
            except KeyError, key:
                print "Unknown element: %s" % (key)
        return mw


    def Hbond_donors(self):
        """Return the number of Hydrogen bond donors in the molecule.
        
        A Hydrogen bond donor is a H bonded to an O, or N.
        """
        if not self.hasBonds:
            self.buildBondsByDistance()
        num_donors = 0
        for a in self.allAtoms:
            if a.element == 'H':
                num_donors += a.bonds[0].neighborAtom(a).element in ('O', 'N')
        return num_donors


    def Hbond_acceptors(self):
        """Return the number of Hydrogen bond acceptors in the molecule.
        
        A Hydrogen bond acceptor is an O or and N.
        """
        num_acceptors = 0
        for a in self.allAtoms:
            num_acceptors += a.element in ('O', 'N')
        return num_acceptors


    def rotatable_bonds(self):
        """Return the number of rotatable bonds in the molecule..
        """
        raise NotImplementedError


    def rigid_rings(self):
        """Return the number of rings in the molecule.
        """
        raise NotImplementedError


    def logP(self):
        """Return the octanol partition coefficient of the molecule.
        """
        raise NotImplementedError


    def mix_in(self, object):
        """Mix this class into the base classes of the given object
        """
        object.__classes.__bases__ += (MolecularDescriptorMixIn,)
        # check to see if you're already there !!

# MolecularDescriptorMixIn

if __name__ == '__main__':
    # read the molecule
    from MolKit import Read
    mol = Read("nfv.pdb")[0]

    # mix in the molecular descriptors
    mol.__class__.__bases__ += (MolecularDescriptorMixIn,)

    print "atomic formula of %s: %s" % (mol.name, repr(mol.atomic_formula()) )
    print "molecular weight of %s: %4.2f" % (mol.name, mol.molecular_weight())
    print "H-bond acceptors of %s: %d" % (mol.name, mol.Hbond_acceptors())
    print "H-bond donors    of %s: %d" % (mol.name, mol.Hbond_donors())
