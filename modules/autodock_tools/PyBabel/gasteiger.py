#############################################################################
#
# Authors: Michel F. SANNER, Ruth HUEY and Garrett M. MORRIS
# Reimplemented from Babel v1.6 from Pat Walters and Math Stahl
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################
#
# $Header: /opt/cvs/python/packages/share1.5/PyBabel/gasteiger.py,v 1.7 2002/09/24 17:32:00 rhuey Exp $
#
# $Id: gasteiger.py,v 1.7 2002/09/24 17:32:00 rhuey Exp $
#
#

"""
This file implements the gasteiger class.

Before this gasteiger object can be used, atoms must have been assigned
a type see (AtomHybridization in types.py).

reimplmentation of Babel1.6 in Python by Michel Sanner April 2000
Original code by W. Patrick Walters and Matthew T. Stahl 
"""



import math
from atomTypes import TypeConverter

###Partial equilization of orbital electronegativity, PEOP
###PEOP method to calculate a, b and c
###  from  J. Gasteiger and M. Marsili, Tetrahedron 36, 3219 (1980)
###
###energy parameters
###
###for neutral atoms: 
###  from J. Hinze and H.H. Jaffe, J. Am. Chem. Soc. 84, 540 (1962)
###     I0 Valence State Ionization Potentials 
###     E0 Valence State Electron Affinities 
###
###for positive ions: 
###  from J. Hinze and H.H. Jaffe, J. Phys. Chem. 67, 1501 (1963)
###     I+ Valence State Ionization Potentials 
###     E+ Valence State Electron Affinities 
###
### orbital ids used:
###     te tetrahedral, sp3
###     tr trigonal, sp2
###     di digonal, sp

energyParms = {}

#elemDict['sp2'] = hybridD = {}
#elemDict['sp2'] = hybridD = {}
#hybridD['I0'] =
#hybridD['E0'] =
#hybridD['I+'] =
#hybridD['E+'] =

#Carbon
energyParms['C'] = {}
energyParms['C']['sp3'] = hybridD ={}     
hybridD['I0'] = 14.61                #tetetete
hybridD['E0'] =  1.34                #tetetete
hybridD['I+'] = 26.71                #tetete, sigma
hybridD['E+'] = 11.37                #tetete, sigma
energyParms['C']['sp2'] = hybridD ={}     
hybridD['I0'] = 15.62                #trtrtr pi, sigma
hybridD['E0'] =  1.95                #trtrtr pi, sigma
hybridD['I+'] = 27.36                #trtr pi, sigma
hybridD['E+'] = 11.91                #trtr pi, sigma
#hybridD['I+'] = 23.68                #trtr pi, pi
#hybridD['E+'] = 10.45                #trtr pi, pi
energyParms['C']['sp'] = hybridD ={}     
hybridD['I0'] = 17.42                #didi pi pi, sigma
hybridD['E0'] =  3.34                #didi pi pi, sigma
hybridD['I+'] = 28.16                #di pi pi, sigma
hybridD['E+'] = 12.96                #di pi pi, sigma
#hybridD['I+'] = 29.85                #didi pi, sigma
#hybridD['E+'] = 13.29                #didi pi, sigma
#hybridD['I+'] = 23.86                #didi pi, pi
#hybridD['E+'] =  9.83                #didi pi, pi

#Nitrogen
energyParms['N'] = {}
energyParms['N']['sp3'] = hybridD ={}     
hybridD['I0'] = 18.93                #te2tetete
hybridD['E0'] =  4.15                #te2tetete
hybridD['I+'] = 33.29                #tetetete, sigma
hybridD['E+'] = 14.14                #tetetete, sigma
energyParms['N']['sp2'] = hybridD ={}     
hybridD['I0'] = 20.60                #tr2trtr pi, sigma
hybridD['E0'] =  5.14                #tr2trtr pi, sigma
hybridD['I+'] = 34.62                #trtrtr pi, sigma
hybridD['E+'] = 15.09                #trtrtr pi, sigma
#hybridD['I+'] = 28.71                #trtrtr pi, pi
#hybridD['E+'] = 11.96                #trtrtr pi, pi
energyParms['N']['sp'] = hybridD ={}     
hybridD['I0'] = 23.91                #di2di pi pi, sigma
hybridD['E0'] =  7.45                #di2di pi pi, sigma
hybridD['I+'] = 37.00                #didi pi pi, sigma
hybridD['E+'] = 17.24                #didi pi pi, sigma
#hybridD['I+'] = 28.70                #didi pi pi, pi
#hybridD['E+'] = 12.06                #didi pi pi, pi

#Oxygen
energyParms['O'] = {}
energyParms['O']['te'] = hybridD ={}     
hybridD['I0'] = 24.39                #te2te2tete
hybridD['E0'] =  6.11                #te2te2tete
hybridD['I+'] = 40.31                #te2tetete, sigma
hybridD['E+'] = 18.70                #te2tetete, sigma
energyParms['O']['sppp'] = hybridD ={}     
hybridD['I0'] = 17.28                #s2 p2 p p
hybridD['E0'] =  2.01                #s2 p2 p p
hybridD['I+'] = 34.15                #s2 p p p, p 
hybridD['E+'] = 14.61                #s2 p p p, p
#hybridization adjusted to an angle of 106deg [methanol]
#106 = 81% of 109.5 [te] + 19% of 90deg [sppp]
energyParms['O']['sp3'] = hybridD ={}     
prop = 0.81
z = energyParms['O']['te']
y = energyParms['O']['sppp']
hybridD['I0'] = prop * z['I0'] + (1.0 - prop) * y['I0']                # methanol
hybridD['E0'] = prop * z['E0'] + (1.0 - prop) * y['E0']                # methanol
hybridD['I+'] = prop * z['I+'] + (1.0 - prop) * y['I+']                # methanol
hybridD['E+'] = prop * z['E+'] + (1.0 - prop) * y['E+']                # methanol
energyParms['O']['sp2'] = hybridD ={}     
hybridD['I0'] = 26.65                #tr2 tr2 tr pi, sigma
hybridD['E0'] =  7.49                #tr2 tr2 tr pi, sigma
hybridD['I+'] = 42.49                #tr2 tr tr pi, sigma
hybridD['E+'] = 20.15                #tr2 tr tr pi, sigma

## new stuff
#Beryllium (2)
energyParms['Be'] = elemDict = {}
elemDict['sp3'] = hybridD = {}
hybridD['I0'] =  7.18           #tete
hybridD['E0'] =  0.51           #tete
hybridD['I+'] = 18.21           #s
hybridD['E+'] =  9.32           #s
#hybridD['I+'] = 14.25           #p
#hybridD['E+'] =  5.32           #p
elemDict['sp2'] = hybridD = {}
hybridD['I0'] =  7.38           #tr pi, sigma
hybridD['E0'] =  0.63           #tr pi, sigma
hybridD['I+'] = 18.21           #s
hybridD['E+'] =  9.32           #s
#hybridD['I+'] = 14.25           #p
#hybridD['E+'] =  5.32           #p

#Sulfur (2)
energyParms['S'] = elemDict = {}
elemDict['sp3'] = hybridD = {}
hybridD['I0'] = 15.50           #te2te2tete
hybridD['E0'] =  4.77           #te2te2tete
hybridD['I+'] = 27.65           #te2tetete, sigma
hybridD['E+'] = 13.64           #te2tetete, sigma
elemDict['sp2'] = hybridD = {}
hybridD['I0'] = 16.27           #tr2 tr tr pi2
hybridD['E0'] =  5.49           #tr2 tr tr pi2
hybridD['I+'] = 28.51           #tr tr tr pi2, sigma
hybridD['E+'] = 14.33           #tr tr tr pi2, sigma

#Phosphorus (3)
energyParms['P'] = elemDict = {}
elemDict['sp3'] = hybridD = {}
hybridD['I0'] = 14.57           #te2tetete
hybridD['E0'] =  3.24           #te2tetete
hybridD['I+'] = 24.10           #tetetete, sigma
hybridD['E+'] = 12.09           #tetetete, sigma
elemDict['sp2'] = hybridD = {}
hybridD['I0'] = 15.59           #tr2 tr tr pi, sigma
hybridD['E0'] =  3.74           #tr2 tr tr pi, sigma
hybridD['I+'] = 25.14           #tr tr tr pi, sigma
hybridD['E+'] = 12.72           #tr tr tr pi, sigma

#Aluminium (3)
energyParms['Al'] = elemDict = {}
elemDict['sp3'] = hybridD = {}
hybridD['I0'] =  8.17           #tetete
hybridD['E0'] =  2.58           #tetete
hybridD['I+'] = 15.75           #tete, sigma
hybridD['E+'] =  6.64           #tete, sigma
elemDict['sp2'] = hybridD = {}
hybridD['I0'] =  8.65           #trtr pi, sigma
hybridD['E0'] =  2.94           #trtr pi, sigma
hybridD['I+'] = 16.28           #tr pi, sigma
hybridD['E+'] =  6.74           #tr pi, sigma

#Boron (3)
energyParms['B'] = elemDict = {}
elemDict['sp3'] = hybridD = {}
hybridD['I0'] = 10.43           #tetete
hybridD['E0'] =  1.53           #tetete
hybridD['I+'] = 20.93           #tete, sigma
hybridD['E+'] =  7.88           #tete, sigma
elemDict['sp2'] = hybridD = {}
hybridD['I0'] = 10.97           #trtr pi, sigma
hybridD['E0'] =  1.87           #trtr pi, sigma
hybridD['I+'] = 21.08           #tr pi, sigma
hybridD['E+'] =  8.02           #tr pi, sigma

#Silicon (4)
energyParms['Si'] = {}
energyParms['Si']['sp3'] = hybridD ={}     
hybridD['I0'] = 11.82                #tetetete
hybridD['E0'] =  2.78                #tetetete
hybridD['I+'] = 18.97                #tetete, sigma
hybridD['E+'] = 10.08                #tetete, sigma
energyParms['Si']['sp2'] = hybridD ={}     
hybridD['I0'] = 12.61                #trtrtr pi, sigma
hybridD['E0'] =  3.20                #trtrtr pi, sigma
hybridD['I+'] = 19.62                #trtr pi, sigma
hybridD['E+'] = 10.57                #trtr pi, sigma
#hybridD['I+'] = 16.53                #trtr pi, pi
#hybridD['E+'] = 9.54                #trtr pi, pi
energyParms['Si']['sp'] = hybridD ={}     
hybridD['I0'] = 14.06                #didi pi pi, sigma
hybridD['E0'] =  4.07                #didi pi pi, sigma
hybridD['I+'] = 20.62                #di pi pi, sigma
hybridD['E+'] = 11.56                #di pi pi, sigma
#hybridD['I0'] =  9.18                #didi pi pi, pi
#hybridD['E0'] =  2.20                #didi pi pi, pi
#hybridD['I+'] = 16.50                #didi pi, pi
#hybridD['E+'] =  8.60                #didi pi, pi

#Magnesium
energyParms['Mg'] = {}
energyParms['Mg']['sp3'] = hybridD ={}     
hybridD['I0'] =  6.28                #tete
hybridD['E0'] =  0.32                #tete
hybridD['I+'] = 15.03                #s, s
hybridD['E+'] =  7.64                #s, s
energyParms['Mg']['sp2'] = hybridD ={}     
hybridD['I0'] =  6.75                #tr pi, sigma
hybridD['E0'] =  0.38                #tr pi, sigma
hybridD['I+'] = 15.03                #s, s
hybridD['E+'] =  7.64                #s, s
energyParms['Mg']['sp'] = hybridD ={}     
hybridD['I0'] =  7.30                #di pi, sigma
hybridD['E0'] =  0.78                #di pi, sigma
hybridD['I+'] = 15.03                #s, s
hybridD['E+'] =  7.64                #s, s


class Gasteiger:
    """ """
    
    def __init__(self):
        """ """
        self.setup_sigma_params()


    def setup_sigma_params(self):
        """ """

        self.par = {
            "H" : [  7.17,  6.24, -0.56, 0.0 ],
            "C3": [  7.98,  9.18,  1.88, 0.0 ],
            "C2": [  8.79,  9.32,  1.51, 0.0 ],
            "C1": [ 10.39,  9.45,  0.73, 0.0 ],
            "N3": [ 11.54, 10.82,  1.36, 0.0 ],
            "N2": [ 12.87, 11.15,  0.85, 0.0 ],
            "N1": [ 15.68, 11.70, -0.27, 0.0 ],
            "O3": [ 14.18, 12.92,  1.39, 0.0 ],
            "O2": [ 17.07, 13.79,  0.47, 0.0 ],
            "F" : [ 14.66, 13.85,  2.31, 0.0 ],
            "Cl": [ 11.00,  9.69,  1.35, 0.0 ],
            "Br": [ 10.08,  8.47,  1.16, 0.0 ],
            "I" : [  9.90,  7.96,  0.96, 0.0 ],
            "S3": [ 10.14,  9.13,  1.38, 0.0 ]
            }

        #add values compute from energyParms information
        #NB: this overwrites specified entries BUT the numbers are
        #    virtually the same...
        for element, elementD in energyParms.items():
            for hybrid, hybridD in elementD.items():
                I0 = hybridD['I0']
                E0 = hybridD['E0']
                IPLUS = hybridD['I+']
                EPLUS = hybridD['E+']
                a = 0.5 * (I0 + E0)
                a = round(a, 3)
                b = 0.25 * (IPLUS + EPLUS - E0)
                b = round(b, 3)
                c = 0.25 * (IPLUS - 2 * I0 + EPLUS - E0)
                c = round(c, 3)
                if hybrid=='sp':
                    key = element + '1'
                elif hybrid not in ['sp2', 'sp3']:
                    continue
                else:
                    key = element + hybrid[-1]
                self.par[key] = [a, b, c, 0.0]
                
        for p in self.par.values():
            p[3] = p[0] + p[1] + p[2]


    def lookup_sigma_params(self, atoms):
        """ """
  
        converter = TypeConverter("MAP")
        for a in atoms:
            type = converter.convert( a.babel_type )
            if type[:2] in self.par.keys():
                a._gast_par = self.par[type]
            else:
                print "Sorry, there are no Gasteiger parameters available for atom %s"%a.full_name()
                a._gast_par = [ 0.0, 0.0, 0.0, 1.0 ]


    def calc_sigma_charges(self, atoms):
        """ """

        # set formal charges
        for a in atoms:
            a.gast_charge = 0.0
            if a.babel_atomic_number==6: # carbons
                if a.babel_type=="C+": a.gast_charge = 1.0
            elif a.babel_atomic_number==7: # nitrogens
                if len(a.bonds)==4: a.gast_charge = 1.0
                if a.babel_type == "N3+": a.gast_charge = 1.0
            elif a.babel_atomic_number==8: # oxygens
                if a.babel_type=="O-": a.gast_charge = -0.5
            a._xx = a._gast_par[0]

        z1 = 1.0
        cycle = 0

        while (1):

            z1 = z1 * 0.5
            d1 = 0.0

            for a in atoms:
                if a._gast_par[3] != 1.0:
                    q1 = a.gast_charge
                    for b in a.bonds:
                        a2 = b.atom1
                        if a2==a: a2 = b.atom2
                        if a2._gast_par[3] != 1.0:
                            hybridD = a2._gast_par[3]
                            if a2._xx > a._xx: hybridD = a._gast_par[3]
                            if a.babel_type == 'H': hybridD = 20.02
                            if a2.babel_type == 'H': hybridD = 20.02
                            a.gast_charge = a.gast_charge + \
                                            (a2._xx - a._xx)/hybridD*z1

                    q1 = math.fabs(a.gast_charge - q1)
                    if q1 > d1: d1 = q1

            if d1 >= 0.001:
                for a in atoms:
                    a._xx = a._gast_par[0] + a._gast_par[1] * a.gast_charge + \
                           a._gast_par[2] * a.gast_charge * a.gast_charge
            cycle = cycle+1
            if d1<=0.001 or cycle>5: break
 

    def compute(self, atoms):
        """compute(atoms) compute gasteiger charges for a set of atoms"""
        
        self.lookup_sigma_params(atoms)
        self.calc_sigma_charges(atoms)

        # cleanup
        for a in atoms:
            delattr(a,'_gast_par')
            delattr(a,'_xx')



if __name__=="__main__":
#    import pdb, sys
#    from MolKit.pdbParser import NewPdbParser

#    parser = NewPdbParser("/tsri/pdb/struct/%s.pdb"%sys.argv[1])
#    mols = parser.parse()
#    mol = mols[0]
#    mol.buildBondsByDistance()
#    allAtoms = mol.chains.residues.atoms

#    print "assigning atom types"
#    from atomTypes import AtomHybridization
#    babel = AtomHybridization()
#    babel.assignHybridization(allAtoms)

#    print "computing charges"
#    Gast = Gasteiger()
#    Gast.compute(allAtoms)

#    print allAtoms.babel_type[:20]
#    print allAtoms.gast_charge[:20]

#    f = open("result.charges", "w")
#    for i in range(len(allAtoms)):
#        a = allAtoms[i]
#        f.write("%2d %4s %10.4f\n"%(i+1, a.babel_type, a.gast_charge))
#    f.close()
    
    for element, elementD in energyParms.items():
        print 'Element = ', element
        for hybrid, hybridD in elementD.items():
            I0 = hybridD['I0']
            E0 = hybridD['E0']
            IPLUS = hybridD['I+']
            EPLUS = hybridD['E+']
            a = 0.5 * (I0 + E0)
            b = 0.25 * (IPLUS + EPLUS - E0)
            c = 0.25 * (IPLUS - 2 * I0 + EPLUS - E0)
            print '    ',hybrid, ' hybridization:'
            print '        a =', a,
            print 'b =', b,
            print 'c =', c
        print '\n',
            


