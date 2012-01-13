# http://users.unimi.it/~ddl/vega/manual/pages/appx_a.htm
#  Type        a        b        c        d
#  =============================================
#  H          7.17     6.24    -0.56    12.85
#  C3         7.98     9.18     1.88    19.04
#  C2         8.79     9.32     1.51    19.62
#  C1        10.39     9.45     0.73    20.57
#  N3        11.54    10.82     1.36    23.72
#  N2        12.87    11.15     0.85    24.87
#  N1        15.68    11.70    -0.27    27.11
#  O3        14.18    12.92     1.39    28.49
#  O2        17.07    13.79     0.47    31.33
#  F         14.66    13.85     2.31    30.82
#  Cl        11.00     9.69     1.35    22.04
#  Br        10.08     8.47     1.16    19.71
#  I          9.90     7.96     0.96    18.82
#  S3        10.14     9.13     1.38    20.65


ground_states = {\
				#      Eg     Ig
'H' : (0.747, 13.595),
'B' : (0.330,  8.296),   
'C' : (1.120, 11.256),
'N' : (1.050, 14.535),
'O' : (1.465, 13.614),
'F' : (3.480, 17.418),
'Si': (1.460,  8.149),
'P' : (0.770, 10.977),
'S' : (2.070, 10.357),
'Cl': (3.690, 12.974),
'Br': (3.550, 11.84),
'I' : (3.210, 10.45) 
} 


Chargeterms = {
'H' :   ( 7.17,       6.24, -0.56,  12.85),
# Carbon
'C.3':  ( 7.98,       9.18,  1.88,  19.04),
'C.cat':( 7.98,       9.18,  1.88,  19.04),
'C.2':  ( 8.79+0.5,   9.32,  1.51,  19.62),
'C.ar': ( 7.98+0.55,  9.18,  1.88,  19.04),
'C.1':  (10.39,       9.45,  0.73,  20.57),
# Nitrogen
'N.3':  (11.54+6.0,  10.28,  1.36,  28.00),
'N.4':  (11.54+6.0,  10.28,  1.36,  28.00),
'N.ar': (12.87-1.29, 11.15,  0.85,  24.87),
'N.2':  (12.87,      11.15,  0.85,  24.87),
'N.pl3':(12.87+0.5,  11.15,  0.85,  24.87),
'N.am': (12.87+3.5,  11.15,  0.85,  24.87),
'N.1':  (15.68,      11.70, -0.27,  27.11),
# Oxygen	
'O.OH': (14.18+0.8,  12.92,  1.39,  28.49),
'O.3':  (14.18-3.1,  12.92,  1.39,  28.49),
'O.2':  (14.18,      12.92,  1.39,  28.49),
'O.co2':(15.25,      13.79,  0.47,  31.33),
# Halogenides etc
'F' :   (12.36,      13.85,  2.31,  30.82),
'Cl':   ( 9.38+1.0,   9.69,  1.35,  22.04),
'Br':   (10.08+0.8,   8.47,  1.16,  19.71),
'I' :   ( 9.90+1.0,   7.96,  0.96,  18.82),
'S.3':  (10.13+0.5,   9.13,  1.38,  20.65),
'S.2':  (10.13+0.5,   9.13,  1.38,  20.65),
'S.o2': (10.13+0.5,   9.13,  1.38,  20.65),
'P.3':  (10.13+0.5,   9.13,  1.38,  20.65),
}


def PEOE( atoms, damp=0.778, k=1.56):
    def calcchi(atom, q):
        if abs(q) > 1.1:
            if q < 0.0: q = -1.1
            else:       q =  1.1
        if (q == 1.0) and (atom.sybylType == 'H'): 
            return 20.02
        else:
            if len(atom.abc) == 4:
                a,b,c,d = atom.abc
                return a + b*q + c*q*q + d*q*q*q
            else:
                a,b,c = atom.abc
                return a + b*q + c*q*q
    abs_qges = 0.0
    counter = 0
    for a in atoms.atoms:
        if not Chargeterms.has_key(a.sybylType):
            raise KeyError, 'PEOE Error: Atomtype <%s> not known, treating atom %s as dummy' % (a.sybylType, a.name)
        if a.sybylType == 'O.3':
            a.chi   = Chargeterms['O.OH'][0]
            a.abc   = Chargeterms['O.OH']
        else:
            a.chi   = Chargeterms[a.sybylType][0]
            a.abc   = Chargeterms[a.sybylType]
        if a.charge != 0.0:      
            a.formal_charge = a.charge*(1/k)
            abs_qges = abs_qges+abs(a.charge)
            counter = counter+1
        else:
            a.formal_charge = 0.0
        a.charge  = 0.0
        a.dq = 0.0
    if abs_qges != 0.0:
        cycles = 7
        for b in range(1,cycles):
            for i in atoms.atoms: # lAtoms
                i.chi = calcchi(i, i.charge)         
                i.dq = 0.0
                for j in i.intrabonds: ### lBondedAtoms
                    for yyy in atoms.atoms:
                        if yyy.name == j:
                            dchi  = (calcchi(yyy, yyy.charge) - i.chi)
                            if dchi > 0.0: i.dq += (dchi / calcchi(i, +1) * (damp**b))
                            else:          i.dq += (dchi / calcchi(yyy, +1) * (damp**b))
            for i in atoms.atoms: # lAtoms
                i.charge += i.dq+(0.166666667*i.formal_charge)  
        for i in atoms.atoms:  # lAtoms
            i.charge = i.charge * k
            i.charge = i.charge
            del i.dq
            del i.abc
        return atoms

    else:                          
    #
        for a in range(1,7):         
            for i in atoms.atoms: # lAtoms
                i.chi = calcchi(i, i.charge)
                i.dq  = 0.0
                for j in i.intrabonds: ### lBondedAtoms
                    for xxx in atoms.atoms:
                        if xxx.name == j:
                            dchi  = (calcchi(xxx, xxx.charge) - i.chi)
                            if dchi > 0.0: i.dq += (dchi / calcchi(i, +1) * (damp**a)) 
                            else:          i.dq += (dchi / calcchi(xxx, +1) * (damp**a))
            for i in atoms.atoms: # lAtoms
                i.charge += i.dq
        for i in atoms.atoms: #lAtoms
            i.charge = i.charge * k             
            i.charge = i.charge
            del i.dq
            del i.abc
        return atoms
