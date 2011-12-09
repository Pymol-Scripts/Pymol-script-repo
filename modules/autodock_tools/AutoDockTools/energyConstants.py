
#############################################################################
#
# Author: Ruth HUEY, William LINDSTROM
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/energyConstants.py,v 1.4 2002/03/12 21:24:55 rhuey Exp $
#
#
# $Id: energyConstants.py,v 1.4 2002/03/12 21:24:55 rhuey Exp $
#
#


#FROM constants.py
# from: gpf3gen 
# Free energy model 140n s coefficients:
#
FE_vdW_coeff   = 0.1485
FE_estat_coeff = 0.1146
FE_hbond_coeff = 0.0656
FE_tors_coeff  = 0.3113
FE_desol_coeff = 0.1711

# AutoDock Binding Free Energy Model 140n
# Uses the C.A.Lig.only solvation model, i.e.
# aliphatic (C) and aromatic (A) carbon atoms only
#
SolVol={}
SolVol[ "C" ] = 12.77
SolVol[ "A" ] = 10.80
SolVol[ "N" ] =  0.
SolVol[ "n" ] =  0.
SolVol[ "O" ] =  0.
SolVol[ "H" ] =  0.

#SolPar[ "C" ] =  0.6844
#SolPar[ "A" ] =  0.1027

SolPar={}
SolPar[ "C" ] =  4.0 * FE_desol_coeff
SolPar[ "A" ] =  0.6 * FE_desol_coeff
SolPar[ "N" ] =  0.
SolPar[ "n" ] =  0.
SolPar[ "O" ] =  0.
SolPar[ "H" ] =  0.

SolCon={}
SolCon[ "C" ] =  0.
SolCon[ "A" ] =  0.
SolCon[ "N" ] =  0.
SolCon[ "n" ] =  0.
SolCon[ "O" ] =  0.236
SolCon[ "H" ] =  0.118

Rij={}
Rij[ "lj" "C" "C" ] = 4.00
Rij[ "lj" "C" "N" ] = 3.75
Rij[ "lj" "C" "n" ] = 3.75
Rij[ "lj" "C" "O" ] = 3.60
Rij[ "lj" "C" "P" ] = 4.10
Rij[ "lj" "C" "S" ] = 4.00
Rij[ "lj" "C" "H" ] = 3.00
Rij[ "lj" "C" "f" ] = 2.65
Rij[ "lj" "C" "F" ] = 3.54
Rij[ "lj" "C" "c" ] = 4.04
Rij[ "lj" "C" "b" ] = 4.17
Rij[ "lj" "C" "I" ] = 4.36
#Rij[ "lj" "C" "M" ] = 2.79
Rij[ "lj" "C" "M" ] = 2.65  #Mg
Rij[ "lj" "C" "Z" ] = 2.74  #Zn
Rij[ "lj" "C" "L" ] = 2.99  #Ca
Rij[ "lj" "C" "X" ] = 3.00  #Xx

Rij[ "lj" "N" "C" ] = 3.75
Rij[ "lj" "N" "N" ] = 3.50
Rij[ "lj" "N" "n" ] = 3.50
Rij[ "lj" "N" "O" ] = 3.35
Rij[ "lj" "N" "P" ] = 3.85
Rij[ "lj" "N" "S" ] = 3.75
Rij[ "lj" "N" "H" ] = 2.75
Rij[ "lj" "N" "f" ] = 2.40
Rij[ "lj" "N" "F" ] = 3.29
Rij[ "lj" "N" "c" ] = 3.79
Rij[ "lj" "N" "b" ] = 3.92
Rij[ "lj" "N" "I" ] = 4.11
#Rij[ "lj" "N" "M" ] = 2.54
Rij[ "lj" "N" "M" ] = 2.40  #Mg
Rij[ "lj" "N" "Z" ] = 2.49  #Zn
Rij[ "lj" "N" "L" ] = 2.74  #Ca
Rij[ "lj" "N" "X" ] = 2.75  #Xx

Rij[ "lj" "n" "C" ] = 3.75
Rij[ "lj" "n" "N" ] = 3.50
Rij[ "lj" "n" "n" ] = 3.50
Rij[ "lj" "n" "O" ] = 3.35
Rij[ "lj" "n" "P" ] = 3.85
Rij[ "lj" "n" "S" ] = 3.75
Rij[ "lj" "n" "H" ] = 2.75
Rij[ "lj" "n" "f" ] = 2.40
Rij[ "lj" "n" "F" ] = 3.29
Rij[ "lj" "n" "c" ] = 3.79
Rij[ "lj" "n" "b" ] = 3.92
Rij[ "lj" "n" "I" ] = 4.11
#Rij[ "lj" "n" "M" ] = 2.54
Rij[ "lj" "n" "M" ] = 2.40  #Mg
Rij[ "lj" "n" "Z" ] = 2.49  #Zn
Rij[ "lj" "n" "L" ] = 2.74  #Ca
Rij[ "lj" "n" "X" ] = 2.75  #Xx

Rij[ "lj" "O" "C" ] = 3.60
Rij[ "lj" "O" "N" ] = 3.35
Rij[ "lj" "O" "n" ] = 3.35
Rij[ "lj" "O" "O" ] = 3.20
Rij[ "lj" "O" "P" ] = 3.70
Rij[ "lj" "O" "S" ] = 3.60
Rij[ "lj" "O" "H" ] = 2.60
Rij[ "lj" "O" "f" ] = 2.25
Rij[ "lj" "O" "F" ] = 3.15
Rij[ "lj" "O" "c" ] = 3.65
Rij[ "lj" "O" "b" ] = 3.77
Rij[ "lj" "O" "I" ] = 3.96
#Rij[ "lj" "O" "M" ] = 2.39
Rij[ "lj" "O" "M" ] = 2.25  #Mg
Rij[ "lj" "O" "Z" ] = 2.34  #Zn
Rij[ "lj" "O" "L" ] = 2.59  #Ca
Rij[ "lj" "O" "X" ] = 2.60  #Xx

Rij[ "lj" "P" "C" ] = 4.10
Rij[ "lj" "P" "N" ] = 3.85
Rij[ "lj" "P" "n" ] = 3.85
Rij[ "lj" "P" "O" ] = 3.70
Rij[ "lj" "P" "P" ] = 4.20
Rij[ "lj" "P" "S" ] = 4.10
Rij[ "lj" "P" "H" ] = 3.10
Rij[ "lj" "P" "f" ] = 2.75
Rij[ "lj" "P" "F" ] = 3.65
Rij[ "lj" "P" "c" ] = 4.14
Rij[ "lj" "P" "b" ] = 4.27
Rij[ "lj" "P" "I" ] = 4.46
#Rij[ "lj" "P" "M" ] = 2.89
Rij[ "lj" "P" "M" ] = 2.75  #Mg
Rij[ "lj" "P" "Z" ] = 2.84  #Zn
Rij[ "lj" "P" "L" ] = 3.09  #Ca
Rij[ "lj" "P" "X" ] = 3.10  #Xx

Rij[ "lj" "S" "C" ] = 4.00
Rij[ "lj" "S" "N" ] = 3.75
Rij[ "lj" "S" "n" ] = 3.75
Rij[ "lj" "S" "O" ] = 3.60
Rij[ "lj" "S" "P" ] = 4.10
Rij[ "lj" "S" "S" ] = 4.00
Rij[ "lj" "S" "H" ] = 3.00
Rij[ "lj" "S" "f" ] = 2.65
Rij[ "lj" "S" "F" ] = 3.54
Rij[ "lj" "S" "c" ] = 4.04
Rij[ "lj" "S" "b" ] = 4.17
Rij[ "lj" "S" "I" ] = 4.36
#Rij[ "lj" "S" "M" ] = 2.79
Rij[ "lj" "S" "M" ] = 2.65  #Mg
Rij[ "lj" "S" "Z" ] = 2.74  #Zn
Rij[ "lj" "S" "L" ] = 2.99  #Ca
Rij[ "lj" "S" "X" ] = 3.00  #Xx

Rij[ "lj" "H" "C" ] = 3.00
Rij[ "lj" "H" "N" ] = 2.75
Rij[ "lj" "H" "n" ] = 2.75
Rij[ "lj" "H" "O" ] = 2.60
Rij[ "lj" "H" "P" ] = 3.10
Rij[ "lj" "H" "S" ] = 3.00
Rij[ "lj" "H" "H" ] = 2.00
Rij[ "lj" "H" "f" ] = 1.65
Rij[ "lj" "H" "F" ] = 2.54
Rij[ "lj" "H" "c" ] = 3.04
Rij[ "lj" "H" "b" ] = 3.17
Rij[ "lj" "H" "I" ] = 3.36
#Rij[ "lj" "H" "M" ] = 1.79
Rij[ "lj" "H" "M" ] = 1.65  #Mg
Rij[ "lj" "H" "Z" ] = 1.74  #Zn
Rij[ "lj" "H" "L" ] = 1.99  #Ca
Rij[ "lj" "H" "X" ] = 2.00  #Xx

Rij[ "lj" "f" "C" ] = 2.65
Rij[ "lj" "f" "N" ] = 2.40
Rij[ "lj" "f" "n" ] = 2.40
Rij[ "lj" "f" "O" ] = 2.25
Rij[ "lj" "f" "P" ] = 2.75
Rij[ "lj" "f" "S" ] = 2.65
Rij[ "lj" "f" "H" ] = 1.65
Rij[ "lj" "f" "f" ] = 1.30
Rij[ "lj" "f" "F" ] = 2.19
Rij[ "lj" "f" "c" ] = 2.69
Rij[ "lj" "f" "b" ] = 2.81
Rij[ "lj" "f" "I" ] = 3.01
#Rij[ "lj" "f" "M" ] = 1.44
Rij[ "lj" "f" "M" ] = 1.30  #Mg
Rij[ "lj" "f" "Z" ] = 1.39  #Zn
Rij[ "lj" "f" "L" ] = 1.64  #Ca
Rij[ "lj" "f" "X" ] = 1.65  #Xx

Rij[ "lj" "F" "C" ] = 3.54
Rij[ "lj" "F" "N" ] = 3.29
Rij[ "lj" "F" "n" ] = 3.29
Rij[ "lj" "F" "O" ] = 3.15
Rij[ "lj" "F" "P" ] = 3.65
Rij[ "lj" "F" "S" ] = 3.54
Rij[ "lj" "F" "H" ] = 2.54
Rij[ "lj" "F" "f" ] = 2.19
Rij[ "lj" "F" "F" ] = 3.09
Rij[ "lj" "F" "c" ] = 3.59
Rij[ "lj" "F" "b" ] = 3.71
Rij[ "lj" "F" "I" ] = 3.90
#Rij[ "lj" "F" "M" ] = 2.33
Rij[ "lj" "F" "M" ] = 2.19  #Mg
Rij[ "lj" "F" "Z" ] = 2.29  #Zn
Rij[ "lj" "F" "L" ] = 2.54  #Ca
Rij[ "lj" "F" "X" ] = 2.54  #Xx

Rij[ "lj" "c" "C" ] = 4.04
Rij[ "lj" "c" "N" ] = 3.79
Rij[ "lj" "c" "n" ] = 3.79
Rij[ "lj" "c" "O" ] = 3.65
Rij[ "lj" "c" "P" ] = 4.14
Rij[ "lj" "c" "S" ] = 4.04
Rij[ "lj" "c" "H" ] = 3.04
Rij[ "lj" "c" "f" ] = 2.69
Rij[ "lj" "c" "F" ] = 3.59
Rij[ "lj" "c" "c" ] = 4.09
Rij[ "lj" "c" "b" ] = 4.21
Rij[ "lj" "c" "I" ] = 4.40
#Rij[ "lj" "c" "M" ] = 2.83
Rij[ "lj" "c" "M" ] = 2.69  #Mg
Rij[ "lj" "c" "Z" ] = 2.79  #Zn
Rij[ "lj" "c" "L" ] = 3.04  #Ca
Rij[ "lj" "c" "X" ] = 3.04  #Xx

Rij[ "lj" "b" "C" ] = 4.17
Rij[ "lj" "b" "N" ] = 3.92
Rij[ "lj" "b" "n" ] = 3.92
Rij[ "lj" "b" "O" ] = 3.77
Rij[ "lj" "b" "P" ] = 4.27
Rij[ "lj" "b" "S" ] = 4.17
Rij[ "lj" "b" "H" ] = 3.17
Rij[ "lj" "b" "f" ] = 2.81
Rij[ "lj" "b" "F" ] = 3.71
Rij[ "lj" "b" "c" ] = 4.21
Rij[ "lj" "b" "b" ] = 4.33
Rij[ "lj" "b" "I" ] = 4.53
#Rij[ "lj" "b" "M" ] = 2.95
Rij[ "lj" "b" "M" ] = 2.81  #Mg
Rij[ "lj" "b" "Z" ] = 2.91  #Zn
Rij[ "lj" "b" "L" ] = 3.16  #Ca
Rij[ "lj" "b" "X" ] = 3.17  #Xx

Rij[ "lj" "I" "C" ] = 4.36
Rij[ "lj" "I" "N" ] = 4.11
Rij[ "lj" "I" "n" ] = 4.11
Rij[ "lj" "I" "O" ] = 3.96
Rij[ "lj" "I" "P" ] = 4.46
Rij[ "lj" "I" "S" ] = 4.36
Rij[ "lj" "I" "H" ] = 3.36
Rij[ "lj" "I" "f" ] = 3.01
Rij[ "lj" "I" "F" ] = 3.90
Rij[ "lj" "I" "c" ] = 4.40
Rij[ "lj" "I" "b" ] = 4.53
Rij[ "lj" "I" "I" ] = 4.72
#Rij[ "lj" "I" "M" ] = 3.15
Rij[ "lj" "I" "M" ] = 3.01  #Mg
Rij[ "lj" "I" "Z" ] = 3.10  #Zn
Rij[ "lj" "I" "L" ] = 3.35  #Ca
Rij[ "lj" "I" "X" ] = 3.36  #Xx

#Rij[ "lj" "M" "C" ] = 2.79
#Rij[ "lj" "M" "N" ] = 2.54
#Rij[ "lj" "M" "n" ] = 2.54
#Rij[ "lj" "M" "O" ] = 2.39
#Rij[ "lj" "M" "P" ] = 2.89
#Rij[ "lj" "M" "S" ] = 2.79
#Rij[ "lj" "M" "H" ] = 1.79
#Rij[ "lj" "M" "f" ] = 1.44
#Rij[ "lj" "M" "F" ] = 2.33
#Rij[ "lj" "M" "c" ] = 2.83
#Rij[ "lj" "M" "b" ] = 2.95
#Rij[ "lj" "M" "I" ] = 3.15
#Rij[ "lj" "M" "M" ] = 1.57

#Mg
Rij[ "lj" "M" "C" ] = 2.65
Rij[ "lj" "M" "N" ] = 2.40
Rij[ "lj" "M" "n" ] = 2.40
Rij[ "lj" "M" "O" ] = 2.25
Rij[ "lj" "M" "P" ] = 2.75
Rij[ "lj" "M" "S" ] = 2.65
Rij[ "lj" "M" "H" ] = 1.65
Rij[ "lj" "M" "f" ] = 1.30
Rij[ "lj" "M" "F" ] = 2.19
Rij[ "lj" "M" "c" ] = 2.69
Rij[ "lj" "M" "b" ] = 2.81
Rij[ "lj" "M" "I" ] = 3.01
Rij[ "lj" "M" "M" ] = 1.30  #Mg
Rij[ "lj" "M" "Z" ] = 1.39  #Zn
Rij[ "lj" "M" "L" ] = 1.64  #Ca
Rij[ "lj" "M" "X" ] = 1.65  #Xx


#Zn
Rij[ "lj" "Z" "C" ] = 2.74
Rij[ "lj" "Z" "N" ] = 2.49
Rij[ "lj" "Z" "n" ] = 2.49
Rij[ "lj" "Z" "O" ] = 2.34
Rij[ "lj" "Z" "P" ] = 2.84
Rij[ "lj" "Z" "S" ] = 2.74
Rij[ "lj" "Z" "H" ] = 1.74
Rij[ "lj" "Z" "f" ] = 1.39
Rij[ "lj" "Z" "F" ] = 2.29
Rij[ "lj" "Z" "c" ] = 2.79
Rij[ "lj" "Z" "b" ] = 2.91
Rij[ "lj" "Z" "I" ] = 3.10
Rij[ "lj" "Z" "M" ] = 1.39  #Mg
Rij[ "lj" "Z" "Z" ] = 1.48  #Zn
Rij[ "lj" "Z" "L" ] = 1.73  #Ca
Rij[ "lj" "Z" "X" ] = 1.74  #Xx

#Ca
Rij[ "lj" "L" "C" ] = 2.99
Rij[ "lj" "L" "N" ] = 2.74
Rij[ "lj" "L" "n" ] = 2.74
Rij[ "lj" "L" "O" ] = 2.59
Rij[ "lj" "L" "P" ] = 3.09
Rij[ "lj" "L" "S" ] = 2.99
Rij[ "lj" "L" "H" ] = 1.99
Rij[ "lj" "L" "f" ] = 1.64
Rij[ "lj" "L" "F" ] = 2.54
Rij[ "lj" "L" "c" ] = 3.04
Rij[ "lj" "L" "b" ] = 3.16
Rij[ "lj" "L" "I" ] = 3.35
Rij[ "lj" "L" "M" ] = 1.64  #Mg
Rij[ "lj" "L" "Z" ] = 1.73  #Zn
Rij[ "lj" "L" "L" ] = 1.98  #Ca
Rij[ "lj" "L" "X" ] = 1.99  #Xx

#Xx
Rij[ "lj" "X" "C" ] = 3.00
Rij[ "lj" "X" "N" ] = 2.75
Rij[ "lj" "X" "n" ] = 2.75
Rij[ "lj" "X" "O" ] = 2.60
Rij[ "lj" "X" "P" ] = 3.10
Rij[ "lj" "X" "S" ] = 3.00
Rij[ "lj" "X" "H" ] = 2.00
Rij[ "lj" "X" "f" ] = 1.65
Rij[ "lj" "X" "F" ] = 2.54
Rij[ "lj" "X" "c" ] = 3.04
Rij[ "lj" "X" "b" ] = 3.17
Rij[ "lj" "X" "I" ] = 3.36
Rij[ "lj" "X" "M" ] = 1.65  #Mg
Rij[ "lj" "X" "Z" ] = 1.74  #Zn
Rij[ "lj" "X" "L" ] = 1.99  #Ca
Rij[ "lj" "X" "X" ] = 2.00  #Xx


#A:aromatic carbon
Rij[ "lj" "C" "A" ] = Rij[ "lj" "C" "C" ]
Rij[ "lj" "A" "C" ] = Rij[ "lj" "C" "A" ]
Rij[ "lj" "A" "A" ] = Rij[ "lj" "C" "C" ]

Rij[ "lj" "A" "N" ] = Rij[ "lj" "C" "N" ]
Rij[ "lj" "A" "n" ] = Rij[ "lj" "C" "n" ]
Rij[ "lj" "A" "O" ] = Rij[ "lj" "C" "O" ]
Rij[ "lj" "A" "S" ] = Rij[ "lj" "C" "S" ]
Rij[ "lj" "A" "H" ] = Rij[ "lj" "C" "H" ]
Rij[ "lj" "A" "P" ] = Rij[ "lj" "C" "P" ]
Rij[ "lj" "A" "f" ] = Rij[ "lj" "C" "f" ]
Rij[ "lj" "A" "F" ] = Rij[ "lj" "C" "F" ]
Rij[ "lj" "A" "c" ] = Rij[ "lj" "C" "c" ]
Rij[ "lj" "A" "b" ] = Rij[ "lj" "C" "b" ]
Rij[ "lj" "A" "I" ] = Rij[ "lj" "C" "I" ]
Rij[ "lj" "A" "M" ] = Rij[ "lj" "C" "M" ]
Rij[ "lj" "A" "Z" ] = Rij[ "lj" "C" "Z" ]
Rij[ "lj" "A" "L" ] = Rij[ "lj" "C" "L" ]
Rij[ "lj" "A" "X" ] = Rij[ "lj" "C" "X" ]


Rij[ "lj" "N" "A" ] = Rij[ "lj" "A" "N" ]
Rij[ "lj" "n" "A" ] = Rij[ "lj" "A" "n" ]
Rij[ "lj" "O" "A" ] = Rij[ "lj" "A" "O" ]
Rij[ "lj" "S" "A" ] = Rij[ "lj" "A" "S" ]
Rij[ "lj" "H" "A" ] = Rij[ "lj" "A" "H" ]
Rij[ "lj" "P" "A" ] = Rij[ "lj" "A" "P" ]
Rij[ "lj" "f" "A" ] = Rij[ "lj" "A" "f" ]
Rij[ "lj" "F" "A" ] = Rij[ "lj" "A" "F" ]
Rij[ "lj" "c" "A" ] = Rij[ "lj" "A" "c" ]
Rij[ "lj" "b" "A" ] = Rij[ "lj" "A" "b" ]
Rij[ "lj" "I" "A" ] = Rij[ "lj" "A" "I" ]
Rij[ "lj" "M" "A" ] = Rij[ "lj" "A" "M" ]
Rij[ "lj" "Z" "A" ] = Rij[ "lj" "A" "Z" ] 
Rij[ "lj" "L" "A" ] = Rij[ "lj" "A" "M" ]
Rij[ "lj" "X" "A" ] = Rij[ "lj" "A" "X" ]

#
# AutoDock Binding Free Energy Model 140n
# New well-depths using new Solvation Model, 
# multiplied by van der Waals coefficient, 0.1485
#

#
# Well-depths prior to Solvation Model, multiplied by FE_vdW_coeff
#
epsij={}
epsij[ "lj" "C" "C" ] = 0.1500 * FE_vdW_coeff
epsij[ "lj" "C" "N" ] = 0.1549 * FE_vdW_coeff
epsij[ "lj" "C" "n" ] = 0.1549 * FE_vdW_coeff
epsij[ "lj" "C" "O" ] = 0.1732 * FE_vdW_coeff
epsij[ "lj" "C" "P" ] = 0.1732 * FE_vdW_coeff
epsij[ "lj" "C" "S" ] = 0.1732 * FE_vdW_coeff
epsij[ "lj" "C" "H" ] = 0.0548 * FE_vdW_coeff
epsij[ "lj" "C" "f" ] = 0.0387 * FE_vdW_coeff
epsij[ "lj" "C" "F" ] = 0.1095 * FE_vdW_coeff
epsij[ "lj" "C" "c" ] = 0.2035 * FE_vdW_coeff
epsij[ "lj" "C" "b" ] = 0.2416 * FE_vdW_coeff
epsij[ "lj" "C" "I" ] = 0.2877 * FE_vdW_coeff
#epsij[ "lj" "C" "M" ] = 0.3623 * FE_vdW_coeff
epsij[ "lj" "C" "M" ] = 0.3623 * FE_vdW_coeff
epsij[ "lj" "C" "Z" ] = 0.2872 * FE_vdW_coeff
epsij[ "lj" "C" "L" ] = 0.2872 * FE_vdW_coeff
epsij[ "lj" "C" "X" ] = 0.0548 * FE_vdW_coeff

epsij[ "lj" "N" "C" ] = 0.1549 * FE_vdW_coeff
epsij[ "lj" "N" "N" ] = 0.1600 * FE_vdW_coeff
epsij[ "lj" "N" "n" ] = 0.1600 * FE_vdW_coeff
epsij[ "lj" "N" "O" ] = 0.1789 * FE_vdW_coeff
epsij[ "lj" "N" "P" ] = 0.1789 * FE_vdW_coeff
epsij[ "lj" "N" "S" ] = 0.1789 * FE_vdW_coeff
epsij[ "lj" "N" "H" ] = 0.0566 * FE_vdW_coeff
epsij[ "lj" "N" "f" ] = 0.0400 * FE_vdW_coeff
epsij[ "lj" "N" "F" ] = 0.1131 * FE_vdW_coeff
epsij[ "lj" "N" "c" ] = 0.2101 * FE_vdW_coeff
epsij[ "lj" "N" "b" ] = 0.2495 * FE_vdW_coeff
epsij[ "lj" "N" "I" ] = 0.2972 * FE_vdW_coeff
#epsij[ "lj" "N" "M" ] = 0.3742 * FE_vdW_coeff
epsij[ "lj" "N" "M" ] = 0.3742 * FE_vdW_coeff
epsij[ "lj" "N" "Z" ] = 0.2966 * FE_vdW_coeff
epsij[ "lj" "N" "L" ] = 0.2966 * FE_vdW_coeff
epsij[ "lj" "N" "X" ] = 0.0566 * FE_vdW_coeff

epsij[ "lj" "n" "C" ] = 0.1549 * FE_vdW_coeff
epsij[ "lj" "n" "N" ] = 0.1600 * FE_vdW_coeff
epsij[ "lj" "n" "n" ] = 0.1600 * FE_vdW_coeff
epsij[ "lj" "n" "O" ] = 0.1789 * FE_vdW_coeff
epsij[ "lj" "n" "P" ] = 0.1789 * FE_vdW_coeff
epsij[ "lj" "n" "S" ] = 0.1789 * FE_vdW_coeff
epsij[ "lj" "n" "H" ] = 0.0566 * FE_vdW_coeff
epsij[ "lj" "n" "f" ] = 0.0400 * FE_vdW_coeff
epsij[ "lj" "n" "F" ] = 0.1131 * FE_vdW_coeff
epsij[ "lj" "n" "c" ] = 0.2101 * FE_vdW_coeff
epsij[ "lj" "n" "b" ] = 0.2495 * FE_vdW_coeff
epsij[ "lj" "n" "I" ] = 0.2972 * FE_vdW_coeff
#epsij[ "lj" "n" "M" ] = 0.3742 * FE_vdW_coeff
epsij[ "lj" "n" "M" ] = 0.3742 * FE_vdW_coeff
epsij[ "lj" "n" "Z" ] = 0.2966 * FE_vdW_coeff
epsij[ "lj" "n" "L" ] = 0.2966 * FE_vdW_coeff
epsij[ "lj" "n" "X" ] = 0.0566 * FE_vdW_coeff

epsij[ "lj" "O" "C" ] = 0.1732 * FE_vdW_coeff
epsij[ "lj" "O" "N" ] = 0.1789 * FE_vdW_coeff
epsij[ "lj" "O" "n" ] = 0.1789 * FE_vdW_coeff
epsij[ "lj" "O" "O" ] = 0.2000 * FE_vdW_coeff
epsij[ "lj" "O" "P" ] = 0.2000 * FE_vdW_coeff
epsij[ "lj" "O" "S" ] = 0.2000 * FE_vdW_coeff
epsij[ "lj" "O" "H" ] = 0.0632 * FE_vdW_coeff
epsij[ "lj" "O" "f" ] = 0.0447 * FE_vdW_coeff
epsij[ "lj" "O" "F" ] = 0.1265 * FE_vdW_coeff
epsij[ "lj" "O" "c" ] = 0.2349 * FE_vdW_coeff
epsij[ "lj" "O" "b" ] = 0.2789 * FE_vdW_coeff
epsij[ "lj" "O" "I" ] = 0.3323 * FE_vdW_coeff
#epsij[ "lj" "O" "M" ] = 0.4183 * FE_vdW_coeff
epsij[ "lj" "O" "M" ] = 0.4183 * FE_vdW_coeff
epsij[ "lj" "O" "Z" ] = 0.3317 * FE_vdW_coeff
epsij[ "lj" "O" "L" ] = 0.3317 * FE_vdW_coeff
epsij[ "lj" "O" "X" ] = 0.0632 * FE_vdW_coeff

epsij[ "lj" "P" "C" ] = 0.1732 * FE_vdW_coeff
epsij[ "lj" "P" "N" ] = 0.1789 * FE_vdW_coeff
epsij[ "lj" "P" "n" ] = 0.1789 * FE_vdW_coeff
epsij[ "lj" "P" "O" ] = 0.2000 * FE_vdW_coeff
epsij[ "lj" "P" "P" ] = 0.2000 * FE_vdW_coeff
epsij[ "lj" "P" "S" ] = 0.2000 * FE_vdW_coeff
epsij[ "lj" "P" "H" ] = 0.0632 * FE_vdW_coeff
epsij[ "lj" "P" "f" ] = 0.0447 * FE_vdW_coeff
epsij[ "lj" "P" "F" ] = 0.1265 * FE_vdW_coeff
epsij[ "lj" "P" "c" ] = 0.2349 * FE_vdW_coeff
epsij[ "lj" "P" "b" ] = 0.2789 * FE_vdW_coeff
epsij[ "lj" "P" "I" ] = 0.3323 * FE_vdW_coeff
#epsij[ "lj" "P" "M" ] = 0.4183 * FE_vdW_coeff
epsij[ "lj" "P" "M" ] = 0.4183 * FE_vdW_coeff       #NB: same as O values
epsij[ "lj" "P" "Z" ] = 0.3317 * FE_vdW_coeff
epsij[ "lj" "P" "L" ] = 0.3317 * FE_vdW_coeff
epsij[ "lj" "P" "X" ] = 0.0632 * FE_vdW_coeff

epsij[ "lj" "S" "C" ] = 0.1732 * FE_vdW_coeff
epsij[ "lj" "S" "N" ] = 0.1789 * FE_vdW_coeff
epsij[ "lj" "S" "n" ] = 0.1789 * FE_vdW_coeff
epsij[ "lj" "S" "O" ] = 0.2000 * FE_vdW_coeff
epsij[ "lj" "S" "P" ] = 0.2000 * FE_vdW_coeff
epsij[ "lj" "S" "S" ] = 0.2000 * FE_vdW_coeff
epsij[ "lj" "S" "H" ] = 0.0632 * FE_vdW_coeff
epsij[ "lj" "S" "f" ] = 0.0447 * FE_vdW_coeff
epsij[ "lj" "S" "F" ] = 0.1265 * FE_vdW_coeff
epsij[ "lj" "S" "c" ] = 0.2349 * FE_vdW_coeff
epsij[ "lj" "S" "b" ] = 0.2789 * FE_vdW_coeff
epsij[ "lj" "S" "I" ] = 0.3323 * FE_vdW_coeff
#epsij[ "lj" "S" "M" ] = 0.4183 * FE_vdW_coeff
epsij[ "lj" "S" "M" ] = 0.4183 * FE_vdW_coeff       #NB: same as O values
epsij[ "lj" "S" "Z" ] = 0.3317 * FE_vdW_coeff
epsij[ "lj" "S" "L" ] = 0.3317 * FE_vdW_coeff
epsij[ "lj" "S" "X" ] = 0.0632 * FE_vdW_coeff

epsij[ "lj" "H" "C" ] = 0.0548 * FE_vdW_coeff
epsij[ "lj" "H" "N" ] = 0.0566 * FE_vdW_coeff
epsij[ "lj" "H" "n" ] = 0.0566 * FE_vdW_coeff
epsij[ "lj" "H" "O" ] = 0.0632 * FE_vdW_coeff
epsij[ "lj" "H" "P" ] = 0.0632 * FE_vdW_coeff
epsij[ "lj" "H" "S" ] = 0.0632 * FE_vdW_coeff
epsij[ "lj" "H" "H" ] = 0.0200 * FE_vdW_coeff
epsij[ "lj" "H" "f" ] = 0.0141 * FE_vdW_coeff
epsij[ "lj" "H" "F" ] = 0.0400 * FE_vdW_coeff
epsij[ "lj" "H" "c" ] = 0.0743 * FE_vdW_coeff
epsij[ "lj" "H" "b" ] = 0.0882 * FE_vdW_coeff
epsij[ "lj" "H" "I" ] = 0.1051 * FE_vdW_coeff
#epsij[ "lj" "H" "M" ] = 0.1323 * FE_vdW_coeff
epsij[ "lj" "H" "M" ] = 0.1323 * FE_vdW_coeff
epsij[ "lj" "H" "Z" ] = 0.1049 * FE_vdW_coeff
epsij[ "lj" "H" "L" ] = 0.1049 * FE_vdW_coeff
epsij[ "lj" "H" "X" ] = 0.0200 * FE_vdW_coeff

epsij[ "lj" "f" "C" ] = 0.0387 * FE_vdW_coeff
epsij[ "lj" "f" "N" ] = 0.0400 * FE_vdW_coeff
epsij[ "lj" "f" "n" ] = 0.0400 * FE_vdW_coeff
epsij[ "lj" "f" "O" ] = 0.0447 * FE_vdW_coeff
epsij[ "lj" "f" "P" ] = 0.0447 * FE_vdW_coeff
epsij[ "lj" "f" "S" ] = 0.0447 * FE_vdW_coeff
epsij[ "lj" "f" "H" ] = 0.0141 * FE_vdW_coeff
epsij[ "lj" "f" "f" ] = 0.0100 * FE_vdW_coeff
epsij[ "lj" "f" "F" ] = 0.0283 * FE_vdW_coeff
epsij[ "lj" "f" "c" ] = 0.0525 * FE_vdW_coeff
epsij[ "lj" "f" "b" ] = 0.0624 * FE_vdW_coeff
epsij[ "lj" "f" "I" ] = 0.0743 * FE_vdW_coeff
#epsij[ "lj" "f" "M" ] = 0.0935 * FE_vdW_coeff
epsij[ "lj" "f" "M" ] = 0.0935 * FE_vdW_coeff
epsij[ "lj" "f" "Z" ] = 0.0742 * FE_vdW_coeff
epsij[ "lj" "f" "L" ] = 0.0742 * FE_vdW_coeff
epsij[ "lj" "f" "X" ] = 0.0141 * FE_vdW_coeff

epsij[ "lj" "F" "C" ] = 0.1095 * FE_vdW_coeff
epsij[ "lj" "F" "N" ] = 0.1131 * FE_vdW_coeff
epsij[ "lj" "F" "n" ] = 0.1131 * FE_vdW_coeff
epsij[ "lj" "F" "O" ] = 0.1265 * FE_vdW_coeff
epsij[ "lj" "F" "P" ] = 0.1265 * FE_vdW_coeff
epsij[ "lj" "F" "S" ] = 0.1265 * FE_vdW_coeff
epsij[ "lj" "F" "H" ] = 0.0400 * FE_vdW_coeff
epsij[ "lj" "F" "f" ] = 0.0283 * FE_vdW_coeff
epsij[ "lj" "F" "F" ] = 0.0800 * FE_vdW_coeff
epsij[ "lj" "F" "c" ] = 0.1486 * FE_vdW_coeff
epsij[ "lj" "F" "b" ] = 0.1764 * FE_vdW_coeff
epsij[ "lj" "F" "I" ] = 0.2101 * FE_vdW_coeff
#epsij[ "lj" "F" "M" ] = 0.2646 * FE_vdW_coeff
epsij[ "lj" "F" "M" ] = 0.2646 * FE_vdW_coeff
epsij[ "lj" "F" "Z" ] = 0.2098 * FE_vdW_coeff
epsij[ "lj" "F" "L" ] = 0.2098 * FE_vdW_coeff
epsij[ "lj" "F" "X" ] = 0.0400 * FE_vdW_coeff

epsij[ "lj" "c" "C" ] = 0.2035 * FE_vdW_coeff
epsij[ "lj" "c" "N" ] = 0.2101 * FE_vdW_coeff
epsij[ "lj" "c" "n" ] = 0.2101 * FE_vdW_coeff
epsij[ "lj" "c" "O" ] = 0.2349 * FE_vdW_coeff
epsij[ "lj" "c" "P" ] = 0.2349 * FE_vdW_coeff
epsij[ "lj" "c" "S" ] = 0.2349 * FE_vdW_coeff
epsij[ "lj" "c" "H" ] = 0.0743 * FE_vdW_coeff
epsij[ "lj" "c" "f" ] = 0.0525 * FE_vdW_coeff
epsij[ "lj" "c" "F" ] = 0.1486 * FE_vdW_coeff
epsij[ "lj" "c" "c" ] = 0.2760 * FE_vdW_coeff
epsij[ "lj" "c" "b" ] = 0.3277 * FE_vdW_coeff
epsij[ "lj" "c" "I" ] = 0.3903 * FE_vdW_coeff
#epsij[ "lj" "c" "M" ] = 0.4914 * FE_vdW_coeff
epsij[ "lj" "c" "M" ] = 0.4914 * FE_vdW_coeff
epsij[ "lj" "c" "Z" ] = 0.3896 * FE_vdW_coeff
epsij[ "lj" "c" "L" ] = 0.3896 * FE_vdW_coeff
epsij[ "lj" "c" "X" ] = 0.0743 * FE_vdW_coeff

epsij[ "lj" "b" "C" ] = 0.2416 * FE_vdW_coeff
epsij[ "lj" "b" "N" ] = 0.2495 * FE_vdW_coeff
epsij[ "lj" "b" "n" ] = 0.2495 * FE_vdW_coeff
epsij[ "lj" "b" "O" ] = 0.2789 * FE_vdW_coeff
epsij[ "lj" "b" "P" ] = 0.2789 * FE_vdW_coeff
epsij[ "lj" "b" "S" ] = 0.2789 * FE_vdW_coeff
epsij[ "lj" "b" "H" ] = 0.0882 * FE_vdW_coeff
epsij[ "lj" "b" "f" ] = 0.0624 * FE_vdW_coeff
epsij[ "lj" "b" "F" ] = 0.1764 * FE_vdW_coeff
epsij[ "lj" "b" "c" ] = 0.3277 * FE_vdW_coeff
epsij[ "lj" "b" "b" ] = 0.3890 * FE_vdW_coeff
epsij[ "lj" "b" "I" ] = 0.4634 * FE_vdW_coeff
#epsij[ "lj" "b" "M" ] = 0.5834 * FE_vdW_coeff
epsij[ "lj" "b" "M" ] = 0.5834 * FE_vdW_coeff
epsij[ "lj" "b" "Z" ] = 0.4625 * FE_vdW_coeff
epsij[ "lj" "b" "L" ] = 0.4625 * FE_vdW_coeff
epsij[ "lj" "b" "X" ] = 0.0882 * FE_vdW_coeff

epsij[ "lj" "I" "C" ] = 0.2877 * FE_vdW_coeff
epsij[ "lj" "I" "N" ] = 0.2972 * FE_vdW_coeff
epsij[ "lj" "I" "n" ] = 0.2972 * FE_vdW_coeff
epsij[ "lj" "I" "O" ] = 0.3323 * FE_vdW_coeff
epsij[ "lj" "I" "P" ] = 0.3323 * FE_vdW_coeff
epsij[ "lj" "I" "S" ] = 0.3323 * FE_vdW_coeff
epsij[ "lj" "I" "H" ] = 0.1051 * FE_vdW_coeff
epsij[ "lj" "I" "f" ] = 0.0743 * FE_vdW_coeff
epsij[ "lj" "I" "F" ] = 0.2101 * FE_vdW_coeff
epsij[ "lj" "I" "c" ] = 0.3903 * FE_vdW_coeff
epsij[ "lj" "I" "b" ] = 0.4634 * FE_vdW_coeff
epsij[ "lj" "I" "I" ] = 0.5520 * FE_vdW_coeff
#epsij[ "lj" "I" "M" ] = 0.6950 * FE_vdW_coeff
epsij[ "lj" "I" "M" ] = 0.6950 * FE_vdW_coeff
epsij[ "lj" "I" "Z" ] = 0.5510 * FE_vdW_coeff
epsij[ "lj" "I" "L" ] = 0.5510 * FE_vdW_coeff
epsij[ "lj" "I" "X" ] = 0.1051 * FE_vdW_coeff

epsij[ "lj" "M" "C" ] = 0.3623 * FE_vdW_coeff
epsij[ "lj" "M" "N" ] = 0.3742 * FE_vdW_coeff
epsij[ "lj" "M" "n" ] = 0.3742 * FE_vdW_coeff
epsij[ "lj" "M" "O" ] = 0.4183 * FE_vdW_coeff
epsij[ "lj" "M" "P" ] = 0.4183 * FE_vdW_coeff
epsij[ "lj" "M" "S" ] = 0.4183 * FE_vdW_coeff
epsij[ "lj" "M" "H" ] = 0.1323 * FE_vdW_coeff
epsij[ "lj" "M" "f" ] = 0.0935 * FE_vdW_coeff
epsij[ "lj" "M" "F" ] = 0.2646 * FE_vdW_coeff
epsij[ "lj" "M" "c" ] = 0.4914 * FE_vdW_coeff
epsij[ "lj" "M" "b" ] = 0.5834 * FE_vdW_coeff
epsij[ "lj" "M" "I" ] = 0.6950 * FE_vdW_coeff
epsij[ "lj" "M" "M" ] = 0.8750 * FE_vdW_coeff
epsij[ "lj" "M" "Z" ] = 0.6937 * FE_vdW_coeff
epsij[ "lj" "M" "L" ] = 0.6937 * FE_vdW_coeff
epsij[ "lj" "M" "X" ] = 0.1323 * FE_vdW_coeff

epsij[ "lj" "Z" "C" ] = 0.2872 * FE_vdW_coeff
epsij[ "lj" "Z" "N" ] = 0.2966 * FE_vdW_coeff
epsij[ "lj" "Z" "n" ] = 0.2966 * FE_vdW_coeff
epsij[ "lj" "Z" "O" ] = 0.3317 * FE_vdW_coeff
epsij[ "lj" "Z" "P" ] = 0.3317 * FE_vdW_coeff
epsij[ "lj" "Z" "S" ] = 0.3317 * FE_vdW_coeff
epsij[ "lj" "Z" "H" ] = 0.1049 * FE_vdW_coeff
epsij[ "lj" "Z" "f" ] = 0.0742 * FE_vdW_coeff
epsij[ "lj" "Z" "F" ] = 0.2098 * FE_vdW_coeff
epsij[ "lj" "Z" "c" ] = 0.3896 * FE_vdW_coeff
epsij[ "lj" "Z" "b" ] = 0.4625 * FE_vdW_coeff
epsij[ "lj" "Z" "I" ] = 0.5510 * FE_vdW_coeff
epsij[ "lj" "Z" "M" ] = 0.6937 * FE_vdW_coeff
epsij[ "lj" "Z" "Z" ] = 0.5500 * FE_vdW_coeff
epsij[ "lj" "Z" "L" ] = 0.5500 * FE_vdW_coeff
epsij[ "lj" "Z" "X" ] = 0.1049 * FE_vdW_coeff

epsij[ "lj" "L" "C" ] = 0.2872 * FE_vdW_coeff   #these are the same as Zn values
epsij[ "lj" "L" "N" ] = 0.2966 * FE_vdW_coeff
epsij[ "lj" "L" "n" ] = 0.2966 * FE_vdW_coeff
epsij[ "lj" "L" "O" ] = 0.3317 * FE_vdW_coeff
epsij[ "lj" "L" "P" ] = 0.3317 * FE_vdW_coeff
epsij[ "lj" "L" "S" ] = 0.3317 * FE_vdW_coeff
epsij[ "lj" "L" "H" ] = 0.1049 * FE_vdW_coeff
epsij[ "lj" "L" "f" ] = 0.0742 * FE_vdW_coeff
epsij[ "lj" "L" "F" ] = 0.2098 * FE_vdW_coeff
epsij[ "lj" "L" "c" ] = 0.3896 * FE_vdW_coeff
epsij[ "lj" "L" "b" ] = 0.4625 * FE_vdW_coeff
epsij[ "lj" "L" "I" ] = 0.5510 * FE_vdW_coeff
epsij[ "lj" "L" "M" ] = 0.6937 * FE_vdW_coeff
epsij[ "lj" "L" "Z" ] = 0.5500 * FE_vdW_coeff
epsij[ "lj" "L" "L" ] = 0.5500 * FE_vdW_coeff
epsij[ "lj" "L" "X" ] = 0.1049 * FE_vdW_coeff

epsij[ "lj" "X" "C" ] = 0.0548 * FE_vdW_coeff
epsij[ "lj" "X" "N" ] = 0.0566 * FE_vdW_coeff
epsij[ "lj" "X" "n" ] = 0.0566 * FE_vdW_coeff
epsij[ "lj" "X" "O" ] = 0.0632 * FE_vdW_coeff
epsij[ "lj" "X" "P" ] = 0.0632 * FE_vdW_coeff
epsij[ "lj" "X" "S" ] = 0.0632 * FE_vdW_coeff
epsij[ "lj" "X" "H" ] = 0.0200 * FE_vdW_coeff
epsij[ "lj" "X" "f" ] = 0.0141 * FE_vdW_coeff
epsij[ "lj" "X" "F" ] = 0.0400 * FE_vdW_coeff
epsij[ "lj" "X" "c" ] = 0.0743 * FE_vdW_coeff
epsij[ "lj" "X" "b" ] = 0.0882 * FE_vdW_coeff
epsij[ "lj" "X" "I" ] = 0.1051 * FE_vdW_coeff
epsij[ "lj" "X" "M" ] = 0.1323 * FE_vdW_coeff
epsij[ "lj" "X" "Z" ] = 0.1049 * FE_vdW_coeff
epsij[ "lj" "X" "L" ] = 0.1049 * FE_vdW_coeff
epsij[ "lj" "X" "X" ] = 0.0200 * FE_vdW_coeff

epsij[ "lj" "C" "A" ] = 0.1500 * FE_vdW_coeff
epsij[ "lj" "A" "C" ] = epsij[ "lj" "C" "A" ]
epsij[ "lj" "A" "A" ] = epsij[ "lj" "C" "C" ]

epsij[ "lj" "A" "N" ] = epsij[ "lj" "C" "N" ]
epsij[ "lj" "A" "n" ] = epsij[ "lj" "C" "n" ]
epsij[ "lj" "A" "O" ] = epsij[ "lj" "C" "O" ]
epsij[ "lj" "A" "S" ] = epsij[ "lj" "C" "S" ]
epsij[ "lj" "A" "H" ] = epsij[ "lj" "C" "H" ]
epsij[ "lj" "A" "P" ] = epsij[ "lj" "C" "P" ]
epsij[ "lj" "A" "f" ] = epsij[ "lj" "C" "f" ]
epsij[ "lj" "A" "F" ] = epsij[ "lj" "C" "F" ]
epsij[ "lj" "A" "c" ] = epsij[ "lj" "C" "c" ]
epsij[ "lj" "A" "b" ] = epsij[ "lj" "C" "b" ]
epsij[ "lj" "A" "I" ] = epsij[ "lj" "C" "I" ]
epsij[ "lj" "A" "M" ] = epsij[ "lj" "C" "M" ]
epsij[ "lj" "A" "Z" ] = epsij[ "lj" "C" "Z" ]
epsij[ "lj" "A" "L" ] = epsij[ "lj" "C" "L" ]
epsij[ "lj" "A" "X" ] = epsij[ "lj" "C" "X" ]


epsij[ "lj" "N" "A" ] = epsij[ "lj" "A" "N" ]
epsij[ "lj" "n" "A" ] = epsij[ "lj" "A" "n" ]
epsij[ "lj" "O" "A" ] = epsij[ "lj" "A" "O" ]
epsij[ "lj" "S" "A" ] = epsij[ "lj" "A" "S" ]
epsij[ "lj" "H" "A" ] = epsij[ "lj" "A" "H" ]
epsij[ "lj" "P" "A" ] = epsij[ "lj" "A" "P" ]
epsij[ "lj" "f" "A" ] = epsij[ "lj" "A" "f" ]
epsij[ "lj" "F" "A" ] = epsij[ "lj" "A" "F" ]
epsij[ "lj" "c" "A" ] = epsij[ "lj" "A" "c" ]
epsij[ "lj" "b" "A" ] = epsij[ "lj" "A" "b" ]
epsij[ "lj" "I" "A" ] = epsij[ "lj" "A" "I" ]
epsij[ "lj" "M" "A" ] = epsij[ "lj" "A" "M" ]
epsij[ "lj" "Z" "A" ] = epsij[ "lj" "A" "Z" ]
epsij[ "lj" "L" "A" ] = epsij[ "lj" "A" "L" ]
epsij[ "lj" "X" "A" ] = epsij[ "lj" "A" "X" ]


#
# Hydrogen-bonding parameters
#
Rij[ "hb" "N" "H" ] = 1.90
Rij[ "hb" "O" "H" ] = 1.90
Rij[ "hb" "S" "H" ] = 2.50
Rij[ "hb" "H" "N" ] = 1.90
Rij[ "hb" "H" "O" ] = 1.90
Rij[ "hb" "H" "S" ] = 2.50

#
# AutoDock Binding Free Energy Model 140n
# New well-depths multiplied by the hbond.difference
# coefficient, 0.0656
#

#
# Original enthalpic well-depths for hydrogen-bonding...
#
epsij[ "hb" "N" "H" ] = 5.0 * FE_hbond_coeff
epsij[ "hb" "O" "H" ] = 5.0 * FE_hbond_coeff
epsij[ "hb" "S" "H" ] = 1.0 * FE_hbond_coeff
epsij[ "hb" "H" "N" ] = 5.0 * FE_hbond_coeff
epsij[ "hb" "H" "O" ] = 5.0 * FE_hbond_coeff
epsij[ "hb" "H" "S" ] = 1.0 * FE_hbond_coeff
